#include <nifti/rsniftiutils.h>
#include "rsgenerateepi_common.h"
#include "rsgenerateepi_ui.h"
#include "utils/rsio.h"
#include <gsl/gsl_sf_exp.h>
#include <sys/time.h>
#include <gsl/gsl_vector_double.h>

void rsGenerateEpiSampleVolume(const rsGenerateEpiParameters *p, const gsl_rng *rng, const gsl_vector *sqrtS, const gsl_matrix *U, const short t);
void rsGenerateEpiAdjustCorrelationForPoint(const rsGenerateEpiParameters *p, const Point3D *point);

void rsGenerateEpiInit(rsGenerateEpiParameters *p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char *[]) {
        (const char *) p->referencePath,
        (const char *) p->correlationPath,
        (const char *) p->meanPath,
        (const char *) p->stdPath,
        (const char *) p->maskPath,
        RSIO_LASTFILE
    });

    BOOL outputsWritable = rsCheckOutputs((const char *[]) {
        (const char *) p->outputPath,
        RSIO_LASTFILE
    });

    if (!inputsReadable || !outputsWritable) {
        return;
    }

    rsSetThreadsNum(p->threads);

    if (p->verbose) {
        fprintf(stdout, "Reference file: %s\n", p->referencePath);
        fprintf(stdout, "Correlation file: %s\n", p->correlationPath);
        fprintf(stdout, "Mean file: %s\n", p->meanPath);
        fprintf(stdout, "Std file: %s\n", p->stdPath);
        fprintf(stdout, "Mask file: %s\n", p->maskPath);
        fprintf(stdout, "Output file: %s\n", p->outputPath);
    }

    // load correlation file
    p->correlationFile = rsOpenNiftiFile(p->correlationPath, RSNIFTI_OPEN_READ);

    if (!p->correlationFile->readable) {
        fprintf(stderr, "\nError: The nifti file containing the correlation values (%s) could not be read.\n", p->correlationPath);
        return;
    }

    if (p->correlationFile->vDim > 1) {
        fprintf(stderr, "\nError: The supplied correlation map must not have more than one volume!\n");
        return;
    }

    // load mask file
    p->maskFile = rsOpenNiftiFile(p->maskPath, RSNIFTI_OPEN_READ);

    if (!p->maskFile->readable) {
        fprintf(stderr, "\nError: The nifti file containing the mask (%s) could not be read.\n", p->maskPath);
        return;
    }

    if (p->maskFile->xDim != p->correlationFile->xDim || p->maskFile->yDim != p->correlationFile->yDim || p->maskFile->zDim != p->correlationFile->zDim) {
        fprintf(stderr, "\nError: The dimensions of both the mask and correlation file must match!\n");
        return;
    }

    rsCloseNiftiFileAndFree(p->maskFile);

    // load reference file
    FILE *referenceFile = fopen(p->referencePath, "r");

    if (referenceFile == NULL) {
        fprintf(stderr, "\nError: The reference txt-file (%s) could not be read.\n", p->referencePath);
        return;
    }

    p->reference = rsReadRegressorFromStream(referenceFile, &p->nReferenceValues);
    fclose(referenceFile);

    if (p->nReferenceValues < 2) {
        fprintf(stderr, "\nError: The reference file (%s) must contain more than one value\n", p->referencePath);
        return;
    }

    if (p->verbose) {
        fprintf(stdout, "Correlation Dim: %d %d %d (%d Volumes)\n", p->correlationFile->xDim, p->correlationFile->yDim, p->correlationFile->zDim, p->correlationFile->vDim);
    }

    // load mean / std files
    if (p->meanPath != NULL) {
        p->meanFile = rsOpenNiftiFile(p->meanPath, RSNIFTI_OPEN_READ);
        if (!p->meanFile->readable) {
            fprintf(stderr, "\nError: The nifti file containing the mean (%s) could not be read.\n", p->meanPath);
            return;
        }
        if (p->meanFile->xDim != p->correlationFile->xDim || p->meanFile->yDim != p->correlationFile->yDim || p->meanFile->zDim != p->correlationFile->zDim) {
            fprintf(stderr, "\nError: The dimensions of both the mean and correlation file must match!\n");
            return;
        }
    }
    if (p->stdPath != NULL) {
        p->stdFile = rsOpenNiftiFile(p->stdPath, RSNIFTI_OPEN_READ);
        if (!p->stdFile->readable) {
            fprintf(stderr, "\nError: The nifti file containing the standard deviation (%s) could not be read.\n", p->stdPath);
            return;
        }
        if (p->stdFile->xDim != p->correlationFile->xDim || p->stdFile->yDim != p->correlationFile->yDim || p->stdFile->zDim != p->correlationFile->zDim) {
            fprintf(stderr, "\nError: The dimensions of both the standard deviation and correlation file must match!\n");
            return;
        }
    }

    // Create output volume
    p->output = rsCloneNiftiFile(p->outputPath, p->correlationFile, RSNIFTI_OPEN_ALLOC, p->nReferenceValues);

    if (!p->output->readable) {
        fprintf(stderr, "\nError: The nifti file that was supplied as output (%s) could not be created.\n", p->outputPath);
        return;
    }

    if (p->verbose) {
        fprintf(stdout, "Output Dim: %d %d %d (%d Volumes)\n", p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim);
    }

    p->parametersValid = TRUE;
    return;
}

void rsGenerateEpiRun(rsGenerateEpiParameters *p)
{
    p->parametersValid = FALSE;

    // de-mean reference timecourse
    double referenceValueMean = 0.0;
    for (unsigned int t=0; t<p->nReferenceValues; t++) {
        referenceValueMean += p->reference[t];
    }
    referenceValueMean /= ((double)p->nReferenceValues);
    for (unsigned int t=0; t<p->nReferenceValues; t++) {
        p->reference[t] -= referenceValueMean;
    }

    // get coordinate matrix
    nifti_image *outputNifti = p->output->fslio->niftiptr;
    mat44 outputWorldMatrix = outputNifti->sform_code == NIFTI_XFORM_UNKNOWN
                              ? outputNifti->qto_xyz
                              : outputNifti->sto_xyz;

    // extract correlation volume from the buffer
    double ***correlationData = d3matrix(p->correlationFile->zDim - 1, p->correlationFile->yDim - 1, p->correlationFile->xDim - 1);
    rsExtractVolumeFromRSNiftiFileBuffer(p->correlationFile, correlationData[0][0], 0);

    // extract list of points that are to be considered
    p->nPoints = 0L;
    p->maskPoints = rsReadMask(p->maskPath, p->correlationFile->xDim, p->correlationFile->yDim, p->correlationFile->zDim, &p->nPoints, NULL, p->correlationFile->fslio, NULL);

    if (p->maskPoints == NULL || p->nPoints < 1) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        return;
    }

    if (p->verbose) {
        fprintf(stdout, "Found %lu points in the mask\n", p->nPoints);
    }

    // create covariance matrix based on the squared exponential kernel (RBF)
    if (p->verbose) {
        fprintf(stdout, "\nCreating covariance matrix for all points based on the squared exponential kernel\n");
    }
    unsigned long i, j;
    const double alpha = 0.01;
    const double cutoff = 90; //gsl_sf_exp(15.0 * alpha); // set covariances between distances larger than 60 to 0.0 for faster computation
    gsl_matrix *cov = gsl_matrix_alloc(p->nPoints, p->nPoints);
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,j) shared(cov,outputWorldMatrix)
    {
        #pragma omp for schedule(guided, 1)
        for (i = 0; i < p->nPoints; i++) {
            gsl_vector *xy = gsl_vector_alloc(3);
            for (j = 0; j < p->nPoints; j++) {
                FloatPoint3D *pointI = rsMakeFloatPoint3D(0, 0, 0);
                FslGetMMCoord(outputWorldMatrix, (&p->maskPoints[i])->x, (&p->maskPoints[i])->y, (&p->maskPoints[i])->z, &pointI->x, &pointI->y, &pointI->z);

                FloatPoint3D *pointJ = rsMakeFloatPoint3D(0, 0, 0);
                FslGetMMCoord(outputWorldMatrix, (&p->maskPoints[j])->x, (&p->maskPoints[j])->y, (&p->maskPoints[j])->z, &pointJ->x, &pointJ->y, &pointJ->z);

                // compute kernel: k = exp(-alpha * ||x-y||^2)

                // compute vector x-y
                gsl_vector_set(xy, 0, pointI->x - pointJ->x);
                gsl_vector_set(xy, 1, pointI->y - pointJ->y);
                gsl_vector_set(xy, 2, pointI->z - pointJ->z);

                // compute norm of x-y
                const double norm = gsl_blas_dnrm2(xy);

                // compute kernel value
                const double k = gsl_sf_exp(-alpha * norm);

                // set corresponding value in covariance matrix
                //gsl_matrix_set(cov, i, j, norm > cutoff ? 0.0 : k);
                gsl_matrix_set(cov, i, j, k);

                // cleanup
                rsFree(pointI);
                rsFree(pointJ);
            }
            gsl_vector_free(xy);
        }
    }

    // factorize covariance matrix
    if (p->verbose) {
        fprintf(stdout, "Factorizing covariance matrix through SVD (%lu² elements, be patient)\n", p->nPoints);
    }

    // decompose: cov = U S V^T (U will be saved in cov!)
    gsl_vector *work = gsl_vector_alloc(p->nPoints);
    gsl_matrix *V = gsl_matrix_alloc(p->nPoints, p->nPoints);
    gsl_vector *S = gsl_vector_alloc(p->nPoints);
    gsl_linalg_SV_decomp(cov, V, S, work);
    gsl_vector *sqrtS = gsl_vector_alloc(p->nPoints);

    for (unsigned long i=0; i<p->nPoints; i++) {
        gsl_vector_set(sqrtS, i, sqrt(gsl_vector_get(S, i)));
    }

    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_matrix_free(V);

    if (p->verbose) {
        fprintf(stdout, "Sampling random timecourses\n");
    }

    //initialize gsl random number generator
    unsigned int seed;
    const unsigned int numThreads = rsGetThreadsNum();
    gsl_rng **rng;
    gsl_rng_env_setup();
    rng = (gsl_rng **)rsMalloc(((size_t)numThreads) * sizeof(gsl_rng *));
    const gsl_rng_type *rng_t = gsl_rng_default;
    struct timeval tv;

    for (unsigned int t = 0; t < numThreads; t++) {
        gettimeofday(&tv, 0);
        seed = tv.tv_sec + tv.tv_usec;
        rng[t] = gsl_rng_alloc(rng_t);
        gsl_rng_set(rng[t], seed * t);
    }

    unsigned int tNum, t;
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(t,tNum) shared(rng,p,sqrtS,cov)
    {
        tNum = omp_get_thread_num();
        #pragma omp for schedule(guided, 1)
        for (t=0; t < p->output->vDim; t++) {
            rsGenerateEpiSampleVolume(p, rng[tNum], sqrtS, cov, t);
        }
    }

    //free random number generator
    for (i = 0; i < numThreads; i++) {
        gsl_rng_free(rng[i]);
    }
    free(rng);

    gsl_matrix_free(cov);

    // align temporal correlation
    if (p->verbose) {
        fprintf(stdout, "Align temporal correlation with given correlation map\n");
    }

    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i) shared(p)
    {
        #pragma omp for schedule(guided, 1)
        for (i = 0; i < p->nPoints; i++) {
            rsGenerateEpiAdjustCorrelationForPoint(p, &(p->maskPoints[i]));
        }
    }

    // apply standard deviation
    if (p->stdFile != NULL) {
        if (p->verbose) {
            fprintf(stdout, "Scaling to the supplied standard deviation\n");
        }
        double ***stdData = d3matrix(p->stdFile->zDim - 1, p->stdFile->yDim - 1, p->stdFile->xDim - 1);
        double ***outputData = d3matrix(p->output->zDim - 1, p->output->yDim - 1, p->output->xDim - 1);
        rsExtractVolumeFromRSNiftiFileBuffer(p->stdFile, stdData[0][0], 0);
        for (short t=0; t<p->output->vDim; t++) {
            rsExtractVolumeFromRSNiftiFileBuffer(p->output, outputData[0][0], t);
            for (short z=0; z<p->output->zDim; z++) {
                for (short y=0; y<p->output->yDim; y++) {
                    for (short x=0; x<p->output->xDim; x++) {
                        outputData[z][y][x] *= stdData[z][y][x];
                    }
                }
            }
            rsWriteVolumeToRSNiftiFileBuffer(p->output, outputData[0][0], t);
        }
        free(stdData[0][0]); free(stdData[0]); free(stdData);
        free(outputData[0][0]); free(outputData[0]); free(outputData);
    }

    // apply mean
    if (p->meanFile != NULL) {
        if (p->verbose) {
            fprintf(stdout, "Adding the supplied mean\n");
        }
        double ***meanData = d3matrix(p->meanFile->zDim - 1, p->meanFile->yDim - 1, p->meanFile->xDim - 1);
        double ***outputData = d3matrix(p->output->zDim - 1, p->output->yDim - 1, p->output->xDim - 1);
        rsExtractVolumeFromRSNiftiFileBuffer(p->meanFile, meanData[0][0], 0);
        for (short t=0; t<p->output->vDim; t++) {
            rsExtractVolumeFromRSNiftiFileBuffer(p->output, outputData[0][0], t);
            for (short z=0; z<p->output->zDim; z++) {
                for (short y=0; y<p->output->yDim; y++) {
                    for (short x=0; x<p->output->xDim; x++) {
                        outputData[z][y][x] += meanData[z][y][x];
                    }
                }
            }
            rsWriteVolumeToRSNiftiFileBuffer(p->output, outputData[0][0], t);
        }
        free(meanData[0][0]); free(meanData[0]); free(meanData);
        free(outputData[0][0]); free(outputData[0]); free(outputData);
    }

    if (p->verbose) {
        fprintf(stdout, "Write out result to: %s\n", p->outputPath);
    }

    rsWriteNiftiHeader(p->output->fslio, p->callString);
    FslWriteVolumes(p->output->fslio, p->output->data, p->output->vDim);

    p->parametersValid = TRUE;
}

void rsGenerateEpiSampleVolume(const rsGenerateEpiParameters *p, const gsl_rng *rng, const gsl_vector *sqrtS, const gsl_matrix *U, const short t) {
    // initialize result with NaN
    double ***data = d3matrix(p->output->zDim - 1, p->output->yDim - 1, p->output->xDim - 1);
    for (short z=0; z<p->output->zDim; z++) {
        for (short y=0; y<p->output->yDim; y++) {
            for (short x=0; x<p->output->xDim; x++) {
                data[z][y][x] = 0.0;
            }
        }
    }

    // draw random samples that follow N(0,1)
    gsl_vector *u = gsl_vector_alloc(p->nPoints);
    for (unsigned long i = 0; i<p->nPoints; i++) {
        const double sample = gsl_rng_uniform(rng);
        gsl_vector_set(u, i, sample);
    }

    // turn them into samples of the gaussian process that follow N(0,cov)
    // see http://stats.stackexchange.com/questions/159313/generating-samples-from-singular-gaussian-distribution#answer-159322
    // z = U * sqrt(S) * u
    gsl_vector *z = gsl_vector_alloc(p->nPoints);
    gsl_vector *tmp = gsl_vector_alloc(p->nPoints);
    for (unsigned long i = 0; i<p->nPoints; i++) {
        const double Si = gsl_vector_get(sqrtS, i);
        const double ui = gsl_vector_get(u, i);
        gsl_vector_set(tmp, i, Si * ui); // <- tmp = sqrt(S) * u
    }
    gsl_blas_dgemv(CblasNoTrans, 1.0, U, tmp, 0.0, z); // <- z = U * tmp

    // write the newly generated samples to the volume
    for (unsigned long i = 0; i < p->nPoints; i++) {
        Point3D *point = &p->maskPoints[i];
        data[point->z][point->y][point->x] = gsl_vector_get(z, i);
    }

    // write volume to file buffer
    rsWriteVolumeToRSNiftiFileBuffer(p->output, data[0][0], t);

    // cleanup
    gsl_vector_free(u);
    gsl_vector_free(z);
    free(data[0][0]); free(data[0]); free(data);
}

void rsGenerateEpiAdjustCorrelationForPoint(const rsGenerateEpiParameters *p, const Point3D *point) {
    // get timeseries for point
    double *series = (double*)rsMalloc(sizeof(double) * ((size_t)p->output->vDim));
    rsExtractTimecourseFromRSNiftiFileBuffer(p->output, series, point);

    // get the correlation coefficient expectation for the given point
    double r = 0.0;
    rsExtractTimecourseFromRSNiftiFileBuffer(p->correlationFile, &r, point);
    double theta = acos(r); // corresponding angle

    // turn into gsl vectors
    gsl_vector *x1 = gsl_vector_alloc(p->output->vDim);
    gsl_vector *x2 = gsl_vector_alloc(p->output->vDim);
    gsl_matrix  *U = gsl_matrix_alloc(p->output->vDim, 1);
    double meanX2 = 0.0;
    for (int i=0; i<p->output->vDim; i++) {
        gsl_vector_set(x1, i, p->reference[i]);
        gsl_matrix_set(U, i, 0, p->reference[i]); // copy of x1 which is needed later on
        gsl_vector_set(x2, i, series[i]);
        meanX2 += series[i] / ((double)p->output->vDim);
    }

    // de-mean the series (should already roughly be the case, but just to make it more numerically stable)
    gsl_vector_add_constant(x2, meanX2 * -1.0);

    // fix the correlation of the series while leaving the reference timecourse untouched.
    // the following is based on:
    // http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable#15040

    // compute the orthogonal basis of x1 (first column of U in the SVD)
    gsl_vector *work = gsl_vector_alloc(1);
    gsl_matrix *V = gsl_matrix_alloc(1, 1);
    gsl_vector *S = gsl_vector_alloc(1);
    gsl_linalg_SV_decomp(U, V, S, work);
    gsl_vector_free(work);
    gsl_matrix_free(V);
    gsl_vector_free(S);

    // P = U * U'
    // projection onto space defined by x1
    gsl_matrix *P = gsl_matrix_alloc(p->output->vDim, p->output->vDim);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, U, U, 0.0, P);
    gsl_matrix_free(U);

    // make x2 orthogonal to x1
    // x2o = (Identity-P) * x2
    gsl_vector *x2o = gsl_vector_alloc(p->output->vDim);
    gsl_matrix *P2 = gsl_matrix_alloc(p->output->vDim, p->output->vDim);
    gsl_matrix_set_identity(P2);
    gsl_matrix_sub(P2, P);
    gsl_blas_dgemv(CblasNoTrans, 1.0, P2, x2, 0.0, x2o);
    gsl_matrix_free(P);
    gsl_matrix_free(P2);

    // X = [x1 x2o]
    gsl_matrix *X = gsl_matrix_alloc(p->output->vDim, 2);
    gsl_matrix_set_col(X, 0, x1);
    gsl_matrix_set_col(X, 1, x2o);
    gsl_vector_free(x1);
    gsl_vector_free(x2o);

    // scale columns to length 1
    // ssX = sum(X.^2, 1)
    // (column sum of squares)
    double ssX[2] = {0.0, 0.0};
    for (int i=0; i<p->output->vDim; i++) {
        ssX[0] += pow(gsl_matrix_get(X, i, 0), 2.0);
        ssX[1] += pow(gsl_matrix_get(X, i, 1), 2.0);
    }

    // diagInvSSX = diag(1 / sqrt(ssX))
    gsl_matrix *diagInvSSX = gsl_matrix_alloc(2,2);
    gsl_matrix_set_zero(diagInvSSX);
    gsl_matrix_set(diagInvSSX, 0, 0, 1 / sqrt(ssX[0]));
    gsl_matrix_set(diagInvSSX, 1, 1, 1 / sqrt(ssX[1]));

    // Y = X * diagInvSSX
    gsl_matrix *Y = gsl_matrix_alloc(p->output->vDim, 2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, diagInvSSX, 0.0, Y);
    gsl_matrix_free(X);
    gsl_matrix_free(diagInvSSX);

    // compute the adjusted vector x2
    // x2 = y2 + (1 / tan(theta)) * y1
    gsl_vector_view y1 = gsl_matrix_column(Y, 0);
    gsl_vector_view y2 = gsl_matrix_column(Y, 1);
    double invTanTheta = 1 / tan(theta);
    gsl_vector_scale(&y1.vector, invTanTheta);
    gsl_vector_memcpy(x2, &y2.vector);
    gsl_vector_add(x2, &y1.vector);

    // finally, set standard deviation to 1
    const double stdX2 = gsl_stats_sd(x2->data, 1, p->output->vDim);
    gsl_vector_scale(x2, 1.0 / stdX2);

    // write result to series c array
    for (int i=0; i<p->output->vDim; i++) {
        series[i] = gsl_vector_get(x2, i);
    }

    gsl_matrix_free(Y);
    gsl_vector_free(x2);

    // write back corrected timeseries
    rsWriteTimecourseToRSNiftiFileBuffer(p->output, series, point);
}

void rsGenerateEpiDestroy(rsGenerateEpiParameters *p)
{
    if (p->stdFile != NULL) {
        rsCloseNiftiFileAndFree(p->stdFile);
        p->stdFile = NULL;
    }

    if (p->meanFile != NULL) {
        rsCloseNiftiFileAndFree(p->meanFile);
        p->meanFile = NULL;
    }

    if (p->correlationFile != NULL) {
        rsCloseNiftiFileAndFree(p->correlationFile);
        p->correlationFile = NULL;
    }

    if (p->output != NULL) {
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }

    rsGenerateEpiFreeParams(p);
}
