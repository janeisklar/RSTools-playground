#include <nifti/rsniftiutils.h>
#include "rsgenerateepi_common.h"
#include "rsgenerateepi_ui.h"
#include "utils/rsio.h"
#include <gsl/gsl_sf_exp.h>
#include <sys/time.h>
#include <gsl/gsl_vector_double.h>
#include <rscommon.h>

void rsGenerateEpiSampleVolume(const rsGenerateEpiParameters *p, const gsl_rng *rng, const short t);
void rsGenerateEpiAdjustCorrelationForPoint(const rsGenerateEpiParameters *p, const Point3D *point);
void rsGenerateEpiEstimateSmoothness(rsGenerateEpiParameters *p);
void rsGenerateEpiApplyNeighbourhoodCorrelationToVolume(const rsGenerateEpiParameters *p, const short t);
void rsGenerateEpiSmoothVoxelTemporally(const rsGenerateEpiParameters *p, const Point3D *point);
void rsConvolveWithKernel(double ***result, double ***input, double ***kernel, short xdim, short ydim, short zdim, short xdimKernel, short ydimKernel, short zdimKernel);
double ***rsCreateGaussianKernel(double sigma, short *xdimkernel, short *ydimkernel, short *zdimkernel, double xvoxsize, double yvoxsize, double zvoxsize);
double ****d4lmatrix(unsigned long th, unsigned long  zh,  unsigned long  yh, unsigned long  xh);

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
        (const char *) p->neighbourhoodCorrelationPath,
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
        fprintf(stdout, "Neighbourhood correlation file: %s\n", p->neighbourhoodCorrelationPath);
        fprintf(stdout, "Timecourse reference file: %s\n", p->referencePath);
        fprintf(stdout, "Correlation file: %s\n", p->correlationPath);
        fprintf(stdout, "Mean file: %s\n", p->meanPath);
        fprintf(stdout, "Std file: %s\n", p->stdPath);
        fprintf(stdout, "Mask file: %s\n", p->maskPath);
        fprintf(stdout, "Output file: %s\n", p->outputPath);
        fprintf(stdout, "TR: %.2f\n", p->TR);
        fprintf(stdout, "Temporal smoothing sigma: %.2f\n", p->sigma);
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

    // load neighbourhood correlation file
    if (p->verbose) {
        fprintf(stdout, "Reading in neighbourhood correlation file..\n");
    }
    p->neighbourhoodCorrelationFile = rsOpenNiftiFile(p->neighbourhoodCorrelationPath, RSNIFTI_OPEN_READ);

    if (!p->neighbourhoodCorrelationFile->readable) {
        fprintf(stderr, "\nError: The nifti file containing the neighbourhood correlation (%s) could not be read.\n", p->neighbourhoodCorrelationPath);
        return;
    }

    if (p->neighbourhoodCorrelationFile->xDim != p->correlationFile->xDim || p->neighbourhoodCorrelationFile->yDim != p->correlationFile->yDim || p->neighbourhoodCorrelationFile->zDim != p->correlationFile->zDim) {
        fprintf(stderr, "\nError: The dimensions x, y and z of both the neighbourhood correlation file and correlation file must match!\n");
        return;
    }

    // derive kernel size from the length of the given neighbourhood correlation file
    p->kernelSize = (int)round(pow(p->neighbourhoodCorrelationFile->vDim, 1.0/3.0));

    if (p->verbose) {
        fprintf(stdout, "Neighbourhood correlation file contains %d volumes. Assuming a correlation matrix width of %d.\n", p->neighbourhoodCorrelationFile->vDim, p->kernelSize);
    }

    if (p->kernelSize % 2 == 0) { // kernel size needs to be an uneven number
        fprintf(stderr, "\nError: The supplied neighbourhood correlation file does not look valid!\n");
        return;
    }

    // load timecourse reference file
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

    // get voxel spacing
    float xSpacing, ySpacing, zSpacing, tr;
    FslGetVoxDim(p->output->fslio, &xSpacing, &ySpacing, &zSpacing, &tr);
    FslSetVoxDim(p->output->fslio, xSpacing, ySpacing, zSpacing, p->TR);

    if (p->verbose) {
        fprintf(stdout, "Neighbourhood Correlation Dim: %d %d %d (%d Volumes)\n", p->neighbourhoodCorrelationFile->xDim, p->neighbourhoodCorrelationFile->yDim, p->neighbourhoodCorrelationFile->zDim, p->neighbourhoodCorrelationFile->vDim);
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
        fprintf(stdout, "Sampling random timecourses\n");
    }

    //initialize gsl random number generator
    {
        unsigned int seed;
        const unsigned int numThreads = rsGetThreadsNum();
        gsl_rng **rng;
        gsl_rng_env_setup();
        rng = (gsl_rng **) rsMalloc(((size_t) numThreads) * sizeof(gsl_rng *));
        const gsl_rng_type *rng_t = gsl_rng_default;
        struct timeval tv;

        for (unsigned int t = 0; t < numThreads; t++) {
            gettimeofday(&tv, 0);
            seed = tv.tv_sec + tv.tv_usec;
            rng[t] = gsl_rng_alloc(rng_t);
            gsl_rng_set(rng[t], seed * t);
        }

        unsigned int tNum, t;
        #pragma omp parallel num_threads(rsGetThreadsNum()) private(t,tNum) shared(rng,p)
        {
            tNum = omp_get_thread_num();
            #pragma omp for schedule(guided, 1)
            for (t = 0; t < p->output->vDim; t++) {
                rsGenerateEpiSampleVolume(p, rng[tNum], t);
            }
        }

        //free random number generator
        for (short i = 0; i < numThreads; i++) {
            gsl_rng_free(rng[i]);
        }
        free(rng);
    }

    // apply neighbourhood correlation
    if (p->verbose) {
        fprintf(stdout, "Re-assemble randomly generated data at every point using the corresponding neighbourhood correlation matrix to achieve an intrinsic structure similar to that of a real EPI\n");
    }

    unsigned int t;
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(t) shared(p)
    {
        #pragma omp for schedule(guided, 1)
        for (t=0; t < p->output->vDim; t++) {
            // apply correlation matrices three times to get more spatially stable results
            rsGenerateEpiApplyNeighbourhoodCorrelationToVolume(p, t);
            rsGenerateEpiApplyNeighbourhoodCorrelationToVolume(p, t);
            rsGenerateEpiApplyNeighbourhoodCorrelationToVolume(p, t);
        }
    }

    if (p->verbose) {
        fprintf(stdout, "Temporal smoothing\n");
    }

    unsigned long i;
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i) shared(p)
    {
        #pragma omp for schedule(guided, 1)
        for (i = 0; i < p->nPoints; i++) {
            rsGenerateEpiSmoothVoxelTemporally(p, &(p->maskPoints[i]));
        }
    }

    // align temporal correlation
    if (p->verbose) {
        fprintf(stdout, "Align temporal correlation with the given correlation map\n");
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

void rsGenerateEpiSampleVolume(const rsGenerateEpiParameters *p, const gsl_rng *rng, const short t)
{

    // initialize volume with NaN
    double ***data = d3matrix(p->output->zDim - 1, p->output->yDim - 1, p->output->xDim - 1);
    for (short z=0; z<p->output->zDim; z++) {
        for (short y=0; y<p->output->yDim; y++) {
            for (short x=0; x<p->output->xDim; x++) {
                data[z][y][x] = log(-1.0);
            }
        }
    }

    // draw random samples that follow N(0,1)
    for (unsigned long i = 0; i < p->nPoints; i++) {
        const Point3D *point = &p->maskPoints[i];
        data[point->z][point->y][point->x] = gsl_rng_uniform(rng);
    }

    // write volume to file buffer
    rsWriteVolumeToRSNiftiFileBuffer(p->output, data[0][0], t);

    // cleanup
    free(data[0][0]); free(data[0]); free(data);
}

void rsGenerateEpiAdjustCorrelationForPoint(const rsGenerateEpiParameters *p, const Point3D *point)
{
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
    rsFree(series);
}

void rsGenerateEpiApplyNeighbourhoodCorrelationToVolume(const rsGenerateEpiParameters *p, const short t)
{
    double ***data = d3matrix(p->output->zDim - 1, p->output->yDim - 1, p->output->xDim - 1);
    rsExtractVolumeFromRSNiftiFileBuffer(p->output, data[0][0], t);
    double ***tmp = d3matrix(p->output->zDim - 1, p->output->yDim - 1, p->output->xDim - 1);
    double ***kernel = d3matrix(p->kernelSize - 1, p->kernelSize - 1, p->kernelSize - 1);

    const short maxX = p->output->xDim - 1;
    const short maxY = p->output->yDim - 1;
    const short maxZ = p->output->zDim - 1;

    const short halfKernelSize = (p->kernelSize - 1) / 2;

    for (short z=0; z < p->output->zDim; z++) {
        for (short y=0; y < p->output->yDim; y++) {
            for (short x=0; x < p->output->xDim; x++) {
                Point3D *point = rsMakePoint3D(x, y, z);

                // extract correlation matrix for current point
                rsExtractTimecourseFromRSNiftiFileBuffer(p->neighbourhoodCorrelationFile, &kernel[0][0][0], point);

                // normalize correlation matrix to be able to use it as kernel
                double sum = 0;
                unsigned int count = 0;
                for (short i = 0; i < p->kernelSize-1; i++) {
                    for (short j = 0; j < p->kernelSize-1; j++) {
                        for (short k = 0; k < p->kernelSize-1; k++) {
                            if (isnan(kernel[i][j][k]) || isinf(kernel[i][j][k])) {
                                continue; // don't consider points for which the kernel does not have any value
                            }
                            short absX = x - halfKernelSize + k;
                            short absY = y - halfKernelSize + j;
                            short absZ = z - halfKernelSize + i;
                            if (absX < 0 || absX > maxX || absY < 0 || absY > maxY || absZ < 0 || absZ > maxZ) {
                                continue; // don't consider points that lie outside of the volume
                            }
                            sum += kernel[i][j][k];
                            count++;
                        }
                    }
                }
                const double norm = sum / count;

                // compute new value for the current data point based on the correlation kernel
                if (isnan(norm) || isinf(norm)) {
                    tmp[z][y][x] = data[z][y][x]; // simply keep the value if we don't have a proper kernel for this point
                } else {
                    double newValue = 0.0;
                    for (short i = 0; i < p->kernelSize-1; i++) {
                        for (short j = 0; j < p->kernelSize - 1; j++) {
                            for (short k = 0; k < p->kernelSize - 1; k++) {
                                if (isnan(kernel[i][j][k]) || isinf(kernel[i][j][k])) {
                                    continue; // don't consider points for which the kernel does not have any value
                                }
                                short absX = x - halfKernelSize + k;
                                short absY = y - halfKernelSize + j;
                                short absZ = z - halfKernelSize + i;
                                if (absX < 0 || absX > maxX || absY < 0 || absY > maxY || absZ < 0 || absZ > maxZ) {
                                    continue; // don't consider points that lie outside of the volume
                                }
                                newValue += data[absZ][absY][absX] * kernel[i][j][k];
                            }
                        }
                    }
                    tmp[z][y][x] = newValue / norm;
                }

                rsFree(point);
            }
        }
    }

    // write back to the buffer
    rsWriteVolumeToRSNiftiFileBuffer(p->output, tmp[0][0], t);

    // free up memory
    free(data[0][0]); free(data[0]); free(data);
    free(tmp[0][0]); free(tmp[0]); free(tmp);
    free(kernel[0][0]); free(kernel[0]); free(kernel);
}

void rsGenerateEpiSmoothVoxelTemporally(const rsGenerateEpiParameters *p, const Point3D *point)
{
    // get timeseries for point
    double *series = (double*)rsMalloc(sizeof(double) * ((size_t)p->output->vDim));
    double *smoothSeries = (double*)rsMalloc(sizeof(double) * ((size_t)p->output->vDim));
    rsExtractTimecourseFromRSNiftiFileBuffer(p->output, series, point);

    // fake 3d volume for convolution
    double **pSeries = (double**)rsMalloc(sizeof(double*));
    double ***ppSeries = (double***)rsMalloc(sizeof(double**));
    pSeries[0] = series;
    ppSeries[0] = pSeries;
    double **pSmoothSeries = (double**)rsMalloc(sizeof(double*));
    double ***ppSmoothSeries = (double***)rsMalloc(sizeof(double**));
    pSmoothSeries[0] = smoothSeries;
    ppSmoothSeries[0] = pSmoothSeries;

    // create gaussian kernel
    short kerneldim[3];
    double ***kernel = rsCreateGaussianKernel(p->sigma, &kerneldim[0], &kerneldim[1], &kerneldim[2], p->TR, 1, 1);
    double ***xKernel = d3matrix(0,0,kerneldim[0]-1);
    const short yMidKernel = (kerneldim[1]-1)/2;
    const short zMidKernel = (kerneldim[2]-1)/2;
    for (short n=0; n<kerneldim[0]; n++) {
        xKernel[0][0][n] = kernel[zMidKernel][yMidKernel][n];
    }

    // convolve
    rsConvolveWithKernel(ppSmoothSeries, ppSeries, xKernel, p->output->vDim, 1, 1, kerneldim[0], 1, 1);

    // write back corrected timeseries
    rsWriteTimecourseToRSNiftiFileBuffer(p->output, smoothSeries, point);

    // cleanup
    rsFree(ppSeries);
    rsFree(pSeries);
    rsFree(series);
    rsFree(ppSmoothSeries);
    rsFree(pSmoothSeries);
    rsFree(smoothSeries);
    free(kernel[0][0]); free(kernel[0]); free(kernel);
    free(xKernel[0][0]); free(xKernel[0]); free(xKernel);
}

void rsConvolveWithKernel(double ***result, double ***input, double ***kernel, short xdim, short ydim, short zdim, short xdimKernel, short ydimKernel, short zdimKernel) {

    const short midx=xdimKernel/2;
    const short midy=ydimKernel/2;
    const short midz=zdimKernel/2;

    for (short z=0; z<zdim; z++) {
        for (short y=0; y<ydim; y++) {
            for (short x=0; x<xdim; x++) {
                size_t nValues = 0;
                double val=0.0,
                    norm= 0.0;
                const short x3=x-midx,
                    y3=y-midy,
                    z3=z-midz;
                for (short mz=0; mz<zdimKernel; mz++) {
                    for (short my=0; my<ydimKernel; my++) {
                        for (short mx=0; mx<xdimKernel; mx++) {
                            const short x2=x3+mx,
                                y2=y3+my,
                                z2=z3+mz;
                            if (x2>=0 && x2<xdim && y2>=0 && y2<ydim && z2>=0 && z2<zdim) {
                                if (isnan(input[z2][y2][x2]) || isinf(input[z2][y2][x2])) {
                                    continue;
                                }
                                val+=input[z2][y2][x2] * kernel[mz][my][mx];
                                norm+=kernel[mz][my][mx];
                                nValues += 1;
                            }
                        }
                    }
                }

                if (nValues > 0) {
                    result[z][y][x] = fabs(norm) > 1e-12 ? val / norm : val;
                } else {
                    // fill with NaN if we excluded all values in the kernel
                    result[z][y][x] = log(-1.0);
                }
            }
        }
    }
}

double ***rsCreateGaussianKernel(double sigma, short *xdimkernel, short *ydimkernel, short *zdimkernel, double xvoxsize, double yvoxsize, double zvoxsize) {
    const double cutoff = 4.0,
        norm   = 2*sigma*sigma,
        dx2    = xvoxsize*xvoxsize,
        dy2    = yvoxsize*yvoxsize,
        dz2    = zvoxsize*zvoxsize;
    const short sx = ((short) ceil(sigma*cutoff/xvoxsize)),
        sy = ((short) ceil(sigma*cutoff/yvoxsize)),
        sz = ((short) ceil(sigma*cutoff/zvoxsize));
    *xdimkernel = 2*sx + 1;
    *ydimkernel = 2*sy + 1;
    *zdimkernel = 2*sz + 1;
    double ***kernel = d3matrix(*zdimkernel-1, *ydimkernel-1, *xdimkernel-1);

    for (short z=0; z<*zdimkernel; z++) {
        for (short y=0; y<*ydimkernel; y++) {
            for (short x=0; x<*xdimkernel; x++) {
                const short x2=x-sx,
                    y2=y-sy,
                    z2=z-sz;
                kernel[z][y][x]=exp(-(x2*x2*dx2+y2*y2*dy2+z2*z2*dz2)/norm);
            }
        }
    }

    return kernel;
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

    if (p->neighbourhoodCorrelationFile != NULL) {
        rsCloseNiftiFileAndFree(p->neighbourhoodCorrelationFile);
        p->neighbourhoodCorrelationFile = NULL;
    }

    if (p->output != NULL) {
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }

    rsGenerateEpiFreeParams(p);
}
