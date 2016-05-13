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
void rsGenerateEpiApplySmoothnessToVolume(const rsGenerateEpiParameters *p, const short t);
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
        (const char *) p->smoothnessReferencePath,
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
        fprintf(stdout, "Smoothness reference file: %s\n", p->smoothnessReferencePath);
        fprintf(stdout, "Timecourse reference file: %s\n", p->referencePath);
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

    // load smoothness reference file
    p->smoothnessReferenceFile = rsOpenNiftiFile(p->smoothnessReferencePath, RSNIFTI_OPEN_NONE);

    if (!p->smoothnessReferenceFile->readable) {
        fprintf(stderr, "\nError: The nifti file containing the smoothness reference (%s) could not be read.\n", p->smoothnessReferencePath);
        return;
    }

    if (p->smoothnessReferenceFile->xDim != p->correlationFile->xDim || p->smoothnessReferenceFile->yDim != p->correlationFile->yDim || p->smoothnessReferenceFile->zDim != p->correlationFile->zDim) {
        fprintf(stderr, "\nError: The dimensions of both the smoothness reference and correlation file must match!\n");
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

    if (p->verbose) {
        fprintf(stdout, "Smoothness Reference Dim: %d %d %d (%d Volumes)\n", p->smoothnessReferenceFile->xDim, p->smoothnessReferenceFile->yDim, p->smoothnessReferenceFile->zDim, p->smoothnessReferenceFile->vDim);
        fprintf(stdout, "Output Dim: %d %d %d (%d Volumes)\n", p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim);
    }

    rsCloseNiftiFile(p->smoothnessReferenceFile, FALSE);

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
        fprintf(stdout, "Estimating smoothness of the supplied smoothness reference epi\n");
    }

    rsGenerateEpiEstimateSmoothness(p);

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
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(t,tNum) shared(rng,p)
    {
        tNum = omp_get_thread_num();
        #pragma omp for schedule(guided, 1)
        for (t=0; t < p->output->vDim; t++) {
            rsGenerateEpiSampleVolume(p, rng[tNum], t);
        }
    }


    //free random number generator
    for (short i = 0; i < numThreads; i++) {
        gsl_rng_free(rng[i]);
    }
    free(rng);

    // smooth data
    if (p->verbose) {
        fprintf(stdout, "Apply spatial smoothness to the randomly generated data\n");
    }

    #pragma omp parallel num_threads(rsGetThreadsNum()) private(t) shared(p)
    {
        #pragma omp for schedule(guided, 1)
        for (t=0; t < p->output->vDim; t++) {
            rsGenerateEpiApplySmoothnessToVolume(p, t);
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
        fprintf(stdout, "Align temporal correlation with given the correlation map\n");
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

void rsGenerateEpiEstimateSmoothness(rsGenerateEpiParameters *p)
{

    p->smoothnessReferenceFile = rsOpenNiftiFile(p->smoothnessReferencePath, RSNIFTI_OPEN_READ);

    if (!p->smoothnessReferenceFile->readable) {
        fprintf(stderr, "\nError: The nifti file containing the smoothness reference (%s) could not be read.\n", p->smoothnessReferencePath);
        return;
    }

    double ****neighbourhoodCorrelations = d4lmatrix(p->nPoints-1, p->kernelSize - 1, p->kernelSize - 1, p->kernelSize - 1);

    // for every point in the mask compute the correlation with its surrounding points
    {
        unsigned long i;
        Point3D *pointA;
        Point3D *pointB;
        const int halfKernelSize = floor(p->kernelSize / 2);
        short j, k, l, x, y, z;
        const short maxX = p->smoothnessReferenceFile->xDim - 1;
        const short maxY = p->smoothnessReferenceFile->yDim - 1;
        const short maxZ = p->smoothnessReferenceFile->zDim - 1;
        double *seriesA;
        double *seriesB;
        double r;
        #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,j,k,l,x,y,z,pointA,pointB,seriesA,seriesB,r) shared(p,neighbourhoodCorrelations)
        {
            #pragma omp for schedule(guided, 1)
            for (i = 0; i < p->nPoints; i=i+1) {
                pointA = &p->maskPoints[i];
                seriesA = (double *) rsMalloc(sizeof(double) * p->smoothnessReferenceFile->vDim);
                seriesB = (double *) rsMalloc(sizeof(double) * p->smoothnessReferenceFile->vDim);
                rsExtractTimecourseFromRSNiftiFileBuffer(p->smoothnessReferenceFile, &seriesA[0], pointA);
                for (j = 0; j < p->kernelSize; j++) {
                    z = pointA->z - halfKernelSize + j;
                    for (k = 0; k < p->kernelSize; k++) {
                        y = pointA->y - halfKernelSize + k;
                        for (l = 0; l < p->kernelSize; l++) {
                            x = pointA->x - halfKernelSize + l;
                            if (x < 0 || x > maxX || y < 0 || y > maxY || z < 0 || z > maxZ) {
                                r = log(-1.0);
                            } else {
                                pointB = rsMakePoint3D(x, y, z);
                                rsExtractTimecourseFromRSNiftiFileBuffer(p->smoothnessReferenceFile, seriesB, pointB);
                                r = rsCorrelation(&seriesA[0], &seriesB[0], p->smoothnessReferenceFile->vDim);
                                rsFree(pointB);
                            }
                            neighbourhoodCorrelations[i][j][k][l] = r;
                        }
                    }

                }
                rsFree(seriesA);
                rsFree(seriesB);
          }
        }
    }

    // mean all correlation maps to build a mean correlation matrix which we later can use as a filter
    p->meanKernel = d3matrix(p->kernelSize - 1, p->kernelSize - 1, p->kernelSize - 1);
    for (short z = 0; z < p->kernelSize; z++) {
        for (short y = 0; y < p->kernelSize; y++) {
            for (short x = 0; x < p->kernelSize; x++) {
                double sum = 0.0;
                unsigned long points = 0;
                for (unsigned long i = 0; i < p->nPoints; i++) {
                    const double r = neighbourhoodCorrelations[i][z][y][x];
                    if (r != r) {
                        continue;
                    }
                    points += 1;
                    sum += r;
                }
                p->meanKernel[z][y][x] = sum / points;
            }
        }
    }


    free(neighbourhoodCorrelations[0][0][0]);
    free(neighbourhoodCorrelations[0][0]);
    free(neighbourhoodCorrelations[0]);
    free(neighbourhoodCorrelations);

    rsCloseNiftiFileAndFree(p->smoothnessReferenceFile);
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
}

void rsGenerateEpiApplySmoothnessToVolume(const rsGenerateEpiParameters *p, const short t)
{
    double ***data = d3matrix(p->output->zDim - 1, p->output->yDim - 1, p->output->xDim - 1);
    rsExtractVolumeFromRSNiftiFileBuffer(p->output, data[0][0], t);

    double ***xKernel = d3matrix(0,0,p->kernelSize-1);
    double ***yKernel = d3matrix(0,p->kernelSize-1,0);
    double ***zKernel = d3matrix(p->kernelSize-1,0,0);

    const short xMidKernel = (p->kernelSize-1)/2;
    const short yMidKernel = (p->kernelSize-1)/2;
    const short zMidKernel = (p->kernelSize-1)/2;

    for (short n=0; n<p->kernelSize; n++) {
        xKernel[0][0][n] = p->meanKernel[zMidKernel][yMidKernel][n];
    }

    for (short n=0; n<p->kernelSize; n++) {
        yKernel[0][n][0] = p->meanKernel[zMidKernel][n][xMidKernel];
    }

    for (short n=0; n<p->kernelSize; n++) {
        zKernel[n][0][0] = p->meanKernel[n][yMidKernel][xMidKernel];
    }

    double ***tmp = d3matrix(p->output->zDim - 1, p->output->yDim - 1, p->output->xDim - 1);

    rsConvolveWithKernel(tmp, data, xKernel, p->output->xDim, p->output->yDim, p->output->zDim, p->kernelSize, 1, 1);
    rsConvolveWithKernel(data, tmp, yKernel, p->output->xDim, p->output->yDim, p->output->zDim, 1, p->kernelSize, 1);
    rsConvolveWithKernel(tmp, data, zKernel, p->output->xDim, p->output->yDim, p->output->zDim, 1, 1, p->kernelSize);

    // write back to the buffer
    rsWriteVolumeToRSNiftiFileBuffer(p->output, tmp[0][0], t);

    // free up memory
    free(data[0][0]); free(data[0]); free(data);
    free(tmp[0][0]); free(tmp[0]); free(tmp);
    free(xKernel[0][0]); free(xKernel[0]); free(xKernel);
    free(yKernel[0][0]); free(yKernel[0]); free(yKernel);
    free(zKernel[0][0]); free(zKernel[0]); free(zKernel);
}

void rsGenerateEpiSmoothVoxelTemporally(const rsGenerateEpiParameters *p, const Point3D *point)
{
    const double sigma = 1.5;

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
    double ***kernel = rsCreateGaussianKernel(sigma, &kerneldim[0], &kerneldim[1], &kerneldim[2], 1, 1, 1);
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

/**
 * Copy of fslio's d4matrix() which uses unsigned long instead of int
 */
double ****d4lmatrix(unsigned long th, unsigned long  zh,  unsigned long  yh, unsigned long  xh)
{

    unsigned long j;
    unsigned long nvol = th+1;
    unsigned long nslice = zh+1;
    unsigned long nrow = yh+1;
    unsigned long ncol = xh+1;
    double ****t;


    /** allocate pointers to vols */
    t=(double ****) malloc((size_t)((nvol)*sizeof(double***)));
    if (!t) RSIOERR("d4lmatrix: allocation failure");

    /** allocate pointers to slices */
    t[0]=(double ***) malloc((size_t)((nvol*nslice)*sizeof(double**)));
    if (!t[0]) RSIOERR("d4lmatrix: allocation failure");

    /** allocate pointers for ydim */
    t[0][0]=(double **) malloc((size_t)((nvol*nslice*nrow)*sizeof(double*)));
    if (!t[0][0]) RSIOERR("d4lmatrix: allocation failure");


    /** allocate the data blob */
    t[0][0][0]=(double *) malloc((size_t)((nvol*nslice*nrow*ncol)*sizeof(double)));
    if (!t[0][0][0]) RSIOERR("d4lmatrix: allocation failure");


    /** point everything to the data blob */
    for(j=1;j<nrow*nslice*nvol;j++) t[0][0][j]=t[0][0][j-1]+ncol;
    for(j=1;j<nslice*nvol;j++) t[0][j]=t[0][j-1]+nrow;
    for(j=1;j<nvol;j++) t[j]=t[j-1]+nslice;

    return t;
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

    if (p->meanKernel != NULL) {
        free(p->meanKernel[0][0]);
        free(p->meanKernel[0]);
        free(p->meanKernel);
    }

    rsGenerateEpiFreeParams(p);
}
