#include <nifti/rsniftiutils.h>
#include "rssmoothing_common.h"
#include "rsmutualinformation_common.h"
#include "rsmutualinformation_ui.h"
#include "utils/rsio.h"
#include <gsl/gsl_sf_exp.h>
#include <sys/time.h>
#include <gsl/gsl_vector_double.h>
#include <rscommon.h>

void rsMutualInformationInit(rsMutualInformationParameters *p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char *[]) {
        (const char *) p->maskPath,
        RSIO_LASTFILE
    });

    inputsReadable = inputsReadable && rsCheckInputs((const char **)p->filePaths);

    if (!inputsReadable) {
        return;
    }

    rsSetThreadsNum(p->threads);

    if (p->verbose) {
        fprintf(stdout, "Mask file: %s\n", p->maskPath);
    }

    // load mask file
    p->maskFile = rsOpenNiftiFile(p->maskPath, RSNIFTI_OPEN_READ);

    if (!p->maskFile->readable) {
        fprintf(stderr, "\nError: The nifti file containing the mask (%s) could not be read.\n", p->maskPath);
        return;
    }

    // extract points from mask
    double ***mask = d3matrix(p->maskFile->zDim-1, p->maskFile->yDim-1, p->maskFile->xDim-1);
    rsExtractVolumeFromRSNiftiFileBuffer(p->maskFile, mask[0][0], 0);

    p->nPoints = 0L;
    unsigned int x=0,y=0,z=0;
    for (x=0; x<p->maskFile->xDim; x=x+1) {
        for (y=0; y<p->maskFile->yDim; y=y+1) {
            for (z=0; z<p->maskFile->zDim; z=z+1) {
                if (mask[z][y][x] > 0.01 ) {
                    p->nPoints += ((unsigned long)1L);
                }
            }
        }
    }

    p->maskPoints = (Point3D*)rsMalloc(p->nPoints*sizeof(Point3D));

    unsigned long i=0;
    for (x=0; x<p->maskFile->xDim; x=x+1) {
        for (y=0; y<p->maskFile->yDim; y=y+1) {
            for (z=0; z<p->maskFile->zDim; z=z+1) {
                if (mask[z][y][x] > 0.01 ) {
                    Point3D *point = rsMakePoint3D(x,y,z);
                    memcpy(&p->maskPoints[i], point, sizeof(Point3D));
                    rsFree(point);
                    i += 1;
                }
            }
        }
    }

    rsFree(mask[0][0]); rsFree(mask[0]); rsFree(mask);

    // load files
    if (p->verbose) {
        fprintf(stdout, "Loading %d input files\n", (int)p->nFiles);
    }

    p->files = (double**)rsMalloc(p->nFiles * sizeof(double*));

    for (unsigned int i=0; i<p->nFiles; i++) {
        if (p->verbose) {
            fprintf(stdout, "%s..\n", p->filePaths[i]);
        }

        rsNiftiFile *file = rsOpenNiftiFile(p->filePaths[i], RSNIFTI_OPEN_READ);

        if (p->maskFile->xDim != file->xDim || p->maskFile->yDim != file->yDim || p->maskFile->zDim != file->zDim) {
            fprintf(stderr, "\nError: The dimensions of both the mask and input file '%s' must match!\n", p->filePaths[i]);
            return;
        }

        p->files[i] = (double*)rsMalloc(p->nPoints * sizeof(double));
        rsExtractPointsFromRSNiftiFileBuffer(file, &(p->files[i][0]), p->maskPoints, p->nPoints, 0);
        rsCloseNiftiFileAndFree(file);
    }

    rsCloseNiftiFileAndFree(p->maskFile);

    p->parametersValid = TRUE;
    return;
}

void rsMutualInformationRun(rsMutualInformationParameters *p)
{
    p->parametersValid = FALSE;

    double mi[p->nFiles][p->nFiles];

    // Initialize diagonal of MI matrix
    for (unsigned int i=0; i<p->nFiles; i++) {
        mi[i][i] = 1.0;
    }

    // Compute minima and maxima of points in the mask for each file
    if (p->verbose) {
        fprintf(stdout, "Computing within-mask min and max intensity values for each file\n");
    }

    p->fileMinima = (double*)rsMalloc(p->nFiles * sizeof(double));
    p->fileMaxima = (double*)rsMalloc(p->nFiles * sizeof(double));

    unsigned int i;
    unsigned long j;
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,j) shared(p)
    {
        #pragma omp for schedule(guided)
        for (i=0; i<p->nFiles; i++) {
            // initialize with first point
            p->fileMinima[i] = p->files[i][0];
            p->fileMaxima[i] = p->fileMinima[i];

            // go through all points
            for (j=0; j<p->nPoints; j++) {
                const double currentValue = p->files[i][j];
                if (isnan(currentValue) || isinf(currentValue)) {
                    continue;
                }
                if (p->fileMinima[i] > currentValue) {
                    p->fileMinima[i] = currentValue;
                }
                if (p->fileMaxima[i] < currentValue) {
                    p->fileMaxima[i] = currentValue;
                }
            }
        }
    }

    if (p->verbose) {
        fprintf(stdout, "Computing mutual information coefficient\n");
    }

    // Iterate over all elements of the resulting matrix and compute the MI
    const double totalComputations = 0.5 * (p->nFiles * p->nFiles + p->nFiles);
    unsigned long currentComputation = 0L;
    for (i=0; i<p->nFiles; i++) {
        for (unsigned int j=0; j<=i; j++) {
            mi[i][j] = rsMutualInformationComputeMI(p, i, j);
            mi[j][i] = mi[i][j];
            currentComputation += 1L;
        }
        if (p->verbose) {
            fprintf(stdout, "% 4.2f%%\n", (100.0*(double)currentComputation)/totalComputations);
        }
    }

    if (p->verbose) {
        fprintf(stdout, "Result:\n");
    }

    // Write out MI matrix
    for (unsigned int i=0; i<p->nFiles; i++) {
        for (unsigned int j=0; j<p->nFiles; j++) {
            fprintf(stdout, "%.10f ", mi[i][j]);
        }
        fprintf(stdout, "\n");
    }

    p->parametersValid = TRUE;
}

double rsMutualInformationComputeMI(const rsMutualInformationParameters *params, const unsigned int indexA, const unsigned int indexB)
{
    // bin the number of points so that each bin can be executed in parallel
    const unsigned int nProcBins = rsGetThreadsNum() * 100;
    const unsigned long procBinWidth = (unsigned long)lround(((double)params->nPoints) / ((double)nProcBins));
    const unsigned int  nHistBins = 256;
    unsigned int b, i, j;
    unsigned long p;

    const double* pointsA = params->files[indexA];
    const double* pointsB = params->files[indexB];

    const double minIntensityA = params->fileMinima[indexA];
    const double minIntensityB = params->fileMinima[indexB];

    const double maxIntensityA = params->fileMaxima[indexA];
    const double maxIntensityB = params->fileMaxima[indexB];

    const double intensityRangeA = maxIntensityA - minIntensityA;
    const double intensityRangeB = maxIntensityB - minIntensityB;

    const double histWidthA = intensityRangeA / nHistBins;
    const double histWidthB = intensityRangeB / nHistBins;

    //unsigned int parallelBinHist[nHistBins][nHistBins][nProcBins];
    unsigned int ***parallelBinHist = (unsigned int***)rsMalloc(nHistBins * sizeof(unsigned int **));
    for (i = 0; i < nHistBins; i++) {
        parallelBinHist[i] = (unsigned int**)rsMalloc(nHistBins * sizeof(unsigned int *));
        for (j = 0; j < nHistBins; j++) {
            parallelBinHist[i][j] = (unsigned int*)rsMalloc(nProcBins * sizeof(unsigned int));
        }
    }

    #pragma omp parallel num_threads(rsGetThreadsNum()) private(b,i,j,p) shared(parallelBinHist,params)
    {
        #pragma omp for schedule(guided)
        for (b=0; b<nProcBins; b++) {
            // initialize histogram for this thread
            for (i=0; i<nHistBins; i++) {
                for (j=0; j<nHistBins; j++) {
                    parallelBinHist[i][j][b] = 0;
                }
            }

            // determine from which point to which point the histogram should be computed
            const unsigned long startPoint = b * procBinWidth;
            const unsigned long endPoint = MIN((b+1) * procBinWidth - 1, params->nPoints);

            // build up histogram for the point range
            for (p=startPoint; p<=endPoint; p++) {
                const double intensityA = pointsA[p];
                const double intensityB = pointsB[p];

                if (isnan(intensityA) || isnan(intensityB) || isinf(intensityA) ||isinf(intensityB)) {
                    continue;
                }

                const unsigned int binA = (unsigned int)floor((intensityA - minIntensityA)/histWidthA);
                const unsigned int binB = (unsigned int)floor((intensityB - minIntensityB)/histWidthB);
                parallelBinHist[MAX(0, MIN(binA, nHistBins-1))][MAX(0, MIN(binB, nHistBins-1))][b] += 1;
            }
        }
    }

    // combine histograms from all threads
    const double epsilon = 1e-16;
    double **binHist = d2matrix(nHistBins-1, nHistBins-1);
    for (i=0; i<nHistBins; i++) {
        for (j=0; j<nHistBins; j++) {
            binHist[i][j] = epsilon;
            for (b=0; b<nProcBins; b++) {
                binHist[i][j] += (double)parallelBinHist[i][j][b];
            }
        }
    }

    // cleanup parallel histograms
    for (i = 0; i < nHistBins; i++)
        for (j = 0; j < nHistBins; j++)
        rsFree(parallelBinHist[i][j]);
    for (i = 0; i < nHistBins; i++)
    rsFree(parallelBinHist[i]);
    rsFree(parallelBinHist);

    if (params->kernelFWHM >= 0.1) {
        // create smoothing kernel
        short kerneldim[3];
        double kernelFWHM = params->kernelFWHM; // in bins
        double kernelSigma = kernelFWHM / (2.0*sqrt(2.0*log(2.0)));
        double ***kernel3D = rsCreateGaussianKernel(kernelSigma, &kerneldim[0], &kerneldim[1], &kerneldim[2], 1, 1, 1);
        double ***xKernel = d3matrix(0,0,kerneldim[0]-1);
        double ***yKernel = d3matrix(0,kerneldim[1]-1,0);
        const short xMidKernel = (kerneldim[0]-1)/2,
                    yMidKernel = (kerneldim[1]-1)/2,
                    zMidKernel = (kerneldim[2]-1)/2;
        for (short n=0; n<kerneldim[0]; n++)
            xKernel[0][0][n] = kernel3D[zMidKernel][yMidKernel][n];
        for (short n=0; n<kerneldim[1]; n++)
            yKernel[0][n][0] = kernel3D[zMidKernel][n][xMidKernel];

        // convolve with histogram
        double **tmpBinHist = d2matrix(nHistBins-1, nHistBins-1);
        rsConvolveWithKernel(&tmpBinHist, &binHist, xKernel, nHistBins, nHistBins, 1, kerneldim[0], 1, 1);
        rsConvolveWithKernel(&binHist, &tmpBinHist, yKernel, nHistBins, nHistBins, 1, 1, kerneldim[0], 1);

        // cleanup smoothing kernel
        rsFree(kernel3D[0][0]); rsFree(kernel3D[0]); rsFree(kernel3D);
        rsFree(xKernel[0][0]);  rsFree(xKernel[0]);  rsFree(xKernel);
        rsFree(yKernel[0][0]);  rsFree(yKernel[0]);  rsFree(yKernel);
        rsFree(tmpBinHist[0]);  rsFree(tmpBinHist);
    }

    // normalize histogram
    double histNorm = 0.0;
    for (i=0; i<nHistBins; i++) {
        for (j=0; j<nHistBins; j++) {
            histNorm += binHist[i][j];
        }
    }
    for (i=0; i<nHistBins; i++) {
        for (j=0; j<nHistBins; j++) {
            binHist[i][j] /= histNorm;
        }
    }

    // compute p(A) and p(B)
    double pA[nHistBins];
    double pB[nHistBins];
    for (i=0; i<nHistBins; i++) {
        pA[i] = 0.0;
        pB[i] = 0.0;
        for (j=0; j<nHistBins; j++) {
            pA[i] += binHist[i][j];
            pB[i] += binHist[j][i];
        }
    }

    // compute H(A) and H(B)
    // where H(A) = - sum_i[p(i) * log2(p(i))]
    double hA = 0.0, hB = 0.0;
    for (i=0; i<nHistBins; i++) {
        hA -= pA[i] * log2(pA[i]);
        hB -= pB[i] * log2(pB[i]);
    }

    // compute H(A,B)
    // where H(A,B) = sum_i[sum_j[ p(i,j) * log2(p(i,j)) ]]
    double hAB = 0.0;
    for (i=0; i<nHistBins; i++) {
        for (j=0; j<nHistBins; j++) {
            hAB += binHist[i][j] * log2(binHist[i][j]);
        }
    }

    // cleanup
    rsFree(binHist[0]); rsFree(binHist);

    // return normalized mutual information
    // -> (H(A) + H(B)) / H(A,B)
    return (hA + hB) / hAB;
}

void rsMutualInformationDestroy(rsMutualInformationParameters *p)
{
    rsMutualInformationFreeParams(p);
}
