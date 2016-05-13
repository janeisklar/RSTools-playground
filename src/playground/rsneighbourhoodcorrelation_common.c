#include <nifti/rsniftiutils.h>
#include "rsneighbourhoodcorrelation_common.h"
#include "rsneighbourhoodcorrelation_ui.h"
#include "utils/rsio.h"
#include <rscommon.h>

BOOL rsNeighbourhoodCorrelationComputeCorrelationForFile(rsNeighbourhoodCorrelationParameters *p, size_t fileIndex);

void rsNeighbourhoodCorrelationInit(rsNeighbourhoodCorrelationParameters *p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char **) p->inputPaths);

    BOOL outputsWritable = rsCheckOutputs((const char *[]) {
        (const char *) p->outputPath,
        RSIO_LASTFILE
    });

    if (!inputsReadable || !outputsWritable) {
        return;
    }

    rsSetThreadsNum(p->threads);

    p->dimLength = (p->distance*2)+1;
    int nNeighbourhood = p->dimLength * p->dimLength * p->dimLength;
    if (p->verbose) {
        fprintf(stdout, "Distance: %d\n", p->distance);
        fprintf(stdout, "Matrix will be of size: %dx%dx%d\n", p->dimLength, p->dimLength, p->dimLength);
        fprintf(stdout, "Output file: %s\n", p->outputPath);
    }

    // load first input as a reference
    p->reference = rsOpenNiftiFile(p->inputPaths[0], RSNIFTI_OPEN_NONE);

    if (!p->reference->readable) {
        fprintf(stderr, "\nError: The first nifti file that was supplied through stdin (%s) could not be read.\n", p->inputPaths[0]);
        return;
    }

    // Create output volume
    p->output = rsCloneNiftiFile(p->outputPath, p->reference, RSNIFTI_OPEN_ALLOC, nNeighbourhood);

    if (!p->output->readable) {
        fprintf(stderr, "\nError: The nifti file that was supplied as output (%s) could not be created.\n", p->outputPath);
        return;
    }

    if (p->verbose) {
        fprintf(stdout, "Output Dim: %d %d %d (%d Volumes)\n", p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim);
    }

    p->meanCount = rsMalloc(rsGetBufferSize(p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim, p->output->dt));

    rsCloseNiftiFile(p->reference, FALSE);

    p->parametersValid = TRUE;
    return;
}

void rsNeighbourhoodCorrelationRun(rsNeighbourhoodCorrelationParameters *p)
{
    p->parametersValid = FALSE;

    // initialize neighbourhood correlation matrix and its count (to compute the mean) with zero
    {
        if (p->verbose) {
            fprintf(stdout, "Initialize mean correlation matrices with zero\n");
        }

        double *zeroVector = (double *) rsMalloc(sizeof(double) * ((size_t) p->output->vDim));
        for (unsigned int t = 0; t < p->output->vDim; t++) {
            zeroVector[t] = 0.0;
        }
        for (short z = 0; z < p->output->zDim; z++) {
            for (short y = 0; y < p->output->yDim; y++) {
                for (short x = 0; x < p->output->xDim; x++) {
                    Point3D *point = rsMakePoint3D(x, y, z);
                    rsWriteTimecourseToRSNiftiFileBuffer(p->output, zeroVector, point);
                    rsWriteTimecourseToBuffer(p->output->dt, zeroVector, p->meanCount, 1.0f, 0.0f, point, p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim);
                    rsFree(point);
                }
            }
        }
        rsFree(zeroVector);
    }

    // compute neighbourhood correlation matrix for every file
    if (p->verbose) {
        fprintf(stdout, "For every file, sum up the correlation matrices in every point\n");
    }
    for (size_t i = 0; i < p->nInputs; i++) {
        rsNeighbourhoodCorrelationComputeCorrelationForFile(p ,i);
        if (p->verbose) {
            fprintf(stdout, "File %d of %d processed.\n", (1 + (int)i), (int)p->nInputs);
        }
    }

    // compute the mean by dividing the sum (currently stored in mean) by the count
    {
        if (p->verbose) {
            fprintf(stdout, "Divide by the total amount of correlation matrices in every point to get the mean\n");
        }
        double *sumVector = (double *) rsMalloc(sizeof(double) * ((size_t) p->output->vDim));
        double *countVector = (double *) rsMalloc(sizeof(double) * ((size_t) p->output->vDim));
        for (short z = 0; z < p->output->zDim; z++) {
            for (short y = 0; y < p->output->yDim; y++) {
                for (short x = 0; x < p->output->xDim; x++) {
                    Point3D *point = rsMakePoint3D(x, y, z);
                    rsExtractTimecourseFromRSNiftiFileBuffer(p->output, sumVector, point);
                    rsExtractTimecourseFromBuffer(p->output->dt, countVector, p->meanCount, 1.0f, 0.0f, point, p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim);
                    for (unsigned int t = 0; t < p->output->vDim; t++) {
                        sumVector[t] /= countVector[t];
                    }
                    rsWriteTimecourseToRSNiftiFileBuffer(p->output, sumVector, point);
                    rsFree(point);
                }
            }
        }
        rsFree(sumVector);
        rsFree(countVector);
        rsFree(p->meanCount);
    }

    if (p->verbose) {
        fprintf(stdout, "Write out result to: %s\n", p->outputPath);
    }

    rsWriteNiftiHeader(p->output->fslio, p->callString);
    FslWriteVolumes(p->output->fslio, p->output->data, p->output->vDim);

    p->parametersValid = TRUE;
}

BOOL rsNeighbourhoodCorrelationComputeCorrelationForFile(rsNeighbourhoodCorrelationParameters *p, size_t fileIndex)
{
    rsNiftiFile *file = rsOpenNiftiFile(p->inputPaths[fileIndex], RSNIFTI_OPEN_READ);

    if (!file->readable) {
        fprintf(stderr, "\nError: The nifti file %s could not be read.\n", p->inputPaths[fileIndex]);
        return FALSE;
    }

    if (file->xDim != p->output->xDim || file->yDim != p->output->yDim || file->zDim != p->output->zDim) {
        fprintf(stderr, "\nError: The nifti file %s was expected to have dimensions %dx%dx%d, but has %dx%dx%d.\n", p->output->xDim, p->output->yDim, p->output->zDim, file->xDim, file->yDim, file->zDim);
        return FALSE;
    }

    // for every point in the mask compute the correlation with its surrounding points
    Point3D *pointA;
    Point3D *pointB;
    short j, k, l, x, y, z, a, b, c;
    const short maxX = file->xDim - 1;
    const short maxY = file->yDim - 1;
    const short maxZ = file->zDim - 1;
    double *seriesA;
    double *seriesB;
    double ***correlationMatrix;
    double ***countMatrix;
    double r;

    #pragma omp parallel num_threads(rsGetThreadsNum()) private(j,k,l,x,y,z,a,b,c,pointA,pointB,seriesA,seriesB,r,correlationMatrix,countMatrix) shared(p,file)
    {
        // iterate over the whole volume
        #pragma omp for schedule(guided, 1)
        for (a = 0; a < file->zDim; a++) {

            seriesA = (double *) rsMalloc(sizeof(double) * file->vDim);
            seriesB = (double *) rsMalloc(sizeof(double) * file->vDim);
            correlationMatrix = d3matrix(p->dimLength - 1, p->dimLength - 1, p->dimLength - 1);
            countMatrix = d3matrix(p->dimLength - 1, p->dimLength - 1, p->dimLength - 1);

            for (b = 0; b < file->yDim; b++) {

                for (c=0; c < file->xDim; c++) {

                    // extract timecourse for current point
                    pointA = rsMakePoint3D(c, b, a);
                    rsExtractTimecourseFromRSNiftiFileBuffer(file, &seriesA[0], pointA);

                    // extract existing correlation matrix for this point
                    rsExtractTimecourseFromRSNiftiFileBuffer(p->output, &correlationMatrix[0][0][0], pointA);

                    // extract existing count matrix for this point
                    rsExtractTimecourseFromBuffer(p->output->dt, &countMatrix[0][0][0], p->meanCount, 1.0f, 0.0f, pointA, p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim);

                    // iterate over all neighbouring points and compute their correlation
                    for (j = 0; j < p->dimLength; j++) {
                        z = pointA->z - p->distance + j;

                        for (k = 0; k < p->dimLength; k++) {
                            y = pointA->y - p->distance + k;

                            for (l = 0; l < p->dimLength; l++) {
                                x = pointA->x - p->distance + l;

                                if (x < 0 || x > maxX || y < 0 || y > maxY || z < 0 || z > maxZ) {
                                    r = log(-1.0); // NaN
                                } else {
                                    pointB = rsMakePoint3D(x, y, z);
                                    rsExtractTimecourseFromRSNiftiFileBuffer(file, seriesB, pointB);
                                    r = rsCorrelation(&seriesA[0], &seriesB[0], file->vDim);
                                    // fix numerical floating point instabilities that might give us correlation coefficients of 1+eps or -1-eps
                                    if (r > 1.0) {
                                        r = 1.0;
                                    } else if (r < -1.0) {
                                        r = -1.0;
                                    }
                                    rsFree(pointB);
                                }

                                // if we're not dealing with NaN values add it to the mean
                                if (!isnan(r) && !isinf(r)) {
                                    countMatrix[j][k][l] += 1.0;
                                    correlationMatrix[j][k][l] += r;
                                }
                            }
                        }
                    }

                    // write back the count matrix and the correlation matrix
                    rsWriteTimecourseToRSNiftiFileBuffer(p->output, &correlationMatrix[0][0][0], pointA);
                    rsWriteTimecourseToBuffer(p->output->dt, &countMatrix[0][0][0], p->meanCount, 1.0f, 0.0f, pointA, p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim);

                    rsFree(pointA);
                }
            }

            rsFree(seriesA);
            rsFree(seriesB);

            free(correlationMatrix[0][0]); free(correlationMatrix[0]); free(correlationMatrix);
            free(countMatrix[0][0]); free(countMatrix[0]); free(countMatrix);
        }
    }

    rsCloseNiftiFileAndFree(file);
    return TRUE;
}

void rsNeighbourhoodCorrelationDestroy(rsNeighbourhoodCorrelationParameters *p)
{
    if (p->meanCount != NULL) {
        rsFree(p->meanCount);
    }

    if (p->output != NULL) {
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }

    rsNeighbourhoodCorrelationFreeParams(p);
}
