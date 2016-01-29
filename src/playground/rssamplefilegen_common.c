#include "rssamplefilegen_common.h"
#include "rssamplefilegen_ui.h"
#include "utils/rsio.h"

void rsSampleFileGenInit(rsSampleFileGenParameters *p)
{
    p->parametersValid = FALSE;
    
    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->referencepath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->outputpath,
        RSIO_LASTFILE
    });

    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }
    
    rsSetThreadsNum(p->threads);

    if ( p->verbose ) {
        fprintf(stdout, "Reference file: %s\n", p->referencepath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
    }
    
    p->reference = rsOpenNiftiFile(p->referencepath, RSNIFTI_OPEN_READ);

    if ( ! p->reference->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied for reference (%s) could not be read.\n", p->referencepath);
        return;
    }
       
    if ( p->verbose ) {
        fprintf(stdout, "Reference Dim: %d %d %d (%d Volumes)\n", p->reference->xDim, p->reference->yDim, p->reference->zDim, p->reference->vDim);
    }

    // Create output volume
    p->output = rsCloneNiftiFile(p->outputpath, p->reference, RSNIFTI_OPEN_ALLOC, -1);
    
    if ( ! p->output->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as output (%s) could not be created.\n", p->outputpath);
        return;
    }
    
    if ( p->verbose ) {
        fprintf(stdout, "Output Dim: %d %d %d (%d Volumes)\n", p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim);
    }
    
    p->parametersValid = TRUE;
    return;
}

void rsSampleFileGenRun(rsSampleFileGenParameters *p)
{
    p->parametersValid = FALSE;

    // Run
    short x,y,z, processedSlices = 0;
    double *signal;
    double *signal2;
    Point3D *point;
        
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(z,y,x,signal,point) shared(processedSlices)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<p->output->zDim; z++) {
            for (y=0; y<p->output->yDim; y++) {
                for (x=0; x<p->output->xDim; x++) {
                    
                    point = rsMakePoint3D(x, y, z);
                    
                    /* create timecourse from input */
                    signal = rsMalloc(p->output->vDim*sizeof(double));
                    if (p->checkerboardBlockLength > 0) {
                        rsSampleFileGenCheckerboard(signal, point, p->output->vDim, p);
                    } else if (p->sphereRepetitionLength > 0) {
                        rsSampleFileGenSpheres(signal, point, p->output->vDim, p);
                    } else if (p->checkerSphereNX > 0) {
                        rsSampleFileGenCheckerSphere(signal, point, p->output->vDim, p);
                    }
                    
                    /* write out filtered data to buffer */
                    rsWriteTimecourseToRSNiftiFileBuffer(p->output, signal, point);

                    rsFree(signal);
                }
            }
            
            /* show progress */
            if (p->verbose) {
                #pragma omp atomic
                processedSlices += 1;
            
                if (processedSlices > 0 && processedSlices % (short)(p->output->zDim / 10) == 0) {
                    fprintf(stdout, "..%.0f%%\n", ceil((float)processedSlices*100.0 / (float)p->output->zDim));
                }
            }
        }
    }
    
    if ( p->verbose ) {
        fprintf(stdout, "Write out result to: %s\n", p->outputpath);
    }
    
    rsWriteNiftiHeader(p->output->fslio, p->callString);
    FslWriteVolumes(p->output->fslio, p->output->data, p->output->vDim);

    p->parametersValid = TRUE;
}

void rsSampleFileGenCheckerboard(double *output, const Point3D *point, const int length, const rsSampleFileGenParameters *p) {
    const int tileIdxX = point->x / p->checkerboardBlockLength;
    const int tileIdxY = point->y / p->checkerboardBlockLength;
    const int tileIdxZ = point->z / p->checkerboardBlockLength;
    
    const BOOL isInBlock = (tileIdxX % 2 == 0) ^ (tileIdxY % 2 == 0) ^ (tileIdxZ % 2 == 0);
    
    for (int i=0; i<length; i++) {
        output[i] = isInBlock ? 1.0 : 0.0;
    }
}

void rsSampleFileGenCheckerSphere(double *output, const Point3D *point, const int length, const rsSampleFileGenParameters *p) {
    const double distance = pow(
        (
            pow((double)p->spherePoint->x - (double)point->x, 2.0) +
            pow((double)p->spherePoint->y - (double)point->y, 2.0) +
            pow((double)p->spherePoint->z - (double)point->z, 2.0)
        ),
        1.0/2.0
    );
    const int intDistance = ceil(distance);
    
    const BOOL isInRadialBlock = (intDistance % (p->checkerSphereTL * 2)) < p->checkerSphereTL;

    const SignedPoint3D* relPoint = rsMakeSignedPoint3D(p->spherePoint->x - point->x, p->spherePoint->y - point->y, p->spherePoint->z - point->z);

    const float normX = sqrt(relPoint->y*relPoint->y + relPoint->z*relPoint->z);
    const float angleX = acos((relPoint->y + relPoint->z) / (normX*sqrt(2.0))); // scalar product with vector [1 1]
    
    const int verticalBlockIndex = round((angleX / M_PI) * (p->checkerSphereNX-1));
    const BOOL isInVerticalBlock = (verticalBlockIndex % 2) == 0;
        
    const float normY = sqrt(relPoint->x*relPoint->x + relPoint->y*relPoint->y);
    const float angleY = acos((relPoint->x + relPoint->y) / (normY*sqrt(2.0))); // scalar product with vector [1 1]
    
    const int horizontalBlockIndex = round((angleY / M_PI) * (p->checkerSphereNY-1));
    const BOOL isInHorizontalBlock = (horizontalBlockIndex % 2) == 0;
    
    const double value = (isInRadialBlock ? 1.0 : -1.0) * ((isInVerticalBlock ? 1.0 : 0.0) + (isInHorizontalBlock ? 1.0 : 0.0));
    
    for (int i=0; i<length; i++) {
        output[i] = value;
    }
}

void rsSampleFileGenSpheres(double *output, const Point3D *point, const int length, const rsSampleFileGenParameters *p) {
    const double distance = pow(
        (
            pow((double)p->spherePoint->x - (double)point->x, 2.0) +
            pow((double)p->spherePoint->y - (double)point->y, 2.0) +
            pow((double)p->spherePoint->z - (double)point->z, 2.0)
        ),
        1.0/2.0
    );
    const int intDistance = ceil(distance);
    
    const BOOL isInBlock = (intDistance % (p->sphereRepetitionLength)) < p->sphereWidth;
    
    for (int i=0; i<length; i++) {
        output[i] = isInBlock ? 1.0 : 0.0;
    }
}

void rsSampleFileGenDestroy(rsSampleFileGenParameters *p)
{
    if ( p->reference != NULL ) {
        rsCloseNiftiFileAndFree(p->reference);
        p->reference = NULL;
    }
    
    if ( p->output != NULL ) {
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }
    
    rsSampleFileGenFreeParams(p);
}
