#include "rsfit_common.h"
#include "rsfit_ui.h"
#include "utils/rsio.h"

void rsFitInit(rsFitParameters *p)
{
    p->parametersValid = FALSE;
    
    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        (const char*)p->targetpath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->betaspath,
        RSIO_LASTFILE
    });

    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }
    
    rsSetThreadsNum(p->threads);

    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Target file: %s\n", p->targetpath);
        fprintf(stdout, "Betas file: %s\n", p->betaspath);
    }
    
    p->input  = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);
    p->target = rsOpenNiftiFile(p->targetpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }
    
    if ( ! p->target->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied for the target (%s) could not be read.\n", p->targetpath);
        return;
    }
       
    if ( p->verbose ) {
        fprintf(stdout, "Input Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
        fprintf(stdout, "Target Dim: %d %d %d (%d Volumes)\n", p->target->xDim, p->target->yDim, p->target->zDim, p->target->vDim);
    }

    // Create output volume
    p->betas = rsCloneNiftiFile(p->betaspath, p->input, RSNIFTI_OPEN_ALLOC, 2);
    
    if ( ! p->betas->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the betas (%s) could not be created.\n", p->betaspath);
        return;
    }
    
    p->parametersValid = TRUE;
    return;
}

void rsFitRun(rsFitParameters *p)
{
    p->parametersValid = FALSE;

    // Run
    short x,y,z, processedSlices = 0;
    double *signal;
    double *signal2;
    Point3D *point;
        
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(z,y,x,signal,signal2,point) shared(processedSlices)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<p->input->zDim; z++) {
            for (y=0; y<p->input->yDim; y++) {
                for (x=0; x<p->input->xDim; x++) {
                    
                    point = rsMakePoint3D(x, y, z);
                    
                    /* read out timecourse from input */
                    signal = rsMalloc(p->input->vDim*sizeof(double));
                    rsExtractTimecourseFromRSNiftiFileBuffer(p->input, signal, point);
                    
                    /* read out timecourse from target */
                    signal2 = rsMalloc(p->target->vDim*sizeof(double));
                    rsExtractTimecourseFromRSNiftiFileBuffer(p->target, signal2, point);
                    
                    /* apply regression */
                    rsFitRegression(signal2, signal, p->input->vDim);
                    
                    /* write out filtered data to buffer */
                    rsWriteTimecourseToRSNiftiFileBuffer(p->betas, signal2, point);

                    rsFree(signal);
                    rsFree(signal2);
                }
            }
            
            /* show progress */
            if (p->verbose) {
                #pragma omp atomic
                processedSlices += 1;
            
                if (processedSlices > 0 && processedSlices % (short)(p->input->zDim / 10) == 0) {
                    fprintf(stdout, "..%.0f%%\n", ceil((float)processedSlices*100.0 / (float)p->input->zDim));
                }
            }
        }
    }
    
    if ( p->verbose ) {
        fprintf(stdout, "Write out result to: %s\n", p->betaspath);
    }
    
    rsWriteNiftiHeader(p->betas->fslio, p->callString);
    FslWriteVolumes(p->betas->fslio, p->betas->data, p->betas->vDim);

    p->parametersValid = TRUE;
}

void rsFitRegression(double *target, const double *input, const int length) {

    const int nRegressors = 2;
    const int nSamples    = length;
    
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nSamples, nRegressors);

    // convert inputs to gsl variables
    gsl_vector *y   = gsl_vector_alloc(nSamples);                 // input/signal
    gsl_matrix *X   = gsl_matrix_alloc(nSamples, nRegressors);    // regressors
    gsl_vector *b   = gsl_vector_alloc(nRegressors);              // betas
    gsl_matrix *cov = gsl_matrix_alloc(nRegressors, nRegressors); // covariance
    //gsl_vector *res = gsl_vector_alloc(nSamples);                 // residuals

    // fill them
    for (int i=0; i < nSamples; i++){
        gsl_vector_set(y, i, target[i]);
    }

    gsl_matrix_set_all(X, 1.0);
    for (int t=0; t<nSamples; t++) {
        gsl_matrix_set(X,t,1,input[t]);
    }

    // remove the regressors' mean
    for (int i=0; i<nRegressors; i=i+1) {
        
        gsl_vector_view regressor = gsl_matrix_column(X, i);
    
        const double meanX = gsl_stats_mean(regressor.vector.data, 1, nSamples);
    
        gsl_vector_add_constant(&regressor.vector, -1.0 * meanX); // remove mean
    }

    // execute linear regression
    double chisq;
    int success = gsl_multifit_linear(X, y, b, cov, &chisq, work);

    // compute residuals
    //gsl_multifit_linear_residuals(X, y, b, res);

    // convert back to basic C variables
    target[0] = gsl_vector_get(b, 0);
    target[1] = gsl_vector_get(b, 1);

    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(b);
    gsl_multifit_linear_free(work);
    
}

void rsFitDestroy(rsFitParameters *p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }
    
    if ( p->target != NULL ) {
        rsCloseNiftiFileAndFree(p->target);
        p->target = NULL;
    }
    
    rsFitFreeParams(p);
}
