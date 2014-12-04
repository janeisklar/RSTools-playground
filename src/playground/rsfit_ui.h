#ifndef rstools_fit_ui_h
#define  rstools_fit_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *inputpath;
    char *targetpath;
    char *betaspath;
    char *callString;
    
    BOOL verbose;
    BOOL parametersValid;
    
    rsNiftiFile *input;
    rsNiftiFile *target;
    rsNiftiFile *betas;

    rsUIInterface *interface;

    int threads;
    size_t wordsize;
    
} rsFitParameters;

rsFitParameters *rsFitParseParams(int argc, char * argv[]);
rsFitParameters *rsFitInitParameters();
void rsFitFreeParams(rsFitParameters *p);
void rsFitBuildInterface(rsFitParameters *p);
    
#ifdef __cplusplus
}
#endif

#endif
