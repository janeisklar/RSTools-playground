#ifndef rstools_fixcorrelation_ui_h
#define  rstools_fixcorrelation_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *inputPath;
    char *outputPath;
    char *correlationPath;
    char *referencePath;
    char *stdPath;
    char *meanPath;
    char *callString;

    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *input;
    rsNiftiFile *correlationFile;
    rsNiftiFile *output;
    rsNiftiFile *stdFile;
    rsNiftiFile *meanFile;

    double *reference;
    unsigned int nReferenceValues;

    rsUIInterface *interface;

    int threads;
    
} rsFixCorrelationParameters;

rsFixCorrelationParameters *rsFixCorrelationParseParams(int argc, char * argv[]);
rsFixCorrelationParameters *rsFixCorrelationInitParameters();
void rsFixCorrelationFreeParams(rsFixCorrelationParameters *p);
void rsFixCorrelationBuildInterface(rsFixCorrelationParameters *p);
    
#ifdef __cplusplus
}
#endif

#endif
