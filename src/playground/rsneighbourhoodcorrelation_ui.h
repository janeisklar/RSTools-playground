#ifndef rstools_neighbourhoodcorrelation_ui_h
#define  rstools_neighbourhoodcorrelation_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *outputPath;
    char *callString;

    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *reference;
    rsNiftiFile *output;

    void *meanCount;

    char **inputPaths;
    size_t nInputs;

    int distance;
    int dimLength;

    rsUIInterface *interface;

    int threads;
    size_t wordsize;
    
} rsNeighbourhoodCorrelationParameters;

rsNeighbourhoodCorrelationParameters *rsGenerateEpiParseParams(int argc, char * argv[]);
rsNeighbourhoodCorrelationParameters *rsGenerateEpiInitParameters();
void rsNeighbourhoodCorrelationFreeParams(rsNeighbourhoodCorrelationParameters *p);
void rsNeighbourhoodCorrelationBuildInterface(rsNeighbourhoodCorrelationParameters *p);
    
#ifdef __cplusplus
}
#endif

#endif
