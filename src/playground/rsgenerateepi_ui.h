#ifndef rstools_generateepi_ui_h
#define  rstools_generateepi_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *outputPath;
    char *stdPath;
    char *meanPath;
    char *maskPath;
    char *correlationPath;
    char *referencePath;
    char *callString;

    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *stdFile;
    rsNiftiFile *meanFile;
    rsNiftiFile *maskFile;
    rsNiftiFile *correlationFile;
    rsNiftiFile *output;

    double *reference;
    unsigned int nReferenceValues;

    Point3D *maskPoints;
    unsigned long nPoints;

    rsUIInterface *interface;

    int threads;
    size_t wordsize;
    
} rsGenerateEpiParameters;

rsGenerateEpiParameters *rsGenerateEpiParseParams(int argc, char * argv[]);
rsGenerateEpiParameters *rsGenerateEpiInitParameters();
void rsGenerateEpiFreeParams(rsGenerateEpiParameters *p);
void rsGenerateEpiBuildInterface(rsGenerateEpiParameters *p);
    
#ifdef __cplusplus
}
#endif

#endif
