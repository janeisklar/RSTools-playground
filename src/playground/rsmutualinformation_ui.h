#ifndef rstools_mutualinformation_ui_h
#define  rstools_mutualinformation_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char **filePaths;
    char *maskPath;
    char *callString;

    BOOL verbose;
    BOOL parametersValid;
    BOOL normalizedMI;

    double **files;
    unsigned int nFiles;

    double *fileMinima;
    double *fileMaxima;

    double kernelFWHM;

    rsNiftiFile *maskFile;
    Point3D *maskPoints;
    unsigned long nPoints;

    rsUIInterface *interface;

    int threads;
    size_t wordsize;
    
} rsMutualInformationParameters;

rsMutualInformationParameters *rsMutualInformationParseParams(int argc, char * argv[]);
rsMutualInformationParameters *rsMutualInformationInitParameters();
void rsMutualInformationFreeParams(rsMutualInformationParameters *p);
void rsMutualInformationBuildInterface(rsMutualInformationParameters *p);
    
#ifdef __cplusplus
}
#endif

#endif
