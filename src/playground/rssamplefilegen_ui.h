#ifndef rstools_samplefilegen_ui_h
#define  rstools_samplefilegen_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *outputpath;
    char *referencepath;
    char *callString;
    
    int checkerboardBlockLength;
    int checkerSphereNX;
    int checkerSphereNY;
    int checkerSphereTL;
    Point3D *spherePoint;
    FloatPoint3D *euclideanCenter;
    int sphereRepetitionLength;
    int sphereWidth;
    
    BOOL verbose;
    BOOL parametersValid;
    
    rsNiftiFile *reference;
    rsNiftiFile *output;

    rsUIInterface *interface;

    int threads;
    size_t wordsize;
    
} rsSampleFileGenParameters;

rsSampleFileGenParameters *rsSampleFileGenParseParams(int argc, char * argv[]);
rsSampleFileGenParameters *rsSampleFileGenInitParameters();
void rsSampleFileGenFreeParams(rsSampleFileGenParameters *p);
void rsSampleFileGenBuildInterface(rsSampleFileGenParameters *p);
    
#ifdef __cplusplus
}
#endif

#endif
