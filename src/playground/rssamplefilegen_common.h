#ifndef rstools_rssamplefilegen_common_h
#define rstools_rssamplefilegen_common_h

#include "rssamplefilegen_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsSampleFileGenInit(rsSampleFileGenParameters *p);
void rsSampleFileGenRun(rsSampleFileGenParameters *p);
void rsSampleFileGenDestroy(rsSampleFileGenParameters *p);
void rsSampleFileGenCheckerboard(double *output, const Point3D *point, const int length, const rsSampleFileGenParameters *p);
void rsSampleFileGenCheckerSphere(double *output, const Point3D *point, const int length, const rsSampleFileGenParameters *p);
void rsSampleFileGenSpheres(double *output, const Point3D *point, const int length, const rsSampleFileGenParameters *p);

#ifdef __cplusplus
}
#endif

#endif
