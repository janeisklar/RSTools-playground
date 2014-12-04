#ifndef rstools_rsfit_common_h
#define rstools_rsfit_common_h

#include "rsfit_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsFitInit(rsFitParameters *p);
void rsFitRun(rsFitParameters *p);
void rsFitDestroy(rsFitParameters *p);
void rsFitRegression(double *target, const double *input, const int length);

#ifdef __cplusplus
}
#endif

#endif
