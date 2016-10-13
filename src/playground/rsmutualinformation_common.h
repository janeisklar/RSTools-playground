#ifndef rstools_rsmutualinformation_common_h
#define rstools_rsmutualinformation_common_h

#include "rsmutualinformation_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsMutualInformationInit(rsMutualInformationParameters *p);
void rsMutualInformationRun(rsMutualInformationParameters *p);
double rsMutualInformationComputeMI(const rsMutualInformationParameters *params, const unsigned int indexA, const unsigned int indexB);
void rsMutualInformationDestroy(rsMutualInformationParameters *p);

#ifdef __cplusplus
}
#endif

#endif
