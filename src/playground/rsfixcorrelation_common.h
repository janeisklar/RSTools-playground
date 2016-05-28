#ifndef rstools_rsfixcorrelation_common_h
#define rstools_rsfixcorrelation_common_h

#include "rsfixcorrelation_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsFixCorrelationInit(rsFixCorrelationParameters *p);
void rsFixCorrelationRun(rsFixCorrelationParameters *p);
void rsFixCorrelationDestroy(rsFixCorrelationParameters *p);

#ifdef __cplusplus
}
#endif

#endif
