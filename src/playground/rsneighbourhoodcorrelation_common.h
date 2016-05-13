#ifndef rstools_rsneighbourhoodcorrelation_common_h
#define rstools_rsneighbourhoodcorrelation_common_h

#include "rsneighbourhoodcorrelation_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsNeighbourhoodCorrelationInit(rsNeighbourhoodCorrelationParameters *p);
void rsNeighbourhoodCorrelationRun(rsNeighbourhoodCorrelationParameters *p);
void rsNeighbourhoodCorrelationDestroy(rsNeighbourhoodCorrelationParameters *p);

#ifdef __cplusplus
}
#endif

#endif
