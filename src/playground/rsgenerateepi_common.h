#ifndef rstools_rsgenerateepi_common_h
#define rstools_rsgenerateepi_common_h

#include "rsgenerateepi_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsGenerateEpiInit(rsGenerateEpiParameters *p);
void rsGenerateEpiRun(rsGenerateEpiParameters *p);
void rsGenerateEpiDestroy(rsGenerateEpiParameters *p);

#ifdef __cplusplus
}
#endif

#endif
