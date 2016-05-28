#include "rsfixcorrelation_common.h"
#include "rsfixcorrelation_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsFixCorrelationParameters *p = rsFixCorrelationParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsFixCorrelationInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsFixCorrelationRun(p);
    }

    // Free memory
    rsFixCorrelationDestroy(p);

    return 0;
}
