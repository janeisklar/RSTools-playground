#include "rsneighbourhoodcorrelation_common.h"
#include "rsneighbourhoodcorrelation_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsNeighbourhoodCorrelationParameters *p = rsGenerateEpiParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsNeighbourhoodCorrelationInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsNeighbourhoodCorrelationRun(p);
    }

    // Free memory
    rsNeighbourhoodCorrelationDestroy(p);

    return 0;
}
