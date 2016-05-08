#include "rsgenerateepi_common.h"
#include "rsgenerateepi_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsGenerateEpiParameters *p = rsGenerateEpiParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsGenerateEpiInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsGenerateEpiRun(p);
    }

    // Free memory
    rsGenerateEpiDestroy(p);

    return 0;
}
