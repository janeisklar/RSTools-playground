#include "rssamplefilegen_common.h"
#include "rssamplefilegen_ui.h"

int main(int argc, char * argv[])
{    
    // Parse run arguments
    rsSampleFileGenParameters *p = rsSampleFileGenParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsSampleFileGenInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsSampleFileGenRun(p);
    }

    // Free memory
    rsSampleFileGenDestroy(p);

    return 0;
}
