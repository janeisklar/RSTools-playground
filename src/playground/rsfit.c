#include "rsfit_common.h"
#include "rsfit_ui.h"

int main(int argc, char * argv[])
{    
    // Parse run arguments
    rsFitParameters *p = rsFitParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsFitInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsFitRun(p);
    }

    // Free memory
    rsFitDestroy(p);

    return 0;
}
