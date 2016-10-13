#include "rsmutualinformation_common.h"
#include "rsmutualinformation_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsMutualInformationParameters *p = rsMutualInformationParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsMutualInformationInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsMutualInformationRun(p);
    }

    // Free memory
    rsMutualInformationDestroy(p);

    return 0;
}
