#include "rsfit_ui.h"

rsFitParameters *rsFitInitParameters()
{
    rsFitParameters *p = (rsFitParameters*)rsMalloc(sizeof(rsFitParameters));
    
    p->inputpath            = NULL;
    p->targetpath           = NULL;
    p->betaspath            = NULL;
    p->callString           = NULL;
    p->verbose              = FALSE;
    p->input                = NULL;;
    p->betas                = NULL;
    p->target               = NULL;
    p->parametersValid      = FALSE;
    p->threads              = 1;
    p->interface            = NULL;
        
    return p;
}

void rsFitFreeParams(rsFitParameters *p)
{
    rsFree(p->inputpath);
    rsFree(p->targetpath);
    rsFree(p->betaspath);
    rsFree(p->input);
    rsFree(p->betas);
    rsFree(p->target);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsFitParameters *rsFitParseParams(int argc, char * argv[])
{

    rsFitParameters *p = rsFitInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsFitBuildInterface(p);
    
    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void*)p);
    
    if ( ! parsingSuccessful ) {
        return p;
    }
    
    if ( p->inputpath == NULL ) {
        fprintf(stderr, "No input volume specified(--input)!\n");
        return p;
    }
    
    if ( p->targetpath == NULL ) {
        fprintf(stderr, "No target volume specified(--target)!\n");
        return p;
    }
    
    if ( p->betaspath == NULL ) {
        fprintf(stderr, "No betas volume specified(--betas)!\n");
        return p;
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

void rsFitBuildInterface(rsFitParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Performs a temporal regression of the following format: target = input * betas. The constant term is added automatically. The resulting betas will contain 2 volumes with the first refering to the constant term and the second to the input.";
    
    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the input volume used for regression";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "target";
    o->shorthand           = 't';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->targetpath;
    o->cli_description     = "the target volume that will be regressed on";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "betas";
    o->shorthand           = 'b';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->betaspath;
    o->cli_description     = "the betas as explained int the program description";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "threads";
    o->shorthand           = 't';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->threads;
    o->cli_description     = "number of threads used for processing";
    o->cli_arg_description = "<n>";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
}
