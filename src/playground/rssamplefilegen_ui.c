#include "rssamplefilegen_ui.h"

rsSampleFileGenParameters *rsSampleFileGenInitParameters()
{
    rsSampleFileGenParameters *p = (rsSampleFileGenParameters*)rsMalloc(sizeof(rsSampleFileGenParameters));
    
    p->referencepath           = NULL;
    p->outputpath              = NULL;
    p->callString              = NULL;
    p->verbose                 = FALSE;
    p->reference               = NULL;;
    p->output                  = NULL;
    p->parametersValid         = FALSE;
    p->checkerboardBlockLength = -1;
    p->spherePoint             = NULL;
    p->sphereRepetitionLength  = -1;
    p->sphereWidth             = -1;
    p->checkerSphereNX         = -1;
    p->checkerSphereNY         = -1;
    p->checkerSphereTL         = -1;
    p->threads                 = 1;
    p->interface               = NULL;
        
    return p;
}

void rsSampleFileGenFreeParams(rsSampleFileGenParameters *p)
{
    rsFree(p->referencepath);
    rsFree(p->outputpath);
    rsFree(p->reference);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsSampleFileGenParameters *rsSampleFileGenParseParams(int argc, char * argv[])
{

    rsSampleFileGenParameters *p = rsSampleFileGenInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsSampleFileGenBuildInterface(p);
    
    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void*)p);
    
    if ( ! parsingSuccessful ) {
        return p;
    }
    
    if ( p->referencepath == NULL ) {
        fprintf(stderr, "No reference volume specified(--reference)!\n");
        return p;
    }
    
    if ( p->outputpath == NULL ) {
        fprintf(stderr, "No output volume specified(--output)!\n");
        return p;
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

gboolean rsSampleFileGenParseCheckerSphere(const gchar *option_name, const gchar *value, gpointer data, GError **error)
{
    rsSampleFileGenParameters *p = (rsSampleFileGenParameters*) data;

    // copy value(const)
    size_t length = strlen(value);
    char v[length+1];
    sprintf(&v[0], "%s", value);

    // parse value
    char *strX;
    char *strY;
    char *strZ;
    char *strNX;
    char *strNY;
    char *strTL;
    strX  = strtok(   v, ",");
    strY  = strtok(NULL, ",");
    strZ  = strtok(NULL, ",");
    strNX = strtok(NULL, ",");
    strNY = strtok(NULL, ",");
    strTL = strtok(NULL, ",");

    int x=-1, y=-1, z=-1, nx=-1, ny=-1, tl=-1;

    // if we were given exactly 5 numbers separated by comma parse them 
    if ( strtok(NULL,",") == NULL && strX != NULL && strY != NULL && strZ != NULL && strNX != NULL && strNY != NULL && strTL != NULL ) {
        x = atoi(strX);
        y = atoi(strY);
        z = atoi(strZ);
        nx = atoi(strNX);
        ny = atoi(strNY);
        tl = atoi(strTL);
    }
    
    // return success if we sucessfully received 6 numbers
    if ( x >= 0 && y >= 0 && z >= 0 && nx > 0 && ny > 0 && tl > 0) {
        p->spherePoint = rsMakePoint3D(x,y,z);
        p->checkerSphereNX = nx;
        p->checkerSphereNY = ny;
        p->checkerSphereTL = tl;
        return TRUE;
    }
    
    // anything else should lead to an error
    g_set_error(
        error,
        G_OPTION_ERROR,
        G_OPTION_ERROR_BAD_VALUE,
        "%s: %s",
        option_name,
        "format should be 'x,y,z,nx,ny,tl'"
    );
    
    return FALSE;
}

gboolean rsSampleFileGenParseSphere(const gchar *option_name, const gchar *value, gpointer data, GError **error)
{
    rsSampleFileGenParameters *p = (rsSampleFileGenParameters*) data;

    // copy value(const)
    size_t length = strlen(value);
    char v[length+1];
    sprintf(&v[0], "%s", value);

    // parse value
    char *strX;
    char *strY;
    char *strZ;
    char *strW;
    char *strR;
    strX = strtok(   v, ",");
    strY = strtok(NULL, ",");
    strZ = strtok(NULL, ",");
    strW = strtok(NULL, ",");
    strR = strtok(NULL, ",");

    int x=-1, y=-1, z=-1, r=-1, w=-1;

    // if we were given exactly 5 numbers separated by comma parse them 
    if ( strtok(NULL,",") == NULL && strX != NULL && strY != NULL && strZ != NULL && strW != NULL && strR != NULL ) {
        x = atoi(strX);
        y = atoi(strY);
        z = atoi(strZ);
        w = atoi(strW);
        r = atoi(strR);
    }
    
    // return success if we sucessfully received 5 numbers
    if ( x >= 0 && y >= 0 && z >= 0 && r > 0 && w > 0) {
        p->spherePoint = rsMakePoint3D(x,y,z);
        p->sphereRepetitionLength = r;
        p->sphereWidth = w;
        return TRUE;
    }
    
    // anything else should lead to an error
    g_set_error(
        error,
        G_OPTION_ERROR,
        G_OPTION_ERROR_BAD_VALUE,
        "%s: %s",
        option_name,
        "format should be 'x,y,z,width,length'"
    );
    
    return FALSE;
}

void rsSampleFileGenBuildInterface(rsSampleFileGenParameters *p)
{
    GOptionArgFunc cbSpheres = (GOptionArgFunc)rsSampleFileGenParseSphere;
    GOptionArgFunc cbCheckerSphere = (GOptionArgFunc)rsSampleFileGenParseCheckerSphere;
    
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Generates sample nifti files (such as a checkerboard pattern)";
    
    o = rsUINewOption();
    o->name                = "reference";
    o->shorthand           = 'r';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->referencepath;
    o->cli_description     = "the reference volume used for specifying the dimensions and orientation";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "output";
    o->shorthand           = 'o';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->outputpath;
    o->cli_description     = "the volume that will be created";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "checkerboard";
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->checkerboardBlockLength;
    o->cli_description     = "creates a checkerboard pattern with a tile length as specified in voxels";
    o->cli_arg_description = "<tileLength>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "checksphere";
    o->type                = G_OPTION_ARG_CALLBACK;
    o->storage             = cbCheckerSphere;
    o->cli_description     = "creates a spherical checkerboard pattern around the point in the specified voxel coordinates (x,y,z) with 'nx' horizontal tiles, 'ny' vertical tiles and a radial tile length 'tl' in voxels";
    o->cli_arg_description = "x,y,z,nx,ny,tl";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "spheres";
    o->type                = G_OPTION_ARG_CALLBACK;
    o->storage             = cbSpheres;
    o->cli_description     = "creates a set of spheres that repeat around a specified point in voxel coordinates (x,y,z) after a specified number of voxels (n) with a line thickness of (w voxels)";
    o->cli_arg_description = "x,y,z,w,n";
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
