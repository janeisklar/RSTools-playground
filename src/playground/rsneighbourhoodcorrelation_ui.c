#include "rsneighbourhoodcorrelation_ui.h"

rsNeighbourhoodCorrelationParameters *rsGenerateEpiInitParameters()
{
    rsNeighbourhoodCorrelationParameters *p = (rsNeighbourhoodCorrelationParameters *) rsMalloc(sizeof(rsNeighbourhoodCorrelationParameters));

    p->outputPath = NULL;
    p->callString = NULL;
    p->verbose = FALSE;
    p->inputPaths = NULL;
    p->nInputs = 0;
    p->output = NULL;
    p->parametersValid = FALSE;
    p->distance = -1;
    p->threads = 1;
    p->interface = NULL;

    return p;
}

void rsNeighbourhoodCorrelationFreeParams(rsNeighbourhoodCorrelationParameters *p)
{
    rsFree(p->outputPath);
    for (size_t i = 0; i<p->nInputs; i++) {
        rsFree(p->inputPaths[i]);
    }
    rsFree(p->inputPaths);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsNeighbourhoodCorrelationParameters *rsGenerateEpiParseParams(int argc, char *argv[])
{

    rsNeighbourhoodCorrelationParameters *p = rsGenerateEpiInitParameters();
    p->callString = rsMergeStringArray(argc, argv);

    rsNeighbourhoodCorrelationBuildInterface(p);

    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void *) p);

    if (!parsingSuccessful) {
        return p;
    }

    if (p->outputPath == NULL) {
        fprintf(stderr, "No output volume specified(--output)!\n");
        return p;
    }

    if (p->distance < 1) {
        fprintf(stderr, "The neighbourhood distance has to be specified and needs to be larger than 0!\n");
        return p;
    }

    char *line = NULL;
    size_t len = 0;
    size_t read;
    int sizeFilesBuffer = 5;
    p->inputPaths = (char**)rsMalloc(sizeof(char*)*sizeFilesBuffer);
    p->nInputs = 0;

    while ((read = getline(&line, &len, stdin)) != -1) {
        p->nInputs++;

        // Check if we're running out of memory and extend the array if necessary
        if ( (p->nInputs+1) >= sizeFilesBuffer ) {
            sizeFilesBuffer = sizeFilesBuffer + 10;
            char **tmpFiles = (char**)realloc(p->inputPaths, sizeFilesBuffer * sizeof(char*));
            if (tmpFiles) {
                p->inputPaths = tmpFiles;
            } else {
                fprintf(stderr, "Could not allocate enough memory to read the file list from stdin.\n");
                exit(EXIT_FAILURE);
            }
        }

        char *str = rsTrimString(rsString(line));
        const size_t strLength = strlen(str);
        if (strLength > 0 && (str[strLength-1] == '\n' || str[strLength-1] == '\r')) {
            str[strLength-1] = '\0';
        }
        p->inputPaths[p->nInputs-1] = str;
    }
    p->inputPaths[p->nInputs] = RSIO_LASTFILE;

    rsFree(line);

    if (p->nInputs < 1) {
        fprintf(stderr, "No inputs specified via stdin!\n");
        return p;
    }

    p->parametersValid = parsingSuccessful;
    return p;
}

void rsNeighbourhoodCorrelationBuildInterface(rsNeighbourhoodCorrelationParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description =
        "This program will go through the list of files that has to be provided through standard input (run ls *.nii | rsneighbourhoodcorrelation -o ...) and creates a correlation matrix for every point of all neighbouring voxels. All the neighbourhood correlation matrices computed for one point are then meaned (by considering all files). The size of the neighbourhood has to be specified using the distance parameter and obviously influences the dimensions of the resulting correlation matrices: (2*distance+1)^3. The first three dimensions of the output will correspond to the location for which the neighourhood correlation matrix has been computed. The correlation matrix itself will be stored in the 4th dimension. The correlation value for a location (x,y,z) within the correlation matrix will be stored in the index i of the 4th volume as given by the following formula: i=x+y*(2*distance+1)+z*(2*distance+1)^2";

    o = rsUINewOption();
    o->name = "output";
    o->shorthand = 'o';
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->outputPath;
    o->cli_description = "the volume that will be created";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "distance";
    o->shorthand = 'd';
    o->type = G_OPTION_ARG_INT;
    o->storage = &p->distance;
    o->cli_description = "number of voxels to consider in each direction. the resulting matrix will be of the size (distance*2 + 1)^3";
    o->cli_arg_description = "<n>";
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
    o->name = "verbose";
    o->shorthand = 'v';
    o->storage = &p->verbose;
    o->cli_description = "show debug information";
    o->showInGUI = FALSE;
    rsUIAddOption(p->interface, o);
}
