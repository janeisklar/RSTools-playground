#include <utils/rsui.h>
#include "rsmutualinformation_ui.h"

rsMutualInformationParameters *rsMutualInformationInitParameters()
{
    rsMutualInformationParameters *p = (rsMutualInformationParameters *) rsMalloc(sizeof(rsMutualInformationParameters));

    p->maskPath = NULL;
    p->callString = NULL;
    p->verbose = FALSE;
    p->maskFile = NULL;
    p->filePaths = NULL;
    p->files = NULL;
    p->nFiles = 0;
    p->kernelFWHM = 0.0;
    p->parametersValid = FALSE;
    p->threads = 1;
    p->interface = NULL;

    return p;
}

void rsMutualInformationFreeParams(rsMutualInformationParameters *p)
{
    rsFree(p->maskPath);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsMutualInformationParameters *rsMutualInformationParseParams(int argc, char *argv[])
{

    rsMutualInformationParameters *p = rsMutualInformationInitParameters();
    p->callString = rsMergeStringArray(argc, argv);

    rsMutualInformationBuildInterface(p);

    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void *) p);

    if (!parsingSuccessful) {
        return p;
    }

    if (p->maskPath == NULL) {
        fprintf(stderr, "No mask specified (--mask)!\n");
        return p;
    }

    p->nFiles = p->interface->nExtraArguments-1;

    if (p->nFiles < 2) {
        fprintf(stderr, "Specify at least two files that are to be compared (got: %d)\n", (int)p->nFiles);
        return p;
    }

    p->filePaths = rsMalloc((p->nFiles+1) * sizeof(char*));

    for (unsigned int i=1; i<p->interface->nExtraArguments; i++) {
        p->filePaths[i-1] = rsString(argv[i]);
    }
    p->filePaths[p->nFiles] = RSIO_LASTFILE;

    p->parametersValid = parsingSuccessful;
    return p;
}

void rsMutualInformationBuildInterface(rsMutualInformationParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description =
        "Computes and outputs a matrix holding the normalized mutual information, i.e. (H(A) + H(B)) / H(A,B) for all images supplied (which need to be of the same size and orientation).\n\nUsage:\nrsmutualinformation [options] -- file1.nii file2.nii ..";

    o = rsUINewOption();
    o->name = "mask";
    o->shorthand = 'm';
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->maskPath;
    o->cli_description = "a 3D-nifti mask-file that specifies the voxels which will be sampled (select as few voxels as possible as the processing is rather time-consuming)";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "smoothing";
    o->shorthand           = 's';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->kernelFWHM;
    o->cli_description     = "The FWHM size of a gaussian kernel in the number of histogram bins. If specified the joint histogram will additionally be smoothed prior to calculating the mutual information. The histogram will be of the size 256x256. The specified kernel width should thus be less than 256 (though probably not more than 10 will be used in practice).";
    o->cli_arg_description = "<float>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "threads";
    o->shorthand = 't';
    o->type = G_OPTION_ARG_INT;
    o->storage = &p->threads;
    o->cli_description = "number of threads used for processing";
    o->cli_arg_description = "<n>";
    o->showInGUI = FALSE;
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "verbose";
    o->shorthand = 'v';
    o->storage = &p->verbose;
    o->cli_description = "show debug information";
    o->showInGUI = FALSE;
    rsUIAddOption(p->interface, o);
}
