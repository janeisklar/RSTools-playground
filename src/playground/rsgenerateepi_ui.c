#include <utils/rsui.h>
#include "rsgenerateepi_ui.h"

rsGenerateEpiParameters *rsGenerateEpiInitParameters()
{
    rsGenerateEpiParameters *p = (rsGenerateEpiParameters *) rsMalloc(sizeof(rsGenerateEpiParameters));

    p->outputPath = NULL;
    p->stdPath = NULL;
    p->meanPath = NULL;
    p->maskPath = NULL;
    p->neighbourhoodCorrelationPath = NULL;
    p->correlationPath = NULL;
    p->referencePath = NULL;
    p->callString = NULL;
    p->verbose = FALSE;
    p->stdFile = NULL;
    p->maskFile = NULL;
    p->meanFile = NULL;
    p->correlationFile = NULL;
    p->output = NULL;
    p->parametersValid = FALSE;
    p->kernelSize = 21;
    p->TR = -1.0;
    p->sigma = -1.0;
    p->threads = 1;
    p->interface = NULL;

    return p;
}

void rsGenerateEpiFreeParams(rsGenerateEpiParameters *p)
{
    rsFree(p->outputPath);
    rsFree(p->stdPath);
    rsFree(p->meanPath);
    rsFree(p->maskPath);
    rsFree(p->correlationPath);
    rsFree(p->referencePath);
    rsFree(p->neighbourhoodCorrelationPath);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsGenerateEpiParameters *rsGenerateEpiParseParams(int argc, char *argv[])
{

    rsGenerateEpiParameters *p = rsGenerateEpiInitParameters();
    p->callString = rsMergeStringArray(argc, argv);

    rsGenerateEpiBuildInterface(p);

    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void *) p);

    if (!parsingSuccessful) {
        return p;
    }

    if (p->neighbourhoodCorrelationPath == NULL) {
        fprintf(stderr, "No correlation volume specified(--neighbourhoodCorrelation)!\n");
        return p;
    }

    if (p->outputPath == NULL) {
        fprintf(stderr, "No output volume specified(--output)!\n");
        return p;
    }

    if (p->referencePath == NULL) {
        fprintf(stderr, "No reference file specified(--reference)!\n");
        return p;
    }

    if (p->maskPath == NULL) {
        fprintf(stderr, "No mask file specified(--mask)!\n");
        return p;
    }

    if (p->TR < 0.01) {
        fprintf(stderr, "TR needs to be specified(--TR)!\n");
        return p;
    }

    if (p->sigma < 0.0001) {
        fprintf(stderr, "The sigma for the temporal smoothing needs to be specified(--sigma)!\n");
        return p;
    }

    p->parametersValid = parsingSuccessful;
    return p;
}

void rsGenerateEpiBuildInterface(rsGenerateEpiParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description =
        "Generates a 4D-nifti file where the correlation between the given reference timecourse and the timecourse of each voxel is sampled such that it matches that from the supplied nifti (--correlation) containing the r-values for each voxel. Furthermore the desired standard deviation and mean can be supplied using --mean and --std which will be incorporated as well. This tool can be used for generating datasets that serve as a ground truth for validating seed voxel regressions.  Please note that for computability/performance reasons the temporal continuity of the dataset is not given as it does not play a role in computing the correlation.";

    o = rsUINewOption();
    o->name = "correlation";
    o->shorthand = 'c';
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->correlationPath;
    o->cli_description = "the reference volume used for specifying the dimensions and orientation";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "output";
    o->shorthand = 'o';
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->outputPath;
    o->cli_description = "the volume that will be created";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "neighbourhoodCorrelation";
    o->shorthand = 'n';
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->neighbourhoodCorrelationPath;
    o->cli_description = "a 4D-nifti created with rsneighbourhoodcorrelation that is applied to the randomly generated data to recover its intrinsic structure";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "mean";
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->meanPath;
    o->cli_description = "a 3D-nifti specifying the mean of the resulting output";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "std";
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->stdPath;
    o->cli_description = "a 3D-nifti specifying the standard deviation of the resulting output";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "mask";
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->maskPath;
    o->cli_description = "a 3D-nifti mask-file that specifies the voxels which will be sampled (select as few voxels as possible as the processing is rather time-consuming)";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "reference";
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->referencePath;
    o->cli_description =
        "a txt-file that contains the reference timecourse. The output timecourses will correlate with this timecourse by the r-values specified using --correlation";
    o->cli_arg_description = "<txt-file>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "sigma";
    o->type = G_OPTION_ARG_DOUBLE;
    o->storage = &p->sigma;
    o->cli_description = "The sigma of the gaussian smoothing kernel that is applied to the timecourses";
    o->cli_arg_description = "<float>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name = "TR";
    o->type = G_OPTION_ARG_DOUBLE;
    o->storage = &p->TR;
    o->cli_description = "The repetition rate (used in conjunction with sigma for the temporal smoothing of the data). Should also be equal to the TR used to create the standard deviation file.";
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
