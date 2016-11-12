#include <utils/rsui.h>
#include "rsfixcorrelation_ui.h"

rsFixCorrelationParameters *rsFixCorrelationInitParameters()
{
    rsFixCorrelationParameters *p = (rsFixCorrelationParameters *) rsMalloc(sizeof(rsFixCorrelationParameters));

    p->inputPath = NULL;
    p->outputPath = NULL;
    p->stdPath = NULL;
    p->meanPath = NULL;
    p->correlationPath = NULL;
    p->referencePath = NULL;
    p->callString = NULL;
    p->verbose = FALSE;
    p->input = NULL;
    p->correlationFile = NULL;
    p->output = NULL;
    p->parametersValid = FALSE;
    p->threads = 1;
    p->interface = NULL;

    return p;
}

void rsFixCorrelationFreeParams(rsFixCorrelationParameters *p)
{
    rsFree(p->input);
    rsFree(p->outputPath);
    rsFree(p->inputPath);
    rsFree(p->correlationPath);
    rsFree(p->referencePath);
    rsFree(p->stdPath);
    rsFree(p->meanPath);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsFixCorrelationParameters *rsFixCorrelationParseParams(int argc, char *argv[])
{

    rsFixCorrelationParameters *p = rsFixCorrelationInitParameters();
    p->callString = rsMergeStringArray(argc, argv);

    rsFixCorrelationBuildInterface(p);

    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void *) p);

    if (!parsingSuccessful) {
        return p;
    }

    if (p->inputPath == NULL) {
        fprintf(stderr, "No input volume specified(--input)!\n");
        return p;
    }

    if (p->correlationPath == NULL) {
        fprintf(stderr, "No correlation volume specified(--correlation)!\n");
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

    p->parametersValid = parsingSuccessful;
    return p;
}

void rsFixCorrelationBuildInterface(rsFixCorrelationParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description =
        "Given a reference timecourse, this program artificially modifies the timecourses in a 4D-nifti file such that they correlate exactly with the pearson correlation coefficient that is specified in a supplied 3D correlation map. This is obviously intended to be used for validation purposes and not for the 'improvement' of study results:)";

    o = rsUINewOption();
    o->name = "input";
    o->shorthand = 'i';
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->inputPath;
    o->cli_description = "the 4D-nifti file that is altered to adhere to have the given correlation coefficients";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

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
    o->name = "reference";
    o->type = G_OPTION_ARG_FILENAME;
    o->storage = &p->referencePath;
    o->cli_description =
        "a txt-file that contains the reference timecourse. The output timecourses will correlate with this timecourse by the r-values specified using --correlation";
    o->cli_arg_description = "<txt-file>";
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
