//
//  rsttest2.c
//  rstools
//
//  Created by André Hoffmann on 9/23/13.
//
//

#include <stdio.h>
#include <string.h>
#include "maths/rsmathutils.h"
#include "nifti/rsniftiutils.h"

#define RSTTEST_OUTPUT_T 0
#define RSTTEST_OUTPUT_P 1
#define RSTTEST_OUTPUT_LN_P 2

void rsTTestPrintHelp() {
    printf(
		RSTOOLS_VERSION_LABEL "\n\n"
        "Takes in a list of niftis with one ore more volumes\n"
        "and performs a one-sample t-test on it.\n"
        "Volumes with the same volume index will be compared.\n"
        "The set of considered niftis are hereby randomly\n"
        "altered for each voxel.\n"
        "basic usage:  stdin | rsvarsubjectsttest2 -output <volume> -n <int>\n"
        "\n"
    );
    
    printf(
        "options:\n"
    );
    
    printf(
        "   -output <volume>       : the volume in which the result will be saved\n"
    );
    
	printf(
    	"   -n <int>               : number of niftis considered\n"
	);
	
	printf(
		"   -repeatedTest <int>    : repeat the randomized t-test n times and\n"
		"                            return the mean p-values\n"
	);
	
	printf(
		"   -p                     : with this option the resulting nifti will be\n"
		"                            containing p-values as opposed to t-values\n"
	);
	
	printf(
		"   -logP                  : much like the -p option, but -ln(p) will be\n"
		"                            used instead\n"
	);

    printf(
        "   -threads <int>         : number of threads used for processing\n"
    );
    
    printf(
        "   -v[erbose]             : show debug information\n"
        "\n"
    );
}

struct rsInputFile {
    char*   path;
	BOOL    readable;
    FSLIO*  fslio;
    short   xDim;
    short   yDim;
    short   zDim;
    short   vDim;
    size_t  dt;
    short   pixtype;
    float   inter;
    float   slope;
    double* data;
};


struct rsInputFile rsOpenInputFile(char* path) {
    struct rsInputFile f;
    
    f.path     = (char*)malloc((strlen(path)+1)*sizeof(char));
    strcpy(f.path, path);
    f.path     = rsTrimString(f.path);
    f.readable = FALSE;
    
    /* open file */
    f.fslio = FslOpen(f.path, "rb");
    
    if (f.fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", f.path);
        return f;
    }
    
    f.readable = TRUE;
    
	/* determine dimensions */
	FslGetDim(f.fslio, &f.xDim, &f.yDim, &f.zDim, &f.vDim);
    
    /* determine scaling */
    f.inter = 0.0;
    f.slope = 1.0;
    
    if (f.fslio->niftiptr->scl_slope != 0) {
        f.slope = f.fslio->niftiptr->scl_slope;
        f.inter = f.fslio->niftiptr->scl_inter;
    }
	
	/* determine datatype */
	f.dt = FslGetDataType(f.fslio, &f.pixtype);
    
    /* read in file */
    f.data = malloc((size_t)f.xDim*(size_t)f.yDim*(size_t)f.zDim*(size_t)f.vDim*(size_t)f.dt/(size_t)8);
    
    FslReadVolumes(f.fslio, f.data, f.vDim);
    
    return f;
}

void rsCloseInputFile(struct rsInputFile* f) {
    free((*f).data);
    FslClose((*f).fslio);
    free((*f).fslio);
    (*f).readable = FALSE;
}

struct rsInputFile *rsReadFileListFromStandardInput(unsigned int *nFiles) {
    char *line = NULL;
    size_t len = 0;
    size_t read;
    int sizeFilesBuffer = 1;
    struct rsInputFile *files = (struct rsInputFile*)malloc(sizeof(struct rsInputFile)*sizeFilesBuffer);
    *nFiles=0;
    
    while ((read = getline(&line, &len, stdin)) != -1) {
        *nFiles = *nFiles + 1;
        
        // Check if we're running out of memory and extend the array if necessary
        if ( *nFiles >= sizeFilesBuffer ) {
            sizeFilesBuffer = sizeFilesBuffer + 10;
            struct rsInputFile* tmpFiles = (struct rsInputFile*)realloc(files, sizeFilesBuffer * sizeof(struct rsInputFile));
            if (tmpFiles) {
                files = tmpFiles;
            } else {
                fprintf(stderr, "Could not allocate enough memory to read the file list from stdin.\n");
                exit(EXIT_FAILURE);
            }
        }
        
        files[*nFiles-1] = rsOpenInputFile(line);
    }
    
    if (line) free(line);
    
    return files;
}

int main(int argc, char * argv[]) {
    
    FSLIO *fslioOutput;
	
	char *outputpath = NULL;
	
	BOOL verbose = FALSE;
    int threads = 1;
	int nNiftis = -1;
	int nRepetitions = 1;
	int statType = RSTTEST_OUTPUT_T;
	
	int ac;
    
	if( argc < 2 ) {
        rsTTestPrintHelp();
        return 1;
    }
    
	/* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			rsTTestPrintHelp();
            return 1;
		} else if ( ! strncmp(argv[ac], "-v", 2) ) {
			verbose = TRUE;
		} else if ( ! strcmp(argv[ac], "-output") ) {
            if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -output\n");
				return 1;
			}
			outputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-threads") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -threads\n");
           		return 1;
           	}
           	threads = atoi(argv[ac]);  /* no string copy, just pointer assignment */
        } else if ( ! strcmp(argv[ac], "-n") ) {
	  		if( ++ac >= argc ) {
	       		fprintf(stderr, "** missing argument for -n\n");
	       		return 1;
	       	}
	       	nNiftis = atoi(argv[ac]);  /* no string copy, just pointer assignment */
	    } else if ( ! strcmp(argv[ac], "-p") ) {
			statType = RSTTEST_OUTPUT_P;
	    } else if ( ! strcmp(argv[ac], "-logP") ) {
   			statType = RSTTEST_OUTPUT_LN_P;
		} else if ( ! strcmp(argv[ac], "-repeatedTest") ) {
	  		if( ++ac >= argc ) {
	       		fprintf(stderr, "** missing argument for -repeatedTest\n");
	       		return 1;
	       	}
	       	nRepetitions = atoi(argv[ac]);  /* no string copy, just pointer assignment */
	    } else {
			fprintf(stderr, "\nError, unrecognized command %s\n", argv[ac]);
		}
	}
    
	if ( outputpath == NULL ) {
		fprintf(stderr, "No output volume specified!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Output file: %s\n", outputpath);
    }
    
    // Load list of files
    unsigned int nFiles = 0;
    struct rsInputFile *files = rsReadFileListFromStandardInput(&nFiles);
	size_t fileListLength = 1;
	
	if ( nNiftis < 0 ) {
		nNiftis = nFiles;
	}
    
    for (int n=0; n<nFiles; n=n+1) {
        const struct rsInputFile file = files[n];
        
        if ( ! file.readable ) {
            fprintf(stderr, "File '%s' is not accessible.\n", file.path);
            return 1;
        }
        
        if (verbose) {
            fprintf(stdout, "File: %s, Volumes: %d\n", file.path, file.vDim);
        }

		fileListLength = fileListLength + strlen(file.path) + 2 + (size_t)fmaxf(rsCountDigits(file.vDim), 4);
    }
    
    if ( nFiles < 1 ) {
        fprintf(stderr, "No files were supplied via standard input!\n");
        return 1;
    }

	char fileList[fileListLength];
	size_t bytesWritten = 0;
    for (int n=0; n<nFiles; n=n+1) {
        const struct rsInputFile file = files[n];
		sprintf(&fileList[bytesWritten], "%s,%04d\n", file.path, file.vDim);
		bytesWritten = bytesWritten + strlen(file.path) + 2 + (size_t)fmaxf(rsCountDigits(file.vDim), 4);
	}
	fileList[fileListLength-1] = '\0';
    
    // Prepare output buffer
    const struct rsInputFile refFile = files[0];
    
    fslioOutput = FslOpen(outputpath, "wb");
    
    if (fslioOutput == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", outputpath);
        return 1;
    }
    
    unsigned int nOutputVolumes = refFile.vDim;
    
    FslCloneHeader(fslioOutput, refFile.fslio);
    FslSetDim(fslioOutput, refFile.xDim, refFile.yDim, refFile.zDim, nOutputVolumes);
    FslSetDimensionality(fslioOutput, 4);
    FslSetDataType(fslioOutput, refFile.pixtype);
	if ( statType == RSTTEST_OUTPUT_P ) {
		FslSetIntent(fslioOutput, NIFTI_INTENT_PVAL, 0, 0, 0);
	} else if ( statType == RSTTEST_OUTPUT_LN_P ) {
		FslSetIntent(fslioOutput, NIFTI_INTENT_LOGPVAL, 0, 0, 0);
	} else {
		FslSetIntent(fslioOutput, NIFTI_INTENT_TTEST, nNiftis-1, 0, 0);
	}
	char *comment1 = rsMergeStringArray(argc, argv);
	char *comment2 = "\nFilelist:\n";
	size_t commentLength = strlen(comment1)+strlen(comment2)+fileListLength+1;
	char comment[commentLength];
	sprintf(&comment[0], "%s%s%s\n", comment1, comment2, fileList);
	comment[commentLength-1] = '\0';	
    rsWriteNiftiHeader(fslioOutput, &comment[0]);
    
    void *buffer = malloc((size_t)refFile.xDim*(size_t)refFile.yDim*(size_t)refFile.zDim*(size_t)nOutputVolumes*(size_t)refFile.dt/(size_t)8);
	gsl_rng *randgen = rsGetRandomNumberGenerator();
	double *tValues;
	double *series;
	double tValue;
	size_t *indices;
	size_t fIndex;
    short t,x,y,z,f,processedSlices = 0;
	int r;
	const int dof=nNiftis-1;
    
    #pragma omp parallel num_threads(threads) private(x,y,t,f,r,tValue,tValues,indices,fIndex,series) shared(randgen,buffer,nFiles,files,nOutputVolumes, nNiftis)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<refFile.zDim; z++) {
			
			tValues = malloc((size_t)nOutputVolumes * sizeof(double));
            indices = malloc((size_t)nFiles         * sizeof(size_t));
			series  = malloc((size_t)nNiftis        * sizeof(double));

			// prepare array with indices that will be shuffled
			for (size_t i = 0; i < (size_t)nFiles; i=i+1) {
				indices[i] = i;
			}
	
            for (y=0; y<refFile.yDim; y=y+1) {
                for (x=0; x<refFile.xDim; x=x+1) {
                    Point3D *point = rsMakePoint3D(x, y, z);

					// shuffle files
					gsl_ran_shuffle(randgen, &indices[0], nFiles, sizeof(size_t));

                    for (t=0; t<refFile.vDim; t=t+1) {
						tValues[t] = 0.0;
						for (r=0; r<nRepetitions; r=r+1) {
	                        for (f=0; f<(int)nNiftis; f=f+1) {
								fIndex = indices[f];
	                            const struct rsInputFile *file = &files[fIndex];

	                            rsExtractPointsFromBuffer(
	                                (*file).pixtype,
	                                &series[f],
	                                (*file).data,
	                                (*file).slope,
	                                (*file).inter,
	                                point,
	                                1L,
	                                t,
	                                (*file).xDim,
	                                (*file).yDim,
	                                (*file).zDim,
	                                (*file).vDim
	                            );
	                        }
                        	
							tValue = rsOneSampleTTest(series, nNiftis, 0.0);

							if ( statType == RSTTEST_OUTPUT_P || statType == RSTTEST_OUTPUT_LN_P ) {
								tValue = rsComputePValueFromTValue(tValue, dof)*2.0;
							}

	                        tValues[t] = tValues[t] + tValue;
	
							if ( nRepetitions > 1 ) {
								gsl_ran_shuffle(randgen, &indices[0], nFiles, sizeof(size_t));
							}
						}
						
						if ( nRepetitions > 1 ) {
							tValues[t] = tValues[t] / (double)nRepetitions;
						}
						
						if ( statType == RSTTEST_OUTPUT_LN_P ) {
							tValues[t] = tValues[t] < 0 ? log(-1.0*tValues[t]) : -log(tValues[t]);
						}
                    }   
                    
                    rsWriteTimecourseToBuffer(refFile.pixtype, tValues, buffer, refFile.slope, refFile.inter, point, refFile.xDim, refFile.yDim, refFile.zDim, nOutputVolumes);
                }
            }
			
			free(series);
			free(indices);
            free(tValues);

			
			/* show progress */
			if (verbose) {
            	#pragma omp atomic
            	processedSlices += 1;
            
            	if (processedSlices > 0 && processedSlices % (short)(refFile.zDim / 10) == 0) {
                	fprintf(stdout, "..%.0f%%\n", ceil((float)processedSlices*100.0 / (float)refFile.zDim));
            	}
			}
        }
    }
    
	rsDestroyRandomNumberGenerator();
    FslWriteVolumes(fslioOutput, buffer, nOutputVolumes);
    
    // Close files
    
    for (int n=0; n<nFiles; n=n+1) {
        struct rsInputFile file = files[n];
        rsCloseInputFile(&file);
    }
    
    FslClose(fslioOutput);
    free(buffer);
    free(files);
}
