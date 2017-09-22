#include "heterogauss.h"
#include <stdio.h>
#include <stdlib.h>
#include "nifti1_io.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <omp.h>

#include "HelpFunctions.cpp"

#define ADD_FILENAME true
#define DONT_ADD_FILENAME true

#define CHECK_EXISTING_FILE true
#define DONT_CHECK_EXISTING_FILE false

int main(int argc, char **argv)
{
	//-----------------------
    // Input
    
    float           *h_Data, *h_Mask; 
  
    float           *h_X_GLM, *h_xtxxt_GLM, *h_X_GLM_Confounds, *h_Contrasts, *h_ctxtxc_GLM, *h_Highres_Regressors, *h_LowpassFiltered_Regressors, *h_Motion_Parameters, *h_Motion_Deriv_Parameters;
                  
    size_t          DATA_W, DATA_H, DATA_D, DATA_T, NUMBER_OF_RUNS;  
	size_t*			DATA_T_PER_RUN;             
    float           VOXEL_SIZE_X, VOXEL_SIZE_Y, VOXEL_SIZE_Z, TR;

	size_t          NUMBER_OF_GLM_REGRESSORS, NUMBER_OF_TOTAL_GLM_REGRESSORS;
    size_t          NUMBER_OF_DETRENDING_REGRESSORS = 4;
    size_t          NUMBER_OF_MOTION_REGRESSORS = 6;	

	int				NUMBER_OF_EVENTS;
	size_t			HIGHRES_FACTOR = 100;

	int				p, q;

    //-----------------------
    // Output

	float			*h_AccPr;    
    float           *h_PPM_Volumes, *h_Beta_Volumes, *h_Beta_Posterior, *h_IBeta_Volumes;
    float           *h_Gamma_Volumes, *h_Gamma_Posterior, *h_IGamma_Volumes;
    float           *h_Rho_Volumes, *h_Rho_Posterior, *h_IRho_Volumes;
	float			*h_Contrast_Volumes, *h_Residuals, *h_Residual_Variances, *h_Statistical_Maps;        
	float           *h_Design_Matrix, *h_Design_Matrix2;
	float			*h_NonStationary_Draws;

	//--------------

    void*           allMemoryPointers[500];
	for (int i = 0; i < 500; i++)
	{
		allMemoryPointers[i] = NULL;
	}
    
	nifti_image*	allNiftiImages[500];
	for (int i = 0; i < 500; i++)
	{
		allNiftiImages[i] = NULL;
	}

    int             numberOfMemoryPointers = 0;
	int				numberOfNiftiImages = 0;

	size_t			allocatedHostMemory = 0;

	//--------------
    
    // Default parameters
        
    bool            DEBUG = false;
    bool            PRINT = true;
	bool			VERBOS = false;
   	bool			CHANGE_OUTPUT_FILENAME = false;    

	int				AR_ORDER = 4;
	int				MCMC_ITERATIONS = 1000;
	int				BURNIN_ITERATIONS = 500;
	double			prcBurnin = (double)BURNIN_ITERATIONS / (double)MCMC_ITERATIONS;

	// Default
	double			tauBeta = 10.0;
	double			tauGamma = 10.0;
	double			tauRho = 1.0;
	double			iota = 1.0;
	double			r = 0.5;

	bool			forceStationarity = false;
	int				nStepsGamma = 2;
	int				hessMethodGamma = 0;
	int 			propDfGamma = 10;
	double			IUpdatePrGamma = 0.6;
	int				linkType = 1;
                   
	int				nOnTrialBeta = 0;
	int				nOnTrialGamma = 0;
	int				nOnTrialRho = 0;

	size_t			NUMBER_OF_CONTRASTS = 1; 
	int	 			NUMBER_OF_STATISTICAL_MAPS = 1;
    float           CLUSTER_DEFINING_THRESHOLD = 2.5f;
	int				STATISTICAL_TEST = 0;
	int				INFERENCE_MODE = 1;
	bool			MASK = false;
	int             USE_TEMPORAL_DERIVATIVES = 0;

	bool REGRESS_GLOBALMEAN = false;
	bool FOUND_DESIGN = false;
	bool FOUND_CONTRASTS = false;
	bool BETAS_ONLY = false;
	bool CONTRASTS_ONLY = false;
	bool RAW_DESIGNMATRIX = false;
	bool REGRESS_MOTION = false;
	bool REGRESS_MOTION_DERIV = false;
	bool REGRESS_GLOBAL_MEAN = false;
	bool REGRESS_CONFOUNDS = false;
	bool RAW_REGRESSORS = false;

	bool GAMMA_REGRESSORS = false;
	bool ONTRIAL_BETA = false;
	bool ONTRIAL_GAMMA = false;
	bool ONTRIAL_RHO = false;
	bool SAVE_FULL_POSTERIOR = false;

	bool updateInclusion = false;

	bool			MULTIPLE_RUNS = false;    
					NUMBER_OF_RUNS = 1;

	bool FIRST_LEVEL = false;
	bool SECOND_LEVEL = false;

	bool WRITE_DESIGNMATRIX = false;
	bool WRITE_ORIGINAL_DESIGNMATRIX = false;

	const char*		MASK_NAME;
	const char*		DESIGN_FILE;        
	const char*		CONTRASTS_FILE;
	const char*		MOTION_PARAMETERS_FILE;
	const char*		MOTION_DERIV_PARAMETERS_FILE;
	const char*		outputFilename;

	const char*		GAMMA_REGRESSORS_NAME;       
	const char*		ONTRIAL_BETA_NAME;
	const char*		ONTRIAL_GAMMA_NAME;
	const char*		ONTRIAL_RHO_NAME;

    //---------------------    
    
    /* Input arguments */
    FILE *fp = NULL; 
    
    // No inputs, so print help text
    if (argc == 1)
    {   
		printf("\nThe function applies the Hetero GLM for single subject analysis.\n\n");     
        printf("Usage first level, design.txt contains all regressors:\n\n");
        printf("HeteroGLM volumes.nii -design design.txt -contrasts contrasts.txt [options]\n\n");
        printf("Usage first level, three runs, design1.txt contains all regressors for run 1:\n\n");
        printf("HeteroGLM -runs 3 volumes1.nii volumes2.nii volumes3.nii -design design1.txt design2.txt design3.txt -contrasts contrasts.txt [options]\n\n");
        printf("Usage first level, regressors.txt contains a file name to each regressor:\n\n");
        printf("HeteroGLM volumes.nii -designfiles regressors.txt -contrasts contrasts.txt [options]\n\n");
        printf("Usage first level, regressors.txt contains a file name to each raw regressor:\n\n");
        printf("HeteroGLM volumes.nii -designfiles regressors.txt -contrasts contrasts.txt -rawregressors [options]\n\n");
        printf("Usage first level, three runs, regressors1.txt contains a file name to each regressor for run 1:\n\n");
        printf("HeteroGLM -runs 3 volumes1.nii volumes2.nii volumes3.nii -designfiles regressors1.txt regressors2.txt regressors3.txt -contrasts contrasts.txt [options]\n\n");
        printf("Options:\n\n");
        printf(" -design                    The design matrix to use, one row per volume, one column per regressor \n");
        printf(" -designfiles               File containing regressor files to use to create design matrix\n");
        printf(" -contrasts                 The contrast vector(s) to apply to the estimated beta values \n\n");
        printf(" -runs                      Number of runs \n");
        printf(" -rawregressors             Use raw regressor files (one per regressor)\n");
        printf(" -detrendingregressors      Set the number of detrending regressors, 1-4 (default 4) \n");
        printf(" -temporalderivatives       Use temporal derivatives for the activity regressors (default no) \n");
        printf(" -regressmotion             Provide file with motion regressors to use in design matrix (default no) \n");
        printf(" -regressmotionderiv        Provide file with derivative of motion regressors to use in design matrix (default no) \n");
        printf(" -regressglobalmean         Include global mean in design matrix (default no) \n\n");

        printf(" -arorder                   Order of auto regressive noise model in each voxel (default 4)\n");
        printf(" -forcestationarity         Force stationary AR process (default false)\n");
        printf(" -taurho                    tau rho parameter for prior variance of rho (default 1.0)\n\n");
        printf(" -iota                      Iota parameter for prior variance of rho (default 1.0)\n\n");

        printf(" -gammaregressors           Text file with regressors for variance (default all regressors, counting from 0) \n");
        printf(" -ontrialbeta               Text file with regressors on trial for beta (default none, counting from 0) \n");
        printf(" -ontrialgamma              Text file with regressors on trial for gamma (default none, counting from 0) \n");
        printf(" -ontrialrho                Text file with regressors on trial for rho (default none, counting from 0) \n");
        printf(" -draws                     Number of draws in the Gibbs sampling (default 1000) \n");
        printf(" -burnin                    Number of burnin draws in the Gibbs sampling that are discarded (default 500) \n");
        printf(" -savefullposterior         Save full posterior for beta, for activity regressors (default no) \n");
        printf(" -updateinclusion           Update inclusion probabilities (default no) \n");
        printf(" \n");
        printf(" \nMisc options \n\n");
        printf(" -mask                      A mask that defines which voxels to run the GLM for (default none) \n");
		printf(" -saveoriginaldesignmatrix  Save the original design matrix used (default no) \n");
        printf(" -savedesignmatrix          Save the total design matrix used (default no) \n");        
		printf(" -output                    Set output filename (default volumes_) \n");
        printf(" -quiet                     Don't print anything to the terminal (default false) \n");
        printf(" -verbose                   Print extra stuff (default false) \n");
        printf("\n\n");
        
        return EXIT_SUCCESS;
    }

	// Check if first argument is -runs
	char *temp = argv[1];
    int i = 2;
    if (strcmp(temp,"-runs") == 0)
    {		
        MULTIPLE_RUNS = true;

		char *p;
		NUMBER_OF_RUNS = (int)strtol(argv[2], &p, 10);
			
		if (!isspace(*p) && *p != 0)
	    {
	        printf("Number of runs must be an integer! You provided %s \n",argv[2]);
			return EXIT_FAILURE;
	    }
        else if (NUMBER_OF_RUNS < 1)
        {
			printf("Number of runs must be > 1!\n");
            return EXIT_FAILURE;
        }
		i = 3 + NUMBER_OF_RUNS;
    }

	DATA_T_PER_RUN = (size_t*)malloc(NUMBER_OF_RUNS*sizeof(size_t));

	/*
    // Try to open file
    else if (argc > 1)
    {        
        fp = fopen(argv[1],"r");
        if (fp == NULL)
        {
            printf("Could not open file %s !\n",argv[1]);
            return EXIT_FAILURE;
        }
        fclose(fp);             
    }
	*/
    
    // Loop over additional inputs

    while (i < argc)
    {
        char *input = argv[i];
        char *p;

        if (strcmp(input,"-design") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -design !\n");
                return EXIT_FAILURE;
			}

            DESIGN_FILE = argv[i+1];
			RAW_DESIGNMATRIX = true;
			FOUND_DESIGN = true;
            i += 1 + NUMBER_OF_RUNS;
        }
        else if (strcmp(input,"-designfiles") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -design !\n");
                return EXIT_FAILURE;
			}

            DESIGN_FILE = argv[i+1];
			RAW_DESIGNMATRIX = false;
			FOUND_DESIGN = true;
            i += 1 + NUMBER_OF_RUNS;
        }
        else if (strcmp(input,"-contrasts") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -contrasts !\n");
                return EXIT_FAILURE;
			}

            CONTRASTS_FILE = argv[i+1];
			FOUND_CONTRASTS = true;
            i += 2;
        }
        else if (strcmp(input,"-rawregressors") == 0)
        {
            RAW_REGRESSORS = true;
            i += 1;
        }
        else if (strcmp(input,"-detrendingregressors") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -detrendingregressors !\n");
                return EXIT_FAILURE;
			}

            NUMBER_OF_DETRENDING_REGRESSORS = (int)strtol(argv[i+1], &p, 10);

			if (!isspace(*p) && *p != 0)
		    {
		        printf("Number of detrending regressors must be an integer! You provided %s \n",argv[i+1]);
				return EXIT_FAILURE;
		    }
            if ((NUMBER_OF_DETRENDING_REGRESSORS < 1) || (NUMBER_OF_DETRENDING_REGRESSORS > 4))
            {
                printf("Number of detrending regressors must be >= 1 & <= 4!\n");
                return EXIT_FAILURE;
            }
            i += 2;
        }


        else if (strcmp(input,"-temporalderivatives") == 0)
        {
            USE_TEMPORAL_DERIVATIVES = 1;
            i += 1;
        }

        else if (strcmp(input,"-regressmotion") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -regressmotion!\n");
                return EXIT_FAILURE;
			}

            MOTION_PARAMETERS_FILE = argv[i+1];
            REGRESS_MOTION = true;
            i += 2;
        }
        else if (strcmp(input,"-regressmotionderiv") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -regressmotionderiv!\n");
                return EXIT_FAILURE;
			}

            MOTION_DERIV_PARAMETERS_FILE = argv[i+1];
            REGRESS_MOTION_DERIV = true;
            i += 2;
        }

        else if (strcmp(input,"-regressglobalmean") == 0)
        {
            REGRESS_GLOBALMEAN = 1;
            i += 1;
        }



        else if (strcmp(input,"-arorder") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -arorder !\n");
                return EXIT_FAILURE;
			}

            AR_ORDER = (int)strtol(argv[i+1], &p, 10);

			if (!isspace(*p) && *p != 0)
		    {
		        printf("AR order must be an integer! You provided %s \n",argv[i+1]);
				return EXIT_FAILURE;
		    }
            if (AR_ORDER <= 0)
            {
                printf("AR order must be > 0!\n");
                return EXIT_FAILURE;
            }
            i += 2;
        }
		else if (strcmp(input,"-forcestationarity") == 0)
        {
            forceStationarity = true;
            i += 1;
        }
		else if (strcmp(input,"-taurho") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -taurho !\n");
                return EXIT_FAILURE;
			}
            
            tauRho = (double)strtod(argv[i+1], &p);
            
			if (!isspace(*p) && *p != 0)
		    {
		        printf("tau rho parameter must be a double! You provided %s \n",argv[i+1]);
				return EXIT_FAILURE;
		    }
  			else if ( tauRho <= 0.0 )
            {
                printf("tau rho must be > 0.0 !\n");
                return EXIT_FAILURE;
            }
            i += 2;
		}

		else if (strcmp(input,"-iota") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -iota !\n");
                return EXIT_FAILURE;
			}
            
            iota = (double)strtod(argv[i+1], &p);
            
			if (!isspace(*p) && *p != 0)
		    {
		        printf("iota parameter must be a double! You provided %s \n",argv[i+1]);
				return EXIT_FAILURE;
		    }
  			else if ( iota <= 0.0 )
            {
                printf("iota must be > 0.0 !\n");
                return EXIT_FAILURE;
            }
            i += 2;
		}


		else if (strcmp(input,"-gammaregressors") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read name after -gammaregressors !\n");
                return EXIT_FAILURE;
			}

			GAMMA_REGRESSORS = true;
            GAMMA_REGRESSORS_NAME = argv[i+1];
            i += 2;
        }
		else if (strcmp(input,"-ontrialbeta") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read name after -ontrialbeta !\n");
                return EXIT_FAILURE;
			}

			ONTRIAL_BETA = true;
            ONTRIAL_BETA_NAME = argv[i+1];
            i += 2;
        }
		else if (strcmp(input,"-ontrialgamma") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read name after -ontrialgamma !\n");
                return EXIT_FAILURE;
			}

			ONTRIAL_GAMMA = true;
            ONTRIAL_GAMMA_NAME = argv[i+1];
            i += 2;
        }
		else if (strcmp(input,"-ontrialrho") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read name after -ontrialrho !\n");
                return EXIT_FAILURE;
			}

			ONTRIAL_RHO = true;
            ONTRIAL_RHO_NAME = argv[i+1];
            i += 2;
        }
        else if (strcmp(input,"-draws") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -draws !\n");
                return EXIT_FAILURE;
			}

            MCMC_ITERATIONS = (int)strtol(argv[i+1], &p, 10);

			if (!isspace(*p) && *p != 0)
		    {
		        printf("Number of draws must be an integer! You provided %s \n",argv[i+1]);
				return EXIT_FAILURE;
		    }
            if (MCMC_ITERATIONS <= 0)
            {
                printf("Number of draws must be > 0!\n");
                return EXIT_FAILURE;
            }
            i += 2;
        }
        else if (strcmp(input,"-burnin") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read value after -burnin !\n");
                return EXIT_FAILURE;
			}

            BURNIN_ITERATIONS = (int)strtol(argv[i+1], &p, 10);

			if (!isspace(*p) && *p != 0)
		    {
		        printf("Number of burnin draws must be an integer! You provided %s \n",argv[i+1]);
				return EXIT_FAILURE;
		    }
            if (BURNIN_ITERATIONS <= 0)
            {
                printf("Number of burnin draws must be > 0!\n");
                return EXIT_FAILURE;
            }
            i += 2;
        }
        else if (strcmp(input,"-savefullposterior") == 0)
        {
            SAVE_FULL_POSTERIOR = true;
            i += 1;
        }
        else if (strcmp(input,"-updateinclusion") == 0)
        {
            updateInclusion = true;
            i += 1;
        }


		else if (strcmp(input,"-mask") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read name after -mask !\n");
                return EXIT_FAILURE;
			}

			MASK = true;
            MASK_NAME = argv[i+1];
            i += 2;
        }
        else if (strcmp(input,"-savedesignmatrix") == 0)
        {
            WRITE_DESIGNMATRIX = true;
            i += 1;
        }
        else if (strcmp(input,"-saveoriginaldesignmatrix") == 0)
        {
            WRITE_ORIGINAL_DESIGNMATRIX = true;
            i += 1;
        }
        else if (strcmp(input,"-debug") == 0)
        {
            DEBUG = true;
            i += 1;
        }
        else if (strcmp(input,"-quiet") == 0)
        {
            PRINT = false;
            i += 1;
        }
        else if (strcmp(input,"-verbose") == 0)
        {
            VERBOS = true;
            i += 1;
        }
        else if (strcmp(input,"-output") == 0)
        {
			if ( (i+1) >= argc  )
			{
			    printf("Unable to read name after -output !\n");
                return EXIT_FAILURE;
			}

			CHANGE_OUTPUT_FILENAME = true;
            outputFilename = argv[i+1];
            i += 2;
        }
        else
        {
            printf("Unrecognized option! %s \n",argv[i]);
            return EXIT_FAILURE;
        }                
    }

	if (!FOUND_DESIGN)
	{
    	printf("No design file detected, aborting! \n");
        return EXIT_FAILURE;
	}

	if (!FOUND_CONTRASTS)
	{
    	printf("No contrasts file detected, aborting! \n");
        return EXIT_FAILURE;
	}

	prcBurnin = (double)BURNIN_ITERATIONS / (double)MCMC_ITERATIONS;

	//------------------------------------------
	// Check for invalid combinations

	if (RAW_DESIGNMATRIX && USE_TEMPORAL_DERIVATIVES)
	{
    	printf("Cannot use temporal derivatives for raw design matrix, aborting! \n");
        return EXIT_FAILURE;
	}

	if (RAW_REGRESSORS && USE_TEMPORAL_DERIVATIVES)
	{
    	printf("Cannot use temporal derivatives for raw regressors, aborting! \n");
        return EXIT_FAILURE;
	}

	// Check if BROCCOLI_DIR variable is set
	if (getenv("BROCCOLI_DIR") == NULL)
	{
        printf("The environment variable BROCCOLI_DIR is not set!\n");
        return EXIT_FAILURE;
	}

	//------------------------------------------
    // Read number of regressors from design matrix file
  	//------------------------------------------

	std::ifstream design;
    std::string tempString;
    int tempNumber;
    std::string NR("NumRegressors");

    design.open(DESIGN_FILE);    
    if (!design.good())
    {
        design.close();
        printf("Unable to open design file %s. Aborting! \n",DESIGN_FILE);
        return EXIT_FAILURE;
    }
    
    // Get number of regressors
    design >> tempString; // NumRegressors as string
    if (tempString.compare(NR) != 0)
    {
        design.close();
        printf("First element of the design file %s should be the string 'NumRegressors', but it is %s. Aborting! \n",DESIGN_FILE,tempString.c_str());
        return EXIT_FAILURE;
    }
    design >> NUMBER_OF_GLM_REGRESSORS;
    
    if (NUMBER_OF_GLM_REGRESSORS <= 0)
    {
        design.close();
        printf("Number of regressors must be > 0 ! You provided %zu regressors in the design file %s. Aborting! \n",NUMBER_OF_GLM_REGRESSORS,DESIGN_FILE);
        return EXIT_FAILURE;
    }
    design.close();

	//------------------------------------------  
    // Read number of contrasts from contrasts file
	//------------------------------------------

   	std::ifstream contrasts;    

    contrasts.open(CONTRASTS_FILE);  
    if (!contrasts.good())
    {
        contrasts.close();
        printf("Unable to open contrasts file %s. Aborting! \n",CONTRASTS_FILE);
        return EXIT_FAILURE;
    }
    
    contrasts >> tempString; // NumRegressors as string
    if (tempString.compare(NR) != 0)
    {
        contrasts.close();
        printf("First element of the contrasts file should be the string 'NumRegressors', but it is %s. Aborting! \n",tempString.c_str());
        return EXIT_FAILURE;
    }
    contrasts >> tempNumber;
    
    // Check for consistency
    if ( tempNumber != NUMBER_OF_GLM_REGRESSORS )
	{
        contrasts.close();
        printf("Design file says that number of regressors is %zu, while contrast file says there are %i regressors. Aborting! \n",NUMBER_OF_GLM_REGRESSORS,tempNumber);
        return EXIT_FAILURE;
    }
    
    contrasts >> tempString; // NumContrasts as string
    std::string NC("NumContrasts");
    if (tempString.compare(NC) != 0)
    {
        contrasts.close();
        printf("Third element of the contrasts file should be the string 'NumContrasts', but it is %s. Aborting! \n",tempString.c_str());
        return EXIT_FAILURE;
    }
    contrasts >> NUMBER_OF_CONTRASTS;
			
    if (NUMBER_OF_CONTRASTS <= 0)
    {
        contrasts.close();
	    printf("Number of contrasts must be > 0 ! You provided %zu in the contrasts file. Aborting! \n",NUMBER_OF_CONTRASTS);
        return EXIT_FAILURE;
    }
    contrasts.close();
  
	//------------------------------------------
    // Read data
	//------------------------------------------

	double startTime = GetWallTime();

	nifti_image *inputData;
	std::vector<nifti_image*> allfMRINiftiImages;

	if (!MULTIPLE_RUNS)
	{
		inputData = nifti_image_read(argv[1],1);
	    allfMRINiftiImages.push_back(inputData);

    	if (inputData == NULL)
    	{
    	    printf("Could not open fMRI data file %s!\n",argv[1]);
    	    return EXIT_FAILURE;
    	}
		allNiftiImages[numberOfNiftiImages] = inputData;
		numberOfNiftiImages++;
	}
	else
	{
		for (int i = 0; i < NUMBER_OF_RUNS; i++)
		{
			inputData = nifti_image_read(argv[3+i],1);
			allfMRINiftiImages.push_back(inputData);    

    		if (inputData == NULL)
    		{
    		    printf("Could not open fMRI data file %s for run %i !\n",argv[3+i],i+1);
    		    return EXIT_FAILURE;
    		}
			allNiftiImages[numberOfNiftiImages] = inputData;
			numberOfNiftiImages++;
		    DATA_T_PER_RUN[i] = inputData->nt;    
		}
	}

	nifti_image *inputMask;
	if (MASK)
	{
	    inputMask = nifti_image_read(MASK_NAME,1);
    
	    if (inputMask == NULL)
	    {
        	printf("Could not open mask volume!\n");
	        return EXIT_FAILURE;
	    }
		allNiftiImages[numberOfNiftiImages] = inputMask;
		numberOfNiftiImages++;
	}
	else
	{
       	printf("\nWarning: No mask being used, doing GLM for all voxels.\n\n");
	}
    	
	double endTime = GetWallTime();

	if (VERBOS)
 	{
		printf("It took %f seconds to read the nifti file(s)\n",(float)(endTime - startTime));
	}

    // Get data dimensions from input data
   	DATA_W = inputData->nx;
    DATA_H = inputData->ny;
    DATA_D = inputData->nz;   

	if (!MULTIPLE_RUNS)
	{   
	    DATA_T = inputData->nt;   
		DATA_T_PER_RUN[0] = DATA_T; 
	}
	else
	{
		DATA_T = 0;
		for (int i = 0; i < NUMBER_OF_RUNS; i++)
		{
		    DATA_T += DATA_T_PER_RUN[i];
		}
	}

    // Get voxel sizes
    VOXEL_SIZE_X = inputData->dx;
    VOXEL_SIZE_Y = inputData->dy;
    VOXEL_SIZE_Z = inputData->dz;

	// Get fMRI repetition time
    TR = inputData->dt;

	// Check if there is more than one volume
	if (DATA_T <= 1)
	{
		printf("Input data is a single volume, cannot run GLM! \n");
		FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
		return EXIT_FAILURE;	
	}

	// Check number of regressors
	if (NUMBER_OF_GLM_REGRESSORS > DATA_T)
	{
		printf("More regressor than data points, cannot run GLM! \n");
		FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
		return EXIT_FAILURE;	
	}
    	
	// Check if mask volume has the same dimensions as the data
	if (MASK)
	{
		size_t TEMP_DATA_W = inputMask->nx;
		size_t TEMP_DATA_H = inputMask->ny;
		size_t TEMP_DATA_D = inputMask->nz;

		if ( (TEMP_DATA_W != DATA_W) || (TEMP_DATA_H != DATA_H) || (TEMP_DATA_D != DATA_D) )
		{
			printf("Input data has the dimensions %zu x %zu x %zu, while the mask volume has the dimensions %zu x %zu x %zu. Aborting! \n",DATA_W,DATA_H,DATA_D,TEMP_DATA_W,TEMP_DATA_H,TEMP_DATA_D);
			FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
			return EXIT_FAILURE;
		}
	}
    

	
    // ------------------------------------------------  
	// Allocate memory
    // ------------------------------------------------  

	startTime = GetWallTime();

    size_t STATISTICAL_MAPS_SIZE, CONTRAST_SCALAR_SIZE;
	CONTRAST_SCALAR_SIZE = NUMBER_OF_CONTRASTS * sizeof(float);
	STATISTICAL_MAPS_SIZE = DATA_W * DATA_H * DATA_D * NUMBER_OF_CONTRASTS * sizeof(float);

	NUMBER_OF_TOTAL_GLM_REGRESSORS = NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + REGRESS_GLOBALMEAN + NUMBER_OF_MOTION_REGRESSORS * REGRESS_MOTION + NUMBER_OF_MOTION_REGRESSORS * REGRESS_MOTION_DERIV;
	p = NUMBER_OF_TOTAL_GLM_REGRESSORS;
	q = NUMBER_OF_TOTAL_GLM_REGRESSORS;

    // ------------------------------------------------

    // Print some info
    if (PRINT)
    {
        printf("Authored by K.A. Eklund \n");
		if (!MULTIPLE_RUNS)
		{
		    printf("Data size: %zu x %zu x %zu x %zu \n", DATA_W, DATA_H, DATA_D, DATA_T);
		}
		else
		{
			for (int i = 0; i < NUMBER_OF_RUNS; i++)
			{
			    printf("Data size for run %i: %zu x %zu x %zu x %zu \n", i+1, DATA_W, DATA_H, DATA_D, DATA_T_PER_RUN[i]);
			}
		    printf("Total data size: %zu x %zu x %zu x %zu \n", DATA_W, DATA_H, DATA_D, DATA_T);
		}
		printf("Number of original regressors: %zu \n",  NUMBER_OF_GLM_REGRESSORS);
		printf("Number of total regressors: %zu \n",  NUMBER_OF_TOTAL_GLM_REGRESSORS);	
        printf("Number of contrasts: %zu \n",  NUMBER_OF_CONTRASTS);
        printf("Number of burnin draws: %i, number of saved draws: %i, ratio: %f \n",  BURNIN_ITERATIONS,MCMC_ITERATIONS,prcBurnin);
    }

    // ------------------------------------------------

    size_t DATA_SIZE = DATA_W * DATA_H * DATA_D * DATA_T * sizeof(float);
    size_t VOLUME_SIZE = DATA_W * DATA_H * DATA_D * sizeof(float);
      
    size_t GLM_SIZE = DATA_T * NUMBER_OF_GLM_REGRESSORS * sizeof(float);
    size_t CONTRAST_SIZE = NUMBER_OF_GLM_REGRESSORS * NUMBER_OF_CONTRASTS * sizeof(float);
    size_t DESIGN_MATRIX_SIZE = NUMBER_OF_TOTAL_GLM_REGRESSORS * DATA_T * sizeof(float);
	size_t HIGHRES_REGRESSORS_SIZE = DATA_T * HIGHRES_FACTOR * NUMBER_OF_GLM_REGRESSORS * sizeof(float);    
	size_t PPM_DATA_SIZE = DATA_W * DATA_H * DATA_D * NUMBER_OF_GLM_REGRESSORS * sizeof(float);
    size_t BETA_DATA_SIZE = DATA_W * DATA_H * DATA_D * NUMBER_OF_TOTAL_GLM_REGRESSORS * sizeof(float);
    size_t GAMMA_DATA_SIZE = DATA_W * DATA_H * DATA_D * NUMBER_OF_TOTAL_GLM_REGRESSORS * sizeof(float);
    size_t RHO_DATA_SIZE = DATA_W * DATA_H * DATA_D * AR_ORDER * sizeof(float);
    size_t RESIDUALS_DATA_SIZE = DATA_W * DATA_H * DATA_D * DATA_T * sizeof(float);
    size_t MOTION_PARAMETERS_SIZE = NUMBER_OF_MOTION_REGRESSORS * DATA_T * sizeof(float);
   
    size_t BETA_FULL_DATA_SIZE = DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS * NUMBER_OF_TOTAL_GLM_REGRESSORS * sizeof(float);
    size_t RHO_FULL_DATA_SIZE = DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS * AR_ORDER * sizeof(float);

	if (!MULTIPLE_RUNS)
	{
		// If the data is in float format, we can just copy the pointer
		if ( inputData->datatype != DT_FLOAT )
		{
			AllocateMemory(h_Data, DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "fMRI_VOLUMES");
		}
		else
		{
			allocatedHostMemory += DATA_SIZE;
		}
	}
	else
	{
		AllocateMemory(h_Data, DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "fMRI_VOLUMES");
	}

    AllocateMemory(h_X_GLM, GLM_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "DESIGN_MATRIX");
    AllocateMemory(h_Highres_Regressors, HIGHRES_REGRESSORS_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "HIGHRES_REGRESSOR");
    AllocateMemory(h_LowpassFiltered_Regressors, HIGHRES_REGRESSORS_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "LOWPASSFILTERED_REGRESSOR");
    AllocateMemory(h_xtxxt_GLM, GLM_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "DESIGN_MATRIX_INVERSE");
    AllocateMemory(h_Contrasts, CONTRAST_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "CONTRASTS");
    AllocateMemory(h_ctxtxc_GLM, CONTRAST_SCALAR_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "CONTRAST_SCALARS");
	AllocateMemory(h_Design_Matrix, DESIGN_MATRIX_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "TOTAL_DESIGN_MATRIX");

	AllocateMemory(h_NonStationary_Draws, VOLUME_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "NONSTATIONARY_DRAWS");

	AllocateMemory(h_Mask, VOLUME_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "MASK");

	if (REGRESS_MOTION)
	{
		AllocateMemory(h_Motion_Parameters, MOTION_PARAMETERS_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "MOTION_PARAMETERS");       
	}
	if (REGRESS_MOTION_DERIV)
	{
		AllocateMemory(h_Motion_Deriv_Parameters, MOTION_PARAMETERS_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "MOTION_DERIV_PARAMETERS");       
	}

    AllocateMemory(h_PPM_Volumes, PPM_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "PPM_VOLUMES");

    AllocateMemory(h_AccPr, VOLUME_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "ACC_PR");

    AllocateMemory(h_Beta_Volumes, BETA_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "BETA_VOLUMES");
    AllocateMemory(h_IBeta_Volumes, BETA_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "IBETA_VOLUMES");

    AllocateMemory(h_Gamma_Volumes, GAMMA_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "GAMMA_VOLUMES");
    AllocateMemory(h_IGamma_Volumes, GAMMA_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "IGAMMA_VOLUMES");

    AllocateMemory(h_Rho_Volumes, RHO_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "RHO_VOLUMES");
    AllocateMemory(h_IRho_Volumes, RHO_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "IRHO_VOLUMES");

	if (SAVE_FULL_POSTERIOR)
	{
	    AllocateMemory(h_Beta_Posterior, BETA_FULL_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "BETA_POSTERIOR");
	    AllocateMemory(h_Gamma_Posterior, BETA_FULL_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "GAMMA_POSTERIOR");
	    AllocateMemory(h_Rho_Posterior, RHO_FULL_DATA_SIZE, allMemoryPointers, numberOfMemoryPointers, allNiftiImages, numberOfNiftiImages, allocatedHostMemory, "RHO_POSTERIOR");
	}

	endTime = GetWallTime();
    
	if (VERBOS)
 	{
		printf("It took %f seconds to allocate memory\n",(float)(endTime - startTime));
	}
   
    // ------------------------------------------------
	// Read events for each regressor    	

	if (RAW_DESIGNMATRIX)
	{
		int designfile;
		if (!MULTIPLE_RUNS)
		{
			designfile = 3;
		}
		else
		{
			designfile = 4 + NUMBER_OF_RUNS;
		}

		int accumulatedTRs = 0;
		for (int run = 0; run < NUMBER_OF_RUNS; run++)
		{
		    // Open design file again
		    design.open(argv[designfile]);

		    if (!design.good())
		    {
		        design.close();
		        printf("Unable to open design file %s. Aborting! \n",argv[designfile]);
				FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
		        FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
		        return EXIT_FAILURE;
		    }

		    // Read first two values again
		    design >> tempString; // NumRegressors as string
		    design >> NUMBER_OF_GLM_REGRESSORS;

			float tempFloat;	
			for (size_t t = 0; t < DATA_T_PER_RUN[run]; t++)
			{
				for (size_t r = 0; r < NUMBER_OF_GLM_REGRESSORS; r++)
				{
					if (! (design >> h_X_GLM[t + accumulatedTRs + r * DATA_T]) )
					{
						design.close();
				        printf("Could not read all values of the design file %s, aborting! Stopped reading at time point %zu for regressor %zu. Please check if the number of regressors and time points are correct. \n",argv[designfile],t,r);      
				        FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
				        FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
				        return EXIT_FAILURE;
					}
				}
			}	
			design.close();
			designfile++;
			accumulatedTRs += DATA_T_PER_RUN[run];
		}
	}
	else
	{
		startTime = GetWallTime();

	    // Each line of the design file is a filename

		int designfile;
		if (!MULTIPLE_RUNS)
		{
			designfile = 3;
		}
		else
		{
			designfile = 4 + NUMBER_OF_RUNS;
		}

		// Reset highres regressors
	    for (size_t t = 0; t < NUMBER_OF_GLM_REGRESSORS * DATA_T * HIGHRES_FACTOR; t++)
    	{
			h_Highres_Regressors[t] = 0.0f;
		}

		int accumulatedTRs = 0;
		for (int run = 0; run < NUMBER_OF_RUNS; run++)
		{  
		    // Open design file again
		    design.open(argv[designfile]);

		    if (!design.good())
		    {
		        design.close();
		        printf("Unable to open design file %s. Aborting! \n",argv[designfile]);
				FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
		        FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
		        return EXIT_FAILURE;
		    }

		    // Read first two values again
		    design >> tempString; // NumRegressors as string
		    design >> NUMBER_OF_GLM_REGRESSORS;

			if (!RAW_REGRESSORS)
			{    
			    // Loop over the number of regressors provided in the design file
			    for (size_t r = 0; r < NUMBER_OF_GLM_REGRESSORS; r++)
		    	{
			        // Each regressor is a filename, so try to open the file
			        std::ifstream regressor;
			        std::string filename;
			        design >> filename;
			        regressor.open(filename.c_str());
			        if (!regressor.good())
			        {
			            regressor.close();
			            printf("Unable to open the regressor file %s . Aborting! \n",filename.c_str());
			            FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
			            FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
			            return EXIT_FAILURE;
			        }
        
			        // Read number of events for current regressor
			        regressor >> tempString; // NumEvents as string
			        std::string NE("NumEvents");
			        if (tempString.compare(NE) != 0)
			        {
    			        design.close();
			            printf("First element of each regressor file should be the string 'NumEvents', but it is %s for the regressor file %s. Aborting! \n",tempString.c_str(),filename.c_str());
			            FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
			            FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
			            return EXIT_FAILURE;
    			    }
			        regressor >> NUMBER_OF_EVENTS;
	
					if (DEBUG)
					{
						printf("Number of events for regressor %zu is %i \n",r,NUMBER_OF_EVENTS);
					}
    	    
    			    // Loop over events
    			    for (int e = 0; e < NUMBER_OF_EVENTS; e++)
    			    {
    			        float onset;
    			        float duration;
    			        float value;
    	        
    			        // Read onset, duration and value for current event
						if (! (regressor >> onset) )
						{
							regressor.close();
    			            design.close();
    			            printf("Unable to read the onset for event %i in regressor file %s, aborting! Check the regressor file. \n",e,filename.c_str());
    			            FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
    			            FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    			            return EXIT_FAILURE;
						}
	
    			        if (! (regressor >> duration) )
						{
							regressor.close();
    			            design.close();
    			            printf("Unable to read the duration for event %i in regressor file %s, aborting! Check the regressor file. \n",e,filename.c_str());
    			            FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
    			            FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    			            return EXIT_FAILURE;
						}

						if (! (regressor >> value) )
						{
							regressor.close();
    			            design.close();
    			            printf("Unable to read the value for event %i in regressor file %s, aborting! Check the regressor file. \n",e,filename.c_str());
    			            FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
    			            FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    			            return EXIT_FAILURE;
						}
    	    
						if (DEBUG)
						{
							printf("Event %i: Onset is %f, duration is %f and value is %f \n",e,onset,duration,value);
						}
    	        
    			        int start = (int)round(onset * (float)HIGHRES_FACTOR / TR);
    			        int activityLength = (int)round(duration * (float)HIGHRES_FACTOR / TR);
    	        
						if (DEBUG)
						{
							printf("Event %i: Start is %i, activity length is %i \n",e,start,activityLength);
						}

    			        // Put values into highres GLM
    			        for (size_t i = 0; i < activityLength; i++)
    			        {
    			            if ((start + i) < (DATA_T_PER_RUN[run] * HIGHRES_FACTOR) )
    			            {
    			                h_Highres_Regressors[start + i + accumulatedTRs*HIGHRES_FACTOR + r * DATA_T*HIGHRES_FACTOR] = value;
    			            }
    			            else
    			            {
								regressor.close();
    			                design.close();
    			                printf("For run %i, the activity start or duration for event %i in regressor file %s is longer than the duration of the fMRI data, aborting! Check the regressor file .\n",run+1,e,filename.c_str());	
    			                FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
    			                FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    			                return EXIT_FAILURE;
    			            }
    			        }            
    			    }
					regressor.close();	
				}
			}
			else if (RAW_REGRESSORS)
			{
				// Loop over the number of regressors provided in the design file
			    for (size_t r = 0; r < NUMBER_OF_GLM_REGRESSORS; r++)
    			{
			        // Each regressor is a filename, so try to open the file
			        std::ifstream regressor;
			        std::string filename;
			        design >> filename;
			        regressor.open(filename.c_str());
			        if (!regressor.good())
			        {
			            regressor.close();
			            printf("Unable to open the regressor file %s . Aborting! \n",filename.c_str());
			            FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
			            FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
			            return EXIT_FAILURE;
			        }

					float value;
					int readValues = 0;
				    for (size_t t = 0; t < DATA_T_PER_RUN[run]; t++)
			    	{
						if (! (regressor >> value) )
						{
							regressor.close();
    			            design.close();
    			            printf("Unable to read the value for TR %zu in regressor file %s, aborting! Check the regressor file. \n",t,filename.c_str());
    			            FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
    			            FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    			            return EXIT_FAILURE;
						}
						h_X_GLM[t + accumulatedTRs + r * DATA_T] = value;
						readValues++;
					}
		
					// Check if number of values is the same as the number of TRs
					if (readValues != DATA_T_PER_RUN[run])
					{
						regressor.close();
    			        design.close();
    			        printf("Number of values in regressor file %s is not the same as the number of TRs in the fMRI data (%i vs %zu), aborting! Check the regressor file. \n",filename.c_str(),readValues,DATA_T_PER_RUN[run]);
    			        FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
    			        FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    			        return EXIT_FAILURE;
					}
		
					regressor.close();
				}
			}
    		design.close();
			designfile++;
			accumulatedTRs += DATA_T_PER_RUN[run];
		}
	}

	if (!RAW_REGRESSORS && !RAW_DESIGNMATRIX)
	{
		// Lowpass filter highres regressors
		LowpassFilterRegressors(h_LowpassFiltered_Regressors,h_Highres_Regressors,DATA_T,HIGHRES_FACTOR,TR,NUMBER_OF_GLM_REGRESSORS);
        
	    // Downsample highres GLM and put values into regular GLM
	    for (size_t t = 0; t < DATA_T; t++)
	    {
			for (size_t r = 0; r < NUMBER_OF_GLM_REGRESSORS; r++)
		    {	
		        h_X_GLM[t + r * DATA_T] = h_LowpassFiltered_Regressors[t*HIGHRES_FACTOR + r * DATA_T * HIGHRES_FACTOR];
			}
	    }
	}



	//------------------------------------------
	// Read the contrasts
	//------------------------------------------

	//if (!BETAS_ONLY)
	//{
		// Open contrast file again
	    contrasts.open(CONTRASTS_FILE);

	    if (!contrasts.good())
	    {
	        contrasts.close();
	        printf("Unable to open contrast file %s. Aborting! \n",CONTRASTS_FILE);
			FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
	        FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
	        return EXIT_FAILURE;
	    }

	    // Read first two values again
		contrasts >> tempString; // NumRegressors as string
	    contrasts >> tempNumber;
	    contrasts >> tempString; // NumContrasts as string
	    contrasts >> tempNumber;
   
		// Read all contrast values
		for (size_t c = 0; c < NUMBER_OF_CONTRASTS; c++)
		{
			for (size_t r = 0; r < NUMBER_OF_GLM_REGRESSORS; r++)
			{
				if (! (contrasts >> h_Contrasts[r + c * NUMBER_OF_GLM_REGRESSORS]) )
				{
				    contrasts.close();
	                printf("Unable to read all the contrast values, aborting! Check the contrasts file. \n");
	                FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
	                FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
	                return EXIT_FAILURE;
				}
			}
		}
		contrasts.close();
	//}

	Eigen::MatrixXd Contrasts(NUMBER_OF_CONTRASTS,NUMBER_OF_GLM_REGRESSORS);


    // Write original design matrix to file
	if (WRITE_ORIGINAL_DESIGNMATRIX)
	{
		std::ofstream designmatrix;

		const char* extension = "_original_designmatrix.txt";
		char* filenameWithExtension;

		CreateFilename(filenameWithExtension, inputData, extension, CHANGE_OUTPUT_FILENAME, outputFilename);
	
   		designmatrix.open(filenameWithExtension);

	    if ( designmatrix.good() )
	    {
    	    for (size_t t = 0; t < DATA_T; t++)
	        {
	    	    for (size_t r = 0; r < NUMBER_OF_GLM_REGRESSORS; r++)
		        {
            		designmatrix << std::setprecision(6) << std::fixed << (double)h_X_GLM[t + r * DATA_T] << "  ";
				}
				designmatrix << std::endl;
			}
		    designmatrix.close();
        } 	
	    else
	    {
			designmatrix.close();
	        printf("Could not open the file for writing the original design matrix!\n");
	    }
		free(filenameWithExtension);
	}
	
    //------------------------------------------
	// Read motion parameters
	//------------------------------------------

	if (REGRESS_MOTION)
	{
	    // Open motion parameters file
		std::ifstream motionparameters;
	    motionparameters.open(MOTION_PARAMETERS_FILE);  

	    if ( motionparameters.good() )
	    {
			for (size_t t = 0; t < DATA_T; t++)
			{
				for (size_t r = 0; r < NUMBER_OF_MOTION_REGRESSORS; r++)
				{
					h_Motion_Parameters[t + r * DATA_T] = 0.0f;

					if (! (motionparameters >> h_Motion_Parameters[t + r * DATA_T]) )
					{
						motionparameters.close();
				        printf("Could not read all values of the motion parameters file %s, aborting! Stopped reading at time point %zu for parameter %zu. Please check the motion parameters file\n",MOTION_PARAMETERS_FILE,t,r);      
				        FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
				        FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
				        return EXIT_FAILURE;
					}
				}
			}	
			motionparameters.close();
		}
		else
		{
			motionparameters.close();
	        printf("Could not open the motion parameters file %s !\n",MOTION_PARAMETERS_FILE);
		}
	}


	if (REGRESS_MOTION_DERIV)
	{
	    // Open motion parameters file
		std::ifstream motionparameters;
	    motionparameters.open(MOTION_DERIV_PARAMETERS_FILE);  

	    if ( motionparameters.good() )
	    {
			for (size_t t = 0; t < DATA_T; t++)
			{
				for (size_t r = 0; r < NUMBER_OF_MOTION_REGRESSORS; r++)
				{
					h_Motion_Deriv_Parameters[t + r * DATA_T] = 0.0f;

					if (! (motionparameters >> h_Motion_Deriv_Parameters[t + r * DATA_T]) )
					{
						motionparameters.close();
				        printf("Could not read all values of the motion deriv parameters file %s, aborting! Stopped reading at time point %zu for parameter %zu. Please check the motion parameters file\n",MOTION_DERIV_PARAMETERS_FILE,t,r);      
				        FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
				        FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
				        return EXIT_FAILURE;
					}
				}
			}	
			motionparameters.close();

			/*
			for (size_t t = 0; t < DATA_T; t++)
			{
				for (size_t r = 0; r < NUMBER_OF_MOTION_REGRESSORS; r++)
				{
					h_Motion_Deriv_Parameters[t + r * DATA_T] = fabs(h_Motion_Deriv_Parameters[t + r * DATA_T]);
				}
			}
			*/
		}
		else
		{
			motionparameters.close();
	        printf("Could not open the motion deriv parameters file %s !\n",MOTION_DERIV_PARAMETERS_FILE);
		}
	}



    //------------------------------------------
	// Read data
	//------------------------------------------

	startTime = GetWallTime();

	// Convert fMRI data to floats
	size_t accumulatedTRs = 0;

	if (MULTIPLE_RUNS)
	{
		for (size_t run = 0; run < NUMBER_OF_RUNS; run++)
		{
			inputData = allfMRINiftiImages[run];
	
		    if ( inputData->datatype == DT_SIGNED_SHORT )
		    {
		        short int *p = (short int*)inputData->data;
    			
		        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D * DATA_T_PER_RUN[run]; i++)
    		    {
    		        h_Data[i + accumulatedTRs * DATA_W * DATA_H * DATA_D] = (float)p[i];
    		    }
			}
		    else if ( inputData->datatype == DT_UINT8 )
		    {
		        unsigned char *p = (unsigned char*)inputData->data;
		    
		        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D * DATA_T_PER_RUN[run]; i++)
		        {
		            h_Data[i + accumulatedTRs * DATA_W * DATA_H * DATA_D] = (float)p[i];
		        }
		    }
		    else if ( inputData->datatype == DT_UINT16 )
		    {
		        unsigned short int *p = (unsigned short int*)inputData->data;
		    
		        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D * DATA_T_PER_RUN[run]; i++)
		        {
		            h_Data[i + accumulatedTRs * DATA_W * DATA_H * DATA_D] = (float)p[i];
		        }
		    }
		    else if ( inputData->datatype == DT_FLOAT )
		    {		
		        float *p = (float*)inputData->data;
		    
		        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D * DATA_T_PER_RUN[run]; i++)
		        {
		            h_Data[i + accumulatedTRs * DATA_W * DATA_H * DATA_D] = (float)p[i];
		        }
		    }
		    else
		    {
		        printf("Unknown data type in fMRI data for run %i, aborting!\n",run+1);
		        FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
				FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
		        return EXIT_FAILURE;
		    }
			accumulatedTRs += DATA_T_PER_RUN[run];
		}
	}
	else
	{
	    if ( inputData->datatype == DT_SIGNED_SHORT )
	    {
	        short int *p = (short int*)inputData->data;
    		
	        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D * DATA_T; i++)
    	    {
    	        h_Data[i] = (float)p[i];
    	    }
		}
	    else if ( inputData->datatype == DT_UINT8 )
	    {
	        unsigned char *p = (unsigned char*)inputData->data;
	    
	        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D * DATA_T; i++)
	        {
	            h_Data[i] = (float)p[i];
	        }
	    }
	    else if ( inputData->datatype == DT_UINT16 )
	    {
	        unsigned short int *p = (unsigned short int*)inputData->data;
	    
	        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D * DATA_T; i++)
	        {
	            h_Data[i] = (float)p[i];
	        }
	    }
		// Correct data type, just copy the pointer
	    else if ( inputData->datatype == DT_FLOAT )
	    {		
			h_Data = (float*)inputData->data;
	
			// Save the pointer in the pointer list
			allMemoryPointers[numberOfMemoryPointers] = (void*)h_Data;
	        numberOfMemoryPointers++;
	    }	
	}

	if (MULTIPLE_RUNS)
	{
		// Free memory, as it is stored in a big array
		for (size_t run = 0; run < NUMBER_OF_RUNS; run++)
		{
			inputData = allfMRINiftiImages[run];
			free(inputData->data);
			inputData->data = NULL;
		}
	}
	else
	{	
		// Free input fMRI data, it has been converted to floats
		if ( inputData->datatype != DT_FLOAT )
		{		
			free(inputData->data);
			inputData->data = NULL;
		}
		// Pointer has been copied to h_Data and pointer list, so set the input data pointer to NULL
		else
		{		
			inputData->data = NULL;
		}
	}


	// Mask is provided by user
	if (MASK)
	{
	    if ( inputMask->datatype == DT_SIGNED_SHORT )
	    {
	        short int *p = (short int*)inputMask->data;
    
	        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D; i++)
	        {
	            h_Mask[i] = (float)p[i];
	        }
	    }
	    else if ( inputMask->datatype == DT_FLOAT )
	    {
	        float *p = (float*)inputMask->data;
    
	        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D; i++)
        	{
	            h_Mask[i] = p[i];
	        }
	    }
	    else if ( inputMask->datatype == DT_UINT8 )
	    {
    	    unsigned char *p = (unsigned char*)inputMask->data;
    
	        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D; i++)
	        {
	            h_Mask[i] = (float)p[i];
	        }
	    }
	    else
	    {
	        printf("Unknown data type in mask volume, aborting!\n");
	        FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
			FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
	        return EXIT_FAILURE;
	    }
	}
	// Mask is NOT provided by user, set all mask voxels to 1
	else
	{
        for (size_t i = 0; i < DATA_W * DATA_H * DATA_D; i++)
        {
            h_Mask[i] = 1.0f;
        }
	}

	endTime = GetWallTime();

	if (VERBOS)
 	{
		printf("It took %f seconds to convert data to floats\n",(float)(endTime - startTime));
	}

	//------------------------
	// Setup the total designmatrix
    //------------------------

	Eigen::MatrixXd X = SetupGLMRegressorsFirstLevel(h_X_GLM, h_Motion_Parameters, h_Motion_Deriv_Parameters, NULL, RAW_REGRESSORS, RAW_DESIGNMATRIX, DATA_T_PER_RUN, REGRESS_MOTION, REGRESS_MOTION_DERIV, REGRESS_GLOBALMEAN, REGRESS_CONFOUNDS, NUMBER_OF_RUNS, USE_TEMPORAL_DERIVATIVES, NUMBER_OF_DETRENDING_REGRESSORS, NUMBER_OF_GLM_REGRESSORS, DATA_T, TR);

	// Use all regressors for variance
	Eigen::MatrixXd Z = X;

	// Motion deriv parameters should always be positive for variance,
	// since variance should always increase with a motion spike
	if (REGRESS_MOTION_DERIV)
	{
		// Last six regressors
		int startRegressor = NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION;
		for (int r = startRegressor; r < NUMBER_OF_TOTAL_GLM_REGRESSORS; r++)
		{
			for (int i = 0; i < DATA_T; i++)
			{
				Z(i,r) = fabs(Z(i,r));
			}
		}
	}

	// Demean regressors and standardize variance

	// Calculate which regressor contains only ones
	int MEAN_REGRESSOR;

	if (NUMBER_OF_RUNS == 1)
	{
		if (!RAW_DESIGNMATRIX)
		{
			MEAN_REGRESSOR = NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1);
		}
		else
		{
			MEAN_REGRESSOR = NUMBER_OF_GLM_REGRESSORS;
		}
			
		for (int r = 0; r < NUMBER_OF_TOTAL_GLM_REGRESSORS; r++)
		{
			if (r != MEAN_REGRESSOR)
			{
				Eigen::VectorXd regressor = X.block(0,r,DATA_T,1);
				DemeanRegressor(regressor,DATA_T);
				NormalizeVariance(regressor,DATA_T);
				X.block(0,r,DATA_T,1) = regressor;

				regressor = Z.block(0,r,DATA_T,1);
				DemeanRegressor(regressor,DATA_T);
				NormalizeVariance(regressor,DATA_T);
				Z.block(0,r,DATA_T,1) = regressor;
			}
		}
	}
	
	for (int i = 0; i < DATA_T; i++)
	{
		for (int r = 0; r < NUMBER_OF_TOTAL_GLM_REGRESSORS; r++)
		{
			h_Design_Matrix[i + r * DATA_T] = X(i,r);
		}
	}

	//-------------------------------------
	// Write total design matrix to file
	//-------------------------------------

	inputData = allfMRINiftiImages[0];

	if (WRITE_DESIGNMATRIX)
	{
		std::ofstream designmatrix;
	    
		const char* extension = "_total_designmatrix.txt";
		char* filenameWithExtension;

		CreateFilename(filenameWithExtension, inputData, extension, CHANGE_OUTPUT_FILENAME, outputFilename);

    	designmatrix.open(filenameWithExtension);    

	    if ( designmatrix.good() )
	    {
    	    for (size_t t = 0; t < DATA_T; t++)
	        {
	    	    for (size_t r = 0; r < NUMBER_OF_TOTAL_GLM_REGRESSORS; r++)
		        {
            		designmatrix << std::setprecision(6) << std::fixed << (double)h_Design_Matrix[t + r * DATA_T] << "  ";
				}
				designmatrix << std::endl;
			}
		    designmatrix.close();
        }  	
	    else
	    {
			designmatrix.close();
	        printf("Could not open the file for writing the total design matrix!\n");
	    }
		free(filenameWithExtension);
	}

	//------------------------------------------  
    // Read regressors for variance
	//------------------------------------------
	
	std::vector<int> gammaRegressorsVector; 
	Eigen::VectorXd gammaRegressorsEigen;

	if (GAMMA_REGRESSORS)
	{
		int regressor;
	   	std::ifstream gammaRegressors;    

    	gammaRegressors.open(GAMMA_REGRESSORS_NAME);  
	    if (!gammaRegressors.good())
	    {
    		gammaRegressors.close();
	        printf("Unable to open file for gamma regressors %s. Aborting! \n",GAMMA_REGRESSORS_NAME);
			FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
			FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    	    return EXIT_FAILURE;
		}
		
		q = 0;
		while (gammaRegressors >> regressor)
		{
			if (regressor >= p)
		    {
	    		gammaRegressors.close();
		        printf("Regressor %i for gamma regressors out of range. Number of regressors for beta is %i and counting starts from 0. Aborting! \n",regressor,p);
				FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
				FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
	    	    return EXIT_FAILURE;
			}
			q++;
			gammaRegressorsVector.push_back(regressor);
		}

		Eigen::MatrixXd Ztemp;
		Ztemp.resize(DATA_T,q);
		gammaRegressorsEigen.resize(q);

		for (int i = 0; i < q; i++)
		{
			Ztemp.col(i) = Z.col(gammaRegressorsVector.at(i));
			gammaRegressorsEigen(i) = gammaRegressorsVector.at(i);
		}
		Z = Ztemp;
    }
	else
	{
		gammaRegressorsEigen.resize(q);
		for (int i = 0; i < q; i++)
		{
			gammaRegressorsEigen(i) = i;
		}
	}

	if (VERBOS)
 	{
		std::cout << "gamma regressors are " << std::endl << std::endl << gammaRegressorsEigen.transpose() << std::endl << std::endl;
	}

	if (WRITE_DESIGNMATRIX)
	{
		std::ofstream designmatrix;
	    
		const char* extension = "_total_designmatrix_variance.txt";
		char* filenameWithExtension;

		CreateFilename(filenameWithExtension, inputData, extension, CHANGE_OUTPUT_FILENAME, outputFilename);

    	designmatrix.open(filenameWithExtension);    

	    if ( designmatrix.good() )
	    {
    	    for (size_t t = 0; t < DATA_T; t++)
	        {
	    	    for (size_t r = 0; r < q; r++)
		        {
            		designmatrix << std::setprecision(6) << std::fixed << (double)Z(t,r) << "  ";
				}
				designmatrix << std::endl;
			}
		    designmatrix.close();
        }  	
	    else
	    {
			designmatrix.close();
	        printf("Could not open the file for writing the total design matrix variance!\n");
	    }
		free(filenameWithExtension);
	}

    //------------------------
	// Setup the parameters for Hetero GLM
    //------------------------

	HeteroGauss HeteroGaussObj;

    Eigen::MatrixXd betaDraws, IbetaDraws;
    Eigen::MatrixXd gammaDraws, IgammaDraws;
    Eigen::MatrixXd rhoDraws, IrhoDraws;    
    Eigen::VectorXd accPrGammaDraws;

	Eigen::VectorXd timeseries;	timeseries.resize(DATA_T);
	
	Eigen::VectorXd muBeta; muBeta.resize(p);
    for (int i = 0; i < p; i++)
    {
        muBeta(i) = 0.0;
    }
	muBeta(NUMBER_OF_GLM_REGRESSORS) = 800; // intercept

	Eigen::VectorXd muGamma; muGamma.resize(q);
    for (int i = 0; i < q; i++)
    {
        muGamma(i) = 0.0;
    }

	if (VERBOS)
 	{
		std::cout << "mu beta is" << std::endl << std::endl << muBeta.transpose() << std::endl << std::endl;
		std::cout << "mu gamma is " << std::endl << std::endl << muGamma.transpose() << std::endl << std::endl;
	}


	Eigen::VectorXd PrInBeta; PrInBeta.resize(p);
    for (int i = 0; i < p; i++)
    {
        PrInBeta(i) = 0.5;
    }
    
    Eigen::VectorXd PrInGamma; PrInGamma.resize(q);
    for (int i = 0; i < q; i++)
    {
        PrInGamma(i) = 0.5;
    }
    
    Eigen::VectorXd PrInRho; PrInRho.resize(AR_ORDER);
    for (int i = 0; i < AR_ORDER; i++)
    {
        PrInRho(i) = 0.5 / sqrt(i+1);
    }

	if (VERBOS)
 	{
		std::cout << "PrInBeta is " << std::endl << std::endl << PrInBeta.transpose() << std::endl << std::endl;
		std::cout << "PrInGamma is " << std::endl << std::endl << PrInGamma.transpose() << std::endl << std::endl;
		std::cout << "PrInRho is " << std::endl << std::endl << PrInRho.transpose() << std::endl << std::endl;
	}

	


	//------------------------------------------  
    // Read on trial settings
	//------------------------------------------

	std::vector<int> betaTrialsVector; 
	std::vector<int> gammaTrialsVector; 
	std::vector<int> rhoTrialsVector; 

	if (ONTRIAL_BETA)
	{
		int regressor;
	   	std::ifstream betaTrials;    

    	betaTrials.open(ONTRIAL_BETA_NAME);  
	    if (!betaTrials.good())
	    {
    		betaTrials.close();
	        printf("Unable to open file for on trial beta %s. Aborting! \n",ONTRIAL_BETA_NAME);
			FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
			FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    	    return EXIT_FAILURE;
		}
		
		while (betaTrials >> regressor)
		{
			if (regressor >= p)
		    {
	    		betaTrials.close();
		        printf("Regressor %i for on trial beta is out of range. Number of regressors for beta is %i and counting starts from 0. Aborting! \n",regressor,p);
				FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
				FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
	    	    return EXIT_FAILURE;
			}
			nOnTrialBeta++;
			betaTrialsVector.push_back(regressor);
		}
    }

	if (ONTRIAL_GAMMA)
	{
		int regressor;
	   	std::ifstream gammaTrials;    

    	gammaTrials.open(ONTRIAL_GAMMA_NAME);  
	    if (!gammaTrials.good())
	    {
    		gammaTrials.close();
	        printf("Unable to open file for on trial gamma %s. Aborting! \n",ONTRIAL_GAMMA_NAME);
			FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
			FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    	    return EXIT_FAILURE;
		}
		
		while (gammaTrials >> regressor)
		{
			if (regressor >= q)
		    {
	    		gammaTrials.close();
		        printf("Regressor %i for on trial gamma is out of range. Number of regressors for gamma is %i and counting starts from 0. Aborting! \n",regressor,q);
				FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
				FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
	    	    return EXIT_FAILURE;
			}
			nOnTrialGamma++;
			gammaTrialsVector.push_back(regressor);
		}
    }

	if (ONTRIAL_RHO)
	{
		int regressor;
	   	std::ifstream rhoTrials;    

    	rhoTrials.open(ONTRIAL_RHO_NAME);  
	    if (!rhoTrials.good())
	    {
    		rhoTrials.close();
	        printf("Unable to open file for on trial rho %s. Aborting! \n",ONTRIAL_RHO_NAME);
			FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
			FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
    	    return EXIT_FAILURE;
		}
		
		while (rhoTrials >> regressor)
		{
			if (regressor >= AR_ORDER)
		    {
	    		rhoTrials.close();
		        printf("Coefficient %i for on trial rho is out of range. Number of AR coefficients is %i and counting starts from 0. Aborting! \n",regressor,AR_ORDER);
				FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
				FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
	    	    return EXIT_FAILURE;
			}
			nOnTrialRho++;
			rhoTrialsVector.push_back(regressor);
		}
    }

    Eigen::VectorXd onTrialBeta; onTrialBeta.resize(nOnTrialBeta);
    Eigen::VectorXd onTrialGamma; onTrialGamma.resize(nOnTrialGamma);   
    Eigen::VectorXd onTrialRho; onTrialRho.resize(nOnTrialRho);
    	
    for (int i = 0; i < nOnTrialBeta; i++)
    {
		onTrialBeta(i) = (double)betaTrialsVector.at(i);
	}
    for (int i = 0; i < nOnTrialGamma; i++)
    {
		onTrialGamma(i) = (double)gammaTrialsVector.at(i);
	}
    for (int i = 0; i < nOnTrialRho; i++)
    {
		onTrialRho(i) = (double)rhoTrialsVector.at(i);
	}

	if (VERBOS)
 	{
		std::cout << "onTrialBeta is " << std::endl << std::endl << onTrialBeta.transpose() << std::endl << std::endl;
		std::cout << "onTrialGamma is " << std::endl << std::endl << onTrialGamma.transpose() << std::endl << std::endl;
		std::cout << "onTrialRho is " << std::endl << std::endl << onTrialRho.transpose() << std::endl << std::endl;
	}

	Eigen::MatrixXd U, Xtilde;
	Eigen::VectorXd beta, u, ytilde, rho;
	Eigen::VectorXd tempVector;
		
    //------------------------
    
	startTime = GetWallTime();

	// Get number of brain voxels
	int numBrainVoxels = 0;
	for (int voxel = 0; voxel < (DATA_W * DATA_H * DATA_D); voxel++)
	{
		if (h_Mask[voxel] == 1.0)
		{
			numBrainVoxels++;
		}
	}

	// Run the HeteroGLM for each voxel in the mask

	int voxel = 0;
	int t = 0;
	int j = 0;
	int it = 0;
	int nonStationaryDraws = 0;
	int analyzedVoxels = 0;
	float analyzedPortion = 0.0f;
	float previousAnalyzedPortion = 0.0f;
	int pp = 0;

	#pragma omp parallel for shared (DATA_W,DATA_H,DATA_D,h_Mask,h_Data,h_Beta_Volumes,h_IBeta_Volumes,h_Gamma_Volumes,h_IGamma_Volumes,h_Rho_Volumes,h_IRho_Volumes,h_Beta_Posterior,h_Gamma_Posterior,h_Rho_Posterior,h_AccPr,analyzedVoxels,analyzedPortion,previousAnalyzedPortion, updateInclusion) private(voxel,t,j,it,pp) firstprivate(HeteroGaussObj,timeseries,betaDraws,IbetaDraws,gammaDraws,IgammaDraws,rhoDraws,IrhoDraws, accPrGammaDraws, beta, u, U, ytilde, Xtilde, rho, X, Z, AR_ORDER, forceStationarity, muBeta, muGamma, tauBeta, tauGamma, tauRho, iota, r, PrInBeta, PrInGamma, PrInRho, onTrialBeta, onTrialGamma, onTrialRho, MCMC_ITERATIONS, prcBurnin, nStepsGamma, hessMethodGamma, linkType, propDfGamma, IUpdatePrGamma, tempVector, nonStationaryDraws) 
	for (voxel = 0; voxel < (DATA_W * DATA_H * DATA_D); voxel++)
	{
		if (h_Mask[voxel] == 1.0)
		{
		    for (t = 0; t < DATA_T; t++)
		    {
		        timeseries(t) = (double)h_Data[voxel + t * DATA_W * DATA_H * DATA_D];
		    }

		    // Run calculations for current voxel
			HeteroGaussObj.GibbsHIGLM(betaDraws,IbetaDraws,gammaDraws,IgammaDraws, rhoDraws, IrhoDraws, accPrGammaDraws,  beta, u, U,  ytilde, Xtilde,  rho, timeseries, X, Z, AR_ORDER, forceStationarity, muBeta, muGamma, tauBeta, tauGamma, tauRho, iota, r, PrInBeta, PrInGamma, PrInRho, onTrialBeta, onTrialGamma, onTrialRho, MCMC_ITERATIONS, prcBurnin, nStepsGamma, hessMethodGamma, linkType, propDfGamma, IUpdatePrGamma, updateInclusion, nonStationaryDraws);
					
			h_AccPr[voxel] = accPrGammaDraws.sum() / accPrGammaDraws.size();

			if (forceStationarity)
			{
				h_NonStationary_Draws[voxel] = nonStationaryDraws;
			}

			if (SAVE_FULL_POSTERIOR)
			{
				for (j = 0; j < p; j++)
				{
					tempVector = betaDraws.col(j);

					for (it = 0; it < MCMC_ITERATIONS; it++)
					{
						// x, y, z, mcmc, regressor, 
						h_Beta_Posterior[voxel + it * DATA_W * DATA_H * DATA_D + j * DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS] = tempVector(it);
					}
				}

				for (j = 0; j < q; j++)
				{
					tempVector = gammaDraws.col(j);

					for (it = 0; it < MCMC_ITERATIONS; it++)
					{
						// x, y, z, mcmc, regressor, 
						h_Gamma_Posterior[voxel + it * DATA_W * DATA_H * DATA_D + j * DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS] = tempVector(it);
					}
				}

				for (j = 0; j < AR_ORDER; j++)
				{
					tempVector = rhoDraws.col(j);

					for (it = 0; it < MCMC_ITERATIONS; it++)
					{
						// x, y, z, mcmc, regressor, 
						h_Rho_Posterior[voxel + it * DATA_W * DATA_H * DATA_D + j * DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS] = tempVector(it);
					}
				}
			}

			// Save PPMs for original regressors
			for (j = 0; j < NUMBER_OF_GLM_REGRESSORS; j++)
			{
				tempVector = betaDraws.col(j);

				pp = 0;
				for (it = 0; it < MCMC_ITERATIONS; it++)
				{
					if (tempVector(it) > 0.0f)
					{
						pp++;
					}
				}

				h_PPM_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = (float)pp / (float)MCMC_ITERATIONS;
			}

			// Save posterior means
			for (j = 0; j < p; j++)
			{
				tempVector = betaDraws.col(j);
				h_Beta_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = tempVector.sum() / tempVector.size();

				tempVector = IbetaDraws.col(j);
				h_IBeta_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = tempVector.sum() / tempVector.size();
			}					

			// Save posterior means
			for (j = 0; j < q; j++)
			{
				tempVector = gammaDraws.col(j);
				h_Gamma_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = tempVector.sum() / tempVector.size();

				tempVector = IgammaDraws.col(j);
				h_IGamma_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = tempVector.sum() / tempVector.size();
			}

			// Save posterior means
			for (j = 0; j < AR_ORDER; j++)
			{
				tempVector = rhoDraws.col(j);
				h_Rho_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = tempVector.sum() / tempVector.size();

				tempVector = IrhoDraws.col(j);
				h_IRho_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = tempVector.sum() / tempVector.size();
			}

			#pragma omp critical
			{				
				double currentTime = GetWallTime();
				analyzedVoxels++;
				analyzedPortion = (float)analyzedVoxels/(float)numBrainVoxels * 100.0f;
				float timeLeft = (float)(currentTime - startTime) / analyzedPortion * (100.0f - analyzedPortion);

				if ((analyzedPortion - previousAnalyzedPortion) > 1.0f)
				{
					printf("Analyzed %f %% of %i voxels in mask in %f minutes, expected time remaining %f minutes \n",analyzedPortion,numBrainVoxels,(float)(currentTime - startTime)/60.0f,timeLeft/60.0f);	
					previousAnalyzedPortion = analyzedPortion;		
				}
			}
		}
		else
		{
			h_AccPr[voxel] = 0.0;

			if (SAVE_FULL_POSTERIOR)
			{
				for (j = 0; j < p; j++)
				{
					for (it = 0; it < MCMC_ITERATIONS; it++)
					{
						// x, y, z, mcmc, regressor, 
						h_Beta_Posterior[voxel + it * DATA_W * DATA_H * DATA_D + j * DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS] = 0.0;
					}
				}

				for (j = 0; j < q; j++)
				{
					for (it = 0; it < MCMC_ITERATIONS; it++)
					{
						// x, y, z, mcmc, regressor, 
						h_Gamma_Posterior[voxel + it * DATA_W * DATA_H * DATA_D + j * DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS] = 0.0;
					}
				}

				for (j = 0; j < AR_ORDER; j++)
				{
					for (it = 0; it < MCMC_ITERATIONS; it++)
					{
						// x, y, z, mcmc, regressor, 
						h_Rho_Posterior[voxel + it * DATA_W * DATA_H * DATA_D + j * DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS] = 0.0;
					}
				}													
			}
					
			// Save PPMs for original regressors
			for (j = 0; j < NUMBER_OF_GLM_REGRESSORS; j++)
			{
				h_PPM_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = 0.0;
			}

			for (j = 0; j < p; j++)
			{
				h_Beta_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = 0.0;
				h_IBeta_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = 0.0;
			}

			for (j = 0; j < q; j++)
			{
				h_Gamma_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = 0.0;
				h_IGamma_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = 0.0;
			}

			for (j = 0; j < AR_ORDER; j++)
			{
				h_Rho_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = 0.0;
				h_IRho_Volumes[voxel + j * DATA_W * DATA_H * DATA_D] = 0.0;
			}					
		}
	} 
	
	endTime = GetWallTime();

	if (VERBOS)
 	{
		printf("\nIt took %f seconds to run the Hetero GLM\n",(float)(endTime - startTime));
	}

	//-------------------------------------
	// Write results to nifti files
	//-------------------------------------

    // Create new nifti image
	nifti_image *outputNifti = nifti_copy_nim_info(inputData);      
	nifti_free_extensions(outputNifti);
    allNiftiImages[numberOfNiftiImages] = outputNifti;
	numberOfNiftiImages++;

	if (CHANGE_OUTPUT_FILENAME)
	{
		nifti_set_filenames(outputNifti, outputFilename, 0, 1);
	}

	outputNifti->nt = 1;
	outputNifti->ndim = 3;
	outputNifti->dim[0] = 3;
    outputNifti->dim[4] = 1;
    outputNifti->nvox = DATA_W * DATA_H * DATA_D;

    std::string ppmString = "_ppm";
    std::string betaString = "_beta";
	std::string IbetaString = "_Ibeta";
    std::string gammaString = "_gamma";
	std::string IgammaString = "_Igamma";
    std::string rhoString = "_rho";
	std::string IrhoString = "_Irho";
	std::string nonStationaryString = "_nonstationary_draws";


	std::string accpr = "_accpr";

	WriteNifti(outputNifti,h_AccPr,accpr.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);

	for (size_t i = 0; i < NUMBER_OF_GLM_REGRESSORS; i++)
	{
		std::string temp = ppmString;

	    std::stringstream ss;
		if ((i+1) < 10)
		{
        	ss << "_regressor000";
		}
		ss << i + 1;

		temp.append(ss.str());

		WriteNifti(outputNifti,&h_PPM_Volumes[i * DATA_W * DATA_H * DATA_D],temp.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
	}

	for (size_t i = 0; i < p; i++)
	{
		std::string temp1 = betaString;
		std::string temp2 = IbetaString;

	    std::stringstream ss;
		if ((i+1) < 10)
		{
        	ss << "_regressor000";
		}
		else if ((i+1) < 100)
		{
			ss << "_regressor00";
		}
		else if ((i+1) < 1000)
		{
			ss << "_regressor0";
		}
		else
		{
			ss << "_regressor";
		}						
		ss << i + 1;

		temp1.append(ss.str());
		temp2.append(ss.str());

		WriteNifti(outputNifti,&h_Beta_Volumes[i * DATA_W * DATA_H * DATA_D],temp1.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
		WriteNifti(outputNifti,&h_IBeta_Volumes[i * DATA_W * DATA_H * DATA_D],temp2.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
	}

	for (size_t i = 0; i < q; i++)
	{
		std::string temp1 = gammaString;
		std::string temp2 = IgammaString;

	    std::stringstream ss;
		if ((i+1) < 10)
		{
        	ss << "_regressor000";
		}
		else if ((i+1) < 100)
		{
			ss << "_regressor00";
		}
		else if ((i+1) < 1000)
		{
			ss << "_regressor0";
		}
		else
		{
			ss << "_regressor";
		}						
		ss << i + 1;

		temp1.append(ss.str());
		temp2.append(ss.str());

		WriteNifti(outputNifti,&h_Gamma_Volumes[i * DATA_W * DATA_H * DATA_D],temp1.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
		WriteNifti(outputNifti,&h_IGamma_Volumes[i * DATA_W * DATA_H * DATA_D],temp2.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
	}

	for (size_t i = 0; i < AR_ORDER; i++)
	{
		std::string temp1 = rhoString;
		std::string temp2 = IrhoString;

	    std::stringstream ss;
		if ((i+1) < 10)
		{
        	ss << "_arcoefficient000";
		}
		else if ((i+1) < 100)
		{
			ss << "_arcoefficient00";
		}
		ss << i + 1;

		temp1.append(ss.str());
		temp2.append(ss.str());

		WriteNifti(outputNifti,&h_Rho_Volumes[i * DATA_W * DATA_H * DATA_D],temp1.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
		WriteNifti(outputNifti,&h_IRho_Volumes[i * DATA_W * DATA_H * DATA_D],temp2.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
	}

	if (forceStationarity)
	{
		std::string temp1 = nonStationaryString;
		WriteNifti(outputNifti,h_NonStationary_Draws,temp1.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
	}

	if (SAVE_FULL_POSTERIOR)
	{
		outputNifti->nt = MCMC_ITERATIONS;
		outputNifti->ndim = 4;
		outputNifti->dim[0] = 4;
	    outputNifti->dim[4] = MCMC_ITERATIONS;
	    outputNifti->nvox = DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS;

		for (size_t i = 0; i < p; i++)
		{	
			std::string temp1 = betaString;
		
		    std::stringstream ss;
			if ((i+1) < 10)
			{
    	    	ss << "_fullposterior_regressor000";
			}
			else if ((i+1) < 100)
			{
				ss << "_fullposterior_regressor00";
			}
			ss << i + 1;

			temp1.append(ss.str());

			WriteNifti(outputNifti,&h_Beta_Posterior[i * DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS],temp1.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
		}

		for (size_t i = 0; i < q; i++)
		{	
			std::string temp1 = gammaString;
		
		    std::stringstream ss;
			if ((i+1) < 10)
			{
    	    	ss << "_fullposterior_regressor000";
			}
			else if ((i+1) < 100)
			{
				ss << "_fullposterior_regressor00";
			}
			ss << i + 1;

			temp1.append(ss.str());

			WriteNifti(outputNifti,&h_Gamma_Posterior[i * DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS],temp1.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
		}

		for (size_t i = 0; i < AR_ORDER; i++)
		{	
			std::string temp1 = rhoString;
		
		    std::stringstream ss;
			if ((i+1) < 10)
			{
    	    	ss << "_fullposterior_arcoefficient000";
			}
			else if ((i+1) < 100)
			{
    	    	ss << "_fullposterior_arcoefficient00";
			}
			ss << i + 1;

			temp1.append(ss.str());

			WriteNifti(outputNifti,&h_Rho_Posterior[i * DATA_W * DATA_H * DATA_D * MCMC_ITERATIONS],temp1.c_str(),ADD_FILENAME,DONT_CHECK_EXISTING_FILE);
		}
	}

    // Free all memory
    FreeAllMemory(allMemoryPointers,numberOfMemoryPointers);
    FreeAllNiftiImages(allNiftiImages,numberOfNiftiImages);
        
    return EXIT_SUCCESS;
}


