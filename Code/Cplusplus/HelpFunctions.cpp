#include <cfloat>
#include <limits>
#include <time.h>
#include <sys/time.h>

void CreateFilename(char *& filenameWithExtension, nifti_image* inputNifti, const char* extension, bool CHANGE_OUTPUT_FILENAME, const char* outputFilename)
{
	// Find the dot in the original filename
	if (!CHANGE_OUTPUT_FILENAME)
	{
	    const char* p = inputNifti->fname;
	    int dotPosition = 0;
	    while ( (p != NULL) && ((*p) != '.') )
	    {
	        p++;
	        dotPosition++;
	    }
    		
	    // Allocate temporary array
	    filenameWithExtension = (char*)malloc(strlen(inputNifti->fname) + strlen(extension) + 1);
		
	    // Copy filename to the dot
	    strncpy(filenameWithExtension,inputNifti->fname,dotPosition);
		filenameWithExtension[dotPosition] = '\0';
	}
	else
	{
	    const char* p = outputFilename;
	    int dotPosition = 0;
		int i = 0;
	    while ( (p != NULL) && ((*p) != '.') && (i < strlen(outputFilename)) )
	    {
	        p++;
	        dotPosition++;
			i++;
	    }
    	
	    // Allocate temporary array
	    filenameWithExtension = (char*)malloc(strlen(outputFilename) + strlen(extension) + 1);
		
	    // Copy filename to the dot
	    strncpy(filenameWithExtension,outputFilename,dotPosition);
		filenameWithExtension[dotPosition] = '\0';
	}
	
	// Add the extension
	strcat(filenameWithExtension,extension);
}

void LowpassFilterRegressor(float* h_LowpassFiltered_Regressor, float* h_Regressor, int DATA_T, int HIGHRES_FACTOR, float TR)
{
	// Allocate memory for lowpass filter
	int FILTER_LENGTH = 151;
	float* h_Filter = (float*)malloc(FILTER_LENGTH * sizeof(float));

	// Create lowpass filter
	int halfSize = (FILTER_LENGTH - 1) / 2;
	double smoothing_FWHM = 150.0;
	double sigma = smoothing_FWHM / 2.354 / (double)TR;
	double sigma2 = 2.0 * sigma * sigma;

	double u;
	float sum = 0.0f;
	for (int i = 0; i < FILTER_LENGTH; i++)
	{
		u = (double)(i - halfSize);
		h_Filter[i] = (float)exp(-pow(u,2.0) / sigma2);
		sum += h_Filter[i];
	}

	// Normalize
	for (int i = 0; i < FILTER_LENGTH; i++)
	{
		h_Filter[i] /= sum;
	}

	// Convolve regressor with filter
	for (int t = 0; t < (DATA_T * HIGHRES_FACTOR); t++)
	{
		h_LowpassFiltered_Regressor[t] = 0.0f;

		// 1D convolution
		//int offset = -(int)(((float)HRF_LENGTH - 1.0f)/2.0f);
		int offset = -(int)(((float)FILTER_LENGTH - 1.0f)/2.0f);
		for (int tt = FILTER_LENGTH - 1; tt >= 0; tt--)
		{
			if ( ((t + offset) >= 0) && ((t + offset) < (DATA_T * HIGHRES_FACTOR)) )
			{
				h_LowpassFiltered_Regressor[t] += h_Regressor[t + offset] * h_Filter[tt];
			}
			offset++;
		}
	}

	free(h_Filter);
}




void LowpassFilterRegressors(float* h_LowpassFiltered_Regressors, float* h_Regressors, int DATA_T, int HIGHRES_FACTOR, float TR, int NUMBER_OF_REGRESSORS)
{
	// Allocate memory for lowpass filter
	int FILTER_LENGTH = 151;
	float* h_Filter = (float*)malloc(FILTER_LENGTH * sizeof(float));

	// Create lowpass filter
	int halfSize = (FILTER_LENGTH - 1) / 2;
	double smoothing_FWHM = 150.0;
	double sigma = smoothing_FWHM / 2.354 / (double)TR;
	double sigma2 = 2.0 * sigma * sigma;

	double u;
	float sum = 0.0f;
	for (int i = 0; i < FILTER_LENGTH; i++)
	{
		u = (double)(i - halfSize);
		h_Filter[i] = (float)exp(-pow(u,2.0) / sigma2);
		sum += h_Filter[i];
	}

	// Normalize
	for (int i = 0; i < FILTER_LENGTH; i++)
	{
		h_Filter[i] /= sum;
	}

	for (int r = 0; r < NUMBER_OF_REGRESSORS; r++)
	{
		// Convolve regressor with filter
		for (int t = 0; t < (DATA_T * HIGHRES_FACTOR); t++)
		{
			h_LowpassFiltered_Regressors[t + r * DATA_T*HIGHRES_FACTOR] = 0.0f;
	
			// 1D convolution
			//int offset = -(int)(((float)HRF_LENGTH - 1.0f)/2.0f);
			int offset = -(int)(((float)FILTER_LENGTH - 1.0f)/2.0f);
			for (int tt = FILTER_LENGTH - 1; tt >= 0; tt--)
			{
				if ( ((t + offset) >= 0) && ((t + offset) < (DATA_T * HIGHRES_FACTOR)) )
				{
					h_LowpassFiltered_Regressors[t + r * DATA_T*HIGHRES_FACTOR] += h_Regressors[t + offset + r * DATA_T*HIGHRES_FACTOR] * h_Filter[tt];
				}
				offset++;
			}
		}
	}

	free(h_Filter);
}


void FreeAllMemory(void **pointers, int N)
{
    for (int i = 0; i < N; i++)
    {
        if (pointers[i] != NULL)
        {
            free(pointers[i]);
        }
    }
}

void FreeAllNiftiImages(nifti_image **niftiImages, int N)
{
    for (int i = 0; i < N; i++)
    {
		if (niftiImages[i] != NULL)
		{
			nifti_image_free(niftiImages[i]);
		}
    }
}

void ReadBinaryFile(float* pointer, int size, const char* filename, void** pointers, int& Npointers, nifti_image** niftiImages, int Nimages)
{
	if (pointer == NULL)
    {
        printf("The provided pointer for file %s is NULL, aborting! \n",filename);
        FreeAllMemory(pointers,Npointers);
		FreeAllNiftiImages(niftiImages,Nimages);
        exit(EXIT_FAILURE);
	}	

	FILE *fp = NULL; 
	fp = fopen(filename,"rb");

    if (fp != NULL)
    {
        fread(pointer,sizeof(float),size,fp);
        fclose(fp);
    }
    else
    {
        printf("Could not open %s , aborting! \n",filename);
        FreeAllMemory(pointers,Npointers);
		FreeAllNiftiImages(niftiImages,Nimages);
        exit(EXIT_FAILURE);
    }
}

void AllocateMemory(float *& pointer, size_t size, void** pointers, int& Npointers, nifti_image** niftiImages, int Nimages, size_t& allocatedMemory, const char* variable)
{
    pointer = (float*)malloc(size);
    if (pointer != NULL)
    {
        pointers[Npointers] = (void*)pointer;
        Npointers++;
		allocatedMemory += size;
    }
    else
    {
   		perror ("The following error occurred");
	    printf("Could not allocate host memory for variable %s ! \n",variable);     
	 	FreeAllMemory(pointers, Npointers);
		FreeAllNiftiImages(niftiImages, Nimages);
		exit(EXIT_FAILURE);        
    }
}

void AllocateMemoryInt(unsigned short int *& pointer, size_t size, void** pointers, int& Npointers, nifti_image** niftiImages, int Nimages, size_t allocatedMemory, const char* variable)
{
    pointer = (unsigned short int*)malloc(size);
    if (pointer != NULL)
    {
        pointers[Npointers] = (void*)pointer;
        Npointers++;
		allocatedMemory += size;
    }
    else
    {
		perror ("The following error occurred");
        printf("Could not allocate host memory for variable %s ! \n",variable);        
		FreeAllMemory(pointers, Npointers);
		FreeAllNiftiImages(niftiImages, Nimages);
		exit(EXIT_FAILURE);        
    }
}


float mymax(float* data, int N)
{
	float max = -100000.0f;
	for (int i = 0; i < N; i++)
	{
		if (data[i] > max)
			max = data[i];
	}

	return max;
}


float mymin(float* data, int N)
{
	float min = 100000.0f;
	for (int i = 0; i < N; i++)
	{
		if (data[i] < min)
			min = data[i];
	}

	return min;
}

bool WriteNifti(nifti_image* inputNifti, float* data, const char* filename, bool addFilename, bool checkFilename)
{       
	if (data == NULL)
    {
        printf("The provided data pointer for file %s is NULL, aborting writing nifti file! \n",filename);
		return false;
	}	
	if (inputNifti == NULL)
    {
        printf("The provided nifti pointer for file %s is NULL, aborting writing nifti file! \n",filename);
		return false;
	}	


    char* filenameWithExtension;
    
    // Add the provided filename extension to the original filename, before the dot
    if (addFilename)
    {
        // Find the dot in the original filename
        const char* p = inputNifti->fname;
        int dotPosition = 0;
        while ( (p != NULL) && ((*p) != '.') )
        {
            p++;
            dotPosition++;
        }
    
        // Allocate temporary array
        filenameWithExtension = (char*)malloc(strlen(inputNifti->fname) + strlen(filename) + 1);
        if (filenameWithExtension == NULL)
        {
            printf("Could not allocate temporary host memory! \n");      
            return false;
        }
    
        // Copy filename to the dot
        strncpy(filenameWithExtension,inputNifti->fname,dotPosition);
        filenameWithExtension[dotPosition] = '\0';
        // Add the extension
        strcat(filenameWithExtension,filename);
        // Add the rest of the original filename
        strcat(filenameWithExtension,inputNifti->fname+dotPosition);    
    }
        
    // Copy information from input data
    nifti_image *outputNifti = nifti_copy_nim_info(inputNifti);    
    // Set data pointer 
    outputNifti->data = (void*)data;        
    // Set data type to float
    outputNifti->datatype = DT_FLOAT;
    outputNifti->nbyper = 4;    
    
	// Change cal_min and cal_max, to get the scaling right in AFNI and FSL
    int N = inputNifti->nx * inputNifti->ny * inputNifti->nz * inputNifti->nt;
	outputNifti->cal_min = mymin(data,N);
	outputNifti->cal_max = mymax(data,N);

    // Change filename and write
    bool written = false;
    if (addFilename)
    {
        if ( nifti_set_filenames(outputNifti, filenameWithExtension, checkFilename, 1) == 0)
        {
            nifti_image_write(outputNifti);
            written = true;
        }
    }
    else if (!addFilename)
    {
        if ( nifti_set_filenames(outputNifti, filename, checkFilename, 1) == 0)
        {
            nifti_image_write(outputNifti);
            written = true;
        }                
    }    
    
    outputNifti->data = NULL;
    nifti_image_free(outputNifti);

    if (addFilename)
    {
        free(filenameWithExtension);
    } 
        
    if (written)
    {      
        return true;
    }
    else
    {
        return false;
    }                        
}

double GetWallTime()
{
    struct timeval time;
    if (gettimeofday(&time,NULL))
    {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

// This code was copied (with permission) from http://www.johndcook.com/Gamma.cpp
// (lgamma seems to exist for Linux compilers, but not for Visual studio)

double LogGamma(double x);

double Gamma(double x);


double LogGamma(double x)
{
    if (x < 12.0)
    {
        return log(fabs(Gamma(x)));
    }

	// Abramowitz and Stegun 6.1.41
    // Asymptotic series should be good to at least 11 or 12 figures
    // For error analysis, see Whittiker and Watson
    // A Course in Modern Analysis (1927), page 252

    static const double c[8] =
    {
		 1.0/12.0,
		-1.0/360.0,
		1.0/1260.0,
		-1.0/1680.0,
		1.0/1188.0,
		-691.0/360360.0,
		1.0/156.0,
		-3617.0/122400.0
    };
    double z = 1.0/(x*x);
    double sum = c[7];
    for (int i=6; i >= 0; i--)
    {
        sum *= z;
        sum += c[i];
    }
    double series = sum/x;

    static const double halfLogTwoPi = 0.91893853320467274178032973640562;
    double logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;
	return logGamma;
}

double Gamma(double x)
{
    // Split the function domain into three intervals:
    // (0, 0.001), [0.001, 12), and (12, infinity)

    ///////////////////////////////////////////////////////////////////////////
    // First interval: (0, 0.001)
	//
	// For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
	// So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
	// The relative error over this interval is less than 6e-7.

	const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

    if (x < 0.001)
        return 1.0/(x*(1.0 + gamma*x));

    ///////////////////////////////////////////////////////////////////////////
    // Second interval: [0.001, 12)

	if (x < 12.0)
    {
        // The algorithm directly approximates gamma over (1,2) and uses
        // reduction identities to reduce other arguments to this interval.

		double y = x;
        int n = 0;
        bool arg_was_less_than_one = (y < 1.0);

        // Add or subtract integers as necessary to bring y into (1,2)
        // Will correct for this below
        if (arg_was_less_than_one)
        {
            y += 1.0;
        }
        else
        {
            n = static_cast<int> (floor(y)) - 1;  // will use n later
            y -= n;
        }

        // numerator coefficients for approximation over the interval (1,2)
        static const double p[] =
        {
            -1.71618513886549492533811E+0,
             2.47656508055759199108314E+1,
            -3.79804256470945635097577E+2,
             6.29331155312818442661052E+2,
             8.66966202790413211295064E+2,
            -3.14512729688483675254357E+4,
            -3.61444134186911729807069E+4,
             6.64561438202405440627855E+4
        };

        // denominator coefficients for approximation over the interval (1,2)
        static const double q[] =
        {
            -3.08402300119738975254353E+1,
             3.15350626979604161529144E+2,
            -1.01515636749021914166146E+3,
            -3.10777167157231109440444E+3,
             2.25381184209801510330112E+4,
             4.75584627752788110767815E+3,
            -1.34659959864969306392456E+5,
            -1.15132259675553483497211E+5
        };

        double num = 0.0;
        double den = 1.0;
        int i;

        double z = y - 1;
        for (i = 0; i < 8; i++)
        {
            num = (num + p[i])*z;
            den = den*z + q[i];
        }
        double result = num/den + 1.0;

        // Apply correction if argument was not initially in (1,2)
        if (arg_was_less_than_one)
        {
            // Use identity gamma(z) = gamma(z+1)/z
            // The variable "result" now holds gamma of the original y + 1
            // Thus we use y-1 to get back the orginal y.
            result /= (y-1.0);
        }
        else
        {
            // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
            for (i = 0; i < n; i++)
                result *= y++;
        }

		return result;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Third interval: [12, infinity)

    if (x > 171.624)
    {
		// Correct answer too large to display. Force +infinity.
		double temp = DBL_MAX;
		return temp*2.0;
    }

    return exp(LogGamma(x));
}



// Gamma probability density function
float Gpdf(double value, double shape, double scale)
{
	return (exp( (shape - 1.0) * log(value) + shape * log(scale) - scale * value - LogGamma(shape) ));
}

// Creates an HRF as in SPM
void CreateHRF(float *& hrf, int& HRF_LENGTH, double TR)
{
	/*
	% p    - parameters of the response function (two gamma functions)
	%
	%							defaults
	%							(seconds)
	%	p(1) - delay of response (relative to onset)	   6
	%	p(2) - delay of undershoot (relative to onset)    16
	%	p(3) - dispersion of response			   1
	%	p(4) - dispersion of undershoot			   1
	%	p(5) - ratio of response to undershoot		   6
	%	p(6) - onset (seconds)				   0
	%	p(7) - length of kernel (seconds)		  32
	*/

	double p[7];
	p[0] = 6.0;
	p[1] = 16.0;
	p[2] = 1.0;
	p[3] = 1.0;
	p[4] = 6.0;
	p[5] = 0.0;
	p[6] = 32.0;
	double fMRI_T = 16.0;
	double dt = ((double)TR)/fMRI_T;

	int length = (int)(p[6]/dt);
	double* highres_hrf = (double*)malloc(sizeof(double) * length);

	for (int i = 0; i < length; i++)
	{
		highres_hrf[i] = (double)i - p[5]/dt;
	}

	for (int i = 0; i < length; i++)
	{
		highres_hrf[i] = Gpdf(highres_hrf[i],p[0]/p[2],dt/p[2]) - 1.0/p[4] * Gpdf(highres_hrf[i],p[1]/p[3],dt/p[3]);
	}

	// Downsample the hrf
	int downsample_factor = 16;
	HRF_LENGTH = length/downsample_factor + 1;

	hrf = (float*)malloc(sizeof(float) * HRF_LENGTH);

	for (int i = 0; i < HRF_LENGTH; i++)
	{
		if ((i * downsample_factor) < length)
		{
			hrf[i] = (float)highres_hrf[i*downsample_factor];
		}
		else
		{
			hrf[i] = 0.0f;
		}
	}

	// Normalize
	float sum = 0.0f;
	for (int i = 0; i < HRF_LENGTH; i++)
	{
		sum += hrf[i];
	}
	for (int i = 0; i < HRF_LENGTH; i++)
	{
		hrf[i] /= sum;
	}

	free(highres_hrf);
}



void GenerateRegressorTemporalDerivatives(float * Regressors_With_Temporal_Derivatives, float* Regressors, int NUMBER_OF_TIMEPOINTS, int NUMBER_OF_REGRESSORS)
{
	int rr = 0;
	for (int r = 0; r < NUMBER_OF_REGRESSORS; r++)
	{
		for (int t = 0; t < NUMBER_OF_TIMEPOINTS; t++)
		{
			// Copy original regressor
			Regressors_With_Temporal_Derivatives[t + rr * NUMBER_OF_TIMEPOINTS] = Regressors[t + r * NUMBER_OF_TIMEPOINTS];

			// Calculate derivative and save as second regressor
			if ( ((t-1) >= 0) )
			{
				Regressors_With_Temporal_Derivatives[t + (rr + 1) * NUMBER_OF_TIMEPOINTS] = Regressors[t + r * NUMBER_OF_TIMEPOINTS] - Regressors[(t-1) + r * NUMBER_OF_TIMEPOINTS];
			}
			else
			{
				Regressors_With_Temporal_Derivatives[t + (rr + 1) * NUMBER_OF_TIMEPOINTS] = 0.0f;
			}
		}
		rr += 2;
	}
}

void ConvolveRegressorsWithHRF(float* Convolved_Regressors, float* Regressors, double TR, int NUMBER_OF_TIMEPOINTS, int NUMBER_OF_REGRESSORS)
{
	float* hrf; 
	int HRF_LENGTH;

	CreateHRF(hrf, HRF_LENGTH, TR);

	for (int r = 0; r < NUMBER_OF_REGRESSORS; r++)
	{
		for (int t = 0; t < NUMBER_OF_TIMEPOINTS; t++)
		{
			Convolved_Regressors[t + r * NUMBER_OF_TIMEPOINTS] = 0.0f;

			// 1D convolution
			//int offset = -(int)(((float)HRF_LENGTH - 1.0f)/2.0f);
			int offset = -(int)(((float)HRF_LENGTH - 1.0f)/1.0f);;
			for (int tt = HRF_LENGTH - 1; tt >= 0; tt--)
			{
				if ( ((t + offset) >= 0) && ((t + offset) < NUMBER_OF_TIMEPOINTS) )
				{
					Convolved_Regressors[t + r * NUMBER_OF_TIMEPOINTS] += Regressors[t + offset + r * NUMBER_OF_TIMEPOINTS] * hrf[tt];
				}
				offset++;
			}
		}
	}

	free(hrf);
}


// Demeans a regressor stored as an Eigen variable
void DemeanRegressor(Eigen::VectorXd& Regressor, int N)
{
	double mean = 0.0;
	for (int t = 0; t < N; t++)
	{
		mean += Regressor(t);
	}
	mean /= (double)N;

	for (int t = 0; t < N; t++)
	{
		Regressor(t) -= mean;
	}
}

void DemeanRegressor(Eigen::VectorXf& Regressor, int N)
{
	float mean = 0.0;
	for (int t = 0; t < N; t++)
	{
		mean += Regressor(t);
	}
	mean /= (float)N;

	for (int t = 0; t < N; t++)
	{
		Regressor(t) -= mean;
	}
}



void NormalizeVariance(Eigen::VectorXd& Regressor, int N)
{
	double mean = Regressor.mean();

	double std = 0.0;
	for (int t = 0; t < N; t++)
	{
		std += (Regressor(t) - mean) * (Regressor(t) - mean);
	}
	std = sqrt( std / (double)(N - 1) );

	for (int t = 0; t < N; t++)
	{
		Regressor(t) = Regressor(t) / std;
	}
}




Eigen::MatrixXd SetupGLMRegressorsFirstLevel(float *h_X_GLM_In, 
											 float *h_Motion_Parameters,
											 float *h_Motion_Deriv_Parameters,
											 float *h_Global_Mean,
                                             bool RAW_REGRESSORS, 
                                             bool RAW_DESIGNMATRIX, 
                                             size_t *DATA_T_PER_RUN, 
                                             bool REGRESS_MOTION, 
                                             bool REGRESS_MOTION_DERIV, 
                                             bool REGRESS_GLOBALMEAN, 
                                             bool REGRESS_CONFOUNDS, 
                                             int NUMBER_OF_RUNS, 
                                             bool USE_TEMPORAL_DERIVATIVES,
											 int NUMBER_OF_DETRENDING_REGRESSORS,
	  										 int NUMBER_OF_GLM_REGRESSORS, 
 											 int DATA_T,
											 double TR)
{
	int NUMBER_OF_TOTAL_GLM_REGRESSORS;
	int NUMBER_OF_MOTION_REGRESSORS = 6;
	int NUMBER_OF_CONFOUND_REGRESSORS = 0;

	// Calculate total number of regressors
	if (!RAW_DESIGNMATRIX)
	{
		NUMBER_OF_TOTAL_GLM_REGRESSORS = NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION_DERIV + REGRESS_GLOBALMEAN + NUMBER_OF_CONFOUND_REGRESSORS*REGRESS_CONFOUNDS;
	}
	else
	{
		NUMBER_OF_TOTAL_GLM_REGRESSORS = NUMBER_OF_GLM_REGRESSORS + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION_DERIV + REGRESS_GLOBALMEAN;
	}

	float *h_X_GLM_Convolved = (float*)malloc(NUMBER_OF_GLM_REGRESSORS * (USE_TEMPORAL_DERIVATIVES+1) * DATA_T * sizeof(float));
	float *h_X_GLM_With_Temporal_Derivatives = (float*)malloc(NUMBER_OF_GLM_REGRESSORS * 2 * DATA_T * sizeof(float));
	float *h_X_GLM_Confounds = NULL;

	std::vector<Eigen::VectorXd> allOneRegressors;
	std::vector<Eigen::VectorXd> allLinearRegressors;
	std::vector<Eigen::VectorXd> allQuadraticRegressors;
	std::vector<Eigen::VectorXd> allCubicRegressors;

	bool meanRegressor[NUMBER_OF_TOTAL_GLM_REGRESSORS];
	bool detrendingRegressor[NUMBER_OF_TOTAL_GLM_REGRESSORS];
	int  detrendingRegressorRun[NUMBER_OF_TOTAL_GLM_REGRESSORS];
	int totalTRs[NUMBER_OF_RUNS];

	for (int r = 0; r < NUMBER_OF_TOTAL_GLM_REGRESSORS; r++)
	{
		meanRegressor[r] = false;
		detrendingRegressor[r] = false;
		detrendingRegressorRun[r] = 0;
	}

	totalTRs[0] = 0;
	for (size_t run = 1; run < NUMBER_OF_RUNS; run++)
	{
		totalTRs[run] = totalTRs[run-1] + DATA_T_PER_RUN[run-1];
	}

	// Create detrending regressors
	for (size_t run = 0; run < NUMBER_OF_RUNS; run++)
	{
		int N = DATA_T_PER_RUN[run];

		Eigen::VectorXd Ones(N,1);
		Eigen::VectorXd Linear(N,1);
		Eigen::VectorXd Quadratic(N,1);
		Eigen::VectorXd Cubic(N,1);

		// Ones and linear trend
		float offset = -((float)N - 1.0f)/2.0f;
		for (int t = 0; t < N; t++)
		{
			Ones(t) = 1.0;
			Linear(t) = offset + (double)t;
		}

		// Calculate quadratic and cubic trends
		Quadratic = Linear.cwiseProduct(Linear);
		Cubic = Linear.cwiseProduct(Linear);
		Cubic = Cubic.cwiseProduct(Linear);

		// Normalize
		Linear = Linear / Linear.maxCoeff();
		double minn = std::abs(Quadratic.minCoeff());
		double maxx = Quadratic.maxCoeff();
		if (maxx > minn)
		{
			Quadratic = Quadratic / maxx;
		}
		else
		{
			Quadratic = Quadratic / minn;
		}
		Cubic = Cubic / Cubic.maxCoeff();

		allOneRegressors.push_back(Ones);
		allLinearRegressors.push_back(Linear);
		allQuadraticRegressors.push_back(Quadratic);
		allCubicRegressors.push_back(Cubic);
	}


	// Create temporal derivatives if requested and then convolve all regressors with HRF
	if (USE_TEMPORAL_DERIVATIVES && !RAW_REGRESSORS && !RAW_DESIGNMATRIX)
	{
		GenerateRegressorTemporalDerivatives(h_X_GLM_With_Temporal_Derivatives, h_X_GLM_In, DATA_T, NUMBER_OF_GLM_REGRESSORS);
		ConvolveRegressorsWithHRF(h_X_GLM_Convolved, h_X_GLM_With_Temporal_Derivatives, TR, DATA_T, NUMBER_OF_GLM_REGRESSORS*2);
	}
	// Convolve regressors with HRF
	else if (!RAW_REGRESSORS && !RAW_DESIGNMATRIX)
	{
		ConvolveRegressorsWithHRF(h_X_GLM_Convolved, h_X_GLM_In, TR, DATA_T, NUMBER_OF_GLM_REGRESSORS);
	}
	// Just copy raw regressors
	else if (RAW_REGRESSORS || RAW_DESIGNMATRIX)
	{
		// Loop over samples
		for (int i = 0; i < DATA_T; i++)
		{
			// Loop over regressors
			for (int r = 0; r < NUMBER_OF_GLM_REGRESSORS; r++)
			{
				h_X_GLM_Convolved[i + r * DATA_T] = h_X_GLM_In[i + r * DATA_T];
			}
		}
	}

	// Setup total design matrix
	Eigen::MatrixXd X(DATA_T,NUMBER_OF_TOTAL_GLM_REGRESSORS);
	for (int i = 0; i < DATA_T; i++)
	{
		for (int r = 0; r < NUMBER_OF_TOTAL_GLM_REGRESSORS; r++)
		{
			X(i,r) = 0.0;		
		}	
	}

	// Detrending regressors
	size_t accumulatedTRs = 0;
	for (int run = 0; run < NUMBER_OF_RUNS; run++)
	{
		Eigen::VectorXd Ones = allOneRegressors[run];
		Eigen::VectorXd Linear = allLinearRegressors[run];
		Eigen::VectorXd Quadratic = allQuadraticRegressors[run];
		Eigen::VectorXd Cubic = allCubicRegressors[run];

		for (int i = 0; i < DATA_T_PER_RUN[run]; i++)
		{	
			X(i+accumulatedTRs,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 0) = Ones(i);

			if (NUMBER_OF_DETRENDING_REGRESSORS >= 2)
			{
				X(i+accumulatedTRs,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 1) = Linear(i);
			}
			if (NUMBER_OF_DETRENDING_REGRESSORS >= 3)
			{
				X(i+accumulatedTRs,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 2) = Quadratic(i);
			}
			if (NUMBER_OF_DETRENDING_REGRESSORS == 4)
			{
				X(i+accumulatedTRs,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 3) = Cubic(i);
			}
		}
		accumulatedTRs += DATA_T_PER_RUN[run];

		meanRegressor[NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 0] = true;
		if (NUMBER_OF_DETRENDING_REGRESSORS >= 2)
		{
			detrendingRegressor[NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 1] = true;
			detrendingRegressorRun[NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 1] = run;
		}
		if (NUMBER_OF_DETRENDING_REGRESSORS >= 3)
		{
			detrendingRegressor[NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 2] = true;
			detrendingRegressorRun[NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 2] = run;
		}
		if (NUMBER_OF_DETRENDING_REGRESSORS == 4)
		{
			detrendingRegressor[NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 3] = true;
			detrendingRegressorRun[NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + run*NUMBER_OF_DETRENDING_REGRESSORS + 3] = run;
		}						
	}

	// Loop over samples
	for (int i = 0; i < DATA_T; i++)
	{
		// Regressors for paradigms (number of regressors is 0 for regressonly)
		for (int r = 0; r < NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1); r++)
		{
			X(i,r) = (double)h_X_GLM_Convolved[i + r * DATA_T];
		}

		if (REGRESS_MOTION)
		{
			// Motion regressors
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + 0) = h_Motion_Parameters[i + 0 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + 1) = h_Motion_Parameters[i + 1 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + 2) = h_Motion_Parameters[i + 2 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + 3) = h_Motion_Parameters[i + 3 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + 4) = h_Motion_Parameters[i + 4 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + 5) = h_Motion_Parameters[i + 5 * DATA_T];
		}


		if (REGRESS_MOTION_DERIV)
		{
			// Motion regressors
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + 0) = h_Motion_Deriv_Parameters[i + 0 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + 1) = h_Motion_Deriv_Parameters[i + 1 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + 2) = h_Motion_Deriv_Parameters[i + 2 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + 3) = h_Motion_Deriv_Parameters[i + 3 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + 4) = h_Motion_Deriv_Parameters[i + 4 * DATA_T];
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + 5) = h_Motion_Deriv_Parameters[i + 5 * DATA_T];
		}


		if (REGRESS_GLOBALMEAN)
		{
			X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION_DERIV) = (double)h_Global_Mean[i];
		}

		if (REGRESS_CONFOUNDS && !RAW_DESIGNMATRIX)
		{
			// Confounding regressors
			for (int r = 0; r < NUMBER_OF_CONFOUND_REGRESSORS; r++)
			{
				X(i,NUMBER_OF_GLM_REGRESSORS*(USE_TEMPORAL_DERIVATIVES+1) + NUMBER_OF_DETRENDING_REGRESSORS*NUMBER_OF_RUNS + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION + NUMBER_OF_MOTION_REGRESSORS*REGRESS_MOTION_DERIV + REGRESS_GLOBALMEAN + r) = (double)h_X_GLM_Confounds[i + r * DATA_T];
			}
		}
	}

	free(h_X_GLM_With_Temporal_Derivatives);
	free(h_X_GLM_Convolved);

	return X;
}


