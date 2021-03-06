
#include "heterogauss.h"
#include <Dense>
#include <Eigen327/Cholesky>
#include <string>
#include <string.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <chrono>
#include <iostream>

#define EIGEN_DONT_PARALLELIZE

HeteroGauss::HeteroGauss()
{
	// Make sure we get different numbers each time we run the code
	unsigned seed_ = std::chrono::system_clock::now().time_since_epoch().count();
	generator.seed(seed_);
}

HeteroGauss::~HeteroGauss()
{

}




double HeteroGauss::RndBeta(double a, double b)
{
	std::gamma_distribution<double> distribution1(a,1.0);
	std::gamma_distribution<double> distribution2(b,1.0);

	double x = distribution1(generator);
	double y = distribution2(generator);

	return (x / (x + y));	
}

double HeteroGauss::RndChi2(int df)
{
	double sum = 0.0;

	for (int i = 0; i < df; i++)
	{
		double number = normalDistribution(generator);
		sum += number*number;
	}

	return sum;
}


void HeteroGauss::RndMultiT(Eigen::VectorXd &draw, bool &error, Eigen::VectorXd &mu, Eigen::MatrixXd &Sigma, int df)
{
	//
	// PURPOSE:     Random number generator for the k-variate student-t distribution with mean mu, 
	//              covariance matrix, Sigma and df degrees of freedom
	//
	// INPUT:       mu       (k-by-1)        
	//              Sigma    (k-by-k)            
	//              df       (Scalar)            
	//
	// OUTPUT:      T        (k-by-1)    A draw from the multivariate student-t.          
	//
	// AUTHOR:      Mattias Villani, Research Department, Sveriges Riksbank and
	//              Department of Statistics, Stockholm University. 
	//
	// C++ AUTHOR:  Anders Eklund, Linköping university
	//              
	//
	// REVISED:     2016-01-19
	//

	int k = mu.size();

	// Cholesky
	Eigen::LLT<Eigen::MatrixXd> lltOfSigma(Sigma);
	Eigen::MatrixXd C = lltOfSigma.matrixU();

	if (lltOfSigma.info() != Eigen::Success)
	{
		error = true;
		draw.resize(0);
	}	
	else
	{
		error = false;

		Eigen::VectorXd z;
		z.resize(k);
	
		for (int i = 0; i < k; i++)
		{
			z(i) = normalDistribution(generator);
		}
	
		double x = RndChi2(df);
		
		draw = mu + C.transpose() * z / sqrt( x / ((double)df - 2.0) ); // Bauwens et al: T=mu+C*z/sqrt(x/(df-2));	
	}
}


void HeteroGauss::RndMultiN(Eigen::VectorXd &draw, bool &error, Eigen::VectorXd &mu, Eigen::MatrixXd &Sigma)
{
	//
	// PURPOSE:     Random number generator for the k-variate normal distribution with mean mu, 
	//              covariance matrix Sigma
	//
	// INPUT:       mu       (k-by-1)        
	//              Sigma    (k-by-k)            
	//
	// OUTPUT:      N        (k-by-1)    A draw from the multivariate normal 
	//
	// AUTHOR:      Mattias Villani, Research Department, Sveriges Riksbank and
	//              Department of Statistics, Stockholm University. 
	//
	// C++ AUTHOR:  Anders Eklund, Linköping university
	//              
	//
	// REVISED:     2016-02-05
	//

	int k = mu.size();

	// Cholesky
	Eigen::LLT<Eigen::MatrixXd> lltOfSigma(Sigma);
	Eigen::MatrixXd C = lltOfSigma.matrixU();

	if (lltOfSigma.info() != Eigen::Success)
	{
		error = true;
		draw.resize(0);
	}	
	else
	{
		error = false;

		Eigen::VectorXd z;
		z.resize(k);
	
		for (int i = 0; i < k; i++)
		{
			z(i) = normalDistribution(generator);
		}
		
		draw = mu + C.transpose() * z;
	}
}


double HeteroGauss::LogPdfMultiN(Eigen::VectorXd &x, Eigen::VectorXd &mu, Eigen::MatrixXd &Sigma)
{
	if ((mu.size() > 0) && (x.size() > 0))
	{
		double LogDetSigma = log(Sigma.determinant());
		Eigen::MatrixXd InvSigma = Sigma.inverse();
		double p = (double)x.size();

		return -(p / 2.0) * log(2.0 * M_PI) - (1.0 / 2.0) * LogDetSigma - 0.5 * (x - mu).transpose() * InvSigma * (x - mu);
	}
	else
	{
		return 0.0;
	}
}

void HeteroGauss::LogPdfMultiT(double &result, bool &error, Eigen::VectorXd &x, Eigen::VectorXd &mu, Eigen::MatrixXd &Sigma, int df)
{
	/*
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PURPOSE:     Density for the k-variate student-t distribution with mean mu, covariance matrix
	%              Sigma and df degrees of freedom
	%
	% INPUT:       mu       (k-by-1)        
	%              Sigma    (k-by-k)            
	%              df       (Scalar)            
	%
	% AUTHOR:      Mattias Villani, Research Department, Sveriges Riksbank and
	%              Department of Statistics, Stockholm University. 
	%
	% C++ AUTHOR:  Anders Eklund, Linköping university
	%
	% REVISED:     2016-01-19
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	*/

	if ((x.size() == 0) || (mu.size() == 0))
	{ 
		result = 0.0;
		error = false;
	}
	else
	{
		// Cholesky
		Eigen::LLT<Eigen::MatrixXd> lltOfSigma(Sigma);

		if (lltOfSigma.info() != Eigen::Success)
		{
			error = true;
			result = nan("");
		}
		else
		{
			error = false;

			double dfd = (double)df;
			double k = (double)mu.size();

			Eigen::MatrixXd P = Sigma.inverse() / (dfd - 2.0);

			double LogCt = lgamma(dfd / 2.0) - lgamma((dfd + k) / 2.0) + (k / 2.0) * log(M_PI) -0.5 * log(P.determinant());
			result = -LogCt - ((dfd + k) / 2.0) * log(1.0 + (x - mu).transpose() * P * (x - mu));
		}
	}
}



Eigen::VectorXd HeteroGauss::invLinkEval(Eigen::VectorXd &param, Eigen::MatrixXd &X, int linkType)
{
	// 
	// PURPOSE
	// --------------
	// Evaluates the inverse link function phi = g^(-1)(eta), where eta = X*param is the vector of linear predictors.
	//
	// CALL
	// --------------
	// phi = LinkEval(param,X,linkType)
	//
	// INPUTS
	// --------------
	// param             p-by-1      Vector with regression coeffs
	// X                 n-by-p      Covariate matrix
	// linkType          string      Name of the link function. 
	//                               Options: identity, log, logit, loglog, cloglog (complimentary log-log)
	//
	// OUTPUTS
	// ---------------
	// phi               n-by-1      Evaluated inverse link function for each observation, 
	//                               i.e. computes the feature from the linear predictor.
	//
	// ORIGINAL AUTHOR
	// ---------------
	// Mattias Villani, Sveriges Riksbank and Stockholm University. e-mail: mattias.villani@riksbank.se
	//
	// C++ AUTHOR
	// ---------------
	// Anders Eklund, Linköping university
	//
	// VERSION DATING
	// ---------------
	// FIRST     2009-01-23
	// CURRENT   2016-01-12
	//
	//
	//

	int n = X.rows();

	Eigen::VectorXd ones;
	ones.resize(n);
	ones.setOnes(n);

	Eigen::VectorXd linPred = X * param; // n-by-1

	switch (linkType)
	{
		case 0: // identity, g(phi) = phi
			break;
		case 1: // log, g(phi) = log(phi)
			// phi = exp(linPred);
			linPred = linPred.array().exp();
			break;
		case 2: // logit, g(phi) = log[phi/(1-phi)]
			// phi = exp(linPred)./(1 + exp(linPred));
			linPred = linPred.array().exp().cwiseQuotient( ones.array() + linPred.array().exp() );
			break;
		case 3: // loglog, g(phi) = log[-log(phi)]
			//phi = exp(-exp(linPred));
			linPred = -1.0 * linPred.array().exp();
			linPred = linPred.array().exp();
			break;
		case 4: // comploglog, g(phi) = log[-log(1-phi)]
			//phi = 1-exp(-exp(linPred));
			linPred = -1.0 * linPred.array().exp();
			linPred = linPred.array().exp();
			linPred = ones - linPred;
			break;
		case 5: // reciprocal, g(phi) = 1/phi
			//phi = 1./linPred;
			linPred = linPred.array().inverse();
			break;
		case 6: // log1, g(phi) = log(phi-1)
			//phi = 1 + exp(linPred);
			linPred = ones.array() + linPred.array().exp();
  			break;
	
		default:
			return linPred;
	}		

	return linPred;
}



void HeteroGauss::LinkDeriv(Eigen::VectorXd &firstDeriv, Eigen::VectorXd &secondDeriv, Eigen::VectorXd &feat, int linkType, bool second)
{
	// 
	// PURPOSE
	// --------------
	// Computes the first and second derivative of several common link functions
	//
	// CALL
	// --------------
	// LinkDeriv(firstDeriv, secondDeriv, feat, linkType, second)
	//
	// INPUTS
	// --------------
	// feat               n-by-1      Vector of parameters, one for each observation.
	// linkType          string      Name of the link function. 
	//                               Options: log, logit, loglog, cloglog, reciprocal or number (power function).
	// second            true, false       If true, second derivative is also computed.
	//
	// OUTPUTS
	// ---------------
	// firstDeriv        n-by-1      The gradients for each observation
	// secondDeriv       n-by-1      The Hessians for each observation
	//	
	// AUTHOR
	// ---------------
	// Mattias Villani, Sveriges Riksbank and Stockholm University. e-mail: mattias.villani@gmail.com
	//
	// C++ AUTHOR
	// ---------------
	// Anders Eklund, Linköping university
	//
	// VERSION DATING
	// ---------------
	// FIRST     2009-01-21
	// CURRENT   2016-01-21
	//

	int n = feat.rows();

	Eigen::VectorXd ones;
	ones.resize(n);
	ones.setOnes(n);

	// First derivative
	switch (linkType)
	{
	    case 0: // 'identity', g(feat) = feat
	        firstDeriv.setOnes(feat.rows(),feat.cols());
			break;
	    case 1: // 'log'   , g(feat) = log(feat)
			//firstDeriv = 1./feat;
			firstDeriv = feat.array().inverse();
			break;
	    case 2: // 'logit'  , g(feat) = log[feat/(1-feat)]
	        //firstDeriv = 1./(feat.*(1-feat));						
			firstDeriv = ( feat.array() * (ones - feat).array() ).array().inverse();
			break;
	    case 3: // 'loglog' , g(feat) = log[-log(feat)]
    	    //firstDeriv = 1./(log(feat).*feat);
			firstDeriv = ( feat.array().log() * feat.array() ).array().inverse();
			break;
	    case 4: // 'comploglog'  g(feat) = log[-log(1-feat)]
	        //firstDeriv = -1./((1-feat).*log(1-feat));
			firstDeriv = -1.0 * ( (ones - feat).array() * (ones - feat).array().log() ).array().inverse();
			break;
	    case 5: // 'reciprocal'
    	    //firstDeriv = -1./feat.^2;
			firstDeriv = -1.0 * (feat.array() * feat.array()).array().inverse();
			break;

		default:
			firstDeriv = feat;

	}

	// Second derivative
	if (second)
	{
	    switch (linkType)
		{
    	    case 0: // 'identity'
    	        secondDeriv.setZero(feat.rows(),feat.cols());
				break;
    	    case 1: // 'log'
    	        //secondDeriv = -1./(feat.^2);
				secondDeriv = -1.0 * ( feat.array() * feat.array() ).array().inverse();
				break;
    	    case 2: // 'logit'
    	        //secondDeriv = (2*feat-1)./(feat.*(1-feat)).^2;				
				secondDeriv = (2.0 * feat.array() - ones.array()).cwiseQuotient( ( feat.array() * (ones - feat).array() ).array().square() );
				break;
    	    case 3: // 'loglog'
    	        //secondDeriv = -(1+log(feat))./(log(feat).*feat).^2;
				secondDeriv = -1.0 * (ones.array() + feat.array().log()).cwiseQuotient( (feat.array().log() * feat.array()).array().square() );
				break;
    	    case 4: // 'comploglog'
    	        //secondDeriv = -(1+log(1-feat))./((1-feat).*log(1-feat)).^2;
				secondDeriv = -1.0 * (ones.array() + (ones - feat).array().log()).cwiseQuotient( ((ones - feat).array() * (ones - feat).array().log()).array().square() );
				break;
    	    case 5: // 'reciprocal'
    	        //secondDeriv = 2./feat.^3;
				secondDeriv = 2.0 * (feat.array() * feat.array() * feat.array()).array().inverse();	
				break;

			default:
				secondDeriv = feat;    	 
		}
	}
	else
	{
	    secondDeriv.setZero(feat.rows(),feat.cols());
	}
}

void HeteroGauss::HeteroGaussGradHess(Eigen::VectorXd &gradObs, 
                                      Eigen::VectorXd &hessObs, 
                                      Eigen::VectorXd &param, 
                                      Eigen::VectorXd &y, 
                                      Eigen::MatrixXd &V, 
                                      int linkType, 
                                      int condDens, 
                                      Eigen::MatrixXd &feat, 
                                      int hessMethod)
{
	// 
	// PURPOSE
	// --------------
	// Computes the gradient and Hessian for a negative binomial model.
	//
	// CALL
	// --------------
	// HeteroGaussGradHess(gradObs, hessObs, param,y,V,linkType,condDens,feat,hessMethod)
	//
	// INPUTS
	// --------------
	// param             p-by-1          Vector with regression parameters in the current feature.
	// y                 n-by-1          Vector with response observations.
	// V                 n-by-p          Matrix with covariates for the current feature
	// linkType          int             Name of the link
	// condDens          integer         The index of the current feature
	// feat              n-by-nFeat      feat(:,j) is an n-by-1 vector of evaluations of the j:th feature (e.g. a mean)
	// hessMethod        integer         0: Expected Hessian, 1: Outer product of gradients, 2: Observed Hessian
	//
	// OUTPUTS
	// ---------------
	// gradObs
	// hessObs
	//
	// ORIGINAL AUTHOR
	// ---------------
	// Mattias Villani, Sveriges Riksbank and Stockholm University. 
	//
	// C++ AUTHOR
	// ---------------
	// Anders Eklund, Linköping university
	//
	// VERSION DATING
	// ---------------
	// FIRST     2009-01-27
	// CURRENT   2016-01-21
	//
	// REFERENCES
	// ---------------
	// Kohn, R. and Villani, M. (2009). 
	// A General Approach to Regression Density Estimation using Smooth Mixtures of Over-Dispersed Models.
	//	

	Eigen::VectorXd mu = feat.col(0);  // mean
	Eigen::VectorXd	phi = feat.col(1); // variance

	switch (condDens)
	{
        case 1: // Full conditional density for mu given phi - compute only mu
		
	        mu = invLinkEval(param,V,linkType);       
        
			// Computing gradient and Hessian with respect to mu for each observation
        	gradObs = (y - mu).cwiseQuotient(phi); 
	
	        switch (hessMethod)
			{
	            case 0: // Expected Hessian
					//hessObs = -1./phi;
					hessObs = -1.0 * phi.array().inverse();
					break;
    	        case 1: // Outer product of gradient
    	            // Do nothing, hessian is not needed
					hessObs.setZero(phi.rows());					
					break;
    	        case 2: // Observed Hessian
					//hessObs = -1./phi;
					hessObs = -1.0 * phi.array().inverse();
					break;

				default: ;
			}
			break;
    	    
	    case 2: // Full conditional density for phi given mu - compute only phi
        	
			phi = invLinkEval(param,V,linkType);  
        
			//gradObs = -1./(2*phi) + ((y-mu).^2)./(2*phi.^2);
			gradObs = -1.0 * (2.0 * phi).array().inverse() + ( (y - mu).array().square() ).cwiseQuotient( 2.0 * phi.array().square() );

	        switch (hessMethod)
			{
	            case 0: // Expected Hessian
					//hessObs = -1./(2*(phi.^2));
					hessObs = -1.0 * ( 2.0 * phi.array().square() ).array().inverse();
					break;
	            case 1: // Outer product of gradient
                	// Do nothing, hessian is not needed
					hessObs.setZero(phi.rows());
					break;
   		        case 2: // Observed Hessian
					//hessObs = 1./(2*(phi.^2)) - ((y-mu).^2)./(phi.^3);
					hessObs = ( 2.0 * phi.array().square() ).array().inverse() - ( (y - mu).array().square() ).cwiseQuotient( phi.array().cube() );
					break;

				default: ;
			}
			break;

		default: ;
	}
}


Eigen::VectorXd HeteroGauss::HeteroGaussLogPost(Eigen::VectorXd &param, 
                                                Eigen::VectorXd &y, 
                                                Eigen::MatrixXd &V, 
                                                Eigen::VectorXd &priorMean, 
                                                Eigen::MatrixXd &priorCov, 
                                                int linkType, 
                                                int condDens, 
                                                Eigen::MatrixXd &feat)
{
	// 
	// PURPOSE
	// --------------
	// Computes the log likelihood/posterior for a Gaussian model. This
	// function does different things and returns different objects depending on condDens
	//
	// CALL
	// --------------
	// logPost = NegBinLogPost(param,y,V,priorMean,priorCov,linkType,condDens,feat)
	//
	// INPUTS
	// --------------
	// param             empty / p-by-1 / nFeat-by-1 cell    Vector with regression parameters.
	// y                 empty / n-by-1 / n-by-1             Vector with response observations.
	// V                 empty / n-ny-p / nFeat-by-1 cell    Matrix with covariates
	// priorMean         empty / empty / p-by-1              Prior mean on the parameters conditional on inclusion.
	// priorCov          empty / empty / p-by-p              Prior precision matrix for the parameters
	// linkType          empty / string / nFeat-by-1 cell    Name of the link
	// condDens          Integer                             condDens = 0, evaluates the log-likelihood for each observation 
	//                                                        for given features. This is used when updating the component allocation.
	//                                                       condDens = i,nfunction used to compute the full conditional for the j:th feature.
	//                                                       condDens = NaN, function used for optimization simultaneously over all features.
	// feat              n-by-nFeat      feat(:,j) is an n-by-1 vector of evaluations of the j:th feature (e.g. a mean)
	//
	// OUTPUTS
	// ---------------
	// logPost           Scalar/Vector   The output depends on how this function was called:
	//                                   condDens = 0 - The log-likelihood of each obs
	//                                   condDens = i, i=1,2,... Returns the log-likelihood summed across all observations.
	//                                   condDens = NaN, Returns the negative log posterior summed across all observations.
	//
	// AUTHOR
	// ---------------
	// Mattias Villani, Sveriges Riksbank and Stockholm University. 	
	//
	// C++ AUTHOR
	// ---------------
	// Anders Eklund, Linköping University. 	
	//
	//
	// VERSION DATING
	// ---------------
	// FIRST     2009-02-27
	// CURRENT   2016-01-12
	//
	// REFERENCES
	// ---------------
	// Kohn, R. and Villani, M. (2009). 
	// A General Approach to Regression Density Estimation using Smooth Mixtures of Over-Dispersed Models.
	//
	
	Eigen::VectorXd mu;
	Eigen::VectorXd phi;

	// Computing the features
	switch (condDens)
	{
	    case 0: // condDens = 0 means the function is used for updating the component allocation. Do not update features.
	        mu = feat.col(0);
	        phi = feat.col(1);
			break;
	    case 1: // Updating the parameters in the first feature
        	mu = invLinkEval(param,V,linkType);
        	phi = feat.col(1);
			break;
	    case 2: // Updating the parameters in the second feature
	        mu = feat.col(0);
	        phi = invLinkEval(param,V,linkType);
			break;
		
		default: ; // We are optimizing the log posterior - param is a vector with all parameters stacked.         
    	
	        //nCovMean = size(V{1},2);
	        //mu = invLinkEval(param.head(nCovMean),V{1},linkType{1});
	        //phi = invLinkEval(param.tail(nCovMean+1:end),V{2},linkType{2});
	}

	// Log likelihood - This is the full likelihood, including normalization constants.
	//logLik = -0.5*log(2*pi*phi) - (1./(2*phi)).*((y-mu).^2)
	Eigen::VectorXd logLik = -0.5 * (2.0 * M_PI * phi).array().log() - ( (2.0 * phi).array().inverse() ).array() * ( (y - mu).array().square() ).array();

	if (condDens == 0)
	{
		return logLik;
	}
	else
	{
		double sum = logLik.sum();
		Eigen::VectorXd result;
		result.resize(1);
		result(0) = sum;
		return result;		
	}

	// Log Prior
	/*
	if (condDens == nanl)  // Computing the log posterior, here we are optimizing. Negative log post is returned.
	{
	    logPrior = LogPdfMultiN(param.head(nCovMean),priorMean{1},priorCov{1}) + LogPdfMultiN(param(nCovMean+1:end),priorMean{2},priorCov{2});
	    logPost = -(logLik + logPrior);
	}
	else
	{
	    logPost = logLik;
	}
	*/
}



void HeteroGauss::NewtonGradHess(Eigen::VectorXd &Grad, 
                                 Eigen::MatrixXd &Hesspp, 
                                 Eigen::MatrixXd &Hesspc, 
                                 Eigen::VectorXd &Beta, 
                                 Eigen::VectorXd &Ic, 
                                 Eigen::VectorXd &Ip, 
                                 Eigen::VectorXd &mu, 
                                 Eigen::MatrixXd &invPsi, 
                                 int hessMethod, 
                                 int gradHessFunc, 
                                 Eigen::VectorXd &y, 
                                 Eigen::MatrixXd &V, 
                                 int linkType, 
                                 int condDens, 
                                 Eigen::MatrixXd &feat)
{	
	// 
	// PURPOSE
	// --------------
	// Computes the Gradient and Hessian of a GLM-type of model with a menu of link functions.
	// This function is called from NewtonIter (generalized Newton iterations)
	// 
	// Beta(k+1) = inv(Hesspp)*(Hesspc*Beta(k) - Grad)
	//
	// CALL
	// --------------
	// NewtonGradHess(Grad, Hesspp, Hesspc,Beta_,Ic,Ip,mu,invPsi,hessMethod,gradHessFunc,y,V,linkType,condDens,feat)
	//
	//
	// INPUTS
	// --------------
	// Beta_             q-by-1      Parameter vector where the gradient and Hessian are evaluated.
	// Ic                q-by-1      Logical vector with variable selection indicator for the current draw.
	// Ip                q-by-1      Logical vector with variable selection indicator for the proposed draw.
	// mu                q-by-1      Prior mean of regression coefficient (Prior: beta ~ N(mu,Psi)
	// invPsi           	q-by-q      Inverse of Prior covariance matrix for the regression coefficients
	// hessMethod        Integer     0: Expected Hessian, 1: Outer product of gradients, 2: Observed Hessian
	// gradHessFunc      string      Name of the function that computes the gradient and Hessian for each observation.
	// y                 n-by-1      Vector of responses
	// V                 n-by-q      Covariate matrix, including column of ones for the intercept
	// linkType          string      Name of the link function, see LinkDeriv.m for valid options.
	// condDens          Integer     Index for the full conditional of the current feature
	// feat              n-by-nFeat  Matrix with features

	//
	// OUTPUTS
	// ---------------
	// Grad              q2-by-1     Gradient of the log pdf at the posterior mode
	// Hesspp            q2-by-q2    Hessian matrix of the log pdf at the posterior mode.
	// Hesspc            q2-by-q1    Hessian matrix of the log pdf at the posterior mode.
	//
	// AUTHOR
	// ---------------
	// Mattias Villani, Sveriges Riksbank and Stockholm University. e-mail: mattias.villani@riksbank.se
	//
	// C++ AUTHOR
	// ---------------
	// Anders Eklund, Linköping university
	//
	// VERSION DATING
	// ---------------
	// FIRST     2008-12-02
	// CURRENT   2016-01-22
	//
	// REFERENCES
	// ---------------
	// Villani, M., Kohn, R. and Giordani, P. (2008). 
	// Regression Density Estimation using Smooth Adaptive Gaussian Mixtures.
	//
	//
	
	// Evaluating the gradient and hessian for each observation

	// Get subsets based on indicator vectors
	int nCurrentIndicators = (int)Ic.sum();
	int nProposedIndicators = (int)Ip.sum();

	Eigen::VectorXd BetaSubsetIc, BetaSubsetIp, muSubsetIp;
	Eigen::MatrixXd VSubsetIc, VSubsetIp, invPsiSubsetIp, invPsiSubsetIc, invPsiSubsetIpIp, invPsiSubsetIpIc;

	BetaSubsetIc.resize(nCurrentIndicators);
	BetaSubsetIp.resize(nProposedIndicators);

	muSubsetIp.resize(nProposedIndicators);
	
	VSubsetIc.resize(V.rows(),nCurrentIndicators);
	VSubsetIp.resize(V.rows(),nProposedIndicators);

	invPsiSubsetIp.resize(Beta.size(),nProposedIndicators);
	invPsiSubsetIc.resize(Beta.size(),nCurrentIndicators);
	invPsiSubsetIpIp.resize(nProposedIndicators,nProposedIndicators);
	invPsiSubsetIpIc.resize(nProposedIndicators,nCurrentIndicators);

	// There is probably a smarter way to do this...

	// Ic
	int current = 0;
	for (int i = 0; i < Ic.size(); i++)
	{
		if (Ic(i) == 1.0)
		{
			BetaSubsetIc(current) = Beta(i);
			VSubsetIc.col(current) = V.col(i);
			invPsiSubsetIc.col(current) = invPsi.col(i);
			current++;
		}
	}

	// Ip
	current = 0;
	for (int i = 0; i < Ip.size(); i++)
	{
		if (Ip(i) == 1.0)
		{
			BetaSubsetIp(current) = Beta(i);
			muSubsetIp(current) = mu(i);
			VSubsetIp.col(current) = V.col(i);
			invPsiSubsetIp.col(current) = invPsi.col(i);
			current++;
		}
	}

	current = 0;
	for (int i = 0; i < Ip.size(); i++)
	{
		if (Ip(i) == 1.0)
		{
			invPsiSubsetIpIp.row(current) = invPsiSubsetIp.row(i);
			invPsiSubsetIpIc.row(current) = invPsiSubsetIc.row(i);
			current++;	
		}	
	}

	//featC = invLinkEval(Beta_(Ic),V(:,Ic),linkType); 
	Eigen::VectorXd featC = invLinkEval(BetaSubsetIc,VSubsetIc,linkType); // The current feature

	//[gradObs, hessObs] = feval(gradHessFunc,Beta_(Ic),y,V(:,Ic),linkType,condDens,feat,hessMethod); // grad is n-by-1 and so is hess        
	Eigen::VectorXd gradObs, hessObs;
	HeteroGaussGradHess(gradObs, hessObs, BetaSubsetIc, y, VSubsetIc, linkType, condDens, feat, hessMethod);

	// Forming the gradient vector ...
	Eigen::VectorXd kPrime, dummy;
	LinkDeriv(kPrime,dummy,featC,linkType,false);

	//d = gradObs./kPrime;
	Eigen::VectorXd d = gradObs.cwiseQuotient(kPrime);
	//Grad = V(:,Ip)'*d - invPsi(Ip,Ip)*(Beta_(Ip)-mu(Ip));
	Grad = VSubsetIp.transpose() * d - invPsiSubsetIpIp * (BetaSubsetIp - muSubsetIp);
 
	// ... and the Hessian matrix
	Eigen::VectorXd diagElements;
	Eigen::VectorXd kPrimeInverse;
	Eigen::VectorXd kBiss;

	switch (hessMethod)
	{
	    case 0: // Expected Hessian 
    	    //diagElements  = hessObs./kPrime.^2;
			diagElements = ( hessObs.array() ).cwiseQuotient( kPrime.array().square() );
			break;
	    case 1: // Outer product of gradient
	        //diagElements  = -d.^2;
			diagElements = -1.0 * d.array().square();
			break;
	    case 2: // Observed Hessian
	        //[~, kBiss] = LinkDeriv(1./kPrime,linkType,1);
    	    //diagElements = hessObs./kPrime.^2   -(gradObs.*kBiss)./kPrime.^2; % d1 + d2			
			kPrimeInverse = kPrime.array().inverse();			
			LinkDeriv(dummy,kBiss,kPrimeInverse,linkType,true);
			diagElements = ( hessObs.array() ).cwiseQuotient( kPrime.array().square() ) - ( gradObs.array() * kBiss.array() ).cwiseQuotient( kPrime.array().square() );
			break;

		default : ;
	}

	Eigen::MatrixXd diagMatrix = diagElements.asDiagonal();

	//Hesspp = MultXDiagX(V(:,Ip),V(:,Ip),diagElements) - invPsi(Ip,Ip); // MultXDiagX computes X1'*diag(w)*X2 fast.
	//Hesspc = MultXDiagX(V(:,Ip),V(:,Ic),diagElements) - invPsi(Ip,Ic);

	Hesspp = VSubsetIp.transpose() * diagMatrix * VSubsetIp - invPsiSubsetIpIp;
	Hesspc = VSubsetIp.transpose() * diagMatrix * VSubsetIc - invPsiSubsetIpIc;
}

Eigen::VectorXd HeteroGauss::GetSubsetVector(Eigen::VectorXd &vector, Eigen::VectorXd &indicators)
{
	int subsetSize = (int)indicators.sum();
	Eigen::VectorXd subset;
	subset.resize(subsetSize);

	int current = 0;
	for (int i = 0; i < indicators.size(); i++)
	{
		if (indicators(i) == 1.0)
		{
			subset(current) = vector(i);
			current++;
		}
	}

	return subset;
}

Eigen::MatrixXd HeteroGauss::GetSubsetMatrix(Eigen::MatrixXd &matrix, Eigen::VectorXd &rowIndicators, Eigen::VectorXd &colIndicators)
{
	int rows = (int)rowIndicators.sum();
	int cols = (int)colIndicators.sum();

	Eigen::MatrixXd subsetTemp, subset;

	subsetTemp.resize(matrix.rows(),cols);  // All rows but only indicator columns
	subset.resize(rows,cols);				// Indicator rows and indicator columns

	// Pick out correct columns
	int current = 0;
	for (int i = 0; i < colIndicators.size(); i++)
	{
		if (colIndicators(i) == 1.0)
		{
			subsetTemp.col(current) = matrix.col(i);
			current++;
		}
	}

	// Pick out correct rows from column subset
	current = 0;
	for (int i = 0; i < rowIndicators.size(); i++)
	{
		if (rowIndicators(i) == 1.0)
		{
			subset.row(current) = subsetTemp.row(i);
			current++;	
		}	
	}

	return subset;
}


Eigen::VectorXd HeteroGauss::SaveFullVector(Eigen::VectorXd &vector, Eigen::VectorXd &indicators)
{
	int fullSize = indicators.size();
	Eigen::VectorXd full;
	full.resize(fullSize);

	int current = 0;
	for (int i = 0; i < indicators.size(); i++)
	{
		if (indicators(i) == 1.0)
		{
			full(i) = vector(current);
			current++;
		}
		else
		{
			full(i) = 0.0;
		}
	}

	return full;
}

Eigen::VectorXd HeteroGauss::FlipIndicators(Eigen::VectorXd &indicators)
{
	Eigen::VectorXd flipped = indicators;

	for (int i = 0; i < indicators.size(); i++)
	{
		flipped(i) = 1.0 - indicators(i);
	}

	return flipped;
}


void HeteroGauss::NewtonIter(Eigen::VectorXd &betaEnd, 
                             Eigen::MatrixXd &hess, 
                             Eigen::VectorXd &betaC, 
                             Eigen::VectorXd &Ic, 
                             Eigen::VectorXd &Ip, 
                             Eigen::VectorXd &mu, 
                             Eigen::MatrixXd &invPsi, 
                             int nSteps, 
                             int hessMethod, 
                             int gradHessFunc, 
                             Eigen::VectorXd &y, 
                             Eigen::MatrixXd &V, 
                             int linkType, 
                             int condDens, 
                             Eigen::MatrixXd &feat)
{
	// 
	// PURPOSE
	// --------------
	// Iterates the generalized (dimension-changing) Newton algorithm. 
	// This function is called from NewtonProp.
	//
	// CALL
	// --------------
	// NewtonIter(betaEnd,hess,betaC,Ic,Ip,mu,invPsi,nSteps,hessMethod,gradHessFunc,y,V,linkType,condDens,feat)
	//
	// INPUTS
	// --------------
	// betaC             nPara-by-1      Starting point of the Newton algorithm.
	// Ic                nPara-by-1      Logical vector with variable selection indicator for the current draw.    
	// Ip                nPara-by-1      Logical vector with variable selection indicator for the proposed draw.
	// mu                nPara-by-1      Prior mean on the parameters conditional on inclusion.
	// invPsi            nPara-by-nPara  Prior precision matrix for the  parameters
	// nSteps            Scalar          The number of iterations of the Newton algorithm. nSteps = Inf is allowed.
	// hessMethod        Integer         0: Expected Hessian, 1: Outer product of gradients, 2: Observed Hessian
	// gradHessFunc      int             The name of the function that computes the gradient and Hessian of the target density.
	// y                 n-by-1          Vector with response observations.
	// V                 n-by-p          Matrix with covariates for the current feature
	// linkType          int             Name of the link
	// condDens          Integer         The index of the current feature
	// feat              n-by-nFeat      feat(:,j) is an n-by-1 vector of evaluations of the j:th feature (e.g. a mean)
    //
	//
	// OUTPUTS
	// ---------------
	// betaEnd           nPara-by-1      The terminal vector of the Newton algorithm.
	// hess              nPara-by-nPara  Hessian of the target density at betaEnd
	// ErrorFlag
	//
	// AUTHOR
	// ---------------
	// Mattias Villani, Sveriges Riksbank and Stockholm University. e-mail: mattias.villani@gmail.com
	//
	// C++ AUTHOR
	// ---------------
	// Anders Eklund, Linköping university
	//
	// VERSION DATING
	// ---------------
	// FIRST     2008-11-20
	// CURRENT   2016-01-22
	//
	// REFERENCES
	// ---------------
	// Villani, M., Kohn, R. and Giordani, P. (2008). Regression Density Estimation using Smooth Adaptive Gaussian Mixtures.
	//

	Eigen::MatrixXd dummymatrix;

	// Prelims
	int nCurrentIndicators = (int)Ic.sum();		
	int nSteps_ = nSteps;

	if (nCurrentIndicators == 0)     // current draw is from zero-dim model (no parameter). Treat this like a dimension-preserving move,
	{              					 // but do twice as many steps because initial value is no good. 
		Ic = Ip;
	    betaC = Eigen::MatrixXd::Zero(Ip.size(),1); // 0.1*randn(length(Ip),1); 
    	nSteps_ *= 2;
	}

	// Finite pre-determined number of Newton steps.

	Eigen::VectorXd grad;
	Eigen::MatrixXd hessPP, hessPC;

	betaEnd	= betaC;

	Eigen::VectorXd betaEndTemp;

	for (int i = 0; i < nSteps_; i++)
	{
		bool allequal = true;
		for (int i = 0; i < Ic.size(); i++)
		{
			if (Ic(i) != Ip(i))
			{
				allequal = false;
				break;
			}
		}		

	    if (!allequal) // Dimension changing move
		{
	        NewtonGradHess(grad,hessPP,hessPC,betaEnd,Ic,Ip,mu,invPsi,hessMethod,gradHessFunc,y,V,linkType,condDens,feat);
    	    //betaEndTemp = (hessPP\hessPC) * betaC(Ic) - hessPP \ grad; // Compact size
			betaEndTemp = ( hessPP.lu().solve(hessPC) ) * GetSubsetVector(betaC,Ic) - hessPP.lu().solve(grad);
			betaEnd = SaveFullVector(betaEndTemp,Ip); // Full size
	        Ic = Ip;
		}      
	    else // Dimension preserving move
		{
	        NewtonGradHess(grad,hess,dummymatrix,betaEnd,Ic,Ip,mu,invPsi,hessMethod,gradHessFunc,y,V,linkType,condDens,feat);
    	    //betaEndTemp = betaC(Ip) - hess\grad;
			betaEndTemp = GetSubsetVector(betaC,Ip) - hess.lu().solve(grad);
			betaEnd = SaveFullVector(betaEndTemp,Ip); // Full size
    	}
    
	    betaC = betaEnd;
	}
    
	// Compute the grad and Hessian at betaEnd, the terminal point.
	NewtonGradHess(grad,hess,dummymatrix,betaEnd,Ip,Ip,mu,invPsi,hessMethod,gradHessFunc,y,V,linkType,condDens,feat);
        
	betaEnd = GetSubsetVector(betaEnd,Ip); // Compact form

	/*	
	if any(any(isnan(hess)))
	{
    	ErrorFlag = 1;
	}
	else
	{
	    ErrorFlag = 0;
	}
	*/		
}



void HeteroGauss::NewtonProp(Eigen::VectorXd &paramCurr, 
                             Eigen::VectorXd &ICurr, 
                             double &AccPr, 
                             bool &errFlagProp, 
                             bool &errFlagRev, 
                             Eigen::VectorXd &PrIn, 
                             Eigen::VectorXd &mu, 
                             Eigen::MatrixXd &invPsi, 
                             Eigen::VectorXd &onTrial, 
                             int logPostFunc, 
                             int gradHessFunc, 
                             Eigen::VectorXd &y, 
                             Eigen::MatrixXd &V, 
                             int linkType, 
                             int condDens, 
                             Eigen::MatrixXd &feat, 
                             int nSteps, 
                             int propDf, 
                             double IUpdatePr, 
                             int hessMethod)
{
	// 
	// PURPOSE
	// --------------
	// Performs a full update of the variable selection indicators and the parameters using variable dimension Newton proposals.
	//
	// CALL
	// --------------
	// NewtonProp(paramCurr,ICurr,PrIn,mu,invPsi,onTrial,logPostFunc,gradHessFunc,y,V,linkType,condDens,feat,nSteps,propDf,IUpdatePr,hessMethod,S,paramCurrAllComp,refComp,currComp)
	//        
	// INPUTS
	// --------------
	// paramCurr         nPara-by-1      Parameter vector for current draw, including zeros for the currently excluded variables.
	// ICurr             nPara-by-1      Logical vector with the current covariate selection indicators.
	// PrIn              Scalar in [0,1] Prior probability of including a covariate in the model.
	// mu                nPara-by-1      Prior mean on the parameters conditional on inclusion.
	// invPsi            nPara-by-nPara  Prior precision matrix for the parameters
	// onTrial           Vector          Vector with indicies for the covariates that are up for inclusion/exclusion.
	//                                   Note: not all of these variable selection indicators will actually be
	//                                   proposed, the have to pass a Bernoulli trial with success prob IUpdatePr.
	// gradHessFunc      String          Name of the .m function containing the gradient and hessian of the full conditional posteriors.
	// logPostFunc       String          Name of the .m function containing the full conditional PDFs
	// y                 nObs-by-1       Vector with observations on the response
	// V                 nObs-by-nPara   Covariate matrix, including a column of ones for the intercept.
	// linkType          String/Integer  Choice of link function. Options: 'identity','log','logit','loglog','comploglog' 
	//                                   or power link (if linkType is an integer).
	// condDens          Integer         Indicating which of the full conditional posteriors that is sampled in the current iteration.
	// feat              nObs-by-nFeat   feat(i,j) contains the values for the j:th feature for the ith observation at the current draw.
	// nSteps            Integer         Number of Newton steps.
	// propDf            Integer         Degrees of freedom in the student-t proposal density. Note that we require propDf>=3, 
	//                                   but this can be generalized if needed.
	// IUpdatePr         Scalar in [0,1] Probability of updating a covariate in the current iteration of the algorithm.
	// hessMethod        Integer         0: Expected Hessian, 1: Outer product of gradients, 2: Observed Hessian    
	//
	// OUTPUTS
	// ---------------
	// paramCurr         nPara-by-1      Updated Parameter vector, including zeros for the currently excluded variables.
	// ICurr             nPara-by-1      Updated Logical vector with covariate selection indicators.
	// AccPr             Scalar          Acceptance probability in the last MH-step where we only updated the parameters.
	//
	// AUTHOR
	// ---------------
	// Mattias Villani, Sveriges Riksbank and Stockholm University. e-mail: mattias.villani@gmail.com
	//
	// C++ AUTHOR
	// --------------
	// Anders Eklund, Linköping university
	//
	// VERSION DATING
	// ---------------
	// FIRST     2008-11-20
	// CURRENT   2016-01-25
	//
	// REFERENCES
	// ---------------
	// Kohn, R. and Villani, M. (2009). A General Approach to Regression Density Estimation.
	// Villani, M., Kohn, R. and Giordani, P. (2008). Regression Density Estimation using Smooth Adaptive Gaussian Mixtures.
	//
	
	errFlagProp = false;
   	errFlagRev = false;	

	double AccPrRaw;

	Eigen::VectorXd propMean, revMean;
	Eigen::MatrixXd propHess, revHess;

	Eigen::VectorXd draw;

	// Prelims
	Eigen::MatrixXd Psi = invPsi.inverse(); // TODO remove this when mvnpdf is replaced with more efficient function

	// Randomly selecting the subset of covariates that we allow to change in the current iteration.

	int nOnTrial = onTrial.size(); // Number of covariates that may be in or out.
	
	Eigen::VectorXd IUpdate, randomElements;
	randomElements.resize(nOnTrial);
	IUpdate.resize(0);

	bool allsmaller = true;
	for (int i = 0; i < PrIn.size(); i++)
	{
		if (PrIn(i) >= 1.0)
		{
			allsmaller = false;
			break;
		}
	}

	//std::cout << "Before all smaller check, IUpdate is " << std::endl << IUpdate << std::endl << std::endl;

	if (allsmaller)
	{
	    //IUpdate = onTrial(rand(1,nOnTrial) < IUpdatePr); // These are the covariates that we update at THIS iteration.

		for (int i = 0; i < nOnTrial; i++)
		{
			double uniformNumber = uniformDistribution(generator);
			if (uniformNumber < IUpdatePr)
			{
				randomElements(i) = 1.0;
			}
			else
			{
				randomElements(i) = 0.0;
			}
		}

		IUpdate = GetSubsetVector(onTrial,randomElements);
	}
	else // No variable selection
	{
	}

	//std::cout << "After all smaller check, IUpdate is " << std::endl << IUpdate << std::endl << std::endl;

	//IUpdate=[IUpdate 0]; // Adding artificial zero for the final step where only the parameters are updated, but not the indicators.
	Eigen::VectorXd oldIUpdate = IUpdate;
	IUpdate.resize(IUpdate.size() + 1);
	IUpdate.head(IUpdate.size() - 1) = oldIUpdate;
	IUpdate(IUpdate.size() - 1) = -1.0; // -1 since C uses 0:n-1, while Matlab uses 1:n

	Eigen::VectorXd IProp;
	Eigen::VectorXd muCond;

	double gProp, gCurr;
	double priorIProp, priorICurr;
	double LogPropProp;
	double LogPriorCurr, logPriorProp, LogPropCurr, logPriorCurr, logPostProp, logPostCurr;

	Eigen::VectorXd logLikProp, logLikCurr;
	
	//std::cout << "Before covariate loop, IUpdate is " << std::endl << IUpdate << std::endl << std::endl;

	// Looping over the selected covariates, jointly proposing variable indicator and parameter
	for (int i = 0; i < IUpdate.size(); i++)
	{
		int update = (int)IUpdate(i);

		// Proposing the indicator. More general than needed because we may want to use Kohn-Nott adaptive proposals later on.
	    IProp = ICurr; // Initializing the indicator proposal.
	
    	if (update > -1) // Variable selection step
		{
        	// Here we use a simple Metropolized proposal of indicators: Out -> In or In -> Out.
        	// TODO: Add adaptive scheme.
                			
        	IProp(update) = 1.0 - IProp(update); // Simple update, if covariate is in, then remove and vice versa.
       
	        gProp = 1.0;      // The proposal probability for the proposed model
	        gCurr = 1.0;      // The proposal probability for the current model
        
			Eigen::VectorXd temp;
			Eigen::VectorXd ones;
			ones.resize(ICurr.size());
			ones.setOnes(ICurr.size());
			
			temp = PrIn.array() * IProp.array() + (ones - PrIn).array() * (ones - IProp).array(); // Prior probability of the proposed model
			priorIProp = temp.prod();
			//std::cout << "prior I prop is " << std::endl << priorIProp << std::endl << std::endl;

			temp = PrIn.array() * ICurr.array() + (ones - PrIn).array() * (ones - ICurr).array(); // Prior probability of the current model
	        priorICurr = temp.prod(); 

			//std::cout << "prior I curr is " << std::endl << priorICurr << std::endl << std::endl;
		}
	    else
		{
        	// Extra step. Only update parameters.
	        priorIProp = 1.0; priorICurr = 1.0; gProp = 1.0; gCurr = 1.0;
		}
    
		// Propose the parameters conditional on indicators using the finite step Newton method
	    Eigen::VectorXd paramProp;
		paramProp.resize(paramCurr.size());
		paramProp.setZero(paramCurr.size()); 
		
		bool allEqual = true;
		for (int i = 0; i < IProp.size(); i++)
		{
			if (IProp(i) != ICurr(i))
			{
				allEqual = false;
				break;
			}
		}

		bool allZero = true;
		for (int i = 0; i < IProp.size(); i++)
		{
			if (IProp(i) != 0.0)
			{
				allZero = false;
				break;
			}
		}

	    if ( ((update > -1) && allEqual) || ((update == -1) && allZero) )
		{
	        // Proposed to stay - do not update parameters (cannot happen with the Metropolized proposals here, but more generally it can) or
        	// Only update the parametes (update==0), but all parameters are out.
	        AccPr = nan("");
		}
	    else
		{
        	// Proposed a move - do update parameters.
        
			bool anyIProp = false;
			for (int i = 0; i < IProp.size(); i++)
			{
				if (IProp(i) != 0.0)
				{
					anyIProp = true;
					break;
				}
			}

        	if (anyIProp) // At least one covariate included in the proposal draw. Do Newton iterations and propose.
			{
            	// Finite-step Newton: paramCurr -> (propMean,propHess),to get the mean and covariance matrix of the proposal.
				Eigen::VectorXd paramCurrTemp = paramCurr;
				Eigen::VectorXd ICurrTemp = ICurr;
            	NewtonIter(propMean,propHess,paramCurrTemp,ICurrTemp,IProp,mu,invPsi,nSteps,hessMethod,gradHessFunc,y,V,linkType,condDens,feat);
            
            	// Generating the proposal draw from t(propMean,-InvHess,dfDelta), where -InvHess is the covariance matrix, NOT the scale matrix.
            	//[paramProp(IProp), errFlagProp] = RndMultiT(propMean,-inv(propHess),propDf);
				Eigen::MatrixXd propHessInverse = -1.0 * propHess.inverse();				
				RndMultiT(draw,errFlagProp,propMean,propHessInverse,propDf);
				if (!errFlagProp)
				{				
					paramProp = SaveFullVector(draw, IProp);
				}
			}          
	        else // No covariate in, no parameters to propose.
			{
	            propMean.resize(0);
    	        propHess.resize(0,0);
			}
                
        	if (!errFlagProp)
			{            
            	// Evaluating the proposal density in the proposal point
				bool dummyError;
				Eigen::MatrixXd propHessInverse = -1.0 * propHess.inverse();
				Eigen::VectorXd paramPropSubset = GetSubsetVector(paramProp,IProp);
            	LogPdfMultiT(LogPropProp,dummyError,paramPropSubset,propMean,propHessInverse,propDf);
            
            	// Computing the proposal density in the current point. Reversing the Newton.
            
            	// Reversing the Newton: paramProp -> paramRev
				bool anyICurr = false;
				for (int i = 0; i < ICurr.size(); i++)
				{
					if (ICurr(i) != 0.0)
					{
						anyICurr = true;
						break;
					}
				}

	            if (anyICurr)
				{
					Eigen::VectorXd paramPropTemp = paramProp;
					Eigen::VectorXd IPropTemp = IProp;
	                NewtonIter(revMean,revHess,paramPropTemp,IPropTemp,ICurr,mu,invPsi,nSteps,hessMethod,gradHessFunc,y,V,linkType,condDens,feat);
					Eigen::VectorXd paramCurrSubset = GetSubsetVector(paramCurr,ICurr);
					Eigen::MatrixXd revHessInverse = -1.0 * revHess.inverse();					
    	            LogPdfMultiT(LogPropCurr,errFlagRev,paramCurrSubset,revMean,revHessInverse,propDf);
				}
	            else
				{
    	            LogPropCurr = 0.0;
    			}
                        
            	if (!errFlagRev)
				{
					Eigen::VectorXd priorMeanDummy;
					Eigen::MatrixXd priorCovDummy;

					Eigen::VectorXd paramPropSubset;
					Eigen::VectorXd paramCurrSubset;

					Eigen::VectorXd NegIProp = FlipIndicators(IProp);
					Eigen::VectorXd NegICurr = FlipIndicators(ICurr);

					Eigen::MatrixXd PsiCond;
					Eigen::MatrixXd PsiSubsetIPropIProp = GetSubsetMatrix(Psi, IProp, IProp);	
					Eigen::MatrixXd PsiSubsetIPropNegIProp = GetSubsetMatrix(Psi, IProp, NegIProp);
					Eigen::MatrixXd PsiSubsetNegIPropIProp = GetSubsetMatrix(Psi, NegIProp, IProp);		
					Eigen::MatrixXd PsiSubsetNegIPropNegIProp = GetSubsetMatrix(Psi, NegIProp, NegIProp);	

					//-------------------------------------------------------
	                // Computing the target density at the proposal point
					//-------------------------------------------------------

                	//PsiCond = Psi(IProp,IProp) - Psi(IProp,~IProp) * inv(Psi(~IProp,~IProp)) * Psi(~IProp,IProp);
					PsiCond = PsiSubsetIPropIProp - PsiSubsetIPropNegIProp * PsiSubsetNegIPropNegIProp.inverse() * PsiSubsetNegIPropIProp;

                	muCond = GetSubsetVector(mu,IProp); // This needs to be changed if the prior mean is different from zero.
					paramPropSubset = GetSubsetVector(paramProp,IProp);
                	logPriorProp = LogPdfMultiN(paramPropSubset,muCond,PsiCond);
                	//logLikProp = feval(logPostFunc,paramProp,y,V,[],[],linkType,condDens,feat);
					logLikProp = HeteroGaussLogPost(paramProp, y, V, priorMeanDummy, priorCovDummy, linkType, condDens, feat);
                	logPostProp = logLikProp(0) + logPriorProp;

                	//feval(logPostFunc,paramProp,y,V,[],[],linkType,condDens,feat); // CHECK: what is this? Should be removed since it is not assigned.

					//-------------------------------------------------------
		            // Computing the target density in the current parameter point
					//-------------------------------------------------------

					Eigen::MatrixXd PsiSubsetICurrICurr = GetSubsetMatrix(Psi, ICurr, ICurr);	
					Eigen::MatrixXd PsiSubsetICurrNegICurr = GetSubsetMatrix(Psi, ICurr, NegICurr);
					Eigen::MatrixXd PsiSubsetNegICurrICurr = GetSubsetMatrix(Psi, NegICurr, ICurr);		
					Eigen::MatrixXd PsiSubsetNegICurrNegICurr = GetSubsetMatrix(Psi, NegICurr, NegICurr);	

	                //PsiCond = Psi(ICurr,ICurr) - Psi(ICurr,~ICurr)*inv(Psi(~ICurr,~ICurr))*Psi(~ICurr,ICurr);
	                PsiCond = PsiSubsetICurrICurr - PsiSubsetICurrNegICurr *  PsiSubsetNegICurrNegICurr.inverse() * PsiSubsetNegICurrICurr;
    	            muCond = GetSubsetVector(mu,ICurr); // This needs to be changed if the prior mean is different from zero.
					paramCurrSubset = GetSubsetVector(paramCurr,ICurr);
	                logPriorCurr = LogPdfMultiN(paramCurrSubset,muCond,PsiCond);
           		    //logLikCurr = feval(logPostFunc,paramCurr,y,V,[],[],linkType,condDens,feat);
					logLikCurr = HeteroGaussLogPost(paramCurr, y, V, priorMeanDummy, priorCovDummy, linkType, condDens, feat);
                	logPostCurr = logLikCurr(0) + logPriorCurr;
                
                	AccPrRaw = exp(logPostProp - logPostCurr + LogPropCurr - LogPropProp) * (priorIProp / priorICurr) * (gCurr / gProp);
					
					if (isnan(AccPrRaw))
					{
	                    AccPrRaw = 0.0;
					}
					
	                AccPr = std::min(1.0,AccPrRaw);
				}
            	else // Bad reversed proposal (Hessian is not ok). Reject.
				{
            	    AccPr = 0.0;
				}  // end if ~errFlagRev  
			}
			else // Bad proposal (Hessian or something else is not valid). Reject.
			{
				AccPr = 0.0;
			}  // end if ~errFlagProp       

			// Decide on acceptance and, possibly, update posterior arguments.
        	if (uniformDistribution(generator) < AccPr)
			{
	            paramCurr = paramProp;
            	ICurr = IProp;
        	}
	    } // End if (update>0 && all(IProp==ICurr)) || (update==0 && all(IProp==0))
    
	} // End Update loop
}











Eigen::VectorXd HeteroGauss::Lag(Eigen::VectorXd &vector, int lag)
{
	// Function that returns the k:th time lag of the input vector.

	int nT = vector.size();
	Eigen::VectorXd result = vector;

	for (int i = 0; i < nT; i++)
	{
		if (i < lag)
		{
			result(i) = 0.0;
		}
		else
		{
			result(i) = vector(i-lag);
		}
	}
	return result;
}

Eigen::MatrixXd HeteroGauss::Lag(Eigen::MatrixXd &matrix, int lag)
{
	// Function that returns the k:th time lag of the input vector.

	Eigen::MatrixXd result = matrix;

	for (int i = 0; i < matrix.cols(); i++)
	{
		Eigen::VectorXd vector = matrix.col(i);
		result.col(i) = Lag(vector,lag);
	}
	return result;
}

Eigen::MatrixXd HeteroGauss::Lags(Eigen::VectorXd &vector, int lag)
{
	// Function that returns a matrix with the first k time lags of the input vector.

	Eigen::MatrixXd result;
	result.resize(vector.size(),lag);

	for (int i = 0; i < lag;  i++)
	{
		result.col(i) = Lag(vector,i+1);
	}
	return result;
}

void HeteroGauss::Prewhitening(Eigen::VectorXd &ytilde, Eigen::MatrixXd &Xtilde, Eigen::VectorXd &y, Eigen::MatrixXd &X, Eigen::VectorXd &rho)
{
	int nLags = rho.size();
	int nT = y.size();

	if (nLags > 0)
	{
		ytilde = y;
		Xtilde = X;

		for (int i = 0; i < nLags; i++)
		{
			ytilde -= Lag(y,i+1) * rho(i);
			Xtilde -= Lag(X,i+1) * rho(i);
		}

		// Remove nLags first values
		Eigen::VectorXd tempV = ytilde.tail(nT-nLags);
		ytilde = tempV;
		Eigen::MatrixXd temp = Xtilde.block(nLags,0,nT-nLags,X.cols());		
		Xtilde = temp;
	}
	else
	{
		ytilde = y;
		Xtilde = X;
	}
}


void HeteroGauss::UpdateLinRegVarSel(Eigen::VectorXd &beta, 
                                     Eigen::VectorXd &I, 
                                     Eigen::VectorXd &y, 
                                     Eigen::MatrixXd &X, 
                                     Eigen::VectorXd &mu, 
                                     Eigen::MatrixXd &Sigma, 
                                     Eigen::VectorXd &PrIn, 
                                     Eigen::VectorXd &onTrialIndex, 
                                     bool AR, 
                                     bool forceStationarity,
									 int &nonStationaryDraws,
									 Eigen::VectorXd &previousBeta)
{
	// Update the indicator vector and the regression coefficient in the linear
	// regression model with unit noise variance
	// Model: y = X*beta_ + eps, eps~N(0,1)
	// Prior: beta_~N(mu,Sigma)
	//
	// INPUT:
	// I             (p-by-1 boolean)    Current variable selection indicators.
	// y             (n-by-1 vector)     Response vector  
	// X             (n-by-p matrix)     Covariate matrix
	// mu            (p-by-1 vector)     Prior mean vector
	// Sigma         (p-by-p diag)       Prior covariance matrix WHICH IS ASSUMED TO BE DIAGONAL (for speed).
	// PrIn          (p-vector)          Prior inclusion probabilities for each of the variables.          
	// onTrialIndex  (pIn-by-1)          Index of the variables that may be deleted from the model
	// AR            (Boolean scalar)    If AR = True, variables in I are only sampled
	//                                   if all previous elements in I are in the
	//                                   model. This is useful for lag order
	//                                   selection in AR processes.
	// forceStationarity (boolean)       Should stationarity be enforced in AR?
	// 
	// OUTPUT:
	// beta_         (p-by-1 vector)     Sampled regression coefficients.
	// I             (p-by-1 binary)     Sampled variable selection indicators.
	//
	// Author: Mattias Villani, Linkoping University. http://mattiasvillani.com
	//
	// C++ Author: Anders Eklund, Linkoping university

	Eigen::MatrixXd preCompMat, CompanionMat, XI, XIXI, SigmaI, IdentityPI, PostCov, invSigmaI;
	Eigen::VectorXd betaHatI, muI, allRows, ICurr, PostMean, draw, betaTildeI;

	allRows.resize(X.rows());
	allRows.setOnes(X.rows());

	double logDensCurr, logDensProp, PrInCurr, PrInProp;

	bool error;

	// Prelims
	int p = I.size();
	bool stationary = true;

	if (AR && forceStationarity) 
	{
		//I.setZero(p);
	    stationary = false;
		//preCompMat = [eye(p-1),zeros(p-1,1)]; // Used for checking stationarity below.
		preCompMat.resize(p-1,p);
		preCompMat.setZero(p-1,p);
		CompanionMat.resize(p,p);
		Eigen::MatrixXd temp;
		temp.setIdentity(p-1,p-1);
		preCompMat.block(0,0,p-1,p-1) = temp;    	
	}
	beta.setZero(p);

	int pI = (int)I.sum();	

	if (pI == 0) // No covariates in this model
	{
	    logDensCurr = -0.5 * y.transpose() * y;
	}
	else
	{
	    XI = GetSubsetMatrix(X, allRows, I);
		XIXI = XI.transpose() * XI;
    	betaHatI = XIXI.inverse() * XI.transpose() * y;     // Pseudo inverse
    	muI = GetSubsetVector(mu,I);
    	SigmaI = GetSubsetMatrix(Sigma,I,I);
    	invSigmaI = SigmaI.inverse();
    	betaTildeI = (XIXI + invSigmaI).lu().solve(XIXI * betaHatI + invSigmaI * muI);
		IdentityPI.setIdentity(pI,pI);
    	logDensCurr = -0.5 * log( ( IdentityPI + SigmaI * XIXI ).determinant() ) -0.5 * (y - XI * betaHatI).transpose() * (y - XI * betaHatI) - 0.5 * (betaTildeI - betaHatI).transpose() * XIXI * (betaTildeI - betaHatI) - 0.5 * (betaTildeI - muI).transpose() * invSigmaI * (betaTildeI - muI);

		//std::cout << "logDensCurr is " << std::endl << logDensCurr << std::endl << std::endl;
	}

	ICurr = I;

	//std::cout << "ICurr is " << std::endl << ICurr << std::endl << std::endl;

	//std::cout << "PrIn is " << std::endl << PrIn << std::endl << std::endl;

	// Gibbs sampling I by iterating through all elements. Maybe TODO: only update a random subset at each sweep.
	for (int j = 0; j < onTrialIndex.size(); j++)
	{
		int i = onTrialIndex(j);
    
		/*
		// Propose a change in I
    	if (AR)
		{
        	if ((i > 0) && (I(i-1) == 0.0)) // Previous lag is out, break out of the loop. Selection is over.
			{
            	break;
			}
	        else
			{
	            I(i) = 1.0;           // For AR lag selection, always add a lag.
	        }
		}
    	else
		{
        	I(i) = 1.0 - I(i);    // Propose a change of I in the i:th position (in->out or out->in)
		}
		*/

		I(i) = 1.0 - I(i);    // Propose a change of I in the i:th position (in->out or out->in)

	    pI = (int)I.sum();
    
	    if (pI == 0) // No covariates in the current model
		{
	         logDensProp = -0.5 * y.transpose() * y;
		}
	    else
		{
		    XI = GetSubsetMatrix(X, allRows, I);
	    	XIXI = XI.transpose() * XI;
			betaHatI = XIXI.inverse() * XI.transpose() * y;     // Pseudo inverse
	    	muI = GetSubsetVector(mu,I);
	    	SigmaI = GetSubsetMatrix(Sigma,I,I);
	    	invSigmaI = SigmaI.inverse();
	    	betaTildeI = (XIXI + invSigmaI).lu().solve(XIXI * betaHatI + invSigmaI * muI);
			IdentityPI.setIdentity(pI,pI);
	    	logDensProp = -0.5 * log( ( IdentityPI + SigmaI * XIXI ).determinant() ) -0.5 * (y - XI * betaHatI).transpose() * (y - XI * betaHatI) - 0.5 * (betaTildeI - betaHatI).transpose() * XIXI * (betaTildeI - betaHatI) - 0.5 * (betaTildeI - muI).transpose() * invSigmaI * (betaTildeI - muI);

			//std::cout << "logDensProp is " << std::endl << logDensProp << std::endl << std::endl;
		}

	    //PrInCurr = prod(PrIn.^ICurr.*(1-PrIn).^(1-ICurr));
	    //PrInProp = prod(PrIn.^I.*(1-PrIn).^(1-I));
		Eigen::VectorXd product1, product2;
		product1.resize(p);
		product2.resize(p);
		for (int i = 0; i < p; i++)
		{
			product1(i) = pow(PrIn(i),ICurr(i)) * pow((1.0 - PrIn(i)),1.0 - ICurr(i));
			product2(i) = pow(PrIn(i),I(i))     * pow((1.0 - PrIn(i)),1.0 - I(i));
		}
		PrInCurr = product1.prod();
		PrInProp = product2.prod();

		//std::cout << "PrInCurr is " << PrInCurr << std::endl << std::endl;
		//std::cout << "PrInProp is " << PrInProp << std::endl << std::endl;

	    //probs = [PrInCurr PrInProp].*exp([logDensCurr logDensProp]);  % TODO may cause over or underflow. Do the cancellation trick.
	    //probs = probs/sum(probs);  % Posterior inclusion probability.

		Eigen::VectorXd probs;
		probs.resize(2);
		probs(0) = PrInCurr * exp(logDensCurr);
		probs(1) = PrInProp * exp(logDensProp);
		probs = probs / probs.sum();
    
		//std::cout << "probs is " << std::endl << probs << std::endl << std::endl;

	    if (uniformDistribution(generator) < probs(0))    // Current wins, Keep ICurr
		{
	        I = ICurr; 
		}
	    else                // Proposed wins. I is correct. Update logDensCurr.
		{
	        logDensCurr = logDensProp;
	        ICurr = I;
	    }    
	}

	// Sample beta | I
	pI = (int)I.sum();
	if (pI > 0.0)
	{
	    XI = GetSubsetMatrix(X,allRows,I);
	    XIXI = XI.transpose() * XI;
		betaHatI = XIXI.inverse() * XI.transpose() * y;     // Pseudo inverse
	    muI = GetSubsetVector(mu,I);
		SigmaI = GetSubsetMatrix(Sigma,I,I);
	    invSigmaI = SigmaI.inverse();
	    PostCov = (XIXI + invSigmaI).inverse();
	    PostMean = PostCov * (XIXI * betaHatI + invSigmaI * muI);
		RndMultiN(draw, error, PostMean, PostCov); 

		// Try to fix "bad" covariance matrices
		double alpha = 0.000000001;
		Eigen::MatrixXd identity;
		identity.setIdentity(I.size(),I.size());
		while(error)
		{
			Eigen::MatrixXd PostCov_ = PostCov + alpha * identity;
			RndMultiN(draw, error, PostMean, PostCov_); 		
			alpha = alpha * 10;	
		}

		/*		
	    RndMultiN(draw, error, PostMean, PostCov); 
		if (error)
		{
			std::cout << "I is " << I << std::endl << std::endl;
			std::cout << "X is " << X << std::endl << std::endl;
			std::cout << "XI is " << XI << std::endl << std::endl;
			std::cout << "XIXI is " << XIXI << std::endl << std::endl;
			std::cout << "betahatI is " << betaHatI << std::endl << std::endl;
			std::cout << "muI is " << muI << std::endl << std::endl;
			std::cout << "SigmaI is " << SigmaI << std::endl << std::endl;
			std::cout << "invSigmaI is " << invSigmaI << std::endl << std::endl;
			std::cout << "PostCov is " << PostCov << std::endl << std::endl;
			std::cout << "PostMean is " << PostMean << std::endl << std::endl;
		}
		*/

		beta = SaveFullVector(draw, I);

		int failedAttempts = 0;

    	while (!stationary && (failedAttempts < 1000)) // TODO: faster check of stationarity?
		{
			//CompanionMat =  [beta_' ; preCompMat ];
			CompanionMat.row(0) = beta.transpose();
			CompanionMat.block(1,0,p-1,p) = preCompMat;
	      
			Eigen::VectorXcd eigs = CompanionMat.eigenvalues();  
			bool allSmallerThanOne = true;

			for (int i = 0; i < eigs.size(); i++)
			{
				double magnitude = sqrt( eigs(i).real() * eigs(i).real() + eigs(i).imag() * eigs(i).imag());
				if ( magnitude >= 1.0 )
				{
					allSmallerThanOne = false;
					break;		
				}
			}

			if (allSmallerThanOne)
			{
				break;
			}
			else
			{
				failedAttempts++;
			}

			RndMultiN(draw, error, PostMean, PostCov); 

			// Try to fix "bad" covariance matrices
			double alpha = 0.000000001;
			Eigen::MatrixXd identity;
			identity.setIdentity(I.size(),I.size());
			while(error)
			{
				Eigen::MatrixXd PostCov_ = PostCov + alpha * identity;
				RndMultiN(draw, error, PostMean, PostCov_); 		
				alpha = alpha * 10;	
			}

			beta = SaveFullVector(draw, I);
		}

		if (failedAttempts >= 1000)
		{
			//beta.setZero(p);
			beta = previousBeta;
			nonStationaryDraws++;
		}
    }
}



void HeteroGauss::UpdateInclusionProbability(Eigen::VectorXd &PrIn, Eigen::VectorXd &Indicators, double a, double b)
{	
	int p = PrIn.size();

	double atot = a + Indicators.sum();
	double btot = b + (double)p - Indicators.sum();

	double value = RndBeta(atot,btot);

	for (int i = 0; i < PrIn.size(); i++)
	{
		PrIn(i) = value;
	}
}



void HeteroGauss::GibbsHIGLM(Eigen::MatrixXd &betaDraws, 
                             Eigen::MatrixXd &IbetaDraws, 
                             Eigen::MatrixXd &gammaDraws, 
                             Eigen::MatrixXd &IgammaDraws, 
                             Eigen::MatrixXd &rhoDraws, 
                             Eigen::MatrixXd &IrhoDraws, 
                             Eigen::VectorXd &accPrGammaDraws, 



							 // For debugging
							 Eigen::VectorXd &beta,
							 Eigen::VectorXd &u,
							 Eigen::MatrixXd &U,
							 Eigen::VectorXd &ytilde,
							 Eigen::MatrixXd &Xtilde,
							 Eigen::VectorXd &rho,

                             Eigen::VectorXd &y, 
                             Eigen::MatrixXd &X, 
                             Eigen::MatrixXd Z, 
                             int ARorder, 
                             bool forceStationarity, 
                             Eigen::VectorXd &muBeta, 
                             Eigen::VectorXd &muGamma, 
                             double tauBeta, 
                             double tauGamma, 
                             double tauRho, 
                             double iota, 
                             double r, 
                             Eigen::VectorXd &PrInBeta, 
                             Eigen::VectorXd &PrInGamma, 
                             Eigen::VectorXd &PrInRho, 
                             Eigen::VectorXd &onTrialBeta, 
                             Eigen::VectorXd &onTrialGamma, 
                             Eigen::VectorXd &onTrialRho, 
                             int nIter, 
                             double prcBurnin, 
                             int nStepsGamma, 
                             int hessMethodGamma, 
							 int linkType,
                             int propDfGamma, 
                             double IUpdatePrGamma,

							 bool updateInclusion,

							 int &nonStationaryDraws)
{

	//
	// Purpose:      Gibbs sampling and variable selection for dynamic GLM with heteroscedastic innovations (HIGLM):
	// 
	// Model:    yt = xt*beta + ut
	//           ut = rho_1*u(t-1) + ... rho_k*u(t-k) + exp(z_t'*gamma)*e(t) 
	//
	//           The algorithm simulates from the joint posterior 
	//           p(beta,gamma,rho, Ibeta, Igamma, Irho | y, X),
	//           where Ibeta, Igamma, Irho are vectors of variable selection indicators.
	// 
	// Priors:   beta  ~ N(0,tauBeta^2*I)  and we also use the notation B0 = tauBeta^2*I
	//           gamma ~ N(0,tauGamma^2*I) and we also use the notation C0 = tauGamma^2*I
	//           rho   ~ N{[r 0 ... 0]' , tauRho^2*Diag(1,1/2^iota,...,1/k^iota)} 
	//                   and we also use the notation R0 = tauRho^2*Diag(1,1/2^iota,...,1/k^iota)
	// 
	// INPUT:
	// y             (T-by-1)    Response variable 
	// X             (T-by-p)    Covariates in the mean
	// Z             (T-by-q)    Covariates in the (log) variance
	// modelOpt      (struct)    Options for the model, e.g. modelOpt.k is the number of lags in the noise model
	// priorOpt      (struct)    Options for the priors, e.g priorOpt.beta0 is the prior mean of beta.
	// algoOpt       (struct)    Options for the algorithm, e.g. algoOpt.nStepsGamma is the number of Newtons steps 
	//                           for the variance parameters.
	// 
	// OUTPUT:
	// paramDraws    (struct)    MCMC draws: paramDraws.gamma, paramDraws.beta, paramDraws.rho, paramDraws.accPrGamma.
	//
	// Original author:      Mattias Villani, mattiasvillani.com
	//
	// C++ author:   		Anders Eklund
	//
	// First ver:    		2015-05-04
	// This ver:     		2016-01-26
	//

	//Eigen::VectorXd ytilde;
	//Eigen::MatrixXd Xtilde;

	double accPrGamma;

	int t = X.rows(); 	  // Number of timepoints
	int p = X.cols();     // Number of regressors for beta
	int q  = Z.cols();    // Number of regressors for gamma

	// Prior for beta
	Eigen::MatrixXd SigmaBeta;
	SigmaBeta.setIdentity(p,p);
	SigmaBeta = tauBeta * tauBeta * SigmaBeta;
	Eigen::MatrixXd invSigmaBeta = SigmaBeta.inverse();

	// Prior for gamma
	Eigen::MatrixXd SigmaGamma;
	SigmaGamma.setIdentity(q,q);
	SigmaGamma = tauGamma * tauGamma * SigmaGamma;
	Eigen::MatrixXd invSigmaGamma = SigmaGamma.inverse();

	// Prior for AR parameters, rho, and pre-whitening
    Eigen::VectorXd muRho; 
	if (ARorder > 0)
	{
		muRho.resize(ARorder);
		muRho(0) = r;
		for (int i = 1; i < ARorder; i++)
		{
			muRho(i) = 0.0;
		}		
	}
	else
	{
		muRho.resize(0);
		ytilde = y;
		Xtilde = X;
	}
	Eigen::VectorXd SigmaRhoVector(ARorder);
	for (int i = 0; i < ARorder; i++)
	{
		SigmaRhoVector(i) = 1.0 / (pow((double)(i + 1),iota));		
	}
	Eigen::MatrixXd SigmaRho = tauRho * tauRho * SigmaRhoVector.asDiagonal();
	Eigen::MatrixXd invSigmaRho = SigmaRho.inverse();

	//std::cout << "SigmaRho is " << std::endl << SigmaRho << std::endl << std::endl;

	betaDraws.resize(nIter,p);
	IbetaDraws.resize(nIter,p);

	gammaDraws.resize(nIter,q);
	IgammaDraws.resize(nIter,q);

	rhoDraws.resize(nIter,ARorder);
	IrhoDraws.resize(nIter,ARorder);

	accPrGammaDraws.resize(nIter);

	// Initial values

	Eigen::VectorXd Ibeta(p);
	Ibeta.setOnes(p);
	// All on trial are out to start with.	
	for (int i = 0; i < onTrialBeta.size(); i++)
	{
		Ibeta(onTrialBeta(i)) = 0.0;
	}	

	//std::cout << "onTrialBeta is " << std::endl << onTrialBeta << std::endl << std::endl;
	//std::cout << "Ibeta is " << std::endl << Ibeta << std::endl << std::endl;

	//Eigen::VectorXd beta = (X.transpose() * X).inverse() * X.transpose() * y; // Pseudo inverse
	beta = (X.transpose() * X).inverse() * X.transpose() * y; // Pseudo inverse

	//std::cout << "beta is " << std::endl << beta << std::endl << std::endl;

	Eigen::VectorXd Igamma(q);
	Igamma.setOnes(q);
	// All on trial are out to start with.
	for (int i = 0; i < onTrialGamma.size(); i++)
	{
		Igamma(onTrialGamma(i)) = 0.0;
	}	

	//gamma_(Igamma_) = Z(:,Igamma_) \ log((y - X*beta_).^2);	
	Eigen::VectorXd Tones(t);
	Tones.setOnes(t);
	Eigen::MatrixXd ZSubsetIgamma = GetSubsetMatrix(Z,Tones,Igamma);
	Eigen::VectorXd tempVector = ((y - X * beta).array().square()).array().log();
	Eigen::VectorXd tempGamma = (ZSubsetIgamma.transpose() * ZSubsetIgamma).inverse() * ZSubsetIgamma.transpose() * tempVector; // pseudo inverse
	Eigen::VectorXd gamma = SaveFullVector(tempGamma,Igamma);
	
	//std::cout << "original gamma is " << std::endl << gamma << std::endl << std::endl;

	Eigen::VectorXd previousRho(ARorder);
	previousRho.setZero(ARorder);

	Eigen::VectorXd previousDummy(100);

	Eigen::VectorXd Irho(ARorder);
	Irho(0) = 1;
	for (int i = 1; i < ARorder; i++)
	{
		Irho(i) = 0.0;
	}
	//Eigen::VectorXd rho = muRho; // Starting Gibbs sampling with rho at its prior mean
	rho = muRho; // Starting Gibbs sampling with rho at its prior mean

	//Eigen::VectorXd u = y - X * beta; // Residuals
	u = y - X * beta; // Residuals

	//Eigen::MatrixXd U = Lags(u,ARorder);
	U = Lags(u,ARorder); 		 		
	u = u.tail(t - ARorder);

	Eigen::MatrixXd temp;
		
	temp = U.block(ARorder,0,t - ARorder,ARorder);
	U = temp;
	
	temp = Z.block(ARorder,0,t - ARorder,q);
	Z = temp;

	Eigen::MatrixXd feat;
	feat.resize(t - ARorder,2);
	if (ARorder == 0)
	{
		Eigen::VectorXd zeros;
		zeros.resize(t - ARorder);
		zeros.setZero(t - ARorder);
		feat.col(0) = zeros;
	}
	else
	{
		feat.col(0) = U * rho;
	}

	Eigen::MatrixXd Utilde;
	Utilde.resize(t - ARorder,ARorder);

	// Gibbs sampling
	int nBurnin = round((float)nIter*(prcBurnin/100.0));
	for (int i = 0; i < (nBurnin + nIter); i++)
	{
		//------------------------------------------------------------
		// Block 1. Update gamma and Igamma (the parameters in the variance)
		//------------------------------------------------------------

		int logPostFunc = 0; int gradHessFunc = 0; int condDens = 2;
		bool errFlagProp, errFlagRev;

	    NewtonProp(gamma, Igamma, accPrGamma, errFlagProp, errFlagRev, PrInGamma, muGamma, invSigmaGamma, onTrialGamma, logPostFunc, gradHessFunc, u, Z, linkType, condDens, feat, nStepsGamma, propDfGamma, IUpdatePrGamma, hessMethodGamma);

		//std::cout << "gamma is " << std::endl << gamma << std::endl << std::endl;
		//std::cout << "I gamma is " << std::endl << Igamma.transpose() << std::endl << std::endl;
		
	    Eigen::VectorXd invSigmas = (-Z * gamma / 2.0).array().exp(); // (n-k)-by-1 vector with noise standard deviations

		//------------------------------------------------------------
		// Update pi-gamma (inclusion probability)
		//------------------------------------------------------------

		if (updateInclusion)
		{
			UpdateInclusionProbability(PrInGamma,Igamma,3.0,3.0);
		}

		//------------------------------------------------------------    
		// Block 2. Update beta and Ibeta
		//------------------------------------------------------------

	    Prewhitening(ytilde, Xtilde, y, X, rho);
	    ytilde = invSigmas.array() * ytilde.array();

	    //Xtilde = repmat(invSigmas,1,p).*Xtilde;
		for (int j = 0; j < p; j++)
		{
			Eigen::VectorXd temp = Xtilde.col(j);
			Xtilde.col(j) = temp.array() * invSigmas.array();
		}

	    UpdateLinRegVarSel(beta, Ibeta, ytilde, Xtilde, muBeta, SigmaBeta, PrInBeta, onTrialBeta, false, false, nonStationaryDraws, previousDummy);
	
		
		//------------------------------------------------------------
		// Update pi-beta (inclusion probability)
		//------------------------------------------------------------

		if (updateInclusion)
		{
			UpdateInclusionProbability(PrInBeta,Ibeta,3.0,3.0);
		}

		//------------------------------------------------------------
		// Block 3. Update rho
		//------------------------------------------------------------

	    u = y - X*beta; // Residuals

		// Do not care about number of non-stationary draws during burnin
		if (i == nBurnin)
		{
			nonStationaryDraws = 0;
		}
		
		if (ARorder > 0)
		{
    	    U = Lags(u,ARorder); 
			Eigen::VectorXd tempv = u.tail(t - ARorder);
			u = tempv;

			temp = U.block(ARorder,0,t - ARorder,ARorder);					
			U = temp;

    	    Eigen::VectorXd utilde = invSigmas.array() * u.array();

			//Utilde = repmat(invSigmas,1,k).*U;			
    	    for (int j = 0; j < ARorder; j++)
			{
				Eigen::VectorXd temp = U.col(j);
				Utilde.col(j) = temp.array() * invSigmas.array();
			}		

    	    UpdateLinRegVarSel(rho, Irho, utilde, Utilde, muRho, SigmaRho, PrInRho, onTrialRho, true, forceStationarity, nonStationaryDraws, previousRho);
			previousRho = rho;			

			//std::cout << "rho is " << std::endl << rho << std::endl << std::endl;
			//std::cout << "Irho is " << std::endl << Irho << std::endl << std::endl;

    	    feat.col(0) = U * rho; // This is the mean in the heteroscedastic regression. 
                           // Note this should NOT be Utilde*rho. This is used for updating gamma.
			
    	}
				    
		// Storing results
	    if (i >= nBurnin)
		{	
			betaDraws.block(i-nBurnin,0,1,p) = beta.transpose();
			IbetaDraws.block(i-nBurnin,0,1,p) = Ibeta.transpose();

			gammaDraws.block(i-nBurnin,0,1,q) = gamma.transpose();
			IgammaDraws.block(i-nBurnin,0,1,q) = Igamma.transpose();

			rhoDraws.block(i-nBurnin,0,1,ARorder) = rho.transpose();
			IrhoDraws.block(i-nBurnin,0,1,ARorder) = Irho.transpose();

    	    accPrGammaDraws(i-nBurnin) = accPrGamma;			
		}    
	}	
}


