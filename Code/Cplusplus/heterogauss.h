#ifndef HETEROGAUSS_H
#define HETEROGAUSS_H

#include <Dense>
#include <Eigen>
#include <random>

class HeteroGauss
{
	public:

		// Constructors & destructor
		HeteroGauss();
		~HeteroGauss();

		// Help functions

		Eigen::VectorXd GetSubsetVector(Eigen::VectorXd &vector, Eigen::VectorXd &indicators);
		Eigen::VectorXd SaveFullVector(Eigen::VectorXd &vector, Eigen::VectorXd &indicators);

		Eigen::MatrixXd GetSubsetMatrix(Eigen::MatrixXd &matrix, Eigen::VectorXd &rowIndicators, Eigen::VectorXd &colIndicators);

		Eigen::VectorXd FlipIndicators(Eigen::VectorXd &indicators);

		double RndBeta(double a, double b);
		double RndChi2(int df);
		void RndMultiT(Eigen::VectorXd &draw, bool &error, Eigen::VectorXd &mu, Eigen::MatrixXd &Sigma, int);
		void RndMultiN(Eigen::VectorXd &draw, bool &error, Eigen::VectorXd &mu, Eigen::MatrixXd &Sigma);

		double LogPdfMultiN(Eigen::VectorXd &x, Eigen::VectorXd &mu, Eigen::MatrixXd &Sigma);
		void LogPdfMultiT(double &result, bool &error, Eigen::VectorXd &x, Eigen::VectorXd &mu, Eigen::MatrixXd &Sigma, int df);

		Eigen::VectorXd invLinkEval(Eigen::VectorXd &param, Eigen::MatrixXd &X, int linkType);
		void LinkDeriv(Eigen::VectorXd &firstDeriv, Eigen::VectorXd &secondDeriv, Eigen::VectorXd &feat, int linkType, bool second);
				
		void HeteroGaussGradHess(Eigen::VectorXd &gradObs, Eigen::VectorXd &hessObs, Eigen::VectorXd &param, Eigen::VectorXd &y, Eigen::MatrixXd &V, int linkType, int condDens, Eigen::MatrixXd &feat, int hessMethod);

		Eigen::VectorXd HeteroGaussLogPost(Eigen::VectorXd &param, Eigen::VectorXd &y, Eigen::MatrixXd &V, Eigen::VectorXd &priorMean, Eigen::MatrixXd &priorCov, int linkType, int condDens, Eigen::MatrixXd &feat);

		void NewtonGradHess(Eigen::VectorXd &Grad, Eigen::MatrixXd &Hesspp, Eigen::MatrixXd &Hesspc, Eigen::VectorXd &Beta_, Eigen::VectorXd &Ic, Eigen::VectorXd &Ip, Eigen::VectorXd &mu, Eigen::MatrixXd &invPsi, int hessMethod, int gradHessFunc, Eigen::VectorXd &y, Eigen::MatrixXd &V, int linkType, int condDens, Eigen::MatrixXd &feat);
		
		void NewtonIter(Eigen::VectorXd &betaEnd, Eigen::MatrixXd &hess, Eigen::VectorXd &betaC, Eigen::VectorXd &Ic, Eigen::VectorXd &Ip, Eigen::VectorXd &mu, Eigen::MatrixXd &invPsi, int nSteps, int hessMethod, int gradHessFunc, Eigen::VectorXd &y, Eigen::MatrixXd &V, int linkType, int condDens, Eigen::MatrixXd &feat);
				
		void NewtonProp(Eigen::VectorXd &paramCurr, Eigen::VectorXd &ICurr, double & AccPr, bool & errFlagProp, bool & errFlagRev, Eigen::VectorXd &PrIn, Eigen::VectorXd &mu, Eigen::MatrixXd &invPsi, Eigen::VectorXd &onTrial, int logPostFunc, int gradHessFunc, Eigen::VectorXd &y, Eigen::MatrixXd &V, int linkType, int condDens, Eigen::MatrixXd &feat, int nSteps, int propDf, double IUpdatePr, int hessMethod);	
		
		void GibbsHIGLM(Eigen::MatrixXd &betaDraws, Eigen::MatrixXd &IbetaDraws, Eigen::MatrixXd &gammaDraws, Eigen::MatrixXd &IgammaDraws, Eigen::MatrixXd &rhoDraws, Eigen::MatrixXd &IrhoDraws, Eigen::VectorXd &accPrGammaDraws,  Eigen::VectorXd &beta, Eigen::VectorXd &u, Eigen::MatrixXd &U, Eigen::VectorXd &ytilde, Eigen::MatrixXd &Xtilde,   Eigen::VectorXd &rho,   Eigen::VectorXd &y, Eigen::MatrixXd &X, Eigen::MatrixXd Z, int ARorder, bool forceStationarity, Eigen::VectorXd &muBeta, Eigen::VectorXd &muGamma, double tauBeta, double tauGamma, double tauRho, double iota, double r, Eigen::VectorXd &PrInBeta, Eigen::VectorXd &PrInGamma, Eigen::VectorXd &PrInRho, Eigen::VectorXd &onTrialBeta, Eigen::VectorXd &onTrialGamma, Eigen::VectorXd &onTrialRho, int nIter, double prcBurnin, int nStepsGamma, int hessMethodGamma, int linkType, int propDfGamma, double IUpdatePrGamma, bool updateInclusion, int &nonStationaryDraws);

		void UpdateLinRegVarSel(Eigen::VectorXd &beta, Eigen::VectorXd &I, Eigen::VectorXd & y, Eigen::MatrixXd & X, Eigen::VectorXd & mu, Eigen::MatrixXd & Sigma, Eigen::VectorXd & PrIn, Eigen::VectorXd & onTrialIndex, bool AR, bool forceStationarity, int &nonStationaryDraws, Eigen::VectorXd & previousBeta);

		Eigen::VectorXd Lag(Eigen::VectorXd &vector, int lag);
		Eigen::MatrixXd Lag(Eigen::MatrixXd &matrix, int lag);
		Eigen::MatrixXd Lags(Eigen::VectorXd &vector, int lag);
		void Prewhitening(Eigen::VectorXd &ytilde, Eigen::MatrixXd &Xtilde, Eigen::VectorXd &y, Eigen::MatrixXd &X, Eigen::VectorXd &rho);

		void UpdateInclusionProbability(Eigen::VectorXd &PrIn, Eigen::VectorXd &Indicators, double a, double b);

    private:

	    std::default_random_engine generator;
	    std::normal_distribution<double> normalDistribution; // default mean 0, std 1
	    std::uniform_real_distribution<double> uniformDistribution; // default min 0, max 1 
};

#endif
