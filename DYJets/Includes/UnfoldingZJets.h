#ifndef _UNFOLDINGZJETS_h_
#define _UNFOLDINGZJETS_h_

#include "RooUnfold.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TH1D.h"
#include "TProfile.h"
#include "SectionedConfig.h"
#include "TGraphAsymmErrors.h"
#include <vector>

//void UnfoldingZJets(TString lepSel, TString algo, TString histoDir, TString unfoldDir, int jetPtMin, int jetEtaMax, TString gen1, TString gen2, TString variable = "", bool noramlized = false);
void UnfoldingZJets(const SectionedConfig& unfCfg, TString lepSel, TString algo,
		    TString histoDir, TString unfoldDir, int jetPtMin, int jetEtaMax, 
		    TString variable = "", bool noramlized = false, int whichSyst = -1);

int UnfoldData(const SectionedConfig& unfCfg, const TString lepSel, const char* variable, RooUnfoldResponse *resp,
	       TH1D *hRecData, TH1D* &hUnfData, TH2D* &hUnfDataStatCov, TH2D* &hUnfMCStatCov,
	       TString name, double integratedLumi, const TString& unfoldDir,
	       TH1* hRecDYJets, TH1* hGenDYJets, bool logy = false,
	       TH1D *hRecDataMinusFakesOdd = 0, TH1D *hRecDataMinusFakesEven = 0, int fixNIterTo = -1,
	       const char* outputFileName = 0);
TH2D* M2H(TMatrixD m);
TH2D* makeCovFromUpAndDown(const TH1D* hUnfDataCentral, const TH1D* hUnfDataUp, const TH1D* hUnfDataDown, TString name,
			   bool noNull = true);
TH1D* foldUnfData(const TH1 *hUnfData, const TMatrixD* cov, const RooUnfoldResponse *resp);
void test();
void createSystPlots(TString outputFileName, TString sysPlotDir, TString lepSel, TString variable, TH1D* hUnf[], bool logy = false);
double MyChi2Test(const TH1 *h1, const TH1 *h2, int nFirstBinsToSkip = 0, int nLastBinsToSkip = 0,
		  Double_t* res = 0, const TH1* herr = 0, bool poisErr = false);
double TH1Chi2Test(const TH1 *h1, const TH1 *h2, int nFirstBinsToSkip = 0, int nLastBinsToSkip = 0, Double_t* res = 0);
void RemoveFakes(TH1* hRecData, TH1* hFakes, TH1* hPurity);
TH1* resample(const TH1* h);
/// uncMode: 0- only measurement stat. unc., 1- both measurement and MC stat. unc, 2- only MC stat. unc.
TH1* unfold(RooUnfold::Algorithm algo, const RooUnfoldResponse* resp,
	    const TH1* hRecDataMinusFakes, int niters, bool smoothPrior, std::vector<TH1*>* hUnfs = 0,
	    int uncMode = 0);
TH1* unfoldWithErr(RooUnfold::Algorithm algo, const RooUnfoldResponse* resp,
		   const TH1* hRecDataMinusFakes, int niters, bool smoothPrior, std::vector<TH1*>* hUnfs = 0,
		   std::vector<TMatrixD>* covs = 0,
		   int uncMode = 0);
std::vector<double> chi2FromToy(RooUnfold::Algorithm algo, bool smoothPrior, const RooUnfoldResponse* resp,
				const TH1* hmeas, int nFirstBinsToSkip, int nLastBinsToSkip,
				int niters, int ntoys,
				std::vector<double>* chi2McErr = 0, std::vector<double>* rms = 0,
				std::vector<TH1*>* hChi2 = 0,  std::vector<TH1*>* hChi2mcErr = 0,
				std::vector<double>* meanII = 0, std::vector<double>* rmsII = 0, 
				std::vector<TH1*>* hChi2II = 0,
				std::vector<TProfile*>* hRes = 0, std::vector<TProfile*>* hResMcErr = 0,
				std::vector<TProfile*>* hResII = 0, bool toyIIrecoResampling = false);
void makeChi2Plot(TH1* hChi2, int thisAlgoNIters, int chosenNIters,
		  double thisAlgoThr, double chosenAlgoThr,
		  const char* outDir, const char* fileBaseName, const char* canvasName, 
		  int chi2Ylog = 0, std::vector<TH1*> hAltChi2 = std::vector<TH1*>(),
		  std::vector<std::string> labels = std::vector<std::string>(), bool limitYmax = false);
void testWithForgesResp(TH1* hUnf, int nIters, const char* outDir, const char* outfileBaseName, int nBinsToSkip = 0);
//int readNiters(const char* lepSel, const char* variable);

TGraphAsymmErrors* statFromToy(RooUnfold::Algorithm algo, bool smoothPrior, const RooUnfoldResponse* resp,
			       const TH1* hmeas, int niters, int ntoys, TH1** hUnfStdUnc = 0);

#endif //_UNFOLDINGZJETS_h_ not defined
