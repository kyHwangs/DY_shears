//-*- mode: c++; c-basic-offset: 4 -*-
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <memory>
#include <math.h>
#include <algorithm>
#include <TStyle.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TLine.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldBinByBin.h>
#include <RooUnfoldSvd.h>
#include <TSVDUnfold.h>
#include <TParameter.h>
#include "fileNamesZJets.h"
#include "getFilesAndHistogramsZJets.h"
#include "variablesOfInterestZJets.h"
#include "UnfoldingZJets.h"
#include "PlotSettings.h"
#include "fixYscale.h"
#include "ConfigVJets.h"
#include "functions.h"
#include "TRandom3.h"
#include <sys/time.h>
#include "TDecompSVD.h"
#include "ResultTables.h"
#include "Uncertainties.h"

extern ConfigVJets cfg;

void createInclusivePlots(bool doNormalized, TString outputFileName, TString lepSel, TH1 *hUnfData, TH2 *hCov[],
			  const std::vector<string>& predictions, double integratedLumi);

//void createInclusivePlots(bool doNormalized, TString outputFileName, TString lepSel, TH1D *hUnfData, TH2D *hCov[], TH1D *hMadGenCrossSection, TH1D *hSheGenCrossSection);
 //void createInclusivePlots(bool doNormalized, TString outputFileName, TString lepSel, TH1D *hUnfData, TH2D *hCov[], TH1D *hMadGenCrossSection);


 const static bool isdatabug = false;
 double pValueToNormChi2(double alpha, int n);
 static int verbosity = 1;

 using namespace std;

 //void UnfoldingZJets(TString lepSel, TString algo, TString histoDir, TString unfoldDir, 
 //        int jetPtMin, int jetEtaMax, TString gen1, TString gen2, TString variable, bool doNormalized)
void UnfoldingZJets(const SectionedConfig& unfCfg, TString lepSel, TString algo, TString histoDir, TString unfoldDir, 
		    int jetPtMin, int jetEtaMax, TString variable, bool doNormalized,
		    int whichSyst)
{
     gStyle->SetOptStat(0);
     //--- create output directory if does not exist ---
     system("mkdir -p " + unfoldDir);

     int start = 0;
     int end = NVAROFINTERESTZJETS;

     TString sysPlotDir = unfoldDir.Strip(TString::kTrailing, '/') + "SysPlots";

     bool cachedNIters = cfg.getB("cachedNIters", false); 
     bool pseudoData = cfg.getB("pseudoData", false);
     std::string pseudoDataSample = cfg.getS("pseudoDataSample", "");
     std::string pseudoDataSamplePrettyName = cfg.getS(TString::Format("prettyName_%s", pseudoDataSample.c_str()), pseudoDataSample);


     if (variable != "") {
	 
	 //ZNGoodJets_Zexc is produced together with ZNGoodJets_Zexc:
	 if(variable=="ZNGoodJets_Zinc") variable = "ZNGoodJets_Zexc";

	 start = findVariable(variable);
	 if (start >= 0) {
	     end = start + 1;
	     std::cout << "Processing only variable: " << variable << std::endl;
	 }
	 else {
	     cerr << "\nError: variable " << variable << " is not interesting." << endl;
	     cerr << "See below the list of interesting variables:" << endl;
	     for (unsigned int i = 0; i < NVAROFINTERESTZJETS; ++i) {
		 cerr << "\t" << i << ": " << VAROFINTERESTZJETS[i].name << "\n" << endl;
	     }
	     return;
	 }
     }

     // Here we declare the different arrays of TFiles. 
     // fData is for the three data files: 
     // 0 - central, 1 - JES up, 2 - JES down
     TFile *fData[3] = {NULL}; 
     // fDYJets is for the five DYJets files:
     // 0 - central, 1 - PU up, 2 - PU down, 3 - JER up, 4 - JER down, 5 - LES up, 6 - LES down, 7 - LER up, 8 - LER down  9 - Alt. Unf.
     TFile *fDYJets[10] = {NULL};
     // fBg is for the NBGDYJETS x 5 systematics files:
     // 0 - central, 1 - PU up, 2 - PU down, 3 - XSEc up, 4 - XSEC down, 5 - LES up, 6 - LES down 
     TFile *fBg[NBGDYJETS][7] = {{NULL}};

     //--- Now run on the different variables ---
     for (int i = start; i < end; ++i) {
	 timeval t0, t1;
	 gettimeofday(&t0,0);

	 //--- Open all files ---------------------------------------------------------------------- 
	 //Note: we close and reopen files for each variable to trigger deletion of histogram of
	 //the previous variables and reduce the memory usage.
	 getAllFiles(histoDir, lepSel, "13TeV", jetPtMin, jetEtaMax, fData, fDYJets, fBg, NBGDYJETS);
	 //----------------------------------------------------------------------------------------- 

	 //reads integrated luminosity
	 double integratedLumi = -1;
	 TH1* Lumi = 0;
	 if(fData[0]) fData[0]->GetObject("Lumi", Lumi);
	 if(Lumi) integratedLumi = Lumi->GetBinContent(1);
	 else {
	     cerr << "Error: Lumi histogram was not found.\n";
	     return;
	 }

	 
	 variable = VAROFINTERESTZJETS[i].name;
	 TString outputFileName = getUnfoldedFileName(unfoldDir, lepSel, variable, algo,
						      jetPtMin, jetEtaMax, "_MGPYTHIA6_", doNormalized);

	 //TFile *outputRootFile = new TFile(outputFileName + ".root", "RECREATE");
	 std::unique_ptr<TFile> outputRootFile(new TFile(outputFileName + ".root", "RECREATE"));


	 //Lumi is stored in the unfolded histo file, so it can be used by the
	 //channel combination code.
	 Lumi->Write();

	 
	 //	TString section = TString::Format("%s_%s", lepSel.Data(), variable.Data());
	 bool withUnfUnc      = cfg.getUnf(lepSel, variable, "withUnfUnc", true);

	 int nFirstBinsToSkip = cfg.getUnf(lepSel, variable, "nFirstBinsToSkip", 0);
	 int nLastBinsToSkip  = cfg.getUnf(lepSel, variable, "nLastBinsToSkip", 0);	
	 int normalizeToBin = cfg.getUnf(lepSel, variable, "normalizeToBin", -1);

	 //FIXME use the names instead of the indices. See Include/Uncertainties.h
	 std::vector<int> disabledUncs   = cfg.getVI("disabledUncs", std::vector<int>());
	 
	 //--- rec Data histograms ---
	 TH1D *hRecData[3] = {NULL};
	 TH1D* hRecDataOdd = 0;
	 TH1D* hRecDataEven = 0;

	 //--- rec DYJets histograms ---
	 TH1D *hRecDYJets[18] = {NULL};
	 //--- fake DYJets histograms ---
	 TH1D *hFakDYJets[18] = {NULL};
	 //--- purity "histograms" ---
	 TH1D *hPurity[18] = {NULL};
	 //--- gen DYJets histograms ---
	 TH1D *hGenDYJets[18] = {NULL};
	 //--- res DYJets histograms ---
	 TH2D *hResDYJets[18] = {NULL};
	 //--- rec Bg histograms ---
	 TH1D *hRecBg[NBGDYJETS][11] = {{NULL}};
	 //--- rec Sum Bg histograms ---
	 TH1D *hRecSumBg[11] = {NULL};
	 //--- response DYJets objects ---
	 RooUnfoldResponse *respDYJets[18] = {NULL};

	 //Get split histogram of data used for cross-validation
	 TString variable_odd = variable + "_Odd";
	 TString variable_even = variable + "_Even";

	 if(fData[0]){
	     fData[0]->cd();
	     hRecDataOdd = (TH1D*) fData[0]->Get(variable_odd);
	     hRecDataEven = (TH1D*) fData[0]->Get(variable_even);
	     if(!hRecDataOdd){ std::cout << "Missing Odd histogram for " << variable << "\n"; }
	     if(!hRecDataEven) std::cout << "Missing Even histogram for " << variable << "\n";
	 }

	 int nbg;
	 //In case of pseudo data, which do not include the background contributions,
	 //we should not load the background histograms otherwise they
	 //will be subtracted from the pseudodata in the fake calculation
	 //within the call to getAllHistos
	 if(!pseudoData) nbg = NBGDYJETS;
	 else nbg= 0;

	 //--- Get all histograms ---
	 getAllHistos(variable, hRecData, fData, 
		      hRecDYJets, hGenDYJets, hResDYJets, fDYJets,
		      hRecBg, hRecSumBg, fBg, nbg, respDYJets, hFakDYJets, hPurity);


	 std::vector<std::string> predictions = cfg.getVS("predictions");
	 std::vector<TH1*> gens = getGenHistos(predictions, lepSel, variable, true, true);
	 gens.resize(3, 0);
	 TH1 *hMadGenCrossSection = gens[0];
	 //hMadGenCrossSection->SetZTitle("aMC@NLO + PY8 (#leq 2j NLO + PS)");
	 TH1 *hGen1CrossSection = gens[1];
	 TH1 *hGen2CrossSection = gens[2];

	 int ipred = -1;
	 for(auto g: gens){
	     ++ipred;
	     if(g) g->SetZTitle(getLegendGen(predictions[ipred].c_str()));
	 }

	 // Here is an array of TH1D to store the various unfolded data:
	 // 0 - Central, 
	 // 1 - JES up, 2 - JES down, 
	 // 3 - PU up, 4 - PU down, 
	 // 5 - JER up, 6 - JER down, 
	 // 7 - XSEC up, 8 - XSEC down
	 // 9 - LES up, 10 - LES down
	 // 11 - LER up, 12 -LER down
	 // 13 - Lumi up, 14 - Lumi down
	 // 15 - SF up, 16 - SF down
	 // 17 - SherpaUnf
	 TString name[] = {"Central", "JESUp", "JESDown", "PUUp", "PUDown", "JERUp", "JERDown", 
			   "XSECUp", "XSECDown", "LESUp", "LESDown", "LERUp", "LERDown",
			   "LumiUp", "LumiDown", "SFUp", "SFDown", "AltUnf"};
	 TH1D *hUnfData[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,};
	 TH2D *hUnfDataStatCov[18] =  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,};
	 TH2D *hUnfMCStatCov[18] =  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,};

	 int nIter[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
 //	int svdKterm(0);
 //	if (lepSel == "DMu")  
 //	    svdKterm = VAROFINTERESTZJETS[i].MuSVDkterm;
 //	else if (lepSel == "DE")   
 //	    svdKterm = VAROFINTERESTZJETS[i].ESVDkterm;
 //	else {
 //	    cerr << "Error: algo " << algo << " or lepSel " << lepSel << "invalid\n"; 
 //	    cerr << "Aborting...\n";
 //	    return;
 //	}

	 bool logy = VAROFINTERESTZJETS[i].logy;

	 int nSysts = fDYJets[9] ? 18 : 17;

	 //--- Unfold the Data histograms for each systematic ---
	 int fixNIterTo = cachedNIters ? -1 : 0;
	 for (unsigned short iSyst = 0; iSyst < nSysts; ++iSyst) {

	     if(iSyst ==  17 && !withUnfUnc) continue;

	     bool isDisabled = false;
	     for(auto disSys: disabledUncs) if (disSys==iSyst) isDisabled = true;
	     
	     if(isDisabled){
		 std::cout << "Uncertainty " << iSyst << " is disabled by configuration.\n";
		 continue;
	     }

	     if(iSyst != 0 && whichSyst >= 0 && iSyst != whichSyst) continue;

	     std::cout << "\n----------------------------------------------------------------------\nProcessing uncertainty variation number " << iSyst << "(0: central)...\n\n";

	     //--- only JES up and down (iSyst = 1 and 2) is applied on data ---
	     unsigned short iData = (iSyst == 1 || iSyst == 2) ? iSyst : 0;
	     unsigned short iBg = 0;
	     if (iSyst == 0 || iSyst == 1 || iSyst == 2 || iSyst == 5 || iSyst == 6 || iSyst == 11 || iSyst == 12 || iSyst == 17) iBg = 0; // Central, JES, JER, LER, Sherpa
	     else if (iSyst == 3 || iSyst == 4) iBg = iSyst - 2; // PU
	     else if (iSyst == 7 || iSyst == 8 || iSyst == 9 || iSyst == 10) iBg = iSyst - 4; // XSec, LES
	     else if (iSyst == 13 || iSyst == 14 || iSyst == 15 || iSyst == 16) iBg = iSyst - 6;  // Lumi, SF

	     TH1D *hRecDataMinusFakes = (TH1D*) hRecData[iData]->Clone(TString::Format("hBkgSubData_%s", hRecData[iData]->GetName()));
	     //	if(hRecDataMinusFakes->GetBinContent(1) > 0){
	     //	  std::cout <<  __FILE__ << __LINE__ << ": "
	     //		    << "hRec: " << hRecDataMinusFakes->GetName()
	     //		    <<  " " << hRecDataMinusFakes->GetBinError(1)
	     //	    / sqrt(hRecDataMinusFakes->GetBinContent(1)) << "\n";
	     //	}
	     //	std::cerr << "DEBUG: hRecData[" << iData << "]->GetEntries() = "
	     //		  << hRecDataMinusFakes->GetEntries()
	     //		  << ", nbins: " << hRecDataMinusFakes->GetNbinsX()
	     //		  <<"\n"
	     //		  << "DEBUG: hRecSumBg[" << iData << "]->GetEntries() = "
	     //		  << hRecSumBg[iData]->GetEntries()
	     //		  << ", nbins: " << hRecSumBg[iData]->GetNbinsX()
	     //		  <<"\n"
	     //		  << "DEBUG: hPurity[" << iData << "]->GetEntries() = "
	     //		  << hPurity[iData]->GetEntries()
	     //		  << ", nbins: " << hPurity[iData]->GetNbinsX()
	     //		  <<"\n"
	     //		  << "hRecDYJets[iData]->GetNbins() = " << hRecDYJets[iData]->GetNbinsX()
	     //		  << std::endl;

	     hRecDataMinusFakes->Add(hRecSumBg[iBg], -1);
	     //	if(hRecDataMinusFakes->GetBinContent(1) > 0){
	     //	  std::cout <<  __FILE__ << __LINE__ << ": "
	     //		    << "hRecMinusBg: " << hRecDataMinusFakes->GetName()
	     //		    << " " <<  hRecDataMinusFakes->GetBinError(1)
	     //	    / sqrt(hRecDataMinusFakes->GetBinContent(1)) << "\n";
	     //	}

	     //Applies purity corrections (aka fake corrections):
	     RemoveFakes(hRecDataMinusFakes, hFakDYJets[iSyst], hPurity[iSyst]);

	     //FIXME: use independent samples for subtraction
	     TH1D *hRecDataMinusFakesOdd;
	     TH1D *hRecDataMinusFakesEven;
	     //FIXME: remove false!
	     if(iData==0 && hRecDataOdd && hRecDataEven && false){
		 hRecDataMinusFakesOdd = (TH1D*) hRecDataOdd->Clone();
		 hRecDataMinusFakesOdd->Add(hRecSumBg[iBg], -0.5);
		 //RemoveFakes(hRecDataMinusFakesOdd, hFakDYJets[iSyst], hPurity[iSyst]);
		 TH1* halfFakes = (TH1*) hFakDYJets[iSyst]->Clone();
		 halfFakes->Scale(0.5);
		 RemoveFakes(hRecDataMinusFakesOdd, halfFakes, hPurity[iSyst]);
		 hRecDataMinusFakesEven = (TH1D*) hRecDataEven->Clone();
		 hRecDataMinusFakesEven->Add(hRecSumBg[iBg], -0.5);
		 //RemoveFakes(hRecDataMinusFakesEven, hFakDYJets[iSyst], hPurity[iSyst]);
		 RemoveFakes(hRecDataMinusFakesEven, halfFakes, hPurity[iSyst]);
		 delete halfFakes;
	     } else{
		 hRecDataMinusFakesEven = hRecDataMinusFakesOdd = 0;
	     }
	     //	if(hRecDataMinusFakes->GetBinContent(1) > 0){
	     //	  std::cout <<  __FILE__ << __LINE__ << ": "
	     //		    << "hRecMinusBgMinusFake: " << hRecDataMinusFakes->GetName()
	     //		    << " " << hRecDataMinusFakes->GetBinError(1)
	     //	    / sqrt(hRecDataMinusFakes->GetBinContent(1))
	     //		    << "\t" << hRecDataMinusFakes->GetBinContent(1) /  hRecData[iData]->GetBinContent(1)
	     //	    	    << "\t" << hRecDataMinusFakes->GetBinError(3)
	     //	    / sqrt(hRecDataMinusFakes->GetBinContent(3))
	     //		    << "\t" << hRecDataMinusFakes->GetBinContent(3) /  hRecData[iData]->GetBinContent(3)
	     //		    << "\n";
	     //	}


	     if (iSyst == 17) cout << "Aternative unfolding" << endl;
	     std::cout << "Starting unfolding of " << variable << " "
		       << name[iSyst] << " for  "<< lepSel << " channel." << "\n";

	     if(hRecDataMinusFakes->GetEntries() == 0){
		 std::cerr << "Warning: histogram " << hRecDataMinusFakes->GetName()
			   << " has no entries. Its unfolding will be skipped.\n";
		 continue;
	     }

	     nIter[iSyst] = UnfoldData(unfCfg, lepSel, variable, respDYJets[iSyst],
				       hRecDataMinusFakes, hUnfData[iSyst], hUnfDataStatCov[iSyst], hUnfMCStatCov[iSyst],
				       name[iSyst], integratedLumi, unfoldDir,
				       hRecDYJets[iSyst], hGenDYJets[iSyst], logy,
				       hRecDataMinusFakesOdd, hRecDataMinusFakesEven, fixNIterTo,
				       outputFileName + "_niters.txt");
	     //The number of unfolding iterations is fixed to the value used for the central value.

	     if(iSyst == 0){
		 fixNIterTo = nIter[0];
		 //	    if(variable=="VisPt_Zinc3jetQun" && lepSel == "DMu"){
		 //		std::cout << "Forcing number of iteration of VisPt_Zinc3jetQun of DMu to 4";
		 //		fixNIterTo = 4;
		 //	    }
		 //	    if(variable=="JetsHT_Zinc1jet" ){
		 //		std::cout << "Forcing number of iteration of JetsHT_Zinc1jet to 4";
		 //		fixNIterTo = 6;
		 //	    }
		 //	    if(variable=="JZB_ptHigh" ){
		 //		std::cout << "Forcing number of iteration of JZB_ptHigh to 4";
		 //		fixNIterTo = 6;
		 //	    }
		 //	    std::ofstream fniter(outputFileName + "_niters.txt");
		 //	    fniter << variable << "\t" << fixNIterTo << "\n";
	     }

	     //--- save the unfolded histograms ---
	     outputRootFile->cd();
	     if(pseudoData){
		 hUnfData[iSyst]->SetZTitle(TString::Format("Pseudodata %s", pseudoDataSamplePrettyName.c_str()).Data());
	     } else{
		 hUnfData[iSyst]->SetZTitle("Measurement");
	     }
	     if(hUnfData[iSyst]) hUnfData[iSyst]->Write();
	 }
	 //----------------------------------------------------------------------------------------- 

	 if (doNormalized) {
	     for(int i = 0; i < nSysts; i++)
		 {
		     if(!hUnfData[i]) continue;
		     double totUnfData = hUnfData[i]->Integral("width"); // normalize to central or itself? 
		     hUnfData[i]->Scale(1.0/totUnfData);
		     if (i == 0) {
			 hUnfDataStatCov[0]->Scale(1.0/pow(totUnfData, 2));
			 hUnfMCStatCov[0]->Scale(1.0/pow(totUnfData, 2));
		     }
		 }
	 }

	 bool noNull = true;
	 if(1 <=  normalizeToBin && normalizeToBin <= hUnfData[0]->GetNbinsX()
	    && hUnfData[0]->GetBinContent(normalizeToBin) > 0){
	     std::cout << "Normalizing " << variable << " to bin "
		       << hUnfData[0]->GetXaxis()->GetBinLowEdge(normalizeToBin)
		       << "..." << hUnfData[0]->GetXaxis()->GetBinUpEdge(normalizeToBin) << "\n";
	     noNull = false;
	     for(int i = 1; i < nSysts; i++)
		 {
		     if(!hUnfData[i]) continue;
		     double s = hUnfData[0]->GetBinContent(normalizeToBin)
			 / hUnfData[i]->GetBinContent(normalizeToBin);
		     hUnfData[i]->Scale(s);
		     hUnfDataStatCov[i]->Scale(pow(s, 2));
		     hUnfMCStatCov[i]->Scale(pow(s, 2));
		 }   
	 }

	 for (unsigned short iSyst = 0; iSyst < nSysts; ++iSyst) {
	     outputRootFile->cd(); 
	     if(hUnfData[iSyst]) hUnfData[iSyst]->Write();
	 }
	 //--- Now create the covariance matrices ---
	 //	 TH2 *hCov[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,};
	 std::vector<TH2*> hCov(kUncCnt, 0);
	 //	 int nCovs = nSysts > 17 ? 11 : 10;

	 //	 if(hUnfDataStatCov[0]) hCov[kStat] = (TH2*) hUnfDataStatCov[0]->Clone("CovDataStat");
	 //	 if(hUnfMCStatCov[0])   hCov[kUnfStat] = (TH2*) hUnfMCStatCov[0]->Clone("CovMCStat");
	 if(hUnfDataStatCov[0]) hCov[kStat] = (TH2*) hUnfDataStatCov[0]->Clone(covName(kStat));
	 if(hUnfMCStatCov[0])   hCov[kUnfStat] = (TH2*) hUnfMCStatCov[0]->Clone(covName(kUnfStat));
	 
	 if(hUnfData[kJESup])   hCov[kJES]    = makeCovFromUpAndDown(hUnfData[0], hUnfData[kJESup], hUnfData[kJESdwn], covName(kJES), noNull);
	 if(hUnfData[kPUup])  hCov[kPU]     = makeCovFromUpAndDown(hUnfData[0], hUnfData[kPUup], hUnfData[kPUdwn], covName(kPU), noNull);
	 if(hUnfData[kJERup])  hCov[kJER]    = makeCovFromUpAndDown(hUnfData[0], hUnfData[kJERup], hUnfData[kJERdwn], covName(kJER), noNull);
	 if(hUnfData[kXsecUp])  hCov[kXsec]   = makeCovFromUpAndDown(hUnfData[0], hUnfData[kXsecUp], hUnfData[kXsecDwn], covName(kXsec), noNull);
	 if(hUnfData[kLESup])  hCov[kLES]    = makeCovFromUpAndDown(hUnfData[0], hUnfData[kLESup], hUnfData[kLESdwn], covName(kLES), noNull);
	 if(hUnfData[kLERup]) hCov[kLER]    = makeCovFromUpAndDown(hUnfData[0], hUnfData[kLERup], hUnfData[kLERdwn], covName(kLER), noNull);
	 if(hUnfData[kLumiUp]) hCov[kLumi]   = makeCovFromUpAndDown(hUnfData[0], hUnfData[kLumiUp], hUnfData[kLumiDwn], covName(kLumi), noNull);
	 if(hUnfData[kSFup]) hCov[kSF]     = makeCovFromUpAndDown(hUnfData[0], hUnfData[kSFup], hUnfData[kSFdwn], covName(kSF), noNull);
	 if(hUnfData[kAltUnf]) hCov[kUnfSys] = makeCovFromUpAndDown(hUnfData[0], hUnfData[kAltUnf], 0, covName(kUnfSys), noNull);

	 if(hUnfMCStatCov[0]) hCov[kTotSys] = (TH2*) hCov[kUnfStat]->Clone(covName(kTotSys));

	 if(hCov[kTotSys]){
	     for (int i = kUnfStat + 1; i < kTotSys; ++i){
		 if(hCov[i]) hCov[kTotSys]->Add(hCov[i]);
	     }
	 }

	 if (doNormalized) {
	     double Madtot = hMadGenCrossSection->Integral("width");
	     hMadGenCrossSection->Scale(1.0/Madtot);
	     if(hGen1CrossSection){
		 double gen1tot = hGen1CrossSection->Integral("width");
		 hGen1CrossSection->Scale(1.0/gen1tot);	
	     }
	     if(hGen2CrossSection){
		 double gen2tot = hGen2CrossSection->Integral("width");
		 hGen2CrossSection->Scale(1.0/gen2tot);	
	     }
	 }

	 int upperBin = hUnfData[0]->GetXaxis()->GetNbins() -  nLastBinsToSkip;
	 int lowerBin = 1 + nFirstBinsToSkip;
	 hUnfData[0]->GetXaxis()->SetRange(lowerBin, upperBin);
	 TCanvas *crossSectionPlot = makeCrossSectionPlot(lepSel, integratedLumi,
							  variable, doNormalized,
							  hUnfData[0], hCov[11],
							  predictions);
	 crossSectionPlot->Draw();

	 crossSectionPlot->SaveAs(outputFileName + ".png");
	 crossSectionPlot->SaveAs(outputFileName + ".pdf");
	 crossSectionPlot->SaveAs(outputFileName + ".eps");
	 crossSectionPlot->SaveAs(outputFileName + ".ps");
	 crossSectionPlot->SaveAs(outputFileName + ".C");
	 crossSectionPlot->SaveAs(outputFileName + "_canvas.root");

	 if(whichSyst < 0){
	     createSystPlots(outputFileName, sysPlotDir, variable, lepSel, hUnfData, logy);
	     //--- print out break down of errors ---
	     //for (int i = 2; i <= nCovs; ++i) {
	     //	 cout << hUnfData[0]->GetBinContent(i);
	     //	 for (int j = 0; j <= 11; ++j) {
	     //	     if(hCov[j]){
	     //		 cout << " +/- " << sqrt(hCov[j]->GetBinContent(i,i))*100./hUnfData[0]->GetBinContent(i) << "%";
	     //	     }
	     //	 }
	     //	 cout << endl;
	     //}
////	     for(int i = 1; i <= hCov[kLumi]->GetNbinsX(); ++i){
////		 std::cout << ">>>> " << i << "\t" << sqrt(hCov[kLumi]->GetBinContent(i,i)) 
////			   << "\t" << sqrt(hCov[kLumi]->GetBinContent(i,i))  / hUnfData[0]->GetBinContent(i)
////			   << "\n";
////	     }
	     createTable(outputFileName + "_withLERS", lepSel, variable,
			 doNormalized, hUnfData[0], hCov, true);
	     createTable(outputFileName, lepSel, variable,
			 doNormalized, hUnfData[0], hCov, false);
	 }

	 if (variable.Index("ZNGoodJets_Zexc") >= 0) {
	     createInclusivePlots(doNormalized, outputFileName, lepSel,
				  hUnfData[0], &hCov[0], predictions, integratedLumi);
	 }
	 //--------------------------------------

	 //--- Save other things --- 
	 outputRootFile->cd();
	 hRecData[0]->Write("hRecDataCentral");
	 hRecSumBg[0]->Write("hRecSumBgCentral");
	 hRecDYJets[0]->Write("hRecDYJetsCentral");
	 hGenDYJets[0]->Write("hGenDYJetsCentral");
	 hMadGenCrossSection->Write("hMadGenDYJetsCrossSection");
	 if(hGen1CrossSection) hGen1CrossSection->Write("hGen1DYJetsCrossSection");
	 if(hGen2CrossSection) hGen2CrossSection->Write("hGen2DYJetsCrossSection");
	 respDYJets[0]->Write("respDYJetsCentral");
	 for (int i = 0; i <= 11; ++i) {
	     if(hCov[i]){
		 hCov[i]->Write();
	     }
	 }
	 TParameter<double> pIntegratedLumi("integratedLumi", integratedLumi);
	 pIntegratedLumi.Write();
	 TParameter<int> pNIter("nIter", nIter[0]);
	 pNIter.Write();
	 std::cout << "number of iterations: " << nIter[0] << "\n";
	 crossSectionPlot->Write();
	 //----------------------------------------------------------------------------------------- 

	 outputRootFile->Close();

	 //if (end == start + 1) system("display " + outputFileName + ".png &");
	 //if (end == start + 1 && variable == "ZNGoodJets_Zexc") system("display " + outputFileName.ReplaceAll("ZNGoodJets_Zexc", "ZNGoodJets_Zinc") + ".png &");

	 gettimeofday(&t1, 0);
	 std::cout << __FILE__ << ":" << __LINE__ << ": unfolding of  variable " << variable << " for channel "
		   << lepSel << " and systematic #" << whichSyst << " took "
		   << ((t1.tv_sec - t0.tv_sec) * 1e3 + (t1.tv_usec - t0.tv_usec) * 1e-3)
		   << " ms.\n";

	 //--- Close all files ----------------------------------------------------------------------
	 closeAllFiles(fData, fDYJets, fBg, nbg);
     }

     //------------------------------------------------------------------------------------------ 

     std::ofstream f(unfoldDir + "/" + "lastUnfConfig.txt");
     f << "#Files automatically genertaed by UnfoldingZJets. The file will be overwritten at next execution.\n\n";
     unfCfg.dumpRetrieved(f);
 }

 void createSystPlots(TString outputFileName, TString sysPlotDir, TString variable, 
		      TString lepSel, TH1D *hUnfData[], bool logy)
 {

     // 0 - Central, 
     // 1 - JES up, 2 - JES down, 
     // 3 - PU up, 4 - PU down, 
     // 5 - JER up, 6 - JER down, 
     // 7 - XSEC up, 8 - XSEC down
     // 9 - LES up, 10 - LES down
     // 11 - LER up, 12 - LER down
     // 13 - Lumi up, 14 - Lumi down
     // 15 - SF up, 16 - SF down
     // 17 - Alernative unfolding
     TString syst[] = {"JES", "PU", "JER", "XSec", "LES", "LER", "Lumi", "S.F."};
     for (int i = 1; i < 16; i += 2) {

	 if(hUnfData[i] == 0) continue;

	 TH1D *hCent = (TH1D*) hUnfData[0]->Clone();
	 hCent->SetMarkerColor(kBlack);
	 hCent->SetMarkerStyle(20);
	 if (variable == "ZNGoodJets_Zexc") hCent->GetXaxis()->SetRangeUser(1, 8);
	 if (variable.Index("JetPt_Zinc") >= 0) hCent->GetXaxis()->SetRangeUser(30, hCent->GetXaxis()->GetXmax());
	 hCent->GetXaxis()->SetLabelSize(0);
	 hCent->GetYaxis()->SetTitle("d#sigma");
	 hCent->GetYaxis()->SetTitleSize(0.05);
	 hCent->GetYaxis()->SetTitleOffset(0.9);

	 TH1D *hUp = (TH1D*) hUnfData[i]->Clone();
	 hUp->SetLineColor(kGreen+2);
	 hUp->SetLineWidth(2);
	 TH1D *hDown = (TH1D*) hUnfData[i+1]->Clone();
	 hDown->SetLineColor(kBlue);
	 hDown->SetLineWidth(2);

	 std::unique_ptr<TCanvas> c(new TCanvas(variable + " - " + syst[i/2], variable + " - " + syst[i/2], 700, 900));
	 c->cd();

	 TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
	 pad1->SetTopMargin(0.11);
	 pad1->SetBottomMargin(0.02);
	 pad1->SetRightMargin(0.03);
	 pad1->SetLeftMargin(0.15);
	 pad1->SetTicks();
	 pad1->SetLogy(logy ? kTRUE : kFALSE);
	 pad1->Draw();
	 pad1->cd();

	 TLegend *leg = new TLegend(0.8, 0.7, 0.95, 0.86);
	 leg->SetBorderSize(0);
	 leg->SetFillStyle(0);
	 hCent->DrawCopy("e");
	 hUp->DrawCopy("samehist");
	 hDown->DrawCopy("samehist");
	 leg->AddEntry(hUp, syst[i/2] + " Up", "l");
	 leg->AddEntry(hCent, "Central", "lp");
	 leg->AddEntry(hDown, syst[i/2] + " Down", "l");
	 leg->Draw();
	 pad1->Draw();
	 fixYscale();
	 c->cd();

	 TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
	 pad2->SetTopMargin(0.02);
	 pad2->SetBottomMargin(0.3);
	 pad2->SetRightMargin(0.03);
	 pad2->SetLeftMargin(0.15);
	 pad2->SetGridy();
	 pad2->SetTicks();
	 pad2->Draw();
	 pad2->cd();
	 hUp->Divide(hCent);
	 hDown->Divide(hCent);
	 hCent->Divide(hCent);

	 double maxRatio = max(max(hCent->GetMaximum(), hUp->GetMaximum()), hDown->GetMaximum());
	 double minRatio = min(min(hCent->GetMinimum(), hUp->GetMinimum()), hDown->GetMinimum());
	 int nBins = hCent->GetNbinsX();
	 for (int j = 1; j <= nBins; j++) {
	     maxRatio = max(maxRatio, 1+hCent->GetBinError(j));
	     minRatio = min(minRatio, 1-hCent->GetBinError(j));
	 }
	 hCent->GetYaxis()->SetRangeUser(minRatio-0.02*(maxRatio-minRatio), maxRatio+0.02*(maxRatio-minRatio));
	 hCent->GetYaxis()->SetRangeUser(0.82, 1.18);
	 if (variable == "ZNGoodJets_Zexc") hCent->GetXaxis()->SetRangeUser(1, 8);
	 hCent->SetTitle(" ");
	 hCent->GetYaxis()->SetTitle("Ratio to Central");
	 hCent->GetYaxis()->SetTitleSize(0.08);
	 hCent->GetYaxis()->SetTitleOffset(0.9);
	 hCent->GetYaxis()->CenterTitle();
	 hCent->GetYaxis()->SetLabelSize(0.08);
	 hCent->GetYaxis()->SetLabelOffset(0.014);
	 hCent->GetXaxis()->SetTitleSize(0.13);
	 //hCent->GetXaxis()->SetLabelSize(0.13);
	 hCent->GetXaxis()->SetLabelSize(0.08);
	 hCent->GetXaxis()->SetLabelOffset(0.012);

	 hCent->DrawCopy("e");
	 hUp->Draw("histsame");
	 hDown->Draw("histsame");
	 c->cd();
	 c->Update();
	 c->Draw();

	 TString systStr = syst[i/2];
	 if (systStr == "S.F.") systStr = "SF";
	 system("mkdir " + sysPlotDir);
	 c->SaveAs(sysPlotDir + "/" + lepSel + "_" + variable + "_" + systStr + ".png");
	 c->SaveAs(sysPlotDir + "/" + lepSel + "_" + variable + "_" + systStr + ".ps");
	 c->SaveAs(sysPlotDir + "/" + lepSel + "_" + variable + "_" + systStr + ".eps");
	 c->SaveAs(sysPlotDir + "/" + lepSel + "_" + variable + "_" + systStr + ".pdf");
	 c->SaveAs(sysPlotDir + "/" + lepSel + "_" + variable + "_" + systStr + ".C");
	 c->SaveAs(sysPlotDir + "/" + lepSel + "_" + variable + "_" + systStr + ".root");	 
     }
 }

void createInclusivePlots(bool doNormalized, TString outputFileName, TString lepSel, TH1 *hUnfData, TH2 *hCov[],
			  const std::vector<std::string>& predictions, double integratedLumi){
    //    std::cerr << "createInclusivePlots disabled!" << __FILE__ << __LINE__ << "\n\n";
    //return;

    TH1D *hInc = (TH1D*) hUnfData->Clone("ZNGoodJets_Zinc");

    std::vector<TH1*> hGens = getGenHistos(predictions, lepSel, "ZNGoodJets_Zinc", true);
    hGens.resize(3, 0);    

    TH1D *hIncMad = hGens[0] ? (TH1D*) hGens[0]->Clone("ZNGoodJets_Zinc_Mad") : 0;
    TH1D *hIncShe = hGens[1] ? (TH1D*) hGens[1]->Clone("ZNGoodJets_Zinc_She") : 0;
    TH1D *hIncPow = hGens[2] ? (TH1D*) hGens[2]->Clone("ZNGoodJets_Zinc_Pow") : 0;

    //    const int kTot = 11;
    //TH2 *hCovInc[kTotSys+1] = {NULL};
    std::vector<TH2*> hCovInc(kTotSys + 1, 0);
    for(unsigned i = 0; i < kUncCnt; ++i){
	//	if(hCov[i]) hCovInc[i] = (TH2*) hCov[i]->Clone(hCov[i[-TString::Format("Cov%s", uncShortNames[i]));
	if(hCov[i]) hCovInc[i] = (TH2*) hCov[i]->Clone(); 
    }
    int nBins = hInc->GetNbinsX();
    for (int i = 1; i <= nBins; i++) {
	double binSum = 0;
	double binSumMad = 0;
	double binSumShe = 0;
	double binSumPow = 0;
	double binStatError2 = 0;
	double binStatMadError2 = 0;
	double binStatSheError2 = 0;
	double binStatPowError2 = 0;
	for (int j = i; j <= nBins; j++) {
	    binSum += hInc->GetBinContent(j);
	    if(hIncMad){
		binSumMad += hIncMad->GetBinContent(j);
		binStatMadError2 += pow(hIncMad->GetBinError(j), 2);
	    }
	    if(hIncShe){ 
		binSumShe += hIncShe->GetBinContent(j);
		binStatSheError2 += pow(hIncShe->GetBinError(j), 2);
	    }
	    if(hIncPow){
		binSumPow += hIncPow->GetBinContent(j);
		binStatPowError2 += pow(hIncPow->GetBinError(j), 2);
	    }
	    binStatError2 += pow(hInc->GetBinError(j), 2);
	}
	hInc->SetBinContent(i, binSum);
	if(hIncMad) hIncMad->SetBinContent(i, binSumMad);
	if(hIncShe) hIncShe->SetBinContent(i, binSumShe);
	if(hIncPow) hIncPow->SetBinContent(i, binSumPow);
	hInc->SetBinError(i, sqrt(binStatError2));
	if(hIncMad) hIncMad->SetBinError(i, sqrt(binStatMadError2));
	if(hIncShe) hIncShe->SetBinError(i, sqrt(binStatSheError2));
	if(hIncPow) hIncPow->SetBinError(i, sqrt(binStatPowError2));
    } //next i (bin)

    //Covariance matrix.
    //We can write:
    //    Y_inc = A * Y_exc, with Y_inc and Y_exc the vector of the respective
    //                       distribution bin contents
    //                       and A_ij = 1 if j >=i, 0 otherwise
    //   => Cov_inc = A * Cov_exc * A^{T}
    for(int m = 0; m < kUncCnt; ++m){
	if(hCov[m]==0) continue;
	for(int i = 1; i <= nBins; ++i){
	    for(int j = 1; j <= nBins; ++j){
		double c = 0;
		for(int k = i; k <= nBins; ++k){
		    for(int l = j; l <= nBins; ++l){
			c += hCov[m]->GetBinContent(k, l);
		    }
		}
		hCovInc[m]->SetBinContent(i, j, c);
	    }
	}
    }

    TCanvas *crossSectionPlot = makeCrossSectionPlot(lepSel, integratedLumi, 
						     TString("ZNGoodJets_Zinc"), doNormalized, hInc, hCovInc[11],
						     predictions); 
    outputFileName.ReplaceAll("ZNGoodJets_Zexc", "ZNGoodJets_Zinc");
    crossSectionPlot->Draw();
    crossSectionPlot->SaveAs(outputFileName + ".png");
    crossSectionPlot->SaveAs(outputFileName + ".pdf");
    crossSectionPlot->SaveAs(outputFileName + ".eps");
    crossSectionPlot->SaveAs(outputFileName + ".ps");
    crossSectionPlot->SaveAs(outputFileName + ".C");
    crossSectionPlot->SaveAs(outputFileName + "_canvas.root");
    createTable(outputFileName + "_withLERS", lepSel, TString("ZNGoodJets_Zinc"), doNormalized,
		hInc, hCovInc, true);
    createTable(outputFileName, lepSel, TString("ZNGoodJets_Zinc"), doNormalized,
		hInc, hCovInc, false);
}


int UnfoldData(const SectionedConfig& unfCfg, const TString lepSel, const char* variable, RooUnfoldResponse *resp,
	       TH1D* hRecDataMinusFakes, TH1D* &hUnfData, TH2D* &hUnfDataStatCov, TH2D* &hUnfMCStatCov,
	       TString name, double integratedLumi, const TString& unfoldDir,
	       TH1* hRecDYJets, TH1* hGenDYJets, bool logy,
	       TH1D *hRecDataMinusFakesOdd, TH1D *hRecDataMinusFakesEven, int fixNIterTo,
	       const char* outputFileName){
    
    //-- sanity check
    //response matrix: x-axis = measured, y-axis = truth
    const TH2& hr = *resp->Hresponse();
    if(!isSameBinning(*hr.GetXaxis(), *hRecDataMinusFakes->GetXaxis())){
	std::cerr << __FILE__ << ":" << __LINE__ << ". Fatal error: axes of measured and hronse matrix histograms are not conistents.\n";
	abort();
    }

    if(!isSameBinning(*hr.GetXaxis(), *hr.GetYaxis())){
	std::cerr << __FILE__ << ":" << __LINE__ << ". Fatal error: x- and y- axes of the response matrices differs. If it's on purpose, please edit the code to suppress this sanity check.\n";
	abort();
    }

    //--- make sure we use OverFlow (should already be set to true) ---
    resp->UseOverflow();

    //    TString variable = TString(hRecDataMinusFakes->GetName());
    TString section = TString::Format("%s_%s", lepSel.Data(), variable);

    const TString algo = unfCfg.get<TString>(section.Data(), "algo");
    int svdKterm       = unfCfg.get<int>(section.Data(), "svdKterm");

    //--- Set the required unfolding algorithm ---
    RooUnfold::Algorithm alg;
    if (algo == "Bayes") {
	alg = RooUnfold::kBayes;
    }
    else if (algo == "SVD") {
	alg = RooUnfold::kSVD;
    }
    else {
	cerr << "Error: the specified algo: " << algo << " is not implemented!" << endl;
	cerr << "       I will proceed with kBayes algo" << endl;
	alg = RooUnfold::kBayes;
    }

    //    std::cout << "-----------------------" << std::endl;
    TString unfoldCheckDir = unfoldDir.Strip(TString::kTrailing, '/') + "Check";
    system(TString("mkdir ") + unfoldCheckDir + "/");
    std::unique_ptr<TFile> f(new TFile(unfoldCheckDir + "/" + lepSel + "_" + variable + "_" + name + ".root", "RECREATE"));
    f->cd();

    TH2D* hresp = (TH2D*) resp->Hresponse()->Clone(TString::Format("hResp%s%s", variable, name.Data()));
    hresp->Write();
    //    std::cout << "Response matrix dimensions: " << resp->Mresponse().GetNrows()
    //	      << "x" << resp->Mresponse().GetNcols() << "\n";
    TDecompSVD svd(resp->Mresponse());
    svd.Decompose();
    double matrixCond = svd.Condition();
    TParameter<double>("respMatCondition", matrixCond).Write();
    TParameter<int>("respMatIsSingular", svd.TestBit(21)).Write();

    double lumiUnc = cfg.getD("lumiUnc", 0.027);


    //bool svd_unfold = cfg.getB("svdUnfold", false);
    //bool tsvd_unfold = cfg.getB("tsvdUnfold", false);
    //bool invert_unfold = cfg.getB("invertUnfold", false);
    //bool binByBin_unfold = cfg.getB("binByBinUnfold", false);
    //int nSkipFirstJetPtBins = cfg.getI("nSkipFirstJetPtBins", 2);
    //verbosity = cfg.getI("unfoldingVerbosity", 1);    
    //int xvalIter = cfg.getI("xvalIter", 0);
    //int minIter = cfg.getI("minIter", 2);
    //int maxIter = cfg.getI("maxIter", 20);
    //bool smoothPrior = cfg.getB("unfSmootPrior", false);
    //int nTestIterMax = std::max(hRecDataMinusFakes->GetNbinsX(), maxIter);
    //bool forceNitersStudy = cfg.getB("forceNitersStudy", false);


    bool svd_unfold          = cfg.getUnf(lepSel, variable, "svdUnfold", false);
    bool tsvd_unfold         = cfg.getUnf(lepSel, variable, "tsvdUnfold", false);
    bool invert_unfold       = cfg.getUnf(lepSel, variable, "invertUnfold", false);
    bool binByBin_unfold     = cfg.getUnf(lepSel, variable, "binByBinUnfold", false);
    //int nSkipFirstJetPtBins  = cfg.getUnf(lepSel, variable, "nSkipFirstJetPtBins", 2);
    verbosity                = cfg.getUnf(lepSel, variable, "unfoldingVerbosity", 1); 
    int xvalIter 	     = cfg.getUnf(lepSel, variable, "xvalIter", 0);
    int nFixedIters          = cfg.getUnf(lepSel, variable, "nFixedIters", 0);
    int minIter  	     = cfg.getUnf(lepSel, variable, "minIter", 2);
    int maxIter  	     = cfg.getUnf(lepSel, variable, "maxIter", 20);
    bool smoothPrior 	     = cfg.getUnf(lepSel, variable, "unfSmootPrior", false);
    bool forceNitersStudy    = cfg.getUnf(lepSel, variable, "forceNitersStudy", false);
    int nFirstBinsToSkip     = cfg.getUnf(lepSel, variable, "nFirstBinsToSkip", 0);
    int nLastBinsToSkip      = cfg.getUnf(lepSel, variable, "nLastBinsToSkip", 0);
    int ntoys                = cfg.getUnf(lepSel, variable, "nToysForChi2", 0);  
    bool toyIIrecoResampling = cfg.getUnf(lepSel, variable, "toyIIrecoResampling", false); 

    //    int nTestIterMax 	     = std::max(hRecDataMinusFakes->GetNbinsX(), maxIter);
    int nTestIterMax = maxIter;
    int chosenIter = -1;

    if(xvalIter == 6 && nFixedIters < 1){
	std::cerr << "Fatal error. Configuration error for " << lepSel << "_" << variable
		  << ". xvalIter set to fixed itertation mode while nFixedIters is either not set or "
		  << "set to a value smaller than 1 for this variable. Aborts\n";
	abort();
    }
     
    if(fixNIterTo > 0){
	chosenIter = fixNIterTo;
    } else if(fixNIterTo < 0){ //use previous number of iterations.
	std::ifstream f(outputFileName);
	if(f.good()){
	    f >> chosenIter;
	}
    }

    if(!forceNitersStudy && (xvalIter == 6)){
	chosenIter = nFixedIters;
	fixNIterTo = nFixedIters;
    }

    if(chosenIter <= 0){
	//test different number of iterations
	int finalNIter = -1;
	int nIter = 99;
	//    int finalNIterXval = -1;
	//    int nIterXval = 99;
	int nIterXval = -1;
	int nIterXvalMin = 0;
	int nIterToy = -1;
	int nIterToyMin = -1;
	int nIterResMax = -1;
	int nIterResMaxXval = -1;
	int nIterResMaxToy = -1;
	int nIterResMaxToyMin = -1;
	double resMaxToyMin = 1.e9;
	int nIterResMaxToyCmp = -1;
	int nIterResMaxToyII = -1;
	int nIterToyIIMin = -1;
	bool nIterToyIIMinFound = false;
	int nIterResMaxToyIIMin = -1;
	int nIterToyII = -1;
	double resMaxToyIIMin = 1.e9;
	int nIterBiasToyIIMin = -1;
	double biasToyIIMin = 1.e9;
	double chi2ToyIIMin = 1.e9;
	double chi2XvalMin = 1.e9;
	double chi2ToyMin = 1.e9;

	TH1D *hchi2 = new TH1D("hchi2", "hchi2", nTestIterMax + 1, -0.5, nTestIterMax + .5);
	hchi2->SetTitle(TString::Format("Reco #chi^{2}/ndf for %s",
					variable));
	hchi2->GetYaxis()->SetTitle("#chi^{2}/ndf");
	hchi2->GetYaxis()->SetTitleOffset(1.40);
	hchi2->GetXaxis()->SetTitle("number of iterations");
	hchi2->GetXaxis()->CenterTitle();
	hchi2->GetXaxis()->SetNdivisions(nTestIterMax, 0, 0);
	hchi2->GetXaxis()->SetLabelSize(0.03);
	hchi2->SetLineWidth(2);


	TH1D* hchi2dataMC = (TH1D*) hchi2->Clone("hchi2dataMC");
	hchi2dataMC->SetTitle(TString::Format("#chi^{2}/ndf test of MC gen and unfolded data for %s",
					      variable));

	TH1D *hchi2Toy = new TH1D("hchi2Toy", "hchi2Toy", nTestIterMax + 1, -0.5, nTestIterMax + .5);
	hchi2Toy->SetTitle(TString::Format("Reco #chi^{2}/ndf for %s",
					   variable));
	hchi2Toy->GetYaxis()->SetTitle("#chi^{2}/ndf");
	hchi2Toy->GetYaxis()->SetTitleOffset(1.40);
	hchi2Toy->GetXaxis()->SetTitle("number of iterations");
	hchi2Toy->GetXaxis()->CenterTitle();
	hchi2Toy->GetXaxis()->SetNdivisions(nTestIterMax, 0, 0);
	hchi2Toy->GetXaxis()->SetLabelSize(0.03);
	hchi2Toy->SetLineWidth(2);    

	TH1D *hchi2ToyMcErr = new TH1D("hchi2ToyMcErr", "hchi2ToyMcErr", nTestIterMax + 1, -0.5, nTestIterMax + .5);
	hchi2ToyMcErr->SetTitle(TString::Format("Reco #chi^{2}/ndf for %s",
						variable));
	hchi2ToyMcErr->GetYaxis()->SetTitle("#chi^{2}/ndf");
	hchi2ToyMcErr->GetYaxis()->SetTitleOffset(1.40);
	hchi2ToyMcErr->GetXaxis()->SetTitle("number of iterations");
	hchi2ToyMcErr->GetXaxis()->CenterTitle();
	hchi2ToyMcErr->GetXaxis()->SetNdivisions(nTestIterMax, 0, 0);
	hchi2ToyMcErr->GetXaxis()->SetLabelSize(0.03);
	hchi2ToyMcErr->SetLineWidth(2);    


	TH1D *hchi2ToyII = new TH1D("hchi2ToyII", "hchi2ToyII", nTestIterMax + 1, -0.5, nTestIterMax + .5);
	hchi2ToyII->SetTitle(TString::Format("Gen. #chi^{2}/ndf for %s",
					     variable));
	hchi2ToyII->GetYaxis()->SetTitle("#chi^{2}/ndf");
	hchi2ToyII->GetYaxis()->SetTitleOffset(1.40);
	hchi2ToyII->GetXaxis()->SetTitle("number of iterations");
	hchi2ToyII->GetXaxis()->CenterTitle();
	hchi2ToyII->GetXaxis()->SetNdivisions(nTestIterMax, 0, 0);
	hchi2ToyII->GetXaxis()->SetLabelSize(0.03);
	hchi2ToyII->SetLineWidth(2);    

	TH1D* hchi2Xval = 0;
	TH1* hResMaxXval = 0;
	if(hRecDataMinusFakesOdd && hRecDataMinusFakesEven){
	    hchi2Xval = (TH1D*) hchi2->Clone("hchi2Xval");
	    hchi2Xval->Reset();
	    hResMaxXval = new TH1D("hResMaxXval",
				   TString::Format("Max(#chi^{2}_{ibin}) for %s - Two samples;number of iterations;max(#chi^{2}_{ibin})",
						   variable),
				   nTestIterMax + 1, -.5, nTestIterMax +  0.5);
	    hResMaxXval->SetLineWidth(2);

	}

	//double chi2Thr = 1./sqrt(2.);
	double chi2Thr = 1.;
	double chi2XvalThr = 1; //pValueToNormChi2(0.05, hRecDataMinusFakes->GetNbinsX() - 1);
	double chi2ToyThr = chi2XvalThr;
	double chi2ToyIIThr = chi2XvalThr;
	double resMaxThr = 1;
	double resMaxXvalThr = 2;
	double resMaxToyThr = 0.4;
	double resMaxToyCmpThr1 = 1.; //applied on resMax
	double resMaxToyCmpThr2 = 0.01;//applied on difference between resMaxToy and resMax
	double resMaxToyIIThr = 1;

	//    chi2Thr= chi2XvalThr = chi2ToyThr = 0.95;
	//
	//    if(hRecDataMinusFakesOdd){
	//      chi2XvalThr = MyChi2Test(hRecDataMinusFakesOdd, hRecDataMinusFakesEven);
	//    } else{
	//      chi2XvalThr = 0;
	//    }

	double chosenAlgoThr = 0.;

	std::cout << "Bayes unfolding, number of first bins to skip: " << nFirstBinsToSkip << "\n";
	std::cout << "Bayes unfolding, number of last bins to skip: " << nLastBinsToSkip << "\n";    

	RooUnfoldResponse *respBis = (RooUnfoldResponse*) resp->Clone();
	TH1D *hRecDataMinusFakesBis = (TH1D*) hRecDataMinusFakes->Clone();

	std::auto_ptr<std::vector<TH1*> > hUnfs (new std::vector<TH1*>);
	std::vector<TMatrixD> covs;

	unfoldWithErr(alg, respBis, hRecDataMinusFakesBis, nTestIterMax, smoothPrior, hUnfs.get(), &covs, 2);
	//unfoldWithErr(alg, respBis, hRecDataMinusFakesBis, nTestIterMax, smoothPrior, hUnfs.get());

	TH1* hResMax = new TH1D("hResMax", TString::Format("Reco Max(#chi^{2}_{ibin}) for %s;Iter;max(#chi^{2}_{ibin})", 
							   variable),
				hUnfs->size(), -.5, hUnfs->size() - 0.5);
	hResMax->SetLineWidth(2);

	TH1* hResMaxToy = new TH1D("hResMaxToy", TString::Format("Reco Max(#chi^{2}_{ibin}) for %s;Iter;max(#chi^{2}_{ibin})",
								 variable),
				   hUnfs->size(), -.5, hUnfs->size() - 0.5);
	hResMaxToy->SetLineWidth(2);

	TH1* hResMaxToyII = new TH1D("hResMaxToyII", TString::Format("Gen. Max(#chi^{2}_{ibin}) for %s;Iter;max(#chi^{2}_{ibin})",
								     variable),
				     hUnfs->size(), -.5, hUnfs->size() - 0.5);
	hResMaxToyII->SetLineWidth(2);

	TH1* hBiasMaxToyII = new TH1D("hBiasMaxToyII", TString::Format("Gen. Max(|pull_{ibin}|) for %s;Iter;max(|pull_{ibin}|)", variable),
				      hUnfs->size(), -.5, hUnfs->size() - 0.5);
	hBiasMaxToyII->SetLineWidth(2);

	
	double chi2dataMCThr = TH1Chi2Test(hRecDataMinusFakes, hRecDYJets, nFirstBinsToSkip, nLastBinsToSkip);

	double prevChi2 = -1;
	
	for (unsigned i = 0; i < hUnfs->size(); ++i) {
	    (*hUnfs)[i]->SetName(TString::Format("hUnf_%d", i));
	    if(i==0){
		(*hUnfs)[i]->SetTitle(TString::Format("Prior for %s %s %d iters", variable, name.Data(), i));
	    } else {
		(*hUnfs)[i]->SetTitle(TString::Format("Unfolded %s %s %d iters", variable, name.Data(), i));
	    }
	    (*hUnfs)[i]->Write();
	    TH1* hfoldUnfData = foldUnfData((*hUnfs)[i], &(covs[i]), respBis);
	    //TH1* hfoldUnfData = foldUnfData((*hUnfs)[i], 0, respBis);
	    hfoldUnfData->SetName(TString::Format("hDatafoldedBack_%d", i));
	    hfoldUnfData->SetTitle(TString::Format("hDataFoldedBack_%d", i));
	    hfoldUnfData->Write();
	    
	    if(i>0){
		TH1* hUnf2 = unfold(alg, resp, hfoldUnfData, i, smoothPrior);
		hUnf2->SetName(TString::Format("hUfu_%d", i));
		hUnf2->SetTitle(TString::Format("hUfu_%d", i));
		hUnf2->Write();
	    }

	    TH1* hRes = (TH1*) hRecDataMinusFakes->Clone(TString::Format("hRes%s_%d", name.Data(), i+1));
	    hRes->Reset();
	    hRes->SetTitle(TString::Format("%s %s residuals - iteration %d",
					   variable, name.Data(), i+1));
	    hRes->GetYaxis()->SetTitle("Residuals");
	    const int nRes = hfoldUnfData->GetNbinsX();
	    Double_t res[nRes];
	    double mychi2 = MyChi2Test(hRecDataMinusFakesBis, hfoldUnfData, nFirstBinsToSkip, nLastBinsToSkip,
				       res, hRecDataMinusFakesBis, false);
	    if(verbosity) std::cout << "Chi2/nbins of data / folded-unfolded distributions: " << mychi2 << "\n"; 
	    if(i > 0) hchi2->SetBinContent(i + 1, mychi2);
	    if (mychi2 < chi2Thr && finalNIter < 0) {
		if(i > 1 && fabs(prevChi2 - chi2Thr) < fabs(mychi2 - chi2Thr)){
		    nIter = i - 1;
		} else{
		    nIter = i;
		}
		finalNIter = nIter;
		std::cout << "Single distribution chi^2 leads to " << nIter << " iterations with a final Chi2/ndf of: " << mychi2 << std::endl;
	    }
	    prevChi2 = mychi2;

	    double dataMcChi2;
	    if(i==0){
		dataMcChi2 =  TH1Chi2Test(hRecDYJets, hRecDataMinusFakes, nFirstBinsToSkip, nLastBinsToSkip);
	    } else{
		dataMcChi2 =  TH1Chi2Test(hGenDYJets, (*hUnfs)[i], nFirstBinsToSkip, nLastBinsToSkip);
	    }
	    hchi2dataMC->SetBinContent(i + 1, dataMcChi2);
	    
	    double maxRes = -1.;
	    for(int j = nFirstBinsToSkip; j < nRes - nLastBinsToSkip; ++j){
		double x = pow((res[j]),2);
		hRes->SetBinContent(1 + j, x);
		if(x > maxRes) maxRes = x;
	    }
	    hRes->Write();
	    if(i>0) hResMax->SetBinContent(hResMax->GetXaxis()->FindBin(i), maxRes);
	    if(nIterResMax < 0  && maxRes < resMaxThr) nIterResMax = i;
	}
	hResMax->Write();

	TH1D *hgen = (TH1D*) respBis->Htruth();
	TString tmpName = "mcGen" + name;
	hgen->SetName(tmpName);
	hgen->SetTitle(tmpName);
	hgen->Write();
    
	TH1D *hfoldgen = foldUnfData(hgen, 0, respBis);
	tmpName = "mcGenFolded" + name;
	hfoldgen->SetName(tmpName);
	hfoldgen->SetTitle(tmpName);
	hfoldgen->Write();

	TH1D *hmes = (TH1D*) respBis->Hmeasured();
	tmpName = "mcReco" + name;
	hmes->SetName(tmpName);
	hmes->SetTitle(tmpName);
	hmes->Write();
    
	//hRecDataMinusFakes->Write("Unf" + name + "_000");
	hRecDataMinusFakes->Write("hBkgSubData");

	if(hRecDataMinusFakesOdd && hRecDataMinusFakesEven){
	    
	    //hchi2Xval = (TH1D*) hchi2->Clone("hchi2Xval");
	    //hchi2Xval->Reset();
      
	    RooUnfoldResponse *respBis2 = (RooUnfoldResponse*) resp->Clone();
	    TH1D *hRecDataMinusFakesBisOdd = (TH1D*) hRecDataMinusFakesOdd->Clone();
	    TH1D *hRecDataMinusFakesBisEven = (TH1D*) hRecDataMinusFakesEven->Clone();

	    std::auto_ptr<std::vector<TH1*> > hUnfs2 (new std::vector<TH1*>);
      
	    //Unfolds the data:
	    //	    unfoldWithErr(alg, respBis2, hRecDataMinusFakesBisOdd, nTestIterMax, smoothPrior, hUnfs2.get(), &covs, 2);
	    unfoldWithErr(alg, respBis2, hRecDataMinusFakesBisOdd, nTestIterMax, smoothPrior, hUnfs2.get());
      
	    for (unsigned i = 0; i < hUnfs2->size(); ++i) {
		(*hUnfs2)[i]->SetName(TString::Format("hUnfXval%d", i));
		if(i==0){
		    (*hUnfs2)[i]->SetTitle(TString::Format("Prior for %s %s", variable, name.Data()));
		} else {
		    (*hUnfs2)[i]->SetTitle(TString::Format("Unfolded %s %s %d iters, cross-validation", variable, name.Data(), i));
		}
		(*hUnfs2)[i]->Write();
		//TH1* hfoldUnfDataOdd = foldUnfData((*hUnfs2)[i], &(covs[i]), respBis);
		TH1* hfoldUnfDataOdd = foldUnfData((*hUnfs2)[i], 0, respBis);
		hfoldUnfDataOdd->SetName(TString::Format("dataOddfoldedBack%s_%d", name.Data(), i));
		hfoldUnfDataOdd->SetTitle(TString::Format("dataOddFoldedBack%s_%d", name.Data(), i));
		hfoldUnfDataOdd->Write();
	
		TH1* hResXval = (TH1*) hRecDataMinusFakes->Clone(TString::Format("hResXval%s_%d", name.Data(), i+1));
		hResXval->Reset();
		hResXval->SetTitle(TString::Format("%s %s residuals - iteration %d, cross validation",
						   variable, name.Data(), i+1));
		hResXval->GetYaxis()->SetTitle("Residuals");

		const int nRes = hfoldUnfDataOdd->GetNbinsX();
		Double_t resXval[nRes];
		double mychi2 = TH1Chi2Test(hRecDataMinusFakesBisEven,
					    i == 0 ? hRecDataMinusFakesBisOdd : hfoldUnfDataOdd,
					    nFirstBinsToSkip, nLastBinsToSkip, resXval);
		
		if(verbosity) std::cout << "Chi2/nbins of data / folded-unfolded distributions: " << mychi2 << "\n"; 

		hchi2Xval->SetBinContent(hchi2Xval->GetXaxis()->FindBin(i), mychi2);
	
		if (i > 0 && mychi2 < chi2XvalThr && nIterXval < 0) {
		    nIterXval = i;
		    //finalNIterXval = i;
		    std::cout << "Cross-validation chi^2 leads to " << nIter << " iterations with a final Chi2/ndf of: " << mychi2 << std::endl;
		}
		if(i > 0 && ((int)i >= minIter) && (mychi2 < chi2XvalMin)){
		    nIterXvalMin =  i;
		    chi2XvalMin = mychi2;
		}
	
		double maxRes = -1.;
		for(int j = nFirstBinsToSkip; j < nRes - nLastBinsToSkip; ++j){
		    double x = pow(resXval[j], 2);
		    hResXval->SetBinContent(1 + j, x);
		    if(x > maxRes) maxRes = x;
		}
		hResXval->Write();
		if(i > 0) hResMaxXval->SetBinContent(hResXval->GetXaxis()->FindBin(i), maxRes);
		if(nIterResMaxXval < 0 && maxRes < resMaxXvalThr){
		    nIterResMaxXval = i;
		}
	    } //next hUnfs2
	    hResMaxXval->Write();

	    TH1D *hgen = (TH1D*) respBis->Htruth();
	    TString tmpName = "mcGen" + name;
	    hgen->SetName(tmpName);
	    hgen->SetTitle(tmpName);
	    hgen->Write();
      
	    TH1D *hfoldgen = foldUnfData(hgen, 0, respBis);
	    tmpName = "mcGenFolded" + name;
	    hfoldgen->SetName(tmpName);
	    hfoldgen->SetTitle(tmpName);
	    hfoldgen->Write();
      
	    TH1D *hmes = (TH1D*) respBis->Hmeasured();
	    tmpName = "mcReco" + name;
	    hmes->SetName(tmpName);
	    hmes->SetTitle(tmpName);
	    hmes->Write();
	    hRecDataMinusFakesOdd->Write("UnfOdd" + name + "_0");
	    hRecDataMinusFakesEven->Write("UnfEven" + name + "_0");
	}//cross-validation
    
	std::vector<double> rms(nTestIterMax+1, 0);
	std::vector<double> chi2ToyMcErrs(nTestIterMax + 1, 0);
	std::vector<TH1*> hChi2s(nTestIterMax+1);
	std::vector<TProfile*> hResToys(nTestIterMax+1);
	std::vector<TProfile*> hResToysII(nTestIterMax+1);
	std::vector<double> rmsII(nTestIterMax+1, 0);
	std::vector<TH1*> hChi2IIs(nTestIterMax+1);
	std::vector<double> chi2ToyII(nTestIterMax+1, 0);


	for(int i = 0; i < nTestIterMax + 1; ++i){
	    hChi2s[i] = new TH1D(TString::Format("hChi2_%d", i),
				 TString::Format("#chi2/ndf distribution for %d iteration", i),
				 100, 0, 10);

	    TString n = TString::Format("hResToy%s_%d", name.Data(), i);
	    TString t = TString::Format("%s %s residuals - iteration %d", variable, name.Data(), i);
	    const TArrayD* bins = hRecDataMinusFakes->GetXaxis()->GetXbins();
	    if(bins){
		hResToys[i] = new TProfile(n, t, hRecDataMinusFakes->GetNbinsX(), bins->GetArray());
	    } else{
		hResToys[i] = new TProfile(n, t, hRecDataMinusFakes->GetNbinsX(),
					   hRecDataMinusFakes->GetXaxis()->GetBinLowEdge(1),
					   hRecDataMinusFakes->GetXaxis()->GetBinUpEdge(hRecDataMinusFakes->GetNbinsX()));
	    }
	    //	    hResToysII[i] = (TProfile*) hResToys[i]->Clone(n+"II");
	    n += "II";
	    if(bins){
		hResToysII[i] = new TProfile(n, t, hRecDataMinusFakes->GetNbinsX(), bins->GetArray(), "S");
	    } else{
		hResToysII[i] = new TProfile(n, t, hRecDataMinusFakes->GetNbinsX(),
					     hRecDataMinusFakes->GetXaxis()->GetBinLowEdge(1),
					     hRecDataMinusFakes->GetXaxis()->GetBinUpEdge(hRecDataMinusFakes->GetNbinsX()),
					     "S");
	    }

	    hChi2IIs[i] = new TH1D(TString::Format("hChi2II_%d", i),
				   TString::Format("#chi2/ndf distribution for %d iteration%s", i, (i>1?"s":"")),
				   100, 0, 10);
      
	}

	std::vector<double> chi2Toy = chi2FromToy(alg, smoothPrior, resp, hRecDataMinusFakes, nFirstBinsToSkip,
						  nLastBinsToSkip,
						  nTestIterMax, ntoys, &chi2ToyMcErrs, &rms,
						  &hChi2s, 0, &chi2ToyII, &rmsII,
						  &hChi2IIs, &hResToys, 0, &hResToysII, toyIIrecoResampling);

	TParameter<int>("ntoys", ntoys).Write();

	for(int iter = 0; iter <= nTestIterMax; ++iter){
	    double x = 0;
	    TH1* h = hResToys[iter];
	    for(int ibin = 1 + nFirstBinsToSkip; ibin <= h->GetNbinsX() - nLastBinsToSkip; ++ibin){
		double c = h->GetBinContent(ibin);
		if(c > x) x = c;
	    }
	    hResMaxToy->SetBinContent(hResMaxToy->GetXaxis()->FindBin(iter), x);
	    if(nIterResMaxToy < 0 && x < resMaxToyThr){
		nIterResMaxToy = iter;
	    }
	    if(iter > 0 && x < resMaxToyMin){
		resMaxToyMin = x;
		nIterResMaxToyMin = iter;
	    }
	    h->Write();
	}

	hResMaxToy->Write();
    
	for(int iter = 1; iter <= nTestIterMax; ++iter){
	    double x = hResMax->GetBinContent(hResMax->GetXaxis()->FindBin(iter));
	    double xtoy = hResMaxToy->GetBinContent(hResMaxToy->GetXaxis()->FindBin(iter));
	    if(nIterResMaxToyCmp < 0 && x < resMaxToyCmpThr1 && (xtoy - x) >  resMaxToyCmpThr2){
		nIterResMaxToyCmp = iter;
	    }
	}

	for(int iter = 1; iter <= nTestIterMax; ++iter){
	    double bias = 0;
	    double ebias = 0;
	    double chi2 = 0;
	    TH1* h = hResToysII[iter];
	    for(int ibin = 1 + nFirstBinsToSkip; ibin <= h->GetNbinsX() - nLastBinsToSkip; ++ibin){
		double c = fabs(h->GetBinContent(ibin));
		if(c > bias){
		    bias = c;
		    if(h->GetBinContent(ibin) >0) ebias = h->GetBinError(ibin) / sqrt(h->GetBinContent(ibin));
		    else ebias = 0;
		}
		c = pow(h->GetBinError(ibin),2) + pow(c,2);
		//c = h->GetBinError(ibin);
		if(c > chi2) chi2 = c;
	    }
	    hResMaxToyII->SetBinContent(hResMaxToyII->GetXaxis()->FindBin(iter), chi2);
	    hBiasMaxToyII->SetBinContent(hBiasMaxToyII->GetXaxis()->FindBin(iter), bias);
	    hBiasMaxToyII->SetBinError(hBiasMaxToyII->GetXaxis()->FindBin(iter), ebias);
	    if(nIterResMaxToyII < 0 && chi2 > resMaxToyIIThr){
		nIterResMaxToyII = iter;
	    }
	    if(chi2 < resMaxToyIIMin){
		nIterResMaxToyIIMin = iter;
		resMaxToyIIMin = chi2;
	    }
	    if(bias < biasToyIIMin){
		nIterBiasToyIIMin = iter;
		biasToyIIMin= bias;
	    }
	    
	    h->Write();
	}
	
	hResMaxToyII->Write();
	hBiasMaxToyII->Write();
	
	
	for(int iter = 0; iter < nTestIterMax + 1; ++iter){
      
	    if(iter> 0 && chi2Toy[iter] < chi2ToyThr && nIterToy < 0) {
		nIterToy = iter;
		std::cout << "Toy method lead to " << nIterToy << " iterations with a final Chi2/ndf of: " << chi2Toy[iter] << std::endl;
	    }

	    if(((int)iter >= minIter) && (chi2Toy[iter] < chi2ToyMin)){
		nIterToyMin =  iter;
		chi2ToyMin = chi2Toy[iter];
	    }

	    if(iter >1 /*iter > 0*/ && !nIterToyIIMinFound && chi2ToyII[iter] < chi2ToyIIMin){
		nIterToyIIMin =  iter;
		chi2ToyIIMin = chi2ToyII[iter];
	    }

	    if(nIterToyIIMin > 0 && chi2ToyIIMin < chi2ToyII[iter]) {//local minimum found
		nIterToyIIMinFound = true;
	    }
	    
	    if(iter > 0 && nIterToyII < 0 && chi2ToyII[iter] < chi2ToyIIThr){
		nIterToyII = iter;
	    }
	    
	    hChi2s[iter]->Write();
	    hChi2IIs[iter]->Write();
	    //if(iter>0){ //skip prior
	    hchi2Toy->SetBinContent(hchi2Toy->FindBin(iter), chi2Toy[iter]);
	    hchi2ToyMcErr->SetBinContent(hchi2ToyMcErr->FindBin(iter), chi2ToyMcErrs[iter]);
	    hchi2Toy->SetBinError(hchi2Toy->FindBin(iter), rms[iter]/ntoys);
	    hchi2ToyII->SetBinContent(hchi2ToyII->FindBin(iter), chi2ToyII[iter]);
	    hchi2ToyII->SetBinError(hchi2ToyII->FindBin(iter), rmsII[iter]/ntoys);
	    //}
	}
     
	if(nIterXval < 0){
	    std::cout << "Chi2/ndf from xval is always larger than 1. Using minimum as number of iterations.\n";
	    nIterXval = nIterXvalMin;
	}

	if(nIterToy < 0){
	    std::cout << "Chi2/ndf from toy is always larger than 1. Using minimum as number of iterations.\n";
	    nIterToy = nIterToyMin;
	}

	if(nIterToyII < 0){
	    std::cout << "Chi2/ndf of Toy II is always larger than 1. Using minimum as number of iterations\n";
	    nIterToyII = nIterToyIIMin;
	}
	
	if(xvalIter == 9){
	    chosenIter = nIterToyII;
	    chosenAlgoThr = -1;
	} else if(xvalIter == 8){
	    if(nIterResMaxToyMin < nTestIterMax) chosenIter = nIterResMaxToyMin;
	    else if(nIterResMax > 0 && nIterResMax < nTestIterMax) chosenIter = nIterResMax;
	    else chosenIter = nIter;
	    chosenAlgoThr = chi2Thr;
	} else if(xvalIter == 7){
	    chosenIter = nIterResMaxToyCmp;
	    chosenAlgoThr = -1;
	} else if(xvalIter == 6){
	    chosenIter = nFixedIters; //readNiters(lepSel.Data(), variable);
	} else if(xvalIter == 5){
	    chosenIter = nIterResMaxToy;
	    chosenAlgoThr = resMaxToyThr;
	} else if(xvalIter == 4){
	    chosenIter = nIterResMaxXval;
	    chosenAlgoThr = resMaxXvalThr;
	} else if(xvalIter == 3){
	    chosenIter = nIterResMax;
	    chosenAlgoThr = resMaxThr;
	} else if(xvalIter == 2){
	    chosenIter = nIterToy;
	    chosenAlgoThr = chi2ToyThr;
	} else if(xvalIter == 1){
	    chosenIter = nIterXval;
	    chosenAlgoThr = chi2XvalThr;
	} else if(xvalIter == 0){
	    chosenIter = nIter;
	    chosenAlgoThr = chi2Thr;
	} else{
	    std::cerr << "Invalid value (" << xvalIter << ") for parameter xvalIter.\n";
	    abort();
	}
	if(chosenIter < minIter){
	    chosenIter = minIter;
	}
	if(chosenIter > maxIter || chosenIter < 0){
	    chosenIter = maxIter;
	}

	const char* mode_names[] = {"chi2", "chi2 xval", "chi2 toy",
				    "max res", "max res xval", "max res toy", 
				    "fixed", "max res cmp", "min of max res toy",
				    "chi2 toy II min"};
	if((unsigned)xvalIter >= sizeof(mode_names)/sizeof(mode_names[0])){
	    std::cerr << "Bug found in " << __FILE__ << ":" << __LINE__ << "\n";
	    abort();
	}
	std::cout << "Bayes unfold, choice of number of iterations:\n"
	    "\tmin: " << minIter << "\tmax: " << maxIter
		  << "\tmode: " << mode_names[xvalIter]
		  << "\n\tnIter(std): " << nIter
		  << "\tnIter(xval): " << nIterXval
		  << "\tnIter(toy): " << nIterToy
		  << "\n\tnIter(res max): " << nIterResMax
		  << "\tnIter(res max, xval): " << nIterResMaxXval
		  << "\tnIter(res max, toy): " << nIterResMaxToy
		  << "\tnIter(chi2 toyII): " << nIterToyII
		  << "\tnIter(min of chi2 toyII): " << nIterToyIIMin
		  << "\n\tSelected value: " << chosenIter << "\n";
	if(outputFileName){
	    std::ofstream fniter(outputFileName);
	    fniter << chosenIter << "\n";
	    fniter << variable << "\t" << nIter << "\t" << nIterXval << "\t" << nIterToy
		   << "\t" << nIterResMax << "\t" << nIterResMaxXval << "\t" << nIterResMaxToy
		   << "\t" << nIterResMaxToyCmp << "\t" << nIterResMaxToyII << "\t" << chosenIter << "\n";
	}

	int chi2Ylog = cfg.getUnf(lepSel, variable, "unfChi2LogScale", 0);
	
	makeChi2Plot(hchi2, xvalIter == 0 ? -1 : nIter, chosenIter,
		     chi2Thr, chosenAlgoThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_chi2", "cChi2", chi2Ylog);

	if(hRecDataMinusFakesOdd && hRecDataMinusFakesEven){
	    hchi2Xval->Write();
	    makeChi2Plot(hchi2Xval, xvalIter == 1 ? -1 : nIterXval, chosenIter,
			 chi2XvalThr, chosenAlgoThr,
			 unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_chi2Xval", "cChi2Xval", chi2Ylog);
	}


	hchi2Toy->Write();
	hchi2ToyMcErr->Write();
	std::cout << nIterToy << ", " << chosenIter << ", " << 		 chi2ToyThr << ", " << chosenAlgoThr << "\n";
	makeChi2Plot(hchi2Toy, xvalIter == 2 ? -1 : nIterToy, chosenIter,
		     chi2ToyThr, chosenAlgoThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_chi2Toy", "cChi2Toy", chi2Ylog);


	makeChi2Plot(hResMax, xvalIter == 3 ? -1 : nIterResMax, chosenIter,
		     resMaxThr, chosenAlgoThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_hResMax", "cResMax", chi2Ylog);

	if(hResMaxXval){
	    makeChi2Plot(hResMaxXval, xvalIter == 4 ? -1 : nIterResMaxXval, chosenIter,
			 resMaxXvalThr, chosenAlgoThr,
			 unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_hResXval", "cResMaxXval", chi2Ylog);
	}

	if(xvalIter==6){
	    makeChi2Plot(hResMaxXval, -1, chosenIter,
			 -1,  -1,
			 unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_hFixed", "cFixed", chi2Ylog);
	}
	makeChi2Plot(hResMaxToy, xvalIter == 5 ? -1 : nIterResMaxToy, chosenIter,
		     resMaxToyThr, chosenAlgoThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_hResMaxToy", "cResMaxToy", chi2Ylog);

	std::vector<TH1*> altH(1);
	std::vector<std::string> labels(2);
	labels[0] = "Toy MC";
	labels[1] = "Same sample";
	altH[0] = hResMax;
	int narr1 = nIterResMaxToyCmp;
	if(xvalIter==7) narr1 = -1;
	if(xvalIter==8) narr1 = nIterResMaxToyMin;
	double thrline = -1;
	if(3 <= xvalIter && xvalIter <= 5) thrline = chosenAlgoThr; //the algo uses a thresold on this curve
	makeChi2Plot(hResMaxToy, narr1 , chosenIter, -1, thrline,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_hResMaxToyCmp", "cResMaxToyCmp",
		     chi2Ylog, altH, labels);
    
	hchi2ToyII->Write();
	makeChi2Plot(hchi2ToyII, nIterToyII, chosenIter,
		     chi2ToyIIThr, chosenAlgoThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_chi2ToyII", "cChi2ToyII", chi2Ylog);

	makeChi2Plot(hResMaxToyII, nIterResMaxToyIIMin, chosenIter,
		     0., chosenAlgoThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_resMaxToyII",
		     "cResMaxToyII", chi2Ylog);

	altH.resize(1);
	altH[0] = hResMaxToyII;
	labels.resize(altH.size()+1);
	labels[0] = "Hist. #chi2/ndf";
	labels[1] = "max(#chi2_{bin})";
	TH1* hchi2ToyIICpy = (TH1*) hchi2ToyII->Clone();
	hchi2ToyIICpy->SetTitle(TString::Format("Gen. #chi^{2} for %s", variable));
	hchi2ToyIICpy->GetYaxis()->SetTitle("#chi^{2}/ndf or max(#chi^{2}_{bin})");
				
	makeChi2Plot(hchi2ToyIICpy, nIterToyII, chosenIter,
		     chi2ToyIIThr, chosenAlgoThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_chi2AndResMaxToyII",
		     "cChi2AndResMaxToyII", chi2Ylog,
		     altH, labels, true);	
	
	altH.resize(2);
	altH[0] = hchi2;
	altH[1] = hchi2Xval;
	labels.resize(altH.size() + 1);
	labels[0] = "Toy MC";
	labels[1] = "Same sample";
	labels[2] = "Two samples";
	makeChi2Plot(hchi2Toy, xvalIter == 2 ? -1 : nIterToy, chosenIter,
		     chi2ToyThr, chosenAlgoThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_chi2All", "cChi2All",
		     chi2Ylog, altH, labels);

	altH.resize(1);
	altH[0] = hchi2Toy;
	labels.resize(altH.size() + 1);
	labels[0] = "Same sample";
	labels[1] = "Toy MC";
	makeChi2Plot(hchi2, nIter, chosenIter,
		     chi2Thr, chi2Thr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_chi2AndToyChi2", "cChi2AndToyChi2",
		     chi2Ylog, altH, labels);


	
	altH.resize(0);
	labels.resize(altH.size() + 1);
	labels[0] = "Data vs MC";
	makeChi2Plot(hchi2dataMC, -1, chosenIter,
		     chi2dataMCThr, chi2ToyThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_chi2DataMC", "cChi2DataMC",
		     chi2Ylog, altH, labels);
	
	makeChi2Plot(hBiasMaxToyII, nIterBiasToyIIMin, chosenIter, -1, chosenAlgoThr,
		     unfoldCheckDir, lepSel + "_" + variable + "_" + name + "_" + algo + "_bias", "bias", chi2Ylog);
    } //end of test of number of iterations

    if(binByBin_unfold){
	std::unique_ptr<RooUnfold> RObjectForDataBinByBin(RooUnfold::New(RooUnfold::kBinByBin, resp, hRecDataMinusFakes));
	RObjectForDataBinByBin->SetVerbose(verbosity);
	TH1D *hUnfDataBinByBin = (TH1D*) RObjectForDataBinByBin->Hreco(RooUnfold::kCovariance);
	hUnfDataBinByBin->SetName("UnfDataBinByBin" + name);
	hUnfDataBinByBin->Write();
    }
    
    if(svd_unfold){
	for (int iter(1); iter <= nTestIterMax; iter++) {
	    std::unique_ptr<RooUnfold> RObjectForDataSVD(RooUnfold::New(RooUnfold::kSVD, resp, hRecDataMinusFakes, iter));
	RObjectForDataSVD->SetVerbose(verbosity);
	TH1D *hUnfDataSVD = (TH1D*) RObjectForDataSVD->Hreco(RooUnfold::kCovariance);
	    hUnfDataSVD->SetName("UnfDataSVD_" + TString::Format("%d", iter) + "_" + name);
	hUnfDataSVD->Write();
	}
    }

    if(tsvd_unfold){
	std::cout << "----- RUN TSVD --- \n";
	std::unique_ptr<TSVDUnfold> unfoldTSVD(new TSVDUnfold(hRecDataMinusFakes, (TH1D*)resp->Htruth(), (TH1D*)resp->Hmeasured(), (TH2D*)resp->Hresponse()));
	TH1D *unfresult = (TH1D*) unfoldTSVD->Unfold(1);
	TH1D *hmodDOriginal = (TH1D*) unfoldTSVD->GetD();
	TH1D *hSV       = (TH1D*) unfoldTSVD->GetSV();
      
	TH1D *hmodD = new TH1D("hmodD", "hmodD", nTestIterMax, 0.5, nTestIterMax+0.5);
	for (int iter(0); iter <= nTestIterMax+1; iter++) {
	    hmodD->SetBinContent(iter, hmodDOriginal->GetBinContent(iter));
	}
	hmodD->SetTitle(hmodDOriginal->GetTitle() + TString(" for ") + TString(hRecDataMinusFakes->GetTitle()));
	hmodD->GetXaxis()->SetNdivisions(nTestIterMax, 0, 0);
	hmodD->GetXaxis()->SetLabelSize(0.03);
	hmodD->GetYaxis()->SetTitle("|d_{i}|");
	hmodD->GetYaxis()->SetTitleOffset(1.40);
	hmodD->GetXaxis()->SetTitle("regularization parameter of the SVD method");
	hmodD->GetXaxis()->CenterTitle();
	TArrow *arrowSvd = new TArrow(svdKterm, 20, svdKterm, 1.1*hmodD->GetBinContent(svdKterm), 0.02, "|>");
	arrowSvd->SetLineColor(kRed);
	arrowSvd->SetFillColor(kRed);
	arrowSvd->SetLineWidth(2);
      
	hmodD->SetName("modD" + name);
	TCanvas *chmodD = new TCanvas("chmodD", "chmodD", 600, 600);
	chmodD->cd();
	chmodD->SetGrid();
	chmodD->SetLogy(logy ? kTRUE : kFALSE);
	hmodD->DrawCopy();
	arrowSvd->Draw();
	chmodD->Write();
	hmodD->Write();
	std::cout << "Wrting tsvd resuslts in " << gDirectory->GetFile()->GetName() << "\n";
	unfresult->Write();
	hSV->Write();
    }



    if(invert_unfold){
	std::cout << "Running unfolding with matrix inversion\n";
	std::unique_ptr<RooUnfold> RObjectForDataInvert(RooUnfold::New(RooUnfold::kInvert, resp, hRecDataMinusFakes));
	RObjectForDataInvert->SetVerbose(verbosity);
	TH1D *hUnfDataInvert = (TH1D*) RObjectForDataInvert->Hreco(RooUnfold::kCovariance);
	hUnfDataInvert->SetName("UnfDataInvert" + name);
	hUnfDataInvert->Write();
    }

    f->Close();
    
    std::cout << "\n---------------------------------------------------------------------------------------------------------------\n-" << std::endl;

    //--- Unfold data minus background ---
    bool useFlatPrior = cfg.getUnf(lepSel, variable, "useFlatPrior", false);
    bool unfSmoothPrior = cfg.getUnf(lepSel, variable, "unfSmoothPrior", false);

    RooUnfold *RObjectForData = RooUnfold::New(alg, resp, hRecDataMinusFakes, chosenIter);
    if(alg==RooUnfold::kBayes) ((RooUnfoldBayes*) RObjectForData)->SetSmoothing(unfSmoothPrior);
    
    RObjectForData->SetVerbose(verbosity);
    RObjectForData->UseFlatPrior(useFlatPrior); //Added on Feb 1, 17
    RObjectForData->IncludeSystematics(0); // new version of RooUnfold: will compute Cov based on Data Statistics only
    hUnfData = (TH1D*) RObjectForData->Hreco(RooUnfold::kCovariance);

    hUnfData->SetName("UnfData" + name);
    hUnfData->SetTitle(hRecDataMinusFakes->GetTitle());

    if (algo == "Bayes") {
	//--- get covariance from statistics on Data ---
	hUnfDataStatCov = M2H(RObjectForData->Ereco(RooUnfold::kCovariance)); // new version of RooUnfold   
	hUnfDataStatCov->SetName("UnfDataStatCov" + name);
	hUnfDataStatCov->SetTitle(hRecDataMinusFakes->GetTitle());

	//--- get covariance from MC stat ---
	RooUnfold *RObjectForMC = RooUnfold::New(alg, resp, hRecDataMinusFakes, chosenIter);
	RObjectForMC->SetVerbose(verbosity);
	RObjectForMC->UseFlatPrior(useFlatPrior); //Added on Feb 9, 17
	if(alg==RooUnfold::kBayes) ((RooUnfoldBayes*) RObjectForMC)->SetSmoothing(smoothPrior);
	RObjectForMC->IncludeSystematics(2); // new version of RooUnfold: will compute Cov based on MC Statistics only
	hUnfMCStatCov = M2H(RObjectForMC->Ereco(RooUnfold::kCovariance)); // new version of RooUnfold
	hUnfMCStatCov->SetName("UnfMCStatCov" + name);
	hUnfMCStatCov->SetTitle(hRecDataMinusFakes->GetTitle());


	TH1* hUnfStdUnc;
	TGraphAsymmErrors* gr = statFromToy(RooUnfold::kBayes, smoothPrior, resp,
					    hRecDataMinusFakes,
					    chosenIter, ntoys, &hUnfStdUnc);
	TFile *f = new TFile(unfoldCheckDir + "/" + lepSel + "_" + variable + "_" + name + ".root", "UPDATE");
	gr->Write(gr->GetName());
	hUnfStdUnc->SetDirectory(f);
	hUnfStdUnc->Write("hUnfStdStat");
	f->Close();
    }

    //--- divide by luminosity ---
    hUnfData->Scale(1./integratedLumi);
    if(hUnfDataStatCov) hUnfDataStatCov->Scale(1./(integratedLumi*integratedLumi));
    if(hUnfDataStatCov) hUnfMCStatCov->Scale(1./(integratedLumi*integratedLumi));
    if ("LumiUp" == name) {
        hUnfData->Scale(1./(1+lumiUnc));
        if(hUnfDataStatCov) hUnfDataStatCov->Scale(1./((1+lumiUnc)*(1+lumiUnc)));
        if(hUnfMCStatCov) hUnfMCStatCov->Scale(1./((1+lumiUnc)*(1+lumiUnc)));
    }
    else if ("LumiDown" == name) {
        hUnfData->Scale(1./(1-lumiUnc));
        if(hUnfDataStatCov) hUnfDataStatCov->Scale(1./((1-lumiUnc)*(1-lumiUnc)));
        if(hUnfMCStatCov) hUnfMCStatCov->Scale(1./((1-lumiUnc)*(1-lumiUnc)));
    }


    //--- divide by bin width to get cross section ---
    int nBins = hUnfData->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
	double binWidth = hUnfData->GetBinWidth(i);
	hUnfData->SetBinContent(i, hUnfData->GetBinContent(i)*1./binWidth);
	hUnfData->SetBinError(i, hUnfData->GetBinError(i)*1./binWidth);
	for (int j = 1; j <= nBins; ++j) {
	    if(hUnfDataStatCov){
		hUnfDataStatCov->SetBinContent(i, j, hUnfDataStatCov->GetBinContent(i, j)*1./(binWidth*binWidth));
		hUnfDataStatCov->SetBinError(i, j, hUnfDataStatCov->GetBinError(i, j)*1./(binWidth*binWidth));
	    }
	    if(hUnfMCStatCov){
		hUnfMCStatCov->SetBinContent(i, j, hUnfMCStatCov->GetBinContent(i, j)*1./(binWidth*binWidth));
		hUnfMCStatCov->SetBinError(i, j, hUnfMCStatCov->GetBinError(i, j)*1./(binWidth*binWidth));
	    }
	}
    }
    return chosenIter;
}


//--- This is needed because we use overflow and the matrices have 2 additional bins ---
TH2D* M2H(const TMatrixD m) 
{
    int nBinsY = m.GetNrows();
    int nBinsX = m.GetNrows();

    TH2D *h = new TH2D(m.GetName(), m.GetName(), nBinsX-2, 0, nBinsX-2, nBinsY-2, 0, nBinsY-2);
    for (int i = 0; i < nBinsX; ++i) {
	for (int j = 0; j < nBinsY; ++j) {
	    h->SetBinContent(i, j, m(i, j));
	}
    }

    return h;
}

TH2D* makeCovFromUpAndDown(const TH1D* hUnfDataCentral, const TH1D* hUnfDataUp,
			   const TH1D* hUnfDataDown, TString name, bool noNull)
{
    int nBins = hUnfDataCentral->GetNbinsX();
    TH2D* h = new TH2D(name, name, nBins, 0, nBins, nBins, 0, nBins);

    double scale = 1.;
    //if no down variations we substitute to the diff. between up and down, twice the
    //diff. between up and central:
    if(hUnfDataDown == 0){
	hUnfDataDown = hUnfDataCentral;
	scale  = 2.;
    }
////    if(name == "CovLumi"){
////	std::cerr << ">>>> " << hUnfDataCentral->GetBinContent(4)
////		  << "\t"  << hUnfDataUp->GetBinContent(4)
////		  << "\t"  << hUnfDataDown->GetBinContent(4)
////		  << "\t"  << 0.5*(hUnfDataUp->GetBinContent(4)-hUnfDataDown->GetBinContent(4))/hUnfDataCentral->GetBinContent(4) <<"\n";
////    }
    
    //1. if noNull is set to true, in case of crossing of up- and down- variation
    //   distribution,  we take for the uncertainty of the bin where the crossing
    //   occur, the average of uncertainties of the two neighbour bins
    //2. For correlation, we assume full correlation and anticorrelation.
    //   correlation sign is determined by comparing the directions of the
    //   variations.
    for (int i = 1; i <= nBins; ++i) {
	double sigma_i = 0.5*scale*fabs(hUnfDataUp->GetBinContent(i) - hUnfDataDown->GetBinContent(i));
	//variation directions used later to determine the sign of the correlation coefficient:
	int sign_i = (hUnfDataUp->GetBinContent(i) - hUnfDataDown->GetBinContent(i) < 0) ? -1 : 1;
	if (noNull && i > 1 && i < nBins) {
	    if ((hUnfDataUp->GetBinContent(i-1) - hUnfDataDown->GetBinContent(i-1))
		* (hUnfDataUp->GetBinContent(i+1) - hUnfDataDown->GetBinContent(i+1)) < 0) {
		sigma_i =   0.5*scale*(  0.5*fabs(hUnfDataUp->GetBinContent(i-1) - hUnfDataDown->GetBinContent(i-1))
			               + 0.5*fabs(hUnfDataUp->GetBinContent(i+1) - hUnfDataDown->GetBinContent(i+1)));
		//		if (name.Index("Sherpa") >= 0) sigma_i *= 2;
	    }
	}
	for (int j = 1; j <= nBins; ++j) {
	    double sigma_j = 0.5*scale*fabs(hUnfDataUp->GetBinContent(j) - hUnfDataDown->GetBinContent(j));
	    //variation directions used together with sign_i to determine the sign of the correlation coefficient:
	    int sign_j = (hUnfDataUp->GetBinContent(j) - hUnfDataDown->GetBinContent(j) < 0) ? -1 : 1;
	    if (noNull && j > 1 && j < nBins) {
		if ((hUnfDataUp->GetBinContent(j-1) - hUnfDataDown->GetBinContent(j-1))
		    * (hUnfDataUp->GetBinContent(j+1) - hUnfDataDown->GetBinContent(j+1)) < 0) {
		    sigma_j = 0.5*scale*(  0.5*fabs(hUnfDataUp->GetBinContent(j-1) - hUnfDataDown->GetBinContent(j-1))
					   + 0.5*fabs(hUnfDataUp->GetBinContent(j+1) - hUnfDataDown->GetBinContent(j+1)));
		    //		    if (name.Index("Sherpa") >= 0) sigma_j *= 2;
		}
	    }
	    double correlation = (i == j) ? 1 : 1;
	    h->SetBinContent(i, j, correlation * sign_i * sign_j * sigma_i * sigma_j);
	}
    }
    return h;
}


TH1D* foldUnfData(const TH1 *hUnfData, const TMatrixD* cov, const RooUnfoldResponse *hresp)
{
    TH1D* hfoldUnfData = (TH1D*) hUnfData->Clone();
    std::unique_ptr<RooUnfoldResponse> resp((RooUnfoldResponse*) hresp->Clone());
    TH2D *hres = (TH2D*) resp->Hresponse();
    TH1D *hgen = (TH1D*) resp->Htruth();


    int nBins = hfoldUnfData->GetNbinsX();

    for (int i = 0; i <= nBins; ++i) {
	double totGen = hgen->GetBinContent(i);
	for (int j = 0; j <= nBins; ++j) {
	    if (totGen != 0.0) {
		hres->SetBinContent(j, i, hres->GetBinContent(j, i)/totGen);
		hres->SetBinError(j, i, hres->GetBinError(j, i)/totGen);
		hres->SetBinError(j, i, 0.0);
	    }
	    else {
		hres->SetBinContent(i, j, 0);
		hres->SetBinError(i, j, hres->GetBinError(i, j));
		hres->SetBinError(i, j, 0.0);
	    }
	}
    }

    for (int i = 0; i <= nBins; ++i) {
	double sum = 0.0;
	double error2 = 0.0;
	for (int j = 0; j <= nBins; ++j) {
	    sum += hres->GetBinContent(i, j) * hUnfData->GetBinContent(j);
	    if(cov){
		for(int k = 0 ; k <= nBins; ++k){
		    for(int l = 0 ; l <= nBins; ++l){
			error2 += hres->GetBinContent(i, k) * hres->GetBinContent(i, l) * (*cov)(l, k);
		    }
		}
	    } else{
		error2 += hres->GetBinContent(i, j) * pow(hUnfData->GetBinError(j), 2);
	    }
	}
	hfoldUnfData->SetBinContent(i, sum);
	hfoldUnfData->SetBinError(i, sqrt(error2));
    }

    //    bool refoldForcePoisson = cfg.getB("refoldForcePoisson", true);
    //    //Forces uncertainty to be poissonian:
    //    if(refoldForcePoisson){
    //	//	std::cout << "Forcing poissonian errors for folded-back histogram.\n";
    //	for(int i = 0; i <= hfoldUnfData->GetNbinsX(); ++i){
    //	    hfoldUnfData->SetBinError(i, sqrt(hfoldUnfData->GetBinContent(i)));
    //	}
    //    }
    
    return hfoldUnfData;
}

void test()
{
    TH1D *hreco = new TH1D("hreco", "hreco", 2, 0, 2);
    TH1D *hgen = new TH1D("hgen", "hgen", 2, 0, 2);
    TH2D *hresp = new TH2D("hresp", "hresp", 2, 0, 2, 2, 0, 2);

    hreco->Sumw2();
    hgen->Sumw2();
    hresp->Sumw2();

    hreco->Fill(0.5, 3.5);
    hreco->Fill(1.5, 5);

    hgen->Fill(0.5, 4);
    hgen->Fill(1.5, 6);

    hresp->Fill(0.5, 0.5, 3);
    hresp->Fill(0.5,1.5,0.5);
    hresp->Fill(1.5,0.5,0.5);
    hresp->Fill(1.5,1.5,4.5);

    RooUnfoldResponse *resp = new RooUnfoldResponse(NULL, hgen, hresp);

    RooUnfold *RObjectForData = RooUnfold::New(RooUnfold::kBayes, resp, hreco, 3);
    TH1D *hUnfData = (TH1D*) RObjectForData->Hreco(RooUnfold::kCovariance);
    TMatrixD  cov = RObjectForData->Ereco(RooUnfold::kCovariance);
    
    TH1D *hfoldUnfData = foldUnfData(hUnfData, &cov, resp);

    std::cout << "reco: " << hreco->GetBinContent(1) << "   " << hreco->GetBinContent(2) << std::endl; 
    std::cout << "gen : " << hgen->GetBinContent(1) << "   " << hgen->GetBinContent(2) << std::endl; 
    std::cout << "unf : " << hUnfData->GetBinContent(1) << "   " << hUnfData->GetBinContent(2) << std::endl; 
    std::cout << "fol : " << hfoldUnfData->GetBinContent(1) << "   " << hfoldUnfData->GetBinContent(2) << std::endl; 
}

double MyChi2Test(const TH1 *h1, const TH1 *h2, int nFirstBinsToSkip, int nLastBinsToSkip, 
		  Double_t* res, const TH1* herr, bool poisErr)
{
    int nbins = h1->GetNbinsX();
    double chi2 = 0;
    if(herr==0) herr = h2;
    int ndof = 0;
    for (int i = 1; i <= nbins; ++i) {
	if(i < nFirstBinsToSkip + 1 || i > nbins - nLastBinsToSkip){
	    if(res) res[i-1] = 0;
	    continue;
	}
	double n1 = h1->GetBinContent(i);
	double n2 = h2->GetBinContent(i);
	double s = poisErr ? sqrt(herr->GetBinContent(i)) : herr->GetBinError(i);
	if(n1 > 0 && n2 > 0){
	    double r = (n2-n1) / s;
	    if(res) res[i-1] = r;
	    chi2 += r*r;
	    ++ndof;
	} else{
	    if(res) res[i-1] = 0.;
	}
    }
    return chi2 / ndof;
}

double TH1Chi2Test(const TH1 *h1, const TH1 *h2, int nFirstBinsToSkip, int nLastBinsToSkip, Double_t* res)
{
    TH1 *h1Copy = (TH1*) h1->Clone();
    TH1 *h2Copy = (TH1*) h2->Clone();
    
    //note: when skipping a first (resp. last) bin we also skip the underflow (resp. oveflow) bin.
    for (int i = 0; i <= nFirstBinsToSkip; ++i) {
	h1Copy->SetBinContent(i, 0.0);
	h2Copy->SetBinContent(i, 0.0);
    }

    for (int i = h1->GetNbinsX() - nLastBinsToSkip + 1; i <= h1->GetNbinsX() + 1; ++i) {
	h1Copy->SetBinContent(i, 0.0);
	h2Copy->SetBinContent(i, 0.0);
    }


    //x =h1Copy->Chi2Test(h2Copy, "WW,P,CHI2/NDF", res);
    //x =h1Copy->Chi2Test(h2Copy, "P,CHI2/NDF", res);
    //x =h1Copy->Chi2Test(h2Copy, "P", res);
    double x = h1Copy->Chi2Test(h2Copy, "CHI2/NDF", res);
    //x =h1Copy->Chi2Test(h2Copy, "", res);

    h1Copy->SetDirectory(0);
    h2Copy->SetDirectory(0);
    delete h1Copy;
    delete h2Copy;
    return x;
}

/** Correct data histogram for signal purity
 * @param hRecData histogram to correct
 * @param hPurity histogram containing the signal purities
 */
void RemoveFakes(TH1* hRecData, TH1* hFakes, TH1* hPurity){
    int fakeMethod = cfg.getI("fakeMethod");
    static bool emitMess = true;
    if(emitMess){
	std::cout << "Fake subtraction method id: " << fakeMethod << "\n";
	emitMess = false;
    }
    switch(fakeMethod){
    case 0: //global rescale to data
	hRecData->Add(hFakes, -1);
	return;
    case 1: //bin-by-bin purity method
	for(int ibin = 0; ibin <= hRecData->GetNbinsX(); ++ibin){
	    double c = hRecData->GetBinContent(ibin);
	    double e = hRecData->GetBinError(ibin);
	    double p = hPurity->GetBinContent(ibin);
	    //	if(ibin<4) std::cout <<  __FILE__ << __LINE__ << ": p(" << ibin
	    //		     << ") = " << p << "\n";
	    hRecData->SetBinContent(ibin, c * p);
	    hRecData->SetBinError(ibin, e * p);
	}
	return;
    default:
	std::cerr << __FILE__ <<  ":" << __LINE__ << ". Value " << fakeMethod
		  << " for parameter fakeMethod " << " is not supported. Exists.\n";
	exit(1);
    }
}

double pValueToNormChi2(double alpha, int n){
    double x =0;
    double dx = 0.01;
    if(alpha < 0 || alpha > 0.9999) return -1.;
    while(TMath::Prob(x, n) > alpha) x += dx;
    double x_corr =  0;
    if(x > dx) x_corr =  (TMath::Prob(x-dx, n) - alpha) /  (TMath::Prob(x, n) - TMath::Prob(x-0.01, n)) * 0.01;
    return (x + x_corr) / n;
}

TH1* unfold(RooUnfold::Algorithm algo, const RooUnfoldResponse* resp,
	    const TH1* hRecDataMinusFakes, int niters, bool smoothPrior, std::vector<TH1*>* hUnfs,
	    int uncMode){
    std::unique_ptr<RooUnfold> rooUnfold(RooUnfold::New(algo, resp, hRecDataMinusFakes, niters));
    if(algo==RooUnfold::kBayes) ((RooUnfoldBayes*) rooUnfold.get())->SetSmoothing(smoothPrior);
    bool verbosity = cfg.getB("unfoldingVerbosity");
    bool useFlatPrior = cfg.getB("useFlatPrior");
    rooUnfold->SetVerbose(verbosity);
    rooUnfold->UseFlatPrior(useFlatPrior);
    rooUnfold->IncludeSystematics(uncMode); // new version of RooUnfold: will compute Cov based on Data Statistics only
    //The following Hreco call will trigger the data unfolding
    TH1* r =  rooUnfold->Hreco(RooUnfold::kCovariance, hUnfs);
    return r;
}

TH1* unfoldWithErr(RooUnfold::Algorithm algo, const RooUnfoldResponse* resp,
		   const TH1* hRecDataMinusFakes, int niters, bool smoothPrior, std::vector<TH1*>* hUnfs,
		   std::vector<TMatrixD>* cov, int uncMode){
    bool verbosity = cfg.getB("unfoldingVerbosity");
    bool useFlatPrior = cfg.getB("useFlatPrior");
    if(hUnfs){
	hUnfs->clear();
	//prior:
	TH1D* hPrior = (TH1D*) resp->Htruth()->Clone(TString::Format("%s_prior", hRecDataMinusFakes->GetName()));
	hPrior->SetTitle (TString::Format("%s prior", hRecDataMinusFakes->GetTitle()));
	if(useFlatPrior){
	    hPrior->Reset();
	    const int nbins = hPrior->GetNbinsX();
	    for(int ibin = 1; ibin <= nbins; ++ibin){
		hPrior->SetBinContent(ibin, 1.);
		hPrior->SetBinError(ibin, 0);
	    }
	}
	double denom = hPrior->Integral();
	if(denom){
	    hPrior->Scale(hRecDataMinusFakes->Integral() / denom);
	}
	hUnfs->push_back(hPrior);
    }
    if(cov){
	cov->clear();
	int nbins = resp->Htruth()->GetNbinsX();
	cov->push_back(TMatrixD(TMatrixD::kZero, TMatrixD(nbins + 2, nbins + 2)));
    }
  
    TH1* r = 0;
    for(int i = 1; i <= niters; ++i){
	if(r && !hUnfs){
	    r->SetDirectory(0);
	    delete r;
	}
	std::unique_ptr<RooUnfold> rooUnfold(RooUnfold::New(algo, resp, hRecDataMinusFakes, i));
	rooUnfold->SetVerbose(verbosity);
	rooUnfold->UseFlatPrior(useFlatPrior);
	if(algo==RooUnfold::kBayes) ((RooUnfoldBayes*) rooUnfold.get())->SetSmoothing(smoothPrior);
	rooUnfold->IncludeSystematics(uncMode); // new version of RooUnfold: will compute Cov based on Data Statistics only
	//The following Hreco call will trigger the data unfolding
	//    r =  rooUnfold->Hreco(RooUnfold::kCovariance, 0);
	//rooUnfold->SetNToys(1000);
	r =  rooUnfold->Hreco(RooUnfold::kCovariance, 0);
	if(hUnfs) hUnfs->push_back(r);
	//	timeval t0, t1;
	//gettimeofday(&t0,0);
	if(cov) cov->push_back(rooUnfold->Ereco(RooUnfold::kCovariance));
	//	gettimeofday(&t1, 0);
	//	std::cout << __FILE__ << ":" << __LINE__ << ": execution of RooUnfold::Ereco took "
	//		  << ((t1.tv_sec - t0.tv_sec) * 1e3 + (t1.tv_usec - t0.tv_usec) * 1e-3)
	//		  << " ms.\n";
    }
    return r;
}


//std::vector<double> chi2FromToy(RooUnfold::Algorithm algo, const RooUnfoldResponse* resp,
//				const TH1* hmeas, int nFirstBinsToSkip, int niters, int ntoys,
//				std::vector<double>* rms, std::vector<TH1*>* hChi2){
//  std::vector<double> acc(niters+1, 0);
//  std::vector<double> acc2(niters+1, 0);
//  if(ntoys==0) ntoys = 1;
//  for(int itoy = 0; itoy < ntoys; ++itoy){
//    std::auto_ptr<std::vector<TH1*> > hUnfs (new std::vector<TH1*>);
//    TH1* h = resample(hmeas);
//    unfold(algo, resp, h, niters, hUnfs.get());
//    
//    for(int iter = 0; iter <= niters; ++iter){
//      TH1* hrefold = foldUnfData(hUnfs->at(iter), resp);
//      double chi2 = MyChi2Test(hmeas, hrefold, nFirstBinsToSkip);
//      acc[iter] += chi2;
//      acc2[iter] += chi2*chi2;
//      if(hChi2) hChi2->at(iter)->Fill(chi2);
//    }
//    delete h;
//  }
//  std::vector<double> mean(niters+1);
//  for(int i = 0; i < niters+1; ++i){
//    mean[i] = acc[i] / ntoys;
//    if(rms) rms->at(i) = 1. / (ntoys-1) * sqrt(acc2[i] - acc[i]*acc[i] / ntoys);
//  }
//  return mean;
//}


std::vector<double> 
chi2FromToy(RooUnfold::Algorithm algo, bool smoothPrior, const RooUnfoldResponse* resp,
	    const TH1* hmeas, int nFirstBinsToSkip, int nLastBinsToSkip, 
	    int niters, int ntoys, std::vector<double>* chi2mcErrs,
	    std::vector<double>* rms, std::vector<TH1*>* hChi2, std::vector<TH1*>* hChi2mcErr,
	    std::vector<double>* meanII, std::vector<double>* rmsII, 
	    std::vector<TH1*>* hChi2II, std::vector<TProfile*>* hRes, std::vector<TProfile*>* hResMcErr,
	    std::vector<TProfile*>* hResII, bool recoResampling){
    
    std::auto_ptr<std::vector<TH1*> > hRefs (new std::vector<TH1*>);

    TH1* hmeas_ = (TH1*) hmeas->Clone("hmeas_");
    hmeas_->SetDirectory(0);

    //Remove bins with less than 10 events for the chi2 calculation:
    //for(int i = 0; i <= hmeas_->GetNbinsX(); ++i){
    //  if(hmeas_->GetBinContent(i) < 10){
    //    hmeas_->SetBinContent(i, 0);
    //    hmeas_->SetBinError(i, 0);
    //  }
    //}

    unfoldWithErr(algo, resp, hmeas_, niters, smoothPrior, hRefs.get(), 0, 0);

    TH1* herr = (TH1*) hmeas_->Clone("herr");
    herr->SetDirectory(0);
    TH1* herrMcErr = (TH1*) hmeas_->Clone("herrMcErr");
    herrMcErr->SetDirectory(0);

    std::vector<double> acc(niters+1, 0);
    std::vector<double> acc2(niters+1, 0);
    std::vector<double> accMcErr(niters+1, 0);
    std::vector<double> accMcErr2(niters+1, 0);
    std::vector<double> accII(niters+1, 0);
    std::vector<double> acc2II(niters+1, 0);


    //Ref. folded histogram used for toy II
    std::vector<TH1*> hFolds (niters + 1);
    for(int iter = 0; iter <= niters; ++iter){
	hFolds[iter] = foldUnfData(hRefs->at(iter), 0, resp);
    }
    
    if(ntoys==0) ntoys = 1;
    for(int itoy = 0; itoy < ntoys; ++itoy){
	std::auto_ptr<std::vector<TH1*> > hUnfs (new std::vector<TH1*>);
	std::vector<TMatrixD> covs;
	TH1* h = resample(hmeas_);
	unfoldWithErr(algo, resp, h, niters, smoothPrior, hUnfs.get(), &covs, 2);

	const int nRes = hmeas_->GetNbinsX();
	Double_t res[nRes];
	Double_t resMcErr[nRes];

	for(int iter = 0; iter <= niters; ++iter){
	    TH1* hrefold;
	    if(iter==0) hrefold = h;
	    else hrefold = foldUnfData(hUnfs->at(iter), &(covs[iter]), resp);

	    //      for(int i = 0; i <= hrefold->GetNbinsX(); ++i){
	    //	  hrefold->SetBinError(i, sqrt(hrefold->GetBinContent(i));
	    //      }

	    for(int i = 0; i <= herr->GetNbinsX(); ++i){
		double err2 = hmeas_->GetBinContent(i);
		herr->SetBinContent(i, hmeas_->GetBinContent(i));
		herr->SetBinError(i, sqrt(err2));
		herrMcErr->SetBinContent(i, hmeas_->GetBinContent(i));
		if(iter > 0) err2 += pow(hrefold->GetBinError(i), 2);
		herrMcErr->SetBinError(i, sqrt(err2));
	    }
	
	    double chi2 = MyChi2Test(hmeas_, hrefold, nFirstBinsToSkip, nLastBinsToSkip, res, herr, false);
	    double chi2mcErr = MyChi2Test(hmeas_, hrefold, nFirstBinsToSkip, nLastBinsToSkip, resMcErr, herrMcErr, false);
	    acc[iter] += chi2;
	    acc2[iter] += chi2*chi2;
	    accMcErr[iter] += chi2mcErr;
	    accMcErr2[iter] += chi2mcErr*chi2mcErr;
	    if(hChi2) hChi2->at(iter)->Fill(chi2);
	    if(hChi2mcErr) hChi2mcErr->at(iter)->Fill(chi2mcErr);
	    for(int iRes = 0; iRes < nRes; ++iRes){
		double x = pow(res[iRes], 2);
		////	  std::cout << "----> iter, hmeas->GetXaxis()->GetFirst() + iRes, x "
		////		    << iter << ", "
		////		    << hmeas->GetXaxis()->GetFirst() + iRes << ", "
		////		    << x
		////		    << "\n";
		hRes->at(iter)->Fill(hRes->at(iter)->GetBinCenter(1 + iRes), x);
		if(hResMcErr) hResMcErr->at(iter)->Fill(hResMcErr->at(iter)->GetBinCenter(1 + iRes),
							pow(resMcErr[iRes], 2));
	    }

	    //toy II: histogram is resampled in the unsmeared space to perform a toy MC
	    //of the folding.
	    if(iter > 0){
#               if 0
		TH1* newFold = resample(hFolds[iter]);
		TH1* newUnf2 = unfold(algo, resp, newFold, iter, smoothPrior);
		double chi2II = MyChi2Test(hRefs->at(iter), newUnf2, nFirstBinsToSkip, 
					   nLastBinsToSkip, res, newUnf2, true);
#               else
		TH1* newUnf = resample(hRefs->at(iter));
		TH1* newFold = foldUnfData(newUnf, 0, resp);
		if(recoResampling) newFold = resample(newFold);
		TH1* newUnf2 = unfold(algo, resp, newFold, iter, smoothPrior);
#                       if 0
		double chi2II = MyChi2Test(hRefs->at(iter), newUnf2, nFirstBinsToSkip,
					   nLastBinsToSkip, res, hRefs->at(iter), true);
#                       else
		double chi2II = MyChi2Test(newUnf, newUnf2, nFirstBinsToSkip,
					   nLastBinsToSkip, res, hRefs->at(iter), true);
#                       endif   
#               endif

		accII[iter] += chi2II;
		acc2II[iter] += pow(chi2II, 2);
		for(int iRes = 0; iRes < nRes; ++iRes){
		    hResII->at(iter)->Fill(hResII->at(iter)->GetBinCenter(1 + iRes), res[iRes]);
		    if(hChi2II) hChi2II->at(iter)->Fill(chi2II);
		}
	    }
	    //	    if(iter>0){
	    //		double chi2II = MyChi2Test(hRefs->at(iter), hUnfs->at(iter), nFirstBinsToSkip, nLastBinsToSkip, res, hRefs->at(iter), false);
	    //		if(hChi2II) hChi2II->at(iter)->Fill(chi2II);
	    //		accII[iter] += chi2II;
	    //		acc2II[iter] += chi2II*chi2II;
	    //		for(int iRes = 0; iRes < nRes; ++iRes){
	    //		    double x = pow(res[iRes], 2);
	    //		    hResII->at(iter)->Fill(hResII->at(iter)->GetBinCenter(1 + iRes), x);
	    //		}
	    //	    }
	}
	delete h;
    }

    if(herr) delete herr;
    if(herrMcErr) delete herrMcErr;
  
    std::vector<double> mean(niters+1);
    for(int i = 0; i < niters+1; ++i){
	mean[i] = acc[i] / ntoys;
	if(chi2mcErrs) chi2mcErrs->at(i) = accMcErr[i]/ntoys;
	if(rms) rms->at(i) = 1. / (ntoys-1) * sqrt(acc2[i] - acc[i]*acc[i] / ntoys);
	if(meanII) meanII->at(i) = accII[i] / ntoys;
	if(rmsII) rmsII->at(i) = 1. / (ntoys-1) * sqrt(acc2II[i] - accII[i]*accII[i] / ntoys);
    }
    double m = 0;
    for(int i = 1; i <= hmeas_->GetNbinsX(); ++i){
	double c = hRes->at(0)->GetBinContent(i);
	if(c>m) m = c;
    }
    return mean;
}


TH1* resample(const TH1* h){
    static TRandom3 rnd;
    TH1* hnew = (TH1*)h->Clone(TString::Format("%s_resampled", h->GetName()));
    for(int i = 0; i < h->GetNbinsX() + 1; ++i){
	double n = rnd.Poisson(h->GetBinContent(i));
	hnew->SetBinContent(i, n);
	hnew->SetBinError(i, sqrt(n));
    }
    return hnew;
}

TGraphAsymmErrors* statFromToy(RooUnfold::Algorithm algo, bool smoothPrior, const RooUnfoldResponse* resp,
			       const TH1* hmeas, int niters, int ntoys, TH1** hUnfStdUnc){
    TH1* hunf = unfold(algo, resp, hmeas, niters, smoothPrior);
    double x[hunf->GetNbinsX()];
    double exl[hunf->GetNbinsX()];
    double exh[hunf->GetNbinsX()];
    double y[hunf->GetNbinsX()];
    double eyl[hunf->GetNbinsX()];
    double eyh[hunf->GetNbinsX()];
    
    //Performs toy MC
    std::vector<TH1*> newUnfs(ntoys);
    for(int itoy = 0; itoy < ntoys; ++itoy){    
	TH1* newMeas = resample(hmeas);
	newUnfs[itoy] = unfold(algo, resp, newMeas, niters, smoothPrior);
    }

    //get stat unc.:
    TAxis* ax = hunf->GetXaxis();
    for(int ibin = 1; ibin <= hunf->GetNbinsX(); ++ibin){
	x[ibin-1] = ax->GetBinCenter(ibin);
	exl[ibin-1] = ax->GetBinCenter(ibin) - ax->GetBinLowEdge(ibin);
	exh[ibin-1] = ax->GetBinUpEdge(ibin) - ax->GetBinCenter(ibin); 
	std::vector<double> sortedVals(ntoys);
	for(int itoy = 0; itoy < ntoys; ++itoy){
	    sortedVals[itoy] = newUnfs[itoy]->GetBinContent(ibin);
	}
	std::sort(sortedVals.begin(), sortedVals.end());
	double quantile_low = sortedVals[unsigned(0.16 * ntoys - 0.5)];
	double quantile_high = sortedVals[unsigned(0.84 * ntoys - 0.5)];
	y[ibin-1] = hunf->GetBinContent(ibin);
	eyl[ibin-1] = quantile_high - y[ibin-1];
	eyh[ibin-1] = y[ibin-1] - quantile_low;
    }

    if(hUnfStdUnc){
	*hUnfStdUnc = hunf;
    } else{
	hunf->SetDirectory(0);
	delete hunf;
    }
    
    //release memory:
    for(int itoy = 0; itoy < ntoys; ++itoy){
	newUnfs[itoy]->SetDirectory(0);
	delete newUnfs[itoy];
    }

    TGraphAsymmErrors* gr = new TGraphAsymmErrors(hunf->GetNbinsX(), x, y, exl, exh, eyl, eyh);
    gr->SetName(TString::Format("grUnfNewStat_%s", hmeas->GetName()));
    return gr;
}

void
makeChi2Plot(TH1* hChi2, int thisAlgoNIters, int chosenNIters,
	     double thisAlgoThr, double chosenAlgoThr,
	     const char* outDir, const char* fileBaseName, const char* canvasName, 
	     int chi2Ylog, std::vector<TH1*> hAltChi2, std::vector<std::string> labels,
	     bool limitYmax){
    
    if(hChi2==0){
	std::cerr << "makeChi2Plot called with a null pointer!. Ignored\n";
	return;
    }
  
    if(hChi2->GetMaximum() > 100 && chi2Ylog==2) chi2Ylog = 1;
    else chi2Ylog = 0;
    
    double ymin = chi2Ylog ? 0.1: 0;
    double ymax = max(1.3, 1.1*hChi2->GetMaximum());
    hChi2->GetYaxis()->SetRangeUser(ymin, ymax);
  
    //get y-coordination of a point moved by a fraction of the graphical y-scale
    auto moveup = [=](double y, double up){
	if(!chi2Ylog) return y + (ymax-ymin) *  up;
	else return y * std::pow(ymax/ymin, up);
    };
    double arrY1 = moveup(hChi2->GetBinContent(hChi2->GetXaxis()->FindBin(thisAlgoNIters)), 0.01);
    //  double arrY1 = min(moveup(thisAlgoThr, 0.01), ymax);
    double arrY0 = moveup(arrY1, 0.05);
  
    TArrow *arrowChi2 = 0;
    if(thisAlgoNIters > 0){
	arrowChi2 = new TArrow(thisAlgoNIters, arrY0, thisAlgoNIters, arrY1, 0.02, "|>");
	arrowChi2->SetLineColor(kRed);
	arrowChi2->SetFillColor(kRed);
	arrowChi2->SetLineWidth(2);
    }

    if(chosenNIters == thisAlgoNIters) arrY1 = moveup(arrY0, 0.01);
    else arrY1 = moveup(hChi2->GetBinContent(hChi2->GetXaxis()->FindBin(chosenNIters)), 0.01);
    arrY0 = moveup(arrY1, 0.05);  
    TArrow *arrowChi2_chosen = new TArrow(chosenNIters+0.01, arrY0, chosenNIters+0.01, arrY1, 0.02, "|>");
    arrowChi2_chosen->SetLineColor(kGreen+2);
    arrowChi2_chosen->SetFillColor(kGreen+2);
    arrowChi2_chosen->SetLineWidth(2);

    TLine* line = 0;
    if(thisAlgoThr > 0){
	TAxis* ax = hChi2->GetXaxis();
	double xmin = ax->GetBinLowEdge(ax->GetFirst());
	double xmax = ax->GetBinUpEdge(ax->GetLast());
	line = new TLine(xmin, thisAlgoThr, xmax, thisAlgoThr);
	line->SetLineColor(kBlack);
	line->SetLineStyle(kDashed);
	line->SetLineWidth(2);
    }

    TString cName;
    if(canvasName) cName =canvasName;
    else cName = TString("c") + hChi2->GetName();

    TCanvas *c = new TCanvas(cName.Data(), cName.Data(), 600, 600);
    c->cd();
    c->SetLogy(chi2Ylog);
    c->SetGrid();
    int icol = 2;
    hChi2->SetLineColor(2);
    hChi2 = hChi2->DrawCopy("HIST");
    if(arrowChi2) arrowChi2->Draw();
    arrowChi2_chosen->Draw();
    if(line) line->Draw();
    TLegend* l = new TLegend(0.689, 0.78, .98, .935);
    if(labels.size() > 0) l->AddEntry(hChi2, labels[0].c_str(), "l");
    unsigned ilabel = 0;
    for(auto h: hAltChi2){
	++ilabel;
	if(!h) continue;
	TH1* h1 = (TH1*) h->Clone();
	h1->SetLineColor(++icol);
	h1->Draw("same,HIST");
	if(labels.size() > ilabel) l->AddEntry(h1, labels[ilabel].c_str(), "l");
    }
  
    if(labels.size() > 0){
	l->SetTextFont(43);
	l->SetTextSize(16);
	l->DrawClone();
    }
    l->Delete();

    fixYscale(1.1, 1.1, 1.e4);	    

    if(limitYmax){
	double maxSecondBinAlt = -std::numeric_limits<double>::max();
	for(auto h: hAltChi2){
	    if(h->GetNbinsX() < 1) continue;
	    maxSecondBinAlt = std::max(h->GetBinContent(2) + h->GetBinError(1), maxSecondBinAlt);
	}
	double maxPrim = -std::numeric_limits<double>::max();
	for(int ibin = 1; ibin < hChi2->GetNbinsX(); ++ibin){
	    maxPrim = std::max(maxPrim, hChi2->GetBinContent(ibin) + hChi2->GetBinError(ibin));
	}
	double ymax = 1.1 * std::max(maxPrim, maxSecondBinAlt);
	hChi2->GetYaxis()->SetRangeUser(gPad->GetUymin(), ymax);
	c->Paint();
    }
    
    hChi2->Write();
    c->Write();
    saveCanvas(c, outDir, fileBaseName);
}


void testWithForgesResp(TH1* hUnf, int nIters, const char* outDir,
			const char* outfileBaseName, int nFirstBinsToSkip, int nLastBinsToSkip){
    bool smoothPrior = false;
    int nbins = hUnf->GetNbinsX();
    double xmin = hUnf->GetXaxis()->GetBinLowEdge(1);
    double xmax = hUnf->GetXaxis()->GetBinUpEdge(nbins);
    TH2* h2resp = new TH2D("resp", "resp",
			   nbins, xmin, xmax,
			   nbins, xmin, xmax);
    TH1* hmeas = (TH1*) hUnf->Clone("hmeasDiagTest");
    hmeas->Reset();
    for(int igen = 1; igen <= nbins; ++igen){
	int ireco = igen;
	double eff = sin(igen*TMath::Pi()/2./nbins*(igen-0.5));
	int ngen = hUnf->GetBinContent(igen);
	int nreco = eff*ngen;
	h2resp->SetBinContent(ireco, igen, nreco);
	h2resp->SetBinError(ireco, igen, sqrt(nreco));
	hmeas->SetBinContent(ireco, nreco);
	hmeas->SetBinError(ireco, sqrt(nreco));
	std::cout << "Eff: " <<  eff << "\n";
    }
  
    RooUnfoldResponse* resp = new RooUnfoldResponse(hmeas, hUnf, h2resp);

    std::auto_ptr<std::vector<TH1*> > hUnfs (new std::vector<TH1*>);
    unfoldWithErr(RooUnfold::kBayes, resp,  hmeas,  nIters, smoothPrior, hUnfs.get());

    TH1D *hchi2 = new TH1D("hchi2Diag", "hchi2Diag", nIters + 1, -0.5, nIters + .5);
    hchi2->SetTitle(TString::Format("Reco #chi^{2}/ndf for %s",
				    hUnf->GetName()));
    hchi2->GetYaxis()->SetTitle("#chi^{2}/ndf");
    hchi2->GetYaxis()->SetTitleOffset(1.40);
    hchi2->GetXaxis()->SetTitle("number of iterations for a diagonal response matrix");
    hchi2->GetXaxis()->CenterTitle();
    hchi2->GetXaxis()->SetNdivisions(nIters, 0, 0);
    hchi2->GetXaxis()->SetLabelSize(0.03);
    hchi2->SetLineWidth(2);

    for(int iter = 1; iter <= nIters; ++iter){
	TH1* hrefold = foldUnfData(hUnfs->at(iter), 0, resp);
	double chi2 = MyChi2Test(hmeas, hrefold, nFirstBinsToSkip, nLastBinsToSkip, 0, hmeas, false);
	hchi2->SetBinContent(hchi2->GetXaxis()->FindBin(iter), chi2);
    }
    hchi2->Write();

    makeChi2Plot(hchi2, 0, 0, 0., 0.,
		 outDir, outfileBaseName, "cChi2");
}


//int readNiters(const char* lepSel, const char* variable){
//    int niters;
//    TString fin_name = TString::Format("%s_niters.txt", lepSel);
//    std::ifstream fin(fin_name);
//    if(!fin.good()){
//	std::cerr << "File " << fin_name << " required for xvalIter = 6 option was not found.\n";
//	abort();
//    }
//    while(fin.good() && !fin.eof()){
//	std::string var; int n;
//	fin >> var;
//	fin >> n;
//	if(var==variable){
//	    niters = n;
//	}		
//    }
//    if(niters < 0){
//	std::cerr << "No entry found for variable " << variable << " in " << fin_name << "\n";
//	abort();
//    }
//    return niters;
//}
