//-*- c-basic-offset: 4; -*-
#define PI 3.14159265359
#define DEBUG 0

#include "ZJets_newformat.h"
#include "ConfigVJets.h"
#include "LHAPDF/LHAPDF.h"
#include "functions.h"
#include "standalone_LumiReWeighting.h"
#include <TCanvas.h>
#include <TDatime.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex.h>
#include <sstream>
#include <sys/time.h>
#include <sys/types.h>
//#include "rochcor2015.h"
#include "TSystem.h"

extern ConfigVJets cfg; // defined in runZJets_newformat.cc

using namespace std;

void ZJets::Loop(bool hasRecoInfo,
                 bool hasGenInfo,
                 int jobNum,
                 int nJobs,
                 TString pdfSet,
                 int pdfMember,
                 double muR,
                 double muF,
                 double yieldScale)
{

    //--- Random generator necessary for BTagging ---
    TRandom3 *RandGen = new TRandom3();
    //--------------------------------------------
    doRochester = true;
    // doRochester = false;
    // rmcor = new rochcor2015();

    //--- Initialize PDF from LHAPDF if needed ---
    if (pdfSet != "") initLHAPDF(pdfSet, pdfMember);
    //--------------------------------------------

    // store job id
    // using unique name for jobinfo to prevent hadd to merge them
    if (jobNum > 0) JobInfo->SetName(TString::Format("%s_%d", JobInfo->GetName(), jobNum));
    JobInfo->SetBinContent(kJobNum, jobNum);
    JobInfo->SetBinContent(kNJobs, nJobs);

    //--- Counters to check the yields ---
    Long64_t nEvents(0);
    unsigned int nEventsVInc0Jets(0), nEventsVInc0JetsNoTrig(0), nEventsVInc1Jets(0),
        nEventsVInc2Jets(0), nEventsVInc3Jets(0);
    unsigned int nGenEventsVInc0Jets(0), nGenEventsVInc1Jets(0), nGenEventsVInc2Jets(0),
        nGenEventsVInc3Jets(0);
    unsigned int nEventsWithTwoGoodLeptons(0), nEventsWithTwoGoodLeptonsWithOppCharge(0),
        nEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass(0);
    unsigned int nGenEventsWithTwoGoodLeptons(0), nGenEventsWithTwoGoodLeptonsWithOppCharge(0),
        nGenEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass(0);

    double nEffEventsVInc0Jets(0), nEffEventsVInc0JetsNoTrig(0), nEffEventsVInc1Jets(0),
        nEffEventsVInc2Jets(0), nEffEventsVInc3Jets(0);
    double nEffGenEventsVInc0Jets(0), nEffGenEventsVInc1Jets(0), nEffGenEventsVInc2Jets(0),
        nEffGenEventsVInc3Jets(0);
    double nEffEventsWithTwoGoodLeptons(0), nEffEventsWithTwoGoodLeptonsWithOppCharge(0),
        nEffEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass(0);
    double nEffGenEventsWithTwoGoodLeptons(0), nEffGenEventsWithTwoGoodLeptonsWithOppCharge(0),
        nEffGenEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass(0);

    unsigned int nEventsPassingTrigger(0);
    double nEffEventsPassingTrigger(0);

    bool UnfoldUnc = cfg.getB("unfoldUnc", false);
    if (UnfoldUnc) {
        std::cout << "Reweighting mode. MC will be reweighted to compute an alternative response "
                     "matrix to be used for the estimation of the model dependency.\n";
    }

    //------------------------------------

    //==========================================================================================================//
    //    int ZMCutLow(71), ZMCutHigh(111);
    double ZMCutLow = cfg.getD("ZMassMin", 71.);
    double ZMCutHigh = cfg.getD("ZMassMax", 111.);
    int MTCutLow(50), METCutLow(0);
    // additional variables
    double ZptRange[6] = {0, 40, 80, 120, 160, 1000};
    int LeptonID(11);
    if (lepSel == "DMu" || lepSel == "SMu") LeptonID = 13;

    muIso_ = cfg.getD("muRelIso");
    eIso_ = cfg.getD("elRelIso");
    bool pogSF = cfg.getB("pogSF", true);

    bool doPuReweight = cfg.getB("doPuReweight", true);

    //==========================================================================================================//
    //         Output file name           //
    //===================================//
    CreateOutputFileName(pdfSet, pdfMember, muR, muF, nJobs == 1 ? 0 : jobNum);
    std::cerr << "Histogram file name " << outputFileName << " for sample " << sampleLabel_ << "\n";
    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    //==========================================================================================================//

    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    //       Load efficiency tables        //
    //====================================//
    // table TableJESunc("EfficiencyTables/JESUnce_FT_53_V21_AN4_Uncertainty_AK5PFchs.txt");
    table TableJESunc("EfficiencyTables/JECUncertainty_Summer15_25nsV5_Data_AK4PF.txt");
    // table TableJESunc("EfficiencyTables/JECUncertainty_Summer15_25nsV6_AK4PF.txt");

    // muon SF
    table TrigMu17Mu8SF("EfficiencyTables/EfficienciesAndSF_ALL_Mu17Mu8OR.txt");

    ////////////////////////////////////////////////////////////////////////
    // 8TeV
    // table LeptIso, LeptID;
    // table
    // SC_RunABCD_TightID("EfficiencyTables/Muon_IDTight_Efficiencies_Run_2012ABCD_53X_Eta_Pt.txt");
    // table
    // SC_RunABCD_LooseIso("EfficiencyTables/Muon_ISOLoose_forTight_Efficiencies_Run_2012ABCD_53X_Eta_Pt.txt");
    // LeptID = SC_RunABCD_TightID;
    // LeptIso = SC_RunABCD_LooseIso;

    // electron SF
    // table Ele_Rec;
    // table Ele_Rec_8TeV("EfficiencyTables/Ele_SF_Reconstruction_2012.txt");
    // table SC_Ele_2012EA("EfficiencyTables/Ele_SF_EA2012.txt");
    ////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // 13TeV 2016 data
    //
    // new for 13 TeV SF for Id+Iso
    table TrigIsoMu24SF("EfficiencyTables/Efficiency_SF_IsoMu24_eta2p1.txt");

    // Muon Aram's numbers (Id+Iso)
    table Iso_TightID13TeV("EfficiencyTables/ratios.txt");
    table MuIdIso = Iso_TightID13TeV;

    // Muon POG numbers
    // table MuIso("EfficiencyTables/Muon_ISOLoose_forTight_Efficiencies_RunD_2015_Eta_Pt.txt");
    // table MuId("EfficiencyTables/Muon_IDTight_Efficiencies_RunD_2015_Eta_Pt.txt");
    table MuIso("EfficiencyTables/EfficienciesAndSF_ISO_2016.txt");
    table MuId("EfficiencyTables/EfficienciesAndSF_ID_2016.txt");
    table MuTrig = TrigMu17Mu8SF;

    // table MuTracking("EfficiencyTables/EfficienciesAndSF_All_Mutracking.txt");

    // table for electron SF
    table ElId("EfficiencyTables/Electron_Id_2015D_SF.txt");
    table ElReco("EfficiencyTables/Electron_Reco_2015D_SF.txt");
    table ElTrig("EfficiencyTables/Trig_Ele17_Ele12_2015D_SF.txt");

    // TTbar SF
    table TTbarSF("EfficiencyTables/sf_multi.txt");
    // Ele_Rec = Ele_Rec_8TeV;
    // if (lepSel == "DE" || lepSel == "SE") LeptID = SC_Ele_2012EA;
    // else if (lepSel == "SMu") LeptTrig = TrigIsoMu24SF;

    // Load Rochester corrections
    rc = new RoccoR("rcdata.2016.v3");
    //==========================================================================================================//
    ////////////////////////////////////////////////////////////////////////

    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    //     Systematics: jec, pu, xsec     //
    //====================================//
    cout << "Lepton Flavor: " << lepSel << "  systematics: " << systematics
         << "  direction: " << direction << endl;

    int mode = (systematics == 1) ? direction : 0;
    // standalone_LumiReWeighting puWeight(2015250, mode);
    int year = cfg.getI("puReweightYear", 2016);
    standalone_LumiReWeighting puWeight(year, mode);

    int scale(0); // 0,+1,-1; (keep 0 for noJEC shift study)
    if (systematics == 2) scale = direction;

    xsecFactor_ = 1.;
    if (systematics == 3) xsecFactor_ = 1. + direction * xsecUnc;

    int smearJet(0);
    if (systematics == 4) smearJet = direction;

    int lepscale(0);
    if (systematics == 5) lepscale = direction;

    int smearlepton(0);
    if (systematics == 6) smearlepton = direction;
    //==========================================================================================================//

    // setting weight when running on MIX of exclusive DY/WJets files to match number of parton
    // events
    double mixingWeightsDY[4] = {
        0.192686,
        0.0718097,
        0.04943495,
        0.0360337}; // here we match all partons, and combine electron and muon side
    double mixingWeightsWJ_SMu[4] = {0.366713, 0.1119323, 0.07641136, 0.03803325};

    //======================================================================
    // additionnal PU weights
    TH1 *addPuWeights = 0;
    string addPuFile = cfg.getS("additionalPuWeightFile");
    if (addPuFile.size() > 0) {
        TFile f(addPuFile.c_str());
        if (f.IsZombie()) {
            std::cerr << "PU reweighting file, " << addPuFile << " was not found!" << endl;
        } else {
            cout << "Event will be reweighting according to their number of vertices using weights "
                    "from file "
                 << addPuFile << "." << endl;
            f.GetObject("hWeights", addPuWeights);
            addPuWeights->SetDirectory(0);
        }
    }
    //======================================================================
    // Data/MC ration  needed for unfolding uncertainties

    // jet pt
    TH1D *FirstJetPt_2_Zinc1jetHratio_fit = 0;
    TH1D *SecondJetPt_2_Zinc2jetHratio_fit = 0;
    TH1D *ThirdJetPt_2_Zinc3jetHratio_fit = 0;
    TH1D *JetsHT_2_Zinc1jetHratio_fit = 0;
    TH1D *JetsHT_2_Zinc2jetHratio_fit = 0;
    TH1D *JetsHT_2_Zinc3jetHratio_fit = 0;
    TH1D *FirstJetAbsRapidity_2_Zinc1jetHratio_fit = 0;
    TH1D *SecondJetAbsRapidity_2_Zinc2jetHratio_fit = 0;
    TH1D *ThirdJetAbsRapidity_2_Zinc3jetHratio_fit = 0;

    TH1D *JZB_2Hratio_fit = 0;
    TH1D *JZB_ptHigh_2Hratio_fit = 0;
    TH1D *JZB_ptLow_2Hratio_fit = 0;

    TH1D *VisPt_2_Zinc1jetQunHratio_fit = 0;
    TH1D *VisPt_2_Zinc2jetQunHratio_fit = 0;
    TH1D *VisPt_2_Zinc3jetQunHratio_fit = 0;

    TH1D *ZNGoodJets_ZexcHratio_fit = 0;

    if (UnfoldUnc) {

        if (lepSel == "DMu") {
            // jet pt
            TFile *fratio1jpt = new TFile("RootRatios/FirstJetPt_2_Zinc1jetHratio_fit.root");
            FirstJetPt_2_Zinc1jetHratio_fit = (TH1D *)fratio1jpt->Get("Hratio");

            TFile *fratio2jpt = new TFile("RootRatios/SecondJetPt_2_Zinc2jetHratio_fit.root");
            SecondJetPt_2_Zinc2jetHratio_fit = (TH1D *)fratio2jpt->Get("Hratio");

            TFile *fratio3jpt = new TFile("RootRatios/ThirdJetPt_2_Zinc3jetHratio_fit.root");
            ThirdJetPt_2_Zinc3jetHratio_fit = (TH1D *)fratio3jpt->Get("Hratio");

            // jet Ht
            TFile *fratioht = new TFile("RootRatios/JetsHT_2_Zinc1jetHratio_fit.root");
            JetsHT_2_Zinc1jetHratio_fit = (TH1D *)fratioht->Get("Hratio");

            TFile *fratioht2 = new TFile("RootRatios/JetsHT_2_Zinc2jetHratio_fit.root");
            JetsHT_2_Zinc2jetHratio_fit = (TH1D *)fratioht2->Get("Hratio");

            TFile *fratioht3 = new TFile("RootRatios/JetsHT_2_Zinc3jetHratio_fit.root");
            JetsHT_2_Zinc3jetHratio_fit = (TH1D *)fratioht3->Get("Hratio");

            // jet rapidity
            TFile *fratio1jrapidity =
                new TFile("RootRatios/FirstJetAbsRapidity_2_Zinc1jetHratio_fit.root");
            FirstJetAbsRapidity_2_Zinc1jetHratio_fit = (TH1D *)fratio1jrapidity->Get("Hratio");

            TFile *fratio2jrapidity =
                new TFile("RootRatios/SecondJetAbsRapidity_2_Zinc2jetHratio_fit.root");
            SecondJetAbsRapidity_2_Zinc2jetHratio_fit = (TH1D *)fratio2jrapidity->Get("Hratio");

            TFile *fratio3jrapidity =
                new TFile("RootRatios/ThirdJetAbsRapidity_2_Zinc3jetHratio_fit.root");
            ThirdJetAbsRapidity_2_Zinc3jetHratio_fit = (TH1D *)fratio3jrapidity->Get("Hratio");

            // pt-balance
            TFile *fratio1jptbal = new TFile("RootRatios/VisPt_2_Zinc1jetQunHratio_fit.root");
            VisPt_2_Zinc1jetQunHratio_fit = (TH1D *)fratio1jptbal->Get("Hratio");

            TFile *fratio2jptbal = new TFile("RootRatios/VisPt_2_Zinc2jetQunHratio_fit.root");
            VisPt_2_Zinc2jetQunHratio_fit = (TH1D *)fratio2jptbal->Get("Hratio");

            TFile *fratio3jptbal = new TFile("RootRatios/VisPt_2_Zinc3jetQunHratio_fit.root");
            VisPt_2_Zinc3jetQunHratio_fit = (TH1D *)fratio3jptbal->Get("Hratio");

            // JZB
            TFile *fratioJZB = new TFile("RootRatios/JZB_2Hratio_fit.root");
            JZB_2Hratio_fit = (TH1D *)fratioJZB->Get("Hratio");

            TFile *fratioJZBlow = new TFile("RootRatios/JZB_ptLow_2Hratio_fit.root");
            JZB_ptLow_2Hratio_fit = (TH1D *)fratioJZBlow->Get("Hratio");

            TFile *fratioJZBhigh = new TFile("RootRatios/JZB_ptHigh_2Hratio_fit.root");
            JZB_ptHigh_2Hratio_fit = (TH1D *)fratioJZBhigh->Get("Hratio");

            //  multiplicity
            TFile *fratioNJexc = new TFile("RootRatios/ZNGoodJets_ZexcHratio_fit.root");
            ZNGoodJets_ZexcHratio_fit = (TH1D *)fratioNJexc->Get("Hratio");
        }

        if (lepSel == "DE") {
            // jet pt
            TFile *fratio1jpt = new TFile("RootRatios/FirstJetPt_2_Zinc1jetHratio_ee_fit.root");
            FirstJetPt_2_Zinc1jetHratio_fit = (TH1D *)fratio1jpt->Get("Hratio");

            TFile *fratio2jpt = new TFile("RootRatios/SecondJetPt_2_Zinc2jetHratio_ee_fit.root");
            SecondJetPt_2_Zinc2jetHratio_fit = (TH1D *)fratio2jpt->Get("Hratio");

            TFile *fratio3jpt = new TFile("RootRatios/ThirdJetPt_2_Zinc3jetHratio_ee_fit.root");
            ThirdJetPt_2_Zinc3jetHratio_fit = (TH1D *)fratio3jpt->Get("Hratio");

            // jet Ht
            TFile *fratioht = new TFile("RootRatios/JetsHT_2_Zinc1jetHratio_ee_fit.root");
            JetsHT_2_Zinc1jetHratio_fit = (TH1D *)fratioht->Get("Hratio");

            TFile *fratioht2 = new TFile("RootRatios/JetsHT_2_Zinc2jetHratio_ee_fit.root");
            JetsHT_2_Zinc2jetHratio_fit = (TH1D *)fratioht2->Get("Hratio");

            TFile *fratioht3 = new TFile("RootRatios/JetsHT_2_Zinc3jetHratio_ee_fit.root");
            JetsHT_2_Zinc3jetHratio_fit = (TH1D *)fratioht3->Get("Hratio");

            // jet rapidity
            TFile *fratio1jrapidity =
                new TFile("RootRatios/FirstJetAbsRapidity_2_Zinc1jetHratio_ee_fit.root");
            FirstJetAbsRapidity_2_Zinc1jetHratio_fit = (TH1D *)fratio1jrapidity->Get("Hratio");

            TFile *fratio2jrapidity =
                new TFile("RootRatios/SecondJetAbsRapidity_2_Zinc2jetHratio_ee_fit.root");
            SecondJetAbsRapidity_2_Zinc2jetHratio_fit = (TH1D *)fratio2jrapidity->Get("Hratio");

            TFile *fratio3jrapidity =
                new TFile("RootRatios/ThirdJetAbsRapidity_2_Zinc3jetHratio_ee_fit.root");
            ThirdJetAbsRapidity_2_Zinc3jetHratio_fit = (TH1D *)fratio3jrapidity->Get("Hratio");

            // pt-balance
            TFile *fratio1jptbal = new TFile("RootRatios/VisPt_2_Zinc1jetQunHratio_ee_fit.root");
            VisPt_2_Zinc1jetQunHratio_fit = (TH1D *)fratio1jptbal->Get("Hratio");

            TFile *fratio2jptbal = new TFile("RootRatios/VisPt_2_Zinc2jetQunHratio_ee_fit.root");
            VisPt_2_Zinc2jetQunHratio_fit = (TH1D *)fratio2jptbal->Get("Hratio");

            TFile *fratio3jptbal = new TFile("RootRatios/VisPt_2_Zinc3jetQunHratio_ee_fit.root");
            VisPt_2_Zinc3jetQunHratio_fit = (TH1D *)fratio3jptbal->Get("Hratio");

            // JZB
            TFile *fratioJZB = new TFile("RootRatios/JZB_2Hratio_ee_fit.root");
            JZB_2Hratio_fit = (TH1D *)fratioJZB->Get("Hratio");

            TFile *fratioJZBlow = new TFile("RootRatios/JZB_ptLow_2Hratio_ee_fit.root");
            JZB_ptLow_2Hratio_fit = (TH1D *)fratioJZBlow->Get("Hratio");

            TFile *fratioJZBhigh = new TFile("RootRatios/JZB_ptHigh_2Hratio_ee_fit.root");
            JZB_ptHigh_2Hratio_fit = (TH1D *)fratioJZBhigh->Get("Hratio");

            //  multiplicity
            TFile *fratioNJexc = new TFile("RootRatios/ZNGoodJets_ZexcHratio_ee_fit.root");
            ZNGoodJets_ZexcHratio_fit = (TH1D *)fratioNJexc->Get("Hratio");
        }
    }

    cout << endl;
    stringstream s;
    cout << "\nProcessing : " << fileName << "\n";
    s << "\t--> " << outputFileName << "\n";
    cout << s.str();
    for (unsigned i = 0; i < s.str().size(); ++i) cout << "-";
    cout << "\n\n" << flush;

    //--- Initialize the tree branches ---
    Init(hasRecoInfo, hasGenInfo);
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntries();

    //------------------------------------

    struct timeval t0;
    int mess_every_n = std::min(1000LL, nentries / 10);
    if (mess_every_n < 1) mess_every_n = 1;

    // ------ Random number for lepton energy resolution smearing -----
    // TRandom* RamMu = new TRandom(10);
    // TRandom* RamEle = new TRandom(20);
    // --------------------------------

    // event yield normalisation for MC
    // norm_ = yieldScale;

    double prev_rate = 0;

    Long64_t nEventsToProcessTot = nentries;
    if (nMaxEvents >= 0 && nEventsToProcessTot > nMaxEvents) nEventsToProcessTot = nMaxEvents;

    Long64_t entry_start = 0;
    Long64_t entry_stop = nEventsToProcessTot;
    if (nJobs > 1) {
        Long64_t eventsPerJob = nEventsToProcessTot / nJobs;
        entry_start = eventsPerJob * (jobNum - 1);
        if (jobNum < nJobs)
            entry_stop = entry_start + eventsPerJob;
        else
            entry_stop = nEventsToProcessTot;
    }

    int nEventsToProcess = entry_stop - entry_start;
    // if(nMaxEvents >= 0 && nEventsToProcess > nMaxEvents) nEventsToProcess = nMaxEvents;
    cout << "We will run on " << nEventsToProcess << " events" << endl;

    processedEventMcWeightSum_ = 0.;
    nEvents = 0;
    double weightSum = 0;
    double genWeightSum = 0;

    if (!GLepBarePrompt) {
        std::cout << "Warning: tau gen veto was not found in the ntuple. It is fine if the samples "
                     "is not DY or if it does not contains Z->\\tau\\tau.\n";
    }

    bool DJALOG = false;
    // entry_start = 19000000;
    // entry_stop = nentries;

    if (DJALOG) printf("{DJA LOG}    nEventsToProcess = %d\n", nEventsToProcess);
    if (DJALOG) printf("{DJA LOG}    nMaxEvents = %ld\n", nMaxEvents);
    if (DJALOG) printf("{DJA LOG}    entry_start = %lld\n", entry_start);
    if (DJALOG) printf("{DJA LOG}    entry_stop = %lld\n", entry_stop);
    if (DJALOG) printf("{DJA LOG}    nentries = %lld\n", nentries);

    //======================================================================
    // Event loop starts here
    //======================================================================
    for (Long64_t jentry = entry_start; jentry < entry_stop; jentry += 1) {
        if (0 <= nMaxEvents && nMaxEvents <= nEvents) break;

        if ((jentry % 1000000) == 0)
            if (DJALOG) printf("{DJA LOG}    Trying Entry = %lld\n", jentry);

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        // cout << "---------------------------------------------------------------------" << endl;

        if (nEvents == mess_every_n) gettimeofday(&t0, 0);
        if (nEvents % mess_every_n == 0 && nEvents > mess_every_n) {
            timeval t1;
            gettimeofday(&t1, 0);
            double rate = ((t1.tv_sec - t0.tv_sec) + 1.e-6 * (t1.tv_usec - t0.tv_usec)) /
                          (nEvents - mess_every_n);
            if (fabs(rate / prev_rate - 1.) > 0.1) {
                prev_rate = 0.5 * (prev_rate + rate);
            }
            double rem = prev_rate * (nEventsToProcess - nEvents);
            int rem_s = int(rem + 0.5);
            int rem_h = int(rem_s) / 3600;
            rem_s -= rem_h * 3600;
            int rem_m = rem_s / 60;
            rem_s -= rem_m * 60;
            cout << "\r" << TString::Format("%4.1f%%", (100. * nEvents) / nEventsToProcess) << " "
                 << std::setw(11) << nEvents << " / " << nEventsToProcess << " " << std::setw(7)
                 << int(prev_rate * 1.e6 + 0.5) << " us/event"
                 << " Remaining time for this dataset: " << std::setw(2) << rem_h << " h "
                 << std::setw(2) << rem_m << " min " << std::setw(2) << rem_s << " s" << std::flush;
        }

        if (fChain->GetEntry(jentry) == 0) {
            std::cerr << "Failed to read Tree entry " << jentry << "!\n";
            continue;
        }

        //	if(nEvents == 0 && !EvtIsRealData){
        //	    if(xsec_ == 0){
        //		std::cerr << "Warning: cross section value for MC sample " << sampleLabel_
        //			  <<  " is null or was not specified. We will assume the event"
        //			  << " weights are normalizes such that the cross section on pb "
        //			  << " is equal to the sum of weights divivided by the numnber of
        // events\n";
        //		norm_ = yieldScale * lumi_ * 1. * xsecFactor_ * skimAccep_[0];
        //	    } else{
        //		norm_ = yieldScale * lumi_ * xsec_ * xsecFactor_ * skimAccep_[0];
        //	    }
        //	    if(norm_ == 0){
        //		std::cerr << "Error: normaliation factor for sample " << fileName
        //			  << " is null! Aborts at " __FILE__ ":"
        //			  << __LINE__ << "." << std::endl;
        //		abort();
        //	    }
        //	}

        //=======================================================================================================//
        //         Continue Statements        //
        //====================================//
        // if (nEvents % 2 == 0) continue;
        // if (EvtVtxCnt <= 14) continue;

        //=======================================================================================================//
        // if (EvtIsRealData && EvtRunNum >= 282037 && EvtRunNum <=283469){

        // if(EvtIsRealData && (EvtRunNum < 282037 || EvtRunNum > 283469))
        // continue;
        const UInt_t runThreshold = 278820; // start of Run G
        // const UInt_t runThreshold = 276811;  //This is the end of Run D

        /*
                if( EvtIsRealData ){
                    if( EvtRunNum > runThreshold ){
                        if( (jentry % 1000000) == 0 )
                            if(DJALOG)    printf("{DJA LOG}    Run = %d\n",EvtRunNum);
                             printf("{DJA LOG} This is from Run GH (run = %d) we also want this\n",
           EvtRunNum);
                        //continue;
                    }else{
                        if(DJALOG)    printf("{DJA LOG} This is from Run B-F (run = %d) we want
           this\n", EvtRunNum);
                    }
                }
        */
        //=======================================================================================================//

        nEvents++;

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

        //=======================================================================================================//
        //         Computing weight            //
        //====================================//
        //        double weight = norm_;
        double weight = 1;

        if (hasRecoInfo && doPuReweight && !EvtIsRealData) {
            // std::cout << "PU weight: " << EvtPuCntTruth << " , " <<
            // puWeight.weight(EvtPuCntTruth) << "\n";
            weight *= puWeight.weight(EvtPuCntTruth);
        }

        if (addPuWeights) {
            double add_w_ =
                addPuWeights->GetBinContent(addPuWeights->GetXaxis()->FindBin(EvtVtxCnt));
            if (add_w_ > 0) weight *= add_w_;
        }

        if (fileName.Index("DYJets") >= 0 && fileName.Index("MIX") >= 0 && GNup > 5)
            weight *= mixingWeightsDY[GNup - 6];
        if (fileName.Index("SMu_8TeV_WJets") >= 0 && fileName.Index("MIX") >= 0 && GNup > 5)
            weight *= mixingWeightsWJ_SMu[GNup - 6];
        // if (fileName.Index("Sherpa") >= 0 && fileName.Index("UNFOLDING") >= 0) {
        //     weight *= mcSherpaWeights_->at(0) / 43597515.;
        // }

        double commonGenWeight = weight; // This variable is used in the fill(x, commonGenWeight,
                                         // EvtWeights) function for taking care of scale/PDF/alphas
                                         // uncertainty graph production if needed (see
                                         // Includes/GenH1D.h).
        if (!EvtIsRealData) {
            if (EvtWeights->size() > 0) {
                weight *= (*EvtWeights)[0];
                processedEventMcWeightSum_ += (*EvtWeights)[0];
            } else {
                processedEventMcWeightSum_ += 1.;
            }
        }

        if (fileName.Index("MLM") >= 0 && EvtWeights->size() > 117) {
            // MLM sample includes weights for several PDFs.
            // Keep only the NNPDF3.0 ones, as current GenH1D implemention
            // won't work with these extra weights.
            bool static weightTrimWarning = true;
            EvtWeights->resize(110);
            if (weightTrimWarning) {
                std::cerr << "The MC sample contains more than 117 event weights. Only the 110 "
                             "first ones will be considered and assumed to have the same defintion "
                             "than for the FxFx sample (nominal + 9 scale variations + 100 NNPF "
                             "replicas.\n";
                weightTrimWarning = false;
            }
        }

        //        if (fileName.Index("Sherpa2") >= 0) {
        //            weight *= EvtWeights->at(0);
        //            weight_amcNLO_sum += EvtWeights->at(1);
        //        }
        //
        //        if (fileName.Index("mcatnlo") >= 0) {
        //            weight *= EvtWeights->at(0);
        //         //  cout << EvtWeights->at(0) << "  , " << EvtWeights->at(0) << "\n";
        //
        //            //if (muR == 0.0 && muF == 0.0 && pdfMember == -1) weight *=
        //            EvtWeights->at(0);
        //            //if (muR == 1.0 && muF == 1.0 && pdfMember == -1) weight *=
        //            EvtWeights->at(0);
        //            // CommentAG: only at(0) available for EvtWeights
        //         /*
        //            if (muR == 1.0 && muF == 2.0 && pdfMember == -1) weight *= EvtWeights->at(2);
        //            if (muR == 1.0 && muF == 0.5 && pdfMember == -1) weight *= EvtWeights->at(3);
        //            if (muR == 2.0 && muF == 1.0 && pdfMember == -1) weight *= EvtWeights->at(4);
        //            if (muR == 2.0 && muF == 2.0 && pdfMember == -1) weight *= EvtWeights->at(5);
        //            if (muR == 2.0 && muF == 0.5 && pdfMember == -1) weight *= EvtWeights->at(6);
        //            if (muR == 0.5 && muF == 1.0 && pdfMember == -1) weight *= EvtWeights->at(7);
        //            if (muR == 0.5 && muF == 2.0 && pdfMember == -1) weight *= EvtWeights->at(8);
        //            if (muR == 0.5 && muF == 0.5 && pdfMember == -1) weight *= EvtWeights->at(9);
        //            if (muR == 0.0 && muF == 0.0 && pdfMember != -1) weight *=
        //            EvtWeights->at(pdfMember+10);
        //            weight_amcNLO_sum += EvtWeights->at(1);
        //	 */
        //            weight_amcNLO_sum += EvtWeights->at(0);
        //        }

        // cout << weight << "\n";
        //==========================================================================================================//
        // Compute the weight for PDF syst    //
        //===================================//
        double wPdf(1);
        if (pdfSet != "") wPdf = computePDFWeight();
        //==========================================================================================================//

        //--- There is no pile-up so no need to reweight for that ---
        double genWeight = weight * wPdf;
        genWeightSum += genWeight;
        //=======================================================================================================//

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

        // Trigger // NEW DA
        // bool passesTrigger = ((*ourTrig_) & triggerMask_);
        bool passesTrigger = false;
        if (EvtIsRealData) {
            if (EvtRunNum < runThreshold) {
                // printf("{DJA LOG}        This is from Runs B-F\n");
                passesTrigger = ((*ourTrig_) & triggerMask_EraBG);
            }
            if (EvtRunNum >= runThreshold) {
                // printf("{DJA LOG}        This is from Run GH\n");
                passesTrigger = ((*ourTrig_) & triggerMask_EraH);
            }
        } else {
            passesTrigger = ((*ourTrig_) & triggerMask_EraH);
        }

        if (passesTrigger) {
            nEventsPassingTrigger++;
            nEffEventsPassingTrigger += weight;
        }

        //=======================================================================================================//
        //         Retrieving leptons          //
        //====================================//
        bool passesLeptonCut(0), passesLeptonCutNoMass(0), passesLeptonChargeCut(0),
            passesTauCut(1);
        // bool passesLeptonMassCut(0);
        unsigned short nLeptons(0), nVetoMuons(0), nVetoElectrons(0);
        vector<leptonStruct> leptons;
        vector<leptonStruct> vetoMuons;
        vector<leptonStruct> vetoElectrons;
        TLorentzVector MET;
        double MT = -99;
        double phi_acop = -99;
        double sinthetastar = -99, costhetastar = -99;
        double phistar = -99;
        TLorentzVector EWKBoson;

        if (hasRecoInfo) {
            //--- get Muons ---
            if (lepSel == "DMu" || lepSel == "SMu") {
                getMuons(leptons, vetoMuons);
            }

            //--- get Electrons ---
            if (lepSel == "DE" || lepSel == "SE") {
                getElectrons(leptons, vetoElectrons);
            }

            //--- get MET ---
            if (lepSel == "SMu" || lepSel == "SE") {
                int whichMET(0); // 0 - slimmedMETs  1- slimmedMETsNoHF 2- slimmedMETsPuppi
                MET.SetXYZM(METPx->at(whichMET), METPy->at(whichMET), 0, 0);
            }

            //--- get the size of the collections ---
            nLeptons = leptons.size();
            nVetoMuons = vetoMuons.size();
            nVetoElectrons = vetoElectrons.size();

            //--- sort leptons by descending pt ---
            sort(leptons.begin(), leptons.end(), LepDescendingOrder);
            sort(vetoMuons.begin(), vetoMuons.end(), LepDescendingOrder);

            sort(vetoElectrons.begin(), vetoElectrons.end(), LepDescendingOrder);

            //-- determine if the event passes the leptons requirements for EWKBoson = Z Boson
            // cout << " nLeptons " << nLeptons << "\n";

            if ((lepSel == "DMu" || lepSel == "DE") && nLeptons >= 2) {
                // --- lepton energy scale and resolution variation ---
                if (lepSel == "DMu") {
                    // --- muon energy scale variation ---
                    leptons[0].v.SetPtEtaPhiE(leptons[0].v.Pt() * (1 + lepscale * 0.002),
                                              leptons[0].v.Eta(),
                                              leptons[0].v.Phi(),
                                              leptons[0].v.E() * (1 + lepscale * 0.002));
                    leptons[1].v.SetPtEtaPhiE(leptons[1].v.Pt() * (1 + lepscale * 0.002),
                                              leptons[1].v.Eta(),
                                              leptons[1].v.Phi(),
                                              leptons[1].v.E() * (1 + lepscale * 0.002));

                } else if (lepSel == "DE") {
                    // --- electron energy scale variation ---
                    if (fabs(leptons[0].v.Eta()) < 1.479) {
                        leptons[0].v.SetPtEtaPhiE(leptons[0].v.Pt() * (1 + lepscale * 0.006),
                                                  leptons[0].v.Eta(),
                                                  leptons[0].v.Phi(),
                                                  leptons[0].v.E() * (1 + lepscale * 0.006));
                    } else {
                        leptons[0].v.SetPtEtaPhiE(leptons[0].v.Pt() * (1 + lepscale * 0.015),
                                                  leptons[0].v.Eta(),
                                                  leptons[0].v.Phi(),
                                                  leptons[0].v.E() * (1 + lepscale * 0.015));
                    }

                    if (fabs(leptons[1].v.Eta()) < 1.479) {
                        leptons[1].v.SetPtEtaPhiE(leptons[1].v.Pt() * (1 + lepscale * 0.006),
                                                  leptons[1].v.Eta(),
                                                  leptons[1].v.Phi(),
                                                  leptons[1].v.E() * (1 + lepscale * 0.006));
                    } else {
                        leptons[1].v.SetPtEtaPhiE(leptons[1].v.Pt() * (1 + lepscale * 0.015),
                                                  leptons[1].v.Eta(),
                                                  leptons[1].v.Phi(),
                                                  leptons[1].v.E() * (1 + lepscale * 0.015));
                    }
                }

                nEventsWithTwoGoodLeptons++;
                nEffEventsWithTwoGoodLeptons += weight;

                //  if(!hasGenInfo){   // CommentAG: comment this out since don't enter 'lepton
                //  energy smearing' block
                // build Electroweak boson candidate: here it is expected to be a Z

                if (!EvtIsRealData) {
                    double effWeight = 1.;
                    if (lepSel == "DMu") {
                        if (pogSF) {
                            effWeight *=
                                MuId.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].v.Eta()));
                            effWeight *=
                                MuId.getEfficiency(leptons[1].v.Pt(), fabs(leptons[1].v.Eta()));
                            effWeight *=
                                MuIso.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].v.Eta()));
                            effWeight *=
                                MuIso.getEfficiency(leptons[1].v.Pt(), fabs(leptons[1].v.Eta()));

                            effWeight *= MuTracking(leptons[0].v.Eta()); // coded in function.cc
                            effWeight *= MuTracking(leptons[1].v.Eta());

                        } else {
                            effWeight *=
                                MuIdIso.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].v.Eta()));
                        }
                        if (useTriggerCorrection)
                            effWeight *= MuTrig.getEfficiency(fabs(leptons[0].v.Eta()),
                                                              fabs(leptons[1].v.Eta()));
                    }
                    if (lepSel == "DE") {
                        effWeight *=
                            ElReco.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].scEta));
                        effWeight *=
                            ElReco.getEfficiency(leptons[1].v.Pt(), fabs(leptons[1].scEta));
                        effWeight *= ElId.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].scEta));
                        effWeight *= ElId.getEfficiency(leptons[1].v.Pt(), fabs(leptons[1].scEta));
                        // cout<< "Reco" <<ElReco.getEfficiency(leptons[1].v.Pt(),
                        // fabs(leptons[1].scEta))<<endl;
                        if (useTriggerCorrection)
                            effWeight *= ElTrig.getEfficiency(fabs(leptons[0].v.Eta()),
                                                              fabs(leptons[1].v.Eta()));
                        // cout<< "Trig "<<ElTrig.getEfficiency(fabs(leptons[0].v.Eta()),
                        // fabs(leptons[1].v.Eta()))<<endl;
                    }
                    weight *= effWeight;
                }

                weightSum += weight;

                EWKBoson = leptons[0].v + leptons[1].v;
                // apply trigger, charge, mass cut
                if (passesTrigger && (leptons[0].charge * leptons[1].charge < 0)) {
                    nEventsWithTwoGoodLeptonsWithOppCharge++;
                    nEffEventsWithTwoGoodLeptonsWithOppCharge += weight;
                    passesLeptonChargeCut = 1;

                    if (leptons[0].v.Pt() > lepPtCutMin && leptons[1].v.Pt() > lepPtCutMin) {
                        passesLeptonCutNoMass = 1;
                        if (EWKBoson.M() > ZMCutLow && EWKBoson.M() < ZMCutHigh) {

                            nEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass++;
                            nEffEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass += weight;
                            passesLeptonCut = 1;
                            // passesLeptonMassCut = 1;

                            //  ----  code Phi* AG
                            phi_acop = PI - deltaPhi(leptons[0].v, leptons[1].v);
                            if (leptons[0].charge < 0)
                                costhetastar = tanh((leptons[0].v.Eta() - leptons[1].v.Eta()) / 2.);
                            else
                                costhetastar = tanh((leptons[1].v.Eta() - leptons[0].v.Eta()) / 2.);
                            sinthetastar = sqrt(1. - costhetastar * costhetastar);
                            //  cout << "phi_acop = " << phi_acop << " , sinthetastar = " <<
                            //  sinthetastar << "\n";
                            phistar = tan(phi_acop / 2.) * sinthetastar;
                            //  cout << "phistar = " << phistar << "\n";
                        }
                    }
                }
                if ((leptons[0].charge * leptons[1].charge < 0) && EWKBoson.M() > ZMCutLow &&
                    EWKBoson.M() < ZMCutHigh && leptons[0].v.Pt() > lepPtCutMin &&
                    leptons[1].v.Pt() > lepPtCutMin) {
                    ++nEventsVInc0JetsNoTrig;
                    nEffEventsVInc0JetsNoTrig += weight;
                }
                //}
            } // end if Z study and nLeptons >= 2

            // determine if the event passes the leptons requirements for EWKBoson = W Boson
            // exactly one muon (or exactly one electron) must be present and no additional
            // charged leptons can be present.
            else if ((lepSel == "SMu" || lepSel == "SE") &&
                     (nLeptons == 1 && nVetoMuons == 0 && nVetoElectrons == 0)) {
                // add the MET to the leptons collection: leptons[1] = MET
                leptons.push_back(leptonStruct(MET.Pt(), 0, MET.Phi(), MET.Pt(), 0, 0, 0, 0, 0));

                // build Electroweak boson candidate: here it is expected to be a W
                EWKBoson = leptons[0].v + MET;
                MT = sqrt(2 * leptons[0].v.Pt() * MET.Pt() *
                          (1 - cos(leptons[0].v.Phi() - MET.Phi())));

                // apply transver mass and MET cut
                if (MT > MTCutLow && MET.Pt() > METCutLow && leptons[0].v.Pt() > lepPtCutMin)
                    passesLeptonCut = 1;

                // apply scale factors only on MC.
                if (!EvtIsRealData) {
                    double effWeight = 1.;
                    if (lepSel == "SMu") {
                        if (pogSF) {
                            effWeight *=
                                MuId.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].v.Eta()));
                            effWeight *=
                                MuId.getEfficiency(leptons[1].v.Pt(), fabs(leptons[1].v.Eta()));
                            effWeight *=
                                MuIso.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].v.Eta()));
                            effWeight *=
                                MuIso.getEfficiency(leptons[1].v.Pt(), fabs(leptons[1].v.Eta()));
                        } else {
                            effWeight *=
                                MuIdIso.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].v.Eta()));
                        }
                        if (useTriggerCorrection)
                            effWeight *= MuTrig.getEfficiency(fabs(leptons[0].v.Eta()),
                                                              fabs(leptons[1].v.Eta()));
                    } else if (lepSel == "SE") {
                        effWeight *=
                            ElReco.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].scEta));
                        effWeight *= ElId.getEfficiency(leptons[0].v.Pt(), fabs(leptons[0].scEta));
                    }
                    weight *= effWeight;
                }
            } // end if W study and nLeptons >= 1
        }     // end has reco info

        // cout << effWeight << "\n";
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //       Retrieving gen leptons        //
        //====================================//
        bool passesgenLeptonCut(0);
        bool passesgenLeptonCutNoMass(0);
        unsigned short nTotgenLeptons(0), ngenLeptons(0), nTotGenPhotons(0);
        vector<leptonStruct> genLeptons;
        vector<int> usedGenPho;
        TLorentzVector genMET;
        double genMT = -99;
        TLorentzVector genEWKBoson;
        int countTauS3 = 0;

        if (hasGenInfo) {
            // CommentAG: this line is commented because can't do countTauS3-- since status 3 is not
            // stored
            // if (hasRecoInfo) countTauS3 = (lepSel == "DMu" || lepSel == "DE") ? 2 : 1; // AG
            nTotGenPhotons = GLepClosePhotEta->size();
            nTotgenLeptons = GLepBareEta->size();
            //-- retriveing generated leptons with status 1
            for (unsigned short i(0); i < nTotgenLeptons; i++) {
                bool lepToBeConsidered(false);
                if ((lepSel == "DMu" || lepSel == "DE") && abs(GLepBareId->at(i)) == LeptonID)
                    lepToBeConsidered = true;
                else if ((lepSel == "SMu" || lepSel == "SE") &&
                         (abs(GLepBareId->at(i)) == LeptonID || abs(GLepBareId->at(i)) == 12 ||
                          abs(GLepBareId->at(i)) == 14))
                    lepToBeConsidered = true;
                // following two lines should give the same result
                if (GLepBareSt->at(i) == 3 && abs(GLepBareId->at(i)) != LeptonID &&
                    (abs(GLepBareId->at(i)) == 15 || abs(GLepBareId->at(i)) == 13 ||
                     abs(GLepBareId->at(i)) == 11))
                    countTauS3++;
                if (GLepBareSt->at(i) == 3 && abs(GLepBareId->at(i)) == LeptonID) countTauS3--;

                bool passesTauVeto = true;
                if (GLepBarePrompt) passesTauVeto = (*GLepBarePrompt)[i];

                if (!passesTauVeto) continue;

                if (!lepToBeConsidered) continue;

                int charge;
                if (abs(GLepBareId->at(i)) == 12 || abs(GLepBareId->at(i)) == 14 ||
                    abs(GLepBareId->at(i)) == 16)
                    charge = 0;
                else if (GLepBareId->at(i) < 0)
                    charge = -1;
                else
                    charge = 1;

                leptonStruct genLep(GLepBarePt->at(i),
                                    GLepBareEta->at(i),
                                    GLepBarePhi->at(i),
                                    GLepBareE->at(i),
                                    charge,
                                    0,
                                    0,
                                    0,
                                    0);
                leptonStruct genLepNoFSR(GLepBarePt->at(i),
                                         GLepBareEta->at(i),
                                         GLepBarePhi->at(i),
                                         GLepBareE->at(i),
                                         charge,
                                         0,
                                         0,
                                         0,
                                         0);

                //-- dress the leptons with photon (cone size = 0.1). Only for status 1 leptons
                //(after FSR)
                if ((GLepBareSt->at(i) == 1 && lepToBeConsidered) ||
                    ((lepSel == "SMu" || lepSel == "SE") && charge == 0)) {

                    for (unsigned short j(0); j < nTotGenPhotons; j++) {
                        TLorentzVector tmpGenPho;
                        tmpGenPho.SetPtEtaPhiM(GLepClosePhotPt->at(j),
                                               GLepClosePhotEta->at(j),
                                               GLepClosePhotPhi->at(j),
                                               0.);
                        int used(0);
                        for (unsigned short k(0); k < usedGenPho.size(); k++) {
                            if (j == usedGenPho[k]) used = 1;
                        }
                        if (deltaR(tmpGenPho.Phi(),
                                   tmpGenPho.Eta(),
                                   genLepNoFSR.v.Phi(),
                                   genLepNoFSR.v.Eta()) < 0.1 &&
                            !used) {
                            genLep.v += tmpGenPho;
                            usedGenPho.push_back(j);
                        }
                    }

                    if ((genLep.v.Pt() >= lepPtCutMin &&
                         fabs(genLep.v.Eta()) <= 0.1 * lepEtaCutMax && genLep.charge != 0) ||
                        ((lepSel == "SMu" || lepSel == "SE") && genLep.charge == 0)) {
                        genLeptons.push_back(genLep);
                    }
                }
            }

            ngenLeptons = genLeptons.size();

            // sort leptons by descending pt
            sort(genLeptons.begin(), genLeptons.end(), LepDescendingOrder);

            // assert(genLeptons.size() < 2 || genLeptons[0].v.Pt() >= genLeptons[1].v.Pt());
            if (genLeptons.size() > 1 && genLeptons[0].v.Pt() < genLeptons[1].v.Pt()) {
                std::cerr << "Problem in Gen jet ordering!\n";
            }

            if (countTauS3 == 0 && fileName.Index("UNFOLDING") >= 0 &&
                fileName.Index("Sherpa") < 0) {
                fill(partonsN, GNup - 5);
                fill(partonsNWeighted, GNup - 5, genWeight);
            }

            //--- if there are taus, but we do not run on the Tau file, thus we run on the DYJets
            // file,
            //    then we don't count the event at reco.
            if (countTauS3 > 0 && fileName.Index("Tau") < 0 && fileName.Index("Sherpa") < 0 &&
                fileName.Index("MG5") < 0) {
                passesTauCut = 0;
                std::cout << "Tau veto changed passesLeptonCut value!\n";
                passesLeptonCut = 0;
            }

            //-- determine if the event passes the leptons requirements
            if ((lepSel == "DMu" || lepSel == "DE") && ngenLeptons >= 2) {

                ++nGenEventsWithTwoGoodLeptons;
                nEffGenEventsWithTwoGoodLeptons += genWeight;

                // build the EWKBoson candidate and the kinematic
                genEWKBoson = genLeptons[0].v + genLeptons[1].v;

                if (genLeptons[0].charge * genLeptons[1].charge < 0) {
                    ++nGenEventsWithTwoGoodLeptonsWithOppCharge;
                    nEffGenEventsWithTwoGoodLeptonsWithOppCharge += genWeight;
                }

                // apply charge, mass and eta cut
                if (genLeptons[0].charge * genLeptons[1].charge < 0) {
                    passesgenLeptonCutNoMass = 1;
                    if (genEWKBoson.M() > ZMCutLow && genEWKBoson.M() < ZMCutHigh) {
                        ++nGenEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass;
                        nEffGenEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass += genWeight;
                        passesgenLeptonCut = 1;
                    }
                }
                //--- if there are taus we don't want the gen level (think about NoMass flag
                // here!!!)
                // if (countTauS3 > 0 && fileName.Index("Bugra") < 0 && fileName.Index("MG5") < 0)
                // passesgenLeptonCut = 0;
            } else if ((lepSel == "SMu" || lepSel == "SE") && (ngenLeptons >= 2)) {

                if (abs(genLeptons[0].charge) > 0 && genLeptons[1].charge == 0) {
                    genMET = genLeptons[1].v;
                } else if (abs(genLeptons[1].charge) > 0 && genLeptons[0].charge == 0) {
                    genMET = genLeptons[0].v;
                    genLeptons[0] = genLeptons[1];
                } else
                    passesgenLeptonCut = 0;

                genEWKBoson = genLeptons[0].v + genMET;
                genMT = sqrt(2 * genLeptons[0].v.Pt() * genMET.Pt() *
                             (1 - cos(genLeptons[0].v.Phi() - genMET.Phi())));

                // apply transverse mass and MET cut
                if (genMT > MTCutLow && genMET.Pt() > METCutLow) passesgenLeptonCut = 1;
                //--- if there are taus we don't want the gen level
                if (countTauS3 > 0) passesgenLeptonCut = 0;
            }
        } // end of hasGenInfo

        //=========== Rochester ==========
        //  cout << "---------" << "\n";
        if (doRochester && hasRecoInfo) {
            if ((lepSel == "DMu" || lepSel == "DE") && nLeptons >= 2) {
                // for each data muon in the loop, use this function to get a scale factor for its
                // momentum:
                if (EvtIsRealData) {
                    double dataSF_lead = 1.;
                    double dataSF_sublead = 1.;
                    dataSF_lead = rc->kScaleDT(leptons[0].charge,
                                               leptons[0].v.Pt(),
                                               leptons[0].v.Eta(),
                                               leptons[0].v.Phi(),
                                               0,
                                               0);
                    dataSF_sublead = rc->kScaleDT(leptons[1].charge,
                                                  leptons[1].v.Pt(),
                                                  leptons[1].v.Eta(),
                                                  leptons[1].v.Phi(),
                                                  0,
                                                  0);
                    // cout << "ptold " <<  leptons[0].v.Pt() << " , "  <<  leptons[1].v.Pt() <<
                    // "\n";
                    // double ptold0 = leptons[0].v.Pt();
                    // double ptold1 = leptons[1].v.Pt();
                    // cout <<  dataSF_lead << " , " << dataSF_sublead << "\n";
                    // cout << "pz old " << leptons[0].v.Pz() << "\n";

                    leptons[0].v.SetPtEtaPhiM(leptons[0].v.Pt() * dataSF_lead,
                                              leptons[0].v.Eta(),
                                              leptons[0].v.Phi(),
                                              .10565);
                    leptons[1].v.SetPtEtaPhiM(leptons[1].v.Pt() * dataSF_sublead,
                                              leptons[1].v.Eta(),
                                              leptons[1].v.Phi(),
                                              .10565);

                    // cout << "dataSF_lead " << dataSF_lead << " dataSF_sublead " << dataSF_sublead
                    // << "\n";
                    // cout << "pz new " << leptons[0].v.Pz() << "\n";
                    // cout << "ptnew " <<  leptons[0].v.Pt() << " , "  <<  leptons[1].v.Pt() <<
                    // "\n";
                    // cout << "M: " << leptons[0].v.E()*leptons[0].v.E()
                    // -(leptons[0].v.Px()*leptons[0].v.Px() + leptons[0].v.Py()*leptons[0].v.Py() +
                    // leptons[0].v.Pz()*leptons[0].v.Pz() ) << "\n";
                }

                if (!EvtIsRealData && hasGenInfo && ngenLeptons >= 2) {

                    // for MC, if matched gen-level muon (genPt) is available, use this function
                    // double u1= gRandom->Rndm();
                    // double u2= gRandom->Rndm();
                    double mcSF_lead = 1.;
                    double mcSF_sublead = 1.;

                    if (deltaR(genLeptons[0].v, leptons[0].v) < 0.1)
                        mcSF_lead = rc->kScaleFromGenMC(leptons[0].charge,
                                                        leptons[0].v.Pt(),
                                                        leptons[0].v.Eta(),
                                                        leptons[0].v.Phi(),
                                                        leptons[0].MuTkLayer,
                                                        genLeptons[0].v.Pt(),
                                                        gRandom->Rndm(),
                                                        0,
                                                        0);
                    else
                        mcSF_lead = rc->kScaleAndSmearMC(leptons[0].charge,
                                                         leptons[0].v.Pt(),
                                                         leptons[0].v.Eta(),
                                                         leptons[0].v.Phi(),
                                                         leptons[0].MuTkLayer,
                                                         gRandom->Rndm(),
                                                         gRandom->Rndm(),
                                                         0,
                                                         0);

                    if (deltaR(genLeptons[1].v, leptons[1].v) < 0.1)
                        mcSF_sublead = rc->kScaleFromGenMC(leptons[1].charge,
                                                           leptons[1].v.Pt(),
                                                           leptons[1].v.Eta(),
                                                           leptons[1].v.Phi(),
                                                           leptons[1].MuTkLayer,
                                                           genLeptons[1].v.Pt(),
                                                           gRandom->Rndm(),
                                                           0,
                                                           0);
                    else
                        mcSF_sublead = rc->kScaleAndSmearMC(leptons[1].charge,
                                                            leptons[1].v.Pt(),
                                                            leptons[1].v.Eta(),
                                                            leptons[1].v.Phi(),
                                                            leptons[1].MuTkLayer,
                                                            gRandom->Rndm(),
                                                            0,
                                                            0);

                    // cout << u1 << " , " << u2 << "\n";
                    // cout << "sf; "<< mcSF_lead << " , " << mcSF_sublead << "\n";
                    leptons[0].v.SetPtEtaPhiM(leptons[0].v.Pt() * mcSF_lead,
                                              leptons[0].v.Eta(),
                                              leptons[0].v.Phi(),
                                              .10565);
                    leptons[1].v.SetPtEtaPhiM(leptons[1].v.Pt() * mcSF_sublead,
                                              leptons[1].v.Eta(),
                                              leptons[1].v.Phi(),
                                              .10565);

                } // if(!EvtIsRealData

                // build Electroweak boson candidate: here it is expected to be a Z
                EWKBoson = leptons[0].v + leptons[1].v;
                // cout << EWKBoson.E()  << " , " << EWKBoson.Pt() << " , " << EWKBoson.Pz() <<
                // "\n";
                // cout << "M_new  "  << sqrt(EWKBoson.E()*EWKBoson.E() -
                // (EWKBoson.Pt()*EWKBoson.Pt()+ EWKBoson.Pz()*EWKBoson.Pz() )) << "\n";
                // cout << "M_new1  "  << sqrt(EWKBoson.E()*EWKBoson.E() - ( (leptons[0].v.Pt() +
                // leptons[1].v.Pt())* (leptons[0].v.Pt() + leptons[1].v.Pt()) +
                // EWKBoson.Pz()*EWKBoson.Pz() )) << "\n";

                // apply charge, mass cut
                if (passesTrigger && leptons[0].charge * leptons[1].charge < 0) {
                    passesLeptonChargeCut = 1;

                    if (leptons[0].v.Pt() > lepPtCutMin && leptons[1].v.Pt() > lepPtCutMin)
                        passesLeptonCutNoMass = 1;
                    else
                        passesLeptonCutNoMass = 0;

                    if (EWKBoson.M() > ZMCutLow && EWKBoson.M() < ZMCutHigh &&
                        leptons[0].v.Pt() > lepPtCutMin && leptons[1].v.Pt() > lepPtCutMin) {
                        // cout <<  "AndGoodMass" << "\n";
                        passesLeptonCut = 1;
                        // passesLeptonMassCut = 1;
                    } else {
                        passesLeptonCut = 0;
                    }
                } // end if re-selection of leptons

            } // if((lepSel == "DMu"
        }     // if (hasRecoInfo && doRochester)
              //============ end of Rochester ======================
        // cout << EWKBoson.M() << "\n";

        //=======================================================================================================//
        //   ------- lepton energy smearing ------
        //==========================================//
        if (hasRecoInfo && hasGenInfo) {
            if ((lepSel == "DMu" || lepSel == "DE") && nLeptons >= 2 && ngenLeptons >= 2) {
                if (smearlepton) {
                    double oldLeptonPt;
                    double genLeptonPt;
                    double newLeptonPt;
                    double smearFactor;

                    smearFactor = (lepSel == "DMu") ? 0.006 : 0.06;

                    oldLeptonPt = leptons[0].v.Pt();
                    genLeptonPt = genLeptons[0].v.Pt();
                    newLeptonPt = SmearLepPt(oldLeptonPt, genLeptonPt, smearlepton, smearFactor);
                    leptons[0].v.SetPtEtaPhiE(newLeptonPt,
                                              leptons[0].v.Eta(),
                                              leptons[0].v.Phi(),
                                              leptons[0].v.E() * newLeptonPt / oldLeptonPt);

                    oldLeptonPt = leptons[1].v.Pt();
                    genLeptonPt = genLeptons[1].v.Pt();
                    newLeptonPt = SmearLepPt(oldLeptonPt, genLeptonPt, smearlepton, smearFactor);
                    leptons[1].v.SetPtEtaPhiE(newLeptonPt,
                                              leptons[1].v.Eta(),
                                              leptons[1].v.Phi(),
                                              leptons[1].v.E() * newLeptonPt / oldLeptonPt);

                    // build Electroweak boson candidate: here it is expected to be a Z
                    EWKBoson = leptons[0].v + leptons[1].v;

                    // apply charge, mass cut
                    if (passesTrigger && leptons[0].charge * leptons[1].charge < 0) {
                        passesLeptonChargeCut = 1;

                        if (leptons[0].v.Pt() > lepPtCutMin && leptons[1].v.Pt() > lepPtCutMin)
                            passesLeptonCutNoMass = 1;
                        else
                            passesLeptonCutNoMass = 0;

                        if (EWKBoson.M() > ZMCutLow && EWKBoson.M() < ZMCutHigh &&
                            leptons[0].v.Pt() > lepPtCutMin && leptons[1].v.Pt() > lepPtCutMin) {
                            // cout <<  "AndGoodMass" << "\n";
                            passesLeptonCut = 1;
                            // passesLeptonMassCut = 1;
                        } else {
                            passesLeptonCut = 0;
                        }
                    } // end if re-selection of leptons
                }
            } // end if DMu/DE channel

        } // end of hasRecoInfo and hasGenInfo

        // cout << MuTkLayerCnt->at(0) << "\n";

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //          Retrieving jets           //
        //====================================//
        unsigned short nGoodJets(0), nGoodJets_20(0), nTotJets(0);
        bool bTagJetFound(false);
        double jetsHT(0);
        TLorentzVector hadronicR(0.0, 0.0, 0.0, 0.0);
        vector<jetStruct> jets;
        vector<jetStruct> jets_20; // additional jet collection with pt threshold of 20 GeV
        TLorentzVector jet1Plus2, jet1Minus2;

        if (hasRecoInfo) {
            nTotJets = JetAk04Eta->size();
            for (unsigned short i(0); i < nTotJets; i++) {
                bool passesBJets = false;
                if (fileName.Index("Sherpa") < 0)
                    passesBJets = (JetAk04BDiscCisvV2->at(i) >= 0.679);

                if (!EvtIsRealData && lepSel == "SMu") {
                    BTagModification(RandGen->Rndm(),
                                     JetAk04Pt->at(i),
                                     JetAk04Eta->at(i),
                                     JetAk04PartFlav->at(i),
                                     passesBJets);
                }

                jetStruct jet(JetAk04Pt->at(i),
                              JetAk04Eta->at(i),
                              JetAk04Phi->at(i),
                              JetAk04E->at(i),
                              i,
                              passesBJets);

                //-- apply jet energy scale uncertainty (need to change the scale when initiating
                // the object)
                double jetEnergyCorr = 0.;
                bool jetPassesPtCut(jet.v.Pt() >= 10);
                jetEnergyCorr = TableJESunc.getEfficiency(jet.v.Pt(), jet.v.Eta());

                jet.v.SetPtEtaPhiE(jet.v.Pt() * (1 + scale * jetEnergyCorr),
                                   jet.v.Eta(),
                                   jet.v.Phi(),
                                   jet.v.E() * (1 + scale * jetEnergyCorr));

                bool jetPassesEtaCut(fabs(jet.v.Rapidity()) <= 0.1 * jetEtaCutMax);
                bool jetPassesIdCut(JetAk04Id->at(i) > 0);
                bool jetPassesMVACut(JetAk04PuMva->at(i) > -0.2);

                bool jetPassesdRCut(1);
                unsigned short nRemovedLep = min(int(nLeptons), 2);
                for (unsigned short j(0); j < nRemovedLep; j++) {
                    if (deltaR(jet.v, leptons[j].v) < 0.4) {
                        jetPassesdRCut = 0;
                    }
                }

                // CAG
                if (passesLeptonCut && jet.v.Pt() >= 30 && jetPassesEtaCut && jetPassesIdCut &&
                    jetPassesdRCut) // no MVA cut
                    fill(puMVA, JetAk04PuMva->at(i), weight);

                if (jetPassesPtCut && jetPassesEtaCut && jetPassesIdCut && jetPassesMVACut &&
                    jetPassesdRCut) {
                    jets.push_back(jet);
                    // as soon as one selected jet is a b-jet, turn bTagJetFound to true
                    bTagJetFound = (bTagJetFound || passesBJets);
                }
            } // End of loop over all the jets

            nGoodJets = jets.size();

            // line below to test reco events that originate from TAU
            if (fileName.Index("Tau") >= 0 && countTauS3 == 0 && hasGenInfo) {
                passesTauCut = 0;
                std::cout << "Tau veto changes passesLeptonCut value!\n";
                passesLeptonCut = 0;
            }
        } // END IF HAS RECO
        //=======================================================================================================//

        // cout << "passesLeptonCut 2 : " << passesLeptonCut << "\n";
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //        Retrieving gen jets          //
        //====================================//
        unsigned short nGoodGenJets(0), nGoodGenJets_20(0), nTotGenJets(0);
        double genJetsHT(0);
        TLorentzVector genHadronicR(0.0, 0.0, 0.0, 0.0);
        vector<jetStruct> genJets;
        vector<TLorentzVector> genLVJets;
        vector<jetStruct> genJets_20;
        TLorentzVector genJet1Plus2, genJet1Minus2;

        if (hasGenInfo) {
            nTotGenJets = GJetAk04Eta->size();
            //-- retrieving generated jets
            for (unsigned short i(0); i < nTotGenJets; i++) {
                jetStruct genJet(GJetAk04Pt->at(i),
                                 GJetAk04Eta->at(i),
                                 GJetAk04Phi->at(i),
                                 GJetAk04E->at(i),
                                 i,
                                 0);
                bool genJetPassesdRCut(1);
                for (unsigned short j(0); j < ngenLeptons; j++) {
                    if (deltaR(genJet.v, genLeptons[j].v) < 0.4) {
                        genJetPassesdRCut = 0;
                    }
                }
                if (genJet.v.Pt() >= 10 && fabs(genJet.v.Eta()) <= 5.0 && genJetPassesdRCut) {
                    genJets.push_back(genJet);
                }
            }
            nGoodGenJets = genJets.size();
        }
        //=======================================================================================================//

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //     Matching gen and reco jets     //
        //====================================//
        vector<int> genJetsIndex(nGoodGenJets, 0);
        vector<vector<int>> matchingTable(nGoodJets, genJetsIndex);
        if (hasRecoInfo && hasGenInfo) {
            for (unsigned short i(0); i < nGoodJets; i++) {
                double mindR(0.15);
                int index(-1);
                double dR(9999);
                for (unsigned short j(0); j < nGoodGenJets; j++) {
                    dR = deltaR(genJets[j].v, jets[i].v);
                    if (dR < mindR) {
                        mindR = dR;
                        index = j;
                    }
                }
                if (index > -1) {
                    matchingTable[i][index] = 1;
                    double oldJetPt = jets[i].v.Pt();
                    double newJetPt =
                        SmearJetPt(oldJetPt, genJets[index].v.Pt(), jets[i].v.Eta(), smearJet);
                    jets[i].v.SetPtEtaPhiE(newJetPt,
                                           jets[i].v.Eta(),
                                           jets[i].v.Phi(),
                                           jets[i].v.E() * newJetPt / oldJetPt);
                }
            }

            /*
                        //-- print the mathcing table
                        //cout << "\n mathcing Table: \n" << endl;
                        for (int i = 0; i < int(matchingTable.size()); i++){
                          for (int j = 0; j < int(matchingTable[i].size()); j++){
                            cout << matchingTable[i][j] << "  ";
                          }
                          cout << endl;
                        }
            */
        }

        //=======================================================================================================//
        // Re-analyze the jets collections and cut on the Pt    //
        // we can do it only now since we needed to smear      //
        // the jet pt distribution for the MC                 //
        //===================================================//

        if (hasRecoInfo) {
            // subsamble the prior jet collection by applying pt cut
            vector<jetStruct> tmpJets;
            for (unsigned short i(0); i < nGoodJets; i++) {
                if (jets[i].v.Pt() >= jetPtCutMin) tmpJets.push_back(jets[i]);
                if (jets[i].v.Pt() >= 20) jets_20.push_back(jets[i]);
            }
            jets.clear();
            jets = tmpJets;
            nGoodJets = jets.size();
            // if(fileName.Index("TT") >= 0) cout <<
            // "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQTTbar"<< nGoodJets << " SF:
            // "<<TTbarSF.getTTbarSF(nGoodJets) << "weight" << weight << endl;
            // if(fileName.Index("TT") >= 0 && (systematics == 3)) {
            //	if(direction > 0 ) weight /= TTbarSF.getTTbarSFHigh(nGoodJets);
            //	else if(direction < 0 ) weight /= TTbarSF.getTTbarSFLow(nGoodJets);
            // }
            // else if(fileName.Index("TT") >= 0) weight /= TTbarSF.getTTbarSF(nGoodJets);
            nGoodJets_20 = jets_20.size();
            sort(jets.begin(), jets.end(), JetDescendingOrder);
            sort(jets_20.begin(), jets_20.end(), JetDescendingOrder);
            if (nGoodJets >= 2) {
                jet1Plus2 = jets[0].v + jets[1].v;
                jet1Minus2 = jets[0].v - jets[1].v;
            }
            jetsHT = 0;
            for (unsigned short i(0); i < nGoodJets; i++) {
                jetsHT += jets[i].v.Pt();
                if (nGoodJets >= 1) hadronicR += jets[i].v;
            }
        }

        if (hasGenInfo) {
            vector<jetStruct> tmpGenJets;
            for (unsigned short i(0); i < nGoodGenJets; i++) {
                if (genJets[i].v.Pt() >= jetPtCutMin &&
                    fabs(genJets[i].v.Rapidity()) <= 0.1 * jetEtaCutMax)
                    tmpGenJets.push_back(genJets[i]);
                if (genJets[i].v.Pt() >= 20 && fabs(genJets[i].v.Rapidity()) <= 0.1 * jetEtaCutMax)
                    genJets_20.push_back(genJets[i]);
            }
            genJets.clear();
            genJets = tmpGenJets;
            nGoodGenJets = genJets.size();
            nGoodGenJets_20 = genJets_20.size();
            sort(genJets.begin(), genJets.end(), JetDescendingOrder);
            if (nGoodGenJets >= 2) {
                genJet1Plus2 = genJets[0].v + genJets[1].v;
                genJet1Minus2 = genJets[0].v - genJets[1].v;
            }
            genJetsHT = 0.;
            for (unsigned short i(0); i < nGoodGenJets; i++) {
                genJetsHT += genJets[i].v.Pt();
                if (nGoodGenJets >= 1) {
                    genHadronicR += genJets[i].v;
                    //   if(genEWKBoson.Pt()>0. && genEWKBoson.Pt()> 50. ) cout << "px = " <<
                    //   genJets[i].v.Px() << " , py = " <<  genJets[i].v.Py() << " , pt = " <<
                    //   genJets[i].v.Pt() << " , eta = " <<  genJets[i].v.Rapidity() << "\n";
                }
            }

            // if(nGoodGenJets>=1 && genEWKBoson.Pt()>0. && genEWKBoson.Pt()> 50. ){
            //    cout << "Zpx = " << genEWKBoson.Px()  <<  " , Zpy = " << genEWKBoson.Py() <<  " ,
            //    Zpt = " << genEWKBoson.Pt() << "\n";
            //    cout << "hadronicR.Pt = " << genHadronicR.Pt()  << ", JZB = " <<
            //    genHadronicR.Pt()-genEWKBoson.Pt() << "\n";
            // }
            // cout << "------------------------" << "\n";

            sort(genJets_20.begin(), genJets_20.end(), JetDescendingOrder);
        }
        //=======================================================================================================//

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        // Select the best pair of jets for DPS  //
        //=======================================//
        pair<TLorentzVector, TLorentzVector> bestTwoJets;
        TLorentzVector bestJet1Plus2, bestJet1Minus2;
        if (hasRecoInfo) {
            bestTwoJetsCandidatesPt(jets, bestTwoJets);
            // bestTwoJetsCandidatesPhi(jets, bestTwoJets);
            if (nGoodJets >= 2) {
                bestJet1Plus2 = bestTwoJets.first + bestTwoJets.second;
                bestJet1Minus2 = bestTwoJets.first - bestTwoJets.second;
            }
        }

        pair<TLorentzVector, TLorentzVector> genBestTwoJets;
        TLorentzVector genBestJet1Plus2, genBestJet1Minus2;
        if (hasGenInfo) {
            bestTwoJetsCandidatesPt(genJets, genBestTwoJets);
            // bestTwoJetsCandidatesPhi(genJets, genBestTwoJets);
            if (nGoodGenJets >= 2) {
                genBestJet1Plus2 = genBestTwoJets.first + genBestTwoJets.second;
                genBestJet1Minus2 = genBestTwoJets.first - genBestTwoJets.second;
            }
        }
        //=======================================================================================================//

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //   Filling gen end parton histos    //
        //====================================//

        double gentau_sum(0), gentau_max(0);
        double gentau_c_sum(0), gentau_c_max(0);
        double gentau_cm_sum(0), gentau_cm_max(0);
        double gentau_c_cm_sum(0), gentau_c_cm_max(0);

        // cout << "passesgenLeptonCut : " << passesgenLeptonCut << endl;

        if (hasGenInfo) {

            if (passesgenLeptonCutNoMass) {
                if (genEWKBoson.M() > 115. && genEWKBoson.M() < 135.)
                    fill(genZPt_Zinc0jetM115_135, genEWKBoson.Pt(), commonGenWeight, EvtWeights);
                if (genEWKBoson.M() >= 111. && genEWKBoson.M() < 130.)
                    fill(genZPt_Zinc0jetM111_130, genEWKBoson.Pt(), commonGenWeight, EvtWeights);
                if (genEWKBoson.M() >= 130. && genEWKBoson.M() <= 3000.)
                    fill(genZPt_Zinc0jetM130_3, genEWKBoson.Pt(), commonGenWeight, EvtWeights);
            }

            if (passesgenLeptonCut) {
                // cout << "Selected at gen level" << endl;
                nGenEventsVInc0Jets++;
                nEffGenEventsVInc0Jets += genWeight;

                double RatioValue = 1.;
                if (UnfoldUnc) {
                    int binNumber = ZNGoodJets_ZexcHratio_fit->GetXaxis()->FindBin(nGoodJets);
                    RatioValue = ZNGoodJets_ZexcHratio_fit->GetBinContent(binNumber);
                }

                fill(genZNGoodJets_Zexc, nGoodGenJets, commonGenWeight * RatioValue, EvtWeights);
                fill(genZNGoodJets_Zinc, 0., commonGenWeight, EvtWeights);
                fill(genZMass_Zinc0jet, genEWKBoson.M(), commonGenWeight, EvtWeights);
                fill(genZPt_Zinc0jet, genEWKBoson.Pt(), commonGenWeight, EvtWeights);
                fill(genZRapidity_Zinc0jet, genEWKBoson.Rapidity(), commonGenWeight, EvtWeights);
                fill(genZEta_Zinc0jet, genEWKBoson.Eta(), commonGenWeight, EvtWeights);
                fill(genlepPt_Zinc0jet, genLeptons[0].v.Pt(), commonGenWeight, EvtWeights);
                fill(genlepPt_Zinc0jet, genLeptons[1].v.Pt(), commonGenWeight, EvtWeights);
                fill(genlepEta_Zinc0jet, genLeptons[0].v.Eta(), commonGenWeight, EvtWeights);
                fill(genlepEta_Zinc0jet, genLeptons[1].v.Eta(), commonGenWeight, EvtWeights);
                fill(genVisPt_Zinc0jetQun, genEWKBoson.Pt(), commonGenWeight, EvtWeights);

                if (nGoodGenJets_20 >= 1) {

                    double RatioValue = 1.;
                    if (nGoodJets_20 >= 1 && UnfoldUnc) {
                        double binNumber =
                            FirstJetPt_2_Zinc1jetHratio_fit->GetXaxis()->FindBin(jets_20[0].v.Pt());
                        RatioValue = FirstJetPt_2_Zinc1jetHratio_fit->GetBinContent(binNumber);
                    }

                    fill(genFirstJetPt_Zinc1jet,
                         genJets_20[0].v.Pt(),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                    fill(genFirstJetPtEta_Zinc1jet,
                         genJets_20[0].v.Pt(),
                         fabs(genJets[0].v.Eta()),
                         genWeight);
                }
                if (nGoodGenJets >= 1) {
                    nGenEventsVInc1Jets++;
                    nEffGenEventsVInc1Jets += genWeight;
                    fill(genAbsFirstJetRapidity_Zinc1jet,
                         fabs(genJets[0].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSumZFirstJetRapidity_Zinc1jet,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);
                    fill(genDifZFirstJetRapidity_Zinc1jet,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);
                    // cross check//////
                    fill(genSumZFirstJetEta_Zinc1jet,
                         fabs(genEWKBoson.Eta() + genJets[0].v.Eta()) / 2.0,
                         commonGenWeight,
                         EvtWeights);
                    fill(genDifZFirstJetEta_Zinc1jet,
                         fabs(genEWKBoson.Eta() - genJets[0].v.Eta()) / 2.0,
                         commonGenWeight,
                         EvtWeights);

                    /// Azimuth cross check//////////////////////////
                    fill(genDPhiZFirstJet_Zinc1jet,
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         commonGenWeight,
                         EvtWeights);

                    if (genEWKBoson.Pt() > 100.) {
                        fill(genAbsZRapidity_ZPt100_Zinc1jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsFirstJetRapidity_ZPt100_Zinc1jet,
                             fabs(genJets[0].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZFirstJetRapidity_ZPt100_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZFirstJetRapidity_ZPt100_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                    }

                    if (genEWKBoson.Pt() > 150.) {
                        fill(genAbsZRapidity_ZPt150_Zinc1jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsFirstJetRapidity_ZPt150_Zinc1jet,
                             fabs(genJets[0].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZFirstJetRapidity_ZPt150_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZFirstJetRapidity_ZPt150_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);

                        fill(genDPhiZFirstJet_ZPt150_Zinc1jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                             commonGenWeight,
                             EvtWeights);
                    }

                    if (genEWKBoson.Pt() > 300.) {
                        fill(genAbsZRapidity_ZPt300_Zinc1jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsFirstJetRapidity_ZPt300_Zinc1jet,
                             fabs(genJets[0].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZFirstJetRapidity_ZPt300_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZFirstJetRapidity_ZPt300_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);

                        fill(genDPhiZFirstJet_ZPt300_Zinc1jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                             commonGenWeight,
                             EvtWeights);
                    }

                    /// different JetPt cuts///////
                    if (genJets[0].v.Pt() > 50.) {
                        fill(genAbsZRapidity_FirstJetPt50_Zinc1jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsFirstJetRapidity_FirstJetPt50_Zinc1jet,
                             fabs(genJets[0].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZFirstJetRapidity_FirstJetPt50_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZFirstJetRapidity_FirstJetPt50_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                    }

                    if (genJets[0].v.Pt() > 80.) {
                        fill(genAbsZRapidity_FirstJetPt80_Zinc1jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsFirstJetRapidity_FirstJetPt80_Zinc1jet,
                             fabs(genJets[0].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZFirstJetRapidity_FirstJetPt80_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZFirstJetRapidity_FirstJetPt80_Zinc1jet,
                             fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                    }

                    for (unsigned short i(0); i < nGoodGenJets; i++) {
                        double trans_mass = genJets[i].v.Mt();

                        gentau_sum += trans_mass *
                                      exp(-fabs(genJets[i].v.Rapidity() - genEWKBoson.Rapidity()));
                        gentau_max =
                            max(gentau_max,
                                trans_mass *
                                    exp(-fabs(genJets[i].v.Rapidity() - genEWKBoson.Rapidity())));

                        gentau_c_sum +=
                            trans_mass /
                            (2 * cosh(genJets[i].v.Rapidity() - genEWKBoson.Rapidity()));
                        gentau_c_max =
                            max(gentau_c_max,
                                trans_mass /
                                    (2 * cosh(genJets[i].v.Rapidity() - genEWKBoson.Rapidity())));

                        gentau_cm_sum += trans_mass * exp(-fabs(genJets[i].v.Rapidity()));
                        gentau_cm_max =
                            max(gentau_cm_max, trans_mass * exp(-fabs(genJets[i].v.Rapidity())));

                        gentau_c_cm_sum += trans_mass / (2 * cosh(genJets[i].v.Rapidity()));
                        gentau_c_cm_max =
                            max(gentau_c_cm_max, trans_mass / (2 * cosh(genJets[i].v.Rapidity())));
                    }

                    for (unsigned short i(0); i < 5; i++) {
                        if (genEWKBoson.Pt() > ZptRange[i] && genEWKBoson.Pt() <= ZptRange[i + 1]) {
                            fill(gentau_sum_Zinc1jet[i], gentau_sum, commonGenWeight, EvtWeights);
                            fill(gentau_max_Zinc1jet[i], gentau_max, commonGenWeight, EvtWeights);
                            fill(gentau_c_sum_Zinc1jet[i],
                                 gentau_c_sum,
                                 commonGenWeight,
                                 EvtWeights);
                            fill(gentau_c_max_Zinc1jet[i],
                                 gentau_c_max,
                                 commonGenWeight,
                                 EvtWeights);
                            fill(gentau_cm_sum_Zinc1jet[i],
                                 gentau_cm_sum,
                                 commonGenWeight,
                                 EvtWeights);
                            fill(gentau_cm_max_Zinc1jet[i],
                                 gentau_cm_max,
                                 commonGenWeight,
                                 EvtWeights);
                            fill(gentau_c_cm_sum_Zinc1jet[i],
                                 gentau_c_cm_sum,
                                 commonGenWeight,
                                 EvtWeights);
                            fill(gentau_c_cm_max_Zinc1jet[i],
                                 gentau_c_cm_max,
                                 commonGenWeight,
                                 EvtWeights);
                        }
                    }

                    fill(genZNGoodJets_Zinc, 1., commonGenWeight, EvtWeights);
                    fill(genZPt_Zinc1jet, genEWKBoson.Pt(), commonGenWeight, EvtWeights);
                    fill(
                        genZRapidity_Zinc1jet, genEWKBoson.Rapidity(), commonGenWeight, EvtWeights);
                    fill(genZAbsRapidity_Zinc1jet,
                         fabs(genEWKBoson.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genZEta_Zinc1jet, genEWKBoson.Eta(), commonGenWeight, EvtWeights);

                    double RatioValue = 1.;
                    double RatioValue1 = 1.;
                    double RatioValue2 = 1.;
                    if (nGoodJets >= 1 && UnfoldUnc) {
                        double binNumber =
                            FirstJetAbsRapidity_2_Zinc1jetHratio_fit->GetXaxis()->FindBin(
                                fabs(jets[0].v.Eta()));
                        RatioValue =
                            FirstJetAbsRapidity_2_Zinc1jetHratio_fit->GetBinContent(binNumber);
                        double binNumber1 =
                            JetsHT_2_Zinc1jetHratio_fit->GetXaxis()->FindBin(jetsHT);
                        RatioValue1 = JetsHT_2_Zinc1jetHratio_fit->GetBinContent(binNumber1);
                        double binNumber2 = VisPt_2_Zinc1jetQunHratio_fit->GetXaxis()->FindBin(
                            fabs((jets[0].v + EWKBoson).Pt()));
                        RatioValue2 = VisPt_2_Zinc1jetQunHratio_fit->GetBinContent(binNumber2);
                    }
                    // cout << RatioValue1 << " , " << RatioValue1 << "\n";

                    fill(genFirstJetEta_Zinc1jet,
                         fabs(genJets[0].v.Eta()),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                    fill(genAbsZRapidity_Zinc1jet,
                         fabs(genEWKBoson.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genFirstJetAbsRapidity_Zinc1jet,
                         fabs(genJets[0].v.Rapidity()),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                    fill(genFirstJetEtaHigh_Zinc1jet,
                         fabs(genJets[0].v.Eta()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genFirstJetRapidityHigh_Zinc1jet,
                         fabs(genJets[0].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genJetsHT_Zinc1jet, genJetsHT, commonGenWeight * RatioValue1, EvtWeights);
                    // fill(genJetsHT_2_Zinc1jet, genJetsHT, commonGenWeight, EvtWeights);
                    fill(genVisPt_Zinc1jetQun,
                         fabs((genJets[0].v + genEWKBoson).Pt()),
                         commonGenWeight * RatioValue2,
                         EvtWeights);
                    fill(genSumZJetRapidity_Zinc1jet,
                         0.5 * fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genDifZJetRapidity_Zinc1jet,
                         0.5 * fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    if (nGoodGenJets == 1) {
                        fill(
                            genFirstJetPt_Zexc1jet, genJets[0].v.Pt(), commonGenWeight, EvtWeights);
                        // Additional Branch
                        fill(genAbsZRapidity_Zexc1jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsJetRapidity_Zexc1jet,
                             fabs(genJets[0].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZJetRapidity_Zexc1jet,
                             fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZJetRapidity_Zexc1jet,
                             fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);

                        if (genEWKBoson.Pt() > 100.) {
                            fill(genAbsZRapidity_ZPt100_Zexc1jet,
                                 fabs(genEWKBoson.Rapidity()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genAbsJetRapidity_ZPt100_Zexc1jet,
                                 fabs(genJets[0].v.Rapidity()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSumZJetRapidity_ZPt100_Zexc1jet,
                                 fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genDifZJetRapidity_ZPt100_Zexc1jet,
                                 fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                                 commonGenWeight,
                                 EvtWeights);
                        }

                        if (genEWKBoson.Pt() > 150.) {
                            fill(genAbsZRapidity_ZPt150_Zexc1jet,
                                 fabs(genEWKBoson.Rapidity()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genAbsJetRapidity_ZPt150_Zexc1jet,
                                 fabs(genJets[0].v.Rapidity()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSumZJetRapidity_ZPt150_Zexc1jet,
                                 fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genDifZJetRapidity_ZPt150_Zexc1jet,
                                 fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                                 commonGenWeight,
                                 EvtWeights);
                        }
                    }
                }
                if (nGoodGenJets_20 >= 2) {

                    double RatioValue = 1.;
                    if (nGoodJets_20 >= 2 && UnfoldUnc) {
                        double binNumber = SecondJetPt_2_Zinc2jetHratio_fit->GetXaxis()->FindBin(
                            jets_20[1].v.Pt());
                        RatioValue = SecondJetPt_2_Zinc2jetHratio_fit->GetBinContent(binNumber);
                    }

                    fill(genSecondJetPt_Zinc2jet,
                         genJets_20[1].v.Pt(),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                }
                if (nGoodGenJets >= 2) {
                    TLorentzVector genJet1Plus2PlusZ = genJet1Plus2 + genEWKBoson;
                    nGenEventsVInc2Jets++;
                    nEffGenEventsVInc2Jets += genWeight;
                    ///////////Special Branch//////////////////
                    fill(genAbsFirstJetRapidity_Zinc2jet,
                         fabs(genJets[0].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSumZFirstJetRapidity_Zinc2jet,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);
                    fill(genDifZFirstJetRapidity_Zinc2jet,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);

                    fill(genAbsZRapidity_Zinc2jet,
                         fabs(genEWKBoson.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genAbsSecondJetRapidity_Zinc2jet,
                         fabs(genJets[1].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSumZSecondJetRapidity_Zinc2jet,
                         fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);
                    fill(genDifZSecondJetRapidity_Zinc2jet,
                         fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);

                    fill(genSumFirstSecondJetRapidity_Zinc2jet,
                         fabs(genJets[0].v.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);
                    fill(genDifFirstSecondJetRapidity_Zinc2jet,
                         fabs(genJets[0].v.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);

                    TLorentzVector genDiJets = genJets[0].v + genJets[1].v;
                    fill(genSumZTwoJetsRapidity_Zinc2jet,
                         fabs(genEWKBoson.Rapidity() + genDiJets.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);
                    fill(genDifZTwoJetsRapidity_Zinc2jet,
                         fabs(genEWKBoson.Rapidity() - genDiJets.Rapidity()) / 2.0,
                         commonGenWeight,
                         EvtWeights);

                    /////Azimuth cross check//////////////////////////////////
                    fill(genDPhiZFirstJet_Zinc2jet,
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         commonGenWeight,
                         EvtWeights);
                    fill(genDPhiZSecondJet_Zinc2jet,
                         fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                         commonGenWeight,
                         EvtWeights);
                    fill(genDPhiFirstSecondJet_Zinc2jet,
                         fabs(genJets[0].v.DeltaPhi(genJets[1].v)),
                         commonGenWeight,
                         EvtWeights);

                    if (genEWKBoson.Pt() > 100.) {
                        fill(genAbsZRapidity_ZPt100_Zinc2jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsSecondJetRapidity_ZPt100_Zinc2jet,
                             fabs(genJets[1].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZSecondJetRapidity_ZPt100_Zinc2jet,
                             fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZSecondJetRapidity_ZPt100_Zinc2jet,
                             fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                    }

                    if (genEWKBoson.Pt() > 150.) {
                        fill(genAbsZRapidity_ZPt150_Zinc2jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsSecondJetRapidity_ZPt150_Zinc2jet,
                             fabs(genJets[1].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZSecondJetRapidity_ZPt150_Zinc2jet,
                             fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZSecondJetRapidity_ZPt150_Zinc2jet,
                             fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);

                        fill(genDPhiZFirstJet_ZPt150_Zinc2jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                             commonGenWeight,
                             EvtWeights);
                    }

                    if (genEWKBoson.Pt() > 300.) {
                        fill(genDPhiZFirstJet_ZPt300_Zinc2jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                             commonGenWeight,
                             EvtWeights);
                    }

                    // Set Jet rapidity discriminator/////

                    if (fabs(genJets[0].v.Rapidity() - genJets[1].v.Rapidity()) > 2) {
                        fill(genAbsZRapidity_DifJetRapidityl2_Zinc2jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet,
                             fabs(genJets[0].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet,
                             fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet,
                             fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                    }

                    if (fabs(genJets[0].v.Rapidity() - genJets[1].v.Rapidity()) < 2) {
                        fill(genAbsZRapidity_DifJetRapiditys2_Zinc2jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet,
                             fabs(genJets[0].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet,
                             fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet,
                             fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                    }

                    fill(genZNGoodJets_Zinc, 2., commonGenWeight, EvtWeights);
                    fill(
                        genTwoJetsPtDiff_Zinc2jet, genJet1Minus2.Pt(), commonGenWeight, EvtWeights);
                    fill(genBestTwoJetsPtDiff_Zinc2jet,
                         genBestJet1Minus2.Pt(),
                         commonGenWeight,
                         EvtWeights);
                    fill(genJetsMass_Zinc2jet, genJet1Plus2.M(), commonGenWeight, EvtWeights);
                    fill(
                        genllJetsMass_Zinc2jet, genJet1Plus2PlusZ.M(), commonGenWeight, EvtWeights);

                    if (genJet1Plus2PlusZ.M() > 450 && genJet1Plus2PlusZ.M() < 600) {
                        if (fabs(genJets[0].v.Eta()) < fabs(genJets[1].v.Eta())) {
                            fill(genCentralJetEta_Zinc2jet,
                                 fabs(genJets[0].v.Eta()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genForwardJetEta_Zinc2jet,
                                 fabs(genJets[1].v.Eta()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genCentralJetPt_Zinc2jet,
                                 genJets[0].v.Pt(),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genForwardJetPt_Zinc2jet,
                                 genJets[1].v.Pt(),
                                 commonGenWeight,
                                 EvtWeights);
                        } else {
                            fill(genCentralJetEta_Zinc2jet,
                                 fabs(genJets[1].v.Eta()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genForwardJetEta_Zinc2jet,
                                 fabs(genJets[0].v.Eta()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genCentralJetPt_Zinc2jet,
                                 genJets[1].v.Pt(),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genForwardJetPt_Zinc2jet,
                                 genJets[0].v.Pt(),
                                 commonGenWeight,
                                 EvtWeights);
                        }
                    }

                    if (EvtVtxCnt < 14)
                        fill(genJetsMassLowPU_Zinc2jet,
                             genJet1Plus2.M(),
                             commonGenWeight,
                             EvtWeights);
                    else if (EvtVtxCnt < 18)
                        fill(genJetsMassMidPU_Zinc2jet,
                             genJet1Plus2.M(),
                             commonGenWeight,
                             EvtWeights);
                    else
                        fill(genJetsMassHigPU_Zinc2jet,
                             genJet1Plus2.M(),
                             commonGenWeight,
                             EvtWeights);
                    fill(genZPt_Zinc2jet, genEWKBoson.Pt(), commonGenWeight, EvtWeights);
                    fill(
                        genZRapidity_Zinc2jet, genEWKBoson.Rapidity(), commonGenWeight, EvtWeights);
                    fill(genZEta_Zinc2jet, genEWKBoson.Eta(), commonGenWeight, EvtWeights);

                    double RatioValue = 1.;
                    double RatioValue1 = 1.;
                    double RatioValue2 = 1.;
                    if (nGoodJets >= 2 && UnfoldUnc) {
                        double binNumber =
                            SecondJetAbsRapidity_2_Zinc2jetHratio_fit->GetXaxis()->FindBin(
                                fabs(jets[1].v.Eta()));
                        RatioValue =
                            SecondJetAbsRapidity_2_Zinc2jetHratio_fit->GetBinContent(binNumber);
                        double binNumber1 =
                            JetsHT_2_Zinc2jetHratio_fit->GetXaxis()->FindBin(jetsHT);
                        RatioValue1 = JetsHT_2_Zinc2jetHratio_fit->GetBinContent(binNumber1);
                        double binNumber2 = VisPt_2_Zinc2jetQunHratio_fit->GetXaxis()->FindBin(
                            fabs((jets[0].v + jets[1].v + EWKBoson).Pt()));
                        RatioValue2 = VisPt_2_Zinc2jetQunHratio_fit->GetBinContent(binNumber2);
                    }

                    fill(genSecondJetEta_Zinc2jet,
                         fabs(genJets[1].v.Eta()),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                    fill(genSecondJetAbsRapidity_Zinc2jet,
                         fabs(genJets[1].v.Rapidity()),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                    fill(genSecondJetEtaHigh_Zinc2jet,
                         fabs(genJets[1].v.Eta()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSecondJetRapidityHigh_Zinc2jet,
                         fabs(genJets[1].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genJetsHT_Zinc2jet, genJetsHT, commonGenWeight * RatioValue1, EvtWeights);
                    fill(genVisPt_Zinc2jetQun,
                         fabs((genJets[0].v + genJets[1].v + genEWKBoson).Pt()),
                         commonGenWeight * RatioValue2,
                         EvtWeights);
                    fill(genptBal_Zinc2jet, genJet1Plus2PlusZ.Pt(), commonGenWeight, EvtWeights);
                    fill(gendPhiJets_Zinc2jet,
                         deltaPhi(genJets[0].v, genJets[1].v),
                         commonGenWeight,
                         EvtWeights);
                    fill(genBestdPhiJets_Zinc2jet,
                         deltaPhi(genBestTwoJets.first, genBestTwoJets.second),
                         commonGenWeight,
                         EvtWeights);
                    fill(gendEtaJets_Zinc2jet,
                         genJets[0].v.Eta() - genJets[1].v.Eta(),
                         commonGenWeight,
                         EvtWeights);
                    fill(gendEtaFirstJetZ_Zinc2jet,
                         genJets[0].v.Eta() - genEWKBoson.Eta(),
                         commonGenWeight,
                         EvtWeights);
                    fill(gendEtaSecondJetZ_Zinc2jet,
                         genJets[1].v.Eta() - genEWKBoson.Eta(),
                         commonGenWeight,
                         EvtWeights);
                    fill(gendEtaJet1Plus2Z_Zinc2jet,
                         genJet1Plus2.Eta() - genEWKBoson.Eta(),
                         commonGenWeight,
                         EvtWeights);
                    fill(genPHI_Zinc2jet,
                         PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                         commonGenWeight,
                         EvtWeights);
                    fill(genBestPHI_Zinc2jet,
                         PHI(genLeptons[0].v,
                             genLeptons[1].v,
                             genBestTwoJets.first,
                             genBestTwoJets.second),
                         commonGenWeight,
                         EvtWeights);
                    fill(genPHI_T_Zinc2jet,
                         PHI_T(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                         commonGenWeight,
                         EvtWeights);
                    fill(genBestPHI_T_Zinc2jet,
                         PHI_T(genLeptons[0].v,
                               genLeptons[1].v,
                               genBestTwoJets.first,
                               genBestTwoJets.second),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSpT_Zinc2jet,
                         SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                         commonGenWeight,
                         EvtWeights);
                    fill(genBestSpT_Zinc2jet,
                         SpT(genLeptons[0].v,
                             genLeptons[1].v,
                             genBestTwoJets.first,
                             genBestTwoJets.second),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSpTJets_Zinc2jet,
                         SpTsub(genJets[0].v, genJets[1].v),
                         commonGenWeight,
                         EvtWeights);
                    fill(genBestSpTJets_Zinc2jet,
                         SpTsub(genBestTwoJets.first, genBestTwoJets.second),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSpTLeptons_Zinc2jet,
                         SpTsub(genLeptons[0].v, genLeptons[1].v),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSPhi_Zinc2jet,
                         SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                         commonGenWeight,
                         EvtWeights);
                    fill(genBestSPhi_Zinc2jet,
                         SPhi(genLeptons[0].v,
                              genLeptons[1].v,
                              genBestTwoJets.first,
                              genBestTwoJets.second),
                         commonGenWeight,
                         EvtWeights);

                    if (genEWKBoson.Pt() < 25) {
                        fill(genptBal_LowPt_Zinc2jet,
                             genJet1Plus2PlusZ.Pt(),
                             commonGenWeight,
                             EvtWeights);
                        fill(gendPhiJets_LowPt_Zinc2jet,
                             deltaPhi(genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genBestdPhiJets_LowPt_Zinc2jet,
                             deltaPhi(genBestTwoJets.first, genBestTwoJets.second),
                             commonGenWeight,
                             EvtWeights);
                        fill(genPHI_T_LowPt_Zinc2jet,
                             PHI_T(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genBestPHI_T_LowPt_Zinc2jet,
                             PHI_T(genLeptons[0].v,
                                   genLeptons[1].v,
                                   genBestTwoJets.first,
                                   genBestTwoJets.second),
                             commonGenWeight,
                             EvtWeights);
                        fill(genPHI_LowPt_Zinc2jet,
                             PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genBestPHI_LowPt_Zinc2jet,
                             PHI(genLeptons[0].v,
                                 genLeptons[1].v,
                                 genBestTwoJets.first,
                                 genBestTwoJets.second),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSpTJets_LowPt_Zinc2jet,
                             SpTsub(genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genBestSpTJets_LowPt_Zinc2jet,
                             SpTsub(genBestTwoJets.first, genBestTwoJets.second),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSpTLeptons_LowPt_Zinc2jet,
                             SpTsub(genLeptons[0].v, genLeptons[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSpT_LowPt_Zinc2jet,
                             SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genBestSpT_LowPt_Zinc2jet,
                             SpT(genLeptons[0].v,
                                 genLeptons[1].v,
                                 genBestTwoJets.first,
                                 genBestTwoJets.second),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSPhi_LowPt_Zinc2jet,
                             SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genBestSPhi_LowPt_Zinc2jet,
                             SPhi(genLeptons[0].v,
                                  genLeptons[1].v,
                                  genBestTwoJets.first,
                                  genBestTwoJets.second),
                             commonGenWeight,
                             EvtWeights);
                        if (SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v) <
                            0.5) {
                            fill(genPHI_LowSpT_LowPt_Zinc2jet,
                                 PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSPhi_LowSpT_LowPt_Zinc2jet,
                                 SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                        } else {
                            fill(genPHI_HighSpT_LowPt_Zinc2jet,
                                 PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSPhi_HighSpT_LowPt_Zinc2jet,
                                 SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                        }
                        if (SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v) <
                            0.5) {
                            fill(genSpT_LowSPhi_LowPt_Zinc2jet,
                                 SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                        } else {
                            fill(genSpT_HighSPhi_LowPt_Zinc2jet,
                                 SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                        }
                    } else {
                        fill(genptBal_HighPt_Zinc2jet,
                             genJet1Plus2PlusZ.Pt(),
                             commonGenWeight,
                             EvtWeights);
                        fill(gendPhiJets_HighPt_Zinc2jet,
                             deltaPhi(genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genPHI_HighPt_Zinc2jet,
                             PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genPHI_T_HighPt_Zinc2jet,
                             PHI_T(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSpTJets_HighPt_Zinc2jet,
                             SpTsub(genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSpTLeptons_HighPt_Zinc2jet,
                             SpTsub(genLeptons[0].v, genLeptons[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSpT_HighPt_Zinc2jet,
                             SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSPhi_HighPt_Zinc2jet,
                             SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                    }
                    if (nGoodGenJets == 2) {
                        //////Special Branch/////////////////////////
                        fill(genAbsZRapidity_Zexc2jet,
                             fabs(genEWKBoson.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genAbsSecondJetRapidity_Zexc2jet,
                             fabs(genJets[1].v.Rapidity()),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSumZSecondJetRapidity_Zexc2jet,
                             fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);
                        fill(genDifZSecondJetRapidity_Zexc2jet,
                             fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                             commonGenWeight,
                             EvtWeights);

                        if (genEWKBoson.Pt() > 100.) {
                            fill(genAbsZRapidity_ZPt100_Zexc2jet,
                                 fabs(genEWKBoson.Rapidity()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genAbsSecondJetRapidity_ZPt100_Zexc2jet,
                                 fabs(genJets[1].v.Rapidity()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSumZSecondJetRapidity_ZPt100_Zexc2jet,
                                 fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genDifZSecondJetRapidity_ZPt100_Zexc2jet,
                                 fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                                 commonGenWeight,
                                 EvtWeights);
                        }

                        if (genEWKBoson.Pt() > 150.) {
                            fill(genAbsZRapidity_ZPt150_Zexc2jet,
                                 fabs(genEWKBoson.Rapidity()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genAbsSecondJetRapidity_ZPt150_Zexc2jet,
                                 fabs(genJets[1].v.Rapidity()),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSumZSecondJetRapidity_ZPt150_Zexc2jet,
                                 fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genDifZSecondJetRapidity_ZPt150_Zexc2jet,
                                 fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                                 commonGenWeight,
                                 EvtWeights);
                        }

                        fill(genTwoJetsPtDiff_Zexc2jet,
                             genJet1Minus2.Pt(),
                             commonGenWeight,
                             EvtWeights);
                        fill(genJetsMass_Zexc2jet, genJet1Plus2.M(), commonGenWeight, EvtWeights);
                        fill(genSecondJetPt_Zexc2jet,
                             genJets[1].v.Pt(),
                             commonGenWeight,
                             EvtWeights);
                        fill(gendPhiJets_Zexc2jet,
                             deltaPhi(genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genPHI_Zexc2jet,
                             PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genPHI_T_Zexc2jet,
                             PHI_T(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(gendEtaJets_Zexc2jet,
                             genJets[0].v.Eta() - genJets[1].v.Eta(),
                             commonGenWeight,
                             EvtWeights);
                        fill(gendEtaFirstJetZ_Zexc2jet,
                             genJets[0].v.Eta() - genEWKBoson.Eta(),
                             commonGenWeight,
                             EvtWeights);
                        fill(gendEtaSecondJetZ_Zexc2jet,
                             genJets[1].v.Eta() - genEWKBoson.Eta(),
                             commonGenWeight,
                             EvtWeights);
                        fill(gendEtaJet1Plus2Z_Zexc2jet,
                             genJet1Plus2.Eta() - genEWKBoson.Eta(),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSpT_Zexc2jet,
                             SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSpTJets_Zexc2jet,
                             SpTsub(genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSpTLeptons_Zexc2jet,
                             SpTsub(genLeptons[0].v, genLeptons[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(genSPhi_Zexc2jet,
                             SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                             commonGenWeight,
                             EvtWeights);
                        fill(
                            genptBal_Zexc2jet, genJet1Plus2PlusZ.Pt(), commonGenWeight, EvtWeights);

                        if (genEWKBoson.Pt() < 25) {
                            fill(genptBal_LowPt_Zexc2jet,
                                 genJet1Plus2PlusZ.Pt(),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(gendPhiJets_LowPt_Zexc2jet,
                                 deltaPhi(genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(
                                genPHI_T_LowPt_Zexc2jet,
                                PHI_T(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                commonGenWeight,
                                EvtWeights);
                            fill(genPHI_LowPt_Zexc2jet,
                                 PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSpTJets_LowPt_Zexc2jet,
                                 SpTsub(genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSpTLeptons_LowPt_Zexc2jet,
                                 SpTsub(genLeptons[0].v, genLeptons[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSpT_LowPt_Zexc2jet,
                                 SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSPhi_LowPt_Zexc2jet,
                                 SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            if (SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v) <
                                0.5) {
                                fill(genPHI_LowSpT_LowPt_Zexc2jet,
                                     PHI(genLeptons[0].v,
                                         genLeptons[1].v,
                                         genJets[0].v,
                                         genJets[1].v),
                                     commonGenWeight,
                                     EvtWeights);
                                fill(genSPhi_LowSpT_LowPt_Zexc2jet,
                                     SPhi(genLeptons[0].v,
                                          genLeptons[1].v,
                                          genJets[0].v,
                                          genJets[1].v),
                                     commonGenWeight,
                                     EvtWeights);
                            } else {
                                fill(genPHI_HighSpT_LowPt_Zexc2jet,
                                     PHI(genLeptons[0].v,
                                         genLeptons[1].v,
                                         genJets[0].v,
                                         genJets[1].v),
                                     commonGenWeight,
                                     EvtWeights);
                                fill(genSPhi_HighSpT_LowPt_Zexc2jet,
                                     SPhi(genLeptons[0].v,
                                          genLeptons[1].v,
                                          genJets[0].v,
                                          genJets[1].v),
                                     commonGenWeight,
                                     EvtWeights);
                            }
                            if (SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v) <
                                0.5) {
                                fill(genSpT_LowSPhi_LowPt_Zexc2jet,
                                     SpT(genLeptons[0].v,
                                         genLeptons[1].v,
                                         genJets[0].v,
                                         genJets[1].v),
                                     commonGenWeight,
                                     EvtWeights);
                            } else {
                                fill(genSpT_HighSPhi_LowPt_Zexc2jet,
                                     SpT(genLeptons[0].v,
                                         genLeptons[1].v,
                                         genJets[0].v,
                                         genJets[1].v),
                                     commonGenWeight,
                                     EvtWeights);
                            }
                        } else {
                            fill(genptBal_HighPt_Zexc2jet,
                                 genJet1Plus2PlusZ.Pt(),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(gendPhiJets_HighPt_Zexc2jet,
                                 deltaPhi(genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genPHI_HighPt_Zexc2jet,
                                 PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(
                                genPHI_T_HighPt_Zexc2jet,
                                PHI_T(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                commonGenWeight,
                                EvtWeights);
                            fill(genSpTJets_HighPt_Zexc2jet,
                                 SpTsub(genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSpTLeptons_HighPt_Zexc2jet,
                                 SpTsub(genLeptons[0].v, genLeptons[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSpT_HighPt_Zexc2jet,
                                 SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                            fill(genSPhi_HighPt_Zexc2jet,
                                 SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                                 commonGenWeight,
                                 EvtWeights);
                        }
                    }
                }
                if (nGoodGenJets_20 >= 3) {

                    double RatioValue = 1.;
                    if (nGoodJets_20 >= 3 && UnfoldUnc) {
                        double binNumber =
                            ThirdJetPt_2_Zinc3jetHratio_fit->GetXaxis()->FindBin(jets_20[2].v.Pt());
                        RatioValue = ThirdJetPt_2_Zinc3jetHratio_fit->GetBinContent(binNumber);
                    }

                    fill(genThirdJetPt_Zinc3jet,
                         genJets_20[2].v.Pt(),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                }
                if (nGoodGenJets >= 3) {
                    nGenEventsVInc3Jets++;
                    nEffGenEventsVInc3Jets += genWeight;
                    fill(genZNGoodJets_Zinc, 3., commonGenWeight, EvtWeights);

                    double RatioValue = 1.;
                    double RatioValue1 = 1.;
                    double RatioValue2 = 1.;
                    if (nGoodJets >= 3 && UnfoldUnc) {
                        double binNumber =
                            ThirdJetAbsRapidity_2_Zinc3jetHratio_fit->GetXaxis()->FindBin(
                                fabs(jets[2].v.Eta()));
                        RatioValue =
                            ThirdJetAbsRapidity_2_Zinc3jetHratio_fit->GetBinContent(binNumber);
                        double binNumber1 =
                            JetsHT_2_Zinc3jetHratio_fit->GetXaxis()->FindBin(jetsHT);
                        // RatioValue1 = JetsHT_2_Zinc3jetHratio_fit->GetBinContent(binNumber1);
                        if (jetsHT >= 90. && jetsHT <= 1200.)
                            RatioValue1 = JetsHT_2_Zinc3jetHratio_fit->GetBinContent(binNumber1);
                        else
                            RatioValue1 = 1.;
                        if (RatioValue1 > 2. || RatioValue1 < 0.5) RatioValue1 = 1.;

                        double binNumber2 = VisPt_2_Zinc3jetQunHratio_fit->GetXaxis()->FindBin(
                            fabs((genJets[0].v + genJets[1].v + genJets[2].v + genEWKBoson).Pt()));
                        // RatioValue2 = VisPt_2_Zinc3jetQunHratio_fit->GetBinContent(binNumber2);
                        if (fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()) >= 0. &&
                            fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()) <= 200.)
                            RatioValue2 = VisPt_2_Zinc3jetQunHratio_fit->GetBinContent(binNumber2);
                        else
                            RatioValue2 = 1.;
                        if (RatioValue2 > 2. || RatioValue2 < 0.5) RatioValue2 = 1.;
                    }
                    // cout << RatioValue2 << "\n";

                    fill(genThirdJetEta_Zinc3jet,
                         fabs(genJets[2].v.Eta()),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                    fill(genThirdJetAbsRapidity_Zinc3jet,
                         fabs(genJets[2].v.Rapidity()),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                    fill(genThirdJetEtaHigh_Zinc3jet,
                         fabs(genJets[2].v.Eta()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genThirdJetRapidityHigh_Zinc3jet,
                         fabs(genJets[2].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genJetsHT_Zinc3jet, genJetsHT, commonGenWeight * RatioValue1, EvtWeights);
                    fill(genVisPt_Zinc3jetQun,
                         fabs((genJets[0].v + genJets[1].v + genJets[2].v + genEWKBoson).Pt()),
                         commonGenWeight * RatioValue2,
                         EvtWeights);

                    /////Azimuth cross check//////////////////////////////////
                    fill(genDPhiZFirstJet_Zinc3jet,
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         commonGenWeight,
                         EvtWeights);
                    fill(genDPhiZSecondJet_Zinc3jet,
                         fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                         commonGenWeight,
                         EvtWeights);
                    fill(genDPhiZThirdJet_Zinc3jet,
                         fabs(genEWKBoson.DeltaPhi(genJets[2].v)),
                         commonGenWeight,
                         EvtWeights);
                    fill(genDPhiFirstSecondJet_Zinc3jet,
                         fabs(genJets[0].v.DeltaPhi(genJets[1].v)),
                         commonGenWeight,
                         EvtWeights);
                    fill(genDPhiFirstThirdJet_Zinc3jet,
                         fabs(genJets[0].v.DeltaPhi(genJets[2].v)),
                         commonGenWeight,
                         EvtWeights);
                    fill(genDPhiSecondThirdJet_Zinc3jet,
                         fabs(genJets[1].v.DeltaPhi(genJets[2].v)),
                         commonGenWeight,
                         EvtWeights);

                    if (genEWKBoson.Pt() > 150.) {
                        fill(genDPhiZFirstJet_ZPt150_Zinc3jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiZSecondJet_ZPt150_Zinc3jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiZThirdJet_ZPt150_Zinc3jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[2].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiFirstSecondJet_ZPt150_Zinc3jet,
                             fabs(genJets[0].v.DeltaPhi(genJets[1].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiFirstThirdJet_ZPt150_Zinc3jet,
                             fabs(genJets[0].v.DeltaPhi(genJets[2].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiSecondThirdJet_ZPt150_Zinc3jet,
                             fabs(genJets[1].v.DeltaPhi(genJets[2].v)),
                             commonGenWeight,
                             EvtWeights);
                    }

                    if (genEWKBoson.Pt() > 300.) {
                        fill(genDPhiZFirstJet_ZPt300_Zinc3jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiZSecondJet_ZPt300_Zinc3jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiZThirdJet_ZPt300_Zinc3jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[2].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiFirstSecondJet_ZPt300_Zinc3jet,
                             fabs(genJets[0].v.DeltaPhi(genJets[1].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiFirstThirdJet_ZPt300_Zinc3jet,
                             fabs(genJets[0].v.DeltaPhi(genJets[2].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiSecondThirdJet_ZPt300_Zinc3jet,
                             fabs(genJets[1].v.DeltaPhi(genJets[2].v)),
                             commonGenWeight,
                             EvtWeights);
                    }

                    if (genEWKBoson.Pt() > 150. &&
                        (genJets[0].v.Pt() + genJets[1].v.Pt() + genJets[2].v.Pt() > 300.)) {
                        fill(genDPhiZFirstJet_ZPt150_HT300_Zinc3jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiZSecondJet_ZPt150_HT300_Zinc3jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                             commonGenWeight,
                             EvtWeights);
                        fill(genDPhiZThirdJet_ZPt150_HT300_Zinc3jet,
                             fabs(genEWKBoson.DeltaPhi(genJets[2].v)),
                             commonGenWeight,
                             EvtWeights);
                    }
                }
                if (nGoodGenJets_20 >= 4)
                    fill(
                        genFourthJetPt_Zinc4jet, genJets_20[3].v.Pt(), commonGenWeight, EvtWeights);
                if (nGoodGenJets >= 4) {
                    fill(genZNGoodJets_Zinc, 4., commonGenWeight, EvtWeights);
                    fill(genFourthJetEta_Zinc4jet,
                         fabs(genJets[3].v.Eta()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genFourthJetAbsRapidity_Zinc4jet,
                         fabs(genJets[3].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genFourthJetEtaHigh_Zinc4jet,
                         fabs(genJets[3].v.Eta()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genFourthJetRapidityHigh_Zinc4jet,
                         fabs(genJets[3].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genJetsHT_Zinc4jet, genJetsHT, commonGenWeight, EvtWeights);
                }
                if (nGoodGenJets_20 >= 5)
                    fill(genFifthJetPt_Zinc5jet, genJets_20[4].v.Pt(), commonGenWeight, EvtWeights);
                if (nGoodGenJets >= 5) {
                    fill(genZNGoodJets_Zinc, 5., commonGenWeight, EvtWeights);
                    fill(genFifthJetEta_Zinc5jet,
                         fabs(genJets[4].v.Eta()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genFifthJetAbsRapidity_Zinc5jet,
                         fabs(genJets[4].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genFifthJetEtaHigh_Zinc5jet,
                         fabs(genJets[4].v.Eta()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genFifthJetRapidityHigh_Zinc5jet,
                         fabs(genJets[4].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genJetsHT_Zinc5jet, genJetsHT, commonGenWeight, EvtWeights);
                }
                if (nGoodGenJets_20 >= 6)
                    fill(genSixthJetPt_Zinc6jet, genJets_20[5].v.Pt(), commonGenWeight, EvtWeights);
                if (nGoodGenJets >= 6) {
                    fill(genZNGoodJets_Zinc, 6., commonGenWeight, EvtWeights);
                    fill(genSixthJetEta_Zinc6jet,
                         fabs(genJets[5].v.Eta()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSixthJetEtaHigh_Zinc6jet,
                         fabs(genJets[5].v.Eta()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genSixthJetRapidityHigh_Zinc6jet,
                         fabs(genJets[5].v.Rapidity()),
                         commonGenWeight,
                         EvtWeights);
                    fill(genJetsHT_Zinc6jet, genJetsHT, commonGenWeight, EvtWeights);
                }
                if (nGoodGenJets >= 7) {
                    fill(genZNGoodJets_Zinc, 7., commonGenWeight, EvtWeights);
                }
                if (nGoodGenJets >= 8) {
                    fill(genZNGoodJets_Zinc, 8., commonGenWeight, EvtWeights);
                }
                if (nGoodGenJets >= 1) {

                    double RatioValue = 1.;
                    double RatioValue1 = 1.;
                    double RatioValue2 = 1.;
                    if (UnfoldUnc) {
                        double binNumber =
                            JZB_2Hratio_fit->GetXaxis()->FindBin(hadronicR.Pt() - EWKBoson.Pt());
                        RatioValue = JZB_2Hratio_fit->GetBinContent(binNumber);

                        // if(EWKBoson.Pt()<= 50){   //check if gen or reco!!!!
                        double binNumber1 = JZB_ptLow_2Hratio_fit->GetXaxis()->FindBin(
                            hadronicR.Pt() - EWKBoson.Pt());
                        if ((hadronicR.Pt() - EWKBoson.Pt()) > -50 &&
                            (hadronicR.Pt() - EWKBoson.Pt()) < 200.)
                            RatioValue1 = JZB_ptLow_2Hratio_fit->GetBinContent(binNumber1);
                        else
                            RatioValue1 = 1.;
                        if (RatioValue1 > 2. || RatioValue1 < 0.5) RatioValue1 = 1.;
                        // cout << RatioValue1 << "\n";

                        // }
                        // if(EWKBoson.Pt()> 50){
                        double binNumber2 = JZB_ptHigh_2Hratio_fit->GetXaxis()->FindBin(
                            hadronicR.Pt() - EWKBoson.Pt());
                        RatioValue2 = JZB_ptHigh_2Hratio_fit->GetBinContent(binNumber2);
                        //}
                    }
                    // cout << hadronicR.Pt()-EWKBoson.Pt() << " , " << RatioValue << " , " <<
                    // RatioValue1 << " , " << RatioValue2 << "\n";

                    fill(genHadRecoil, genHadronicR.Pt(), commonGenWeight, EvtWeights);
                    fill(genJZB,
                         genHadronicR.Pt() - genEWKBoson.Pt(),
                         commonGenWeight * RatioValue,
                         EvtWeights);
                    if (genEWKBoson.Pt() <= 50)
                        fill(genJZB_ptLow,
                             genHadronicR.Pt() - genEWKBoson.Pt(),
                             commonGenWeight * RatioValue1,
                             EvtWeights);
                    else {
                        fill(genJZB_ptHigh,
                             genHadronicR.Pt() - genEWKBoson.Pt(),
                             commonGenWeight * RatioValue2,
                             EvtWeights);
                    }
                }
            }
        }
        //=======================================================================================================//
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

        //=======================================================================================================//
        //      Selection for Reco Histos      //
        //====================================//
        double tau_sum(0), tau_max(0);
        double tau_c_sum(0), tau_c_max(0);
        double tau_cm_sum(0), tau_cm_max(0);
        double tau_c_cm_sum(0), tau_c_cm_max(0);

        if (hasRecoInfo && passesLeptonChargeCut && passesTauCut) {
            fill(ZMassFrom60_Zinc0jet, (leptons[0].v + leptons[1].v).M(), weight);
        }
        //  cout << passesLeptonCut << " , " << (!bTagJetFound || !rejectBTagEvents)  << "\n";
        if (hasRecoInfo && passesLeptonCutNoMass &&
            (!bTagJetFound || !rejectBTagEvents)) { // no Mass cut, fill the whole range
            fill(Mass_Zinc0jet, EWKBoson.M(), weight);
            if (EWKBoson.M() > 115. && EWKBoson.M() < 135.)
                fill(ZPt_Zinc0jetM115_135, EWKBoson.Pt(), weight);
            if (EWKBoson.M() >= 111. && EWKBoson.M() < 130.)
                fill(ZPt_Zinc0jetM111_130, EWKBoson.Pt(), weight);
            if (EWKBoson.M() >= 130. && EWKBoson.M() <= 3000.)
                fill(ZPt_Zinc0jetM130_3, EWKBoson.Pt(), weight);
        }

        if (hasRecoInfo && passesLeptonCut && (!bTagJetFound || !rejectBTagEvents)) {
            //=======================================================================================================//

            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            //=======================================================================================================//
            //      Start filling histograms      //
            //====================================//

            // cout << "Selected at reco level" << endl;
            fill(NVtx, EvtVtxCnt, weight);

            double weightNoPUweight(weight);
            if (hasRecoInfo && !EvtIsRealData)
                weightNoPUweight /= puWeight.weight(int(EvtPuCntTruth));
            fill(NVtx_NoPUweight, EvtVtxCnt, weightNoPUweight);

            nEventsVInc0Jets++;
            nEffEventsVInc0Jets += weight;
            fill(ZNGoodJetsNVtx_Zexc, nGoodJets, EvtVtxCnt, weight);
            fill(ZNGoodJets_Zinc, 0., weight);

            double RatioValue = 1;
            if (UnfoldUnc) {
                double binNumber = ZNGoodJets_ZexcHratio_fit->GetXaxis()->FindBin(nGoodJets);
                RatioValue = ZNGoodJets_ZexcHratio_fit->GetBinContent(binNumber);
            }

            fill(ZNGoodJets_Zexc, nGoodJets, weight * RatioValue);
            fill(ZNGoodJets_Zinc_NoWeight, 0.);
            fill(ZMass_Zinc0jet, EWKBoson.M(), weight);
            fill(ZPt_Zinc0jet, EWKBoson.Pt(), weight);
            fill(ZRapidity_Zinc0jet, EWKBoson.Rapidity(), weight);
            fill(ZEta_Zinc0jet, EWKBoson.Eta(), weight);
            fill(ZEtaUpTo5_Zinc0jet, EWKBoson.Eta(), weight);
            fill(lepPt_Zinc0jet, leptons[0].v.Pt(), weight);
            fill(lepEta_Zinc0jet, leptons[0].v.Eta(), weight);
            fill(lepPhi_Zinc0jet, leptons[0].v.Phi(), weight);
            fill(lepPt_Zinc0jet, leptons[1].v.Pt(), weight);
            fill(lepEta_Zinc0jet, leptons[1].v.Eta(), weight);
            fill(lepPhi_Zinc0jet, leptons[1].v.Phi(), weight);
            fill(dPhiLeptons_Zinc0jet, deltaPhi(leptons[0].v, leptons[1].v), weight);
            fill(dEtaLeptons_Zinc0jet, leptons[0].v.Eta() - leptons[1].v.Eta(), weight);
            fill(dRLeptons_Zinc0jet, deltaR(leptons[0].v, leptons[1].v), weight);
            fill(SpTLeptons_Zinc0jet, SpTsub(leptons[0].v, leptons[1].v), weight);
            fill(VisPt_Zinc0jetQun, EWKBoson.Pt(), weight);
            fill(VisPt_2_Zinc0jetQun, EWKBoson.Pt(), weight);
            if (nEvents % 2)
                fill(VisPt_Zinc0jetQun_Odd, EWKBoson.Pt(), weight);
            else
                fill(VisPt_Zinc0jetQun_Even, EWKBoson.Pt(), weight);
            // 2016
            fill(Phistar_Zinc0jet, phistar, weight);

            if (nGoodJets == 0) {
                // fill(TruePU_0, EvtPuCntTruth, weight);
                // fill(PU_0, EvtPuCnt, weight);
                fill(PU_0, EvtVtxCnt, weight);
                fill(ZNGoodJets_Zexc_NoWeight, 0.);
                fill(ZPt_Zexc0jet, EWKBoson.Pt(), weight);
                fill(ZRapidity_Zexc0jet, EWKBoson.Rapidity(), weight);
                fill(ZEta_Zexc0jet, EWKBoson.Eta(), weight);
                fill(lepPt_Zexc0jet, leptons[0].v.Pt(), weight);
                fill(lepEta_Zexc0jet, leptons[0].v.Eta(), weight);
                fill(lepPt_Zexc0jet, leptons[1].v.Pt(), weight);
                fill(lepEta_Zexc0jet, leptons[1].v.Eta(), weight);
                fill(dPhiLeptons_Zexc0jet, deltaPhi(leptons[0].v, leptons[1].v), weight);
                fill(dEtaLeptons_Zexc0jet, leptons[0].v.Eta() - leptons[1].v.Eta(), weight);
                fill(SpTLeptons_Zexc0jet, SpTsub(leptons[0].v, leptons[0].v), weight);
            }

            if (nGoodJets_20 >= 1) {

                double RatioValue = 1;
                if (UnfoldUnc) {
                    double binNumber =
                        FirstJetPt_2_Zinc1jetHratio_fit->GetXaxis()->FindBin(jets_20[0].v.Pt());
                    RatioValue = FirstJetPt_2_Zinc1jetHratio_fit->GetBinContent(binNumber);
                }

                fill(FirstJetPt_Zinc1jet, jets_20[0].v.Pt(), weight * RatioValue);
                if (nEvents % 2)
                    fill(FirstJetPt_Zinc1jet_Odd, jets_20[0].v.Pt(), weight * RatioValue);
                else
                    fill(FirstJetPt_Zinc1jet_Even, jets_20[0].v.Pt(), weight * RatioValue);
                fill(FirstJetPt_2_Zinc1jet, jets_20[0].v.Pt(), weight);
                fill(FirstJetPt_Zinc1jet_NVtx, jets_20[0].v.Pt(), EvtVtxCnt, weight);
                fill(FirstJetPtEta_Zinc1jet, jets_20[0].v.Pt(), fabs(jets[0].v.Eta()), weight);
            }

            if (nGoodJets >= 1) {
                nEventsVInc1Jets++;
                nEffEventsVInc1Jets += weight;
                fill(AbsZRapidity_Zinc1jet, fabs(EWKBoson.Rapidity()), weight);
                fill(AbsFirstJetRapidity_Zinc1jet, fabs(jets[0].v.Rapidity()), weight);
                fill(SumZFirstJetRapidity_Zinc1jet,
                     fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                     weight);
                fill(DifZFirstJetRapidity_Zinc1jet,
                     fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                     weight);

                // cross check//////
                fill(
                    SumZFirstJetEta_Zinc1jet, fabs(EWKBoson.Eta() + jets[0].v.Eta()) / 2.0, weight);
                fill(
                    DifZFirstJetEta_Zinc1jet, fabs(EWKBoson.Eta() - jets[0].v.Eta()) / 2.0, weight);

                /// Azimuth cross check/////////////////////////////
                fill(DPhiZFirstJet_Zinc1jet, fabs(EWKBoson.DeltaPhi(jets[0].v)), weight);

                if (EWKBoson.Pt() > 100.) {
                    fill(AbsZRapidity_ZPt100_Zinc1jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsFirstJetRapidity_ZPt100_Zinc1jet, fabs(jets[0].v.Rapidity()), weight);
                    fill(SumZFirstJetRapidity_ZPt100_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZFirstJetRapidity_ZPt100_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                if (EWKBoson.Pt() > 150.) {
                    fill(AbsZRapidity_ZPt150_Zinc1jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsFirstJetRapidity_ZPt150_Zinc1jet, fabs(jets[0].v.Rapidity()), weight);
                    fill(SumZFirstJetRapidity_ZPt150_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZFirstJetRapidity_ZPt150_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         weight);

                    fill(DPhiZFirstJet_ZPt150_Zinc1jet, fabs(EWKBoson.DeltaPhi(jets[0].v)), weight);
                }

                if (EWKBoson.Pt() > 300.) {
                    fill(AbsZRapidity_ZPt300_Zinc1jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsFirstJetRapidity_ZPt300_Zinc1jet, fabs(jets[0].v.Rapidity()), weight);
                    fill(SumZFirstJetRapidity_ZPt300_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZFirstJetRapidity_ZPt300_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         weight);

                    fill(DPhiZFirstJet_ZPt300_Zinc1jet, fabs(EWKBoson.DeltaPhi(jets[0].v)), weight);
                }

                /// different JetPt Cuts//////
                if (jets[0].v.Pt() > 50.) {
                    fill(AbsZRapidity_FirstJetPt50_Zinc1jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsFirstJetRapidity_FirstJetPt50_Zinc1jet,
                         fabs(jets[0].v.Rapidity()),
                         weight);
                    fill(SumZFirstJetRapidity_FirstJetPt50_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZFirstJetRapidity_FirstJetPt50_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                if (jets[0].v.Pt() > 80.) {
                    fill(AbsZRapidity_FirstJetPt80_Zinc1jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsFirstJetRapidity_FirstJetPt80_Zinc1jet,
                         fabs(jets[0].v.Rapidity()),
                         weight);
                    fill(SumZFirstJetRapidity_FirstJetPt80_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZFirstJetRapidity_FirstJetPt80_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                fill(ZNGoodJets_Zinc, 1., weight);
                fill(ZNGoodJets_Zinc_NoWeight, 1.);
                fill(ZPt_Zinc1jet, EWKBoson.Pt(), weight);
                fill(lepPt_Zinc1jet, leptons[0].v.Pt(), weight);
                fill(lepEta_Zinc1jet, leptons[0].v.Eta(), weight);
                fill(lepPt_Zinc1jet, leptons[1].v.Pt(), weight);
                fill(lepEta_Zinc1jet, leptons[1].v.Eta(), weight);
                fill(dPhiLeptons_Zinc1jet, deltaPhi(leptons[0].v, leptons[1].v), weight);
                fill(dRLeptons_Zinc1jet, deltaPhi(leptons[0].v, leptons[1].v), weight);
                fill(dEtaLeptons_Zinc1jet, leptons[0].v.Eta() - leptons[1].v.Eta(), weight);
                fill(ZMass_Zinc1jet, EWKBoson.M(), weight);
                fill(ZRapidity_Zinc1jet, EWKBoson.Rapidity(), weight);
                fill(ZAbsRapidity_Zinc1jet, fabs(EWKBoson.Rapidity()), weight);
                fill(ZEta_Zinc1jet, EWKBoson.Eta(), weight);
                fill(ZEtaUpTo5_Zinc1jet, EWKBoson.Eta(), weight);
                fill(SpTLeptons_Zinc1jet, SpTsub(leptons[0].v, leptons[1].v), weight);
                fill(Phistar_Zinc1jet, phistar, weight);

                double RatioValue = 1.;
                double RatioValue1 = 1.;
                double RatioValue2 = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        FirstJetAbsRapidity_2_Zinc1jetHratio_fit->GetXaxis()->FindBin(
                            fabs(jets[0].v.Eta()));
                    RatioValue = FirstJetAbsRapidity_2_Zinc1jetHratio_fit->GetBinContent(binNumber);
                    double binNumber1 = JetsHT_2_Zinc1jetHratio_fit->GetXaxis()->FindBin(jetsHT);
                    RatioValue1 = JetsHT_2_Zinc1jetHratio_fit->GetBinContent(binNumber1);
                    double binNumber2 = VisPt_2_Zinc1jetQunHratio_fit->GetXaxis()->FindBin(
                        fabs((jets[0].v + EWKBoson).Pt()));
                    RatioValue2 = VisPt_2_Zinc1jetQunHratio_fit->GetBinContent(binNumber2);
                }

                fill(FirstJetEta_Zinc1jet, fabs(jets[0].v.Eta()), weight * RatioValue);
                fill(FirstJetEta_2_Zinc1jet, fabs(jets[0].v.Eta()), weight);
                fill(FirstJetAbsRapidity_Zinc1jet, fabs(jets[0].v.Rapidity()), weight * RatioValue);
                fill(FirstJetAbsRapidity_2_Zinc1jet,
                     fabs(jets[0].v.Rapidity()),
                     weight * RatioValue);
                if (nEvents % 2)
                    fill(FirstJetAbsRapidity_Zinc1jet_Odd,
                         fabs(jets[0].v.Rapidity()),
                         weight * RatioValue);
                else
                    fill(FirstJetAbsRapidity_Zinc1jet_Even,
                         fabs(jets[0].v.Rapidity()),
                         weight * RatioValue);
                fill(FirstJetEtaHigh_Zinc1jet, fabs(jets[0].v.Eta()), weight);
                fill(FirstJetRapidityHigh_Zinc1jet, fabs(jets[0].v.Rapidity()), weight);
                fill(FirstJetEtaFull_Zinc1jet, jets[0].v.Eta(), weight);
                fill(FirstJetPhi_Zinc1jet, jets[0].v.Phi(), weight);
                fill(JetsHT_Zinc1jet, jetsHT, weight * RatioValue1);
                fill(JetsHT_2_Zinc1jet, jetsHT, weight);
                fill(dEtaBosonJet_Zinc1jet, fabs(jets[0].v.Eta() - EWKBoson.Eta()), weight);
                fill(SumZJetRapidity_Zinc1jet,
                     0.5 * fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()),
                     weight);
                fill(DifZJetRapidity_Zinc1jet,
                     0.5 * fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()),
                     weight);
                fill(VisPt_Zinc1jetQun, fabs((jets[0].v + EWKBoson).Pt()), weight * RatioValue2);
                fill(VisPt_2_Zinc1jetQun, fabs((jets[0].v + EWKBoson).Pt()), weight * RatioValue2);
                if (nEvents % 2)
                    fill(VisPt_Zinc1jetQun_Odd,
                         fabs((jets[0].v + EWKBoson).Pt()),
                         weight * RatioValue2);
                else
                    fill(VisPt_Zinc1jetQun_Even,
                         fabs((jets[0].v + EWKBoson).Pt()),
                         weight * RatioValue2);

                for (unsigned short i(0); i < nGoodJets; i++) {
                    double trans_mass = jets[i].v.Mt();

                    tau_sum += trans_mass * exp(-fabs(jets[i].v.Rapidity() - EWKBoson.Rapidity()));
                    tau_max =
                        max(tau_max,
                            trans_mass * exp(-fabs(jets[i].v.Rapidity() - EWKBoson.Rapidity())));

                    tau_c_sum +=
                        trans_mass / (2 * cosh(jets[i].v.Rapidity() - EWKBoson.Rapidity()));
                    tau_c_max =
                        max(tau_c_max,
                            trans_mass / (2 * cosh(jets[i].v.Rapidity() - EWKBoson.Rapidity())));

                    tau_cm_sum += trans_mass * exp(-fabs(jets[i].v.Rapidity()));
                    tau_cm_max = max(tau_cm_max, trans_mass * exp(-fabs(jets[i].v.Rapidity())));

                    tau_c_cm_sum += trans_mass / (2 * cosh(jets[i].v.Rapidity()));
                    tau_c_cm_max = max(tau_c_cm_max, trans_mass / (2 * cosh(jets[i].v.Rapidity())));
                }

                for (unsigned short i(0); i < 5; i++) {
                    if (EWKBoson.Pt() > ZptRange[i] && EWKBoson.Pt() <= ZptRange[i + 1]) {
                        fill(tau_sum_Zinc1jet[i], tau_sum, weight);
                        fill(tau_max_Zinc1jet[i], tau_max, weight);
                        fill(tau_c_sum_Zinc1jet[i], tau_c_sum, weight);
                        fill(tau_c_max_Zinc1jet[i], tau_c_max, weight);
                        fill(tau_cm_sum_Zinc1jet[i], tau_cm_sum, weight);
                        fill(tau_cm_max_Zinc1jet[i], tau_cm_max, weight);
                        fill(tau_c_cm_sum_Zinc1jet[i], tau_c_cm_sum, weight);
                        fill(tau_c_cm_max_Zinc1jet[i], tau_c_cm_max, weight);
                    }
                }

                if (nGoodJets == 1) {
                    // fill(TruePU_1, EvtPuCntTruth, weight);
                    // fill(PU_1, EvtPuCnt, weight);
                    fill(PU_1, EvtVtxCnt, weight);
                    fill(ZNGoodJets_Zexc_NoWeight, 1.);
                    fill(ZPt_Zexc1jet, EWKBoson.Pt(), weight);
                    fill(ZRapidity_Zexc1jet, EWKBoson.Rapidity(), weight);
                    fill(ZEta_Zexc1jet, EWKBoson.Eta(), weight);
                    fill(SpTLeptons_Zexc1jet, SpTsub(leptons[0].v, leptons[1].v), weight);
                    fill(FirstJetPt_Zexc1jet, jets[0].v.Pt(), weight);
                    fill(FirstJetEta_Zexc1jet, jets[0].v.Eta(), weight);
                    fill(FirstJetPhi_Zexc1jet, jets[0].v.Phi(), weight);
                    fill(dEtaBosonJet_Zexc1jet, fabs(jets[0].v.Eta() - EWKBoson.Eta()), weight);
                    // Additional Branch
                    fill(AbsZRapidity_Zexc1jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsJetRapidity_Zexc1jet, fabs(jets[0].v.Rapidity()), weight);
                    fill(SumZJetRapidity_Zexc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZJetRapidity_Zexc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         weight);

                    if (EWKBoson.Pt() > 100.) {
                        fill(AbsZRapidity_ZPt100_Zexc1jet, fabs(EWKBoson.Rapidity()), weight);
                        fill(AbsJetRapidity_ZPt100_Zexc1jet, fabs(jets[0].v.Rapidity()), weight);
                        fill(SumZJetRapidity_ZPt100_Zexc1jet,
                             fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                             weight);
                        fill(DifZJetRapidity_ZPt100_Zexc1jet,
                             fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                             weight);
                    }

                    if (EWKBoson.Pt() > 150.) {
                        fill(AbsZRapidity_ZPt150_Zexc1jet, fabs(EWKBoson.Rapidity()), weight);
                        fill(AbsJetRapidity_ZPt150_Zexc1jet, fabs(jets[0].v.Rapidity()), weight);
                        fill(SumZJetRapidity_ZPt150_Zexc1jet,
                             fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                             weight);
                        fill(DifZJetRapidity_ZPt150_Zexc1jet,
                             fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                             weight);
                    }
                }
            }
            if (nGoodJets_20 >= 2) {
                double RatioValue = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        SecondJetPt_2_Zinc2jetHratio_fit->GetXaxis()->FindBin(jets_20[1].v.Pt());
                    RatioValue = SecondJetPt_2_Zinc2jetHratio_fit->GetBinContent(binNumber);
                }

                fill(SecondJetPt_Zinc2jet, jets_20[1].v.Pt(), weight * RatioValue);
                if (nEvents % 2)
                    fill(SecondJetPt_Zinc2jet_Odd, jets_20[1].v.Pt(), weight * RatioValue);
                else
                    fill(SecondJetPt_Zinc2jet_Even, jets_20[1].v.Pt(), weight * RatioValue);
                fill(SecondJetPt_2_Zinc2jet, jets_20[1].v.Pt(), weight * RatioValue);
            }
            if (nGoodJets >= 2) {
                nEventsVInc2Jets++;
                nEffEventsVInc2Jets += weight;
                //////////////////Special Branch////////////////////////
                fill(AbsFirstJetRapidity_Zinc2jet, fabs(jets[0].v.Rapidity()), weight);
                fill(SumZFirstJetRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                     weight);
                fill(DifZFirstJetRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                     weight);

                fill(AbsZRapidity_Zinc2jet, fabs(EWKBoson.Rapidity()), weight);
                fill(AbsSecondJetRapidity_Zinc2jet, fabs(jets[1].v.Rapidity()), weight);
                fill(SumZSecondJetRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                     weight);
                fill(DifZSecondJetRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                     weight);

                fill(SumFirstSecondJetRapidity_Zinc2jet,
                     fabs(jets[0].v.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                     weight);
                fill(DifFirstSecondJetRapidity_Zinc2jet,
                     fabs(jets[0].v.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                     weight);

                TLorentzVector DiJets = jets[0].v + jets[1].v;
                fill(SumZTwoJetsRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() + DiJets.Rapidity()) / 2.0,
                     weight);
                fill(DifZTwoJetsRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() - DiJets.Rapidity()) / 2.0,
                     weight);

                /// Azimuth cross check
                fill(DPhiZFirstJet_Zinc2jet, fabs(EWKBoson.DeltaPhi(jets[0].v)), weight);
                fill(DPhiZSecondJet_Zinc2jet, fabs(EWKBoson.DeltaPhi(jets[1].v)), weight);
                fill(DPhiFirstSecondJet_Zinc2jet, fabs(jets[0].v.DeltaPhi(jets[1].v)), weight);

                if (EWKBoson.Pt() > 100.) {
                    fill(AbsZRapidity_ZPt100_Zinc2jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsSecondJetRapidity_ZPt100_Zinc2jet, fabs(jets[1].v.Rapidity()), weight);
                    fill(SumZSecondJetRapidity_ZPt100_Zinc2jet,
                         fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZSecondJetRapidity_ZPt100_Zinc2jet,
                         fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                         weight);
                }

                if (EWKBoson.Pt() > 150.) {
                    fill(AbsZRapidity_ZPt150_Zinc2jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsSecondJetRapidity_ZPt150_Zinc2jet, fabs(jets[1].v.Rapidity()), weight);
                    fill(SumZSecondJetRapidity_ZPt150_Zinc2jet,
                         fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZSecondJetRapidity_ZPt150_Zinc2jet,
                         fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                         weight);

                    fill(DPhiZFirstJet_ZPt150_Zinc2jet, fabs(EWKBoson.DeltaPhi(jets[0].v)), weight);
                }

                if (EWKBoson.Pt() > 300.) {
                    fill(DPhiZFirstJet_ZPt300_Zinc2jet, fabs(EWKBoson.DeltaPhi(jets[0].v)), weight);
                }

                // set jet rapidity discriminator/////

                if (fabs(jets[0].v.Rapidity() - jets[1].v.Rapidity()) > 2) {
                    fill(AbsZRapidity_DifJetRapidityl2_Zinc2jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet,
                         fabs(jets[0].v.Rapidity()),
                         weight);
                    fill(SumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                if (fabs(jets[0].v.Rapidity() - jets[1].v.Rapidity()) < 2) {
                    fill(AbsZRapidity_DifJetRapiditys2_Zinc2jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet,
                         fabs(jets[0].v.Rapidity()),
                         weight);
                    fill(SumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                TLorentzVector jet1Plus2PlusZ = jet1Plus2 + EWKBoson;
                fill(ZNGoodJets_Zinc, 2., weight);
                fill(ZNGoodJets_Zinc_NoWeight, 2.);
                fill(TwoJetsPtDiff_Zinc2jet, jet1Minus2.Pt(), weight);
                fill(BestTwoJetsPtDiff_Zinc2jet, bestJet1Minus2.Pt(), weight);
                fill(JetsMass_Zinc2jet, jet1Plus2.M(), weight);
                fill(llJetsMass_Zinc2jet, jet1Plus2PlusZ.M(), genWeight);

                if (jet1Plus2PlusZ.M() > 450 && jet1Plus2PlusZ.M() < 600) {
                    if (fabs(jets[0].v.Eta()) < fabs(jets[1].v.Eta())) {
                        fill(CentralJetEta_Zinc2jet, fabs(jets[0].v.Eta()), weight);
                        fill(ForwardJetEta_Zinc2jet, fabs(jets[1].v.Eta()), weight);
                        fill(CentralJetPt_Zinc2jet, jets[0].v.Pt(), weight);
                        fill(ForwardJetPt_Zinc2jet, jets[1].v.Pt(), weight);
                    } else {
                        fill(CentralJetEta_Zinc2jet, fabs(jets[1].v.Eta()), weight);
                        fill(ForwardJetEta_Zinc2jet, fabs(jets[0].v.Eta()), weight);
                        fill(CentralJetPt_Zinc2jet, jets[1].v.Pt(), weight);
                        fill(ForwardJetPt_Zinc2jet, jets[0].v.Pt(), weight);
                    }
                }

                double RatioValue = 1.;
                double RatioValue1 = 1.;
                double RatioValue2 = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        SecondJetAbsRapidity_2_Zinc2jetHratio_fit->GetXaxis()->FindBin(
                            fabs(jets[1].v.Eta()));
                    RatioValue =
                        SecondJetAbsRapidity_2_Zinc2jetHratio_fit->GetBinContent(binNumber);
                    double binNumber1 = JetsHT_2_Zinc2jetHratio_fit->GetXaxis()->FindBin(jetsHT);
                    RatioValue1 = JetsHT_2_Zinc2jetHratio_fit->GetBinContent(binNumber1);
                    double binNumber2 = VisPt_2_Zinc2jetQunHratio_fit->GetXaxis()->FindBin(
                        fabs((jets[0].v + jets[1].v + EWKBoson).Pt()));
                    RatioValue2 = VisPt_2_Zinc2jetQunHratio_fit->GetBinContent(binNumber2);
                }

                if (EvtVtxCnt < 14)
                    fill(JetsMassLowPU_Zinc2jet, jet1Plus2.M(), weight);
                else if (EvtVtxCnt < 18)
                    fill(JetsMassMidPU_Zinc2jet, jet1Plus2.M(), weight);
                else
                    fill(JetsMassHigPU_Zinc2jet, jet1Plus2.M(), weight);
                fill(ZPt_Zinc2jet, EWKBoson.Pt(), weight);
                fill(VisPt_Zinc2jetQun,
                     fabs((jets[0].v + jets[1].v + EWKBoson).Pt()),
                     weight * RatioValue2);
                fill(VisPt_2_Zinc2jetQun,
                     fabs((jets[0].v + jets[1].v + EWKBoson).Pt()),
                     weight * RatioValue2);
                if (nEvents % 2)
                    fill(VisPt_Zinc2jetQun_Odd,
                         fabs((jets[0].v + jets[1].v + EWKBoson).Pt()),
                         weight * RatioValue2);
                else
                    fill(VisPt_Zinc2jetQun_Even,
                         fabs((jets[0].v + jets[1].v + EWKBoson).Pt()),
                         weight * RatioValue2);
                fill(ZRapidity_Zinc2jet, EWKBoson.Rapidity(), weight);
                fill(ZEta_Zinc2jet, EWKBoson.Eta(), weight);
                fill(SpTLeptons_Zinc2jet, SpTsub(leptons[0].v, leptons[1].v), weight);

                fill(SecondJetEta_Zinc2jet, fabs(jets[1].v.Eta()), weight * RatioValue);
                fill(SecondJetEta_2_Zinc2jet, fabs(jets[1].v.Eta()), weight);
                fill(
                    SecondJetAbsRapidity_Zinc2jet, fabs(jets[1].v.Rapidity()), weight * RatioValue);
                fill(SecondJetAbsRapidity_2_Zinc2jet,
                     fabs(jets[1].v.Rapidity()),
                     weight * RatioValue);
                if (nEvents % 2)
                    fill(SecondJetAbsRapidity_Zinc2jet_Odd,
                         fabs(jets[1].v.Rapidity()),
                         weight * RatioValue);
                else
                    fill(SecondJetAbsRapidity_Zinc2jet_Even,
                         fabs(jets[1].v.Rapidity()),
                         weight * RatioValue);
                fill(SecondJetEtaHigh_Zinc2jet, fabs(jets[1].v.Eta()), weight);
                fill(SecondJetRapidityHigh_Zinc2jet, fabs(jets[1].v.Rapidity()), weight);
                fill(SecondJetEtaFull_Zinc2jet, jets[1].v.Eta(), weight);
                fill(SecondJetPhi_Zinc2jet, jets[1].v.Phi(), weight);
                fill(JetsHT_Zinc2jet, jetsHT, weight * RatioValue1);
                if (nEvents % 2)
                    fill(JetsHT_Zinc2jet_Odd, jetsHT, weight * RatioValue1);
                else
                    fill(JetsHT_Zinc2jet_Even, jetsHT, weight * RatioValue1);
                fill(JetsHT_2_Zinc2jet, jetsHT, weight);
                fill(ptBal_Zinc2jet, jet1Plus2PlusZ.Pt(), weight);
                fill(dPhiJets_Zinc2jet, deltaPhi(jets[0].v, jets[1].v), weight);
                fill(
                    BestdPhiJets_Zinc2jet, deltaPhi(bestTwoJets.first, bestTwoJets.second), weight);
                fill(dEtaJets_Zinc2jet, jets[0].v.Eta() - jets[1].v.Eta(), weight);
                fill(dEtaFirstJetZ_Zinc2jet, jets[0].v.Eta() - EWKBoson.Eta(), weight);
                fill(dEtaSecondJetZ_Zinc2jet, jets[1].v.Eta() - EWKBoson.Eta(), weight);
                fill(dEtaJet1Plus2Z_Zinc2jet, jet1Plus2.Eta() - EWKBoson.Eta(), weight);
                fill(PHI_Zinc2jet, PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v), weight);
                fill(BestPHI_Zinc2jet,
                     PHI(leptons[0].v, leptons[1].v, bestTwoJets.first, bestTwoJets.second),
                     weight);
                fill(PHI_T_Zinc2jet,
                     PHI_T(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                     weight);
                fill(BestPHI_T_Zinc2jet,
                     PHI_T(leptons[0].v, leptons[1].v, bestTwoJets.first, bestTwoJets.second),
                     weight);
                fill(SpT_Zinc2jet, SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v), weight);
                fill(BestSpT_Zinc2jet,
                     SpT(leptons[0].v, leptons[1].v, bestTwoJets.first, bestTwoJets.second),
                     weight);
                fill(SpTJets_Zinc2jet, SpTsub(jets[0].v, jets[1].v), weight);
                fill(BestSpTJets_Zinc2jet, SpTsub(bestTwoJets.first, bestTwoJets.second), weight);
                fill(SPhi_Zinc2jet, SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v), weight);
                fill(BestSPhi_Zinc2jet,
                     SPhi(leptons[0].v, leptons[1].v, bestTwoJets.first, bestTwoJets.second),
                     weight);

                if (EWKBoson.Pt() < 25) {
                    fill(ptBal_LowPt_Zinc2jet, jet1Plus2PlusZ.Pt(), weight);
                    fill(dPhiJets_LowPt_Zinc2jet, deltaPhi(jets[0].v, jets[1].v), weight);
                    fill(BestdPhiJets_LowPt_Zinc2jet,
                         deltaPhi(bestTwoJets.first, bestTwoJets.second),
                         weight);
                    fill(dPhiLeptons_LowPt_Zinc2jet, deltaPhi(leptons[0].v, leptons[1].v), weight);
                    fill(PHI_LowPt_Zinc2jet,
                         PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(BestPHI_LowPt_Zinc2jet,
                         PHI(leptons[0].v, leptons[1].v, bestTwoJets.first, bestTwoJets.second),
                         weight);
                    fill(PHI_T_LowPt_Zinc2jet,
                         PHI_T(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(BestPHI_T_LowPt_Zinc2jet,
                         PHI_T(leptons[0].v, leptons[1].v, bestTwoJets.first, bestTwoJets.second),
                         weight);
                    fill(SpT_LowPt_Zinc2jet,
                         SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(BestSpT_LowPt_Zinc2jet,
                         SpT(leptons[0].v, leptons[1].v, bestTwoJets.first, bestTwoJets.second),
                         weight);
                    fill(SpTJets_LowPt_Zinc2jet, SpTsub(jets[0].v, jets[1].v), weight);
                    fill(BestSpTJets_LowPt_Zinc2jet,
                         SpTsub(bestTwoJets.first, bestTwoJets.second),
                         weight);
                    fill(SpTLeptons_LowPt_Zinc2jet, SpTsub(leptons[0].v, leptons[1].v), weight);
                    fill(SPhi_LowPt_Zinc2jet,
                         SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(BestSPhi_LowPt_Zinc2jet,
                         SPhi(leptons[0].v, leptons[1].v, bestTwoJets.first, bestTwoJets.second),
                         weight);
                    if (SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                        fill(PHI_LowSpT_LowPt_Zinc2jet,
                             PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SPhi_LowSpT_LowPt_Zinc2jet,
                             SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    } else {
                        fill(PHI_HighSpT_LowPt_Zinc2jet,
                             PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SPhi_HighSpT_LowPt_Zinc2jet,
                             SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    }
                    if (SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                        fill(SpT_LowSPhi_LowPt_Zinc2jet,
                             SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    } else {
                        fill(SpT_HighSPhi_LowPt_Zinc2jet,
                             SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    }
                } else {
                    fill(ptBal_HighPt_Zinc2jet, jet1Plus2PlusZ.Pt(), weight);
                    fill(dPhiJets_HighPt_Zinc2jet, deltaPhi(jets[0].v, jets[1].v), weight);
                    fill(dPhiLeptons_HighPt_Zinc2jet, deltaPhi(leptons[0].v, leptons[1].v), weight);
                    fill(PHI_HighPt_Zinc2jet,
                         PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(PHI_T_HighPt_Zinc2jet,
                         PHI_T(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(SpT_HighPt_Zinc2jet,
                         SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(SpTJets_HighPt_Zinc2jet, SpTsub(jets[0].v, jets[1].v), weight);
                    fill(SpTLeptons_HighPt_Zinc2jet, SpTsub(leptons[0].v, leptons[1].v), weight);
                    fill(SPhi_HighPt_Zinc2jet,
                         SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    if (SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                        fill(PHI_LowSpT_HighPt_Zinc2jet,
                             PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SPhi_LowSpT_HighPt_Zinc2jet,
                             SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    } else {
                        fill(PHI_HighSpT_HighPt_Zinc2jet,
                             PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SPhi_HighSpT_HighPt_Zinc2jet,
                             SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    }
                    if (SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                        fill(SpT_LowSPhi_HighPt_Zinc2jet,
                             SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    } else {
                        fill(SpT_HighSPhi_HighPt_Zinc2jet,
                             SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    }
                }
                if (SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                    fill(SpT_LowSPhi_Zinc2jet,
                         SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                } else {
                    fill(SpT_HighSPhi_Zinc2jet,
                         SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                }
                if (SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                    fill(PHI_LowSpT_Zinc2jet,
                         PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(SPhi_LowSpT_Zinc2jet,
                         SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                } else {
                    fill(PHI_HighSpT_Zinc2jet,
                         PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(SPhi_HighSpT_Zinc2jet,
                         SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                }
                if (nGoodJets == 2) {
                    ////////////Special Branch/////////////////////////////////
                    fill(AbsZRapidity_Zexc2jet, fabs(EWKBoson.Rapidity()), weight);
                    fill(AbsSecondJetRapidity_Zexc2jet, fabs(jets[1].v.Rapidity()), weight);
                    fill(SumZSecondJetRapidity_Zexc2jet,
                         fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                         weight);
                    fill(DifZSecondJetRapidity_Zexc2jet,
                         fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                         weight);
                    fill(ZNGoodJetsNVtx_Zexc, 2., EvtVtxCnt, weight);

                    if (EWKBoson.Pt() > 100.) {
                        fill(AbsZRapidity_ZPt100_Zexc2jet, fabs(EWKBoson.Rapidity()), weight);
                        fill(AbsSecondJetRapidity_ZPt100_Zexc2jet,
                             fabs(jets[1].v.Rapidity()),
                             weight);
                        fill(SumZSecondJetRapidity_ZPt100_Zexc2jet,
                             fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                             weight);
                        fill(DifZSecondJetRapidity_ZPt100_Zexc2jet,
                             fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                             weight);
                    }

                    if (EWKBoson.Pt() > 150.) {
                        fill(AbsZRapidity_ZPt150_Zexc2jet, fabs(EWKBoson.Rapidity()), weight);
                        fill(AbsSecondJetRapidity_ZPt150_Zexc2jet,
                             fabs(jets[1].v.Rapidity()),
                             weight);
                        fill(SumZSecondJetRapidity_ZPt150_Zexc2jet,
                             fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                             weight);
                        fill(DifZSecondJetRapidity_ZPt150_Zexc2jet,
                             fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                             weight);
                    }

                    // fill(TruePU_2, EvtPuCntTruth, weight);
                    // fill(PU_2, EvtPuCnt, weight);
                    fill(PU_2, EvtVtxCnt, weight);
                    fill(ZNGoodJets_Zexc_NoWeight, 2.);
                    fill(ZPt_Zexc2jet, EWKBoson.Pt(), weight);
                    fill(ZRapidity_Zexc2jet, EWKBoson.Rapidity(), weight);
                    fill(ZEta_Zexc2jet, EWKBoson.Eta(), weight);
                    fill(SpTLeptons_Zexc2jet, SpTsub(leptons[0].v, leptons[1].v), weight);
                    fill(SecondJetPt_Zexc2jet, jets[1].v.Pt(), weight);
                    fill(SecondJetEta_Zexc2jet, jets[1].v.Eta(), weight);
                    fill(SecondJetPhi_Zexc2jet, jets[1].v.Phi(), weight);

                    //-- DPS Histograms
                    fill(TwoJetsPtDiff_Zexc2jet, jet1Minus2.Pt(), weight);
                    fill(JetsMass_Zexc2jet, jet1Plus2.M(), weight);
                    fill(ptBal_Zexc2jet, jet1Plus2PlusZ.Pt(), weight);
                    fill(dPhiJets_Zexc2jet, deltaPhi(jets[0].v, jets[1].v), weight);
                    fill(dEtaJets_Zexc2jet, jets[0].v.Eta() - jets[1].v.Eta(), weight);
                    fill(dEtaFirstJetZ_Zexc2jet, jets[0].v.Eta() - EWKBoson.Eta(), weight);
                    fill(dEtaSecondJetZ_Zexc2jet, jets[1].v.Eta() - EWKBoson.Eta(), weight);
                    fill(dEtaJet1Plus2Z_Zexc2jet, jet1Plus2.Eta() - EWKBoson.Eta(), weight);
                    fill(PHI_Zexc2jet,
                         PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(PHI_T_Zexc2jet,
                         PHI_T(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(SpT_Zexc2jet,
                         SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    fill(SpTJets_Zexc2jet, SpTsub(jets[0].v, jets[1].v), weight);
                    fill(SPhi_Zexc2jet,
                         SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                         weight);
                    if (EWKBoson.Pt() < 25) {
                        fill(ptBal_LowPt_Zexc2jet, jet1Plus2PlusZ.Pt(), weight);
                        fill(dPhiJets_LowPt_Zexc2jet, deltaPhi(jets[0].v, jets[1].v), weight);
                        fill(PHI_LowPt_Zexc2jet,
                             PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(PHI_T_LowPt_Zexc2jet,
                             PHI_T(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SpT_LowPt_Zexc2jet,
                             SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SpTJets_LowPt_Zexc2jet, SpTsub(jets[0].v, jets[1].v), weight);
                        fill(SpTLeptons_LowPt_Zexc2jet, SpTsub(leptons[0].v, leptons[1].v), weight);
                        fill(SPhi_LowPt_Zexc2jet,
                             SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        if (SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                            fill(PHI_LowSpT_LowPt_Zexc2jet,
                                 PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                            fill(SPhi_LowSpT_LowPt_Zexc2jet,
                                 SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                        } else {
                            fill(PHI_HighSpT_LowPt_Zexc2jet,
                                 PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                            fill(SPhi_HighSpT_LowPt_Zexc2jet,
                                 SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                        }
                        if (SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                            fill(SpT_LowSPhi_LowPt_Zexc2jet,
                                 SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                        } else {
                            fill(SpT_HighSPhi_LowPt_Zexc2jet,
                                 SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                        }
                    } else {
                        fill(ptBal_HighPt_Zexc2jet, jet1Plus2PlusZ.Pt(), weight);
                        fill(dPhiJets_HighPt_Zexc2jet, deltaPhi(jets[0].v, jets[1].v), weight);
                        fill(PHI_HighPt_Zexc2jet,
                             PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(PHI_T_HighPt_Zexc2jet,
                             PHI_T(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SpT_HighPt_Zexc2jet,
                             SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SpTJets_HighPt_Zexc2jet, SpTsub(jets[0].v, jets[1].v), weight);
                        fill(
                            SpTLeptons_HighPt_Zexc2jet, SpTsub(leptons[0].v, leptons[1].v), weight);
                        fill(SPhi_HighPt_Zexc2jet,
                             SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        if (SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                            fill(PHI_LowSpT_HighPt_Zexc2jet,
                                 PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                            fill(SPhi_LowSpT_HighPt_Zexc2jet,
                                 SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                        } else {
                            fill(PHI_HighSpT_HighPt_Zexc2jet,
                                 PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                            fill(SPhi_HighSpT_HighPt_Zexc2jet,
                                 SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                        }
                        if (SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                            fill(SpT_LowSPhi_HighPt_Zexc2jet,
                                 SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                        } else {
                            fill(SpT_HighSPhi_HighPt_Zexc2jet,
                                 SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                                 weight);
                        }
                    }
                    if (SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                        fill(SpT_LowSPhi_Zexc2jet,
                             SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    } else {
                        fill(SpT_HighSPhi_Zexc2jet,
                             SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    }
                    if (SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v) < 0.5) {
                        fill(PHI_LowSpT_Zexc2jet,
                             PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SPhi_LowSpT_Zexc2jet,
                             SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    } else {
                        fill(PHI_HighSpT_Zexc2jet,
                             PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                        fill(SPhi_HighSpT_Zexc2jet,
                             SPhi(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                             weight);
                    }
                }
            }
            if (nGoodJets_20 >= 3) {
                double RatioValue = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        ThirdJetPt_2_Zinc3jetHratio_fit->GetXaxis()->FindBin(jets_20[2].v.Pt());
                    RatioValue = ThirdJetPt_2_Zinc3jetHratio_fit->GetBinContent(binNumber);
                }

                fill(ThirdJetPt_Zinc3jet, jets_20[2].v.Pt(), weight * RatioValue);
                if (nEvents % 2)
                    fill(ThirdJetPt_Zinc3jet_Odd, jets_20[2].v.Pt(), weight * RatioValue);
                else
                    fill(ThirdJetPt_Zinc3jet_Even, jets_20[2].v.Pt(), weight * RatioValue);
                fill(ThirdJetPt_2_Zinc3jet, jets_20[2].v.Pt(), weight);
            }
            if (nGoodJets >= 3) {
                nEventsVInc3Jets++;
                nEffEventsVInc3Jets += weight;
                fill(ZNGoodJets_Zinc, 3., weight);
                fill(ZNGoodJets_Zinc_NoWeight, 3.);

                double RatioValue = 1.;
                double RatioValue1 = 1.;
                double RatioValue2 = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        ThirdJetAbsRapidity_2_Zinc3jetHratio_fit->GetXaxis()->FindBin(
                            fabs(jets[2].v.Eta()));
                    RatioValue = ThirdJetAbsRapidity_2_Zinc3jetHratio_fit->GetBinContent(binNumber);
                    double binNumber1 = JetsHT_2_Zinc3jetHratio_fit->GetXaxis()->FindBin(jetsHT);

                    if (jetsHT >= 90. && jetsHT <= 1200.)
                        RatioValue1 = JetsHT_2_Zinc3jetHratio_fit->GetBinContent(binNumber1);
                    else
                        RatioValue1 = 1.;
                    if (RatioValue1 > 2. || RatioValue1 < 0.5) RatioValue1 = 1.;

                    double binNumber2 = VisPt_2_Zinc3jetQunHratio_fit->GetXaxis()->FindBin(
                        fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()));
                    // RatioValue2 = VisPt_2_Zinc3jetQunHratio_fit->GetBinContent(binNumber2);
                    if (fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()) >= 0. &&
                        fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()) <= 200.)
                        RatioValue2 = VisPt_2_Zinc3jetQunHratio_fit->GetBinContent(binNumber2);
                    else
                        RatioValue2 = 1.;
                    if (RatioValue2 > 2. || RatioValue2 < 0.5) RatioValue2 = 1.;
                }

                // cout << RatioValue1 << "\n";

                fill(ThirdJetEta_Zinc3jet, fabs(jets[2].v.Eta()), weight * RatioValue);
                fill(ThirdJetEta_2_Zinc3jet, fabs(jets[2].v.Eta()), weight);
                fill(ThirdJetAbsRapidity_Zinc3jet, fabs(jets[2].v.Rapidity()), weight * RatioValue);
                fill(ThirdJetAbsRapidity_2_Zinc3jet,
                     fabs(jets[2].v.Rapidity()),
                     weight * RatioValue);

                if (nEvents % 2)
                    fill(ThirdJetAbsRapidity_Zinc3jet_Odd,
                         fabs(jets[2].v.Rapidity()),
                         weight * RatioValue);
                else
                    fill(ThirdJetAbsRapidity_Zinc3jet_Even,
                         fabs(jets[2].v.Rapidity()),
                         weight * RatioValue);
                fill(ThirdJetEtaHigh_Zinc3jet, fabs(jets[2].v.Eta()), weight);
                fill(ThirdJetRapidityHigh_Zinc3jet, fabs(jets[2].v.Rapidity()), weight);
                fill(ThirdJetEtaFull_Zinc3jet, jets[2].v.Eta(), weight);
                fill(ThirdJetPhi_Zinc3jet, jets[2].v.Phi(), weight);
                fill(JetsHT_Zinc3jet, jetsHT, weight * RatioValue1);
                if (nEvents % 2)
                    fill(JetsHT_Zinc3jet_Odd, jetsHT, weight * RatioValue1);
                else
                    fill(JetsHT_Zinc3jet_Even, jetsHT, weight * RatioValue1);
                fill(JetsHT_2_Zinc3jet, jetsHT, weight);
                fill(VisPt_Zinc3jetQun,
                     fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()),
                     weight * RatioValue2);
                fill(VisPt_2_Zinc3jetQun,
                     fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()),
                     weight * RatioValue2);
                if (nEvents % 2)
                    fill(VisPt_Zinc3jetQun_Odd,
                         fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()),
                         weight * RatioValue2);
                else
                    fill(VisPt_Zinc3jetQun_Even,
                         fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()),
                         weight * RatioValue2);

                /// Azimuth cross check
                fill(DPhiZFirstJet_Zinc3jet, fabs(EWKBoson.DeltaPhi(jets[0].v)), weight);
                fill(DPhiZSecondJet_Zinc3jet, fabs(EWKBoson.DeltaPhi(jets[1].v)), weight);
                fill(DPhiZThirdJet_Zinc3jet, fabs(EWKBoson.DeltaPhi(jets[2].v)), weight);
                fill(DPhiFirstSecondJet_Zinc3jet, fabs(jets[0].v.DeltaPhi(jets[1].v)), weight);
                fill(DPhiFirstThirdJet_Zinc3jet, fabs(jets[0].v.DeltaPhi(jets[2].v)), weight);
                fill(DPhiSecondThirdJet_Zinc3jet, fabs(jets[1].v.DeltaPhi(jets[2].v)), weight);

                if (EWKBoson.Pt() > 150.) {
                    fill(DPhiZFirstJet_ZPt150_Zinc3jet, fabs(EWKBoson.DeltaPhi(jets[0].v)), weight);
                    fill(
                        DPhiZSecondJet_ZPt150_Zinc3jet, fabs(EWKBoson.DeltaPhi(jets[1].v)), weight);
                    fill(DPhiZThirdJet_ZPt150_Zinc3jet, fabs(EWKBoson.DeltaPhi(jets[2].v)), weight);
                    fill(DPhiFirstSecondJet_ZPt150_Zinc3jet,
                         fabs(jets[0].v.DeltaPhi(jets[1].v)),
                         weight);
                    fill(DPhiFirstThirdJet_ZPt150_Zinc3jet,
                         fabs(jets[0].v.DeltaPhi(jets[2].v)),
                         weight);
                    fill(DPhiSecondThirdJet_ZPt150_Zinc3jet,
                         fabs(jets[1].v.DeltaPhi(jets[2].v)),
                         weight);
                }

                if (EWKBoson.Pt() > 300.) {
                    fill(DPhiZFirstJet_ZPt300_Zinc3jet, fabs(EWKBoson.DeltaPhi(jets[0].v)), weight);
                    fill(
                        DPhiZSecondJet_ZPt300_Zinc3jet, fabs(EWKBoson.DeltaPhi(jets[1].v)), weight);
                    fill(DPhiZThirdJet_ZPt300_Zinc3jet, fabs(EWKBoson.DeltaPhi(jets[2].v)), weight);
                    fill(DPhiFirstSecondJet_ZPt300_Zinc3jet,
                         fabs(jets[0].v.DeltaPhi(jets[1].v)),
                         weight);
                    fill(DPhiFirstThirdJet_ZPt300_Zinc3jet,
                         fabs(jets[0].v.DeltaPhi(jets[2].v)),
                         weight);
                    fill(DPhiSecondThirdJet_ZPt300_Zinc3jet,
                         fabs(jets[1].v.DeltaPhi(jets[2].v)),
                         weight);
                }

                if (EWKBoson.Pt() > 150. &&
                    (jets[0].v.Pt() + jets[1].v.Pt() + jets[2].v.Pt() > 300.)) {
                    fill(DPhiZFirstJet_ZPt150_HT300_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[0].v)),
                         weight);
                    fill(DPhiZSecondJet_ZPt150_HT300_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[1].v)),
                         weight);
                    fill(DPhiZThirdJet_ZPt150_HT300_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[2].v)),
                         weight);
                }

                if (nGoodJets == 3) {
                    // fill(TruePU_3, EvtPuCntTruth, weight);
                    // fill(PU_3, EvtPuCnt, weight);
                    fill(PU_3, EvtVtxCnt, weight);
                    fill(ZNGoodJets_Zexc_NoWeight, 3.);
                }
            }
            if (nGoodJets_20 >= 4) fill(FourthJetPt_Zinc4jet, jets_20[3].v.Pt(), weight);
            if (nGoodJets >= 4) {
                fill(ZNGoodJets_Zinc, 4., weight);
                fill(ZNGoodJets_Zinc_NoWeight, 4.);
                fill(FourthJetEta_Zinc4jet, fabs(jets[3].v.Eta()), weight);
                fill(FourthJetAbsRapidity_Zinc4jet, fabs(jets[3].v.Rapidity()), weight);
                fill(FourthJetEtaHigh_Zinc4jet, fabs(jets[3].v.Eta()), weight);
                fill(FourthJetRapidityHigh_Zinc4jet, fabs(jets[3].v.Rapidity()), weight);
                fill(FourthJetEtaFull_Zinc4jet, jets[3].v.Eta(), weight);
                fill(FourthJetPhi_Zinc4jet, jets[3].v.Phi(), weight);
                fill(JetsHT_Zinc4jet, jetsHT, weight);
                if (nGoodJets == 4) {
                    // fill(TruePU_4, EvtPuCntTruth, weight);
                    // fill(PU_4, EvtPuCnt, weight);
                    fill(PU_4, EvtVtxCnt, weight);
                    fill(ZNGoodJets_Zexc_NoWeight, 4.);
                }
            }
            if (nGoodJets_20 >= 5) fill(FifthJetPt_Zinc5jet, jets_20[4].v.Pt(), weight);
            if (nGoodJets >= 5) {
                fill(ZNGoodJets_Zinc, 5., weight);
                fill(ZNGoodJets_Zinc_NoWeight, 5.);
                fill(FifthJetEta_Zinc5jet, fabs(jets[4].v.Eta()), weight);
                fill(FifthJetAbsRapidity_Zinc5jet, fabs(jets[4].v.Rapidity()), weight);
                fill(FifthJetEtaHigh_Zinc5jet, fabs(jets[4].v.Eta()), weight);
                fill(FifthJetRapidityHigh_Zinc5jet, fabs(jets[4].v.Rapidity()), weight);
                fill(FifthJetEtaFull_Zinc5jet, jets[4].v.Eta(), weight);
                fill(FifthJetPhi_Zinc5jet, jets[4].v.Phi(), weight);
                fill(JetsHT_Zinc5jet, jetsHT, weight);
                if (nGoodJets == 5) {
                    // fill(TruePU_5, EvtPuCntTruth, weight);
                    // fill(PU_5, EvtPuCnt, weight);
                    fill(PU_5, EvtVtxCnt, weight);
                    fill(ZNGoodJets_Zexc_NoWeight, 5.);
                }
            }
            if (nGoodJets_20 >= 6) fill(SixthJetPt_Zinc6jet, jets_20[5].v.Pt(), weight);
            if (nGoodJets >= 6) {
                fill(ZNGoodJets_Zinc, 6., weight);
                fill(ZNGoodJets_Zinc_NoWeight, 6.);
                fill(SixthJetEta_Zinc6jet, fabs(jets[5].v.Eta()), weight);
                fill(SixthJetEtaHigh_Zinc6jet, fabs(jets[5].v.Eta()), weight);
                fill(SixthJetRapidityHigh_Zinc6jet, fabs(jets[5].v.Rapidity()), weight);
                fill(SixthJetEtaFull_Zinc6jet, jets[5].v.Eta(), weight);
                fill(SixthJetPhi_Zinc6jet, jets[5].v.Phi(), weight);
                fill(JetsHT_Zinc6jet, jetsHT, weight);
                if (nGoodJets == 6) {
                    // fill(TruePU_6, EvtPuCntTruth, weight);
                    // fill(PU_6, EvtPuCnt, weight);
                    fill(PU_6, EvtVtxCnt, weight);
                    fill(ZNGoodJets_Zexc_NoWeight, 6.);
                }
            }
            if (nGoodJets >= 7) {
                fill(ZNGoodJets_Zinc, 7., weight);
                if (nGoodJets == 7) {
                    // fill(TruePU_7, EvtPuCntTruth, weight);
                    // fill(PU_7, EvtPuCnt, weight);
                    fill(PU_7, EvtVtxCnt, weight);
                }
            }
            if (nGoodJets >= 8) {
                fill(ZNGoodJets_Zinc, 8., weight);
            }

            if (nGoodJets >= 1) {

                double RatioValue = 1.;
                double RatioValue1 = 1.;
                double RatioValue2 = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        JZB_2Hratio_fit->GetXaxis()->FindBin(hadronicR.Pt() - EWKBoson.Pt());
                    RatioValue = JZB_2Hratio_fit->GetBinContent(binNumber);
                    double binNumber1 =
                        JZB_ptLow_2Hratio_fit->GetXaxis()->FindBin(hadronicR.Pt() - EWKBoson.Pt());
                    // RatioValue1 =  JZB_ptLow_2Hratio_fit->GetBinContent(binNumber1);
                    if ((hadronicR.Pt() - EWKBoson.Pt()) > -50 &&
                        (hadronicR.Pt() - EWKBoson.Pt()) < 200.)
                        RatioValue1 = JZB_ptLow_2Hratio_fit->GetBinContent(binNumber1);
                    else
                        RatioValue1 = 1.;
                    if (RatioValue1 > 2. || RatioValue1 < 0.5) RatioValue1 = 1.;

                    double binNumber2 =
                        JZB_ptHigh_2Hratio_fit->GetXaxis()->FindBin(hadronicR.Pt() - EWKBoson.Pt());
                    RatioValue2 = JZB_ptHigh_2Hratio_fit->GetBinContent(binNumber2);
                }

                fill(HadRecoil, hadronicR.Pt(), weight);
                fill(JZB, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue);
                fill(JZB_2, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue);
                if (nEvents % 2)
                    fill(JZB_Odd, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue);
                else
                    fill(JZB_Even, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue);

                if (EWKBoson.Pt() <= 50) {
                    fill(JZB_ptLow, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue1);
                    fill(JZB_ptLow_2, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue1);
                    if (nEvents % 2)
                        fill(JZB_ptLow_Odd, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue1);
                    else
                        fill(JZB_ptLow_Even, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue1);

                } else {
                    fill(JZB_ptHigh, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue2);
                    fill(JZB_ptHigh_2, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue2);
                    if (nEvents % 2)
                        fill(JZB_ptHigh_Odd, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue2);
                    else
                        fill(JZB_ptHigh, hadronicR.Pt() - EWKBoson.Pt(), weight * RatioValue2);
                }
            }
            //=======================================================================================================//
        }

        //=======================================================================================================//
        //             Unfolding               //
        //====================================//

        if (hasRecoInfo && hasGenInfo && (!bTagJetFound || !rejectBTagEvents)) {
            if (nGoodJets_20 >= 1 && nGoodGenJets_20 >= 1) {
                // looks for matching gen jet:
                int igen = 0;
                double dr2 = 999.;
                for (int i = 0; i < nGoodGenJets_20; ++i) {
                    double dr2_ = std::pow(jets_20[0].v.Eta() - genJets_20[i].v.Eta(), 2) +
                                  std::pow(jets_20[0].v.Pt() - genJets_20[i].v.Pt(), 2);
                    if (dr2_ < dr2) {
                        dr2 = dr2_;
                        igen = i;
                    }
                }
                fill(FirstJetPtRecoOvGen_Zinc1jet_NVtx,
                     jets_20[0].v.Pt() / genJets_20[igen].v.Pt(),
                     EvtVtxCnt,
                     weight);
            }

            //-- EWKBoson Mass and jet multiplicity

            if (passesgenLeptonCutNoMass && passesLeptonCutNoMass) {

                if (genEWKBoson.M() > 115. && genEWKBoson.M() < 135. && EWKBoson.M() > 115. &&
                    EWKBoson.M() < 135.)
                    fill(hresponseZPt_Zinc0jetM115_135, EWKBoson.Pt(), genEWKBoson.Pt(), weight);

                if (genEWKBoson.M() >= 111. && genEWKBoson.M() < 130. && EWKBoson.M() >= 111. &&
                    EWKBoson.M() < 130.)
                    fill(hresponseZPt_Zinc0jetM111_130, EWKBoson.Pt(), genEWKBoson.Pt(), weight);

                if (genEWKBoson.M() >= 130. && genEWKBoson.M() <= 3000. && EWKBoson.M() >= 130. &&
                    EWKBoson.M() <= 3000.)
                    fill(hresponseZPt_Zinc0jetM130_3, EWKBoson.Pt(), genEWKBoson.Pt(), weight);
            }

            if (passesgenLeptonCut && passesLeptonCut) {

                double RatioValue = 1.;
                if (UnfoldUnc) {
                    double binNumber = ZNGoodJets_ZexcHratio_fit->GetXaxis()->FindBin(nGoodJets);
                    RatioValue = ZNGoodJets_ZexcHratio_fit->GetBinContent(binNumber);
                }

                fill(hresponseZNGoodJets_Zexc, nGoodJets, nGoodGenJets, weight * RatioValue);
                fill(hresponseZPt_Zinc0jet, EWKBoson.Pt(), genEWKBoson.Pt(), weight);
                fill(hresponseVisPt_Zinc0jetQun, EWKBoson.Pt(), genEWKBoson.Pt(), weight);
            }

            //-- First Jet Pt
            if (nGoodGenJets >= 1 && passesgenLeptonCut && nGoodJets >= 1 && passesLeptonCut) {

                fill(hresponseZPt_Zinc1jet, EWKBoson.Pt(), genEWKBoson.Pt(), weight);
                fill(hresponseZAbsRapidity_Zinc1jet,
                     fabs(EWKBoson.Rapidity()),
                     fabs(genEWKBoson.Rapidity()),
                     weight);

                double RatioValue = 1.;
                double RatioValue1 = 1.;
                double RatioValue2 = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        FirstJetAbsRapidity_2_Zinc1jetHratio_fit->GetXaxis()->FindBin(
                            fabs(jets[0].v.Eta()));
                    RatioValue = FirstJetAbsRapidity_2_Zinc1jetHratio_fit->GetBinContent(binNumber);
                    double binNumber1 = JetsHT_2_Zinc1jetHratio_fit->GetXaxis()->FindBin(jetsHT);
                    RatioValue1 = JetsHT_2_Zinc1jetHratio_fit->GetBinContent(binNumber1);
                    double binNumber2 = VisPt_2_Zinc1jetQunHratio_fit->GetXaxis()->FindBin(
                        fabs((jets[0].v + EWKBoson).Pt()));
                    RatioValue2 = VisPt_2_Zinc1jetQunHratio_fit->GetBinContent(binNumber2);
                }

                fill(hresponseFirstJetEta_Zinc1jet,
                     fabs(jets[0].v.Eta()),
                     fabs(genJets[0].v.Eta()),
                     weight * RatioValue);
                fill(hresponseFirstJetAbsRapidity_Zinc1jet,
                     fabs(jets[0].v.Rapidity()),
                     fabs(genJets[0].v.Rapidity()),
                     weight * RatioValue);
                fill(hresponseFirstJetEtaHigh_Zinc1jet,
                     fabs(jets[0].v.Eta()),
                     fabs(genJets[0].v.Eta()),
                     weight);
                fill(hresponseFirstJetRapidityHigh_Zinc1jet,
                     fabs(jets[0].v.Rapidity()),
                     fabs(genJets[0].v.Rapidity()),
                     weight);
                fill(hresponseJetsHT_Zinc1jet, jetsHT, genJetsHT, weight * RatioValue1);
                fill(hresponseVisPt_Zinc1jetQun,
                     fabs((jets[0].v + EWKBoson).Pt()),
                     fabs((genJets[0].v + genEWKBoson).Pt()),
                     weight * RatioValue2);

                // Additional Abs responses of variables
                fill(hresponseAbsZRapidity_Zinc1jet,
                     fabs(EWKBoson.Rapidity()),
                     fabs(genEWKBoson.Rapidity()),
                     weight);
                fill(hresponseAbsFirstJetRapidity_Zinc1jet,
                     fabs(jets[0].v.Rapidity()),
                     fabs(genJets[0].v.Rapidity()),
                     weight);
                fill(hresponseSumZFirstJetRapidity_Zinc1jet,
                     fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                     weight);
                fill(hresponseDifZFirstJetRapidity_Zinc1jet,
                     fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                     weight);

                // cross check//////
                fill(hresponseSumZFirstJetEta_Zinc1jet,
                     fabs(EWKBoson.Eta() + jets[0].v.Eta()) / 2.0,
                     fabs(genEWKBoson.Eta() + genJets[0].v.Eta()) / 2.0,
                     weight);
                fill(hresponseDifZFirstJetEta_Zinc1jet,
                     fabs(EWKBoson.Eta() - jets[0].v.Eta()) / 2.0,
                     fabs(genEWKBoson.Eta() - genJets[0].v.Eta()) / 2.0,
                     weight);

                /////Azimuthal cross check////////////////////////////
                fill(hresponseDPhiZFirstJet_Zinc1jet,
                     fabs(EWKBoson.DeltaPhi(jets[0].v)),
                     fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                     weight);

                if (genEWKBoson.Pt() > 100. && EWKBoson.Pt() > 100.) {
                    fill(hresponseAbsZRapidity_ZPt100_Zinc1jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsFirstJetRapidity_ZPt100_Zinc1jet,
                         fabs(jets[0].v.Rapidity()),
                         fabs(genJets[0].v.Rapidity()),
                         weight);
                    fill(hresponseSumZFirstJetRapidity_ZPt100_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZFirstJetRapidity_ZPt100_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                if (genEWKBoson.Pt() > 150. && EWKBoson.Pt() > 150.) {
                    fill(hresponseAbsZRapidity_ZPt150_Zinc1jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsFirstJetRapidity_ZPt150_Zinc1jet,
                         fabs(jets[0].v.Rapidity()),
                         fabs(genJets[0].v.Rapidity()),
                         weight);
                    fill(hresponseSumZFirstJetRapidity_ZPt150_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZFirstJetRapidity_ZPt150_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         weight);

                    fill(hresponseDPhiZFirstJet_ZPt150_Zinc1jet,
                         fabs(EWKBoson.DeltaPhi(jets[0].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         weight);
                }

                if (genEWKBoson.Pt() > 300. && EWKBoson.Pt() > 300.) {
                    fill(hresponseAbsZRapidity_ZPt300_Zinc1jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsFirstJetRapidity_ZPt300_Zinc1jet,
                         fabs(jets[0].v.Rapidity()),
                         fabs(genJets[0].v.Rapidity()),
                         weight);
                    fill(hresponseSumZFirstJetRapidity_ZPt300_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZFirstJetRapidity_ZPt300_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         weight);

                    fill(hresponseDPhiZFirstJet_ZPt300_Zinc1jet,
                         fabs(EWKBoson.DeltaPhi(jets[0].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         weight);
                }

                /// different JetPt Cuts///////

                if (genJets[0].v.Pt() > 50. && jets[0].v.Pt() > 50.) {
                    fill(hresponseAbsZRapidity_FirstJetPt50_Zinc1jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsFirstJetRapidity_FirstJetPt50_Zinc1jet,
                         fabs(jets[0].v.Rapidity()),
                         fabs(genJets[0].v.Rapidity()),
                         weight);
                    fill(hresponseSumZFirstJetRapidity_FirstJetPt50_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZFirstJetRapidity_FirstJetPt50_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                if (genJets[0].v.Pt() > 80. && jets[0].v.Pt() > 80.) {
                    fill(hresponseAbsZRapidity_FirstJetPt80_Zinc1jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsFirstJetRapidity_FirstJetPt80_Zinc1jet,
                         fabs(jets[0].v.Rapidity()),
                         fabs(genJets[0].v.Rapidity()),
                         weight);
                    fill(hresponseSumZFirstJetRapidity_FirstJetPt80_Zinc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZFirstJetRapidity_FirstJetPt80_Zinc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                fill(hresponseSumZJetRapidity_Zinc1jet,
                     0.5 * fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()),
                     0.5 * fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()),
                     weight);
                fill(hresponseDifZJetRapidity_Zinc1jet,
                     0.5 * fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()),
                     0.5 * fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()),
                     weight);

                for (unsigned short i(0); i < 5; i++) {
                    if (EWKBoson.Pt() > ZptRange[i] && EWKBoson.Pt() <= ZptRange[i + 1]) {
                        fill(hresponsetau_sum_Zinc1jet[i], tau_sum, gentau_sum, weight);
                        fill(hresponsetau_max_Zinc1jet[i], tau_max, gentau_max, weight);
                        fill(hresponsetau_c_sum_Zinc1jet[i], tau_c_sum, gentau_c_sum, weight);
                        fill(hresponsetau_c_max_Zinc1jet[i], tau_c_max, gentau_c_max, weight);
                        fill(hresponsetau_cm_sum_Zinc1jet[i], tau_cm_sum, gentau_cm_sum, weight);
                        fill(hresponsetau_cm_max_Zinc1jet[i], tau_cm_max, gentau_cm_max, weight);
                        fill(hresponsetau_c_cm_sum_Zinc1jet[i],
                             tau_c_cm_sum,
                             gentau_c_cm_sum,
                             weight);
                        fill(hresponsetau_c_cm_max_Zinc1jet[i],
                             tau_c_cm_max,
                             gentau_c_cm_max,
                             weight);
                    }
                }
            }

            if (nGoodGenJets_20 >= 1 && passesgenLeptonCut && nGoodJets_20 >= 1 &&
                passesLeptonCut) {

                double RatioValue = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        FirstJetPt_2_Zinc1jetHratio_fit->GetXaxis()->FindBin(jets_20[0].v.Pt());
                    RatioValue = FirstJetPt_2_Zinc1jetHratio_fit->GetBinContent(binNumber);
                }

                fill(hresponseFirstJetPt_Zinc1jet,
                     jets_20[0].v.Pt(),
                     genJets_20[0].v.Pt(),
                     weight * RatioValue);
                if (hresponseFirstJetPt_Zinc1jet) {
                    fill(hresponseFirstJetPtEta_Zinc1jet,
                         0.5 +
                             FirstJetPtEta_Zinc1jet->FindBin(jets_20[0].v.Pt(),
                                                             fabs(jets_20[0].v.Eta())),
                         0.5 +
                             FirstJetPtEta_Zinc1jet->FindBin(genJets_20[0].v.Pt(),
                                                             fabs(genJets_20[0].v.Eta())),
                         weight);
                }
            }

            // exclusive one jet case
            if (nGoodGenJets == 1 && passesgenLeptonCut && nGoodJets == 1 && passesLeptonCut) {
                fill(hresponseAbsZRapidity_Zexc1jet,
                     fabs(EWKBoson.Rapidity()),
                     fabs(genEWKBoson.Rapidity()),
                     weight);
                fill(hresponseAbsJetRapidity_Zexc1jet,
                     fabs(jets[0].v.Rapidity()),
                     fabs(genJets[0].v.Rapidity()),
                     weight);
                fill(hresponseSumZJetRapidity_Zexc1jet,
                     fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                     weight);
                fill(hresponseDifZJetRapidity_Zexc1jet,
                     fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                     weight);

                if (EWKBoson.Pt() > 100. && genEWKBoson.Pt() > 100.) {
                    fill(hresponseAbsZRapidity_ZPt100_Zexc1jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsJetRapidity_ZPt100_Zexc1jet,
                         fabs(jets[0].v.Rapidity()),
                         fabs(genJets[0].v.Rapidity()),
                         weight);
                    fill(hresponseSumZJetRapidity_ZPt100_Zexc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZJetRapidity_ZPt100_Zexc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                if (EWKBoson.Pt() > 150. && genEWKBoson.Pt() > 150.) {
                    fill(hresponseAbsZRapidity_ZPt150_Zexc1jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsJetRapidity_ZPt150_Zexc1jet,
                         fabs(jets[0].v.Rapidity()),
                         fabs(genJets[0].v.Rapidity()),
                         weight);
                    fill(hresponseSumZJetRapidity_ZPt150_Zexc1jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZJetRapidity_ZPt150_Zexc1jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         weight);
                }
            }

            //-- Second Jet Pt inclusive
            if (nGoodGenJets >= 2 && passesgenLeptonCut && nGoodJets >= 2 && passesLeptonCut) {
                ////////////////////////Special Branch//////////////////////
                fill(hresponseZPt_Zinc2jet, EWKBoson.Pt(), genEWKBoson.Pt(), weight);

                fill(hresponseAbsFirstJetRapidity_Zinc2jet,
                     fabs(jets[0].v.Rapidity()),
                     fabs(genJets[0].v.Rapidity()),
                     weight);
                fill(hresponseSumZFirstJetRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                     weight);
                fill(hresponseDifZFirstJetRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                     weight);

                fill(hresponseAbsZRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity()),
                     fabs(genEWKBoson.Rapidity()),
                     weight);
                fill(hresponseAbsSecondJetRapidity_Zinc2jet,
                     fabs(jets[1].v.Rapidity()),
                     fabs(genJets[1].v.Rapidity()),
                     weight);
                fill(hresponseSumZSecondJetRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                     weight);
                fill(hresponseDifZSecondJetRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                     weight);

                fill(hresponseSumFirstSecondJetRapidity_Zinc2jet,
                     fabs(jets[0].v.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                     fabs(genJets[0].v.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                     weight);
                fill(hresponseDifFirstSecondJetRapidity_Zinc2jet,
                     fabs(jets[0].v.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                     fabs(genJets[0].v.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                     weight);

                TLorentzVector genDiJets = genJets[0].v + genJets[1].v;
                TLorentzVector DiJets = jets[0].v + jets[1].v;
                fill(hresponseSumZTwoJetsRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() + DiJets.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() + genDiJets.Rapidity()) / 2.0,
                     weight);
                fill(hresponseDifZTwoJetsRapidity_Zinc2jet,
                     fabs(EWKBoson.Rapidity() - DiJets.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() - genDiJets.Rapidity()) / 2.0,
                     weight);

                /////Azimuthal cross check//////////////////////////////
                fill(hresponseDPhiZFirstJet_Zinc2jet,
                     fabs(EWKBoson.DeltaPhi(jets[0].v)),
                     fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                     weight);
                fill(hresponseDPhiZSecondJet_Zinc2jet,
                     fabs(EWKBoson.DeltaPhi(jets[1].v)),
                     fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                     weight);
                fill(hresponseDPhiFirstSecondJet_Zinc2jet,
                     fabs(jets[0].v.DeltaPhi(jets[1].v)),
                     fabs(genJets[0].v.DeltaPhi(genJets[1].v)),
                     weight);

                if (EWKBoson.Pt() > 100. && genEWKBoson.Pt() > 100.) {
                    fill(hresponseAbsZRapidity_ZPt100_Zinc2jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsSecondJetRapidity_ZPt100_Zinc2jet,
                         fabs(jets[1].v.Rapidity()),
                         fabs(genJets[1].v.Rapidity()),
                         weight);
                    fill(hresponseSumZSecondJetRapidity_ZPt100_Zinc2jet,
                         fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZSecondJetRapidity_ZPt100_Zinc2jet,
                         fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                         weight);
                }

                if (EWKBoson.Pt() > 150. && genEWKBoson.Pt() > 150.) {
                    fill(hresponseAbsZRapidity_ZPt150_Zinc2jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsSecondJetRapidity_ZPt150_Zinc2jet,
                         fabs(jets[1].v.Rapidity()),
                         fabs(genJets[1].v.Rapidity()),
                         weight);
                    fill(hresponseSumZSecondJetRapidity_ZPt150_Zinc2jet,
                         fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZSecondJetRapidity_ZPt150_Zinc2jet,
                         fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                         weight);

                    fill(hresponseDPhiZFirstJet_ZPt150_Zinc2jet,
                         fabs(EWKBoson.DeltaPhi(jets[0].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         weight);
                }

                if (EWKBoson.Pt() > 300. && genEWKBoson.Pt() > 300.) {
                    fill(hresponseDPhiZFirstJet_ZPt300_Zinc2jet,
                         fabs(EWKBoson.DeltaPhi(jets[0].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         weight);
                }

                // set jet rapidity discriminator////

                if (fabs(genJets[0].v.Rapidity() - genJets[1].v.Rapidity()) > 2 &&
                    fabs(jets[0].v.Rapidity() - jets[1].v.Rapidity()) > 2) {
                    fill(hresponseAbsZRapidity_DifJetRapidityl2_Zinc2jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet,
                         fabs(jets[0].v.Rapidity()),
                         fabs(genJets[0].v.Rapidity()),
                         weight);
                    fill(hresponseSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                if (fabs(genJets[0].v.Rapidity() - genJets[1].v.Rapidity()) < 2 &&
                    fabs(jets[0].v.Rapidity() - jets[1].v.Rapidity()) < 2) {
                    fill(hresponseAbsZRapidity_DifJetRapiditys2_Zinc2jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet,
                         fabs(jets[0].v.Rapidity()),
                         fabs(genJets[0].v.Rapidity()),
                         weight);
                    fill(hresponseSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet,
                         fabs(EWKBoson.Rapidity() + jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[0].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet,
                         fabs(EWKBoson.Rapidity() - jets[0].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[0].v.Rapidity()) / 2.0,
                         weight);
                }

                double RatioValue = 1.;
                double RatioValue1 = 1.;
                double RatioValue2 = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        SecondJetAbsRapidity_2_Zinc2jetHratio_fit->GetXaxis()->FindBin(
                            fabs(jets[1].v.Eta()));
                    RatioValue =
                        SecondJetAbsRapidity_2_Zinc2jetHratio_fit->GetBinContent(binNumber);
                    double binNumber1 = JetsHT_2_Zinc2jetHratio_fit->GetXaxis()->FindBin(jetsHT);
                    RatioValue1 = JetsHT_2_Zinc2jetHratio_fit->GetBinContent(binNumber1);
                    double binNumber2 = VisPt_2_Zinc2jetQunHratio_fit->GetXaxis()->FindBin(
                        fabs((jets[0].v + jets[1].v + EWKBoson).Pt()));
                    RatioValue2 = VisPt_2_Zinc2jetQunHratio_fit->GetBinContent(binNumber2);
                }

                fill(hresponseSecondJetEta_Zinc2jet,
                     fabs(jets[1].v.Eta()),
                     fabs(genJets[1].v.Eta()),
                     weight * RatioValue);
                fill(hresponseSecondJetAbsRapidity_Zinc2jet,
                     fabs(jets[1].v.Rapidity()),
                     fabs(genJets[1].v.Rapidity()),
                     weight * RatioValue);
                fill(hresponseSecondJetEtaHigh_Zinc2jet,
                     fabs(jets[1].v.Eta()),
                     fabs(genJets[1].v.Eta()),
                     weight);
                fill(hresponseSecondJetRapidityHigh_Zinc2jet,
                     fabs(jets[1].v.Rapidity()),
                     fabs(genJets[1].v.Rapidity()),
                     weight);
                fill(hresponseJetsHT_Zinc2jet, jetsHT, genJetsHT, weight * RatioValue1);
                fill(hresponseVisPt_Zinc2jetQun,
                     fabs((jets[0].v + jets[1].v + EWKBoson).Pt()),
                     fabs((genJets[0].v + genJets[1].v + genEWKBoson).Pt()),
                     weight * RatioValue2);
                // fill(hresponseJetsHT_2_Zinc2jet, jetsHT, genJetsHT, weight);
                // fill(responseTwoJetsPtDiffInc, jet1Minus2.Pt(), genJet1Minus2.Pt(), weight);
                // fill(responseBestTwoJetsPtDiffInc, bestJet1Minus2.Pt(), genBestJet1Minus2.Pt(),
                // weight);
                // fill(responseJetsMassInc, jet1Plus2.M(), genJet1Plus2.M(), weight);
                fill(hresponseJetsMass_Zinc2jet, jet1Plus2.M(), genJet1Plus2.M(), weight);

                if (EvtVtxCnt < 14)
                    fill(hresponseJetsMassLowPU_Zinc2jet, jet1Plus2.M(), genJet1Plus2.M(), weight);
                else if (EvtVtxCnt < 18)
                    fill(hresponseJetsMassMidPU_Zinc2jet, jet1Plus2.M(), genJet1Plus2.M(), weight);
                else
                    fill(hresponseJetsMassHigPU_Zinc2jet, jet1Plus2.M(), genJet1Plus2.M(), weight);

                // fill(responseBestJetsMassInc, bestJet1Plus2.M(), genBestJet1Plus2.M(), weight);
                // fill(responseSpTJets_Zinc2jet, SpTsub(jets[0].v, jets[1].v), SpTsub(genJets[0].v,
                // genJets[1].v), weight);
                // fill(responseBestSpTJets_Zinc2jet, SpTsub(bestTwoJets.first, bestTwoJets.second),
                // SpTsub(genBestTwoJets.first, genBestTwoJets.second), weight);
                // fill(responseSpT_Zinc2jet, SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                // SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v), weight);
                // fill(responseBestSpT_Zinc2jet, SpT(leptons[0].v, leptons[1].v, bestTwoJets.first,
                // bestTwoJets.second), SpT(genLeptons[0].v, genLeptons[1].v, genBestTwoJets.first,
                // genBestTwoJets.second), weight);
                // fill(responsedPhiJets_Zinc2jet, deltaPhi(jets[0].v, jets[1].v),
                // deltaPhi(genJets[0].v, genJets[1].v), weight);
                // fill(responseBestdPhiJets_Zinc2jet, deltaPhi(bestTwoJets.first,
                // bestTwoJets.second), deltaPhi(genBestTwoJets.first, genBestTwoJets.second),
                // weight);
                // fill(responsePHI_Zinc2jet, PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                // PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v), weight);
                // fill(responseBestPHI_Zinc2jet, PHI(leptons[0].v, leptons[1].v, bestTwoJets.first,
                // bestTwoJets.second), PHI(genLeptons[0].v, genLeptons[1].v, genBestTwoJets.first,
                // genBestTwoJets.second), weight);
                // fill(responsePHI_T_Zinc2jet, PHI_T(leptons[0].v, leptons[1].v, jets[0].v,
                // jets[1].v), PHI_T(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                // weight);
                // fill(responseBestPHI_T_Zinc2jet, PHI_T(leptons[0].v, leptons[1].v,
                // bestTwoJets.first, bestTwoJets.second), PHI_T(genLeptons[0].v, genLeptons[1].v,
                // genBestTwoJets.first, genBestTwoJets.second), weight);
                // fill(responseSPhi_Zinc2jet, SPhi(leptons[0].v, leptons[1].v, jets[0].v,
                // jets[1].v), SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                // weight);
                // fill(responsedEtaJets_Zinc2jet,
                // fabs(genJets[0].eta-genJets[1].eta),fabs(jets[0].eta-jets[1].eta), weight);
                // fill(responseBestSPhi_Zinc2jet, SPhi(leptons[0].v, leptons[1].v,
                // bestTwoJets.first, bestTwoJets.second), SPhi(genLeptons[0].v, genLeptons[1].v,
                // genBestTwoJets.first, genBestTwoJets.second), weight);
            }

            if (nGoodGenJets_20 >= 2 && passesgenLeptonCut && nGoodJets_20 >= 2 &&
                passesLeptonCut) {
                double RatioValue = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        SecondJetPt_2_Zinc2jetHratio_fit->GetXaxis()->FindBin(jets_20[1].v.Pt());
                    RatioValue = SecondJetPt_2_Zinc2jetHratio_fit->GetBinContent(binNumber);
                }

                fill(hresponseSecondJetPt_Zinc2jet,
                     jets_20[1].v.Pt(),
                     genJets_20[1].v.Pt(),
                     weight * RatioValue);
            }

            //-- Second Jet Pt exclusive
            if (nGoodGenJets == 2 && passesgenLeptonCut && nGoodJets == 2 && passesLeptonCut) {
                //////////////////////////////Special Branch/////////////////
                fill(hresponseAbsZRapidity_Zexc2jet,
                     fabs(EWKBoson.Rapidity()),
                     fabs(genEWKBoson.Rapidity()),
                     weight);
                fill(hresponseAbsSecondJetRapidity_Zexc2jet,
                     fabs(jets[1].v.Rapidity()),
                     fabs(genJets[1].v.Rapidity()),
                     weight);
                fill(hresponseSumZSecondJetRapidity_Zexc2jet,
                     fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                     weight);
                fill(hresponseDifZSecondJetRapidity_Zexc2jet,
                     fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                     fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                     weight);

                if (EWKBoson.Pt() > 100. && genEWKBoson.Pt() > 100.) {
                    fill(hresponseAbsZRapidity_ZPt100_Zexc2jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsSecondJetRapidity_ZPt100_Zexc2jet,
                         fabs(jets[1].v.Rapidity()),
                         fabs(genJets[1].v.Rapidity()),
                         weight);
                    fill(hresponseSumZSecondJetRapidity_ZPt100_Zexc2jet,
                         fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZSecondJetRapidity_ZPt100_Zexc2jet,
                         fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                         weight);
                }

                if (EWKBoson.Pt() > 150. && genEWKBoson.Pt() > 150.) {
                    fill(hresponseAbsZRapidity_ZPt150_Zexc2jet,
                         fabs(EWKBoson.Rapidity()),
                         fabs(genEWKBoson.Rapidity()),
                         weight);
                    fill(hresponseAbsSecondJetRapidity_ZPt150_Zexc2jet,
                         fabs(jets[1].v.Rapidity()),
                         fabs(genJets[1].v.Rapidity()),
                         weight);
                    fill(hresponseSumZSecondJetRapidity_ZPt150_Zexc2jet,
                         fabs(EWKBoson.Rapidity() + jets[1].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() + genJets[1].v.Rapidity()) / 2.0,
                         weight);
                    fill(hresponseDifZSecondJetRapidity_ZPt150_Zexc2jet,
                         fabs(EWKBoson.Rapidity() - jets[1].v.Rapidity()) / 2.0,
                         fabs(genEWKBoson.Rapidity() - genJets[1].v.Rapidity()) / 2.0,
                         weight);
                }

                //

                // fill(responseTwoJetsPtDiffExc, jet1Minus2.Pt(), genJet1Minus2.Pt(), weight);
                // fill(responseJetsMassExc, jet1Plus2.M(), genJet1Plus2.M(), weight);
                // fill(responseSpTJets_Zexc2jet, SpTsub(jets[0].v, jets[1].v), SpTsub(genJets[0].v,
                // genJets[1].v), weight);
                // fill(responseSpT_Zexc2jet, SpT(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                // SpT(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v), weight);
                // fill(responsedPhiJets_Zexc2jet, deltaPhi(jets[0].v, jets[1].v),
                // deltaPhi(genJets[0].v, genJets[1].v), weight);
                // fill(responsePHI_Zexc2jet, PHI(leptons[0].v, leptons[1].v, jets[0].v, jets[1].v),
                // PHI(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v), weight);
                // fill(responsePHI_T_Zexc2jet, PHI_T(leptons[0].v, leptons[1].v, jets[0].v,
                // jets[1].v), PHI_T(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                // weight);
                // fill(responsedEtaJets_Zexc2jet,
                // fabs(genJets[0].eta-genJets[1].eta),fabs(jets[0].eta-jets[1].eta), weight);
                // fill(responseSPhi_Zexc2jet, SPhi(leptons[0].v, leptons[1].v, jets[0].v,
                // jets[1].v), SPhi(genLeptons[0].v, genLeptons[1].v, genJets[0].v, genJets[1].v),
                // weight);
            }

            //-- Third Jet Pt
            if (nGoodGenJets >= 3 && passesgenLeptonCut && nGoodJets >= 3 && passesLeptonCut) {

                double RatioValue = 1.;
                double RatioValue1 = 1.;
                double RatioValue2 = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        ThirdJetAbsRapidity_2_Zinc3jetHratio_fit->GetXaxis()->FindBin(
                            fabs(jets[2].v.Eta()));
                    RatioValue = ThirdJetAbsRapidity_2_Zinc3jetHratio_fit->GetBinContent(binNumber);
                    double binNumber1 = JetsHT_2_Zinc3jetHratio_fit->GetXaxis()->FindBin(jetsHT);
                    // RatioValue1 = JetsHT_2_Zinc3jetHratio_fit->GetBinContent(binNumber1);
                    if (jetsHT >= 90. && jetsHT <= 1200.)
                        RatioValue1 = JetsHT_2_Zinc3jetHratio_fit->GetBinContent(binNumber1);
                    else
                        RatioValue1 = 1.;
                    if (RatioValue1 > 2. || RatioValue1 < 0.5) RatioValue1 = 1.;

                    double binNumber2 = VisPt_2_Zinc3jetQunHratio_fit->GetXaxis()->FindBin(
                        fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()));
                    // RatioValue2 = VisPt_2_Zinc3jetQunHratio_fit->GetBinContent(binNumber2);
                    if (fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()) >= 0. &&
                        fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()) <= 200.)
                        RatioValue2 = VisPt_2_Zinc3jetQunHratio_fit->GetBinContent(binNumber2);
                    else
                        RatioValue2 = 1.;
                    if (RatioValue2 > 2. || RatioValue2 < 0.5) RatioValue2 = 1.;
                }

                fill(hresponseThirdJetEta_Zinc3jet,
                     fabs(jets[2].v.Eta()),
                     fabs(genJets[2].v.Eta()),
                     weight * RatioValue);
                fill(hresponseThirdJetAbsRapidity_Zinc3jet,
                     fabs(jets[2].v.Rapidity()),
                     fabs(genJets[2].v.Rapidity()),
                     weight * RatioValue);
                fill(hresponseThirdJetEtaHigh_Zinc3jet,
                     fabs(jets[2].v.Eta()),
                     fabs(genJets[2].v.Eta()),
                     weight);
                fill(hresponseThirdJetRapidityHigh_Zinc3jet,
                     fabs(jets[2].v.Rapidity()),
                     fabs(genJets[2].v.Rapidity()),
                     weight);
                fill(hresponseJetsHT_Zinc3jet, jetsHT, genJetsHT, weight * RatioValue1);
                fill(hresponseVisPt_Zinc3jetQun,
                     fabs((jets[0].v + jets[1].v + jets[2].v + EWKBoson).Pt()),
                     fabs((genJets[0].v + genJets[1].v + genJets[2].v + genEWKBoson).Pt()),
                     weight * RatioValue2);

                /////Azimuthal cross check//////////////////////////////
                fill(hresponseDPhiZFirstJet_Zinc3jet,
                     fabs(EWKBoson.DeltaPhi(jets[0].v)),
                     fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                     weight);
                fill(hresponseDPhiZSecondJet_Zinc3jet,
                     fabs(EWKBoson.DeltaPhi(jets[1].v)),
                     fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                     weight);
                fill(hresponseDPhiZThirdJet_Zinc3jet,
                     fabs(EWKBoson.DeltaPhi(jets[2].v)),
                     fabs(genEWKBoson.DeltaPhi(genJets[2].v)),
                     weight);
                fill(hresponseDPhiFirstSecondJet_Zinc3jet,
                     fabs(jets[0].v.DeltaPhi(jets[1].v)),
                     fabs(genJets[0].v.DeltaPhi(genJets[1].v)),
                     weight);
                fill(hresponseDPhiFirstThirdJet_Zinc3jet,
                     fabs(jets[0].v.DeltaPhi(jets[2].v)),
                     fabs(genJets[0].v.DeltaPhi(genJets[2].v)),
                     weight);
                fill(hresponseDPhiSecondThirdJet_Zinc3jet,
                     fabs(jets[1].v.DeltaPhi(jets[2].v)),
                     fabs(genJets[1].v.DeltaPhi(genJets[2].v)),
                     weight);

                if (EWKBoson.Pt() > 150. && genEWKBoson.Pt() > 150.) {
                    fill(hresponseDPhiZFirstJet_ZPt150_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[0].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         weight);
                    fill(hresponseDPhiZSecondJet_ZPt150_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[1].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                         weight);
                    fill(hresponseDPhiZThirdJet_ZPt150_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[2].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[2].v)),
                         weight);
                    fill(hresponseDPhiFirstSecondJet_ZPt150_Zinc3jet,
                         fabs(jets[0].v.DeltaPhi(jets[1].v)),
                         fabs(genJets[0].v.DeltaPhi(genJets[1].v)),
                         weight);
                    fill(hresponseDPhiFirstThirdJet_ZPt150_Zinc3jet,
                         fabs(jets[0].v.DeltaPhi(jets[2].v)),
                         fabs(genJets[0].v.DeltaPhi(genJets[2].v)),
                         weight);
                    fill(hresponseDPhiSecondThirdJet_ZPt150_Zinc3jet,
                         fabs(jets[1].v.DeltaPhi(jets[2].v)),
                         fabs(genJets[1].v.DeltaPhi(genJets[2].v)),
                         weight);
                }

                if (EWKBoson.Pt() > 300. && genEWKBoson.Pt() > 300.) {
                    fill(hresponseDPhiZFirstJet_ZPt300_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[0].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         weight);
                    fill(hresponseDPhiZSecondJet_ZPt300_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[1].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                         weight);
                    fill(hresponseDPhiZThirdJet_ZPt300_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[2].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[2].v)),
                         weight);
                    fill(hresponseDPhiFirstSecondJet_ZPt300_Zinc3jet,
                         fabs(jets[0].v.DeltaPhi(jets[1].v)),
                         fabs(genJets[0].v.DeltaPhi(genJets[1].v)),
                         weight);
                    fill(hresponseDPhiFirstThirdJet_ZPt300_Zinc3jet,
                         fabs(jets[0].v.DeltaPhi(jets[2].v)),
                         fabs(genJets[0].v.DeltaPhi(genJets[2].v)),
                         weight);
                    fill(hresponseDPhiSecondThirdJet_ZPt300_Zinc3jet,
                         fabs(jets[1].v.DeltaPhi(jets[2].v)),
                         fabs(genJets[1].v.DeltaPhi(genJets[2].v)),
                         weight);
                }

                if (EWKBoson.Pt() > 150. && genEWKBoson.Pt() > 150. &&
                    (jets[0].v.Pt() + jets[1].v.Pt() + jets[2].v.Pt() > 300.) &&
                    (genJets[0].v.Pt() + genJets[1].v.Pt() + genJets[2].v.Pt() > 300.)) {
                    fill(hresponseDPhiZFirstJet_ZPt150_HT300_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[0].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[0].v)),
                         weight);
                    fill(hresponseDPhiZSecondJet_ZPt150_HT300_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[1].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[1].v)),
                         weight);
                    fill(hresponseDPhiZThirdJet_ZPt150_HT300_Zinc3jet,
                         fabs(EWKBoson.DeltaPhi(jets[2].v)),
                         fabs(genEWKBoson.DeltaPhi(genJets[2].v)),
                         weight);
                }
            }

            if (nGoodGenJets_20 >= 3 && passesgenLeptonCut && nGoodJets_20 >= 3 &&
                passesLeptonCut) {

                double RatioValue = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        ThirdJetPt_2_Zinc3jetHratio_fit->GetXaxis()->FindBin(jets_20[2].v.Pt());
                    RatioValue = ThirdJetPt_2_Zinc3jetHratio_fit->GetBinContent(binNumber);
                }

                fill(hresponseThirdJetPt_Zinc3jet,
                     jets_20[2].v.Pt(),
                     genJets_20[2].v.Pt(),
                     weight * RatioValue);
            }

            //-- Fourth Jet Pt
            if (nGoodGenJets >= 4 && passesgenLeptonCut && nGoodJets >= 4 && passesLeptonCut) {

                fill(hresponseFourthJetEta_Zinc4jet,
                     fabs(jets[3].v.Eta()),
                     fabs(genJets[3].v.Eta()),
                     weight);
                fill(hresponseFourthJetAbsRapidity_Zinc4jet,
                     fabs(jets[3].v.Rapidity()),
                     fabs(genJets[3].v.Rapidity()),
                     weight);
                fill(hresponseFourthJetEtaHigh_Zinc4jet,
                     fabs(jets[3].v.Eta()),
                     fabs(genJets[3].v.Eta()),
                     weight);
                fill(hresponseFourthJetRapidityHigh_Zinc4jet,
                     fabs(jets[3].v.Rapidity()),
                     fabs(genJets[3].v.Rapidity()),
                     weight);
                fill(hresponseJetsHT_Zinc4jet, jetsHT, genJetsHT, weight);
            }

            if (nGoodGenJets_20 >= 4 && passesgenLeptonCut && nGoodJets_20 >= 4 &&
                passesLeptonCut) {
                fill(
                    hresponseFourthJetPt_Zinc4jet, jets_20[3].v.Pt(), genJets_20[3].v.Pt(), weight);
            }

            //-- Fifth Jet Pt
            if (nGoodGenJets >= 5 && passesgenLeptonCut && nGoodJets >= 5 && passesLeptonCut) {

                fill(hresponseFifthJetEta_Zinc5jet,
                     fabs(jets[4].v.Eta()),
                     fabs(genJets[4].v.Eta()),
                     weight);
                fill(hresponseFifthJetAbsRapidity_Zinc5jet,
                     fabs(jets[4].v.Rapidity()),
                     fabs(genJets[4].v.Rapidity()),
                     weight);
                fill(hresponseFifthJetEtaHigh_Zinc5jet,
                     fabs(jets[4].v.Eta()),
                     fabs(genJets[4].v.Eta()),
                     weight);
                fill(hresponseFifthJetRapidityHigh_Zinc5jet,
                     fabs(jets[4].v.Rapidity()),
                     fabs(genJets[4].v.Rapidity()),
                     weight);
                fill(hresponseJetsHT_Zinc5jet, jetsHT, genJetsHT, weight);
            }

            if (nGoodGenJets_20 >= 5 && passesgenLeptonCut && nGoodJets_20 >= 5 &&
                passesLeptonCut) {
                fill(hresponseFifthJetPt_Zinc5jet, jets_20[4].v.Pt(), genJets_20[4].v.Pt(), weight);
            }

            if (nGoodGenJets >= 1 && passesgenLeptonCut && nGoodJets >= 1 && passesLeptonCut) {

                double RatioValue = 1.;
                double RatioValue1 = 1.;
                double RatioValue2 = 1.;
                if (UnfoldUnc) {
                    double binNumber =
                        JZB_2Hratio_fit->GetXaxis()->FindBin(hadronicR.Pt() - EWKBoson.Pt());
                    RatioValue = JZB_2Hratio_fit->GetBinContent(binNumber);
                    double binNumber1 =
                        JZB_ptLow_2Hratio_fit->GetXaxis()->FindBin(hadronicR.Pt() - EWKBoson.Pt());
                    // RatioValue1 =  JZB_ptLow_2Hratio_fit->GetBinContent(binNumber1);
                    if ((hadronicR.Pt() - EWKBoson.Pt()) > -50 &&
                        (hadronicR.Pt() - EWKBoson.Pt()) < 200.)
                        RatioValue1 = JZB_ptLow_2Hratio_fit->GetBinContent(binNumber1);
                    else
                        RatioValue1 = 1.;
                    if (RatioValue1 > 2. || RatioValue1 < 0.5) RatioValue1 = 1.;

                    double binNumber2 =
                        JZB_ptHigh_2Hratio_fit->GetXaxis()->FindBin(hadronicR.Pt() - EWKBoson.Pt());
                    RatioValue2 = JZB_ptHigh_2Hratio_fit->GetBinContent(binNumber2);
                }

                fill(hresponseHadRecoil, hadronicR.Pt(), genHadronicR.Pt(), weight);
                fill(hresponseJZB,
                     (hadronicR.Pt() - EWKBoson.Pt()),
                     (genHadronicR.Pt() - genEWKBoson.Pt()),
                     weight * RatioValue);
                if (EWKBoson.Pt() <= 50 && genEWKBoson.Pt() <= 50) {
                    fill(hresponseJZB_ptLow,
                         (hadronicR.Pt() - EWKBoson.Pt()),
                         (genHadronicR.Pt() - genEWKBoson.Pt()),
                         weight * RatioValue1);
                }
                if (EWKBoson.Pt() > 50 && genEWKBoson.Pt() > 50) {
                    fill(hresponseJZB_ptHigh,
                         (hadronicR.Pt() - EWKBoson.Pt()),
                         (genHadronicR.Pt() - genEWKBoson.Pt()),
                         weight * RatioValue2);
                }
            }
        }
        //=======================================================================================================//

        //	if(passesLeptonCut && !passesgenLeptonCut){
        //	    std::cout << "Event " << nEvents << " passes Reco cut but not gen cuts: reco -
        // gen\n"
        //		      << "Nleptons: " << nLeptons << " - " << ngenLeptons << "\n"
        //		      << "Pt l1:    "
        //		      << (nLeptons > 0 ? leptons[0].v.Pt() : -1) << " - "
        //		      << (ngenLeptons > 0 ? genLeptons[0].v.Pt() : -1) << "\n"
        //		      << "Pt l2:    "
        //		      << (nLeptons > 1 ? leptons[1].v.Pt() : -1) << " - "
        //		      << (ngenLeptons > 1 ? genLeptons[1].v.Pt() : -1) << "\n"
        //		      << "Eta l1:    "
        //		      << (nLeptons > 0 ? leptons[0].v.Eta() : -1) << " - "
        //		      << (ngenLeptons > 0 ? genLeptons[0].v.Eta() : -1) << "\n"
        //		      << "Eta l2:    "
        //		      << (nLeptons > 1 ? leptons[1].v.Eta() : -1) << " - "
        //		      << (ngenLeptons > 1 ? genLeptons[1].v.Eta() : -1) << "\n";
        //
        //	    if(nLeptons >= 2 && ngenLeptons >= 2){
        //		std::cout << "Charges: " << (leptons[0].charge ? '+' : '-' )
        //                          << (leptons[1].charge ? '+' : '-' )
        //                          << "-" << (genLeptons[0].charge ? '+' : '-')
        //			  << (genLeptons[1].charge ? '+' : '-');
        //	    }
        //	}
    } // End of loop over all the events
    cout << endl;
    //==========================================================================================================//

    JobInfo->SetBinContent(kNEvts, nEventsToProcess);
    JobInfo->SetBinContent(kNEvtsSample, nentries);
    JobInfo->SetBinContent(kNEvtsAllJobs, nEventsToProcessTot);
    double data_frac = EvtIsRealData ? (nEventsToProcessTot / double(nentries)) : yieldScale;
    JobInfo->SetBinContent(kLumi, lumi_ * data_frac);
    // store integrated luminosity in Lumi histogram:
    Lumi->SetBinContent(1., lumi_ * data_frac);
    JobInfo->SetBinContent(kXsec, xsec_);

    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    //         Writing file              //
    //==================================//

    outputFile->cd();

    //--- Save all the histograms ---
    unsigned short numbOfHistograms = listOfHistograms.size();

    // factor the histogram contents must be scaled down to
    // matches the data integrated luminosity times, in case of
    // of multi MC jobs, the fraction of MC events processed by this job
    if (!EvtIsRealData) {
        double xsec = xsec_;
        if (xsec == 0) {
            std::cerr << "Warning: cross section value for MC sample " << sampleLabel_
                      << " is null or was not specified. We will assume the event"
                      << " weights are normalizes such that the cross section on pb "
                      << " is equal to the sum of weights divivided by the numnber of events\n";
            if (InEvtCount_ > 0) {
                norm_ = data_frac * lumi_ * xsecFactor_ / InEvtCount_ *
                        (nentries / nEventsToProcessTot);
                std::cout << "Used norm_: data_frac * lumi_ * xsecFactor_  / InEvtCount_ "
                          << "* (nentries / nEventsToProcessTot)" << data_frac << " * " << lumi_
                          << " * " << xsecFactor_ << " / " << InEvtCount_ << " * (" << nentries
                          << " / " << nEventsToProcessTot << ")"
                          << " = " << norm_ << "\n";

            } else {
                norm_ = data_frac * lumi_ * xsecFactor_ / nEventsToProcessTot;
                std::cout << "Used norm_: data_frac * lumi_ * xsecFactor_  / nEventsToProcessTot = "
                          << data_frac << " * " << lumi_ << " * " << xsecFactor_ << " / "
                          << nEventsToProcessTot << " = " << norm_ << "\n";
            }
        } else {
            // sum of weights before any cut over the full dataset:
            if (InEvtWeightSums_.size() > 0) {
                // normalisation is defined to get perfect normalisation when running
                // on the full dataset statistics by just adding up the histograms.
                // In case of partial dataset processing, direct sum will include an
                // approximation (*) which can be removed by using the JobWeight information
                // stored in the JobInfo histograms.
                //
                //(*) sum of processed event weights equals to the sum over all the events times
                // the fraction of processed events.
                //
                norm_ = data_frac * lumi_ * xsec * xsecFactor_ / InEvtWeightSums_[0] * nentries /
                        nEventsToProcessTot;
                std::cout
                    << "Used norm_: data_frac * lumi_ * xsec * xsecFactor_  / nEvtWeightSums_[0] "
                       "* nentries / nEventsToProcessTot = "
                    << data_frac << " * " << lumi_ << " * " << xsec << " * " << xsecFactor_ << " / "
                    << InEvtWeightSums_[0] << " * " << nentries << " / " << nEventsToProcessTot
                    << "=" << norm_ << "\n";
            } else {
                norm_ = data_frac * lumi_ * xsec * xsecFactor_ / processedEventMcWeightSum_ *
                        nEventsToProcess / nEventsToProcessTot;
                std::cout << "Used norm_: data_frac * lumi_ * xsec * xsecFactor_  / "
                             "processedEventMcWeightSum_"
                             "* nEventsToProcess / nEventsToProcessTot = "
                          << data_frac << " * " << lumi_ << " * " << xsec << " * " << xsecFactor_
                          << " / " << processedEventMcWeightSum_ << " * " << nEventsToProcess
                          << " / " << nEventsToProcessTot << "=" << norm_ << "\n";
            }
        }

        if (norm_ == 0) {
            std::cerr << "Error: normaliation factor for sample " << fileName
                      << " is null! Aborts at " __FILE__ ":" << __LINE__ << "." << std::endl;
            abort();
        }
    } else {
        norm_ = 1;
    }

    if (!EvtIsRealData) {
        JobInfo->SetBinContent(kXsec, xsec_ * xsecFactor_);
        double a = 1.;
        if (EvtWeightSums_.size()) a = EvtWeightSums_[0];
        JobInfo->SetBinContent(kJobWeight, processedEventMcWeightSum_ / a);
    }

    for (unsigned short i(0); i < numbOfHistograms; i++) {
        string hName = listOfHistograms[i]->GetName();
        if ((!hasGenInfo && hName.find("gen") != string::npos) ||
            (!hasRecoInfo && hName.find("gen") == string::npos))
            continue;
        // finalize normalisation of MC histograms:
        if (!EvtIsRealData && !listOfHistograms[i]->TestBit(TH1::kIsAverage)) {
            listOfHistograms[i]->Scale(norm_);
        }

        listOfHistograms[i]->Write();
    }

    //--- let's delete all histograms ---
    for (unsigned short i(0); i < numbOfHistograms; i++) {
        delete listOfHistograms[i];
    }

    outputFile->Write();
    outputFile->Close();

    //==========================================================================================================//

    nEffEventsPassingTrigger *= nEvents / genWeightSum;
    nEffEventsVInc0JetsNoTrig *= nEvents / genWeightSum;
    nEffEventsVInc0Jets *= nEvents / genWeightSum;
    nEffEventsVInc1Jets *= nEvents / genWeightSum;
    nEffEventsVInc2Jets *= nEvents / genWeightSum;
    nEffEventsVInc3Jets *= nEvents / genWeightSum;
    nEffEventsWithTwoGoodLeptons *= nEvents / genWeightSum;
    nEffEventsWithTwoGoodLeptonsWithOppCharge *= nEvents / genWeightSum;
    nEffEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass *= nEvents / genWeightSum;

    nEffGenEventsVInc0Jets *= nEvents / genWeightSum;
    nEffGenEventsVInc1Jets *= nEvents / genWeightSum;
    nEffGenEventsVInc2Jets *= nEvents / genWeightSum;
    nEffGenEventsVInc3Jets *= nEvents / genWeightSum;
    nEffGenEventsWithTwoGoodLeptons *= nEvents / genWeightSum;
    nEffGenEventsWithTwoGoodLeptonsWithOppCharge *= nEvents / genWeightSum;
    nEffGenEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass *= nEvents / genWeightSum;

    cout << "Number of processed events                                : " << nEvents << endl;
    if (maxFiles_ < 0) {
        cout << "Fraction of processed events from dataset                 : " << nEvents << " / "
             << EvtCount_ << " = " << (nEvents / double(EvtCount_)) << endl;
        if (EvtIsRealData) {
            cout << "\tvalue stored in file .mcYieldScale for the '--mcYieldScale -1' auto "
                    "normalisation option.\n";
            std::ofstream f(outputDirectory + "/.mcYieldScale#");
            f << data_frac << "\n";
            f.close();
            rename(outputDirectory + "/.mcYieldScale#", outputDirectory + "/.mcYieldScale");
        }
    }
    cout << "Number of events passing the trigger                      : " << nEventsPassingTrigger
         << "\n";
    cout << "Number with two good leptons (gen)                        : "
         << nEventsWithTwoGoodLeptons << " (" << nGenEventsWithTwoGoodLeptons << ")" << endl;
    cout << "Number with two good leptons of opp. charge (gen)         : "
         << nEventsWithTwoGoodLeptonsWithOppCharge << " ("
         << nGenEventsWithTwoGoodLeptonsWithOppCharge << ")" << endl;
    cout << "Number with two good leptons of opp. charge and good mass : "
         << nEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass << " ("
         << nGenEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass << ")" << endl;
    cout << "Number Reco Inclusif V + 0 jets                           : " << nEventsVInc0Jets
         << endl;
    cout << "Number Reco Inclusif V + 0 jets no trigger                : " << nEventsVInc0JetsNoTrig
         << endl;
    cout << "Number Reco Inclusif V + 1 jets                           : " << nEventsVInc1Jets
         << endl;
    cout << "Number Reco Inclusif V + 2 jets                           : " << nEventsVInc2Jets
         << endl;
    cout << "Number Reco Inclusif V + 3 jets                           : " << nEventsVInc3Jets
         << endl;
    cout << "Number GEN Inclusif V + 0 jets                            : " << nGenEventsVInc0Jets
         << endl;
    cout << "Number GEN Inclusif V + 1 jets                            : " << nGenEventsVInc1Jets
         << endl;
    cout << "Number GEN Inclusif V + 2 jets                            : " << nGenEventsVInc2Jets
         << endl;
    cout << "Number GEN Inclusif V + 3 jets                            : " << nGenEventsVInc3Jets
         << endl;

    cout << "Eff. number with two good leptons (gen)                        : "
         << nEffEventsWithTwoGoodLeptons << " (" << nEffGenEventsWithTwoGoodLeptons << ")" << endl;
    cout << "Eff. number with two good leptons of opp. charge (gem)         : "
         << nEffEventsWithTwoGoodLeptonsWithOppCharge << " ("
         << nEffGenEventsWithTwoGoodLeptonsWithOppCharge << ")" << endl;
    cout << "Eff. number with two good leptons of opp. charge and good mass : "
         << nEffEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass << " ("
         << nEffGenEventsWithTwoGoodLeptonsWithOppChargeAndGoodMass << ")" << endl;
    cout << "Eff. number Reco Inclusif V + 0 jets                           : "
         << nEffEventsVInc0Jets << endl;
    cout << "Eff. number Reco Inclusif V + 0 jets no trigger                : "
         << nEffEventsVInc0JetsNoTrig << endl;
    cout << "Eff. number Reco Inclusif V + 1 jets                           : "
         << nEffEventsVInc1Jets << endl;
    cout << "Eff. number Reco Inclusif V + 2 jets                           : "
         << nEffEventsVInc2Jets << endl;
    cout << "Eff. number Reco Inclusif V + 3 jets                           : "
         << nEffEventsVInc3Jets << endl;
    cout << "Eff. number GEN Inclusif V + 0 jets                            : "
         << nEffGenEventsVInc0Jets << endl;
    cout << "Eff. number GEN Inclusif V + 1 jets                            : "
         << nEffGenEventsVInc1Jets << endl;
    cout << "Eff. number GEN Inclusif V + 2 jets                            : "
         << nEffGenEventsVInc2Jets << endl;
    cout << "Eff. number GEN Inclusif V + 3 jets                            : "
         << nEffGenEventsVInc3Jets << endl;
    cout << "Sum of MC event weights                                        : "
         << processedEventMcWeightSum_ << endl;
    if (!EvtIsRealData) {
        if (xsec_ > 0) {
            cout << "MC norm., yield_scale*lumi*xsec*skim_accep/sum_weights*unc_var. : "
                 << yieldScale << " * " << lumi_ << " * " << xsec_ << " * " << skimAccep_[0]
                 << " / " << processedEventMcWeightSum_ << " * " << xsecFactor_ << " = "
                 << norm_ / processedEventMcWeightSum_ << endl;
        } else {
            cout << "MC norm., yield_scale*lumi*skim_accep/n_events*unc_var. : " << yieldScale
                 << " * " << lumi_ << " * " << skimAccep_[0] << " / " << nEvents << " * "
                 << xsecFactor_ << " = " << norm_ / nEvents << endl;
        }
    }
}

void ZJets::initLHAPDF(TString pdfSet, int pdfMember)
{
    LHAPDF::initPDFSet(1, pdfSet.Data(), pdfMember);
    // LHAPDF::initPDFSet(1, "CT10.LHgrid");
    LHAPDF::initPDFSet(2, "cteq6ll.LHpdf");
    if (pdfMember > LHAPDF::numberPDF(1) + 1) {
        cout << "Warning pdfMember to high" << endl;
        return;
    }
}

double ZJets::computePDFWeight()
{

    //-- get the pdgId of the two colliding partons
    double wPdf(1.);
    /*

        int id1 = pdfInfo_->at(0);
        int id2 = pdfInfo_->at(1);
        if (id1 == 21) id1 = 0; // 21 is Pythia convention for gluon, but needs to be 0 for LHAPDF
        if (id2 == 21) id2 = 0;

        double pdf1  = LHAPDF::xfx(1, pdfInfo_->at(2), pdfInfo_->at(4), id1);
        double pdf2  = LHAPDF::xfx(1, pdfInfo_->at(3), pdfInfo_->at(4), id2);
        double pdf01 = LHAPDF::xfx(2, pdfInfo_->at(2), pdfInfo_->at(4), id1);
        double pdf02 = LHAPDF::xfx(2, pdfInfo_->at(3), pdfInfo_->at(4), id2);

        if (pdfInfo_->at(2) * pdfInfo_->at(3) > 0) {
            wPdf = pdf1 * pdf2;
            if (pdf01*pdf02 <= 0 || pdf1*pdf2 <= 0) {
                cout << "Small problem" << endl;
                wPdf = 1;
            }
            else {
                wPdf /= (pdf01 * pdf02);
            }
        }
    */

    return wPdf;
}

void ZJets::getMuons(vector<leptonStruct> &leptons, vector<leptonStruct> &vetoMuons)
{

    //--- get the number of Muon candidates from the vector size ---
    unsigned short nTotLeptons(MuEta->size());

    //    bool eventTrigger = false;
    // we also have event trigger variables --> we should at least match one of the leptons to
    // trigger
    /* // CommentAG: check patMuonTrig
        for (unsigned short i(0); i < nTotLeptons; i++) {
            int whichTrigger(patMuonTrig_->at(i));
            if (lepSel == "SMu" && (whichTrigger & 0x1)) eventTrigger = true;
        }
      */

    for (unsigned short i(0); i < nTotLeptons; i++) {
        //        double muonId = 0;

        // CommentAG: don't have patMuonCombId
        /*
                if(fileName.Index("mcatnlo") >= 0 || fileName.Index("MG-MLM") >= 0) {
                    muonId = (double) patMuonCombId_Int->at(i);
                }
                else {
                    muonId = (double) patMuonCombId_Double->at(i);
                }
        */
        leptonStruct mu(MuPt->at(i),
                        MuEta->at(i),
                        MuPhi->at(i),
                        MuE->at(i),
                        MuCh->at(i),
                        MuIdTight->at(i), //  CommentAG: muonId; Tight muons are selected in Bonzai
                        MuPfIso->at(i),
                        MuEta->at(i),
                        MuTkLayerCnt->at(i)); // CommentAG

        //  cout << mu.MuTkLayer << "\n";

        // float qter = 1.0;
        /*
                if (doRochester) {
                    if (!EvtIsRealData) {
                         //rmcor->momcor_mc(mu.v, (float)mu.charge, 0, qter);
                         //cout << MuTkLayerCnt->at(i) << "\n";
                         //TRandom *rochester;
                         //float aRand = rochester->Uniform(1.0);
                         //float aRand2 = rochester->Uniform(1.0);
                          double mcSF = rc->kScaleAndSmearMC((float)mu.charge, mu.v.Pt(),
           mu.v.Eta(), mu.v.Phi(), 6, gRandom->Rndm(), gRandom->Rndm());
                         // cout << mu.v.Pt() << " , " << mu.v.Eta() << " , " <<  mu.v.Phi() << " ,
           " << mu.v.E() << "\n";
                          mu.v.SetPtEtaPhiE(mu.v.Pt()*mcSF, mu.v.Eta(), mu.v.Phi(), mu.v.E());

                         // cout << mcSF << "\n";
                    }
                    else {
                         //rmcor->momcor_data(mu.v, (float)mu.charge, 0, qter);
                          double dataSF = rc->kScaleDT((float)mu.charge, mu.v.Pt(), mu.v.Eta(),
           mu.v.Phi());
                          mu.v.SetPtEtaPhiE(mu.v.Pt()*dataSF, mu.v.Eta(), mu.v.Phi(),mu.v.E());
                          //cout << dataSF << "\n";
                    }
                }
        */
        bool muPassesPtCut(mu.v.Pt() >= (lepPtCutMin * 0.8));
        bool muPassesEtaCut(fabs(mu.v.Eta()) <= 0.1 * lepEtaCutMax);
        // bool muPassesIdCut(mu.id & 0x1);  //CommentAG: Tight muons Id are selected in the Bonzai
        // Maker
        bool muPassesIdCut(mu.id & 1);

        bool muPassesIsoCut(0);
        ///        if (lepSel == "DMu" && mu.iso < 0.25) muPassesIsoCut = 1;
        //        else if (lepSel == "SMu" && mu.iso < 0.15) muPassesIsoCut = 1;
        muPassesIsoCut = (mu.iso < muIso_);
        bool muPassesTrig(1); // no matching with leptons offtrigger
        // if (lepSel == "DMu" && (mu.trigger & 0x4)) muPassesTrig = 1;       // HLT_Mu17_Mu8 !!!!
        // changed from 0x8 to 0x4
        //	const Long64_t dimutrigmask = 19 | 21;
        //	const Long64_t mutrigmask = 18 | 22;
        //        if (lepSel == "DMu" && (TrigHltDiMu  & dimutrigmask)) muPassesTrig = 1;
        //	else if (lepSel == "SMu" && (mu.trigger & mutrigmask)) muPassesTrig = 1;  //
        // HLT_IsoMu24_eta2p1_v

        //--- veto muons ---
        bool muPassesVetoPtCut(mu.v.Pt() >= 15);
        bool muPassesVetoEtaCut(fabs(mu.v.Eta()) <= 2.4);
        // bool muPassesVetoIdCut(mu.id > 0); // muon Id  // CommentAG: not sure what does it mean

        /// for files obtained form bugra
        if (fileName.Index("Sherpa_Bugra_1_13_UNFOLDING") >= 0 && mu.trigger > 0)
            muPassesTrig = 1; // Bugra only keeps the double electron trigger !!!!!

        // select the good muons only
        // cout << i << " pt=" << mu.v.Pt() << " eta=" << mu.v.Eta() << " phi=" << mu.v.Phi() <<
        // endl;
        // cout << "id=" << muPassesIdCut << " iso=" << muPassesIsoCut << " trig=" << muPassesTrig
        // << endl;

        // if (muPassesPtCut && muPassesEtaCut && muPassesIdCut && muPassesIsoCut &&
        // (!useTriggerCorrection || muPassesTrig || eventTrigger)) {   // CommentAG: this is
        // original line which is replaced by:

        if (muPassesPtCut && muPassesEtaCut && muPassesIdCut && muPassesIsoCut && muPassesTrig) {
            leptons.push_back(mu);
        }
        // select the veto muons
        else if (lepSel == "SMu" && muPassesVetoPtCut &&
                 muPassesVetoEtaCut) { // CommentAG:  need to check muPassesVetoIdCut!
            vetoMuons.push_back(mu);
        }

    } // End of loop over all the muons
}

void ZJets::getElectrons(vector<leptonStruct> &leptons, vector<leptonStruct> &vetoElectrons)
{
    //--- get the number of Electron candidates from the vector size ---
    unsigned short nTotLeptons(ElEta->size());

    // if we don't really care to match both leptons to trigger
    // bool eventTrigger = false;
    for (unsigned short i(0); i < nTotLeptons; i++) {
        /* // CommentAG
                int whichTrigger(patElecTrig_->at(i));
                if (lepSel == "DE" && (whichTrigger & 0x2)) eventTrigger = true;      //
           HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v
                else if (lepSel == "SE" && (whichTrigger & 0x1)) eventTrigger = true; //
           HLT_Ele27_WP80_v
        */
    }
    for (unsigned short i(0); i < nTotLeptons; i++) {

        leptonStruct ele(ElPt->at(i),
                         ElEta->at(i),
                         ElPhi->at(i),
                         ElE->at(i),
                         ElCh->at(i),
                         ElId->at(i),
                         ElPfIsoRho->at(i),
                         ElEtaSc->at(i),
                         0.); // CommentAG patElecTrig_->at(i)

        //--- good electrons ---
        bool elePassesPtCut(ele.v.Pt() >= (lepPtCutMin * 0.8));
        bool elePassesEtaCut(fabs(ele.scEta) <= min(1.4442, 0.1 * lepEtaCutMax) ||
                             (fabs(ele.scEta) >= 1.566 && fabs(ele.scEta) <= 0.1 * lepEtaCutMax));
        bool elePassesIdCut(ele.id >= 4); /// 4 is medium ID, 2 is Loose ID
        bool elePassesIsoCut(ele.iso < eIso_);
        // bool elePassesAnyTrig(ele.trigger & 0x2);
        bool elePassesAnyTrig(true); // no matching with lepton from trigger.
        if (fileName.Index("Sherpa_Bugra_1_13_UNFOLDING") > 0) elePassesAnyTrig = true;

        //--- veto electrons ---
        bool elePassesVetoPtCut(ele.v.Pt() >= 15);
        bool elePassesVetoEtaCut(fabs(ele.v.Eta()) <= 2.4);
        bool elePassesVetoIdCut(ele.id >= 2); // ele Id
        bool elePassesVetoIsoCut(ele.iso < 0.25);

        // select the good electrons only
        // if (elePassesPtCut && elePassesEtaCut && elePassesIdCut && elePassesIsoCut &&
        // (!useTriggerCorrection || elePassesAnyTrig || eventTrigger))
        if (elePassesPtCut && elePassesEtaCut && elePassesIdCut && elePassesIsoCut &&
            elePassesAnyTrig) {
            leptons.push_back(ele);
        }
        // select the veto electrons
        else if (elePassesVetoPtCut && elePassesVetoEtaCut && elePassesVetoIdCut &&
                 elePassesVetoIsoCut) {
            vetoElectrons.push_back(ele);
        }

    } // End of loop over all the electrons
}

ZJets::ZJets(const TString &lepSel_,
             TString sampleLabel,
             TString fileName_,
             float lumi,
             bool useTriggerCorrection_,
             int systematics_,
             int direction_,
             float xsecUnc_,
             int lepPtCutMin_,
             int lepEtaCutMax_,
             int jetPtCutMin_,
             int jetEtaCutMax_,
             Long_t maxEvents_,
             TString outDir_,
             TString bonzaiDir,
             int maxFiles)
    : HistoSetZJets(lepSel_),
      outputDirectory(outDir_),
      fileName(fileName_),
      lumi_(lumi),
      useTriggerCorrection(useTriggerCorrection_),
      systematics(systematics_),
      direction(direction_),
      xsecUnc(xsecUnc_),
      lepPtCutMin(lepPtCutMin_),
      lepEtaCutMax(lepEtaCutMax_),
      jetPtCutMin(jetPtCutMin_),
      jetEtaCutMax(jetEtaCutMax_),
      nMaxEvents(maxEvents_),
      lepSel(lepSel_),
      xsec_(0.),
      sampleLabel_(sampleLabel),
      maxFiles_(maxFiles),
      triggerMask_EraBG(0),
      triggerMask_EraH(0),
      triggerMaskSet_(false),
      muIso_(0),
      eIso_(0)
{
    //--- Create output directory if necessary ---
    TString command = "mkdir -p " + outputDirectory;
    system(command);

    //--------------------------------------------

    rejectBTagEvents = lepSel.BeginsWith("S");

    fChain = new TChain("", "");

    TString fullFileName;
    TString baseName;

    canonizeInputFilePath(bonzaiDir, fileName, &fullFileName, &baseName);

    fileName = baseName;

    Input->SetTitle(fullFileName);

    readCatalog(fullFileName,
                bonzaiDir,
                maxFiles,
                &lumi_,
                &xsec_,
                fChain,
                &fBonzaiHeaderChain,
                &fBitFieldsChain);

    getMcNorm();

    if (!setTriggerMask()) {
        std::cerr << "Failed to set trigger mask. Aborts at " __FILE__ ":" << __LINE__ << "."
                  << std::endl;
        abort();
    }
    // NEW DA
    // std::cout << "Trigger mask: " << std::hex << "0x" << triggerMask_ << std::dec << "\n";
    std::cout << "Trigger mask EraBF: " << std::hex << "0x" << triggerMask_EraBG << std::dec
              << "\n";
    std::cout << "Trigger mask EraGH: " << std::hex << "0x" << triggerMask_EraH << std::dec << "\n";
}

void ZJets::canonizeInputFilePath(const TString &bonzaiDir,
                                  const TString &fileName,
                                  TString *fullFileName,
                                  TString *baseName,
                                  TString *ext)
{

    if (fileName.BeginsWith("/")) { // absolute path
        *fullFileName = fileName;
    } else {
        *fullFileName = bonzaiDir + "/" + fileName;
    }

    if (fullFileName->BeginsWith("/store/")) {
        fullFileName->Insert(0, "root://eoscms.cern.ch//eos/cms");
    }

    // fileName is expected to contain only the basename without extension
    // remove the .root, .txt extensions:
    if (baseName) {
        *baseName = gSystem->BaseName(fileName);
        if (baseName->EndsWith(".root")) {
            baseName->Remove(baseName->Length() - 5, 5);
            if (ext) *ext = ".root";
        }
        if (baseName->EndsWith(".txt")) {
            baseName->Remove(baseName->Length() - 4, 4);
            if (ext) *ext = ".txt";
        }
    }
}

void ZJets::readCatalog(const TString &fullFileName,
                        const TString &bonzaiDir,
                        int maxFiles,
                        double *pLumi,
                        double *pXsec,
                        TChain *pEventTreeChain,
                        TChain *pBonzaiHeaderChain,
                        TChain *pBitFieldsChain)
{
    regex_t xsecLine;
    int rc = regcomp(
        &xsecLine, "[#*][[:space:]]*sample xsec[[:space:]:=]\\+\\([[:digit:].eE+-]\\+\\)", 0);
    if (rc) {
        char buffer[256];
        regerror(rc, &xsecLine, buffer, sizeof(buffer));
        buffer[sizeof(buffer) - 1] = 0;
        std::cerr << "Bug found in " << __FILE__ << ":" << __LINE__ << ": " << buffer << "\n";
    }

    regex_t lumiLine;
    rc = regcomp(&lumiLine, "[#*][[:space:]]*lumi[[:space:]:=]\\+\\([[:digit:].eE+-]\\+\\)", 0);
    if (rc) {
        char buffer[256];
        regerror(rc, &lumiLine, buffer, sizeof(buffer));
        buffer[sizeof(buffer) - 1] = 0;
        std::cerr << "Bug found in " << __FILE__ << ":" << __LINE__ << ": " << buffer << "\n";
    }

    if (isRootFile(fullFileName)) {
        TString treePath = fullFileName + "/tupel/EventTree";
        TString bonzaiHeaderPath = fullFileName + "/tupel/BonzaiHeader";
        TString bonzaiBitFieldsPath = fullFileName + "/tupel/BitFields";
        cout << "Loading file: " << fullFileName << endl;
        if (pEventTreeChain) pEventTreeChain->Add(treePath);
        if (pBonzaiHeaderChain) {
            // check presence of the BonzaiHeader tree:
            TFile *f = TFile::Open(fullFileName);
            if (f && !f->IsZombie()) {
                if (f->GetDirectory("tupel")->FindKey("BonzaiHeader")) {
                    pBonzaiHeaderChain->Add(bonzaiHeaderPath);
                } else {
                    std::cerr
                        << "Warning: the tree BonzaiHeader was not found in file " << fullFileName
                        << ". We will assume we run on a boabab file and not Baobab->Bonzai "
                        << "acceptance correction will be considered. This message can be ignored "
                        << "if for this sample Boabab ntuples are usd as input.\n";
                }
            }
        }
        if (pBitFieldsChain) pBitFieldsChain->Add(bonzaiBitFieldsPath);
    } else {
        int (*closeFunc)(FILE *);
        FILE *f = eosOpen(fullFileName, &closeFunc);
        if (!f) {
            std::cerr << "Failed to  open file " << fullFileName << ".\n";
        } else {
            std::cout << "Reading input files from catalog file " << fullFileName << "\n";
            string line;
            char *buffer = 0;
            size_t buffer_size = 0;
            int ifile = 0;
            // if maxFiles = 0 only catalog header is read.
            enum { False = 0, True, Unknown } isBonzai = Unknown;
            while (!feof(f)) {
                ssize_t len = getline(&buffer, &buffer_size, f);
                if (len < 0) break;
                char *line = buffer;
                size_t n = len - 1;
                // trim white spaces:
                while (line[0] == ' ' || line[0] == '\t') {
                    ++line;
                    --n;
                }
                while (n >= 0 &&
                       (line[n] == ' ' || line[n] == '\t' || line[n] == '\r' || line[n] == '\n')) {
                    line[n] = 0;
                    --n;
                }

                regmatch_t pmatch[2];
                if (pXsec &&
                    !regexec(&xsecLine, line, sizeof(pmatch) / sizeof(pmatch[0]), pmatch, 0)) {
                    line[pmatch[1].rm_eo] = 0;
                    *pXsec = strtod(line + pmatch[1].rm_so, 0);
                    if (*pXsec == 0) {
                        std::cerr << "Value of parameter 'sample xsec', " << line + pmatch[1].rm_so
                                  << " found in file " << fullFileName << " is not valid.\n";
                    }
                } else if (pLumi &&
                           !regexec(
                               &lumiLine, line, sizeof(pmatch) / sizeof(pmatch[0]), pmatch, 0)) {
                    line[pmatch[1].rm_eo] = 0;
                    *pLumi = strtod(line + pmatch[1].rm_so, 0);
                    if (*pLumi == 0) {
                        std::cerr << "Integrated luminosity parameter value, "
                                  << line + pmatch[1].rm_so << " found in file " << fullFileName
                                  << " is not valid.\n";
                    }
                }

                // skip empty lines,  comment lines and metadata lines:
                if (line[0] == 0 || line[0] == '#' || line[0] == '*') continue;

                if (maxFiles == 0 || (pEventTreeChain == 0 && pEventTreeChain == 0)) break;

                // keep content of first column only:
                char *p = line;
                while (*p != 0 && *p != ' ' && *p != '\t' && *p != '\r' && *p != '\n') ++p;
                *p = 0;

                // following check is done having read the header
                // such that maxFiles = 0 can be used to read only
                // the header.
                if (maxFiles >= 0 && ifile >= maxFiles) break;

                TString filePath = TString(line);

                if (filePath[0] != '/') {
                    filePath.Insert(0, TString(bonzaiDir) + "/");
                }
                if (filePath.BeginsWith("/store/")) {
                    filePath.Insert(0, "root://eoscms.cern.ch//eos/cms");
                }

                TString treePath = filePath + "/tupel";
                TString bonzaiHeaderPath = treePath + "/BonzaiHeader";
                TString bonzaiBitFieldsPath = treePath + "/BitFields";
                treePath += "/EventTree";
                // std::cout << "Adding path " << treePath << " to the tree chain.\n";
                if (pEventTreeChain) pEventTreeChain->Add(treePath);

                if (pBonzaiHeaderChain && (isBonzai == Unknown)) {
                    // check presence of the BonzaiHeader tree. It is checked
                    // only on the first file which can be succesfully opened
                    // assuming that all files of the catalog are the same.
                    TFile *f = TFile::Open(filePath);
                    if (f && !f->IsZombie()) {
                        if (f->GetDirectory("tupel")->FindKey("BonzaiHeader")) {
                            isBonzai = True;
                        } else {
                            isBonzai = False;
                            std::cerr << "Warning: the tree BonzaiHeader was not found in file "
                                      << fullFileName << ". We will assume we run on a boabab file "
                                                         "and not Baobab->Bonzai "
                                      << "acceptance correction will be considered. This message "
                                         "can be ignored "
                                      << "if for this sample Boabab ntuples are usd as input.\n";
                        }
                    }
                }
                if (isBonzai) pBonzaiHeaderChain->Add(bonzaiHeaderPath);
                if (pBitFieldsChain) pBitFieldsChain->Add(bonzaiBitFieldsPath);
                ++ifile;
            } // next line
            std::cout << "Closing catalog file " << fullFileName << "\n";
            if (buffer) free(buffer);
            closeFunc(f);
        } // file opening succeeded
    }     // is root file
    regfree(&xsecLine);
    regfree(&lumiLine);
}

//#define weight_bug //to read bonzai version 2.

void ZJets::getMcNorm()
{
    Int_t InEvtCount = 0;
    //#ifdef weight_bug
    //    std::vector<Double_t> InEvtWeightSums(1,0);
    //    std::vector<Double_t> EvtWeightSums(1,0);
    //    fBonzaiHeaderChain.SetBranchAddress("InEvtWeightSums", &InEvtWeightSums[0]);
    //    fBonzaiHeaderChain.SetBranchAddress("EvtWeightSums", &EvtWeightSums[0]);
    //#else
    std::vector<Double_t> *InEvtWeightSums = 0;
    std::vector<Double_t> *EvtWeightSums = 0;
    EvtCount_ = fChain->GetEntries();

    if (fBonzaiHeaderChain.GetListOfFiles()->IsEmpty()) {
        std::cerr << "Running on a boabab file, skim acceptance = 1\n";
        skimAccep_ = std::vector<double>(1, 1.);
        InEvtWeightSums_ = std::vector<Double_t>(InEvtWeightSums->size(), 0);
        EvtWeightSums_ = std::vector<Double_t>(EvtWeightSums->size(), 0);
    } else {
        fBonzaiHeaderChain.SetBranchAddress("InEvtWeightSums", &InEvtWeightSums);
        fBonzaiHeaderChain.SetBranchAddress("EvtWeightSums", &EvtWeightSums);
        //#endif
        // for(Long64_t i = 0; i < nfiles; ++ i){
        int nheaders =
            fBonzaiHeaderChain.GetEntries(); // can be several in case files were merged with haddd
        for (int ientry = 0; ientry < nheaders; ++ientry) {
            fBonzaiHeaderChain.GetEntry(ientry);
            if (ientry == 0) {
                InEvtWeightSums_ = std::vector<Double_t>(InEvtWeightSums->size(), 0);
                EvtWeightSums_ = std::vector<Double_t>(EvtWeightSums->size(), 0);
                if (InEvtWeightSums->size() != EvtWeightSums->size()) {
                    std::cerr << "InEvtWeightSums and EvtWeightSums branches "
                                 "of input BonzaiHeader tree have different size ("
                                 "resp. "
                              << InEvtWeightSums->size() << " and " << EvtWeightSums->size()
                              << "). Aborts at " __FILE__ ":" << __LINE__ << ".\n";
                    abort();
                }
            }
            if (InEvtWeightSums->size() != InEvtWeightSums_.size()) {
                std::cerr << "Inconsistency in number of elements of "
                          << " InEvtWeightSums branch of input files! Aborts at " __FILE__ ":"
                          << __LINE__ << ".\n";
                abort();
            }
            if (EvtWeightSums->size() != EvtWeightSums_.size()) {
                std::cerr << "Inconsistency in number of elements of EvtWeightSums "
                             "branch of input files! Aborts at " __FILE__ ":"
                          << __LINE__ << ".\n";
                abort();
            }
            for (size_t i = 0; i < InEvtWeightSums_.size(); ++i) {
                InEvtWeightSums_[i] += (*InEvtWeightSums)[i];
            }
            for (size_t i = 0; i < EvtWeightSums_.size(); ++i) {
                EvtWeightSums_[i] += (*EvtWeightSums)[i];
            }
            InEvtCount_ += InEvtCount;
        }
        //}

        if (InEvtWeightSums_.size() > 0) {
            std::cerr << "InEvtWeightSums_[0] = " << InEvtWeightSums_[0] << "\n";
        }

        if (EvtWeightSums_.size() == 0 || InEvtWeightSums_.size() == 0 ||
            InEvtWeightSums_[0] == 0) {
            if (InEvtCount_) {
                skimAccep_ = std::vector<double>(1, EvtCount_ / InEvtCount_);
            } else {
                std::cout << "Warning: InEvtCount is equal to 0. Event yield normalization might "
                             "be wrong!"
                          << std::endl;
            }
        } else {
            skimAccep_ = std::vector<double>(InEvtWeightSums_.size());
            for (unsigned i = 0; i < InEvtWeightSums_.size() && i < EvtWeightSums_.size(); ++i) {
                skimAccep_[i] = EvtWeightSums_[i] / InEvtWeightSums_[i];
            }
        }
    }
    // cout << "we are here" << "\n";
    //    delete InEvtWeightSums;
    //    delete EvtWeightSums;
}

ZJets::~ZJets()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

void ZJets::CreateOutputFileName(
    const TString &pdfSet, int pdfMember, double muR, double muF, int iJob)
{
    outputFileName = CreateOutputFileName(pdfSet,
                                          pdfMember,
                                          muR,
                                          muF,
                                          iJob,
                                          lepSel,
                                          sampleLabel_,
                                          useTriggerCorrection,
                                          systematics,
                                          direction,
                                          jetPtCutMin,
                                          jetEtaCutMax,
                                          outputDirectory);
}

string ZJets::CreateOutputFileName(const TString &pdfSet,
                                   int pdfMember,
                                   double muR,
                                   double muF,
                                   int iJob,
                                   const TString &lepSel,
                                   const TString &sampleLabel_,
                                   bool useTriggerCorrection,
                                   int systematics,
                                   int direction,
                                   int jetPtCutMin,
                                   int jetEtaCutMax,
                                   TString outputDirectory)
{
    ostringstream result;
    //    result << outputDirectory << fileName;
    result << outputDirectory << lepSel << "_"
           << "13TeV"
           << "_" << sampleLabel_;
    result << "_TrigCorr_" << useTriggerCorrection;
    result << "_Syst_" << systematics;
    if (direction == 1)
        result << "_Up";
    else if (direction == -1)
        result << "_Down";
    result << "_JetPtMin_" << jetPtCutMin;
    result << "_JetEtaMax_" << jetEtaCutMax;
    if (iJob > 0) result << "_" << iJob;

    if (muR != 0 && muF != 0 && muR != 1 && muF != 1) result << "_muR_" << muR << "_muF_" << muF;
    if (pdfSet != "") result << "_PDF_" << pdfSet << "_" << pdfMember;
    if (pdfSet == "" && pdfMember != -1) result << "_NNPDF_" << pdfMember;
    //--- Add your test names here ---
    // result << "_NoPUCut";
    // result << "_LooseID";

    result << ".root";
    return result.str();
}

Int_t ZJets::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t ZJets::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void ZJets::Init(bool hasRecoInfo, bool hasGenInfo)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    //  pdfInfo_ = 0;
    GLepBarePt = 0;
    GLepBareEta = 0;
    GLepBarePhi = 0;
    GLepBareE = 0;
    // genLepQ_ = 0;
    GLepBareId = 0;
    GLepBareSt = 0;
    GLepBarePrompt = 0;
    GLepClosePhotPt = 0;
    GLepClosePhotEta = 0;
    GLepClosePhotPhi = 0;
    GJetAk04Pt = 0;
    GJetAk04Eta = 0;
    GJetAk04Phi = 0;
    GJetAk04E = 0;

    ElPt = 0;
    ElEta = 0;
    ElPhi = 0;
    ElE = 0;
    ElCh = 0;
    ElId = 0;
    // patElecTrig_ = 0;
    ElPfIsoRho = 0;
    ElEtaSc = 0;

    MuPt = 0;
    MuEta = 0;
    MuPhi = 0;
    MuE = 0;
    MuIdTight = 0;
    MuCh = 0;
    //    MuId = 0;
    // patMuonCombId_Double = 0;
    // patMuonTrig_ = 0;
    MuPfIso = 0;
    MuTkLayerCnt = 0;

    JetAk04E = 0;
    JetAk04Pt = 0;
    JetAk04Eta = 0;
    JetAk04Phi = 0;
    JetAk04Id = 0;
    JetAk04PuMva = 0;
    JetAk04BDiscCisvV2 = 0;
    JetAk04PartFlav = 0;

    METPt = 0;
    METPx = 0;
    METPy = 0;
    METsig = 0;
    // mcSherpaWeights_ = 0;
    // weight_amcNLO_ = 0;
    // weight_amcNLO_sum_ = 0;
    EvtWeights = 0;

    //    TrigHlt = 0;
    TrigHltPhot = 0;
    TrigHltMu = 0;
    TrigHltDiMu = 0;
    TrigHltEl = 0;
    TrigHltDiEl = 0;

    // Set branch addresses and branch pointers
    fCurrent = -1;
    fChain->SetMakeClass(1);
    fChain->SetBranchAddress("EvtIsRealData", &EvtIsRealData, &b_EvtIsRealData);
    //    if (fileName.Index("Data") < 0) {
    fChain->SetBranchAddress("EvtPuCntTruth", &EvtPuCntTruth, &b_EvtPuCntTruth);
    fChain->SetBranchAddress("EvtPuCnt", &EvtPuCnt, &b_EvtPuCnt);
    // }
    if (hasRecoInfo) {
        fChain->SetBranchAddress("EvtVtxCnt", &EvtVtxCnt, &b_EvtVtxCnt);
        fChain->SetBranchAddress("EvtRunNum", &EvtRunNum, &b_EvtRunNum);
        fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
        fChain->SetBranchAddress("EvtLumiNum", &EvtLumiNum, &b_EvtLumiNum);
        fChain->SetBranchAddress("JetAk04E", &JetAk04E, &b_JetAk04E);
        fChain->SetBranchAddress("JetAk04Pt", &JetAk04Pt, &b_JetAk04Pt);
        fChain->SetBranchAddress("JetAk04Eta", &JetAk04Eta, &b_JetAk04Eta);
        fChain->SetBranchAddress("JetAk04Phi", &JetAk04Phi, &b_JetAk04Phi);
        fChain->SetBranchAddress("JetAk04Id", &JetAk04Id, &b_JetAk04Id);
        fChain->SetBranchAddress("JetAk04PuMva", &JetAk04PuMva, &b_JetAk04PuMva);
        fChain->SetBranchAddress("JetAk04BDiscCisvV2", &JetAk04BDiscCisvV2, &b_JetAk04BDiscCisvV2);
        // fChain->SetBranchAddress("JetAk04PartFlav", &JetAk04PartFlav, &b_JetAk04PartFlav);
        fChain->SetBranchAddress("METPt", &METPt, &b_METPt);
        fChain->SetBranchAddress("METPx", &METPx, &b_METPx);
        fChain->SetBranchAddress("METPy", &METPy, &b_METPy);
        // fChain->SetBranchAddress("METsig", &METsig, &b_METsig); // not used
        //        fChain->SetBranchAddress("TrigHlt", &TrigHlt, &b_TrigHlt);
        fChain->SetBranchAddress("TrigHltMu", &TrigHltMu, &b_TrigHltMu);
        fChain->SetBranchAddress("TrigHltDiMu", &TrigHltDiMu, &b_TrigHltDiMu);
        fChain->SetBranchAddress("TrigHltEl", &TrigHltEl, &b_TrigHltEl);
        fChain->SetBranchAddress("TrigHltDiEl", &TrigHltDiEl, &b_TrigHltDiEl);

        if (lepSel == "DE" || lepSel == "SE") {
            fChain->SetBranchAddress("ElPt", &ElPt, &b_ElPt);
            fChain->SetBranchAddress("ElEta", &ElEta, &b_ElEta);
            fChain->SetBranchAddress("ElPhi", &ElPhi, &b_ElPhi);
            fChain->SetBranchAddress("ElE", &ElE, &b_ElE);
            fChain->SetBranchAddress("ElCh", &ElCh, &b_ElCh);
            fChain->SetBranchAddress("ElId", &ElId, &b_ElId);
            // fChain->SetBranchAddress("patElecTrig_", &patElecTrig_, &b_patElecTrig_);
            fChain->SetBranchAddress("ElPfIsoRho", &ElPfIsoRho, &b_ElPfIsoRho);
            fChain->SetBranchAddress("ElEtaSc", &ElEtaSc, &b_ElEtaSc);
        }
        if (lepSel == "DMu" || lepSel == "SMu") {
            fChain->SetBranchAddress("MuPt", &MuPt, &b_MuPt);
            fChain->SetBranchAddress("MuEta", &MuEta, &b_MuEta);
            fChain->SetBranchAddress("MuPhi", &MuPhi, &b_MuPhi);
            fChain->SetBranchAddress("MuE", &MuE, &b_MuE);
            fChain->SetBranchAddress("MuCh", &MuCh, &b_MuCh);
            fChain->SetBranchAddress("MuIdTight", &MuIdTight, &b_MuIdTight);
            //            if(fileName.Index("mcatnlo") >= 0 || fileName.Index("MG-MLM") >= 0) {
            //                // CommentAG: before: patMuonCombId_Int
            //                fChain->SetBranchAddress("MuId", &MuId, &b_MuId);
            //            }
            //  else {  // CommentAG: check this
            //      fChain->SetBranchAddress("patMuonCombId_", &patMuonCombId_Double,
            //      &b_patMuonCombId_Double);
            //  }
            //   fChain->SetBranchAddress("patMuonTrig_", &patMuonTrig_, &b_patMuonTrig_);
            fChain->SetBranchAddress("MuPfIso", &MuPfIso, &b_MuPfIso);
            fChain->SetBranchAddress("MuTkLayerCnt", &MuTkLayerCnt, &b_MuTkLayerCnt);
        }
    }
    if (hasGenInfo) {
        fChain->SetBranchAddress("GLepBarePt", &GLepBarePt, &b_GLepBarePt);
        fChain->SetBranchAddress("GLepBareEta", &GLepBareEta, &b_GLepBareEta);
        fChain->SetBranchAddress("GLepBarePhi", &GLepBarePhi, &b_GLepBarePhi);
        fChain->SetBranchAddress("GLepBareE", &GLepBareE, &b_GLepBareE);
        // fChain->SetBranchAddress("genLepQ_", &genLepQ_, &b_genLepQ_);
        fChain->SetBranchAddress("GLepBareId", &GLepBareId, &b_GLepBareId);
        fChain->SetBranchAddress("GLepBareSt", &GLepBareSt, &b_GLepBareSt);
        fChain->SetBranchAddress("GJetAk04Pt", &GJetAk04Pt, &b_GJetAk04Pt);
        fChain->SetBranchAddress("GJetAk04Eta", &GJetAk04Eta, &b_GJetAk04Eta);
        fChain->SetBranchAddress("GJetAk04Phi", &GJetAk04Phi, &b_GJetAk04Phi);
        fChain->SetBranchAddress("GJetAk04E", &GJetAk04E, &b_GJetAk04E);
        fChain->SetBranchAddress("GLepClosePhotPt", &GLepClosePhotPt, &b_GLepClosePhotPt);
        fChain->SetBranchAddress("GLepClosePhotEta", &GLepClosePhotEta, &b_GLepClosePhotEta);
        fChain->SetBranchAddress("GLepClosePhotPhi", &GLepClosePhotPhi, &b_GLepClosePhotPhi);
        if (fChain->GetBranch("GLepBarePrompt")) {
            fChain->SetBranchAddress("GLepBarePrompt", &GLepBarePrompt, &b_GLepClosePhotPhi);
        }
        if (fileName.Index("MIX") >= 0 && fileName.Index("UNFOLDING") >= 0) {
            // fChain->SetBranchAddress("pdfInfo_", &pdfInfo_, &b_pdfInfo_);
            fChain->SetBranchAddress("GNup", &GNup, &b_GNup);
        }
        //  if (fileName.Index("Sherpa") >= 0 && fileName.Index("UNFOLDING") >= 0) {
        //      fChain->SetBranchAddress("mcSherpaWeights_", &mcSherpaWeights_,
        //      &b_mcSherpaWeights_);
        //  }
        // if (fileName.Index("mcatnlo") >= 0) {
        //    fChain->SetBranchAddress("EvtWeights", &EvtWeights, &b_EvtWeights);
        //}
        // if(fileName.Index("Sherpa2") >= 0){
        //    fChain->SetBranchAddress("EvtWeights", &EvtWeights, &b_EvtWeights);
        //}
        // if((fileName.Index("Sherpa") >= 0 && fileName.Index("UNFOLDING") >= 0) ||
        // fileName.Index("mcatnlo") >= 0) {
        //    fChain->SetBranchAddress("weight_amcNLO_", &weight_amcNLO_, &b_weight_amcNLO_);
        //    fChain->SetBranchAddress("weight_amcNLO_sum_", &weight_amcNLO_sum_,
        //    &b_weight_amcNLO_sum_);
        //}
    }

    fChain->SetBranchAddress("EvtWeights", &EvtWeights, &b_EvtWeights);

    Notify();
    cout << "Branches are properly initialized." << endl;
}

Bool_t ZJets::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void ZJets::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

Int_t ZJets::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    printf("entry %lld", entry);
    return 1;
}

bool ZJets::setTriggerMask()
{
    // if (triggerMaskSet_) return triggerMask_;
    if (triggerMaskSet_) return true;
    // NEW DA
    // std::vector<std::string> triggers = cfg.getVS("triggers");
    std::vector<std::string> triggers_EraBF = cfg.getVS("triggers_EraBF");
    std::vector<std::string> triggers_EraGH = cfg.getVS("triggers_EraGH");

    if (triggers_EraBF.size()) {
        std::cout << "An or of the following trigger paths will be used:";
        for (unsigned i = 0; i < triggers_EraBF.size(); ++i) {
            std::cout << " " << triggers_EraBF[i];
        }
        std::cout << "\n";
    } else {
        std::cout << "No trigger requirement will be applied.\n\n";
    }

    if (triggers_EraGH.size()) {
        std::cout << "An or of the following trigger paths will be used:";
        for (unsigned i = 0; i < triggers_EraGH.size(); ++i) {
            std::cout << " " << triggers_EraGH[i];
        }
        std::cout << "\n";
    } else {
        std::cout << "No trigger requirement will be applied.\n\n";
    }

    TString branchName;
    std::vector<std::string> *TrigHlt = 0;

    if (lepSel == "DMu") {
        branchName = "TrigHltDiMu";
        ourTrig_ = &TrigHltDiMu;
    }
    // if(lepSel == "DMu") { branchName = "TrigHltMu"; ourTrig_ = &TrigHltMu;  }
    else if (lepSel == "DE") {
        branchName = "TrigHltDiEl";
        ourTrig_ = &TrigHltDiEl;
    } else if (lepSel == "SMu") {
        branchName = "TrigHltMu";
        ourTrig_ = &TrigHltMu;
    } else if (lepSel == "SE") {
        branchName = "TrigHltEl";
        ourTrig_ = &TrigHltEl;
    } else {
        std::cerr << __FILE__ << ":" << __LINE__ << ". "
                  << "lepSel value, '" << lepSel
                  << "', was  not recognzed. We cannot set the trigger bits.\n\n";
        return false;
    }

    // NEW DA
    if (triggers_EraBF.size() == 0) {
        triggerMask_EraBG = (ULong64_t)-1; // all bits set.
        return true;
    }
    if (triggers_EraGH.size() == 0) {
        triggerMask_EraH = (ULong64_t)-1; // all bits set.
        return true;
    }

    if (fBitFieldsChain.GetBranch(branchName)) {
        fBitFieldsChain.SetBranchAddress(branchName, &TrigHlt);
    } else {
        std::cerr << "Cannot set the trigger bits, because the Branch " << branchName
                  << " was not found in the tree BitFields.\n\n";
        return false;
    }

    // FIXME: check all files?
    if (fBitFieldsChain.GetEntry(0) <= 0) {
        std::cerr << "Failed to read BitFields tree. Is the tree empty? "
                  << "Cannot set the trigger bits.\n\n";
        return false;
    }

    std::cout << "Available triggers for " << branchName << ":";
    for (unsigned i = 0; i < TrigHlt->size(); ++i) {
        if ((*TrigHlt)[i].size()) std::cout << " " << (*TrigHlt)[i];
    }
    std::cout << "\n";

    // NEW DA
    bool rc = true;
    triggerMask_EraBG = 0;
    triggerMask_EraH = 0;
    for (unsigned i = 0; i < triggers_EraBF.size(); ++i) {
        int found = 0;
        for (unsigned bit = 0; bit < TrigHlt->size(); ++bit) {
            if (triggers_EraBF[i] == (*TrigHlt)[bit]) {
                triggerMask_EraBG |= (1 << bit);
                ++found;
            }
        }
        if (found == 0) {
            std::cerr << "Trigger " << triggers_EraBF[i]
                      << " was not found in the input sample!\n\n";
            rc = rc && false;
        }
        if (found > 1) {
            std::cerr << "Trigger " << triggers_EraBF[i] << " is assigned to several bits!\n\n";
            rc = rc && true;
        }
    }

    for (unsigned i = 0; i < triggers_EraGH.size(); ++i) {
        int found = 0;
        for (unsigned bit = 0; bit < TrigHlt->size(); ++bit) {
            if (triggers_EraGH[i] == (*TrigHlt)[bit]) {
                triggerMask_EraH |= (1 << bit);
                ++found;
            }
        }
        if (found == 0) {
            std::cerr << "Trigger " << triggers_EraGH[i]
                      << " was not found in the input sample!\n\n";
            rc = rc && false;
        }
        if (found > 1) {
            std::cerr << "Trigger " << triggers_EraGH[i] << " is assigned to several bits!\n\n";
            rc = rc && true;
        }
    }

    return rc;
}

// bool ZJets::compTrigger(const char* a, const char* bv) const{
//  int i = 0;
//  for(;a[i]!=0 && bv[i]!=0; ++i){
//    if(a[i]!= bv[i]) return false;
//  }
//  if(a[i]) return false;
//  if(bv[i]==0) return true;
//  if(bv[i] != '_') return false;
//  if(bv[++i]!='v') return false;
//  for(;;){
//    if(bv[++i]==0) return true;
//    if(!isdigit(bv[i])) return false;
//  }
//  return true;
//};
