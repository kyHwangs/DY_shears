// -*- c-basic-offset:8; -*-
#ifndef _HistoSetZJets_h_
#define _HistoSetZJets_h_

#include "TArray.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TNamed.h"

#include "GenH1D.h"
#include <cstdarg>
#include <iostream>
#include <set>
#include <vector>
#define JETPT_ZINC1_BINS 9
#define JETPT_20_ZINC1_BINS 11

using namespace std;

class HistoSetZJets
{
  public:
    vector<TH1 *> listOfHistograms;
    HistoSetZJets(TString leptonFlavor = "DMu");
    ~HistoSetZJets();

    vector<double> makeVector(int num, ...);
    void insertVector(vector<double> &veca, int num, ...);
    vector<double> buildVecFineBin(int nStdBin, double arrStdBin[], int factChop);

    // List of actually used histograms;
    std::set<std::string> varList;

    GenH1D *newTH1D(string, string, string, int, double *);
    GenH1D *newTH1D(string, string, string, int, double, double);
    GenH1D *newTH1D(string, string, string, vector<double> &);
    TH2D *newTH2D(string, string, int, double *, int, double *);
    TH2D *newTH2D(string, string, vector<double> &, vector<double> &);
    TH2D *newTH2D(string, string, int, double *, int, double, double);
    TH2D *newTH2D(string, string, int, double, double, int, double *);
    TH2D *newTH2D(string, string, int, double, double, int, double, double);

    void readHistList();
    void writeHistList() const;

    bool filterHist(const char *name) const;

    // bins indices of JobInfo histogram:
    enum { kJobNum = 1, kNJobs, kJobWeight, kNEvts, kNEvtsAllJobs, kNEvtsSample, kXsec, kLumi };

    TH1D *Input;
    TH1D *JobInfo;
    TH1D *Lumi;
    TH1D *NumberPFcandidates;
    TH1D *NumberOfEvents;
    TH1D *ZMass_lowDeltaR;
    TH1D *ZMass_Zinc0jet;
    TH1D *ZMass_Zinc1jet;
    TH1D *ZMassFrom60_Zinc0jet;
    GenH1D *genZMass_Zinc0jet;
    TH1D *Mass_Zinc0jet;
    TH1D *Phistar_Zinc0jetM50_76;
    TH1D *Phistar_Zinc0jetM76_106;
    TH1D *Phistar_Zinc0jetM106_170;
    TH1D *Phistar_Zinc0jetM170_350;
    TH1D *Phistar_Zinc0jetM170_inf;
    TH1D *Phistar_Zinc0jet;
    TH1D *Phistar_Zinc0jetM111_130;
    TH1D *Phistar_Zinc0jetM130_170;
    TH1D *Phistar_Zinc0jetM170_250;
    TH1D *Phistar_Zinc0jetM250_3;
    GenH1D *genPhistar_Zinc0jetM50_76;
    GenH1D *genPhistar_Zinc0jetM76_106;
    GenH1D *genPhistar_Zinc0jetM106_170;
    GenH1D *genPhistar_Zinc0jetM170_350;
    GenH1D *genPhistar_Zinc0jetM170_inf;
    GenH1D *genPhistar_Zinc0jet;
    GenH1D *genPhistar_Zinc0jetM111_130;
    GenH1D *genPhistar_Zinc0jetM130_170;
    GenH1D *genPhistar_Zinc0jetM170_250;
    GenH1D *genPhistar_Zinc0jetM250_3;
    TH2D *hresponsePhistar_Zinc0jetM50_76;
    TH2D *hresponsePhistar_Zinc0jetM76_106;
    TH2D *hresponsePhistar_Zinc0jetM50_106;
    TH2D *hresponsePhistar_Zinc0jetM106_170;
    TH2D *hresponsePhistar_Zinc0jetM170_350;
    TH2D *hresponsePhistar_Zinc0jetM170_inf;
    TH2D *hresponsePhistar_Zinc0jet;
    TH2D *hresponsePhistar_Zinc0jetM111_130;
    TH2D *hresponsePhistar_Zinc0jetM130_170;
    TH2D *hresponsePhistar_Zinc0jetM170_250;
    TH2D *hresponsePhistar_Zinc0jetM250_3;
    TH1D *Phistar_Zinc1jetM50_76;
    TH1D *Phistar_Zinc1jet;
    TH1D *Phistar_Zinc1jetM111_130;
    TH1D *Phistar_Zinc1jetM130_170;
    TH1D *Phistar_Zinc1jetM170_250;
    TH1D *Phistar_Zinc1jetM250_3;
    GenH1D *genPhistar_Zinc1jetM50_76;
    GenH1D *genPhistar_Zinc1jet;
    GenH1D *genPhistar_Zinc1jetM111_130;
    GenH1D *genPhistar_Zinc1jetM130_170;
    GenH1D *genPhistar_Zinc1jetM170_250;
    GenH1D *genPhistar_Zinc1jetM250_3;

    GenH1D *genZPt_Zinc0jetM15_50;
    GenH1D *genZPt_Zinc0jetM50_76;
    GenH1D *genZPt_Zinc0jetM111_130;
    GenH1D *genZPt_Zinc0jetM130_170;
    GenH1D *genZPt_Zinc0jetM250_3;
    GenH1D *genZPt_2_Zinc0jetM250_3;
    GenH1D *genZPt_Zinc0jetM320_3;
    GenH1D *genZPt_Zinc0jetM170_250;
    GenH1D *genZPt_2_Zinc0jetM170_250;
    GenH1D *genZPt_Zinc0jetM170_250_new;

    TH2D *hresponsePhistar_Zinc1jetM50_76;
    TH2D *hresponsePhistar_Zinc1jet;
    TH2D *hresponsePhistar_Zinc1jetM111_130;
    TH2D *hresponsePhistar_Zinc1jetM130_170;
    TH2D *hresponsePhistar_Zinc1jetM170_250;
    TH2D *hresponsePhistar_Zinc1jetM250_3;

    TH2D *hresponseZPt_Zinc0jetM15_50;
    TH2D *hresponseZPt_Zinc0jetM50_76;
    TH2D *hresponseZPt_Zinc0jetM111_130;
    TH2D *hresponseZPt_Zinc0jetM130_170;
    TH2D *hresponseZPt_Zinc0jetM250_3;
    TH2D *hresponseZPt_2_Zinc0jetM250_3;
    TH2D *hresponseZPt_Zinc0jetM320_3;
    TH2D *hresponseZPt_Zinc0jetM170_250;
    TH2D *hresponseZPt_2_Zinc0jetM170_250;
    TH2D *hresponseZPt_Zinc0jetM170_250_new;
    TH2D *hresponseZPt_Zinc0jetM250_3_dR;
    TH2D *hresponseZPt_Zinc0jetM250_3_dPt;
    TH2D *hresponseZPt_Zinc0jetM170_inf;

    TH2D *hresponselep0Pt_Zinc0jet;
    TH2D *hresponselep1Pt_Zinc0jet;

    TH1D *ZPt_Zinc0jet;
    TH1D *ZPt_Zinc0jet_new;
    TH1D *ZPt_Zinc0jetM115_135;
    TH1D *ZPt_Zinc0jetM130_3;
    TH1D *ZPt_Zinc0jetM170_inf;

    TH1D *ZPt_Zinc0jetM15_50;
    TH1D *ZPt_Zinc0jetM50_76;
    TH1D *ZPt_Zinc0jetM106_170;
    TH1D *ZPt_Zinc0jetM76_106;
    TH1D *ZPt_Zinc0jetM76_106_Mbin;
    TH1D *ZPt_Zinc0jetM170_350;
    TH1D *ZPt_Zinc0jetM111_130;
    TH1D *ZPt_Zinc0jetM130_170;
    TH1D *ZPt_Zinc0jetM250_3;
    TH1D *ZPt_2_Zinc0jetM250_3;
    TH1D *ZPt_Zinc0jetM320_3;
    TH1D *ZPt_Zinc0jetM170_250;
    TH1D *ZPt_2_Zinc0jetM170_250;
    TH1D *ZPt_Zinc0jetM170_250_new;

    TH1D *ZPt_Zinc0jet_Odd;
    TH1D *ZPt_Zinc0jet_Even;
    TH1D *ZPt_Zinc1jetM50_76;
    TH1D *ZPt_Zinc1jetM76_106;
    TH1D *ZPt_Zinc1jetM76_106_Mbin;
    TH1D *ZPt_Zinc1jetM106_170;
    TH1D *ZPt_Zinc1jetM170_350;
    TH1D *ZPt_Zinc1jetM170_inf;
    TH1D *ZPt_Zinc1jet;
    TH1D *ZPt_Zinc1jetM111_130;
    TH1D *ZPt_Zinc1jetM130_170;
    TH1D *ZPt_Zinc1jetM170_250;
    TH1D *ZPt_Zinc1jetM250_3;

    TH1D *ZPt_Zinc1jet_Odd;
    ;
    TH1D *ZPt_Zinc1jet_Even;
    TH1D *HadRecoil;
    TH1D *JZB;
    TH1D *JZB_Odd;
    ;
    TH1D *JZB_Even;
    TH1D *JZB_2;
    TH1D *JZB_ptLow;
    TH1D *JZB_ptLow_Odd;
    ;
    TH1D *JZB_ptLow_Even;
    TH1D *JZB_ptLow_2;
    TH1D *JZB_ptHigh;
    TH1D *JZB_ptHigh_Odd;
    ;
    TH1D *JZB_ptHigh_Even;
    TH1D *JZB_ptHigh_2;
    TH1D *ZPt_Zinc2jet;
    TH1D *VisPt_Zinc0jetQun;
    TH1D *VisPt_Zinc0jetQun_Odd;
    ;
    TH1D *VisPt_Zinc0jetQun_Even;
    TH1D *VisPt_2_Zinc0jetQun;
    GenH1D *genVisPt_Zinc0jetQun;
    TH1D *VisPt_Zinc1jetQun;
    TH1D *VisPt_Zinc1jetQun_Odd;
    ;
    TH1D *VisPt_Zinc1jetQun_Even;
    TH1D *VisPt_2_Zinc1jetQun;
    GenH1D *genVisPt_Zinc1jetQun;
    TH1D *VisPt_Zinc2jetQun;
    TH1D *VisPt_Zinc2jetQun_Odd;
    ;
    TH1D *VisPt_Zinc2jetQun_Even;
    TH1D *VisPt_2_Zinc2jetQun;
    GenH1D *genVisPt_Zinc2jetQun;
    TH1D *VisPt_Zinc3jetQun;
    TH1D *VisPt_Zinc3jetQun_Odd;
    ;
    TH1D *VisPt_Zinc3jetQun_Even;
    TH1D *VisPt_2_Zinc3jetQun;
    GenH1D *genVisPt_Zinc3jetQun;
    GenH1D *genZPt_Zinc0jet;
    GenH1D *genZPt_Zinc0jetM106_170;
    GenH1D *genZPt_Zinc0jetM76_106;
    GenH1D *genZPt_Zinc0jetM76_106_Mbin;
    GenH1D *genZPt_Zinc0jetM170_inf;
    GenH1D *genZPt_Zinc0jetM170_350;
    GenH1D *genZPt_Zinc0jet_new;
    GenH1D *genZPt_Zinc0jetM115_135;
    GenH1D *genZPt_Zinc0jetM130_3;
    TH2D *hresponseZPt_Zinc0jet;
    TH2D *hresponseZPt_Zinc0jetM76_106;
    TH2D *hresponseZPt_Zinc0jetM76_106_Mbin;
    TH2D *hresponseZPt_Zinc0jetM106_170;
    TH2D *hresponseZPt_Zinc0jetM170_350;
    TH2D *hresponseZPt_Zinc0jet_new;
    TH2D *hresponseZPt_Zinc0jet_lowNVtx;
    TH2D *hresponseZPt_Zinc0jet_highNVtx;
    TH2D *hresponseZPt_Zinc0jetM115_135;
    TH2D *hresponseZPt_Zinc0jetM130_3;
    GenH1D *genZPt_Zinc1jetM50_76;
    GenH1D *genZPt_Zinc1jetM76_106;
    GenH1D *genZPt_Zinc1jetM76_106_Mbin;
    GenH1D *genZPt_Zinc1jetM106_170;
    GenH1D *genZPt_Zinc1jetM170_350;
    GenH1D *genZPt_Zinc1jetM170_inf;
    GenH1D *genZPt_Zinc1jet;
    GenH1D *genZPt_Zinc1jetM111_130;
    GenH1D *genZPt_Zinc1jetM130_170;
    GenH1D *genZPt_Zinc1jetM170_250;
    GenH1D *genZPt_Zinc1jetM250_3;
    TH2D *hresponseZPt_Zinc1jetM50_76;
    TH2D *hresponseZPt_Zinc1jetM76_106;
    TH2D *hresponseZPt_Zinc1jetM76_106_Mbin;
    TH2D *hresponseZPt_Zinc1jetM106_170;
    TH2D *hresponseZPt_Zinc1jetM170_350;
    TH2D *hresponseZPt_Zinc1jetM170_inf;
    TH2D *hresponseZPt_Zinc1jet;
    TH2D *hresponseZPt_Zinc1jetM111_130;
    TH2D *hresponseZPt_Zinc1jetM130_170;
    TH2D *hresponseZPt_Zinc1jetM170_250;
    TH2D *hresponseZPt_Zinc1jetM250_3;
    GenH1D *genHadRecoil;
    GenH1D *genJZB;
    GenH1D *genJZB_ptLow;
    GenH1D *genJZB_ptHigh;
    GenH1D *genZPt_Zinc2jet;
    TH1D *ZPt_Zexc0jet;
    TH1D *ZPt_Zexc1jet;
    TH1D *ZPt_Zexc2jet;
    TH1D *ZRapidity_Zinc0jet;
    TH1D *ZRapidity_Zinc1jet;
    TH1D *ZAbsRapidity_Zinc1jet;
    TH1D *ZRapidity_Zinc2jet;
    GenH1D *genZRapidity_Zinc0jet;
    GenH1D *genZRapidity_Zinc1jet;
    GenH1D *genZAbsRapidity_Zinc1jet;
    GenH1D *genZRapidity_Zinc2jet;
    TH1D *ZRapidity_Zexc0jet;
    TH1D *ZRapidity_Zexc1jet;
    TH1D *ZRapidity_Zexc2jet;

    TH1D *ZEta_Zinc0jet;
    TH1D *ZEtaUpTo5_Zinc0jet;
    TH1D *ZEta_Zinc1jet;
    TH1D *ZEtaUpTo5_Zinc1jet;
    TH1D *ZEta_Zinc2jet;

    GenH1D *genZEta_Zinc0jet;
    GenH1D *genZEta_Zinc1jet;
    GenH1D *genZEta_Zinc2jet;

    TH1D *ZEta_Zexc0jet;
    TH1D *ZEta_Zexc1jet;
    TH1D *ZEta_Zexc2jet;

    TH1D *lepEta_Zinc0jet;
    TH1D *lepLEta_Zinc0jet;
    TH1D *lepSEta_Zinc0jet;
    TH1D *lepEta_Zinc1jet;
    TH1D *lepEtaUpTo4_Zinc0jet;
    GenH1D *genlepEta_Zinc0jet;

    TH1D *lepEta_Zexc0jet;

    TH1D *lepPhi_Zinc0jet;
    TH1D *lepPhi_Zexc0jet;

    TH1D *FirstJetEta_Zinc1jet;
    TH1D *FirstJetEta_Zinc1jet_res08;
    TH1D *FirstJetEta_Zinc1jet_res16;
    TH1D *FirstJetEta_Zinc1jet_res24;
    TH2D *FirstJetEta_Zinc1jet_res2D;

    TH1D *FirstJetEta_2_Zinc1jet;
    TH1D *SecondJetEta_Zinc2jet;
    TH1D *SecondJetEta_2_Zinc2jet;
    TH1D *ThirdJetEta_Zinc3jet;
    TH1D *ThirdJetEta_2_Zinc3jet;
    TH1D *FourthJetEta_Zinc4jet;
    TH1D *FifthJetEta_Zinc5jet;
    TH1D *SixthJetEta_Zinc6jet;

    TH1D *FirstJetEtaHigh_Zinc1jet;
    TH1D *SecondJetEtaHigh_Zinc2jet;
    TH1D *ThirdJetEtaHigh_Zinc3jet;
    TH1D *FourthJetEtaHigh_Zinc4jet;
    TH1D *FifthJetEtaHigh_Zinc5jet;
    TH1D *SixthJetEtaHigh_Zinc6jet;

    TH2D *FirstJetAbsYvsAbsEta_Zinc1jet;

    TH1D *FirstJetAbsRapidity_Zinc1jet;
    TH1D *FirstJetAbsRapidity_SmearMatch_Zinc1jet; // DJALOG
    TH1D *FirstJetAbsRapidity_SmearGauss_Zinc1jet; // DJALOG
    TH1D *FirstJetAbsRapidity_Zinc1jet_Odd;
    TH1D *FirstJetAbsRapidity_Zinc1jet_Even;
    TH1D *FirstJetAbsRapidity_2_Zinc1jet;
    TH1D *FirstJetRapidityHigh_Zinc1jet;

    TH1D *SecondJetAbsRapidity_Zinc2jet;
    TH1D *SecondJetAbsRapidity_Zinc2jet_Odd;
    TH1D *SecondJetAbsRapidity_Zinc2jet_Even;
    TH1D *SecondJetAbsRapidity_2_Zinc2jet;
    TH1D *SecondJetRapidityHigh_Zinc2jet;

    TH1D *ThirdJetAbsRapidity_Zinc3jet;
    TH1D *ThirdJetAbsRapidity_Zinc3jet_Odd;
    TH1D *ThirdJetAbsRapidity_Zinc3jet_Even;
    TH1D *ThirdJetAbsRapidity_2_Zinc3jet;
    TH1D *ThirdJetRapidityHigh_Zinc3jet;

    TH1D *FourthJetAbsRapidity_Zinc4jet;
    TH1D *FourthJetRapidityHigh_Zinc4jet;
    TH1D *FifthJetAbsRapidity_Zinc5jet;
    TH1D *FifthJetRapidityHigh_Zinc5jet;
    TH1D *SixthJetAbsRapidity_Zinc6jet;
    TH1D *SixthJetRapidityHigh_Zinc6jet;

    TH1D *FirstJetEtaFull_Zinc1jet;
    TH1D *SecondJetEtaFull_Zinc2jet;
    TH1D *ThirdJetEtaFull_Zinc3jet;
    TH1D *FourthJetEtaFull_Zinc4jet;
    TH1D *FifthJetEtaFull_Zinc5jet;
    TH1D *SixthJetEtaFull_Zinc6jet;

    GenH1D *genFirstJetEta_Zinc1jet;
    GenH1D *genFirstJetEta_2_Zinc1jet;
    GenH1D *genSecondJetEta_Zinc2jet;
    GenH1D *genSecondJetEta_2_Zinc2jet;
    GenH1D *genThirdJetEta_Zinc3jet;
    GenH1D *genThirdJetEta_2_Zinc3jet;
    GenH1D *genFourthJetEta_Zinc4jet;
    GenH1D *genFifthJetEta_Zinc5jet;
    GenH1D *genSixthJetEta_Zinc6jet;

    GenH1D *genFirstJetEtaHigh_Zinc1jet;
    GenH1D *genSecondJetEtaHigh_Zinc2jet;
    GenH1D *genThirdJetEtaHigh_Zinc3jet;
    GenH1D *genFourthJetEtaHigh_Zinc4jet;
    GenH1D *genFifthJetEtaHigh_Zinc5jet;
    GenH1D *genSixthJetEtaHigh_Zinc6jet;

    GenH1D *genFirstJetRapidityHigh_Zinc1jet;
    GenH1D *genFirstJetAbsRapidity_Zinc1jet;
    GenH1D *genSecondJetRapidityHigh_Zinc2jet;
    GenH1D *genSecondJetAbsRapidity_Zinc2jet;
    GenH1D *genThirdJetRapidityHigh_Zinc3jet;
    GenH1D *genThirdJetAbsRapidity_Zinc3jet;
    GenH1D *genFourthJetRapidityHigh_Zinc4jet;
    GenH1D *genFourthJetAbsRapidity_Zinc4jet;
    GenH1D *genFifthJetRapidityHigh_Zinc5jet;
    GenH1D *genFifthJetAbsRapidity_Zinc5jet;
    GenH1D *genSixthJetRapidityHigh_Zinc6jet;
    GenH1D *genSixthJetAbsRapidity_Zinc6jet;

    TH1D *FirstJetEta_Zexc1jet;
    TH1D *SecondJetEta_Zexc2jet;

    TH1D *FirstJetPhi_Zinc1jet;
    TH2D *FirstJetPt_Zinc1jet_NVtx;
    TH2D *FirstJetPtRecoOvGen_Zinc1jet_NVtx;
    TH2D *lepResolution_pt;
    TH2D *lepResolution_nvtx;
    TH2D *lepResolution_pt_rel;
    TH2D *Phistar_Zpt;
    TH2D *Phistar_Zpt_test;

    TH1D *SecondJetPhi_Zinc2jet;
    TH1D *ThirdJetPhi_Zinc3jet;
    TH1D *FourthJetPhi_Zinc4jet;
    TH1D *FifthJetPhi_Zinc5jet;
    TH1D *SixthJetPhi_Zinc6jet;
    TH1D *FirstJetPhi_Zexc1jet;
    TH1D *SecondJetPhi_Zexc2jet;
    TH1D *lepPt_Zinc0jet;
    TH1D *lepLPt_Zinc0jet;
    TH1D *lepSPt_Zinc0jet;
    TH1D *lepPt_Zinc1jet;
    TH1D *lepLPt_Zinc1jet;
    TH1D *lepSPt_Zinc1jet;
    TH1D *lepPtFrom15_Zinc0jet;
    GenH1D *genlepPt_Zinc0jet;
    TH1D *lepPt_Zexc0jet;
    TH1D *lepLPt_Zexc0jet;
    TH1D *lepSPt_Zexc0jet;
    TH1D *dPhiLeptons_Zexc0jet;
    TH1D *dPhiLeptons_Zinc0jet;
    TH1D *dPhiLeptons_Zinc1jet;
    TH1D *dEtaLeptons_Zexc0jet;
    TH1D *dEtaLeptons_Zinc0jet;
    TH1D *dEtaLeptons_Zinc1jet;
    TH1D *dRLeptons_Zinc0jet;
    TH1D *dRLeptons_Zinc1jet;
    TH1D *SpTLeptons_Zexc0jet;
    TH1D *SpTLeptons_Zexc1jet;
    TH1D *SpTLeptons_Zexc2jet;
    GenH1D *genSpTLeptons_Zexc2jet;
    TH1D *SpTLeptons_Zinc0jet;
    TH1D *SpTLeptons_Zinc1jet;
    TH1D *SpTLeptons_Zinc2jet;
    GenH1D *genSpTLeptons_Zinc2jet;

    TH1D *JetsHT_Zinc1jet;
    TH1D *JetsHT_Zinc1jet_Odd;
    TH1D *JetsHT_Zinc1jet_Even;
    TH1D *JetsHT_2_Zinc1jet;
    TH1D *JetsHT_Zinc2jet;
    TH1D *JetsHT_Zinc2jet_Odd;
    TH1D *JetsHT_Zinc2jet_Even;
    TH1D *JetsHT_2_Zinc2jet;
    TH1D *JetsHT_Zinc3jet;
    TH1D *JetsHT_Zinc3jet_Odd;
    TH1D *JetsHT_Zinc3jet_Even;
    TH1D *JetsHT_2_Zinc3jet;
    TH1D *JetsHT_Zinc4jet;
    TH1D *JetsHT_Zinc5jet;
    TH1D *JetsHT_Zinc6jet;

    GenH1D *genJetsHT_Zinc1jet;
    GenH1D *genJetsHT_2_Zinc1jet;
    GenH1D *genJetsHT_Zinc2jet;
    GenH1D *genJetsHT_2_Zinc2jet;
    GenH1D *genJetsHT_Zinc3jet;
    GenH1D *genJetsHT_2_Zinc3jet;
    GenH1D *genJetsHT_Zinc4jet;
    GenH1D *genJetsHT_Zinc5jet;
    GenH1D *genJetsHT_Zinc6jet;

    TH1D *FirstJetPt_Zinc1jet;
    TH1D *FirstJetPt_SmearMatch_Zinc1jet; // DJALOG
    TH1D *FirstJetPt_SmearGauss_Zinc1jet; // DJALOG
    TH1D *FirstJetPt_Zinc1jet_Odd;
    TH1D *FirstJetPt_Zinc1jet_Even;
    TH1D *FirstJetPt_2_Zinc1jet;

    TH1D *SecondJetPt_Zinc2jet;
    TH1D *SecondJetPt_Zinc2jet_Odd;
    TH1D *SecondJetPt_Zinc2jet_Even;
    TH1D *SecondJetPt_2_Zinc2jet;

    TH1D *ThirdJetPt_Zinc3jet;
    TH1D *ThirdJetPt_Zinc3jet_Odd;
    TH1D *ThirdJetPt_Zinc3jet_Even;

    TH1D *ThirdJetPt_2_Zinc3jet;
    TH1D *FourthJetPt_Zinc4jet;
    TH1D *FifthJetPt_Zinc5jet;
    TH1D *SixthJetPt_Zinc6jet;

    GenH1D *genFirstJetPt_Zinc1jet;
    GenH1D *genFirstJetPt_2_Zinc1jet;
    GenH1D *genSecondJetPt_Zinc2jet;
    GenH1D *genSecondJetPt_2_Zinc2jet;
    GenH1D *genThirdJetPt_Zinc3jet;
    GenH1D *genThirdJetPt_2_Zinc3jet;
    GenH1D *genFourthJetPt_Zinc4jet;
    GenH1D *genFifthJetPt_Zinc5jet;
    GenH1D *genSixthJetPt_Zinc6jet;
    TH1D *FirstJetPt_Zexc1jet;
    TH1D *SecondJetPt_Zexc2jet;
    GenH1D *genFirstJetPt_Zexc1jet;
    GenH1D *genSecondJetPt_Zexc2jet;

    TH1D *ZNGoodJets_Zexc;
    TH1D *ZNGoodJets_Zexc_Odd;
    TH1D *ZNGoodJets_Zexc_Even;
    TH2D *ZNGoodJetsNVtx_Zexc;
    TH2D *ZNGoodJets20NVtx_Zexc;
    TH1D *ZNGoodJets_Zinc;
    TH1D *ZNGoodJets_Zinc_nvtx10;
    TH1D *ZNGoodJets_Zinc_nvtx20;
    TH1D *ZNGoodJets_Zinc_nvtx30;
    TH1D *ZNGoodJets_Zinc_nvtx45;
    TH1D *ZNGoodJets_Zinc_5;
    TH1D *ZNGoodJets_Zinc_10;
    TH1D *ZNGoodJets_Zinc_15;
    TH1D *ZNGoodJets_Zinc_20;
    TH1D *ZNGoodJets_Zinc_25;
    TH1D *ZNGoodJets_Zexc_NoWeight;
    TH1D *ZNGoodJets_Zinc_NoWeight;
    TH1D *ZNGoodJets_SameChargePair_Zinc;

    TH1D *ptBal_Zinc1jet;
    TH1D *ptBal_Zinc3jet;
    TH1D *TwoJetsPtDiff_Zexc2jet;
    GenH1D *genTwoJetsPtDiff_Zexc2jet;
    TH1D *JetsMass_Zexc2jet;
    GenH1D *genJetsMass_Zexc2jet;
    TH1D *ptBal_Zexc2jet;
    GenH1D *genptBal_Zexc2jet;
    TH1D *dPhiJets_Zexc2jet;
    GenH1D *gendPhiJets_Zexc2jet;
    TH1D *dEtaJets_Zexc2jet;
    GenH1D *gendEtaJets_Zexc2jet;
    TH1D *dEtaFirstJetZ_Zexc2jet;
    GenH1D *gendEtaFirstJetZ_Zexc2jet;
    TH1D *dEtaSecondJetZ_Zexc2jet;
    GenH1D *gendEtaSecondJetZ_Zexc2jet;
    TH1D *dEtaJet1Plus2Z_Zexc2jet;
    GenH1D *gendEtaJet1Plus2Z_Zexc2jet;
    TH1D *PHI_Zexc2jet;
    GenH1D *genPHI_Zexc2jet;
    TH1D *PHI_T_Zexc2jet;
    GenH1D *genPHI_T_Zexc2jet;
    TH1D *SpT_Zexc2jet;
    GenH1D *genSpT_Zexc2jet;
    TH1D *SpTJets_Zexc2jet;
    GenH1D *genSpTJets_Zexc2jet;
    TH1D *SPhi_Zexc2jet;
    GenH1D *genSPhi_Zexc2jet;

    TH1D *TwoJetsPtDiff_Zinc2jet;
    GenH1D *genTwoJetsPtDiff_Zinc2jet;
    TH1D *BestTwoJetsPtDiff_Zinc2jet;
    GenH1D *genBestTwoJetsPtDiff_Zinc2jet;

    TH1D *llJetsMass_Zinc2jet;
    GenH1D *genllJetsMass_Zinc2jet;

    TH1D *JetsMass_Zinc2jet;
    GenH1D *genJetsMass_Zinc2jet;

    TH1D *JetsMassLowPU_Zinc2jet;
    GenH1D *genJetsMassLowPU_Zinc2jet;

    TH1D *JetsMassMidPU_Zinc2jet;
    GenH1D *genJetsMassMidPU_Zinc2jet;

    TH1D *JetsMassHigPU_Zinc2jet;
    GenH1D *genJetsMassHigPU_Zinc2jet;

    TH1D *BestJetsMass_Zinc2jet;
    GenH1D *genBestJetsMass_Zinc2jet;
    TH1D *ptBal_Zinc2jet;
    GenH1D *genptBal_Zinc2jet;
    TH1D *dPhiJets_Zinc2jet;
    GenH1D *gendPhiJets_Zinc2jet;
    TH1D *BestdPhiJets_Zinc2jet;
    GenH1D *genBestdPhiJets_Zinc2jet;
    TH1D *dEtaJets_Zinc2jet;
    GenH1D *gendEtaJets_Zinc2jet;
    TH1D *dEtaFirstJetZ_Zinc2jet;
    GenH1D *gendEtaFirstJetZ_Zinc2jet;
    TH1D *dEtaSecondJetZ_Zinc2jet;
    GenH1D *gendEtaSecondJetZ_Zinc2jet;
    TH1D *dEtaJet1Plus2Z_Zinc2jet;
    GenH1D *gendEtaJet1Plus2Z_Zinc2jet;
    TH1D *PHI_Zinc2jet;
    GenH1D *genPHI_Zinc2jet;
    TH1D *BestPHI_Zinc2jet;
    GenH1D *genBestPHI_Zinc2jet;
    TH1D *PHI_T_Zinc2jet;
    GenH1D *genPHI_T_Zinc2jet;
    TH1D *BestPHI_T_Zinc2jet;
    GenH1D *genBestPHI_T_Zinc2jet;
    TH1D *SpT_Zinc2jet;
    GenH1D *genSpT_Zinc2jet;
    TH1D *BestSpT_Zinc2jet;
    GenH1D *genBestSpT_Zinc2jet;
    TH1D *SpTJets_Zinc2jet;
    GenH1D *genSpTJets_Zinc2jet;
    TH1D *BestSpTJets_Zinc2jet;
    GenH1D *genBestSpTJets_Zinc2jet;
    TH1D *SPhi_Zinc2jet;
    GenH1D *genSPhi_Zinc2jet;
    TH1D *BestSPhi_Zinc2jet;
    GenH1D *genBestSPhi_Zinc2jet;

    //-- low Z pT;
    TH1D *TwoJetsPtDiff_LowPt_Zexc2jet;
    GenH1D *genTwoJetsPtDiff_LowPt_Zexc2jet;
    TH1D *JetsMass_LowPt_Zexc2jet;
    GenH1D *genJetsMass_LowPt_Zexc2jet;
    TH1D *ptBal_LowPt_Zexc2jet;
    GenH1D *genptBal_LowPt_Zexc2jet;
    TH1D *dPhiJets_LowPt_Zexc2jet;
    GenH1D *gendPhiJets_LowPt_Zexc2jet;
    TH1D *dPhiLeptons_LowPt_Zexc2jet;
    GenH1D *gendPhiLeptons_LowPt_Zexc2jet;
    TH1D *PHI_LowPt_Zexc2jet;
    GenH1D *genPHI_LowPt_Zexc2jet;
    TH1D *PHI_T_LowPt_Zexc2jet;
    GenH1D *genPHI_T_LowPt_Zexc2jet;
    TH1D *SpT_LowPt_Zexc2jet;
    GenH1D *genSpT_LowPt_Zexc2jet;
    TH1D *SpTJets_LowPt_Zexc2jet;
    GenH1D *genSpTJets_LowPt_Zexc2jet;
    TH1D *SpTLeptons_LowPt_Zexc2jet;
    GenH1D *genSpTLeptons_LowPt_Zexc2jet;
    TH1D *SPhi_LowPt_Zexc2jet;
    GenH1D *genSPhi_LowPt_Zexc2jet;

    TH1D *TwoJetsPtDiff_LowPt_Zinc2jet;
    GenH1D *genTwoJetsPtDiff_LowPt_Zinc2jet;
    TH1D *BestTwoJetsPtDiff_LowPt_Zinc2jet;
    GenH1D *genBestTwoJetsPtDiff_LowPt_Zinc2jet;

    TH1D *JetsMass_LowPt_Zinc2jet;
    GenH1D *genJetsMass_LowPt_Zinc2jet;
    TH1D *BestJetsMass_LowPt_Zinc2jet;
    GenH1D *genBestJetsMass_LowPt_Zinc2jet;
    TH1D *ptBal_LowPt_Zinc2jet;
    GenH1D *genptBal_LowPt_Zinc2jet;
    TH1D *dPhiJets_LowPt_Zinc2jet;
    GenH1D *gendPhiJets_LowPt_Zinc2jet;
    TH1D *BestdPhiJets_LowPt_Zinc2jet;
    GenH1D *genBestdPhiJets_LowPt_Zinc2jet;
    TH1D *dPhiLeptons_LowPt_Zinc2jet;
    GenH1D *gendPhiLeptons_LowPt_Zinc2jet;
    TH1D *PHI_LowPt_Zinc2jet;
    GenH1D *genPHI_LowPt_Zinc2jet;
    TH1D *BestPHI_LowPt_Zinc2jet;
    GenH1D *genBestPHI_LowPt_Zinc2jet;
    TH1D *PHI_T_LowPt_Zinc2jet;
    GenH1D *genPHI_T_LowPt_Zinc2jet;
    TH1D *BestPHI_T_LowPt_Zinc2jet;
    GenH1D *genBestPHI_T_LowPt_Zinc2jet;
    TH1D *SpT_LowPt_Zinc2jet;
    GenH1D *genSpT_LowPt_Zinc2jet;
    TH1D *BestSpT_LowPt_Zinc2jet;
    GenH1D *genBestSpT_LowPt_Zinc2jet;
    TH1D *SpTJets_LowPt_Zinc2jet;
    GenH1D *genSpTJets_LowPt_Zinc2jet;
    TH1D *BestSpTJets_LowPt_Zinc2jet;
    GenH1D *genBestSpTJets_LowPt_Zinc2jet;
    TH1D *SpTLeptons_LowPt_Zinc2jet;
    GenH1D *genSpTLeptons_LowPt_Zinc2jet;
    TH1D *SPhi_LowPt_Zinc2jet;
    GenH1D *genSPhi_LowPt_Zinc2jet;
    TH1D *BestSPhi_LowPt_Zinc2jet;
    GenH1D *genBestSPhi_LowPt_Zinc2jet;

    //-- low Z pT and low SpT;
    TH1D *PHI_LowSpT_LowPt_Zexc2jet;
    GenH1D *genPHI_LowSpT_LowPt_Zexc2jet;
    TH1D *SPhi_LowSpT_LowPt_Zexc2jet;
    GenH1D *genSPhi_LowSpT_LowPt_Zexc2jet;

    TH1D *PHI_LowSpT_LowPt_Zinc2jet;
    GenH1D *genPHI_LowSpT_LowPt_Zinc2jet;
    TH1D *SPhi_LowSpT_LowPt_Zinc2jet;
    GenH1D *genSPhi_LowSpT_LowPt_Zinc2jet;

    //-- low Z pT and high SpT;
    TH1D *PHI_HighSpT_LowPt_Zexc2jet;
    GenH1D *genPHI_HighSpT_LowPt_Zexc2jet;
    TH1D *SPhi_HighSpT_LowPt_Zexc2jet;
    GenH1D *genSPhi_HighSpT_LowPt_Zexc2jet;

    TH1D *PHI_HighSpT_LowPt_Zinc2jet;
    GenH1D *genPHI_HighSpT_LowPt_Zinc2jet;
    TH1D *SPhi_HighSpT_LowPt_Zinc2jet;
    GenH1D *genSPhi_HighSpT_LowPt_Zinc2jet;

    //-- low Z pT and low SPhi;
    TH1D *SpT_LowSPhi_LowPt_Zexc2jet;
    GenH1D *genSpT_LowSPhi_LowPt_Zexc2jet;

    TH1D *SpT_LowSPhi_LowPt_Zinc2jet;
    GenH1D *genSpT_LowSPhi_LowPt_Zinc2jet;

    //-- low Z pT and high SPhi;
    TH1D *SpT_HighSPhi_LowPt_Zexc2jet;
    GenH1D *genSpT_HighSPhi_LowPt_Zexc2jet;
    ;
    TH1D *SpT_HighSPhi_LowPt_Zinc2jet;
    GenH1D *genSpT_HighSPhi_LowPt_Zinc2jet;

    //-- high Z pT;
    TH1D *ptBal_HighPt_Zexc2jet;
    GenH1D *genptBal_HighPt_Zexc2jet;
    TH1D *dPhiJets_HighPt_Zexc2jet;
    GenH1D *gendPhiJets_HighPt_Zexc2jet;
    TH1D *dPhiLeptons_HighPt_Zexc2jet;
    GenH1D *gendPhiLeptons_HighPt_Zexc2jet;
    TH1D *PHI_HighPt_Zexc2jet;
    GenH1D *genPHI_HighPt_Zexc2jet;
    TH1D *PHI_T_HighPt_Zexc2jet;
    GenH1D *genPHI_T_HighPt_Zexc2jet;
    TH1D *SpT_HighPt_Zexc2jet;
    GenH1D *genSpT_HighPt_Zexc2jet;
    TH1D *SpTJets_HighPt_Zexc2jet;
    GenH1D *genSpTJets_HighPt_Zexc2jet;
    TH1D *SpTLeptons_HighPt_Zexc2jet;
    GenH1D *genSpTLeptons_HighPt_Zexc2jet;
    TH1D *SPhi_HighPt_Zexc2jet;
    GenH1D *genSPhi_HighPt_Zexc2jet;

    TH1D *ptBal_HighPt_Zinc2jet;
    GenH1D *genptBal_HighPt_Zinc2jet;
    TH1D *dPhiJets_HighPt_Zinc2jet;
    GenH1D *gendPhiJets_HighPt_Zinc2jet;
    TH1D *dPhiLeptons_HighPt_Zinc2jet;
    GenH1D *gendPhiLeptons_HighPt_Zinc2jet;
    TH1D *PHI_HighPt_Zinc2jet;
    GenH1D *genPHI_HighPt_Zinc2jet;
    TH1D *PHI_T_HighPt_Zinc2jet;
    GenH1D *genPHI_T_HighPt_Zinc2jet;
    TH1D *SpT_HighPt_Zinc2jet;
    GenH1D *genSpT_HighPt_Zinc2jet;
    TH1D *SpTJets_HighPt_Zinc2jet;
    GenH1D *genSpTJets_HighPt_Zinc2jet;
    TH1D *SpTLeptons_HighPt_Zinc2jet;
    GenH1D *genSpTLeptons_HighPt_Zinc2jet;
    TH1D *SPhi_HighPt_Zinc2jet;
    GenH1D *genSPhi_HighPt_Zinc2jet;

    //-- high Z pT and low SpT
    TH1D *PHI_LowSpT_HighPt_Zexc2jet;
    TH1D *SPhi_LowSpT_HighPt_Zexc2jet;

    TH1D *PHI_LowSpT_HighPt_Zinc2jet;
    TH1D *SPhi_LowSpT_HighPt_Zinc2jet;

    //-- high Z pT and high SpT
    TH1D *PHI_HighSpT_HighPt_Zexc2jet;
    TH1D *SPhi_HighSpT_HighPt_Zexc2jet;

    TH1D *PHI_HighSpT_HighPt_Zinc2jet;
    TH1D *SPhi_HighSpT_HighPt_Zinc2jet;

    //-- high Z pT and low SPhi
    TH1D *SpT_LowSPhi_HighPt_Zexc2jet;
    TH1D *SpT_LowSPhi_HighPt_Zinc2jet;

    //-- high Z pT and high SPhi
    TH1D *SpT_HighSPhi_HighPt_Zexc2jet;
    TH1D *SpT_HighSPhi_HighPt_Zinc2jet;

    //-- low SPhi
    TH1D *SpT_LowSPhi_Zexc2jet;
    TH1D *SpT_LowSPhi_Zinc2jet;

    //-- high SPhi
    TH1D *SpT_HighSPhi_Zexc2jet;
    TH1D *SpT_HighSPhi_Zinc2jet;

    //-- low SpT
    TH1D *PHI_LowSpT_Zexc2jet;
    TH1D *SPhi_LowSpT_Zexc2jet;

    TH1D *PHI_LowSpT_Zinc2jet;
    TH1D *SPhi_LowSpT_Zinc2jet;

    //-- high SpT
    TH1D *PHI_HighSpT_Zexc2jet;
    TH1D *SPhi_HighSpT_Zexc2jet;

    TH1D *PHI_HighSpT_Zinc2jet;
    TH1D *SPhi_HighSpT_Zinc2jet;

    //-- gen stuff
    GenH1D *gendPhiJetsDeltaR_Zexc2jet;
    TH1D *resdPhiJetsDeltaR_Zexc2jet;
    GenH1D *genPHI_TDeltaR_Zexc2jet;
    TH1D *resPHI_TDeltaR_Zexc2jet;
    GenH1D *genSpTJetsDeltaR_Zexc2jet;
    TH1D *resSpTJetsDeltaR_Zexc2jet;
    GenH1D *genSpTDeltaR_Zexc2jet;
    TH1D *resSpTDeltaR_Zexc2jet;

    GenH1D *gendPhiJetsDPS_Zexc2jet;
    GenH1D *gendPhiJetsDPSDeltaR_Zexc2jet;
    GenH1D *genPHI_TDPS_Zexc2jet;
    GenH1D *genPHI_TDPSDeltaR_Zexc2jet;
    GenH1D *genSpTJetsDPS_Zexc2jet;
    GenH1D *genSpTJetsDPSDeltaR_Zexc2jet;
    GenH1D *genSpTDPS_Zexc2jet;
    GenH1D *genSpTDPSDeltaR_Zexc2jet;
    GenH1D *genSpTDPSPartons_Zexc2jet;
    GenH1D *genZNGoodJets_Zinc;
    GenH1D *genZNGoodJets_Zexc;
    TH2D *gendPhiJetsDPSDeltaR_ZpT_Zexc2jet;
    TH2D *partonX2D;

    GenH1D *gendeltaRjetMu;

    /// additional information
    // Muoisolation

    TH1D *MuDetIsoRhoCorr;
    TH1D *MuPFIsoDBetaCorr;
    TH1D *MuPFIsoDBetaCorrj1;
    TH1D *MuPFIsoDBetaCorrj2;
    TH1D *MuPFIsoDBetaCorrj3;

    TH1D *deltaRMuRecGen_lead_lowM;
    TH1D *deltaRMuRecGen_sublead_lowM;
    TH1D *deltaRMuRecGen_lead_highM;
    TH1D *deltaRMuRecGen_sublead_highM;
    TH1D *deltaPtMuRecGen_lead_highM;
    TH1D *deltaPtMuRecGen_sublead_highM;
    TH1D *deltaPtMuRecGen_lead_lowM;
    TH1D *deltaPtMuRecGen_sublead_lowM;

    TH1D *deltaRjetMu;
    TH1D *deltaPtjetMu;
    TH1D *TruePU_0;
    TH1D *TruePU_1;
    TH1D *TruePU_2;
    TH1D *TruePU_3;
    TH1D *TruePU_4;
    TH1D *TruePU_5;
    TH1D *TruePU_6;
    TH1D *TruePU_7;

    TH1D *NVtx_Zinc0jet;
    TH1D *NVtx_NoPUweight_Zinc0jet;
    TH1D *NVtx_Zexc0jet;
    TH1D *NVtx_Zexc1jet;
    TH1D *NVtx_Zexc2jet;
    TH1D *NVtx_Zexc3jet;
    TH1D *NVtx_Zexc4jet;
    TH1D *NVtx_Zexc5jet;
    TH1D *NVtx_Zexc6jet;
    TH1D *NVtx_Zexc7jet;

    TH1D *Beta;
    TH1D *BetaStar;
    TH2D *ZNGoodJetsBeta_Zexc;
    TH1D *puBeta_JetsMatchGenJets;
    TH1D *puBetaStar_JetsMatchGenJets;
    TH1D *puBeta_JetsNoMatchGenJets;
    TH1D *puBetaStar_JetsNoMatchGenJets;
    TH1D *puMVA;
    TH1D *puMVA_JetsMatchGenJets;
    TH1D *puMVA_JetsNoMatchGenJets;
    TH1D *jetsEta_JetsMatchGenJets;
    TH1D *jetsEta_JetsNoMatchGenJets;
    TH2D *puMVAvsBeta;
    TH1D *PUWeight;
    TH1D *PUWeight0;
    TH1D *PUWeight1;

    //--- number of addtion genrated parton
    TH1D *partonsN;
    TH1D *partonsNWeighted;
    TH1D *partonsNAfterGenCut;
    TH1D *partonsNAfterGenCutWeighted;

    // -- Vector boson jet properties
    TH1D *dEtaBosonJet_Zexc1jet;
    TH1D *dEtaBosonJet_Zinc1jet;
    GenH1D *gendEtaBosonJet_Zexc1jet;
    GenH1D *gendEtaBosonJet_Zinc1jet;

    TH1D *trMass;
    TH1D *trE;
    TH1D *E;
    TH1D *eta;
    TH1D *rapidity;

    TH1D *tau_sum_Zinc1jet[5];
    TH1D *tau_max_Zinc1jet[5];
    TH1D *tau_c_sum_Zinc1jet[5];
    TH1D *tau_c_max_Zinc1jet[5];
    TH1D *tau_cm_sum_Zinc1jet[5];
    TH1D *tau_cm_max_Zinc1jet[5];
    TH1D *tau_c_cm_sum_Zinc1jet[5];
    TH1D *tau_c_cm_max_Zinc1jet[5];

    GenH1D *gentau_sum_Zinc1jet[5];
    GenH1D *gentau_max_Zinc1jet[5];
    GenH1D *gentau_c_sum_Zinc1jet[5];
    GenH1D *gentau_c_max_Zinc1jet[5];
    GenH1D *gentau_cm_sum_Zinc1jet[5];
    GenH1D *gentau_cm_max_Zinc1jet[5];
    GenH1D *gentau_c_cm_sum_Zinc1jet[5];
    GenH1D *gentau_c_cm_max_Zinc1jet[5];

    // Additional Histograms
    TH1D *AbsZRapidity_Zinc1jet;
    GenH1D *genAbsZRapidity_Zinc1jet;
    TH1D *AbsFirstJetRapidity_Zinc1jet;
    GenH1D *genAbsFirstJetRapidity_Zinc1jet;
    TH1D *SumZFirstJetRapidity_Zinc1jet;
    GenH1D *genSumZFirstJetRapidity_Zinc1jet;
    TH1D *DifZFirstJetRapidity_Zinc1jet;
    GenH1D *genDifZFirstJetRapidity_Zinc1jet;

    TH1D *SumZFirstJetEta_Zinc1jet;
    GenH1D *genSumZFirstJetEta_Zinc1jet;
    TH1D *DifZFirstJetEta_Zinc1jet;
    GenH1D *genDifZFirstJetEta_Zinc1jet;

    TH1D *DPhiZFirstJet_Zinc1jet;
    GenH1D *genDPhiZFirstJet_Zinc1jet;

    TH1D *AbsZRapidity_Zexc1jet;
    GenH1D *genAbsZRapidity_Zexc1jet;
    TH1D *AbsJetRapidity_Zexc1jet;
    GenH1D *genAbsJetRapidity_Zexc1jet;
    TH1D *SumZJetRapidity_Zexc1jet;
    GenH1D *genSumZJetRapidity_Zexc1jet;
    TH1D *DifZJetRapidity_Zexc1jet;
    GenH1D *genDifZJetRapidity_Zexc1jet;
    TH1D *AbsFirstJetRapidity_Zinc2jet;
    GenH1D *genAbsFirstJetRapidity_Zinc2jet;
    TH1D *SumZFirstJetRapidity_Zinc2jet;
    GenH1D *genSumZFirstJetRapidity_Zinc2jet;
    TH1D *DifZFirstJetRapidity_Zinc2jet;
    GenH1D *genDifZFirstJetRapidity_Zinc2jet;

    TH1D *SumZTwoJetsRapidity_Zinc2jet;
    GenH1D *genSumZTwoJetsRapidity_Zinc2jet;
    TH1D *DifZTwoJetsRapidity_Zinc2jet;
    GenH1D *genDifZTwoJetsRapidity_Zinc2jet;

    TH1D *AbsZRapidity_Zinc2jet;
    GenH1D *genAbsZRapidity_Zinc2jet;
    TH1D *AbsSecondJetRapidity_Zinc2jet;
    GenH1D *genAbsSecondJetRapidity_Zinc2jet;
    TH1D *SumZSecondJetRapidity_Zinc2jet;
    GenH1D *genSumZSecondJetRapidity_Zinc2jet;
    TH1D *DifZSecondJetRapidity_Zinc2jet;
    GenH1D *genDifZSecondJetRapidity_Zinc2jet;

    TH1D *SumFirstSecondJetRapidity_Zinc2jet;
    GenH1D *genSumFirstSecondJetRapidity_Zinc2jet;
    TH1D *DifFirstSecondJetRapidity_Zinc2jet;
    GenH1D *genDifFirstSecondJetRapidity_Zinc2jet;

    TH1D *DPhiZFirstJet_Zinc2jet;
    GenH1D *genDPhiZFirstJet_Zinc2jet;
    TH1D *DPhiZSecondJet_Zinc2jet;
    GenH1D *genDPhiZSecondJet_Zinc2jet;
    TH1D *DPhiFirstSecondJet_Zinc2jet;
    GenH1D *genDPhiFirstSecondJet_Zinc2jet;

    TH1D *AbsZRapidity_Zexc2jet;
    GenH1D *genAbsZRapidity_Zexc2jet;
    TH1D *AbsSecondJetRapidity_Zexc2jet;
    GenH1D *genAbsSecondJetRapidity_Zexc2jet;
    TH1D *SumZSecondJetRapidity_Zexc2jet;
    GenH1D *genSumZSecondJetRapidity_Zexc2jet;
    TH1D *DifZSecondJetRapidity_Zexc2jet;
    GenH1D *genDifZSecondJetRapidity_Zexc2jet;

    TH1D *DPhiZFirstJet_Zinc3jet;
    GenH1D *genDPhiZFirstJet_Zinc3jet;
    TH1D *DPhiZSecondJet_Zinc3jet;
    GenH1D *genDPhiZSecondJet_Zinc3jet;
    TH1D *DPhiZThirdJet_Zinc3jet;
    GenH1D *genDPhiZThirdJet_Zinc3jet;
    TH1D *DPhiFirstSecondJet_Zinc3jet;
    GenH1D *genDPhiFirstSecondJet_Zinc3jet;
    TH1D *DPhiFirstThirdJet_Zinc3jet;
    GenH1D *genDPhiFirstThirdJet_Zinc3jet;
    TH1D *DPhiSecondThirdJet_Zinc3jet;
    GenH1D *genDPhiSecondThirdJet_Zinc3jet;

    ////Different Z boson Pt cut branches///////////////////
    TH1D *AbsZRapidity_ZPt100_Zinc1jet;
    GenH1D *genAbsZRapidity_ZPt100_Zinc1jet;
    TH1D *AbsFirstJetRapidity_ZPt100_Zinc1jet;
    GenH1D *genAbsFirstJetRapidity_ZPt100_Zinc1jet;
    TH1D *SumZFirstJetRapidity_ZPt100_Zinc1jet;
    GenH1D *genSumZFirstJetRapidity_ZPt100_Zinc1jet;
    TH1D *DifZFirstJetRapidity_ZPt100_Zinc1jet;
    GenH1D *genDifZFirstJetRapidity_ZPt100_Zinc1jet;

    TH1D *AbsZRapidity_ZPt100_Zexc1jet;
    GenH1D *genAbsZRapidity_ZPt100_Zexc1jet;
    TH1D *AbsJetRapidity_ZPt100_Zexc1jet;
    GenH1D *genAbsJetRapidity_ZPt100_Zexc1jet;
    TH1D *SumZJetRapidity_ZPt100_Zexc1jet;
    GenH1D *genSumZJetRapidity_ZPt100_Zexc1jet;
    TH1D *DifZJetRapidity_ZPt100_Zexc1jet;
    GenH1D *genDifZJetRapidity_ZPt100_Zexc1jet;

    TH1D *AbsZRapidity_ZPt100_Zinc2jet;
    GenH1D *genAbsZRapidity_ZPt100_Zinc2jet;
    TH1D *AbsSecondJetRapidity_ZPt100_Zinc2jet;
    GenH1D *genAbsSecondJetRapidity_ZPt100_Zinc2jet;
    TH1D *SumZSecondJetRapidity_ZPt100_Zinc2jet;
    GenH1D *genSumZSecondJetRapidity_ZPt100_Zinc2jet;
    TH1D *DifZSecondJetRapidity_ZPt100_Zinc2jet;
    GenH1D *genDifZSecondJetRapidity_ZPt100_Zinc2jet;

    TH1D *AbsZRapidity_ZPt100_Zexc2jet;
    GenH1D *genAbsZRapidity_ZPt100_Zexc2jet;
    TH1D *AbsSecondJetRapidity_ZPt100_Zexc2jet;
    GenH1D *genAbsSecondJetRapidity_ZPt100_Zexc2jet;
    TH1D *SumZSecondJetRapidity_ZPt100_Zexc2jet;
    GenH1D *genSumZSecondJetRapidity_ZPt100_Zexc2jet;
    TH1D *DifZSecondJetRapidity_ZPt100_Zexc2jet;
    GenH1D *genDifZSecondJetRapidity_ZPt100_Zexc2jet;

    TH1D *AbsZRapidity_ZPt150_Zinc1jet;
    GenH1D *genAbsZRapidity_ZPt150_Zinc1jet;
    TH1D *AbsFirstJetRapidity_ZPt150_Zinc1jet;
    GenH1D *genAbsFirstJetRapidity_ZPt150_Zinc1jet;
    TH1D *SumZFirstJetRapidity_ZPt150_Zinc1jet;
    GenH1D *genSumZFirstJetRapidity_ZPt150_Zinc1jet;
    TH1D *DifZFirstJetRapidity_ZPt150_Zinc1jet;
    GenH1D *genDifZFirstJetRapidity_ZPt150_Zinc1jet;

    TH1D *AbsZRapidity_ZPt150_Zexc1jet;
    GenH1D *genAbsZRapidity_ZPt150_Zexc1jet;
    TH1D *AbsJetRapidity_ZPt150_Zexc1jet;
    GenH1D *genAbsJetRapidity_ZPt150_Zexc1jet;
    TH1D *SumZJetRapidity_ZPt150_Zexc1jet;
    GenH1D *genSumZJetRapidity_ZPt150_Zexc1jet;
    TH1D *DifZJetRapidity_ZPt150_Zexc1jet;
    GenH1D *genDifZJetRapidity_ZPt150_Zexc1jet;

    TH1D *AbsZRapidity_ZPt150_Zinc2jet;
    GenH1D *genAbsZRapidity_ZPt150_Zinc2jet;
    TH1D *AbsSecondJetRapidity_ZPt150_Zinc2jet;
    GenH1D *genAbsSecondJetRapidity_ZPt150_Zinc2jet;
    TH1D *SumZSecondJetRapidity_ZPt150_Zinc2jet;
    GenH1D *genSumZSecondJetRapidity_ZPt150_Zinc2jet;
    TH1D *DifZSecondJetRapidity_ZPt150_Zinc2jet;
    GenH1D *genDifZSecondJetRapidity_ZPt150_Zinc2jet;

    TH1D *AbsZRapidity_ZPt150_Zexc2jet;
    GenH1D *genAbsZRapidity_ZPt150_Zexc2jet;
    TH1D *AbsSecondJetRapidity_ZPt150_Zexc2jet;
    GenH1D *genAbsSecondJetRapidity_ZPt150_Zexc2jet;
    TH1D *SumZSecondJetRapidity_ZPt150_Zexc2jet;
    GenH1D *genSumZSecondJetRapidity_ZPt150_Zexc2jet;
    TH1D *DifZSecondJetRapidity_ZPt150_Zexc2jet;
    GenH1D *genDifZSecondJetRapidity_ZPt150_Zexc2jet;

    TH1D *AbsZRapidity_ZPt300_Zinc1jet;
    GenH1D *genAbsZRapidity_ZPt300_Zinc1jet;
    TH1D *AbsFirstJetRapidity_ZPt300_Zinc1jet;
    GenH1D *genAbsFirstJetRapidity_ZPt300_Zinc1jet;
    TH1D *SumZFirstJetRapidity_ZPt300_Zinc1jet;
    GenH1D *genSumZFirstJetRapidity_ZPt300_Zinc1jet;
    TH1D *DifZFirstJetRapidity_ZPt300_Zinc1jet;
    GenH1D *genDifZFirstJetRapidity_ZPt300_Zinc1jet;

    ////Azimuthal cross check
    TH1D *DPhiZFirstJet_ZPt150_Zinc1jet;
    GenH1D *genDPhiZFirstJet_ZPt150_Zinc1jet;
    TH1D *DPhiZFirstJet_ZPt150_Zinc2jet;
    GenH1D *genDPhiZFirstJet_ZPt150_Zinc2jet;
    TH1D *DPhiZFirstJet_ZPt150_Zinc3jet;
    GenH1D *genDPhiZFirstJet_ZPt150_Zinc3jet;
    TH1D *DPhiZSecondJet_ZPt150_Zinc3jet;
    GenH1D *genDPhiZSecondJet_ZPt150_Zinc3jet;
    TH1D *DPhiZThirdJet_ZPt150_Zinc3jet;
    GenH1D *genDPhiZThirdJet_ZPt150_Zinc3jet;

    TH1D *DPhiFirstSecondJet_ZPt150_Zinc3jet;
    GenH1D *genDPhiFirstSecondJet_ZPt150_Zinc3jet;
    TH1D *DPhiFirstThirdJet_ZPt150_Zinc3jet;
    GenH1D *genDPhiFirstThirdJet_ZPt150_Zinc3jet;
    TH1D *DPhiSecondThirdJet_ZPt150_Zinc3jet;
    GenH1D *genDPhiSecondThirdJet_ZPt150_Zinc3jet;

    TH1D *DPhiZFirstJet_ZPt300_Zinc1jet;
    GenH1D *genDPhiZFirstJet_ZPt300_Zinc1jet;
    TH1D *DPhiZFirstJet_ZPt300_Zinc2jet;
    GenH1D *genDPhiZFirstJet_ZPt300_Zinc2jet;
    TH1D *DPhiZFirstJet_ZPt300_Zinc3jet;
    GenH1D *genDPhiZFirstJet_ZPt300_Zinc3jet;
    TH1D *DPhiZSecondJet_ZPt300_Zinc3jet;
    GenH1D *genDPhiZSecondJet_ZPt300_Zinc3jet;
    TH1D *DPhiZThirdJet_ZPt300_Zinc3jet;
    GenH1D *genDPhiZThirdJet_ZPt300_Zinc3jet;

    TH1D *DPhiFirstSecondJet_ZPt300_Zinc3jet;
    GenH1D *genDPhiFirstSecondJet_ZPt300_Zinc3jet;
    TH1D *DPhiFirstThirdJet_ZPt300_Zinc3jet;
    GenH1D *genDPhiFirstThirdJet_ZPt300_Zinc3jet;
    TH1D *DPhiSecondThirdJet_ZPt300_Zinc3jet;
    GenH1D *genDPhiSecondThirdJet_ZPt300_Zinc3jet;

    TH1D *DPhiZFirstJet_ZPt150_HT300_Zinc3jet;
    GenH1D *genDPhiZFirstJet_ZPt150_HT300_Zinc3jet;
    TH1D *DPhiZSecondJet_ZPt150_HT300_Zinc3jet;
    GenH1D *genDPhiZSecondJet_ZPt150_HT300_Zinc3jet;
    TH1D *DPhiZThirdJet_ZPt150_HT300_Zinc3jet;
    GenH1D *genDPhiZThirdJet_ZPt150_HT300_Zinc3jet;

    /// Branches with different JetPt cut///////////

    TH1D *AbsZRapidity_FirstJetPt50_Zinc1jet;
    GenH1D *genAbsZRapidity_FirstJetPt50_Zinc1jet;
    TH1D *AbsFirstJetRapidity_FirstJetPt50_Zinc1jet;
    GenH1D *genAbsFirstJetRapidity_FirstJetPt50_Zinc1jet;
    TH1D *SumZFirstJetRapidity_FirstJetPt50_Zinc1jet;
    GenH1D *genSumZFirstJetRapidity_FirstJetPt50_Zinc1jet;
    TH1D *DifZFirstJetRapidity_FirstJetPt50_Zinc1jet;
    GenH1D *genDifZFirstJetRapidity_FirstJetPt50_Zinc1jet;

    TH1D *AbsZRapidity_FirstJetPt80_Zinc1jet;
    GenH1D *genAbsZRapidity_FirstJetPt80_Zinc1jet;
    TH1D *AbsFirstJetRapidity_FirstJetPt80_Zinc1jet;
    GenH1D *genAbsFirstJetRapidity_FirstJetPt80_Zinc1jet;
    TH1D *SumZFirstJetRapidity_FirstJetPt80_Zinc1jet;
    GenH1D *genSumZFirstJetRapidity_FirstJetPt80_Zinc1jet;
    TH1D *DifZFirstJetRapidity_FirstJetPt80_Zinc1jet;
    GenH1D *genDifZFirstJetRapidity_FirstJetPt80_Zinc1jet;

    // Set Jet rapidity difference discriminator/////
    TH1D *AbsZRapidity_DifJetRapidityl2_Zinc2jet;
    GenH1D *genAbsZRapidity_DifJetRapidityl2_Zinc2jet;
    TH1D *AbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet;
    GenH1D *genAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet;
    TH1D *SumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet;
    GenH1D *genSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet;
    TH1D *DifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet;
    GenH1D *genDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet;

    TH1D *AbsZRapidity_DifJetRapiditys2_Zinc2jet;
    GenH1D *genAbsZRapidity_DifJetRapiditys2_Zinc2jet;
    TH1D *AbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet;
    GenH1D *genAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet;
    TH1D *SumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet;
    GenH1D *genSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet;
    TH1D *DifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet;
    GenH1D *genDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet;
    TH1D *SumZJetRapidity_Zinc1jet;
    GenH1D *genSumZJetRapidity_Zinc1jet;
    TH1D *DifZJetRapidity_Zinc1jet;
    GenH1D *genDifZJetRapidity_Zinc1jet;

    TH1D *CentralJetPt_Zinc2jet;
    TH1D *ForwardJetPt_Zinc2jet;
    TH1D *CentralJetEta_Zinc2jet;
    TH1D *ForwardJetEta_Zinc2jet;
    GenH1D *genCentralJetPt_Zinc2jet;
    GenH1D *genForwardJetPt_Zinc2jet;
    GenH1D *genCentralJetEta_Zinc2jet;
    GenH1D *genForwardJetEta_Zinc2jet;

    // DJALOG TH2D
    TH1D *hjetResolution;
    TH1D *hEvtFastJetRho;
    TH2D *hhjetResolutionPt;
    TH2D *hhjetResolutionEta;
    TH1D *hjetScaleFactor;
    TH2D *hhJetMatching_StatMatches;
    TH2D *hhJetMatching_StatFail;
    TH2D *hhJetMatching_StatMatches_Lep;
    TH2D *hhJetMatching_StatFail_Lep;
    TH2D *hhRecoGenJetMatchingMatrix;
    TH2D *hhRecoGenJetMatchingMatrix_Lep;
    TH2D *hhJetMatching_FirstJetRap_Match;
    TH2D *hhJetMatching_FirstJetRap_NoMatch;
    TH2D *hhJetMatching_FirstJetNoMatch_Matrix;
    TH1D *hJetMatching_FirstJetdPT_Match;
    TH1D *hJetMatching_FirstJetdPT_Other;
    TH1D *hJetMatching_dPT_SecondMatch;
    TH2D *hhRecoGenLeptonMatching;
    TH1D *hJetMatching_dRFirstJet;
    TH1D *hJetMathing_FailedRecoJet_AllEta;
    TH1D *hJetMathing_FailedRecoJet_FirstEta;
    TH1D *hJetMatching_dRFailedFirstLeadLep;
    TH1D *hJetMatching_dRFailedFirstSubLep;
    TH1D *hJetMatching_FirstJetPt_FaileddR;
    TH1D *hJetMatching_FirstJetAbsRapidity_FaileddR;
    TH1D *hJetMatching_FirstJetPt_FailednGen;
    TH1D *hJetMatching_FirstJetAbsRapidity_FailednGen;

    // TH2D for reponse objects

    TH2D *hresponsetau_sum_Zinc1jet[5];
    TH2D *hresponsetau_max_Zinc1jet[5];
    TH2D *hresponsetau_c_sum_Zinc1jet[5];
    TH2D *hresponsetau_c_max_Zinc1jet[5];
    TH2D *hresponsetau_cm_sum_Zinc1jet[5];
    TH2D *hresponsetau_cm_max_Zinc1jet[5];
    TH2D *hresponsetau_c_cm_sum_Zinc1jet[5];
    TH2D *hresponsetau_c_cm_max_Zinc1jet[5];

    TH2D *hresponseZNGoodJets_Zexc;

    TH2D *hresponseHadRecoil;
    TH2D *hresponseJZB;
    TH2D *hresponseJZB_ptLow;
    TH2D *hresponseJZB_ptHigh;
    TH2D *hresponseZPt_Zinc2jet;
    TH2D *hresponseVisPt_Zinc0jetQun;
    TH2D *hresponseVisPt_Zinc1jetQun;
    TH2D *hresponseVisPt_Zinc2jetQun;
    TH2D *hresponseVisPt_Zinc3jetQun;
    TH2D *hresponseZAbsRapidity_Zinc1jet;
    TH2D *hresponseSumZJetRapidity_Zinc1jet;
    TH2D *hresponseDifZJetRapidity_Zinc1jet;

    TH2D *hresponseFirstJetPt_Zinc1jet;
    TH2D *hresponseFirstJetPtMatch_Zinc1jet;
    TH2D *hresponseFirstJetPt_2_Zinc1jet;
    TH2D *hresponseSecondJetPt_Zinc2jet;
    TH2D *hresponseSecondJetPtMatch_Zinc2jet; // DJALOG
    TH2D *hresponseSecondJetPt_2_Zinc2jet;
    TH2D *hresponseThirdJetPt_Zinc3jet;
    TH2D *hresponseThirdJetPtMatch_Zinc3jet; // DJALOG
    TH2D *hresponseThirdJetPt_2_Zinc3jet;
    TH2D *hresponseFourthJetPt_Zinc4jet;
    TH2D *hresponseFourthJetPtMatch_Zinc4jet; // DJALOG
    TH2D *hresponseFifthJetPt_Zinc5jet;

    TH2D *hresponseJetsHT_Zinc1jet;
    TH2D *hresponseJetsHT_2_Zinc1jet;
    TH2D *hresponseJetsHT_Zinc2jet;
    TH2D *hresponseJetsHT_2_Zinc2jet;
    TH2D *hresponseJetsHT_Zinc3jet;
    TH2D *hresponseJetsHT_2_Zinc3jet;
    TH2D *hresponseJetsHT_Zinc4jet;
    TH2D *hresponseJetsHT_Zinc5jet;

    TH2D *hresponseFirstJetEta_Zinc1jet;
    TH2D *hresponseFirstJetEtaMatch_Zinc1jet; // DJALOG
    TH2D *hresponseFirstJetEta_2_Zinc1jet;

    TH2D *hresponseSecondJetEta_Zinc2jet;
    TH2D *hresponseSecondJetEtaMatch_Zinc2jet; // DJALOG
    TH2D *hresponseSecondJetEta_2_Zinc2jet;
    TH2D *hresponseThirdJetEta_Zinc3jet;
    TH2D *hresponseThirdJetEtaMatch_Zinc3jet; // DJALOG
    TH2D *hresponseThirdJetEta_2_Zinc3jet;
    TH2D *hresponseFourthJetEta_Zinc4jet;
    TH2D *hresponseFourthJetEtaMatch_Zinc4jet; // DJALOG
    TH2D *hresponseFifthJetEta_Zinc5jet;

    TH2D *hresponseFirstJetEtaHigh_Zinc1jet;
    TH2D *hresponseFirstJetEtaHighMatch_Zinc1jet; // DJALOG
    TH2D *hresponseSecondJetEtaHigh_Zinc2jet;
    TH2D *hresponseSecondJetEtaHighMatch_Zinc2jet; // DJALOG
    TH2D *hresponseThirdJetEtaHigh_Zinc3jet;
    TH2D *hresponseThirdJetEtaHighMatch_Zinc3jet; // DJALOG
    TH2D *hresponseFourthJetEtaHigh_Zinc4jet;
    TH2D *hresponseFourthJetEtaHighMatch_Zinc4jet; // DJALOG
    TH2D *hresponseFifthJetEtaHigh_Zinc5jet;

    TH2D *hresponseFirstJetRapidityHigh_Zinc1jet;
    TH2D *hresponseFirstJetRapidityHighMatch_Zinc1jet; // DJALOG
    TH2D *hresponseFirstJetAbsRapidity_Zinc1jet;
    TH2D *hresponseFirstJetAbsRapidityMatch_Zinc1jet; // DJALOG
    TH2D *hresponseSecondJetRapidityHigh_Zinc2jet;
    TH2D *hresponseSecondJetAbsRapidity_Zinc2jet;
    TH2D *hresponseSecondJetRapidityHighMatch_Zinc2jet; // DJALOG
    TH2D *hresponseSecondJetAbsRapidityMatch_Zinc2jet;  // DJALOG
    TH2D *hresponseThirdJetRapidityHigh_Zinc3jet;
    TH2D *hresponseThirdJetAbsRapidity_Zinc3jet;
    TH2D *hresponseThirdJetRapidityHighMatch_Zinc3jet; // DJALOG
    TH2D *hresponseThirdJetAbsRapidityMatch_Zinc3jet;  // DJALOG
    TH2D *hresponseFourthJetRapidityHigh_Zinc4jet;
    TH2D *hresponseFourthJetAbsRapidity_Zinc4jet;
    TH2D *hresponseFourthJetRapidityHighMatch_Zinc4jet; // DJALOG
    TH2D *hresponseFourthJetAbsRapidityMatch_Zinc4jet;  // DJALOG
    TH2D *hresponseFifthJetRapidityHigh_Zinc5jet;
    TH2D *hresponseFifthJetAbsRapidity_Zinc5jet;

    TH2D *hresponseJetsMass_Zinc2jet;
    TH2D *hresponseJetsMassLowPU_Zinc2jet;
    TH2D *hresponseJetsMassMidPU_Zinc2jet;
    TH2D *hresponseJetsMassHigPU_Zinc2jet;

    TH2D *hresponseFirstJetPtEta_Zinc1jet;
    TH2D *hresponseFirstJetPtEtaMatch_Zinc1jet; // DJALOG
    TH2D *FirstJetPtEta_Zinc1jet;
    TH2D *genFirstJetPtEta_Zinc1jet;

    // Additional response variables
    TH2D *hresponseAbsZRapidity_Zinc1jet;
    TH2D *hresponseAbsFirstJetRapidity_Zinc1jet;
    TH2D *hresponseSumZFirstJetRapidity_Zinc1jet;
    TH2D *hresponseDifZFirstJetRapidity_Zinc1jet;

    TH2D *hresponseSumZFirstJetEta_Zinc1jet;
    TH2D *hresponseDifZFirstJetEta_Zinc1jet;

    TH2D *hresponseDPhiZFirstJet_Zinc1jet;

    TH2D *hresponseAbsZRapidity_Zexc1jet;
    TH2D *hresponseAbsJetRapidity_Zexc1jet;
    TH2D *hresponseSumZJetRapidity_Zexc1jet;
    TH2D *hresponseDifZJetRapidity_Zexc1jet;

    TH2D *hresponseAbsFirstJetRapidity_Zinc2jet;
    TH2D *hresponseSumZFirstJetRapidity_Zinc2jet;
    TH2D *hresponseDifZFirstJetRapidity_Zinc2jet;

    TH2D *hresponseAbsZRapidity_Zinc2jet;
    TH2D *hresponseAbsSecondJetRapidity_Zinc2jet;
    TH2D *hresponseSumZSecondJetRapidity_Zinc2jet;
    TH2D *hresponseDifZSecondJetRapidity_Zinc2jet;

    TH2D *hresponseSumFirstSecondJetRapidity_Zinc2jet;
    TH2D *hresponseDifFirstSecondJetRapidity_Zinc2jet;

    TH2D *hresponseSumZTwoJetsRapidity_Zinc2jet;
    TH2D *hresponseDifZTwoJetsRapidity_Zinc2jet;

    TH2D *hresponseDPhiZFirstJet_Zinc2jet;
    TH2D *hresponseDPhiZSecondJet_Zinc2jet;
    TH2D *hresponseDPhiFirstSecondJet_Zinc2jet;

    TH2D *hresponseAbsZRapidity_Zexc2jet;
    TH2D *hresponseAbsSecondJetRapidity_Zexc2jet;
    TH2D *hresponseSumZSecondJetRapidity_Zexc2jet;
    TH2D *hresponseDifZSecondJetRapidity_Zexc2jet;

    TH2D *hresponseDPhiZFirstJet_Zinc3jet;
    TH2D *hresponseDPhiZSecondJet_Zinc3jet;
    TH2D *hresponseDPhiZThirdJet_Zinc3jet;
    TH2D *hresponseDPhiFirstSecondJet_Zinc3jet;
    TH2D *hresponseDPhiFirstThirdJet_Zinc3jet;
    TH2D *hresponseDPhiSecondThirdJet_Zinc3jet;

    /// Different Z boson Pt cut////////////////////
    TH2D *hresponseAbsZRapidity_ZPt100_Zinc1jet;
    TH2D *hresponseAbsFirstJetRapidity_ZPt100_Zinc1jet;
    TH2D *hresponseSumZFirstJetRapidity_ZPt100_Zinc1jet;
    TH2D *hresponseDifZFirstJetRapidity_ZPt100_Zinc1jet;

    TH2D *hresponseAbsZRapidity_ZPt100_Zexc1jet;
    TH2D *hresponseAbsJetRapidity_ZPt100_Zexc1jet;
    TH2D *hresponseSumZJetRapidity_ZPt100_Zexc1jet;
    TH2D *hresponseDifZJetRapidity_ZPt100_Zexc1jet;

    TH2D *hresponseAbsZRapidity_ZPt100_Zinc2jet;
    TH2D *hresponseAbsSecondJetRapidity_ZPt100_Zinc2jet;
    TH2D *hresponseSumZSecondJetRapidity_ZPt100_Zinc2jet;
    TH2D *hresponseDifZSecondJetRapidity_ZPt100_Zinc2jet;

    TH2D *hresponseAbsZRapidity_ZPt100_Zexc2jet;
    TH2D *hresponseAbsSecondJetRapidity_ZPt100_Zexc2jet;
    TH2D *hresponseSumZSecondJetRapidity_ZPt100_Zexc2jet;
    TH2D *hresponseDifZSecondJetRapidity_ZPt100_Zexc2jet;

    TH2D *hresponseAbsZRapidity_ZPt150_Zinc1jet;
    TH2D *hresponseAbsFirstJetRapidity_ZPt150_Zinc1jet;
    TH2D *hresponseSumZFirstJetRapidity_ZPt150_Zinc1jet;
    TH2D *hresponseDifZFirstJetRapidity_ZPt150_Zinc1jet;

    TH2D *hresponseAbsZRapidity_ZPt150_Zexc1jet;
    TH2D *hresponseAbsJetRapidity_ZPt150_Zexc1jet;
    TH2D *hresponseSumZJetRapidity_ZPt150_Zexc1jet;
    TH2D *hresponseDifZJetRapidity_ZPt150_Zexc1jet;

    TH2D *hresponseAbsZRapidity_ZPt150_Zinc2jet;
    TH2D *hresponseAbsSecondJetRapidity_ZPt150_Zinc2jet;
    TH2D *hresponseSumZSecondJetRapidity_ZPt150_Zinc2jet;
    TH2D *hresponseDifZSecondJetRapidity_ZPt150_Zinc2jet;

    TH2D *hresponseAbsZRapidity_ZPt150_Zexc2jet;
    TH2D *hresponseAbsSecondJetRapidity_ZPt150_Zexc2jet;
    TH2D *hresponseSumZSecondJetRapidity_ZPt150_Zexc2jet;
    TH2D *hresponseDifZSecondJetRapidity_ZPt150_Zexc2jet;

    TH2D *hresponseAbsZRapidity_ZPt300_Zinc1jet;
    TH2D *hresponseAbsFirstJetRapidity_ZPt300_Zinc1jet;
    TH2D *hresponseSumZFirstJetRapidity_ZPt300_Zinc1jet;
    TH2D *hresponseDifZFirstJetRapidity_ZPt300_Zinc1jet;

    /// Azimuthal cross check/////////////////
    TH2D *hresponseDPhiZFirstJet_ZPt150_Zinc1jet;
    TH2D *hresponseDPhiZFirstJet_ZPt150_Zinc2jet;
    TH2D *hresponseDPhiZFirstJet_ZPt150_Zinc3jet;
    TH2D *hresponseDPhiZSecondJet_ZPt150_Zinc3jet;
    TH2D *hresponseDPhiZThirdJet_ZPt150_Zinc3jet;
    TH2D *hresponseDPhiFirstSecondJet_ZPt150_Zinc3jet;
    TH2D *hresponseDPhiFirstThirdJet_ZPt150_Zinc3jet;
    TH2D *hresponseDPhiSecondThirdJet_ZPt150_Zinc3jet;

    TH2D *hresponseDPhiZFirstJet_ZPt300_Zinc1jet;
    TH2D *hresponseDPhiZFirstJet_ZPt300_Zinc2jet;
    TH2D *hresponseDPhiZFirstJet_ZPt300_Zinc3jet;
    TH2D *hresponseDPhiZSecondJet_ZPt300_Zinc3jet;
    TH2D *hresponseDPhiZThirdJet_ZPt300_Zinc3jet;
    TH2D *hresponseDPhiFirstSecondJet_ZPt300_Zinc3jet;
    TH2D *hresponseDPhiFirstThirdJet_ZPt300_Zinc3jet;
    TH2D *hresponseDPhiSecondThirdJet_ZPt300_Zinc3jet;

    TH2D *hresponseDPhiZFirstJet_ZPt150_HT300_Zinc3jet;
    TH2D *hresponseDPhiZSecondJet_ZPt150_HT300_Zinc3jet;
    TH2D *hresponseDPhiZThirdJet_ZPt150_HT300_Zinc3jet;

    /// Branches with different JetPt Cuts////////

    TH2D *hresponseAbsZRapidity_FirstJetPt50_Zinc1jet;
    TH2D *hresponseAbsFirstJetRapidity_FirstJetPt50_Zinc1jet;
    TH2D *hresponseSumZFirstJetRapidity_FirstJetPt50_Zinc1jet;
    TH2D *hresponseDifZFirstJetRapidity_FirstJetPt50_Zinc1jet;

    TH2D *hresponseAbsZRapidity_FirstJetPt80_Zinc1jet;
    TH2D *hresponseAbsFirstJetRapidity_FirstJetPt80_Zinc1jet;
    TH2D *hresponseSumZFirstJetRapidity_FirstJetPt80_Zinc1jet;
    TH2D *hresponseDifZFirstJetRapidity_FirstJetPt80_Zinc1jet;

    // Set Jet rapidity discriminator//////

    TH2D *hresponseAbsZRapidity_DifJetRapidityl2_Zinc2jet;
    TH2D *hresponseAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet;
    TH2D *hresponseSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet;
    TH2D *hresponseDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet;

    TH2D *hresponseAbsZRapidity_DifJetRapiditys2_Zinc2jet;
    TH2D *hresponseAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet;
    TH2D *hresponseSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet;
    TH2D *hresponseDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet;
};

#endif
