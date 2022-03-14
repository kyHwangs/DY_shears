#include "functions.h"
#include "ConfigVJets.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TString.h"
#include "TSystem.h"
#include <algorithm>
#include <cstdarg>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdlib.h>
#include <vector>

#ifndef DYJETS_NEW_API
#   include "RooUnfoldResponse.h"
#endif // DYJETS_NEW_API

extern ConfigVJets cfg;

using namespace std;

bool LepDescendingOrder(leptonStruct l1, leptonStruct l2) { return (l1.v.Pt() > l2.v.Pt()); }

bool JetDescendingOrder(jetStruct j1, jetStruct j2) { return (j1.v.Pt() > j2.v.Pt()); }

//--- for WJets ---
bool JetYDescendingOrder(TLorentzVector tj1, TLorentzVector tj2)
{
    return (tj1.Rapidity() > tj2.Rapidity());
}

double deltaRYPhi(TLorentzVector j1, TLorentzVector j2)
{
    double dY = j1.Rapidity() - j2.Rapidity();
    double dPhi = deltaPhi(j1, j2);
    return sqrt(dY * dY + dPhi * dPhi);
}
//-----------------

vector<double> makeVector(int num, ...)
{
    va_list list;
    va_start(list, num);
    vector<double> vec;
    for (int i(0); i < num; i++) {
        double next = va_arg(list, double);
        vec.push_back(next);
    }
    va_end(list);
    return vec;
}

void insertVector(vector<double> &veca, int num, ...)
{
    va_list list;
    va_start(list, num);
    vector<double> vecb;
    for (int i(0); i < num; i++) {
        double next = va_arg(list, double);
        vecb.push_back(next);
    }
    va_end(list);
    veca.insert(veca.end(), vecb.begin(), vecb.end());
}

double phi0to2pi(double phi)
{
    double pi = 3.141592653589793238;
    while (phi >= 2. * pi) phi -= 2. * pi;
    while (phi < 0.) phi += 2. * pi;
    return phi;
}

double ZPtviaPhistar(double phistar)
{
    return log(-3.36739e+00 / (-5.14753e-01 - phistar) - 1.) * (-6.48782e+01) +
           1.11148e+02; // from fit "[0] - [1]/(1+exp((x-[2])/[3]))"
}

double deltaPhi(TLorentzVector v1, TLorentzVector v2)
{
    // build the delta Phi angle between the two vectors
    double pi = 3.141592653589793238;
    double phi1 = phi0to2pi(v1.Phi());
    double phi2 = phi0to2pi(v2.Phi());
    double dPhi = phi0to2pi(phi1 - phi2);
    dPhi = (dPhi > (2 * pi - dPhi)) ? 2 * pi - dPhi : dPhi;
    return dPhi;
}

double deltaPhi(double Phi1, double Phi2)
{
    // build the delta Phi angle between the two vectors
    double pi = 3.141592653589793238;
    double phi1 = phi0to2pi(Phi1);
    double phi2 = phi0to2pi(Phi2);
    double dPhi = phi0to2pi(phi1 - phi2);
    dPhi = (dPhi > (2 * pi - dPhi)) ? 2 * pi - dPhi : dPhi;
    // cout << "      DeltaPhi: " << endl;
    // cout << "      phi1 = " << phi1 << "  phi2 = " << phi2 << endl;
    // cout << "      DeltaPhi = " << dPhi << endl;
    return dPhi;
}

double deltaR(TLorentzVector v1, TLorentzVector v2)
{
    double dEta = v1.Eta() - v2.Eta();
    double dPhi = deltaPhi(v1, v2);
    return sqrt(dEta * dEta + dPhi * dPhi);
}

double deltaR(double Phi1, double Eta1, double Phi2, double Eta2)
{
    double dEta = Eta1 - Eta2;
    double dPhi = deltaPhi(Phi1, Phi2);
    return sqrt(dEta * dEta + dPhi * dPhi);
}

double PHI(TLorentzVector l1, TLorentzVector l2, TLorentzVector j1, TLorentzVector j2)
{
    // build the angle PHI between the two subsytems (l1+l2, j1+j2) vectors
    double lPx = (l1.Px() + l2.Px());
    double lPy = (l1.Py() + l2.Py());
    double lPz = (l1.Pz() + l2.Pz());
    double lNorm = sqrt(lPx * lPx + lPy * lPy + lPz * lPz);
    double jPx = (j1.Px() + j2.Px());
    double jPy = (j1.Py() + j2.Py());
    double jPz = (j1.Pz() + j2.Pz());
    double jNorm = sqrt(jPx * jPx + jPy * jPy + jPz * jPz);
    return acos((jPx * lPx + jPy * lPy + jPz * lPz) / (jNorm * lNorm));
}

double PHI_T(TLorentzVector l1, TLorentzVector l2, TLorentzVector j1, TLorentzVector j2)
{
    // build the angle PHI between the two subsytems (l1+l2, j1+j2) vectors in the transverse plane
    double lPx = (l1.Px() + l2.Px());
    double lPy = (l1.Py() + l2.Py());
    double lNorm = sqrt(lPx * lPx + lPy * lPy);
    double jPx = (j1.Px() + j2.Px());
    double jPy = (j1.Py() + j2.Py());
    double jNorm = sqrt(jPx * jPx + jPy * jPy);
    return acos((jPx * lPx + jPy * lPy) / (jNorm * lNorm));
}

double SpTsub(TLorentzVector v1, TLorentzVector v2)
{
    return sqrt(pow(v1.Px() + v2.Px(), 2) + pow(v1.Py() + v2.Py(), 2)) / (v1.Pt() + v2.Pt());
}

double SpT(TLorentzVector l1, TLorentzVector l2, TLorentzVector j1, TLorentzVector j2)
{
    return sqrt(pow(SpTsub(l1, l2), 2) + pow(SpTsub(j1, j2), 2)) / sqrt(2.);
}

double SPhi(TLorentzVector l1, TLorentzVector l2, TLorentzVector j1, TLorentzVector j2)
{
    return sqrt(deltaPhi(l1, l2) * deltaPhi(l1, l2) + deltaPhi(j1, j2) * deltaPhi(j1, j2)) /
           sqrt(2.);
}

char GetRunData(int runNumber)
{
    if ((RUNB_2016 <= runNumber) && (runNumber < RUNC_2016)) {
        return 'B';
    }
    if ((RUNC_2016 <= runNumber) && (runNumber < RUND_2016)) {
        return 'C';
    }
    if ((RUND_2016 <= runNumber) && (runNumber < RUNE_2016)) {
        return 'D';
    }
    if ((RUNE_2016 <= runNumber) && (runNumber < RUNF_2016)) {
        return 'E';
    }
    if ((RUNF_2016 <= runNumber) && (runNumber < RUNG_2016)) {
        return 'F';
    }
    if ((RUNG_2016 <= runNumber) && (runNumber < RUNH_2016)) {
        return 'G';
    }
    if (RUNH_2016 <= runNumber) {
        return 'H';
    }
    return 'Z';
}

char GetRunMC(Long64_t *mcEraBoundary, Long64_t eventNumber)
{
    const size_t _nLetters = 7;
    const std::string _runLettersAll = "BCDEFGH";
    for (size_t iLetter = 0; iLetter < _nLetters; iLetter++) {
        if (eventNumber <= mcEraBoundary[iLetter]) return _runLettersAll[iLetter];
    }
    return 'Z';
}

double SmearLepPt(double recoPt, double genPt, int smearlepton, double smearFactor)
{

    double smearedPt(0);

    if (smearlepton == 0) {
        smearedPt = std::max(0., recoPt);
    }

    else if (smearlepton == 1) {
        smearedPt = std::max(0., genPt + (1.0 + smearFactor) * (recoPt - genPt));
    }

    else if (smearlepton == -1) {
        smearedPt = std::max(0., genPt + (1.0 - smearFactor) * (recoPt - genPt));
    }

    return smearedPt;
}

// Helper function for SmearsJetPt
double GetJetSF(double eta, int direction)
{
    // Fall 2015 resolution scale factor
    // twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    size_t year = 2016;
    double centralSF(1.00);
    double upSF(1.00);
    double downSF(1.00);
    if (year == 2015) {
        if (fabs(eta) < 0.5)
            centralSF = 1.095;
        else if (fabs(eta) < 0.8)
            centralSF = 1.120;
        else if (fabs(eta) < 1.1)
            centralSF = 1.097;
        else if (fabs(eta) < 1.3)
            centralSF = 1.103;
        else if (fabs(eta) < 1.7)
            centralSF = 1.118;
        else if (fabs(eta) < 1.9)
            centralSF = 1.100;
        else if (fabs(eta) < 2.1)
            centralSF = 1.162;
        else if (fabs(eta) < 2.3)
            centralSF = 1.160;
        else if (fabs(eta) < 2.5)
            centralSF = 1.161;
        else if (fabs(eta) < 2.8)
            centralSF = 1.209;
        else if (fabs(eta) < 3.0)
            centralSF = 1.564;
        else if (fabs(eta) < 3.2)
            centralSF = 1.384;
        else if (fabs(eta) < 5.0)
            centralSF = 1.216;
        else
            centralSF = 1.320;

        if (fabs(eta) < 0.5)
            upSF = 1.095 + 0.018;
        else if (fabs(eta) < 0.8)
            upSF = 1.120 + 0.028;
        else if (fabs(eta) < 1.1)
            upSF = 1.097 + 0.017;
        else if (fabs(eta) < 1.3)
            upSF = 1.103 + 0.033;
        else if (fabs(eta) < 1.7)
            upSF = 1.118 + 0.014;
        else if (fabs(eta) < 1.9)
            upSF = 1.100 + 0.033;
        else if (fabs(eta) < 2.1)
            upSF = 1.162 + 0.044;
        else if (fabs(eta) < 2.3)
            upSF = 1.160 + 0.048;
        else if (fabs(eta) < 2.5)
            upSF = 1.161 + 0.060;
        else if (fabs(eta) < 2.8)
            upSF = 1.209 + 0.059;
        else if (fabs(eta) < 3.0)
            upSF = 1.564 + 0.321;
        else if (fabs(eta) < 3.2)
            upSF = 1.384 + 0.033;
        else if (fabs(eta) < 5.0)
            upSF = 1.216 + 0.050;
        else
            upSF = 1.606;

        if (fabs(eta) < 0.5)
            downSF = 1.095 - 0.018;
        else if (fabs(eta) < 0.8)
            downSF = 1.120 - 0.028;
        else if (fabs(eta) < 1.1)
            downSF = 1.097 - 0.017;
        else if (fabs(eta) < 1.3)
            downSF = 1.103 - 0.033;
        else if (fabs(eta) < 1.7)
            downSF = 1.118 - 0.014;
        else if (fabs(eta) < 1.9)
            downSF = 1.100 - 0.033;
        else if (fabs(eta) < 2.1)
            downSF = 1.162 - 0.044;
        else if (fabs(eta) < 2.3)
            downSF = 1.160 - 0.048;
        else if (fabs(eta) < 2.5)
            downSF = 1.161 - 0.060;
        else if (fabs(eta) < 2.8)
            downSF = 1.209 - 0.059;
        else if (fabs(eta) < 3.0)
            downSF = 1.564 - 0.321;
        else if (fabs(eta) < 3.2)
            downSF = 1.384 - 0.033;
        else if (fabs(eta) < 5.0)
            downSF = 1.216 - 0.050;
        else
            downSF = 1.034;
    }
    if (year == 2016) {
        if (fabs(eta) < 0.5)
            centralSF = 1.109;
        else if (fabs(eta) < 0.8)
            centralSF = 1.138;
        else if (fabs(eta) < 1.1)
            centralSF = 1.114;
        else if (fabs(eta) < 1.3)
            centralSF = 1.123;
        else if (fabs(eta) < 1.7)
            centralSF = 1.084;
        else if (fabs(eta) < 1.9)
            centralSF = 1.082;
        else if (fabs(eta) < 2.1)
            centralSF = 1.140;
        else if (fabs(eta) < 2.3)
            centralSF = 1.067;
        else if (fabs(eta) < 2.5)
            centralSF = 1.177;
        else if (fabs(eta) < 2.8)
            centralSF = 1.364;
        else if (fabs(eta) < 3.0)
            centralSF = 1.857;
        else if (fabs(eta) < 3.2)
            centralSF = 1.328;
        else if (fabs(eta) < 5.0)
            centralSF = 1.160;
        else
            centralSF = 1.320;

        if (fabs(eta) < 0.5)
            centralSF = 1.109 + 0.008;
        else if (fabs(eta) < 0.8)
            centralSF = 1.138 + 0.013;
        else if (fabs(eta) < 1.1)
            centralSF = 1.114 + 0.013;
        else if (fabs(eta) < 1.3)
            centralSF = 1.123 + 0.024;
        else if (fabs(eta) < 1.7)
            centralSF = 1.084 + 0.011;
        else if (fabs(eta) < 1.9)
            centralSF = 1.082 + 0.035;
        else if (fabs(eta) < 2.1)
            centralSF = 1.140 + 0.047;
        else if (fabs(eta) < 2.3)
            centralSF = 1.067 + 0.053;
        else if (fabs(eta) < 2.5)
            centralSF = 1.177 + 0.041;
        else if (fabs(eta) < 2.8)
            centralSF = 1.364 + 0.039;
        else if (fabs(eta) < 3.0)
            centralSF = 1.857 + 0.071;
        else if (fabs(eta) < 3.2)
            centralSF = 1.328 + 0.022;
        else if (fabs(eta) < 5.0)
            centralSF = 1.160 + 0.029;
        else
            upSF = 1.606;

        if (fabs(eta) < 0.5)
            centralSF = 1.109 - 0.008;
        else if (fabs(eta) < 0.8)
            centralSF = 1.138 - 0.013;
        else if (fabs(eta) < 1.1)
            centralSF = 1.114 - 0.013;
        else if (fabs(eta) < 1.3)
            centralSF = 1.123 - 0.024;
        else if (fabs(eta) < 1.7)
            centralSF = 1.084 - 0.011;
        else if (fabs(eta) < 1.9)
            centralSF = 1.082 - 0.035;
        else if (fabs(eta) < 2.1)
            centralSF = 1.140 - 0.047;
        else if (fabs(eta) < 2.3)
            centralSF = 1.067 - 0.053;
        else if (fabs(eta) < 2.5)
            centralSF = 1.177 - 0.041;
        else if (fabs(eta) < 2.8)
            centralSF = 1.364 - 0.039;
        else if (fabs(eta) < 3.0)
            centralSF = 1.857 - 0.071;
        else if (fabs(eta) < 3.2)
            centralSF = 1.328 - 0.022;
        else if (fabs(eta) < 5.0)
            centralSF = 1.160 - 0.029;
        else
            downSF = 1.034;
    }

    if (direction == 0)
        return centralSF;
    else if (direction == 1)
        return upSF;
    else if (direction == -1)
        return downSF;
    else {
        std::cerr
            << "function.cc::GetJetSF: Something went wrong getting the SF for jet smearing\n";
        abort();
    }
}

// Helper function for SmearsJetPt
double GetJetResolution() { return 0.10; }

// Smearing when a gen jet match is found:
double SmearJetPt(double recoPt, double genPt, double recoJetEta, int direction)
{
    double smearedPt = recoPt;
    smearedPt = std::max(0.0, genPt + GetJetSF(recoJetEta, direction) * (recoPt - genPt));
    return smearedPt;
}
// Smearing with stochastic method (no match found)
double SmearJetPt(double recoPt, double recoJetEta, int direction)
{
    TRandom3 *random = new TRandom3();
    double smearedPt = recoPt;
    double smearFactor = 1 +
                         random->Gaus(0, GetJetResolution()) *
                             sqrt(std::max(pow(GetJetSF(recoJetEta, direction), 2) - 1.0, 0.0));
    smearedPt = smearedPt * smearFactor;
    // printf("SmearFactor = %F\n",smearFactor);
    return smearedPt;
}

void bestTwoJetsCandidatesPt(vector<jetStruct> jets,
                             pair<TLorentzVector, TLorentzVector> &bestTwoJets)
{
    int nGoodJets(jets.size());
    if (nGoodJets >= 2) {
        // cout << "\nMore than 2 jets, selecting best pair" << endl;
        double minPt(999999.);
        for (int i(0); i < nGoodJets - 1; i++) {
            TLorentzVector jeti = jets[i].v;
            for (int j(i + 1); j < nGoodJets; j++) {
                TLorentzVector jetj = jets[j].v;
                TLorentzVector jetij = jeti + jetj;
                // cout << i << " " << j << ": Pair pt = " << jetij.Pt() << endl;
                if (jetij.Pt() < minPt) {
                    bestTwoJets.first = jeti;
                    bestTwoJets.second = jetj;
                    minPt = jetij.Pt();
                    // cout << "Smallest pt = " << jetij.Pt() << endl;
                }
            }
        }
    }
}

void bestTwoJetsCandidatesPhi(vector<jetStruct> jets,
                              pair<TLorentzVector, TLorentzVector> &bestTwoJets)
{
    int nGoodJets(jets.size());
    if (nGoodJets >= 2) {
        // cout << "\nMore than 2 jets, selecting best pair" << endl;
        double maxdPhi(-0.0001);
        for (int i(0); i < nGoodJets - 1; i++) {
            TLorentzVector jeti = jets[i].v;
            for (int j(i + 1); j < nGoodJets; j++) {
                TLorentzVector jetj = jets[j].v;
                double dPhi = deltaPhi(jeti, jetj);
                // cout << i << " " << j << ": dPhi = " << dPhi << endl;
                if (dPhi > maxdPhi) {
                    bestTwoJets.first = jeti;
                    bestTwoJets.second = jetj;
                    maxdPhi = dPhi;
                    // cout << "Biggest dPhi = " << dPhi << endl;
                }
            }
        }
    }
}

void BTagModification(double randNumber, double pt, double eta, int jetFlavour, bool &passBJets)
{
    double x = 0.679; /// discrim_cut;
    if (abs(jetFlavour) == 5) {

        float effb = -1.73338329789 * x * x * x * x + 1.26161794785 * x * x * x +
                     0.784721653518 * x * x + -1.03328577451 * x + 1.04305075822;
        float SFb = (0.938887 + (0.00017124 * pt)) + (-2.76366e-07 * (pt * pt));

        float SFb_error = 0;

        if (pt >= 20 && pt < 30)
            SFb_error = 0.0415707;
        else if (pt >= 30 && pt < 40)
            SFb_error = 0.0204209;
        else if (pt >= 40 && pt < 50)
            SFb_error = 0.0223227;
        else if (pt >= 50 && pt < 60)
            SFb_error = 0.0206655;
        else if (pt >= 60 && pt < 70)
            SFb_error = 0.0199325;
        else if (pt >= 70 && pt < 80)
            SFb_error = 0.0174121;
        else if (pt >= 80 && pt < 100)
            SFb_error = 0.0202332;
        else if (pt >= 100 && pt < 120)
            SFb_error = 0.0182446;
        else if (pt >= 120 && pt < 160)
            SFb_error = 0.0159777;
        else if (pt >= 160 && pt < 210)
            SFb_error = 0.0218531;
        else if (pt >= 210 && pt < 260)
            SFb_error = 0.0204688;
        else if (pt >= 260 && pt < 320)
            SFb_error = 0.0265191;
        else if (pt >= 320 && pt < 400)
            SFb_error = 0.0313175;
        else if (pt >= 400 && pt < 500)
            SFb_error = 0.0415417;
        else if (pt >= 500 && pt < 600)
            SFb_error = 0.0740446;
        else if (pt >= 600)
            SFb_error = 0.0596716;

        float SFb_up = SFb + SFb_error;
        float SFb_down = SFb - SFb_error;

        // F values for rand comparison
        float f = 0.0;
        float f_up = 0.0;
        float f_down = 0.0;

        if (SFb < 1.0) f = (1.0 - SFb);
        if (SFb_up < 1.0) f_up = (1.0 - SFb_up);
        if (SFb_down < 1.0) f_down = (1.0 - SFb_down);

        if (SFb > 1.0) f = (1.0 - SFb) / (1.0 - 1.0 / effb);
        if (SFb_up > 1.0) f_up = (1.0 - SFb_up) / (1.0 - 1.0 / effb);
        if (SFb_down > 1.0) f_down = (1.0 - SFb_down) / (1.0 - 1.0 / effb);

        bool passBJets_SFB_sys_up = passBJets; // Initialize the systematic_up as the central value
        bool passBJets_SFB_sys_down =
            passBJets; // Initialize the systematic_down as the central value

        // Untag a tagged jet

        if (passBJets && SFb < 1.0 && randNumber < f) passBJets = false; // for central value
        if (passBJets_SFB_sys_up && SFb < 1.0 && randNumber < f_up)
            passBJets_SFB_sys_up = false; // for systematic_up
        if (passBJets_SFB_sys_down && SFb < 1.0 && randNumber < f_down)
            passBJets_SFB_sys_down = false; // for sytematic_down

        // Tag an untagged jet
        if (!passBJets && SFb > 1.0 && randNumber < f) passBJets = true; // for central value
        if (!passBJets_SFB_sys_up && SFb > 1.0 && randNumber < f_up)
            passBJets_SFB_sys_up = true; // for systematic_up
        if (!passBJets_SFB_sys_down && SFb > 1.0 && randNumber < f_down)
            passBJets_SFB_sys_down = true; // for sytematic_down

    }
    // ---------------- For Real C-jets--------------- //
    else if (abs(jetFlavour) == 4) {

        float effc = -1.5734604211 * x * x * x * x + 1.52798999269 * x * x * x +
                     0.866697059943 * x * x + -1.66657942274 * x + 0.780639301724;
        float SFc = (0.938887 + (0.00017124 * pt)) + (-2.76366e-07 * (pt * pt));

        float SFc_error = 0.;

        if (pt >= 20 && pt < 30)
            SFc_error = 0.0415707;
        else if (pt >= 30 && pt < 40)
            SFc_error = 0.0204209;
        else if (pt >= 40 && pt < 50)
            SFc_error = 0.0223227;
        else if (pt >= 50 && pt < 60)
            SFc_error = 0.0206655;
        else if (pt >= 60 && pt < 70)
            SFc_error = 0.0199325;
        else if (pt >= 70 && pt < 80)
            SFc_error = 0.0174121;
        else if (pt >= 80 && pt < 100)
            SFc_error = 0.0202332;
        else if (pt >= 100 && pt < 120)
            SFc_error = 0.0182446;
        else if (pt >= 120 && pt < 160)
            SFc_error = 0.0159777;
        else if (pt >= 160 && pt < 210)
            SFc_error = 0.0218531;
        else if (pt >= 210 && pt < 260)
            SFc_error = 0.0204688;
        else if (pt >= 260 && pt < 320)
            SFc_error = 0.0265191;
        else if (pt >= 320 && pt < 400)
            SFc_error = 0.0313175;
        else if (pt >= 400 && pt < 500)
            SFc_error = 0.0415417;
        else if (pt >= 500 && pt < 600)
            SFc_error = 0.0740446;
        else if (pt >= 600)
            SFc_error = 0.0596716;

        float SFc_up = SFc + 2 * SFc_error;
        float SFc_down = SFc - 2 * SFc_error;

        // F values for rand comparison
        float f = 0.0;
        float f_up = 0.0;
        float f_down = 0.0;

        if (SFc < 1.0) f = (1.0 - SFc);
        if (SFc_up < 1.0) f_up = (1.0 - SFc_up);
        if (SFc_down < 1.0) f_down = (1.0 - SFc_down);

        if (SFc > 1.0) f = (1.0 - SFc) / (1.0 - 1.0 / effc);
        if (SFc_up > 1.0) f_up = (1.0 - SFc_up) / (1.0 - 1.0 / effc);
        if (SFc_down > 1.0) f_down = (1.0 - SFc_down) / (1.0 - 1.0 / effc);

        bool passBJets_SFB_sys_up = passBJets; // Initialize the systematic_up as the central value
        bool passBJets_SFB_sys_down =
            passBJets; // Initialize the systematic_down as the central value

        // Untag a tagged jet

        if (passBJets && SFc < 1.0 && randNumber < f) passBJets = false; // for central value
        if (passBJets_SFB_sys_up && SFc < 1.0 && randNumber < f_up)
            passBJets_SFB_sys_up = false; // for systematic_up
        if (passBJets_SFB_sys_down && SFc < 1.0 && randNumber < f_down)
            passBJets_SFB_sys_down = false; // for sytematic_down

        // Tag an untagged jet
        if (!passBJets && SFc > 1.0 && randNumber < f) passBJets = true; // for central value
        if (!passBJets_SFB_sys_up && SFc > 1.0 && randNumber < f_up)
            passBJets_SFB_sys_up = true; // for systematic_up
        if (!passBJets_SFB_sys_down && SFc > 1.0 && randNumber < f_down)
            passBJets_SFB_sys_down = true; // for sytematic_down

    }
    // ---------------- For REAL Light-jets --------------- //
    else if (abs(jetFlavour) < 4) {

        float SFlight = 1.0;
        float SFlight_up = 1.0;
        float SFlight_down = 1.0;
        float eff_l = 0.0;

        if (fabs(eta) <= 0.8) {
            SFlight = (((1.07541 + (0.00231827 * pt)) + (-4.74249e-06 * (pt * pt))) +
                       (2.70862e-09 * (pt * (pt * pt))));
            SFlight_up = (((1.18638 + (0.00314148 * pt)) + (-6.68993e-06 * (pt * pt))) +
                          (3.89288e-09 * (pt * (pt * pt))));
            SFlight_down = (((0.964527 + (0.00149055 * pt)) + (-2.78338e-06 * (pt * pt))) +
                            (1.51771e-09 * (pt * (pt * pt))));
            eff_l = ((0.00967751 + (2.54564e-05 * pt)) + (-6.92256e-10 * (pt * pt)));
        } else if (fabs(eta) <= 1.6) {
            SFlight = (((1.05613 + (0.00114031 * pt)) + (-2.56066e-06 * (pt * pt))) +
                       (1.67792e-09 * (pt * (pt * pt))));
            SFlight_up = (((1.16624 + (0.00151884 * pt)) + (-3.59041e-06 * (pt * pt))) +
                          (2.38681e-09 * (pt * (pt * pt))));
            SFlight_down = (((0.946051 + (0.000759584 * pt)) + (-1.52491e-06 * (pt * pt))) +
                            (9.65822e-10 * (pt * (pt * pt))));
            eff_l = ((0.00974141 + (5.09503e-05 * pt)) + (2.0641e-08 * (pt * pt)));
        } else if (fabs(eta) <= 2.4) {
            SFlight = (((1.05625 + (0.000487231 * pt)) + (-2.22792e-06 * (pt * pt))) +
                       (1.70262e-09 * (pt * (pt * pt))));
            SFlight_up = (((1.15575 + (0.000693344 * pt)) + (-3.02661e-06 * (pt * pt))) +
                          (2.39752e-09 * (pt * (pt * pt))));
            SFlight_down = (((0.956736 + (0.000280197 * pt)) + (-1.42739e-06 * (pt * pt))) +
                            (1.0085e-09 * (pt * (pt * pt))));
            eff_l = ((0.013595 + (0.000104538 * pt)) + (-1.36087e-08 * (pt * pt)));
        }

        // F values for rand comparison
        float f = 0.0;
        float f_up = 0.0;
        float f_down = 0.0;

        if (SFlight < 1.0) f = (1.0 - SFlight);
        if (SFlight_up < 1.0) f_up = (1.0 - SFlight_up);
        if (SFlight_down < 1.0) f_down = (1.0 - SFlight_down);

        if (SFlight > 1.0) f = (1.0 - SFlight) / (1.0 - 1.0 / eff_l);
        if (SFlight_up > 1.0) f_up = (1.0 - SFlight_up) / (1.0 - 1.0 / eff_l);
        if (SFlight_down > 1.0) f_down = (1.0 - SFlight_down) / (1.0 - 1.0 / eff_l);

        bool passBJets_SFB_sys_up = passBJets; // Initialize the systematic_up as the central value
        bool passBJets_SFB_sys_down =
            passBJets; // Initialize the systematic_down as the central value

        // Untag a tagged jet

        if (passBJets && SFlight < 1.0 && randNumber < f) passBJets = false; // for central value
        if (passBJets_SFB_sys_up && SFlight < 1.0 && randNumber < f_up)
            passBJets_SFB_sys_up = false; // for systematic_up
        if (passBJets_SFB_sys_down && SFlight < 1.0 && randNumber < f_down)
            passBJets_SFB_sys_down = false; // for sytematic_down

        // Tag an untagged jet
        if (!passBJets && SFlight > 1.0 && randNumber < f) passBJets = true; // for central value
        if (!passBJets_SFB_sys_up && SFlight > 1.0 && randNumber < f_up)
            passBJets_SFB_sys_up = true; // for systematic_up
        if (!passBJets_SFB_sys_down && SFlight > 1.0 && randNumber < f_down)
            passBJets_SFB_sys_down = true; // for sytematic_down

    } ////////flavour lop
}
/*
double BTagweight(double randNumber, double pt, double eta, int jetFlavour, bool passBJets) {
   std::string tg="";
   BTagEntry::JetFlavor flavor;
    if (std::abs(jetFlavour) == 5) {
        flavor = BTagEntry::FLAV_B;
        tg="bjet";
    } else if (std::abs(jetFlavour) == 4) {
        flavor = BTagEntry::FLAV_C;
        tg="cjet";
    } else {
        flavor = BTagEntry::FLAV_UDSG;
        tg="udsgjet";
    }

    double sf = _btag_calibration_reader.eval_auto_bounds(
        "central", flavor, std::abs(j.v.Eta()), j.v.Pt());
   // double eff = _bjet_tag_eff[flavor];
   //
}*/
FILE *eosOpen(const char *path, int (**closeFunc)(FILE *))
{
    TString tspath(path);
    if (tspath.BeginsWith("root://")) {
        *closeFunc = pclose;
        tspath.Remove(0, strlen("root://"));
        Ssiz_t p = tspath.First("/");
        if (p == TString::kNPOS) return 0;
        TString server(tspath(0, p));
        TString filepath(tspath(p + 1, tspath.Length() - p));
        //   std::cout << ">>> server: " << server << ", path: " << filepath << "\n";
        return popen(TString::Format("xrdfs %s cat %s", server.Data(), filepath.Data()), "r");
    } else {
        *closeFunc = fclose;
        return fopen(path, "r");
    }
}

bool isRootFile(const char *path)
{
    int (*funcClose)(FILE *);
    FILE *f = eosOpen(path, &funcClose);

    char buffer[5];
    if (f) {
        bool rc = (fread(buffer, 5, 1, f) == 1) && (memcmp(buffer, "root\0", 5) == 0);
        funcClose(f);
        return rc;
    } else {
        std::cerr << "Warning: failed to read file " << path
                  << ". File type determined from its extension.\n";
        return TString(path).EndsWith(".root");
    }
}


#ifndef DYJETS_NEW_API

bool mergeHistFiles(const std::vector<std::string> &src, const std::string &dest)
{
    TH1::SetDefaultSumw2();

    std::cerr << "Creating file " << dest << std::endl;
    TFile fout(dest.c_str(), "RECREATE");
    if (fout.IsZombie()) {
        std::cerr << "Error: failed to create file " << dest << " (" << __FILE__ << ":" << __LINE__
                  << ").\n\n";
        return false;
    }

    union u
    {
        TObject *o;
        TH1 *h;
        RooUnfoldResponse *r;
    } obj;
    std::map<std::string, u> objs;
    for (unsigned iFile = 0; iFile < src.size(); ++iFile) {
        TFile in(src[iFile].c_str());
        if (in.IsZombie()) continue;
        TIter next(in.GetListOfKeys());
        TKey *key;
        while ((key = (TKey *)next())) {
            obj.o = key->ReadObj();
            enum { kHist, kResp, kOther } type = kOther;
            if (obj.o->InheritsFrom("TH1")) type = kHist;
            if (obj.o->InheritsFrom("RooUnfoldResponse")) type = kResp;
            if (type == kOther) continue;
            if (iFile == 0) {
                if (type == kHist) obj.h->SetDirectory(0);
                objs[obj.o->GetName()].o = obj.o;
            } else {
                std::map<std::string, u>::iterator it = objs.find(obj.o->GetName());
                if (it != objs.end()) {
                    if (type == kHist) (*it).second.h->Add(obj.h);
                    if (type == kResp) (*it).second.r->Add(*(obj.r));
                }
            }
        } // next Object
        in.Close();
    } // next File

    fout.cd();
    for (std::map<std::string, u>::iterator it = objs.begin(); it != objs.end(); ++it) {
        it->second.o->Write();
        delete it->second.o;
    }
    fout.Close();
    return true;
}

#endif // DYJETS_NEW_API

// Check that two Root TAxis have indentical boudaries and binning:
bool isSameBinning(const TAxis &ax1, const TAxis &ax2)
{
    const static bool verbose = true;
    // check number of bins
    if (ax1.GetNbins() != ax2.GetNbins()) return false;

    for (int ibin = 1; ibin <= ax1.GetNbins(); ++ibin) {
        if (verbose)
            std::cout << ax1.GetBinLowEdge(ibin) << " (" << ax2.GetBinLowEdge(ibin) << ")\t";
        if (ax1.GetBinLowEdge(ibin) != ax2.GetBinLowEdge(ibin)) return false;
    }
    if (verbose)
        std::cout << ax1.GetBinUpEdge(ax1.GetNbins()) << " (" << ax2.GetBinUpEdge(ax2.GetNbins())
                  << ")\n";
    if (ax1.GetBinUpEdge(ax1.GetNbins()) != ax2.GetBinUpEdge(ax2.GetNbins())) return false;

    return true;
}

#ifndef DYJETS_NEW_API

void saveCanvas(TCanvas *c, const char *outputDir, const char *baseName)
{
    std::string mainFormat = cfg.getS("mainFormat", "pdf");
    std::vector<std::string> extraFormats_ =
        cfg.getVS("extraFormats", std::vector<std::string>(1, "root"));
    // remove duplicates if any:
    std::set<std::string> extraFormats;
    for (auto s : extraFormats_) {
        extraFormats.insert(s);
    }
    c->SaveAs(TString(outputDir) + "/" + baseName + "." + mainFormat.c_str());
    for (auto ext : extraFormats) {
        gSystem->mkdir(TString(outputDir) + "/" + ext);
        // extra tag for .root file to distinguish with the file containing all the histograms:
        const char *canvas = (ext == "root" ? "_canvas" : "");
        c->SaveAs(TString(outputDir) + "/" + ext + "/" + baseName + canvas + "." + ext.c_str());
    }
}

#endif

/** Round a number to n digits before or after the decimal point.
 * @param val: input and value.
 * @param l10: log10 of the precisison, i.e. position of the least
 * significant digit position:
 *     1-> tens, 0->unit, -1-> tenth, -2 -> hundredth, etc.
 * @param sVal: output value as a string with proper formatting.
 * @return rounded value
 */
double ndec_round(double val, int l10, std::string &sVal)
{
    // following rounding step is needed for l10 > 0
    //(precisions of unit, tens, hundreds,...)
    double a = pow(10, l10);
    double rounded = round(val / a) * a;
    int ndecimals = std::max(0, -l10);
    //  std::cout << "(" << ndecimals << " decimals)\n";
    sVal = TString::Format("%#.*f", ndecimals, rounded);
    return strtod(sVal.c_str(), 0);
}

double nsignif_round(double val, int nsignif, std::string &sVal)
{
    sVal = TString::Format("%#.*g", nsignif, val);
    return strtod(sVal.c_str(), 0);
}

/** Rounds figures of a measurement according to CMS convention
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/Internal/PubGuidelines#Significant_figures_for_measurem
 * rev. 188 and matching the precision of the central value to the precision of the largest
 * uncertainty.
 */
void pground(double val,
             const std::vector<double> &unc,
             std::string &sVal,
             std::vector<std::string> &sUnc,
             bool matchUncPrecOnCentralValue)
{

    const bool verbose = false;

    int nsignif = 2;

    double maxUnc = 0;
    for (auto u : unc) {
        if (u > maxUnc) maxUnc = u;
    }

    // Precision to keep (see ndec_round), 2 significant digits on the largest uncertainty:
    int l10 = floor(log10(maxUnc)) - nsignif + 1;

    // round central value:
    ndec_round(val, l10, sVal);

    sUnc.resize(unc.size());

    std::string s;
    int i = 0;
    for (auto &u : unc) {
        if (matchUncPrecOnCentralValue) {
            ndec_round(u, l10, sUnc[i]);
        } else {
            nsignif_round(u, nsignif, sUnc[i]);
        }
        ++i;
    }

    if (verbose) {
        std::cout << val;
        for (auto u : unc) std::cout << "\t\\pm " << u;
        std::cout << "\nrounded to:\n";
        std::cout << sVal;
        for (auto u : sUnc) std::cout << "\t\\pm " << u;
        std::cout << "\n";
    }
}

void pground(double val, double unc, std::string &sVal, std::string &sUnc)
{
    std::vector<double> uncs(1, unc);
    std::vector<std::string> sUncs;
    pground(val, uncs, sVal, sUncs, false);
    sUnc = sUncs[0];
}
