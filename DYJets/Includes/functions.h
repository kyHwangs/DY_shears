#ifndef _functions_h_
#define _functions_h_

#include <TLorentzVector.h>
#include <cstdarg>
#include <iostream>
#include <vector>
#include "BTagCalibrationStandalone.h"


#include "lepton.h"
#include "tables.h"

// Beginning of run for 2016
#define RUNB_2016 273150
#define RUNC_2016 275656
#define RUND_2016 276315
#define RUNE_2016 276831
#define RUNF_2016 277932
#define RUNG_2016 278820
#define RUNH_2016 281613

// Integrated Lumi for each run 2016
#define LUMI_RUNB_2016 5.748
#define LUMI_RUNC_2016 2.573
#define LUMI_RUND_2016 4.248
#define LUMI_RUNE_2016 4.009
#define LUMI_RUNF_2016 3.102
#define LUMI_RUNG_2016 7.540
#define LUMI_RUNH_2016 8.606

class TAxis;
class TCanvas;

using namespace std;

void barre_de_progression(int);

struct leptonStruct : public physics::lepton
{

    leptonStruct();
    leptonStruct(float pt_,
                 float eta_,
                 float phi_,
                 float en_,
                 float charge_,
                 unsigned id_,
                 float iso_,
                 float scEta_,
                 int trigger_,
                 char flavor_,
                 int TkLayerCnt_)
    {
        v.SetPtEtaPhiE(pt_, eta_, phi_, en_);
        charge = charge_;
        id = id_;
        iso = iso_;
        scEta = scEta_;
        trigger = trigger_;
        flavor = flavor_;
        TkLayerCnt = TkLayerCnt_;
    }

    TLorentzVector v;
    double scEta;
    int TkLayerCnt;
    int trigger;
    char flavor;
};

struct jetStruct
{
    jetStruct();
    jetStruct(double pt_, double eta_, double phi_, double en_, int patIndex_, float jetw_)
    {
        v.SetPtEtaPhiE(pt_, eta_, phi_, en_);
        patIndex = patIndex_;
        jetw = jetw_;
    }

    void setSmearMatch(bool match) { smearMatch = match; }
    void setGenMatchIndex(size_t match) { genMatchIndex = match; }

    TLorentzVector v;
    int patIndex;
    float jetw;
    bool smearMatch; // True implies a gen match was found for smearing, vs false is guassian smear.
    size_t genMatchIndex;
};

bool LepDescendingOrder(leptonStruct, leptonStruct);
bool JetDescendingOrder(jetStruct, jetStruct);
//--- for WJets ---
bool JetYDescendingOrder(TLorentzVector, TLorentzVector);
double deltaRYPhi(TLorentzVector, TLorentzVector);
//-----------------

vector<double> makeVector(int num, ...);
void insertVector(vector<double> &veca, int num, ...);

double phi0to2pi(double);
double ZPtviaPhistar(double);

double deltaPhi(TLorentzVector, TLorentzVector);

double deltaPhi(double, double);
double deltaR(TLorentzVector, TLorentzVector);
double deltaR(double, double, double, double);
double PHI(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
double PHI_T(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
double SpTsub(TLorentzVector, TLorentzVector);
double SpT(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
double SPhi(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);

using table = util::table;

char GetRunData(int runNumber);
char GetRunMC(Long64_t *mcEraBoundary, Long64_t eventNumber);
double SmearLepPt(double recoPt, double genPt, int smearlepton, double smearFactor);
double GetJetSF(double, int);
double GetJetResolution();
double SmearJetPt(double, double, double, int);
double SmearJetPt(double, double, int);
void bestTwoJetsCandidatesPt(vector<jetStruct>, pair<TLorentzVector, TLorentzVector> &);
void bestTwoJetsCandidatesPhi(vector<jetStruct>, pair<TLorentzVector, TLorentzVector> &);
void BTagModification(double randNumber, double pt, double eta, int jetFlavour, bool &passBJets);
double BTagweight(double randNumber, double pt, double eta, int jetFlavour, bool passBJets);

//    BTagCalibrationReader _btag_calibration_reader;

/** Open in read mode a file on eos or a local file. For EOS file,
 * performance might be limited and it should be kept to small file:
 * it executes the "xrdfs cat" command (see man xrdfs).
 * @param path file path. It is assumes an EOS file if it starts with
 * root:// string.
 * @param closeFunc pointer where to store the function to call for close
 * the file instead of fclose.
 * @return pointer to the file stream.
 */
FILE *eosOpen(const char *path, int (**closeFunc)(FILE *));

/** Test if a file is a root file. The test is based on the "magic number"
 * contained in the file which identifies its type. A ROOT file starts with
 * the sequence r,o,o,t,\0
 * @param path path of the file to test.
 * @return true iff the file is a ROOT file.
 */
bool isRootFile(const char *path);

/** Adds histograms and RooUnfoldResponse objects  with same name
 * and definition read from different files and writes the result
 * in a new file.
 * @param src list of input files
 * @param dest output files
 * @return true on success, false on failure
 */
bool mergeHistFiles(const std::vector<std::string> &src, const std::string &dest);

/** Check that two Root TAxis have indentical boudaries and binning:
 * @param ax1 first axis to compare
 * @param ax2 second axis to compare
 * @return true iff the test succeeds
 */

bool isSameBinning(const TAxis &ax1, const TAxis &ax2);

#ifndef DYJETS_NEW_API

// void saveCanvas(const char* fileBaseName, const TCanvas* c = 0);

/** Save a root canvas in the file formats defined in the configuration
 * parameters mainFormat and extraFormats
 */
void saveCanvas(TCanvas *c, const char *outputDir, const char *baseName);

#endif // DYJETS_NEW_API

///@{
/** Rounds figures of a measurement according to CMS convention
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/Internal/PubGuidelines#Significant_figures_for_measurem
 * rev. 188 and matching the precision of the central value to the precision of the largest
 * uncertainty.
 */
void pground(double val,
             const std::vector<double> &unc,
             std::string &sVal,
             std::vector<std::string> &sUnc,
             bool matchUncPrecOnCentralValue);

void pground(double val, double unc, std::string &sVal, std::string &sUnc);

///@}
#endif
