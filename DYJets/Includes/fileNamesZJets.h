//-*- c-basic-offset: 4; -*-
#ifndef _FILENAMESZJETS_h_
#define _FILENAMESZJETS_h_

#include <TString.h>

//-- directory of input root files --------------------
// Commented out, as it should be obtained from the configuration parameter histoDir
// const TString FILESDIRECTORY("HistoFiles/");
//---------- lets add basic information on samples inot common struct
//-------------------------------------------
struct processInfoStruct
{
    TString name;
    char merge; //'+' indicates the sample be merged to the next one
                // until '=' is found. Line with '=' refers to the sample
    // resulting from the merge.
    double NEvents, xsec, xsecFrac, xsecError;
    int colorReco, colorGen;
    TString legendReco, legendGen;
};

//--- first element must point to the data
//--- last element must point to the MC Signal(s)
// 8TeV colours: MLM: kBlue-10, Sherpa: kOrange-2, FXFX: kGreen-8
const processInfoStruct Samples[] = {
    //--  Name    --- merge - #events  -- xsec   - BR - xsec rel. unc. - colorReco - colorGen -
    //legendReco - legendGen
    /* 0*/ {"Data", ' ', 1., 1., 1, 1, kBlack, kBlack, " Data", " Data"},
    /* 1*/ {"TT", ' ', 1., 1., 1, 0.00, kBlue, kBlue, " t#bar{t}", " t#bar{t}"},
    /* 2*/ {"ST_sch",
            '+',
            1.,
            1.,
            1,
            0.06,
            kBlue + 2,
            kBlue + 2,
            " Single Top s-ch",
            " Single Top s-ch"},
    /* 3*/ {"ST_tch",
            '+',
            1.,
            1.,
            1,
            0.06,
            kBlue + 4,
            kBlue + 4,
            " Single Top t-ch",
            " Single Top t-ch"},
    /* 4*/ {"STbar_tW", '+', 1., 1., 1, 0.06, kBlue + 6, kBlue + 6, " #bar{t}W", " #bar{t}W"},
    /* 5*/ {"ST_tW", '+', 1., 1., 1, 0.06, kBlue + 8, kBlue + 8, " tW", " tW"},
    /* 6*/ {"Top", '=', 1., 1., 1, 0.06, kMagenta, kMagenta, " Single top", "Single top"},
    /* 7*/ {"WToLNu", ' ', 1., 1., 1, 0.06, kOrange, kOrange, " W", " W"},
    //    /* 7*/{"TauTau", 	' ',     1.,       1.,      1,  0.06,     kOrange,   kOrange,   "
    //    Z/#gamma^{*} #rightarrow #tau#tau", "Z/#gamma^{*} #rightarrow  #tau#tau"},
    /* 8*/ {"ZZ", '+', 1., 1., 1, 0.06, kOrange, kOrange, " ZZ", " ZZ"},
    /* 9*/ {"WWTo2L2Nu", '+', 1., 1., 1, 0.06, kViolet + 5, kViolet + 5, " WW", " WW"},
    /*10*/ {"WZ", '+', 1., 1., 1, 0.06, kRed + 1, kRed + 1, " WZ", " WZ"},
    /*11*/ {"VV", '=', 1., 1., 1, 0.06, kRed + 1, kRed + 1, " VV", "VV"},
    /*12*/ {"DYJets_UNFOLDING",
            ' ',
            1.,
            1.,
            1,
            0.06,
            kGreen - 8,
            kGreen - 8,
            " Z/#gamma^{*} #rightarrow ll",
            "MG5_aMC + PY8 (#leq 2j NLO + PS)"},
    /*13*/ {"DYJets_MLM",
            ' ',
            1.,
            1.,
            1,
            0.06,
            kBlue - 10,
            kBlue - 10,
            " Z/#gamma^{*} #rightarrow ll",
            "MG5_aMC + PY8 (#leq 4j LO + PS)"},
    /*14*/ {"DYJets_GE",
            ' ',
            1.,
            1.,
            1,
            0.00,
            kBlue - 5,
            kBlue - 5,
            " Z/#gamma^{*} #rightarrow ll",
            "GE + PY8 (NNLL'_{#tau}+NNLO_{0}) #alpha_{s}=0.118"},
    /*15*/ {"DYJets_GEas1135",
            ' ',
            1.,
            1.,
            1,
            0.00,
            kBlue - 2,
            kBlue - 2,
            " Z/#gamma^{*} #rightarrow ll",
            "GE + PY8 (NNLL'_{#tau}+NNLO_{0}) #alpha_{s}=0.1135"},
    /*16*/ {"DYJets_ZjNNLO",
            ' ',
            1.,
            1.,
            1,
            0.00,
            kOrange,
            kOrange,
            " Z/#gamma^{*} #rightarrow ll",
            "N_{jetti} NNLO (1j NNLO)"},
    /*17*/ {"DYJets_GE10as1135",
            ' ',
            1.,
            1.,
            1,
            0.00,
            kBlue - 2,
            kBlue - 2,
            " Z/#gamma^{*} #rightarrow ll",
            "GE + PY8 (NNLL'_{#tau}+NNLO_{0}) #alpha_{s}=0.1135"},
    /*18*/ {"DYJets_GE10as118",
            ' ',
            1.,
            1.,
            1,
            0.00,
            kBlue - 2,
            kBlue - 2,
            " Z/#gamma^{*} #rightarrow ll",
            "GE + PY8 (NNLL'_{#tau}+NNLO_{0})"},
    /* 19*/ {"Data_SMu", ' ', 1., 1., 1, 1, kBlack, kBlack, " Data_SMu", " Data_SMu"},
      /*20*/{"DYJets_FxFx", ' ', 1.,    1.,      1,  0.06,     kGreen-8,  kGreen-8,  " Z/#gamma^{*} #rightarrow ll", "MG5_aMC + PY8 (#leq 2j NLO + PS)"},	
           /*21*/{"DYJets_UNFOLDING_UNC", ' ', 1.,    1.,      1,  0.06,     kGreen-8,  kGreen-8,  " Z/#gamma^{*} #rightarrow ll", "MG5_aMC + PY8 (#leq 2j NLO + PS) Reweighted"},
};

const int NSamples = sizeof(Samples) / sizeof(Samples[0]);
const int DATA = 0;
const int DATA_SMu = 19;
const int DYJETS = 12; // Signal MC

/** Total number of samples after sample grouping, including real data, background MC, and signal MC
 */
const unsigned int NFILESDYJETS = 6;

/** Number of background MC samples: all minus data and the signal MC
 */
const unsigned int NBGDYJETS = NFILESDYJETS - 2;

/** List of indices of samples from Samples to be used for ZJets analysis
 * When samples are grouped, only the merged sample is included in this list.
 * 1st element must be data, last one MC signal and the one in between the MC background
 * samples. In the reco comparison plots the backgound samples are stacked from bottom
 * to top in the order they appear in this list
 */
//
// const unsigned int FilesDYJets[NFILESDYJETS] = {0, 1, 6, 7, 11, 12};
const unsigned int FilesDYJets[NFILESDYJETS] = {0, 7, 6, 11, 1, 12};
// const unsigned int FilesDYJets[NFILESDYJETS] = {0, 1, 6, 7, 11, 14};
// const unsigned int FilesDYJets[NFILESDYJETS] = {0, 1, 2};

// AG
// const TString
// DYPOWHEGFILENAME("DYJetsToLL_M-50_TuneCUETP8M1_8TeV-amcatnloFXFX-Bonzai_fixed_allWeights");
const TString DYAMCATNLOFILENAME("DYJets_UNFOLDING");
const TString DYAMCATNLOLEGEND("MG5_aMC + PY8 (#leq 2j NLO + PS)");
// const TString DYSHERPA2FILENAME("DYJets_Sherpa2_0_16000");
// const TString DYSHERPA2LEGEND("SHERPA 2 (#leq 2j NLO 3,4j LO + PS)");
const TString DYMLM2FILENAME("DYJets_MLM");
const TString DYMLM2LEGEND("MG5_aMC + PY8 (#leq 4j LO + PS)");

const TString DYGEFILENAME("DYJets_GE");
const TString DYGELEGEND("GE + PY8 (NNLL'_{#tau}+NNLO_{0}) #alpha_{s}=0.118");

const TString DYGEAS1135FILENAME("DYJets_GEas1135");
const TString DYGEAS1135LEGEND("GE + PY8 (NNLL'_{#tau}+NNLO_{0}) #alpha_{s}=0.1135");

// ALT_UNFOLDING_FILENAME: alternate signal sample to use to estimate
// unfolding systematic uncertainties. Use empty string to disable
// the calculaiton
const TString ALT_UNFOLDING_FILENAME("DYJets_UNFOLDING_UNC");
// const TString ALT_UNFOLDING_FILENAME("DYJets_Sherpa_Bugra_1_13_UNFOLDING");

const TString DYSHERPA14LEGEND("Sherpa1.4 LO");
// const TString DYMGPYTHIA8FILENAME("DYJetsToLL_M-50_TuneCUETP8M1_8TeV-MG-MLM-Bonzai");
// const TString DYMGPYTHIA8LEGEND("MG+PYthia8 legend");

#endif
