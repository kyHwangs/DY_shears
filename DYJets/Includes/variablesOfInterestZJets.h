#ifndef __VARIABLESOFINTERESTZJETS_H_
#define __VARIABLESOFINTERESTZJETS_H_

#include <TList.h>
#include <TObjString.h>
#include <TObject.h>

struct variableStruct
{
    TString name;
    bool logx;
    bool logy;
    bool decrease;
    int MuSVDkterm, MuBayeskterm, ESVDkterm, EBayeskterm, BayesktermCombine;
};

//--- Z Jets ------------------------------------

const variableStruct VAROFINTERESTZJETS[] = {
    //--  Name  ------------------ log x - log y - decrease - Mu SVD kterm - Mu Bayes kterm - E SVD
    //kterm - E Bayes kterm  - E Bayes kterm COMBINE
    //-- jet Pt distribution
    // PAPER//{"FirstJetPt_Zinc1jet",          1, 1,      1,          11,              5,
    // 9,            9 ,               6    }, // 1
    //{"FirstJetPtMatch_Zinc1jet",          1, 1,      1,          11,              5,
    //9,            9 ,               6    }, // 1
    // PAPER//{"SecondJetPt_Zinc2jet",         1, 1,      1,          14,              8,
    // 11,            7 ,               7    }, // 2
    //{"SecondJetPtMatch_Zinc2jet",         1, 1,      1,          14,              8,
    //11,            7 ,               7    }, // 2
    // PAPER//{"ThirdJetPt_Zinc3jet",          1, 1,      1,           6,              9,
    // 8,            7 ,               7    }, // 3
    //{"ThirdJetPtMatch_Zinc3jet",          1, 1,      1,           6,              9,
    //8,            7 ,               7    }, // 3
    // PAPER//{"FourthJetPt_Zinc4jet",         1, 1,      1,           4,              4,
    // 6,            5 ,               4    }, // 4
    //{"FourthJetPtMatch_Zinc4jet",         1, 1,      1,           4,              4,
    //6,            5 ,               4    }, // 4
    // PAPER//{"FifthJetPt_Zinc5jet",          1, 1,      1,           3,              2,
    // 4,            3 ,               2    }, // 5
    //-- jets HT distribution
    // PAPER//{"JetsHT_Zinc1jet",              1, 1,      1,           9,              4,
    // 7,            3 ,               3    }, // 6
    // PAPER//{"JetsHT_Zinc2jet",              1, 1,      1,           5,              3,
    // 9,            4 ,               2    }, // 7
    // PAPER//{"JetsHT_Zinc3jet",              1, 1,      1,           6,              3,
    // 6,            6 ,               3    }, // 8
    //    {"JetsHT_Zinc4jet",              1, 1,      1,           5,              5,             4,
    //    4 ,               3    }, // 9
    //    {"JetsHT_Zinc5jet",              1, 1,      1,           3,              2,             4,
    //    3 ,               1    }, // 10
    //-- jet Eta distribution
    //    {"FirstJetEta_Zinc1jet",         0, 0,      1,          10,              4,            18,
    //    4 ,               5    }, // 11
    //    {"SecondJetEta_Zinc2jet",        0, 0,      1,           6,              5,            16,
    //    5 ,               5    }, // 12
    //    {"ThirdJetEta_Zinc3jet",         0, 0,      1,           4,              4,             8,
    //    4 ,               5    }, // 13
    //    {"FourthJetEta_Zinc4jet",        0, 0,      1,           4,              4,             4,
    //    5 ,               4    }, // 14
    //    {"FifthJetEta_Zinc5jet",         0, 0,      1,           3,              2,             2,
    //    2 ,               1    }, // 15
    //-- jet Rapidity distribution
    // PAPER//{"FirstJetAbsRapidity_Zinc1jet",         0, 0,      1,          10,              4,
    // 18,            4 ,               5    }, // 11
    //{"FirstJetAbsRapidityMatch_Zinc1jet",         0, 0,      1,          10,              4,
    //18,            4 ,               5    }, // 11
    // PAPER//{"SecondJetAbsRapidity_Zinc2jet",        0, 0,      1,           6,              5,
    // 16,            5 ,               5    }, // 12
    //{"SecondJetAbsRapidityMatch_Zinc2jet",        0, 0,      1,           6,              5,
    //16,            5 ,               5    }, // 12
    // PAPER//{"ThirdJetAbsRapidity_Zinc3jet",         0, 0,      1,           4,              4,
    // 8,            4 ,               5    }, // 13
    //{"ThirdJetAbsRapidityMatch_Zinc3jet",         0, 0,      1,           4,              4,
    //8,            4 ,               5    }, // 13
    // PAPER//{"FourthJetAbsRapidity_Zinc4jet",        0, 0,      1,           4,              4,
    // 4,            5 ,               4    }, // 14
    //{"FourthJetAbsRapidityMatch_Zinc4jet",        0, 0,      1,           4,              4,
    //4,            5 ,               4    }, // 14
    // PAPER//{"FifthJetAbsRapidity_Zinc5jet",         0, 0,      1,           3,              2,
    // 2,            2 ,               1    }, // 15
    //-- Jets Mass distribution
    //{"JetsMass_Zinc2jet",            0, 0,      1,           9,              5,            13,
    //11 ,               5    }, // 16
  
    {"Phistar_Zinc0jetM50_76", 1, 0, 1, 9, 5, 13, 11, 5}, // 17 
    {"Phistar_Zinc0jet", 1, 0, 1, 9, 5, 13, 11, 5},         // 17
    {"Phistar_Zinc0jetM76_106", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"Phistar_Zinc0jetM106_170", 1, 0, 1, 9, 5, 13, 11, 5}, // 17

    {"Phistar_Zinc0jetM170_350", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"Phistar_Zinc0jetM170_inf", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"Phistar_Zinc0jetM111_130", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"Phistar_Zinc0jetM130_170", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"Phistar_Zinc0jetM170_250", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"Phistar_Zinc0jetM250_3", 1, 0, 1, 9, 5, 13, 11, 5},   // 17

   //
   // {"Phistar_Zinc1jetM50_71", 1, 0, 1, 9, 5, 13, 11, 5},         // 17
   // {"Phistar_Zinc1jet", 1, 0, 1, 9, 5, 13, 11, 5},         // 17
   // {"Phistar_Zinc1jetM111_130", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
   // {"Phistar_Zinc1jetM130_170", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
   // {"Phistar_Zinc1jetM170_250", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
   // {"Phistar_Zinc1jetM250_3", 1, 0, 1, 9, 5, 13, 11, 5},   // 17
    //  {"ZPt_Zinc0jetM15_50", 1, 0, 1, 9, 5, 13, 11, 5},       // 17
 //
    {"ZPt_Zinc0jetM50_76", 1, 0, 1, 9, 5, 13, 11, 5},       // 17
    {"ZPt_Zinc0jet", 1, 0, 1, 9, 5, 13, 11, 5},             // 17
    {"ZPt_Zinc0jetM76_106", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM76_106_Mbin", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM106_170", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM170_350", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM170_inf", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM111_130", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM115_135", 1, 0,  1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM111_130", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM130_170", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM170_250", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc0jetM250_3", 1, 0, 1, 9, 5, 13, 11, 5},   // 17
//
    {"ZPt_Zinc1jetM50_76", 1, 0, 1, 9, 5, 13, 11, 5},       // 17
    {"ZPt_Zinc1jet", 1, 0, 1, 9, 5, 13, 11, 5},             // 17
  //   {"ZPt_Zinc0jetM115_135", 1, 0,      1,      9,   5, 13,  11 ,      5    }, // 17
    {"ZPt_Zinc1jetM76_106", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc1jetM76_106_Mbin", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc1jetM106_170", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc1jetM170_350", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc1jetM170_inf", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc1jetM111_130", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc1jetM130_170", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc1jetM170_250", 1, 0, 1, 9, 5, 13, 11, 5}, // 17
    {"ZPt_Zinc1jetM250_3", 1, 0, 1, 9, 5, 13, 11, 5},   // 17
    //  {"ZPt_Zinc0jetM135_170",                 1, 0,      1,           9,              5,
    //  13,           11 ,               5    }, // 17
    //  {"ZPt_Zinc0jetM170_3",                   1, 0,      1,           9,              5,
    //  13,           11 ,               5    }, // 17
    //{"ZAbsRapidity_Zinc1jet",        0, 0,      1,           9,              5,            13,
    //11 ,               5    }, // 18
    //    {"SumZJetRapidity_Zinc1jet",     0, 0,      1,           9,              5,            13,
    //    11 ,               5    }, // 19
    //    {"DifZJetRapidity_Zinc1jet",     0, 0,      1,           9,              5,            13,
    //    11 ,               5    }, // 20
    //    {"FirstJetAbsRapidity_Zinc1jet", 0, 0,      1,           9,              5,            13,
    //    11 ,               5    }, // 21
    //-- jet Pt distribution
    //    {"AbsZRapidity_Zinc1jet",         0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"AbsFirstJetRapidity_Zinc1jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"SumZFirstJetRapidity_Zinc1jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DifZFirstJetRapidity_Zinc1jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"SumZFirstJetEta_Zinc1jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DifZFirstJetEta_Zinc1jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //   {"AbsZRapidity_Zexc1jet",          0, 1,      1,           3,              3,
    //   3,            3 ,               3    }, // 1
    //   {"AbsJetRapidity_Zexc1jet",          0, 1,      1,           3,              3,
    //   3,            3 ,               3    }, // 1
    //   {"SumZJetRapidity_Zexc1jet",          0, 1,      1,           3,              3,
    //   3,            3 ,               3    }, // 1
    //   {"DifZJetRapidity_Zexc1jet",          0, 1,      1,           3,              3,
    //   3,            3 ,               3    }, // 1
    //{"AbsZRapidity_Zinc2jet",         0, 1,      1,           3,              3,             3,
    //3 ,               3    }, // 1
    //    {"AbsFirstJetRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"SumZFirstJetRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DifZFirstJetRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"AbsSecondJetRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"SumZSecondJetRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DifZSecondJetRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"SumFirstSecondJetRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DifFirstSecondJetRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"SumZTwoJetsRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DifZTwoJetsRapidity_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_Zinc1jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZSecondJet_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZThirdJet_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiFirstSecondJet_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiFirstThirdJet_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiSecondThirdJet_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_ZPt150_Zinc1jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_ZPt150_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_ZPt150_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZSecondJet_ZPt150_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZThirdJet_ZPt150_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiFirstSecondJet_ZPt150_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiFirstThirdJet_ZPt150_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiSecondThirdJet_ZPt150_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_ZPt300_Zinc1jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_ZPt300_Zinc2jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_ZPt300_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZSecondJet_ZPt300_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZThirdJet_ZPt300_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiFirstSecondJet_ZPt300_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiFirstThirdJet_ZPt300_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiSecondThirdJet_ZPt300_Zinc3jet",          0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"DPhiZFirstJet_ZPt150_HT300_Zinc3jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DPhiZSecondJet_ZPt150_HT300_Zinc3jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DPhiZThirdJet_ZPt150_HT300_Zinc3jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    // different Z boson PT cuts/////////
    //    {"AbsZRapidity_ZPt100_Zinc1jet",         0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"AbsFirstJetRapidity_ZPt100_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"SumZFirstJetRapidity_ZPt100_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DifZFirstJetRapidity_ZPt100_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"AbsZRapidity_ZPt100_Zinc2jet",         0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"AbsSecondJetRapidity_ZPt100_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"SumZSecondJetRapidity_ZPt100_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DifZSecondJetRapidity_ZPt100_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"AbsZRapidity_ZPt150_Zinc1jet",         0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"AbsFirstJetRapidity_ZPt150_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"SumZFirstJetRapidity_ZPt150_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DifZFirstJetRapidity_ZPt150_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"AbsZRapidity_ZPt150_Zinc2jet",         0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"AbsSecondJetRapidity_ZPt150_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"SumZSecondJetRapidity_ZPt150_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DifZSecondJetRapidity_ZPt150_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"AbsZRapidity_ZPt300_Zinc1jet",         0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"AbsFirstJetRapidity_ZPt300_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"SumZFirstJetRapidity_ZPt300_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DifZFirstJetRapidity_ZPt300_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    {"ZPt_Zinc1jet", 1, 1, 1, 3, 3, 3, 3, 3}, // 1
    //  {"Phistar_Zinc1jet",          1, 1,      1,           3,              3,             3,
    //  3 ,               3    }, // 1
    //  {"Phistar_Zinc1jetM111_130",          1, 1,      1,           3,              3,
    //  3,            3 ,               3    }, // 1
    //  {"Phistar_Zinc1jetM130_3",          1, 1,      1,           3,              3,
    //  3,            3 ,               3    }, // 1
    //  {"ZPt_Zinc1jet",          1, 1,      1,           3,              3,             3,
    //  3 ,               3    }, // 1
    //  {"ZPt_Zinc1jetM111_130",          1, 1,      1,           3,              3,             3,
    //  3 ,               3    }, // 1
    // {"ZPt_Zinc1jetM130_3",          1, 1,      1,           3,              3,             3,
    // 3 ,               3    }, // 1
    //    {"ZPt_Zinc2jet",          0, 1,      1,           3,              3,             3,
    //    3 ,               3    }, // 1
    //  {"VisPt_Zinc0jetQun",     0, 0,      1,           3,              3,             3,
    //  3 ,               3    }, // 1
    //  {"VisPt_Zinc1jetQun",     0, 0,      1,           3,              3,             3,
    //  3 ,               3    }, // 1
    //  {"VisPt_Zinc2jetQun",     0, 0,      1,           3,              3,             3,
    //  3 ,               3    }, // 1
    //  {"VisPt_Zinc3jetQun",     0, 0,      1,           3,              3,             3,
    //  3 ,               3    }, // 1
    //    {"HadRecoil",             0, 1,      1,           3,              3,             3,
    //    3 ,               3    }, // 1
    // PAPER//{"JZB",                   0, 0,      1,           3,              3,             3,
    // 3 ,               3    }, // 1
    // PAPER//{"JZB_ptHigh",                   0, 0,      1,           3,              3,
    // 3,            3 ,               3    }, // 1
    // PAPER//{"JZB_ptLow",                   0, 0,      1,           3,              3,
    // 3,            3 ,               3    }, // 1
    // different JetPt cuts////////
    //    {"AbsZRapidity_FirstJetPt50_Zinc1jet",         0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"AbsFirstJetRapidity_FirstJetPt50_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"SumZFirstJetRapidity_FirstJetPt50_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DifZFirstJetRapidity_FirstJetPt50_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"AbsZRapidity_FirstJetPt80_Zinc1jet",         0, 1,      1,           3,              3,
    //    3,            3 ,               3    }, // 1
    //    {"AbsFirstJetRapidity_FirstJetPt80_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"SumZFirstJetRapidity_FirstJetPt80_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DifZFirstJetRapidity_FirstJetPt80_Zinc1jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    // Set jet rapidity discriminator////
    //    {"AbsZRapidity_DifJetRapidityl2_Zinc2jet",         0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"AbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"SumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"AbsZRapidity_DifJetRapiditys2_Zinc2jet",         0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"AbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"SumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    //    {"DifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",          0, 1,      1,           3,
    //    3,             3,            3 ,               3    }, // 1
    {"ZNGoodJets_Zexc", 0, 0, 1, 7, 4, 3, 3, 3}, // 0
};

const unsigned int NVAROFINTERESTZJETS = sizeof(VAROFINTERESTZJETS) / sizeof(VAROFINTERESTZJETS[0]);

/** Gets the index of a variable for the VAROFINTERESTZJETS array.
 * @param variable: variable name
 * @return the index if the variable was found, -1 otherwise
 */
int findVariable(const TString &variable);

/** Interprets a variable name and decompose it in observable, jet multiplicity,
 * jet multiplicity inclusiveness. For the output parameter pointers, 0 can be passed
 * if the corresponding output is not needed.
 * @param name [in] variable name in the form Observable_ZxxxNjet, with: xxx = inc,exc; N = 0,1,2,..
 * @param observable [out] to store observable name
 * @param jetMultStr [out] to store ZincNjet/ZexcNjet string
 * @param nJets [out] to store the number of jets
 * @param isInc [out] to store a flag which is true
 * @return true in case of success of name interpretation, false otherwise.
 */
bool decomposeVarName(const std::string &name,
                      std::string *observable = 0,
                      std::string *jetMultStr = 0,
                      int *nJets = 0,
                      bool *isInc = 0);
#endif
