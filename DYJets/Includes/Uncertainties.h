#ifndef UNCERTAINTIES_H
#define UNCERTAINTIES_H

#include <TString.h>

enum VarIndices {
    kCentral = 0,
    kJESup = 1,
    kJESdwn = 2,
    kPUup = 3,
    kPUdwn = 4,
    kJERup = 5,
    kJERdwn = 6,
    kXsecUp = 7,
    kXsecDwn = 8,
    kLESup = 9,
    kLESdwn = 10,
    kLERup = 11,
    kLERdwn = 12,
    kLumiUp = 13,
    kLumiDwn = 14,
    kSFup = 15,
    kSFdwn = 16,
    kAltUnf = 17
};
enum UncIndices {
    kStat = 0,
    kUnfStat = 1,
    kJES = 2,
    kPU = 3,
    kJER = 4,
    kXsec = 5,
    kLES = 6,
    kLER = 7,
    kLumi = 8,
    kSF = 9,
    kUnfSys = 10,
    kTotSys = 11,
    kUncCnt = 12 // number of uncertainty indices
};

extern const char *uncShortNames[kUncCnt]; // defined in Uncertainties.cc

TString covName(unsigned i, const char *lepSel = "");

#endif // UNCERTAINTIES_H not defined
