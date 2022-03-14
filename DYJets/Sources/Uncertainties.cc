#include "Uncertainties.h"

#include <assert.h>

const char *uncShortNames[kUncCnt] = {
    "Stat", "UnfStat", "JES", "PU", "JER", "Xsec", "LES", "LER", "Lumi", "SF", "AltUnf", "Tot"};

TString covName(unsigned i, const char *lepSel)
{
    assert(i < kUncCnt);
    return TString(lepSel) + "Cov" + uncShortNames[i] + lepSel;
}
