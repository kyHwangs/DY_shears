#ifndef RESULTTABLES_H
#define RESULTTABLES_H

#include "TH2.h"
#include "TString.h"
#include <vector>

// Version used by Unfolding
// void createTable(TString outputFileName, TString lepSel, TString variable,
//		 bool doNormalized, TH1 *hUnfData, TH2 *hCov[]);

// Version used by Combination
void createTable(TString outputFileName,
                 TString lepSel,
                 TString variable,
                 bool doNormalized,
                 TH1 *hMeas,
                 std::vector<TH2 *> &covuxaxb,
                 bool withLERS = false,
                 bool forceNoLumiSep = false);

#endif // RESULTTABLES_H not defined
