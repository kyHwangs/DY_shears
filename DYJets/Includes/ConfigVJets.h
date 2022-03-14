#ifndef CONFIGVJETS_H
#define CONFIGVJETS_H

#include "Config.h"
#include "SectionedConfig.h"
#include "TObject.h"
#include "TString.h"
#include "stdlib.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

class ConfigVJets : public Config
{
  public:
    ConfigVJets(const char *filename = 0);

    bool read(const char *filename = 0);

    template <typename T>
    bool getUnf1(const char *lepSel, const char *variable, const char *param, T &val) const;

    template <typename T>
    T getUnf(const char *lepSel,
             const char *variable,
             const char *param,
             const T &defaultVal) const;

    template <typename T>
    T getUnf(const char *lepSel, const char *variable, const char *param) const;

    TString getUnfSec(const char *lepSel, const char *variable) const
    {
        return TString::Format("%s_%s", lepSel, variable);
    }

  private:
    SectionedConfig unfCfg;
};

extern ConfigVJets cfg;

template <typename T>
bool ConfigVJets::getUnf1(const char *lepSel, const char *variable, const char *param, T &val) const
{
    return unfCfg.get1<T>(getUnfSec(lepSel, variable).Data(), param, val);
}

template <typename T>
T ConfigVJets::getUnf(const char *lepSel,
                      const char *variable,
                      const char *param,
                      const T &defaultVal) const
{
    return unfCfg.get<T>(getUnfSec(lepSel, variable).Data(), param, defaultVal);
}

template <typename T>
T ConfigVJets::getUnf(const char *lepSel, const char *variable, const char *param) const
{
    return unfCfg.get<T>(getUnfSec(lepSel, variable).Data(), param);
}

#endif // CONFIGVJETS_H not defined
