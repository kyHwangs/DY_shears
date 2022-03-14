#include "ConfigVJets.h"
#include "TObject.h"
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

ConfigVJets cfg = ConfigVJets();

bool ConfigVJets::read(const char *filename)
{
    bool rc = Config::read(filename);
    TString unfCfgFile = getS("unfConf");
    if (unfCfgFile.Length() > 0) rc &= unfCfg.read(unfCfgFile);
    return rc;
}

ConfigVJets::ConfigVJets(const char *filename) : Config(filename)
{
    TString unfCfgFile = getS("unfConf");
    if (unfCfgFile.Length() > 0) unfCfg.read(unfCfgFile);
}
