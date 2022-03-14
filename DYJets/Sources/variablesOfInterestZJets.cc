#include "variablesOfInterestZJets.h"
#include <cstdlib>
#include <iostream>
#include <regex.h>

int findVariable(const TString &variable)
{
    int index = -1;
    for (unsigned int i = 0; i < NVAROFINTERESTZJETS; ++i) {
        if (VAROFINTERESTZJETS[i].name == variable) {
            index = i;
            break;
        }
    }

    return index;
}

bool decomposeVarName(const std::string &name,
                      std::string *observable,
                      std::string *jetMultStr,
                      int *nJets,
                      bool *isInc)
{

    // variable name example: FirstJetPt_Zinc1jet
    regex_t reg;
    regmatch_t match[6];
    const char *p = "\\(.*\\)_\\(\\(Zinc\\|Zexc\\)\\(\\([[:digit:]]\\+\\)jet\\)\\?\\)";
    int rc = regcomp(&reg, p, 0);
    if (rc != 0) {
        std::cerr << "Bug found in " __FILE__ << ":" << __LINE__
                  << ": invalid regular expression.\n";
        char buffer[256];
        regerror(rc, &reg, buffer, sizeof(buffer));
        buffer[sizeof(buffer) - 1] = 0;
        std::cerr << buffer << ".\n";
        abort();
    }
    rc = regexec(&reg, name.c_str(), sizeof(match) / sizeof(match[0]), match, 0);
    if (rc != 0) {
        std::cerr << "Failed to interpret the variable name '" << name
                  << "'. It does not follow the common pattern '" << p << "'.\n";
    } else {

        if (observable) *observable = name.substr(match[1].rm_so, match[1].rm_eo - match[1].rm_so);
        if (jetMultStr) *jetMultStr = name.substr(match[2].rm_so, match[2].rm_eo - match[2].rm_so);
        if (isInc)
            *isInc = (name.compare(match[3].rm_so, match[3].rm_eo - match[3].rm_so, "Zinc") == 0);
        if (nJets) {
            if (match[4].rm_so < 0)
                *nJets = 0; // Zinc or Zexc without n jet specification => 0 jet
            else
                *nJets = strtol(
                    name.substr(match[5].rm_so, match[5].rm_eo - match[5].rm_so).c_str(), 0, 10);
        }
    }
    regfree(&reg);
    return (rc == 0);
}
