#include "variablesOfInterestZJets.h"
#include <iostream>
#include <string>
#include <vector>

/** Test decomposeVarName function and compliance of variable names
 * to common convention.
 */

int main()
{
    std::vector<std::string> invalidVars;
    for (unsigned ivar = 0; ivar < NVAROFINTERESTZJETS; ++ivar) {
        std::string vname = VAROFINTERESTZJETS[ivar].name.Data();
        std::string obs;
        std::string jetMultStr;
        int nJets;
        bool isInc;
        std::cout << vname << ":\n\t";
        if (decomposeVarName(vname, &obs, &jetMultStr, &nJets, &isInc)) {
            std::cout << "Observable: " << obs << "\n\tMultiplicity label: " << jetMultStr
                      << "\n\tMultiplicity: " << (isInc ? " >= " : " = ") << nJets << " jet(s)\n";
        } else {
            invalidVars.push_back(vname);
        }
        std::cout << "\n";
    }
    if (invalidVars.size() > 0) {
        std::cout << "Following variables found in Includes/variablesOfInterestZjets.h do not "
                     " follow the expected naming scheme:\n\n";
        for (unsigned i = 0; i < invalidVars.size(); ++i)
            std::cout << "\t" << invalidVars[i] << "\n";
        std::cout << "\n";
    }
}
