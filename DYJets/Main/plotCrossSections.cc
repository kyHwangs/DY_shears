#include "ArgParser.h"
#include "ConfigVJets.h"
#include "PlotSettings.h"

int main(int argc, char **argv)
{
    TString variable("");
    TString lepSel;
    TString ref = cfg.getS("ref");
    for (int i = 1; i < argc; ++i) {
        TString currentArg = argv[i];
        //--- possible options ---
        if (currentArg.BeginsWith("variable=")) {
            getArg(currentArg, variable);
        }
        if (currentArg.BeginsWith("lepSel=")) {
            getArg(currentArg, lepSel);
            cfg.set("lepSel", lepSel);
        }
        if (currentArg.BeginsWith("ref=")) {
            getArg(currentArg, ref);
            cfg.set("ref", ref);
        }
    }
    makeCrossSectionPlot(variable.Length() > 0 ? variable.Data() : 0, ref);
    return 0;
}
