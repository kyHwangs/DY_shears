#include "ArgParser.h"
#include "Combination.h"
#include "ConfigVJets.h"
#include <iostream>

int main(int argc, char **argv)
{

    //--- Load configuration ---
    //-----------------------------------------------------------------------

    //--- Settings ---
    // doWhat = "DATA", "BACKGROUND", "TAU", "DYJETS",
    //          "WJETS", "ALL", "PDF", "SHERPA"

    TString unfoldDir = cfg.getS("unfoldDir");
    TString combDir = cfg.getS("combDir");
    TString histoDir = cfg.getS("histoDir");
    TString algo = cfg.getS("algo");
    int jetPtMin = cfg.getI("jetPtMin", 30);
    int jetEtaMax = cfg.getI("jetEtaMax", 24);
    bool diagXChanCov = cfg.getB("diagXChanCov", true);
    bool fullXChanCov = cfg.getB("fullXChanCov", true);
    bool fullSChanCov = cfg.getB("fullSChanCov", true);
    bool modifiedSWA = cfg.getB("modifiedSWA", true);

    TString variable = "";
    TString genList;

    bool doNormalized(false);
    // bool doNormband(false);

    //--- Parse the arguments -----------------------------------------------------
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            TString currentArg = argv[i];
            //--- possible options ---
            if (currentArg.BeginsWith("unfoldDir=")) {
                getArg(currentArg, unfoldDir);
                cfg.set("unfoldDir", unfoldDir);
            } else if (currentArg.BeginsWith("combDir=")) {
                getArg(currentArg, combDir);
                cfg.set("combDir", combDir);
            } else if (currentArg.BeginsWith("histoDir")) {
                getArg(currentArg, histoDir);
                cfg.set("histoDir", histoDir);
            } else if (currentArg.BeginsWith("algo=")) {
                getArg(currentArg, algo);
                cfg.set("algo", algo);
            } else if (currentArg.BeginsWith("predictions=")) {
                getArg(currentArg, genList);
                cfg.set("predictions", genList);
            } else if (currentArg.BeginsWith("jetPtMin=")) {
                getArg(currentArg, jetPtMin);
                cfg.set("jetPtMin", jetPtMin);
            } else if (currentArg.BeginsWith("jetEtaMax=")) {
                getArg(currentArg, jetEtaMax);
                cfg.set("jetEtaMax", jetEtaMax);
            } else if (currentArg.BeginsWith("diagXChanCov=")) {
                getArg(currentArg, diagXChanCov);
                cfg.set("diagXChanCov", diagXChanCov);
            } else if (currentArg.BeginsWith("fullXChanCov=")) {
                getArg(currentArg, fullXChanCov);
                cfg.set("fullXChanCov", fullXChanCov);
            } else if (currentArg.BeginsWith("fullSChanCov=")) {
                getArg(currentArg, fullSChanCov);
                cfg.set("fullSChanCov", fullSChanCov);
            } else if (currentArg.BeginsWith("modifiedSWA=")) {
                getArg(currentArg, modifiedSWA);
                cfg.set("modifiedSWA", modifiedSWA);
            } else if (currentArg.BeginsWith("variable=")) {
                getArg(currentArg, variable);
                cfg.set("variable", variable);
            } else if (currentArg.BeginsWith("doNormalized=")) {
                getArg(currentArg, doNormalized);
                cfg.set("doNormalized", doNormalized);
            }
            /*      else if (currentArg.BeginsWith("doNormband=")) {
                      getArg(currentArg, doNormband);
                  }*/
            //--- asking for help ---
            else if (currentArg.Contains("help") || currentArg.BeginsWith("-h")) {
                std::cout << "\nUsage: ./runCombination [unfoldDir=(path)] [combDir=(path)] "
                             "[histoDir=(path)] [algo=(Bayes, SVD)] [jetPtMin=(int)] "
                             "[jetEtaMax=(int*10)]";
                std::cout << "[diagXChanCov=(1,0)] [fullXChanCov=(1,0)] [fullSChanCov=(1,0)] "
                             "[modifiedSWA=(1,0)] [variable=(variableName)] [doNormalized=(0, 1)]  "
                             "[predictions=(comma-separated list)] [--help]"
                          << std::endl;
                std::cout << "eg: ./runCombination fullXChanCov=0 jetEtaMax=24" << std::endl;
                std::cout << "unspecified options will be read from vjets.cfg\n" << std::endl;
                return 0;
            }
            //--- bad option ---
            else {
                std::cerr << "Warning: unknown option \"" << currentArg << "\"" << std::endl;
                std::cerr << "Please issue ./runCombination --help for more information on "
                             "possible options"
                          << std::endl;
                return 0;
            }
        }
    }

    if (!unfoldDir.EndsWith("/")) unfoldDir += "/";
    if (!combDir.EndsWith("/")) combDir += "/";

    Combination(unfoldDir,
                combDir,
                algo,
                jetPtMin,
                jetEtaMax,
                diagXChanCov,
                fullXChanCov,
                fullSChanCov,
                modifiedSWA,
                variable,
                doNormalized);
    return 0;
}
