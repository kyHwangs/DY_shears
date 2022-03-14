#include "ArgParser.h"
#include "ConfigVJets.h"
#include "UnfoldingZJets.h"
#include <TString.h>
#include <iostream>

int main(int argc, char **argv)
{
    //--- Loads configuration -----------------------------------------------------

    TString histoDir = cfg.getS("histoDir");
    TString unfoldDir = cfg.getS("unfoldDir");
    TString lepSel = cfg.getS("lepSel");
    TString algo = cfg.getS("algo");
    int jetPtMin = cfg.getI("jetPtMin");
    int jetEtaMax = cfg.getI("jetEtaMax");
    int whichSyst = cfg.getI("whichSyst");
    TString generator1 = cfg.getS("generator1", "mgpythia8");
    TString generator2 = cfg.getS("generator2", "");

    TString variable = "";
    // TString variable = "ZNGoodJets_Zinc";
    // TString variable = "ZPt_Zinc0jet";
    // TString variable =
    // "FirstJetPt_Zinc1jet,FirstJetAbsRapidity_Zinc1jet,SecondJetPt_Zinc2jet,SecondJetAbsRapidity_Zinc2jet,ThirdJetPt_Zinc3jet,ThirdJetAbsRapidity_Zinc3jet,FourthJetPt_Zinc4jet,FourthJetAbsRapidity_Zinc4jet";

    bool doNormalized(false);

    int nIters = 0;

    //-----------------------------------------------------------------------------

    //--- Parse the arguments -----------------------------------------------------
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            TString currentArg = argv[i];
            //--- possible options ---
            if (currentArg.BeginsWith("histoDir=")) {
                getArg(currentArg, histoDir);
                cfg.set("histoDir", histoDir);
            } else if (currentArg.BeginsWith("unfoldDir=")) {
                getArg(currentArg, unfoldDir);
                cfg.set("unfoldDir", unfoldDir);
            } else if (currentArg.BeginsWith("lepSel=")) {
                getArg(currentArg, lepSel);
                cfg.set("lepSel", lepSel);
            } else if (currentArg.BeginsWith("algo=")) {
                getArg(currentArg, algo);
                cfg.set("algo", algo);
            }

            else if (currentArg.BeginsWith("generator1=")) {
                getArg(currentArg, generator1);
                cfg.set("generator1", generator1);
            }

            /*          else if (currentArg.BeginsWith("generator2=")) {
                         getArg(currentArg, generator2);
                         cfg.set("generator2", generator2);
                     }
         */
            else if (currentArg.BeginsWith("jetPtMin=")) {
                getArg(currentArg, jetPtMin);
                cfg.set("jetPtMin", jetPtMin);
            } else if (currentArg.BeginsWith("jetEtaMax=")) {
                getArg(currentArg, jetEtaMax);
                cfg.set("jetEtaMax", jetEtaMax);
            } else if (currentArg.BeginsWith("variable=")) {
                getArg(currentArg, variable);
                cfg.set("variable", variable);
            } else if (currentArg.BeginsWith("doNormalized=")) {
                getArg(currentArg, doNormalized);
                cfg.set("doNormalized", doNormalized);
            } else if (currentArg.BeginsWith("whichSyst=")) {
                getArg(currentArg, whichSyst);
                cfg.set("whichSyst", whichSyst);
            } else if (currentArg.BeginsWith("nIters=")) {
                getArg(currentArg, nIters);
                cfg.set("minIter", nIters);
                cfg.set("maxIter", nIters);
            }
            //--- asking for help ---
            else if (currentArg.Contains("help") || currentArg.BeginsWith("-h")) {
                std::cout << "\nUsage: \n\t./runUnfolding [lepSel=(DMu, DE)] [algo=(Bayes, SVD)] "
                             "[jetPtMin=(int)] [jetEtaMax=(int*10)] [histoDir=(path)] "
                             "[unfoldDir=(path)] [variable=(var1),(var2),...] [doNormalized=(0, "
                             "1)] [whichSyst=(-1)] [--help]"
                          << std::endl;
                std::cout << "\neg: ./runUnfolding lepSel=DMu jetEtaMax=24" << std::endl;
                std::cout << "\nunspecified options will be read from vjets.cfg\n" << std::endl;
                return 0;
            }
            //--- bad option ---
            else {
                std::cerr << "Warning: unknown option \"" << currentArg << "\"" << std::endl;
                std::cerr
                    << "Please issue ./runUnfolding --help for more information on possible options"
                    << std::endl;
                return 0;
            }
        }
    }

    if (!histoDir.EndsWith("/")) histoDir += "/";
    if (!unfoldDir.EndsWith("/")) unfoldDir += "/";

    TString unfCfgFile = cfg.getS("unfConf");
    static SectionedConfig unfCfg;
    unfCfg.read(unfCfgFile);

    // UnfoldingZJets(lepSel, algo, histoDir, unfoldDir, jetPtMin, jetEtaMax, generator1,
    // generator2, variable, doNormalized);
    // UnfoldingZJets(lepSel, algo, histoDir, unfoldDir, jetPtMin, jetEtaMax, variable,
    // doNormalized, whichSyst);

    // support lepSel list DE,DMu:
    TString lepSel_;
    Ssiz_t from0 = 0;
    while (lepSel.Tokenize(lepSel_, from0, "[ \\t]*[, ][ \\t]*")) {
        if (variable.Length() > 0) {
            TString variable_;
            Ssiz_t from1 = 0;
            // support for variable list like JetPt_Zinc2jet,JetAbsRapidity_Zinc2jet
            while (variable.Tokenize(variable_, from1, "[ \\t]*[, ][ \\t]*")) {
                UnfoldingZJets(unfCfg,
                               lepSel_,
                               algo,
                               histoDir,
                               unfoldDir,
                               jetPtMin,
                               jetEtaMax,
                               variable_,
                               doNormalized,
                               whichSyst);
            }
        } else {
            UnfoldingZJets(unfCfg,
                           lepSel_,
                           algo,
                           histoDir,
                           unfoldDir,
                           jetPtMin,
                           jetEtaMax,
                           "",
                           doNormalized,
                           whichSyst);
        }
    }

    return 0;
}
