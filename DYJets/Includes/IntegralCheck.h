#include "ConfigVJets.h"
#include "getFilesAndHistogramsZJets.h"
#include "variablesOfInterestZJets.h"

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "printTable.h"
#include <cassert>
#include <string>
//#include "PlotSettings.h"
#include "SectionedConfig.h"
#include "TextTools.h"

/** Computes the integral of the different measured differental cross sections
 * for cross check. A latex table with the result is produced.
 * The list of considered differential cross sections is taken from the file
 * Inclues/variableOfInterest.h
 */

class IntegralCheck
{
  public:
    double integrate(TH1 *h, int istart, const char *opt)
    {
        double r = 0;
        if (istart == 0) {
            r += h->GetBinContent(0);
            ++istart;
        }
        return h->Integral(istart, h->GetNbinsX(), opt) + h->GetBinContent(h->GetNbinsX() + 1);
    }

    void run()
    {
        std::string histoDir = cfg.getS("histoDir");
        std::string unfoldDir = cfg.getS("unfoldDir");
        std::string lepSel = cfg.getS("lepSel");
        int jetEtaMax = cfg.getI("jetEtaMax");
        int jetPtMin = cfg.getI("jetPtMin");
        std::string algo = cfg.getS("algo");

        TString unfCfgFile = cfg.getS("unfConf");
        static SectionedConfig unfCfg;
        unfCfg.read(unfCfgFile, true);

        const char *outDir = "IntegralCheck";
        gSystem->mkdir(outDir);

        TFile *recoFile =
            getFile(histoDir, lepSel, getEnergy(), "Data", jetPtMin, jetEtaMax, "", "0");
        if (!recoFile || recoFile->IsZombie()) {
            std::cerr << "The histo file for Data lepSel =  " << lepSel << " was not found in"
                      << histoDir << "\n";
        }

        // vector index corresponds to the jet multiplicity.
        std::map<std::string, std::vector<double>> recoInt[2]; // 0->incl. 1->excl.
        std::map<std::string, std::vector<double>> unfInt[2];  // 0->incl. 1->excl.
        std::vector<bool> mults[2]; // flags recording  for each jet multiplicity
        //                            if at least one measurement is available
        std::vector<double> refRecoXsec[2];
        std::vector<double> refUnfXsec[2];
        std::map<std::string, std::string> obsLabels;
        // Gets cross sections from jet multipliciy plots
        const char *hnames[] = {"ZNGoodJets_Zinc", "ZNGoodJets_Zexc"};
        if (recoFile) {
            TH1 *hReco;
            for (int iInc = 0; iInc < 2; ++iInc) {
                recoFile->GetObject(hnames[iInc], hReco);
                if (!hReco) {
                    std::cerr << "Histogram " << hnames[iInc] << " was not found in file "
                              << recoFile->GetName() << "\n";
                } else {
                    for (int i = 1; i <= hReco->GetNbinsX(); ++i) {
                        int x = (int)floor(hReco->GetXaxis()->GetBinLowEdge(i)) + 1;
                        if (x >= (int)refRecoXsec[iInc].size()) refRecoXsec[iInc].resize(x + 1);
                        refRecoXsec[iInc][x] = hReco->GetBinContent(i);
                    }
                }
            }
        }

        TFile *unfFile = TFile::Open(
            getUnfoldedFileName(
                unfoldDir, lepSel, hnames[1], algo, jetPtMin, jetEtaMax, "_MGPYTHIA6_", false) +
            ".root");
        if (!unfFile || unfFile->IsZombie()) {
            std::cerr << "Warning. File of the unfolded distribution of " << hnames[1]
                      << "was not found.\n";
        } else {
            TH1 *hUnfDataCentral;
            unfFile->GetObject("UnfDataCentral", hUnfDataCentral);
            if (!hUnfDataCentral) {
                std::cerr << "Histogram "
                          << "UnfDataCentral"
                          << " was not found in file " << unfFile->GetName() << "\n";
            } else {
                int xmin = hUnfDataCentral->GetXaxis()->GetBinLowEdge(1);
                int xmax = hUnfDataCentral->GetXaxis()->GetBinLowEdge(hUnfDataCentral->GetNbinsX());
                refUnfXsec[1] = refUnfXsec[0] = std::vector<double>(xmax - xmin + 1);
                for (int i = 1; i <= hUnfDataCentral->GetNbinsX(); ++i) {
                    int x = (int)floor(hUnfDataCentral->GetXaxis()->GetBinLowEdge(i)) + 1;
                    refUnfXsec[1][x] = hUnfDataCentral->GetBinContent(i);
                }
                double acc = hUnfDataCentral->GetBinContent(hUnfDataCentral->GetNbinsX() + 1);
                for (int x = xmax; x >= xmin; --x) {
                    acc += refUnfXsec[1][x];
                    refUnfXsec[0][x] = acc;
                }
            }
        }
        if (unfFile) delete unfFile;

        for (unsigned i = 0; i < NVAROFINTERESTZJETS; ++i) {
            std::string obs;
            int nJets;
            bool isInc;
            if (!decomposeVarName(VAROFINTERESTZJETS[i].name.Data(), &obs, 0, &nJets, &isInc)) {
                continue;
            }
            obsLabels[obs] = obs;
            TString section =
                TString::Format("%s_%s", lepSel.c_str(), VAROFINTERESTZJETS[i].name.Data());
            int nFirstBinsToSkip = unfCfg.get(section.Data(), "nFirstBinsToSkip", 0);

            // Jet pt distributions starts at a lower value than the threshold used for the
            // per-jet-multiplicity cross sections.
            int istart =
                (VAROFINTERESTZJETS[i].name.Index("JetPt_Zinc") > 0) ? nFirstBinsToSkip + 1 : 0;

            assert(nJets < 100);
            int iInc = isInc ? 0 : 1;
            if (nJets >= (int)mults[iInc].size()) {
                mults[iInc].resize(nJets + 1, false);
                for (auto &&r : recoInt[iInc]) {
                    r.second.resize(nJets + 1, -1);
                }
                for (auto &&r : unfInt[iInc]) {
                    r.second.resize(nJets + 1, -1);
                }
            }
            mults[iInc][nJets] = true;

            // create entries for obs if they don't exist yet:
            recoInt[iInc].insert(std::make_pair(obs, std::vector<double>(mults[iInc].size(), -1)));
            unfInt[iInc].insert(std::make_pair(obs, std::vector<double>(mults[iInc].size(), -1)));

            if (recoInt[iInc][obs][nJets] > -1. || unfInt[iInc][obs][nJets] > -1.) {
                std::cerr << "Fatal Error. Several histograms found for the differential "
                          << " cross-section of " << obs << " for "
                          << (isInc ? "inclusive" : "exclusive") << "jet multiplicity of " << nJets
                          << ".\n";
                abort();
            }

            TH1 *hReco;
            if (recoFile) {
                recoFile->GetObject(VAROFINTERESTZJETS[i].name.Data(), hReco);
                if (!hReco) {
                    std::cerr << "Histogram " << VAROFINTERESTZJETS[i].name
                              << " was not found in file " << recoFile->GetName() << "\n";
                } else {
                    recoInt[iInc][obs][nJets] = integrate(hReco, istart, "");
                }
            }
            TFile *unfFile = TFile::Open(getUnfoldedFileName(unfoldDir,
                                                             lepSel,
                                                             VAROFINTERESTZJETS[i].name,
                                                             algo,
                                                             jetPtMin,
                                                             jetEtaMax,
                                                             "_MGPYTHIA6_",
                                                             false) +
                                         ".root");

            if (!unfFile || unfFile->IsZombie()) {
                std::cerr << "Unfold file not found for " << VAROFINTERESTZJETS[i].name << "\n";
            } else {
                TH1 *hUnfDataCentral;
                unfFile->GetObject("UnfDataCentral", hUnfDataCentral);
                if (!hUnfDataCentral) {
                    std::cerr << "Histogram "
                              << "UnfDataCentral"
                              << " was not found in file " << unfFile->GetName() << "\n";
                } else {
                    TString xtitle = hUnfDataCentral->GetXaxis()->GetTitle();
                    TString title;
                    TString var, varUnit;
                    TString dSigma, dSigmaUnit;
                    bool sepLumiUnc = false;
                    createTitleVariableAnddSigma(VAROFINTERESTZJETS[i].name,
                                                 false,
                                                 xtitle,
                                                 title,
                                                 var,
                                                 varUnit,
                                                 dSigma,
                                                 dSigmaUnit,
                                                 sepLumiUnc);
                    if (var.Length() > 0) obsLabels[obs] = var.Data();
                    obsLabels["JZB_ptLow"] = "JZB, low $p_\\text{T}$";
                    obsLabels["JZB_ptHigh"] = "JZB, high $p_\\text{T}$";
                    obsLabels["JetsHT"] = "$H_{\\text{T}}$";

                    //	  std::cout << "xtitle = " << xtitle << "\n"
                    //		    << "title = " << title << "\n"
                    //		    << "var = " << var << "\n"
                    //		    << "dSigma = " << dSigma << "\n";

                    unfInt[iInc][obs][nJets] = integrate(hUnfDataCentral, istart, "width");
                }
            }
            if (unfFile) delete unfFile;
        }
        if (recoFile) delete recoFile;

        const char *multLabel[2] = {"= ", "$\\ge$ "};

        int ixsec = -1;
        const char *refLabels[] = {"N_{\\text{tot}}", "\\sigma_{\\text{tot}}"};
        for (auto xsecs : {recoInt, unfInt}) {
            ++ixsec;
            for (int iInc : {0, 1}) {
                std::vector<double> *ref[2] = {refRecoXsec + iInc, refUnfXsec + iInc};
                std::vector<std::string> lineHeaders(xsecs[iInc].size() + 1);
                std::vector<std::vector<std::string>> vals(xsecs[iInc].size() + 1);
                std::vector<std::string> colHeaders;
                colHeaders.push_back("Observable");
                for (unsigned nJets = 0; nJets < mults[iInc].size(); ++nJets) {
                    if (mults[iInc][nJets]) {
                        colHeaders.push_back(
                            TString::Format(
                                "%s %d jet%s", multLabel[iInc], nJets, (nJets > 0 ? "" : "s"))
                                .Data());
                    }
                } // next nJets

                lineHeaders[0] = refLabels[ixsec];
                vals[0] = std::vector<std::string>(mults[iInc].size(), "-");
                for (unsigned nJets = 0; nJets < mults[iInc].size(); ++nJets) {
                    if (nJets < ref[ixsec]->size()) {
                        vals[0][nJets] = TString::Format("%#.3g", (*ref[ixsec])[nJets]).Data();
                    }
                }
                int irow = 1;
                for (auto rec : xsecs[iInc]) {
                    assert(rec.second.size() == mults[iInc].size());
                    lineHeaders[irow] = obsLabels[rec.first];
                    int ival = -1;
                    vals[irow] = std::vector<std::string>(colHeaders.size() - 1, "-");
                    for (unsigned nJets = 0; nJets < mults[iInc].size(); ++nJets) {
                        if (!mults[iInc][nJets]) continue;
                        ++ival;
                        if (rec.second[nJets] > -1.) {
                            vals[irow][ival] = TString::Format("%#.3g", rec.second[nJets]).Data();
                        }
                    }
                    ++irow;
                } // next row

                const char *outFileNames[2][2] = {
                    {"RecoXsecIntegralTable_Zinc.tex", "RecoXsecIntegralTable_Zexc.tex"},
                    {"UnfXsecIntegralTable_Zinc.tex", "UnfXsecIntegralTable_Zexc.tex"}};

                std::ofstream f(TString::Format("%s/%s", outDir, outFileNames[ixsec][iInc]).Data());
                const char *captions[2][2] = {
                    {"Check of the integrals of the reconstruction-level histograms of the "
                     "measured observables inclusive in jet multiplicity.",
                     "Check of the integrals of the reconstruction-level  histograms of the "
                     "measured observables excluive in jetx multiplicity."},
                    {"Check of the integrals of the differential cross sections inclusive in jet "
                     "multiplicity.",
                     "Check of the integrals of the differential cross sections exclusive in jet "
                     "multiplicity."}};
                const char *labels[2][2] = {{"tab:tot_exc_xsec_reco", "tab:tot_inc_xsec_reco"},
                                            {"tab:tot_exc_xsec_unf", "tab:tot_inc_xsec_unf"}};
                printTable(
                    f, colHeaders, lineHeaders, vals, captions[ixsec][iInc], labels[ixsec][iInc]);
            }
        }
    }
};
