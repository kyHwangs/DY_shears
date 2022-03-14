#include "Combination.h"
#include "BLUEMeth.h"
#include "ConfigVJets.h"
#include "PlotSettings.h"
#include "ResultTables.h"
#include "SectionedConfig.h"
#include "Uncertainties.h"
#include "assert.h"
#include "functions.h"
#include "getFilesAndHistogramsZJets.h"
#include "variablesOfInterestZJets.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include <fstream>
#include <iostream>
#include <vector>

extern ConfigVJets cfg; // defined in runCombination.cc

using namespace std;

void createInclusivePlots(bool doNormalized,
                          double lumi,
                          TString outputFileDir,
                          TString outputFileName,
                          TH1 *hUnfData,
                          vector<TH2 *> hCov,
                          TH2 *hCovSyst,
                          const std::vector<std::string> &predictions);

void Combination(TString unfoldDir,
                 TString combDir,
                 TString algo,
                 int jetPtMin,
                 int jetEtaMax,
                 bool diagXChanCov,
                 bool fullXChanCov,
                 bool fullSChanCov,
                 bool modifiedSWA,
                 TString variable,
                 bool doNormalized)
{
    //--- create output directory if does not exist ---
    system("mkdir -p " + combDir);

    int start = 0;
    int end = NVAROFINTERESTZJETS;

    std::vector<std::string> predictions = cfg.getVS("predictions");
    int verbosity = cfg.getI("verbosity");

    if (variable != "") {
        start = findVariable(variable);
        if (start >= 0) {
            end = start + 1;
        } else {
            std::cerr << "\nError: variable " << variable << " is not interesting." << std::endl;
            std::cerr << "See below the list of interesting variables:" << std::endl;
            for (unsigned int i = 0; i < NVAROFINTERESTZJETS; ++i) {
                std::cerr << "\t" << i << ": " << VAROFINTERESTZJETS[i].name << "\n" << std::endl;
            }
            return;
        }
    }

    //--- loop over the variable of interest ---
    for (int i = start; i < end; ++i) {
        variable = VAROFINTERESTZJETS[i].name;

        //--- fetch the electron and muon unfolded files ---
        TString commonName = "_unfolded_" + variable + "_" + algo;
        commonName += "_JetPtMin_";
        commonName += jetPtMin;
        commonName += "_JetEtaMax_";
        commonName += jetEtaMax;
        commonName += "_MGPYTHIA6_"; //+ gen1 + "_" + gen2;
        // TString commonName2;
        // if(doNormband) commonName2 = commonName + "_normalized.root";
        commonName += doNormalized ? "_normalized" : "";
        commonName += ".root";
        TFile *fDE = new TFile(unfoldDir + "DE" + commonName);
        if (!fDE->IsOpen()) {
            std::cerr << "\nError: file " << unfoldDir + "DE" + commonName << " does not exist."
                      << std::endl;
            std::cerr
                << "       You can create it using\n\t ./runUnfoldingZJets lepSel=DE variable="
                << variable << std::endl;
            std::cerr << "Skipping variable " << variable << ".\n" << std::endl;
            continue;
        }

        TFile *fDMu = new TFile(unfoldDir + "DMu" + commonName);
        if (!fDMu->IsOpen()) {
            std::cerr << "\nError: file " << unfoldDir + "DMu" + commonName << " does not exist."
                      << std::endl;
            std::cerr
                << "       You can create it using\n\t ./runUnfoldingZJets lepSel=DMu variable="
                << variable << std::endl;
            std::cerr << "Skipping variable " << variable << ".\n" << std::endl;
            continue;
        }
        //---------------------------------------------------------------------
        //--- fetch the cross section histogram and the covariance matrices ---
        fDE->cd();
        TH1 *hUnfDE = (TH1 *)fDE->Get("UnfDataCentral");
        //        TH1 *hMadGenDE = (TH1*) fDE->Get("hMadGenDYJetsCrossSection");
        TH1 *hLumiDE = (TH1 *)fDE->Get("Lumi");
        if (!hLumiDE) {
            std::cerr << "Error. Luminosity histogram was not found in the file " << fDE->GetName()
                      << ".\n";
            continue;
        }
        double lumiDE = hLumiDE->GetBinContent(1);

        fDMu->cd();
        TH1 *hUnfDMu = (TH1 *)fDMu->Get("UnfDataCentral");
        TH1 *hLumiDMu = (TH1 *)fDMu->Get("Lumi");
        if (!hLumiDMu) {
            std::cerr << "Error. Luminosity histogram was not found in the file " << fDMu->GetName()
                      << ".\n";
            continue;
        }
        double lumiDMu = hLumiDMu->GetBinContent(1);

        double lumi = 0.5 * (lumiDE + lumiDMu);

        //---------------------------------------------------------------------
        TString genList;
        for (auto g : predictions) {
            genList += TString("_") + g;
        }
        genList.ReplaceAll("DYJets_", "");
        genList.ReplaceAll("UNFOLDING", "FXFX");

        //--- create the output root file ---
        TString outputFileName = variable + "_" + algo;
        outputFileName += "_diagXChanCov_";
        outputFileName += (int)diagXChanCov;
        outputFileName += "_fullXChanCov_";
        outputFileName += (int)fullXChanCov;
        outputFileName += "_fullSChanCov_";
        outputFileName += (int)fullSChanCov;
        outputFileName += "_modifiedSWA_";
        outputFileName += (int)modifiedSWA;
        outputFileName += "_JetPtMin_";
        outputFileName += jetPtMin;
        outputFileName += "_JetEtaMax_";
        outputFileName += jetEtaMax;
        // outputFileName += genList;
        outputFileName += doNormalized ? "_normalized" : "";
        TString outputFilePath = combDir + outputFileName;

        TFile *outputRootFile = new TFile(outputFilePath + ".root", "RECREATE");

        //--- fill in the vector of measurements ---
        vector<TH1 *> measurements{hUnfDE, hUnfDMu};

        //--- fill in the vector of vector of covariances ---
        vector<vector<TH2 *>> covariances(2, std::vector<TH2 *>(kUncCnt, 0));
        TFile *f[2] = {fDE, fDMu};
        const char *lepLabels[2] = {"DE", "DMu"};

        for (unsigned ilep = 0; ilep < 2; ++ilep) {
            f[ilep]->cd(); //<-- Is this line required?
            for (unsigned iunc = 0; iunc < kUncCnt; ++iunc) {
                covariances[ilep][iunc] = (TH2 *)f[ilep]->Get(covName(iunc));
                if (!covariances[ilep][iunc]) {
                    std::cerr << "Error. TH2 " << covName(iunc) << " was not found in file "
                              << f[ilep]->GetName() << ". This uncertainty will be missing.\n";
                }
            }
        }
        //--- create objects to be filled with output of combination ---
        vector<TH2 *> covuxaxb(covariances[0].size(), NULL); // each covariance matrix
        assert(covuxaxb.size() == kUncCnt);
        TH2 *covxaxb = NULL;      // total covariance matrix
        TH1 *hTotComUnc = NULL;   // total uncertainty if combined cross section
        TH1 *hCombination = NULL; // combined cross section
        //---------------------------------------------------------------------

        //--- create the BLUEMeth object to compute the covariance ---
        BLUEMeth *blueXSec = new BLUEMeth(measurements, covariances, variable);
        // BLUEMeth* blueXSecNorm;
        //---------------------------------------------------------------------

        //--- set up how you want to combine the channels and do the combination ---
        // one could rerun first do the combination with simple weighted average
        //
        // hCombination = blue->GetCombination(false, false, false, covuxaxb, covxaxb);
        //
        // and then do it with full covariance to get the proper error without fetching
        // the output of the measurement:
        //
        // blue->GetCombination(true, true, true, covuxaxb, covxaxb);
        //
        hCombination = blueXSec->GetCombination(
            diagXChanCov, fullXChanCov, fullSChanCov, modifiedSWA, covuxaxb, covxaxb);

        covuxaxb[kTotSys] = (TH2 *)covxaxb->Clone();
        covuxaxb[kTotSys]->Reset();
        for (int i = kStat + 1; i < kTotSys; ++i) {
            covuxaxb[kTotSys]->Add(covuxaxb[i]);
        }
        covuxaxb[kTotSys]->SetName(covName(kTotSys, "Comb"));

        if (verbosity > 0) {
            int nbins = covxaxb->GetNbinsX();
            for (int i = 1; i <= nbins; i++) {
                std::cout << i << "  " << sqrt(covxaxb->GetBinContent(i, i)) << std::endl;
            }
        }
        //---------------------------------------------------------------------
        TString unfCfgFile = cfg.getS("unfConf");
        static SectionedConfig unfCfg;
        unfCfg.read(unfCfgFile, true);

        TString sectionDE = TString::Format("DE_%s", variable.Data());
        TString sectionDMu = TString::Format("DMu_%s", variable.Data());
        int nFirstBinsToSkip = unfCfg.get(sectionDE.Data(), "nFirstBinsToSkip", 0);
        nFirstBinsToSkip =
            std::max(nFirstBinsToSkip, unfCfg.get(sectionDMu.Data(), "nFirstBinsToSkip", 0));
        int nLastBinsToSkip = unfCfg.get(sectionDE.Data(), "nLastBinsToSkip", 0);
        nLastBinsToSkip =
            std::max(nLastBinsToSkip, unfCfg.get(sectionDMu.Data(), "nLastBinsToSkip", 0));

        TCanvas *crossSectionPlot;
        int upperBin = hCombination->GetXaxis()->GetNbins() - nLastBinsToSkip;
        int lowerBin = 1 + nFirstBinsToSkip;

        int verbosity = cfg.getI("verbosity");

        hCombination->GetXaxis()->SetRange(lowerBin, upperBin);
        crossSectionPlot = makeCrossSectionPlot(
            "", lumi, variable, doNormalized, hCombination, covuxaxb[kTotSys], predictions);
        crossSectionPlot->Draw();
        saveCanvas(crossSectionPlot, combDir, outputFileName);

        //--- print out the combined cross section measurement and fill the total uncertainty ---
        double tempunc = 0;
        hTotComUnc = (TH1 *)hCombination->Clone();
        for (int i = 1; i <= hCombination->GetNbinsX(); ++i) {
            if (verbosity > 0) std::cout << hCombination->GetBinContent(i) << std::endl;
            tempunc =
                sqrt(covuxaxb[0]->GetBinContent(i, i) + covuxaxb[kTotSys]->GetBinContent(i, i));
            hTotComUnc->SetBinContent(i, tempunc);
        }
        //--- print out break down of errors ---
        if (verbosity > 0) {
            for (int i = 2; i <= 11; ++i) {
                std::cout << hCombination->GetBinContent(i);
                for (unsigned j = 0; j < covuxaxb.size(); ++j) {
                    if (!covuxaxb[j]) continue;
                    std::cout << " +/- "
                              << sqrt(covuxaxb[j]->GetBinContent(i, i)) * 100. /
                                     hCombination->GetBinContent(i)
                              << "%";
                }
                std::cout << std::endl;
            }
        }

        createTable(outputFilePath + "_withLERS",
                    "",
                    variable,
                    doNormalized,
                    hCombination,
                    covuxaxb,
                    true,
                    true);

        createTable(
            outputFilePath, "", variable, doNormalized, hCombination, covuxaxb, false, true);

        if (variable.Index("ZNGoodJets_Zexc") >= 0) {
            createInclusivePlots(doNormalized,
                                 lumi,
                                 combDir,
                                 outputFileName,
                                 hCombination,
                                 covuxaxb,
                                 covuxaxb[kTotSys],
                                 predictions);
        }

        //--- save results and inputs to root file ---
        outputRootFile->cd();
        crossSectionPlot->Write();
        hCombination->Write("CombDataCentral");
        hTotComUnc->Write("CombTotUnc");
        covxaxb->Write("CombCovTot");
        covuxaxb[kTotSys]->Write("CombCovTotSyst");
        for (unsigned int i = 0; i < covuxaxb.size(); ++i) {
            if (covuxaxb[i]) covuxaxb[i]->Write();
        }

        hUnfDE->Write("DEUnfDataCentral");
        hUnfDMu->Write("DMuUnfDataCentral");

        for (unsigned ilep = 0; ilep < 2; ++ilep) {
            for (unsigned iunc = 0; iunc < kUncCnt; ++iunc) {
                if (covariances[ilep][iunc])
                    covariances[ilep][iunc]->Write(covName(iunc, lepLabels[ilep]));
            }
        }
        //--- Close all files ---
        outputRootFile->Close();
        fDE->Close();
        fDMu->Close();

        // if (end == start + 1) system("display " + outputFilePath + ".png &");
        // if (end == start + 1 && variable == "ZNGoodJets_Zexc") system("display " +
        // outputFilePath.ReplaceAll("ZNGoodJets_Zexc", "ZNGoodJets_Zinc") + ".png &");
    }
}

void createInclusivePlots(bool doNormalized,
                          double lumi,
                          TString outputFileDir,
                          TString outputFileName,
                          TH1 *hUnfData,
                          vector<TH2 *> hCov,
                          TH2 *hCovSyst,
                          const std::vector<std::string> &predictions)
{
    TH1 *hInc = (TH1 *)hUnfData->Clone("ZNGoodJets_Zinc");

    std::vector<TH1 *> hGens = getGenHistos(predictions, "", "ZNGoodJets_Zinc", true, true);
    hGens.resize(3, 0);

    TH1D *hIncMad = hGens[0] ? (TH1D *)hGens[0]->Clone("ZNGoodJets_Zinc_Mad") : 0;
    TH1D *hIncShe = hGens[1] ? (TH1D *)hGens[1]->Clone("ZNGoodJets_Zinc_She") : 0;
    TH1D *hIncPow = hGens[2] ? (TH1D *)hGens[2]->Clone("ZNGoodJets_Zinc_Pow") : 0;

    //    const int kTot = 11;
    if (hCov.size() != kUncCnt) {
        std::cerr << "Bug found in " << __FILE__ << ":" << __LINE__
                  << ": hCov vector has an unexpected size: " << hCov.size() << " instead of "
                  << kUncCnt << "! Aborts.\n";
        abort();
    }
    hCov.push_back(hCovSyst);

    std::vector<TH2 *> hCovInc(kUncCnt, 0);
    for (unsigned i = 0; i < kUncCnt; ++i) {
        if (hCov[i]) hCovInc[i] = (TH2 *)hCov[i]->Clone();
    }

    int nBins = hInc->GetNbinsX();
    for (int i = 1; i <= nBins; i++) {
        double binSum = 0;
        double binSumMad = 0;
        double binSumShe = 0;
        double binSumPow = 0;
        double binStatError2 = 0;
        double binStatMadError2 = 0;
        double binStatSheError2 = 0;
        double binStatPowError2 = 0;
        for (int j = i; j <= nBins; j++) {
            binSum += hInc->GetBinContent(j);
            if (hIncMad) binSumMad += hIncMad->GetBinContent(j);
            if (hIncShe) binSumShe += hIncShe->GetBinContent(j);
            if (hIncPow) binSumPow += hIncPow->GetBinContent(j);
            binStatError2 += pow(hInc->GetBinError(j), 2);
            if (hIncMad) binStatMadError2 += pow(hIncMad->GetBinError(j), 2);
            if (hIncShe) binStatSheError2 += pow(hIncShe->GetBinError(j), 2);
            if (hIncPow) binStatPowError2 += pow(hIncPow->GetBinError(j), 2);
        }
        hInc->SetBinContent(i, binSum);
        if (hIncMad) hIncMad->SetBinContent(i, binSumMad);
        if (hIncShe) hIncShe->SetBinContent(i, binSumShe);
        if (hIncPow) hIncPow->SetBinContent(i, binSumPow);
        hInc->SetBinError(i, sqrt(binStatError2));
        if (hIncMad) hIncMad->SetBinError(i, sqrt(binStatMadError2));
        if (hIncShe) hIncShe->SetBinError(i, sqrt(binStatSheError2));
        if (hIncPow) hIncPow->SetBinError(i, sqrt(binStatPowError2));
    }

    // Covariance matrix.
    // We can write:
    //    Y_inc = A * Y_exc, with Y_inc and Y_exc the vector of the respective
    //                       distribution bin contents
    //                       and A_ij = 1 if j >=i, 0 otherwise
    //   => Cov_inc = A * Cov_exc * A^{T}
    for (int m = 0; m < kUncCnt; ++m) {
        if (hCov[m] == 0) continue;
        for (int i = 1; i <= nBins; ++i) {
            for (int j = 1; j <= nBins; ++j) {
                double c = 0;
                for (int k = i; k <= nBins; ++k) {
                    for (int l = j; l <= nBins; ++l) {
                        c += hCov[m]->GetBinContent(k, l);
                    }
                }
                hCovInc[m]->SetBinContent(i, j, c);
            }
        }
    }

    TCanvas *crossSectionPlot = makeCrossSectionPlot(TString(""),
                                                     lumi,
                                                     TString("ZNGoodJets_Zinc"),
                                                     doNormalized,
                                                     hInc,
                                                     hCovInc[kTotSys],
                                                     predictions);
    outputFileName.ReplaceAll("ZNGoodJets_Zexc", "ZNGoodJets_Zinc");
    crossSectionPlot->Draw();
    saveCanvas(crossSectionPlot, outputFileDir, outputFileName);
    createTable(outputFileDir + "/" + outputFileName + "_withLERS",
                "",
                TString("ZNGoodJets_Zinc"),
                doNormalized,
                hInc,
                hCovInc,
                true,
                true);
    createTable(outputFileDir + "/" + outputFileName,
                "",
                TString("ZNGoodJets_Zinc"),
                doNormalized,
                hInc,
                hCovInc,
                false,
                true);
}
