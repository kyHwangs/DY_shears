#include "ConfigVJets.h"
#include "functions.h"
#include "getFilesAndHistogramsZJets.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <iostream>
#include <sstream>
#include <vector>

//--  Setting global variables --------------------------------------------------------------
#include "fileNamesZJets.h"
//-------------------------------------------------------------------------------------------

using namespace std;

namespace /* anonymous */ {
    /**
     * \brief Tunes an axis to be used on a log scale.
     *
     * If the first bin extends to 0, it is modified to start from a higher
     * value. This is needed because otherwise ROOT will make it span half of
     * the plots.
     *
     * \note Histograms with a uniform binning are not supported.
     */
    void prepare_axis_for_log(TAxis &axis)
    {
        if (axis.GetXmin() > 0 || axis.GetNbins() < 2) {
            // Nothing to do
            return;
        }

        // Get the bin edges in a safe container (std::vector over TArray).
        const double *edges_ptr = axis.GetXbins()->GetArray();
        if (edges_ptr == nullptr) {
            // Happens when the histogram has a uniform binning.
            return;
        }
        std::vector<double> edges(edges_ptr, edges_ptr + axis.GetNbins() + 1);

        // First, we get the extent of the second bin. It's proportional to the
        // ratio of the edges.
        double ratio = edges[2] / edges[1];

        // Then, we modify the first bin boundaries so that it gets the same
        // displayed size.
        edges[0] = edges[1] / ratio;

        // We want the lower bound to be a round number, so we round it. It's
        // more complicated that a simple floor() because we want 0.25 to become
        // 0.2 and not 0.0.
        double logfactor = std::pow(10, std::ceil(std::log10(edges[0])) - 1);
        edges[0] = logfactor * std::floor(edges[0] / logfactor);
        edges[0] *= 1.001; // Avoid tick labels.

        // Modify the axis to use the new bin edges.
        axis.Set(edges.size() - 1, edges.data());
    }
} // namespace anonymous

/** Draw data/MC comparison plots from the histograms of the individual contributions.
 * The list of histograms to superimposed is taken from Samples array defined in fileNamesZJets.h.
 * Histogram colours and labels are defined in the same array.
 * The first element is used as signal data and the last one as signal MC.
 * @param lepSel DE, DMu. SE, SMu for Z+jet electron channel, Z+jet muon channel, W+jet...
 * @param histDir location of histogram to use as input.
 * @param recoCompDir directory where the produced plots should be stored
 * @param jetPtMin lower bound for jet pt histograms
 * @param jetEtaMax upper bound in jet |eta| for jet eta histograms
 */
void RecoComparison(
    TString lepSel, TString histoDir, TString recoCompDir, int jetPtMin, int jetEtaMax)
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    //    TString energy = "13TeV";
    ConfigVJets cfg;
    TString energy = TString::Format("%gTeV", cfg.getD("energy"));

    bool isPrel = cfg.getB("preliminaryTag", true);

    int Colors[NFILESDYJETS];
    TString legendNames[NFILESDYJETS];

    //-- get the files, legend names and colors
    //-----------------------------------------------------------
    TFile *fSamples[NFILESDYJETS];
    for (unsigned short i = 0; i < NFILESDYJETS; ++i) {

        int iSample = FilesDYJets[i];

        if (iSample < 0) continue;

        //--- get the file ---
        TString syst = "0";
        if (iSample != 0) syst = "0";
        fSamples[i] =
            getFile(histoDir, lepSel, energy, Samples[iSample].name, jetPtMin, jetEtaMax, "", syst);
        if (!fSamples[i]) return;

        //-- set the legend name for the current file ---
        if (iSample == 0) {
            if (lepSel == "DMu")
                legendNames[i] = " #mu#mu Data";
            else if (lepSel == "DE")
                legendNames[i] = " ee Data";
            else if (lepSel == "SMu")
                legendNames[i] = " #mu Data";
            else if (lepSel == "SE")
                legendNames[i] = " e Data";
            else
                legendNames[i] = " Data";
        }
        //        else if (i == NFILESDYJETS-1)
        //    legendNames[i] = (lepSel == "DMu") ? " Z/#gamma^{*} #rightarrow #mu#mu" :
        //    "Z/#gamma^{*} #rightarrow ee";
        else
            legendNames[i] = Samples[iSample].legendReco;
        //--- set the legend color for the current file ---
        Colors[i] = Samples[iSample].colorReco;
    }
    //-----------------------------------------------------------------------------------------------------

    TString outputFileName = recoCompDir;
    system("mkdir -p " + recoCompDir);
    outputFileName += "Comparison_" + lepSel + "_" + energy + "_Data_All_MC";
    // outputFileName += "_JetPtMin_";
    // outputFileName += jetPtMin;
    // outputFileName += "_JetEtaMax_";
    // outputFileName += jetEtaMax;
    //--- create the directory if it doesn't exist ---
    system("mkdir -p " + outputFileName);
    TString outputFileRoot = outputFileName + ".root";
    cout << "Output directory is: " << outputFileName << endl;
    cout << "Output root file is: " << outputFileRoot << endl;
    TFile *outputFile = new TFile(outputFileRoot, "RECREATE");
    outputFile->cd();

    //--- vector to contain the names and titles of the reco histograms ---
    vector<TString> vhNames;
    vector<TString> vhTitles;
    //---------------------------------------------------------------------

    //--- total number of histograms inside the data file
    //    (included gen and TH2 that we dont want for comparison) ---
    unsigned short nTotHist = fSamples[0]->GetListOfKeys()->GetEntries();
    for (unsigned short i = 0; i < nTotHist; ++i) {
        TString hName = fSamples[0]->GetListOfKeys()->At(i)->GetName();
        TH1D *hTemp = (TH1D *)fSamples[0]->Get(hName);
        TString hTitle = hTemp->GetName();
        //--- skip histogram if it is gen or has no entry or is not a TH1 ---
        if (hName.Index("gen") >= 0 || hTemp->GetEntries() < 1 ||
            !hTemp->InheritsFrom(TH1D::Class()))
            continue;

        //--- store the histograme name  and title ---
        vhNames.push_back(hName);
        vhTitles.push_back(hTitle);
    }

    //--- vhNames size gives us the number of reco histograms
    //    interesting for comparison at reco level
    int nHist = vhNames.size();

    TH1D *hist[NFILESDYJETS][nHist];
    THStack *hSumMC[nHist];
    TLegend *legend[nHist];

    TLatex *cmsColl = new TLatex();
    cmsColl->SetTextSize(0.04);
    cmsColl->SetTextFont(61);
    cmsColl->SetLineWidth(2);
    cmsColl->SetTextColor(kBlack);
    cmsColl->SetNDC();
    cmsColl->SetTextAlign(11);

    TLatex *cmsPrel = new TLatex();
    cmsPrel->SetTextSize(0.04);
    cmsPrel->SetTextFont(52);
    cmsPrel->SetLineWidth(2);
    cmsPrel->SetTextColor(kBlack);
    cmsPrel->SetNDC();
    cmsPrel->SetTextAlign(11);

    TLatex *jetAlgo = new TLatex();
    jetAlgo->SetTextSize(0.035);
    jetAlgo->SetTextFont(42);
    jetAlgo->SetLineWidth(2);
    jetAlgo->SetTextColor(kBlack);
    jetAlgo->SetNDC();
    jetAlgo->SetTextAlign(11);

    TLatex *jetCuts = new TLatex();
    jetCuts->SetTextSize(0.035);
    jetCuts->SetTextFont(42);
    jetCuts->SetLineWidth(2);
    jetCuts->SetTextColor(kBlack);
    jetCuts->SetNDC();
    jetCuts->SetTextAlign(11);

    TLatex *intLumi = new TLatex();
    intLumi->SetTextSize(0.03);
    intLumi->SetTextFont(42);
    intLumi->SetLineWidth(2);
    intLumi->SetTextColor(kBlack);
    intLumi->SetNDC();
    intLumi->SetTextAlign(31);

    for (unsigned int i = 0; i < NFILESDYJETS; ++i) {
        double scale = cfg.getD(TString("scale_") + Samples[FilesDYJets[i]].name, 1.);
        if (scale != 1.) {
            std::cout << "Info. Scale factor " << scale
                      << " will be applied on the event yield of sample "
                      << Samples[FilesDYJets[i]].name << "\n";
        }

        for (int j = 0; j < nHist; ++j) {
            hist[i][j] = getHisto(fSamples[i], vhNames[j]);
            if (!hist[i][j]) {
                std::cerr << "Histogram " << vhNames[j] << " was not found for sample "
                          << Samples[FilesDYJets[i]].name << "\n";
                continue;
            }
            hist[i][j]->Scale(scale);
            hist[i][j]->SetTitle(vhTitles[j]);
            if (i == 0) {
                hist[0][j]->SetMarkerStyle(20);
                hist[0][j]->SetMarkerColor(Colors[0]);
                hist[0][j]->SetLineColor(Colors[0]);
                hSumMC[j] = new THStack(vhNames[j], vhTitles[j]);

                // if (!doPASPlots) {
                // legend[j] = new TLegend(0.72, 0.5, 0.76, 0.86);
                // legend[j]->SetTextSize(0.032);
                //}
                // else {
                legend[j] = new TLegend(0.63, 0.60, 0.81, 0.87);
                legend[j]->SetTextSize(0.042);
                //}
                legend[j]->SetFillStyle(0);
                legend[j]->SetBorderSize(0);
                legend[j]->SetTextFont(42);
            } else {
                hist[i][j]->SetFillStyle(1001);
                hist[i][j]->SetFillColor(Colors[i]);
                hist[i][j]->SetLineColor(Colors[i]);
                hSumMC[j]->Add(hist[i][j]);
                // if (!doPASPlots || i == 1 || i == 3 || i == 5 || i == 11)
                // legend[j]->AddEntry(hist[i][j], legendNames[i], "f");
            }
        } // next histo j
    }     // next file i

    // Fill the legend in reverse order of drawing in order
    // that the legend lines order matches with stacked histogram one.
    for (int j = 0; j < nHist; ++j) {
        if (NFILESDYJETS > 0) legend[j]->AddEntry(hist[0][j], legendNames[0], "ep");
        for (int i = NFILESDYJETS - 1; i > 0; --i) {
            legend[j]->AddEntry(hist[i][j], legendNames[i], "f");
        }
    }
    // reads integrated luminosity
    double lumi = -1;
    TH1 *Lumi;
    fSamples[0]->GetObject("Lumi", Lumi);
    if (Lumi)
        lumi = Lumi->GetBinContent(1);
    else
        cerr << "Warning: Lumi histogram was not found. The integrated luminosity indicaion will "
                "be missing from the plots.\n";

    cout << "Now creating the pdf files ..." << endl;
    double minRatioY = cfg.getD("minRatioYReco", 0.51);
    ;
    double maxRatioY = cfg.getD("maxRatioYReco", 1.49);

    for (unsigned short i = 0; i < nHist; ++i) {

        TCanvas *canvas = new TCanvas(vhNames[i], vhNames[i], 700, 900);
        canvas->cd();

        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
        pad1->SetTopMargin(0.11);
        pad1->SetBottomMargin(0.);
        pad1->SetRightMargin(0.03);
        pad1->SetTicks();
        pad1->SetLogy();
        // pad1->SetLogx(); //DJALOG
        pad1->Draw();
        pad1->cd();

        if (vhNames[i].Index("Phistar") >= 0 || vhNames[i].Index("Mass_Zinc0jet") >= 0 ||
            vhNames[i].Index("ZPt_Zinc") >= 0)
            pad1->SetLogx();
        if (vhNames[i].Index("ZPt_Zinc0jetM115_135") >= 0) pad1->SetLogx(0);

        TH1D *hRatio = (TH1D *)hSumMC[i]->GetStack()->Last()->Clone();
        // Need to draw MC Stack first other wise
        // cannot access Xaxis !!!

        if (vhNames[i].Index("ZNGoodJets_Zexc") >= 0) {
            std::cout << __FILE__ << ":" << __LINE__
                      << ". Range of ZNGoodJets_Zexc x-axis is being modified.!\n";
            // hSumMC[i]->GetXaxis()->Set(maxX-minX, minX, maxX);
            //	hSumMC[i]->GetXaxis()->SetRangeUser(-0.5, 4.5);
            hRatio->GetXaxis()->SetBinLabel(1, "= 0");
            hRatio->GetXaxis()->SetBinLabel(2, "= 1");
            hRatio->GetXaxis()->SetBinLabel(3, "= 2");
            hRatio->GetXaxis()->SetBinLabel(4, "= 3");
            hRatio->GetXaxis()->SetBinLabel(5, "= 4");
            hRatio->GetXaxis()->SetBinLabel(6, "= 5");
            hRatio->GetXaxis()->SetBinLabel(7, "= 6");
            // hSumMC[i]->GetXaxis()->SetBinLabel(8, "= 7");
            //     hSumMC[i]->GetXaxis()->SetBinLabel(9, "= 8");
            //  hSumMC[i]->GetXaxis()->SetLabelSize(0.18);
            //  hSumMC[i]->GetXaxis()->SetLabelOffset(0.01);
        }

        if (vhNames[i].Index("ZNGoodJets_Zinc") >= 0) {
            std::cout << __FILE__ << ":" << __LINE__
                      << ". Range of ZNGoodJets_Zexc x-axis is being modified.!\n";
            // hSumMC[i]->GetXaxis()->Set(maxX-minX, minX, maxX);
            //	  if (vhNames[i].Index("ZNGoodJets_Zinc") >= 0) {
            // hSumMC[i]->GetXaxis()->SetRangeUser(-0.5, 4.5);
            hRatio->GetXaxis()->SetBinLabel(1, "#geq 0");
            hRatio->GetXaxis()->SetBinLabel(2, "#geq 1");
            hRatio->GetXaxis()->SetBinLabel(3, "#geq 2");
            hRatio->GetXaxis()->SetBinLabel(4, "#geq 3");
            hRatio->GetXaxis()->SetBinLabel(5, "#geq 4");
            hRatio->GetXaxis()->SetBinLabel(6, "#geq 5");
            hRatio->GetXaxis()->SetBinLabel(7, "#geq 6");
            // hSumMC[i]->GetXaxis()->SetBinLabel(8, "= 7");
            //     hSumMC[i]->GetXaxis()->SetBinLabel(9, "= 8");
            //  hSumMC[i]->GetXaxis()->SetLabelSize(0.18);
            //  hSumMC[i]->GetXaxis()->SetLabelOffset(0.01);

            // DJALOG
            printf("Printing ZNGoodJets_Zinc bin values and signal to background ratios...\n");
            printf("Entries in stack=%d\n", hSumMC[i]->GetStack()->GetEntries());

            printf("         Bin|     Signal| Background|      Ratio\n");
            for (int iBin = 1; iBin < 8; iBin++) {
                double signal = ((TH1D *)hSumMC[i]->GetStack()->At(4))->GetBinContent(iBin);
                double background = ((TH1D *)hSumMC[i]->GetStack()->At(3))->GetBinContent(iBin);
                signal = signal - background;
                printf(
                    "%12d|%11.2F|%11.2F|%11.2F\n", iBin, signal, background, signal / background);
            }
            printf("\n");

            for (int iBin = 0; iBin < 8; iBin++) {
                for (int iHisto = 1; iHisto < 6; iHisto++) {
                    if (iBin == 0) {
                        printf("%30s|", legendNames[iHisto].Data());
                    } else {
                        printf("%30F|", hist[iHisto][i]->GetBinContent(iBin));
                    }
                }
                printf("\n");
            }

            printf("\n");
            printf("WJets Integral = %F", hist[3][i]->Integral());

            printf("\n");
        }
        if (vhNames[i].Index("ZNGoodJets_SameChargePair_Zinc") >= 0) {
            // DJALOG
            printf("Printing ZNGoodJets_SameSignCharge_Zinc bin values and signal to background "
                   "ratios...\n");
            printf("Entries in stack=%d\n", hSumMC[i]->GetStack()->GetEntries());

            printf("         Bin|       Data|     Signal| Background|      Ratio|        QCD\n");
            for (int iBin = 1; iBin < 8; iBin++) {
                double signal = ((TH1D *)hSumMC[i]->GetStack()->At(4))->GetBinContent(iBin);
                double background = ((TH1D *)hSumMC[i]->GetStack()->At(3))->GetBinContent(iBin);
                signal = signal - background;
                printf("%12d|%11.2F|%11.2F|%11.2F|%11.2F|%11.2F\n",
                       iBin,
                       hist[0][i]->GetBinContent(iBin),
                       signal,
                       background,
                       signal / background,
                       hist[0][i]->GetBinContent(iBin) - signal - background);
            }
        }

        hSumMC[i]->Draw("HIST");
        if (vhNames[i].Index("ZMass_Z") >= 0) {
            hist[0][i]->GetXaxis()->SetRangeUser(71, 110.9);
            hSumMC[i]->GetXaxis()->SetRangeUser(71, 110.9);
            hRatio->GetXaxis()->SetRangeUser(71, 110.9);
        }
        // if (vhNames[i].Index("JetEta") >= 0){
        //    hist[0][i]->GetXaxis()->SetRangeUser(-2.4,2.4);
        //    hSumMC[i]->GetXaxis()->SetRangeUser(-2.4,2.4);
        //    hRatio->GetXaxis()->SetRangeUser(-2.4,2.4);

        //}

        hSumMC[i]->SetTitle("");
        hSumMC[i]->GetYaxis()->SetLabelSize(0.04);
        hSumMC[i]->GetYaxis()->SetLabelOffset(0.002);
        hSumMC[i]->GetYaxis()->SetTitle("# Events");
        hSumMC[i]->GetYaxis()->SetTitleSize(0.04);
        hSumMC[i]->GetYaxis()->SetTitleOffset(1.32);
        hSumMC[i]->SetMinimum(8);
        hSumMC[i]->SetMaximum(100 * hSumMC[i]->GetMaximum());

        // DJALOG
        // Change the axes of the jet pt to be the cut value.
        if (vhNames[i].Index("JetPt") >= 0) {
            printf("Pt Histo\n");
            hSumMC[i]->GetXaxis()->SetRangeUser(
                jetPtMin, hSumMC[i]->GetXaxis()->GetBinUpEdge(hSumMC[i]->GetXaxis()->GetLast()));
        }

        // first pad plots
        hist[0][i]->DrawCopy("e same");
        legend[i]->Draw();

        cmsColl->DrawLatex(0.17, 0.83, isPrel ? "CMS Preliminary" : "CMS");
        // if (energy == "13TeV") intLumi->DrawLatex(0.5,0.77, "#int L dt = 2.25 fb^{-1},  #sqrt{s}
        // = 13 TeV");
        if (energy == "13TeV") {
            if (lumi > 0)
                intLumi->DrawLatex(
                    0.5,
                    0.77,
                    TString::Format("#int L dt = %.3g fb^{-1},  #sqrt{s} = 13 TeV", lumi / 1000.));
                    //  TString::Format("#int L dt = %.3g fb^{-1},  #sqrt{s} = 13 TeV", 35.9));
            else
                intLumi->DrawLatex(0.5, 0.77, "#sqrt{s} = 13 TeV");
        }

        if (vhNames[i].Index("inc0") < 0) {
            ostringstream ptLegend;

            if (vhNames[i].Index("ZB") > 0)
                ptLegend << "p_{T}^{jet} > " << jetPtMin
                         << " GeV,  |y^{jet}| < 2.4, N_{jets} #geq 1";
            else if (vhNames[i].Index("inc1") > 0)
                ptLegend << "p_{T}^{jet} > " << jetPtMin
                         << " GeV,  |y^{jet}| < 2.4, N_{jets} #geq 1";
            else if (vhNames[i].Index("inc2") > 0)
                ptLegend << "p_{T}^{jet} > " << jetPtMin
                         << " GeV,  |y^{jet}| < 2.4, N_{jets} #geq 2";
            else if (vhNames[i].Index("inc3") > 0)
                ptLegend << "p_{T}^{jet} > " << jetPtMin
                         << " GeV,  |y^{jet}| < 2.4, N_{jets} #geq 3";
            else
                ptLegend << "p_{T}^{jet} > " << jetPtMin << " GeV,  |y^{jet}| < 2.4";

            jetCuts->DrawLatex(0.17, 0.66, ptLegend.str().c_str());
            jetAlgo->DrawLatex(0.17, 0.715, "anti-k_{t} jets,  R = 0.4");
            if (vhNames[i].Index("ptLow") > 0)
                jetCuts->DrawLatex(0.17, 0.61, " p_{T}(Z) #leq 50 GeV");
            if (vhNames[i].Index("ptHigh") > 0)
                jetCuts->DrawLatex(0.17, 0.61, " p_{T}(Z) > 50 GeV");

            pad1->Draw();
        }
        //-------------------------
        // cmsColl->DrawLatex(0.13,0.82, "CMS Preliminary");
        //        cmsPrel->DrawLatex(0.13,0.78, "Preliminary");
        //        if (energy == "7TeV")      intLumi->DrawLatex(0.97,0.9, "5.05 fb^{-1} (7 TeV)");
        //        else if (energy == "8TeV") intLumi->DrawLatex(0.97,0.9, "19.6 fb^{-1} (8 TeV)");
        // if(lumi >= 0) intLumi->DrawLatex(0.97, 0.9, TString::Format("%.3g fb^{-1} (%s)", lumi,
        // energy.Data()));
        if (vhNames[i].Index("inc0") < 0) {
            //            if (!doPASPlots) {
            //                ostringstream ptLegend;
            //                if (vhNames[i].Index("JetPt_Zinc") > 0) {
            //                    ptLegend << "p_{T}^{jet} > 20 GeV,  |y^{jet}| < " <<
            //                    (0.1*jetEtaMax);
            //                }
            //                else {
            //                    ptLegend << "p_{T}^{jet} > " << jetPtMin << "GeV,  |y^{jet}| < "
            //                    << (0.1*jetEtaMax);
            //                }
            //                jetAlgo->DrawLatex(0.13,0.68, "anti-k_{t} jets,  R = 0.4");
            //                jetCuts->DrawLatex(0.13,0.63, ptLegend.str().c_str());
            //            }
            pad1->Draw();
        }
        canvas->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
        pad2->SetTopMargin(0.);
        pad2->SetBottomMargin(0.3);
        pad2->SetRightMargin(0.03);
        pad2->SetGridy();
        // pad2->SetLogx(); //DJALOG
        pad2->SetTicks();
        pad2->Draw();
        pad2->cd();

        bool logx = false;

        if (vhNames[i].Index("Phistar") >= 0 || vhNames[i].Index("Mass_Zinc0jet") >= 0 ||
            vhNames[i].Index("ZPt_Zinc") >= 0) {
            logx = true;
            pad2->SetLogx();
        }
        if (vhNames[i].Index("ZPt_Zinc0jetM115_135") >= 0) {
            logx = false;
            pad2->SetLogx(0);
        }



        hRatio->SetStats(0);
        hRatio->SetTitle("");

        hRatio->SetMarkerStyle(20);
        hRatio->SetMarkerColor(Colors[0]);
        hRatio->SetLineColor(Colors[0]);

        hRatio->GetXaxis()->SetTickLength(0.03);
        hRatio->GetXaxis()->SetTitleSize(0.1);
        hRatio->GetXaxis()->SetTitleOffset(1.2);
        hRatio->GetXaxis()->SetLabelSize(0.10);
        hRatio->GetXaxis()->SetLabelOffset(0.017);

        if (logx) {
            prepare_axis_for_log(*hRatio->GetXaxis());
        }
        // DJALOG
        // Change the axes of the jet pt to be the cut value.
        if (vhNames[i].Index("JetPt") >= 0) {
            printf("Pt Histo\n");
            hRatio->GetXaxis()->SetRangeUser(
                jetPtMin, hRatio->GetXaxis()->GetBinUpEdge(hRatio->GetXaxis()->GetLast()));
        }

        //        hRatio->GetYaxis()->SetRangeUser(0.51,1.49);
        hRatio->GetYaxis()->SetRangeUser(minRatioY, maxRatioY);

        hRatio->GetYaxis()->SetNdivisions(5, 5, 0);
        hRatio->GetYaxis()->SetTitle("Simulation/Data");
        hRatio->GetYaxis()->SetTitleSize(0.1);
        hRatio->GetYaxis()->SetTitleOffset(0.5);
        hRatio->GetYaxis()->CenterTitle();
        hRatio->GetYaxis()->SetLabelSize(0.08);

/*
      if(vhNames[i].Index("ZMass_Zinc0jet") >= 0){
          for(int j =1 ; j<=hRatio->GetNbinsX();j++)    
                  std::cout<<hRatio->GetBinContent(j) << ","; //<< hist[0][i] ->GetBinContent(j) <<"\n ";
                  std::cout<<std::endl;
         }
*/

        hRatio->Divide(hist[0][i]);


     //  if(vhNames[i].Index("ZMass_Zinc0jet") >= 0){
     //       for(int i =1 ; i<=hRatio->GetNbinsX();i++)    
      //            std::cout<<hRatio->GetBinContent(i)<<" , ";
     //             std::cout<<std::endl;
     //    }
        hRatio->DrawCopy("EP");

        canvas->cd();
        canvas->Update();

        // TString outputFilePDF = outputFileName + "/" + vhNames[i] + ".pdf";
        // canvas->Print(outputFilePDF);
        // outputFile->cd();
        // canvas->Write();

        // TString outputFileBase = outputFileName + "/" + vhNames[i];
        // canvas->SaveAs(outputFileBase + ".root");
        // canvas->SaveAs(outputFileBase + ".C");
        // canvas->SaveAs(outputFileBase + ".png");
        saveCanvas(canvas, outputFileName, vhNames[i]);

        hSumMC[i]->SetMaximum(1.5 * hSumMC[i]->GetMaximum());
        TCanvas *tmpCanvas = (TCanvas *)canvas->Clone();
        tmpCanvas->cd();
        tmpCanvas->Draw();
        TPad *tmpPad = (TPad *)tmpCanvas->GetPrimitive("pad1");
        tmpPad->SetLogy(0);
        vhNames[i] += "_Lin";
        tmpCanvas->SetTitle(vhNames[i]);
        tmpCanvas->SetName(vhNames[i]);
        tmpCanvas->Update();
        // TString outputFileLinPDF = outputFileName + "/" + vhNames[i] + ".pdf";
        // tmpCanvas->Print(outputFileLinPDF);
        // outputFile->cd();
        // tmpCanvas->Write();

        // TString outputFileLinBase = outputFileName + "/" + vhNames[i];
        // tmpCanvas->SaveAs(outputFileLinBase + ".root");
        // tmpCanvas->SaveAs(outputFileLinBase + ".C");
        saveCanvas(tmpCanvas, outputFileName, vhNames[i]);
    }

    outputFile->cd();
    outputFile->Close();

    //-- Close all the files ------------------------------
    for (unsigned short i(0); i < NFILESDYJETS; ++i) closeFile(fSamples[i]);
    //-----------------------------------------------------
}
