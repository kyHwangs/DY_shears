#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    gROOT->SetBatch();

    if (argc < 3) {
        cout << "You need to provide lepSel and variable." << endl;
        cout << "ex:  ./plotRespMatrix DMu ZNGoodJets_Zexc" << endl;
        return 1;
    }

    TApplication *myApp = new TApplication("myApp", &argc, argv);

    TString lepSel(argv[1]);
    TString variable(argv[2]);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.0f");

    TString outputFileName = "ResponseMatrix/" + lepSel + "_" + variable + "_ResponseMatrix";
    TString respName = "hresponse" + variable;

    TFile *fMad = new TFile("HistoFilesBarrel/" + lepSel + "_13TeV_DYJets_UNFOLDING_Syst_0.root");
    TH1D *hMad = (TH1D *)fMad->Get(variable);
    TH2D *hrespMad = (TH2D *)fMad->Get(respName);
    TH2D *hrespNormMad = (TH2D *)hrespMad->Clone();

    /*    TFile *fShe = new TFile("../HistoFiles/" + lepSel +
       "_13TeV_DYJets_Sherpa_Bugra_1_13_UNFOLDING_dR_TrigCorr_0_Syst_0_JetPtMin_30_JetEtaMax_24.root");
        TH1D *hShe = (TH1D*) fShe->Get(variable);
        TH2D *hrespShe = (TH2D*) fShe->Get(respName);
        TH2D *hrespNormShe = (TH2D*) hrespShe->Clone();
    */
    TString title(hMad->GetTitle());
    TString xTitle(hMad->GetXaxis()->GetTitle());
    TString yTitle("gen " + xTitle);

    hrespNormMad->SetTitle("aMC@NLO+Pythia8 Resp. Matrix for " + title);
    hrespNormMad->GetXaxis()->SetTitle(xTitle);
    hrespNormMad->GetXaxis()->SetTitleOffset(1.4);
    hrespNormMad->GetYaxis()->SetTitle(yTitle);
    hrespNormMad->GetYaxis()->SetTitleOffset(1.6);
    hrespNormMad->GetZaxis()->SetTitle("");
    hrespNormMad->GetZaxis()->SetRangeUser(0, 100);
    hrespNormMad->SetMarkerSize(1.7);
    /*
        hrespNormShe->SetTitle("Sherpa Resp. Matrix for " + title);
        hrespNormShe->GetXaxis()->SetTitle(xTitle);
        hrespNormShe->GetXaxis()->SetTitleOffset(1.4);
        hrespNormShe->GetYaxis()->SetTitle(yTitle);
        hrespNormShe->GetYaxis()->SetTitleOffset(1.6);
        hrespNormShe->GetZaxis()->SetTitle("");
        hrespNormShe->GetZaxis()->SetRangeUser(0,100);
    */

    int nBinsX = hrespMad->GetNbinsX();
    int nBinsY = hrespMad->GetNbinsY();

    for (int i(0); i <= nBinsY + 1; i++) {
        double totRowMad(0);
        //        double totRowShe(0);
        for (int j(0); j <= nBinsX + 1; j++) {
            totRowMad += hrespMad->GetBinContent(j, i);
            //            totRowShe += hrespShe->GetBinContent(j, i);
        }
        for (int j(0); j <= nBinsX + 1; j++) {
            double binContentMad = hrespMad->GetBinContent(j, i);
            //            double binContentShe = hrespShe->GetBinContent(j, i);
            hrespNormMad->SetBinContent(j, i, 100 * binContentMad / totRowMad);
            //            hrespNormShe->SetBinContent(j, i, 100*binContentShe/totRowShe);
        }
    }

    TCanvas *cMad = new TCanvas("cMad", "MadGraph", 800, 800);
    cMad->Connect("Closed()", "TApplication", myApp, "Terminate()");
    cMad->cd();

    TPad *padMad = new TPad("padMad", "padMad", 0, 0, 1, 1);
    padMad->SetRightMargin(0.12);
    padMad->SetLeftMargin(0.12);
    padMad->SetBottomMargin(0.12);
    if (variable.Index("Pt") >= 0 || variable.Index("HT") >= 0) {
        padMad->SetLogy();
        padMad->SetLogx();
    }
    padMad->Draw();
    padMad->cd();
    hrespNormMad->DrawCopy("colztext");
    padMad->Draw();
    cMad->Update();
    cMad->cd();
    cMad->SaveAs(outputFileName + "_MadGraph.png");
    cMad->SaveAs(outputFileName + "_MadGraph.pdf");
    cMad->SaveAs(outputFileName + "_MadGraph.ps");
    cMad->SaveAs(outputFileName + "_MadGraph.eps");
    cMad->SaveAs(outputFileName + "_MadGraph.root");
    /*
        TCanvas *cShe = new TCanvas("cShe", "Sherpa", 800, 800);
        cShe->Connect("Closed()", "TApplication", myApp,  "Terminate()");
        cShe->cd();

        TPad *padShe = new TPad("padShe", "padShe", 0, 0, 1, 1);
        padShe->SetRightMargin(0.12);
        padShe->SetLeftMargin(0.12);
        padShe->SetBottomMargin(0.12);
        padShe->Draw();
        padShe->cd();
        hrespNormShe->DrawCopy("colztext");
        padShe->Draw();
        cShe->Update();
        cShe->cd();
        cShe->SaveAs(outputFileName + "_Sherpa.png");
        cShe->SaveAs(outputFileName + "_Sherpa.pdf");
        cShe->SaveAs(outputFileName + "_Sherpa.root");
    */

    // myApp->Run();
    fMad->Close();
    //    fShe->Close();
    return 0;
}
