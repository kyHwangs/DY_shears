//-*- mode: c++; c-basic-offset: 4 -*-
#include "PlotSettings.h"
#include "TStyle.h"
#include "TFile.h"
#include "TMath.h"
#include "TMathText.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TParameter.h"
#include <sstream>
#include <string>
#include <set>
#include <limits>
#include "rootFunctions.h"
#include "functions.h"
#include "variablesOfInterestZJets.h"
#include "getFilesAndHistogramsZJets.h"
#include "fixYscale.h"
#include "SectionedConfig.h"
#include "Uncertainties.h"

using namespace std;
#include "ConfigVJets.h"

extern ConfigVJets cfg;

#define NEW_PLOTS

const static double padHeightRatio = 0.3; //height ration of a ratio frame to the top frame.

static double ratioYaxisExtend = 1.;

static void extendAxis(double& minY, double& maxY){
    double mid = 0.5 * (minY + maxY);
    minY = mid - ratioYaxisExtend * (mid - minY);
    maxY = mid + ratioYaxisExtend * (maxY - mid);
}

static struct shearsTextStyles {
    const int defaultFont = 43;    
    const int xLabelSize  = 20;
    const int xTitleSize = 20;
    const float xTitleOffset = 1.2;
    //    const float xTitleOffset = 2.; 
    const int yLabelSize  = 20;
    const float mainYTitleOffset = 1.1;
    //    const int mainYTitleSize = 42;
    const int mainYTitleSize = 24;
    //--//    const int ratioYTitleSize = 15;
    //--//    const float ratioYTitleOffset = 0.66; //this number will be devided by hbottomratio
    const int ratioYTitleSize = 18;
    const float ratioYTitleOffset = 1.3; //this number will be devided by the number of ratio plots
    const int cmsLabelFont = 63;
    const int cmsLabelSize = 30;
    const int prelimLabelFont = 53;
    const int prelimLabelSize = 20;
    const int lumiLabelFont = 43;
    const int lumiLabelSize = 20;
    const int mainLegendTextSize = 15;//
    const int ratioLegendTextSize = 20;
    const int descTextSize = 16;
} ts;

#ifdef NEW_PLOTS
static double hratio = 1.;
static double hbottomratio = 1.;
static double ratioYTitleOffsetCorr = 1.;

TCanvas* setAndDrawTPad(const TString& canvasName, const TString& canvasTitle, TPad*& plot, int plotNumber, int numbOfGenerator){

    if(numbOfGenerator<=1) numbOfGenerator = 1;

    double absHeightTopFrameTopMargin = 44;
    double absHeightTopFrameHisto = 356;
    double absHeightTopFrameTot = absHeightTopFrameTopMargin + absHeightTopFrameHisto;
    double absHeightGap = 4;
    double absHeightMidFrame = 120;
    double absHeightBottomMargin = 65;
    double absHeightBottomFrame = absHeightMidFrame * ratioYaxisExtend;
    double absHeightBottomFrameTot = absHeightBottomFrame + absHeightBottomMargin;
    double absHeightCanvas = absHeightTopFrameTot + absHeightGap + absHeightMidFrame * (numbOfGenerator - 1) + absHeightBottomFrameTot;
    
    hratio       = absHeightMidFrame / absHeightCanvas;;
    hbottomratio = absHeightBottomFrameTot / absHeightCanvas;;
    ratioYTitleOffsetCorr = absHeightCanvas / absHeightTopFrameTot;

    //ad-hoc correction
    //    absHeightCanvas -= 30 * (numbOfGenerator - 4);

//--//    double margin0 = absHeightTopFrameTopMargin / absHeightTopFrameHeight;
//--//    double margin1 = absHeightGap / absHeightCanvas; 
//--//    double margin2 = 

//--//    double actualPadHeightRatio = padHeightRatio * ratioYaxisExtend;
//--//    double margin0 = 0.11;   //margin on top of main subframe
//--//    double margin1 = 0.005;  //margin between main subframe and ratio plot frames
//--//    double margin2 = 0.35;    //margin at bottom of all frames (bottom margin of bottom subframe)    
//--//    double htop = 1. * (1. + margin0 + margin1);
//--//    hratio = htop * actualPadHeightRatio;
//--//    hbottomratio = hratio / (1. - margin2);
//--//    double htot = htop + hratio * (numbOfGenerator -1) + hbottomratio;

//--//    //rescale to fit within htot = 0.98:
//--//    double fac = 0.98 / htot;
//--//    htop *=  fac;
//--//    hratio *= fac;
//--//    hbottomratio *= fac;
    //---

    //    std::cout << "fac, htop, hratio, hborromratio, sum: "
    //	      << fac << ", " << htop << ", " << hratio << ", " 
    //	      << hbottomratio << ", "
    //	      << (htop + (numbOfGenerator-1) * hratio + hbottomratio) << "\n";

    TCanvas* c;
    if(plotNumber==1){
	//	c = new TCanvas(canvasName, canvasTitle, 600, 400*(1 + padHeightRatio * ratioYaxisExtend * ratioFrames.size()));
	//hratio =        fac * htop  * actualPadHeightRatio
	//hbottomratio = (fac * htop) / (1.-magin2) * actualPadHeightRatio
	//c = new TCanvas(canvasName, canvasTitle, 600, 400*(1 + actualPadHeightRatio * ( 1. / (1. - margin2) + 1. * (numbOfGenerator-1))));
	//--//	c = new TCanvas(canvasName, canvasTitle, 600, 400*(1 + (hbottomratio / htop) +  (hratio / htop)  * (numbOfGenerator-1)));
	c = newTCanvas(canvasName, canvasTitle, 600, absHeightCanvas);
	c->SetTopMargin(0);
	c->SetBottomMargin(0);
    } else if (plot != 0){
	c = plot->GetCanvas();
    }
    TString pname = TString::Format("plot%d", plotNumber);
    if(plot==0) plot = new TPad(pname, pname, 0, 0, 0, 0); 

    plot->SetCanvas(c);
    plot->SetNumber(plotNumber);
    double y0 = 0;
    double y1 = 0;
    if(plotNumber == 1){
	//--//	 plot->SetTopMargin(margin0);
	//--//	 plot->SetBottomMargin(margin1);
	plot->SetTopMargin(absHeightTopFrameTopMargin / absHeightTopFrameTot);
	plot->SetBottomMargin(0);
	//--//y0 = 0.99;
	//--// y1 = 0.99 - htop;
	y0 = absHeightCanvas;
	y1 = y0 - absHeightTopFrameTot;
    } else if(plotNumber < numbOfGenerator + 1){
	 plot->SetTopMargin(0.0);
	 plot->SetBottomMargin(0.0);
	 //	 y0 = 0.99 - htop * (1 + margin0 + margin1) - hratio * (plotNumber - 2);
	 //--// y0 = 0.99 - htop - hratio * (plotNumber - 2);
	 //--//y1 = y0 - hratio;
	 y0 = absHeightCanvas - absHeightGap - absHeightTopFrameTot - absHeightMidFrame * (plotNumber - 2);
	 y1 = y0 - absHeightMidFrame;
    } else{
	 plot->SetTopMargin(0.0);
	 //--//plot->SetBottomMargin(margin2);
	 plot->SetBottomMargin(absHeightBottomMargin / absHeightBottomFrameTot); 
	 //--// y0 = 0.99 - htop - hratio * (plotNumber - 2);
	 //--// y1 = y0 - hbottomratio;
	 y0 = absHeightCanvas - absHeightGap - absHeightTopFrameTot - absHeightMidFrame * (plotNumber - 2);
	 y1 = y0 - absHeightBottomFrameTot;
    }
    //--//plot->SetPad(0.01, y1, 0.99, y0);
    plot->SetPad(0, y1 / absHeightCanvas, 1, y0 / absHeightCanvas);
    if (plotNumber == 1 && (canvasName.Index("Eta") < 0 && canvasName.Index("AbsRapidity") < 0 && canvasName.Index("DPhi") < 0)) plot->SetLogy();
    if (plotNumber == 1 && canvasName.Index("DPhiZFirstJet") > 0) plot->SetLogy();
    if (canvasName.Index("ZPt_") > 0){
	plot->SetLogx();
	plot->SetLogy(0);
    }
    plot->SetLeftMargin(0.13);
    plot->SetRightMargin(0.07);
    plot->SetFillStyle(0);
    plot->Draw();
    plot->cd();

//    std::cout << ">> " << canvasName << "\t" << plotNumber << ":\t" << plot->VtoAbsPixel(0) - plot->VtoAbsPixel(1) 
//	      << "\t" << ( plotNumber == 1 ? absHeightTopFrameTot : 
//			   ((plotNumber-1) == numbOfGenerator ? absHeightBottomFrameTot : absHeightMidFrame))
//	      << "\t" << c->VtoAbsPixel(0) - c->VtoAbsPixel(1) << "\t" << absHeightCanvas
//	      << "\n";

    return c;
}
#else
TCanvas* setAndDrawTPad(TString canvasName, TPad *plot, int plotNumber, int numbOfGenerator)
{
    plot->SetNumber(plotNumber);
    if (numbOfGenerator == 1) {
        if (plotNumber == 1) {
            plot->SetPad(0.01, 0.35, 0.99, 0.99);
            plot->SetTopMargin(0.11);
            plot->SetBottomMargin(0.005);
        }
        else if (plotNumber == 2) {
            plot->SetPad(0.01, 0.01, 0.99, 0.35);
            plot->SetTopMargin(0.0);
            plot->SetBottomMargin(0.3);
        }
    }
    else if (numbOfGenerator == 2) {
        if (plotNumber == 1) {
            plot->SetPad(0.01, 0.45, 0.99, 0.99);
            plot->SetTopMargin(0.11);
            plot->SetBottomMargin(0.005);
        }
        else if (plotNumber == 2) {
            plot->SetPad(0.01, 0.27, 0.99, 0.45);
            plot->SetTopMargin(0.0);
            plot->SetBottomMargin(0.0);
        }
        else if (plotNumber == 3) {
            plot->SetPad(0.01, 0.01, 0.99, 0.27);
            plot->SetTopMargin(0.0);
            plot->SetBottomMargin(0.3);
        }

    }
    else if (numbOfGenerator == 3) {
        if (plotNumber == 1) {
            plot->SetPad(0.01, 0.55, 0.99, 0.99);
            plot->SetTopMargin(0.11);
            plot->SetBottomMargin(0.005);
        }
        else if (plotNumber == 2) {
            plot->SetPad(0.01, 0.39, 0.99, 0.55);
            plot->SetTopMargin(0.0);
            plot->SetBottomMargin(0.0);
        }
        else if (plotNumber == 3) {
            plot->SetPad(0.01, 0.23, 0.99, 0.39);
            plot->SetTopMargin(0.0);
            plot->SetBottomMargin(0.0);
        }
        else if (plotNumber == 4) {
            plot->SetPad(0.01, 0.01, 0.99, 0.23);
            plot->SetTopMargin(0.0);
            plot->SetBottomMargin(0.3);
        }
    }
    if (plotNumber == 1 && (canvasName.Index("Eta") < 0 && canvasName.Index("AbsRapidity") < 0 && canvasName.Index("DPhi") < 0)) plot->SetLogy();
    if (plotNumber == 1 && canvasName.Index("DPhiZFirstJet") > 0) plot->SetLogy();
    plot->SetLeftMargin(0.13);
    plot->SetRightMargin(0.07);
    plot->SetFillStyle(0);
    plot->Draw();
    plot->cd();
    return plot->GetCanvas();
}
#endif //defined NEW_PLOTS


void customizeMainLegend(TString canvasName, TLegend *legend, int numbOfGenerator)
{
    legend->SetFillColor(0);
    legend->SetFillStyle(1001);
    legend->SetBorderSize(1);
    legend->SetMargin(0.15);
    legend->SetTextFont(ts.defaultFont);
    legend->SetTextSize(ts.mainLegendTextSize);
    legend->SetX1(0.39);
    //    legend->SetY1(std::max(0., 0.91 - numbOfGenerator*0.07));
    legend->SetY1(std::max(0., 0.63));
    legend->SetX2(0.96);
    legend->SetY2(0.98);
}

void customizeRatioLegend(TString canvasName, TLegend *legend, int genNumb, int numbOfGenerator, int pos)
{

    legend->SetFillColor(0);
    legend->SetFillStyle(ZJetsFillStyle);
    legend->SetBorderSize(0);
    legend->SetMargin(0.3);
    legend->SetTextAlign(12);

    if(pos >= 1){
	legend->SetY2(1-0.02);
	legend->SetY1(1-0.19);
    } else{
	legend->SetY1(0.03);
	legend->SetY2(0.20);
    }
    
    if(pos == 2){
	legend->SetX1(0.2);
    } else{
	//legend->SetX1(0.15);
	legend->SetX1(0.16);
    }
    //    legend->SetX2(0.64);
    legend->SetX2(legend->GetX1() + 0.49);

    legend->SetTextFont(ts.defaultFont);
    legend->SetTextSize(ts.mainLegendTextSize);

//   if (canvasName.Index("JZB") > 0 || canvasName.Index("VisPt") > 0){
//       legend->SetY1(0.88);
//       legend->SetX1(0.16);
//       legend->SetX2(0.43);
//       legend->SetY2(0.97);
//       //legend->SetTextSize(0.06);
//       legend->SetTextSize(0.08);
//    }


   if (genNumb == numbOfGenerator) {
       //legend->SetY1(legend->GetY1() + 0.3);
       //legend->SetY2(legend->GetY2() + 0.3);
       legend->SetY1(1 - (1 - legend->GetY1()) * hratio / hbottomratio);
       legend->SetY2(1 - (1 - legend->GetY2()) * hratio / hbottomratio);
   }
}

void customizeCentral(TGraphAsymmErrors *grCentral, bool ratio)
{
    customizeCentral(grCentral, (TLegend*) NULL);
    if (ratio) grCentral->SetMarkerSize(0);
}

void customizeCentral(TGraphAsymmErrors *grCentral, TLegend *legend, TString legText)
{
    grCentral->SetLineColor(kBlack);
    grCentral->SetLineWidth(2);
    grCentral->SetMarkerStyle(20);
    grCentral->SetFillColor(12);
    gStyle->SetHatchesSpacing(1.5);
    gStyle->SetHatchesLineWidth(2);
    grCentral->SetFillStyle(3354);
    grCentral->SetMarkerColor(kBlack);

    grCentral->GetXaxis()->SetTitleOffset(ts.xTitleOffset);
    grCentral->GetXaxis()->SetTitleFont(ts.defaultFont);
    grCentral->GetXaxis()->SetTitleSize(ts.xTitleSize);
    grCentral->GetXaxis()->SetLabelFont(ts.defaultFont);
    grCentral->GetXaxis()->SetLabelSize(ts.xLabelSize);
    grCentral->GetXaxis()->SetLabelFont(ts.defaultFont);

    grCentral->GetYaxis()->SetTitleOffset(ts.mainYTitleOffset);
    grCentral->GetYaxis()->SetTitleFont(ts.defaultFont);
    grCentral->GetYaxis()->SetTitleSize(ts.mainYTitleSize);
    grCentral->GetYaxis()->SetLabelFont(ts.defaultFont);
    grCentral->GetYaxis()->SetLabelSize(ts.yLabelSize);

    grCentral->SetTitle();
    grCentral->GetXaxis()->SetTitle();
    //    if (legend) legend->AddEntry(grCentral, legText, "PLEF");
}

TGraphAsymmErrors* createGrFromHist(const TH1 *h)
{
    const TAxis* ax = h->GetXaxis();
    int nPoints = std::max(0, ax->GetLast() - ax->GetFirst() + 1);

    std::vector<double> xCoor(nPoints);
    std::vector<double> yCoor(nPoints);
    std::vector<double> xErr(nPoints);
    std::vector<double> yErr(nPoints);

    for(int i = 0; i < nPoints;  ++i){
        xCoor[i] = h->GetBinCenter(ax->GetFirst() + i);
        xErr[i]  = 0.5*h->GetBinWidth(ax->GetFirst() + i);
        yCoor[i] = h->GetBinContent(ax->GetFirst() + i);
        yErr[i] = h->GetBinError(ax->GetFirst() + i);
    }

    TGraphAsymmErrors *gr = new TGraphAsymmErrors(nPoints, &xCoor[0], &yCoor[0],
						  &xErr[0], &xErr[0], &yErr[0], &yErr[0]);

    return gr;
}

//======================================================================
// This function creates a TGraphAsymmErrors from a TGraphAsymmErrors.
// The output TGraph has all yCoor set to 1. and y low and high errors
// representing the relative errors of the input TGraph.
//======================================================================
TGraphAsymmErrors* createRatioGraph(const TGraphAsymmErrors* grCentral)
{       
    int nPoints = grCentral->GetN();
    double *xCoor = new double[nPoints];
    double *yCoor = new double[nPoints];
    double *xErr  = new double[nPoints];
    double *yErrL = new double[nPoints];
    double *yErrH = new double[nPoints];

    for (int i(0); i < nPoints; i++) {
        grCentral->GetPoint(i, xCoor[i], yCoor[i]);
        xErr[i] = grCentral->GetErrorXlow(i);
        yErrL[i] = grCentral->GetErrorYlow(i)/yCoor[i];
        yErrH[i] = grCentral->GetErrorYhigh(i)/yCoor[i];
        yCoor[i] = 1.;
    }   

    TGraphAsymmErrors *grCentralRatio = new TGraphAsymmErrors(nPoints, xCoor, yCoor, xErr, xErr, yErrL, yErrH);

    delete [] xCoor; delete [] yCoor; delete [] xErr; delete [] yErrL; delete [] yErrH; 
    return grCentralRatio;
}

//======================================================================
// This function creates a TGraphErrors from a TGraphErrors.
// The output TGraph has all yCoor set to 1. and y errors
// representing the relative errors of the input TGraph.
//======================================================================
TGraphErrors* createRatioGraph(const TGraphErrors* grCentral)
{       
    int nPoints = grCentral->GetN();
    double *xCoor = new double[nPoints];
    double *yCoor = new double[nPoints];
    double *xErr  = new double[nPoints];
    double *yErr = new double[nPoints];

    for (int i(0); i < nPoints; i++) {
        grCentral->GetPoint(i, xCoor[i], yCoor[i]);
        xErr[i] = grCentral->GetErrorX(i);
        yErr[i] = grCentral->GetErrorY(i)/yCoor[i];
        yCoor[i] = 1.;
    }   

    TGraphErrors *grCentralRatio = new TGraphErrors(nPoints, xCoor, yCoor, xErr, yErr);

    delete [] xCoor; delete [] yCoor; delete [] xErr; delete [] yErr; 
    return grCentralRatio;
}

//================================================================================
// This function creates a TGraphAsymmErrors from a TH1 and a TGraphAsymmErrors.
// The output TGraph is the ratio of the TH1 by the TGraphAsymmErrors.
//================================================================================
TGraphAsymmErrors *createGenToCentral(const TH1 *gen, const TGraphAsymmErrors *grCentral){
    if(!gen) return 0;

    const TAxis* ax = gen->GetXaxis();
    int nPoints = grCentral->GetN();
    int nPoints2 = ax->GetLast() - ax->GetFirst() + 1;
    if (nPoints != nPoints2) {
	std::cerr << "createGenToCentral() function called with inconsistent inputs. Aborts. ("
		  << __FILE__ << ":" << __LINE__ << ").\n";
	abort();
    }

    std::vector<double> xCoor(nPoints);
    std::vector<double> yCoor(nPoints);
    std::vector<double> xErr(nPoints);
    std::vector<double> yErr(nPoints);


    for (int i(0); i < nPoints; i++) {
        grCentral->GetPoint(i, xCoor[i], yCoor[i]);
        xErr[i] = grCentral->GetErrorXlow(i);
        yErr[i] = 0.;
        if (yCoor[i] != 0) {
            yErr[i]  = gen->GetBinError(ax->GetFirst() + i) / yCoor[i];
            yCoor[i] = gen->GetBinContent(ax->GetFirst() + i) / yCoor[i];
        }
    }
    TGraphAsymmErrors *grGenToCentral = new TGraphAsymmErrors(nPoints, &xCoor[0], &yCoor[0],
							      &xErr[0], &xErr[0], &yErr[0], &yErr[0]);
    return grGenToCentral;
}

TGraphAsymmErrors* createScaleSystGraph(TString sample, TString lepSel, TString variable,
					const TGraphAsymmErrors *grGenToCentral)
{
    int nPoints = grGenToCentral->GetN();
    double *xCoor    = new double[nPoints];
    double *yCoor    = new double[nPoints];
    double *xErr     = new double[nPoints];
    double *yErrUp   = new double[nPoints];
    double *yErrDown = new double[nPoints];

    TString histoDir = cfg.getS("histoDir");
    TFile *fDE;
    if (lepSel == "DE" || lepSel == "") {
        fDE = new TFile(histoDir + "/DE_13TeV_" + sample + "_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
	if(!fDE || fDE->IsZombie()){
	    std::cerr << "Fatal error. Failed to open file  " << fDE->GetName()  << ".\n";
	    abort();
	}
    }

    TFile *fDMu;
    if (lepSel == "DMu" || lepSel == "") {
        fDMu = new TFile(histoDir + "/DMu_13TeV_" + sample + "_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
	if(!fDMu || fDMu->IsZombie()){
	    std::cerr << "Fatal error. Failed to open file  " << fDMu->GetName()  << ".\n";
	    abort();
	}
    }

    TGraphAsymmErrors *grDE, *grDMu;
    if (lepSel == "DE" || lepSel == "") {
        grDE = (TGraphAsymmErrors*) fDE->Get("gen" + variable + "_scaleUnc");
	//FIXME mem. leak
	//if(!grDE) return 0;
    }
    if (lepSel == "DMu" || lepSel == "") {
        grDMu = (TGraphAsymmErrors*) fDMu->Get("gen" + variable + "_scaleUnc");
	//FIXME mem. leak
	//if(!grDMu) return 0;
    }

    // ---- this variable is used to fetch the TGraph of scale uncertainty from input file ----
    double *xMeanDMu  = new double[nPoints];
    double *yMeanDMu  = new double[nPoints];
    double *xMeanDE   = new double[nPoints];
    double *yMeanDE   = new double[nPoints];

    for (int i(0); i < nPoints; i++) {
        grGenToCentral->GetPoint(i, xCoor[i], yCoor[i]);

        xErr[i] = grGenToCentral->GetErrorXlow(i);

        yErrUp[i] = pow(grGenToCentral->GetErrorYhigh(i), 2);
        yErrDown[i] = pow(grGenToCentral->GetErrorYlow(i), 2);

        if ((lepSel == "DMu" && grDMu) || (lepSel == "DMu" && grDMu && !grDE)) {
            grDMu->GetPoint(i, xMeanDMu[i], yMeanDMu[i]);
            yErrUp[i] += pow((grDMu->GetErrorYhigh(i)/yMeanDMu[i]) * yCoor[i], 2);
            yErrDown[i] += pow((grDMu->GetErrorYlow(i)/yMeanDMu[i]) * yCoor[i], 2);
        }

        if ((lepSel == "DE" && grDE) || (lepSel == "DMu" && !grDMu && grDE)) {
            grDE->GetPoint(i, xMeanDE[i], yMeanDE[i]);
            yErrUp[i] += pow((grDE->GetErrorYhigh(i)/yMeanDE[i]) * yCoor[i], 2);
            yErrDown[i] += pow((grDE->GetErrorYlow(i)/yMeanDE[i]) * yCoor[i], 2);
        }

        if (lepSel == "" && grDE && grDMu) {
            grDMu->GetPoint(i, xMeanDMu[i], yMeanDMu[i]);
            grDE->GetPoint(i, xMeanDE[i], yMeanDE[i]);
            yErrUp[i] += pow(((grDMu->GetErrorYhigh(i) + grDE->GetErrorYhigh(i)) / (yMeanDMu[i] + yMeanDE[i])) * yCoor[i], 2);
            //yErrDown[i] += pow(((grDMu->GetErrorYhigh(i) + grDE->GetErrorYhigh(i)) / (yMeanDMu[i] + yMeanDE[i])) * yCoor[i], 2);
	    yErrDown[i] += pow(((grDMu->GetErrorYlow(i) + grDE->GetErrorYlow(i)) / (yMeanDMu[i] + yMeanDE[i])) * yCoor[i], 2);
        }

        yErrUp[i] = sqrt(yErrUp[i]);
        yErrDown[i] = sqrt(yErrDown[i]);

    }

    TGraphAsymmErrors *grScaleSyst = new TGraphAsymmErrors(nPoints, xCoor, yCoor, xErr, xErr, yErrDown, yErrUp);
    delete [] xCoor; delete [] yCoor; delete [] xErr; delete [] yErrDown; delete [] yErrUp;
    delete [] xMeanDMu; delete [] yMeanDMu; delete [] xMeanDE; delete [] yMeanDE;
    if (lepSel == "DE" || lepSel == "") {
        fDE->Close();
    }
    if (lepSel == "DMu" || lepSel == "") {
        fDMu->Close();
    }
    return grScaleSyst;

}

// --- This function is dedicated for Geneva theoretical prediction (inclusive) ---
TGraphAsymmErrors* createGenevaIncScaleSystGraph(TString sample, TString lepSel, TString variable,
						 const TGraphAsymmErrors *grGenToCentral)
{
    int nPoints = grGenToCentral->GetN();
    double *xCoor    = new double[nPoints];
    double *yCoor    = new double[nPoints];
    double *xErr     = new double[nPoints];
    double *yErrUp   = new double[nPoints];
    double *yErrDown = new double[nPoints];

    TString histoDir = cfg.getS("histoDir");
    TFile *fDE;
    if (lepSel == "DE" || lepSel == "") {
        fDE = new TFile(histoDir + "/DE_13TeV_" + sample + "_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
	if(!fDE || fDE->IsZombie()){
	    std::cerr << "Fatal error. Failed to open file  " << fDE->GetName()  << ".\n";
	    abort();
	}
    }

    TFile *fDMu;
    if (lepSel == "DMu" || lepSel == "") {
        fDMu = new TFile(histoDir + "/DMu_13TeV_" + sample + "_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
	if(!fDMu || fDMu->IsZombie()){
	    std::cerr << "Fatal error. Failed to open file  " << fDMu->GetName()  << ".\n";
	    abort();
	}
    }

    TGraphAsymmErrors *grDE, *grDMu;
    if (lepSel == "DE" || lepSel == "") {
        grDE = (TGraphAsymmErrors*) fDE->Get("gen" + variable + "_scaleUncInc");
	//FIXME mem. leak
	if(!grDE) return 0;
    }
    if (lepSel == "DMu" || lepSel == "") {
        grDMu = (TGraphAsymmErrors*) fDMu->Get("gen" + variable + "_scaleUncInc");
	//FIXME mem. leak
	if(!grDMu) return 0;
    }

    // ---- this variable is used to fetch the TGraph of scale uncertainty from input file ----
    double *xMeanDMu  = new double[nPoints];
    double *yMeanDMu  = new double[nPoints];
    double *xMeanDE   = new double[nPoints];
    double *yMeanDE   = new double[nPoints];

    for (int i(0); i < nPoints; i++) {
        grGenToCentral->GetPoint(i, xCoor[i], yCoor[i]);

        xErr[i] = grGenToCentral->GetErrorXlow(i);

        yErrUp[i] = pow(grGenToCentral->GetErrorYhigh(i), 2);
        yErrDown[i] = pow(grGenToCentral->GetErrorYlow(i), 2);

        if ((lepSel == "DMu" && grDMu) || (lepSel == "DMu" && grDMu && !grDE)) {
            grDMu->GetPoint(i, xMeanDMu[i], yMeanDMu[i]);
            yErrUp[i] += pow((grDMu->GetErrorYhigh(i)/yMeanDMu[i]) * yCoor[i], 2);
            yErrDown[i] += pow((grDMu->GetErrorYlow(i)/yMeanDMu[i]) * yCoor[i], 2);
        }

        if ((lepSel == "DE" && grDE) || (lepSel == "DMu" && !grDMu && grDE)) {
            grDE->GetPoint(i, xMeanDE[i], yMeanDE[i]);
            yErrUp[i] += pow((grDE->GetErrorYhigh(i)/yMeanDE[i]) * yCoor[i], 2);
            yErrDown[i] += pow((grDE->GetErrorYlow(i)/yMeanDE[i]) * yCoor[i], 2);
        }

        if (lepSel == "" && grDE && grDMu) {
            grDMu->GetPoint(i, xMeanDMu[i], yMeanDMu[i]);
            grDE->GetPoint(i, xMeanDE[i], yMeanDE[i]);
            yErrUp[i] += pow(((grDMu->GetErrorYhigh(i) + grDE->GetErrorYhigh(i)) / (yMeanDMu[i] + yMeanDE[i])) * yCoor[i], 2);
            yErrDown[i] += pow(((grDMu->GetErrorYhigh(i) + grDE->GetErrorYhigh(i)) / (yMeanDMu[i] + yMeanDE[i])) * yCoor[i], 2);
        }

        yErrUp[i] = sqrt(yErrUp[i]);
        yErrDown[i] = sqrt(yErrDown[i]);

    }

    TGraphAsymmErrors *grScaleSyst = new TGraphAsymmErrors(nPoints, xCoor, yCoor, xErr, xErr, yErrDown, yErrUp);
    delete [] xCoor; delete [] yCoor; delete [] xErr; delete [] yErrDown; delete [] yErrUp;
    delete [] xMeanDMu; delete [] yMeanDMu; delete [] xMeanDE; delete [] yMeanDE;
    if (lepSel == "DE" || lepSel == "") {
        fDE->Close();
    }
    if (lepSel == "DMu" || lepSel == "") {
        fDMu->Close();
    }
    return grScaleSyst;

}

// --- This function is dedicated for NNLO theoretical prediction ---
TGraphAsymmErrors* createNNLOScaleSystGraph(TString lepSel, TString variable, const TGraphAsymmErrors *grGenToCentral)
{
    int nPoints = grGenToCentral->GetN();
    double *xCoor    = new double[nPoints];
    double *yCoor    = new double[nPoints];
    double *xErr     = new double[nPoints];
    double *yErrUp   = new double[nPoints];
    double *yErrDown = new double[nPoints];

    TString histoDir = cfg.getS("histoDir");

    TFile *fnnlo[3];
    if (lepSel == "DE" || lepSel == "DMu" || lepSel == "") {
        fnnlo[0] = new TFile(histoDir + "/DMu_13TeV_DYJets_ZjNNLO_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
	if(!fnnlo[0] || fnnlo[0]->IsZombie()){
	    std::cerr << "Fatal error. Histo file for ZjNNLO prediction, "
		      << histoDir + "NNLO.root"
		      << " was not found.\n";
	}
        fnnlo[1] = new TFile(histoDir + "/DMu_13TeV_DYJets_ZjNNLO_Scaleup_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
	if(!fnnlo[1] || fnnlo[1]->IsZombie()){
	    std::cerr << "Fatal error. Histo file for scale up uncertainties of ZjNNLO prediction, "
		      << histoDir + "NNLO.root"
		      << " was not found.\n";
	}
	
        fnnlo[2] = new TFile(histoDir + "/DMu_13TeV_DYJets_ZjNNLO_Scaledn_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
	if(!fnnlo[2] || fnnlo[2]->IsZombie()){
	    std::cerr << "Fatal error. Histo file for scale down uncertainties of ZjNNLO prediction, "
		      << histoDir + "NNLO.root"
		      << " was not found.\n";
	}

    }

    TH1D *hnnlo[3];
    if (lepSel == "DE" || lepSel == "DMu" || lepSel == "") {
        hnnlo[0] = (TH1D*) fnnlo[0]->Get("gen" + variable);
        hnnlo[1] = (TH1D*) fnnlo[1]->Get("gen" + variable);
        hnnlo[2] = (TH1D*) fnnlo[2]->Get("gen" + variable);
    }

    // ---- this variable is used to fetch the TGraph of scale uncertainty from input file ----

    for (int i(0); i < nPoints; i++) {
        grGenToCentral->GetPoint(i, xCoor[i], yCoor[i]);

        xErr[i] = grGenToCentral->GetErrorXlow(i);

        yErrUp[i] = pow(grGenToCentral->GetErrorYhigh(i), 2);
        yErrDown[i] = pow(grGenToCentral->GetErrorYlow(i), 2);

        if (lepSel == "DMu" || lepSel == "DE" || lepSel == "") {
            yErrUp[i] += pow((hnnlo[1]->GetBinContent(i+1) - hnnlo[0]->GetBinContent(i+1))/hnnlo[0]->GetBinContent(i+1) * yCoor[i], 2);
            yErrDown[i] += pow((hnnlo[0]->GetBinContent(i+1) - hnnlo[2]->GetBinContent(i+1))/hnnlo[0]->GetBinContent(i+1) * yCoor[i], 2);
        }

        yErrUp[i] = sqrt(yErrUp[i]);
        yErrDown[i] = sqrt(yErrDown[i]);

    }

    TGraphAsymmErrors *grScaleSyst = new TGraphAsymmErrors(nPoints, xCoor, yCoor, xErr, xErr, yErrDown, yErrUp);
    delete [] xCoor; delete [] yCoor; delete [] xErr; delete [] yErrDown; delete [] yErrUp;
    if (lepSel == "DE" || lepSel == "DMu" || lepSel == "") {
        fnnlo[0]->Close();
        fnnlo[1]->Close();
        fnnlo[2]->Close();
    }
    return grScaleSyst;

}

// ---- create PDF systematic graph with other uncertainty in quadrature ----
TGraphAsymmErrors* createPDFSystGraph(TString sample, TString lepSel, TString variable,
				      const TGraphAsymmErrors *grGenToCentral, const TGraphAsymmErrors *grGen3ScaleSyst)
{
    int nPoints = grGenToCentral->GetN();
    double *xCoor    = new double[nPoints];
    double *yCoor    = new double[nPoints];
    double *xErr     = new double[nPoints];
    double *yErrUp   = new double[nPoints];
    double *yErrDown = new double[nPoints];

    TString histoDir = cfg.getS("histoDir");
    TFile *fDE;
    if (lepSel == "DE" || lepSel == "") {
        fDE = new TFile(histoDir + "/DE_13TeV_" + sample + "_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
    }

    TFile *fDMu;
    if (lepSel == "DMu" || lepSel == "") {
        fDMu = new TFile(histoDir + "/DMu_13TeV_" + sample + "_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
    }

    TGraphAsymmErrors *grDE, *grDMu;
    if (lepSel == "DE" || lepSel == "") {
        grDE = (TGraphAsymmErrors*) fDE->Get("gen" + variable + "_pdfUncPlain");
	//FIXME mem. leak
	if(!grDE) return 0;
    }
    if (lepSel == "DMu" || lepSel == "") {
        grDMu = (TGraphAsymmErrors*) fDMu->Get("gen" + variable + "_pdfUncPlain");
	//FIXME mem. leak
	if(!grDMu) return 0;
    }

    // ---- this variable is used to fetch the TGraph of scale uncertainty from input file ----
    double *xMeanDMu  = new double[nPoints];
    double *yMeanDMu  = new double[nPoints];
    double *xMeanDE   = new double[nPoints];
    double *yMeanDE   = new double[nPoints];

    for (int i(0); i < nPoints; i++) {
        grGenToCentral->GetPoint(i, xCoor[i], yCoor[i]);

        xErr[i] = grGenToCentral->GetErrorXlow(i);
        yErrUp[i] = pow(grGenToCentral->GetErrorYhigh(i), 2);
        yErrDown[i] = pow(grGenToCentral->GetErrorYlow(i), 2);

	if(grGen3ScaleSyst){
	    yErrUp[i]   += pow(grGen3ScaleSyst->GetErrorYhigh(i), 2);
	    yErrDown[i] += pow(grGen3ScaleSyst->GetErrorYlow(i), 2);
	}
	
        if ((lepSel == "DMu" && grDMu) || (lepSel == "" && grDMu && !grDE)) {
            grDMu->GetPoint(i, xMeanDMu[i], yMeanDMu[i]);
            yErrUp[i] += pow((grDMu->GetErrorYhigh(i)/yMeanDMu[i]) * yCoor[i], 2);
            yErrDown[i] += pow((grDMu->GetErrorYlow(i)/yMeanDMu[i]) * yCoor[i], 2);
        }

        if ((lepSel == "DE" && grDE) || (lepSel == "" && !grDMu && grDE)) {
            grDE->GetPoint(i, xMeanDE[i], yMeanDE[i]);
            yErrUp[i] += pow((grDE->GetErrorYhigh(i)/yMeanDE[i]) * yCoor[i], 2);
            yErrDown[i] += pow((grDE->GetErrorYlow(i)/yMeanDE[i]) * yCoor[i], 2);
        }

        if (lepSel == "" && grDE && grDMu) {
            grDMu->GetPoint(i, xMeanDMu[i], yMeanDMu[i]);
            grDE->GetPoint(i, xMeanDE[i], yMeanDE[i]);
            yErrUp[i] += pow(((grDMu->GetErrorYhigh(i) + grDE->GetErrorYhigh(i)) / (yMeanDMu[i] + yMeanDE[i])) * yCoor[i], 2);
            yErrDown[i] += pow(((grDMu->GetErrorYhigh(i) + grDE->GetErrorYhigh(i)) / (yMeanDMu[i] + yMeanDE[i])) * yCoor[i], 2);
        }

        yErrUp[i] = sqrt(yErrUp[i]);
        yErrDown[i] = sqrt(yErrDown[i]);

    }

    TGraphAsymmErrors *grPDFSyst = new TGraphAsymmErrors(nPoints, xCoor, yCoor, xErr, xErr, yErrDown, yErrUp);
    delete [] xCoor; delete [] yCoor; delete [] xErr; delete [] yErrDown; delete [] yErrUp;
    delete [] xMeanDMu; delete [] yMeanDMu; delete [] xMeanDE; delete [] yMeanDE;

    if (lepSel == "DE" || lepSel == "") {
        fDE->Close();
    }
    if (lepSel == "DMu" || lepSel == "") {
        fDMu->Close();
    }

    return grPDFSyst;
}

//// ---- create PDF systematical graph without other uncertainty in quadrature ----
//TGraphAsymmErrors* createPDFSystGraph(TString lepSel, TString variable, const TGraphAsymmErrors *grGenToCentral)
//{
//    int nPoints = grGenToCentral->GetN();
//    double *xCoor    = new double[nPoints];
//    double *yCoor    = new double[nPoints];
//    double *xErr     = new double[nPoints];
//    double *yErrUp   = new double[nPoints];
//    double *yErrDown = new double[nPoints];
//
//    TString histoDir = cfg.getS("histoDir");
//    TFile *fDE;
//    if (lepSel == "DE" || lepSel == "") {
//      //        fDE = new TFile("HistoFilesUnc/DE_13TeV_DYJets_UNFOLDING_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
//        fDE = new TFile(histoDir + "/DE_13TeV_DYJets_UNFOLDING_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
//    }
//
//    TFile *fDMu;
//    if (lepSel == "DMu" || lepSel == "") {
//        fDMu = new TFile(histoDir + "/DMu_13TeV_DYJets_UNFOLDING_TrigCorr_1_Syst_0_JetPtMin_30_JetEtaMax_24.root");
//    }
//
//    TGraphAsymmErrors *grDE, *grDMu;
//    if (lepSel == "DE" || lepSel == "") {
//        grDE = (TGraphAsymmErrors*) fDE->Get("gen" + variable + "_pdfUncPlain");
//    }
//    if (lepSel == "DMu" || lepSel == "") {
//        grDMu = (TGraphAsymmErrors*) fDMu->Get("gen" + variable + "_pdfUncPlain");
//    }
//
//    // ---- this variable is used to fetch the TGraph of scale uncertainty from input file ----
//    double *xMeanDMu  = new double[nPoints];
//    double *yMeanDMu  = new double[nPoints];
//    double *xMeanDE   = new double[nPoints];
//    double *yMeanDE   = new double[nPoints];
//
//    for (int i(0); i < nPoints; i++) {
//        grGenToCentral->GetPoint(i, xCoor[i], yCoor[i]);
//
//        xErr[i] = grGenToCentral->GetErrorXlow(i);
//
//        if (lepSel == "DMu") {
//            grDMu->GetPoint(i, xMeanDMu[i], yMeanDMu[i]);
//            yErrUp[i] += pow((grDMu->GetErrorYhigh(i)/yMeanDMu[i]) * yCoor[i], 2);
//            yErrDown[i] += pow((grDMu->GetErrorYlow(i)/yMeanDMu[i]) * yCoor[i], 2);
//        }
//
//        if (lepSel == "DE") {
//            grDE->GetPoint(i, xMeanDE[i], yMeanDE[i]);
//            yErrUp[i] += pow((grDE->GetErrorYhigh(i)/yMeanDE[i]) * yCoor[i], 2);
//            yErrDown[i] += pow((grDE->GetErrorYlow(i)/yMeanDE[i]) * yCoor[i], 2);
//        }
//
//        if (lepSel == "") {
//            grDMu->GetPoint(i, xMeanDMu[i], yMeanDMu[i]);
//            grDE->GetPoint(i, xMeanDE[i], yMeanDE[i]);
//            yErrUp[i] += pow(((grDMu->GetErrorYhigh(i) + grDE->GetErrorYhigh(i)) / (yMeanDMu[i] + yMeanDE[i])) * yCoor[i], 2);
//            yErrDown[i] += pow(((grDMu->GetErrorYhigh(i) + grDE->GetErrorYhigh(i)) / (yMeanDMu[i] + yMeanDE[i])) * yCoor[i], 2);
//        }
//
//    }
//
//    TGraphAsymmErrors *grPDFSyst = new TGraphAsymmErrors(nPoints, xCoor, yCoor, xErr, xErr, yErrDown, yErrUp);
//    delete [] xCoor; delete [] yCoor; delete [] xErr; delete [] yErrDown; delete [] yErrUp;
//    delete [] xMeanDMu; delete [] yMeanDMu; delete [] xMeanDE; delete [] yMeanDE;
//
//    if (lepSel == "DE" || lepSel == "") {
//        fDE->Close();
//    }
//    if (lepSel == "DMu" || lepSel == "") {
//        fDMu->Close();
//    }
//
//    return grPDFSyst;
//}

void customizeRatioGraph(TH1 *hAxis, TGraphAsymmErrors *gen,
			 TGraphAsymmErrors *gScale, TGraphAsymmErrors *gPDF,
			 int genNum, TString yTitle, int numbOfGenerator,
			 TLegend *legend){

    double minRatioY = cfg.getD("minRatioYUnf", 0.2);
    double maxRatioY = cfg.getD("maxRatioYUnf", 1.8);
    if(genNum == 3) extendAxis(minRatioY, maxRatioY); //only the 3rd gen (GE) needs y-axis extension

    if(hAxis){
	hAxis->GetYaxis()->SetRangeUser(minRatioY, maxRatioY);
	hAxis->GetYaxis()->SetNdivisions(507);
	hAxis->GetYaxis()->SetLabelFont(ts.defaultFont);
	hAxis->GetYaxis()->SetLabelSize(ts.yLabelSize);
	hAxis->GetYaxis()->SetTitle(yTitle);
	hAxis->GetYaxis()->SetTitleFont(ts.defaultFont);
	hAxis->GetYaxis()->SetTitleSize(ts.ratioYTitleSize);
       	hAxis->GetYaxis()->SetTitleOffset(ts.ratioYTitleOffset * ratioYTitleOffsetCorr);
	hAxis->GetYaxis()->CenterTitle();
	hAxis->SetTitle("");
	
	if (genNum == numbOfGenerator) {
	    hAxis->GetXaxis()->SetLabelFont(ts.defaultFont);
	    hAxis->GetXaxis()->SetLabelSize(ts.xLabelSize);
	    hAxis->GetXaxis()->SetLabelOffset(0.02);
	    hAxis->GetXaxis()->SetTitleFont(ts.defaultFont);
	    hAxis->GetXaxis()->SetTitleSize(ts.xTitleSize);
	    hAxis->GetXaxis()->SetTitleOffset(ts.xTitleOffset / hbottomratio);
	    //	    hAxis->GetYaxis()->SetTitleOffset(ts.ratioYTitleOffset * bottomYTitleOffsetCorr);
       	} else{ //disable display of other axes than the bottom plot:
	    hAxis->GetXaxis()->SetLabelSize(0);
	    hAxis->GetXaxis()->SetTitleSize(0);
	}
    }

    if(gen){
	gen->SetFillColor(ZJetsFillColor[genNum-1]);
	gen->SetFillStyle(ZJetsFillStyle);
	gen->SetLineColor(ZJetsLineColor[genNum-1]);
	gen->SetLineWidth(2);
	gen->SetMarkerColor(ZJetsMarkerColor[genNum-1]);
	gen->SetMarkerStyle(ZJetsMarkerStyle[genNum-1]);
    }

    if(gScale){
	gScale->SetFillStyle(ZJetsFillStyle);
	gScale->SetLineColor(ZJetsLineColor[genNum-1]);
	gScale->SetLineWidth(2);
	gScale->SetFillColor(ZJetsScaleFillColor[genNum-1]);
    }
    
    if(gPDF){
	gPDF->SetFillStyle(0);
	gPDF->SetLineColor(ZJetsLineColor[genNum-1]);
	gPDF->SetLineWidth(2);
    }
    
    if (legend) {
        TLegendEntry *leEntry;
        TLegendEntry *statEntry;
        TLegendEntry *pdfEntry;
	int nentries = 1;
	if(gScale) ++nentries;
	if(gPDF) ++nentries;
	
	//   if(/*genNum == 3 ||*/ genNum == 1) {
	//	legend->SetX2(0.64);
	legend->SetNColumns(3);
	//statEntry = legend->AddEntry(gen, "Stat", "f");
	TString l = "Stat";
	if(!gScale && !gPDF) l += " unc.";
	statEntry = legend->AddEntry((TObject*)0, l, "f");
	statEntry->SetFillStyle(ZJetsFillStyle);
	statEntry->SetFillColor(ZJetsFillColor[genNum-1]);
	statEntry->SetLineColor(ZJetsFillColor[genNum-1]);

	l = "#oplus theo";
	if(gScale){
	    //leEntry = legend->AddEntry(gScale, "#oplus Theory", "f");
	    if(!gPDF) l += " unc.                            ";
	    leEntry = legend->AddEntry((TObject*)0, l, "f");
	    leEntry->SetFillColor(ZJetsScaleFillColor[genNum-1]);
	    leEntry->SetFillStyle(ZJetsFillStyle);
	    leEntry->SetLineColor(ZJetsScaleFillColor[genNum-1]);
	} 

	if(gPDF){
            pdfEntry = legend->AddEntry(gPDF, "#oplus PDF #oplus #alpha_{s} unc.   ", "f");
            pdfEntry->SetFillStyle(0);
	}   
	
	
	//        }
	//        else {
	//            //leEntry = legend->AddEntry(gen, "Stat unc.", "f");
	//            leEntry = legend->AddEntry((TObject*)0, "Stat unc.", "f");
	//            leEntry->SetFillColor(ZJetsFillColor[genNum-1]);
	//            leEntry->SetFillStyle(ZJetsFillStyle);
	//            leEntry->SetLineColor(ZJetsFillColor[genNum-1]);
	//        }
    }
}


//void customizeRatioGraph(TH1 *hAxis, TGraphAsymmErrors *gen, TGraphAsymmErrors *gPDF,
//			 int genNum, TString yTitle, int numbOfGenerator, TLegend *legend)
//{
//    double minRatioY = cfg.getD("minRatioYUnf", 0.2);
//    double maxRatioY = cfg.getD("maxRatioYUnf", 1.8);
//
//    if(hAxis){
//	hAxis->GetYaxis()->SetRangeUser(minRatioY, maxRatioY);
//	hAxis->GetYaxis()->SetNdivisions(507);
//	hAxis->GetYaxis()->SetLabelFont(ts.defaultFont);
//	hAxis->GetYaxis()->SetLabelSize(ts.yLabelSize);
//	hAxis->GetYaxis()->SetTitle(yTitle);
//	hAxis->GetYaxis()->SetTitleFont(ts.defaultFont);
//	hAxis->GetYaxis()->SetTitleSize(ts.ratioYTitleSize);
//	hAxis->GetYaxis()->SetTitleOffset(3.);
//	hAxis->GetYaxis()->CenterTitle();
//    }
//
//    if(gen){
//	gen->SetFillColor(ZJetsFillColor[genNum-1]);
//	gen->SetFillStyle(ZJetsFillStyle);
//	gen->SetLineColor(ZJetsLineColor[genNum-1]);
//	gen->SetLineWidth(2);
//	gen->SetMarkerColor(ZJetsLineColor[genNum-1]);
//	gen->SetMarkerStyle(ZJetsMarkerStyle[genNum-1]);
//    }
//
//    if(gPDF){
//	gPDF->SetFillStyle(0);
//	gPDF->SetLineColor(ZJetsLineColor[genNum-1]);
//	gPDF->SetLineWidth(2);
//    }
//    
//    if (genNum == numbOfGenerator && hAxis) {
//        hAxis->GetXaxis()->SetLabelFont(ts.defaultFont);
//        hAxis->GetXaxis()->SetLabelSize(ts.xLabelSize);
//        hAxis->GetXaxis()->SetTitleFont(ts.defaultFont);
//        hAxis->GetXaxis()->SetTitleSize(ts.xTitleSize);
//        hAxis->GetXaxis()->SetTitleOffset(3.0);
//    }
//    else if(hAxis){
//        hAxis->GetXaxis()->SetTitle();
//    }
//
//    if (legend) {
//        TLegendEntry *leEntry;
//        //leEntry = legend->AddEntry(gen, "Stat unc.", "f");
//        leEntry = legend->AddEntry((TObject*)0, "Stat unc.", "f");
//        leEntry->SetFillColor(ZJetsFillColor[genNum-1]);
//        leEntry->SetFillStyle(ZJetsFillStyle);
//        leEntry->SetLineColor(ZJetsFillColor[genNum-1]);
//    }
//}


void customizeGenHist(TH1 *gen, int genNumb, TLegend *legend, TString legText)
{

    //--- Customize gen Sherpa ---
    gen->SetFillColor(ZJetsFillColor[genNumb-1]);
    gen->SetFillStyle(ZJetsFillStyle);
    gen->SetLineColor(ZJetsLineColor[genNumb-1]);
    gen->SetLineWidth(2);
    gen->SetMarkerColor(ZJetsMarkerColor[genNumb-1]);
    gen->SetMarkerStyle(ZJetsMarkerStyle[genNumb-1]);
    //TLegendEntry *le = legend->AddEntry(gen, legText, "pefl");
    TLegendEntry *le = legend->AddEntry(gen, legText, "pfl");
    le->SetFillColor(ZJetsFillColor[genNumb-1]);
    le->SetFillStyle(ZJetsFillStyle);
    le->SetLineColor(ZJetsLineColor[genNumb-1]);
    le->SetMarkerColor(ZJetsMarkerColor[genNumb-1]);
    le->SetMarkerStyle(ZJetsMarkerStyle[genNumb-1]);
}


void configYaxis(TH1 *grCentralSyst, TH1 *gen1, TH1 *gen2, TH1 *gen3)
{
    //--- Configure Y axis of the plot ---
    double minimumToPlot = std::numeric_limits<double>::max();
    if(grCentralSyst) minimumToPlot = TMath::Min(minimumToPlot, grCentralSyst->GetMinimum());
    if (gen1) minimumToPlot = TMath::Min(minimumToPlot, gen1->GetBinContent(gen1->GetMinimumBin()));
    if (gen2) minimumToPlot = TMath::Min(minimumToPlot, gen2->GetBinContent(gen2->GetMinimumBin()));
    if (gen3) minimumToPlot = TMath::Min(minimumToPlot, gen3->GetBinContent(gen3->GetMinimumBin()));

    double maximumToPlot = -std::numeric_limits<double>::max();
    if(grCentralSyst) maximumToPlot = TMath::Max(maximumToPlot, grCentralSyst->GetMaximum());
    if (gen1) maximumToPlot = TMath::Max(maximumToPlot, gen1->GetBinContent(gen1->GetMaximumBin()));
    if (gen2) maximumToPlot = TMath::Max(maximumToPlot, gen2->GetBinContent(gen2->GetMaximumBin()));
    if (gen3) maximumToPlot = TMath::Max(maximumToPlot, gen3->GetBinContent(gen3->GetMaximumBin()));

    //    if(grCentralSyst) grCentralSyst->GetYaxis()->SetRangeUser(0.2*minimumToPlot, 5*maximumToPlot);
    //if (TString(grCentralSyst->GetName()).Contains("Eta")) {
    //    grCentralSyst->GetYaxis()->SetRangeUser(0.001, 1.4*maximumToPlot);
    //}
}

//void configXaxis(TGraphAsymmErrors *grCentralSyst, TH1D *gen1)
void configXaxis(TH1 *grCentralSyst, TH1 *gen1, TString variable)
{
    //--- Configure X axis of the plot ---
    //double minX, tmp;
    //double maxX;
    //TString variable = gen1->GetName();
    //grCentralSyst->GetPoint(firstBin, minX, tmp);
    //grCentralSyst->GetPoint(grCentralSyst->GetN()-1, maxX, tmp);
    //minX -= grCentralSyst->GetErrorXlow(firstBin); 
    //maxX += grCentralSyst->GetErrorXhigh(grCentralSyst->GetN()-1);
    if(grCentralSyst){	
	if (variable.Index("ZNGoodJets_Zexc") >= 0) {	
	    std::cout << __FILE__ << ":" << __LINE__ 
		      << ". Range of ZNGoodJets_Zexc x-axis is being modified.!\n";
	    //grCentralSyst->GetXaxis()->Set(maxX-minX, minX, maxX);
	    //	grCentralSyst->GetXaxis()->SetRangeUser(-0.5, 4.5);
	    grCentralSyst->GetXaxis()->SetBinLabel(1, "= 0");
	    grCentralSyst->GetXaxis()->SetBinLabel(2, "= 1");
	    grCentralSyst->GetXaxis()->SetBinLabel(3, "= 2");
	    grCentralSyst->GetXaxis()->SetBinLabel(4, "= 3");
	    grCentralSyst->GetXaxis()->SetBinLabel(5, "= 4");
	    grCentralSyst->GetXaxis()->SetBinLabel(6, "= 5");
	    grCentralSyst->GetXaxis()->SetBinLabel(7, "= 6");
	    //grCentralSyst->GetXaxis()->SetBinLabel(8, "= 7");
	    //     grCentralSyst->GetXaxis()->SetBinLabel(9, "= 8");
	    //-->
	    //        grCentralSyst->GetXaxis()->SetLabelSize(0.18);
	    //        grCentralSyst->GetXaxis()->SetLabelOffset(0.01);
	}  else if (variable.Index("ZNGoodJets_Zinc") >= 0) {
	    std::cout << __FILE__ << ":" << __LINE__ 
		      << ". Range of ZNGoodJets_Zexc x-axis is being modified.!\n";
	    //	grCentralSyst->GetXaxis()->SetRangeUser(-0.5, 4.5);
	    //grCentralSyst->GetXaxis()->Set(maxX-minX, minX, maxX);
	    grCentralSyst->GetXaxis()->SetBinLabel(1, "#geq 0");
	    grCentralSyst->GetXaxis()->SetBinLabel(2, "#geq 1");
	    grCentralSyst->GetXaxis()->SetBinLabel(3, "#geq 2");
	    grCentralSyst->GetXaxis()->SetBinLabel(4, "#geq 3");
	    grCentralSyst->GetXaxis()->SetBinLabel(5, "#geq 4");
	    grCentralSyst->GetXaxis()->SetBinLabel(6, "#geq 5");
	    grCentralSyst->GetXaxis()->SetBinLabel(7, "#geq 6");
	    //        grCentralSyst->GetXaxis()->SetBinLabel(8, "#geq 7");
	    // -->
	    //	    grCentralSyst->GetXaxis()->SetLabelSize(0.18);
	    //	    grCentralSyst->GetXaxis()->SetLabelOffset(0.01);
	}
	//grCentralSyst->GetXaxis()->SetRangeUser(minX, maxX);
    }
    TString xtitle;
    if(gen1){
	xtitle = gen1->GetXaxis()->GetTitle();
    } else if(grCentralSyst){	
	xtitle = grCentralSyst->GetXaxis()->GetTitle();
    }

    if (xtitle.Index("^{gen}") >= 0) xtitle = xtitle.ReplaceAll("^{gen}","");
    if (xtitle.Index("H_{T}") >= 0) {
	//	    TString njets;
	//	    if (variable.Index("Zinc1jet") >= 0) njets = "1";
	//	    else if (variable.Index("Zinc2jet") >= 0) njets = "2";
	//	    else if (variable.Index("Zinc3jet") >= 0) njets = "3";
	//	    else if (variable.Index("Zinc4jet") >= 0) njets = "4";
	//	    else if (variable.Index("Zinc5jet") >= 0) njets = "5";
	//	    else if (variable.Index("Zinc6jet") >= 0) njets = "6";
	//	    else if (variable.Index("Zinc7jet") >= 0) njets = "7";
	//	    else if (variable.Index("Zinc8jet") >= 0) njets = "8";
	xtitle = "H_{T} [GeV]";
    }

    xtitle.ReplaceAll("p_{T} balance",  "p_{T}^{bal}");
    if(xtitle.BeginsWith("JZB")) xtitle = "JZB [GeV]";

    if(gen1) gen1->GetXaxis()->SetTitle(xtitle);
    if(grCentralSyst) grCentralSyst->GetXaxis()->SetTitle(xtitle);
   
//   if(variable.Index("ZPt_") >= 0){
//       TAxis* a = grCentralSyst->GetXaxis();
//       a->SetRangeUser(10., a->GetBinUpEdge(a->GetNbins()));
//   }

//   if(grCentralSyst){
//       if(gen1) grCentralSyst->GetXaxis()->SetTitle(xtitle);
//       grCentralSyst->GetXaxis()->SetTitleFont(ts.defaultFont);
//       grCentralSyst->GetXaxis()->SetTitleSize(ts.xTitleSize);
//       grCentralSyst->GetXaxis()->SetTitleOffset(1.2/hbottomratio);
//       grCentralSyst->GetXaxis()->SetLabelSize(ts.xLabelSize);
//       grCentralSyst->GetXaxis()->SetTitleFont(ts.defaultFont);
//   }
    //    if(grCentralSyst) grCentralSyst->GetXaxis()->SetTitleSize(0.12);
    
    //-----------------------------------------

}

std::string getYaxisTitle(bool doNormalized, const TH1 *gen1)
{
    std::string title = "";
    std::string xtitle = gen1->GetXaxis()->GetTitle();
    std::string shortVar = xtitle.substr(0, xtitle.find(" "));
    std::string unit = "";
    if (xtitle.find("^{gen}") != std::string::npos) {
        xtitle.replace(xtitle.find("^{gen}"),6,"");
    }
    if (xtitle.find("[") != std::string::npos){
        size_t begin = xtitle.find("[") + 1;
        unit = xtitle.substr(begin);
        unit = unit.substr(0, unit.find("]"));
    }
    title = "d#sigma/d" + shortVar;
    if (doNormalized) {
        title = "1/#sigma " + title;
    }
    else {
        title += "  [pb";
        if (unit != "" ) title += "/" + unit;
        title += "]";
    }
    return title;
}




void makeCrossSectionPlot(const char* variable, const char* ref){
    if(!variable){
	for(unsigned i = 0 ; i < NVAROFINTERESTZJETS; ++i){
	    const char* v = VAROFINTERESTZJETS[i].name.Data();
	    if(v && v[0]){
		makeCrossSectionPlot(v, ref);
	    }
	}
	return;
    }


    std::string lepSel = cfg.getS("lepSel");
    std::string algo = cfg.getS("algo");
    bool doNormalized = false;
    int jetPtMin = cfg.getI("jetPtMin");
    int jetEtaMax = cfg.getI("jetEtaMax");

    std::vector<std::string> predictions = cfg.getVS("predictions");

    std::vector<std::string> tmp;
    for(auto v: predictions){
	std::cerr << "====> " << v << "\n";
	if(v.size() != 0) tmp.push_back(v);
    }
    std::swap(tmp, predictions);

    TString unfCfgFile = cfg.getS("unfConf");
    static SectionedConfig unfCfg;
    unfCfg.read(unfCfgFile, true);
  
    TString section = TString::Format("%s_%s", lepSel.size() > 0 ? lepSel.c_str() : "DMu", variable);
    int nFirstBinsToSkip = unfCfg.get(section.Data(), "nFirstBinsToSkip", 0);
    int nLastBinsToSkip = unfCfg.get(section.Data(), "nLastBinsToSkip", 0);

    TString path;
    TString dataHistName;
    TString covPrefix; //prefix in covariance matrix name
    if(lepSel.size() > 0){ //single channel
	std::string unfoldDir = cfg.getS("unfoldDir");
	path = getUnfoldedFileName(unfoldDir, lepSel, variable, algo,
				   jetPtMin, jetEtaMax, "_MGPYTHIA6_", doNormalized);
	dataHistName = "UnfDataCentral";
	covPrefix = "";
    } else{//combined result
	std::string combDir = cfg.getS("combDir");
	bool diagXChanCov  = cfg.getB("diagXChanCov", true);
	bool fullXChanCov  = cfg.getB("fullXChanCov", true);
	bool fullSChanCov  = cfg.getB("fullSChanCov", true);
	bool modifiedSWA   = cfg.getB("modifiedSWA", true);
	path = getCombinedFileName(combDir, variable, algo,
				   diagXChanCov, fullXChanCov, fullSChanCov, modifiedSWA,
				   jetPtMin,jetEtaMax, "_MGPYTHIA6_",
				   doNormalized);
	dataHistName = "CombDataCentral";
	covPrefix = "Comb";
    }
    
    Ssiz_t p = path.Index(TRegexp("[^/]*$"));
    TString dir = path(0, p-1);
    TString fname = path(p, path.Length());

    TH1* hRef = 0;
    TH2* hCov = 0;
    
    double lumi = 0;

    if(TString(ref).CompareTo("data", TString::kIgnoreCase) == 0){
	TFile* fdata = TFile::Open(path + ".root");
	if(fdata && fdata->IsZombie()){
	    delete fdata;
	    fdata = 0;
	}
	
	if(!fdata){
	    std::cerr << "Fatal error. Failed to open unfolded data file "
		      << (path + ".root") << "\n";
	    abort();
	}
	
	hRef = (TH1*) fdata->Get(dataHistName);
	if(!hRef){
	    std::cerr << "Fatal error. Histogram UnfDataCentral  was not found in file "
		      << fdata->GetName() << "\n";
	    abort();
	}
	hRef->SetZTitle("Data");
		
	hCov = (TH2*) fdata->Get(covName(kTotSys, covPrefix));
	if(!hCov){
	    std::cerr << "Fatal error. Histogram " << covName(kTotSys, covPrefix) 
		      << "  was not found in file "
		      << fdata->GetName() << "\n";
	    abort();
	}
	hRef->SetDirectory(0);
	hCov->SetDirectory(0);

	 TH1* Lumi = 0;
	 if(fdata) fdata->GetObject("Lumi", Lumi);
	 if(Lumi) lumi = Lumi->GetBinContent(1);
	 else {
	     cerr << "Error: Lumi histogram was not found.\n";
	     return;
	 }
	delete fdata;
    } else{
	hRef = (TH1*)getGenHistos(std::vector<std::string>(1, ref), lepSel.c_str(), variable, true, true)[0];
	if(hRef){
	    //to avoid double rescale (see Scale(fac) below and in makeCrossSectionPlot method)
	    //FIXME: use a smart pointer
	    hRef = (TH1*) hRef->Clone();
	    hRef->SetZTitle(getLegendGen(ref));
	} 
    }

    if(hRef){
	double fac = cfg.getF(TString("scale_") + ref, 1.);
	if(fac != 1.){
	    std::cout << "Scaling " << ref << " by factor " << fac << std::endl;
	    hRef->Scale(fac);
	}
	
	int upperBin = hRef->GetXaxis()->GetNbins() -  nLastBinsToSkip;
	int lowerBin = 1 + nFirstBinsToSkip;
	hRef->GetXaxis()->SetRange(lowerBin, upperBin);
	TCanvas *crossSectionPlot = makeCrossSectionPlot(lepSel, lumi, variable, doNormalized,
							 hRef, hCov,
							 predictions);
	saveCanvas(crossSectionPlot, dir, fname);
    } else{
	std::cerr << "Warning: no reference histogram was found for variable " << variable << ", channel "
		  << lepSel << "\n";
    }
    delete hRef;
    delete hCov;
}

TH1* makeCrossSectionHist(TH1* hGenDYJets, double integratedLumi)
{
    TH1 *hGenCrossSection = (TH1*) hGenDYJets->Clone();
    //--- divide by luminosity ---
    hGenCrossSection->Scale(1./integratedLumi);

    int nBins = hGenCrossSection->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
	double binWidth = hGenCrossSection->GetBinWidth(i);
	hGenCrossSection->SetBinContent(i, hGenCrossSection->GetBinContent(i)*1./binWidth);
	hGenCrossSection->SetBinError(i, hGenCrossSection->GetBinError(i)*1./binWidth);
    }

    return hGenCrossSection;
}


//======


TCanvas* makeCrossSectionPlot(TString lepSel, double lumi, TString variable,
			      bool doNormalized, TH1* hStat, TH2* hCovSyst,
			      std::vector<std::string> gens){
//    if(gens.size() > 3){
//	std::cerr << "Warning. Maxium three generator comparison is supported. Only the first three will be considered."
//		  << " (" << __FILE__ << ":" << __LINE__ << ").\n";
//	gens.resize(3);
//    }


    double savedRatioYaxisExtend = ratioYaxisExtend;
    if(variable.BeginsWith("JZB")) ratioYaxisExtend = 1.67;
    if(variable == "VisPt_Zinc1jetQun" || variable == "VisPt_Zinc2jetQun") ratioYaxisExtend = 1.67;
    
    gStyle->SetOptStat(0);

    bool isPrel = cfg.getB("preliminaryTag", true);    
    std::vector<TH1*> hGens = getGenHistos(gens, lepSel, variable);

    std::vector<TH1*> tmp1;
    std::vector<std::string> tmp2;
    
    struct GenRcd {
	GenRcd(const std::string& n, TH1* h, TGraphAsymmErrors* ratio = 0,
	       TGraphAsymmErrors* pdf = 0, TGraphAsymmErrors* scale = 0):
	    name(n), histo(h), dataRatio(ratio), dataRatioPDFSyst(pdf),
	    dataRatioScaleSyst(scale){}
	GenRcd(const GenRcd& a){
	    this->name = a.name;
	    this->histo = a.histo;
	    this->dataRatio = a.dataRatio;
	    this->dataRatioPDFSyst = a.dataRatioPDFSyst;
	    this->dataRatioScaleSyst = a.dataRatioScaleSyst;
	}
	std::string name;
	TH1* histo;
	TGraphAsymmErrors* dataRatio;	
	TGraphAsymmErrors* dataRatioPDFSyst;
	TGraphAsymmErrors* dataRatioScaleSyst;
    };

    //Predictions to show in the main frame
    std::vector<GenRcd> mainFrame;
    //Predictions to show on ratio plots. Several
    //predictions can be superimposed in the same frame
    // [iframe][ipred]
    std::vector<std::vector<GenRcd> > ratioFrames;

    std::map<std::string, std::vector<GenRcd> > genRatiosToSuperimposed;

    int ipred = -1;
    for(auto g: hGens){
	++ipred;
	if(!g) continue;
	if(!alignRanges(hStat->GetXaxis(), g->GetXaxis())){
	    std::cerr << "Fatal error. X-axis range of " << hStat->GetName() << " and "
		      << g->GetName() << " are inconsistent. Aborts. ("
		      << __FILE__ << ":" << __LINE__ << ")\n";
	    std::cerr << "\t" << hStat->GetName() << ": " << hStat->GetXaxis()->GetNbins() << "bins:";
	    for(int i = 1; i <= hStat->GetNbinsX() + 1; ++i){
		std::cerr << "\t" << hStat->GetXaxis()->GetBinLowEdge(i);
	    }
	    std::cerr << "\n\t" << g->GetName() << ": " << g->GetXaxis()->GetNbins() << "bins:";
	    for(int i = 1; i <= g->GetXaxis()->GetNbins() + 1; ++i){
		std::cerr << "\t" << g->GetXaxis()->GetBinLowEdge(i);
	    }

	    abort();
	}
	    
	tmp1.push_back(g);
	tmp2.push_back(gens[ipred]);

	g->SetName(TString("gen") + gens[ipred]);

	if(g) g->SetZTitle(getLegendGen(gens[ipred].c_str()));
	//if(TString(gens[ipred]).BeginsWith("DYJets_GE")){
	//    g->Scale(2.);
	//}
	double fac = cfg.getF(TString("scale_") + gens[ipred], 1.);
	if(fac != 1.){
	    std::cout << "Scaling " << gens[ipred] << " by factor " << fac << std::endl;
	    g->Scale(fac);
	}
	std::string display = cfg.getS(TString(gens[ipred]) + "_display");
	if(TString(display.c_str()).BeginsWith("on-top-of_")){
	    std::string altGen = display.substr(strlen("on-top-of_"));
	    if(std::find(gens.begin(), gens.end(), altGen) == gens.end()){
		std::cerr << "Warning: prediction " << altGen 
			  << " indicated in the parameter " << display
			  << " is not listed by the parameter predictions. "
			  << "None of the two predictions will be displayed.\n";
	    } else {
		genRatiosToSuperimposed[altGen].push_back(GenRcd(gens[ipred], g));
	    }
		    
	} else{
	    mainFrame.push_back(GenRcd(gens[ipred], g));
	}
    }

    for(auto g: mainFrame){
	//FIXME: should go in the configuration file
//	if(TString(g.name).BeginsWith("DYJets_GE") 
//	 && (variable == "VisPt_Zinc2jetQun" || variable == "VisPt_Zinc3jetQun")){
//	    continue;
//	}

	ratioFrames.push_back(std::vector<GenRcd>(1, g));
	auto it = genRatiosToSuperimposed.find(g.name);
	if(it != genRatiosToSuperimposed.end()){
	    for(auto r: it->second){
		ratioFrames.back().push_back(r);
	    }
	}
    }

    //drops missing predictions (we have ZjNNLO for only few distributions):
    hGens = tmp1;
    gens = tmp2;

    //    int nRatioPlots = ratioFrames.size(); 
    
    //    hGens.resize(3, 0);

    //--- TGraph for data central value and stat unc.:
    TGraphAsymmErrors *grCentralStat = createGrFromHist(hStat);
    grCentralStat->SetName("gr" + variable + "CentralStatError");
    TGraphAsymmErrors *grCentralStatRatio = createRatioGraph(grCentralStat);
    
    TGraphAsymmErrors *grCentralSyst = 0;
    TGraphAsymmErrors *grCentralSystRatio = 0;
    TH1* hSyst = (TH1*) hStat->Clone("hSyst");
    if(hCovSyst){
	int nBins = hSyst->GetNbinsX();
	for (int i = 1; i <= nBins; ++i) {
	    hSyst->SetBinError(i, sqrt(pow(hStat->GetBinError(i), 2) + hCovSyst->GetBinContent(i, i)));
	}
    }
    grCentralSyst = createGrFromHist(hSyst);
    grCentralSystRatio = createRatioGraph(grCentralSyst);
    grCentralSyst->SetName("gr" + variable + "CentralTotError"); 
  
    std::vector<TGraphAsymmErrors*> grGen1ToCentral(hGens.size(), 0);
    std::vector<TGraphAsymmErrors*> grGen1ScaleSyst(hGens.size(), 0);
    std::vector<TGraphAsymmErrors*> grGen1PDFSyst(hGens.size(), 0);

    //create ratio plots
    for(auto& f: ratioFrames){
	for(auto& g: f){
	    g.dataRatio = createGenToCentral(g.histo, grCentralStat);
	    int showSys = cfg.getI(TString::Format("%s_unc", g.name.c_str()));
	    if(showSys == 1){
		g.dataRatioScaleSyst = createScaleSystGraph(g.name, lepSel, variable, g.dataRatio);
		g.dataRatioPDFSyst   = createPDFSystGraph(g.name, lepSel, variable, g.dataRatio, g.dataRatioScaleSyst);
		
	    } else if(showSys == 2){
		g.dataRatioScaleSyst = createNNLOScaleSystGraph(lepSel, variable, g.dataRatio);
	    } else if(showSys == 3){
		g.dataRatioScaleSyst = createScaleSystGraph(g.name, lepSel, variable, g.dataRatio);
	    } else if(showSys == 4){
		g.dataRatioScaleSyst = createGenevaIncScaleSystGraph(g.name, lepSel, variable, g.dataRatio);
	    }
	}
    }

    //--- Main Canvas ---
    //double maximum = hGen1->GetMaximum();
    //double minimum = hGen1->GetMinimum();
    double maximum = hStat->GetMaximum();
    double minimum = hStat->GetMinimum();
    TString canvasName = "canvas" + variable;
#ifdef NEW_PLOTS
    //    TCanvas *plots = new TCanvas(canvasName, hStat->GetTitle(), 600, 400*(1 + padHeightRatio * ratioYaxisExtend * ratioFrames.size()));
#else
    TCanvas *plots = new TCanvas(canvasName, hStat->GetTitle(), 600, 800);
#endif
    //-------------------

    //--- First Pad ---
    TPad *plot1 = 0;
    TString canvasTitle  = hStat->GetTitle();
#ifdef NEW_PLOTS
    TCanvas *plots = setAndDrawTPad(canvasName, canvasTitle, plot1, 1, ratioFrames.size());
#else
    setAndDrawTPad(canvasName, plot1, 1, ratioFrames.size());
#endif

    //    plots->cd();

    //--- TLegend ---
    TLegend *legend = new TLegend(0.7, 0.74, 0.99, 0.98);
    customizeMainLegend(canvasName, legend, ratioFrames.size());
    if(grCentralSyst){
	//	legend->AddEntry(grCentralSyst, hStat->GetZaxis()->GetTitle(), "PLEF");
	legend->AddEntry(grCentralSyst, hStat->GetZaxis()->GetTitle(), "plf");
    } else{
	//	legend->AddEntry(grCentralStat, hStat->GetZaxis()->GetTitle(), "PLE");	
	legend->AddEntry(grCentralStat, hStat->GetZaxis()->GetTitle(), "pl");	
    }
    //------------------
    if(grCentralSyst) customizeCentral(grCentralSyst, 0, hStat->GetZaxis()->GetTitle());
    customizeCentral(grCentralStat, false);
    //customizeCentral(grCentralStat, legend, hStat->GetZaxis()->GetTitle());
    if(grCentralSystRatio) customizeCentral(grCentralSystRatio, true);
    if(grCentralStatRatio) customizeCentral(grCentralStatRatio, true);
    if(hSyst){
	hSyst->SetLineColor(kWhite);
	hSyst->SetMarkerColor(kWhite);
	hSyst->SetTitle("");
	hSyst->GetXaxis()->SetLabelSize(0);
	hSyst->GetYaxis()->SetTitle("");
	//-->
	//	hSyst->GetYaxis()->SetLabelSize(0.055);
	hSyst->GetYaxis()->SetLabelFont(ts.defaultFont);
	hSyst->GetYaxis()->SetLabelSize(ts.yLabelSize);

	if (canvasName.Contains("Eta") || canvasName.Contains("AbsRapidity")) {
	    hSyst->GetYaxis()->SetRangeUser(0.001, 1.4*maximum);
	}
	if (canvasName.Contains("DPhi")) {
	    hSyst->GetYaxis()->SetRangeUser(0.2*minimum, 1.5*maximum);
	}
	if (canvasName.Contains("Vis")) {
	    hSyst->GetYaxis()->SetRangeUser(0.2*minimum, 1.3*maximum);
	}
	hSyst->SetStats(0);
	configXaxis(hSyst, 0, variable);
	hSyst->GetXaxis()->SetLabelSize(0.);
	hSyst->GetXaxis()->SetTitleSize(0.);
	//hSyst->DrawCopy("e");
	if(ratioFrames.size() > 0){
	    //if ratio plots are drawn bellow the main frame,
	    //disable the x-axis drawing:
	    hSyst->GetXaxis()->SetTitleSize(0);
	    hSyst->GetXaxis()->SetLabelSize(0);
	}
	
	hSyst->Draw("e");
	if(grCentralSyst){
	    grCentralSyst->SetName("grCentralSyst");
	    grCentralSyst->Draw("2");
	}
    }

    //    for(auto hGen: hGens){
    int igen = -1;
    for(auto& g: mainFrame){
	++igen;
	if(!g.histo) continue;
	customizeGenHist(g.histo, igen + 1, legend, TString::Format("%s", g.histo->GetZaxis()->GetTitle()));
    	g.histo->SetStats(0);
	//	g.histo->DrawCopy("ESAME");
	g.histo->Draw("ESAME");
    }

    grCentralStat->SetName("grCentralStat");
    grCentralStat->Draw("p");

    legend->SetName("mainLegend");
    legend->Draw("same");
    legend->SetTextSize(15);

    //    if (canvasName.Contains("JZB")) fixYscale(1.2, 1.2);
    //if(canvasName.Contains("ZPt")>=0) fixYscale(1.5, 1.5);
    //--- TLatex stuff ---
    //TLatex *latexLabel = new TLatex(); 
    //latexLabel->SetNDC();
    //
    //latexLabel->SetLineWidth(2);
    //latexLabel->SetTextFont(63);
    //cmsLabel->SetTextSize(16);
    TLatex cmsLabel;
    cmsLabel.SetName("cmsLabel");
    cmsLabel.SetNDC();
    cmsLabel.SetTextFont(ts.cmsLabelFont);
    cmsLabel.SetTextSize(ts.cmsLabelSize);
    
    TLatex prelimLabel;
    prelimLabel.SetNDC();
    prelimLabel.SetTextFont(ts.prelimLabelFont);
    prelimLabel.SetTextSize(ts.prelimLabelSize);

    TLatex lumiLabel;
    lumiLabel.SetNDC();
    lumiLabel.SetTextFont(ts.lumiLabelFont);
    lumiLabel.SetTextSize(ts.lumiLabelSize);

    if(TString(hStat->GetZaxis()->GetTitle()).BeginsWith("data", TString::kIgnoreCase)
       || TString(hStat->GetZaxis()->GetTitle()).BeginsWith("meas", TString::kIgnoreCase)){
	cmsLabel.SetTextAlign(12);
	cmsLabel.DrawLatex(0.17,0.833,"CMS");
	if(isPrel) prelimLabel.DrawLatex(0.17,0.833-0.09,"Preliminary");
	if(lumi > 0) lumiLabel.DrawLatex(0.13,0.905, TString::Format("%.3g fb^{-1} (13 TeV)", lumi/1000.));
    } else{
	cmsLabel.DrawLatex(0.17,0.82,"MC study");	
	if(isPrel) prelimLabel.DrawLatex(0.17,0.82-0.09,"Preliminary");
    }

    TLatex descLabel;
    descLabel.SetNDC();
    descLabel.SetTextFont(ts.defaultFont);
    descLabel.SetTextSize(ts.descTextSize);

    double xlabel;
    double ylabel;

    if (canvasName.Contains("JZB_ptLow")){
	xlabel = 0.25;
	ylabel = 0.25;
    }  else if (canvasName.Contains("JZB_ptHigh")){
	xlabel = 0.3;
	ylabel = 0.23;
    } else if (canvasName.Contains("JZB")){
	xlabel = 0.30;
	ylabel = 0.23;
    } else if (canvasName.Contains("ZPt_Zinc1jet")){
	xlabel = 0.18;
	ylabel = 0.60;
    }  else if (canvasName.Contains("ZPt_Zinc0jet")){
	xlabel = 0.60;
	ylabel = 0.4;
    } else{
	xlabel = 0.18;
	ylabel = 0.23;
    }
    if (!canvasName.EndsWith("Zinc0jet")) descLabel.DrawLatex(xlabel, ylabel - 0.05, "Anti-#it{k}_{T} (R = 0.4) jets");

    if (canvasName.Contains("FirstJetPt50")){
        descLabel.DrawLatex(xlabel,ylabel-0.11,"p_{T}^{jet} > 50 GeV, |y^{jet}| < 2.4 ");
    }else if (canvasName.Contains("FirstJetPt80")){
        descLabel.DrawLatex(xlabel,ylabel-0.11,"p_{T}^{jet} > 80 GeV, |y^{jet}| < 2.4 ");
    } else if (canvasName.Contains("ZPt150")){
        descLabel.DrawLatex(xlabel,ylabel-0.11,"p_{T}^{Z} > 150 GeV, p_{T}^{jet} > 30 GeV, |y^{jet}| < 2.4 ");
    } else if (canvasName.Contains("ZPt300")){
        descLabel.DrawLatex(xlabel,ylabel-0.11,"p_{T}^{Z} > 300 GeV, p_{T}^{jet} > 30 GeV, |y^{jet}| < 2.4 ");
    } else if (canvasName.Contains("DifJetRapidityl2")){
        descLabel.DrawLatex(xlabel,ylabel-0.11,"p_{T}^{jet} > 30 GeV, |y^{jet}| < 2.4, |y_{jet1}-y_{jet2}| > 2 ");
    } else if (canvasName.Contains("DifJetRapiditys2")){
        descLabel.DrawLatex(xlabel,ylabel-0.11,"p_{T}^{jet} > 30 GeV, |y^{jet}| < 2.4, |y_{jet1}-y_{jet2}| < 2 ");
    } else if (canvasName.Contains("ZPt150_HT300")){
        descLabel.DrawLatex(xlabel,ylabel-0.11,"p_{T}^{Z} > 150 GeV, p_{T}^{jet} > 30 GeV, |y^{jet}| < 2.4, H_{T}^{jet} > 300 GeV ");
    //} else if (canvasName.Contains("Vis")){
    //	descLabel.DrawLatex(xlabel, 0.7-0.06,"p_{T}^{#mu} > 20 GeV, p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.4 ");
    } else if (!canvasName.EndsWith("Zinc0jet")){
        descLabel.DrawLatex(xlabel, ylabel - 0.11,"p_{T}^{jet} > 30 GeV, |y^{jet}| < 2.4 ");
    }

    if (lepSel == "") { 
	TLatex ll;
	ll.SetTextFont(42);
	ll.SetTextSize(0.05);
	ll.SetNDC();
	if(canvasName.Contains("inc1"))  ll.DrawLatex(xlabel,ylabel-0.17,"\\text{Z}/\\gamma^{*} \\rightarrow \\ell^{+} \\ell^{-}, \\text{N}_{\\text{jets}} \\geq 1");
        else if(canvasName.Contains("inc2"))  ll.DrawLatex(xlabel,ylabel-0.17,"\\text{Z}/\\gamma^{*} \\rightarrow \\ell^{+} \\ell^{-}, \\text{N}_{\\text{jets}} \\geq 2");
        else if(canvasName.Contains("inc3"))  ll.DrawLatex(xlabel,ylabel-0.17,"\\text{Z}/\\gamma^{*} \\rightarrow \\ell^{+} \\ell^{-}, \\text{N}_{\\text{jets}} \\geq 3");
        else if(canvasName.Contains("JZB_ptHigh")) ll.DrawLatex(xlabel,ylabel-0.17,"\\text{Z}/\\gamma^{*} \\rightarrow \\ell^{+} \\ell^{-}, \\text{N}_{\\text{jets}} \\geq 1, \\text{p}_{\\text{T}}(\\text{Z}) > 50\\,\\text{GeV}");
        else if(canvasName.Contains("JZB_ptLow")) ll.DrawLatex(xlabel,ylabel-0.17,"\\text{Z}/\\gamma^{*} \\rightarrow \\ell^{+} \\ell^{-}, \\text{N}_{\\text{jets}} \\geq 1, \\text{p}_{\\text{T}}(\\text{Z}) \\leq 50\\,\\text{GeV}");
        else if(canvasName.Contains("JZB")) ll.DrawLatex(xlabel,ylabel-0.17,"\\text{Z}/\\gamma^{*} \\rightarrow \\ell^{+} \\ell^{-}, \\text{N}_{\\text{jets}} \\geq 1");
	else ll.DrawLatex(xlabel,ylabel-0.17,"\\text{Z}/\\gamma^{*} \\rightarrow \\ell^{+} \\ell^{-}");
    }

    else if (lepSel == "DMu"){
         if(canvasName.Contains("inc1")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow #mu#mu channel, N_{jets} #geq 1");
         else if(canvasName.Contains("inc2")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow #mu#mu channel, N_{jets} #geq 2");
         else if(canvasName.Contains("inc3")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow #mu#mu channel, N_{jets} #geq 3");
         else if(canvasName.Contains("JZB_ptHigh")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow #mu#mu channel, N_{jets} #geq 1, p_{T}(Z) > 50 GeV");
         else if(canvasName.Contains("JZB_ptLow")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow #mu#mu channel, N_{jets} #geq 1, p_{T}(Z) #leq 50 GeV");
         else if(canvasName.Contains("JZB")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow #mu#mu channel, N_{jets} #geq 1");
	 else descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow #mu#mu channel");

     }

    else if (lepSel == "DE") {
         if(canvasName.Contains("inc1")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow ee channel, N_{jets} #geq 1");
         else if(canvasName.Contains("inc2")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow ee channel, N_{jets} #geq 2");
         else if(canvasName.Contains("inc3")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow ee channel, N_{jets} #geq 3");
         else if(canvasName.Contains("JZB_ptHigh")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow ee channel, N_{jets} #geq 1, p_{T}(Z) > 50 GeV");
         else if(canvasName.Contains("JZB_ptLow")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow ee channel, N_{jets} #geq 1, p_{T}(Z) < 50 GeV");
         else if(canvasName.Contains("JZB")) descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow ee channel, N_{jets} #geq 1");
	 else descLabel.DrawLatex(xlabel,ylabel-0.17,"Z/#gamma*#rightarrow ee channel");
    }


    //    descLabel.SetName("descLabel");
    //    descLabel.Draw("same");

    TLatex *ytitle = new TLatex();
    ytitle->SetName("ytitle");
    //    if (gens.size() == 2) ytitle->SetTextSize(0.05);
    //    if (gens.size() == 3) ytitle->SetTextSize(0.06);
    ytitle->SetTextFont(ts.defaultFont);
    ytitle->SetTextSize(ts.mainYTitleSize);
    ytitle->SetLineWidth(2);
    ytitle->SetTextColor(kBlack);
    ytitle->SetNDC();
    ytitle->SetTextAlign(33);
    ytitle->SetTextAngle(90);
    //TODO: check if we can use hStat anytime
    //    std::string strYtitle = hGen1 ?  getYaxisTitle(doNormalized, hGen1) :  getYaxisTitle(doNormalized, hStat);
    std::string strYtitle = getYaxisTitle(doNormalized, hStat);
    if (strYtitle.find("eta") != std::string::npos) {
        size_t first = strYtitle.find("#eta");
        std::string tmp1 = strYtitle.substr(0, first);
        std::string tmp2 = strYtitle.substr(first);
        strYtitle = tmp1 + "|" + tmp2;
        size_t second = strYtitle.find(")");
        tmp1 = strYtitle.substr(0, second+1);
        tmp2 = strYtitle.substr(second+1);
        strYtitle = tmp1 + "|" + tmp2;
    }
    ytitle->DrawLatex(0.008,0.91,strYtitle.c_str());

    fixYscale(1.2, 1.1);

    if(variable.Contains("JZB_ptHigh")) hSyst->GetYaxis()->SetRangeUser(1.1e-4, 1.e2);

    //--- End Of first Pad ---

    //--- Ratio Pads ---
    //    igen = -1;
    int ipad = 1;
    igen = -1;
    double minRatioY = cfg.getD("minRatioYUnf", 0.2);
    double maxRatioY = cfg.getD("maxRatioYUnf", 1.8);
    extendAxis(minRatioY, maxRatioY);
    
    for(auto& f: ratioFrames){
	++ipad;
	    //	    ++igen;
	    //	    int ipad = igen + 2;
	plots->cd();
	TString padName = TString::Format("plot%d", ipad);
	TPad *pad = new TPad(padName, padName, 0., 0., 0., 0.);
	setAndDrawTPad(canvasName, canvasTitle, pad, ipad, ratioFrames.size());

	//--- TLegend ---
	TLegend *legend = new TLegend(0.16, 0.05, 0.42, 0.20);
        int pos = 0;
	//std::cout << "XXX " << variable << "\n";
	if(variable == "JetsHT_Zinc1jet") pos = 1;
	if(variable == "JetsHT_Zinc2jet") pos = 1;
	if(variable == "JetsHT_Zinc3jet") pos = 1;
	if(variable.BeginsWith("ZPt")) pos = 1;


	if(variable.BeginsWith("VisPt")){
	    if(variable.BeginsWith("VisPt_Zinc2jet") && TString(f[0].name).BeginsWith("DYJets_GE")){
		pos = 2;
	    } else{
		pos = 1;
	    }
	}


	if(variable.BeginsWith("JZB") && f[0].name != "DYJets_UNFOLDING")  pos = 1;

	customizeRatioLegend(canvasName, legend, ipad - 1, ratioFrames.size(), pos);
	//TString generator = hGen->GetZaxis()->GetTitle();
	//generator = generator(0, generator.Index(" "));
	TString ref_shortname = hStat->GetZaxis()->GetTitle();
	if(ref_shortname.Contains(" ")){
	    ref_shortname = ref_shortname(0, ref_shortname.Index(" "));
	}
	if(ref_shortname.Length()==0) ref_shortname = "Measurement";
//	customizeRatioGraph(hSyst, dataRatioToCentral[igen], dataRatioScaleSyst[igen], dataRatioPDFSyst[igen], igen + 1,
//			    //TString("#frac{") + generator + "}{" + ref_shortname + "}", ratioFrames.size(), legend);
//			    TString::Format("#frac{Prediction}{%s}", ref_shortname.Data()), ratioFrames.size(), legend);
		
	//Histogram to draw axis of ratio plots
	TH1* hAxis = (TH1*) hStat->Clone(TString("hAxis%d", ipad));
	hAxis->Reset();
	configXaxis(hAxis, 0, variable);
//	hAxis->GetXaxis()->SetLabelSize(0);
//	hAxis->GetXaxis()->SetTitleSize(0);
	customizeRatioGraph(hAxis, 0, 0, 0, ipad - 1,
			    TString::Format("#frac{Prediction}{%s}",
					    ref_shortname.Data()),
			    ratioFrames.size(), 0);
	hAxis->DrawClone("");
	// hSyst->DrawCopy("e");
	bool mainGen = true;
	std::vector<TGraphAsymmErrors*> stairs_to_draw;
	std::vector<int> stairs_to_draw_maxPoints;
	for(auto& g: f){
	    ++igen;
	    customizeRatioGraph(mainGen ? hAxis : 0, g.dataRatio, g.dataRatioScaleSyst,
				g.dataRatioPDFSyst, ipad - 1,
				TString::Format("#frac{Prediction}{%s}",
						ref_shortname.Data()),
				ratioFrames.size(), mainGen ? legend : 0);
	    int showSys = cfg.getI(TString::Format("%s_unc", g.name.c_str()));
	    if(g.dataRatioScaleSyst) g.dataRatioScaleSyst->Draw("2");
	    if(g.dataRatioPDFSyst)   g.dataRatioPDFSyst->Draw("2");

	    //FIXME: should go in configuration file
	    if(!mainGen && g.dataRatio && TString(g.name).Contains("as1135")){
		TLegend* l = new TLegend(0.65, 0.05, 0.8, 0.15);
		l->SetFillColor(legend->GetFillColor());
		l->SetBorderSize(legend->GetBorderSize());
		l->SetMargin(legend->GetMargin());
		l->SetTextFont(legend->GetTextFont());
		l->SetTextSize(legend->GetTextSize()*1.2);
		l->SetTextAlign(legend->GetTextAlign());
		l->SetY1(legend->GetY1());
		l->SetY2(legend->GetY2());

		l->AddEntry(g.dataRatio, "#alpha_{s} = 0.1135", "l");
		l->Draw();
	    }
	    
	    if(g.dataRatio){
		g.dataRatio->SetName("grGen1ToCentral");
		if(showSys >= 0){
		    g.dataRatio->Draw("2p");
		} else{
		    g.dataRatio->SetLineStyle(kDashed);
		    g.dataRatio->SetLineWidth(4.);
		    int maxPoints = 99999;
		    if(variable.BeginsWith("ZNGoodJet")){
			maxPoints = 3;
		    }
		    stairs_to_draw.push_back(g.dataRatio);
		    stairs_to_draw_maxPoints.push_back(maxPoints);
		}
		//FIXME: replace by one call to draw with the combined X2p option?
		//g.dataRatio->Draw("Xp");
	    }
	    //	    legend->Draw("same");
	    legend->Draw();
	
	//if (canvasName.Contains("JetPt_Zinc")) {
	//    dataRatioToCentral->GetXaxis()->SetRangeUser(30, x + ex);
	//}

	//	if(variable.Contains("ZPt_") && igen == 0){
	//    draw_axis_labels(hcopy->GetXaxis());
	//	    plots->Update();
	//	} 
	
	    mainGen = false;
	}//next gen of the ratio frame
	
	//draw data/data:
	if(grCentralSystRatio){
	    grCentralSystRatio->SetName("grCentralSystRatio");
	    //TODO: replace with one draw with "2p" option?
	    grCentralSystRatio->Draw("2");
	    grCentralStatRatio->Draw("p");
	}

	int i = 0;
	for(auto g: stairs_to_draw){
	    // graph_draw_stairs(g, stairs_to_draw_maxPoints[i]);
	    graph_draw_stairs(g, minRatioY, maxRatioY, false);
	    ++i;
	}

	
	pad->RedrawAxis();
    }
    //--- End of Ratio Pads ---
    
    plots->Update();

    ratioYaxisExtend = savedRatioYaxisExtend;
    
    return plots;
}

