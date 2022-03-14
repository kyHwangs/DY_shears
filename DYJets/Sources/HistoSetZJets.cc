#define PI 3.14159265359
#include "HistoSetZJets.h"
#include "ConfigVJets.h"
#include <RooUnfoldResponse.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <iostream>
#include <sstream>

extern ConfigVJets cfg; // defined in runZJets_newformat.cc

using namespace std;

HistoSetZJets::~HistoSetZJets() {}

bool HistoSetZJets::filterHist(const char *name) const
{
    if (varList.size() == 0) return true;
    std::string n(name);
    // DJALOG
    // printf("HistoSetZJets::filterHist name=%s\n",name);
    for (std::set<std::string>::const_iterator it = varList.begin(); it != varList.end(); ++it) {
        if ((*it) == n) {
            // printf("Got it\n");
            return true;
        }
        if ((*it) + "_Odd" == n) return true;
        if ((*it) + "_Even" == n) return true;
        if ((*it) + "_2" == n) return true;
        if (std::string("gen") + (*it) == n) return true;
        if (std::string("hresponse") + (*it) == n) return true;
    }
    return false;
}

void HistoSetZJets::readHistList()
{
    std::string fname = cfg.getS("histList");
    if (fname.empty()) return;
    std::ifstream f(fname);
    if (!f.good()) {
        std::cerr << "Failed to read file " << fname << " defined by parameter varList\n";
        abort();
    }
    while (!f.eof()) {
        std::string v;
        f >> v;
        varList.insert(v);
    }
    // Histogram to always include:
    varList.insert("JobInfo");
    varList.insert("lumi");
    varList.insert("input");
}

vector<double> HistoSetZJets::makeVector(int num, ...)
{
    va_list list;
    va_start(list, num);
    vector<double> vec;
    for (int i(0); i < num; i++) {
        double next = va_arg(list, double);
        vec.push_back(next);
    }
    va_end(list);
    return vec;
}

void HistoSetZJets::insertVector(vector<double> &veca, int num, ...)
{
    va_list list;
    va_start(list, num);
    vector<double> vecb;
    for (int i(0); i < num; i++) {
        double next = va_arg(list, double);
        vecb.push_back(next);
    }
    va_end(list);
    veca.insert(veca.end(), vecb.begin(), vecb.end());
}

vector<double> HistoSetZJets::buildVecFineBin(int nStdBin, double arrStdBin[], int factChop)
{
    vector<double> vecTemp;
    for (int i = 0; i < nStdBin; i++) {
        double binWidth = (arrStdBin[i + 1] - arrStdBin[i]) / 5;
        for (int j = 0; j < factChop; j++) {
            double element(0.);
            element = arrStdBin[i] + (j * binWidth);
            vecTemp.push_back(element);
        }
    }
    vecTemp.push_back(arrStdBin[nStdBin]);
    return vecTemp;
}

GenH1D *HistoSetZJets::newTH1D(string name, string title, string xTitle, int nBins, double *xBins)
{
    if (!filterHist(name.c_str())) return 0;
    GenH1D *hist = new GenH1D(name.c_str(), title.c_str(), nBins, xBins);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
    listOfHistograms.push_back(hist);
    return hist;
}

GenH1D *HistoSetZJets::newTH1D(string name, string title, string xTitle, vector<double> &xBinsVect)
{
    if (!filterHist(name.c_str())) return 0;
    int nBins = xBinsVect.size() - 1;
    double *xBins = new double[xBinsVect.size()];
    std::copy(xBinsVect.begin(), xBinsVect.end(), xBins);
    GenH1D *hist = new GenH1D(name.c_str(), title.c_str(), nBins, xBins);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
    delete[] xBins;
    listOfHistograms.push_back(hist);
    return hist;
}

GenH1D *
HistoSetZJets::newTH1D(string name, string title, string xTitle, int nBins, double xLow, double xUp)
{
    if (!filterHist(name.c_str())) return 0;
    GenH1D *hist = new GenH1D(name.c_str(), title.c_str(), nBins, xLow, xUp);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
    hist->SetOption("HIST");
    listOfHistograms.push_back(hist);
    return hist;
}

TH2D *HistoSetZJets::newTH2D(
    string name, string title, int nBinsX, double *xBins, int nBinsY, double *yBinsY)
{
    if (!filterHist(name.c_str())) return 0;
    TH2D *hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xBins, nBinsY, yBinsY);
    hist->GetZaxis()->SetTitle("# Events");
    listOfHistograms.push_back(hist);
    return hist;
}

TH2D *HistoSetZJets::newTH2D(
    string name, string title, int nBinsX, double *xBins, int nBinsY, double yLow, double yUp)
{
    if (!filterHist(name.c_str())) return 0;
    TH2D *hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xBins, nBinsY, yLow, yUp);
    hist->GetZaxis()->SetTitle("# Events");
    listOfHistograms.push_back(hist);
    return hist;
}

TH2D *HistoSetZJets::newTH2D(
    string name, string title, int nBinsX, double xLow, double xUp, int nBinsY, double *yBins)
{
    if (!filterHist(name.c_str())) return 0;
    TH2D *hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xLow, xUp, nBinsY, yBins);
    hist->GetZaxis()->SetTitle("# Events");
    listOfHistograms.push_back(hist);
    return hist;
}

TH2D *HistoSetZJets::newTH2D(string name,
                             string title,
                             int nBinsX,
                             double xLow,
                             double xUp,
                             int nBinsY,
                             double yLow,
                             double yUp)
{
    if (!filterHist(name.c_str())) return 0;
    TH2D *hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xLow, xUp, nBinsY, yLow, yUp);
    hist->GetZaxis()->SetTitle("# Events");
    hist->SetOption("HIST");
    listOfHistograms.push_back(hist);
    return hist;
}

TH2D *HistoSetZJets::newTH2D(string name,
                             string title,
                             vector<double> &xBinsVect,
                             vector<double> &yBinsVect)
{
    if (!filterHist(name.c_str())) return 0;
    int nBins_x = xBinsVect.size() - 1;
    int nBins_y = yBinsVect.size() - 1;
    double *xBins = new double[xBinsVect.size()];
    double *yBins = new double[yBinsVect.size()];
    std::copy(xBinsVect.begin(), xBinsVect.end(), xBins);
    std::copy(yBinsVect.begin(), yBinsVect.end(), yBins);
    TH2D *hist = new TH2D(name.c_str(), title.c_str(), nBins_x, xBins, nBins_y, yBins);
    hist->GetZaxis()->SetTitle("# Events");
    delete[] xBins;
    delete[] yBins;
    listOfHistograms.push_back(hist);
    return hist;
}

HistoSetZJets::HistoSetZJets(TString leptonFlavor)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    readHistList();

    string ZpT = "p_{T}(Z) [GeV]", Zrap = "y(Z)", Zeta = "#eta(Z)";
    string phistar = "#phi_{#eta}^{*}";
    string ZpTVis = "p_{T} balance [GeV]";
    string HRecoi = "Hadronic Recoil [GeV]";
    string JZb = "JZB [GeV]";
    string JZb_low = "JZB [GeV]";
    string JZb_high = "JZB [GeV]";
    string HT = "H_{T}(jets) [GeV]", Mjj = "M_{j_{1}j_{2}} [GeV]",
           jSpt = "#Delta_{pT}^{rel}(j_{1}j_{2})", jdPhi = "#Delta#phi(j_{1}j_{2})",
           jdEta = "#Delta#eta(j_{1}j_{2})";
    string Mll = "M_{#mu#mu} [GeV]", leta = "#eta(#mu)", lphi = "#phi(#mu)",
           lpT = "p_{T}(#mu) [GeV]", ldPhi = "#Delta#phi(#mu_{1}#mu_{2})",
           ldEta = "#Delta#eta(#mu_{1}#mu_{2})", ldR = "#DeltaR(#mu_{1}#mu_{2})";
    string lSpt = "#Delta_{pT}^{rel}(#mu_{1}#mu_{2})";
    string Spt = "#Delta_{pT}^{rel}(j_{1}j_{2}#mu_{1}#mu_{2})";
    string Sphi = "Sphi(j_{1}j_{2}#mu_{1}#mu_{2})";
    string lJetdEta = "#Delta#eta(#mu_{1}#mu_{2},j_{1})";

    bool doWJets = false;
    if (leptonFlavor == "Electrons" || leptonFlavor == "DE" || leptonFlavor == "DE_") {
        Mll = "M_{ee} [GeV]";
        leta = "#eta(e)";
        lpT = "p_{T}(e) [GeV]";
        ldPhi = "#Delta#phi(e_{1}e_{2})";
        ldEta = "#Delta#eta(e_{1}e_{2})";
        ldR = "#DeltaR(e_{1}e_{2})";
        lSpt = "#Delta_{pT}^{rel}(e_{1}e_{2})";
        Spt = "#Delta_{pT}^{rel}(j_{1}j_{2}e_{1}e_{2})";
        Sphi = "Sphi(j_{1}j_{2}e_{1}e_{2})";
        lJetdEta = "#Delta#eta(e_{1}e_{2},j_{1})";
    } else if (leptonFlavor == "Electron" || leptonFlavor == "SE" || leptonFlavor == "SE_") {
        doWJets = true;
        Mll = "M_{e#nu} [GeV]";
        leta = "#eta(e)";
        lpT = "p_{T}(e) [GeV]";
        ldPhi = "#Delta#phi(e_{1}#nu_{2})";
        ldEta = "#Delta#eta(e_{1}#nu_{2})";
        ldR = "#DeltaR(e_{1}#nu_{2})";
        lSpt = "#Delta_{pT}^{rel}(e_{1}#nu_{2})";
        Spt = "#Delta_{pT}^{rel}(j_{1}j_{2}e_{1}#nu_{2})";
        Sphi = "Sphi(j_{1}j_{2}e_{1}#nu_{2})";
        lJetdEta = "#Delta#eta(e,j_{1})";

    } else if (leptonFlavor == "Muon" || leptonFlavor == "SMu" || leptonFlavor == "SMu_") {
        doWJets = true;
        Mll = "M_{#mu#nu} [GeV]";
        ldPhi = "#Delta#phi(#mu_{1}#nu_{2})";
        ldEta = "#Delta#eta(#mu_{1}#nu_{2})";
        ldR = "#DeltaR(e_{1}#nu_{2})";
        lSpt = "#Delta_{pT}^{rel}(#mu_{1}#nu_{2})";
        Spt = "#Delta_{pT}^{rel}(j_{1}j_{2}#mu_{1}#nu_{2})";
        Sphi = "Sphi(j_{1}j_{2}#mu_{1}#nu_{2})";
        lJetdEta = "#Delta#eta(#mu,j_{1})";
    }

    //----------- 2016 analysis ----------------

    int nPhistar_Zinc0jet(34);
    double phistar_Zinc0jet[35] = {
        0.001, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.051,
        0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219,
        0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277}; // first bin


    int nPhistar_Zinc1jet(17);
    double phistar_Zinc1jet[18] = {
        0.001, 0.008,  0.016,  0.024,  0.034,  0.045, 
        0.057,  0.072, 0.091,  0.114,  0.145,  0.189, 
        0.258,  0.391,  0.695, 1.153,  1.947,  3.277}; // first bin

    int nPhistar_Zinc1jetMbins(8);
    double phistar_Zinc1jetMbins[9] = {
        0.001,  0.016,   0.034,  
        0.057,   0.091,   0.145,  
        0.258,   0.695,   1.947}; // first bin

    // double mass_Zinc0jet[] = {116,130,150,175,200,230,260,300,380,500,700,1000,1500,2000}
    int nMass_Zinc0jet(81);
    double mass_Zinc0jet[82] = {20,   35,   50,   55,   60,   65,   70,   75,   80,   85,   90,
                                95,   100,  105,  110,  115,  120,  125,  130,  135,  140,  145,
                                150,  160,  170,  180,  190,  200,  220,  240,  260,  280,  300,
                                320,  340,  360,  380,  400,  420,  440,  460,  480,  500,  520,
                                540,  560,  580,  600,  630,  660,  690,  720,  750,  780,  810,
                                840,  870,  900,  940,  980,  1020, 1060, 1100, 1140, 1180, 1220,
                                1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1660, 1720, 1780,
                                1840, 1900, 1980, 2060, 2140}; //  , 2220, 2300, 2380, 2460, 2500,
                                                               //  2600, 2700, 2800, 2900, 3000,
                                                               //  3100, 3200, 3300, 3400, 3500,
                                                               //  3600, 3700, 3800, 3900, 4000

    // int nZPt_Zinc0jet(25);
    // double zPt_Zinc0jet[26] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
    // 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300};

    // int nZPt_Zinc0jet(33);
    // double zPt_Zinc0jet[34] = {0., 1.25, 2.5, 3.75, 5, 6.25, 7.5, 8.75, 10, 11.25, 12.5, 15,
    // 17.5, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 220, 250,
    // 400, 1000};
    // new 2016 binning
    int nZPt_Zinc0jet(37);
    double zPt_Zinc0jet[38] = {0.1,  1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,
                               10.,  11.,  12.,  13.,  14.,  16.,  18.,  20.,  22.,  25.,
                               28.,  32.,  37.,  43.,  52.,  65.,  85.,  120., 160., 190.,
                               220., 250., 300., 350., 400., 450., 500., 1000.};

    int nZPt_Zinc0jetM115_135(7);
    double zPt_Zinc0jetM115_135[8] = {0., 15., 30., 45., 85., 125., 200., 350.};

    int nZPt_Zinc0jetMbins(16);
    double zPt_Zinc0jetMbins[17] = {
        0.1, 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 37., 52., 85., 160., 240., 1000.};

    vector<double> zPt_2_Zinc0jetMbins;
    zPt_2_Zinc0jetMbins = buildVecFineBin(nZPt_Zinc0jetMbins, zPt_Zinc0jetMbins, 5);


    int nZPt_Zinc1jetMbins(12);
    double zPt_Zinc1jetMbins[13] = {
        0.1, 4., 8., 12., 18., 22., 28., 37., 52., 85., 160., 240., 1000.};

    int nZPt_Zinc1jetHighMbins(10);
    double zPt_Zinc1jetHighMbins[11] = {
        0.1,  8.,  18., 22., 28., 37., 52., 85., 160., 240., 1000.};

    /*

        int nZPt_Zinc0jetMbins170(16);
       // double zPt_Zinc0jetMbins170[16] = {0.5, 2., 4., 6.,  8., 10., 12.,14., 18., 22.,  28.,
       37.,  52.,  85.,  200., 1000.};
        double zPt_Zinc0jetMbins170[17] = {0.1, 2., 4., 6.,  8., 10., 12.,14., 18., 22.,  28., 37.,
       52.,  85.,  160.,  240., 1000.};
    */

    int nZPt_Zinc0jetMbins250(12);
    double zPt_Zinc0jetMbins250[13] = {
        0.1, 2.5, 5., 7.5, 10., 13., 18., 22., 28., 40., 90., 200., 1000.};

    vector<double> zPt_2_Zinc0jetMbins250;
    zPt_2_Zinc0jetMbins250 = buildVecFineBin(nZPt_Zinc0jetMbins250, zPt_Zinc0jetMbins250, 5);
    // double zPt_Zinc0jetMbins250[11] = {0.1, 4., 8., 13., 18., 22.,  28., 40., 90., 200., 1000.};

    int nZPt_Zinc0jetMbins320(12);
    // double zPt_Zinc0jetMbins250[11] = {0.5, 4., 8., 12., 18., 22.,  28., 40., 90., 200., 1000.};
    //  double zPt_Zinc0jetMbins320[11] = {0.1, 4., 8., 13., 18., 22.,  28., 40., 90., 200., 1000.};
    double zPt_Zinc0jetMbins320[13] = {
        0.1, 2.5, 5., 7.5, 10., 13., 18., 22., 28., 40., 90., 200., 1000.};

    // int nZPt_Zinc0jetMbins1(15);
    // double zPt_Zinc0jetMbins1[16] = {0.5, 2.5, 5., 7.5, 10, 12.5, 15., 18., 22.,  28., 37.,  52.,
    // 85.,  160.,  240., 1000.};

    int nZPt_Zinc1jet(33); // 0,2.5,5,7.5,10,12.5,17.5,25,35,45,60,80,100,130,170,220,400,600
    double zPt_Zinc1jet[34] = {0.1,  1.25, 2.5, 3.75, 5,   6.25, 7.5, 8.75, 10,  11.25, 12.5, 15,
                               17.5, 20,   25,  30,   35,  40,   45,  50,   60,  70,    80,   90,
                               100,  110,  130, 150,  170, 190,  220, 250,  400, 1000};
    // double zPt_Zinc1jet[25] = {0., 2.5, 5, 7.5, 10, 12.5, 17.5, 25, 35, 45, 50, 60, 70, 80, 90,
    // 100, 110, 130, 150, 170, 190, 220, 250, 400, 1000};
    vector<double> zPt_2_Zinc1jet = buildVecFineBin(nZPt_Zinc1jet, zPt_Zinc1jet, 5);

    int nZPt_Zinc2jet(24);
    double zPt_Zinc2jet[25] = {0., 2.5, 5,   7.5, 10,  12.5, 17.5, 25,  35,  45,  50,  60,  70,
                               80, 90,  100, 110, 130, 150,  170,  190, 220, 250, 400, 1000};
    vector<double> zPt_2_Zinc2jet = buildVecFineBin(nZPt_Zinc2jet, zPt_Zinc2jet, 5);

    //----------------------
    int nZPt_Zinc2jetQun(11);
    double zPt_Zinc2jetQun[12] = {0, 10, 20, 35, 50, 65, 80, 100, 125, 150, 175, 200}; //, 260,
                                                                                       //400};//
                                                                                       //240, 260,
                                                                                       //280, 300,
                                                                                       //400, 800};
    vector<double> zPt_2_Zinc2jetQun = buildVecFineBin(nZPt_Zinc2jetQun, zPt_Zinc2jetQun, 5);

    int nZPt_Zinc2JetQun(10);
    double zPt_Zinc2JetQun[11] = {0, 15, 30, 45, 60, 80, 100, 125, 150, 175, 200};
    vector<double> zPt_2_Zinc2JetQun = buildVecFineBin(nZPt_Zinc2JetQun, zPt_Zinc2JetQun, 5);

    int nZPt_Zinc3jetQun(8);
    double zPt_Zinc3jetQun[9] = {0, 20, 40, 65, 90, 120, 150, 175, 200};
    vector<double> zPt_2_Zinc3jetQun = buildVecFineBin(nZPt_Zinc3jetQun, zPt_Zinc3jetQun, 5);

    int nZPt_Zinc2jetQunJZB(16);
    double zPt_Zinc2jetQunJZB[17] = {
        -200, -165, -140, -105, -80, -60, -40, -20, 0, 20, 40, 60, 85, 110, 140, 165, 200};
    vector<double> zPt_2_Zinc2jetQunJZB =
        buildVecFineBin(nZPt_Zinc2jetQunJZB, zPt_Zinc2jetQunJZB, 5);

    int nZPt_Zinc2jetQunJZBptHigh(11);
    double zPt_Zinc2jetQunJZBptHigh[12] = {-165, -125, -95, -70, -45, -20, 0, 25, 55, 85, 120, 150};
    vector<double> zPt_2_Zinc2jetQunJZBptHigh =
        buildVecFineBin(nZPt_Zinc2jetQunJZBptHigh, zPt_Zinc2jetQunJZBptHigh, 5);

    // int nZPt_Zinc2jetQunJZB(30);
    // double zPt_Zinc2jetQunJZB[31] = {-200, -180, -160, -140, -120, -100, -90, -80, -70, -60, -50,
    // -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200};//,
    // 260, 400};// 240, 260, 280, 300, 400, 800};

    int nZPt_Zinc1jetQunJZB(10);
    double zPt_Zinc1jetQunJZB[11] = {-50, -30, -15, 0, 15, 30, 50, 75, 105, 150, 200};
    vector<double> zPt_2_Zinc1jetQunJZB =
        buildVecFineBin(nZPt_Zinc1jetQunJZB, zPt_Zinc1jetQunJZB, 5);
    // int nZPt_Zinc1jetQunJZB(20);
    // double zPt_Zinc1jetQunJZB[21] = {-200, -180, -160, -140, -120, -100, -90, -80, -70, -60, -50,
    // -40, -30, -20, -10, 0, 10, 20, 30, 40, 50};//, 260, 400};// 240, 260, 280, 300, 400, 800};
    //    int nZPt_Zinc4jetQun(8);
    //    double zPt_Zinc4jetQun[9] = {0, 15, 30, 45, 60, 80, 120, 160, 200};//, 260, 400};// 240,
    //    260, 280, 300, 400, 800};

    // int nJetPt_Zinc1jet(22);
    // double jetPt_Zinc1jet[23] = {20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220,
    // 258, 300, 350, 400, 450, 500, 590, 700, 1000};

    // int nJetPt_Zinc1jet(9);
    // double jetPt_Zinc1jet[10] = {20, 30, 41, 59, 83, 118, 168, 220, 300, 400};
    double jetPt_Zinc1jet[11] = {20, 24, 30, 41, 59, 83, 118, 168, 220, 300, 400};
    int nJetPt_Zinc1jet(sizeof(jetPt_Zinc1jet) / sizeof(jetPt_Zinc1jet[0]) - 1);
    vector<double> jetPt_2_Zinc1jet;
    jetPt_2_Zinc1jet = buildVecFineBin(nJetPt_Zinc1jet, jetPt_Zinc1jet, 5);

    //  int nJetPt_Zinc2jet(21);
    //  double jetPt_Zinc2jet[22] = {20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220,
    //  258, 300, 350, 400, 450, 500, 590, 800};
    // int nJetPt_Zinc2jet(7);
    // double jetPt_Zinc2jet[8] = {20, 30, 41, 59, 83, 118, 168, 250};
    double jetPt_Zinc2jet[9] = {20, 24, 30, 41, 59, 83, 118, 168, 250};
    int nJetPt_Zinc2jet(sizeof(jetPt_Zinc2jet) / sizeof(jetPt_Zinc2jet[0]) - 1);
    vector<double> jetPt_2_Zinc2jet;
    jetPt_2_Zinc2jet = buildVecFineBin(nJetPt_Zinc2jet, jetPt_Zinc2jet, 5);

    // int nJetPt_Zinc3jet(11);
    // double jetPt_Zinc3jet[12] = {20, 24, 30, 39, 49, 62, 78, 105, 142, 185, 235, 300};
    //  double jetPt_Zinc2jet[22] = {20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220,
    //  258, 300, 350, 400, 450, 500, 590, 800};
    // int nJetPt_Zinc3jet(7);
    // double jetPt_Zinc3jet[8] = {20, 30, 41, 59, 83, 118, 168, 250};
    double jetPt_Zinc3jet[9] = {20, 24, 30, 41, 59, 83, 118, 168, 250};
    int nJetPt_Zinc3jet(sizeof(jetPt_Zinc3jet) / sizeof(jetPt_Zinc3jet[0]) - 1);
    vector<double> jetPt_2_Zinc3jet;
    jetPt_2_Zinc3jet = buildVecFineBin(nJetPt_Zinc3jet, jetPt_Zinc3jet, 5);

    // DJALOG Making 4 and 5 jet histogram binning the same as 3 jet
    double jetPt_Zinc4jet[9] = {20, 24, 30, 41, 59, 83, 118, 168, 250};
    int nJetPt_Zinc4jet(sizeof(jetPt_Zinc4jet) / sizeof(jetPt_Zinc4jet[0]) - 1);
    double jetPt_Zinc5jet[9] = {20, 24, 30, 41, 59, 83, 118, 168, 250};
    int nJetPt_Zinc5jet(sizeof(jetPt_Zinc5jet) / sizeof(jetPt_Zinc5jet[0]) - 1);

    // double jetPt_Zinc4jet[9]  = {20, 24, 30, 39, 49, 62, 78, 96, 150};
    // int nJetPt_Zinc4jet(8);
    // int nJetPt_Zinc5jet(6);
    // double jetPt_Zinc5jet[7]  =  {20, 24, 30, 39, 49, 62, 100};

    // double jetPt_Zinc1jet[23] = {39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350,
    // 400, 450, 500, 590, 700, 1000};
    // double jetPt_Zinc1jet[9] = {20, 29, 41, 59, 83, 118, 168, 250, 350};

    // int nJetHT_Zinc1jet(17);
    // double jetHT_Zinc1jet[18] = {30, 39, 49, 62, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540,
    // 650, 800, 1000, 1500};
    int nJetHT_Zinc1jet(11);
    // double jetHT_Zinc1jet[10] = {20, 30, 41, 59, 83, 118, 168, 250, 350, 450, 540};
    double jetHT_Zinc1jet[12] = {30, 41, 59, 83, 118, 168, 220, 300, 400, 550, 780, 1100};

    vector<double> jetHT_2_Zinc1jet;
    jetHT_2_Zinc1jet = buildVecFineBin(nJetHT_Zinc1jet, jetHT_Zinc1jet, 5);

    // int nJetHT_Zinc2jet(13);
    // double jetHT_Zinc2jet[14] = {60, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800,
    // 1200};
    int nJetHT_Zinc2jet(9);
    double jetHT_Zinc2jet[10] = {60, 83, 118, 168, 220, 300, 400, 550, 780, 1100};

    vector<double> jetHT_2_Zinc2jet;
    jetHT_2_Zinc2jet = buildVecFineBin(nJetHT_Zinc2jet, jetHT_Zinc2jet, 5);

    int nJetHT_Zinc3jet(8);
    double jetHT_Zinc3jet[9] = {90, 130, 168, 220, 300, 400, 550, 780, 1100};
    //  double jetHT_Zinc3jet[10] = {30, 41, 59, 83, 118, 168, 220, 300, 400, 550};

    vector<double> jetHT_2_Zinc3jet;
    jetHT_2_Zinc3jet = buildVecFineBin(nJetHT_Zinc3jet, jetHT_Zinc3jet, 5);

    int nJetHT_Zinc4jet(9);
    double jetHT_Zinc4jet[10] = {120, 140, 167, 203, 253, 320, 410, 530, 690, 910};

    int nJetHT_Zinc5jet(6);
    double jetHT_Zinc5jet[7] = {150, 200, 282, 365, 485, 650, 880};

    int nJetsMass_Zinc2jet(19);
    double jetsMass_Zinc2jet[20] = {0,   25,  52,  81,  112, 145, 180, 217, 256, 297,
                                    340, 385, 432, 481, 532, 585, 640, 700, 830, 1000};

    Input = newTH1D("input", "", "", 1, 0, 1);
    const char *JobInfo_Labels[] = {
        "JobNum", "nJobs", "JobWeight", "nEvtsJob", "nEvtsAllJobs", "nEvtsSample", "Xsec", "Lumi"};
    unsigned nBins = sizeof(JobInfo_Labels) / sizeof(JobInfo_Labels[0]);
    JobInfo = newTH1D("JobInfo", "", "", nBins, -0.5, -0.5 + nBins);
    JobInfo->SetTitle("Information to merge histograms produced in multi job mode");
    JobInfo->SetBit(TH1::kIsAverage);
    for (unsigned i = 0; i < nBins; ++i) {
        JobInfo->GetXaxis()->SetBinLabel(i + 1, JobInfo_Labels[i]);
    }

    Lumi = newTH1D("Lumi", "Integrated luminosity (fb^{-1})", "", 1, 0, 1);
    Lumi->SetBit(TH1::kIsAverage);

    NumberPFcandidates = newTH1D("NumberPFcandidates",
                                 "NumberPFcandidates",
                                 "Number of lepton PF candidates",
                                 20,
                                 -0.5,
                                 19.5);

    NumberOfEvents = newTH1D("NumberOfEvents",
                             "NumberOfEvents",
                             "Number of events after various selection",
                             10,
                             -0.5,
                             9.5);

    ZMass_lowDeltaR = newTH1D("ZMass_lowDeltaR", "ZMass_lowDeltaR", Mll, 120, 50, 169);
    ZMass_Zinc0jet =
        newTH1D("ZMass_Zinc0jet", "Z Invariant Mass (N_{jets} #geq 0)", Mll, 210, 50, 260);
    ZMass_Zinc1jet =
        newTH1D("ZMass_Zinc1jet", "Z Invariant Mass (N_{jets} #geq 1)", Mll, 210, 50, 260);
    ZMassFrom60_Zinc0jet =
        newTH1D("ZMassFrom60_Zinc0jet", "Z Invariant Mass (N_{jets} #geq 0)", Mll, 300, 60, 660);

    genZMass_Zinc0jet =
        newTH1D("genZMass_Zinc0jet", "Z Invariant Mass (N_{jets} #geq 0)", Mll, 111, 50, 260);

    Mass_Zinc0jet = newTH1D(
        "Mass_Zinc0jet", "Invariant Mass (N_{jets} #geq 0)", Mll, nMass_Zinc0jet, mass_Zinc0jet);

    //--------------- phistar -------------

   Phistar_Zinc0jetM50_76 = newTH1D("Phistar_Zinc0jetM50_76",
                               "#phi^{*}_{#eta} (N_{jets} #geq 0, 50 < M < 71 GeV)",
                               phistar,
                               nPhistar_Zinc0jet,
                               phistar_Zinc0jet);
    genPhistar_Zinc0jetM50_76 = newTH1D("genPhistar_Zinc0jetM50_76",
                                  " gen #phi^{*}_{#eta} (N_{jets} #geq 0, 50 < M < 71 GeV)",
                                  phistar,
                                  nPhistar_Zinc0jet,
                                  phistar_Zinc0jet);
    hresponsePhistar_Zinc0jetM50_76 = newTH2D("hresponsePhistar_Zinc0jetM50_76",
                                        "response #phi^{*}_{#eta} (N_{jets} #geq 0, 50 < M < 71 GeV)",
                                        nPhistar_Zinc0jet,
                                        phistar_Zinc0jet,
                                        nPhistar_Zinc0jet,
                                        phistar_Zinc0jet);


  // M71_111
    Phistar_Zinc0jet = newTH1D("Phistar_Zinc0jet",
                               "#phi^{*}_{#eta} (N_{jets} #geq 0)",
                               phistar,
                               nPhistar_Zinc0jet,
                               phistar_Zinc0jet);
    genPhistar_Zinc0jet = newTH1D("genPhistar_Zinc0jet",
                                  " gen #phi^{*}_{#eta} (N_{jets} #geq 0)",
                                  phistar,
                                  nPhistar_Zinc0jet,
                                  phistar_Zinc0jet);
    hresponsePhistar_Zinc0jet = newTH2D("hresponsePhistar_Zinc0jet",
                                        "response #phi^{*}_{#eta} (N_{jets} #geq 0)",
                                        nPhistar_Zinc0jet,
                                        phistar_Zinc0jet,
                                        nPhistar_Zinc0jet,
                                        phistar_Zinc0jet);

 //  M76_106
    Phistar_Zinc0jetM76_106 = newTH1D("Phistar_Zinc0jetM76_106",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 0, 76 < M < 106 GeV)",
                                       phistar,
                                       nPhistar_Zinc0jet,
                                       phistar_Zinc0jet);
    genPhistar_Zinc0jetM76_106 =
        newTH1D("genPhistar_Zinc0jetM76_106",
                "gen #phi^{*}_{#eta} (N_{jets} #geq 0, 76 < M < 106 GeV)",
                phistar,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);
    hresponsePhistar_Zinc0jetM76_106 =
        newTH2D("hresponsePhistar_Zinc0jetM76_106",
                "response #phi^{*}_{#eta} (N_{jets} #geq 0, 76 < M < 106 GeV)",
                nPhistar_Zinc0jet,
                phistar_Zinc0jet,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);

  //  M111_130
    Phistar_Zinc0jetM111_130 = newTH1D("Phistar_Zinc0jetM111_130",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 0, 111 < M < 130 GeV)",
                                       phistar,
                                       nPhistar_Zinc0jet,
                                       phistar_Zinc0jet);
    genPhistar_Zinc0jetM111_130 =
        newTH1D("genPhistar_Zinc0jetM111_130",
                "gen #phi^{*}_{#eta} (N_{jets} #geq 0, 111 < M < 130 GeV)",
                phistar,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);
    hresponsePhistar_Zinc0jetM111_130 =
        newTH2D("hresponsePhistar_Zinc0jetM111_130",
                "response #phi^{*}_{#eta} (N_{jets} #geq 0, 111 < M < 130 GeV)",
                nPhistar_Zinc0jet,
                phistar_Zinc0jet,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);

  //  M106_170
    Phistar_Zinc0jetM106_170 = newTH1D("Phistar_Zinc0jetM106_170",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 0, 106 < M < 170 GeV)",
                                       phistar,
                                       nPhistar_Zinc0jet,
                                       phistar_Zinc0jet);
    genPhistar_Zinc0jetM106_170 =
        newTH1D("genPhistar_Zinc0jetM106_170",
                "gen #phi^{*}_{#eta} (N_{jets} #geq 0, 106 < M < 170 GeV)",
                phistar,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);

    hresponsePhistar_Zinc0jetM106_170 =
        newTH2D("hresponsePhistar_Zinc0jetM106_170",
                "response #phi^{*}_{#eta} (N_{jets} #geq 0, 106 < M < 170 GeV)",
                nPhistar_Zinc0jet,
                phistar_Zinc0jet,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);

    //  M130_170
    Phistar_Zinc0jetM130_170 = newTH1D("Phistar_Zinc0jetM130_170",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 0, 130 < M < 170 GeV)",
                                       phistar,
                                       nPhistar_Zinc0jet,
                                       phistar_Zinc0jet);
    genPhistar_Zinc0jetM130_170 =
        newTH1D("genPhistar_Zinc0jetM130_170",
                "gen #phi^{*}_{#eta} (N_{jets} #geq 0, 130 < M < 170 GeV)",
                phistar,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);
    hresponsePhistar_Zinc0jetM130_170 =
        newTH2D("hresponsePhistar_Zinc0jetM130_170",
                "response #phi^{*}_{#eta} (N_{jets} #geq 0, 130 < M < 170 GeV)",
                nPhistar_Zinc0jet,
                phistar_Zinc0jet,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);


  //  M170_250
    Phistar_Zinc0jetM170_250 = newTH1D("Phistar_Zinc0jetM170_250",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 0, 170 < M < 250 GeV)",
                                       phistar,
                                       nPhistar_Zinc0jet,
                                       phistar_Zinc0jet);

    genPhistar_Zinc0jetM170_250 =
        newTH1D("genPhistar_Zinc0jetM170_250",
                " gen #phi^{*}_{#eta} (N_{jets} #geq 0, 170 < M < 250 GeV)",
                phistar,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);
    hresponsePhistar_Zinc0jetM170_250 =
        newTH2D("hresponsePhistar_Zinc0jetM170_250",
                "response #phi^{*}_{#eta} (N_{jets} #geq 0,  170 < M < 250 GeV)",
                nPhistar_Zinc0jet,
                phistar_Zinc0jet,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);

   //  M170_350

    Phistar_Zinc0jetM170_350 = newTH1D("Phistar_Zinc0jetM170_350",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 0, 170 < M < 350 GeV)",
                                       phistar,
                                       nPhistar_Zinc0jet,
                                       phistar_Zinc0jet);

    genPhistar_Zinc0jetM170_350 =
        newTH1D("genPhistar_Zinc0jetM170_350",
                " gen #phi^{*}_{#eta} (N_{jets} #geq 0, 170 < M < 350 GeV)",
                phistar,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);

    hresponsePhistar_Zinc0jetM170_350 =
        newTH2D("hresponsePhistar_Zinc0jetM170_350",
                "response #phi^{*}_{#eta} (N_{jets} #geq 0,  170 < M < 350 GeV)",
                nPhistar_Zinc0jet,
                phistar_Zinc0jet,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);

   //  M170_inf

    Phistar_Zinc0jetM170_inf = newTH1D("Phistar_Zinc0jetM170_inf",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 0, M > 170 GeV)",
                                       phistar,
                                       nPhistar_Zinc0jet,
                                       phistar_Zinc0jet);

    genPhistar_Zinc0jetM170_inf =
        newTH1D("genPhistar_Zinc0jetM170_inf",
                " gen #phi^{*}_{#eta} (N_{jets} #geq 0, M > 170 GeV)",
                phistar,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);

    hresponsePhistar_Zinc0jetM170_inf =
        newTH2D("hresponsePhistar_Zinc0jetM170_inf",
                "response #phi^{*}_{#eta} (N_{jets} #geq 0,  M > 170 GeV)",
                nPhistar_Zinc0jet,
                phistar_Zinc0jet,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);


// M250_320

    Phistar_Zinc0jetM250_3 = newTH1D("Phistar_Zinc0jetM250_3",
                                     "#phi^{*}_{#eta} (N_{jets} #geq 0, 250 < M < 3000 GeV)",
                                     phistar,
                                     nPhistar_Zinc0jet,
                                     phistar_Zinc0jet);
    genPhistar_Zinc0jetM250_3 =
        newTH1D("genPhistar_Zinc0jetM250_3",
                " gen #phi^{*}_{#eta} (N_{jets} #geq 0, 250 < M < 3000 GeV)",
                phistar,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);
    hresponsePhistar_Zinc0jetM250_3 =
        newTH2D("hresponsePhistar_Zinc0jetM250_3",
                "response #phi^{*}_{#eta} (N_{jets} #geq 0,  250 < M < 3000 GeV)",
                nPhistar_Zinc0jet,
                phistar_Zinc0jet,
                nPhistar_Zinc0jet,
                phistar_Zinc0jet);

    // -- Pjistar 1 jet

   Phistar_Zinc1jetM50_76 = newTH1D("Phistar_Zinc1jetM50_76",
                               "#phi^{*}_{#eta} (N_{jets} #geq 1, 50 < M < 71 GeV)",
                               phistar,
                               nPhistar_Zinc1jetMbins,
                               phistar_Zinc1jetMbins);
    genPhistar_Zinc1jetM50_76 = newTH1D("genPhistar_Zinc1jetM50_76",
                                  " gen #phi^{*}_{#eta} (N_{jets} #geq 1, 50 < M < 71 GeV)",
                                  phistar,
                                  nPhistar_Zinc1jetMbins,
                                  phistar_Zinc1jetMbins);
    hresponsePhistar_Zinc1jetM50_76 = newTH2D("hresponsePhistar_Zinc1jetM50_76",
                                        "response #phi^{*}_{#eta} (N_{jets} #geq 1, 50 < M < 71 GeV)",
                                        nPhistar_Zinc1jetMbins,
                                        phistar_Zinc1jetMbins,
                                        nPhistar_Zinc1jetMbins,
                                        phistar_Zinc1jetMbins);




    Phistar_Zinc1jet = newTH1D("Phistar_Zinc1jet",
                               "#phi^{*}_{#eta} (N_{jets} #geq 1)",
                               phistar,
                               nPhistar_Zinc1jet,
                               phistar_Zinc1jet);

    genPhistar_Zinc1jet = newTH1D("genPhistar_Zinc1jet",
                                  " gen #phi^{*}_{#eta} (N_{jets} #geq 1)",
                                  phistar,
                                  nPhistar_Zinc1jet,
                                  phistar_Zinc1jet);
    hresponsePhistar_Zinc1jet = newTH2D("hresponsePhistar_Zinc1jet",
                                        "response #phi^{*}_{#eta} (N_{jets} #geq 1)",
                                        nPhistar_Zinc1jet,
                                        phistar_Zinc1jet,
                                        nPhistar_Zinc1jet,
                                        phistar_Zinc1jet);

    // M111_130
    Phistar_Zinc1jetM111_130 = newTH1D("Phistar_Zinc1jetM111_130",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 0, 111 < M < 130 GeV)",
                                       phistar,
                                       nPhistar_Zinc1jetMbins,
                                       phistar_Zinc1jetMbins);
    genPhistar_Zinc1jetM111_130 =
        newTH1D("genPhistar_Zinc1jetM111_130",
                "gen #phi^{*}_{#eta} (N_{jets} #geq 0, 111 < M < 130 GeV)",
                phistar,
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins);
    hresponsePhistar_Zinc1jetM111_130 =
        newTH2D("hresponsePhistar_Zinc1jetM111_130",
                "response #phi^{*}_{#eta} (N_{jets} #geq 0, 111 < M < 130 GeV)",
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins,
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins);

    // M130_170
    Phistar_Zinc1jetM130_170 = newTH1D("Phistar_Zinc1jetM130_170",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                       phistar,
                                       nPhistar_Zinc1jetMbins,
                                       phistar_Zinc1jetMbins);
    genPhistar_Zinc1jetM130_170 =
        newTH1D("genPhistar_Zinc1jetM130_170",
                " gen #phi^{*}_{#eta} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                phistar,
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins);
    hresponsePhistar_Zinc1jetM130_170 =
        newTH2D("hresponsePhistar_Zinc1jetM130_170",
                "response #phi^{*}_{#eta} (N_{jets} #geq 1,  130 < M < 170 GeV)",
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins,
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins);

    // M170_250
    Phistar_Zinc1jetM170_250 = newTH1D("Phistar_Zinc1jetM170_250",
                                       "#phi^{*}_{#eta} (N_{jets} #geq 1, 170 < M < 250 GeV)",
                                       phistar,
                                       nPhistar_Zinc1jetMbins,
                                       phistar_Zinc1jetMbins);
    genPhistar_Zinc1jetM170_250 =
        newTH1D("genPhistar_Zinc1jetM170_250",
                " gen #phi^{*}_{#eta} (N_{jets} #geq 1, 170 < M < 250 GeV)",
                phistar,
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins);
    hresponsePhistar_Zinc1jetM170_250 =
        newTH2D("hresponsePhistar_Zinc1jetM170_250",
                "response #phi^{*}_{#eta} (N_{jets} #geq 1,  170 < M < 250 GeV)",
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins,
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins);

    // M250_3
    Phistar_Zinc1jetM250_3 = newTH1D("Phistar_Zinc1jetM250_3",
                                     "#phi^{*}_{#eta} (N_{jets} #geq 1, 250 < M < 3000 GeV)",
                                     phistar,
                                     nPhistar_Zinc1jetMbins,
                                     phistar_Zinc1jetMbins);
    genPhistar_Zinc1jetM250_3 =
        newTH1D("genPhistar_Zinc1jetM250_3",
                " gen #phi^{*}_{#eta} (N_{jets} #geq 1, 250 < M < 3000 GeV)",
                phistar,
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins);
    hresponsePhistar_Zinc1jetM250_3 =
        newTH2D("hresponsePhistar_Zinc1jetM250_3",
                "response #phi^{*}_{#eta} (N_{jets} #geq 1,  250 < M < 3000 GeV)",
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins,
                nPhistar_Zinc1jetMbins,
                phistar_Zinc1jetMbins);

    //-- 2D
    Phistar_Zpt = newTH2D("Phistar_Zpt",
                          "Phistar vs Zpt",
                          nZPt_Zinc0jet,
                          zPt_Zinc0jet, // x-axis
                          nPhistar_Zinc0jet,
                          phistar_Zinc0jet); // y-axis

    Phistar_Zpt_test = newTH2D("Phistar_Zpt_test",
                               "Phistar vs Zpt(phistar)",
                               nZPt_Zinc0jet,
                               zPt_Zinc0jet, // x-axis
                               nPhistar_Zinc0jet,
                               phistar_Zinc0jet); // y-axis
    //--------------- Zpt --------
 
    // nZPt_Zinc0jetM15_50
    ZPt_Zinc0jetM15_50 = newTH1D("ZPt_Zinc0jetM15_50",
                                 "Z p_{T} (N_{jets} #geq 0 15 < M < 50 GeV)",
                                 ZpT,
                                 nZPt_Zinc0jetMbins,
                                 zPt_Zinc0jetMbins);

    genZPt_Zinc0jetM15_50 = newTH1D("genZPt_Zinc0jetM15_50",
                                    "gen Z p_{T} (N_{jets} #geq 0 15 < M < 50 GeV)",
                                    ZpT,
                                    nZPt_Zinc0jetMbins,
                                    zPt_Zinc0jetMbins);
    hresponseZPt_Zinc0jetM15_50 = newTH2D("hresponseZPt_Zinc0jetM15_50",
                                          "response Z p_{T} (N_{jets} #geq 0 15 < M < 50 GeV)",
                                          nZPt_Zinc0jetMbins,
                                          zPt_Zinc0jetMbins,
                                          nZPt_Zinc0jetMbins,
                                          zPt_Zinc0jetMbins);

    // nZPt_Zinc0jetM50_76
    ZPt_Zinc0jetM50_76 = newTH1D("ZPt_Zinc0jetM50_76",
                                 "Z p_{T} (N_{jets} #geq 0 50 < M < 71 GeV)",
                                 ZpT,
                                 nZPt_Zinc0jetMbins,
                                 zPt_Zinc0jetMbins);

    genZPt_Zinc0jetM50_76 = newTH1D("genZPt_Zinc0jetM50_76",
                                    "gen Z p_{T} (N_{jets} #geq 0 50 < M < 71 GeV)",
                                    ZpT,
                                    nZPt_Zinc0jetMbins,
                                    zPt_Zinc0jetMbins);
    hresponseZPt_Zinc0jetM50_76 = newTH2D("hresponseZPt_Zinc0jetM50_76",
                                          "response Z p_{T} (N_{jets} #geq 0 50 < M < 71 GeV)",
                                          nZPt_Zinc0jetMbins,
                                          zPt_Zinc0jetMbins,
                                          nZPt_Zinc0jetMbins,
                                          zPt_Zinc0jetMbins);

    // Z peak 71 -111

   ZPt_Zinc0jet =
        newTH1D("ZPt_Zinc0jet", "Z p_{T} (N_{jets} #geq 0)", ZpT, nZPt_Zinc0jet, zPt_Zinc0jet);
    ZPt_Zinc0jet_new = newTH1D("ZPt_Zinc0jet_new",
                               "Z p_{T} via #phi^{*} (N_{jets} #geq 0)",
                               ZpT,
                               nZPt_Zinc0jet,
                               zPt_Zinc0jet);
    genZPt_Zinc0jet = newTH1D(
        "genZPt_Zinc0jet", "gen Z p_{T} (N_{jets} #geq 0)", ZpT, nZPt_Zinc0jet, zPt_Zinc0jet);
    genZPt_Zinc0jet_new = newTH1D(
        "genZPt_Zinc0jet_new", "gen Z p_{T} (N_{jets} #geq 0)", ZpT, nZPt_Zinc0jet, zPt_Zinc0jet);


    hresponseZPt_Zinc0jet = newTH2D("hresponseZPt_Zinc0jet",
                                    "response Z p_{T} (N_{jets} #geq 0)",
                                    nZPt_Zinc0jet,
                                    zPt_Zinc0jet,
                                    nZPt_Zinc0jet,
                                    zPt_Zinc0jet);

    hresponseZPt_Zinc0jet_lowNVtx = newTH2D("hresponseZPt_Zinc0jet_lowNVtx",
                                            "response Z p_{T} (N_{jets} #geq 0) lowNVtx",
                                            nZPt_Zinc0jet,
                                            zPt_Zinc0jet,
                                            nZPt_Zinc0jet,
                                            zPt_Zinc0jet);

    hresponseZPt_Zinc0jet_highNVtx = newTH2D("hresponseZPt_Zinc0jet_highNVtx",
                                             "response Z p_{T} (N_{jets} #geq 0) highNVtx",
                                             nZPt_Zinc0jet,
                                             zPt_Zinc0jet,
                                             nZPt_Zinc0jet,
                                             zPt_Zinc0jet);

    hresponseZPt_Zinc0jet_new = newTH2D("hresponseZPt_Zinc0jet_new",
                                        "response Z p_{T} (N_{jets} #geq 0)",
                                        nPhistar_Zinc0jet,
                                        phistar_Zinc0jet,
                                        nZPt_Zinc0jet,
                                        zPt_Zinc0jet);

    // Z peak 76 -106

     ZPt_Zinc0jetM76_106 =
        newTH1D("ZPt_Zinc0jetM76_106", "Z p_{T} (N_{jets} #geq 0)", ZpT, nZPt_Zinc0jet, zPt_Zinc0jet);

     genZPt_Zinc0jetM76_106 = newTH1D(
        "genZPt_Zinc0jetM76_106", "gen Z p_{T} (N_{jets} #geq 0)", ZpT, nZPt_Zinc0jet, zPt_Zinc0jet);

     hresponseZPt_Zinc0jetM76_106 = newTH2D("hresponseZPt_Zinc0jetM76_106",
                                    "response Z p_{T} (N_{jets} #geq 0)",
                                    nZPt_Zinc0jet,
                                    zPt_Zinc0jet,
                                    nZPt_Zinc0jet,
                                    zPt_Zinc0jet);

    // Z peak 76 -106 bins as for high Mass 

     ZPt_Zinc0jetM76_106_Mbin =
        newTH1D("ZPt_Zinc0jetM76_106_Mbin", "Z p_{T} (N_{jets} #geq 0)", ZpT,
                                   nZPt_Zinc0jetMbins,
                                   zPt_Zinc0jetMbins);

     genZPt_Zinc0jetM76_106_Mbin = newTH1D(
        "genZPt_Zinc0jetM76_106_Mbin", "gen Z p_{T} (N_{jets} #geq 0)", ZpT,
                                   nZPt_Zinc0jetMbins,
                                   zPt_Zinc0jetMbins);

     hresponseZPt_Zinc0jetM76_106_Mbin = newTH2D("hresponseZPt_Zinc0jetM76_106_Mbin",
                                    "response Z p_{T} (N_{jets} #geq 0)",
                                           nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins,
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins);

    // M 106 -170

    ZPt_Zinc0jetM106_170 = newTH1D("ZPt_Zinc0jetM106_170",
                                   "Z p_{T} (N_{jets} #geq 0 106 < M < 170 GeV)",
                                   ZpT,
                                   nZPt_Zinc0jetMbins,
                                   zPt_Zinc0jetMbins);

    genZPt_Zinc0jetM106_170 = newTH1D("genZPt_Zinc0jetM106_170",
                                      "gen Z p_{T} (N_{jets} #geq 0 106 < M < 170 GeV)",
                                      ZpT,
                                      nZPt_Zinc0jetMbins,
                                      zPt_Zinc0jetMbins);
    hresponseZPt_Zinc0jetM106_170 = newTH2D("hresponseZPt_Zinc0jetM106_170",
                                            "response Z p_{T} (N_{jets} #geq 0 106 < M < 170 GeV)",
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins,
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins);

    // nZPt_Zinc0jetM111_130
    ZPt_Zinc0jetM111_130 = newTH1D("ZPt_Zinc0jetM111_130",
                                   "Z p_{T} (N_{jets} #geq 0 111 < M < 130 GeV)",
                                   ZpT,
                                   nZPt_Zinc0jetMbins,
                                   zPt_Zinc0jetMbins);

    genZPt_Zinc0jetM111_130 = newTH1D("genZPt_Zinc0jetM111_130",
                                      "gen Z p_{T} (N_{jets} #geq 0 111 < M < 130 GeV)",
                                      ZpT,
                                      nZPt_Zinc0jetMbins,
                                      zPt_Zinc0jetMbins);
    hresponseZPt_Zinc0jetM111_130 = newTH2D("hresponseZPt_Zinc0jetM111_130",
                                            "response Z p_{T} (N_{jets} #geq 0 111 < M < 130 GeV)",
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins,
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins);

    // M 170 -350

    ZPt_Zinc0jetM170_350 = newTH1D("ZPt_Zinc0jetM170_350",
                                   "Z p_{T} (N_{jets} #geq 0 170 < M < 350 GeV)",
                                   ZpT,
                                   nZPt_Zinc0jetMbins,
                                   zPt_Zinc0jetMbins);

    genZPt_Zinc0jetM170_350 = newTH1D("genZPt_Zinc0jetM170_350",
                                      "gen Z p_{T} (N_{jets} #geq 0 170 < M < 350 GeV)",
                                      ZpT,
                                      nZPt_Zinc0jetMbins,
                                      zPt_Zinc0jetMbins);
    hresponseZPt_Zinc0jetM170_350 = newTH2D("hresponseZPt_Zinc0jetM170_350",
                                            "response Z p_{T} (N_{jets} #geq 0 170 < M < 350 GeV)",
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins,
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins);

    // M 170 inf

    ZPt_Zinc0jetM170_inf = newTH1D("ZPt_Zinc0jetM170_inf",
                                   "Z p_{T} (N_{jets} #geq 0, M > 170 GeV)",
                                   ZpT,
                                   nZPt_Zinc0jetMbins,
                                   zPt_Zinc0jetMbins);

    genZPt_Zinc0jetM170_inf = newTH1D("genZPt_Zinc0jetM170_inf",
                                      "gen Z p_{T} (N_{jets} #geq 0, M > 170 GeV)",
                                      ZpT,
                                      nZPt_Zinc0jetMbins,
                                      zPt_Zinc0jetMbins);
    hresponseZPt_Zinc0jetM170_inf = newTH2D("hresponseZPt_Zinc0jetM170_inf",
                                            "response Z p_{T} (N_{jets} #geq 0, M > 170 GeV)",
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins,
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins);


    // nZPt_Zinc0jetM130_170
    ZPt_Zinc0jetM130_170 = newTH1D("ZPt_Zinc0jetM130_170",
                                   "Z p_{T} (N_{jets} #geq 0 130 < M < 170 GeV)",
                                   ZpT,
                                   nZPt_Zinc0jetMbins,
                                   zPt_Zinc0jetMbins);

    genZPt_Zinc0jetM130_170 = newTH1D("genZPt_Zinc0jetM130_170",
                                      "gen Z p_{T} (N_{jets} #geq 0 130 < M < 170 GeV)",
                                      ZpT,
                                      nZPt_Zinc0jetMbins,
                                      zPt_Zinc0jetMbins);
    hresponseZPt_Zinc0jetM130_170 = newTH2D("hresponseZPt_Zinc0jetM130_170",
                                            "response Z p_{T} (N_{jets} #geq 0 130 < M < 170 GeV)",
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins,
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins);

    // nZPt_Zinc0jetM170_250
    ZPt_Zinc0jetM170_250 = newTH1D("ZPt_Zinc0jetM170_250",
                                   "Z p_{T} (N_{jets} #geq 0 170 < M < 250 GeV)",
                                   ZpT,
                                   nZPt_Zinc0jetMbins,
                                   zPt_Zinc0jetMbins);
    genZPt_Zinc0jetM170_250 = newTH1D("genZPt_Zinc0jetM170_250",
                                      "gen Z p_{T} (N_{jets} #geq 0 170 < M < 250 GeV)",
                                      ZpT,
                                      nZPt_Zinc0jetMbins,
                                      zPt_Zinc0jetMbins);
    ZPt_Zinc0jetM170_250_new = newTH1D("ZPt_Zinc0jetM170_250_new",
                                       "Z p_{T} (N_{jets} #geq 0 170 < M < 250 GeV)",
                                       ZpT,
                                       nZPt_Zinc0jetMbins,
                                       zPt_Zinc0jetMbins);
    genZPt_Zinc0jetM170_250_new = newTH1D("genZPt_Zinc0jetM170_250_new",
                                          "gen Z p_{T} (N_{jets} #geq 0 170 < M < 250 GeV)",
                                          ZpT,
                                          nZPt_Zinc0jetMbins,
                                          zPt_Zinc0jetMbins);
    hresponseZPt_Zinc0jetM170_250 = newTH2D("hresponseZPt_Zinc0jetM170_250",
                                            "response Z p_{T} (N_{jets} #geq 0 170 < M < 250 GeV)",
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins,
                                            nZPt_Zinc0jetMbins,
                                            zPt_Zinc0jetMbins);

    hresponseZPt_Zinc0jetM170_250_new =
        newTH2D("hresponseZPt_Zinc0jetM170_250_new",
                "response Z p_{T} (N_{jets} #geq 0 170 < M < 250 GeV)",
                nZPt_Zinc0jetMbins,
                zPt_Zinc0jetMbins,
                nZPt_Zinc0jetMbins,
                zPt_Zinc0jetMbins);

//----------- narrow binning ----------------

    ZPt_2_Zinc0jetM170_250 = newTH1D("ZPt_2_Zinc0jetM170_250",
                                   "Z p_{T} (N_{jets} #geq 0 170 < M < 250 GeV)",
                                   ZpT,
                                   zPt_2_Zinc0jetMbins);

    genZPt_2_Zinc0jetM170_250 = newTH1D("genZPt_2_Zinc0jetM170_250",
                                      "gen Z p_{T} (N_{jets} #geq 0 170 < M < 250 GeV)",
                                      ZpT,
                                      zPt_2_Zinc0jetMbins);


    hresponseZPt_2_Zinc0jetM170_250 = newTH2D("hresponseZPt_2_Zinc0jetM170_250",
                                            "response Z p_{T} (N_{jets} #geq 0 170 < M < 250 GeV)",
                                            zPt_2_Zinc0jetMbins,
                                            zPt_2_Zinc0jetMbins);

//-----------------------------

    // nZPt_Zinc0jetM250_3
    ZPt_Zinc0jetM250_3 = newTH1D("ZPt_Zinc0jetM250_3",
                                 "Z p_{T} (N_{jets} #geq 0 250 < M < 320 GeV)",
                                 ZpT,
                                 nZPt_Zinc0jetMbins250,
                                 zPt_Zinc0jetMbins250);


    genZPt_Zinc0jetM250_3 = newTH1D("genZPt_Zinc0jetM250_3",
                                    "gen Z p_{T} (N_{jets} #geq 0 250 < M < 320 GeV)",
                                    ZpT,
                                    nZPt_Zinc0jetMbins250,
                                    zPt_Zinc0jetMbins250);
    hresponseZPt_Zinc0jetM250_3 = newTH2D("hresponseZPt_Zinc0jetM250_3",
                                          "response Z p_{T} (N_{jets} #geq 0 250 < M < 320 GeV)",
                                          nZPt_Zinc0jetMbins250,
                                          zPt_Zinc0jetMbins250,
                                          nZPt_Zinc0jetMbins250,
                                          zPt_Zinc0jetMbins250);
//-------narrow binning----------------

    ZPt_2_Zinc0jetM250_3 = newTH1D("ZPt_2_Zinc0jetM250_3",
                                 "Z p_{T} (N_{jets} #geq 0 250 < M < 320 GeV)",
                                 ZpT,
                                 zPt_2_Zinc0jetMbins250);


    genZPt_2_Zinc0jetM250_3 = newTH1D("genZPt_2_Zinc0jetM250_3",
                                    "gen Z p_{T} (N_{jets} #geq 0 250 < M < 320 GeV)",
                                    ZpT,
                                    zPt_2_Zinc0jetMbins250);
    hresponseZPt_2_Zinc0jetM250_3 = newTH2D("hresponseZPt_2_Zinc0jetM250_3",
                                          "response Z p_{T} (N_{jets} #geq 0 250 < M < 320 GeV)",
                                          zPt_2_Zinc0jetMbins250,
                                          zPt_2_Zinc0jetMbins250);

//--------------------


    hresponseZPt_Zinc0jetM250_3_dR =
        newTH2D("hresponseZPt_Zinc0jetM250_3_dR",
                "response Z p_{T} (N_{jets} #geq 0 250 < M < 320 GeV, dR(muRec,muGen)<0.1)",
                nZPt_Zinc0jetMbins250,
                zPt_Zinc0jetMbins250,
                nZPt_Zinc0jetMbins250,
                zPt_Zinc0jetMbins250);
    hresponseZPt_Zinc0jetM250_3_dPt =
        newTH2D("hresponseZPt_Zinc0jetM250_3_dPt",
                "response Z p_{T} (N_{jets} #geq 0 250 < M < 320 GeV, dPt(murec-mugen)<0.1)",
                nZPt_Zinc0jetMbins250,
                zPt_Zinc0jetMbins250,
                nZPt_Zinc0jetMbins250,
                zPt_Zinc0jetMbins250);

    ZPt_Zinc0jetM320_3 = newTH1D("ZPt_Zinc0jetM320_3",
                                 "Z p_{T} (N_{jets} #geq 0 320 < M < 3000 GeV)",
                                 ZpT,
                                 nZPt_Zinc0jetMbins320,
                                 zPt_Zinc0jetMbins320);

    genZPt_Zinc0jetM320_3 = newTH1D("genZPt_Zinc0jetM320_3",
                                    "gen Z p_{T} (N_{jets} #geq 0 320 < M < 3000 GeV)",
                                    ZpT,
                                    nZPt_Zinc0jetMbins320,
                                    zPt_Zinc0jetMbins320);
    hresponseZPt_Zinc0jetM320_3 = newTH2D("hresponseZPt_Zinc0jetM320_3",
                                          "response Z p_{T} (N_{jets} #geq 0 320 < M < 3000 GeV)",
                                          nZPt_Zinc0jetMbins320,
                                          zPt_Zinc0jetMbins320,
                                          nZPt_Zinc0jetMbins320,
                                          zPt_Zinc0jetMbins320);

    // nZPt_Zinc0jetM115_135 (for Higgs comparison)
    ZPt_Zinc0jetM115_135 = newTH1D("ZPt_Zinc0jetM115_135",
                                   "Z p_{T} (N_{jets} #geq 0 115 < M < 135 GeV)",
                                   ZpT,
                                   nZPt_Zinc0jetM115_135,
                                   zPt_Zinc0jetM115_135);
    genZPt_Zinc0jetM115_135 = newTH1D("genZPt_Zinc0jetM115_135",
                                      "gen Z p_{T} (N_{jets} #geq 0 115 < M < 135 GeV)",
                                      ZpT,
                                      nZPt_Zinc0jetM115_135,
                                      zPt_Zinc0jetM115_135);
    hresponseZPt_Zinc0jetM115_135 = newTH2D("hresponseZPt_Zinc0jetM115_135",
                                            "response Z p_{T} (N_{jets} #geq 0 115 < M < 135 GeV)",
                                            nZPt_Zinc0jetM115_135,
                                            zPt_Zinc0jetM115_135,
                                            nZPt_Zinc0jetM115_135,
                                            zPt_Zinc0jetM115_135);

    //--------------- 1 jet ----------------------------------------



    // nZPt_Zinc1jetM50_76
    ZPt_Zinc1jetM50_76 = newTH1D("ZPt_Zinc1jetM50_76",
                                 "Z p_{T} (N_{jets} #geq 1, 50 < M < 71 GeV)",
                                 ZpT,
                                 nZPt_Zinc1jetMbins,
                                 zPt_Zinc1jetMbins);

    genZPt_Zinc1jetM50_76 = newTH1D("genZPt_Zinc1jetM50_76",
                                    "gen Z p_{T} (N_{jets} #geq 1, 50 < M < 71 GeV)",
                                    ZpT,
                                    nZPt_Zinc1jetMbins,
                                    zPt_Zinc1jetMbins);

    hresponseZPt_Zinc1jetM50_76 = newTH2D("hresponseZPt_Zinc1jetM50_76",
                                          "response Z p_{T} (N_{jets} #geq 1, 50 < M < 71 GeV)",
                                          nZPt_Zinc1jetMbins,
                                          zPt_Zinc1jetMbins,
                                          nZPt_Zinc1jetMbins,
                                          zPt_Zinc1jetMbins);



    // Z peak 71-111

    ZPt_Zinc1jet =
        newTH1D("ZPt_Zinc1jet", "Z p_{T} (N_{jets} #geq 1)", ZpT, nZPt_Zinc0jet, zPt_Zinc0jet);
    genZPt_Zinc1jet = newTH1D(
        "genZPt_Zinc1jet", "gen Z p_{T} (N_{jets} #geq 1)", ZpT, nZPt_Zinc0jet, zPt_Zinc0jet);

    hresponseZPt_Zinc1jet = newTH2D("hresponseZPt_Zinc1jet",
                                    "response Z p_{T} (N_{jets} #geq 1)",
                                    nZPt_Zinc0jet,
                                    zPt_Zinc0jet,
                                    nZPt_Zinc0jet,
                                    zPt_Zinc0jet);

    // Z peak 76-106

    ZPt_Zinc1jetM76_106 =
        newTH1D("ZPt_Zinc1jetM76_106", "Z p_{T} (N_{jets} #geq 1)", ZpT, nZPt_Zinc0jet, zPt_Zinc0jet);
    genZPt_Zinc1jetM76_106 = newTH1D(
        "genZPt_Zinc1jetM76_106", "gen Z p_{T} (N_{jets} #geq 1)", ZpT, nZPt_Zinc0jet, zPt_Zinc0jet);

    hresponseZPt_Zinc1jetM76_106 = newTH2D("hresponseZPt_Zinc1jetM76_106",
                                    "response Z p_{T} (N_{jets} #geq 1)",
                                    nZPt_Zinc0jet,
                                    zPt_Zinc0jet,
                                    nZPt_Zinc0jet,
                                    zPt_Zinc0jet);

    // Z peak 76-106 bins as in high Mass region

    ZPt_Zinc1jetM76_106_Mbin =
        newTH1D("ZPt_Zinc1jetM76_106_Mbin", "Z p_{T} (N_{jets} #geq 1)",  ZpT,
                                   nZPt_Zinc1jetMbins,
                                   zPt_Zinc1jetMbins);
    genZPt_Zinc1jetM76_106_Mbin = newTH1D(
        "genZPt_Zinc1jetM76_106_Mbin", "gen Z p_{T} (N_{jets} #geq 1)",  ZpT,
                                   nZPt_Zinc1jetMbins,
                                   zPt_Zinc1jetMbins);

    hresponseZPt_Zinc1jetM76_106_Mbin = newTH2D("hresponseZPt_Zinc1jetM76_106_Mbin",
                                    "response Z p_{T} (N_{jets} #geq 1)",
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins,
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins);


    // M106_170
    ZPt_Zinc1jetM106_170 = newTH1D("ZPt_Zinc1jetM106_170",
                                   "Z p_{T} (N_{jets} #geq 1, 106 < M < 170 GeV)",
                                   ZpT,
                                   nZPt_Zinc1jetMbins,
                                   zPt_Zinc1jetMbins);
    genZPt_Zinc1jetM106_170 = newTH1D("genZPt_Zinc1jetM106_170",
                                      "gen Z p_{T} (N_{jets} #geq 1, 106 < M < 170 GeV)",
                                      ZpT,
                                      nZPt_Zinc1jetMbins,
                                      zPt_Zinc1jetMbins);
    hresponseZPt_Zinc1jetM106_170 = newTH2D("hresponseZPt_Zinc1jetM106_170",
                                            "response Z p_{T} (N_{jets} #geq 1, 106 < M < 170 GeV)",
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins,
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins);




    // M111_130
    ZPt_Zinc1jetM111_130 = newTH1D("ZPt_Zinc1jetM111_130",
                                   "Z p_{T} (N_{jets} #geq 1, 111 < M < 130 GeV)",
                                   ZpT,
                                   nZPt_Zinc1jetMbins,
                                   zPt_Zinc1jetMbins);
    genZPt_Zinc1jetM111_130 = newTH1D("genZPt_Zinc1jetM111_130",
                                      "gen Z p_{T} (N_{jets} #geq 1, 111 < M < 130 GeV)",
                                      ZpT,
                                      nZPt_Zinc1jetMbins,
                                      zPt_Zinc1jetMbins);
    hresponseZPt_Zinc1jetM111_130 = newTH2D("hresponseZPt_Zinc1jetM111_130",
                                            "response Z p_{T} (N_{jets} #geq 1, 111 < M < 130 GeV)",
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins,
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins);

    // M130_170
    ZPt_Zinc1jetM130_170 = newTH1D("ZPt_Zinc1jetM130_170",
                                   "Z p_{T} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                   ZpT,
                                   nZPt_Zinc1jetHighMbins,
                                   zPt_Zinc1jetHighMbins);
    genZPt_Zinc1jetM130_170 = newTH1D("genZPt_Zinc1jetM130_170",
                                      "gen Z p_{T} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                      ZpT,
                                      nZPt_Zinc1jetHighMbins,
                                      zPt_Zinc1jetHighMbins);
    hresponseZPt_Zinc1jetM130_170 = newTH2D("hresponseZPt_Zinc1jetM130_170",
                                            "response Z p_{T} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                            nZPt_Zinc1jetHighMbins,
                                            zPt_Zinc1jetHighMbins,
                                            nZPt_Zinc1jetHighMbins,
                                            zPt_Zinc1jetHighMbins);

    // M170_250
    ZPt_Zinc1jetM170_250 = newTH1D("ZPt_Zinc1jetM170_250",
                                   "Z p_{T} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                   ZpT,
                                   nZPt_Zinc1jetMbins,
                                   zPt_Zinc1jetMbins);
    genZPt_Zinc1jetM170_250 = newTH1D("genZPt_Zinc1jetM170_250",
                                      "gen Z p_{T} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                      ZpT,
                                      nZPt_Zinc1jetMbins,
                                      zPt_Zinc1jetMbins);
    hresponseZPt_Zinc1jetM170_250 = newTH2D("hresponseZPt_Zinc1jetM170_250",
                                            "response Z p_{T} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins,
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins);


   // M170_350
    ZPt_Zinc1jetM170_350 = newTH1D("ZPt_Zinc1jetM170_350",
                                   "Z p_{T} (N_{jets} #geq 1, 170 < M < 350 GeV)",
                                   ZpT,
                                   nZPt_Zinc1jetMbins,
                                   zPt_Zinc1jetMbins);
    genZPt_Zinc1jetM170_350 = newTH1D("genZPt_Zinc1jetM170_350",
                                      "gen Z p_{T} (N_{jets} #geq 1, 170 < M < 350 GeV)",
                                      ZpT,
                                      nZPt_Zinc1jetMbins,
                                      zPt_Zinc1jetMbins);
    hresponseZPt_Zinc1jetM170_350 = newTH2D("hresponseZPt_Zinc1jetM170_350",
                                            "response Z p_{T} (N_{jets} #geq 1, 170 < M < 350 GeV)",
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins,
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins);

   // M170_inf
    ZPt_Zinc1jetM170_inf = newTH1D("ZPt_Zinc1jetM170_inf",
                                   "Z p_{T} (N_{jets} #geq 1, M > 170 GeV)",
                                   ZpT,
                                   nZPt_Zinc1jetMbins,
                                   zPt_Zinc1jetMbins);
    genZPt_Zinc1jetM170_inf = newTH1D("genZPt_Zinc1jetM170_inf",
                                      "gen Z p_{T} (N_{jets} #geq 1, M > 170 GeV)",
                                      ZpT,
                                      nZPt_Zinc1jetMbins,
                                      zPt_Zinc1jetMbins);
    hresponseZPt_Zinc1jetM170_inf = newTH2D("hresponseZPt_Zinc1jetM170_inf",
                                            "response Z p_{T} (N_{jets} #geq 1, M > 170 GeV)",
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins,
                                            nZPt_Zinc1jetMbins,
                                            zPt_Zinc1jetMbins);
//---------------- narrow binning
/*
    ZPt_2_Zinc1jetM170_250 = newTH1D("ZPt_2_Zinc1jetM170_250",
                                   "Z p_{T} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                   ZpT,
                                   nZPt_Zinc1jetMbins,
                                   zPt_2_Zinc1jetMbins);
    genZPt_2_Zinc1jetM170_250 = newTH1D("genZPt_2_Zinc1jetM170_250",
                                      "gen Z p_{T} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                      ZpT,
                                      nZPt_Zinc1jetMbins,
                                      zPt_2_Zinc1jetMbins);
    hresponseZPt_2_Zinc1jetM170_250 = newTH2D("hresponseZPt_2_Zinc1jetM170_250",
                                            "response Z p_{T} (N_{jets} #geq 1, 130 < M < 170 GeV)",
                                            nZPt_Zinc1jetMbins,
                                            zPt_2_Zinc1jetMbins,
                                            nZPt_Zinc1jetMbins,
                                            zPt_2_Zinc1jetMbins);
*/
//--------------------------


    // M250_3
    ZPt_Zinc1jetM250_3 = newTH1D("ZPt_Zinc1jetM250_3",
                                 "Z p_{T} (N_{jets} #geq 1, 250 < M < 3000 GeV)",
                                 ZpT,
                                 nZPt_Zinc1jetHighMbins,
                                 zPt_Zinc1jetHighMbins);
    genZPt_Zinc1jetM250_3 = newTH1D("genZPt_Zinc1jetM250_3",
                                    "gen Z p_{T} (N_{jets} #geq 1, 250 < M < 3000 GeV)",
                                    ZpT,
                                    nZPt_Zinc1jetHighMbins,
                                    zPt_Zinc1jetHighMbins);
    hresponseZPt_Zinc1jetM250_3 = newTH2D("hresponseZPt_Zinc1jetM250_3",
                                          "response Z p_{T} (N_{jets} #geq 1, 250 < M < 3000 GeV)",
                                          nZPt_Zinc1jetHighMbins,
                                          zPt_Zinc1jetHighMbins,
                                          nZPt_Zinc1jetHighMbins,
                                          zPt_Zinc1jetHighMbins);

    //--------------------------------------------------
    HadRecoil = newTH1D(
        "HadRecoil", "hadRecoil p_{T} (N_{jets} #geq 1)", HRecoi, nZPt_Zinc1jet, zPt_Zinc1jet);
    genHadRecoil = newTH1D("genHadRecoil",
                           "genHadRecoil p_{T} (N_{jets} #geq 1)",
                           HRecoi,
                           nZPt_Zinc1jet,
                           zPt_Zinc1jet);
    JZB =
        newTH1D("JZB", "JZB p_{T} (N_{jets} #geq 1)", JZb, nZPt_Zinc2jetQunJZB, zPt_Zinc2jetQunJZB);
    JZB_Odd = newTH1D(
        "JZB_Odd", "JZB p_{T} (N_{jets} #geq 1)", JZb, nZPt_Zinc2jetQunJZB, zPt_Zinc2jetQunJZB);
    JZB_Even = newTH1D(
        "JZB_Even", "JZB p_{T} (N_{jets} #geq 1)", JZb, nZPt_Zinc2jetQunJZB, zPt_Zinc2jetQunJZB);
    JZB_2 = newTH1D("JZB_2", "JZB p_{T} (N_{jets} #geq 1)", JZb, zPt_2_Zinc2jetQunJZB);
    genJZB = newTH1D(
        "genJZB", "gen JZB p_{T} (N_{jets} #geq 1)", JZb, nZPt_Zinc2jetQunJZB, zPt_Zinc2jetQunJZB);
    JZB_ptLow = newTH1D("JZB_ptLow",
                        "JZB p_{T} (N_{jets} #geq 1)",
                        JZb_low,
                        nZPt_Zinc1jetQunJZB,
                        zPt_Zinc1jetQunJZB);
    JZB_ptLow_Odd = newTH1D("JZB_ptLow_Odd",
                            "JZB p_{T} (N_{jets} #geq 1)",
                            JZb_low,
                            nZPt_Zinc1jetQunJZB,
                            zPt_Zinc1jetQunJZB);
    JZB_ptLow_Even = newTH1D("JZB_ptLow_Even",
                             "JZB p_{T} (N_{jets} #geq 1)",
                             JZb_low,
                             nZPt_Zinc1jetQunJZB,
                             zPt_Zinc1jetQunJZB);
    JZB_ptLow_2 =
        newTH1D("JZB_ptLow_2", "JZB p_{T} (N_{jets} #geq 1)", JZb_low, zPt_2_Zinc1jetQunJZB);
    genJZB_ptLow = newTH1D("genJZB_ptLow",
                           "gen JZB p_{T} (N_{jets} #geq 1)",
                           JZb_low,
                           nZPt_Zinc1jetQunJZB,
                           zPt_Zinc1jetQunJZB);
    JZB_ptHigh = newTH1D("JZB_ptHigh",
                         "JZB p_{T} (N_{jets} #geq 1)",
                         JZb_high,
                         nZPt_Zinc2jetQunJZBptHigh,
                         zPt_Zinc2jetQunJZBptHigh);
    JZB_ptHigh_Odd = newTH1D("JZB_ptHigh_Odd",
                             "JZB p_{T} (N_{jets} #geq 1)",
                             JZb_high,
                             nZPt_Zinc2jetQunJZBptHigh,
                             zPt_Zinc2jetQunJZBptHigh);
    JZB_ptHigh_Even = newTH1D("JZB_ptHigh_Even",
                              "JZB p_{T} (N_{jets} #geq 1)",
                              JZb_high,
                              nZPt_Zinc2jetQunJZBptHigh,
                              zPt_Zinc2jetQunJZBptHigh);
    JZB_ptHigh_2 = newTH1D(
        "JZB_ptHigh_2", "JZB p_{T} (N_{jets} #geq 1)", JZb_high, zPt_2_Zinc2jetQunJZBptHigh);
    genJZB_ptHigh = newTH1D("genJZB_ptHigh",
                            "gen JZB p_{T} (N_{jets} #geq 1)",
                            JZb_high,
                            nZPt_Zinc2jetQunJZBptHigh,
                            zPt_Zinc2jetQunJZBptHigh);

    ZPt_Zinc2jet =
        newTH1D("ZPt_Zinc2jet", "Z p_{T} (N_{jets} #geq 2)", ZpT, nZPt_Zinc2jet, zPt_Zinc2jet);
    genZPt_Zinc2jet = newTH1D(
        "genZPt_Zinc2jet", "gen Z p_{T} (N_{jets} #geq 2)", ZpT, nZPt_Zinc2jet, zPt_Zinc2jet);
    VisPt_Zinc0jetQun = newTH1D("VisPt_Zinc0jetQun",
                                "visible p_{T} (N_{jets} #geq 0)",
                                ZpTVis,
                                nZPt_Zinc1jet,
                                zPt_Zinc1jet);
    VisPt_Zinc0jetQun_Odd = newTH1D("VisPt_Zinc0jetQun_Odd",
                                    "visible p_{T} (N_{jets} #geq 0)",
                                    ZpTVis,
                                    nZPt_Zinc1jet,
                                    zPt_Zinc1jet);
    VisPt_Zinc0jetQun_Even = newTH1D("VisPt_Zinc0jetQun_Even",
                                     "visible p_{T} (N_{jets} #geq 0)",
                                     ZpTVis,
                                     nZPt_Zinc1jet,
                                     zPt_Zinc1jet);
    VisPt_2_Zinc0jetQun =
        newTH1D("VisPt_2_Zinc0jetQun", "visible p_{T} (N_{jets} #geq 0)", ZpTVis, zPt_2_Zinc1jet);
    genVisPt_Zinc0jetQun = newTH1D("genVisPt_Zinc0jetQun",
                                   "gen vis p_{T} (N_{jets} #geq 0)",
                                   ZpTVis,
                                   nZPt_Zinc1jet,
                                   zPt_Zinc1jet);
    VisPt_Zinc1jetQun = newTH1D("VisPt_Zinc1jetQun",
                                "visible p_{T} (N_{jets} #geq 1)",
                                ZpTVis,
                                nZPt_Zinc2jetQun,
                                zPt_Zinc2jetQun);
    VisPt_Zinc1jetQun_Odd = newTH1D("VisPt_Zinc1jetQun_Odd",
                                    "visible p_{T} (N_{jets} #geq 1)",
                                    ZpTVis,
                                    nZPt_Zinc2jetQun,
                                    zPt_Zinc2jetQun);
    VisPt_Zinc1jetQun_Even = newTH1D("VisPt_Zinc1jetQun_Even",
                                     "visible p_{T} (N_{jets} #geq 1)",
                                     ZpTVis,
                                     nZPt_Zinc2jetQun,
                                     zPt_Zinc2jetQun);
    VisPt_2_Zinc1jetQun = newTH1D(
        "VisPt_2_Zinc1jetQun", "visible p_{T} (N_{jets} #geq 1)", ZpTVis, zPt_2_Zinc2jetQun);
    genVisPt_Zinc1jetQun = newTH1D("genVisPt_Zinc1jetQun",
                                   "gen vis p_{T} (N_{jets} #geq 1)",
                                   ZpTVis,
                                   nZPt_Zinc2jetQun,
                                   zPt_Zinc2jetQun);
    VisPt_Zinc2jetQun = newTH1D("VisPt_Zinc2jetQun",
                                "visible p_{T} (N_{jets} #geq 2)",
                                ZpTVis,
                                nZPt_Zinc2JetQun,
                                zPt_Zinc2JetQun);
    VisPt_Zinc2jetQun_Odd = newTH1D("VisPt_Zinc2jetQun_Odd",
                                    "visible p_{T} (N_{jets} #geq 2)",
                                    ZpTVis,
                                    nZPt_Zinc2JetQun,
                                    zPt_Zinc2JetQun);
    VisPt_Zinc2jetQun_Even = newTH1D("VisPt_Zinc2jetQun_Even",
                                     "visible p_{T} (N_{jets} #geq 2)",
                                     ZpTVis,
                                     nZPt_Zinc2JetQun,
                                     zPt_Zinc2JetQun);
    VisPt_2_Zinc2jetQun = newTH1D(
        "VisPt_2_Zinc2jetQun", "visible p_{T} (N_{jets} #geq 2)", ZpTVis, zPt_2_Zinc2JetQun);
    genVisPt_Zinc2jetQun = newTH1D("genVisPt_Zinc2jetQun",
                                   "gen vis p_{T} (N_{jets} #geq 2)",
                                   ZpTVis,
                                   nZPt_Zinc2JetQun,
                                   zPt_Zinc2JetQun);
    VisPt_Zinc3jetQun = newTH1D("VisPt_Zinc3jetQun",
                                "visible p_{T} (N_{jets} #geq 3)",
                                ZpTVis,
                                nZPt_Zinc3jetQun,
                                zPt_Zinc3jetQun);
    VisPt_Zinc3jetQun_Odd = newTH1D("VisPt_Zinc3jetQun_Odd",
                                    "visible p_{T} (N_{jets} #geq 3)",
                                    ZpTVis,
                                    nZPt_Zinc3jetQun,
                                    zPt_Zinc3jetQun);
    VisPt_Zinc3jetQun_Even = newTH1D("VisPt_Zinc3jetQun_Even",
                                     "visible p_{T} (N_{jets} #geq 3)",
                                     ZpTVis,
                                     nZPt_Zinc3jetQun,
                                     zPt_Zinc3jetQun);
    VisPt_2_Zinc3jetQun = newTH1D(
        "VisPt_2_Zinc3jetQun", "visible p_{T} (N_{jets} #geq 3)", ZpTVis, zPt_2_Zinc3jetQun);
    genVisPt_Zinc3jetQun = newTH1D("genVisPt_Zinc3jetQun",
                                   "gen vis p_{T} (N_{jets} #geq 3)",
                                   ZpTVis,
                                   nZPt_Zinc3jetQun,
                                   zPt_Zinc3jetQun);

    ZPt_Zexc0jet = newTH1D("ZPt_Zexc0jet", "Z p_{T} (N_{jets} = 0)", ZpT, 40, 0, 400);
    ZPt_Zexc1jet = newTH1D("ZPt_Zexc1jet", "Z p_{T} (N_{jets} = 1)", ZpT, 40, 0, 400);
    ZPt_Zexc2jet = newTH1D("ZPt_Zexc2jet", "Z p_{T} (N_{jets} = 2)", ZpT, 40, 0, 400);

    ZRapidity_Zinc0jet =
        newTH1D("ZRapidity_Zinc0jet", "Z Rapidity (N_{jets} #geq 0)", Zrap, 30, -3, 3);
    ZRapidity_Zinc1jet =
        newTH1D("ZRapidity_Zinc1jet", "Z Rapidity (N_{jets} #geq 1)", Zrap, 30, -3, 3);
    ZAbsRapidity_Zinc1jet = newTH1D(
        "ZAbsRapidity_Zinc1jet", "Z Absolute Rapidity (N_{jets} #geq 1)", "|y_{Z}|", 12, 0, 2.4);
    ZRapidity_Zinc2jet =
        newTH1D("ZRapidity_Zinc2jet", "Z Rapidity (N_{jets} #geq 2)", Zrap, 30, -3, 3);

    genZRapidity_Zinc0jet =
        newTH1D("genZRapidity_Zinc0jet", "gen Z Rapidity (N_{jets} #geq 0)", Zrap, 30, -3, 3);
    genZRapidity_Zinc1jet =
        newTH1D("genZRapidity_Zinc1jet", "gen Z Rapidity (N_{jets} #geq 1)", Zrap, 30, -3, 3);
    genZAbsRapidity_Zinc1jet = newTH1D("genZAbsRapidity_Zinc1jet",
                                       "gen Z Absolute Rapidity (N_{jets} #geq 1)",
                                       "|y_{Z}|",
                                       12,
                                       0,
                                       2.4);
    genZRapidity_Zinc2jet =
        newTH1D("genZRapidity_Zinc2jet", "gen Z Rapidity (N_{jets} #geq 2)", Zrap, 30, -3, 3);

    ZRapidity_Zexc0jet =
        newTH1D("ZRapidity_Zexc0jet", "Z Rapidity (N_{jets} = 0)", Zrap, 30, -3, 3);
    ZRapidity_Zexc1jet =
        newTH1D("ZRapidity_Zexc1jet", "Z Rapidity (N_{jets} = 1)", Zrap, 30, -3, 3);
    ZRapidity_Zexc2jet =
        newTH1D("ZRapidity_Zexc2jet", "Z Rapidity (N_{jets} = 2)", Zrap, 30, -3, 3);

    ZEta_Zinc0jet = newTH1D("ZEta_Zinc0jet", "Z #eta (N_{jets} #geq 0)", Zeta, 30, -3, 3);
    ZEtaUpTo5_Zinc0jet =
        newTH1D("ZEtaUpTo5_Zinc0jet", "Z #eta (N_{jets} #geq 0)", Zeta, 120, -6, 6);
    ZEta_Zinc1jet = newTH1D("ZEta_Zinc1jet", "Z #eta (N_{jets} #geq 1)", Zeta, 30, -3, 3);
    ZEtaUpTo5_Zinc1jet =
        newTH1D("ZEtaUpTo5_Zinc1jet", "Z #eta (N_{jets} #geq 1)", Zeta, 120, -6, 6);
    ZEta_Zinc2jet = newTH1D("ZEta_Zinc2jet", "Z #eta (N_{jets} #geq 2)", Zeta, 30, -3, 3);

    genZEta_Zinc0jet = newTH1D("genZEta_Zinc0jet", "gen Z #eta (N_{jets} #geq 0)", Zeta, 30, -3, 3);
    genZEta_Zinc1jet = newTH1D("genZEta_Zinc1jet", "gen Z #eta (N_{jets} #geq 1)", Zeta, 30, -3, 3);
    genZEta_Zinc2jet = newTH1D("genZEta_Zinc2jet", "gen Z #eta (N_{jets} #geq 2)", Zeta, 30, -3, 3);

    ZEta_Zexc0jet = newTH1D("ZEta_Zexc0jet", "Z #eta (N_{jets} = 0)", Zeta, 30, -3, 3);
    ZEta_Zexc1jet = newTH1D("ZEta_Zexc1jet", "Z #eta (N_{jets} = 1)", Zeta, 30, -3, 3);
    ZEta_Zexc2jet = newTH1D("ZEta_Zexc2jet", "Z #eta (N_{jets} = 2)", Zeta, 30, -3, 3);

    lepEta_Zinc0jet =
        newTH1D("lepEta_Zinc0jet", "1st & 2nd lep #eta (N_{jets} #geq 0)", leta, 24, -2.4, 2.4);
    lepLEta_Zinc0jet =
        newTH1D("lepLEta_Zinc0jet", "1st lep #eta (N_{jets} #geq 0)", leta, 24, -2.4, 2.4);
    lepSEta_Zinc0jet =
        newTH1D("lepSEta_Zinc0jet", "2st lep #eta (N_{jets} #geq 0)", leta, 24, -2.4, 2.4);
    lepEta_Zinc1jet =
        newTH1D("lepEta_Zinc1jet", "1st & 2nd lep #eta (N_{jets} #geq 1)", leta, 24, -2.4, 2.4);
    lepEtaUpTo4_Zinc0jet =
        newTH1D("lepEtaUpTo4_Zinc0jet", "1st & 2nd lep #eta (N_{jets} #geq 0)", leta, 80, -4, 4);

    lepPhi_Zinc0jet =
        newTH1D("lepPhi_Zinc0jet", "1st & 2nd lep #phi (N_{jets} #geq 0)", lphi, 24, -PI, PI);

    genlepEta_Zinc0jet =
        newTH1D("genlepEta_Zinc0jet", "1st & 2nd lep #eta (N_{jets} #geq 0)", leta, 24, -2.4, 2.4);

    lepEta_Zexc0jet =
        newTH1D("lepEta_Zexc0jet", "1st & 2nd lep #eta (N_{jets} = 0)", leta, 24, -2.4, 2.4);

    lepPhi_Zexc0jet =
        newTH1D("lepPhi_Zexc0jet", "1st & 2nd lep #phi (N_{jets} = 0)", lphi, 24, -PI, PI);

    FirstJetEtaFull_Zinc1jet = newTH1D(
        "FirstJetEtaFull_Zinc1jet", "1st jet #eta (N_{jets} #geq 1)", "#eta(j_{1})", 48, -2.4, 2.4);
    SecondJetEtaFull_Zinc2jet = newTH1D("SecondJetEtaFull_Zinc2jet",
                                        "2nd jet #eta (N_{jets} #geq 2)",
                                        "#eta(j_{2})",
                                        48,
                                        -2.4,
                                        2.4);
    ThirdJetEtaFull_Zinc3jet = newTH1D(
        "ThirdJetEtaFull_Zinc3jet", "3rd jet #eta (N_{jets} #geq 3)", "#eta(j_{3})", 16, -2.4, 2.4);
    FourthJetEtaFull_Zinc4jet = newTH1D(
        "FourthJetEtaFull_Zinc4jet", "4th jet #eta (N_{jets} #geq 4)", "#eta(j_{4})", 8, -2.4, 2.4);
    FifthJetEtaFull_Zinc5jet = newTH1D(
        "FifthJetEtaFull_Zinc5jet", "5th jet #eta (N_{jets} #geq 5)", "#eta(j_{5})", 4, -2.4, 2.4);
    SixthJetEtaFull_Zinc6jet = newTH1D("SixthJetEtaFull_Zinc6jet",
                                       "#geq 6th jets #eta (N_{jets} #geq 6)",
                                       "#eta(j_{6})",
                                       4,
                                       -2.4,
                                       2.4);

    FirstJetEta_Zinc1jet = newTH1D(
        "FirstJetEta_Zinc1jet", "1st jet |#eta| (N_{jets} #geq 1)", "|#eta(j_{1})|", 32, 0., 2.4);
    FirstJetEta_Zinc1jet_res08 = newTH1D("FirstJetEta_Zinc1jet_res08",
                                         "1st jet |#eta| (N_{jets} #geq 1) res08",
                                         "|#eta(j_{1})|",
                                         40,
                                         -0.1,
                                         0.1);
    FirstJetEta_Zinc1jet_res16 = newTH1D("FirstJetEta_Zinc1jet_res16",
                                         "1st jet |#eta| (N_{jets} #geq 1) res16",
                                         "|#eta(j_{1})|",
                                         40,
                                         -0.1,
                                         0.1);
    FirstJetEta_Zinc1jet_res24 = newTH1D("FirstJetEta_Zinc1jet_res24",
                                         "1st jet |#eta| (N_{jets} #geq 1) res24",
                                         "|#eta(j_{1})|",
                                         40,
                                         -0.1,
                                         0.1);

    FirstJetEta_Zinc1jet_res2D = newTH2D("FirstJetEta_Zinc1jet_res2D",
                                         "1st jet p_{T} (N_{jets} #geq 1) res2D}",
                                         32,
                                         0.,
                                         2.4, // x-axis
                                         40,
                                         -0.1,
                                         0.1); // y-axis

    FirstJetEta_2_Zinc1jet = newTH1D("FirstJetEta_2_Zinc1jet",
                                     "1st jet |#eta| (N_{jets} #geq 1)2",
                                     "|#eta(j_{1})|",
                                     160,
                                     0.,
                                     2.4);
    SecondJetEta_Zinc2jet = newTH1D(
        "SecondJetEta_Zinc2jet", "2nd jet |#eta| (N_{jets} #geq 2)", "|#eta(j_{2})|", 24, 0., 2.4);
    SecondJetEta_2_Zinc2jet = newTH1D("SecondJetEta_2_Zinc2jet",
                                      "2nd jet |#eta| (N_{jets} #geq 2)2",
                                      "|#eta(j_{2})|",
                                      120,
                                      0.,
                                      2.4);
    ThirdJetEta_Zinc3jet = newTH1D(
        "ThirdJetEta_Zinc3jet", "3rd jet |#eta| (N_{jets} #geq 3)", "|#eta(j_{3})|", 12, 0., 2.4);
    ThirdJetEta_2_Zinc3jet = newTH1D("ThirdJetEta_2_Zinc3jet",
                                     "3rd jet |#eta| (N_{jets} #geq 3)2",
                                     "|#eta(j_{3})|",
                                     60,
                                     0.,
                                     2.4);
    FourthJetEta_Zinc4jet = newTH1D(
        "FourthJetEta_Zinc4jet", "4th jet |#eta| (N_{jets} #geq 4)", "|#eta(j_{4})|", 12, 0., 2.4);
    FifthJetEta_Zinc5jet = newTH1D(
        "FifthJetEta_Zinc5jet", "5th jet |#eta| (N_{jets} #geq 5)", "|#eta(j_{5})|", 12, 0., 2.4);
    SixthJetEta_Zinc6jet = newTH1D(
        "SixthJetEta_Zinc6jet", "6th jet |#eta| (N_{jets} #geq 6)", "|#eta(j_{6})|", 12, 0., 2.4);

    FirstJetEtaHigh_Zinc1jet = newTH1D("FirstJetEtaHigh_Zinc1jet",
                                       "1st jet |#eta| (N_{jets} #geq 1)",
                                       "|#eta(j_{1})|",
                                       47,
                                       0,
                                       4.7);
    SecondJetEtaHigh_Zinc2jet = newTH1D("SecondJetEtaHigh_Zinc2jet",
                                        "2nd jet |#eta| (N_{jets} #geq 2)",
                                        "|#eta(j_{2})|",
                                        47,
                                        0,
                                        4.7);
    ThirdJetEtaHigh_Zinc3jet = newTH1D("ThirdJetEtaHigh_Zinc3jet",
                                       "3rd jet |#eta| (N_{jets} #geq 3)",
                                       "|#eta(j_{3})|",
                                       24,
                                       0.,
                                       4.7);
    FourthJetEtaHigh_Zinc4jet = newTH1D("FourthJetEtaHigh_Zinc4jet",
                                        "4th jet |#eta| (N_{jets} #geq 4)",
                                        "|#eta(j_{4})|",
                                        12,
                                        0.,
                                        4.7);
    FifthJetEtaHigh_Zinc5jet = newTH1D("FifthJetEtaHigh_Zinc5jet",
                                       "5th jet |#eta| (N_{jets} #geq 5)",
                                       "|#eta(j_{5})|",
                                       6,
                                       0.,
                                       4.7);
    SixthJetEtaHigh_Zinc6jet = newTH1D("SixthJetEtaHigh_Zinc6jet",
                                       "6th jet |#eta| (N_{jets} #geq 6)",
                                       "|#eta(j_{6})|",
                                       6,
                                       0.,
                                       4.7);

    // DJALOG
    FirstJetAbsYvsAbsEta_Zinc1jet = newTH2D(
        "FirstJetAbsYvsAbsEta_Zinc1jet", "FirstJetAbsYvsAbsEta_Zinc1jet", 12, 0, 2.4, 12, 0, 2.4);

    FirstJetAbsRapidity_Zinc1jet = newTH1D(
        "FirstJetAbsRapidity_Zinc1jet", "1st jet |y| (N_{jets} #geq 1)", "|y(j_{1})|", 12, 0, 2.4);
    FirstJetAbsRapidity_SmearMatch_Zinc1jet = newTH1D("FirstJetAbsRapidity_SmearMatch_Zinc1jet",
                                                      "1st jet |y| Smear Match (N_{jets} #geq 1)",
                                                      "|y(j_{1})|",
                                                      12,
                                                      0,
                                                      2.4); // DJALOG
    FirstJetAbsRapidity_SmearGauss_Zinc1jet = newTH1D("FirstJetAbsRapidity_SmearGauss_Zinc1jet",
                                                      "1st jet |y| Smear Gauss (N_{jets} #geq 1)",
                                                      "|y(j_{1})|",
                                                      12,
                                                      0,
                                                      2.4); // DJALOG
    FirstJetAbsRapidity_Zinc1jet_Odd = newTH1D("FirstJetAbsRapidity_Zinc1jet_Odd",
                                               "1st jet |y| (N_{jets} #geq 1)",
                                               "|y(j_{1})|",
                                               12,
                                               0,
                                               2.4);
    FirstJetAbsRapidity_Zinc1jet_Even = newTH1D("FirstJetAbsRapidity_Zinc1jet_Even",
                                                "1st jet |y| (N_{jets} #geq 1)",
                                                "|y(j_{1})|",
                                                12,
                                                0,
                                                2.4);
    FirstJetAbsRapidity_2_Zinc1jet = newTH1D("FirstJetAbsRapidity_2_Zinc1jet",
                                             "1st jet |y| (N_{jets} #geq 1)",
                                             "|y(j_{1})|",
                                             5 * 12,
                                             0,
                                             2.4);

    FirstJetRapidityHigh_Zinc1jet = newTH1D(
        "FirstJetRapidityHigh_Zinc1jet", "1st jet |y| (N_{jets} #geq 1)", "|y(j_{1})|", 47, 0, 4.7);

    SecondJetAbsRapidity_Zinc2jet = newTH1D(
        "SecondJetAbsRapidity_Zinc2jet", "2nd jet |y| (N_{jets} #geq 2)", "|y(j_{2})|", 12, 0, 2.4);
    SecondJetAbsRapidity_Zinc2jet_Odd = newTH1D("SecondJetAbsRapidity_Zinc2jet_Odd",
                                                "2nd jet |y| (N_{jets} #geq 2)",
                                                "|y(j_{2})|",
                                                12,
                                                0,
                                                2.4);
    SecondJetAbsRapidity_Zinc2jet_Even = newTH1D("SecondJetAbsRapidity_Zinc2jet_Even",
                                                 "2nd jet |y| (N_{jets} #geq 2)",
                                                 "|y(j_{2})|",
                                                 12,
                                                 0,
                                                 2.4);
    SecondJetAbsRapidity_2_Zinc2jet = newTH1D("SecondJetAbsRapidity_2_Zinc2jet",
                                              "2nd jet |y| (N_{jets} #geq 2)",
                                              "|y(j_{2})|",
                                              5 * 12,
                                              0,
                                              2.4);

    SecondJetRapidityHigh_Zinc2jet = newTH1D("SecondJetRapidityHigh_Zinc2jet",
                                             "2nd jet |y| (N_{jets} #geq 2)",
                                             "|y(j_{2})|",
                                             47,
                                             0,
                                             4.7);

    ThirdJetAbsRapidity_Zinc3jet = newTH1D(
        "ThirdJetAbsRapidity_Zinc3jet", "3rd jet |y| (N_{jets} #geq 3)", "|y(j_{3})|", 8, 0., 2.4);
    ThirdJetAbsRapidity_Zinc3jet_Odd = newTH1D("ThirdJetAbsRapidity_Zinc3jet_Odd",
                                               "3rd jet |y| (N_{jets} #geq 3)",
                                               "|y(j_{3})|",
                                               8,
                                               0.,
                                               2.4);
    ThirdJetAbsRapidity_Zinc3jet_Even = newTH1D("ThirdJetAbsRapidity_Zinc3jet_Even",
                                                "3rd jet |y| (N_{jets} #geq 3)",
                                                "|y(j_{3})|",
                                                8,
                                                0.,
                                                2.4);
    ThirdJetAbsRapidity_2_Zinc3jet = newTH1D("ThirdJetAbsRapidity_2_Zinc3jet",
                                             "3rd jet |y| (N_{jets} #geq 3)",
                                             "|y(j_{3})|",
                                             5 * 8,
                                             0.,
                                             2.4);

    // to match with binning of current response matrices..
    //    ThirdJetAbsRapidity_Zinc3jet        = newTH1D("ThirdJetAbsRapidity_Zinc3jet",        "3rd
    //    jet |y| (N_{jets} #geq 3)",              "|y(j_{3})|",  12, 0., 2.4);
    //    ThirdJetAbsRapidity_Zinc3jet_Odd    = newTH1D("ThirdJetAbsRapidity_Zinc3jet_Odd",     "3rd
    //    jet |y| (N_{jets} #geq 3)",              "|y(j_{3})|",  12, 0., 2.4);
    //    ThirdJetAbsRapidity_Zinc3jet_Even   = newTH1D("ThirdJetAbsRapidity_Zinc3jet_Even",    "3rd
    //    jet |y| (N_{jets} #geq 3)",              "|y(j_{3})|",  12, 0., 2.4);
    //    ThirdJetAbsRapidity_2_Zinc3jet        = newTH1D("ThirdJetAbsRapidity_2_Zinc3jet",
    //    "3rd jet |y| (N_{jets} #geq 3)",              "|y(j_{3})|",  5*12, 0., 2.4);

    ThirdJetRapidityHigh_Zinc3jet = newTH1D("ThirdJetRapidityHigh_Zinc3jet",
                                            "3rd jet |y| (N_{jets} #geq 3)",
                                            "|y(j_{3})|",
                                            24,
                                            0.,
                                            4.7);
    FourthJetAbsRapidity_Zinc4jet = newTH1D(
        "FourthJetAbsRapidity_Zinc4jet", "4th jet |y| (N_{jets} #geq 4)", "|y(j_{4})|", 8, 0., 2.4);
    FourthJetRapidityHigh_Zinc4jet = newTH1D("FourthJetRapidityHigh_Zinc4jet",
                                             "4th jet |y| (N_{jets} #geq 4)",
                                             "|y(j_{4})|",
                                             12,
                                             0.,
                                             4.7);
    FifthJetAbsRapidity_Zinc5jet = newTH1D(
        "FifthJetAbsRapidity_Zinc5jet", "5th jet |y| (N_{jets} #geq 5)", "|y(j_{5})|", 6, 0., 2.4);
    FifthJetRapidityHigh_Zinc5jet = newTH1D(
        "FifthJetRapidityHigh_Zinc5jet", "5th jet |y| (N_{jets} #geq 5)", "|y(j_{5})|", 6, 0., 4.7);
    SixthJetAbsRapidity_Zinc6jet = newTH1D(
        "SixthJetAbsRapidity_Zinc6jet", "6th jet |y| (N_{jets} #geq 6)", "|y(j_{6})|", 6, 0., 2.4);
    SixthJetRapidityHigh_Zinc6jet = newTH1D(
        "SixthJetRapidityHigh_Zinc6jet", "6th jet |y| (N_{jets} #geq 6)", "|y(j_{6})|", 6, 0., 4.7);

    genFirstJetEta_Zinc1jet = newTH1D("genFirstJetEta_Zinc1jet",
                                      "gen 1st jet #eta (N_{jets} #geq 1)",
                                      "|#eta(j_{1})|",
                                      32,
                                      0.,
                                      2.4);
    genFirstJetEta_2_Zinc1jet = newTH1D("genFirstJetEta_2_Zinc1jet",
                                        "gen 1st jet #eta (N_{jets} #geq 1)2",
                                        "|#eta(j_{1})|",
                                        160,
                                        0.,
                                        2.4);
    genSecondJetEta_Zinc2jet = newTH1D("genSecondJetEta_Zinc2jet",
                                       "gen 2nd jet #eta (N_{jets} #geq 2)",
                                       "|#eta(j_{2})|",
                                       24,
                                       0.,
                                       2.4);
    genSecondJetEta_2_Zinc2jet = newTH1D("genSecondJetEta_2_Zinc2jet",
                                         "gen 2nd jet #eta (N_{jets} #geq 2)2",
                                         "|#eta(j_{2})|",
                                         120,
                                         0.,
                                         2.4);
    genThirdJetEta_Zinc3jet = newTH1D("genThirdJetEta_Zinc3jet",
                                      "gen 3rd jet #eta (N_{jets} #geq 3)",
                                      "|#eta(j_{3})|",
                                      12,
                                      0.,
                                      2.4);
    genThirdJetEta_2_Zinc3jet = newTH1D("genThirdJetEta_2_Zinc3jet",
                                        "gen 3rd jet #eta (N_{jets} #geq 3)2",
                                        "|#eta(j_{3})|",
                                        60,
                                        0.,
                                        2.4);
    genFourthJetEta_Zinc4jet = newTH1D("genFourthJetEta_Zinc4jet",
                                       "gen 4th jet #eta (N_{jets} #geq 4)",
                                       "|#eta(j_{4})|",
                                       8,
                                       0.,
                                       2.4);
    genFifthJetEta_Zinc5jet = newTH1D("genFifthJetEta_Zinc5jet",
                                      "gen 5th jet #eta (N_{jets} #geq 5)",
                                      "|#eta(j_{5})|",
                                      6,
                                      0.,
                                      2.4);
    genSixthJetEta_Zinc6jet = newTH1D("genSixthJetEta_Zinc6jet",
                                      "gen 6th jet #eta (N_{jets} #geq 6)",
                                      "|#eta(j_{6})|",
                                      6,
                                      0.,
                                      2.4);

    genFirstJetEtaHigh_Zinc1jet = newTH1D("genFirstJetEtaHigh_Zinc1jet",
                                          "gen 1st jet #eta (N_{jets} #geq 1)",
                                          "|#eta(j_{1})|",
                                          47,
                                          0,
                                          4.7);
    genSecondJetEtaHigh_Zinc2jet = newTH1D("genSecondJetEtaHigh_Zinc2jet",
                                           "gen 2nd jet #eta (N_{jets} #geq 2)",
                                           "|#eta(j_{2})|",
                                           47,
                                           0,
                                           4.7);
    genThirdJetEtaHigh_Zinc3jet = newTH1D("genThirdJetEtaHigh_Zinc3jet",
                                          "gen 3rd jet #eta (N_{jets} #geq 3)",
                                          "|#eta(j_{3})|",
                                          24,
                                          0.,
                                          4.7);
    genFourthJetEtaHigh_Zinc4jet = newTH1D("genFourthJetEtaHigh_Zinc4jet",
                                           "gen 4th jet #eta (N_{jets} #geq 4)",
                                           "|#eta(j_{4})|",
                                           12,
                                           0.,
                                           4.7);
    genFifthJetEtaHigh_Zinc5jet = newTH1D("genFifthJetEtaHigh_Zinc5jet",
                                          "gen 5th jet #eta (N_{jets} #geq 5)",
                                          "|#eta(j_{5})|",
                                          6,
                                          0.,
                                          4.7);
    genSixthJetEtaHigh_Zinc6jet = newTH1D("genSixthJetEtaHigh_Zinc6jet",
                                          "gen 6th jet #eta (N_{jets} #geq 6)",
                                          "|#eta(j_{6})|",
                                          6,
                                          0.,
                                          4.7);

    genFirstJetAbsRapidity_Zinc1jet = newTH1D("genFirstJetAbsRapidity_Zinc1jet",
                                              "gen 1st jet |y| (N_{jets} #geq 1)",
                                              "|y(j_{1})|",
                                              12,
                                              0,
                                              2.4);
    genFirstJetRapidityHigh_Zinc1jet = newTH1D("genFirstJetRapidityHigh_Zinc1jet",
                                               "gen 1st jet |y| (N_{jets} #geq 1)",
                                               "|y(j_{1})|",
                                               47,
                                               0,
                                               4.7);
    genSecondJetAbsRapidity_Zinc2jet = newTH1D("genSecondJetAbsRapidity_Zinc2jet",
                                               "gen 2nd jet |y| (N_{jets} #geq 2)",
                                               "|y(j_{2})|",
                                               12,
                                               0,
                                               2.4);
    genSecondJetRapidityHigh_Zinc2jet = newTH1D("genSecondJetRapidityHigh_Zinc2jet",
                                                "gen 2nd jet |y| (N_{jets} #geq 2)",
                                                "|y(j_{2})|",
                                                47,
                                                0,
                                                4.7);
    genThirdJetAbsRapidity_Zinc3jet = newTH1D("genThirdJetAbsRapidity_Zinc3jet",
                                              "gen 3rd jet |y| (N_{jets} #geq 3)",
                                              "|y(j_{3})|",
                                              8,
                                              0,
                                              2.4);
    genThirdJetRapidityHigh_Zinc3jet = newTH1D("genThirdJetRapidityHigh_Zinc3jet",
                                               "gen 3rd jet |y| (N_{jets} #geq 3)",
                                               "|y(j_{3})|",
                                               24,
                                               0,
                                               4.7);
    genFourthJetAbsRapidity_Zinc4jet = newTH1D("genFourthJetAbsRapidity_Zinc4jet",
                                               "gen 4th jet |y| (N_{jets} #geq 4)",
                                               "|y(j_{4})|",
                                               8,
                                               0,
                                               2.4);
    genFourthJetRapidityHigh_Zinc4jet = newTH1D("genFourthJetRapidityHigh_Zinc4jet",
                                                "gen 4th jet |y| (N_{jets} #geq 4)",
                                                "|y(j_{4})|",
                                                12,
                                                0,
                                                4.7);
    genFifthJetAbsRapidity_Zinc5jet = newTH1D("genFifthJetAbsRapidity_Zinc5jet",
                                              "gen 5th jet |y| (N_{jets} #geq 5)",
                                              "|y(j_{5})|",
                                              6,
                                              0,
                                              2.4);
    genFifthJetRapidityHigh_Zinc5jet = newTH1D("genFifthJetRapidityHigh_Zinc5jet",
                                               "gen 5th jet |y| (N_{jets} #geq 5)",
                                               "|y(j_{5})|",
                                               6,
                                               0,
                                               4.7);
    genSixthJetAbsRapidity_Zinc6jet = newTH1D("genSixthJetAbsRapidity_Zinc6jet",
                                              "gen 6th jet |y| (N_{jets} #geq 6)",
                                              "|y(j_{6})|",
                                              6,
                                              0,
                                              2.4);
    genSixthJetRapidityHigh_Zinc6jet = newTH1D("genSixthJetRapidityHigh_Zinc6jet",
                                               "gen 6th jet |y| (N_{jets} #geq 6)",
                                               "|y(j_{6})|",
                                               6,
                                               0,
                                               4.7);

    FirstJetEta_Zexc1jet = newTH1D(
        "FirstJetEta_Zexc1jet", "1st jet #eta (N_{jets} = 1)", "#eta(j_{1})", 47, -4.7, 4.7);
    SecondJetEta_Zexc2jet = newTH1D(
        "SecondJetEta_Zexc2jet", "2nd jet #eta (N_{jets} = 2)", "#eta(j_{2})", 47, -4.7, 4.7);

    FirstJetPhi_Zinc1jet = newTH1D(
        "FirstJetPhi_Zinc1jet", "1st jet #phi (N_{jets} #geq 1)", "#phi(j_{1})", 30, -PI, PI);
    SecondJetPhi_Zinc2jet = newTH1D(
        "SecondJetPhi_Zinc2jet", "2nd jet #phi (N_{jets} #geq 2)", "#phi(j_{2})", 30, -PI, PI);
    ThirdJetPhi_Zinc3jet = newTH1D(
        "ThirdJetPhi_Zinc3jet", "3rd jet #phi (N_{jets} #geq 3)", "#phi(j_{3})", 30, -PI, PI);
    FourthJetPhi_Zinc4jet = newTH1D(
        "FourthJetPhi_Zinc4jet", "4th jet #phi (N_{jets} #geq 4)", "#phi(j_{4})", 30, -PI, PI);
    FifthJetPhi_Zinc5jet = newTH1D(
        "FifthJetPhi_Zinc5jet", "5th jet #phi (N_{jets} #geq 5)", "#phi(j_{5})", 30, -PI, PI);
    SixthJetPhi_Zinc6jet = newTH1D(
        "SixthJetPhi_Zinc6jet", "6th jet #phi (N_{jets} #geq 6)", "#phi(j_{6})", 30, -PI, PI);

    FirstJetPhi_Zexc1jet =
        newTH1D("FirstJetPhi_Zexc1jet", "1st jet #phi (N_{jets} = 1)", "#phi(j_{1})", 30, -PI, PI);
    SecondJetPhi_Zexc2jet =
        newTH1D("SecondJetPhi_Zexc2jet", "2nd jet #phi (N_{jets} = 2)", "#phi(j_{2})", 30, -PI, PI);

    lepPt_Zinc0jet =
        newTH1D("lepPt_Zinc0jet", "1st & 2nd lep p_{T} (N_{jets} #geq 0)", lpT, 40, 0, 200);

    hresponselep0Pt_Zinc0jet =
        newTH2D("hresponselep0Pt_Zinc0jet", "hresp lep lead pt", 40, 20, 200, 40, 20, 200);
    lepLPt_Zinc0jet =
        newTH1D("lepLPt_Zinc0jet", "1st lep p_{T} (N_{jets} #geq 0)", lpT, 40, 0, 200);
    lepSPt_Zinc0jet =
        newTH1D("lepSPt_Zinc0jet", "2nd lep p_{T} (N_{jets} #geq 0)", lpT, 40, 0, 200);
    hresponselep1Pt_Zinc0jet =
        newTH2D("hresponselep1Pt_Zinc0jet", "hresp lep sublead pt", 40, 20, 200, 40, 20, 200);

    lepPt_Zinc1jet =
        newTH1D("lepPt_Zinc1jet", "1st & 2nd lep p_{T} (N_{jets} #geq 1)", lpT, 40, 0, 200);
    lepLPt_Zinc1jet =
        newTH1D("lepLPt_Zinc1jet", "1st lep p_{T} (N_{jets} #geq 1)", lpT, 40, 0, 200);
    lepSPt_Zinc1jet =
        newTH1D("lepSPt_Zinc1jet", "2nd lep p_{T} (N_{jets} #geq 1)", lpT, 40, 0, 200);
    lepPtFrom15_Zinc0jet =
        newTH1D("lepPtFrom15_Zinc0jet", "1st & 2nd lep p_{T} (N_{jets} #geq 0)", lpT, 100, 0, 200);
    genlepPt_Zinc0jet =
        newTH1D("genlepPt_Zinc0jet", "gen 1st & 2nd lep p_{T} (N_{jets} #geq 0)", lpT, 40, 0, 200);
    lepPt_Zexc0jet =
        newTH1D("lepPt_Zexc0jet", "1st & 2nd lep p_{T} (N_{jets} = 0)", lpT, 40, 0, 200);
    lepLPt_Zexc0jet = newTH1D("lepLPt_Zexc0jet", "1st lep p_{T} (N_{jets} = 0)", lpT, 40, 0, 200);
    lepSPt_Zexc0jet = newTH1D("lepSPt_Zexc0jet", "2nd lep p_{T} (N_{jets} = 0)", lpT, 40, 0, 200);

    dPhiLeptons_Zexc0jet =
        newTH1D("dPhiLeptons_Zexc0jet", "#Delta #phi btw lep (N_{jets} = 0)", ldPhi, 50, 0, PI);

    dPhiLeptons_Zinc0jet =
        newTH1D("dPhiLeptons_Zinc0jet", "#Delta #phi btw lep (N_{jets} #geq 0)", ldPhi, 50, 0, PI);
    dPhiLeptons_Zinc1jet =
        newTH1D("dPhiLeptons_Zinc1jet", "#Delta #phi btw lep (N_{jets} #geq 1)", ldPhi, 50, 0, PI);

    dEtaLeptons_Zexc0jet =
        newTH1D("dEtaLeptons_Zexc0jet", "#Delta #eta btw lep (N_{jets} = 0)", ldEta, 50, -5, 5);

    dEtaLeptons_Zinc0jet =
        newTH1D("dEtaLeptons_Zinc0jet", "#Delta #eta btw lep (N_{jets} #geq 0)", ldEta, 50, -5, 5);
    dEtaLeptons_Zinc1jet =
        newTH1D("dEtaLeptons_Zinc1jet", "#Delta #eta btw lep (N_{jets} #geq 1)", ldEta, 50, -5, 5);

    dRLeptons_Zinc0jet =
        newTH1D("dRLeptons_Zinc0jet", "#Delta R btw lep (N_{jets} #geq 0)", ldR, 50, 0, 5);

    dRLeptons_Zinc1jet =
        newTH1D("dRLeptons_Zinc1jet", "#Delta R btw lep (N_{jets} #geq 1)", ldR, 50, 0, 5);

    SpTLeptons_Zexc0jet =
        newTH1D("SpTLeptons_Zexc0jet", "#Delta_{pT}^{rel} lep (N_{jets} = 0)", lSpt, 50, 0, 1);
    SpTLeptons_Zexc1jet =
        newTH1D("SpTLeptons_Zexc1jet", "#Delta_{pT}^{rel} lep (N_{jets} = 1)", lSpt, 50, 0, 1);
    SpTLeptons_Zexc2jet =
        newTH1D("SpTLeptons_Zexc2jet", "#Delta_{pT}^{rel} lep (N_{jets} = 2)", lSpt, 50, 0, 1);

    genSpTLeptons_Zexc2jet = newTH1D(
        "genSpTLeptons_Zexc2jet", "gen #Delta_{pT}^{rel} lep (N_{jets} = 2)", lSpt, 50, 0., 1.);

    SpTLeptons_Zinc0jet =
        newTH1D("SpTLeptons_Zinc0jet", "#Delta_{pT}^{rel} lep (N_{jets} #geq 0)", lSpt, 50, 0, 1);
    SpTLeptons_Zinc1jet =
        newTH1D("SpTLeptons_Zinc1jet", "#Delta_{pT}^{rel} lep (N_{jets} #geq 1)", lSpt, 50, 0, 1);
    SpTLeptons_Zinc2jet =
        newTH1D("SpTLeptons_Zinc2jet", "#Delta_{pT}^{rel} lep (N_{jets} #geq 2)", lSpt, 50, 0, 1);

    genSpTLeptons_Zinc2jet = newTH1D(
        "genSpTLeptons_Zinc2jet", "gen #Delta_{pT}^{rel} lep (N_{jets} #geq 2)", lSpt, 50, 0, 1);

    JetsHT_Zinc1jet = newTH1D("JetsHT_Zinc1jet",
                              "Scalar sum jets p_{T} (N_{jets} #geq 1)",
                              HT,
                              nJetHT_Zinc1jet,
                              jetHT_Zinc1jet);
    JetsHT_Zinc1jet_Odd = newTH1D("JetsHT_Zinc1jet_Odd",
                                  "Scalar sum jets p_{T} (N_{jets} #geq 1)",
                                  HT,
                                  nJetHT_Zinc1jet,
                                  jetHT_Zinc1jet);
    JetsHT_Zinc1jet_Even = newTH1D("JetsHT_Zinc1jet_Even",
                                   "Scalar sum jets p_{T} (N_{jets} #geq 1)",
                                   HT,
                                   nJetHT_Zinc1jet,
                                   jetHT_Zinc1jet);

    JetsHT_2_Zinc1jet = newTH1D(
        "JetsHT_2_Zinc1jet", "Scalar sum jets p_{T} (N_{jets} #geq 1)2", HT, jetHT_2_Zinc1jet);

    JetsHT_Zinc2jet = newTH1D("JetsHT_Zinc2jet",
                              "Scalar sum jets p_{T} (N_{jets} #geq 2)",
                              HT,
                              nJetHT_Zinc2jet,
                              jetHT_Zinc2jet);
    JetsHT_Zinc2jet_Odd = newTH1D("JetsHT_Zinc2jet_Odd",
                                  "Scalar sum jets p_{T} (N_{jets} #geq 2)",
                                  HT,
                                  nJetHT_Zinc2jet,
                                  jetHT_Zinc2jet);
    JetsHT_Zinc2jet_Even = newTH1D("JetsHT_Zinc2jet_Even",
                                   "Scalar sum jets p_{T} (N_{jets} #geq 2)",
                                   HT,
                                   nJetHT_Zinc2jet,
                                   jetHT_Zinc2jet);

    JetsHT_2_Zinc2jet = newTH1D(
        "JetsHT_2_Zinc2jet", "Scalar sum jets p_{T} (N_{jets} #geq 2)2", HT, jetHT_2_Zinc2jet);

    JetsHT_Zinc3jet = newTH1D("JetsHT_Zinc3jet",
                              "Scalar sum jets p_{T} (N_{jets} #geq 3)",
                              HT,
                              nJetHT_Zinc3jet,
                              jetHT_Zinc3jet);
    JetsHT_Zinc3jet_Odd = newTH1D("JetsHT_Zinc3jet_Odd",
                                  "Scalar sum jets p_{T} (N_{jets} #geq 3)",
                                  HT,
                                  nJetHT_Zinc3jet,
                                  jetHT_Zinc3jet);
    JetsHT_Zinc3jet_Even = newTH1D("JetsHT_Zinc3jet_Even",
                                   "Scalar sum jets p_{T} (N_{jets} #geq 3)",
                                   HT,
                                   nJetHT_Zinc3jet,
                                   jetHT_Zinc3jet);

    JetsHT_2_Zinc3jet = newTH1D(
        "JetsHT_2_Zinc3jet", "Scalar sum jets p_{T} (N_{jets} #geq 3)2", HT, jetHT_2_Zinc3jet);

    JetsHT_Zinc4jet = newTH1D("JetsHT_Zinc4jet",
                              "Scalar sum jets p_{T} (N_{jets} #geq 4)",
                              HT,
                              nJetHT_Zinc4jet,
                              jetHT_Zinc4jet);
    JetsHT_Zinc5jet = newTH1D("JetsHT_Zinc5jet",
                              "Scalar sum jets p_{T} (N_{jets} #geq 5)",
                              HT,
                              nJetHT_Zinc5jet,
                              jetHT_Zinc5jet);
    JetsHT_Zinc6jet = newTH1D("JetsHT_Zinc6jet",
                              "Scalar sum jets p_{T} (N_{jets} #geq 6)",
                              HT,
                              nJetHT_Zinc5jet,
                              jetHT_Zinc5jet);

    genJetsHT_Zinc1jet = newTH1D("genJetsHT_Zinc1jet",
                                 "gen Scalar sum jets p_{T} (N_{jets} #geq 1)",
                                 HT,
                                 nJetHT_Zinc1jet,
                                 jetHT_Zinc1jet);

    genJetsHT_2_Zinc1jet = newTH1D("genJetsHT_2_Zinc1jet",
                                   "gen Scalar sum jets p_{T} (N_{jets} #geq 1)2",
                                   HT,
                                   jetHT_2_Zinc1jet);

    genJetsHT_Zinc2jet = newTH1D("genJetsHT_Zinc2jet",
                                 "gen Scalar sum jets p_{T} (N_{jets} #geq 2)",
                                 HT,
                                 nJetHT_Zinc2jet,
                                 jetHT_Zinc2jet);

    genJetsHT_2_Zinc2jet = newTH1D("genJetsHT_2_Zinc2jet",
                                   "gen Scalar sum jets p_{T} (N_{jets} #geq 2)2",
                                   HT,
                                   jetHT_2_Zinc2jet);

    genJetsHT_Zinc3jet = newTH1D("genJetsHT_Zinc3jet",
                                 "gen Scalar sum jets p_{T} (N_{jets} #geq 3)",
                                 HT,
                                 nJetHT_Zinc3jet,
                                 jetHT_Zinc3jet);

    genJetsHT_2_Zinc3jet = newTH1D("genJetsHT_2_Zinc3jet",
                                   "gen Scalar sum jets p_{T} (N_{jets} #geq 3)2",
                                   HT,
                                   jetHT_2_Zinc3jet);

    genJetsHT_Zinc4jet = newTH1D("genJetsHT_Zinc4jet",
                                 "gen Scalar sum jets p_{T} (N_{jets} #geq 4)",
                                 HT,
                                 nJetHT_Zinc4jet,
                                 jetHT_Zinc4jet);
    genJetsHT_Zinc5jet = newTH1D("genJetsHT_Zinc5jet",
                                 "gen Scalar sum jets p_{T} (N_{jets} #geq 5)",
                                 HT,
                                 nJetHT_Zinc5jet,
                                 jetHT_Zinc5jet);
    genJetsHT_Zinc6jet = newTH1D("genJetsHT_Zinc6jet",
                                 "gen Scalar sum jets p_{T} (N_{jets} #geq 6)",
                                 HT,
                                 nJetHT_Zinc5jet,
                                 jetHT_Zinc5jet);

    FirstJetPt_Zinc1jet_NVtx = newTH2D(
        "FirstJetPt_Zinc1jet_NVtx",
        "1st jet p_{T} (N_{jets} #geq 1) and vertex multiplicity;p_{T}(j_{1}) [GeV];N_{vertices}",
        nJetPt_Zinc1jet,
        jetPt_Zinc1jet, // x-axis
        45,
        0.5,
        45.5); // y-axis

    FirstJetPtRecoOvGen_Zinc1jet_NVtx =
        newTH2D("FirstJetPtRecoOvGen_Zinc1jet_NVtx",
                "1st jet p_{T}^{reco}/p_{T}^{gen} (N_{jets} #geq 1) and vertex "
                "multiplicity;p_{T}(j_{1})^{reco}/p_{T}(j_{1})^{gen};N_{vertices}",
                100,
                0,
                2, // x-axis
                45,
                0.5,
                45.5); // y-axis

    FirstJetPt_Zinc1jet = newTH1D("FirstJetPt_Zinc1jet",
                                  "1st jet p_{T} (N_{jets} #geq 1)",
                                  "p_{T}(j_{1}) [GeV]",
                                  nJetPt_Zinc1jet,
                                  jetPt_Zinc1jet);
    FirstJetPt_SmearMatch_Zinc1jet = newTH1D("FirstJetPt_SmearMatch_Zinc1jet",
                                             "1st jet p_{T} SmearMatch (N_{jets} #geq 1)",
                                             "p_{T}(j_{1}) [GeV]",
                                             nJetPt_Zinc1jet,
                                             jetPt_Zinc1jet); // DJALOG
    FirstJetPt_SmearGauss_Zinc1jet = newTH1D("FirstJetPt_SmearGauss_Zinc1jet",
                                             "1st jet p_{T} SmearGauss (N_{jets} #geq 1)",
                                             "p_{T}(j_{1}) [GeV]",
                                             nJetPt_Zinc1jet,
                                             jetPt_Zinc1jet); // DJALOG
    FirstJetPt_Zinc1jet_Odd = newTH1D("FirstJetPt_Zinc1jet_Odd",
                                      "1st jet p_{T} (N_{jets} #geq 1)",
                                      "p_{T}(j_{1}) [GeV]",
                                      nJetPt_Zinc1jet,
                                      jetPt_Zinc1jet);
    FirstJetPt_Zinc1jet_Even = newTH1D("FirstJetPt_Zinc1jet_Even",
                                       "1st jet p_{T} (N_{jets} #geq 1)",
                                       "p_{T}(j_{1}) [GeV]",
                                       nJetPt_Zinc1jet,
                                       jetPt_Zinc1jet);
    // FirstJetPt_2_Zinc1jet               = newTH1D("FirstJetPt_2_Zinc1jet",                 "1st
    // jet p_{T} (N_{jets} #geq 1)",             "p_{T}(j_{1}) [GeV]",     nJetPt_Zinc1jet,
    // jetPt_Zinc1jet);
    FirstJetPt_2_Zinc1jet = newTH1D("FirstJetPt_2_Zinc1jet",
                                    "1st jet p_{T} (N_{jets} #geq 1)2",
                                    "p_{T}(j_{1}) [GeV]",
                                    jetPt_2_Zinc1jet);

    SecondJetPt_Zinc2jet = newTH1D("SecondJetPt_Zinc2jet",
                                   "2nd jet p_{T} (N_{jets} #geq 2)",
                                   "p_{T}(j_{2}) [GeV]",
                                   nJetPt_Zinc2jet,
                                   jetPt_Zinc2jet);
    SecondJetPt_Zinc2jet_Odd = newTH1D("SecondJetPt_Zinc2jet_Odd",
                                       "2nd jet p_{T} (N_{jets} #geq 2)",
                                       "p_{T}(j_{2}) [GeV]",
                                       nJetPt_Zinc2jet,
                                       jetPt_Zinc2jet);
    SecondJetPt_Zinc2jet_Even = newTH1D("SecondJetPt_Zinc2jet_Even",
                                        "2nd jet p_{T} (N_{jets} #geq 2)",
                                        "p_{T}(j_{2}) [GeV]",
                                        nJetPt_Zinc2jet,
                                        jetPt_Zinc2jet);

    SecondJetPt_2_Zinc2jet = newTH1D("SecondJetPt_2_Zinc2jet",
                                     "2nd jet p_{T} (N_{jets} #geq 2)2",
                                     "p_{T}(j_{2}) [GeV]",
                                     jetPt_2_Zinc2jet);

    ThirdJetPt_Zinc3jet = newTH1D("ThirdJetPt_Zinc3jet",
                                  "3rd jet p_{T} (N_{jets} #geq 3)",
                                  "p_{T}(j_{3}) [GeV]",
                                  nJetPt_Zinc3jet,
                                  jetPt_Zinc3jet);
    ThirdJetPt_Zinc3jet_Odd = newTH1D("ThirdJetPt_Zinc3jet_Odd",
                                      "3rd jet p_{T} (N_{jets} #geq 3)",
                                      "p_{T}(j_{3}) [GeV]",
                                      nJetPt_Zinc3jet,
                                      jetPt_Zinc3jet);
    ThirdJetPt_Zinc3jet_Even = newTH1D("ThirdJetPt_Zinc3jet_Even",
                                       "3rd jet p_{T} (N_{jets} #geq 3)",
                                       "p_{T}(j_{3}) [GeV]",
                                       nJetPt_Zinc3jet,
                                       jetPt_Zinc3jet);

    ThirdJetPt_2_Zinc3jet = newTH1D("ThirdJetPt_2_Zinc3jet",
                                    "3rd jet p_{T} (N_{jets} #geq 3)2",
                                    "p_{T}(j_{3}) [GeV]",
                                    jetPt_2_Zinc3jet);
    FourthJetPt_Zinc4jet = newTH1D("FourthJetPt_Zinc4jet",
                                   "4th jet p_{T} (N_{jets} #geq 4)",
                                   "p_{T}(j_{4}) [GeV]",
                                   nJetPt_Zinc4jet,
                                   jetPt_Zinc4jet);
    FifthJetPt_Zinc5jet = newTH1D("FifthJetPt_Zinc5jet",
                                  "5th jet p_{T} (N_{jets} #geq 5)",
                                  "p_{T}(j_{5}) [GeV]",
                                  nJetPt_Zinc5jet,
                                  jetPt_Zinc5jet);
    SixthJetPt_Zinc6jet = newTH1D("SixthJetPt_Zinc6jet",
                                  "6th jet p_{T} (N_{jets} #geq 6)",
                                  "p_{T}(j_{6}) [GeV]",
                                  nJetPt_Zinc5jet,
                                  jetPt_Zinc5jet);

    FirstJetPtEta_Zinc1jet =
        newTH2D("FirstJetPtEta_Zinc1jet", "FirstJetPtEta_Zinc1jet", 10, 20, 420, 6, 0, 2.4);
    genFirstJetPtEta_Zinc1jet =
        newTH2D("genFirstJetPtEta_Zinc1jet", "genFirstJetPtEta_Zinc1jet", 10, 20, 420, 6, 0, 2.4);

    genFirstJetPt_Zinc1jet = newTH1D("genFirstJetPt_Zinc1jet",
                                     "gen 1st jet p_{T} (N_{jets} #geq 1)",
                                     "p_{T}(j_{1}) [GeV]",
                                     nJetPt_Zinc1jet,
                                     jetPt_Zinc1jet);
    genFirstJetPt_2_Zinc1jet = newTH1D("genFirstJetPt_2_Zinc1jet",
                                       "gen 1st jet p_{T} (N_{jets} #geq 1) 2",
                                       "p_{T}(j_{1}) [GeV]",
                                       jetPt_2_Zinc1jet);
    genSecondJetPt_Zinc2jet = newTH1D("genSecondJetPt_Zinc2jet",
                                      "gen 2nd jet p_{T} (N_{jets} #geq 2)",
                                      "p_{T}(j_{2}) [GeV]",
                                      nJetPt_Zinc2jet,
                                      jetPt_Zinc2jet);
    genSecondJetPt_2_Zinc2jet = newTH1D("genSecondJetPt_2_Zinc2jet",
                                        "gen 2nd jet p_{T} (N_{jets} #geq 2) 2",
                                        "p_{T}(j_{2}) [GeV]",
                                        jetPt_2_Zinc2jet);
    genThirdJetPt_Zinc3jet = newTH1D("genThirdJetPt_Zinc3jet",
                                     "gen 3rd jet p_{T} (N_{jets} #geq 3)",
                                     "p_{T}(j_{3}) [GeV]",
                                     nJetPt_Zinc3jet,
                                     jetPt_Zinc3jet);
    genThirdJetPt_2_Zinc3jet = newTH1D("genThirdJetPt_2_Zinc3jet",
                                       "gen 3rd jet p_{T} (N_{jets} #geq 3)2",
                                       "p_{T}(j_{3}) [GeV]",
                                       jetPt_2_Zinc3jet);

    genFourthJetPt_Zinc4jet = newTH1D("genFourthJetPt_Zinc4jet",
                                      "gen 4th jet p_{T} (N_{jets} #geq 4)",
                                      "p_{T}(j_{4}) [GeV]",
                                      nJetPt_Zinc4jet,
                                      jetPt_Zinc4jet);
    genFifthJetPt_Zinc5jet = newTH1D("genFifthJetPt_Zinc5jet",
                                     "gen 5th jet p_{T} (N_{jets} #geq 5)",
                                     "p_{T}(j_{5}) [GeV]",
                                     nJetPt_Zinc5jet,
                                     jetPt_Zinc5jet);
    genSixthJetPt_Zinc6jet = newTH1D("genSixthJetPt_Zinc6jet",
                                     "gen 6th jet p_{T} (N_{jets} #geq 6)",
                                     "p_{T}(j_{6}) [GeV]",
                                     nJetPt_Zinc5jet,
                                     jetPt_Zinc5jet);

    FirstJetPt_Zexc1jet = newTH1D("FirstJetPt_Zexc1jet",
                                  "1st jet p_{T} (N_{jets} = 1)",
                                  "p_{T}(j_{1}) [GeV]",
                                  nJetPt_Zinc1jet,
                                  jetPt_Zinc1jet);
    SecondJetPt_Zexc2jet = newTH1D("SecondJetPt_Zexc2jet",
                                   "2nd jet p_{T} (N_{jets} = 2)",
                                   "p_{T}(j_{2}) [GeV]",
                                   nJetPt_Zinc2jet,
                                   jetPt_Zinc2jet);

    genFirstJetPt_Zexc1jet = newTH1D("genFirstJetPt_Zexc1jet",
                                     "gen 1st jet p_{T} (N_{jets} = 1)",
                                     "p_{T}(j_{1}) [GeV]",
                                     nJetPt_Zinc1jet,
                                     jetPt_Zinc1jet);
    genSecondJetPt_Zexc2jet = newTH1D("genSecondJetPt_Zexc2jet",
                                      "gen 2nd jet p_{T} (N_{jets} = 2)",
                                      "p_{T}(j_{2}) [GeV]",
                                      nJetPt_Zinc2jet,
                                      jetPt_Zinc2jet);

    //    ZNGoodJets_Zexc = newTH1D("ZNGoodJets_Zexc","Jet Multiplicity (excl.)", "N_{jets}", 8,
    //    -0.5, 7.5);
    ZNGoodJets_Zexc =
        newTH1D("ZNGoodJets_Zexc", "Jet Multiplicity (excl.)", "N_{jets}", 7, -0.5, 6.5);
    if (ZNGoodJets_Zexc) {
        for (int ibin = 1; ibin < ZNGoodJets_Zexc->GetNbinsX(); ++ibin) {
            ZNGoodJets_Zexc->GetXaxis()->SetBinLabel(ibin, TString::Format("= %d", ibin - 1));
        }

        //	ZNGoodJets_Zexc->GetXaxis()->SetBinLabel(1, "= 0");
        //	ZNGoodJets_Zexc->GetXaxis()->SetBinLabel(2, "= 1");
        //	ZNGoodJets_Zexc->GetXaxis()->SetBinLabel(3, "= 2");
        //	ZNGoodJets_Zexc->GetXaxis()->SetBinLabel(4, "= 3");
        //	ZNGoodJets_Zexc->GetXaxis()->SetBinLabel(5, "= 4");
        //	ZNGoodJets_Zexc->GetXaxis()->SetBinLabel(6, "= 5");
        //	ZNGoodJets_Zexc->GetXaxis()->SetBinLabel(7, "= 6");
        //     //ZNGoodJets_Zexc->GetXaxis()->SetBinLabel(8, "= 7");
        ZNGoodJets_Zexc_Odd = (TH1D *)ZNGoodJets_Zexc->Clone("ZNGoodJets_Zexc_Odd");
        ZNGoodJets_Zexc_Odd->SetName("ZNGoodJets_Zexc_Odd");
        ZNGoodJets_Zexc_Even = (TH1D *)ZNGoodJets_Zexc->Clone("ZNGoodJets_Zexc_Even");
        ZNGoodJets_Zexc_Even->SetName("ZNGoodJets_Zexc_Even");
    }

    SumZJetRapidity_Zinc1jet =
        newTH1D("SumZJetRapidity_Zinc1jet", "SumZJetRapidity_Zinc1jet", "y_{sum}", 12, 0, 2.4);
    genSumZJetRapidity_Zinc1jet = newTH1D(
        "genSumZJetRapidity_Zinc1jet", "genSumZJetRapidity_Zinc1jet", "y_{sum}", 12, 0, 2.4);
    DifZJetRapidity_Zinc1jet =
        newTH1D("DifZJetRapidity_Zinc1jet", "DifZJetRapidity_Zinc1jet", "y_{dif}", 12, 0, 2.4);
    genDifZJetRapidity_Zinc1jet = newTH1D(
        "genDifZJetRapidity_Zinc1jet", "genDifZJetRapidity_Zinc1jet", "y_{dif}", 12, 0, 2.4);
    hresponseSumZJetRapidity_Zinc1jet = newTH2D("hresponseSumZJetRapidity_Zinc1jet",
                                                "hresponseSumZJetRapidity_Zinc1jet",
                                                12,
                                                0.,
                                                2.4,
                                                12,
                                                0.,
                                                2.4);
    hresponseDifZJetRapidity_Zinc1jet = newTH2D("hresponseDifZJetRapidity_Zinc1jet",
                                                "hresponseDifZJetRapidity_Zinc1jet",
                                                12,
                                                0.,
                                                2.4,
                                                12,
                                                0.,
                                                2.4);

    CentralJetPt_Zinc2jet =
        newTH1D("CentralJetPt_Zinc2jet", "Central Jet Pt Zinc2jet", "P_{T}", 40, 30, 830);
    ForwardJetPt_Zinc2jet =
        newTH1D("ForwardJetPt_Zinc2jet", "Forward Jet Pt Zinc2jet", "P_{T}", 40, 30, 830);
    genCentralJetPt_Zinc2jet =
        newTH1D("genCentralJetPt_Zinc2jet", "genCentral Jet Pt Zinc2jet", "P_{T}", 40, 30, 830);
    genForwardJetPt_Zinc2jet =
        newTH1D("genForwardJetPt_Zinc2jet", "genForward Jet Pt Zinc2jet", "P_{T}", 40, 30, 830);
    CentralJetEta_Zinc2jet =
        newTH1D("CentralJetEta_Zinc2jet", "Central Jet Eta Zinc2jet", "Eta", 12, 0, 2.4);
    ForwardJetEta_Zinc2jet =
        newTH1D("ForwardJetEta_Zinc2jet", "Forward Jet Eta Zinc2jet", "Eta", 12, 0, 2.4);
    genCentralJetEta_Zinc2jet =
        newTH1D("genCentralJetEta_Zinc2jet", "genCentral Jet Eta Zinc2jet", "Eta", 12, 0, 2.4);
    genForwardJetEta_Zinc2jet =
        newTH1D("genForwardJetEta_Zinc2jet", "genForward Jet Eta Zinc2jet", "Eta", 12, 0, 2.4);

    // DJALOG
    hjetResolution = newTH1D("JetResolution", "Jet Resolution", "Resolution", 100, 0.0, 1.0);
    hEvtFastJetRho = newTH1D("hEvtFastJetRho", "hEvtFastJetRho", "Rho", 100, 0.0, 100.0);
    hhjetResolutionPt =
        newTH2D("JetResolutionPt", "Jet Resolution vs Pt", 10, 0.0, 200.0, 30, 0.05, 0.35);
    hhjetResolutionEta =
        newTH2D("JetResolutionEta", "Jet Resolution vs Eta", 12, 0.0, 2.4, 30, 0.05, 0.35);
    hjetScaleFactor = newTH1D("JetScaleFactor", "Jet Scale Factor", "Scale Factor", 50, 1.0, 1.5);
    hhJetMatching_StatMatches =
        newTH2D("hhJetMatching_StatMatches", "Number of Matches", 7, 0.5, 7.5, 6, -0.5, 5.5);
    hhJetMatching_StatFail =
        newTH2D("hhJetMatching_StatFail", "Reason for No Match", 7, 0.5, 7.5, 3, -0.5, 2.5);
    hhJetMatching_StatMatches_Lep = newTH2D("hhJetMatching_StatMatches_Lep",
                                            "Number of Matches with Lepton",
                                            7,
                                            0.5,
                                            7.5,
                                            6,
                                            -0.5,
                                            5.5);
    hhJetMatching_StatFail_Lep = newTH2D(
        "hhJetMatching_StatFail_Lep", "Reason for No Match with Lepton", 7, 0.5, 7.5, 3, -0.5, 2.5);

    hhRecoGenJetMatchingMatrix = newTH2D(
        "hhRecoGenJetMatchingMatrix", "Reco Gen Jet dR Matching Stat", 7, 0.5, 7.5, 7, 0.5, 7.5);
    hhRecoGenJetMatchingMatrix_Lep = newTH2D("hhRecoGenJetMatchingMatrixLep",
                                             "Reco Gen Jet dR Matching Stat with Lepton",
                                             7,
                                             0.5,
                                             7.5,
                                             7,
                                             0.5,
                                             7.5);

    hhJetMatching_FirstJetRap_Match = newTH2D("hhJetMatching_FirstJetRap_Match",
                                              "Reco + Gen Jet Rapidity with match",
                                              12,
                                              0,
                                              2.4,
                                              12,
                                              0,
                                              2.4);
    hhJetMatching_FirstJetRap_NoMatch = newTH2D("hhJetMatching_FirstJetRap_NoMatch",
                                                "Reco + Gen Jet Rapidity with no match",
                                                12,
                                                0,
                                                2.4,
                                                12,
                                                0,
                                                2.4);
    hJetMatching_FirstJetdPT_Match = newTH1D("hJetMatching_dPT_FirstMatch",
                                             "dPT first reco jet and first match",
                                             "GeV",
                                             20,
                                             -200.0,
                                             200.0);
    hJetMatching_FirstJetdPT_Other = newTH1D("hJetMatching_dPT_FirstMatch",
                                             "dPT first reco jet and first match",
                                             "GeV",
                                             20,
                                             -200.0,
                                             200.0);
    hJetMatching_dPT_SecondMatch = newTH1D("hJetMatching_dPT_SecondMatch",
                                           "dPT first reco jet and second match",
                                           "GeV",
                                           20,
                                           -200.0,
                                           200.0);
    hJetMatching_dRFirstJet =
        newTH1D("hJetMatching_dRFirstJet", "hJetMatching_dRFirstJet", "dR", 200, 0.0, 6.0);
    hJetMathing_FailedRecoJet_AllEta = newTH1D(
        "hJetMathing_FailedRecoJet_AllEta", "hJetMathing_FailedRecoJet_AllEta", "eta", 12, 0, 2.4);
    hJetMathing_FailedRecoJet_FirstEta = newTH1D("hJetMathing_FailedRecoJet_FirstEta",
                                                 "hJetMathing_FailedRecoJet_FirstEta",
                                                 "eta",
                                                 12,
                                                 0,
                                                 2.4);
    hJetMatching_dRFailedFirstLeadLep = newTH1D("hJetMatching_dRFailedFirstLeadLep",
                                                "hJetMatching_dRFailedFirstLeadLep",
                                                "dR",
                                                200,
                                                0,
                                                6.0);
    hJetMatching_dRFailedFirstSubLep = newTH1D(
        "hJetMatching_dRFailedFirstSubLep", "hJetMatching_dRFailedFirstSubLep", "dR", 200, 0, 6.0);

    hJetMatching_FirstJetPt_FaileddR = newTH1D("hJetMatching_FirstJetPt_FaileddR",
                                               "hJetMatching_FirstJetPt_FaileddR",
                                               "GeV",
                                               28,
                                               20.0,
                                               300.0);
    hJetMatching_FirstJetAbsRapidity_FaileddR = newTH1D("hJetMatching_FirstJetAbsRapidity_FaileddR",
                                                        "hJetMatching_FirstJetAbsRapidity_FaileddR",
                                                        "dR",
                                                        12,
                                                        0,
                                                        2.4);
    hJetMatching_FirstJetPt_FailednGen = newTH1D("hJetMatching_FirstJetPt_FailednGen",
                                                 "hJetMatching_FirstJetPt_FailednGen",
                                                 "GeV",
                                                 28,
                                                 20.0,
                                                 300.0);
    hJetMatching_FirstJetAbsRapidity_FailednGen =
        newTH1D("hJetMatching_FirstJetAbsRapidity_FailednGen",
                "hJetMatching_FirstJetAbsRapidity_FailednGen",
                "dR",
                12,
                0,
                2.4);

    hhRecoGenLeptonMatching =
        newTH2D("hhRecoGenLeptonMatching", "Reco Gen Lepton Matching", 4, 0.0, 4.0, 4, 0.0, 4.0);

    hresponseZNGoodJets_Zexc =
        newTH2D("hresponseZNGoodJets_Zexc", "hresp ZNGoodJets_Zexc", 7, -0.5, 6.5, 7, -0.5, 6.5);

    // hresponseZPt_Zinc1jet =
    // newTH2D("hresponseZPt_Zinc1jet","hresponseZPt_Zinc1jet",nZPt_Zinc1jet, zPt_Zinc1jet,
    // nZPt_Zinc1jet, zPt_Zinc1jet);
    hresponseHadRecoil = newTH2D("hresponseHadRecoil",
                                 "hresponseHadRecoil",
                                 nZPt_Zinc1jet,
                                 zPt_Zinc1jet,
                                 nZPt_Zinc1jet,
                                 zPt_Zinc1jet);
    hresponseJZB = newTH2D("hresponseJZB",
                           "hresp Scalar JZB  p_{T} (N_{jets} #geq 1)",
                           nZPt_Zinc2jetQunJZB,
                           zPt_Zinc2jetQunJZB,
                           nZPt_Zinc2jetQunJZB,
                           zPt_Zinc2jetQunJZB);
    hresponseJZB_ptLow = newTH2D("hresponseJZB_ptLow",
                                 "hresp Scalar JZB  p_{T} (N_{jets} #geq 1)",
                                 nZPt_Zinc1jetQunJZB,
                                 zPt_Zinc1jetQunJZB,
                                 nZPt_Zinc1jetQunJZB,
                                 zPt_Zinc1jetQunJZB);
    hresponseJZB_ptHigh = newTH2D("hresponseJZB_ptHigh",
                                  "hresp Scalar JZB  p_{T} (N_{jets} #geq 1)",
                                  nZPt_Zinc2jetQunJZBptHigh,
                                  zPt_Zinc2jetQunJZBptHigh,
                                  nZPt_Zinc2jetQunJZBptHigh,
                                  zPt_Zinc2jetQunJZBptHigh);
    hresponseZPt_Zinc2jet = newTH2D("hresponseZPt_Zinc2jet",
                                    "hresponseZPt_Zinc2jet",
                                    nZPt_Zinc2jet,
                                    zPt_Zinc2jet,
                                    nZPt_Zinc2jet,
                                    zPt_Zinc2jet);
    hresponseZAbsRapidity_Zinc1jet = newTH2D(
        "hresponseZAbsRapidity_Zinc1jet", "hresp ZAbsRapidity_Zinc1jet", 12, 0., 2.4, 12, 0., 2.4);

    hresponseFirstJetPt_Zinc1jet = newTH2D("hresponseFirstJetPt_Zinc1jet",
                                           "hresp 1st jet pt",
                                           nJetPt_Zinc1jet,
                                           jetPt_Zinc1jet,
                                           nJetPt_Zinc1jet,
                                           jetPt_Zinc1jet);
    hresponseFirstJetPtMatch_Zinc1jet = newTH2D("hresponseFirstJetPtMatch_Zinc1jet",
                                                "hresp 1st jet pt match",
                                                nJetPt_Zinc1jet,
                                                jetPt_Zinc1jet,
                                                nJetPt_Zinc1jet,
                                                jetPt_Zinc1jet);

    hresponseFirstJetPt_2_Zinc1jet = newTH2D("hresponseFirstJetPt_2_Zinc1jet",
                                             "hresp 1st jet pt (2)",
                                             jetPt_2_Zinc1jet,
                                             jetPt_2_Zinc1jet);
    hresponseSecondJetPt_Zinc2jet = newTH2D("hresponseSecondJetPt_Zinc2jet",
                                            "hresp 2nd jet pt",
                                            nJetPt_Zinc2jet,
                                            jetPt_Zinc2jet,
                                            nJetPt_Zinc2jet,
                                            jetPt_Zinc2jet);
    hresponseSecondJetPtMatch_Zinc2jet = newTH2D("hresponseSecondJetPtMatch_Zinc2jet",
                                                 "hresp 2nd jet pt Match",
                                                 nJetPt_Zinc2jet,
                                                 jetPt_Zinc2jet,
                                                 nJetPt_Zinc2jet,
                                                 jetPt_Zinc2jet); // DJALOG
    hresponseSecondJetPt_2_Zinc2jet = newTH2D("hresponseSecondJetPt_2_Zinc2jet",
                                              "hresp 2nd jet pt (2)",
                                              jetPt_2_Zinc2jet,
                                              jetPt_2_Zinc2jet);
    hresponseThirdJetPt_Zinc3jet = newTH2D("hresponseThirdJetPt_Zinc3jet",
                                           "hresp 3rd jet pt",
                                           nJetPt_Zinc3jet,
                                           jetPt_Zinc3jet,
                                           nJetPt_Zinc3jet,
                                           jetPt_Zinc3jet);
    hresponseThirdJetPtMatch_Zinc3jet = newTH2D("hresponseThirdJetPtMatch_Zinc3jet",
                                                "hresp 3rd jet pt match",
                                                nJetPt_Zinc3jet,
                                                jetPt_Zinc3jet,
                                                nJetPt_Zinc3jet,
                                                jetPt_Zinc3jet); // DJALOG
    hresponseThirdJetPt_2_Zinc3jet = newTH2D("hresponseThirdJetPt_2_Zinc3jet",
                                             "hresp 3rd jet pt (2)",
                                             jetPt_2_Zinc3jet,
                                             jetPt_2_Zinc3jet);
    hresponseFourthJetPt_Zinc4jet = newTH2D("hresponseFourthJetPt_Zinc4jet",
                                            "hresp 4th jet pt",
                                            nJetPt_Zinc4jet,
                                            jetPt_Zinc4jet,
                                            nJetPt_Zinc4jet,
                                            jetPt_Zinc4jet);
    hresponseFourthJetPtMatch_Zinc4jet = newTH2D("hresponseFourthJetPtMatch_Zinc4jet",
                                                 "hresp 4th jet pt match",
                                                 nJetPt_Zinc4jet,
                                                 jetPt_Zinc4jet,
                                                 nJetPt_Zinc4jet,
                                                 jetPt_Zinc4jet); // DJALOG
    hresponseFifthJetPt_Zinc5jet = newTH2D("hresponseFifthJetPt_Zinc5jet",
                                           "hresp 5th jet pt",
                                           nJetPt_Zinc5jet,
                                           jetPt_Zinc5jet,
                                           nJetPt_Zinc5jet,
                                           jetPt_Zinc5jet);
    hresponseJetsHT_Zinc1jet = newTH2D("hresponseJetsHT_Zinc1jet",
                                       "hresp Scalar sum jets p_{T} (N_{jets} #geq 1)",
                                       nJetHT_Zinc1jet,
                                       jetHT_Zinc1jet,
                                       nJetHT_Zinc1jet,
                                       jetHT_Zinc1jet);
    hresponseJetsHT_2_Zinc1jet = newTH2D("hresponseJetsHT_2_Zinc1jet",
                                         "hresp Scalar sum jets p_{T} (N_{jets} #geq 1)2",
                                         jetHT_2_Zinc1jet,
                                         jetHT_2_Zinc1jet);
    hresponseJetsHT_Zinc2jet = newTH2D("hresponseJetsHT_Zinc2jet",
                                       "hresp Scalar sum jets p_{T} (N_{jets} #geq 2)",
                                       nJetHT_Zinc2jet,
                                       jetHT_Zinc2jet,
                                       nJetHT_Zinc2jet,
                                       jetHT_Zinc2jet);
    hresponseVisPt_Zinc0jetQun = newTH2D("hresponseVisPt_Zinc0jetQun",
                                         "hresp Scalar visible  p_{T} (N_{jets} #geq 0)",
                                         nZPt_Zinc1jet,
                                         zPt_Zinc1jet,
                                         nZPt_Zinc1jet,
                                         zPt_Zinc1jet);
    hresponseVisPt_Zinc1jetQun = newTH2D("hresponseVisPt_Zinc1jetQun",
                                         "hresp Scalar visible  p_{T} (N_{jets} #geq 1)",
                                         nZPt_Zinc2jetQun,
                                         zPt_Zinc2jetQun,
                                         nZPt_Zinc2jetQun,
                                         zPt_Zinc2jetQun);
    hresponseVisPt_Zinc2jetQun = newTH2D("hresponseVisPt_Zinc2jetQun",
                                         "hresp Scalar visible  p_{T} (N_{jets} #geq 2)",
                                         nZPt_Zinc2JetQun,
                                         zPt_Zinc2JetQun,
                                         nZPt_Zinc2JetQun,
                                         zPt_Zinc2JetQun);
    hresponseVisPt_Zinc3jetQun = newTH2D("hresponseVisPt_Zinc3jetQun",
                                         "hresp Scalar visible  p_{T} (N_{jets} #geq 3)",
                                         nZPt_Zinc3jetQun,
                                         zPt_Zinc3jetQun,
                                         nZPt_Zinc3jetQun,
                                         zPt_Zinc3jetQun);

    hresponseJetsHT_2_Zinc2jet = newTH2D("hresponseJetsHT_2_Zinc2jet",
                                         "hresp Scalar sum jets p_{T} (N_{jets} #geq 2)2",
                                         jetHT_2_Zinc2jet,
                                         jetHT_2_Zinc2jet);

    hresponseJetsHT_Zinc3jet = newTH2D("hresponseJetsHT_Zinc3jet",
                                       "hresp Scalar sum jets p_{T} (N_{jets} #geq 3)",
                                       nJetHT_Zinc3jet,
                                       jetHT_Zinc3jet,
                                       nJetHT_Zinc3jet,
                                       jetHT_Zinc3jet);

    hresponseJetsHT_2_Zinc3jet = newTH2D("hresponseJetsHT_2_Zinc3jet",
                                         "hresp Scalar sum jets p_{T} (N_{jets} #geq 3)2",
                                         jetHT_2_Zinc3jet,
                                         jetHT_2_Zinc3jet);

    hresponseJetsHT_Zinc4jet = newTH2D("hresponseJetsHT_Zinc4jet",
                                       "hresp Scalar sum jets p_{T} (N_{jets} #geq 4)",
                                       nJetHT_Zinc4jet,
                                       jetHT_Zinc4jet,
                                       nJetHT_Zinc4jet,
                                       jetHT_Zinc4jet);
    hresponseJetsHT_Zinc5jet = newTH2D("hresponseJetsHT_Zinc5jet",
                                       "hresp Scalar sum jets p_{T} (N_{jets} #geq 5)",
                                       nJetHT_Zinc5jet,
                                       jetHT_Zinc5jet,
                                       nJetHT_Zinc5jet,
                                       jetHT_Zinc5jet);

    hresponseFirstJetEta_Zinc1jet = newTH2D("hresponseFirstJetEta_Zinc1jet",
                                            "hresp 1st jet #eta (N_{jets} #geq 1)",
                                            32,
                                            0,
                                            2.4,
                                            32,
                                            0,
                                            2.4);
    hresponseFirstJetEta_2_Zinc1jet = newTH2D("hresponseFirstJetEta_2_Zinc1jet",
                                              "hresp 1st jet #eta (N_{jets} #geq 1)2",
                                              160,
                                              0,
                                              2.4,
                                              160,
                                              0,
                                              2.4);
    // DJALOG
    hresponseFirstJetEtaMatch_Zinc1jet = newTH2D("hresponseFirstJetEtaMatch_Zinc1jet",
                                                 "hresp 1st jet Match #eta (N_{jets} #geq 1)",
                                                 32,
                                                 0,
                                                 2.4,
                                                 32,
                                                 0,
                                                 2.4);
    hresponseSecondJetEta_Zinc2jet = newTH2D("hresponseSecondJetEta_Zinc2jet",
                                             "hresp 2nd jet #eta (N_{jets} #geq 2)",
                                             24,
                                             0,
                                             2.4,
                                             24,
                                             0,
                                             2.4);
    hresponseSecondJetEtaMatch_Zinc2jet = newTH2D("hresponseSecondJetEtaMatch_Zinc2jet",
                                                  "hresp 2nd jet match #eta (N_{jets} #geq 2)",
                                                  24,
                                                  0,
                                                  2.4,
                                                  24,
                                                  0,
                                                  2.4); // DJALOG
    hresponseSecondJetEta_2_Zinc2jet = newTH2D("hresponseSecondJetEta_2_Zinc2jet",
                                               "hresp 2nd jet #eta (N_{jets} #geq 2)2",
                                               120,
                                               0,
                                               2.4,
                                               120,
                                               0,
                                               2.4);
    hresponseThirdJetEta_Zinc3jet = newTH2D("hresponseThirdJetEta_Zinc3jet",
                                            "hresp 3rd jet #eta (N_{jets} #geq 3)",
                                            12,
                                            0,
                                            2.4,
                                            12,
                                            0,
                                            2.4);
    hresponseThirdJetEtaMatch_Zinc3jet = newTH2D("hresponseThirdJetEtaMatch_Zinc3jet",
                                                 "hresp 3rd jet match #eta (N_{jets} #geq 3)",
                                                 12,
                                                 0,
                                                 2.4,
                                                 12,
                                                 0,
                                                 2.4); // DJALOG
    hresponseThirdJetEta_2_Zinc3jet = newTH2D("hresponseThirdJetEta_2_Zinc3jet",
                                              "hresp 3rd jet #eta (N_{jets} #geq 3)2",
                                              60,
                                              0,
                                              2.4,
                                              60,
                                              0,
                                              2.4);
    hresponseFourthJetEta_Zinc4jet = newTH2D("hresponseFourthJetEta_Zinc4jet",
                                             "hresp 4th jet #eta (N_{jets} #geq 4)",
                                             8,
                                             0,
                                             2.4,
                                             8,
                                             0,
                                             2.4);
    hresponseFourthJetEtaMatch_Zinc4jet = newTH2D("hresponseFourthJetEtaMatch_Zinc4jet",
                                                  "hresp 4th jet match #eta (N_{jets} #geq 4)",
                                                  8,
                                                  0,
                                                  2.4,
                                                  8,
                                                  0,
                                                  2.4); // DJALOG
    hresponseFifthJetEta_Zinc5jet = newTH2D("hresponseFifthJetEta_Zinc5jet",
                                            "hresp 5th jet #eta (N_{jets} #geq 5)",
                                            6,
                                            0,
                                            2.4,
                                            6,
                                            0,
                                            2.4);

    hresponseFirstJetEtaHigh_Zinc1jet = newTH2D("hresponseFirstJetEtaHigh_Zinc1jet",
                                                "hresp 1st jet |#eta| (N_{jets} #geq 1)",
                                                47,
                                                0,
                                                4.7,
                                                47,
                                                0,
                                                4.7);
    hresponseFirstJetEtaHighMatch_Zinc1jet = newTH2D("hresponseFirstJetEtaHighMatch_Zinc1jet",
                                                     "hresp 1st jet match |#eta| (N_{jets} #geq 1)",
                                                     47,
                                                     0,
                                                     4.7,
                                                     47,
                                                     0,
                                                     4.7); // DJALOG
    hresponseSecondJetEtaHigh_Zinc2jet = newTH2D("hresponseSecondJetEtaHigh_Zinc2jet",
                                                 "hresp 2nd jet |#eta| (N_{jets} #geq 2)",
                                                 47,
                                                 0,
                                                 4.7,
                                                 47,
                                                 0,
                                                 4.7);
    hresponseSecondJetEtaHighMatch_Zinc2jet =
        newTH2D("hresponseSecondJetEtaHighMatch_Zinc2jet",
                "hresp 2nd jet match |#eta| (N_{jets} #geq 2)",
                47,
                0,
                4.7,
                47,
                0,
                4.7); // DJALOG
    hresponseThirdJetEtaHigh_Zinc3jet = newTH2D("hresponseThirdJetEtaHigh_Zinc3jet",
                                                "hresp 3rd jet |#eta| (N_{jets} #geq 3)",
                                                24,
                                                0,
                                                4.7,
                                                24,
                                                0,
                                                4.7);
    hresponseThirdJetEtaHighMatch_Zinc3jet = newTH2D("hresponseThirdJetEtaHighMatch_Zinc3jet",
                                                     "hresp 3rd jet match |#eta| (N_{jets} #geq 3)",
                                                     24,
                                                     0,
                                                     4.7,
                                                     24,
                                                     0,
                                                     4.7); // DJALOG
    hresponseFourthJetEtaHigh_Zinc4jet = newTH2D("hresponseFourthJetEtaHigh_Zinc4jet",
                                                 "hresp 4th jet |#eta| (N_{jets} #geq 4)",
                                                 12,
                                                 0,
                                                 4.7,
                                                 12,
                                                 0,
                                                 4.7);
    hresponseFourthJetEtaHighMatch_Zinc4jet =
        newTH2D("hresponseFourthJetEtaHighMatch_Zinc4jet",
                "hresp 4th jet match |#eta| (N_{jets} #geq 4)",
                12,
                0,
                4.7,
                12,
                0,
                4.7); // DJALOG
    hresponseFifthJetEtaHigh_Zinc5jet = newTH2D("hresponseFifthJetEtaHigh_Zinc5jet",
                                                "hresp 5th jet |#eta| (N_{jets} #geq 5)",
                                                6,
                                                0,
                                                4.7,
                                                6,
                                                0,
                                                4.7);

    hresponseFirstJetAbsRapidity_Zinc1jet = newTH2D("hresponseFirstJetAbsRapidity_Zinc1jet",
                                                    "hresp 1st jet |y| (N_{jets} #geq 1)",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseFirstJetAbsRapidityMatch_Zinc1jet =
        newTH2D("hresponseFirstJetAbsRapidityMatch_Zinc1jet",
                "hresp 1st jet match |y| (N_{jets} #geq 1)",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    // DJALOG
    hresponseFirstJetAbsRapidityMatch_Zinc1jet =
        newTH2D("hresponseFirstJetAbsRapidityMatch_Zinc1jet",
                "hresp 1st jet |y| Match (N_{jets} #geq 1)",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseFirstJetRapidityHigh_Zinc1jet = newTH2D("hresponseFirstJetRapidityHigh_Zinc1jet",
                                                     "hresp 1st jet |y| (N_{jets} #geq 1)",
                                                     47,
                                                     0,
                                                     4.7,
                                                     47,
                                                     0,
                                                     4.7);
    hresponseFirstJetRapidityHighMatch_Zinc1jet =
        newTH2D("hresponseFirstJetRapidityHighMatch_Zinc1jet",
                "hresp 1st jet match |y| (N_{jets} #geq 1)",
                47,
                0,
                4.7,
                47,
                0,
                4.7); // DJALOG
    hresponseSecondJetAbsRapidity_Zinc2jet = newTH2D("hresponseSecondJetAbsRapidity_Zinc2jet",
                                                     "hresp 2nd jet |y| (N_{jets} #geq 2)",
                                                     12,
                                                     0,
                                                     2.4,
                                                     12,
                                                     0,
                                                     2.4);
    hresponseSecondJetRapidityHigh_Zinc2jet = newTH2D("hresponseSecondJetRapidityHigh_Zinc2jet",
                                                      "hresp 2nd jet |y| (N_{jets} #geq 2)",
                                                      47,
                                                      0,
                                                      4.7,
                                                      47,
                                                      0,
                                                      4.7);
    hresponseSecondJetAbsRapidityMatch_Zinc2jet =
        newTH2D("hresponseSecondJetAbsRapidityMatch_Zinc2jet",
                "hresp 2nd jet match |y| (N_{jets} #geq 2)",
                12,
                0,
                2.4,
                12,
                0,
                2.4); // DJALOG
    hresponseSecondJetRapidityHighMatch_Zinc2jet =
        newTH2D("hresponseSecondJetRapidityHighMatch_Zinc2jet",
                "hresp 2nd jet match |y| (N_{jets} #geq 2)",
                47,
                0,
                4.7,
                47,
                0,
                4.7); // DJALOG

    hresponseThirdJetAbsRapidity_Zinc3jet = newTH2D("hresponseThirdJetAbsRapidity_Zinc3jet",
                                                    "hresp 3rd jet |y| (N_{jets} #geq 3)",
                                                    8,
                                                    0,
                                                    2.4,
                                                    8,
                                                    0,
                                                    2.4);
    hresponseThirdJetRapidityHigh_Zinc3jet = newTH2D("hresponseThirdJetRapidityHigh_Zinc3jet",
                                                     "hresp 3rd jet |y| (N_{jets} #geq 3)",
                                                     24,
                                                     0,
                                                     4.7,
                                                     24,
                                                     0,
                                                     4.7);
    hresponseThirdJetAbsRapidityMatch_Zinc3jet =
        newTH2D("hresponseThirdJetAbsRapidityMatch_Zinc3jet",
                "hresp 3rd jet match |y| (N_{jets} #geq 3)",
                8,
                0,
                2.4,
                8,
                0,
                2.4); // DJALOG
    hresponseThirdJetRapidityHighMatch_Zinc3jet =
        newTH2D("hresponseThirdJetRapidityHighMatch_Zinc3jet",
                "hresp 3rd jet match |y| (N_{jets} #geq 3)",
                24,
                0,
                4.7,
                24,
                0,
                4.7); // DJALOG
    hresponseFourthJetAbsRapidity_Zinc4jet = newTH2D("hresponseFourthJetAbsRapidity_Zinc4jet",
                                                     "hresp 4th jet |y| (N_{jets} #geq 4)",
                                                     8,
                                                     0,
                                                     2.4,
                                                     8,
                                                     0,
                                                     2.4);
    hresponseFourthJetRapidityHigh_Zinc4jet = newTH2D("hresponseFourthJetRapidityHigh_Zinc4jet",
                                                      "hresp 4th jet |y| (N_{jets} #geq 4)",
                                                      12,
                                                      0,
                                                      4.7,
                                                      12,
                                                      0,
                                                      4.7);
    hresponseFourthJetAbsRapidityMatch_Zinc4jet =
        newTH2D("hresponseFourthJetAbsRapidityMatch_Zinc4jet",
                "hresp 4th jet match |y| (N_{jets} #geq 4)",
                8,
                0,
                2.4,
                8,
                0,
                2.4); // DJALOG
    hresponseFourthJetRapidityHighMatch_Zinc4jet =
        newTH2D("hresponseFourthJetRapidityHighMatch_Zinc4jet",
                "hresp 4th jet match |y| (N_{jets} #geq 4)",
                12,
                0,
                4.7,
                12,
                0,
                4.7); // DJALOG
    hresponseFifthJetAbsRapidity_Zinc5jet = newTH2D("hresponseFifthJetAbsRapidity_Zinc5jet",
                                                    "hresp 5th jet |y| (N_{jets} #geq 5)",
                                                    6,
                                                    0,
                                                    2.4,
                                                    6,
                                                    0,
                                                    2.4);
    hresponseFifthJetRapidityHigh_Zinc5jet = newTH2D("hresponseFifthJetRapidityHigh_Zinc5jet",
                                                     "hresp 5th jet |y| (N_{jets} #geq 5)",
                                                     6,
                                                     0,
                                                     4.7,
                                                     6,
                                                     0,
                                                     4.7);

    hresponseFirstJetPtEta_Zinc1jet = newTH2D("hresponseFirstJetPtEta_Zinc1jet",
                                              "hresponseFirstJetPtEta_Zinc1jet",
                                              10 * 6,
                                              0,
                                              10 * 6,
                                              10 * 6,
                                              0,
                                              10 * 6);

    ZNGoodJetsNVtx_Zexc = newTH2D(
        "ZNGoodJetsNVtx_Zexc", "NVtx vs Jet Counter (excl.)", 11, -0.5, 10.5, 45, 0.5, 45.5);
    if (ZNGoodJetsNVtx_Zexc) {
        for (int ibin = 1; ibin < ZNGoodJetsNVtx_Zexc->GetNbinsX(); ++ibin) {
            ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(ibin, TString::Format("= %d", ibin - 1));
        }
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(1, "= 0");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(2, "= 1");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(3, "= 2");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(4, "= 3");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(5, "= 4");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(6, "= 5");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(7, "= 6");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(8, "= 7");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(9, "= 8");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(10,"= 9");
        //	ZNGoodJetsNVtx_Zexc->GetXaxis()->SetBinLabel(11,"= 10");
    }

    ZNGoodJets20NVtx_Zexc = newTH2D(
        "ZNGoodJets20NVtx_Zexc", "NVtx vs Jet20 Counter (excl.)", 11, -0.5, 10.5, 45, 0.5, 45.5);
    if (ZNGoodJetsNVtx_Zexc) {
        for (int ibin = 1; ibin < ZNGoodJets20NVtx_Zexc->GetNbinsX(); ++ibin) {
            ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(ibin, TString::Format("= %d", ibin - 1));
        }
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(1, "= 0");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(2, "= 1");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(3, "= 2");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(4, "= 3");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(5, "= 4");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(6, "= 5");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(7, "= 6");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(8, "= 7");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(9, "= 8");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(10,"= 9");
        //	ZNGoodJets20NVtx_Zexc->GetXaxis()->SetBinLabel(11,"= 10");	}
    }

    ZNGoodJets_Zinc = newTH1D("ZNGoodJets_Zinc", "Jet Counter (incl.)", "N_{jets}", 7, -0.5, 6.5);
    ZNGoodJets_SameChargePair_Zinc = newTH1D("ZNGoodJets_SameChargePair_Zinc",
                                             "Jet Counter Same Sign Leptons(incl.)",
                                             "N_{jets}",
                                             7,
                                             -0.5,
                                             6.5);
    if (ZNGoodJets_Zinc) {
        for (int ibin = 1; ibin < ZNGoodJets_Zinc->GetNbinsX(); ++ibin) {
            ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(ibin, TString::Format("#geq %d", ibin));
            ZNGoodJets_SameChargePair_Zinc->GetXaxis()->SetBinLabel(
                ibin, TString::Format("#geq %d", ibin));
        }
        //    ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(1, "#geq 0");
        //    ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(2, "#geq 1");
        //    ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(3, "#geq 2");
        //    ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(4, "#geq 3");
        //    ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(5, "#geq 4");
        //    ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(6, "#geq 5");
        //    ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(7, "#geq 6");
        ////  ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(8, "#geq 7");
        ////  ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(9, "#geq 8");
        ////  ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(10,"#geq 9");
        ////  ZNGoodJets_Zinc->GetXaxis()->SetBinLabel(11,"#geq 10");
    }

    ZNGoodJets_Zinc_nvtx10 =
        newTH1D("ZNGoodJets_Zinc_nvtx10", "1 Jet Counter (incl.)", "N_{jets}", 7, -0.5, 6.5);
    ZNGoodJets_Zinc_nvtx20 =
        newTH1D("ZNGoodJets_Zinc_nvtx20", "2 Jet Counter (incl.)", "N_{jets}", 7, -0.5, 6.5);
    ZNGoodJets_Zinc_nvtx30 =
        newTH1D("ZNGoodJets_Zinc_nvtx30", "3 Jet Counter (incl.)", "N_{jets}", 7, -0.5, 6.5);
    ZNGoodJets_Zinc_nvtx45 =
        newTH1D("ZNGoodJets_Zinc_nvtx45", "4 Jet Counter (incl.)", "N_{jets}", 7, -0.5, 6.5);

    ZNGoodJets_Zexc_NoWeight = newTH1D(
        "ZNGoodJets_Zexc_NoWeight", "Unweighted jet Counter (excl.)", "N_{jets}", 8, -0.5, 7.5);
    if (ZNGoodJets_Zexc_NoWeight) {
        for (int ibin = 1; ibin < ZNGoodJets_Zexc_NoWeight->GetNbinsX(); ++ibin) {
            ZNGoodJets_Zexc_NoWeight->GetXaxis()->SetBinLabel(ibin,
                                                              TString::Format("= %d", ibin - 1));
        }

        //	ZNGoodJets_Zexc_NoWeight->GetXaxis()->SetBinLabel(1,"= 0");
        //	ZNGoodJets_Zexc_NoWeight->GetXaxis()->SetBinLabel(2,"= 1");
        //	ZNGoodJets_Zexc_NoWeight->GetXaxis()->SetBinLabel(3,"= 2");
        //	ZNGoodJets_Zexc_NoWeight->GetXaxis()->SetBinLabel(4,"= 3");
        //	ZNGoodJets_Zexc_NoWeight->GetXaxis()->SetBinLabel(5,"= 4");
        //	ZNGoodJets_Zexc_NoWeight->GetXaxis()->SetBinLabel(6,"= 5");
        //	ZNGoodJets_Zexc_NoWeight->GetXaxis()->SetBinLabel(7,"= 6");
        //	ZNGoodJets_Zexc_NoWeight->GetXaxis()->SetBinLabel(8,"= 7");
    }

    ZNGoodJets_Zinc_NoWeight = newTH1D(
        "ZNGoodJets_Zinc_NoWeight", "Unweighted jet Counter (incl.)", "N_{jets}", 8, -0.5, 7.5);
    if (ZNGoodJets_Zinc_NoWeight) {
        for (int ibin = 1; ibin < ZNGoodJets_Zinc_NoWeight->GetNbinsX(); ++ibin) {
            ZNGoodJets_Zinc_NoWeight->GetXaxis()->SetBinLabel(ibin,
                                                              TString::Format("#geq %d", ibin));
        }
        //	ZNGoodJets_Zinc_NoWeight->GetXaxis()->SetBinLabel(1,"#geq 0");
        //	ZNGoodJets_Zinc_NoWeight->GetXaxis()->SetBinLabel(2,"#geq 1");
        //	ZNGoodJets_Zinc_NoWeight->GetXaxis()->SetBinLabel(3,"#geq 2");
        //	ZNGoodJets_Zinc_NoWeight->GetXaxis()->SetBinLabel(4,"#geq 3");
        //	ZNGoodJets_Zinc_NoWeight->GetXaxis()->SetBinLabel(5,"#geq 4");
        //	ZNGoodJets_Zinc_NoWeight->GetXaxis()->SetBinLabel(6,"#geq 5");
        //	ZNGoodJets_Zinc_NoWeight->GetXaxis()->SetBinLabel(7,"#geq 6");
        //	ZNGoodJets_Zinc_NoWeight->GetXaxis()->SetBinLabel(8,"#geq 7");
    }

    // DPS histograms
    // binning
    int nbinSpt = 21;
    double binSpt[22] = {0,   0.05, 0.1, 0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  .45,  .5,
                         .55, .6,   .65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.94, 0.98, 1};

    //-- jets and Z
    TwoJetsPtDiff_Zexc2jet = newTH1D("TwoJetsPtDiff_Zexc2jet",
                                     "pT diff of the two highest jet (N_{jets} = 2)",
                                     "#Delta pT(j_{1}j_{2}) [GeV]",
                                     10,
                                     0,
                                     100);
    genTwoJetsPtDiff_Zexc2jet = newTH1D("genTwoJetsPtDiff_Zexc2jet",
                                        "gen pT diff of the two highest jet (N_{jets} = 2)",
                                        "#Delta pT(j_{1}j_{2}) [GeV]",
                                        10,
                                        0,
                                        100);
    JetsMass_Zexc2jet =
        newTH1D("JetsMass_Zexc2jet", "2Jets Invariant Mass (N_{jets} = 2)", Mjj, 25, 0, 625);
    genJetsMass_Zexc2jet =
        newTH1D("genJetsMass_Zexc2jet", "gen 2Jets Invariant Mass (N_{jets} = 2)", Mjj, 25, 0, 625);
    ptBal_Zexc2jet = newTH1D("ptBal_Zexc2jet",
                             "Vectorial pT sum: Z_{pT} + DiJet_{pT} (N_{jets} = 2)",
                             "#Sigma pT [GeV]",
                             50,
                             0,
                             100);

    genptBal_Zexc2jet = newTH1D("genptBal_Zexc2jet",
                                "gen Vectorial pT sum: Z_{pT} + DiJet_{pT} (N_{jets} = 2)",
                                "#Sigma pT [GeV]",
                                50,
                                0,
                                100);
    dPhiJets_Zexc2jet =
        newTH1D("dPhiJets_Zexc2jet", "#Delta#phi btwn jets (N_{jets} = 2)", jdPhi, 20, 0, PI);
    gendPhiJets_Zexc2jet = newTH1D(
        "gendPhiJets_Zexc2jet", "gen #Delta#phi btwn jets (N_{jets} = 2)", jdPhi, 20, 0, PI);
    dEtaJets_Zexc2jet =
        newTH1D("dEtaJets_Zexc2jet", "#Delta#eta btwn jets (N_{jets} = 2)", jdEta, 48, 0, 4.8);
    gendEtaJets_Zexc2jet = newTH1D(
        "gendEtaJets_Zexc2jet", "gen #Delta#eta btwn jets (N_{jets} = 2)", jdEta, 48, 0, 4.8);
    dEtaFirstJetZ_Zexc2jet = newTH1D("dEtaFirstJetZ_Zexc2jet",
                                     "#Delta#eta btwn Jet_{1} and Z (N_{jets} = 2)",
                                     "#Delta#eta(j_{1}Z)",
                                     50,
                                     -6,
                                     6);
    gendEtaFirstJetZ_Zexc2jet = newTH1D("gendEtaFirstJetZ_Zexc2jet",
                                        "gen #Delta#eta btwn Jet_{1} and Z (N_{jets} = 2)",
                                        "#Delta#eta(j_{1}Z)",
                                        50,
                                        -6,
                                        6);
    dEtaSecondJetZ_Zexc2jet = newTH1D("dEtaSecondJetZ_Zexc2jet",
                                      "#Delta#eta btwn Jet_{2} and Z (N_{jets} = 2)",
                                      "#Delta#eta(j_{2}Z)",
                                      50,
                                      -6,
                                      6);
    gendEtaSecondJetZ_Zexc2jet = newTH1D("gendEtaSecondJetZ_Zexc2jet",
                                         "gen #Delta#eta btwn Jet_{2} and Z (N_{jets} = 2)",
                                         "#Delta#eta(j_{2}Z)",
                                         50,
                                         -6,
                                         6);
    dEtaJet1Plus2Z_Zexc2jet = newTH1D("dEtaJet1Plus2Z_Zexc2jet",
                                      "#Delta#eta btwn jets and Z (N_{jets} = 2)",
                                      "#Delta#eta(j_{12}Z)",
                                      120,
                                      -6,
                                      6);
    gendEtaJet1Plus2Z_Zexc2jet = newTH1D("gendEtaJet1Plus2Z_Zexc2jet",
                                         "gen #Delta#eta btwn jets and Z (N_{jets} = 2)",
                                         "#Delta#eta(j_{12}Z)",
                                         120,
                                         -6,
                                         6);
    PHI_Zexc2jet = newTH1D("PHI_Zexc2jet",
                           "#phi: Angle btwn the two subsystems planes (N_{jets} = 2)",
                           "#phi(j_{12}Z)",
                           25,
                           0,
                           PI);
    genPHI_Zexc2jet = newTH1D("genPHI_Zexc2jet",
                              "gen #phi: Angle btwn the two subsystems planes (N_{jets} = 2)",
                              "#phi(j_{12}Z)",
                              25,
                              0,
                              PI);
    PHI_T_Zexc2jet = newTH1D("PHI_T_Zexc2jet",
                             "#Delta S Angle btwn lep and jet pair in T-plane (N_{jets} = 2)",
                             "#Delta S(j_{12}Z)",
                             10,
                             0,
                             PI);
    genPHI_T_Zexc2jet =
        newTH1D("genPHI_T_Zexc2jet",
                "gen #Delta S Angle btwn lep and jet pair in T-plane (N_{jets} = 2)",
                "#Delta S(j_{12}Z)",
                10,
                0,
                PI);
    SpT_Zexc2jet = newTH1D(
        "SpT_Zexc2jet", "#Delta_{pT}^{rel} lep and jets combined (N_{jets} = 2)", Spt, 20, 0, 1);
    genSpT_Zexc2jet = newTH1D("genSpT_Zexc2jet",
                              "gen #Delta_{pT}^{rel} lep and jets combined (N_{jets} = 2)",
                              Spt,
                              20,
                              0,
                              1);
    SpTJets_Zexc2jet =
        newTH1D("SpTJets_Zexc2jet", "#Delta_{pT}^{rel} jets (N_{jets} = 2)", jSpt, 20, 0, 1);
    genSpTJets_Zexc2jet =
        newTH1D("genSpTJets_Zexc2jet", "gen #Delta_{pT}^{rel} jets (N_{jets} = 2)", jSpt, 20, 0, 1);
    SPhi_Zexc2jet =
        newTH1D("SPhi_Zexc2jet", "S_{#phi} lep and jets combined (N_{jets} = 2)", Sphi, 50, 0, PI);
    genSPhi_Zexc2jet = newTH1D(
        "genSPhi_Zexc2jet", "gen S_{#phi} lep and jets combined (N_{jets} = 2)", Sphi, 50, 0, PI);

    TwoJetsPtDiff_Zinc2jet = newTH1D("TwoJetsPtDiff_Zinc2jet",
                                     "pT diff of the two highest jet (N_{jets} #geq 2)",
                                     "#Delta pT(j_{1}j_{2}) [GeV]",
                                     10,
                                     0,
                                     100);
    genTwoJetsPtDiff_Zinc2jet = newTH1D("genTwoJetsPtDiff_Zinc2jet",
                                        "gen pT diff of the two highest jet (N_{jets} #geq 2)",
                                        "#Delta pT(j_{1}j_{2}) [GeV]",
                                        10,
                                        0,
                                        100);
    BestTwoJetsPtDiff_Zinc2jet = newTH1D("BestTwoJetsPtDiff_Zinc2jet",
                                         "Best pT diff of the two highest jet (N_{jets} #geq 2)",
                                         "#Delta pT(j_{1}j_{2}) [GeV]",
                                         10,
                                         0,
                                         100);
    genBestTwoJetsPtDiff_Zinc2jet =
        newTH1D("genBestTwoJetsPtDiff_Zinc2jet",
                "gen Best pT diff of the two highest jet (N_{jets} #geq 2)",
                "#Delta pT(j_{1}j_{2}) [GeV]",
                10,
                0,
                100);

    llJetsMass_Zinc2jet = newTH1D(
        "llJetsMass_Zinc2jet", "ll+2Jets Invariant Mass (N_{jets} #geq 2)", "M(lljj)", 50, 0, 1000);
    genllJetsMass_Zinc2jet = newTH1D("genllJetsMass_Zinc2jet",
                                     "gen ll+2Jets Invariant Mass (N_{jets} #geq 2)",
                                     "M(lljj)",
                                     50,
                                     0,
                                     1000);

    JetsMass_Zinc2jet = newTH1D("JetsMass_Zinc2jet",
                                "2Jets Invariant Mass (N_{jets} #geq 2)",
                                Mjj,
                                nJetsMass_Zinc2jet,
                                jetsMass_Zinc2jet);
    genJetsMass_Zinc2jet = newTH1D("genJetsMass_Zinc2jet",
                                   "gen 2Jets Invariant Mass (N_{jets} #geq 2)",
                                   Mjj,
                                   nJetsMass_Zinc2jet,
                                   jetsMass_Zinc2jet);
    hresponseJetsMass_Zinc2jet = newTH2D("hresponseJetsMass_Zinc2jet",
                                         "hresponseJetsMass_Zinc2jet",
                                         nJetsMass_Zinc2jet,
                                         jetsMass_Zinc2jet,
                                         nJetsMass_Zinc2jet,
                                         jetsMass_Zinc2jet);

    JetsMassLowPU_Zinc2jet = newTH1D("JetsMassLowPU_Zinc2jet",
                                     "2Jets Invariant Mass LowPU(N_{jets} #geq 2)",
                                     Mjj,
                                     nJetsMass_Zinc2jet,
                                     jetsMass_Zinc2jet);
    genJetsMassLowPU_Zinc2jet = newTH1D("genJetsMassLowPU_Zinc2jet",
                                        "gen 2Jets Invariant Mass LowPU(N_{jets} #geq 2)",
                                        Mjj,
                                        nJetsMass_Zinc2jet,
                                        jetsMass_Zinc2jet);
    hresponseJetsMassLowPU_Zinc2jet = newTH2D("hresponseJetsMassLowPU_Zinc2jet",
                                              "hresponseJetsMassLowPU_Zinc2jet",
                                              nJetsMass_Zinc2jet,
                                              jetsMass_Zinc2jet,
                                              nJetsMass_Zinc2jet,
                                              jetsMass_Zinc2jet);

    JetsMassMidPU_Zinc2jet = newTH1D("JetsMassMidPU_Zinc2jet",
                                     "2Jets Invariant Mass MidPU(N_{jets} #geq 2)",
                                     Mjj,
                                     nJetsMass_Zinc2jet,
                                     jetsMass_Zinc2jet);
    genJetsMassMidPU_Zinc2jet = newTH1D("genJetsMassMidPU_Zinc2jet",
                                        "gen 2Jets Invariant Mass MidPU(N_{jets} #geq 2)",
                                        Mjj,
                                        nJetsMass_Zinc2jet,
                                        jetsMass_Zinc2jet);
    hresponseJetsMassMidPU_Zinc2jet = newTH2D("hresponseJetsMassMidPU_Zinc2jet",
                                              "hresponseJetsMassMidPU_Zinc2jet",
                                              nJetsMass_Zinc2jet,
                                              jetsMass_Zinc2jet,
                                              nJetsMass_Zinc2jet,
                                              jetsMass_Zinc2jet);

    JetsMassHigPU_Zinc2jet = newTH1D("JetsMassHigPU_Zinc2jet",
                                     "2Jets Invariant Mass HigPU (N_{jets} #geq 2)",
                                     Mjj,
                                     nJetsMass_Zinc2jet,
                                     jetsMass_Zinc2jet);
    genJetsMassHigPU_Zinc2jet = newTH1D("genJetsMassHigPU_Zinc2jet",
                                        "gen 2Jets Invariant Mass HigPU (N_{jets} #geq 2)",
                                        Mjj,
                                        nJetsMass_Zinc2jet,
                                        jetsMass_Zinc2jet);
    hresponseJetsMassHigPU_Zinc2jet = newTH2D("hresponseJetsMassHigPU_Zinc2jet",
                                              "hresponseJetsMassHigPU_Zinc2jet",
                                              nJetsMass_Zinc2jet,
                                              jetsMass_Zinc2jet,
                                              nJetsMass_Zinc2jet,
                                              jetsMass_Zinc2jet);

    BestJetsMass_Zinc2jet = newTH1D("BestJetsMass_Zinc2jet",
                                    "Best 2Jets Invariant Mass (N_{jets} #geq 2)",
                                    Mjj,
                                    nJetsMass_Zinc2jet,
                                    jetsMass_Zinc2jet);
    genBestJetsMass_Zinc2jet = newTH1D("genBestJetsMass_Zinc2jet",
                                       "gen Best 2Jets Invariant Mass (N_{jets} #geq 2)",
                                       Mjj,
                                       nJetsMass_Zinc2jet,
                                       jetsMass_Zinc2jet);
    ptBal_Zinc1jet = newTH1D("ptBal_Zinc1jet",
                             "Vectorial pT sum: Z_{pT} + 1Jet_{pT} (N_{jets} #geq 1)",
                             "#Sigma pT [GeV]",
                             50,
                             0,
                             100);
    ptBal_Zinc2jet = newTH1D("ptBal_Zinc2jet",
                             "Vectorial pT sum: Z_{pT} + 2Jet_{pT} (N_{jets} #geq 2)",
                             "#Sigma pT [GeV]",
                             50,
                             0,
                             100);
    ptBal_Zinc3jet = newTH1D("ptBal_Zinc3jet",
                             "Vectorial pT sum: Z_{pT} + 3Jet_{pT} (N_{jets} #geq 3)",
                             "#Sigma pT [GeV]",
                             50,
                             0,
                             100);
    genptBal_Zinc2jet = newTH1D("genptBal_Zinc2jet",
                                "gen Vectorial pT sum: Z_{pT} + DiJet_{pT} (N_{jets} #geq 2)",
                                "#Sigma pT [GeV]",
                                50,
                                0,
                                100);
    dPhiJets_Zinc2jet =
        newTH1D("dPhiJets_Zinc2jet", "#Delta#phi btwn jets (N_{jets} #geq 2)", jdPhi, 20, 0, PI);
    gendPhiJets_Zinc2jet = newTH1D(
        "gendPhiJets_Zinc2jet", "gen #Delta#phi btwn jets (N_{jets} #geq 2)", jdPhi, 20, 0, PI);
    BestdPhiJets_Zinc2jet = newTH1D(
        "BestdPhiJets_Zinc2jet", "Best #Delta#phi btwn jets (N_{jets} #geq 2)", jdPhi, 20, 0, PI);
    genBestdPhiJets_Zinc2jet = newTH1D("genBestdPhiJets_Zinc2jet",
                                       "gen Best #Delta#phi btwn jets (N_{jets} #geq 2)",
                                       jdPhi,
                                       20,
                                       0,
                                       PI);
    dEtaJets_Zinc2jet =
        newTH1D("dEtaJets_Zinc2jet", "#Delta#eta btwn jets (N_{jets} #geq 2)", jdEta, 48, 0, 4.8);
    gendEtaJets_Zinc2jet = newTH1D(
        "gendEtaJets_Zinc2jet", "gen #Delta#eta btwn jets (N_{jets} #geq 2)", jdEta, 48, 0, 4.8);
    dEtaFirstJetZ_Zinc2jet = newTH1D("dEtaFirstJetZ_Zinc2jet",
                                     "#Delta#eta btwn Jet_{1} and Z (N_{jets} #geq 2)",
                                     "#Delta#eta(j_{1}Z)",
                                     50,
                                     -6,
                                     6);
    gendEtaFirstJetZ_Zinc2jet = newTH1D("gendEtaFirstJetZ_Zinc2jet",
                                        "gen #Delta#eta btwn Jet_{1} and Z (N_{jets} #geq 2)",
                                        "#Delta#eta(j_{1}Z)",
                                        50,
                                        -6,
                                        6);
    dEtaSecondJetZ_Zinc2jet = newTH1D("dEtaSecondJetZ_Zinc2jet",
                                      "#Delta#eta btwn Jet_{2} and Z (N_{jets} #geq 2)",
                                      "#Delta#eta(j_{2}Z)",
                                      50,
                                      -6,
                                      6);
    gendEtaSecondJetZ_Zinc2jet = newTH1D("gendEtaSecondJetZ_Zinc2jet",
                                         "gen #Delta#eta btwn Jet_{2} and Z (N_{jets} #geq 2)",
                                         "#Delta#eta(j_{2}Z)",
                                         120,
                                         -6,
                                         6);
    dEtaJet1Plus2Z_Zinc2jet = newTH1D("dEtaJet1Plus2Z_Zinc2jet",
                                      "#Delta#eta btwn jets and Z (N_{jets} #geq 2)",
                                      "#Delta#eta(j_{12}Z)",
                                      120,
                                      -6,
                                      6);
    gendEtaJet1Plus2Z_Zinc2jet = newTH1D("gendEtaJet1Plus2Z_Zinc2jet",
                                         "gen #Delta#eta btwn jets and Z (N_{jets} #geq 2)",
                                         "#Delta#eta(j_{12}Z)",
                                         120,
                                         -6,
                                         6);
    PHI_Zinc2jet = newTH1D("PHI_Zinc2jet",
                           "#phi: Angle btwn the two subsystems planes (N_{jets} #geq 2)",
                           "#phi(j_{12}Z)",
                           25,
                           0,
                           PI);
    genPHI_Zinc2jet = newTH1D("genPHI_Zinc2jet",
                              "gen #phi: Angle btwn the two subsystems planes (N_{jets} #geq 2)",
                              "#phi(j_{12}Z)",
                              25,
                              0,
                              PI);
    BestPHI_Zinc2jet = newTH1D("BestPHI_Zinc2jet",
                               "Best #phi: Angle btwn the two subsystems planes (N_{jets} #geq 2)",
                               "#phi(j_{12}Z)",
                               25,
                               0,
                               PI);
    genBestPHI_Zinc2jet =
        newTH1D("genBestPHI_Zinc2jet",
                "gen Best #phi: Angle btwn the two subsystems planes (N_{jets} #geq 2)",
                "#phi(j_{12}Z)",
                25,
                0,
                PI);
    PHI_T_Zinc2jet = newTH1D("PHI_T_Zinc2jet",
                             "#Delta S Angle btwn lep and jet pair in T-plane (N_{jets} #geq 2)",
                             "#Delta S(j_{12}Z)",
                             10,
                             0,
                             PI);
    genPHI_T_Zinc2jet =
        newTH1D("genPHI_T_Zinc2jet",
                "gen #Delta S Angle btwn lep and jet pair in T-plane (N_{jets} #geq 2)",
                "#Delta S(j_{12}Z)",
                10,
                0,
                PI);
    BestPHI_T_Zinc2jet =
        newTH1D("BestPHI_T_Zinc2jet",
                "Best #Delta S Angle btwn lep and jet pair in T-plane (N_{jets} #geq 2)",
                "#Delta S(j_{12}Z)",
                10,
                0,
                PI);
    genBestPHI_T_Zinc2jet =
        newTH1D("genBestPHI_T_Zinc2jet",
                "gen Best #Delta S Angle btwn lep and jet pair in T-plane (N_{jets} #geq 2)",
                "#Delta S(j_{12}Z)",
                10,
                0,
                PI);
    SpT_Zinc2jet = newTH1D(
        "SpT_Zinc2jet", "#Delta_{pT}^{rel} lep and jets combined (N_{jets} #geq 2)", Spt, 20, 0, 1);
    genSpT_Zinc2jet = newTH1D("genSpT_Zinc2jet",
                              "gen #Delta_{pT}^{rel} lep and jets combined (N_{jets} #geq 2)",
                              Spt,
                              20,
                              0,
                              1);
    BestSpT_Zinc2jet = newTH1D("BestSpT_Zinc2jet",
                               "Best #Delta_{pT}^{rel} lep and jets combined (N_{jets} #geq 2)",
                               Spt,
                               20,
                               0,
                               1);
    genBestSpT_Zinc2jet =
        newTH1D("genBestSpT_Zinc2jet",
                "gen Best #Delta_{pT}^{rel} lep and jets combined (N_{jets} #geq 2)",
                Spt,
                20,
                0,
                1);
    SpTJets_Zinc2jet =
        newTH1D("SpTJets_Zinc2jet", "#Delta_{pT}^{rel} jets (N_{jets} #geq 2)", jSpt, 20, 0, 1);
    genSpTJets_Zinc2jet = newTH1D(
        "genSpTJets_Zinc2jet", "gen #Delta_{pT}^{rel} jets (N_{jets} #geq 2)", jSpt, 20, 0, 1);
    BestSpTJets_Zinc2jet = newTH1D(
        "BestSpTJets_Zinc2jet", "Best #Delta_{pT}^{rel} jets (N_{jets} #geq 2)", jSpt, 20, 0, 1);
    genBestSpTJets_Zinc2jet = newTH1D("genBestSpTJets_Zinc2jet",
                                      "gen Best #Delta_{pT}^{rel} jets (N_{jets} #geq 2)",
                                      jSpt,
                                      20,
                                      0,
                                      1);
    SPhi_Zinc2jet = newTH1D(
        "SPhi_Zinc2jet", "S_{#phi} lep and jets combined (N_{jets} #geq 2)", Sphi, 50, 0, PI);
    genSPhi_Zinc2jet = newTH1D("genSPhi_Zinc2jet",
                               "gen S_{#phi} lep and jets combined (N_{jets} #geq 2)",
                               Sphi,
                               50,
                               0,
                               PI);
    BestSPhi_Zinc2jet = newTH1D("BestSPhi_Zinc2jet",
                                "Best S_{#phi} lep and jets combined (N_{jets} #geq 2)",
                                Sphi,
                                50,
                                0,
                                PI);
    genBestSPhi_Zinc2jet = newTH1D("genBestSPhi_Zinc2jet",
                                   "gen Best S_{#phi} lep and jets combined (N_{jets} #geq 2)",
                                   Sphi,
                                   50,
                                   0,
                                   PI);

    //-- low Z pT
    TwoJetsPtDiff_LowPt_Zexc2jet =
        newTH1D("TwoJetsPtDiff_LowPt_Zexc2jet",
                "pT diff of the two highest jet at low Z_{pT} (N_{jets} = 2)",
                "#Delta pT(j_{1}j_{2}) [GeV]",
                10,
                0,
                100);
    genTwoJetsPtDiff_LowPt_Zexc2jet =
        newTH1D("genTwoJetsPtDiff_LowPt_Zexc2jet",
                "gen pT diff of the two highest jet at low Z_{pT} (N_{jets} = 2)",
                "#Delta pT(j_{1}j_{2}) [GeV]",
                10,
                0,
                100);
    JetsMass_LowPt_Zexc2jet = newTH1D("JetsMass_LowPt_Zexc2jet",
                                      "2Jets Invariant Mass at low Z_{pT} (N_{jets} = 2)",
                                      Mjj,
                                      25,
                                      0,
                                      625);
    genJetsMass_LowPt_Zexc2jet = newTH1D("genJetsMass_LowPt_Zexc2jet",
                                         "gen 2Jets Invariant Mass at low Z_{pT} (N_{jets} = 2)",
                                         Mjj,
                                         25,
                                         0,
                                         625);
    ptBal_LowPt_Zexc2jet =
        newTH1D("ptBal_LowPt_Zexc2jet",
                "Vectorial pT sum: Z_{pT} + DiJet_{pT} at low Z_{pT} (N_{jets} = 2)",
                "#Sigma pT [GeV]",
                50,
                0,
                100);
    genptBal_LowPt_Zexc2jet =
        newTH1D("genptBal_LowPt_Zexc2jet",
                "gen Vectorial pT sum: Z_{pT} + DiJet_{pT} at low Z_{pT} (N_{jets} = 2)",
                "#Sigma pT [GeV]",
                50,
                0,
                100);
    dPhiJets_LowPt_Zexc2jet = newTH1D("dPhiJets_LowPt_Zexc2jet",
                                      "#Delta#phi btwn jets at low Z_{pT} (N_{jets} = 2)",
                                      jdPhi,
                                      15,
                                      0,
                                      PI);
    gendPhiJets_LowPt_Zexc2jet = newTH1D("gendPhiJets_LowPt_Zexc2jet",
                                         "gen #Delta#phi btwn jets at low Z_{pT} (N_{jets} = 2)",
                                         jdPhi,
                                         15,
                                         0,
                                         PI);
    dPhiLeptons_LowPt_Zexc2jet = newTH1D("dPhiLeptons_LowPt_Zexc2jet",
                                         "#Delta#phi btwn leptons at low Z_{pT} (N_{jets} = 2)",
                                         ldPhi,
                                         50,
                                         0,
                                         PI);
    gendPhiLeptons_LowPt_Zexc2jet =
        newTH1D("gendPhiLeptons_LowPt_Zexc2jet",
                "gen #Delta#phi btwn leptons at low Z_{pT} (N_{jets} = 2)",
                ldPhi,
                50,
                0,
                PI);
    PHI_LowPt_Zexc2jet =
        newTH1D("PHI_LowPt_Zexc2jet",
                "#phi: Angle btwn the two subsystems planes at low Z_{pT} (N_{jets} = 2)",
                "#phi(j_{12}Z)",
                25,
                0,
                PI);
    genPHI_LowPt_Zexc2jet =
        newTH1D("genPHI_LowPt_Zexc2jet",
                "gen #phi: Angle btwn the two subsystems planes at low Z_{pT} (N_{jets} = 2)",
                "#phi(j_{12}Z)",
                25,
                0,
                PI);
    PHI_T_LowPt_Zexc2jet =
        newTH1D("PHI_T_LowPt_Zexc2jet",
                "#Delta S Angle btwn lepton and jet pair in T-plane at low Z_{pT} (N_{jets} = 2)",
                "#Delta S(j_{12}Z)",
                10,
                0,
                PI);
    genPHI_T_LowPt_Zexc2jet = newTH1D(
        "genPHI_T_LowPt_Zexc2jet",
        "gen #Delta S Angle btwn lepton and jet pair in T-plane at low Z_{pT} (N_{jets} = 2)",
        "#Delta S(j_{12}Z)",
        10,
        0,
        PI);
    SpT_LowPt_Zexc2jet =
        newTH1D("SpT_LowPt_Zexc2jet",
                "#Delta_{pT}^{rel} leptons and jets combined at low Z_{pT} (N_{jets} = 2)",
                Spt,
                25,
                0,
                1);
    genSpT_LowPt_Zexc2jet =
        newTH1D("genSpT_LowPt_Zexc2jet",
                "gen #Delta_{pT}^{rel} leptons and jets combined at low Z_{pT} (N_{jets} = 2)",
                Spt,
                25,
                0,
                1);
    SpTJets_LowPt_Zexc2jet = newTH1D("SpTJets_LowPt_Zexc2jet",
                                     "#Delta_{pT}^{rel} jets at low Z_{pT} (N_{jets} = 2)",
                                     jSpt,
                                     15,
                                     0,
                                     1);
    genSpTJets_LowPt_Zexc2jet = newTH1D("genSpTJets_LowPt_Zexc2jet",
                                        "gen #Delta_{pT}^{rel} jets at low Z_{pT} (N_{jets} = 2)",
                                        jSpt,
                                        15,
                                        0,
                                        1);
    SpTLeptons_LowPt_Zexc2jet = newTH1D("SpTLeptons_LowPt_Zexc2jet",
                                        "#Delta_{pT}^{rel} leptons at low Z_{pT} (N_{jets} = 2)",
                                        lSpt,
                                        50,
                                        0,
                                        1);
    genSpTLeptons_LowPt_Zexc2jet =
        newTH1D("genSpTLeptons_LowPt_Zexc2jet",
                "gen #Delta_{pT}^{rel} leptons at low Z_{pT} (N_{jets} = 2)",
                lSpt,
                50,
                0,
                1);
    SPhi_LowPt_Zexc2jet =
        newTH1D("SPhi_LowPt_Zexc2jet",
                "S_{#phi}: leptons and jets combined at low Z_{pT} (N_{jets} = 2)",
                Sphi,
                50,
                0,
                PI);
    genSPhi_LowPt_Zexc2jet =
        newTH1D("genSPhi_LowPt_Zexc2jet",
                "gen S_{#phi}: leptons and jets combined at low Z_{pT} (N_{jets} = 2)",
                Sphi,
                50,
                0,
                PI);

    TwoJetsPtDiff_LowPt_Zinc2jet =
        newTH1D("TwoJetsPtDiff_LowPt_Zinc2jet",
                "pT diff of the two highest jet at low Z_{pT} (N_{jets} #geq 2)",
                "#Delta pT(j_{1}j_{2}) [GeV]",
                10,
                0,
                100);
    genTwoJetsPtDiff_LowPt_Zinc2jet =
        newTH1D("genTwoJetsPtDiff_LowPt_Zinc2jet",
                "gen pT diff of the two highest jet at low Z_{pT}  (N_{jets} #geq 2)",
                "#Delta pT(j_{1}j_{2}) [GeV]",
                10,
                0,
                100);
    BestTwoJetsPtDiff_LowPt_Zinc2jet =
        newTH1D("BestTwoJetsPtDiff_LowPt_Zinc2jet",
                "Best pT diff of the two highest jet at low Z_{pT} (N_{jets} #geq 2)",
                "#Delta pT(j_{1}j_{2}) [GeV]",
                10,
                0,
                100);
    genBestTwoJetsPtDiff_LowPt_Zinc2jet =
        newTH1D("genBestTwoJetsPtDiff_LowPt_Zinc2jet",
                "gen Best pT diff of the two highest jet at low Z_{pT} (N_{jets} #geq 2)",
                "#Delta pT(j_{1}j_{2}) [GeV]",
                10,
                0,
                100);
    JetsMass_LowPt_Zinc2jet = newTH1D("JetsMass_LowPt_Zinc2jet",
                                      "2Jets Invariant Mass at low Z_{pT} (N_{jets} #geq 2)",
                                      Mjj,
                                      25,
                                      0,
                                      625);
    genJetsMass_LowPt_Zinc2jet = newTH1D("genJetsMass_LowPt_Zinc2jet",
                                         "gen 2Jets Invariant Mass at low Z_{pT} (N_{jets} #geq 2)",
                                         Mjj,
                                         25,
                                         0,
                                         625);
    BestJetsMass_LowPt_Zinc2jet =
        newTH1D("BestJetsMass_LowPt_Zinc2jet",
                "Best 2Jets Invariant Mass at low Z_{pT} (N_{jets} #geq 2)",
                Mjj,
                25,
                0,
                625);
    genBestJetsMass_LowPt_Zinc2jet =
        newTH1D("genBestJetsMass_LowPt_Zinc2jet",
                "gen Best 2Jets Invariant Mass at low Z_{pT} (N_{jets} #geq 2)",
                Mjj,
                25,
                0,
                625);
    ptBal_LowPt_Zinc2jet =
        newTH1D("ptBal_LowPt_Zinc2jet",
                "Vectorial pT sum: Z_{pT} + DiJet_{pT} at low Z_{pT} (N_{jets} #geq 2)",
                "#Sigma pT [GeV]",
                50,
                0,
                100);
    genptBal_LowPt_Zinc2jet =
        newTH1D("genptBal_LowPt_Zinc2jet",
                "gen Vectorial pT sum: Z_{pT} + DiJet_{pT} at low Z_{pT} (N_{jets} #geq 2)",
                "#Sigma pT [GeV]",
                50,
                0,
                100);
    dPhiJets_LowPt_Zinc2jet = newTH1D("dPhiJets_LowPt_Zinc2jet",
                                      "#Delta#phi btwn jets at low Z_{pT} (N_{jets} #geq 2)",
                                      jdPhi,
                                      15,
                                      0,
                                      PI);
    gendPhiJets_LowPt_Zinc2jet = newTH1D("gendPhiJets_LowPt_Zinc2jet",
                                         "gen#Delta#phi btwn jets at low Z_{pT} (N_{jets} #geq 2)",
                                         jdPhi,
                                         15,
                                         0,
                                         PI);
    BestdPhiJets_LowPt_Zinc2jet =
        newTH1D("BestdPhiJets_LowPt_Zinc2jet",
                "Best #Delta#phi btwn jets at low Z_{pT} (N_{jets} #geq 2)",
                jdPhi,
                15,
                0,
                PI);
    genBestdPhiJets_LowPt_Zinc2jet =
        newTH1D("genBestdPhiJets_LowPt_Zinc2jet",
                "gen Best #Delta#phi btwn jets at low Z_{pT} (N_{jets} #geq 2)",
                jdPhi,
                15,
                0,
                PI);
    dPhiLeptons_LowPt_Zinc2jet = newTH1D("dPhiLeptons_LowPt_Zinc2jet",
                                         "#Delta #phi btwn leptons at low Z_{pT} (N_{jets} #geq 2)",
                                         ldPhi,
                                         50,
                                         0,
                                         PI);
    gendPhiLeptons_LowPt_Zinc2jet =
        newTH1D("gendPhiLeptons_LowPt_Zinc2jet",
                "gen #Delta #phi btwn leptons at low Z_{pT} (N_{jets} #geq 2)",
                ldPhi,
                50,
                0,
                PI);
    PHI_LowPt_Zinc2jet =
        newTH1D("PHI_LowPt_Zinc2jet",
                "#phi: Angle btwn the two subsystems planes at low Z_{pT} (N_{jets} #geq 2)",
                "#phi(j_{12}Z)",
                25,
                0,
                PI);
    genPHI_LowPt_Zinc2jet =
        newTH1D("genPHI_LowPt_Zinc2jet",
                "gen #phi: Angle btwn the two subsystems planes at low Z_{pT} (N_{jets} #geq 2)",
                "#phi(j_{12}Z)",
                25,
                0,
                PI);
    BestPHI_LowPt_Zinc2jet =
        newTH1D("BestPHI_LowPt_Zinc2jet",
                "Best #phi: Angle btwn the two subsystems planes at low Z_{pT} (N_{jets} #geq 2)",
                "#phi(j_{12}Z)",
                25,
                0,
                PI);
    genBestPHI_LowPt_Zinc2jet = newTH1D(
        "genBestPHI_LowPt_Zinc2jet",
        "gen Best #phi: Angle btwn the two subsystems planes at low Z_{pT} (N_{jets} #geq 2)",
        "#phi(j_{12}Z)",
        25,
        0,
        PI);
    PHI_T_LowPt_Zinc2jet = newTH1D(
        "PHI_T_LowPt_Zinc2jet",
        "#Delta S Angle btwn lepton and jet pair in T-plane at low Z_{pT} (N_{jets} #geq 2)",
        "#Delta S(j_{12}Z)",
        10,
        0,
        PI);
    genPHI_T_LowPt_Zinc2jet = newTH1D(
        "genPHI_T_LowPt_Zinc2jet",
        "gen #Delta S Angle btwn lepton and jet pair in T-plane at low Z_{pT} (N_{jets} #geq 2)",
        "#Delta S(j_{12}Z)",
        10,
        0,
        PI);
    BestPHI_T_LowPt_Zinc2jet = newTH1D(
        "BestPHI_T_LowPt_Zinc2jet",
        "Best #Delta S Angle btwn lepton and jet pair in T-plane at low Z_{pT} (N_{jets} #geq 2)",
        "#Delta S(j_{12}Z)",
        10,
        0,
        PI);
    genBestPHI_T_LowPt_Zinc2jet = newTH1D("genBestPHI_T_LowPt_Zinc2jet",
                                          "gen Best #Delta S Angle btwn lepton and jet pair in "
                                          "T-plane at low Z_{pT} (N_{jets} #geq 2)",
                                          "#Delta S(j_{12}Z)",
                                          10,
                                          0,
                                          PI);
    SpT_LowPt_Zinc2jet =
        newTH1D("SpT_LowPt_Zinc2jet",
                "#Delta_{pT}^{rel} leptons and jets combined at low Z_{pT} (N_{jets} #geq 2)",
                Spt,
                25,
                0,
                1);
    genSpT_LowPt_Zinc2jet =
        newTH1D("genSpT_LowPt_Zinc2jet",
                "gen #Delta_{pT}^{rel} leptons and jets combined at low Z_{pT} (N_{jets} #geq 2)",
                Spt,
                25,
                0,
                1);
    BestSpT_LowPt_Zinc2jet =
        newTH1D("BestSpT_LowPt_Zinc2jet",
                "Best #Delta_{pT}^{rel} leptons and jets combined at low Z_{pT} (N_{jets} #geq 2)",
                Spt,
                25,
                0,
                1);
    genBestSpT_LowPt_Zinc2jet = newTH1D(
        "genBestSpT_LowPt_Zinc2jet",
        "gen Best #Delta_{pT}^{rel} leptons and jets combined at low Z_{pT} (N_{jets} #geq 2)",
        Spt,
        25,
        0,
        1);
    SpTJets_LowPt_Zinc2jet = newTH1D("SpTJets_LowPt_Zinc2jet",
                                     "#Delta_{pT}^{rel} jets at low Z_{pT} (N_{jets} #geq 2)",
                                     jSpt,
                                     15,
                                     0,
                                     1);
    genSpTJets_LowPt_Zinc2jet =
        newTH1D("genSpTJets_LowPt_Zinc2jet",
                "gen #Delta_{pT}^{rel} jets at low Z_{pT} (N_{jets} #geq 2)",
                jSpt,
                15,
                0,
                1);
    BestSpTJets_LowPt_Zinc2jet =
        newTH1D("BestSpTJets_LowPt_Zinc2jet",
                "Best #Delta_{pT}^{rel} jets at low Z_{pT} (N_{jets} #geq 2)",
                jSpt,
                15,
                0,
                1);
    genBestSpTJets_LowPt_Zinc2jet =
        newTH1D("genBestSpTJets_LowPt_Zinc2jet",
                "gen Best #Delta_{pT}^{rel} jets at low Z_{pT} (N_{jets} #geq 2)",
                jSpt,
                15,
                0,
                1);
    SpTLeptons_LowPt_Zinc2jet = newTH1D("SpTLeptons_LowPt_Zinc2jet",
                                        "#Delta_{pT}^{rel} leptons at low Z_{pT} (N_{jets} #geq 2)",
                                        lSpt,
                                        50,
                                        0,
                                        1);
    genSpTLeptons_LowPt_Zinc2jet =
        newTH1D("genSpTLeptons_LowPt_Zinc2jet",
                "gen #Delta_{pT}^{rel} leptons at low Z_{pT} (N_{jets} #geq 2)",
                lSpt,
                50,
                0,
                1);
    SPhi_LowPt_Zinc2jet =
        newTH1D("SPhi_LowPt_Zinc2jet",
                "S_{#phi}: leptons and jets combined at low Z_{pT} (N_{jets} #geq 2)",
                Sphi,
                50,
                0,
                PI);
    genSPhi_LowPt_Zinc2jet =
        newTH1D("genSPhi_LowPt_Zinc2jet",
                "gen S_{#phi}: leptons and jets combined at low Z_{pT} (N_{jets} #geq 2)",
                Sphi,
                50,
                0,
                PI);
    BestSPhi_LowPt_Zinc2jet =
        newTH1D("BestSPhi_LowPt_Zinc2jet",
                "Best S_{#phi}: leptons and jets combined at low Z_{pT} (N_{jets} #geq 2)",
                Sphi,
                50,
                0,
                PI);
    genBestSPhi_LowPt_Zinc2jet =
        newTH1D("genBestSPhi_LowPt_Zinc2jet",
                "gen Best S_{#phi}: leptons and jets combined at low Z_{pT} (N_{jets} #geq 2)",
                Sphi,
                50,
                0,
                PI);

    //-- low Z pT and low SpT
    PHI_LowSpT_LowPt_Zexc2jet = newTH1D("PHI_LowSpT_LowPt_Zexc2jet",
                                        "#phi: Angle btwn the two subsystems planes at low "
                                        "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} = 2)",
                                        "#phi",
                                        25,
                                        0.,
                                        PI);
    genPHI_LowSpT_LowPt_Zexc2jet = newTH1D("genPHI_LowSpT_LowPt_Zexc2jet",
                                           "gen #phi: Angle btwn the two subsystems planes at low "
                                           "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} = 2)",
                                           "#phi",
                                           25,
                                           0.,
                                           PI);
    SPhi_LowSpT_LowPt_Zexc2jet = newTH1D("SPhi_LowSpT_LowPt_Zexc2jet",
                                         "S_{#phi}: leptons and jets combined at low "
                                         "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} = 2)",
                                         "S_{#phi}",
                                         50,
                                         2.5,
                                         PI);
    genSPhi_LowSpT_LowPt_Zexc2jet = newTH1D("genSPhi_LowSpT_LowPt_Zexc2jet",
                                            "gen S_{#phi}: leptons and jets combined at low "
                                            "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} = 2)",
                                            "S_{#phi}",
                                            50,
                                            2.5,
                                            PI);

    PHI_LowSpT_LowPt_Zinc2jet = newTH1D("PHI_LowSpT_LowPt_Zinc2jet",
                                        "#phi: Angle btwn the two subsystems planes at low "
                                        "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} #geq 2)",
                                        "#phi",
                                        25,
                                        0.,
                                        PI);
    genPHI_LowSpT_LowPt_Zinc2jet = newTH1D("genPHI_LowSpT_LowPt_Zinc2jet",
                                           "gen #phi: Angle btwn the two subsystems planes at low "
                                           "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} #geq 2)",
                                           "#phi",
                                           25,
                                           0.,
                                           PI);
    SPhi_LowSpT_LowPt_Zinc2jet = newTH1D("SPhi_LowSpT_LowPt_Zinc2jet",
                                         "S_{#phi}: leptons and jets combined at low "
                                         "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} #geq 2)",
                                         "S_{#phi}",
                                         50,
                                         2.5,
                                         PI);
    genSPhi_LowSpT_LowPt_Zinc2jet = newTH1D("genSPhi_LowSpT_LowPt_Zinc2jet",
                                            "gen S_{#phi}: leptons and jets combined at low "
                                            "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} #geq 2)",
                                            "S_{#phi}",
                                            50,
                                            2.5,
                                            PI);

    //-- low Z pT and high SpT
    PHI_HighSpT_LowPt_Zexc2jet = newTH1D("PHI_HighSpT_LowPt_Zexc2jet",
                                         "#phi: Angle btwn the two subsystems planes at high "
                                         "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} = 2)",
                                         "#phi",
                                         50,
                                         0.,
                                         PI);
    genPHI_HighSpT_LowPt_Zexc2jet = newTH1D("genPHI_HighSpT_LowPt_Zexc2jet",
                                            "gen #phi: Angle btwn the two subsystems planes at "
                                            "high #Delta_{pT}^{rel} and low Z_{pT} (N_{jets} = 2)",
                                            "#phi",
                                            50,
                                            0.,
                                            PI);
    SPhi_HighSpT_LowPt_Zexc2jet = newTH1D("SPhi_HighSpT_LowPt_Zexc2jet",
                                          "S_{#phi}: leptons and jets combined at high "
                                          "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} = 2)",
                                          "S_{#phi}",
                                          50,
                                          0.,
                                          PI);
    genSPhi_HighSpT_LowPt_Zexc2jet = newTH1D("genSPhi_HighSpT_LowPt_Zexc2jet",
                                             "gen S_{#phi}: leptons and jets combined at high "
                                             "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} = 2)",
                                             "S_{#phi}",
                                             50,
                                             0.,
                                             PI);

    PHI_HighSpT_LowPt_Zinc2jet = newTH1D("PHI_HighSpT_LowPt_Zinc2jet",
                                         "#phi: Angle btwn the two subsystems planes at high "
                                         "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} #geq 2)",
                                         "#phi",
                                         50,
                                         0.,
                                         PI);
    genPHI_HighSpT_LowPt_Zinc2jet =
        newTH1D("genPHI_HighSpT_LowPt_Zinc2jet",
                "gen #phi: Angle btwn the two subsystems planes at high #Delta_{pT}^{rel} and low "
                "Z_{pT} (N_{jets} #geq 2)",
                "#phi",
                50,
                0.,
                PI);
    SPhi_HighSpT_LowPt_Zinc2jet = newTH1D("SPhi_HighSpT_LowPt_Zinc2jet",
                                          "S_{#phi}: leptons and jets combined at high "
                                          "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} #geq 2)",
                                          "S_{#phi}",
                                          50,
                                          0.,
                                          PI);
    genSPhi_HighSpT_LowPt_Zinc2jet = newTH1D("genSPhi_HighSpT_LowPt_Zinc2jet",
                                             "gen S_{#phi}: leptons and jets combined at high "
                                             "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} #geq 2)",
                                             "S_{#phi}",
                                             50,
                                             0.,
                                             PI);

    //-- low Z pT and low SPhi
    SpT_LowSPhi_LowPt_Zexc2jet = newTH1D(
        "SpT_LowSPhi_LowPt_Zexc2jet",
        "#Delta_{pT}^{rel} leptons and jets combined at low S_{#phi} and low Z_{pT} (N_{jets} = 2)",
        "#Delta_{pT}^{rel}",
        50,
        0.,
        1.);
    genSpT_LowSPhi_LowPt_Zexc2jet = newTH1D("genSpT_LowSPhi_LowPt_Zexc2jet",
                                            "gen #Delta_{pT}^{rel} leptons and jets combined at "
                                            "low S_{#phi} and low Z_{pT} (N_{jets} = 2)",
                                            "#Delta_{pT}^{rel}",
                                            50,
                                            0.,
                                            1.);

    SpT_LowSPhi_LowPt_Zinc2jet = newTH1D("SpT_LowSPhi_LowPt_Zinc2jet",
                                         "#Delta_{pT}^{rel} leptons and jets combined at low "
                                         "S_{#phi} and low Z_{pT} (N_{jets} #geq 2)",
                                         "#Delta_{pT}^{rel}",
                                         50,
                                         0.,
                                         1.);
    genSpT_LowSPhi_LowPt_Zinc2jet = newTH1D("genSpT_LowSPhi_LowPt_Zinc2jet",
                                            "gen #Delta_{pT}^{rel} leptons and jets combined at "
                                            "low S_{#phi} and low Z_{pT} (N_{jets} #geq 2)",
                                            "#Delta_{pT}^{rel}",
                                            50,
                                            0.,
                                            1.);

    //-- low Z pT and high SPhi
    SpT_HighSPhi_LowPt_Zexc2jet = newTH1D("SpT_HighSPhi_LowPt_Zexc2jet",
                                          "#Delta_{pT}^{rel} leptons and jets combined at high "
                                          "S_{#phi} and low Z_{pT} (N_{jets} = 2)",
                                          "#Delta_{pT}^{rel}",
                                          50,
                                          0.,
                                          1.);
    genSpT_HighSPhi_LowPt_Zexc2jet = newTH1D("genSpT_HighSPhi_LowPt_Zexc2jet",
                                             "gen #Delta_{pT}^{rel} leptons and jets combined at "
                                             "high S_{#phi} and low Z_{pT} (N_{jets} = 2)",
                                             "#Delta_{pT}^{rel}",
                                             50,
                                             0.,
                                             1.);

    SpT_HighSPhi_LowPt_Zinc2jet = newTH1D("SpT_HighSPhi_LowPt_Zinc2jet",
                                          "#Delta_{pT}^{rel} leptons and jets combined at high "
                                          "S_{#phi} and low Z_{pT} (N_{jets} #geq 2)",
                                          "#Delta_{pT}^{rel}",
                                          50,
                                          0.,
                                          1.);
    genSpT_HighSPhi_LowPt_Zinc2jet = newTH1D("genSpT_HighSPhi_LowPt_Zinc2jet",
                                             "gen #Delta_{pT}^{rel} leptons and jets combined at "
                                             "high S_{#phi} and low Z_{pT} (N_{jets} #geq 2)",
                                             "#Delta_{pT}^{rel}",
                                             50,
                                             0.,
                                             1.);

    //-- high Z pT
    ptBal_HighPt_Zexc2jet =
        newTH1D("ptBal_HighPt_Zexc2jet",
                "Vectorial pT sum: Z_{pT} + DiJet_{pT} at high Z_{pT} (N_{jets} = 2)",
                "#Sigma pT [GeV]",
                50,
                0.,
                100.);
    genptBal_HighPt_Zexc2jet =
        newTH1D("genptBal_HighPt_Zexc2jet",
                "gen Vectorial pT sum: Z_{pT} + DiJet_{pT} at high Z_{pT} (N_{jets} = 2)",
                "#Sigma pT [GeV]",
                50,
                0.,
                100.);
    dPhiJets_HighPt_Zexc2jet = newTH1D("dPhiJets_HighPt_Zexc2jet",
                                       "#Delta#phi btwn jets at high Z_{pT} (N_{jets} = 2)",
                                       jdPhi,
                                       15,
                                       0,
                                       PI);
    gendPhiJets_HighPt_Zexc2jet = newTH1D("gendPhiJets_HighPt_Zexc2jet",
                                          "gen #Delta#phi btwn jets at high Z_{pT} (N_{jets} = 2)",
                                          jdPhi,
                                          15,
                                          0,
                                          PI);
    dPhiLeptons_HighPt_Zexc2jet = newTH1D("dPhiLeptons_HighPt_Zexc2jet",
                                          "#Delta#phi btwn leptons at high Z_{pT} (N_{jets} = 2)",
                                          ldPhi,
                                          50,
                                          0.,
                                          PI);
    gendPhiLeptons_HighPt_Zexc2jet =
        newTH1D("gendPhiLeptons_HighPt_Zexc2jet",
                "gen #Delta#phi btwn leptons at high Z_{pT} (N_{jets} = 2)",
                ldPhi,
                50,
                0.,
                PI);
    PHI_HighPt_Zexc2jet =
        newTH1D("PHI_HighPt_Zexc2jet",
                "#phi: Angle btwn the two subsystems planes at high Z_{pT} (N_{jets} = 2)",
                "#phi",
                50,
                0.,
                PI);
    genPHI_HighPt_Zexc2jet =
        newTH1D("genPHI_HighPt_Zexc2jet",
                "gen #phi: Angle btwn the two subsystems planes at high Z_{pT} (N_{jets} = 2)",
                "#phi",
                50,
                0.,
                PI);
    PHI_T_HighPt_Zexc2jet =
        newTH1D("PHI_T_HighPt_Zexc2jet",
                "#Delta S Angle btwn lepton and jet pair in T-plane at high Z_{pT} (N_{jets} = 2)",
                "#Delta S",
                10,
                0.,
                PI);
    genPHI_T_HighPt_Zexc2jet = newTH1D(
        "genPHI_T_HighPt_Zexc2jet",
        "gen #Delta S Angle btwn lepton and jet pair in T-plane at high Z_{pT} (N_{jets} = 2)",
        "#Delta S",
        10,
        0.,
        PI);
    SpT_HighPt_Zexc2jet =
        newTH1D("SpT_HighPt_Zexc2jet",
                "#Delta_{pT}^{rel} leptons and jets combined at high Z_{pT} (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);
    genSpT_HighPt_Zexc2jet =
        newTH1D("genSpT_HighPt_Zexc2jet",
                "gen #Delta_{pT}^{rel} leptons and jets combined at high Z_{pT} (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);
    SpTJets_HighPt_Zexc2jet = newTH1D("SpTJets_HighPt_Zexc2jet",
                                      "#Delta_{pT}^{rel} jets at high Z_{pT} (N_{jets} = 2)",
                                      "#Delta_{pT}^{rel}",
                                      15,
                                      0.,
                                      1.);
    genSpTJets_HighPt_Zexc2jet = newTH1D("genSpTJets_HighPt_Zexc2jet",
                                         "gen #Delta_{pT}^{rel} jets at high Z_{pT} (N_{jets} = 2)",
                                         "#Delta_{pT}^{rel}",
                                         15,
                                         0.,
                                         1.);
    SpTLeptons_HighPt_Zexc2jet = newTH1D("SpTLeptons_HighPt_Zexc2jet",
                                         "#Delta_{pT}^{rel} leptons at high Z_{pT} (N_{jets} = 2)",
                                         "#Delta_{pT}^{rel}",
                                         50,
                                         0.,
                                         1.);
    genSpTLeptons_HighPt_Zexc2jet =
        newTH1D("genSpTLeptons_HighPt_Zexc2jet",
                "gen #Delta_{pT}^{rel} leptons at high Z_{pT} (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);
    SPhi_HighPt_Zexc2jet =
        newTH1D("SPhi_HighPt_Zexc2jet",
                "S_{#phi}: leptons and jets combined at high Z_{pT} (N_{jets} = 2)",
                "S_{#phi}",
                50,
                0.,
                PI);
    genSPhi_HighPt_Zexc2jet =
        newTH1D("genSPhi_HighPt_Zexc2jet",
                "gen S_{#phi}: leptons and jets combined at high Z_{pT} (N_{jets} = 2)",
                "S_{#phi}",
                50,
                0.,
                PI);

    ptBal_HighPt_Zinc2jet =
        newTH1D("ptBal_HighPt_Zinc2jet",
                "Vectorial pT sum: Z_{pT} + DiJet_{pT} at high Z_{pT} (N_{jets} #geq 2)",
                "#Sigma pT [GeV]",
                50,
                0.,
                100.);
    genptBal_HighPt_Zinc2jet =
        newTH1D("genptBal_HighPt_Zinc2jet",
                "gen Vectorial pT sum: Z_{pT} + DiJet_{pT} at high Z_{pT} (N_{jets} #geq 2)",
                "#Sigma pT [GeV]",
                50,
                0.,
                100.);
    dPhiJets_HighPt_Zinc2jet = newTH1D("dPhiJets_HighPt_Zinc2jet",
                                       "#Delta#phi btwn jets at high Z_{pT} (N_{jets} #geq 2)",
                                       jdPhi,
                                       15,
                                       0,
                                       PI);
    gendPhiJets_HighPt_Zinc2jet =
        newTH1D("gendPhiJets_HighPt_Zinc2jet",
                "gen #Delta#phi btwn jets at high Z_{pT} (N_{jets} #geq 2)",
                jdPhi,
                15,
                0,
                PI);
    dPhiLeptons_HighPt_Zinc2jet =
        newTH1D("dPhiLeptons_HighPt_Zinc2jet",
                "#Delta#phi btwn leptons at high Z_{pT (N_{jets} #geq 2)}",
                ldPhi,
                50,
                0.,
                PI);
    gendPhiLeptons_HighPt_Zinc2jet =
        newTH1D("gendPhiLeptons_HighPt_Zinc2jet",
                "gen #Delta#phi btwn leptons at high Z_{pT} (N_{jets} #geq 2)",
                ldPhi,
                50,
                0.,
                PI);
    PHI_HighPt_Zinc2jet =
        newTH1D("PHI_HighPt_Zinc2jet",
                "#phi: Angle btwn the two subsystems planes at high Z_{pT} (N_{jets} #geq 2)",
                "#phi",
                50,
                0.,
                PI);
    genPHI_HighPt_Zinc2jet =
        newTH1D("genPHI_HighPt_Zinc2jet",
                "gen #phi: Angle btwn the two subsystems planes at high Z_{pT} (N_{jets} #geq 2)",
                "#phi",
                50,
                0.,
                PI);
    PHI_T_HighPt_Zinc2jet = newTH1D(
        "PHI_T_HighPt_Zinc2jet",
        "#Delta S Angle btwn lepton and jet pair in T-plane at high Z_{pT} (N_{jets} #geq 2)",
        "#Delta S",
        10,
        0.,
        PI);
    genPHI_T_HighPt_Zinc2jet = newTH1D(
        "genPHI_T_HighPt_Zinc2jet",
        "gen#Delta S Angle btwn lepton and jet pair in T-plane at high Z_{pT} (N_{jets} #geq 2)",
        "#Delta S",
        10,
        0.,
        PI);
    SpT_HighPt_Zinc2jet =
        newTH1D("SpT_HighPt_Zinc2jet",
                "#Delta_{pT}^{rel} leptons and jets combined at high Z_{pT} (N_{jets} #geq 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);
    genSpT_HighPt_Zinc2jet =
        newTH1D("genSpT_HighPt_Zinc2jet",
                "gen #Delta_{pT}^{rel} leptons and jets combined at high Z_{pT} (N_{jets} #geq 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);
    SpTJets_HighPt_Zinc2jet = newTH1D("SpTJets_HighPt_Zinc2jet",
                                      "#Delta_{pT}^{rel} jets at high Z_{pT} (N_{jets} #geq 2)",
                                      "#Delta_{pT}^{rel}",
                                      15,
                                      0.,
                                      1.);
    genSpTJets_HighPt_Zinc2jet =
        newTH1D("genSpTJets_HighPt_Zinc2jet",
                "gen #Delta_{pT}^{rel} jets at high Z_{pT} (N_{jets} #geq 2)",
                "#Delta_{pT}^{rel}",
                15,
                0.,
                1.);
    SpTLeptons_HighPt_Zinc2jet =
        newTH1D("SpTLeptons_HighPt_Zinc2jet",
                "#Delta_{pT}^{rel} leptons at high Z_{pT} (N_{jets} #geq 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);
    genSpTLeptons_HighPt_Zinc2jet =
        newTH1D("genSpTLeptons_HighPt_Zinc2jet",
                "gen #Delta_{pT}^{rel} leptons at high Z_{pT} (N_{jets} #geq 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);
    SPhi_HighPt_Zinc2jet =
        newTH1D("SPhi_HighPt_Zinc2jet",
                "S_{#phi}: leptons and jets combined at high Z_{pT} (N_{jets} #geq 2)",
                "S_{#phi}",
                50,
                0.,
                PI);
    genSPhi_HighPt_Zinc2jet =
        newTH1D("genSPhi_HighPt_Zinc2jet",
                "gen S_{#phi}: leptons and jets combined at high Z_{pT} (N_{jets} #geq 2)",
                "S_{#phi}",
                50,
                0.,
                PI);

    //-- high Z pT and low SpT
    PHI_LowSpT_HighPt_Zexc2jet = newTH1D("PHI_LowSpT_HighPt_Zexc2jet",
                                         "#phi: Angle btwn the two subsystems planes at low "
                                         "#Delta_{pT}^{rel} and high Z_{pT} (N_{jets} = 2)",
                                         "#Phi",
                                         50,
                                         0.,
                                         PI);
    SPhi_LowSpT_HighPt_Zexc2jet = newTH1D("SPhi_LowSpT_HighPt_Zexc2jet",
                                          "S_{#phi}: leptons and jets combined at low "
                                          "#Delta_{pT}^{rel} and high Z_{pT} (N_{jets} = 2)",
                                          "S_{#phi}",
                                          50,
                                          2.5,
                                          PI);

    PHI_LowSpT_HighPt_Zinc2jet = newTH1D("PHI_LowSpT_HighPt_Zinc2jet",
                                         "#phi: Angle btwn the two subsystems planes at low "
                                         "#Delta_{pT}^{rel} and high Z_{pT} (N_{jets} #geq 2)",
                                         "#Phi",
                                         50,
                                         0.,
                                         PI);
    SPhi_LowSpT_HighPt_Zinc2jet = newTH1D("SPhi_LowSpT_HighPt_Zinc2jet",
                                          "S_{#phi}: leptons and jets combined at low "
                                          "#Delta_{pT}^{rel} and high Z_{pT} (N_{jets} #geq 2)",
                                          "S_{#phi}",
                                          50,
                                          2.5,
                                          PI);

    //-- high Z pT and high SpT
    PHI_HighSpT_HighPt_Zexc2jet = newTH1D("PHI_HighSpT_HighPt_Zexc2jet",
                                          "#phi: Angle btwn the two subsystems planes at high "
                                          "#Delta_{pT}^{rel} and high Z_{pT} (N_{jets} = 2)",
                                          "#phi",
                                          50,
                                          0.,
                                          PI);
    SPhi_HighSpT_HighPt_Zexc2jet = newTH1D("SPhiHighSpT_HighPt_Zexc2jet",
                                           "S_{#phi}: leptons and jets combined at high "
                                           "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} = 2)",
                                           "S_{#phi}",
                                           50,
                                           0.,
                                           PI);

    PHI_HighSpT_HighPt_Zinc2jet = newTH1D("PHI_HighSpT_HighPt_Zinc2jet",
                                          "#phi: Angle btwn the two subsystems planes at high "
                                          "#Delta_{pT}^{rel} and high Z_{pT} (N_{jets} #geq 2)",
                                          "#phi",
                                          50,
                                          0.,
                                          PI);
    SPhi_HighSpT_HighPt_Zinc2jet = newTH1D("SPhiHighSpT_HighPt_Zinc2jet",
                                           "S_{#phi}: leptons and jets combined at high "
                                           "#Delta_{pT}^{rel} and low Z_{pT} (N_{jets} #geq 2)",
                                           "S_{#phi}",
                                           50,
                                           0.,
                                           PI);

    //-- high Z pT and low SPhi
    SpT_LowSPhi_HighPt_Zexc2jet = newTH1D("SpT_LowSPhi_HighPt_Zexc2jet",
                                          "#Delta_{pT}^{rel} leptons and jets combined at low "
                                          "S_{#phi} and high Z_{pT} (N_{jets} = 2)",
                                          "#Delta_{pT}^{rel}",
                                          50,
                                          0.,
                                          1.);

    SpT_LowSPhi_HighPt_Zinc2jet = newTH1D("SpT_LowSPhi_HighPt_Zinc2jet",
                                          "#Delta_{pT}^{rel} leptons and jets combined at low "
                                          "S_{#phi} and high Z_{pT} (N_{jets} #geq 2)",
                                          "#Delta_{pT}^{rel}",
                                          50,
                                          0.,
                                          1.);

    //-- high Z pT and high SPhi
    SpT_HighSPhi_HighPt_Zexc2jet = newTH1D("SpT_HighSPhi_HighPt_Zexc2jet",
                                           "#Delta_{pT}^{rel} leptons and jets combined at high "
                                           "S_{#phi} and high Z_{pT} (N_{jets} = 2)",
                                           "#Delta_{pT}^{rel}",
                                           50,
                                           0.,
                                           1.);

    SpT_HighSPhi_HighPt_Zinc2jet = newTH1D("SpT_HighSPhi_HighPt_Zinc2jet",
                                           "#Delta_{pT}^{rel} leptons and jets combined at high "
                                           "S_{#phi} and high Z_{pT} (N_{jets} #geq 2)",
                                           "#Delta_{pT}^{rel}",
                                           50,
                                           0.,
                                           1.);

    //-- low SPhi
    SpT_LowSPhi_Zexc2jet =
        newTH1D("SpT_LowSPhi_Zexc2jet",
                "#Delta_{pT}^{rel} leptons and jets combined at low S_{#phi} (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);

    SpT_LowSPhi_Zinc2jet =
        newTH1D("SpT_LowSPhi_Zinc2jet",
                "#Delta_{pT}^{rel} leptons and jets combined at low S_{#phi} (N_{jets} #geq 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);

    //-- high SPhi
    SpT_HighSPhi_Zexc2jet =
        newTH1D("SpT_HighSPhi_Zexc2jet",
                "#Delta_{pT}^{rel} leptons and jets combined at high S_{#phi} (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);

    SpT_HighSPhi_Zinc2jet =
        newTH1D("SpT_HighSPhi_Zinc2jet",
                "#Delta_{pT}^{rel} leptons and jets combined at high S_{#phi} (N_{jets} #geq 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);

    //-- low SpT
    PHI_LowSpT_Zexc2jet = newTH1D(
        "PHI_LowSpT_Zexc2jet",
        "#phi: Angle btwn the two subsystems planes at low #Delta_{pT}^{rel} (N_{jets} = 2)",
        "#Phi",
        50,
        0.,
        PI);
    SPhi_LowSpT_Zexc2jet =
        newTH1D("SPhi_LowSpT_Zexc2jet",
                "S_{#phi}: leptons and jets combined at low #Delta_{pT}^{rel} (N_{jets} = 2)",
                "S_{#phi}",
                50,
                2.5,
                PI);

    PHI_LowSpT_Zinc2jet = newTH1D(
        "PHI_LowSpT_Zinc2jet",
        "#phi: Angle btwn the two subsystems planes at low #Delta_{pT}^{rel} (N_{jets} #geq 2)",
        "#Phi",
        50,
        0.,
        PI);
    SPhi_LowSpT_Zinc2jet =
        newTH1D("SPhi_LowSpT_Zinc2jet",
                "S_{#phi}: leptons and jets combined at low #Delta_{pT}^{rel} (N_{jets} #geq 2)",
                "S_{#phi}",
                50,
                2.5,
                PI);

    //-- high SpT
    PHI_HighSpT_Zexc2jet = newTH1D(
        "PHI_HighSpT_Zexc2jet",
        "#phi: Angle btwn the two subsystems planes at high #Delta_{pT}^{rel} (N_{jets} = 2)",
        "#Phi",
        50,
        0.,
        PI);
    SPhi_HighSpT_Zexc2jet =
        newTH1D("SPhi_HighSpT_Zexc2jet",
                "S_{#phi}: leptons and jets combined at high #Delta_{pT}^{rel} (N_{jets} = 2)",
                "S_{#phi}",
                50,
                0.,
                PI);

    PHI_HighSpT_Zinc2jet = newTH1D(
        "PHI_HighSpT_Zinc2jet",
        "#phi: Angle btwn the two subsystems planes at high #Delta_{pT}^{rel} (N_{jets} #geq 2)",
        "#Phi",
        50,
        0.,
        PI);
    SPhi_HighSpT_Zinc2jet =
        newTH1D("SPhi_HighSpT_Zinc2jet",
                "S_{#phi}: leptons and jets combined at high #Delta_{pT}^{rel} (N_{jets} #geq 2)",
                "S_{#phi}",
                50,
                0.,
                PI);

    //-- gen stuff
    gendPhiJetsDeltaR_Zexc2jet =
        newTH1D("gendPhiJetsDeltaR_Zexc2jet",
                "#Delta #phi btwn gen jets with #Delta R < 0.5 (N_{jets} = 2)",
                "#Delta#phi",
                50,
                0.,
                PI);
    resdPhiJetsDeltaR_Zexc2jet =
        newTH1D("resdPhiJetsDeltaR_Zexc2jet",
                "#Delta #phi btwn gen jets with #Delta R < 0.5 (N_{jets} = 2)",
                "#Delta#phi",
                50,
                -2.5,
                2.5);
    genPHI_TDeltaR_Zexc2jet = newTH1D("genPHI_TDeltaR_Zexc2jet",
                                      "#Delta S Angle btwn gen lep and gen jet pair in T-plane "
                                      "with #Delta R < 0.5 (N_{jets} = 2)",
                                      "#Delta S",
                                      50,
                                      0.,
                                      PI);
    resPHI_TDeltaR_Zexc2jet = newTH1D("resPHI_TDeltaR_Zexc2jet",
                                      "#Delta S Angle btwn gen lep and gen jet pair in T-plane "
                                      "with #Delta R < 0.5 (N_{jets} = 2)",
                                      "#Delta S",
                                      50,
                                      -2.5,
                                      2.5);
    genSpTJetsDeltaR_Zexc2jet =
        newTH1D("genSpTJetsDeltaR_Zexc2jet",
                "#Delta_{pT}^{rel} Gen jets with #Delta R < 0.5 (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                50,
                0.,
                1.);
    resSpTJetsDeltaR_Zexc2jet =
        newTH1D("resSpTJetsDeltaR_Zexc2jet",
                "#Delta_{pT}^{rel} Gen jets with #Delta R < 0.5 (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                50,
                -2.5,
                2.5);
    genSpTDeltaR_Zexc2jet = newTH1D("genSpTDeltaR_Zexc2jet",
                                    "#Delta_{pT}^{rel} with #Delta R < 0.5 (N_{jets} = 2)",
                                    "#Delta_{pT}^{rel}",
                                    50,
                                    0.,
                                    1.);
    resSpTDeltaR_Zexc2jet = newTH1D("resSpTDeltaR_Zexc2jet",
                                    "#Delta_{pT}^{rel} with #Delta R < 0.5 (N_{jets} = 2)",
                                    "#Delta_{pT}^{rel}",
                                    50,
                                    -2.5,
                                    2.5);

    gendPhiJetsDPS_Zexc2jet =
        newTH1D("gendPhiJetsDPS_Zexc2jet",
                "#Delta #phi btwn gen jets matching DPS parton (N_{jets} = 2)",
                "#Delta#phi_{j_{1}j_{2}}",
                50,
                0.,
                PI);
    gendPhiJetsDPSDeltaR_Zexc2jet =
        newTH1D("gendPhiJetsDPSDeltaR_Zexc2jet",
                "#Delta #phi btwn gen jets matching DPS parton with #Delta R < 0.5 (N_{jets} = 2)",
                "#Delta#phi",
                50,
                0.,
                PI);
    genPHI_TDPS_Zexc2jet =
        newTH1D("genPHI_TDPS_Zexc2jet",
                "#Delta S Angle btwn gen lepton and jet pair in T-plane (N_{jets} = 2)",
                "#Delta S",
                50,
                0.,
                PI);
    genPHI_TDPSDeltaR_Zexc2jet = newTH1D(
        "genPHI_TDPSDeltaR_Zexc2jet",
        "#Delta S Angle btwn gen lepton and jet pair in T-plane with #Delta R < 0.5 (N_{jets} = 2)",
        "#Delta S",
        50,
        0.,
        PI);
    genSpTJetsDPS_Zexc2jet =
        newTH1D("genSpTJetsDPS_Zexc2jet",
                "#Delta_{pT}^{rel} Gen jets matching DPS parton (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                nbinSpt,
                binSpt);
    genSpTJetsDPSDeltaR_Zexc2jet =
        newTH1D("genSpTJetsDPSDeltaR_Zexc2jet",
                "#Delta_{pT}^{rel} Gen jets matching DPS parton with #Delta R < 0.5 (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                nbinSpt,
                binSpt);
    genSpTDPS_Zexc2jet =
        newTH1D("genSpTDPS_Zexc2jet",
                "#Delta_{pT}^{rel} with gen jets matching DPS parton (N_{jets} = 2)",
                "#Delta_{pT}^{rel}",
                nbinSpt,
                binSpt);
    genSpTDPSDeltaR_Zexc2jet = newTH1D(
        "genSpTDPSDeltaR_Zexc2jet",
        "#Delta_{pT}^{rel} with gen jets matching DPS parton with #Delta R < 0.5 (N_{jets} = 2)",
        "#Delta_{pT}^{rel}",
        nbinSpt,
        binSpt);
    genSpTDPSPartons_Zexc2jet = newTH1D("genSpTDPSPartons_Zexc2jet",
                                        "#Delta_{pT}^{rel} DPS partons (N_{jets} = 2)",
                                        "#Delta_{pT}^{rel}",
                                        nbinSpt,
                                        binSpt);

    genZNGoodJets_Zinc =
        newTH1D("genZNGoodJets_Zinc", "Jet Counter (incl.)", "N_{jets}", 7, -0.5, 6.5);
    if (genZNGoodJets_Zinc) {
        for (int ibin = 1; ibin < genZNGoodJets_Zinc->GetNbinsX(); ++ibin) {
            genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(ibin, TString::Format("#geq %d", ibin));
        }

        //	genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(1,"#geq 0");
        //	genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(2,"#geq 1");
        //	genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(3,"#geq 2");
        //	genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(4,"#geq 3");
        //	genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(5,"#geq 4");
        //	genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(6,"#geq 5");
        //	genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(7,"#geq 6");
        //    //genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(8,"#geq 7");
        //    //genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(9,"#geq 8");
        //    //genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(10,"#geq 9");
        //    //genZNGoodJets_Zinc->GetXaxis()->SetBinLabel(11,"#geq 10");
    }
    if (doWJets)
        genZNGoodJets_Zexc =
            newTH1D("genZNGoodJets_Zexc", "Jet Counter (excl.)", "N_{jets}", 11, -0.5, 10.5);
    else
        genZNGoodJets_Zexc =
            newTH1D("genZNGoodJets_Zexc", "Jet Counter (excl.)", "N_{jets}", 7, -0.5, 6.5);

    if (genZNGoodJets_Zexc) {
        for (int ibin = 1; ibin < genZNGoodJets_Zexc->GetNbinsX(); ++ibin) {
            genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(ibin, TString::Format("= %d", ibin - 1));
        }

        //    genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(1,"= 0");
        //    genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(2,"= 1");
        //    genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(3,"= 2");
        //    genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(4,"= 3");
        //    genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(5,"= 4");
        //    genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(6,"= 5");
        //    genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(7,"= 6");
        //    if ( doWJets ){
        //        genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(8,"= 7");
        //        genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(9,"#geq 8");
        //        genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(10,"#geq 9");
        //        genZNGoodJets_Zexc->GetXaxis()->SetBinLabel(11,"#geq 10");
        //    }
    }
    // Correlations

    gendPhiJetsDPSDeltaR_ZpT_Zexc2jet = newTH2D("gendPhiJetsDPSDeltaR_ZpT_Zexc2jet",
                                                "gendPhiJetsDPSDeltaR_ZpT_Zexc2jet",
                                                50,
                                                0,
                                                PI,
                                                100,
                                                0,
                                                100);
    partonX2D = newTH2D("partonX2D", "parton X: x1 vs x2", 100, 0.0001, 0.2, 100, 0.0001, 0.2);

    gendeltaRjetMu = newTH1D("gendeltaRjetMu", "gen delta R btwn jet and muon", "#R", 50, 0., 2.5);

    lepResolution_pt = newTH2D("lepResolution_pt",
                               "muon pt resolution",
                               40,
                               0,
                               200, // x-axis
                               40,
                               0.,
                               0.1); // y-axis

    lepResolution_nvtx = newTH2D("lepResolution_nvtx",
                                 "muon pt resolution vs nvtx",
                                 55,
                                 0.5,
                                 55.5, // x-axis
                                 40,
                                 0.,
                                 0.1); // y-axis

    lepResolution_pt_rel = newTH2D("lepResolution_pt_rel",
                                   "muon pt resolution rel",
                                   40,
                                   0,
                                   200, // x-axis
                                   40,
                                   -0.1,
                                   0.1); // y-axis

    deltaRMuRecGen_lead_lowM = newTH1D("deltaRMuRecGen_lead_lowM",
                                       "gen delta R btwn reco and gen leading muon low M",
                                       "#R",
                                       40,
                                       0.,
                                       0.2);
    deltaRMuRecGen_sublead_lowM = newTH1D("deltaRMuRecGen_sublead_lowM",
                                          "gen delta R btwn reco and gen subleading muon low M",
                                          "#R",
                                          40,
                                          0.,
                                          0.2);

    deltaRMuRecGen_lead_highM = newTH1D("deltaRMuRecGen_lead_highM",
                                        "gen delta R btwn reco and gen leading muon high M",
                                        "#R",
                                        40,
                                        0.,
                                        0.2);
    deltaRMuRecGen_sublead_highM = newTH1D("deltaRMuRecGen_sublead_highM",
                                           "gen delta R btwn reco and gen subleading muon high M",
                                           "#R",
                                           40,
                                           0.,
                                           0.2);

    //   deltaPtMuRecGen_lead_highM                = newTH1D("deltaPtMuRecGen_lead_highM", "gen
    //   delta Pt btwn reco and gen leading muon high M", "#R", 100, 0., 30);
    //   deltaPtMuRecGen_sublead_highM             = newTH1D("deltaPtMuRecGen_sublead_highM", "gen
    //   delta Pt btwn reco and gen subleading muon high M", "#R", 100, 0., 30);

    //   deltaPtMuRecGen_lead_lowM                = newTH1D("deltaPtMuRecGen_lead_lowM", "gen delta
    //   Pt btwn reco and gen leading muon low M", "#R", 100, 0., 30);
    //   deltaPtMuRecGen_sublead_lowM             = newTH1D("deltaPtMuRecGen_sublead_lowM", "gen
    //   delta Pt btwn reco and gen subleading muon low M", "#R", 100, 0., 30);

    deltaPtMuRecGen_lead_highM = newTH1D("deltaPtMuRecGen_lead_highM",
                                         "gen delta Pt btwn reco and gen leading muon high M",
                                         "#R",
                                         50,
                                         -0.05,
                                         0.05);
    deltaPtMuRecGen_sublead_highM = newTH1D("deltaPtMuRecGen_sublead_highM",
                                            "gen delta Pt btwn reco and gen subleading muon high M",
                                            "#R",
                                            50,
                                            -0.05,
                                            0.05);

    deltaPtMuRecGen_lead_lowM = newTH1D("deltaPtMuRecGen_lead_lowM",
                                        "gen delta Pt btwn reco and gen leading muon low M",
                                        "#R",
                                        50,
                                        -0.05,
                                        0.05);
    deltaPtMuRecGen_sublead_lowM = newTH1D("deltaPtMuRecGen_sublead_lowM",
                                           "gen delta Pt btwn reco and gen subleading muon low M",
                                           "#R",
                                           50,
                                           -0.05,
                                           0.05);

    /// additional information
    // Muoisolation

    MuDetIsoRhoCorr =
        newTH1D("MuDetIsoRhoCorr", "Muon Detect. Iso #rho corr.", "l_{Iso}^{Det.}", 30, 0, 1.5);
    MuPFIsoDBetaCorr =
        newTH1D("MuPFIsoDBetaCorr", "Muon PF Iso DBeta corr.", "l_{Iso}^{PF}", 30, 0, 1.5);

    MuPFIsoDBetaCorrj1 =
        newTH1D("MuPFIsoDBetaCorrj1", "Muon PF Iso DBeta corr j1.", "l_{Iso}^{PF}", 30, 0, 1.5);
    MuPFIsoDBetaCorrj2 =
        newTH1D("MuPFIsoDBetaCorrj2", "Muon PF Iso DBeta corr. j2", "l_{Iso}^{PF}", 30, 0, 1.5);
    MuPFIsoDBetaCorrj3 =
        newTH1D("MuPFIsoDBetaCorrj3", "Muon PF Iso DBeta corr. j3", "l_{Iso}^{PF}", 30, 0, 1.5);

    deltaRjetMu = newTH1D("deltaRjetMu", "delta R btwn jet and muon", "#R", 50, 0., 2.5);
    deltaPtjetMu =
        newTH1D("deltaPtjetMu", "delta Pt btwn jet and muon if dR<0.5", "#R", 150, -75., 75.);

    // TH2D* jecVspt=newTH1D("jecVspt","jec Vs pt","jec","pt",80,0.,400,100,0,0.5);
    NVtx_Zinc0jet = newTH1D("NVtx_Zinc0jet", "Number of vertices 0 jet inc", "#Vtx", 45, 0.5, 45.5);
    NVtx_NoPUweight_Zinc0jet =
        newTH1D("NVtx_NoPUweight_Zinc0jet", "Number of vertices", "#Vtx", 45, 0.5, 45.5);
    TruePU_0 = newTH1D("TruePU_0", "True pile-up 0 jet", "#pu", 45, 0.5, 45.5);
    TruePU_1 = newTH1D("TruePU_1", "True pile-up 1 jet", "#pu", 45, 0.5, 45.5);
    TruePU_2 = newTH1D("TruePU_2", "True pile-up 2 jets", "#pu", 45, 0.5, 45.5);
    TruePU_3 = newTH1D("TruePU_3", "True pile-up 3 jets", "#pu", 45, 0.5, 45.5);
    TruePU_4 = newTH1D("TruePU_4", "True pile-up 4 jets", "#pu", 45, 0.5, 45.5);
    TruePU_5 = newTH1D("TruePU_5", "True pile-up 5 jets", "#pu", 45, 0.5, 45.5);
    TruePU_6 = newTH1D("TruePU_6", "True pile-up 6 jets", "#pu", 45, 0.5, 45.5);
    TruePU_7 = newTH1D("TruePU_7", "True pile-up 7 jets", "#pu", 45, 0.5, 45.5);
    NVtx_Zexc0jet =
        newTH1D("NVtx_Zexc0jet", "Number of vertices 0 jet exc", "#NVtx", 45, 0.5, 45.5);
    NVtx_Zexc1jet =
        newTH1D("NVtx_Zexc1jet", "Number of vertices 1 jet exc", "#NVtx", 45, 0.5, 45.5);
    NVtx_Zexc2jet =
        newTH1D("NVtx_Zexc2jet", "Number of vertices 2 jet exc", "#NVtx", 45, 0.5, 45.5);
    NVtx_Zexc3jet =
        newTH1D("NVtx_Zexc3jet", "Number of vertices 3 jet exc", "#NVtx", 45, 0.5, 45.5);
    NVtx_Zexc4jet =
        newTH1D("NVtx_Zexc4jet", "Number of vertices 4 jet exc", "#NVtx", 45, 0.5, 45.5);
    NVtx_Zexc5jet =
        newTH1D("NVtx_Zexc5jet", "Number of vertices 5 jet exc", "#NVtx", 45, 0.5, 45.5);
    NVtx_Zexc6jet =
        newTH1D("NVtx_Zexc6jet", "Number of vertices 6 jet exc", "#NVtx", 45, 0.5, 45.5);
    NVtx_Zexc7jet =
        newTH1D("NVtx_Zexc7jet", "Number of vertices 7 jet exc", "#NVtx", 45, 0.5, 45.5);

    ZNGoodJetsBeta_Zexc = newTH2D(
        "ZNGoodJetsBeta_Zexc", "Beta cut vs Jet Counter (excl.) ", 11, -0.5, 10.5, 10, -0.5, 9.5);
    if (ZNGoodJetsBeta_Zexc) {
        for (int ibin = 1; ibin < ZNGoodJetsBeta_Zexc->GetNbinsX(); ++ibin) {
            ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(ibin, TString::Format("= %d", ibin - 1));
        }
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(1, "= 0");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(2, "= 1");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(3, "= 2");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(4, "= 3");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(5, "= 4");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(6, "= 5");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(7, "= 6");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(8, "= 7");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(9, "= 8");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(10,"= 9");
        //    ZNGoodJetsBeta_Zexc->GetXaxis()->SetBinLabel(11,"= 10");
    }

    Beta = newTH1D("Beta", "Jet PU variable Beta", "Beta", 50, 0., 1.);
    BetaStar = newTH1D("BetaStar", "Jet PU variable BetaStar", "BetaStar", 50, 0., 1.);
    puBeta_JetsMatchGenJets =
        newTH1D("puBeta_JetsMatchGenJets", "puBeta_JetsMatchGenJets", "Beta", 50, 0, 1);
    puBetaStar_JetsMatchGenJets =
        newTH1D("puBetaStar_JetsMatchGenJets", "puBetaStar_JetsMatchGenJets", "Beta", 50, 0, 1);
    puBeta_JetsNoMatchGenJets =
        newTH1D("puBeta_JetsNoMatchGenJets", "puBeta_JetsNoMatchGenJets", "Beta", 50, 0, 1);
    puBetaStar_JetsNoMatchGenJets =
        newTH1D("puBetaStar_JetsNoMatchGenJets", "puBetaStar_JetsNoMatchGenJets", "Beta", 50, 0, 1);
    puMVA = newTH1D("puMVA", "Jet PU variable from MVA", "puMVA", 40, -1., 1.);
    puMVA_JetsMatchGenJets = newTH1D("puMVA_JetsMatchGenJets",
                                     "Jet PU variable from MVA for matching jets",
                                     "puMVA",
                                     40,
                                     -1.,
                                     1.);
    puMVA_JetsNoMatchGenJets = newTH1D("puMVA_JetsNoMatchGenJets",
                                       "Jet PU variable from MVA for non matching jets",
                                       "puMVA",
                                       40,
                                       -1.,
                                       1.);
    jetsEta_JetsMatchGenJets =
        newTH1D("jetsEta_JetsMatchGenJets", "Jet Eta for matching jets", "puMVA", 48, -2.4, 2.4);
    jetsEta_JetsNoMatchGenJets = newTH1D(
        "jetsEta_JetsNoMatchGenJets", "Jet Eta for non matching jets", "puMVA", 48, -2.4, 2.4);
    puMVAvsBeta =
        newTH2D("puMVA vs beta", "Jet PU variable from MVA vs Beta", 50, -1., 1., 50, 0., 1.);

    PUWeight = newTH1D("PUWeight", "PU weight Z all", "PU weight Z all", 500, 0., 14.);
    PUWeight0 = newTH1D("PUWeight0", "PU weight Z+0jet", "PU weight Z+0jet", 500, 0., 14.);
    PUWeight1 = newTH1D("PUWeigh1", "PU weight Z+jet>0 ", "PU weight Z+jet>0", 500, 0., 14.);

    partonsN = newTH1D("partonsN", "Number of ME partons ", "N_{partons}", 16, -0.5, 15.5);
    partonsNWeighted = newTH1D(
        "partonsNWeighted", "Number of ME partons: weighted ", "N_{partons}", 16, -0.5, 15.5);
    partonsNAfterGenCut = newTH1D("partonsNAfterGenCut",
                                  "Number of ME partons passing the gen cut",
                                  "N_{partons}",
                                  16,
                                  -0.5,
                                  15.5);
    partonsNAfterGenCutWeighted = newTH1D("partonsNAfterGenCutWeighted",
                                          "Number of ME partons passing the gen cut:weighted",
                                          "N_{partons}",
                                          16,
                                          -0.5,
                                          15.5);

    // vector boson and single jet
    dEtaBosonJet_Zexc1jet = newTH1D("dEtaBosonJet_Zexc1",
                                    "#Delta#eta btwn leading jet and V (N_{jets} #eq 1)",
                                    lJetdEta,
                                    72,
                                    0,
                                    4.8);
    gendEtaBosonJet_Zexc1jet = newTH1D("gendEtaBosonJet_Zexc1",
                                       "gen #Delta#eta btwn leading jet and V (N_{jets} #eq 1)",
                                       lJetdEta,
                                       72,
                                       0,
                                       4.8);
    dEtaBosonJet_Zinc1jet = newTH1D("dEtaBosonJet_Zinc1",
                                    "#Delta#eta btwn leading jet and V (N_{jets} #geq 1)",
                                    lJetdEta,
                                    72,
                                    0,
                                    4.8);
    gendEtaBosonJet_Zinc1jet = newTH1D("gendEtaBosonJet_Zinc1",
                                       "gen #Delta#eta btwn leading jet and V (N_{jets} #geq 1)",
                                       lJetdEta,
                                       72,
                                       0,
                                       4.8);

    // Additional sum and difference of Z+jet rapidity
    AbsZRapidity_Zinc1jet =
        newTH1D("AbsZRapidity_Zinc1jet", "AbsZRapidity_Zinc1jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_Zinc1jet =
        newTH1D("genAbsZRapidity_Zinc1jet", "genAbsZRapidity_Zinc1jet", "|y_{Z}|", 12, 0, 2.4);

    AbsFirstJetRapidity_Zinc1jet = newTH1D(
        "AbsFirstJetRapidity_Zinc1jet", "AbsFirstJetRapidity_Zinc1jet", "|y_{jet1}|", 12, 0, 2.4);
    genAbsFirstJetRapidity_Zinc1jet = newTH1D("genAbsFirstJetRapidity_Zinc1jet",
                                              "genAbsFirstJetRapidity_Zinc1jet",
                                              "|y_{jet1}|",
                                              12,
                                              0,
                                              2.4);
    SumZFirstJetRapidity_Zinc1jet = newTH1D("SumZFirstJetRapidity_Zinc1jet",
                                            "SumZFirstJetRapidity_Zinc1jet",
                                            "y_{sum(Z,jet1)}",
                                            12,
                                            0,
                                            2.4);
    genSumZFirstJetRapidity_Zinc1jet = newTH1D("genSumZFirstJetRapidity_Zinc1jet",
                                               "genSumZFirstJetRapidity_Zinc1jet",
                                               "y_{sum(Z,jet1)}",
                                               12,
                                               0,
                                               2.4);
    DifZFirstJetRapidity_Zinc1jet = newTH1D("DifZFirstJetRapidity_Zinc1jet",
                                            "DifZFirstJetRapidity_Zinc1jet",
                                            "y_{diff(Z,jet1)}",
                                            12,
                                            0,
                                            2.4);
    genDifZFirstJetRapidity_Zinc1jet = newTH1D("genDifZFirstJetRapidity_Zinc1jet",
                                               "genDifZFirstJetRapidity_Zinc1jet",
                                               "y_{diff(Z,jet1)}",
                                               12,
                                               0,
                                               2.4);

    /// Cross check//////
    SumZFirstJetEta_Zinc1jet =
        newTH1D("SumZFirstJetEta_Zinc1jet", "SumZFirstJetEta_Zinc1jet", "#eta_{sum}", 12, 0, 2.4);
    genSumZFirstJetEta_Zinc1jet = newTH1D(
        "genSumZFirstJetEta_Zinc1jet", "genSumZFirstJetEta_Zinc1jet", "#eta_{sum}", 12, 0, 2.4);
    DifZFirstJetEta_Zinc1jet =
        newTH1D("DifZFirstJetEta_Zinc1jet", "DifZFirstJetEta_Zinc1jet", "#eta_{diff}", 12, 0, 2.4);
    genDifZFirstJetEta_Zinc1jet = newTH1D(
        "genDifZFirstJetEta_Zinc1jet", "genDifZFirstJetEta_Zinc1jet", "#eta_{diff}", 12, 0, 2.4);

    ////Azimuth correlation cross check/////////////////////////////////////
    DPhiZFirstJet_Zinc1jet = newTH1D(
        "DPhiZFirstJet_Zinc1jet", "DPhiZFirstJet_Zinc1jet", "#Delta#phi_{Z,Jet1}", 25, 0, 3.14);
    genDPhiZFirstJet_Zinc1jet = newTH1D("genDPhiZFirstJet_Zinc1jet",
                                        "genDPhiZFirstJet_Zinc1jet",
                                        "#Delta#phi_{Z,Jet1}",
                                        25,
                                        0,
                                        3.14);

    AbsZRapidity_Zexc1jet =
        newTH1D("AbsZRapidity_Zexc1jet", "AbsZRapidity_Zexc1jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_Zexc1jet =
        newTH1D("genAbsZRapidity_Zexc1jet", "genAbsZRapidity_Zexc1jet", "|y_{Z}|", 12, 0, 2.4);
    AbsJetRapidity_Zexc1jet =
        newTH1D("AbsJetRapidity_Zexc1jet", "AbsJetRapidity_Zexc1jet", "|y_{jet}|", 12, 0, 2.4);
    genAbsJetRapidity_Zexc1jet = newTH1D(
        "genAbsJetRapidity_Zexc1jet", "genAbsJetRapidity_Zexc1jet", "|y_{jet}|", 12, 0, 2.4);
    SumZJetRapidity_Zexc1jet =
        newTH1D("SumZJetRapidity_Zexc1jet", "SumZJetRapidity_Zexc1jet", "y_{sum}", 12, 0, 2.4);
    genSumZJetRapidity_Zexc1jet = newTH1D(
        "genSumZJetRapidity_Zexc1jet", "genSumZJetRapidity_Zexc1jet", "y_{sum}", 12, 0, 2.4);
    DifZJetRapidity_Zexc1jet =
        newTH1D("DifZJetRapidity_Zexc1jet", "DifZJetRapidity_Zexc1jet", "y_{diff}", 12, 0, 2.4);
    genDifZJetRapidity_Zexc1jet = newTH1D(
        "genDifZJetRapidity_Zexc1jet", "genDifZJetRapidity_Zexc1jet", "y_{diff}", 12, 0, 2.4);

    AbsFirstJetRapidity_Zinc2jet = newTH1D(
        "AbsFirstJetRapidity_Zinc2jet", "AbsFirstJetRapidity_Zinc2jet", "|y_{jet1}|", 12, 0, 2.4);
    genAbsFirstJetRapidity_Zinc2jet = newTH1D("genAbsFirstJetRapidity_Zinc2jet",
                                              "genAbsFirstJetRapidity_Zinc2jet",
                                              "|y_{jet1}|",
                                              12,
                                              0,
                                              2.4);
    SumZFirstJetRapidity_Zinc2jet = newTH1D("SumZFirstJetRapidity_Zinc2jet",
                                            "SumZFirstJetRapidity_Zinc2jet",
                                            "y_{sum(Z,jet1)}",
                                            12,
                                            0,
                                            2.4);
    genSumZFirstJetRapidity_Zinc2jet = newTH1D("genSumZFirstJetRapidity_Zinc2jet",
                                               "genSumZFirstJetRapidity_Zinc2jet",
                                               "y_{sum(Z,jet1)}",
                                               12,
                                               0,
                                               2.4);
    DifZFirstJetRapidity_Zinc2jet = newTH1D("DifZFirstJetRapidity_Zinc2jet",
                                            "DifZFirstJetRapidity_Zinc2jet",
                                            "y_{diff(Z,jet1)}",
                                            12,
                                            0,
                                            2.4);
    genDifZFirstJetRapidity_Zinc2jet = newTH1D("genDifZFirstJetRapidity_Zinc2jet",
                                               "genDifZFirstJetRapidity_Zinc2jet",
                                               "y_{diff(Z,jet1)}",
                                               12,
                                               0,
                                               2.4);

    AbsZRapidity_Zinc2jet =
        newTH1D("AbsZRapidity_Zinc2jet", "AbsZRapidity_Zinc2jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_Zinc2jet =
        newTH1D("genAbsZRapidity_Zinc2jet", "genAbsZRapidity_Zinc2jet", "|y_{Z}|", 12, 0, 2.4);

    AbsSecondJetRapidity_Zinc2jet = newTH1D(
        "AbsSecondJetRapidity_Zinc2jet", "AbsSecondJetRapidity_Zinc2jet", "|y_{jet2}|", 12, 0, 2.4);
    genAbsSecondJetRapidity_Zinc2jet = newTH1D("genAbsSecondJetRapidity_Zinc2jet",
                                               "genAbsSecondJetRapidity_Zinc2jet",
                                               "|y_{jet2}|",
                                               12,
                                               0,
                                               2.4);
    SumZSecondJetRapidity_Zinc2jet = newTH1D("SumZSecondJetRapidity_Zinc2jet",
                                             "SumZSecondJetRapidity_Zinc2jet",
                                             "y_{sum(Z,jet2)}",
                                             12,
                                             0,
                                             2.4);
    genSumZSecondJetRapidity_Zinc2jet = newTH1D("genSumZSecondJetRapidity_Zinc2jet",
                                                "genSumZSecondJetRapidity_Zinc2jet",
                                                "y_{sum(Z,jet2)}",
                                                12,
                                                0,
                                                2.4);
    DifZSecondJetRapidity_Zinc2jet = newTH1D("DifZSecondJetRapidity_Zinc2jet",
                                             "DifZSecondJetRapidity_Zinc2jet",
                                             "y_{diff(Z,jet2)}",
                                             12,
                                             0,
                                             2.4);
    genDifZSecondJetRapidity_Zinc2jet = newTH1D("genDifZSecondJetRapidity_Zinc2jet",
                                                "genDifZSecondJetRapidity_Zinc2jet",
                                                "y_{diff(Z,jet2)}",
                                                12,
                                                0,
                                                2.4);

    SumFirstSecondJetRapidity_Zinc2jet = newTH1D("SumFirstSecondJetRapidity_Zinc2jet",
                                                 "SumFirstSecondJetRapidity_Zinc2jet",
                                                 "y_{sum(jet1,jet2)}",
                                                 12,
                                                 0,
                                                 2.4);
    genSumFirstSecondJetRapidity_Zinc2jet = newTH1D("genSumFirstSecondJetRapidity_Zinc2jet",
                                                    "genSumFirstSecondJetRapidity_Zinc2jet",
                                                    "y_{sum(jet1,jet2)}",
                                                    12,
                                                    0,
                                                    2.4);
    DifFirstSecondJetRapidity_Zinc2jet = newTH1D("DifFirstSecondJetRapidity_Zinc2jet",
                                                 "DifFirstSecondJetRapidity_Zinc2jet",
                                                 "y_{diff(jet1,jet2)}",
                                                 12,
                                                 0,
                                                 2.4);
    genDifFirstSecondJetRapidity_Zinc2jet = newTH1D("genDifFirstSecondJetRapidity_Zinc2jet",
                                                    "genDifFirstSecondJetRapidity_Zinc2jet",
                                                    "y_{diff(jet1,jet2)}",
                                                    12,
                                                    0,
                                                    2.4);

    SumZTwoJetsRapidity_Zinc2jet = newTH1D("SumZTwoJetsRapidity_Zinc2jet",
                                           "SumZTwoJetsRapidity_Zinc2jet",
                                           "y_{sum(Z,dijet)}",
                                           12,
                                           0,
                                           2.4);
    genSumZTwoJetsRapidity_Zinc2jet = newTH1D("genSumZTwoJetsRapidity_Zinc2jet",
                                              "genSumZTwoJetsRapidity_Zinc2jet",
                                              "y_{sum(Z,dijet)}",
                                              12,
                                              0,
                                              2.4);
    DifZTwoJetsRapidity_Zinc2jet = newTH1D("DifZTwoJetsRapidity_Zinc2jet",
                                           "DifZTwoJetsRapidity_Zinc2jet",
                                           "y_{diff(Z,dijet)}",
                                           12,
                                           0,
                                           2.4);
    genDifZTwoJetsRapidity_Zinc2jet = newTH1D("genDifZTwoJetsRapidity_Zinc2jet",
                                              "genDifZTwoJetsRapidity_Zinc2jet",
                                              "y_{diff(Z,dijet)}",
                                              12,
                                              0,
                                              2.4);

    ////////Azimuth cross check///////////////////////////////////////
    DPhiZFirstJet_Zinc2jet = newTH1D(
        "DPhiZFirstJet_Zinc2jet", "DPhiZFirstJet_Zinc2jet", "#Delta#phi_{Z,Jet1}", 25, 0, 3.14);
    genDPhiZFirstJet_Zinc2jet = newTH1D("genDPhiZFirstJet_Zinc2jet",
                                        "genDPhiZFirstJet_Zinc2jet",
                                        "#Delta#phi_{Z,Jet1}",
                                        25,
                                        0,
                                        3.14);
    DPhiZSecondJet_Zinc2jet = newTH1D(
        "DPhiZSecondJet_Zinc2jet", "DPhiZSecondJet_Zinc2jet", "#Delta#phi_{Z,Jet2}", 25, 0, 3.14);
    genDPhiZSecondJet_Zinc2jet = newTH1D("genDPhiZSecondJet_Zinc2jet",
                                         "genDPhiZSecondJet_Zinc2jet",
                                         "#Delta#phi_{Z,Jet2}",
                                         25,
                                         0,
                                         3.14);
    DPhiFirstSecondJet_Zinc2jet = newTH1D("DPhiFirstSecondJet_Zinc2jet",
                                          "DPhiFirstSecondJet_Zinc2jet",
                                          "#Delta#phi_{Jet1,Jet2}",
                                          25,
                                          0,
                                          3.14);
    genDPhiFirstSecondJet_Zinc2jet = newTH1D("genDPhiFirstSecondJet_Zinc2jet",
                                             "genDPhiFirstSecondJet_Zinc2jet",
                                             "#Delta#phi_{Jet1,Jet2}",
                                             25,
                                             0,
                                             3.14);

    AbsZRapidity_Zexc2jet =
        newTH1D("AbsZRapidity_Zexc2jet", "AbsZRapidity_Zexc2jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_Zexc2jet =
        newTH1D("genAbsZRapidity_Zexc2jet", "genAbsZRapidity_Zexc2jet", "|y_{Z}|", 12, 0, 2.4);
    AbsSecondJetRapidity_Zexc2jet = newTH1D(
        "AbsSecondJetRapidity_Zexc2jet", "AbsSecondJetRapidity_Zexc2jet", "|y_{jet}|", 12, 0, 2.4);
    genAbsSecondJetRapidity_Zexc2jet = newTH1D("genAbsSecondJetRapidity_Zexc2jet",
                                               "genAbsSecondJetRapidity_Zexc2jet",
                                               "|y_{jet}|",
                                               12,
                                               0,
                                               2.4);
    SumZSecondJetRapidity_Zexc2jet = newTH1D(
        "SumZSecondJetRapidity_Zexc2jet", "SumZSecondJetRapidity_Zexc2jet", "y_{sum}", 12, 0, 2.4);
    genSumZSecondJetRapidity_Zexc2jet = newTH1D("genSumZSecondJetRapidity_Zexc2jet",
                                                "genSumZSecondJetRapidity_Zexc2jet",
                                                "y_{sum}",
                                                12,
                                                0,
                                                2.4);
    DifZSecondJetRapidity_Zexc2jet = newTH1D(
        "DifZSecondJetRapidity_Zexc2jet", "DifZSecondJetRapidity_Zexc2jet", "y_{diff}", 12, 0, 2.4);
    genDifZSecondJetRapidity_Zexc2jet = newTH1D("genDifZSecondJetRapidity_Zexc2jet",
                                                "genDifZSecondJetRapidity_Zexc2jet",
                                                "y_{diff}",
                                                12,
                                                0,
                                                2.4);

    /////////Azimuth cross check//////////////////////////
    DPhiZFirstJet_Zinc3jet = newTH1D(
        "DPhiZFirstJet_Zinc3jet", "DPhiZFirstJet_Zinc3jet", "#Delta#phi_{Z,Jet1}", 25, 0, 3.14);
    genDPhiZFirstJet_Zinc3jet = newTH1D("genDPhiZFirstJet_Zinc3jet",
                                        "genDPhiZFirstJet_Zinc3jet",
                                        "#Delta#phi_{Z,Jet1}",
                                        25,
                                        0,
                                        3.14);
    DPhiZSecondJet_Zinc3jet = newTH1D(
        "DPhiZSecondJet_Zinc3jet", "DPhiZSecondJet_Zinc3jet", "#Delta#phi_{Z,Jet2}", 25, 0, 3.14);
    genDPhiZSecondJet_Zinc3jet = newTH1D("genDPhiZSecondJet_Zinc3jet",
                                         "genDPhiZSecondJet_Zinc3jet",
                                         "#Delta#phi_{Z,Jet2}",
                                         25,
                                         0,
                                         3.14);
    DPhiZThirdJet_Zinc3jet = newTH1D(
        "DPhiZThirdJet_Zinc3jet", "DPhiZThirdJet_Zinc3jet", "#Delta#phi_{Z,Jet3}", 25, 0, 3.14);
    genDPhiZThirdJet_Zinc3jet = newTH1D("genDPhiZThirdJet_Zinc3jet",
                                        "genDPhiZThirdJet_Zinc3jet",
                                        "#Delta#phi_{Z,Jet3}",
                                        25,
                                        0,
                                        3.14);
    DPhiFirstSecondJet_Zinc3jet = newTH1D("DPhiFirstSecondJet_Zinc3jet",
                                          "DPhiFirstSecondJet_Zinc3jet",
                                          "#Delta#phi_{Jet1,Jet2}",
                                          25,
                                          0,
                                          3.14);
    genDPhiFirstSecondJet_Zinc3jet = newTH1D("genDPhiFirstSecondJet_Zinc3jet",
                                             "genDPhiFirstSecondJet_Zinc3jet",
                                             "#Delta#phi_{Jet1,Jet2}",
                                             25,
                                             0,
                                             3.14);
    DPhiFirstThirdJet_Zinc3jet = newTH1D("DPhiFirstThirdJet_Zinc3jet",
                                         "DPhiFirstThirdJet_Zinc3jet",
                                         "#Delta#phi_{Jet1,Jet3}",
                                         25,
                                         0,
                                         3.14);
    genDPhiFirstThirdJet_Zinc3jet = newTH1D("genDPhiFirstThirdJet_Zinc3jet",
                                            "genDPhiFirstThirdJet_Zinc3jet",
                                            "#Delta#phi_{Jet1,Jet3}",
                                            25,
                                            0,
                                            3.14);
    DPhiSecondThirdJet_Zinc3jet = newTH1D("DPhiSecondThirdJet_Zinc3jet",
                                          "DPhiSecondThirdJet_Zinc3jet",
                                          "#Delta#phi_{Jet2,Jet3}",
                                          25,
                                          0,
                                          3.14);
    genDPhiSecondThirdJet_Zinc3jet = newTH1D("genDPhiSecondThirdJet_Zinc3jet",
                                             "genDPhiSecondThirdJet_Zinc3jet",
                                             "#Delta#phi_{Jet2,Jet3}",
                                             25,
                                             0,
                                             3.14);

    /////////Different Z boson Pt cuts/////////////////////////////
    AbsZRapidity_ZPt100_Zinc1jet = newTH1D(
        "AbsZRapidity_ZPt100_Zinc1jet", "AbsZRapidity_ZPt100_Zinc1jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_ZPt100_Zinc1jet = newTH1D("genAbsZRapidity_ZPt100_Zinc1jet",
                                              "genAbsZRapidity_ZPt100_Zinc1jet",
                                              "|y_{Z}|",
                                              12,
                                              0,
                                              2.4);
    AbsFirstJetRapidity_ZPt100_Zinc1jet = newTH1D("AbsFirstJetRapidity_ZPt100_Zinc1jet",
                                                  "AbsFirstJetRapidity_ZPt100_Zinc1jet",
                                                  "|y_{jet}|",
                                                  12,
                                                  0,
                                                  2.4);
    genAbsFirstJetRapidity_ZPt100_Zinc1jet = newTH1D("genAbsFirstJetRapidity_ZPt100_Zinc1jet",
                                                     "genAbsFirstJetRapidity_ZPt100_Zinc1jet",
                                                     "|y_{jet}|",
                                                     12,
                                                     0,
                                                     2.4);
    SumZFirstJetRapidity_ZPt100_Zinc1jet = newTH1D("SumZFirstJetRapidity_ZPt100_Zinc1jet",
                                                   "SumZFirstJetRapidity_ZPt100_Zinc1jet",
                                                   "y_{sum}",
                                                   12,
                                                   0,
                                                   2.4);
    genSumZFirstJetRapidity_ZPt100_Zinc1jet = newTH1D("genSumZFirstJetRapidity_ZPt100_Zinc1jet",
                                                      "genSumZFirstJetRapidity_ZPt100_Zinc1jet",
                                                      "y_{sum}",
                                                      12,
                                                      0,
                                                      2.4);
    DifZFirstJetRapidity_ZPt100_Zinc1jet = newTH1D("DifZFirstJetRapidity_ZPt100_Zinc1jet",
                                                   "DifZFirstJetRapidity_ZPt100_Zinc1jet",
                                                   "y_{diff}",
                                                   12,
                                                   0,
                                                   2.4);
    genDifZFirstJetRapidity_ZPt100_Zinc1jet = newTH1D("genDifZFirstJetRapidity_ZPt100_Zinc1jet",
                                                      "genDifZFirstJetRapidity_ZPt100_Zinc1jet",
                                                      "y_{diff}",
                                                      12,
                                                      0,
                                                      2.4);

    AbsZRapidity_ZPt100_Zexc1jet = newTH1D(
        "AbsZRapidity_ZPt100_Zexc1jet", "AbsZRapidity_ZPt100_Zexc1jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_ZPt100_Zexc1jet = newTH1D("genAbsZRapidity_ZPt100_Zexc1jet",
                                              "genAbsZRapidity_ZPt100_Zexc1jet",
                                              "|y_{Z}|",
                                              12,
                                              0,
                                              2.4);
    AbsJetRapidity_ZPt100_Zexc1jet = newTH1D("AbsJetRapidity_ZPt100_Zexc1jet",
                                             "AbsJetRapidity_ZPt100_Zexc1jet",
                                             "|y_{jet}|",
                                             12,
                                             0,
                                             2.4);
    genAbsJetRapidity_ZPt100_Zexc1jet = newTH1D("genAbsJetRapidity_ZPt100_Zexc1jet",
                                                "genAbsJetRapidity_ZPt100_Zexc1jet",
                                                "|y_{jet}|",
                                                12,
                                                0,
                                                2.4);
    SumZJetRapidity_ZPt100_Zexc1jet = newTH1D("SumZJetRapidity_ZPt100_Zexc1jet",
                                              "SumZJetRapidity_ZPt100_Zexc1jet",
                                              "y_{sum}",
                                              12,
                                              0,
                                              2.4);
    genSumZJetRapidity_ZPt100_Zexc1jet = newTH1D("genSumZJetRapidity_ZPt100_Zexc1jet",
                                                 "genSumZJetRapidity_ZPt100_Zexc1jet",
                                                 "y_{sum}",
                                                 12,
                                                 0,
                                                 2.4);
    DifZJetRapidity_ZPt100_Zexc1jet = newTH1D("DifZJetRapidity_ZPt100_Zexc1jet",
                                              "DifZJetRapidity_ZPt100_Zexc1jet",
                                              "y_{diff}",
                                              12,
                                              0,
                                              2.4);
    genDifZJetRapidity_ZPt100_Zexc1jet = newTH1D("genDifZJetRapidity_ZPt100_Zexc1jet",
                                                 "genDifZJetRapidity_ZPt100_Zexc1jet",
                                                 "y_{diff}",
                                                 12,
                                                 0,
                                                 2.4);

    AbsZRapidity_ZPt100_Zinc2jet = newTH1D(
        "AbsZRapidity_ZPt100_Zinc2jet", "AbsZRapidity_ZPt100_Zinc2jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_ZPt100_Zinc2jet = newTH1D("genAbsZRapidity_ZPt100_Zinc2jet",
                                              "genAbsZRapidity_ZPt100_Zinc2jet",
                                              "|y_{Z}|",
                                              12,
                                              0,
                                              2.4);
    AbsSecondJetRapidity_ZPt100_Zinc2jet = newTH1D("AbsSecondJetRapidity_ZPt100_Zinc2jet",
                                                   "AbsSecondJetRapidity_ZPt100_Zinc2jet",
                                                   "|y_{jet}|",
                                                   12,
                                                   0,
                                                   2.4);
    genAbsSecondJetRapidity_ZPt100_Zinc2jet = newTH1D("genAbsSecondJetRapidity_ZPt100_Zinc2jet",
                                                      "genAbsSecondJetRapidity_ZPt100_Zinc2jet",
                                                      "|y_{jet}|",
                                                      12,
                                                      0,
                                                      2.4);
    SumZSecondJetRapidity_ZPt100_Zinc2jet = newTH1D("SumZSecondJetRapidity_ZPt100_Zinc2jet",
                                                    "SumZSecondJetRapidity_ZPt100_Zinc2jet",
                                                    "y_{sum}",
                                                    12,
                                                    0,
                                                    2.4);
    genSumZSecondJetRapidity_ZPt100_Zinc2jet = newTH1D("genSumZSecondJetRapidity_ZPt100_Zinc2jet",
                                                       "genSumZSecondJetRapidity_ZPt100_Zinc2jet",
                                                       "y_{sum}",
                                                       12,
                                                       0,
                                                       2.4);
    DifZSecondJetRapidity_ZPt100_Zinc2jet = newTH1D("DifZSecondJetRapidity_ZPt100_Zinc2jet",
                                                    "DifZSecondJetRapidity_ZPt100_Zinc2jet",
                                                    "y_{diff}",
                                                    12,
                                                    0,
                                                    2.4);
    genDifZSecondJetRapidity_ZPt100_Zinc2jet = newTH1D("genDifZSecondJetRapidity_ZPt100_Zinc2jet",
                                                       "genDifZSecondJetRapidity_ZPt100_Zinc2jet",
                                                       "y_{diff}",
                                                       12,
                                                       0,
                                                       2.4);

    AbsZRapidity_ZPt100_Zexc2jet = newTH1D(
        "AbsZRapidity_ZPt100_Zexc2jet", "AbsZRapidity_ZPt100_Zexc2jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_ZPt100_Zexc2jet = newTH1D("genAbsZRapidity_ZPt100_Zexc2jet",
                                              "genAbsZRapidity_ZPt100_Zexc2jet",
                                              "|y_{Z}|",
                                              12,
                                              0,
                                              2.4);
    AbsSecondJetRapidity_ZPt100_Zexc2jet = newTH1D("AbsSecondJetRapidity_ZPt100_Zexc2jet",
                                                   "AbsSecondJetRapidity_ZPt100_Zexc2jet",
                                                   "|y_{jet}|",
                                                   12,
                                                   0,
                                                   2.4);
    genAbsSecondJetRapidity_ZPt100_Zexc2jet = newTH1D("genAbsSecondJetRapidity_ZPt100_Zexc2jet",
                                                      "genAbsSecondJetRapidity_ZPt100_Zexc2jet",
                                                      "|y_{jet}|",
                                                      12,
                                                      0,
                                                      2.4);
    SumZSecondJetRapidity_ZPt100_Zexc2jet = newTH1D("SumZSecondJetRapidity_ZPt100_Zexc2jet",
                                                    "SumZSecondJetRapidity_ZPt100_Zexc2jet",
                                                    "y_{sum}",
                                                    12,
                                                    0,
                                                    2.4);
    genSumZSecondJetRapidity_ZPt100_Zexc2jet = newTH1D("genSumZSecondJetRapidity_ZPt100_Zexc2jet",
                                                       "genSumZSecondJetRapidity_ZPt100_Zexc2jet",
                                                       "y_{sum}",
                                                       12,
                                                       0,
                                                       2.4);
    DifZSecondJetRapidity_ZPt100_Zexc2jet = newTH1D("DifZSecondJetRapidity_ZPt100_Zexc2jet",
                                                    "DifZSecondJetRapidity_ZPt100_Zexc2jet",
                                                    "y_{diff}",
                                                    12,
                                                    0,
                                                    2.4);
    genDifZSecondJetRapidity_ZPt100_Zexc2jet = newTH1D("genDifZSecondJetRapidity_ZPt100_Zexc2jet",
                                                       "genDifZSecondJetRapidity_ZPt100_Zexc2jet",
                                                       "y_{diff}",
                                                       12,
                                                       0,
                                                       2.4);

    AbsZRapidity_ZPt150_Zinc1jet = newTH1D(
        "AbsZRapidity_ZPt150_Zinc1jet", "AbsZRapidity_ZPt150_Zinc1jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_ZPt150_Zinc1jet = newTH1D("genAbsZRapidity_ZPt150_Zinc1jet",
                                              "genAbsZRapidity_ZPt150_Zinc1jet",
                                              "|y_{Z}|",
                                              12,
                                              0,
                                              2.4);
    AbsFirstJetRapidity_ZPt150_Zinc1jet = newTH1D("AbsFirstJetRapidity_ZPt150_Zinc1jet",
                                                  "AbsFirstJetRapidity_ZPt150_Zinc1jet",
                                                  "|y_{jet1}|",
                                                  12,
                                                  0,
                                                  2.4);
    genAbsFirstJetRapidity_ZPt150_Zinc1jet = newTH1D("genAbsFirstJetRapidity_ZPt150_Zinc1jet",
                                                     "genAbsFirstJetRapidity_ZPt150_Zinc1jet",
                                                     "|y_{jet1}|",
                                                     12,
                                                     0,
                                                     2.4);
    SumZFirstJetRapidity_ZPt150_Zinc1jet = newTH1D("SumZFirstJetRapidity_ZPt150_Zinc1jet",
                                                   "SumZFirstJetRapidity_ZPt150_Zinc1jet",
                                                   "y_{sum(Z,jet1)}",
                                                   12,
                                                   0,
                                                   2.4);
    genSumZFirstJetRapidity_ZPt150_Zinc1jet = newTH1D("genSumZFirstJetRapidity_ZPt150_Zinc1jet",
                                                      "genSumZFirstJetRapidity_ZPt150_Zinc1jet",
                                                      "y_{sum(Z,jet1)}",
                                                      12,
                                                      0,
                                                      2.4);
    DifZFirstJetRapidity_ZPt150_Zinc1jet = newTH1D("DifZFirstJetRapidity_ZPt150_Zinc1jet",
                                                   "DifZFirstJetRapidity_ZPt150_Zinc1jet",
                                                   "y_{diff(Z,jet1)}",
                                                   12,
                                                   0,
                                                   2.4);
    genDifZFirstJetRapidity_ZPt150_Zinc1jet = newTH1D("genDifZFirstJetRapidity_ZPt150_Zinc1jet",
                                                      "genDifZFirstJetRapidity_ZPt150_Zinc1jet",
                                                      "y_{diff(Z,jet1)}",
                                                      12,
                                                      0,
                                                      2.4);

    AbsZRapidity_ZPt150_Zexc1jet = newTH1D(
        "AbsZRapidity_ZPt150_Zexc1jet", "AbsZRapidity_ZPt150_Zexc1jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_ZPt150_Zexc1jet = newTH1D("genAbsZRapidity_ZPt150_Zexc1jet",
                                              "genAbsZRapidity_ZPt150_Zexc1jet",
                                              "|y_{Z}|",
                                              12,
                                              0,
                                              2.4);
    AbsJetRapidity_ZPt150_Zexc1jet = newTH1D("AbsJetRapidity_ZPt150_Zexc1jet",
                                             "AbsJetRapidity_ZPt150_Zexc1jet",
                                             "|y_{jet}|",
                                             12,
                                             0,
                                             2.4);
    genAbsJetRapidity_ZPt150_Zexc1jet = newTH1D("genAbsJetRapidity_ZPt150_Zexc1jet",
                                                "genAbsJetRapidity_ZPt150_Zexc1jet",
                                                "|y_{jet}|",
                                                12,
                                                0,
                                                2.4);
    SumZJetRapidity_ZPt150_Zexc1jet = newTH1D("SumZJetRapidity_ZPt150_Zexc1jet",
                                              "SumZJetRapidity_ZPt150_Zexc1jet",
                                              "y_{sum}",
                                              12,
                                              0,
                                              2.4);
    genSumZJetRapidity_ZPt150_Zexc1jet = newTH1D("genSumZJetRapidity_ZPt150_Zexc1jet",
                                                 "genSumZJetRapidity_ZPt150_Zexc1jet",
                                                 "y_{sum}",
                                                 12,
                                                 0,
                                                 2.4);
    DifZJetRapidity_ZPt150_Zexc1jet = newTH1D("DifZJetRapidity_ZPt150_Zexc1jet",
                                              "DifZJetRapidity_ZPt150_Zexc1jet",
                                              "y_{diff}",
                                              12,
                                              0,
                                              2.4);
    genDifZJetRapidity_ZPt150_Zexc1jet = newTH1D("genDifZJetRapidity_ZPt150_Zexc1jet",
                                                 "genDifZJetRapidity_ZPt150_Zexc1jet",
                                                 "y_{diff}",
                                                 12,
                                                 0,
                                                 2.4);

    AbsZRapidity_ZPt150_Zinc2jet = newTH1D(
        "AbsZRapidity_ZPt150_Zinc2jet", "AbsZRapidity_ZPt150_Zinc2jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_ZPt150_Zinc2jet = newTH1D("genAbsZRapidity_ZPt150_Zinc2jet",
                                              "genAbsZRapidity_ZPt150_Zinc2jet",
                                              "|y_{Z}|",
                                              12,
                                              0,
                                              2.4);
    AbsSecondJetRapidity_ZPt150_Zinc2jet = newTH1D("AbsSecondJetRapidity_ZPt150_Zinc2jet",
                                                   "AbsSecondJetRapidity_ZPt150_Zinc2jet",
                                                   "|y_{jet}|",
                                                   12,
                                                   0,
                                                   2.4);
    genAbsSecondJetRapidity_ZPt150_Zinc2jet = newTH1D("genAbsSecondJetRapidity_ZPt150_Zinc2jet",
                                                      "genAbsSecondJetRapidity_ZPt150_Zinc2jet",
                                                      "|y_{jet}|",
                                                      12,
                                                      0,
                                                      2.4);
    SumZSecondJetRapidity_ZPt150_Zinc2jet = newTH1D("SumZSecondJetRapidity_ZPt150_Zinc2jet",
                                                    "SumZSecondJetRapidity_ZPt150_Zinc2jet",
                                                    "y_{sum}",
                                                    12,
                                                    0,
                                                    2.4);
    genSumZSecondJetRapidity_ZPt150_Zinc2jet = newTH1D("genSumZSecondJetRapidity_ZPt150_Zinc2jet",
                                                       "genSumZSecondJetRapidity_ZPt150_Zinc2jet",
                                                       "y_{sum}",
                                                       12,
                                                       0,
                                                       2.4);
    DifZSecondJetRapidity_ZPt150_Zinc2jet = newTH1D("DifZSecondJetRapidity_ZPt150_Zinc2jet",
                                                    "DifZSecondJetRapidity_ZPt150_Zinc2jet",
                                                    "y_{diff}",
                                                    12,
                                                    0,
                                                    2.4);
    genDifZSecondJetRapidity_ZPt150_Zinc2jet = newTH1D("genDifZSecondJetRapidity_ZPt150_Zinc2jet",
                                                       "genDifZSecondJetRapidity_ZPt150_Zinc2jet",
                                                       "y_{diff}",
                                                       12,
                                                       0,
                                                       2.4);

    AbsZRapidity_ZPt150_Zexc2jet = newTH1D(
        "AbsZRapidity_ZPt150_Zexc2jet", "AbsZRapidity_ZPt150_Zexc2jet", "|y_{Z}|", 12, 0, 2.4);
    genAbsZRapidity_ZPt150_Zexc2jet = newTH1D("genAbsZRapidity_ZPt150_Zexc2jet",
                                              "genAbsZRapidity_ZPt150_Zexc2jet",
                                              "|y_{Z}|",
                                              12,
                                              0,
                                              2.4);
    AbsSecondJetRapidity_ZPt150_Zexc2jet = newTH1D("AbsSecondJetRapidity_ZPt150_Zexc2jet",
                                                   "AbsSecondJetRapidity_ZPt150_Zexc2jet",
                                                   "|y_{jet}|",
                                                   12,
                                                   0,
                                                   2.4);
    genAbsSecondJetRapidity_ZPt150_Zexc2jet = newTH1D("genAbsSecondJetRapidity_ZPt150_Zexc2jet",
                                                      "genAbsSecondJetRapidity_ZPt150_Zexc2jet",
                                                      "|y_{jet}|",
                                                      12,
                                                      0,
                                                      2.4);
    SumZSecondJetRapidity_ZPt150_Zexc2jet = newTH1D("SumZSecondJetRapidity_ZPt150_Zexc2jet",
                                                    "SumZSecondJetRapidity_ZPt150_Zexc2jet",
                                                    "y_{sum}",
                                                    12,
                                                    0,
                                                    2.4);
    genSumZSecondJetRapidity_ZPt150_Zexc2jet = newTH1D("genSumZSecondJetRapidity_ZPt150_Zexc2jet",
                                                       "genSumZSecondJetRapidity_ZPt150_Zexc2jet",
                                                       "y_{sum}",
                                                       12,
                                                       0,
                                                       2.4);
    DifZSecondJetRapidity_ZPt150_Zexc2jet = newTH1D("DifZSecondJetRapidity_ZPt150_Zexc2jet",
                                                    "DifZSecondJetRapidity_ZPt150_Zexc2jet",
                                                    "y_{diff}",
                                                    12,
                                                    0,
                                                    2.4);
    genDifZSecondJetRapidity_ZPt150_Zexc2jet = newTH1D("genDifZSecondJetRapidity_ZPt150_Zexc2jet",
                                                       "genDifZSecondJetRapidity_ZPt150_Zexc2jet",
                                                       "y_{diff}",
                                                       12,
                                                       0,
                                                       2.4);

    int nRapidity_ZPt300_Zinc1jet(7);
    double Rapidity_ZPt300_Zinc1jet[8] = {0, 0.4, 0.6, 0.9, 1.2, 1.5, 1.9, 2.4};

    AbsZRapidity_ZPt300_Zinc1jet = newTH1D("AbsZRapidity_ZPt300_Zinc1jet",
                                           "AbsZRapidity_ZPt300_Zinc1jet",
                                           "|y_{Z}|",
                                           nRapidity_ZPt300_Zinc1jet,
                                           Rapidity_ZPt300_Zinc1jet);
    genAbsZRapidity_ZPt300_Zinc1jet = newTH1D("genAbsZRapidity_ZPt300_Zinc1jet",
                                              "genAbsZRapidity_ZPt300_Zinc1jet",
                                              "|y_{Z}|",
                                              nRapidity_ZPt300_Zinc1jet,
                                              Rapidity_ZPt300_Zinc1jet);
    AbsFirstJetRapidity_ZPt300_Zinc1jet = newTH1D("AbsFirstJetRapidity_ZPt300_Zinc1jet",
                                                  "AbsFirstJetRapidity_ZPt300_Zinc1jet",
                                                  "|y_{jet1}|",
                                                  nRapidity_ZPt300_Zinc1jet,
                                                  Rapidity_ZPt300_Zinc1jet);
    genAbsFirstJetRapidity_ZPt300_Zinc1jet = newTH1D("genAbsFirstJetRapidity_ZPt300_Zinc1jet",
                                                     "genAbsFirstJetRapidity_ZPt300_Zinc1jet",
                                                     "|y_{jet1}|",
                                                     nRapidity_ZPt300_Zinc1jet,
                                                     Rapidity_ZPt300_Zinc1jet);
    SumZFirstJetRapidity_ZPt300_Zinc1jet = newTH1D("SumZFirstJetRapidity_ZPt300_Zinc1jet",
                                                   "SumZFirstJetRapidity_ZPt300_Zinc1jet",
                                                   "y_{sum(Z,jet1)}",
                                                   nRapidity_ZPt300_Zinc1jet,
                                                   Rapidity_ZPt300_Zinc1jet);
    genSumZFirstJetRapidity_ZPt300_Zinc1jet = newTH1D("genSumZFirstJetRapidity_ZPt300_Zinc1jet",
                                                      "genSumZFirstJetRapidity_ZPt300_Zinc1jet",
                                                      "y_{sum(Z,jet1)}",
                                                      nRapidity_ZPt300_Zinc1jet,
                                                      Rapidity_ZPt300_Zinc1jet);
    DifZFirstJetRapidity_ZPt300_Zinc1jet = newTH1D("DifZFirstJetRapidity_ZPt300_Zinc1jet",
                                                   "DifZFirstJetRapidity_ZPt300_Zinc1jet",
                                                   "y_{diff(Z,jet1)}",
                                                   nRapidity_ZPt300_Zinc1jet,
                                                   Rapidity_ZPt300_Zinc1jet);
    genDifZFirstJetRapidity_ZPt300_Zinc1jet = newTH1D("genDifZFirstJetRapidity_ZPt300_Zinc1jet",
                                                      "genDifZFirstJetRapidity_ZPt300_Zinc1jet",
                                                      "y_{diff(Z,jet1)}",
                                                      nRapidity_ZPt300_Zinc1jet,
                                                      Rapidity_ZPt300_Zinc1jet);

    /////Azimuthal cross check/////////////////

    int nDPhiZFirstJet_ZPt150(13);
    double DPhiZFirstJet_ZPt150[14] = {0,
                                       0.16 * PI,
                                       0.28 * PI,
                                       0.40 * PI,
                                       0.52 * PI,
                                       0.64 * PI,
                                       0.72 * PI,
                                       0.76 * PI,
                                       0.8 * PI,
                                       0.84 * PI,
                                       0.88 * PI,
                                       0.92 * PI,
                                       0.96 * PI,
                                       1.0 * PI};

    int nDPhiZSecondJet_ZPt150(15);
    double DPhiZSecondJet_ZPt150[16] = {0,
                                        0.12 * PI,
                                        0.2 * PI,
                                        0.28 * PI,
                                        0.36 * PI,
                                        0.44 * PI,
                                        0.52 * PI,
                                        0.6 * PI,
                                        0.68 * PI,
                                        0.76 * PI,
                                        0.8 * PI,
                                        0.84 * PI,
                                        0.88 * PI,
                                        0.92 * PI,
                                        0.96 * PI,
                                        1.0 * PI};

    int nDPhiZThirdJet_ZPt150(23);
    double DPhiZThirdJet_ZPt150[24] = {
        0,         0.12 * PI, 0.16 * PI, 0.2 * PI,  0.24 * PI, 0.28 * PI, 0.32 * PI, 0.36 * PI,
        0.4 * PI,  0.44 * PI, 0.48 * PI, 0.52 * PI, 0.56 * PI, 0.6 * PI,  0.64 * PI, 0.68 * PI,
        0.72 * PI, 0.76 * PI, 0.8 * PI,  0.84 * PI, 0.88 * PI, 0.92 * PI, 0.96 * PI, 1.0 * PI};

    DPhiZFirstJet_ZPt150_Zinc1jet = newTH1D("DPhiZFirstJet_ZPt150_Zinc1jet",
                                            "DPhiZFirstJet_ZPt150_Zinc1jet",
                                            "#Delta#phi_{Z,Jet1}",
                                            nDPhiZFirstJet_ZPt150,
                                            DPhiZFirstJet_ZPt150);
    genDPhiZFirstJet_ZPt150_Zinc1jet = newTH1D("genDPhiZFirstJet_ZPt150_Zinc1jet",
                                               "genDPhiZFirstJet_ZPt150_Zinc1jet",
                                               "#Delta#phi_{Z,Jet1}",
                                               nDPhiZFirstJet_ZPt150,
                                               DPhiZFirstJet_ZPt150);
    DPhiZFirstJet_ZPt150_Zinc2jet = newTH1D("DPhiZFirstJet_ZPt150_Zinc2jet",
                                            "DPhiZFirstJet_ZPt150_Zinc2jet",
                                            "#Delta#phi_{Z,Jet1}",
                                            nDPhiZFirstJet_ZPt150,
                                            DPhiZFirstJet_ZPt150);
    genDPhiZFirstJet_ZPt150_Zinc2jet = newTH1D("genDPhiZFirstJet_ZPt150_Zinc2jet",
                                               "genDPhiZFirstJet_ZPt150_Zinc2jet",
                                               "#Delta#phi_{Z,Jet1}",
                                               nDPhiZFirstJet_ZPt150,
                                               DPhiZFirstJet_ZPt150);
    DPhiZFirstJet_ZPt150_Zinc3jet = newTH1D("DPhiZFirstJet_ZPt150_Zinc3jet",
                                            "DPhiZFirstJet_ZPt150_Zinc3jet",
                                            "#Delta#phi_{Z,Jet1}",
                                            nDPhiZFirstJet_ZPt150,
                                            DPhiZFirstJet_ZPt150);
    genDPhiZFirstJet_ZPt150_Zinc3jet = newTH1D("genDPhiZFirstJet_ZPt150_Zinc3jet",
                                               "genDPhiZFirstJet_ZPt150_Zinc3jet",
                                               "#Delta#phi_{Z,Jet1}",
                                               nDPhiZFirstJet_ZPt150,
                                               DPhiZFirstJet_ZPt150);

    DPhiZSecondJet_ZPt150_Zinc3jet = newTH1D("DPhiZSecondJet_ZPt150_Zinc3jet",
                                             "DPhiZSecondJet_ZPt150_Zinc3jet",
                                             "#Delta#phi_{Z,Jet2}",
                                             nDPhiZSecondJet_ZPt150,
                                             DPhiZSecondJet_ZPt150);
    genDPhiZSecondJet_ZPt150_Zinc3jet = newTH1D("genDPhiZSecondJet_ZPt150_Zinc3jet",
                                                "genDPhiZSecondJet_ZPt150_Zinc3jet",
                                                "#Delta#phi_{Z,Jet2}",
                                                nDPhiZSecondJet_ZPt150,
                                                DPhiZSecondJet_ZPt150);
    DPhiZThirdJet_ZPt150_Zinc3jet = newTH1D("DPhiZThirdJet_ZPt150_Zinc3jet",
                                            "DPhiZThirdJet_ZPt150_Zinc3jet",
                                            "#Delta#phi_{Z,Jet3}",
                                            nDPhiZThirdJet_ZPt150,
                                            DPhiZThirdJet_ZPt150);
    genDPhiZThirdJet_ZPt150_Zinc3jet = newTH1D("genDPhiZThirdJet_ZPt150_Zinc3jet",
                                               "genDPhiZThirdJet_ZPt150_Zinc3jet",
                                               "#Delta#phi_{Z,Jet3}",
                                               nDPhiZThirdJet_ZPt150,
                                               DPhiZThirdJet_ZPt150);
    DPhiFirstSecondJet_ZPt150_Zinc3jet = newTH1D("DPhiFirstSecondJet_ZPt150_Zinc3jet",
                                                 "DPhiFirstSecondJet_ZPt150_Zinc3jet",
                                                 "#Delta#phi_{Jet1,Jet2}",
                                                 25,
                                                 0,
                                                 3.14);
    genDPhiFirstSecondJet_ZPt150_Zinc3jet = newTH1D("genDPhiFirstSecondJet_ZPt150_Zinc3jet",
                                                    "genDPhiFirstSecondJet_ZPt150_Zinc3jet",
                                                    "#Delta#phi_{Jet1,Jet2}",
                                                    25,
                                                    0,
                                                    3.14);
    DPhiFirstThirdJet_ZPt150_Zinc3jet = newTH1D("DPhiFirstThirdJet_ZPt150_Zinc3jet",
                                                "DPhiFirstThirdJet_ZPt150_Zinc3jet",
                                                "#Delta#phi_{Jet1,Jet3}",
                                                25,
                                                0,
                                                3.14);
    genDPhiFirstThirdJet_ZPt150_Zinc3jet = newTH1D("genDPhiFirstThirdJet_ZPt150_Zinc3jet",
                                                   "genDPhiFirstThirdJet_ZPt150_Zinc3jet",
                                                   "#Delta#phi_{Jet1,Jet3}",
                                                   25,
                                                   0,
                                                   3.14);
    DPhiSecondThirdJet_ZPt150_Zinc3jet = newTH1D("DPhiSecondThirdJet_ZPt150_Zinc3jet",
                                                 "DPhiSecondThirdJet_ZPt150_Zinc3jet",
                                                 "#Delta#phi_{Jet2,Jet3}",
                                                 25,
                                                 0,
                                                 3.14);
    genDPhiSecondThirdJet_ZPt150_Zinc3jet = newTH1D("genDPhiSecondThirdJet_ZPt150_Zinc3jet",
                                                    "genDPhiSecondThirdJet_ZPt150_Zinc3jet",
                                                    "#Delta#phi_{Jet2,Jet3}",
                                                    25,
                                                    0,
                                                    3.14);

    int nDPhiZFirstJet_ZPt300(8);
    double DPhiZFirstJet_ZPt300[9] = {0.52 * PI,
                                      0.68 * PI,
                                      0.76 * PI,
                                      0.8 * PI,
                                      0.84 * PI,
                                      0.88 * PI,
                                      0.92 * PI,
                                      0.96 * PI,
                                      1.0 * PI};

    int nDPhiZSecondJet_ZPt300(12);
    double DPhiZSecondJet_ZPt300[13] = {0,
                                        0.154 * PI,
                                        0.231 * PI,
                                        0.308 * PI,
                                        0.384 * PI,
                                        0.461 * PI,
                                        0.538 * PI,
                                        0.615 * PI,
                                        0.692 * PI,
                                        0.769 * PI,
                                        0.846 * PI,
                                        0.923 * PI,
                                        1.0 * PI};

    int nDPhiZThirdJet_ZPt300(8);
    double DPhiZThirdJet_ZPt300[9] = {0,
                                      0.231 * PI,
                                      0.384 * PI,
                                      0.538 * PI,
                                      0.692 * PI,
                                      0.769 * PI,
                                      0.846 * PI,
                                      0.923 * PI,
                                      1.0 * PI};

    int nDPhiJets_ZPt300(9);
    double DPhiJets_ZPt300[10] = {0,
                                  0.08 * PI,
                                  0.16 * PI,
                                  0.24 * PI,
                                  0.32 * PI,
                                  0.40 * PI,
                                  0.52 * PI,
                                  0.64 * PI,
                                  0.76 * PI,
                                  1.0 * PI};

    DPhiZFirstJet_ZPt300_Zinc1jet = newTH1D("DPhiZFirstJet_ZPt300_Zinc1jet",
                                            "DPhiZFirstJet_ZPt300_Zinc1jet",
                                            "#Delta#phi_{Z,Jet1}",
                                            nDPhiZFirstJet_ZPt300,
                                            DPhiZFirstJet_ZPt300);
    genDPhiZFirstJet_ZPt300_Zinc1jet = newTH1D("genDPhiZFirstJet_ZPt300_Zinc1jet",
                                               "genDPhiZFirstJet_ZPt300_Zinc1jet",
                                               "#Delta#phi_{Z,Jet1}",
                                               nDPhiZFirstJet_ZPt300,
                                               DPhiZFirstJet_ZPt300);
    DPhiZFirstJet_ZPt300_Zinc2jet = newTH1D("DPhiZFirstJet_ZPt300_Zinc2jet",
                                            "DPhiZFirstJet_ZPt300_Zinc2jet",
                                            "#Delta#phi_{Z,Jet1}",
                                            nDPhiZFirstJet_ZPt300,
                                            DPhiZFirstJet_ZPt300);
    genDPhiZFirstJet_ZPt300_Zinc2jet = newTH1D("genDPhiZFirstJet_ZPt300_Zinc2jet",
                                               "genDPhiZFirstJet_ZPt300_Zinc2jet",
                                               "#Delta#phi_{Z,Jet1}",
                                               nDPhiZFirstJet_ZPt300,
                                               DPhiZFirstJet_ZPt300);
    DPhiZFirstJet_ZPt300_Zinc3jet = newTH1D("DPhiZFirstJet_ZPt300_Zinc3jet",
                                            "DPhiZFirstJet_ZPt300_Zinc3jet",
                                            "#Delta#phi_{Z,Jet1}",
                                            nDPhiZFirstJet_ZPt300,
                                            DPhiZFirstJet_ZPt300);
    genDPhiZFirstJet_ZPt300_Zinc3jet = newTH1D("genDPhiZFirstJet_ZPt300_Zinc3jet",
                                               "genDPhiZFirstJet_ZPt300_Zinc3jet",
                                               "#Delta#phi_{Z,Jet1}",
                                               nDPhiZFirstJet_ZPt300,
                                               DPhiZFirstJet_ZPt300);

    DPhiZSecondJet_ZPt300_Zinc3jet = newTH1D("DPhiZSecondJet_ZPt300_Zinc3jet",
                                             "DPhiZSecondJet_ZPt300_Zinc3jet",
                                             "#Delta#phi_{Z,Jet2}",
                                             nDPhiZSecondJet_ZPt300,
                                             DPhiZSecondJet_ZPt300);
    genDPhiZSecondJet_ZPt300_Zinc3jet = newTH1D("genDPhiZSecondJet_ZPt300_Zinc3jet",
                                                "genDPhiZSecondJet_ZPt300_Zinc3jet",
                                                "#Delta#phi_{Z,Jet2}",
                                                nDPhiZSecondJet_ZPt300,
                                                DPhiZSecondJet_ZPt300);
    DPhiZThirdJet_ZPt300_Zinc3jet = newTH1D("DPhiZThirdJet_ZPt300_Zinc3jet",
                                            "DPhiZThirdJet_ZPt300_Zinc3jet",
                                            "#Delta#phi_{Z,Jet3}",
                                            nDPhiZThirdJet_ZPt300,
                                            DPhiZThirdJet_ZPt300);
    genDPhiZThirdJet_ZPt300_Zinc3jet = newTH1D("genDPhiZThirdJet_ZPt300_Zinc3jet",
                                               "genDPhiZThirdJet_ZPt300_Zinc3jet",
                                               "#Delta#phi_{Z,Jet3}",
                                               nDPhiZThirdJet_ZPt300,
                                               DPhiZThirdJet_ZPt300);
    DPhiFirstSecondJet_ZPt300_Zinc3jet = newTH1D("DPhiFirstSecondJet_ZPt300_Zinc3jet",
                                                 "DPhiFirstSecondJet_ZPt300_Zinc3jet",
                                                 "#Delta#phi_{Jet1,Jet2}",
                                                 nDPhiJets_ZPt300,
                                                 DPhiJets_ZPt300);
    genDPhiFirstSecondJet_ZPt300_Zinc3jet = newTH1D("genDPhiFirstSecondJet_ZPt300_Zinc3jet",
                                                    "genDPhiFirstSecondJet_ZPt300_Zinc3jet",
                                                    "#Delta#phi_{Jet1,Jet2}",
                                                    nDPhiJets_ZPt300,
                                                    DPhiJets_ZPt300);
    DPhiFirstThirdJet_ZPt300_Zinc3jet = newTH1D("DPhiFirstThirdJet_ZPt300_Zinc3jet",
                                                "DPhiFirstThirdJet_ZPt300_Zinc3jet",
                                                "#Delta#phi_{Jet1,Jet3}",
                                                nDPhiJets_ZPt300,
                                                DPhiJets_ZPt300);
    genDPhiFirstThirdJet_ZPt300_Zinc3jet = newTH1D("genDPhiFirstThirdJet_ZPt300_Zinc3jet",
                                                   "genDPhiFirstThirdJet_ZPt300_Zinc3jet",
                                                   "#Delta#phi_{Jet1,Jet3}",
                                                   nDPhiJets_ZPt300,
                                                   DPhiJets_ZPt300);
    DPhiSecondThirdJet_ZPt300_Zinc3jet = newTH1D("DPhiSecondThirdJet_ZPt300_Zinc3jet",
                                                 "DPhiSecondThirdJet_ZPt300_Zinc3jet",
                                                 "#Delta#phi_{Jet2,Jet3}",
                                                 nDPhiJets_ZPt300,
                                                 DPhiJets_ZPt300);
    genDPhiSecondThirdJet_ZPt300_Zinc3jet = newTH1D("genDPhiSecondThirdJet_ZPt300_Zinc3jet",
                                                    "genDPhiSecondThirdJet_ZPt300_Zinc3jet",
                                                    "#Delta#phi_{Jet2,Jet3}",
                                                    nDPhiJets_ZPt300,
                                                    DPhiJets_ZPt300);

    int nDPhiZFirstJet_ZPt150_HT300(11);
    double DPhiZFirstJet_ZPt150_HT300[12] = {0,
                                             0.48 * PI,
                                             0.64 * PI,
                                             0.68 * PI,
                                             0.72 * PI,
                                             0.76 * PI,
                                             0.8 * PI,
                                             0.84 * PI,
                                             0.88 * PI,
                                             0.92 * PI,
                                             0.96 * PI,
                                             1.0 * PI};

    int nDPhiZSecondJet_ZPt150_HT300(12);
    double DPhiZSecondJet_ZPt150_HT300[13] = {0,
                                              0.12 * PI,
                                              0.2 * PI,
                                              0.28 * PI,
                                              0.36 * PI,
                                              0.44 * PI,
                                              0.52 * PI,
                                              0.6 * PI,
                                              0.68 * PI,
                                              0.76 * PI,
                                              0.84 * PI,
                                              0.92 * PI,
                                              1.0 * PI};

    int nDPhiZThirdJet_ZPt150_HT300(12);
    double DPhiZThirdJet_ZPt150_HT300[13] = {0,
                                             0.12 * PI,
                                             0.2 * PI,
                                             0.28 * PI,
                                             0.36 * PI,
                                             0.44 * PI,
                                             0.52 * PI,
                                             0.6 * PI,
                                             0.68 * PI,
                                             0.76 * PI,
                                             0.84 * PI,
                                             0.92 * PI,
                                             1.0 * PI};

    DPhiZFirstJet_ZPt150_HT300_Zinc3jet = newTH1D("DPhiZFirstJet_ZPt150_HT300_Zinc3jet",
                                                  "DPhiZFirstJet_ZPt150_HT300_Zinc3jet",
                                                  "#Delta#phi_{Z,Jet1}",
                                                  nDPhiZFirstJet_ZPt150_HT300,
                                                  DPhiZFirstJet_ZPt150_HT300);
    genDPhiZFirstJet_ZPt150_HT300_Zinc3jet = newTH1D("genDPhiZFirstJet_ZPt150_HT300_Zinc3jet",
                                                     "genDPhiZFirstJet_ZPt150_HT300_Zinc3jet",
                                                     "#Delta#phi_{Z,Jet1}",
                                                     nDPhiZFirstJet_ZPt150_HT300,
                                                     DPhiZFirstJet_ZPt150_HT300);
    DPhiZSecondJet_ZPt150_HT300_Zinc3jet = newTH1D("DPhiZSecondJet_ZPt150_HT300_Zinc3jet",
                                                   "DPhiZSecondJet_ZPt150_HT300_Zinc3jet",
                                                   "#Delta#phi_{Z,Jet2}",
                                                   nDPhiZSecondJet_ZPt150_HT300,
                                                   DPhiZSecondJet_ZPt150_HT300);
    genDPhiZSecondJet_ZPt150_HT300_Zinc3jet = newTH1D("genDPhiZSecondJet_ZPt150_HT300_Zinc3jet",
                                                      "genDPhiZSecondJet_ZPt150_HT300_Zinc3jet",
                                                      "#Delta#phi_{Z,Jet2}",
                                                      nDPhiZSecondJet_ZPt150_HT300,
                                                      DPhiZSecondJet_ZPt150_HT300);
    DPhiZThirdJet_ZPt150_HT300_Zinc3jet = newTH1D("DPhiZThirdJet_ZPt150_HT300_Zinc3jet",
                                                  "DPhiZThirdJet_ZPt150_HT300_Zinc3jet",
                                                  "#Delta#phi_{Z,Jet3}",
                                                  nDPhiZThirdJet_ZPt150_HT300,
                                                  DPhiZThirdJet_ZPt150_HT300);
    genDPhiZThirdJet_ZPt150_HT300_Zinc3jet = newTH1D("genDPhiZThirdJet_ZPt150_HT300_Zinc3jet",
                                                     "genDPhiZThirdJet_ZPt150_HT300_Zinc3jet",
                                                     "#Delta#phi_{Z,Jet3}",
                                                     nDPhiZThirdJet_ZPt150_HT300,
                                                     DPhiZThirdJet_ZPt150_HT300);

    // Branches with different JetPt cuts///////

    AbsZRapidity_FirstJetPt50_Zinc1jet = newTH1D("AbsZRapidity_FirstJetPt50_Zinc1jet",
                                                 "AbsZRapidity_FirstJetPt50_Zinc1jet",
                                                 "|y_{Z}|",
                                                 12,
                                                 0,
                                                 2.4);
    genAbsZRapidity_FirstJetPt50_Zinc1jet = newTH1D("genAbsZRapidity_FirstJetPt50_Zinc1jet",
                                                    "genAbsZRapidity_FirstJetPt50_Zinc1jet",
                                                    "|y_{Z}|",
                                                    12,
                                                    0,
                                                    2.4);
    AbsFirstJetRapidity_FirstJetPt50_Zinc1jet = newTH1D("AbsFirstJetRapidity_FirstJetPt50_Zinc1jet",
                                                        "AbsFirstJetRapidity_FirstJetPt50_Zinc1jet",
                                                        "|y_{jet1}|",
                                                        12,
                                                        0,
                                                        2.4);
    genAbsFirstJetRapidity_FirstJetPt50_Zinc1jet =
        newTH1D("genAbsFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "genAbsFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "|y_{jet1}|",
                12,
                0,
                2.4);
    SumZFirstJetRapidity_FirstJetPt50_Zinc1jet =
        newTH1D("SumZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "SumZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "y_{sum(Z,jet1)}",
                12,
                0,
                2.4);
    genSumZFirstJetRapidity_FirstJetPt50_Zinc1jet =
        newTH1D("genSumZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "genSumZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "y_{sum(Z,jet1)}",
                12,
                0,
                2.4);
    DifZFirstJetRapidity_FirstJetPt50_Zinc1jet =
        newTH1D("DifZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "DifZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "y_{diff(Z,jet1)}",
                12,
                0,
                2.4);
    genDifZFirstJetRapidity_FirstJetPt50_Zinc1jet =
        newTH1D("genDifZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "genDifZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "y_{diff(Z,jet1)}",
                12,
                0,
                2.4);

    AbsZRapidity_FirstJetPt80_Zinc1jet = newTH1D("AbsZRapidity_FirstJetPt80_Zinc1jet",
                                                 "AbsZRapidity_FirstJetPt80_Zinc1jet",
                                                 "|y_{Z}|",
                                                 12,
                                                 0,
                                                 2.4);
    genAbsZRapidity_FirstJetPt80_Zinc1jet = newTH1D("genAbsZRapidity_FirstJetPt80_Zinc1jet",
                                                    "genAbsZRapidity_FirstJetPt80_Zinc1jet",
                                                    "|y_{Z}|",
                                                    12,
                                                    0,
                                                    2.4);
    AbsFirstJetRapidity_FirstJetPt80_Zinc1jet = newTH1D("AbsFirstJetRapidity_FirstJetPt80_Zinc1jet",
                                                        "AbsFirstJetRapidity_FirstJetPt80_Zinc1jet",
                                                        "|y_{jet1}|",
                                                        12,
                                                        0,
                                                        2.4);
    genAbsFirstJetRapidity_FirstJetPt80_Zinc1jet =
        newTH1D("genAbsFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "genAbsFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "|y_{jet1}|",
                12,
                0,
                2.4);
    SumZFirstJetRapidity_FirstJetPt80_Zinc1jet =
        newTH1D("SumZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "SumZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "y_{sum(Z,jet1)}",
                12,
                0,
                2.4);
    genSumZFirstJetRapidity_FirstJetPt80_Zinc1jet =
        newTH1D("genSumZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "genSumZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "y_{sum(Z,jet1)}",
                12,
                0,
                2.4);
    DifZFirstJetRapidity_FirstJetPt80_Zinc1jet =
        newTH1D("DifZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "DifZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "y_{diff(Z,jet1)}",
                12,
                0,
                2.4);
    genDifZFirstJetRapidity_FirstJetPt80_Zinc1jet =
        newTH1D("genDifZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "genDifZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "y_{diff(Z,jet1)}",
                12,
                0,
                2.4);

    // Set jet rapidity discriminator/////
    AbsZRapidity_DifJetRapidityl2_Zinc2jet = newTH1D("AbsZRapidity_DifJetRapidityl2_Zinc2jet",
                                                     "AbsZRapidity_DifJetRapidityl2_Zinc2jet",
                                                     "|y_{Z}|",
                                                     12,
                                                     0,
                                                     2.4);
    genAbsZRapidity_DifJetRapidityl2_Zinc2jet = newTH1D("genAbsZRapidity_DifJetRapidityl2_Zinc2jet",
                                                        "genAbsZRapidity_DifJetRapidityl2_Zinc2jet",
                                                        "|y_{Z}|",
                                                        12,
                                                        0,
                                                        2.4);
    AbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet =
        newTH1D("AbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "AbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "|y_{jet}|",
                12,
                0,
                2.4);
    genAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet =
        newTH1D("genAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "genAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "|y_{jet}|",
                12,
                0,
                2.4);
    SumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet =
        newTH1D("SumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "SumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "y_{sum}",
                12,
                0,
                2.4);
    genSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet =
        newTH1D("genSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "genSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "y_{sum}",
                12,
                0,
                2.4);
    DifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet =
        newTH1D("DifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "DifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "y_{diff}",
                12,
                0,
                2.4);
    genDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet =
        newTH1D("genDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "genDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "y_{diff}",
                12,
                0,
                2.4);

    AbsZRapidity_DifJetRapiditys2_Zinc2jet = newTH1D("AbsZRapidity_DifJetRapiditys2_Zinc2jet",
                                                     "AbsZRapidity_DifJetRapiditys2_Zinc2jet",
                                                     "|y_{Z}|",
                                                     12,
                                                     0,
                                                     2.4);
    genAbsZRapidity_DifJetRapiditys2_Zinc2jet = newTH1D("genAbsZRapidity_DifJetRapiditys2_Zinc2jet",
                                                        "genAbsZRapidity_DifJetRapiditys2_Zinc2jet",
                                                        "|y_{Z}|",
                                                        12,
                                                        0,
                                                        2.4);
    AbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet =
        newTH1D("AbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "AbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "|y_{jet}|",
                12,
                0,
                2.4);
    genAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet =
        newTH1D("genAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "genAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "|y_{jet}|",
                12,
                0,
                2.4);
    SumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet =
        newTH1D("SumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "SumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "y_{sum}",
                12,
                0,
                2.4);
    genSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet =
        newTH1D("genSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "genSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "y_{sum}",
                12,
                0,
                2.4);
    DifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet =
        newTH1D("DifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "DifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "y_{diff}",
                12,
                0,
                2.4);
    genDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet =
        newTH1D("genDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "genDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "y_{diff}",
                12,
                0,
                2.4);

    //--Additional response TH2D
    hresponseAbsZRapidity_Zinc1jet = newTH2D(
        "hresponseAbsZRapidity_Zinc1jet", "hresponseAbsZRapidity_Zinc1jet", 12, 0, 2.4, 12, 0, 2.4);
    hresponseAbsFirstJetRapidity_Zinc1jet = newTH2D("hresponseAbsFirstJetRapidity_Zinc1jet",
                                                    "hresponseAbsFirstJetRapidity_Zinc1jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseSumZFirstJetRapidity_Zinc1jet = newTH2D("hresponseSumZFirstJetRapidity_Zinc1jet",
                                                     "hresponseSumZFirstJetRapidity_Zinc1jet",
                                                     12,
                                                     0,
                                                     2.4,
                                                     12,
                                                     0,
                                                     2.4);
    hresponseDifZFirstJetRapidity_Zinc1jet = newTH2D("hresponseDifZFirstJetRapidity_Zinc1jet",
                                                     "hresponseDifZFirstJetRapidity_Zinc1jet",
                                                     12,
                                                     0,
                                                     2.4,
                                                     12,
                                                     0,
                                                     2.4);

    /// cross check//////
    hresponseSumZFirstJetEta_Zinc1jet = newTH2D("hresponseSumZFirstJetEta_Zinc1jet",
                                                "hresponseSumZFirstJetEta_Zinc1jet",
                                                12,
                                                0,
                                                2.4,
                                                12,
                                                0,
                                                2.4);
    hresponseDifZFirstJetEta_Zinc1jet = newTH2D("hresponseDifZFirstJetEta_Zinc1jet",
                                                "hresponseDifZFirstJetEta_Zinc1jet",
                                                12,
                                                0,
                                                2.4,
                                                12,
                                                0,
                                                2.4);

    ////Azimuth cross check////////////////
    hresponseDPhiZFirstJet_Zinc1jet = newTH2D("hresponseDPhiZFirstJet_Zinc1jet",
                                              "hresponseDPhiZFirstJet_Zinc1jet",
                                              25,
                                              0,
                                              3.14,
                                              25,
                                              0,
                                              3.14);

    hresponseAbsZRapidity_Zexc1jet = newTH2D(
        "hresponseAbsZRapidity_Zexc1jet", "hresponseAbsZRapidity_Zexc1jet", 12, 0, 2.4, 12, 0, 2.4);
    hresponseAbsJetRapidity_Zexc1jet = newTH2D("hresponseAbsJetRapidity_Zexc1jet",
                                               "hresponseAbsJetRapidity_Zexc1jet",
                                               12,
                                               0,
                                               2.4,
                                               12,
                                               0,
                                               2.4);
    hresponseSumZJetRapidity_Zexc1jet = newTH2D("hresponseSumZJetRapidity_Zexc1jet",
                                                "hresponseSumZJetRapidity_Zexc1jet",
                                                12,
                                                0,
                                                2.4,
                                                12,
                                                0,
                                                2.4);
    hresponseDifZJetRapidity_Zexc1jet = newTH2D("hresponseDifZJetRapidity_Zexc1jet",
                                                "hresponseDifZJetRapidity_Zexc1jet",
                                                12,
                                                0,
                                                2.4,
                                                12,
                                                0,
                                                2.4);

    hresponseAbsFirstJetRapidity_Zinc2jet = newTH2D("hresponseAbsFirstJetRapidity_Zinc2jet",
                                                    "hresponseAbsFirstJetRapidity_Zinc2jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseSumZFirstJetRapidity_Zinc2jet = newTH2D("hresponseSumZFirstJetRapidity_Zinc2jet",
                                                     "hresponseSumZFirstJetRapidity_Zinc2jet",
                                                     12,
                                                     0,
                                                     2.4,
                                                     12,
                                                     0,
                                                     2.4);
    hresponseDifZFirstJetRapidity_Zinc2jet = newTH2D("hresponseDifZFirstJetRapidity_Zinc2jet",
                                                     "hresponseDifZFirstJetRapidity_Zinc2jet",
                                                     12,
                                                     0,
                                                     2.4,
                                                     12,
                                                     0,
                                                     2.4);

    hresponseAbsZRapidity_Zinc2jet = newTH2D(
        "hresponseAbsZRapidity_Zinc2jet", "hresponseAbsZRapidity_Zinc2jet", 12, 0, 2.4, 12, 0, 2.4);
    hresponseAbsSecondJetRapidity_Zinc2jet = newTH2D("hresponseAbsSecondJetRapidity_Zinc2jet",
                                                     "hresponseAbsSecondJetRapidity_Zinc2jet",
                                                     12,
                                                     0,
                                                     2.4,
                                                     12,
                                                     0,
                                                     2.4);
    hresponseSumZSecondJetRapidity_Zinc2jet = newTH2D("hresponseSumZSecondJetRapidity_Zinc2jet",
                                                      "hresponseSumZSecondJetRapidity_Zinc2jet",
                                                      12,
                                                      0,
                                                      2.4,
                                                      12,
                                                      0,
                                                      2.4);
    hresponseDifZSecondJetRapidity_Zinc2jet = newTH2D("hresponseDifZSecondJetRapidity_Zinc2jet",
                                                      "hresponseDifZSecondJetRapidity_Zinc2jet",
                                                      12,
                                                      0,
                                                      2.4,
                                                      12,
                                                      0,
                                                      2.4);

    hresponseSumFirstSecondJetRapidity_Zinc2jet =
        newTH2D("hresponseSumFirstSecondJetRapidity_Zinc2jet",
                "hresponseSumFirstSecondJetRapidity_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifFirstSecondJetRapidity_Zinc2jet =
        newTH2D("hresponseDifFirstSecondJetRapidity_Zinc2jet",
                "hresponseDifFirstSecondJetRapidity_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    hresponseSumZTwoJetsRapidity_Zinc2jet = newTH2D("hresponseSumZTwoJetsRapidity_Zinc2jet",
                                                    "hresponseSumZTwoJetsRapidity_Zinc2jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseDifZTwoJetsRapidity_Zinc2jet = newTH2D("hresponseDifZTwoJetsRapidity_Zinc2jet",
                                                    "hresponseDifZTwoJetsRapidity_Zinc2jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);

    hresponseDPhiZFirstJet_Zinc2jet = newTH2D("hresponseDPhiZFirstJet_Zinc2jet",
                                              "hresponseDPhiZFirstJet_Zinc2jet",
                                              25,
                                              0,
                                              3.14,
                                              25,
                                              0,
                                              3.14);
    hresponseDPhiZSecondJet_Zinc2jet = newTH2D("hresponseDPhiZSecondJet_Zinc2jet",
                                               "hresponseDPhiZSecondJet_Zinc2jet",
                                               25,
                                               0,
                                               3.14,
                                               25,
                                               0,
                                               3.14);
    hresponseDPhiFirstSecondJet_Zinc2jet = newTH2D("hresponseDPhiFirstSecondJet_Zinc2jet",
                                                   "hresponseDPhiFirstSecondJet_Zinc2jet",
                                                   25,
                                                   0,
                                                   3.14,
                                                   25,
                                                   0,
                                                   3.14);

    hresponseAbsZRapidity_Zexc2jet = newTH2D(
        "hresponseAbsZRapidity_Zexc2jet", "hresponseAbsZRapidity_Zexc2jet", 12, 0, 2.4, 12, 0, 2.4);
    hresponseAbsSecondJetRapidity_Zexc2jet = newTH2D("hresponseAbsSecondJetRapidity_Zexc2jet",
                                                     "hresponseAbsSecondJetRapidity_Zexc2jet",
                                                     12,
                                                     0,
                                                     2.4,
                                                     12,
                                                     0,
                                                     2.4);
    hresponseSumZSecondJetRapidity_Zexc2jet = newTH2D("hresponseSumZSecondJetRapidity_Zexc2jet",
                                                      "hresponseSumZSecondJetRapidity_Zexc2jet",
                                                      12,
                                                      0,
                                                      2.4,
                                                      12,
                                                      0,
                                                      2.4);
    hresponseDifZSecondJetRapidity_Zexc2jet = newTH2D("hresponseDifZSecondJetRapidity_Zexc2jet",
                                                      "hresponseDifZSecondJetRapidity_Zexc2jet",
                                                      12,
                                                      0,
                                                      2.4,
                                                      12,
                                                      0,
                                                      2.4);

    hresponseDPhiZFirstJet_Zinc3jet = newTH2D("hresponseDPhiZFirstJet_Zinc3jet",
                                              "hresponseDPhiZFirstJet_Zinc3jet",
                                              25,
                                              0,
                                              3.14,
                                              25,
                                              0,
                                              3.14);
    hresponseDPhiZSecondJet_Zinc3jet = newTH2D("hresponseDPhiZSecondJet_Zinc3jet",
                                               "hresponseDPhiZSecondJet_Zinc3jet",
                                               25,
                                               0,
                                               3.14,
                                               25,
                                               0,
                                               3.14);
    hresponseDPhiZThirdJet_Zinc3jet = newTH2D("hresponseDPhiZThirdJet_Zinc3jet",
                                              "hresponseDPhiZThirdJet_Zinc3jet",
                                              25,
                                              0,
                                              3.14,
                                              25,
                                              0,
                                              3.14);
    hresponseDPhiFirstSecondJet_Zinc3jet = newTH2D("hresponseDPhiFirstSecondJet_Zinc3jet",
                                                   "hresponseDPhiFirstSecondJet_Zinc3jet",
                                                   25,
                                                   0,
                                                   3.14,
                                                   25,
                                                   0,
                                                   3.14);
    hresponseDPhiFirstThirdJet_Zinc3jet = newTH2D("hresponseDPhiFirstThirdJet_Zinc3jet",
                                                  "hresponseDPhiFirstThirdJet_Zinc3jet",
                                                  25,
                                                  0,
                                                  3.14,
                                                  25,
                                                  0,
                                                  3.14);
    hresponseDPhiSecondThirdJet_Zinc3jet = newTH2D("hresponseDPhiSecondThirdJet_Zinc3jet",
                                                   "hresponseDPhiSecondThirdJet_Zinc3jet",
                                                   25,
                                                   0,
                                                   3.14,
                                                   25,
                                                   0,
                                                   3.14);

    /// different Z boson Pt cuts///////////////////
    hresponseAbsZRapidity_ZPt100_Zinc1jet = newTH2D("hresponseAbsZRapidity_ZPt100_Zinc1jet",
                                                    "hresponseAbsZRapidity_ZPt100_Zinc1jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseAbsFirstJetRapidity_ZPt100_Zinc1jet =
        newTH2D("hresponseAbsFirstJetRapidity_ZPt100_Zinc1jet",
                "hresponseAbsFirstJetRapidity_ZPt100_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZFirstJetRapidity_ZPt100_Zinc1jet =
        newTH2D("hresponseSumZFirstJetRapidity_ZPt100_Zinc1jet",
                "hresponseSumZFirstJetRapidity_ZPt100_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZFirstJetRapidity_ZPt100_Zinc1jet =
        newTH2D("hresponseDifZFirstJetRapidity_ZPt100_Zinc1jet",
                "hresponseDifZFirstJetRapidity_ZPt100_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    hresponseAbsZRapidity_ZPt100_Zexc1jet = newTH2D("hresponseAbsZRapidity_ZPt100_Zexc1jet",
                                                    "hresponseAbsZRapidity_ZPt100_Zexc1jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseAbsJetRapidity_ZPt100_Zexc1jet = newTH2D("hresponseAbsJetRapidity_ZPt100_Zexc1jet",
                                                      "hresponseAbsJetRapidity_ZPt100_Zexc1jet",
                                                      12,
                                                      0,
                                                      2.4,
                                                      12,
                                                      0,
                                                      2.4);
    hresponseSumZJetRapidity_ZPt100_Zexc1jet = newTH2D("hresponseSumZJetRapidity_ZPt100_Zexc1jet",
                                                       "hresponseSumZJetRapidity_ZPt100_Zexc1jet",
                                                       12,
                                                       0,
                                                       2.4,
                                                       12,
                                                       0,
                                                       2.4);
    hresponseDifZJetRapidity_ZPt100_Zexc1jet = newTH2D("hresponseDifZJetRapidity_ZPt100_Zexc1jet",
                                                       "hresponseDifZJetRapidity_ZPt100_Zexc1jet",
                                                       12,
                                                       0,
                                                       2.4,
                                                       12,
                                                       0,
                                                       2.4);

    hresponseAbsZRapidity_ZPt100_Zinc2jet = newTH2D("hresponseAbsZRapidity_ZPt100_Zinc2jet",
                                                    "hresponseAbsZRapidity_ZPt100_Zinc2jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseAbsSecondJetRapidity_ZPt100_Zinc2jet =
        newTH2D("hresponseAbsSecondJetRapidity_ZPt100_Zinc2jet",
                "hresponseAbsSecondJetRapidity_ZPt100_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZSecondJetRapidity_ZPt100_Zinc2jet =
        newTH2D("hresponseSumZSecondJetRapidity_ZPt100_Zinc2jet",
                "hresponseSumZSecondJetRapidity_ZPt100_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZSecondJetRapidity_ZPt100_Zinc2jet =
        newTH2D("hresponseDifZSecondJetRapidity_ZPt100_Zinc2jet",
                "hresponseDifZSecondJetRapidity_ZPt100_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    hresponseAbsZRapidity_ZPt100_Zexc2jet = newTH2D("hresponseAbsZRapidity_ZPt100_Zexc2jet",
                                                    "hresponseAbsZRapidity_ZPt100_Zexc2jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseAbsSecondJetRapidity_ZPt100_Zexc2jet =
        newTH2D("hresponseAbsSecondJetRapidity_ZPt100_Zexc2jet",
                "hresponseAbsSecondJetRapidity_ZPt100_Zexc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZSecondJetRapidity_ZPt100_Zexc2jet =
        newTH2D("hresponseSumZSecondJetRapidity_ZPt100_Zexc2jet",
                "hresponseSumZSecondJetRapidity_ZPt100_Zexc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZSecondJetRapidity_ZPt100_Zexc2jet =
        newTH2D("hresponseDifZSecondJetRapidity_ZPt100_Zexc2jet",
                "hresponseDifZSecondJetRapidity_ZPt100_Zexc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    hresponseAbsZRapidity_ZPt150_Zinc1jet = newTH2D("hresponseAbsZRapidity_ZPt150_Zinc1jet",
                                                    "hresponseAbsZRapidity_ZPt150_Zinc1jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseAbsFirstJetRapidity_ZPt150_Zinc1jet =
        newTH2D("hresponseAbsFirstJetRapidity_ZPt150_Zinc1jet",
                "hresponseAbsFirstJetRapidity_ZPt150_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZFirstJetRapidity_ZPt150_Zinc1jet =
        newTH2D("hresponseSumZFirstJetRapidity_ZPt150_Zinc1jet",
                "hresponseSumZFirstJetRapidity_ZPt150_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZFirstJetRapidity_ZPt150_Zinc1jet =
        newTH2D("hresponseDifZFirstJetRapidity_ZPt150_Zinc1jet",
                "hresponseDifZFirstJetRapidity_ZPt150_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    hresponseAbsZRapidity_ZPt150_Zexc1jet = newTH2D("hresponseAbsZRapidity_ZPt150_Zexc1jet",
                                                    "hresponseAbsZRapidity_ZPt150_Zexc1jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseAbsJetRapidity_ZPt150_Zexc1jet = newTH2D("hresponseAbsJetRapidity_ZPt150_Zexc1jet",
                                                      "hresponseAbsJetRapidity_ZPt150_Zexc1jet",
                                                      12,
                                                      0,
                                                      2.4,
                                                      12,
                                                      0,
                                                      2.4);
    hresponseSumZJetRapidity_ZPt150_Zexc1jet = newTH2D("hresponseSumZJetRapidity_ZPt150_Zexc1jet",
                                                       "hresponseSumZJetRapidity_ZPt150_Zexc1jet",
                                                       12,
                                                       0,
                                                       2.4,
                                                       12,
                                                       0,
                                                       2.4);
    hresponseDifZJetRapidity_ZPt150_Zexc1jet = newTH2D("hresponseDifZJetRapidity_ZPt150_Zexc1jet",
                                                       "hresponseDifZJetRapidity_ZPt150_Zexc1jet",
                                                       12,
                                                       0,
                                                       2.4,
                                                       12,
                                                       0,
                                                       2.4);

    hresponseAbsZRapidity_ZPt150_Zinc2jet = newTH2D("hresponseAbsZRapidity_ZPt150_Zinc2jet",
                                                    "hresponseAbsZRapidity_ZPt150_Zinc2jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseAbsSecondJetRapidity_ZPt150_Zinc2jet =
        newTH2D("hresponseAbsSecondJetRapidity_ZPt150_Zinc2jet",
                "hresponseAbsSecondJetRapidity_ZPt150_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZSecondJetRapidity_ZPt150_Zinc2jet =
        newTH2D("hresponseSumZSecondJetRapidity_ZPt150_Zinc2jet",
                "hresponseSumZSecondJetRapidity_ZPt150_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZSecondJetRapidity_ZPt150_Zinc2jet =
        newTH2D("hresponseDifZSecondJetRapidity_ZPt150_Zinc2jet",
                "hresponseDifZSecondJetRapidity_ZPt150_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    hresponseAbsZRapidity_ZPt150_Zexc2jet = newTH2D("hresponseAbsZRapidity_ZPt150_Zexc2jet",
                                                    "hresponseAbsZRapidity_ZPt150_Zexc2jet",
                                                    12,
                                                    0,
                                                    2.4,
                                                    12,
                                                    0,
                                                    2.4);
    hresponseAbsSecondJetRapidity_ZPt150_Zexc2jet =
        newTH2D("hresponseAbsSecondJetRapidity_ZPt150_Zexc2jet",
                "hresponseAbsSecondJetRapidity_ZPt150_Zexc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZSecondJetRapidity_ZPt150_Zexc2jet =
        newTH2D("hresponseSumZSecondJetRapidity_ZPt150_Zexc2jet",
                "hresponseSumZSecondJetRapidity_ZPt150_Zexc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZSecondJetRapidity_ZPt150_Zexc2jet =
        newTH2D("hresponseDifZSecondJetRapidity_ZPt150_Zexc2jet",
                "hresponseDifZSecondJetRapidity_ZPt150_Zexc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    hresponseAbsZRapidity_ZPt300_Zinc1jet = newTH2D("hresponseAbsZRapidity_ZPt300_Zinc1jet",
                                                    "hresponseAbsZRapidity_ZPt300_Zinc1jet",
                                                    nRapidity_ZPt300_Zinc1jet,
                                                    Rapidity_ZPt300_Zinc1jet,
                                                    nRapidity_ZPt300_Zinc1jet,
                                                    Rapidity_ZPt300_Zinc1jet);
    hresponseAbsFirstJetRapidity_ZPt300_Zinc1jet =
        newTH2D("hresponseAbsFirstJetRapidity_ZPt300_Zinc1jet",
                "hresponseAbsFirstJetRapidity_ZPt300_Zinc1jet",
                nRapidity_ZPt300_Zinc1jet,
                Rapidity_ZPt300_Zinc1jet,
                nRapidity_ZPt300_Zinc1jet,
                Rapidity_ZPt300_Zinc1jet);
    hresponseSumZFirstJetRapidity_ZPt300_Zinc1jet =
        newTH2D("hresponseSumZFirstJetRapidity_ZPt300_Zinc1jet",
                "hresponseSumZFirstJetRapidity_ZPt300_Zinc1jet",
                nRapidity_ZPt300_Zinc1jet,
                Rapidity_ZPt300_Zinc1jet,
                nRapidity_ZPt300_Zinc1jet,
                Rapidity_ZPt300_Zinc1jet);
    hresponseDifZFirstJetRapidity_ZPt300_Zinc1jet =
        newTH2D("hresponseDifZFirstJetRapidity_ZPt300_Zinc1jet",
                "hresponseDifZFirstJetRapidity_ZPt300_Zinc1jet",
                nRapidity_ZPt300_Zinc1jet,
                Rapidity_ZPt300_Zinc1jet,
                nRapidity_ZPt300_Zinc1jet,
                Rapidity_ZPt300_Zinc1jet);

    /// Azimuthal cross check//////////////////////////
    hresponseDPhiZFirstJet_ZPt150_Zinc1jet = newTH2D("hresponseDPhiZFirstJet_ZPt150_Zinc1jet",
                                                     "hresponseDPhiZFirstJet_ZPt150_Zinc1jet",
                                                     nDPhiZFirstJet_ZPt150,
                                                     DPhiZFirstJet_ZPt150,
                                                     nDPhiZFirstJet_ZPt150,
                                                     DPhiZFirstJet_ZPt150);
    hresponseDPhiZFirstJet_ZPt150_Zinc2jet = newTH2D("hresponseDPhiZFirstJet_ZPt150_Zinc2jet",
                                                     "hresponseDPhiZFirstJet_ZPt150_Zinc2jet",
                                                     nDPhiZFirstJet_ZPt150,
                                                     DPhiZFirstJet_ZPt150,
                                                     nDPhiZFirstJet_ZPt150,
                                                     DPhiZFirstJet_ZPt150);
    hresponseDPhiZFirstJet_ZPt150_Zinc3jet = newTH2D("hresponseDPhiZFirstJet_ZPt150_Zinc3jet",
                                                     "hresponseDPhiZFirstJet_ZPt150_Zinc3jet",
                                                     nDPhiZFirstJet_ZPt150,
                                                     DPhiZFirstJet_ZPt150,
                                                     nDPhiZFirstJet_ZPt150,
                                                     DPhiZFirstJet_ZPt150);
    hresponseDPhiZSecondJet_ZPt150_Zinc3jet = newTH2D("hresponseDPhiZSecondJet_ZPt150_Zinc3jet",
                                                      "hresponseDPhiZSecondJet_ZPt150_Zinc3jet",
                                                      nDPhiZSecondJet_ZPt150,
                                                      DPhiZSecondJet_ZPt150,
                                                      nDPhiZSecondJet_ZPt150,
                                                      DPhiZSecondJet_ZPt150);
    hresponseDPhiZThirdJet_ZPt150_Zinc3jet = newTH2D("hresponseDPhiZThirdJet_ZPt150_Zinc3jet",
                                                     "hresponseDPhiZThirdJet_ZPt150_Zinc3jet",
                                                     nDPhiZThirdJet_ZPt150,
                                                     DPhiZThirdJet_ZPt150,
                                                     nDPhiZThirdJet_ZPt150,
                                                     DPhiZThirdJet_ZPt150);
    hresponseDPhiFirstSecondJet_ZPt150_Zinc3jet =
        newTH2D("hresponseDPhiFirstSecondJet_ZPt150_Zinc3jet",
                "hresponseDPhiFirstSecondJet_ZPt150_Zinc3jet",
                25,
                0,
                3.14,
                25,
                0,
                3.14);
    hresponseDPhiFirstThirdJet_ZPt150_Zinc3jet =
        newTH2D("hresponseDPhiFirstThirdJet_ZPt150_Zinc3jet",
                "hresponseDPhiFirstThirdJet_ZPt150_Zinc3jet",
                25,
                0,
                3.14,
                25,
                0,
                3.14);
    hresponseDPhiSecondThirdJet_ZPt150_Zinc3jet =
        newTH2D("hresponseDPhiSecondThirdJet_ZPt150_Zinc3jet",
                "hresponseDPhiSecondThirdJet_ZPt150_Zinc3jet",
                25,
                0,
                3.14,
                25,
                0,
                3.14);

    hresponseDPhiZFirstJet_ZPt300_Zinc1jet = newTH2D("hresponseDPhiZFirstJet_ZPt300_Zinc1jet",
                                                     "hresponseDPhiZFirstJet_ZPt300_Zinc1jet",
                                                     nDPhiZFirstJet_ZPt300,
                                                     DPhiZFirstJet_ZPt300,
                                                     nDPhiZFirstJet_ZPt300,
                                                     DPhiZFirstJet_ZPt300);
    hresponseDPhiZFirstJet_ZPt300_Zinc2jet = newTH2D("hresponseDPhiZFirstJet_ZPt300_Zinc2jet",
                                                     "hresponseDPhiZFirstJet_ZPt300_Zinc2jet",
                                                     nDPhiZFirstJet_ZPt300,
                                                     DPhiZFirstJet_ZPt300,
                                                     nDPhiZFirstJet_ZPt300,
                                                     DPhiZFirstJet_ZPt300);
    hresponseDPhiZFirstJet_ZPt300_Zinc3jet = newTH2D("hresponseDPhiZFirstJet_ZPt300_Zinc3jet",
                                                     "hresponseDPhiZFirstJet_ZPt300_Zinc3jet",
                                                     nDPhiZFirstJet_ZPt300,
                                                     DPhiZFirstJet_ZPt300,
                                                     nDPhiZFirstJet_ZPt300,
                                                     DPhiZFirstJet_ZPt300);
    hresponseDPhiZSecondJet_ZPt300_Zinc3jet = newTH2D("hresponseDPhiZSecondJet_ZPt300_Zinc3jet",
                                                      "hresponseDPhiZSecondJet_ZPt300_Zinc3jet",
                                                      nDPhiZSecondJet_ZPt300,
                                                      DPhiZSecondJet_ZPt300,
                                                      nDPhiZSecondJet_ZPt300,
                                                      DPhiZSecondJet_ZPt300);
    hresponseDPhiZThirdJet_ZPt300_Zinc3jet = newTH2D("hresponseDPhiZThirdJet_ZPt300_Zinc3jet",
                                                     "hresponseDPhiZThirdJet_ZPt300_Zinc3jet",
                                                     nDPhiZThirdJet_ZPt300,
                                                     DPhiZThirdJet_ZPt300,
                                                     nDPhiZThirdJet_ZPt300,
                                                     DPhiZThirdJet_ZPt300);
    hresponseDPhiFirstSecondJet_ZPt300_Zinc3jet =
        newTH2D("hresponseDPhiFirstSecondJet_ZPt300_Zinc3jet",
                "hresponseDPhiFirstSecondJet_ZPt300_Zinc3jet",
                nDPhiJets_ZPt300,
                DPhiJets_ZPt300,
                nDPhiJets_ZPt300,
                DPhiJets_ZPt300);
    hresponseDPhiFirstThirdJet_ZPt300_Zinc3jet =
        newTH2D("hresponseDPhiFirstThirdJet_ZPt300_Zinc3jet",
                "hresponseDPhiFirstThirdJet_ZPt300_Zinc3jet",
                nDPhiJets_ZPt300,
                DPhiJets_ZPt300,
                nDPhiJets_ZPt300,
                DPhiJets_ZPt300);
    hresponseDPhiSecondThirdJet_ZPt300_Zinc3jet =
        newTH2D("hresponseDPhiSecondThirdJet_ZPt300_Zinc3jet",
                "hresponseDPhiSecondThirdJet_ZPt300_Zinc3jet",
                nDPhiJets_ZPt300,
                DPhiJets_ZPt300,
                nDPhiJets_ZPt300,
                DPhiJets_ZPt300);

    hresponseDPhiZFirstJet_ZPt150_HT300_Zinc3jet =
        newTH2D("hresponseDPhiZFirstJet_ZPt150_HT300_Zinc3jet",
                "hresponseDPhiZFirstJet_ZPt150_HT300_Zinc3jet",
                nDPhiZFirstJet_ZPt150_HT300,
                DPhiZFirstJet_ZPt150_HT300,
                nDPhiZFirstJet_ZPt150_HT300,
                DPhiZFirstJet_ZPt150_HT300);
    hresponseDPhiZSecondJet_ZPt150_HT300_Zinc3jet =
        newTH2D("hresponseDPhiZSecondJet_ZPt150_HT300_Zinc3jet",
                "hresponseDPhiZSecondJet_ZPt150_HT300_Zinc3jet",
                nDPhiZSecondJet_ZPt150_HT300,
                DPhiZSecondJet_ZPt150_HT300,
                nDPhiZSecondJet_ZPt150_HT300,
                DPhiZSecondJet_ZPt150_HT300);
    hresponseDPhiZThirdJet_ZPt150_HT300_Zinc3jet =
        newTH2D("hresponseDPhiZThirdJet_ZPt150_HT300_Zinc3jet",
                "hresponseDPhiZThirdJet_ZPt150_HT300_Zinc3jet",
                nDPhiZThirdJet_ZPt150_HT300,
                DPhiZThirdJet_ZPt150_HT300,
                nDPhiZThirdJet_ZPt150_HT300,
                DPhiZThirdJet_ZPt150_HT300);

    // different JetPt Cuts//////

    hresponseAbsZRapidity_FirstJetPt50_Zinc1jet =
        newTH2D("hresponseAbsZRapidity_FirstJetPt50_Zinc1jet",
                "hresponseAbsZRapidity_FirstJetPt50_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseAbsFirstJetRapidity_FirstJetPt50_Zinc1jet =
        newTH2D("hresponseAbsFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "hresponseAbsFirstJetRapidity_FirstJetPt50_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZFirstJetRapidity_FirstJetPt50_Zinc1jet =
        newTH2D("hresponseSumZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "hresponseSumZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZFirstJetRapidity_FirstJetPt50_Zinc1jet =
        newTH2D("hresponseDifZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                "hresponseDifZFirstJetRapidity_FirstJetPt50_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    hresponseAbsZRapidity_FirstJetPt80_Zinc1jet =
        newTH2D("hresponseAbsZRapidity_FirstJetPt80_Zinc1jet",
                "hresponseAbsZRapidity_FirstJetPt80_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseAbsFirstJetRapidity_FirstJetPt80_Zinc1jet =
        newTH2D("hresponseAbsFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "hresponseAbsFirstJetRapidity_FirstJetPt80_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZFirstJetRapidity_FirstJetPt80_Zinc1jet =
        newTH2D("hresponseSumZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "hresponseSumZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZFirstJetRapidity_FirstJetPt80_Zinc1jet =
        newTH2D("hresponseDifZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                "hresponseDifZFirstJetRapidity_FirstJetPt80_Zinc1jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    // Set Jet rapidity discriminator/////

    hresponseAbsZRapidity_DifJetRapidityl2_Zinc2jet =
        newTH2D("hresponseAbsZRapidity_DifJetRapidityl2_Zinc2jet",
                "hresponseAbsZRapidity_DifJetRapidityl2_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet =
        newTH2D("hresponseAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "hresponseAbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet =
        newTH2D("hresponseSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "hresponseSumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet =
        newTH2D("hresponseDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                "hresponseDifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    hresponseAbsZRapidity_DifJetRapiditys2_Zinc2jet =
        newTH2D("hresponseAbsZRapidity_DifJetRapiditys2_Zinc2jet",
                "hresponseAbsZRapidity_DifJetRapiditys2_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet =
        newTH2D("hresponseAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "hresponseAbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet =
        newTH2D("hresponseSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "hresponseSumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);
    hresponseDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet =
        newTH2D("hresponseDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                "hresponseDifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet",
                12,
                0,
                2.4,
                12,
                0,
                2.4);

    //--- tackmann variable ---
    for (int i(0); i < 5; i++) {
        stringstream i_str;
        i_str << "_" << i;
        tau_sum_Zinc1jet[i] = newTH1D(
            string("tau_sum_Zinc1jet" + i_str.str()).c_str(), "#Sigma#tau", "", 100, 0, 100);
        gentau_sum_Zinc1jet[i] = newTH1D(
            string("gentau_sum_Zinc1jet" + i_str.str()).c_str(), "#Sigma#tau", "", 100, 0, 100);
        hresponsetau_sum_Zinc1jet[i] =
            newTH2D(string("hresponsetau_sum_Zinc1jet" + i_str.str()).c_str(),
                    "#Sigma#tau",
                    100,
                    0,
                    100,
                    100,
                    0,
                    100);

        tau_max_Zinc1jet[i] =
            newTH1D(string("tau_max_Zinc1jet" + i_str.str()).c_str(), "max#tau", "", 100, 0, 100);
        gentau_max_Zinc1jet[i] = newTH1D(
            string("gentau_max_Zinc1jet" + i_str.str()).c_str(), "max#tau", "", 100, 0, 100);
        hresponsetau_max_Zinc1jet[i] =
            newTH2D(string("hresponsetau_max_Zinc1jet" + i_str.str()).c_str(),
                    "max#tau",
                    100,
                    0,
                    100,
                    100,
                    0,
                    100);

        tau_c_sum_Zinc1jet[i] = newTH1D(
            string("tau_c_sum_Zinc1jet" + i_str.str()).c_str(), "#Sigma#tau^{c}", "", 100, 0, 100);
        gentau_c_sum_Zinc1jet[i] = newTH1D(string("gentau_c_sum_Zinc1jet" + i_str.str()).c_str(),
                                           "#Sigma#tau^{c}",
                                           "",
                                           100,
                                           0,
                                           100);
        hresponsetau_c_sum_Zinc1jet[i] =
            newTH2D(string("hresponsetau_c_sum_Zinc1jet" + i_str.str()).c_str(),
                    "#Sigma#tau^{c}",
                    100,
                    0,
                    100,
                    100,
                    0,
                    100);

        tau_c_max_Zinc1jet[i] = newTH1D(
            string("tau_c_max_Zinc1jet" + i_str.str()).c_str(), "max#tau^{c}", "", 100, 0, 100);
        gentau_c_max_Zinc1jet[i] = newTH1D(
            string("gentau_c_max_Zinc1jet" + i_str.str()).c_str(), "max#tau^{c}", "", 100, 0, 100);
        hresponsetau_c_max_Zinc1jet[i] =
            newTH2D(string("hresponsetau_c_max_Zinc1jet" + i_str.str()).c_str(),
                    "max#tau^{c}",
                    100,
                    0,
                    100,
                    100,
                    0,
                    100);

        tau_cm_sum_Zinc1jet[i] = newTH1D(string("tau_cm_sum_Zinc1jet" + i_str.str()).c_str(),
                                         "#Sigma#tau_{cm}",
                                         "",
                                         100,
                                         0,
                                         100);
        gentau_cm_sum_Zinc1jet[i] = newTH1D(string("gentau_cm_sum_Zinc1jet" + i_str.str()).c_str(),
                                            "#Sigma#tau_{cm}",
                                            "",
                                            100,
                                            0,
                                            100);
        hresponsetau_cm_sum_Zinc1jet[i] =
            newTH2D(string("hresponsetau_cm_sum_Zinc1jet" + i_str.str()).c_str(),
                    "#Sigma#tau_{cm}",
                    100,
                    0,
                    100,
                    100,
                    0,
                    100);

        tau_cm_max_Zinc1jet[i] = newTH1D(
            string("tau_cm_max_Zinc1jet" + i_str.str()).c_str(), "max#tau_{cm}", "", 100, 0, 100);
        gentau_cm_max_Zinc1jet[i] = newTH1D(string("gentau_cm_max_Zinc1jet" + i_str.str()).c_str(),
                                            "max#tau_{cm}",
                                            "",
                                            100,
                                            0,
                                            100);
        hresponsetau_cm_max_Zinc1jet[i] =
            newTH2D(string("hresponsetau_cm_max_Zinc1jet" + i_str.str()).c_str(),
                    "max#tau_{cm}",
                    100,
                    0,
                    100,
                    100,
                    0,
                    100);

        tau_c_cm_sum_Zinc1jet[i] = newTH1D(string("tau_c_cm_sum_Zinc1jet" + i_str.str()).c_str(),
                                           "#Sigma#tau^{c}_{cm}",
                                           "",
                                           100,
                                           0,
                                           100);
        gentau_c_cm_sum_Zinc1jet[i] =
            newTH1D(string("gentau_c_cm_sum_Zinc1jet" + i_str.str()).c_str(),
                    "#Sigma#tau^{c}_{cm}",
                    "",
                    100,
                    0,
                    100);
        hresponsetau_c_cm_sum_Zinc1jet[i] =
            newTH2D(string("hresponsetau_c_cm_sum_Zinc1jet" + i_str.str()).c_str(),
                    "#Sigma#tau^{c}_{cm}",
                    100,
                    0,
                    100,
                    100,
                    0,
                    100);

        tau_c_cm_max_Zinc1jet[i] = newTH1D(string("tau_c_cm_max_Zinc1jet" + i_str.str()).c_str(),
                                           "max#tau^{c}_{cm}",
                                           "",
                                           100,
                                           0,
                                           100);
        gentau_c_cm_max_Zinc1jet[i] =
            newTH1D(string("gentau_c_cm_max_Zinc1jet" + i_str.str()).c_str(),
                    "max#tau^{c}_{cm}",
                    "",
                    100,
                    0,
                    100);
        hresponsetau_c_cm_max_Zinc1jet[i] =
            newTH2D(string("hresponsetau_c_cm_max_Zinc1jet" + i_str.str()).c_str(),
                    "max#tau^{c}_{cm}",
                    100,
                    0,
                    100,
                    100,
                    0,
                    100);
    }

    writeHistList();
}

void HistoSetZJets::writeHistList() const
{
    std::ofstream f(".histList");
    f << "# File generated by runZJets_newformat. The file .histList will be overwritten at each "
         "execution.\n\n";
    for (std::vector<TH1 *>::const_iterator it = listOfHistograms.begin();
         it != listOfHistograms.end();
         ++it) {
        f << (*it)->GetName() << "\n";
    }
}
