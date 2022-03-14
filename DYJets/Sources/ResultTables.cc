#include "ResultTables.h"
#include "ConfigVJets.h"
#include "TextTools.h"
#include "Uncertainties.h"
#include "functions.h"
#include <assert.h>
#include <fstream>
#include <iostream>

extern ConfigVJets cfg;

TString cellWrapper(const TString &id, const TString &value)
{
    bool debug = false;
    if (debug) {
        return id;
    } else {
        return value;
    }
}

void createTable(TString outputFilePath,
                 TString lepSel,
                 TString variable,
                 bool doNormalized,
                 TH1 *hCombination,
                 vector<TH2 *> &covuxaxb,
                 bool withLERS,
                 bool forceNoLumiSep)
{

    assert(covuxaxb.size() == kUncCnt);

    //--- print out break down of errors ---

    TString title = hCombination->GetTitle();
    int nBins = hCombination->GetNbinsX();
    TString var;
    TString varUnit;
    TString dSigma;
    TString dSigmaUnit;
    TString xtitle = hCombination->GetXaxis()->GetTitle();
    bool sepLumiUnc = false;
    createTitleVariableAnddSigma(
        variable, doNormalized, xtitle, title, var, varUnit, dSigma, dSigmaUnit, sepLumiUnc);

    bool withJEC = true;
    if (variable.EndsWith("Zinc0jet")) withJEC = 0;

    if (forceNoLumiSep) sepLumiUnc = false;

    title.ReplaceAll("Exclusive", "exclusive");
    title.ReplaceAll("Inclusive", "inclusive");
    title.ReplaceAll("Hadronic", "hadronic");

    int verbosity = cfg.getI("verbosity");

    if (verbosity > 0) {
        std::cout << "Title: " << title << std::endl;
        std::cout << "Var: " << var << std::endl;
        std::cout << "dSig: " << dSigma << std::endl;
        std::cout << "nBins: " << nBins << std::endl;
    }

    TString table = "\\begin{table}[htb!]\n\\begin{center}\n";
    table += "\\caption{Differential cross section in " + title;
    table += " and break down of the systematic uncertainties for the ";
    if (lepSel == "DMu") table += "muon decay channel.}\n";
    if (lepSel == "DE") table += "electron decay channel.}\n";
    if (lepSel == "") table += "combination of both decay channels.}\n";

    table += "\\scriptsize{\n";
    if (sepLumiUnc) {
        table += "\\begin{tabular}{c|ccc|cccccc";
    } else {
        table += "\\begin{tabular}{c|cc|ccccccc";
    }
    if (withLERS) table += "cc";
    if (withJEC) table += "cc";
    table += "}\n";
    TString ul;
    table += var + " & " + dSigma + " & ";
    // FIXME: need to separtate unit.....
    ul += varUnit + " & " + dSigmaUnit + " & ";
    if (sepLumiUnc) {
        table += "\\tiny{Subtot. Unc} & \\tiny{Lumi} & ";
        ul += "\\tiny{[\\%]} & \\tiny{[\\%]} & ";
    } else {
        table += "\\tiny{Tot. Unc} & ";
        ul += "\\tiny{[\\%]}& ";
    }
    table += "\\tiny{Stat} & ";
    ul += "\\tiny{[\\%]} & ";
    if (withJEC) {
        table += "\\tiny{JES} & \\tiny{JER} & ";
        ul += "\\tiny{[\\%]} & \\tiny{[\\%]} & ";
    }
    table += "\\tiny{Eff} & ";
    ul += "\\tiny{[\\%]} & ";
    if (!sepLumiUnc) {
        table += "\\tiny{Lumi} & ";
        ul += "\\tiny{[\\%]} & ";
    }
    table += "\\tiny{Bkg} & ";
    ul += "\\tiny{[\\%]} & ";
    if (withLERS) {
        table += "\\tiny{LES} & \\tiny{LER} & ";
        ul += "\\tiny{[\\%]} & \\tiny{[\\%]} & ";
    }
    table += "\\tiny{PU} & \\tiny{Unf model} & \\tiny{Unf stat}\\\\\n";
    ul += "\\tiny{[\\%]} & \\tiny{[\\%]} & \\tiny{[\\%]}";

    table += ul + "\\\\\n\\hline\n";

    int start = 1;
    /*if (title.Index("multiplicity", 0, TString::ECaseCompare::kIgnoreCase) >= 0) {
        start = 2;
        //nBins--;
    }*/
    if (title.Index("jet $p_{\\text{T}}$", 0, TString::ECaseCompare::kIgnoreCase) >= 0) start = 3;

    for (int i = start; i <= nBins; ++i) {
        double xs = hCombination->GetBinContent(i);
        TString numbers;
        if (title.Index("exclusive jet multiplicity", 0, TString::ECaseCompare::kIgnoreCase) >= 0) {
            numbers.Form("= %d", i - 1);
        } else if (title.Index(
                       "inclusive jet multiplicity", 0, TString::ECaseCompare::kIgnoreCase) >= 0) {
            numbers.Form("$\\geq$ %d", i - 1);
        } else {
            numbers.Form("$%g \\ -\\ %g$",
                         hCombination->GetBinLowEdge(i),
                         hCombination->GetBinLowEdge(i + 1));
        }
        table += cellWrapper("bin", numbers) + " & ";
        // numbers.Form("%#.3g", xs);
        // table += numbers + " & ";
        //	double totUnc = sqrt(covuxaxb[0]->GetBinContent(i,i) +
        //covxaxbSyst->GetBinContent(i,i));
        double totUnc =
            sqrt(covuxaxb[0]->GetBinContent(i, i) + covuxaxb[kTotSys]->GetBinContent(i, i));
        if (sepLumiUnc) {
            totUnc = sqrt(pow(totUnc, 2) - pow(covuxaxb[kLumi]->GetBinContent(i, i), 2));
        }
        std::string sXs, dummy;
        pground(xs, totUnc, sXs, dummy);
        table += cellWrapper("$\\sigma$", sXs) + " & ";
        // total uncertainty
        numbers.Form("%#.2g", totUnc * 100. / xs);
        // numbers.Form("%#.2g", sqrt(covuxaxb[0]->GetBinContent(i,i) +
        //                           covuxaxb[1]->GetBinContent(i,i) +
        //                           covuxaxb[2]->GetBinContent(i,i) +
        //                           covuxaxb[3]->GetBinContent(i,i) +
        //                           covuxaxb[4]->GetBinContent(i,i) +
        //                           covuxaxb[5]->GetBinContent(i,i) +
        //                           covuxaxb[6]->GetBinContent(i,i) +
        //                           covuxaxb[7]->GetBinContent(i,i) +
        //                           covuxaxb[8]->GetBinContent(i,i) )*100./xs);
        table += cellWrapper("totUnc", numbers) + " & ";
        if (sepLumiUnc) {
            numbers.Form("%#.2g", sqrt(covuxaxb[kLumi]->GetBinContent(i, i)) * 100. / xs);
            table += cellWrapper("lumi", numbers) + " & ";
        }
        // stat uncertainty
        numbers.Form("%#.2g", sqrt(covuxaxb[kStat]->GetBinContent(i, i)) * 100. / xs);
        table += cellWrapper("statUnc", numbers) + " & ";
        if (withJEC) {
            // JES uncertainty
            numbers.Form("%#.2g", sqrt(covuxaxb[kJES]->GetBinContent(i, i)) * 100. / xs);
            table += cellWrapper("jes", numbers) + " & ";
            // JER uncertainty
            numbers.Form("%#.2g", sqrt(covuxaxb[kJER]->GetBinContent(i, i)) * 100. / xs);
            table += cellWrapper("jer", numbers) + " & ";
        }
        // SF uncertinaty
        numbers.Form("%#.2g", sqrt(covuxaxb[kSF]->GetBinContent(i, i)) * 100. / xs);
        table += cellWrapper("sf", numbers) + " & ";
        // Lumi uncertainty
        if (!sepLumiUnc) {
            numbers.Form("%#.2g", sqrt(covuxaxb[kLumi]->GetBinContent(i, i)) * 100. / xs);
            table += cellWrapper("lumi", numbers) + " & ";
        }
        // XSec (Bgnd) uncertainty
        numbers.Form("%#.2g", sqrt(covuxaxb[kXsec]->GetBinContent(i, i)) * 100. / xs);
        table += cellWrapper("bgnd", numbers) + " & ";
        if (withLERS) {
            // LES uncertainty
            numbers.Form("%#.2g", sqrt(covuxaxb[kLES]->GetBinContent(i, i)) * 100. / xs);
            table += cellWrapper("les", numbers) + " & ";
            // LER uncertainty
            numbers.Form("%#.2g", sqrt(covuxaxb[kLER]->GetBinContent(i, i)) * 100. / xs);
            table += cellWrapper("ler", numbers) + " & ";
        }
        // PU uncertainty
        numbers.Form("%#.2g", sqrt(covuxaxb[kPU]->GetBinContent(i, i)) * 100. / xs);
        table += cellWrapper("pu", numbers) + " & ";
        // Unf uncertainty
        if (covuxaxb[kUnfSys])
            numbers.Form("%#.2g", sqrt(covuxaxb[kUnfSys]->GetBinContent(i, i)) * 100. / xs);
        else
            numbers = "-";
        table += cellWrapper("unf. sys.", numbers) + " & ";
        // MC stat uncertainty
        numbers.Form("%#.2g", sqrt(covuxaxb[kUnfStat]->GetBinContent(i, i)) * 100. / xs);
        table += cellWrapper("unf. stat", numbers) + " \\\\\n";
    }

    table += "\\end{tabular}}\n";
    table += "\\label{tab:" + (lepSel.Length() ? lepSel : TString("comb")) + variable + "}\n";
    table += "\\end{center}\\end{table}\n";
    std::ofstream out(outputFilePath + ".tex");
    out << table;
    out.close();
    if (verbosity > 0) std::cout << table << std::endl;
}
