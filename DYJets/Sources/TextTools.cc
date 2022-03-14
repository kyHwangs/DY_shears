#include "TextTools.h"
#include "TString.h"

void createTitleVariableAnddSigma(TString variable,
                                  bool doNormalized,
                                  TString xtitle,
                                  TString &title,
                                  TString &var,
                                  TString &varUnit,
                                  TString &dSigma,
                                  TString &dSigmaUnit,
                                  bool &sepLumUnc)
{
    sepLumUnc = false;
    varUnit = "";
    dSigmaUnit = "";

    // jet multiplicity
    if (variable.Index("ZNGoodJets_Zexc") >= 0) {
        title = "Exclusive jet multiplicity";
        var = "$N_{\\text{jets}}$";
        dSigma = "$\\frac{d\\sigma}{dN_{\\text{jets}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (variable.Index("ZNGoodJets_Zinc") >= 0) {
        title = "Inclusive jet multiplicity";
        var = "$N_{\\text{jets}}$";
        dSigma = "$\\frac{d\\sigma}{dN_{\\text{jets}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    //    if (xtitle.Index("p_{T}(Z)") >= 0) {
    //        title = "$p_{\\text{T}}(\\text{Z})$";
    //        var = "$p_{\\text{T}}(Z)$ \\tiny{[GeV]}";
    //        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(Z)}$ ${\\scriptstyle
    //        [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    //    }

    if (xtitle.Index("p_{T} balance") >= 0) {
        title = "$p_{\\text{T}}^{\\text{bal}}$";
        var = "$p_{\\text{T}}^{\\text{bal}}$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}^{\\text{bal}}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }

    if (xtitle.Index("Recoil") >= 0) {
        title = "$Hadronic recoil$";
        var = "Hadronic recoil";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(hadronic recoil)}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }

    if (xtitle.Index("JZB_ptHigh") >= 0) {
        title = "JZB ($p_{\\text{T}}(\\text{Z}) > 50\\,\\text{GeV}$)";
        var = "JZB";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{d\\text{JZB}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    } else if (xtitle.Index("JZB_ptLow") >= 0) {
        title = "JZB ($p_{\\text{T}}(\\text{Z}) < 50\\,\\text{GeV}$)";
        var = "JZB";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{d\\text{JZB}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    } else if (xtitle.Index("JZB") >= 0) {
        title = "JZB";
        var = "JZB";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{d\\text{JZB}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }

    // Z pt
    if (variable.Index("ZPt_Zinc0jet") >= 0) {
        title = "$p_{\\text{T}}(\\text{Z})$  ($N_{\\text{jets}} \\geq 0$)";
        var = "$p_{\\text{T}}(\\text{Z})$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(\\text{Z})}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
        sepLumUnc = true;
    } else if (variable.Index("ZPt_Zinc1jet") >= 0) {
        title = "$p_{\\text{T}}(\\text{Z})$  ($N_{\\text{jets}} \\geq 1$)";
        var = "$p_{\\text{T}}(\\text{Z})$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(\\text{Z})}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
        sepLumUnc = true;
    }

    // pt balance
    if (xtitle.Index("VisPt_Zinc1jetQun") >= 0) {
        title = "$p_{\\text{T}}^{\\text{bal}}  ($N_{\\text{jets}} \\geq 1$)";
        var = "p_{\\text{T}}^{\\text{bal}}";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}^{\\text{bal}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("VisPt_Zinc2jetQun") >= 0) {
        title = "$p_{\\text{T}}^{\\text{bal}}  ($N_{\\text{jets}} \\geq 2$)";
        var = "p_{\\text{T}}^{\\text{bal}}";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}^{\\text{bal}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("VisPt_Zinc3jetQun") >= 0) {
        title = "$p_{\\text{T}}^{\\text{bal}}  ($N_{\\text{jets}} \\geq 3$)";
        var = "p_{\\text{T}}^{\\text{bal}}";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}^{\\text{bal}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }

    // Z Boson and Jet rapidity
    if (variable.Index("AbsZRapidity_Zexc1jet") >= 0) {
        title = "$|y_\\text{Z}|$ ($N_{\\text{jets}} = 1$)";
        var = "$|y_{\\text{Z}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{Z}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsZRapidity_Zinc1jet") >= 0) {
        title = "$|y_\\text{Z}|$ ($N_{\\text{jets}} \\geq 1$)";
        var = "$|y_{\\text{Z}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{Z}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsJetRapidity_Zexc1jet") >= 0) {
        title = "$|y_\\text{jet}|$ ($N_{\\text{jets}} = 1$)";
        var = "$|y_{\\text{jet}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{jet}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsFirstJetRapidity_Zinc1jet") >= 0) {
        title = "$|y_\\text{jet1}|$ ($N_{\\text{jets}} \\geq 1$)";
        var = "$|y_{\\text{jet1}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{jet1}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZJetRapidity_Zexc1jet") >= 0) {
        title = "$y_\\text{sum}$ ($N_{\\text{jets}} = 1$)";
        var = "$y_{\\text{sum}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZFirstJetRapidity_Zinc1jet") >= 0) {
        title = "$y_\\text{sum(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1$)";
        var = "$y_{\\text{sum(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZJetRapidity_Zexc1jet") >= 0) {
        title = "$y_\\text{diff}$ ($N_{\\text{jets}} = 1$)";
        var = "$y_{\\text{diff}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZFirstJetRapidity_Zinc1jet") >= 0) {
        title = "$y_\\text{diff(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1$)";
        var = "$y_{\\text{diff(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsZRapidity_Zinc2jet") >= 0) {
        title = "$|y_\\text{Z}|$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$|y_{\\text{Z}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{Z}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsFirstJetRapidity_Zinc2jet") >= 0) {
        title = "$|y_\\text{jet1}|$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$|y_{\\text{jet1}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{jet1}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsSecondJetRapidity_Zinc2jet") >= 0) {
        title = "$|y_\\text{jet2}|$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$|y_{\\text{jet2}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{jet2}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZFirstJetRapidity_Zinc2jet") >= 0) {
        title = "$y_\\text{sum(Z,jet1)}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$y_{\\text{sum(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZSecondJetRapidity_Zinc2jet") >= 0) {
        title = "$y_\\text{sum(Z,jet2)}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$y_{\\text{sum(Z,jet2)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet2)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumFirstSecondJetRapidity_Zinc2jet") >= 0) {
        title = "$y_\\text{sum(jet1,jet2)}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$y_{\\text{sum(jet1,jet2)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(jet1,jet2)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZFirstJetRapidity_Zinc2jet") >= 0) {
        title = "$y_\\text{diff(Z,jet1)}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$y_{\\text{diff(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZSecondJetRapidity_Zinc2jet") >= 0) {
        title = "$y_\\text{diff(Z,jet2)}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$y_{\\text{diff(Z,jet2)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet2)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifFirstSecondJetRapidity_Zinc2jet") >= 0) {
        title = "$y_\\text{diff(jet1,jet2)}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$y_{\\text{diff(jet1,jet2)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(jet1,jet2)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZTwoJetsRapidity_Zinc2jet") >= 0) {
        title = "$y_\\text{sum(Z,jet1+jet2)}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$y_{\\text{sum(Z,jet1+jet2)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet1+jet2)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZTwoJetsRapidity_Zinc2jet") >= 0) {
        title = "$y_\\text{diff(Z,jet1+jet2)}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$y_{\\text{diff(Z,jet1+jet2)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet1+jet2)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsZRapidity_ZPt150_Zinc1jet") >= 0) {
        title = "$|y_\\text{Z}|$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{Z}} \\geq 150 "
                "\\text{GeV}$)";
        var = "$|y_{\\text{Z}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{Z}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsZRapidity_ZPt300_Zinc1jet") >= 0) {
        title = "$|y_\\text{Z}|$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{Z}} \\geq 300 "
                "\\text{GeV}$)";
        var = "$|y_{\\text{Z}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{Z}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsFirstJetRapidity_ZPt150_Zinc1jet") >= 0) {
        title = "$|y_\\text{jet1}|$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{Z}} \\geq "
                "150 \\text{GeV}$)";
        var = "$|y_{\\text{jet1}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{jet1}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsFirstJetRapidity_ZPt300_Zinc1jet") >= 0) {
        title = "$|y_\\text{jet1}|$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{Z}} \\geq "
                "300 \\text{GeV}$)";
        var = "$|y_{\\text{jet1}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{jet1}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZFirstJetRapidity_ZPt150_Zinc1jet") >= 0) {
        title = "$y_\\text{sum(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{Z}} "
                "\\geq 150 \\text{GeV}$)";
        var = "$y_{\\text{sum(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZFirstJetRapidity_ZPt300_Zinc1jet") >= 0) {
        title = "$y_\\text{sum(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{Z}} "
                "\\geq 300 \\text{GeV}$)";
        var = "$y_{\\text{sum(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZFirstJetRapidity_ZPt150_Zinc1jet") >= 0) {
        title = "$y_\\text{diff(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{Z}} "
                "\\geq 150 \\text{GeV}$)";
        var = "$y_{\\text{diff(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZFirstJetRapidity_ZPt300_Zinc1jet") >= 0) {
        title = "$y_\\text{diff(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{Z}} "
                "\\geq 300 \\text{GeV}$)";
        var = "$y_{\\text{diff(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZFirstJetRapidity_FirstJetPt50_Zinc1jet") >= 0) {
        title = "$y_\\text{sum(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{jet1}} "
                "\\geq 50 \\text{GeV}$)";
        var = "$y_{\\text{sum(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZFirstJetRapidity_FirstJetPt80_Zinc1jet") >= 0) {
        title = "$y_\\text{sum(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{jet1}} "
                "\\geq 80 \\text{GeV}$)";
        var = "$y_{\\text{sum(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZFirstJetRapidity_FirstJetPt50_Zinc1jet") >= 0) {
        title = "$y_\\text{diff(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{jet1}} "
                "\\geq 50 \\text{GeV}$)";
        var = "$y_{\\text{diff(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZFirstJetRapidity_FirstJetPt80_Zinc1jet") >= 0) {
        title = "$y_\\text{diff(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1, p_{\\text{T}}^{\\text{jet1}} "
                "\\geq 80 \\text{GeV}$)";
        var = "$y_{\\text{diff(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsZRapidity_DifJetRapiditys2_Zinc2jet") >= 0) {
        title = "$|y_\\text{Z}|$ ($N_{\\text{jets}} \\geq 2, |y_{jet1} - y_{jet2}| \\leq 2 $)";
        var = "$|y_{\\text{Z}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{Z}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsZRapidity_DifJetRapidityl2_Zinc2jet") >= 0) {
        title = "$|y_\\text{Z}|$ ($N_{\\text{jets}} \\geq 2, |y_{jet1} - y_{jet2}| \\geq 2 $)";
        var = "$|y_{\\text{Z}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{Z}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsFirstJetRapidity_DifJetRapiditys2_Zinc2jet") >= 0) {
        title = "$|y_\\text{jet1}|$ ($N_{\\text{jets}} \\geq 2, |y_{jet1} - y_{jet2}| \\leq 2 $)";
        var = "$|y_{\\text{jet1}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{jet1}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("AbsFirstJetRapidity_DifJetRapidityl2_Zinc2jet") >= 0) {
        title = "$|y_\\text{jet1}|$ ($N_{\\text{jets}} \\geq 2, |y_{jet1} - y_{jet2}| \\geq 2 $)";
        var = "$|y_{\\text{jet1}}|$";
        dSigma = "$\\frac{d\\sigma}{d|y_{\\text{jet1}}|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZFirstJetRapidity_DifJetRapiditys2_Zinc2jet") >= 0) {
        title =
            "$y_\\text{sum(Z,jet1)}$ ($N_{\\text{jets}} \\geq 2, |y_{jet1} - y_{jet2}| \\leq 2 $)";
        var = "$y_{\\text{sum(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZFirstJetRapidity_DifJetRapidityl2_Zinc2jet") >= 0) {
        title =
            "$y_\\text{sum(Z,jet1)}$ ($N_{\\text{jets}} \\geq 2, |y_{jet1} - y_{jet2}| \\geq 2 $)";
        var = "$y_{\\text{sum(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{sum(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZFirstJetRapidity_DifJetRapiditys2_Zinc2jet") >= 0) {
        title =
            "$y_\\text{diff(Z,jet1)}$ ($N_{\\text{jets}} \\geq 2, |y_{jet1} - y_{jet2}| \\leq 2 $)";
        var = "$y_{\\text{diff(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZFirstJetRapidity_DifJetRapidityl2_Zinc2jet") >= 0) {
        title =
            "$y_\\text{diff(Z,jet1)}$ ($N_{\\text{jets}} \\geq 2, |y_{jet1} - y_{jet2}| \\geq 2 $)";
        var = "$y_{\\text{diff(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{dy_{\\text{diff(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("SumZFirstJetEta_Zinc1jet") >= 0) {
        title = "$\\eta_\\text{sum(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1 $)";
        var = "$\\eta_{\\text{sum(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{d\\eta_{\\text{sum(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DifZFirstJetEta_Zinc1jet") >= 0) {
        title = "$\\eta_\\text{diff(Z,jet1)}$ ($N_{\\text{jets}} \\geq 1 $)";
        var = "$\\eta_{\\text{diff(Z,jet1)}}$";
        dSigma = "$\\frac{d\\sigma}{d\\eta_{\\text{diff(Z,jet1)}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_Zinc1jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 1 $)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_Zinc2jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 2 $)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 3 $)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_ZPt150_Zinc1jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 1, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_ZPt150_Zinc2jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 2, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_ZPt150_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_ZPt300_Zinc1jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 1, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 300 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_ZPt300_Zinc2jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 2, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 300 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_ZPt300_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 300 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZSecondJet_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet2}$ ($N_{\\text{jets}} \\geq 3 $)";
        var = "$\\Delta\\phi_{\\text{Z,jet2}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet2}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZThirdJet_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet3}$ ($N_{\\text{jets}} \\geq 3 $)";
        var = "$\\Delta\\phi_{\\text{Z,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZSecondJet_ZPt150_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet2}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet2}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet2}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZThirdJet_ZPt150_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet3}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZSecondJet_ZPt300_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet2}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 300 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet2}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet2}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZThirdJet_ZPt300_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet3}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 300 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiFirstSecondJet_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{jet1,jet2}$ ($N_{\\text{jets}} \\geq 3 $)";
        var = "$\\Delta\\phi_{\\text{jet1,jet2}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{jet1,jet2}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiFirstThirdJet_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{jet1,jet3}$ ($N_{\\text{jets}} \\geq 3 $)";
        var = "$\\Delta\\phi_{\\text{jet1,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{jet1,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiSecondThirdJet_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{jet2,jet3}$ ($N_{\\text{jets}} \\geq 3 $)";
        var = "$\\Delta\\phi_{\\text{jet2,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{jet2,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiFirstSecondJet_ZPt150_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{jet1,jet2}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{jet1,jet2}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{jet1,jet2}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiFirstThirdJet_ZPt150_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{jet1,jet3}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{jet1,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{jet1,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiSecondThirdJet_ZPt150_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{jet2,jet3}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{jet2,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{jet2,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiFirstSecondJet_ZPt300_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{jet1,jet2}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 300 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{jet1,jet2}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{jet1,jet2}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiFirstThirdJet_ZPt300_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{jet1,jet3}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 300 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{jet1,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{jet1,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiSecondThirdJet_ZPt300_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{jet2,jet3}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 300 \\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{jet2,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{jet2,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZFirstJet_ZPt150_HT300_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet1}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}, H_{\\text{T}} \\geq 300 "
                "\\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet1}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet1}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZSecondJet_ZPt150_HT300_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet2}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}, H_{\\text{T}} \\geq 300 "
                "\\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet2}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet2}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    if (variable.Index("DPhiZThirdJet_ZPt150_HT300_Zinc3jet") >= 0) {
        title = "$\\Delta\\phi_\\text{Z,jet3}$ ($N_{\\text{jets}} \\geq 3, "
                "p_{\\text{T}}^{\\text{Z}} \\geq 150 \\text{GeV}, H_{\\text{T}} \\geq 300 "
                "\\text{GeV}$)";
        var = "$\\Delta\\phi_{\\text{Z,jet3}}$";
        dSigma = "$\\frac{d\\sigma}{d\\Delta\\phi_{\\text{Z,jet3}}}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    // jet pt distributions
    if (xtitle.Index("p_{T}(j_{1})") >= 0) {
        title = "$1^{\\text{st}}$ jet $p_{\\text{T}}$ ($N_{\\text{jets}} \\geq 1$)";
        var = "$p_{\\text{T}}(j_{1})$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(j_{1})}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("p_{T}(j_{2})") >= 0) {
        title = "$2^{\\text{nd}}$ jet $p_{\\text{T}}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$p_{\\text{T}}(j_{2})$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(j_{2})}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("p_{T}(j_{3})") >= 0) {
        title = "$3^{\\text{rd}}$ jet $p_{\\text{T}}$ ($N_{\\text{jets}} \\geq 3$)";
        var = "$p_{\\text{T}}(j_{3})$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(j_{3})}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("p_{T}(j_{4})") >= 0) {
        title = "$4^{\\text{th}}$ jet $p_{\\text{T}}$ ($N_{\\text{jets}} \\geq 4$)";
        var = "$p_{\\text{T}}(j_{4})$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(j_{4})}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("p_{T}(j_{5})") >= 0) {
        title = "$5^{\\text{th}}$ jet $p_{\\text{T}}$ ($N_{\\text{jets}} \\geq 5$)";
        var = "$p_{\\text{T}}(j_{5})$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(j_{5})}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("p_{T}(j_{6})") >= 0) {
        title = "$6^{\\text{th}}$ jet $p_{\\text{T}}$ ($N_{\\text{jets}} \\geq 6$)";
        var = "$p_{\\text{T}}(j_{6})$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dp_{\\text{T}}(j_{6})}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }

    // jet HT distributions
    if (xtitle.Index("H_{T}") >= 0 && title.Index("N_{jets} #geq 1") >= 0) {
        title = "$H_{\\text{T}}$ ($N_{\\text{jets}} \\geq 1$)";
        var = "$H_{\\text{T}}$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dH_{\\text{T}}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("H_{T}") >= 0 && title.Index("N_{jets} #geq 2") >= 0) {
        title = "$H_{\\text{T}}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$H_{\\text{T}}$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dH_{\\text{T}}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("H_{T}") >= 0 && title.Index("N_{jets} #geq 3") >= 0) {
        title = "$H_{\\text{T}}$ ($N_{\\text{jets}} \\geq 3$)";
        var = "$H_{\\text{T}}$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dH_{\\text{T}}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("H_{T}") >= 0 && title.Index("N_{jets} #geq 4") >= 0) {
        title = "$H_{\\text{T}}$ ($N_{\\text{jets}} \\geq 4$)";
        var = "$H_{\\text{T}}$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dH_{\\text{T}}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("H_{T}") >= 0 && title.Index("N_{jets} #geq 5") >= 0) {
        title = "$H_{\\text{T}}$ ($N_{\\text{jets}} \\geq 5$)";
        var = "$H_{\\text{T}}$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dH_{\\text{T}}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }
    if (xtitle.Index("H_{T}") >= 0 && title.Index("N_{jets} #geq 6") >= 0) {
        title = "$H_{\\text{T}}$ ($N_{\\text{jets}} \\geq 6$)";
        var = "$H_{\\text{T}}$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dH_{\\text{T}}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }

    // jet eta distributions
    if (xtitle.Index("eta(j_{1})") >= 0) {
        title = "$1^{\\text{st}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 1$)";
        var = "$\\eta(j_{1})$";
        dSigma = "$\\frac{d\\sigma}{d\\eta(j_{1})}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("eta(j_{2})") >= 0) {
        title = "$2^{\\text{nd}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$\\eta(j_{2})$";
        dSigma = "$\\frac{d\\sigma}{d\\eta(j_{2})}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("eta(j_{3})") >= 0) {
        title = "$3^{\\text{rd}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 3$)";
        var = "$\\eta(j_{3})$";
        dSigma = "$\\frac{d\\sigma}{d\\eta(j_{3})}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("eta(j_{4})") >= 0) {
        title = "$4^{\\text{th}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 4$)";
        var = "$\\eta(j_{4})$";
        dSigma = "$\\frac{d\\sigma}{d\\eta(j_{4})}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("eta(j_{5})") >= 0) {
        title = "$5^{\\text{th}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 5$)";
        var = "$\\eta(j_{5})$";
        dSigma = "$\\frac{d\\sigma}{d\\eta(j_{5})}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("eta(j_{6})") >= 0) {
        title = "$6^{\\text{th}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 6$)";
        var = "$\\eta(j_{6})$";
        dSigma = "$\\frac{d\\sigma}{d\\eta(j_{6})}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    // abs rapidity distributions
    if (xtitle.Index("|y(j_{1})|") >= 0) {
        title = "$1^{\\text{st}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 1$)";
        var = "$|y(j_{1})|$";
        dSigma = "$\\frac{d\\sigma}{d|y(j_{1})|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("|y(j_{2})|") >= 0) {
        title = "$2^{\\text{nd}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$|y(j_{2})|$";
        dSigma = "$\\frac{d\\sigma}{d|y(j_{2})|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("|y(j_{3})|") >= 0) {
        title = "$3^{\\text{rd}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 3$)";
        var = "$|y(j_{3})|$";
        dSigma = "$\\frac{d\\sigma}{d|y(j_{3})|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("|y(j_{4})|") >= 0) {
        title = "$4^{\\text{th}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 4$)";
        var = "$|y(j_{4})|$";
        dSigma = "$\\frac{d\\sigma}{d|y(j_{4})|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("|y(j_{5})|") >= 0) {
        title = "$5^{\\text{th}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 5$)";
        var = "$|y(j_{5})|$";
        dSigma = "$\\frac{d\\sigma}{d|y(j_{5})|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }
    if (xtitle.Index("|y(j_{6})|") >= 0) {
        title = "$6^{\\text{th}}$ jet $\\vert\\eta\\vert$ ($N_{\\text{jets}} \\geq 6$)";
        var = "$|y(j_{6})|$";
        dSigma = "$\\frac{d\\sigma}{d|y(j_{6})|}$";
        dSigmaUnit = "\\tiny{[\\text{pb}]}";
    }

    // dijet mass distribution
    if (xtitle.Index("M_{j_{1}j_{2}}") >= 0) {
        title = "dijet mass $M_{jj}$ ($N_{\\text{jets}} \\geq 2$)";
        var = "$M_{jj}$";
        varUnit = "\\tiny{[GeV]}";
        dSigma = "$\\frac{d\\sigma}{dM_{jj}}$";
        dSigmaUnit = "${\\scriptstyle [\\frac{\\text{pb}}{\\text{GeV}}]}$";
    }

    if (doNormalized) {
        dSigma = "$\\frac{1}{\\sigma}$ " + dSigma;
        dSigma.ReplaceAll("\\frac{\\text{pb}}", "\\frac{1}");
        dSigma.ReplaceAll("\\tiny{\\left[\\text{pb}\\right]}", "");
    }
}
