#include "jets.h"

#include <algorithm>

#include <boost/math/constants/constants.hpp>

#include "functions.h"

namespace physics
{

jets::jets(util::job::info &info, const util::options &opt)
    : Jet_pt(info.reader, "Jet_pt"),
      Jet_eta(info.reader, "Jet_eta"),
      Jet_phi(info.reader, "Jet_phi"),
      Jet_mass(info.reader, "Jet_mass"),
      Jet_jetId(info.reader, "Jet_jetId"),
      Jet_puIdDisc(info.reader, "Jet_puIdDisc"),//Iti:check
      Jet_btagCSVV2(info.reader, "Jet_btagCSVV2"),
      //Jet_hadronFlavour(info.reader, "Jet_hadronFlavour"),
      //EvtFastJetRho(info.reader, "EvtFastJetRho"),
      fixedGridRhoFastjetAll(info.reader, "fixedGridRhoFastjetAll"), //Iti:check
      GenJet_pt(info.init_optional_branch<decltype(GenJet_pt)>("GenJet_pt")),
      GenJet_eta(info.init_optional_branch<decltype(GenJet_eta)>("GenJet_eta")),
      GenJet_phi(info.init_optional_branch<decltype(GenJet_eta)>("GenJet_phi")),
      GenJet_mass(info.init_optional_branch<decltype(GenJet_eta)>("GenJet_mass"))
{
    configure(opt);
    m_JetResolution =
        new JME::JetResolution("EfficiencyTables/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt");
    m_JetResolutionScaleFactor =
        new JME::JetResolutionScaleFactor("EfficiencyTables/Spring16_25nsV10_MC_SF_AK4PFchs.txt");
    m_JetParameters = new JME::JetParameters();
}

void jets::configure(const util::options &opt)
{

    const YAML::Node node = opt.config["jets"];
    util::set_value_safe(node, _pt_cut, "pt", "jet pt cut", [](double val) { return val >= 0; });
    util::set_value_safe(node, _y_cut, "rapidity", "jet rapidity cut", [](double val) { return val > 0; });
    util::set_value_safe(node, _pumva_cut, "pu mva", "jet PU MVA cut", [](double val) {
        return val >= -1 && val < 1;
    });
    util::set_value_safe(
        node, _deltar_cut, "lepton delta r", "jet-lepton Delta R cut", [](double val) {
            return val > 0;
        });
    util::set_value_safe(node, _jer_smearing, "JER smearing", "JER smearing toggle");
}

void jets::declare_histograms(util::histo_set &h)
{
    const double pi = boost::math::constants::pi<double>();

    h.declare("nJets", "Jet multiplicity (excl.)", 7, -0.5, 6.5);
    h.declare("nJetsIncl", "Jet multiplicity (incl.)", 7, -0.5, 6.5);
    h.declare("jetPuMva", "Jet PU variable from MVA", 40, -1, 1);
    h.declare("jetbdisc", "Jet bdisc variable ", 40, -1, 1);
    h.declare("jetPt", "Jet pt", 40, 0, 200);
    h.declare("jetEta", "Jet eta", 24, -2.4, 2.4);
    h.declare("jetPhi", "Jet phi", 24, -pi, pi);
    h.declare("jetDr_gen", "Jet DeltaR", 24, -pi, pi);
}

std::vector<jet> jets::getGen()
{
    return getGen(_pt_cut, _y_cut);
}

std::vector<jet> jets::getGen(double ptmin, double rapmax)
{
    std::vector<jet> Gjets;
    for (unsigned i = 0; i < GenJet_pt->GetSize(); ++i) {
        jet j;
        if (GenJet_pt->At(i) < ptmin) {
            continue;
        }
        j.v.SetPtEtaPhiM(GenJet_pt->At(i), GenJet_eta->At(i), GenJet_phi->At(i), GenJet_mass->At(i));
        if (std::abs(j.v.Rapidity()) > rapmax) {
            continue;
        }
        Gjets.push_back(j);
    }
    return Gjets;
}

std::vector<jet> jets::get(bool isdata, const std::vector<lepton> &leptons, double ptcut)
{
    if (ptcut < 0) {
        ptcut = _pt_cut;
    }

    std::vector<jet> jets;
    for (unsigned i = 0; i < Jet_pt.GetSize(); ++i) {
        jet j;
        if (Jet_puIdDisc[i] < _pumva_cut || Jet_jetId[i] <= 0) {
            continue;
        }
        j.v.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
        m_JetParameters->setJetPt(j.v.Pt());
        m_JetParameters->setJetEta(j.v.Eta());
        m_JetParameters->setRho(*fixedGridRhoFastjetAll);
        jetResolution = m_JetResolution->getResolution(*m_JetParameters);
        jetSF = m_JetResolutionScaleFactor->getScaleFactor(*m_JetParameters, m_Variation);

        double deltarjjmin =0.2;
        bool matched =false;
        if (!isdata && _jer_smearing) {
            float smearFactor = 1.0;
            std::vector<jet> gjet = getGen(10 , 5);
            veto(gjet, leptons);
            for (const auto &gj: gjet) {
                float deltarjj = gj.v.DeltaR(j.v);
                float dPt = abs(j.v.Pt() - gj.v.Pt());
                if (deltarjj < deltarjjmin && dPt < (3 * j.v.Pt() * jetResolution)) {
                    deltarjjmin = deltarjj;
                    smearFactor = 1.0 + (jetSF - 1.0) * (j.v.Pt() - gj.v.Pt()) / j.v.Pt();
                    matched = true;
                }
           }
           if (!matched) {
                TRandom3 *random = new TRandom3(0);
                smearFactor =
                    1.0 +
                    random->Gaus(0.0, jetResolution) * sqrt(std::max(pow(jetSF, 2) - 1.0, 0.0));
                delete random;
            }
            float oldJetPt = j.v.Pt();
            float newJetPt = oldJetPt * smearFactor;
            j.v.SetPtEtaPhiM(newJetPt, j.v.Eta(), j.v.Phi(), j.v.M() * newJetPt / oldJetPt);
        }

        if (std::abs(j.v.Rapidity()) > _y_cut) {
            continue;
        }
        j.raw_v = j.v;
        j.id = Jet_jetId[i];
        j.puMva = Jet_puIdDisc[i];
        j.bdisc = Jet_btagCSVV2[i];
        j.hadflav = 0;

        if (j.v.Pt() < ptcut) continue;

        jets.push_back(j);
    }
        std::sort(jets.begin(), jets.end(), [](const jet &lhs, const jet &rhs)
            {
                return lhs.v.Pt() > rhs.v.Pt();
            }
);
    return jets;
}

void jets::veto(std::vector<jet> &jets, const std::vector<lepton> &leptons) const
{
    jets.erase(std::remove_if(jets.begin(),
                              jets.end(),
                              [&](const jet &j) {
                                  for (const lepton &l : leptons) {
                                      if (deltaR(j.v, l.v) < _deltar_cut) {
                                          return true;
                                      }
                                  }
                                  return false;
                              }),
               jets.end());
}

void jets::fill(util::histo_set &h,
                const std::string &tag,
                const std::vector<jet> &jets,
                const weights &w)
{
    h.fill("nJets", tag, jets.size(), w.global_weight());
    for (std::size_t njets = 0; njets <= jets.size(); ++njets) {
        h.fill("nJetsIncl", tag, njets, w.global_weight());
    }
    for (const jet &j : jets) {
        h.fill("jetPt", tag, j.v.Pt(), w.global_weight());
        h.fill("jetEta", tag, j.v.Eta(), w.global_weight());
        h.fill("jetPhi", tag, j.v.Phi(), w.global_weight());
        h.fill("jetPuMva", tag, j.puMva, w.global_weight());
        h.fill("jetbdisc", tag, j.bdisc, w.global_weight());
    }
    if (jets.size() > 0) {
        const jet &j = jets[0];
        h.fill("jetPt", "leading_" + tag, j.v.Pt(), w.global_weight());
        h.fill("jetEta", "leading_" + tag, j.v.Eta(), w.global_weight());
        h.fill("jetPhi", "leading_" + tag, j.v.Phi(), w.global_weight());
    }
    if (jets.size() > 1) {
        const jet &j = jets[1];
        h.fill("jetPt", "subleading_" + tag, j.v.Pt(), w.global_weight());
        h.fill("jetEta", "subleading_" + tag, j.v.Eta(), w.global_weight());
        h.fill("jetPhi", "subleading_" + tag, j.v.Phi(), w.global_weight());
    }
    if (jets.size() > 2) {
        const jet &j = jets[2];
        h.fill("jetPt", "third_" + tag, j.v.Pt(), w.global_weight());
        h.fill("jetEta", "third_" + tag, j.v.Eta(), w.global_weight());
        h.fill("jetPhi", "third_" + tag, j.v.Phi(), w.global_weight());
    }
    if (jets.size() > 3) {
        const jet &j = jets[3];
        h.fill("jetPt", "fourth_" + tag, j.v.Pt(), w.global_weight());
        h.fill("jetEta", "fourth_" + tag, j.v.Eta(), w.global_weight());
        h.fill("jetPhi", "fourth_" + tag, j.v.Phi(), w.global_weight());
    }
}
}
