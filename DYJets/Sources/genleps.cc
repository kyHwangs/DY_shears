#include "genleps.h"

#include <boost/math/constants/constants.hpp>

#include "logging.h"
//#include "RoccoR.h"

namespace physics
{

genleps::genleps(util::job::info &info, const util::options &opt, util::histo_set &h)
    : GenDressedLepton_pt(info.init_optional_branch<decltype(GenDressedLepton_pt)>("GenDressedLepton_pt")),
      GenDressedLepton_eta(info.init_optional_branch<decltype(GenDressedLepton_eta)>("GenDressedLepton_eta")),
      GenDressedLepton_phi(info.init_optional_branch<decltype(GenDressedLepton_phi)>("GenDressedLepton_phi")),
      //GenLepE(info.reader, "GLepDr01E"),
      GenDressedLepton_mass(info.init_optional_branch<decltype(GenDressedLepton_mass)>("GenDressedLepton_mass")),
      GenDressedLepton_pdgId(info.init_optional_branch<decltype(GenDressedLepton_pdgId)>("GenDressedLepton_pdgId")),
      //GenLepPrompt(info.reader,"GenLepPrompt"),
      GenPart_statusFlags(info.init_optional_branch<decltype(GenPart_statusFlags)>("GenPart_statusFlags")), //Iti: check
      GenDressedLepton_hasTauAnc(info.init_optional_branch<decltype(GenDressedLepton_hasTauAnc)>("GenDressedLepton_hasTauAnc")),
      //GenLepSt(info.reader, "GLepDr01St"),
      GenPart_status(info.init_optional_branch<decltype(GenPart_status)>("GenPart_status")) //Iti: check
      //LHEZChild1Id(info.reader, "LHEZChild1id"),
      //LHEZChild2Id(info.reader, "LHEZChild2id"),
      //LHEZChild1Px(info.reader, "LHEChild1Px"),
      //LHEZChild2Px(info.reader, "LHEChild2Px"),
      //LHEZChild1Py(info.reader, "LHEChild1Py"),
      //LHEZChild2Py(info.reader, "LHEChild2Py")


{
    configure(opt);

    const double pi = boost::math::constants::pi<double>();

    h.declare("genLepPt", "Gen lepton pt;Gen lepton p_{T} [GeV]", 50, 0, 200);
    h.declare("genLepEta", "Gen lepton eta;Gen lepton #eta", 24, -2.4, 2.4);
    h.declare("genLepPhi", "Gen lepton phi;Gen lepton #phi", 24, -pi, pi);
}

void genleps::configure(const util::options &opt)
{
    const YAML::Node node = opt.config["muons"];
    util::set_value_safe(node, _pt_cut, "pt", "muon pt cut", [](double val) { return val >= 0; });
    util::set_value_safe(node, _eta_cut, "eta", "muon eta cut", [](double val) { return val > 0; });
    
}

std::vector<lepton> genleps::get()
{
//    std::cout<<"________"<<std::endl;
    std::vector<lepton> genleps;

    for (unsigned i = 0; i < GenDressedLepton_pt->GetSize(); ++i) {
        lepton l;

        if (std::abs(GenDressedLepton_eta->At(i)) > _eta_cut || GenDressedLepton_hasTauAnc->At(i) ) {
            continue;
        }

        l.v.SetPtEtaPhiM(GenDressedLepton_pt->At(i), GenDressedLepton_eta->At(i), GenDressedLepton_phi->At(i), GenDressedLepton_mass->At(i));
        l.raw_v = l.v;
        l.charge = GenDressedLepton_pdgId->At(i)/std::abs(GenDressedLepton_pdgId->At(i));
        l.pdgid = GenDressedLepton_pdgId->At(i);

        if (l.v.Pt() < _pt_cut) {
           continue;
        }
        genleps.push_back(l);
    }
    std::sort(genleps.begin(), genleps.end(), [](const lepton &lhs, const lepton &rhs)
        {
            return lhs.v.Pt() > rhs.v.Pt();
        }
    );

    return genleps;
}
void genleps::fill(util::histo_set &h,
                 const std::string &tag,
                 const std::vector<lepton> &genleps,
                 const weights &w)
{
    for (const lepton &mu : genleps) {
        h.fill("genLepPt", tag, mu.v.Pt(), w.gen_weight());
        h.fill("genLepEta", tag, mu.v.Eta(), w.gen_weight());
        h.fill("genLepPhi", tag, mu.v.Phi(), w.gen_weight());
    }
    if (genleps.size() > 0) {
        const lepton &mu = genleps[0];
        h.fill("genLepPt", "leading_" + tag, mu.v.Pt(), w.gen_weight());
        h.fill("genLepEta", "leading_" + tag, mu.v.Eta(), w.gen_weight());
        h.fill("genLepPhi", "leading_" + tag, mu.v.Phi(), w.gen_weight());
    }
    if (genleps.size() > 1) {
        const lepton &mu = genleps[1];
        h.fill("genLepPt", "subleading_" + tag, mu.v.Pt(), w.gen_weight());
        h.fill("genLepEta", "subleading_" + tag, mu.v.Eta(), w.gen_weight());
        h.fill("genLepPhi", "subleading_" + tag, mu.v.Phi(), w.gen_weight());
    }
    if (genleps.size() > 2) {
        const lepton &mu = genleps[2];
        h.fill("genLepPt", "third_" + tag, mu.v.Pt(), w.gen_weight());
        h.fill("genLepEta", "third_" + tag, mu.v.Eta(), w.gen_weight());
        h.fill("genLepPhi", "third_" + tag, mu.v.Phi(), w.gen_weight());
    }
    if (genleps.size() > 3) {
        const lepton &mu = genleps[3];
        h.fill("genLepPt", "fourth_" + tag, mu.v.Pt(), w.gen_weight());
        h.fill("genLepEta", "fourth_" + tag, mu.v.Eta(), w.gen_weight());
        h.fill("genLepPhi", "fourth_" + tag, mu.v.Phi(), w.gen_weight());
    }
}
} // namespace physics
