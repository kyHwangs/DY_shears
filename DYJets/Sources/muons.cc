#include "muons.h"

#include <boost/math/constants/constants.hpp>

#include "logging.h"
#include "RoccoR.h"
#include "TRandom.h"

namespace physics
{

muons::muons(util::job::info &info, const util::options &opt, util::histo_set &h)
    : Muon_pt(info.reader, "Muon_pt"),
      Muon_eta(info.reader, "Muon_eta"),
      Muon_phi(info.reader, "Muon_phi"),
      Muon_mass(info.reader, "Muon_mass"),
      Muon_charge(info.reader, "Muon_charge"),
      Muon_pfRelIso04_all(info.reader, "Muon_pfRelIso04_all"), //Iti: check
      Muon_nTrackerLayers(info.reader, "Muon_nTrackerLayers"),
      Muon_pfIsoId(info.reader, "Muon_pfIsoId"), //Iti: check
      Muon_tightId(info.reader, "Muon_tightId"),
      Muon_mediumId(info.reader, "Muon_mediumId"),
      Muon_looseId(info.reader, "Muon_looseId")
{
    configure(opt);

    const double pi = boost::math::constants::pi<double>();

    h.declare("muPt", "Muon pt;Muon p_{T} [GeV]", 50, 0, 200);
    h.declare("muEta", "Muon eta;Muon #eta", 24, -2.4, 2.4);
    h.declare("muPhi", "Muon phi;Muon #phi", 24, -pi, pi);
    h.declare("muIso", "Muon relative isolation;Muon relative isolation", 40, 0, 1);
}

void muons::configure(const util::options &opt)
{
    const YAML::Node node = opt.config["muons"];
    util::set_value_safe(node, _pt_cut, "pt", "muon pt cut", [](double val) { return val >= 0; });
    util::set_value_safe(node, _eta_cut, "eta", "muon eta cut", [](double val) { return val > 0; });
//    util::set_value_safe(
//        node, _iso_cut, "isolation", "muon isolation cut", [](double val) { return val >= 0; });
    util::set_value_safe(node, _id_sf_enabled, "use id scale factors", "id scale factors toggle");
    util::set_value_safe(
        node, _iso_sf_enabled, "use isolation scale factors", "isolation scale factors toggle");
    util::set_value_safe(
        node, _trk_sf_enabled, "use tracking scale factors", "tracking scale factors toggle");
    util::set_value_safe(
        node, _roccor_enabled, "use rochester correction", "rochester correction toggle");
    if (_roccor_enabled) {
        std::string roccor_dir = "rcdata.2016.v3";
        if (node["rochester correction path"]) {
            roccor_dir = node["rochester correction path"].as<std::string>();
        }
        roccor_dir = "EfficiencyTables/" + roccor_dir;
        _roccor = std::make_shared<RoccoR>(roccor_dir);//hardcoded now, since there is a problem with yml
    }

    if (node["id"]) {
        std::string id = node["id"].as<std::string>();
        if (id == "loose") {
            _id_cut = muons::id::loose;
        } else if (id == "medium") {
            _id_cut = muons::id::medium;
        } else if (id == "tight") {
            _id_cut = muons::id::tight;
        } else {
            throw std::invalid_argument("Unknown muon id: \"" + id + "\"");
        }
    }
    if (node["iso"]) {
        std::string iso = node["iso"].as<std::string>();
        if (iso == "loose") {
            _iso_cut = muons::iso::loose;
        } else if (iso == "medium") {
            _iso_cut = muons::iso::medium;
        } else if (iso == "tight") {
            _iso_cut = muons::iso::tight;
        } else if (iso == "veryloose") {
            _iso_cut = muons::iso::veryloose;
        } else {
            throw std::invalid_argument("Unknown muon Isoid: \"" + iso + "\"");
        }
    }
}

std::vector<lepton> muons::get(bool isdata, std::mt19937 &rng, std::vector<lepton> gl, int &nVetoMuons)
{
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    nVetoMuons =0 ;
    std::vector<lepton> muons;
    for (unsigned i = 0; i < Muon_pt.GetSize(); ++i) {
       lepton l;
        if (std::abs(Muon_eta[i]) > _eta_cut) {
            continue;
        }
        l.v.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
        l.raw_v = l.v;
        l.charge = Muon_charge[i];
        l.iso = Muon_pfIsoId[i];
        l.id = Muon_looseId[i];
        l.pdgid = 13;
        
#ifdef DEBUG_PRINTOUT
        l.tkLayerCnt = Muon_nTrackerLayers[i];
#endif // DEBUG_PRINTOUT

        switch (_id_cut) {
        case id::loose:
            l.passes_id = Muon_looseId[i];
            break;
        case id::medium:
            l.passes_id = Muon_mediumId[i];
            break;
        case id::tight:
            l.passes_id = Muon_tightId[i];
            break;
        }
        if (!l.passes_id) {
            continue;
        }

        switch (_iso_cut) {
        case iso::veryloose:
            l.passes_iso = (Muon_pfIsoId[i] >= 1);
            break;
        case iso::loose:
            l.passes_iso = (Muon_pfIsoId[i] >= 2);
            break;
        case iso::medium:
            l.passes_iso = (Muon_pfIsoId[i] >= 3);
            break;
        case iso::tight:
            l.passes_iso = (Muon_pfIsoId[i] >= 4);
            break;
        }
        if (!l.passes_iso) {
            continue;
        }
        
        if (_roccor_enabled) {
            if (isdata) {
                l.v *= _roccor->kScaleDT(l.charge, l.v.Pt(), l.v.Eta(), l.v.Phi(), 0, 0);
#ifdef DEBUG_PRINTOUT
                l.fnUsed = 2;
                l.gllPt = 0;
#endif // DEBUG_PRINTOUT
            } else {
                lepton gll;double drmin=99.;bool match =false;
                for (auto &v : gl) {
                if(fabs(v.pdgid)==13 &&v.v.DeltaR(l.v)<0.1 && v.v.DeltaR(l.v)<drmin){
                gll=v;match=true;
                drmin=v.v.DeltaR(l.v);
                }}
                if(match){
                l.v *= _roccor->kScaleFromGenMC(l.charge,
                                                 l.v.Pt(),
                                                 l.v.Eta(),
                                                 l.v.Phi(),
                                                 Muon_nTrackerLayers[i],
                                                 gll.v.Pt(),
//                                                 uniform(rng),
						gRandom->Rndm(),
  //0.75,
                                                 0,
                                                 0);
#ifdef DEBUG_PRINTOUT
                l.fnUsed = 0;
                l.gllPt = gll.v.Pt();
#endif // DEBUG_PRINTOUT
                }
                else{
                l.v *= _roccor->kScaleAndSmearMC(l.charge,
                                                 l.v.Pt(),
                                                 l.v.Eta(),
                                                 l.v.Phi(),
                                                 Muon_nTrackerLayers[i],
//  0.75,                                          
//  0.75,                                          

//                                                 uniform(rng),
//                                                 uniform(rng),
						gRandom->Rndm(),
						gRandom->Rndm(),
                                                 0,
                                                 0);
#ifdef DEBUG_PRINTOUT
                l.fnUsed = 1;
#endif // DEBUG_PRINTOUT
		}
            }
        }
        if (l.v.Pt() < _pt_cut) {
            continue;
        }
        muons.push_back(l);
    }
    if (_roccor_enabled) {
        std::sort(muons.begin(), muons.end(), [](const lepton &lhs, const lepton &rhs)
            {
                return lhs.v.Pt() > rhs.v.Pt();
            }
        );
    }
    return muons;
}

void muons::apply_sf(weights &w, const std::vector<lepton> &muons, const util::tables &tab) const
{
    if (w.ismc()) {
        for (const lepton &mu : muons) {
            if (_id_sf_enabled) {
                w.use_weight(tab.at("muon id").getEfficiency(mu.v.Pt(), std::abs(mu.v.Eta())));
            }
            if (_iso_sf_enabled) {
                w.use_weight(tab.at("muon isolation")
                                .getEfficiency(mu.v.Pt(), std::abs(mu.v.Eta())));
            }
            if (_trk_sf_enabled) {
                w.use_weight(
                    tab.at("muon tracking").getEfficiency(mu.raw_v.Pt(), mu.raw_v.Eta()));
            }
        }
    }
}

void muons::fill(util::histo_set &h,
                 const std::string &tag,
                 const std::vector<lepton> &muons,
                 const weights &w)
{
    for (const lepton &mu : muons) {
        h.fill("muPt", tag, mu.v.Pt(), w.global_weight());
        h.fill("muEta", tag, mu.v.Eta(), w.global_weight());
        h.fill("muPhi", tag, mu.v.Phi(), w.global_weight());
        h.fill("muIso", tag, mu.iso, w.global_weight());
    }
    if (muons.size() > 0) {
        const lepton &mu = muons[0];
        h.fill("muPt", "leading_" + tag, mu.v.Pt(), w.global_weight());
        h.fill("muEta", "leading_" + tag, mu.v.Eta(), w.global_weight());
        h.fill("muPhi", "leading_" + tag, mu.v.Phi(), w.global_weight());
        h.fill("muIso", "leading_" + tag, mu.iso, w.global_weight());
    }
    if (muons.size() > 1) {
        const lepton &mu = muons[1];
        h.fill("muPt", "subleading_" + tag, mu.v.Pt(), w.global_weight());
        h.fill("muEta", "subleading_" + tag, mu.v.Eta(), w.global_weight());
        h.fill("muPhi", "subleading_" + tag, mu.v.Phi(), w.global_weight());
        h.fill("muIso", "subleading_" + tag, mu.iso, w.global_weight());
    }
    if (muons.size() > 2) {
        const lepton &mu = muons[2];
        h.fill("muPt", "third_" + tag, mu.v.Pt(), w.global_weight());
        h.fill("muEta", "third_" + tag, mu.v.Eta(), w.global_weight());
        h.fill("muPhi", "third_" + tag, mu.v.Phi(), w.global_weight());
        h.fill("muIso", "third_" + tag, mu.iso, w.global_weight());
    }
    if (muons.size() > 3) {
        const lepton &mu = muons[3];
        h.fill("muPt", "fourth_" + tag, mu.v.Pt(), w.global_weight());
        h.fill("muEta", "fourth_" + tag, mu.v.Eta(), w.global_weight());
        h.fill("muPhi", "fourth_" + tag, mu.v.Phi(), w.global_weight());
        h.fill("muIso", "fourth_" + tag, mu.iso, w.global_weight());
    }
}
} // namespace physics
