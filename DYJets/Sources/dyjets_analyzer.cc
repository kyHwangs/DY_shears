#include "dyjets_analyzer.h"

#include <algorithm>
#include <array>
#include <random>

#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "cmake_config.h" // DEBUG_PRINTOUT
#include "functions.h"
#include "lepton.h"

dyjets_analyzer::dyjets_analyzer(util::job::info &info, const util::options &opt)
    : boson_jets_analyzer(info, opt), _zfinder(opt, "Z"), event(info.reader, "event"),run(info.reader, "run")
{
    counter.declare("With two good leptons");
    counter.declare("With two good electrons");
    counter.declare("With two good muons");
    counter.declare("With a good Z boson");
    counter.declare("With two good gen leptons");
    counter.declare("With two good gen electrons");
    counter.declare("With two good gen muons");
    counter.declare("With a good gen Z boson");

    std::string binning_file = "dyjets-binnings.yml";
    if (opt.config["binning file"]) {
        binning_file = opt.config["binning file"].as<std::string>();
    }
    util::logging::info << "Taking binnings from file " << binning_file << std::endl;
    histo_set.set_style(util::style_list(YAML::LoadFile(binning_file)));
    histo_set2D.set_style(util::style_list(YAML::LoadFile(binning_file)));
}

namespace /* anonymous */
{

void apply_dimu_trigger_sf(physics::weights &w,
                           const physics::lepton &mu1,
                           const physics::lepton &mu2,
                           const util::tables &tab)
{
    if (w.ismc()) {
        if (mu1.raw_v.Pt() > mu2.raw_v.Pt()) {
            w.use_weight(tab.at("dimu trigger")
                             .getEfficiency(std::abs(mu1.raw_v.Eta()), std::abs(mu2.raw_v.Eta())));
        } else {
            w.use_weight(tab.at("dimu trigger")
                             .getEfficiency(std::abs(mu2.raw_v.Eta()), std::abs(mu1.raw_v.Eta())));
        }
    }
}

void apply_diel_trigger_sf(physics::weights &w,
                           const physics::lepton &e1,
                           const physics::lepton &e2,
                           const util::tables &tab)
{
    if (w.ismc()) {
        w.use_weight(tab.at("diel trigger leg1")
                        .getEfficiency(e1.v.Pt(), e1.raw_v.Eta()));
        w.use_weight(tab.at("diel trigger leg2")
                        .getEfficiency(e2.v.Pt(), e2.raw_v.Eta()));
    }
}

void apply_emu_trigger_sf(physics::weights &w,
                          const physics::lepton &l1,
                          const physics::lepton &l2,
                          const util::tables &tab)
{
    if (w.ismc()) {
        auto mu = l1.pdgid == 13 ? l1 : l2;
        w.use_weight(tab.at("emu trigger")
                        .getEfficiency(l1.v.Pt(), std::abs(l1.raw_v.Eta())));
    }
}
} // namespace anonymous

void dyjets_analyzer::apply_trigger_sf(physics::weights &weights,
                                       const std::vector<physics::lepton> &leptons)
{
    using physics::zfinder;

    switch (_zfinder.get_flavor_mode()) {
    case zfinder::flavor_mode::mumu:
        apply_dimu_trigger_sf(weights, leptons[0], leptons[1], tables());
        break;
    case zfinder::flavor_mode::ee:
        apply_diel_trigger_sf(weights, leptons[0], leptons[1], tables());
        break;
    case zfinder::flavor_mode::emu:
        apply_emu_trigger_sf(weights, leptons[0], leptons[1], tables());
        break;
    case zfinder::flavor_mode::same:
    case zfinder::flavor_mode::none:
        if (leptons[0].pdgid == 13) {
            apply_dimu_trigger_sf(weights, leptons[0], leptons[1], tables());
        } else {
            apply_diel_trigger_sf(weights, leptons[0], leptons[1], tables());
        }
        break;
    }
}

void dyjets_analyzer::fill(const util::matched<std::string> &tags,
                           const util::matched<event_contents> &evt)
{
    boson_jets_analyzer::fill(tags, evt);

    auto mass = evt.apply(&event_contents::get_boson_p).apply(&TLorentzVector::M);

#ifdef DEBUG_PRINTOUT
    if (mass.rec && tags.rec && *tags.rec == "inc0jet_mass76_106") {
        std::cout << "----\n";
        for (const auto &mu : (*evt.rec).leptons) {
            double musf = tables().at("muon id").getEfficiency(mu.v.Pt(), mu.v.Eta());
            double elsf = tables().at("electron reco").getEfficiency(mu.v.Pt(), mu.raw_v.Eta());
            std::cout << "Mu eta, pt, sf  "
                      << setprecision(5)
                      << mu.v.Eta()
                      << ", "
                      << mu.v.Pt()
                      << ", "
                      << (weights().ismc() ? (std::abs(mu.pdgid) == 11 ? elsf : musf) : 1.0)
                      << "\n";
        }
        std::cout << "event "
                  << *event
                  << ", weight  = "
                  << weights().global_weight()
                  << "\n";
        if (std::abs(evt.rec->leptons[0].pdgid) == 13) {
            // Sort wrt raw_v for Rochester input
            auto lep_copy = (*evt.rec).leptons;
            std::sort(lep_copy.begin(),
                    lep_copy.end(),
                    [](const physics::lepton &a, const physics::lepton &b) {
                        return a.raw_v.Pt() > b.raw_v.Pt();
                    });
            for (const auto &l : lep_copy) {
                std::cout << "Rochester input eta, pt, phi, charge, function name : "
                    << setprecision(5)
                    << l.raw_v.Eta()
                    << ", "
                    << l.raw_v.Pt()
                    << ", "
                    << l.raw_v.Phi()
                    << ", "
                    << l.charge;
                switch (l.fnUsed) {
                    case 0:
                        cout << ", kScaleFromGenMC\n";
                        break;
                    case 1:
                        cout << ", kScaleAndSmearMC\n";
                        break;
                    case 2:
                        cout << ", kScaleDT\n";
                        break;
                    default:
                        cout << ", UNKNOWN\n";
                        break;
                }
            }
        }
        auto jets = evt.apply(&event_contents::get_jets);
        if (jets.rec) {
            std::cout << jets.rec->size() << std::endl;
            for (const auto &j : *jets.rec) {
                std::cout <<setprecision(5)<< "Jet " << j.v.Eta() << ", " << j.v.Pt() << std::endl;
            }
        }
        if (mass.gen && tags.gen
                     && *tags.gen == "inc0jet_mass76_106"
                     && evt.gen->leptons.size() >= 2
                     && std::abs(evt.gen->leptons[0].v.Eta()) < 2.4
                     && std::abs(evt.gen->leptons[1].v.Eta()) < 2.4
                     && std::abs(evt.gen->leptons[0].v.Pt()) > 25
                     && std::abs(evt.gen->leptons[1].v.Pt()) > 20) {
            for (const auto &l : (*evt.gen).leptons) {
                std::cout << "Gen mu eta, pt, phi "
                          << setprecision(5)
                          << l.v.Eta()
                          << ", "
                          << l.v.Pt()
                          << ", "
                          << l.v.Phi()
                          << "\n";
            }
            std::cout << "event " << *event << "\n";
        }
    }
#endif // DEBUG_PRINTOUT

    fill_unfolded("mass", tags, mass);
    fill_unfolded("mass_wide_range", tags, mass);

    auto pt = evt.apply(&event_contents::get_boson_p)
                 .apply((double (TLorentzVector::*)() const) &TLorentzVector::Pt);
    fill_unfolded("pt", tags, pt);

    if (!tags.rec || !evt.rec) {
        return;
    }

    const auto boson = evt.rec->leptons;
    const auto jets = evt.rec->jets;

    physics::dilepton Z(boson[0], boson[1]);

    histo_set.fill("phistar", *tags.rec, Z.phistar(), weights().global_weight());

    /*
     * Variables in Z rest frame
     */

    TLorentzVector pZ = Z.v;
    TLorentzVector pLep1 = Z.a.v;
    TLorentzVector pLep2 = Z.b.v;
    if (Z.a.charge > 0) {
        std::swap(pLep1, pLep2);
    }

    // Rotate to decay frame with Z axis parallel to boson momentum
    double phi = Z.v.Phi();
    double theta = Z.v.Theta();

    pZ.RotateZ(-phi);
    pZ.RotateY(-theta);
    pLep1.RotateZ(-phi);
    pLep1.RotateY(-theta);
    pLep2.RotateZ(-phi);
    pLep2.RotateY(-theta);

    // Boost to decay frame with boson at rest
    pLep1.Boost(-pZ.BoostVector());
    pLep2.Boost(-pZ.BoostVector());
    pZ.Boost(-pZ.BoostVector());

    histo_set.fill("decay_costheta", *tags.rec,
                   std::cos(pLep1.Theta()), weights().global_weight());
}

std::vector<physics::lepton>
dyjets_analyzer::find_boson(const std::vector<physics::lepton> &muons,
                            const std::vector<physics::lepton> &electrons)
{
    using physics::dilepton;
    using physics::lepton;
    using physics::zfinder;

    if (_zfinder.get_flavor_mode() == zfinder::flavor_mode::emu) {
        // Special treatment for same flavor mode
        if (muons.size() + electrons.size() < 2) {
            return {};
        }
    }

    if (muons.size() >= 2) {
        counter.count("With two good muons", weights().global_weight());
        
    }
    //cout << "mu size" <<muons.size() << endl;
    if (electrons.size() >= 2) {
        counter.count("With two good electrons", weights().global_weight());
    }

    std::vector<lepton> leptons;
    switch (_zfinder.get_flavor_mode()) {
    case zfinder::flavor_mode::mumu:
        leptons = muons;
        break;
    case zfinder::flavor_mode::ee:
        leptons = electrons;
        break;
    case zfinder::flavor_mode::emu:
    case zfinder::flavor_mode::same:
    case zfinder::flavor_mode::none:
        leptons.insert(leptons.end(), muons.begin(), muons.end());
        leptons.insert(leptons.end(), electrons.begin(), electrons.end());
        // Sort the merged collection
        std::sort(leptons.begin(), leptons.end(), [](const lepton &lhs, const lepton &rhs)
            {
                return lhs.v.Pt() > rhs.v.Pt();
            }
        );
        break;
    }

    if (leptons.size() < 2) {
        return {};
    }
    
//    std::cout << "Event no, Mu pt, Mu eta "
//                      << setprecision(5)
//                      << *event
//                      << ", "
//                      << leptons[0].v.Pt()
//                      << ", "
//                      << leptons[1].v.Pt()
//                      << ", "
//                      << leptons[0].v.Eta()
//                      << ", "
//                      << leptons[1].v.Eta()
//                      << "\n";

    if (_zfinder.get_flavor_mode() != zfinder::flavor_mode::emu && leptons[0].v.Pt() < 25) {
        return {};
    } else if (_zfinder.get_flavor_mode() == zfinder::flavor_mode::emu) {
        // Since we use an SMu trigger, always cut on the muon
        auto mu = leptons[0].pdgid == 13 ? leptons[0] : leptons[1];
        if (mu.v.Pt() < 25) {
            return {};
        }
    }

    counter.count("With two good leptons", weights().global_weight());

    std::vector<dilepton> candidates = _zfinder.find({leptons[0], leptons[1]});
    if (candidates.size() == 0) {
        return {};
    }
    counter.count("With a good Z boson", weights().global_weight());

    std::sort(candidates.begin(), candidates.end(), dilepton::zmass_ordering);
    dilepton Z = candidates[0];
    return {Z.a, Z.b};
}

std::vector<physics::lepton>
dyjets_analyzer::find_gen_boson(const std::vector<physics::lepton> &genleps)
{
    if (genleps.size() < 2) {
        return {};
    }
    counter.count("With two good gen leptons", weights().global_weight());

    if (std::abs(genleps[0].pdgid) == 13) {
        counter.count("With two good gen muons", weights().global_weight());
    } else {
        counter.count("With two good gen electrons", weights().global_weight());
    }

    if (genleps[0].v.Pt() < 25) {
        return {};
    }

    std::vector<physics::dilepton> candidates = _zfinder.find({genleps[0], genleps[1]});
    if (candidates.size() == 0) {
        return {};
    }
    counter.count("With a good gen Z boson", weights().global_weight());

    std::sort(candidates.begin(), candidates.end(), physics::dilepton::zmass_ordering);
    physics::dilepton Z = candidates[0];
    return {Z.a, Z.b};
}

po::options_description dyjets_analyzer::options()
{
    return po::options_description("Physics options");
}
