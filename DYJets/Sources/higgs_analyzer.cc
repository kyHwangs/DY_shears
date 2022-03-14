#include "higgs_analyzer.h"

#include <algorithm>
#include <array>
#include <utility>

#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "functions.h"

higgs_analyzer::higgs_analyzer(util::job::info &info, const util::options &opt)
    : boson_jets_analyzer(info, opt),
      _zfinder_good(opt, "good Z"),
      _zfinder_bad(opt, "bad Z")
{
    counter.declare("With four good leptons");
    counter.declare("With four good electrons");
    counter.declare("With four good muons");
    counter.declare("ZZ->4e");
    counter.declare("ZZ->4mu");
    counter.declare("ZZ->2e2mu");

    std::string binning_file = "higgs-binnings.yml";
    if (opt.config["binning file"]) {
        binning_file = opt.config["binning file"].as<std::string>();
    }
    util::logging::info << "Taking binnings from file " << binning_file << std::endl;
    histo_set.set_style(util::style_list(YAML::LoadFile(binning_file)));
    histo_set2D.set_style(util::style_list(YAML::LoadFile(binning_file)));
}

void higgs_analyzer::fill(const util::matched<std::string> &tags,
                          const util::matched<event_contents> &evt)
{
    boson_jets_analyzer::fill(tags, evt);

    if (!tags.rec || !evt.rec) {
        return;
    }

    const auto boson = evt.rec->leptons;
    const auto jets = evt.rec->jets;

    physics::dilepton good_Z(boson[0], boson[1]);
    physics::dilepton bad_Z(boson[2], boson[3]);

    TLorentzVector higgs = good_Z.v + bad_Z.v;

    histo_set.fill("mass", *tags.rec, higgs.M(), weights().global_weight());
//     FIXME
//     histo_set.fill("phistar", *tags.rec, Z.phistar(), weights().global_weight());
    histo_set.fill("pt", *tags.rec, higgs.Pt(), weights().global_weight());
}

void higgs_analyzer::apply_trigger_sf(
    physics::weights &weights, const std::vector<physics::lepton> &leptons)
{
    // FIXME
}

std::vector<physics::lepton>
higgs_analyzer::find_boson(const std::vector<physics::lepton> &muons,
                           const std::vector<physics::lepton> &electrons)
{
    using physics::dilepton;
    using physics::lepton;

    if (muons.size() + electrons.size() < 4) {
        return {};
    }
    counter.count("With four good leptons", weights().global_weight());

    if (electrons.size() >= 4) {
        counter.count("With four good electrons", weights().global_weight());
    }
    if (muons.size() >= 4) {
        counter.count("With four good muons", weights().global_weight());
    }

    std::vector<lepton> leptons = muons;
    leptons.insert(leptons.end(), electrons.begin(), electrons.end());

    bool found = false;
    std::pair<dilepton, dilepton> chosen_pair;

    std::vector<dilepton> good = _zfinder_good.find(leptons);
    for (const dilepton &z1 : good) {
        // Build a lepton list without the ones used to build the good Z candidate
        std::vector<lepton> vetoed_leptons = leptons;
        vetoed_leptons.erase(
            std::find(vetoed_leptons.begin(), vetoed_leptons.end(), z1.a));
        vetoed_leptons.erase(
            std::find(vetoed_leptons.begin(), vetoed_leptons.end(), z1.b));

        // Loop on all bad Z candidates
        std::vector<dilepton> bad = _zfinder_bad.find(vetoed_leptons);
        for (const dilepton &z2 : bad) {
            if (!dilepton::zmass_ordering(z1, z2)) {
                continue;
            }
            // List of leptons, sorted by ascending pt
            std::array<lepton, 4> used{ z1.a, z1.b, z2.a, z2.b };
            std::sort(used.begin(),
                      used.end(),
                      [](const lepton &a, const lepton &b) { return a.v.Pt() < b.v.Pt(); });

            // Lepton pt cuts
            if (used[3].v.Pt() < 20 || used[2].v.Pt() < 10) {
                continue;
            }

            // All M_ll > 4 GeV
            bool pass = true;
            for (auto it1 = used.begin(); it1 != used.end(); ++it1) {
                for (auto it2 = it1 + 1; it2 != used.end(); ++it2) {
                    const TLorentzVector &p1 = it1->v;
                    const TLorentzVector &p2 = it2->v;
                    if ((p1 + p2).M() < 4 || deltaR(p1, p2) < 0.02) {
                        pass = false;
                        break;
                    }
                }
                if (!pass) {
                    break;
                }
            }
            if (!pass) {
                continue;
            }

            // Other Z candidate
            /// @todo If 4e/4mu
            /*if (TODO)*/ {
                dilepton other1(z1.a, z1.a.charge == z2.a.charge ? z2.b : z2.a);
                dilepton other2(z1.b, z1.b.charge == z2.b.charge ? z2.a : z2.b);
                if (other2.v.M() < 12 && other1.distance_to_z() < z1.distance_to_z()) {
                    continue;
                }
            }

            chosen_pair = {z1, z2};
            found = true;

            break;
        }
        if (found) {
            break;
        }
    }

    if (!found) {
        return {};
    }

    if (std::abs(chosen_pair.first.a.pdgid) == 11 &&
        std::abs(chosen_pair.second.a.pdgid) == 11) {
        counter.count("ZZ->4e", weights().global_weight());
    } else if (std::abs(chosen_pair.first.a.pdgid) == 13 &&
               std::abs(chosen_pair.second.a.pdgid) == 13) {
        counter.count("ZZ->4mu", weights().global_weight());
    } else {
        counter.count("ZZ->2e2mu", weights().global_weight());
    }

    return {
        chosen_pair.first.a,
        chosen_pair.first.b,
        chosen_pair.second.a,
        chosen_pair.second.b};
}

std::vector<physics::lepton>
higgs_analyzer::find_gen_boson(const std::vector<physics::lepton> &genleps)
{
    // TODO
    return {};
}

po::options_description higgs_analyzer::options()
{
    return po::options_description("Physics options");
}
