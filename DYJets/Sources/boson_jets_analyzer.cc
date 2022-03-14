#include "boson_jets_analyzer.h"

#include <algorithm>

namespace physics
{

namespace /* anonymous */
{

/**
 * \brief Returns a tag of the form @c low_high where low and high are the bounds of the
 *        bin that contains @c value.
 *
 * Returns @c "" if the value doesn't fit in a bin.
 */
std::string make_tag(double value, const std::vector<double> &bins)
{
    auto high = std::lower_bound(bins.begin(), bins.end(), value);
    if (high == bins.begin() || high == bins.end()) {
        // Out of bounds
        return "";
    }
    auto low = std::prev(high);
    return std::to_string(int(*low)) + "_" + std::to_string(int(*high));
}

/**
 * \brief Retrives the trigger list for the current @ref sample.
 */
std::string get_triggers(util::job::info &info, const util::options &opt, int era)
{
    if (era == 0 && info.sample.has_triggers().first) {
        return info.sample.triggers().first;
    } else if (era == 1 && info.sample.has_triggers().second) {
        return info.sample.triggers().second;
    }
    std::string name = (era == 0) ? "triggers B-F" : "triggers G-H";
    if (opt.config[name]) {
        return opt.config[name].as<std::string>();
    }
}

} // anonymous namespace

boson_jets_analyzer::boson_jets_analyzer(util::job::info &info,
                                         const util::options &opt) :
    run(info.reader, "run"),
    event(info.reader, "event"),
    L1PreFiringWeight_Nom(info.reader, "L1PreFiringWeight_Nom"), //Iti: check
    _rng(0 /*std::random_device()()*/),
    _genleps(info, opt, histo_set),
    _mask_eraBG(info, get_triggers(info, opt, 0)),
    _mask_eraH(info, get_triggers(info, opt, 1)),
    _muons(info, opt, histo_set),
    _electrons(info, opt, histo_set),
    _jets(info, opt),
    _pileup(info, opt),
    _btagger(opt,histo_set2D),
    _reweighing(info, opt),
    _weights(info)
{
    if (opt.config["tables B-F"]) {
        _tables_eraBF = opt.config["tables B-F"].as<util::tables>();
    }
    if (opt.config["tables G-H"]) {
        _tables_eraGH = opt.config["tables G-H"].as<util::tables>();
    }

    if (!opt.config["b jet veto"]) {
        throw std::runtime_error("Missing mandatory section in config file: \"b jet veto\"");
    }

    util::set_value_safe(opt.config["b jet veto"], _bjet_veto, "use", "use b jet veto");

   if (!opt.config["prefiring weights"]) {
        throw std::runtime_error("Missing mandatory section in config file: \"prefiring weights\"");
    }
       util::set_value_safe(opt.config["prefiring weights"], _pref, "use", "use apply prefiring weights");

 if (opt.config["mass bins"]) {
        _mass_bins = opt.config["mass bins"].as<std::vector<double>>();
        std::sort(_mass_bins.begin(), _mass_bins.end());
    }

    _jets.declare_histograms(histo_set);
    _pileup.declare_histograms(histo_set);

    counter.declare("Total");
    counter.declare("Passing the trigger");
}

boson_jets_analyzer::~boson_jets_analyzer()
{}

void boson_jets_analyzer::operator()()
{
    _weights.process_event();
    _reweighing.reweigh(_weights);
    counter.count("Total", weights().global_weight());

    /*
     * Choose the right era for this event
     */
    const unsigned run_threshold = 278820u; // start of Run G
    const double run_lumi_fraction = 0.5493217216546642; // lumi fraction before run G

    std::uniform_real_distribution<> uniform(0.0, 1.0);

    if (weights().isdata()) {
        // Run number based era selection
        if (*run < run_threshold) {
            _era = 0;
        } else {
            _era = 1;
        }
    } else {
        // Monte-Carlo based era selection
        if (uniform(rng()) < run_lumi_fraction) {
            _era = 0;
        } else {
            _era = 1;
        }
    }

    // Event
    util::matched<event_contents> evt;
    std::vector<lepton> genleps = {};

    // fill the generator level histograms only when the sample is MC
    if( !weights().isdata() ) {
        genleps = _genleps.get();
        std::vector<lepton> genleptons = find_gen_boson(genleps);

        if (!genleptons.empty()) {
            // Gen boson found
            evt.gen = event_contents();
            evt.gen->leptons = genleptons;
            evt.gen->boson_p = std::accumulate(
                genleptons.begin(),
                genleptons.end(),
                TLorentzVector(),
                [](const TLorentzVector &p, const lepton &lep) { return p + lep.v; });

            _genleps.fill(histo_set, "geninc0jet_noweight", genleptons, weights());
        }
    }

    /*
     * Handle the trigger
     */
    bool triggered = passes_trigger();
    if (!evt.gen && !triggered) {
        // End early if nothing to do
        return;
    } else if (triggered) {
        counter.count("Passing the trigger", weights().global_weight());
    }

    /*
     * Pre-firing weight
     */
    //if (EvtPrefiringweight.GetSize() > 0 && _pref) {
      //  _weights.use_weight(EvtPrefiringweight[0]);
    //}

    /*
     * Read leptons and find the boson
     */
    if (triggered) {
        int nVetoMuons=0;
        int nVetoElecs=0;

        std::vector<lepton> muons = _muons.get(weights().isdata(), rng(), genleps,nVetoMuons);
        std::vector<lepton> electrons = _electrons.get(nVetoElecs);
        std::vector<lepton> leptons = find_boson(muons, electrons);
        if (!leptons.empty()&&(nVetoMuons+nVetoElecs)<=2) {
            // Rec boson found
            evt.rec = event_contents();
            evt.rec->leptons = leptons;
            evt.rec->boson_p = std::accumulate(
                leptons.begin(),
                leptons.end(),
                TLorentzVector(),
                [](const TLorentzVector &p, const lepton &lep) { return p + lep.v; });
        }
    }

    if (!evt.gen && !evt.rec) {
        // End early if nothing to do
        return;
    }

    /*
     * At this point at least one boson was found, either gen or rec.
     * Apply rec scale factors and load rec jets, then apply b tagging.
     */

    // Create lists of chosen muons and electrons to use for scale factors
    if (evt.rec) {
        std::vector<lepton> chosen_muons, chosen_electrons;
        std::copy_if(evt.rec->leptons.begin(),
                     evt.rec->leptons.end(),
                     std::back_inserter(chosen_muons),
                     [](const lepton &lep) { return lep.pdgid == 13; });
        std::copy_if(evt.rec->leptons.begin(),
                     evt.rec->leptons.end(),
                     std::back_inserter(chosen_electrons),
                     [](const lepton &lep) { return lep.pdgid == 11; });

        // Apply lepton scale factors
        _muons.apply_sf(_weights, chosen_muons, tables());
        _electrons.apply_sf(_weights, chosen_electrons, tables());

        // Fill lepton control plots
        _muons.fill(histo_set, "inc0jet_noweight", chosen_muons, weights());
        _electrons.fill(histo_set, "inc0jet_noweight", chosen_electrons, weights());
    }

    /*
     * Handle jets and pileup
     */
    if (evt.rec) {
       std::vector<lepton> l ;
        if (evt.gen) l = evt.gen->leptons;
        evt.rec->jets = _jets.get(weights().isdata(),l);
        evt.rec->jets20 = _jets.get(weights().isdata(), l, 20); // For b veto
        _jets.veto(evt.rec->jets, evt.rec->leptons);
        _jets.veto(evt.rec->jets20, evt.rec->leptons);

        // Calculate b efficiencies and apply scale factors
        if (_bjet_veto && _btagger.any(evt.rec->jets20, _weights, histo_set2D, tables())) {
            evt.rec = boost::none;
            if (!evt.gen && !evt.rec) {
                // End early if vetoed and no gen boson
                return;
            }
        }

        if (evt.rec) { // May have been zero'ed by b veto
            _jets.fill(histo_set, "inc0jet_noweight", evt.rec->jets, weights());
            _pileup.fill(histo_set, "inc0jet_noweight", weights());

            _pileup.reweight(_weights);
        }
    }
    if (evt.gen) {
        evt.gen->jets = _jets.getGen();
        _jets.veto(evt.gen->jets, evt.gen->leptons);
    }

    /*
     * Apply lepton trigger scale factors
     */
    if (evt.rec) {
        apply_trigger_sf(_weights, evt.rec->leptons);
    }

    /*
     * Fill histograms w.r.t. N_jets and invariant mass
     */
    auto mass_tags = evt.apply(&event_contents::get_boson_p)
                        .apply(&TLorentzVector::M)
                        .apply(make_tag, _mass_bins);
    if (mass_tags.rec && mass_tags.rec->empty()) {
        mass_tags.rec = boost::none;
    } else if (mass_tags.rec) {
        mass_tags.rec = "_mass" + *mass_tags.rec;
    }
    if (mass_tags.gen && mass_tags.gen->empty()) {
        mass_tags.gen = boost::none;
    } else if (mass_tags.gen) {
        mass_tags.gen = "_mass" + *mass_tags.gen;
    }

    auto njets = evt.apply(&event_contents::get_jets)
                    .apply(&std::vector<jet>::size);

    {
        // Exclusive
        util::matched<std::string> tags; // eg "exc1jet"
        util::matched<std::string> tags_mass; // eg "exc1jet_mass50_71"

        if (njets.rec && *njets.rec < 3) {
            tags.rec = "exc" + std::to_string(*njets.rec) + "jet";
            tags_mass.rec = *tags.rec + *mass_tags.rec;
        }
        if (njets.gen && *njets.gen < 3) {
            tags.gen = "exc" + std::to_string(*njets.gen) + "jet";
            tags_mass.gen = *tags.gen + *mass_tags.gen;
        }

        if (tags.gen || tags.rec) {
            fill(tags, evt);
            fill(tags_mass, evt);
        }
    }

    // Inclusive
    for (std::size_t nj = 0; nj < 3; ++nj) {
        // Exclusive
        util::matched<std::string> tags; // eg "inc1jet"
        util::matched<std::string> tags_mass; // eg "inc1jet_mass50_71"

        if (njets.rec && *njets.rec >= nj) {
            tags.rec = "inc" + std::to_string(nj) + "jet";
            tags_mass.rec = *tags.rec + *mass_tags.rec;
        }
        if (njets.gen && *njets.gen >= nj) {
            tags.gen = "inc" + std::to_string(nj) + "jet";
            tags_mass.gen = *tags.gen + *mass_tags.gen;
        }

        if (tags.gen || tags.rec) {
            fill(tags, evt);
            fill(tags_mass, evt);
        }
    }
}

void boson_jets_analyzer::fill(const util::matched<std::string> &tags,
                               const util::matched<event_contents> &evt)
{
    if (tags.rec && evt.rec) {
        _jets.fill(histo_set, *tags.rec, evt.rec->jets, weights());
        _pileup.fill(histo_set, *tags.rec, weights());

        // Create lists of chosen muons and electrons
        std::vector<lepton> chosen_muons, chosen_electrons;
        std::copy_if(evt.rec->leptons.begin(), evt.rec->leptons.end(), std::back_inserter(chosen_muons),
                    [](const lepton &lep) { return lep.pdgid == 13; });
        std::copy_if(evt.rec->leptons.begin(), evt.rec->leptons.end(), std::back_inserter(chosen_electrons),
                    [](const lepton &lep) { return lep.pdgid == 11; });

        // Fill lepton control plots
        _muons.fill(histo_set, *tags.rec, chosen_muons, weights());
        _electrons.fill(histo_set, *tags.rec, chosen_electrons, weights());
    }
    if (tags.gen && evt.gen) {
        _genleps.fill(histo_set, *tags.gen, evt.gen->leptons, weights());
    }
}

void boson_jets_analyzer::fill_unfolded(const std::string &name,
                                        const util::matched<std::string> &tags,
                                        const util::matched<double> &value)
{
    // Fill 1D distributions
    if (tags.rec && value.rec) {
        histo_set.fill(name, *tags.rec, *value.rec, weights().global_weight());
    }
    if (tags.gen && value.gen) {
        histo_set.fill(name, *tags.gen + "-gen", *value.gen, weights().gen_weight());
    }
    // Fill response matrix
    if (tags.rec && tags.gen && value.rec && value.gen && tags.rec == tags.gen) {
        histo_set2D.fill(name,
                         *tags.rec + "-matrix",
                         *value.rec,
                         *value.gen, weights().global_weight());
    }
}

bool boson_jets_analyzer::passes_trigger()
{
    if (_weights.ismc()) {
        return _mask_eraH.passes();
    }
    return era_select(_mask_eraBG, _mask_eraH).passes();
}

void boson_jets_analyzer::write()
{
    counter.print();
    _weights.write(&histo_set);
    histo_set.write();
    histo_set2D.write();

}

} // namespace physics
