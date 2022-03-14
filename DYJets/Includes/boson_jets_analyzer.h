#ifndef BOSON_JETS_ANALYZER_H
#define BOSON_JETS_ANALYZER_H

#include <random>
#include <string>
#include <vector>

#include "BTagCalibrationStandalone.h"
#include "btagger.h"
#include "electrons.h"
#include "event_counter.h"
#include "genleps.h"
#include "histo_set.h"
#include "histo_set2D.h"
#include "jets.h"
#include "job.h"
#include "lepton.h"
#include "matched.h"
#include "muons.h"
#include "pileup.h"
#include "reweighing.h"
#include "tables.h"
#include "triggers.h"
#include "weights.h"

namespace physics
{

/// \brief Base class for boson-jets analyzers
class boson_jets_analyzer
{
    TTreeReaderValue<unsigned> run;
    TTreeReaderValue<unsigned long long> event;
    TTreeReaderArray<float> L1PreFiringWeight_Nom; //check: Iti
protected:
    util::event_counter counter;
    util::histo_set histo_set;
    util::histo_set2D histo_set2D;

private:
    std::string _short_name, _long_name;

    std::mt19937 _rng;
    int _era;

    util::tables _tables_eraBF;
    util::tables _tables_eraGH;

    trigger_mask _mask_eraBG;
    trigger_mask _mask_eraH;

    genleps _genleps;
    muons _muons;
    electrons _electrons;
    jets _jets;
    pileup _pileup;
    btagger _btagger;

    bool _bjet_veto = true;
    bool _pref =false;
    std::vector<double> _mass_bins;

    physics::reweighing _reweighing;
    physics::weights _weights;

public:
    /// \brief Groups together the contents of a boson-jets event
    class event_contents
    {
      public:
        std::vector<lepton> leptons; ///< \brief Leptons making up the boson
        TLorentzVector boson_p;      ///< \brief Reconstructed boson 4-momentum
        std::vector<jet> jets;       ///< \brief Jets in the event
        std::vector<jet> jets20;     ///< \brief Jets with 20 GeV cut

        /// \brief Returns the list of leptons making up the boson
        const std::vector<lepton> &get_leptons() const { return leptons; }

        /// \brief Returns the reconstructed boson 4-momentum
        TLorentzVector get_boson_p() const { return boson_p; }

        /// \brief Returns the list of jets
        const std::vector<jet> &get_jets() const { return jets; }
    };

    /// \brief Constructor
    boson_jets_analyzer(util::job::info &info, const util::options &opt);

    /// \brief Destructor
    virtual ~boson_jets_analyzer();

    /// \brief Entry point, called for every event
    virtual void operator()() final;

    /// \brief Function called at the end of the processing.
    virtual void write();

protected:
    /// \brief Applies trigger scale factors
    virtual void apply_trigger_sf(class weights &weights,
                                  const std::vector<physics::lepton> &leptons) = 0;

    /**
     * \brief Returns the current era (0 for eraBG, 1 for eraGH).
     *
     * For MC events, the era is chosen randomly.
     */
    int era() const { return _era; }

    /**
     * \brief Returns a reference to the random number generator used by the
     *        analyzer.
     */
    std::mt19937 &rng() { return _rng; }

    /// \brief Selects one of two values based on the current era.
    template <class T> T &era_select(T &eraBG, T &eraGH) const
    {
        return era() == 0 ? eraBG : eraGH;
    }

    /// \brief Fills histograms
    virtual void fill(const util::matched<std::string> &tags,
                      const util::matched<event_contents> &evt);

    /**
     * \brief Reconstructs the boson candidate and returns its constituents.
     *
     * This method must be implemented in derived classes.
     */
    virtual std::vector<lepton> find_boson(const std::vector<lepton> &muons,
                                           const std::vector<lepton> &electrons) = 0;

    virtual std::vector<lepton> find_gen_boson(const std::vector<lepton> &genleps) = 0;

    /// \brief Checks whether the current event passes the trigger.
    virtual bool passes_trigger();

    /// \brief Fills histograms for an unfolded variable
    void fill_unfolded(const std::string &name,
                       const util::matched<std::string> &tags,
                       const util::matched<double> &value);

    util::tables tables() const
    {
        return era_select(_tables_eraBF, _tables_eraGH);
    }

    /// \brief Retrieves the weight information for the current event.
    const physics::weights &weights() const { return _weights; }
};

} // namespace physics

#endif // BOSON_JETS_ANALYZER_H
