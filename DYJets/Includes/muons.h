#ifndef MUONS_H
#define MUONS_H

#include <memory>
#include <random>
#include <vector>

#include <TTreeReaderArray.h>

#include "histo_set.h"
#include "job.h"
#include "lepton.h"
#include "options.h"
#include "tables.h"
#include "weights.h"

class RoccoR;

namespace physics
{

/// \brief Handles muons.
class muons
{
  public:
    /// \brief Represents muon IDs
    enum class id
    {
        loose,  ///< \brief Loose ID
        medium, ///< \brief Medium ID
        tight,  ///< \brief Tight ID
    };

    enum class iso
    {
        veryloose,
        loose,  ///< \brief Loose ID
        medium, ///< \brief Medium ID
        tight,  ///< \brief Tight ID
    };


  private:
    TTreeReaderArray<float> Muon_pt;
    TTreeReaderArray<float> Muon_eta;
    TTreeReaderArray<float> Muon_phi;
    TTreeReaderArray<float> Muon_mass;
    TTreeReaderArray<int> Muon_charge;
    TTreeReaderArray<unsigned char> Muon_pfIsoId;
    TTreeReaderArray<int> Muon_nTrackerLayers;
    TTreeReaderArray<float> Muon_pfRelIso04_all;
    TTreeReaderArray<bool> Muon_tightId;
    TTreeReaderArray<bool> Muon_mediumId;
    TTreeReaderArray<bool> Muon_looseId;

    double _pt_cut = 20;
    double _eta_cut = 2.4;
    iso _iso_cut = iso::tight;
    id _id_cut = id::tight;

    bool _id_sf_enabled = true;
    bool _iso_sf_enabled = true;
    bool _trk_sf_enabled = true;
    bool _roccor_enabled = true;

    std::shared_ptr<RoccoR> _roccor = nullptr;

  public:
    /// \brief Constructor.
    explicit muons(util::job::info &info, const util::options &opt, util::histo_set &h);

    /// \brief Configures the analyzer from user input.
    void configure(const util::options &opt);

    /**
     * \brief Retrieves a list of all muons in the current event.
     *
     * The list is already filtered according to config file options.
     */
    std::vector<lepton> get(bool isdata, std::mt19937 &rng, std::vector<lepton> gl,int &nVetoMuons);

    /**
     * \brief Reweighs an event to take scale factors into account.
     *
     * If enabled in the configuration file, the following tables will be used:
     *
     * - `muon id`
     * - `muon isolation`
     * - `muon tracking`
     *
     * It is the user responsibility to load them. An exception is thrown if they're not present.
     *
     * \param w Weights to be reweighed
     * \param muons List of muons to take into account
     * \param tab Tables to get the scale factors from
     *
     * \throws std::out_of_range if a table is enabled and not present.
     */
    void apply_sf(weights &w, const std::vector<lepton> &muons, const util::tables &tab) const;

    /**
     * \brief Fills muon control plots.
     * \param muons The list of muons in the event.
     * \param tag   A tag to pass to \ref util::histo_set
     */
    void fill(util::histo_set &h,
              const std::string &tag,
              const std::vector<lepton> &muons,
              const weights &w);

    /// \brief Writes histograms to the current directory.
    void write();
};
} // namespace physics

#endif // MUONS_H
