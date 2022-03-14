#ifndef ELECTRONS_H
#define ELECTRONS_H

#include <memory>
#include <vector>

#include <TTreeReaderArray.h>

#include "histo_set.h"
#include "job.h"
#include "lepton.h"
#include "options.h"
#include "tables.h"
#include "weights.h"

namespace physics
{

/// \brief Handles electrons.
class electrons
{
  public:
    /// \brief Represents electron IDs
    enum class id
    {
        veto,   ///< \brief veto ID
        loose,  ///< \brief Loose ID
        medium, ///< \brief Medium ID
        tight,  ///< \brief Tight ID
    };

  private:
    TTreeReaderArray<float> Electron_pt;
    TTreeReaderArray<float> Electron_eta;
    TTreeReaderArray<float> Electron_phi;
    TTreeReaderArray<float> Electron_mass;
    TTreeReaderArray<int> Electron_charge;
    TTreeReaderArray<float> Electron_deltaEtaSC;
    TTreeReaderArray<float> Electron_miniPFRelIso_all;
    TTreeReaderArray<int> Electron_cutBased;

    double _pt_cut = 20;
    double _eta_cut = 2.4;
    double _iso_cut = 0.25;
    id _id_cut = id::tight;

    bool _id_sf_enabled = true;
    bool _reco_sf_enabled = true;

  public:
    /// \brief Constructor.
    explicit electrons(util::job::info &info, const util::options &opt, util::histo_set &h);

    /// \brief Configures the analyzer from user input.
    void configure(const util::options &opt);

    /**
     * \brief Retrieves a list of all electrons in the current event.
     *
     * The list is already filtered according to config file options.
     */
    std::vector<lepton> get( int & nVetoElecs);

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
     * \param electrons List of electrons to take into account
     * \param tab Tables to get the scale factors from
     *
     * \throws std::out_of_range if a table is enabled and not present.
     */
    void apply_sf(weights &w, const std::vector<lepton> &electrons, const util::tables &tab) const;

    /**
     * \brief Fills muon control plots.
     * \param electrons The list of electrons in the event.
     * \param tag       A tag to pass to \ref util::histo_set
     */
    void fill(util::histo_set &h,
              const std::string &tag,
              const std::vector<lepton> &electrons,
              const weights &w);

    /// \brief Writes histograms to the current directory.
    void write();
};
} // namespace physics

#endif // ELECTRONS_H
