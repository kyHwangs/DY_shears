#ifndef HIGGS_ANALYZER_H
#define HIGGS_ANALYZER_H

#include <boost/program_options/options_description.hpp>

#include <TH1D.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "boson_jets_analyzer.h"
#include "jets.h"
#include "job.h"
#include "zfinder.h"

namespace po = boost::program_options;

/// \brief Implements a \f$ H \to 4l \f$ analysis.
class higgs_analyzer : public physics::boson_jets_analyzer
{
    physics::zfinder _zfinder_good, _zfinder_bad;

  public:
    /// \brief Constructor.
    explicit higgs_analyzer(util::job::info &info, const util::options &opt);

    /// \brief Destructor.
    virtual ~higgs_analyzer() = default;

    // Overridden from base class
    void apply_trigger_sf(physics::weights &weights,
                          const std::vector<physics::lepton> &leptons) override;

    /// \brief Fills histograms
    void fill(const util::matched<std::string> &tags,
              const util::matched<event_contents> &evt) override;

    // Overridden from base class
    std::vector<physics::lepton> find_boson(
        const std::vector<physics::lepton> &muons,
        const std::vector<physics::lepton> &electrons) override;

    // Overridden from base class
    std::vector<physics::lepton> find_gen_boson(
        const std::vector<physics::lepton> &genleps) override;

    /// \brief Returns the list of options supported by the analyzer.
    static po::options_description options();
};

#endif // HIGGS_ANALYZER_H
