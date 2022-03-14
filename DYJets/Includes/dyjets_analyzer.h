#ifndef DYJETS_ANALYZER_H
#define DYJETS_ANALYZER_H

#include <boost/program_options/options_description.hpp>

#include <TH1D.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "boson_jets_analyzer.h"
#include "jets.h"
#include "job.h"
#include "pileup.h"
#include "zfinder.h"

namespace po = boost::program_options;

/// \brief Implements a \f$ Z \to 2l \f$ analysis.
class dyjets_analyzer : public physics::boson_jets_analyzer
{
    physics::zfinder _zfinder;
    TTreeReaderValue<unsigned long long> event;
    TTreeReaderValue<unsigned> run;

  public:
    /// \brief Constructor.
    explicit dyjets_analyzer(util::job::info &info, const util::options &opt);

    /// \brief Destructor.
    virtual ~dyjets_analyzer() = default;

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

    std::vector<physics::lepton> find_gen_boson(
    const std::vector<physics::lepton> &genleps) override;
    /// \brief Returns the list of options supported by the analyzer.
    static po::options_description options();
};

#endif // DYJETS_ANALYZER_H
