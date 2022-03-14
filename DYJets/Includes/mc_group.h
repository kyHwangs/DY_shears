#ifndef MC_GROUP_H
#define MC_GROUP_H

#include <memory>
#include <set>
#include <string>

#include "sample.h"

class TFile;
class TH1;

namespace data
{

class comparison_entry;

/// \brief Describes a group of MC samples that will be added together.
class mc_group
{
  private:
    friend struct YAML::convert<data::mc_group>;

    struct sample_data
    {
        data::sample sample;
        std::shared_ptr<comparison_entry> centry = nullptr;
    };

    bool _required = false;
    double _scale_factor = 1;
    int _color;
    bool _is_signal = false;
    std::string _legend;
    std::vector<sample_data> _sample_data;
    std::vector<std::string> _sample_names;

  public:
    virtual ~mc_group() = default;

    /// \brief Returns whether the group is required.
    bool required() const { return _required; }

    /// \brief Returns the group color (using \c ROOT conventions).
    int color() const { return _color; }

    /// \brief Returns whether the group is tagged as 'signal'
    bool is_signal() const { return _is_signal; }

    /// \brief Returns the group legend.
    std::string legend() const { return _legend; }

    /// \brief Returns the list of samples in this group
    std::vector<sample> samples() const;

    /// \brief Adds the name of all available histograms to \c histos
    void add_histograms(std::set<std::string> &histos);

    /**
     * \brief Retrieves the histogram named \c name, or \c nullptr if it's not available
     *
     * The histogram will be ready for drawing.
     */
    std::unique_ptr<TH1> get(const std::string &name);

    /**
     * \brief Loads the list of groups from the configuration file.
     *
     * \param opt Configuration
     * \param analyzer_name The name of the analyzer
     * \param input_dir The directory where files are stored
     * \param all_samples The list of all available samples
     * \param open_files Whether to open files (not everything will work if set to \c false)
     */
    static std::vector<mc_group> load(const util::options &opt,
                                      const std::string &analyzer_name,
                                      const std::string &input_dir,
                                      const std::vector<data::sample> &all_samples,
                                      bool open_files = true);

  private:
    void init(const std::vector<sample> &all_samples,
              const std::string &analyzer_name,
              const std::string &input_dir,
              bool open_files);
};
}

#endif // MC_GROUP_H
