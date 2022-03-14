#ifndef SAMPLE_H
#define SAMPLE_H

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "catalog.h"
#include "options.h"

class TFile;

namespace data
{

/**
 * \brief Represents a data sample.
 *
 * Objects of this class can be loaded from the configuration file using \ref load.
 */
class sample
{
  public:
    /// \brief Identifies data, MC and background samples
    enum type {
        data,      ///< \brief Data sample
        mc,        ///< \brief MC sample
        background ///< \brief Background sample
    };

  private:
    friend struct YAML::convert<data::sample>;

    type _type;
    std::string _nano_dir;
    std::string _catalog;
    std::string _name;
    unsigned _jobs = 1;
    double _xsec = -1;

    std::pair<bool, bool> _has_triggers;
    std::pair<std::string, std::string> _triggers;

    std::vector<std::string> _flags;

  public:
    /// \brief Retrives the catalog for this sample.
    data::catalog catalog() const;

    /**
     * \brief Finds and opens the file containing histograms for this sample.
     *
     * If \c directory doesn't exist, it is created first (except when \c mode is \c "READ").
     *
     * If the file couldn't be opened, returns a null pointer.
     *
     * \param analyzer_name A unique prefix for the analyzer
     * \param directory     The directory in which files are stored
     * \param mode          Passed to \c TFile::Open.
     * \param job_id        The job id If the sample is divided in jobs, else \c -1.
     */
    std::shared_ptr<TFile> histogram_file(const std::string &analyzer_name,
                                          const std::string &directory,
                                          const std::string &mode = "READ",
                                          int job_id = -1) const;

    /// \brief Retrieves the number of jobs to use for the sample.
    unsigned jobs() const { return _jobs; }

    /// \brief Retrieves the name of the sample.
    std::string name() const { return _name; }

    /// \brief Retrieves the name of the sample.
    double xsec() const { return _xsec; }

    /// \brief Returns the list of triggers to be used in this sample
    std::pair<std::string, std::string> triggers() const { return _triggers; }

    /// \brief Returns the list of triggers to be used in this sample
    std::pair<bool, bool> has_triggers() const { return _has_triggers; }

    /// \brief Retrieves the \ref type of the sample.
    enum type type() const { return _type; }

    std::vector<std::string> flags() const { return _flags; }

    /// \brief Loads the list of samples from the configuration file.
    static std::vector<sample> load(const util::options &opt);
};
} // namespace data

#endif // SAMPLE_H
