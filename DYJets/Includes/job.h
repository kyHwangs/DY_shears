#ifndef JOB_H
#define JOB_H

#include <optional>
#include <string>
#include <vector>

#include <boost/program_options/errors.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/ref.hpp>

#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>

#include "ansi_seq.h"
#include "catalog.h"
#include "chains.h"
#include "logging.h"
#include "options.h"
#include "sample.h"
#include "timer.h"

namespace po = boost::program_options;

namespace util
{

/**
 * \brief The core class to run analysis codes.
 *
 * This class handles all the complexity of looping on events. Its main features
 * are:
 *
 *   * Command line options handling
 *   * Running on a fraction of the data only
 *   * Division of files in a number of jobs (file-based)
 *   * Progress display (through \ref timer)
 *   * Graceful handling of exceptions thrown by the analysis code
 *   * Graceful handling of SIGINT
 */
class job
{
  public:
    struct info
    {
        explicit info(const data::catalog &catalog, const util::chains &chains);

        data::catalog catalog;
        data::sample sample;
        util::chains chains;
        TTreeReader reader;

        /**
         * \brief Initalizes a TTreeReaderValue if the branch is present, 
         *        returns nullopt otherwise
         *
         * This function can be used to initialize a variable of type 
         * @c std::optional<TTreeReaderValue<...>>. It is used as follows in
         * constructors:
         *
         *     genWeight(info.init_optional_branch<decltype(genWeight)>("genWeight")),
         *
         * The same function can be used with @c TTreeReaderArray.
         */
        template<typename T>
        auto init_optional_branch(const char *name)
        {
            // We test whether the branch is present. If it is, we create the 
            // reader (using std::move because TTreeReaderX classes aren't 
            // copyable). If not, we return an empty optional.
            return reader.GetTree()->GetBranch(name) != nullptr
                ? std::move(
                    std::make_optional<typename T::value_type>(reader, name))
                : std::nullopt;
        }
    };

  private:
    int _job_id = 0;
    int _job_count = 1;
    int _max_files = std::numeric_limits<int>::max();
    long long _max_events = std::numeric_limits<long long>::max();
    bool _fatal_exceptions = false;
    bool _graceful_sigint;
    std::atomic<bool> _sigint_caught;
    std::vector<data::sample> _samples;
    std::string _sample_name;
    std::string _analyzer_name;
    std::string _output_dir;

  public:
    /// \brief Constructor.
    explicit job(const std::string &analyzer_name);

    /// \brief Retrieves the sample that will be processed.
    data::sample sample() const;

    /// \brief Retrieves the name of the output file.
    std::string output_dirname() const;

    /**
     * \brief Loops on data.
     *
     * This function creates an instance of the \c Analyzer class, passing a \ref job_info as the
     * first argument, followed by any arguments passed to this function:
     *
     * ~~~
     * new Analyzer(<job_info>, <args...>);
     * ~~~
     *
     * The \c Analyzer class should have the following methods:
     *
     * ~~~{.cpp}
     * void operator()();
     * void write();
     * ~~~
     *
     * The `operator()` method will be called for every event. `write()` will be called once at
     * the end of the processing.
     *
     * \tparam Analyzer An analyzer class.
     * \param args Arguments to pass to the analyzer constructor.
     */
    template <class Analyzer, class... Args> inline void run(Args &... args);

    /**
     * \brief Sets whether the even loop should exit gracefully on \c SIGINT (ie \c Ctrl+C).
     *
     * This setting defaults to \c true if \c stdin, \c stdout or \c stderr is a terminal, and to
     * \c false otherwise.
     *
     * \note This function must be called before \ref run.
     */
    void set_graceful_sigint(bool enable) { _graceful_sigint = enable; }

    /// \brief Configures the job from user input.
    void configure(const options &opt);

    /// \brief Retrieves the list of command-line options supported by this class.
    static po::options_description options();

  private:
    /// \brief Retrieves the list of files that will be processed.
    std::vector<std::string> files(const data::catalog &catalog) const;

    /// \brief Configures the job from a configuration file.
    void configure(const YAML::Node &node);

    /// \brief Configures the job from command-line options.
    void configure(const po::variables_map &varmap);

    /**
     * \brief Handles the subset of options available on the command line and in the configuration
     *        file.
     */
    template <class Container> void configure_common(const Container &);

    /// \brief Sets up the \c SIGINT handler.
    static void setup_sigint_handler(class job *job);

    /// \brief \c SIGINT handler.
    static void sigint_handler(int);
};

template <class Analyzer, class... Args> void job::run(Args &... args)
{
    using namespace logging;

    data::catalog catalog = sample().catalog();

    std::vector<std::string> files = job::files(catalog);
    if (files.empty()) {
        throw std::runtime_error("No file set for input.");
    } else {
        logging::info << "Initializing reader (this can take a while)..." << std::endl;

        util::chains chains(files);
        struct info info(catalog, chains);
        info.sample = sample();

        Analyzer ana(info, args...);

        long long count = info.reader.GetEntries(true);
        if (count == 0) {
            throw std::runtime_error(
                "Input files don't appear to contain data. Is your proxy valid?");
        }
        count = std::min(_max_events, count);

        if (count <= 0) {
            warn << "Running on zero event. Skipping event loop." << std::endl;
        } else {
            logging::info << "We will run on " << count << " events." << std::endl;

            if (_graceful_sigint) {
                setup_sigint_handler(this);
            }

            bool had_exception = false;

            timer time(count);
            time.start();
            for (long long entry = 0; entry < count; ++entry) {
                TTreeReader::EEntryStatus status = info.reader.SetEntry(entry);
                if (status != TTreeReader::kEntryValid) {
                    error << "Reader status code not valid: " << status << std::endl;
                    if (_fatal_exceptions) {
                        throw std::runtime_error("Reader status code not valid: " +
                                                 std::to_string(status));
                    }
                    continue;
                }

                try {
                    ana();
                } catch (std::exception &e) {
                    had_exception = true;
                    error << "Caught exception while processing events: " << e.what() << std::endl;
                    if (_fatal_exceptions) {
                        throw;
                    }
                }
                time.next();

                if (_sigint_caught) {
                    warn << "Caught SIGINT, exiting..." << std::endl;
                    break;
                }
            }
            time.stop();

            if (_graceful_sigint) {
                setup_sigint_handler(nullptr);
            }

            if (had_exception) {
                warn << "Caught exceptions while processing events. Output may not be complete."
                     << std::endl;
            }
        }

        logging::info << "Writing output into: " << output_dirname() << std::endl;
        std::shared_ptr<TFile> out = sample().histogram_file(
            _analyzer_name, output_dirname(), "RECREATE", _job_count > 1 ? _job_id : -1);
        out->cd();
        ana.write();
        out->Close();
        logging::info << "Done writing output." << std::endl;
    }
}
} // namespace util

#endif // JOB_H
