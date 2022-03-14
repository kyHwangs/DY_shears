#include "job.h"

#include <algorithm>
#include <csignal>
#include <string>
#include <vector>

namespace util
{

job::info::info(const data::catalog &catalog, const util::chains &chains)
    : catalog(catalog),
      chains(chains),
      reader(this->chains.events().get())
{
}

job::job(const std::string &analyzer_name)
    : _graceful_sigint(isatty(fileno(stdin)) || isatty(fileno(stdout)) || isatty(fileno(stderr))),
      _sigint_caught(false),
      _analyzer_name(analyzer_name)
{
}

std::vector<std::string> job::files(const data::catalog &catalog) const
{
    std::vector<std::string> files = catalog.files();

    std::size_t total = files.size();
    std::size_t begin = total * _job_id / _job_count;
    std::size_t end = std::min(begin + _max_files, total * (_job_id + 1) / _job_count);

    assert(end <= files.size());

    return std::vector<std::string>(files.begin() + begin, files.begin() + end);
}

data::sample job::sample() const
{
    auto it = std::find_if(_samples.begin(), _samples.end(), [this](const data::sample &s) {
        return s.name() == _sample_name;
    });
    if (it == _samples.end()) {
        throw std::runtime_error("Sample \"" + _sample_name + "\" not found");
    } else {
        return *it;
    }
}

std::string job::output_dirname() const
{
    if (!_output_dir.empty()) {
        return _output_dir;
    }

    std::string name = _analyzer_name + "-histograms";
    if (_max_files < std::numeric_limits<int>::max()) {
        name += "-max-files-" + std::to_string(_max_files);
    }
    if (_max_events < std::numeric_limits<long long>::max()) {
        name += "-max-events-" + std::to_string(_max_events);
    }
    return name;
}

void job::configure(const class options &opt)
{
    _samples = data::sample::load(opt);

    if (opt.config["job"]) {
        configure(opt.config["job"]);
    }
    configure(opt.map);
}

void job::configure(const YAML::Node &node) { configure_common(node); }

void job::configure(const po::variables_map &varmap)
{
    configure_common(varmap);
    util::set_value_safe(
        varmap, _job_count, "job-count", "number of jobs", [](int val) { return val > 0; });
    util::set_value_safe(varmap, _job_id, "job-id", "job id", [&](int val) -> bool {
        return val >= 0 && val < _job_count;
    });
    _sample_name = varmap["sample"].as<std::string>();
    if (varmap.count("output-dir") > 0) {
        _output_dir = varmap["output-dir"].as<std::string>();
    }
}

template <class Container> void job::configure_common(const Container &container)
{
    util::set_value_safe(container, _max_files, "max files", "maximum number of files");
    if (_max_files < 0) {
        _max_files = std::numeric_limits<int>::max();
    }
    util::set_value_safe(container, _max_events, "max events", "maximum number of events");
    if (_max_events < 0) {
        _max_events = std::numeric_limits<long long>::max();
    }
}

po::options_description job::options()
{
    po::options_description options("Job control options");
    options.add_options()(
        "sample,s", po::value<std::string>()->default_value("data"), "Name of the sample to use");
    options.add_options()("output-dir,o", po::value<std::string>(), "Output directory");
    options.add_options()("max-events",
                          po::value<long long>()->default_value(-1),
                          "Maximum number of events to read (-1 for no limit)");
    options.add_options()("max-files",
                          po::value<int>()->default_value(-1),
                          "Maximum number of files to read (-1 for no limit)");
    options.add_options()(
        "job-id", po::value<int>()->default_value(0), "Job id (useful when running on batch)");
    options.add_options()("job-count",
                          po::value<int>()->default_value(1),
                          "Number of jobs (useful when running on batch)");
    return options;
}

namespace
{
static job *sigint_handler_target = nullptr;
}

void job::setup_sigint_handler(class job *job)
{
    if (job == nullptr) {
        logging::debug << "Unsetting SIGINT handler for " << sigint_handler_target << std::endl;
        sigint_handler_target = nullptr;
        std::signal(SIGINT, SIG_DFL);
        return;
    }
    logging::debug << "Setting up SIGINT handler for " << job << std::endl;
    if (sigint_handler_target != nullptr) {
        logging::error << "SIGINT handler already defined." << std::endl;
        return;
    }
    sigint_handler_target = job;
    if (job != nullptr) {
        std::signal(SIGINT, sigint_handler);
    }
}

void job::sigint_handler(int)
{
    if (sigint_handler_target != nullptr) {
        sigint_handler_target->_sigint_caught = true;
    }
}
} // namespace util
