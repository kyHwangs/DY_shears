#include "options.h"
#include "sample.h"

/// \brief Retrieves the list of command-line options supported by this program.
po::options_description options();

int main(int argc, char **argv)
{
    util::options opt;
    try {
        opt.default_init(argc, argv, "higgs.yml", {options()});

        const std::vector<data::sample> samples = data::sample::load(opt);

        for (const data::sample &sample : samples) {
            if (opt.map.count("disable-data") > 0 && sample.type() == data::sample::data) {
                continue;
            }
            if (opt.map.count("disable-mc") > 0 && sample.type() == data::sample::mc) {
                continue;
            }
            if (opt.map.count("disable-background") > 0 &&
                sample.type() == data::sample::background) {
                continue;
            }

            for (unsigned job = 0; job < sample.jobs(); ++job) {
                std::cout << "bin/higgs-loop -s " << sample.name() << " -o {results}";
                if (opt.map.count("config") > 0) {
                    std::cout << " -c " << opt.config_file();
                }
                if (sample.jobs() > 1) {
                    std::cout << " --job-id " << job;
                    std::cout << " --job-count " << sample.jobs();
                }
                std::cout << std::endl;
            }
        }
    } catch (std::exception &e) {
        util::logging::fatal << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

po::options_description options()
{
    po::options_description options("Sample selection options");
    options.add_options()("disable-data", "Don't include jobs for data");
    options.add_options()("disable-mc", "Don't include jobs for MC");
    options.add_options()("disable-background", "Don't include jobs for background");
    return options;
}
