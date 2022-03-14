
#include "logging.h"
#include "reco_compare_builder.h"
#include "signal_bg_compare_builder.h"
#include "ss_yield_compare_builder.h"

void usage(const std::string &program_name);

int main(int argc, char **argv)
{
    try {
        std::unique_ptr<util::compare_builder_base> builder;

        std::string tool = (argc >= 2 ? argv[1] : "reco-level-agreement");
        if (tool == "-h" || tool == "--help") {
            usage(argv[0]);
            return EXIT_SUCCESS;
        } else if (tool[0] == '-' /* option */ || tool == "reco-level-agreement") {
            builder = std::make_unique<util::reco_compare_builder>("dyjets");
        } else if (tool == "signal-over-background") {
            builder = std::make_unique<util::signal_bg_compare_builder>("dyjets");
        } else if (tool == "ss-yield") {
            builder = std::make_unique<util::ss_yield_compare_builder>("dyjets");
        } else {
            throw std::runtime_error("Unknown tool '" + tool + "'");
        }

        builder->parse_options(argc, argv);
        builder->build();

    } catch (std::exception &e) {
        util::logging::fatal << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

void usage(const std::string &program_name)
{
    std::cerr << "Usage: " << program_name << " [tool] [options...]" << std::endl
              << std::endl
              << "Available tools:" << std::endl
              << "\treco-level-agreement (default)" << std::endl
              << "\tsignal-over-background" << std::endl
              << "\tss-yield" << std::endl
              << std::endl
              << "Use " << program_name
              << " <tool> --help for the corresponding list of options." << std::endl;
}
