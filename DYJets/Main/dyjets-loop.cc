#include "catalog.h"
#include "dyjets_analyzer.h"
#include "job.h"
#include "logging.h"
#include "options.h"
#include "sample.h"
#include "timer.h"

int main(int argc, char **argv)
{
    util::options opt;
    try {
        opt.default_init(argc, argv, "dyjets.yml", {dyjets_analyzer::options(), util::job::options()});

        util::job j("dyjets");
        j.configure(opt);
        j.run<dyjets_analyzer>(opt);

    } catch (std::exception &e) {
        util::logging::fatal << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
