#include "reweighing.h"

#include "logging.h"

namespace physics {

reweighing::reweighing(util::job::info &info, const util::options &opt)
{
    const YAML::Node node = opt.config["reweighing"];
    if (!node) {
        // Not mandatory
        return;
    }

    // Get the name of the sample we are running on
    auto sample_name = info.sample.name();

    // npNLO reweighing
    if (node["npNLO"]) {
        auto tables = node["npNLO"].as<util::tables>();
        if (tables.count(sample_name) > 0) {
            util::logging::info << "This sample will be reweighed using table "
                                << node["npNLO"][sample_name].as<std::string>()
                                << std::endl;
            npNLO.reset(new TTreeReaderValue<int>(info.reader, "npNLO"));
            _table_npNLO = tables.at(sample_name);
        }
    }
}

void reweighing::reweigh(weights &w) const
{
    if (npNLO != nullptr) {
        // npNLO reweighing
        int npartons = **npNLO;
        // Get the weight
        double weight = _table_npNLO.getEfficiency(npartons, 0);
        // Use it
        w.use_gen_weight(weight);
        w.use_weight(weight);
    }
}

} // namespace physics
