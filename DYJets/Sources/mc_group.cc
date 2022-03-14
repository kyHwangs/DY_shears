#include "mc_group.h"

#include <stdexcept>

#include <TH1.h>

#include "comparison_entry.h"
#include "logging.h"

namespace data
{

void mc_group::add_histograms(std::set<std::string> &histos)
{
    for (sample_data &sd : _sample_data) {
        if (sd.centry == nullptr) {
            continue;
        }
        sd.centry->add_histograms(histos);
    }
}

std::unique_ptr<TH1> mc_group::get(const std::string &name)
{
    std::unique_ptr<TH1> res = nullptr;
    for (sample_data &sd : _sample_data) {
        if (sd.centry == nullptr) {
            continue;
        }

        std::unique_ptr<TH1> histo = sd.centry->get(name, 1);
        if (histo == nullptr) {
            continue;
        }

        if (res == nullptr) {
            res = std::move(histo);
        } else {
            res->Add(histo.get());
        }

        sd.centry->reset_drawing_state();
    }
    if (res != nullptr) {
        res->SetStats(0);
        res->SetFillStyle(1001);
        res->SetFillColor(_color);
        res->SetLineColor(_color);
        res->Scale(_scale_factor);
    }
    return res;
}

void mc_group::init(const std::vector<sample> &all_samples,
                    const std::string &analyzer_name,
                    const std::string &input_dir,
                    bool open_files)
{
    bool all_found = true;
    for (const std::string &name : _sample_names) {
        sample_data sd;

        // Find sample
        bool found = false;
        for (const sample &sample : all_samples) {
            if (name == sample.name()) {
                sd.sample = sample;
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("Sample " + name + " (required by MC group " + _legend +
                                     ") doesn't exist.");
        }

        if (open_files) {
            // Open file
            try {
                sd.centry = std::make_shared<data_comparison_entry>(analyzer_name, sd.sample, input_dir);
            } catch (std::runtime_error e) {
                util::logging::warn << "File not found for sample " << name << std::endl;
                all_found = false;
            }
        }

        _sample_data.push_back(sd);
    }
    if (!all_found) {
        if (_required) {
            throw std::runtime_error("Missing files for group " + _legend);
        } else {
            util::logging::warn << "Sample " << _legend << " is not complete" << std::endl;
        }
    }
}

std::vector<sample> mc_group::samples() const
{
    std::vector<sample> ret;
    for (const sample_data &sd : _sample_data) {
        ret.push_back(sd.sample);
    }
    return ret;
}

std::vector<mc_group> mc_group::load(const util::options &opt,
                                     const std::string &analyzer_name,
                                     const std::string &input_dir,
                                     const std::vector<sample> &all_samples,
                                     bool open_files)
{
    std::vector<mc_group> groups = opt.config["MC grouping"].as<std::vector<mc_group>>();
    for (mc_group &g : groups) {
        g.init(all_samples, analyzer_name, input_dir, open_files);
    }
    return groups;
}
} // namespace data

/// \cond
namespace YAML
{

template <> struct convert<data::mc_group>
{
    static bool decode(const Node &node, data::mc_group &group)
    {
        if (!node["legend"]) {
            throw std::runtime_error("MC group legend is not set");
        }
        group._legend = node["legend"].as<std::string>();

        if (!node["color"]) {
            throw std::runtime_error("MC group color is not set for " + group._legend);
        }
        group._color = node["color"].as<int>();

        if (node["required"]) {
            group._required = node["required"].as<bool>();
        }

        if (node["scale factor"]) {
            group._scale_factor = node["scale factor"].as<double>();
        }

        if (!node["samples"]) {
            throw std::runtime_error("MC group list of samples is not set for " + group._legend);
        }
        group._sample_names = node["samples"].as<std::vector<std::string>>();

        if (node["is signal"] && node["is signal"].as<bool>()) {
            group._is_signal = true;
        }

        return true;
    }
};
} // namespace YAML
/// \endcond
