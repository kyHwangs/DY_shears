#include "histo_set2D.h"

#include "logging.h"

namespace util
{

histo_set2D::histogram_type &histo_set2D::get(const std::string &name, const std::string &tag)
{
    auto it = _histograms.find({name, tag});
    if (it != _histograms.end()) {
        return it->second;
    } else {
        // Histogram doesn't exist (yet), create it
        std::string hname = combined_name(name, tag);
        // Try to create from style
        auto res = create_from_style(name, tag);
        if (res != _histograms.end()) {
            logging::debug << "Properties of histogram '" << hname << "' taken from stylesheet"
                           << std::endl;
            return res->second;
        } else {
            // Try to find a model
            auto model = _models.find(name);
            if (model != _models.end()) {
                // We rely on the model's copy constructor
                auto res = _histograms.emplace(std::make_pair(name, tag), model->second);
                return res.first->second;
            } else {
                throw std::logic_error("Histogram not declared: " + name);
            }
        }
    }
}

void histo_set2D::write()
{
    logging::debug << "Writing histograms..." << std::endl;
    for (auto &pair : _histograms) {
        pair.second.Write(combined_name(pair.first.first, pair.first.second).c_str());
    }
    logging::debug << _histograms.size() << " histograms written." << std::endl;
}

std::string histo_set2D::combined_name(const std::string &name, const std::string &tag)
{
    if (tag.empty()) {
        return name;
    }
    return name + "_" + tag;
}

auto histo_set2D::create_from_style(const std::string &name, const std::string &tag)
    -> decltype(_histograms)::iterator
{
    std::string fullname = combined_name(name, tag);
    std::string method = _style.get<std::string>("binning", fullname, "ignored");
    if (method == "ignored") {
        return _histograms.end();
    } else if (method == "uniform") {
        int bin_count = _style.get<int>("bin count", fullname, 100);
        double binning_min = _style.get<double>("binning min", fullname);
        double binning_max = _style.get<double>("binning max", fullname);
        auto res = _histograms.emplace(
            std::make_pair(name, tag),
            histogram_type(
                fullname.c_str(), fullname.c_str(), bin_count, binning_min, binning_max, bin_count, binning_min, binning_max));
        return res.first;
    } else if (method == "custom") {
        auto bins = _style.get<std::vector<double>>("bin edges", fullname);
        auto res = _histograms.emplace(
            std::make_pair(name, tag),
            histogram_type(fullname.c_str(), fullname.c_str(), bins.size() - 1, bins.data(), bins.size() - 1, bins.data()));
        return res.first;
    } else {
        throw std::runtime_error("Undefined binning method '" + method +
                                 "' matched by histogram '" + fullname + "'");
    }
}
} // namespace util
