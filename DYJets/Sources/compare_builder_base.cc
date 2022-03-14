#include "compare_builder_base.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <TAxis.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>

namespace util
{

compare_builder_base::compare_builder_base(const std::string &analyzer_name,
                                           const std::string &default_config_file) :
    _analyzer_name(analyzer_name),
    _default_config_file(default_config_file),
    _preliminary(true)
{
}

void compare_builder_base::parse_options(int argc, char **argv)
{
    auto options = this->options();

    options.add_options()("output,o", po::value<std::string>(), "Sets the output directory");
    options.add_options()("format,f",
                          po::value<std::string>()->default_value("png"),
                          "Sets the output format (png, pdf, ...)");
    options.add_options()("histogram-name,n",
                          po::value<std::vector<std::string>>(),
                          "Produce the given histogram (can be used several times)");
    options.add_options()("lin", "Use a linear scale for the y axis (the default is a log scale)");

    _opt.default_init(argc, argv, _default_config_file, {options});

    if (_opt.map.count("output") > 0) {
        _output_dir_name = _opt.map["output"].as<std::string>();
    }
    _output_format = _opt.map["format"].as<std::string>();

    if (_opt.config["plots"]) {
        _style = util::style_list(_opt.config["plots"]);
    }

    if (_opt.config["preliminary"]) {
        _preliminary = _opt.config["preliminary"].as<bool>();
    }
}

void compare_builder_base::build()
{
    load();
    create_output_dir();
    filter_histogram_names();

    bool log = (_opt.map.count("lin") == 0);

    for (const std::string &name : _histogram_names) {
        util::logging::debug << "Producing histogram: " << name << std::endl;

        // Apply style
        _current_histo_name = name;
        _logx = _style.get<bool>("log x", name, false);

        TCanvas canvas(name.c_str(), "", 700, 900);

        // Upper panel
        TPad upper("upper", "upper", 0, 0.3, 1, 1);
        upper.SetTopMargin(0.11);
        upper.SetRightMargin(0.03);
        upper.SetTicks();
        if (log) {
            upper.SetLogy();
        }
        upper.Draw();
        upper.cd();

        fill_upper_panel(name);

        // CMS label
        TLatex cms;
        cms.SetTextSize(0.04);
        cms.SetTextFont(42);
        cms.SetTextAlign(kHAlignLeft + kVAlignBottom);
        cms.SetNDC();
        cms.SetText(0.1,
                    0.9,
                    _preliminary ? "#bf{CMS} #it{Preliminary}" : "#bf{CMS}");
        cms.Draw();

        // Lumi label
        double lumi = get_lumi();
        TLatex label;
        if (lumi > 0) {
            label.SetTextSize(0.04);
            label.SetTextFont(42);
            label.SetTextAlign(kHAlignRight + kVAlignBottom);
            label.SetNDC();

            std::stringstream ss;
            ss << std::setprecision(3) << (lumi / 1000);
            ss << " fb^{-1} (13 TeV)";
            label.SetText(0.97, 0.9, ss.str().c_str());
            label.Draw();
        }

        // Legend
        TLegend legend(0.63, 0.60, 0.81, 0.87);
        legend.SetTextSize(0.042);
        legend.SetFillStyle(0);
        legend.SetBorderSize(0);
        legend.SetTextFont(42);
        legend.Draw();

        fill_legend(legend, name);

        // Get back to the canvas
        canvas.cd();

        // Lower panel
        TPad lower("lower", "lower", 0, 0.05, 1, 0.3);
        lower.SetTopMargin(0.);
        lower.SetBottomMargin(0.3);
        lower.SetRightMargin(0.03);
        lower.SetGridy();
        lower.SetTicks();
        lower.Draw();
        lower.cd();

        override_lower_panel_settings(lower);

        if (fill_lower_panel(name)) {
            upper.SetBottomMargin(0.);
        }

        // Apply style
        if (_logx) {
            upper.SetLogx();
            lower.SetLogx();
        }
        // Write file
        canvas.Print(
            (_output_dir_name + "/" + name + "." + _output_format).c_str());

        // Cleanup
        reset_drawing_state();

    }

}

std::unique_ptr<data::data_comparison_entry> compare_builder_base::load_data(
        const std::string &input_dir)
{
    data::sample data;
    std::vector<data::sample> samples = data::sample::load(_opt);
    for (data::sample &s : samples) {
        if (s.name() == "data") {
            data = s;
        }
    }

    auto ptr = std::make_unique<data::data_comparison_entry>(_analyzer_name, data, input_dir);
    ptr->add_histograms(_histogram_names);
    return ptr;
}

std::unique_ptr<data::mc_comparison_entry> compare_builder_base::load_mc(
        const std::string &input_dir,
        bool keep_signal,
        bool keep_background)
{
    auto ptr = std::make_unique<data::mc_comparison_entry>(_opt,
                                                           _analyzer_name,
                                                           input_dir,
                                                           keep_signal,
                                                           keep_background);
    ptr->add_histograms(_histogram_names);
    return ptr;
}

namespace /* anonymous */ {
    /**
     * \brief Tunes an axis to be used on a log scale.
     *
     * If the first bin extends to 0, it is modified to start from a higher
     * value. This is needed because otherwise ROOT will make it span half of
     * the plots.
     *
     * \note Histograms with a uniform binning are not supported.
     */
    void prepare_axis_for_log(TAxis &axis)
    {
        if (axis.GetXmin() > 0 || axis.GetNbins() < 2) {
            // Nothing to do
            return;
        }

        // Get the bin edges in a safe container (std::vector over TArray).
        const double *edges_ptr = axis.GetXbins()->GetArray();
        if (edges_ptr == nullptr) {
            // Happens when the histogram has a uniform binning.
            return;
        }
        std::vector<double> edges(edges_ptr, edges_ptr + axis.GetNbins() + 1);

        // First, we get the extent of the second bin. It's proportional to the
        // ratio of the edges.
        double ratio = edges[2] / edges[1];

        // Then, we modify the first bin boundaries so that it gets the same
        // displayed size.
        edges[0] = edges[1] / ratio;

        // We want the lower bound to be a round number, so we round it. It's
        // more complicated that a simple floor() because we want 0.25 to become
        // 0.2 and not 0.0.
        double logfactor = std::pow(10, std::ceil(std::log10(edges[0])) - 1);
        edges[0] = logfactor * std::floor(edges[0] / logfactor);
        edges[0] *= 1.001; // Avoid tick labels.

        // Modify the axis to use the new bin edges.
        axis.Set(edges.size() - 1, edges.data());
    }
} // namespace anonymous

void compare_builder_base::format_upper_x_axis(TAxis &axis) const
{
    if (_logx) {
        prepare_axis_for_log(axis);
    }
}

void compare_builder_base::format_upper_y_axis(TAxis &axis, const std::string &title) const
{
    axis.SetLabelSize(0.04);
    axis.SetLabelOffset(0.002);
    axis.SetTitle(title.c_str());
    axis.SetTitleSize(0.04);
    axis.SetTitleOffset(1.32);
}

void compare_builder_base::format_lower_x_axis(TAxis &axis) const
{
    axis.SetTickLength(0.03);
    axis.SetTitleSize(0.1);
    axis.SetTitleOffset(1.2);
    axis.SetLabelSize(0.10);
    axis.SetLabelOffset(0.017);

    if (_logx) {
        prepare_axis_for_log(axis);
    }

    std::string label = _style.get_formatted("x axis label", _current_histo_name, "");
    if (!label.empty()) {
        auto precisions = _style.get_formatted_all("x axis detail", _current_histo_name);
        if (!precisions.empty()) {
            label += " (" + boost::algorithm::join(precisions, ", ") + ")";
        }
        auto unit = _style.get_formatted("x axis unit", _current_histo_name, "");
        if (!unit.empty()) {
            label += " [" + unit + "]";
        }
        axis.SetTitle(label.c_str());
    } else if (std::strlen(axis.GetTitle()) == 0) {
        axis.SetTitle(_current_histo_name.c_str());
    }
}

void compare_builder_base::format_lower_y_axis(TAxis &axis, const std::string &title) const
{
    axis.SetNdivisions(5, 5, 0);
    axis.SetTitle(title.c_str());
    axis.SetTitleSize(0.1);
    axis.SetTitleOffset(0.5);
    axis.CenterTitle();
    axis.SetLabelSize(0.08);
}

void compare_builder_base::create_output_dir() const
{
    using namespace boost::filesystem;

    // Create output directory if it doesn't exist
    if (!is_directory(_output_dir_name)) {
        if (exists(_output_dir_name)) {
            // "_output_dir_name" exists and is not a directory...
            throw std::runtime_error("Path " + _output_dir_name +
                                     " exists and is not a directory");
        } else {
            util::logging::info << "Creating directory " << _output_dir_name << std::endl;
            create_directories(_output_dir_name);
        }
    }
}

void compare_builder_base::filter_histogram_names()
{
    if (_opt.map.count("histogram-name") > 0) {
        // Read from command line
        std::vector<std::string> names =
            _opt.map["histogram-name"].as<std::vector<std::string>>();

        // Check that names exist
        for (auto name : names) {
            if (std::find(_histogram_names.begin(),
                          _histogram_names.end(),
                          name) == _histogram_names.end()) {
                throw std::runtime_error("Histogram " + name + " doesn't exist");
            }
        }

        _histogram_names.clear();
        std::copy(names.begin(),
                  names.end(),
                  std::inserter(_histogram_names, _histogram_names.begin()));
    } else {
        // Produce everything
        auto before_filter = _histogram_names.size();

        util::logging::info << "Found " << _histogram_names.size() << " histograms."
                            << std::endl;

        // Remove histograms vetoed by style
        for (auto it = _histogram_names.begin(); it != _histogram_names.end(); ) {
            if (_style.get<bool>("produce", *it, true)) {
                ++it;
            } else {
                util::logging::debug << "Not producing histogram " << *it
                                     << " due to plot rules." << std::endl;
                it = _histogram_names.erase(it);
            }
        }

        auto removed = before_filter - _histogram_names.size();
        util::logging::info << removed << " histograms skipped due to style rules."
                            << std::endl;
    }
}

} // namespace util
