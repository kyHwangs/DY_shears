#include "signal_bg_compare_builder.h"

#include <TPad.h>

namespace util
{

signal_bg_compare_builder::signal_bg_compare_builder(const std::string &analyzer_name) :
    compare_builder_base(analyzer_name, analyzer_name + ".yml")
{
}

po::options_description signal_bg_compare_builder::options() const
{
    po::options_description options = po::options_description("Comparison options");
    options.add_options()("input,i",
                          po::value<std::string>()->default_value("dyjets-histograms"),
                          "Sets the directory to search for histogram files");
    return options;
}

void signal_bg_compare_builder::load()
{
    std::string input_dir = parsed_options().map["input"].as<std::string>();
    set_default_output_dir(input_dir + "/signal-over-background");

    _data_entry = load_data(input_dir);
    _signal_entry = load_mc(input_dir, true, false);
    _background_entry = load_mc(input_dir, false, true);

    _lumi = _data_entry->lumi();
    util::logging::info << "Normalizing MC to " << (_lumi / 1000) << " fb^-1" << std::endl;
}

void signal_bg_compare_builder::fill_upper_panel(const std::string &name)
{
    _signal = _signal_entry->get(name, _lumi);

    if (_signal != nullptr) {
        format_upper_y_axis(*_signal->GetYaxis());
        _signal->SetMarkerStyle(20);
        _signal->SetMarkerColor(kGreen + 3);
        _signal->SetLineColor(kGreen + 3);
        _signal->SetStats(0);

        _signal->Draw("e2");
        _background_entry->draw(name, _lumi, true);
        _signal->Draw("epsame");

        if (auto axis = _background_entry->get_x_axis(name, _lumi)) {
            format_upper_x_axis(*axis);
        }
        format_upper_x_axis(*_signal->GetXaxis());
    } else {
        _background_entry->draw(name, _lumi);

        if (auto axis = _background_entry->get_x_axis(name, _lumi)) {
            format_upper_x_axis(*axis);
        }
    }
}

void signal_bg_compare_builder::fill_legend(TLegend &legend, const std::string &name)
{
    _signal_entry->add_to_legend(legend, name, _lumi);
    _background_entry->add_to_legend(legend, name, _lumi);
}

bool signal_bg_compare_builder::fill_lower_panel(const std::string &name)
{
    _ratio = _signal_entry->get(name, _lumi);
    std::unique_ptr<TH1> den = _background_entry->get(name, _lumi);

    if (_ratio == nullptr || den == nullptr) {
        return false;
    }

    _ratio->Divide(den.get());

    format_lower_x_axis(*_ratio->GetXaxis());
    format_lower_y_axis(*_ratio->GetYaxis(), "Signal/Background");

    _ratio->SetMarkerStyle(20);
    _ratio->SetMarkerColor(kGreen + 3);
    _ratio->SetLineColor(kGreen + 3);
    _ratio->SetStats(0);
    _ratio->SetTitle("");
    _ratio->Draw("ep");

    return true;
}

void signal_bg_compare_builder::override_lower_panel_settings(TPad &pad)
{
    pad.SetLogy();
}

void signal_bg_compare_builder::reset_drawing_state()
{
    _signal_entry->reset_drawing_state();
    _background_entry->reset_drawing_state();
    _ratio.reset();
    _signal.reset();
}

} // namespace util
