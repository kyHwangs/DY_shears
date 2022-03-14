#include "ss_yield_compare_builder.h"

#include <TPad.h>

namespace util
{

ss_yield_compare_builder::ss_yield_compare_builder(const std::string &analyzer_name) :
    compare_builder_base(analyzer_name, analyzer_name + ".yml")
{
}

po::options_description ss_yield_compare_builder::options() const
{
    po::options_description options = po::options_description("Comparison options");
    options.add_options()("os",
                          po::value<std::string>(),
                          "Sets the directory to search for OS histogram files");
    options.add_options()("ss",
                          po::value<std::string>(),
                          "Sets the directory to search for SS histogram files");
    return options;
}

void ss_yield_compare_builder::load()
{
    const auto &opt = parsed_options();

    if (opt.map.count("os") != 1) {
        throw std::runtime_error("'os' argument must be given exactly once");
    }
    if (opt.map.count("ss") != 1) {
        throw std::runtime_error("'ss' argument must be given exactly once");
    }

    std::string os_input_dir = opt.map["os"].as<std::string>();
    std::string ss_input_dir = opt.map["ss"].as<std::string>();

    set_default_output_dir(os_input_dir + "/ss-yield");

    _os_data_entry = load_data(os_input_dir);
    _ss_data_entry = load_data(ss_input_dir);
    _ss_mc_entry = load_mc(ss_input_dir);

    _lumi = _os_data_entry->lumi();
    if (_lumi != _ss_data_entry->lumi()) {
        throw std::runtime_error("Trying to compare folders with different luminsities");
    }
    util::logging::info << "Normalizing MC to " << (_lumi / 1000) << " fb^-1" << std::endl;
}

void ss_yield_compare_builder::fill_upper_panel(const std::string &name)
{
    _ss_data_entry->draw(name, _lumi);
    _ss_mc_entry->draw(name, _lumi, true);
    _ss_data_entry->draw(name, _lumi, true);

    if (auto axis = _ss_data_entry->get_x_axis(name, _lumi)) {
        format_upper_x_axis(*axis);
    }
    if (auto axis = _ss_mc_entry->get_x_axis(name, _lumi)) {
        format_upper_x_axis(*axis);
    }
}

void ss_yield_compare_builder::fill_legend(TLegend &legend, const std::string &name)
{
    _ss_data_entry->add_to_legend(legend, name, _lumi);
    _ss_mc_entry->add_to_legend(legend, name, _lumi);
}

bool ss_yield_compare_builder::fill_lower_panel(const std::string &name)
{
    _yield = _ss_data_entry->get(name, _lumi);
    std::unique_ptr<TH1> ss_mc = _ss_mc_entry->get(name, _lumi);
    std::unique_ptr<TH1> os_data = _os_data_entry->get(name, _lumi);

    if (_yield == nullptr || os_data == nullptr) {
        return false;
    }
    if (ss_mc != nullptr) {
        _yield->Add(ss_mc.get(), -1);
    }
    _yield->Divide(os_data.get());

    format_lower_x_axis(*_yield->GetXaxis());
    format_lower_y_axis(*_yield->GetYaxis(), "(SS Data - SS MC)/OS Data");

    _yield->SetMarkerStyle(20);
    _yield->SetMarkerColor(kBlack);
    _yield->SetLineColor(kBlack);
    _yield->SetStats(0);
    _yield->SetTitle("");
    _yield->Draw("ep");

    return true;
}

void ss_yield_compare_builder::override_lower_panel_settings(TPad &pad)
{
    pad.SetLogy();
}

void ss_yield_compare_builder::reset_drawing_state()
{
    _ss_data_entry->reset_drawing_state();
    _ss_mc_entry->reset_drawing_state();
    _os_data_entry->reset_drawing_state();
    _yield.reset();
}

} // namespace util
