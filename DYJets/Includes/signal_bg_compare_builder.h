#ifndef SIGNAL_BG_COMPARE_BUILDER_H
#define SIGNAL_BG_COMPARE_BUILDER_H

#include "compare_builder_base.h"

class TH1;

namespace util
{

class signal_bg_compare_builder : public compare_builder_base
{
    std::unique_ptr<data::data_comparison_entry> _data_entry;
    std::unique_ptr<data::mc_comparison_entry> _signal_entry;
    std::unique_ptr<data::mc_comparison_entry> _background_entry;

    std::unique_ptr<TH1> _ratio, _signal;

    double _lumi;

  public:
    /// \brief Constructor
    explicit signal_bg_compare_builder(const std::string &analyzer_name);

  protected:
    // Overriden from base class
    po::options_description options() const override;

    // Overriden from base class
    void load() override;

    // Overriden from base class
    void fill_upper_panel(const std::string &name) override;

    // Overriden from base class
    void fill_legend(TLegend &legend, const std::string &name) override;

    // Overriden from base class
    bool fill_lower_panel(const std::string &name) override;

    // Overriden from base class
    void override_lower_panel_settings(TPad &pad) override;

    // Overriden from base class
    void reset_drawing_state() override;
};

} // namespace util

#endif // SIGNAL_BG_COMPARE_BUILDER_H
