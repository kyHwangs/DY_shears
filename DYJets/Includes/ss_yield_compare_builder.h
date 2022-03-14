#ifndef SS_YIELD_COMPARE_BUILDER_H
#define SS_YIELD_COMPARE_BUILDER_H

#include "compare_builder_base.h"

class TH1;

namespace util
{

class ss_yield_compare_builder : public compare_builder_base
{
    std::unique_ptr<data::data_comparison_entry> _os_data_entry;
    std::unique_ptr<data::data_comparison_entry> _ss_data_entry;
    std::unique_ptr<data::mc_comparison_entry> _ss_mc_entry;

    std::unique_ptr<TH1> _yield;

    double _lumi;

  public:
    /// \brief Constructor
    explicit ss_yield_compare_builder(const std::string &analyzer_name);

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

#endif // SS_YIELD_COMPARE_BUILDER_H
