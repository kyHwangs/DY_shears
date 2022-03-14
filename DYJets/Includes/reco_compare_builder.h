#ifndef RECO_COMPARE_BUILDER_H
#define RECO_COMPARE_BUILDER_H

#include "compare_builder_base.h"

class TH1;

namespace util
{

class reco_compare_builder : public compare_builder_base
{
    std::unique_ptr<data::data_comparison_entry> _data_entry;
    std::unique_ptr<data::mc_comparison_entry> _mc_entry;

    std::unique_ptr<TH1> _ratio;

    bool _reversed;

    double _lumi;

  public:
    /// \brief Constructor
    explicit reco_compare_builder(const std::string &analyzer_name);

  protected:
    // Overriden from base class
    po::options_description options() const override;

    // Overriden from base class
    void load() override;

    // Overriden from base class
    double get_lumi() const override { return _lumi; }

    // Overriden from base class
    void fill_upper_panel(const std::string &name) override;

    // Overriden from base class
    void fill_legend(TLegend &legend, const std::string &name) override;

    // Overriden from base class
    bool fill_lower_panel(const std::string &name) override;

    // Overriden from base class
    void reset_drawing_state() override;
};

} // namespace util

#endif // RECO_COMPARE_BUILDER_H
