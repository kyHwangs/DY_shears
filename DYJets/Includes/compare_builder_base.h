#ifndef COMPARE_BUILDER_BASE_h
#define COMPARE_BUILDER_BASE_h

#include <memory>
#include <set>
#include <string>

#include "comparison_entry.h"
#include "style_list.h"

class TAxis;
class TLegend;
class TPad;

namespace po = boost::program_options;

namespace util
{

class compare_builder_base
{
    std::string _analyzer_name;
    std::string _default_config_file;
    std::string _output_dir_name;
    std::string _output_format;
    std::set<std::string> _histogram_names;

    util::options _opt;
    util::style_list _style;
    bool _preliminary;

    std::string _current_histo_name;
    bool _logx;

  public:
    /// \brief Constructor
    explicit compare_builder_base(const std::string &analyzer_name,
                                  const std::string &default_config_file);

    /// \brief Destructor
    virtual ~compare_builder_base() = default;

    /// \brief Parses options passed to the program
    void parse_options(int argc, char **argv);

    /**
     * \brief Sets the default output directory
     *
     * It can be overriden using command line option \c -o.
     */
    void set_default_output_dir(const std::string &output)
    {
        if (_output_dir_name.empty()) {
            _output_dir_name = output;
        }
    }

    /// \brief Make all plots
    void build();

    /// \brief Gives access to the options object
    const util::options &parsed_options() const { return _opt; }

  protected:
    /// \brief Loads data samples
    std::unique_ptr<data::data_comparison_entry> load_data(const std::string &input_dir);

    /// \brief Loads MC samples
    std::unique_ptr<data::mc_comparison_entry> load_mc(const std::string &input_dir,
                                                       bool keep_signal = true,
                                                       bool keep_background = true);

    /// \brief Gives access to the style
    const util::style_list &style() const { return _style; }

    /// \brief Formats the upper panel's x axis
    void format_upper_x_axis(TAxis &axis) const;

    /// \brief Formats the upper panel's y axis
    void format_upper_y_axis(TAxis &axis, const std::string &title = "# Events") const;

    /// \brief Formats the lower panel's x axis
    void format_lower_x_axis(TAxis &axis) const;

    /// \brief Formats the lower panel's y axis
    void format_lower_y_axis(TAxis &axis, const std::string &title = "") const;

    /// \brief Returns the options supported by this builder
    virtual po::options_description options() const = 0;

    /// \brief Loads all needed \ref comparison_entry
    virtual void load() = 0;

    /// \brief Gets the luminosity
    virtual double get_lumi() const { return 0; }

    /// \brief Fills the upper panel with plots
    virtual void fill_upper_panel(const std::string &name) = 0;

    /// \brief Fills the legend with plots
    virtual void fill_legend(TLegend &legend, const std::string &name) = 0;

    /// \brief Fills the lower panel with plots
    virtual bool fill_lower_panel(const std::string &name) = 0;

    /// \brief Can be used to override lower \c TPad settings
    virtual void override_lower_panel_settings(TPad &) {}

    /// \brief Resets drawing state after a plot was made
    virtual void reset_drawing_state() = 0;

  private:
    /// \brief Creates the directory results are written to
    void create_output_dir() const;

    /// \brief Filters the \c _histogram_names set to only contain plots that should be created
    void filter_histogram_names();
};

} // namespace util

#endif // COMPARE_BUILDER_BASE_h
