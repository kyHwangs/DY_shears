#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>
#include <TKey.h>
#include <TStyle.h>
#include <TVectorD.h>

#include "logging.h"
#include "options.h"
#include "sample.h"
#include "style_list.h"

struct config
{
    std::string input_dir_name;
    std::string output_dir_name;
    std::string output_format = "png";
    util::style_list style;
    bool preliminary = true;
    bool verticalNormalization = false;
};

config parse_options(int argc, char **argv, util::options &options);
void create_output_dir(const config &c);
std::vector<data::sample> load(const util::options &options);
void filter_histogram_names(const util::options &options,
                            const config &c,
                            std::set<std::string> &histogram_names);
double get_xsec_wsum(const data::sample &sample, std::shared_ptr<TFile> &file);
void prepare_axis_for_log(TAxis &axis);

int main(int argc, char **argv)
{
    try {
        // Parse options
        util::options options;
        config c = parse_options(argc, argv, options);

        // Load files
        auto samples = load(options);
        std::vector<std::shared_ptr<TFile>> files;
        for (auto &s : samples) {
            auto file = s.histogram_file("dyjets", c.input_dir_name);
            if (file == nullptr) {
                throw std::runtime_error("Could not open files");
            }
            files.push_back(file);
        }

        create_output_dir(c);

        // Style options
        gStyle->SetOptStat(0);
        gStyle->SetPaintTextFormat("4.0f");

        // Find response matrices
        std::set<std::string> matrices;
        for (auto &file : files) {
            for (const auto &&obj : *file->GetListOfKeys()) {
                if (auto key = dynamic_cast<const TKey *>(obj)) {
                    std::string name = key->GetName();
                    if (boost::ends_with(name, "-matrix")) {
                        matrices.insert(name);
                    }
                }
            }
        }

        filter_histogram_names(options, c, matrices);

        // Plot 'em all!
        for (const auto &name : matrices) {
            util::logging::debug << "Plotting matrix " << name << std::endl;

            // Get matrix
            std::unique_ptr<TH2> matrix = nullptr;
            for (std::size_t i = 0; i < files.size(); ++i) {
                auto &file = files[i];

                TH2 *m = nullptr;
                file->GetObject(name.c_str(), m);
                if (m == nullptr) {
                    util::logging::warn << "Could not read matrix "
                                        << name
                                        << " for sample "
                                        << samples[i].name()
                                        << std::endl;
                    continue;
                }

                std::cout << samples[i].name() << " "
                          << get_xsec_wsum(samples[i], file)
                          << " " << m
                          << " " << matrix.get()
                          << std::endl;
                m->Scale(get_xsec_wsum(samples[i], file));
                if (matrix == nullptr) {
                    matrix.reset(m); // Can't clone, seems to triggers an issue.
                } else {
                    matrix->Add(m);
                }
            }
            if (matrix == nullptr) {
                util::logging::warn << "Could not read any matrix: "
                                    << name
                                    << std::endl;
                continue;
            }

            // Normalize
            int nBinsX = matrix->GetNbinsX();
            int nBinsY = matrix->GetNbinsY();

            if (c.verticalNormalization) {
                for (int x = 0; x <= nBinsX + 1; x++) {
                    double row_total = 0;
                    for (int y = 0; y <= nBinsY + 1; y++) {
                        row_total += matrix->GetBinContent(x, y);
                    }
                    if (row_total != 0) {
                        for (int y = 0; y <= nBinsY + 1; y++) {
                            double contents = matrix->GetBinContent(x, y);
                            matrix->SetBinContent(x, y, 100 * contents / row_total);
                        }
                    }
                }
            } else {
                for (int y = 0; y <= nBinsY + 1; y++) {
                    double row_total = 0;
                    for (int x = 0; x <= nBinsX + 1; x++) {
                        row_total += matrix->GetBinContent(x, y);
                    }
                    if (row_total != 0) {
                        for (int x = 0; x <= nBinsX + 1; x++) {
                            double contents = matrix->GetBinContent(x, y);
                            matrix->SetBinContent(x, y, 100 * contents / row_total);
                        }
                    }
                }
            }

            // Draw
            TCanvas canvas("canvas", "", 800, 800);
            matrix->GetZaxis()->SetRangeUser(0, 100);
            matrix->Draw("colz text");

            // Styling
            bool logx = c.style.get<bool>("log x", name, false);
            if (logx) {
                canvas.SetLogx(true);
                canvas.SetLogy(true);
                prepare_axis_for_log(*matrix->GetXaxis());
                prepare_axis_for_log(*matrix->GetYaxis());
            }

            // Title
            std::string title = c.style.get_formatted("x axis label", name, "");
            if (!title.empty()) {
                auto precisions = c.style.get_formatted_all("x axis detail", name);
                if (!precisions.empty()) {
                    title += " (" + boost::algorithm::join(precisions, ", ") + ")";
                }
                matrix->SetTitle(("Response matrix for " + title).c_str());
            } else if (std::strlen(matrix->GetTitle()) == 0) {
                matrix->SetTitle(("Response matrix for " + name).c_str());
            }

            // Axis labels
            auto unit = c.style.get_formatted("x axis unit", name, "");
            if (!unit.empty()) {
                matrix->GetXaxis()->SetTitle(("Reconstructed level [" + unit + "]").c_str());
                matrix->GetYaxis()->SetTitle(("Generated level [" + unit + "]").c_str());
            } else {
                matrix->GetXaxis()->SetTitle("Reconstructed level");
                matrix->GetYaxis()->SetTitle("Generated level");
            }

            canvas.Print((c.output_dir_name + "/" + name + "." + c.output_format).c_str());
        }
    } catch (std::exception &e) {
        util::logging::fatal << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

config parse_options(int argc, char **argv, util::options &options)
{
    po::options_description optdesc;

    optdesc.add_options()("input,i",
                          po::value<std::string>()->default_value("dyjets-histograms"),
                          "Sets the directory to search for histogram files");
    optdesc.add_options()("output,o", po::value<std::string>(), "Sets the output directory");
    optdesc.add_options()("format,f",
                          po::value<std::string>()->default_value("png"),
                          "Sets the output format (png, pdf, ...)");
    optdesc.add_options()("histogram-name,n",
                          po::value<std::vector<std::string>>(),
                          "Produce the given histogram (can be used several times)");
    optdesc.add_options()("vertical", "Normalize on columns instead of rows");

    options.default_init(argc, argv, "dyjets.yml", {optdesc});

    config c;

    if (options.map.count("input") > 0) {
        c.input_dir_name = options.map["input"].as<std::string>();
    } else {
        c.input_dir_name = "dyjets-histograms";
    }

    if (options.map.count("output") > 0) {
        c.output_dir_name = options.map["output"].as<std::string>();
    } else {
        c.output_dir_name = c.input_dir_name + "/unfolding";
    }
    c.output_format = options.map["format"].as<std::string>();

    if (options.config["plots"]) {
        c.style = util::style_list(options.config["plots"]);
    }

    if (options.config["preliminary"]) {
        c.preliminary = options.config["preliminary"].as<bool>();
    }

    c.verticalNormalization = (options.map.count("vertical") > 0);

    return c;
}

std::vector<data::sample> load(const util::options &options)
{
    std::vector<data::sample> samples = data::sample::load(options);
    // Keep only MC
    samples.erase(std::remove_if(
            samples.begin(),
            samples.end(),
            [](const data::sample &s) { return s.type() != data::sample::mc; }
        ), samples.end());
    return samples;
}

double get_xsec_wsum(const data::sample &sample, std::shared_ptr<TFile> &file)
{
    // Read job info histograms
    TH1 *job_info = nullptr;
    file->GetObject("_job_info", job_info);
    if (job_info == nullptr) {
        throw std::runtime_error("File " + std::string(file->GetName()) +
                                 " doesn't have the _job_info histogram.");
    }

    TVectorD *job_info_average = nullptr;
    file->GetObject("_job_info_average", job_info_average);
    if (job_info_average == nullptr) {
        throw std::runtime_error("File " + std::string(file->GetName()) +
                                 " doesn't have the _job_info_average vector.");
    }

    double wsum = job_info->GetBinContent(2);
    double xsec = sample.xsec() > 0 ? sample.xsec() : (*job_info_average)[1];
    return xsec / wsum;
}

void create_output_dir(const config &c)
{
    using namespace boost::filesystem;

    // Create output directory if it doesn't exist
    if (!is_directory(c.output_dir_name)) {
        if (exists(c.output_dir_name)) {
            // "c.output_dir_name" exists and is not a directory...
            throw std::runtime_error("Path " + c.output_dir_name +
                                     " exists and is not a directory");
        } else {
            util::logging::info << "Creating directory " << c.output_dir_name << std::endl;
            create_directories(c.output_dir_name);
        }
    }
}

void filter_histogram_names(const util::options &options,
                            const config &c,
                            std::set<std::string> &histogram_names)
{
    if (options.map.count("histogram-name") > 0) {
        // Read from command line
        std::vector<std::string> names =
            options.map["histogram-name"].as<std::vector<std::string>>();

        // Check that names exist
        for (auto name : names) {
            if (std::find(histogram_names.begin(), histogram_names.end(), name) ==
                histogram_names.end()) {
                throw std::runtime_error("Histogram " + name + " doesn't exist");
            }
        }

        histogram_names.clear();
        std::copy(
            names.begin(), names.end(), std::inserter(histogram_names, histogram_names.begin()));
    } else {
        // Produce everything
        auto before_filter = histogram_names.size();

        util::logging::info << "Found " << histogram_names.size() << " histograms." << std::endl;

        // Remove histograms vetoed by style
        for (auto it = histogram_names.begin(); it != histogram_names.end();) {
            if (c.style.get<bool>("produce", *it, true)) {
                ++it;
            } else {
                util::logging::debug << "Not producing histogram " << *it << " due to plot rules."
                                     << std::endl;
                it = histogram_names.erase(it);
            }
        }

        auto removed = before_filter - histogram_names.size();
        util::logging::info << removed << " histograms skipped due to style rules." << std::endl;
    }
}

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
