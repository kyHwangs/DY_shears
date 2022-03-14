#include "comparison_entry.h"

#include <cmath>

#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TKey.h>
#include <TLegend.h>
#include <TList.h>
#include <TVectorD.h>

#include "logging.h"
#include "sample.h"

namespace data
{

mc_comparison_entry::mc_comparison_entry(const util::options &opt,
                                         const std::string &analyzer_name,
                                         const std::string &input_dir,
                                         bool keep_signal,
                                         bool keep_background)
{
    std::vector<data::sample> samples = data::sample::load(opt);
    _groups = data::mc_group::load(opt, analyzer_name, input_dir, samples);
    if (!keep_signal) {
        _groups.erase(std::remove_if(_groups.begin(),
                                     _groups.end(),
                                     [](const mc_group &g) { return g.is_signal(); }),
                      _groups.end());
    }
    if (!keep_background) {
        _groups.erase(std::remove_if(_groups.begin(),
                                     _groups.end(),
                                     [](const mc_group &g) {
                                         return !g.is_signal(); }),
                      _groups.end());
    }
}

void mc_comparison_entry::add_histograms(std::set<std::string> &histos)
{
    for (data::mc_group &group : _groups) {
        group.add_histograms(histos);
    }
}

void mc_comparison_entry::add_to_legend(TLegend &legend, const std::string &plotname, double lumi)
{
    if (_stack == nullptr) {
        create_stack(plotname, lumi);
        if (_stack == nullptr) {
            return;
        }
    }
    if (_stack->GetNhists() == 0) {
        return;
    }
    TList *histograms = _stack->GetHists(); // Not owned
    if ((unsigned) histograms->GetSize() != _legend.size()) {
        throw std::logic_error("Inconsistent sizes in mc_comparison_entry::add_to_legend");
    }
    for (int i = _legend.size() - 1; i >= 0; --i) {
        legend.AddEntry(dynamic_cast<TH1 *>(histograms->At(i)), _legend[i].c_str());
    }
}

void mc_comparison_entry::draw(const std::string &name, double lumi, bool same)
{
    if (_stack == nullptr) {
        create_stack(name, lumi);
        if (_stack == nullptr) {
            return;
        }
    }
    _stack->Draw(same ? "hist same" : "hist");
    _stack->GetYaxis()->SetLabelSize(0.04);
    _stack->GetYaxis()->SetLabelOffset(0.002);
    _stack->GetYaxis()->SetTitle("# Events");
    _stack->GetYaxis()->SetTitleSize(0.04);
    _stack->GetYaxis()->SetTitleOffset(1.32);
    _stack->SetMinimum(8);
}

std::unique_ptr<TH1> mc_comparison_entry::get(const std::string &name, double lumi)
{
    if (_stack == nullptr) {
        create_stack(name, lumi);
        if (_stack == nullptr) {
            return nullptr;
        }
    }
    if (_stack->GetNhists() == 0) {
        return nullptr;
    }
    TObjArray *partial_sums = _stack->GetStack(); // Not owned
    std::unique_ptr<TH1> res(dynamic_cast<TH1 *>(partial_sums->Last()->Clone()));
    return res;
}

TAxis *mc_comparison_entry::get_x_axis(const std::string &name, double lumi)
{
    if (_stack == nullptr) {
        create_stack(name, lumi);
        if (_stack == nullptr) {
            return nullptr;
        }
    }
    if (_stack->GetNhists() == 0) {
        return nullptr;
    }
    return _stack->GetXaxis();
}

void mc_comparison_entry::reset_drawing_state()
{
    _stack.reset();
}

void mc_comparison_entry::create_stack(const std::string &name, double lumi)
{
    _stack = std::make_unique<THStack>("stack", "");
    _legend.clear();
    bool had_histo = false;
    for (auto it = _groups.rbegin(); it != _groups.rend(); ++it) {
        // Get the histogram
        std::unique_ptr<TH1> histo = it->get(name);

        // Add it to the stack
        if (histo != nullptr) {
            histo->Scale(lumi);
            _stack->Add(dynamic_cast<TH1 *>(histo.get()->Clone()));
            _legend.push_back(it->legend());
            had_histo = true;
        }
    }
    if (!had_histo) {
        _stack.reset(); // Prevents crashes
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

data_comparison_entry::data_comparison_entry(const std::string &analyzer_name,
                                             const sample &sample,
                                             const std::string &input_dir)
    : _file(sample.histogram_file(analyzer_name, input_dir)),
      _histo(nullptr)
{
    if (_file == nullptr) {
        throw std::runtime_error("Could not open file for sample " + sample.name());
    }

    // Read job info histograms
    TH1 *job_info = nullptr;
    _file->GetObject("_job_info", job_info);
    if (job_info == nullptr) {
        throw std::runtime_error("File " + std::string(_file->GetName()) +
                                 " doesn't have the _job_info histogram.");
    }

    TVectorD *job_info_average = nullptr;
    _file->GetObject("_job_info_average", job_info_average);
    if (job_info_average == nullptr) {
        throw std::runtime_error("File " + std::string(_file->GetName()) +
                                 " doesn't have the _job_info_average vector.");
    }

    _frac = job_info->GetBinContent(1);
    if (std::isnan(_frac) || _frac < 0 || _frac > 1) {
        util::logging::warn << "Inconsistent processed fraction of sample: "
                            << _frac << ". Setting it to 1." << std::endl;
        _frac = 1;
    }
    _wsum = job_info->GetBinContent(2);
    _lumi = (*job_info_average)[0];
    _xsec = sample.xsec() > 0 ? sample.xsec() : (*job_info_average)[1];
}

void data_comparison_entry::add_histograms(std::set<std::string> &histos)
{
    // Iterate on keys
    for (const auto &&obj : *_file->GetListOfKeys()) {
        if (auto key = dynamic_cast<const TKey *>(obj)) {
            histos.insert(key->GetName());
        }
    }
}

void data_comparison_entry::add_to_legend(TLegend &legend, const std::string &plotname, double lumi)
{
    if (_histo == nullptr) {
        create_histo(plotname, lumi);
        if (_histo == nullptr) {
            return;
        }
    }
    legend.AddEntry(_histo.get(), "Data");
}

void data_comparison_entry::draw(const std::string &name, double lumi, bool same)
{
    if (_histo == nullptr) {
        create_histo(name, lumi);
        if (_histo == nullptr) {
            return;
        }
    }
    _histo->Draw(same ? "e same" : "e");
    _histo->SetMarkerStyle(20);
    _histo->SetMarkerColor(kBlack);
    _histo->SetLineColor(kBlack);
}

std::unique_ptr<TH1> data_comparison_entry::get(const std::string &name, double lumi)
{
    if (_histo == nullptr) {
        create_histo(name, lumi);
        if (_histo == nullptr) {
            return nullptr;
        }
    }
    std::unique_ptr<TH1> res(dynamic_cast<TH1 *>(_histo->Clone()));
    return res;
}

TAxis *data_comparison_entry::get_x_axis(const std::string &name, double lumi)
{
    if (_histo == nullptr) {
        create_histo(name, lumi);
        if (_histo == nullptr) {
            return nullptr;
        }
    }
    return _histo->GetXaxis();
}

void data_comparison_entry::reset_drawing_state()
{
    _histo = nullptr;
}

void data_comparison_entry::create_histo(const std::string &name, double lumi)
{
    TH1 *histo = nullptr;
    _file->cd();
    _file->GetObject(name.c_str(), histo);
    if (histo != nullptr) {
        _histo.reset(dynamic_cast<TH1 *>(histo->Clone()));
        if (_lumi > 0) { // Data
            _histo->Scale(lumi / _lumi / _frac);
        } else { // MC
            _histo->Scale(lumi * _xsec / _wsum);
        }
    } else {
        _histo.reset();
    }
}
} // namespace data
