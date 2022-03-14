#include <numeric>

#include <TChain.h>
#include <TH1F.h>
#include <TTreeReader.h>

#include "catalog.h"
#include "chains.h"
#include "logging.h"
#include "options.h"
#include "sample.h"
#include "timer.h"

void process_files(const data::catalog &catalog,
                   std::unique_ptr<TH1F> &njets,
                   std::unique_ptr<TH1F> &pt);

int main(int argc, char **argv)
{
    util::options opt;
    try {
        opt.default_init(argc, argv, "dyjets.yml", {});

        std::string bonzai_dir =
            "/store/group/phys_smp/AnalysisFramework/Bonzai/MC/v9.5/Catalogs/";

        std::vector<std::string> catalogs = {
            "Bonzais-DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-DMuUnf.txt",
            "Bonzais-DYToLL_0J_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-DMuUnf.txt",
            "Bonzais-DYToLL_1J_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-DMuUnf.txt",
            "Bonzais-DYToLL_2J_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-DMuUnf.txt",
        };

        std::unique_ptr<TH1F> njets, pt;
        for (const auto &path : catalogs) {
            process_files(data::catalog(path, bonzai_dir), njets, pt);
        }

        {
            // double inc_xs = 5941.0; // Not used
            std::array<double, 5> exc_xs{ 0, 4620.52, 859.59, 338.26, 0 };

            std::unique_ptr<TH1F> sf;
            sf.reset(dynamic_cast<TH1F *>(njets->Clone()));

            for (int bin = 1; bin <= njets->GetNbinsX(); ++bin) {
                sf->SetBinContent(bin, exc_xs[bin] / njets->GetBinContent(bin));
            }

            // Scale the factors to keep the sum of weights constant
            // This is important because reweighing.[h,cc] doesn't modify the
            // sum of weights
            sf->Scale(njets->Integral(0, njets->GetNbinsX() + 1) /
                      std::accumulate(exc_xs.begin(), exc_xs.end(), 0.));

            double sum;

            util::logging::info << "Njets table follows:" << std::endl;
            for (int bin = 1; bin <= njets->GetNbinsX(); ++bin) {
                if (sf->GetBinContent(bin) > 0)
                    sum += sf->GetBinContent(bin) * njets->GetBinContent(bin);
                // We use the first dimension only
                std::cout << "-1e10\t1e10\t"
                          << sf->GetBinLowEdge(bin) << "\t"
                          << sf->GetBinLowEdge(bin + 1) << "\t"
                          << sf->GetBinContent(bin) << "\t"
                          << 0 << "\t"
                          << 0 << std::endl;
            }
        }
    } catch (std::exception &e) {
        util::logging::fatal << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


void process_files(const data::catalog &catalog,
                   std::unique_ptr<TH1F> &njets,
                   std::unique_ptr<TH1F> &pt)
{
    auto chains = util::chains(catalog.files());
    auto headers = chains.bonzai_header();

    TTreeReader reader(headers.get());
    auto count = reader.GetEntries(true);
    if (count == 0) {
        throw std::runtime_error(
            "Input files don't appear to contain data. Is your proxy valid?");
    }

    auto njets_reader = TTreeReaderValue<TH1F>(reader, "hWeightednpNLO");
    auto pt_reader = TTreeReaderValue<TH1F>(reader, "hWeightedLHEZPt");

    for (long long entry = 0; entry < count; ++entry) {
        TTreeReader::EEntryStatus status = reader.SetEntry(entry);
        if (status != TTreeReader::kEntryValid) {
            throw std::runtime_error("Reader status code not valid: " +
                                        std::to_string(status));
            continue;
        }

        if (njets == nullptr) {
            njets.reset(dynamic_cast<TH1F *>(njets_reader->Clone()));
        } else {
            njets->Add(&*njets_reader);
        }

        if (pt == nullptr) {
            pt.reset(dynamic_cast<TH1F *>(pt_reader->Clone()));
        } else {
            pt->Add(&*pt_reader);
        }
    }
}
