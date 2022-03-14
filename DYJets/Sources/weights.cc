#include "weights.h"

#include <TVectorD.h>

namespace physics
{

weights::weights(util::job::info &info)
    : genWeight(info.init_optional_branch<decltype(genWeight)>("genWeight")),
//I: from tuppel: EvtWeights_->push_back(genEventInfoProd->weight()); EvtWeights_->push_back(lheEvent->weights()[iw].wgt)
//available in nano: genWeight, LHEPdfWeight, LHEReweightingWeight, LHEScaleWeight, LHEWeight, PSWeight (w_var / w_nominal), )
      _ismc(info.catalog.ismc()),
      
      _events_in_chain(info.reader.GetEntries(true)),
      _xsec(info.catalog.xsec()),
      _lumi(info.catalog.lumi())
{
}
 
void weights::process_event()
{
    if (weights_count() > 0) {
        /*if (!_is_weights_sum_divided) {
            // This needs to be done in case several samples are merged (think
            // exclusive 0, 1, 2 jets). Their weights will in general not be on
            // the same scale.

            // We assume that all abs(weight) are equal. Additional support in
            // the pruner is needed if this assumption is not correct.

            _weights_sum /= std::abs(weight_at(0));
            _is_weights_sum_divided = true;
        }*/

// normalizing the sum weights. will need to revisit once we have exclusive samples.
        _gen_weight = weight_at(0)/ std::abs(weight_at(0));
        _global_weight = weight_at(0)/ std::abs(weight_at(0));
        _processed_weights_sum += weight_at(0)/std::abs(weight_at(0));
    } else {
        _gen_weight = 1;
        _global_weight = 1;
    }
    _processed_events++;
}

void weights::write(util::histo_set *histos)
{
    double fraction_processed = 1;

    // Declare and fill summed info histogram
    histos->declare("_job_info", "Job information", 4, 0, 4);
    util::histo_set::histogram_type &job_info = histos->get("_job_info");

    job_info.GetXaxis()->SetBinLabel(1, "fraction_processed"); // For data
    job_info.SetBinContent(1, fraction_processed);

    job_info.GetXaxis()->SetBinLabel(2, "weights_sum"); // For MC
    job_info.SetBinContent(2, _processed_weights_sum);

    // Declare and fill info vector
    // We use TVectorD because hadd won't sum them, as it should be for sample lumi and xsec
    // (TH1::kSetAverage is broken in hadd)
    TVectorD job_info_average(2);
    job_info_average[0] = _lumi; // For data
    job_info_average[1] = _xsec; // For MC
    job_info_average.Write("_job_info_average");

    if (isdata()) {
        util::logging::info << "Processed fraction of sample: " << fraction_processed << std::endl;
    }
}
} // namespace physics
