#include "pileup.h"

#include "histo_set.h"
#include "options.h"
#include "weights.h"

namespace physics
{

pileup::pileup(util::job::info &info, const util::options &opt)
    : /*Pileup_nTrueInt(info.reader, "Pileup_nTrueInt"),*/
      PV_npvsGood(info.reader, "PV_npvsGood"),
      _standalone_lrw(opt.config["year"].as<int>(), 0)
{
    util::set_value_safe(
        opt.config, _reweighing_enabled, "use pileup reweighing", "pileup reweighing toggle");
}

void pileup::reweight(weights &w)
{
    if (_reweighing_enabled && w.ismc()) {
        //w.use_weight(_standalone_lrw.weight(*Pileup_nTrueInt));
        //w.use_gen_weight(_standalone_lrw.weight(*Pileup_nTrueInt));
    }
}

void pileup::declare_histograms(util::histo_set &h) const
{
    h.declare("nvtx", "Number of vertices;#Vtx", 60, 0.5, 60.5);
}

void pileup::fill(util::histo_set &h, const std::string &tag, const weights &w)
{
    h.fill("nvtx", tag, *PV_npvsGood, w.global_weight());
}
} // namespace physics
