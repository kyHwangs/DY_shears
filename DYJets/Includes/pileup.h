#ifndef PILEUP_H
#define PILEUP_H

#include <TTreeReaderValue.h>

#include "job.h"
#include "standalone_LumiReWeighting.h"

namespace util
{
class histo_set;
class options;
} // namespace util

namespace physics
{

class weights;

/// \brief Handles pileup
class pileup
{
    //TTreeReaderValue<int> Pileup_nTrueInt;
    TTreeReaderValue<int> PV_npvsGood;

    bool _reweighing_enabled = true;
    standalone_LumiReWeighting _standalone_lrw;

  public:
    /// \brief Constructor
    explicit pileup(util::job::info &info, const util::options &opt);

    /// \brief Returns the number of pileup vertices
    int nvtx() { return *PV_npvsGood; }

    /// \brief Reweights \c weights to take PU into account
    void reweight(weights &w);

    /// \brief Declares PU control histograms
    void declare_histograms(util::histo_set &h) const;

    /// \brief Fills PU control histograms
    void fill(util::histo_set &h, const std::string &tag, const weights &w);
};
} // namespace physics

#endif // PILEUP_H
