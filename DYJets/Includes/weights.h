#ifndef WEIGHTS_H
#define WEIGHTS_H

#include "chains.h"
#include "histo_set.h"
#include "job.h"

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <optional>

namespace physics
{

/**
 * \brief Analyzer for event weights.
 *
 * This class is responsible for the recording of processed event weights, lumiosity and cross
 * section. It does *not* mean that all histograms will be weighted automatically.
 */
class weights
{
    std::optional<TTreeReaderArray<float>> genWeight;

    bool _ismc;

    long long _processed_events = 0;
    long long _events_in_chain = 0;

    double _processed_weights_sum = 0;

    double _xsec;
    double _lumi;

    double _gen_weight;
    double _global_weight;

  public:
    /// \brief Constructor.
    explicit weights(util::job::info &info);

    /// \brief Call this for every processed event.
    void process_event();

    /// \brief Writes results to the current (ROOT) directory.
    void write(util::histo_set *histos);

    /// \brief Returns the size of the current weight vector.
    std::size_t weights_count() { return genWeight ? genWeight->GetSize() : 0; }

    /// \brief Returns the contents of the current weight vector at index \c i (checked).
    double weight_at(std::size_t i) { return genWeight->At(i); }

    /// \brief Returns the gen-level weight of the event.
    double gen_weight() const { return _gen_weight; }

    /// \brief Returns the global weight of the event.
    double global_weight() const { return _global_weight; }

    /// \brief Edits the global weight of the event.
    void use_weight(double weight) { _global_weight *= weight; }

    /// \brief Edits the global weight of the event.
    void use_gen_weight(double weight) { _gen_weight *= weight; }

    /// \brief Checks whether the current event is from MC.
    bool ismc() const { return _ismc; }

    /// \brief Checks whether the current event is real data.
    bool isdata() const { return !_ismc; }
};
}

#endif // WEIGHTS_H
