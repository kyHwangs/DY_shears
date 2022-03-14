#ifndef REWEIGHING_H
#define REWEIGHING_H

#include <memory>

#include <TTreeReaderValue.h>

#include "tables.h"
#include "weights.h"

namespace physics {

/**
 * \brief Supports reweighing events to merge inclusive and exclusive samples.
 *
 * The following binned samples are supported:
 *
 *   * `npNLO` (number of partons)
 */
class reweighing
{
    mutable std::unique_ptr<TTreeReaderValue<int>> npNLO;
    util::table _table_npNLO;

  public:
    /// \brief Constructor.
    explicit reweighing(util::job::info &info, const util::options &opt);

    /// \brief Reweighs the current event.
    void reweigh(weights &w) const;
};

} // namespace physics

#endif // REWEIGHING_H
