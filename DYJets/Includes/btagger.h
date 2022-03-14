#ifndef BTAGGER_H
#define BTAGGER_H

#include "BTagCalibrationStandalone.h"
#include "jets.h"
#include "job.h"
#include "options.h"
#include "tables.h"
#include "weights.h"
#include "histo_set2D.h"

namespace physics
{

/// \brief Wrapper around b-tagging utilities
class btagger
{
    double _bjet_cut = 0.5426;
    std::string bjet_cut; 
    BTagCalibrationReader _btag_calibration_reader;

  public:
    /// \brief Constructor
    explicit btagger(const util::options &opt, util::histo_set2D &h);

    /// \brief Checks whether any jet is a b jet and apply the b-tagging scale factors
    bool any(const std::vector<jet> &jets, weights &w,util::histo_set2D &h, const util::tables &tab) const;

  private:
    /// \brief Applies the b-tagging scale factors for the given jet
    void apply_sf(const jet &j, weights &w,util::histo_set2D &h, const util::tables &tab, double &wu) const;
    void fill_eff(const jet &j, weights &w, util::histo_set2D &h, double wu) const;

};

} // namespace physics

#endif // BTAGGER_H
