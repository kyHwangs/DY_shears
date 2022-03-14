#ifndef ZFINDER_H
#define ZFINDER_H

#include <string>

#include "lepton.h"
#include "options.h"

namespace physics
{

/**
 * \brief Represents a pair of leptons, interpreted as originating from a single boson.
 */
class dilepton
{
  public:
    TLorentzVector v;     ///< \brief The four-momentum of the dilepton
    float charge_product; ///< \brief The product of the charges of \c a and \c b
    lepton a, b;          ///< \brief The leptons used to create the dilepton

    /**
     * \brief Constructs a dilepton with two default-constructed leptons.
     */
    explicit dilepton() = default;

    /**
     * \brief Constructs a dilepton from two leptons.
     */
    explicit dilepton(const lepton &a, const lepton &b);

    /**
     * \brief Returns the distance of this dilepton to the Z mass
     */
    double distance_to_z() const { return std::abs(v.M() - 91.1876); }

    /**
     * \brief Returns the \f$ \phi^* \f$ variable
     */
    double phistar() const;

    /**
     * \brief Returns \c true if \c a is closer to the Z mass than \c b.
     *
     * This function is intended to be used in functions requiring a \c Compare predicate (e.g.
     * \c std::sort).
     */
    static bool zmass_ordering(const dilepton &a, const dilepton &b)
    {
        return a.distance_to_z() < b.distance_to_z();
    }
};

/**
 * \brief Finds Z bosons (dileptons) in a list of leptons.
 */
class zfinder
{
  public:
    enum class charge_mode { none, same_sign, neutral };

    /// \brief Modes for the flavor check
    enum class flavor_mode
    {
        none, ///< No check
        same, ///< Same flavor
        ee,   ///< Two electrons
        mumu, ///< Two muons
        emu   ///< Electron and muon
    };

  private:
    double _mass_low = 0;
    double _mass_high = std::numeric_limits<double>::max();
    charge_mode _charge_mode = charge_mode::neutral;
    flavor_mode _flavor_mode = flavor_mode::same;

  public:
    /**
     * \brief Constructs a Z finder.
     *
     * The finder's options will be taken from the configuration file, under the section given by
     * \c name (if it is found).
     */
    explicit zfinder(const util::options &opt, const std::string &name);

    /// \brief Destructor.
    virtual ~zfinder() = default;

    /**
     * \brief Finds dileptons pairs in the list of input leptons
     *
     * Criteria from the configuration are applied. The output can be empty.
     */
    std::vector<dilepton> find(const std::vector<lepton> &inputs) const;

    /**
     * \brief Checks if a dilepton candidate meets this finder's criteria.
     */
    bool valid(const dilepton &candidate) const;

    /**
     * \brief Returns the flavor mode for this finder.
     */
    flavor_mode get_flavor_mode() const { return _flavor_mode; }
};

} // namespace physics

#endif // ZFINDER_H
