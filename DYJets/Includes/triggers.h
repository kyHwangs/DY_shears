#ifndef TRIGGERS_H
#define TRIGGERS_H

#include <memory>
#include <string>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "options.h"
#include "job.h"

class TTree;

namespace physics
{

/**
 * \brief Holds the list of triggers used in an analysis.
 *
 * This class can be used when one needs to accept events that pass any of a
 * given set of triggers, and reject others. To do this, one first passes all
 * accepted triggers to \ref accept, and then calls \ref passes for every event.
 */
class trigger_mask
{
    bool _accepts_any_trigger;
    std::vector<TTreeReaderValue<bool>> _accepted_triggers;
    std::vector<TTreeReaderValue<bool>> _vetoed_triggers;
  public:
    /**
     * \brief Constructor.
     * \throws std::invalid_argument if \c bitFieldsChain is empty.
     */
    explicit trigger_mask(util::job::info &info);

    /**
     * \brief Constructs a trigger mask from user input.
     *
     * The definition provided must be a comma- or space-separated list of (possibly quoted)
     * trigger names. A trigger can be specified as a veto by putting a \c ^ before its name.
     *
     * For sample, the definition `HLT_IsoMu24, HLT_IsoTkMu24, ^HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL`
     * will produce a mask accepting events passing either of the `HLT_IsoMu24` and `HLT_IsoTkMu24`
     * triggers, and rejecting events passing the `HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL`.
     *
     * If the definition doesn't have any accepted trigger, anything will be accepted (as per
     * \ref set_accepts_any_trigger).
     *
     * If \c verbose is enabled, the constructor will log some information to \c cout.
     *
     * \throws std::invalid_argument if \c bitFieldsChain is empty, or a trigger is not found.
     */
    explicit trigger_mask(util::job::info &info,
                          const std::string &definition,
                          bool verbose = false);

    /**
     * \brief Sets the given trigger path to be accepted by \ref passes.
     * \returns \c true if the path was found.
     */
    bool accept(TTreeReader &tr, const std::string &name);

    /**
     * \brief Sets all of the given trigger paths to be accepted by \ref passes.
     *
     * If some paths aren't found, this function will set all others to be
     * accepted and return \c false.
     *
     * \returns \c true if all paths were found.
     */
    bool accept(TTreeReader &tr, const std::vector<std::string> &name);

    /// \brief Sets whether all triggers should be accepted.
    void set_accepts_any_trigger(bool enable) { _accepts_any_trigger = enable; }

    /// \brief Returns whether all triggers are accepted.
    bool accepts_any_trigger() { return _accepts_any_trigger; }

    /**
     * \brief Sets the given trigger path to cause \ref passes to return \c false.
     * \returns \c true if the path was found.
     */
    bool veto(TTreeReader &tr, const std::string &name);

    /**
     * \brief Sets all of the given trigger paths cause \ref passes to return \c false.
     *
     * If some paths aren't found, this function will set all others to be
     * vetoed and return \c false.
     *
     * \returns \c true if all paths were found.
     */
    bool veto(TTreeReader &tr, const std::vector<std::string> &name);

    /**
     * \brief Checks whether the event passes the trigger requirements.
     *
     * An event passes the trigger requirement if no veto'ed trigger fired, and at least one
     * accepted trigger fired (or \c accepts_any_trigger was set).
     */
    bool passes();
};
} // namespace physics

#endif // TRIGGERS_H
