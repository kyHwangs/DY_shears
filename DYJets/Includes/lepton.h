#ifndef LEPTON_H
#define LEPTON_H

#include <TLorentzVector.h>

#include "cmake_config.h"

/// \brief Namespace for physics objects.
namespace physics
{

/// \brief Represents a lepton.
class lepton
{
  public:
    TLorentzVector v; ///< Four-momentum
    TLorentzVector raw_v; ///< Four-momentum as provided by CMSSW (ie before applying corrections)
    float charge;     ///< Charge
    unsigned iso;        ///< Isolation ID
    bool passes_id;   ///< Does the lepton pass the Id cut?
    bool passes_iso;   ///< Does the lepton pass the Iso cut?
    unsigned id;      ///< Id
    int pdgid;        ///< PDG ID (absolute value): electron = 11, muon = 13

#ifdef DEBUG_PRINTOUT
    int tkLayerCnt;
    int fnUsed;
    int gllPt;
#endif // DEBUG_PRINTOUT

    bool operator== (const lepton &other) const
    {
        return v == other.v
            && raw_v == other.raw_v
            && charge == other.charge
            && iso == other.iso
            && passes_id == other.passes_id
            && passes_iso == other.passes_iso
            && id == other.id
            && pdgid == other.pdgid;
    }
};

} // namespace physics

#endif // LEPTON_H
