#ifndef JETS_H
#define JETS_H

#include <vector>

#include <TLorentzVector.h>
#include <TTreeReaderArray.h>

#include "JetResolution.h"
#include "JetResolutionObject.h"
#include "TRandom3.h"
#include "histo_set.h"
#include "job.h"
#include "lepton.h"
#include "options.h"
#include "weights.h"
namespace physics
{

/// \brief Represents a jet.
class jet
{
  public:
    TLorentzVector raw_v; ///< Unsmeared four-momentum
    TLorentzVector v; ///< Four-momentum
    float id;         ///< Jet ID
    float puMva;      ///< Result of the pileup MVA
    float bdisc;      ///< b-tag ID score
    float hadflav;    ///< jet hadron flavor
};

/// \brief Handles jets.
class jets
{
    TTreeReaderArray<float> Jet_pt;
    TTreeReaderArray<float> Jet_eta;
    TTreeReaderArray<float> Jet_phi;
    TTreeReaderArray<float> Jet_mass;
    TTreeReaderArray<int> Jet_jetId;
    TTreeReaderArray<float> Jet_puIdDisc;
    TTreeReaderArray<float> Jet_btagCSVV2;
    //TTreeReaderArray<float> Jet_hadronFlavour;
    TTreeReaderValue<float> fixedGridRhoFastjetAll;

    std::optional<TTreeReaderArray<float>> GenJet_pt;
    std::optional<TTreeReaderArray<float>> GenJet_eta;
    std::optional<TTreeReaderArray<float>> GenJet_phi;
    std::optional<TTreeReaderArray<float>> GenJet_mass;

    JME::JetResolution *m_JetResolution;
    JME::JetResolutionScaleFactor *m_JetResolutionScaleFactor;
    JME::JetParameters *m_JetParameters;
    float jetSF;
    double jetResolution;
    Variation m_Variation = Variation::NOMINAL;

    double _pt_cut = 30;
    double _y_cut = 2.4;
    double _pumva_cut = -0.2;
    double _deltar_cut = 0.4;
    bool _jer_smearing = true;

  public:
    /// \brief Constructor.
    explicit jets(util::job::info &info, const util::options &opt);

    /// \brief Configures the analyzer from user input.
    void configure(const util::options &opt);

    /**
     * \brief Retrieves a list of all jets in the current event.
     *
     * The list is already filtered according to config file options.
     */
    std::vector<jet> get(bool isdata, const std::vector<lepton> &leptons, double ptcut = -1);

    std::vector<jet> getGen();
    std::vector<jet> getGen(double ptmin, double rapmax);

    /// \brief Vetoes \c jets too close to one of the given \c leptons.
    void veto(std::vector<jet> &jets, const std::vector<lepton> &leptons) const;

    /// \brief Declares histograms filled by this class.
    void declare_histograms(util::histo_set &h);

    /**
     * \brief Fills muon control plots.
     * \param jets The list of jets in the event.
     * \param tag   A tag to pass to \ref util::histo_set
     */
    void fill(util::histo_set &h,
              const std::string &tag,
              const std::vector<jet> &jets,
              const weights &w);

    /// \brief Writes histograms to the current directory.
    void write();
};
} // namespace physics

#endif // JETS_H
