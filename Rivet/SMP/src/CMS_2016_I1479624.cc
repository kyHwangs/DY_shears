// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"

#include <cstdlib>
#include <cstring>

#include <iostream>

//// Credits: code was copied from CMS_2015_I1410737.cc and adapted by Ph. Gras to SMP-16-005

namespace Rivet {

  class CMS_2016_I1479624: public Analysis {
  public:
  int ievt = 0;
  

    /// Constructor
    CMS_2016_I1479624()
      : Analysis("CMS_2016_I1479624"), ievt(0)
    {  }


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs; ///< @todo No cuts?
      VisibleFinalState visfs(fs);
      
      WFinder wFinder(fs, /*lepton cuts*/ Cuts::abseta < 2.4 && Cuts::pT > 25*GeV, PID::MUON, 
		      /*mT window*/ 50*GeV, sqrtS(),
		      /*missET cut*/ 0, 
		      /*dR dressing*/0.1, /*lepton dressing mode*/ WFinder::CLUSTERALL,
		      /*photon tracking*/ WFinder::NOTRACK, 
		      /*mass cut type*/ WFinder::TRANSMASS);

      addProjection(wFinder, "wmunuFinder");

      VetoedFinalState jetConstits(visfs);
      jetConstits.addVetoOnThisFinalState(wFinder);

      FastJets jets(jetConstits, FastJets::ANTIKT, 0.4);
      addProjection(jets, "AntiKt04Jets");

      _h_excmult_jets_tot   = bookHisto1D(1, 1, 1);
      _h_incmult_jets_tot   = bookHisto1D(2, 1, 1);
      _h_leading_jet_pt_tot = bookHisto1D(3, 1, 1);
      _h_second_jet_pt_tot  = bookHisto1D(4, 1, 1);
      _h_third_jet_pt_tot   = bookHisto1D(5, 1, 1);
      _h_leading_jet_y_tot  = bookHisto1D(6, 1, 1);
      _h_second_jet_y_tot   = bookHisto1D(7, 1, 1);
      _h_third_jet_y_tot    = bookHisto1D(8, 1, 1);
      _h_ht1_tot 	    = bookHisto1D(9, 1, 1);
      _h_ht2_tot 	    = bookHisto1D(10, 1, 1);
      _h_ht3_tot 	    = bookHisto1D(11, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {;
      ++ievt;
      const WFinder& wFS = applyProjection<WFinder>(event, "wmunuFinder");

      //const Particles& zees = zeeFS.bosons();
      const Particles& wmunu = wFS.bosons();

      // We did not find exactly one Z. No good.
      //if (zees.size() + zmumus.size() != 1) {
      if (wmunu.size() != 1) {
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }

      // Find the (dressed!) leptons
      const Particles& dressedLeptons =  wFS.constituentLeptons();

      // Cluster jets
      // NB. Veto has already been applied on leptons and photons used for dressing
      const FastJets& fj = applyProjection<FastJets>(event, "AntiKt04Jets");
      const Jets& jets = fj.jetsByPt(Cuts::absrap < 2.4 && Cuts::pT > 30*GeV);

      // Perform lepton-jet overlap and HT calculation
      double ht = 0;
      Jets goodjets;
      foreach (const Jet& j, jets) {
        // Decide if this jet is "good", i.e. isolated from the leptons
        /// @todo Nice use-case for any() and a C++11 lambda
        bool overlap = false;
        foreach (const Particle& l, dressedLeptons) {
          if (Rivet::deltaR(j, l) < 0.4) {
            overlap = true;
            break;
          }
        }

        // Fill HT and good-jets collection
        if (overlap) continue;
        goodjets.push_back(j);
        ht += j.pT();
      }

      // We don't care about events with no isolated jets
      //      if (goodjets.empty()) {
      //        MSG_DEBUG("No jets in event");
      //  vetoEvent;
      //}


      /////////////////


      // Weight to be used for histo filling
      const double w = event.weight();

      // Fill jet number integral histograms
      _h_excmult_jets_tot->fill(goodjets.size(), w);
      /// @todo Could be better computed by toIntegral transform on exclusive histo
      for (size_t iJet = 0; iJet <= goodjets.size(); iJet++ )
        _h_incmult_jets_tot->fill(iJet, w);


      if (goodjets.size() < 1) return;
      
      // Fill leading jet histograms
      const Jet& j1 = goodjets[0];
      _h_leading_jet_pt_tot->fill(j1.pT()/GeV, w);
      _h_leading_jet_y_tot->fill(j1.absrapidity(), w);
      _h_ht1_tot->fill(ht/GeV, w);

      // Fill 2nd jet histograms
      if (goodjets.size() < 2) return;
      const Jet& j2 = goodjets[1];
      _h_second_jet_pt_tot->fill(j2.pT()/GeV, w);
      _h_second_jet_y_tot->fill(j2.absrapidity(), w);
      _h_ht2_tot->fill(ht/GeV, w);

      // Fill 3rd jet histograms
      if (goodjets.size() < 3) return;
      const Jet& j3 = goodjets[2];
      _h_third_jet_pt_tot->fill(j3.pT()/GeV, w);
      _h_third_jet_y_tot->fill(j3.absrapidity(), w);
      _h_ht3_tot->fill(ht/GeV, w);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = (sumOfWeights() != 0) ? crossSection()/sumOfWeights() : 1.0;

      MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
      MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
      MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());
      MSG_INFO("Norm factor   = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(6) << norm);

      MSG_INFO("My event count = " << std::setfill(' ') << std::setw(14) << std::fixed << ievt);
      
      scale(_h_excmult_jets_tot, norm );
      scale(_h_incmult_jets_tot, norm );
      scale(_h_leading_jet_pt_tot, norm );
      scale(_h_second_jet_pt_tot, norm );
      scale(_h_third_jet_pt_tot, norm );
      scale(_h_leading_jet_y_tot, norm );
      scale(_h_second_jet_y_tot, norm );
      scale(_h_third_jet_y_tot, norm );
      scale(_h_ht1_tot, norm );
      scale(_h_ht2_tot, norm );
      scale(_h_ht3_tot, norm );
    }


  private:

    /// @name Histograms

    Histo1DPtr _h_excmult_jets_tot,  _h_incmult_jets_tot;
    Histo1DPtr _h_leading_jet_pt_tot, _h_second_jet_pt_tot, _h_third_jet_pt_tot, _h_fourth_jet_pt_tot;
    Histo1DPtr _h_leading_jet_y_tot, _h_second_jet_y_tot, _h_third_jet_y_tot, _h_fourth_jet_y_tot;
    Histo1DPtr _h_ht1_tot, _h_ht2_tot, _h_ht3_tot, _h_ht4_tot;

  };


  DECLARE_RIVET_PLUGIN(CMS_2016_I1479624);


}
