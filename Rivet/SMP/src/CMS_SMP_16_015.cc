// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

#include <cstdlib>
#include <cstring>

#include <iostream>

//// Credits: code was copied from CMS_2015_I1410737.cc and adapted by Ph. Gras to SMP-15-010

namespace Rivet {

  class CMS_SMP_16_015 : public Analysis {
  public:

    /// Constructor
    CMS_SMP_16_015()
      : Analysis("CMS_SMP_16_015")
    {  }


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs; ///< @todo No cuts?
      VisibleFinalState visfs(fs);

      VetoedFinalState fs_notaudeday(fs);
      fs_notaudeday.addDecayProductsVeto(PID::TAU);
      fs_notaudeday.addDecayProductsVeto(-PID::TAU);
      
      ZFinder zeeFinder(fs_notaudeday, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::ELECTRON, 71.0*GeV, 111.0*GeV);
      addProjection(zeeFinder, "ZeeFinder");

      ZFinder zmumuFinder(fs_notaudeday, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::MUON, 71.0*GeV, 111.0*GeV);
      addProjection(zmumuFinder, "ZmumuFinder");


      FastJets jets(visfs, FastJets::ANTIKT, 0.4);
      addProjection(jets, "jets");


      _h_excmult_jets_tot   = bookHisto1D(1, 1, 1);
      _h_incmult_jets_tot   = bookHisto1D(2, 1, 1);
      _h_leading_jet_pt_tot = bookHisto1D(3, 1, 1);
      _h_second_jet_pt_tot  = bookHisto1D(4, 1, 1);
      _h_third_jet_pt_tot   = bookHisto1D(5, 1, 1);
      _h_leading_jet_y_tot  = bookHisto1D(6, 1, 1);
      _h_second_jet_y_tot   = bookHisto1D(7, 1, 1);
      _h_third_jet_y_tot    = bookHisto1D(8, 1, 1);
      _h_ht1_tot            = bookHisto1D(9, 1, 1);
      _h_ht2_tot 	    = bookHisto1D(10, 1, 1);
      _h_ht3_tot 	    = bookHisto1D(11, 1, 1);
      _h_ptbal1 	    = bookHisto1D(12, 1, 1);
      _h_ptbal2 	    = bookHisto1D(13, 1, 1);
      _h_ptbal3 	    = bookHisto1D(14, 1, 1);
      _h_jzb 		    = bookHisto1D(15, 1, 1);
      _h_jzb_ptHigh 	    = bookHisto1D(16, 1, 1);
      _h_jzb_ptLow  	    = bookHisto1D(17, 1, 1);
      _h_zpt0  		    = bookHisto1D(18, 1, 1);
      _h_zpt1  		    = bookHisto1D(19, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {;

      const ZFinder& zeeFS = applyProjection<ZFinder>(event, "ZeeFinder");
      const ZFinder& zmumuFS = applyProjection<ZFinder>(event, "ZmumuFinder");

      const Particles& zees = zeeFS.bosons();
      const Particles& zmumus = zmumuFS.bosons();


      // We did not find exactly one Z. No good.
      if (zees.size() + zmumus.size() != 1) {
	MSG_DEBUG("Did not find exactly one good Z candidate");
	vetoEvent;
      }

      
      // Find the (dressed!) leptons
      const ZFinder& z = zees.size() ? zeeFS : zmumuFS;
      const Particles& dressedLeptons =  z.constituents();

      // Cluster jets
      // NB. Veto has already been applied on leptons and photons used for dressing
      const FastJets& fj = applyProjection<FastJets>(event, "jets");
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


      // Fill inclusize histogram
      _h_zpt0->fill(z.boson().pT(), w);

      if (goodjets.size() < 1) return;

      //hadronic recoil:
      FourMomentum recoil;
      foreach (const Jet& j, goodjets) {
	recoil += j.momentum();
      }
      //Jet-Z balance = |recoil_T| - |pt(Z)|
      double jzb  = recoil.pT() - z.boson().pT();
      //Pt balance:
      double ptbal = (recoil + z.boson().momentum()).pT();

      // Fill leading jet histograms
      _h_zpt1->fill(z.boson().pT(), w);
      const Jet& j1 = goodjets[0];
      _h_leading_jet_pt_tot->fill(j1.pT()/GeV, w);
      _h_leading_jet_y_tot->fill(j1.absrapidity(), w);
      _h_ht1_tot->fill(ht/GeV, w);
      _h_jzb->fill(jzb/GeV, w);
      if(z.boson().pT() > 50*GeV){
	_h_jzb_ptHigh->fill(jzb/GeV, w);
      } else{
	_h_jzb_ptLow->fill(jzb/GeV, w);
      }
      _h_ptbal1->fill(ptbal/GeV, w);
      
      // Fill 2nd jet histograms
      if (goodjets.size() < 2) return;
      const Jet& j2 = goodjets[1];
      _h_second_jet_pt_tot->fill(j2.pT()/GeV, w);
      _h_second_jet_y_tot->fill(j2.absrapidity(), w);
      _h_ht2_tot->fill(ht/GeV, w);
      _h_ptbal2->fill(ptbal/GeV, w);

      // Fill 3rd jet histograms
      if (goodjets.size() < 3) return;
      const Jet& j3 = goodjets[2];
      _h_third_jet_pt_tot->fill(j3.pT()/GeV, w);
      _h_third_jet_y_tot->fill(j3.absrapidity(), w);
      _h_ht3_tot->fill(ht/GeV, w);
      _h_ptbal3->fill(ptbal/GeV, w);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = (sumOfWeights() != 0) ? crossSection()/sumOfWeights() : 1.0;

      norm /= 2.; //to allow use of same MC cross section than CMS_I1409817, which keeps Z->\mu\mu only

      MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
      MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
      MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());
      MSG_INFO("Norm factor   = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(6) << norm);

      scale(_h_excmult_jets_tot, norm);
      scale(_h_incmult_jets_tot, norm);
      scale(_h_leading_jet_pt_tot, norm);
      scale(_h_second_jet_pt_tot, norm);
      scale(_h_third_jet_pt_tot, norm);
      scale(_h_leading_jet_y_tot, norm);
      scale(_h_second_jet_y_tot, norm);
      scale(_h_third_jet_y_tot, norm);
      scale(_h_ht1_tot, norm);
      scale(_h_ht2_tot, norm);
      scale(_h_ht3_tot, norm);
      scale(_h_ptbal1, norm);
      scale(_h_ptbal2, norm);
      scale(_h_ptbal3, norm);
      scale(_h_jzb, norm);
      scale(_h_jzb_ptHigh, norm);
      scale(_h_jzb_ptLow, norm);
      scale(_h_zpt0, norm);
      scale(_h_zpt1, norm);
    }


  private:

    /// @name Histograms

    Histo1DPtr _h_excmult_jets_tot,  _h_incmult_jets_tot;
    Histo1DPtr _h_leading_jet_pt_tot, _h_second_jet_pt_tot, _h_third_jet_pt_tot, _h_fourth_jet_pt_tot;
    Histo1DPtr _h_leading_jet_y_tot, _h_second_jet_y_tot, _h_third_jet_y_tot, _h_fourth_jet_y_tot;
    Histo1DPtr _h_ht1_tot, _h_ht2_tot, _h_ht3_tot, _h_ht4_tot;
    Histo1DPtr _h_ptbal1, _h_ptbal2, _h_ptbal3;
    Histo1DPtr _h_jzb, _h_jzb_ptHigh, _h_jzb_ptLow;
    Histo1DPtr _h_zpt0, _h_zpt1;
  };

  DECLARE_RIVET_PLUGIN(CMS_SMP_16_015);
}
