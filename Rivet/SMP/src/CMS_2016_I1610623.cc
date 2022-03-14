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
#include <vector>
#include <iostream>
//// Credits: code was copied from CMS_2015_I1410737.cc and adapted by Ph. Gras to SMP-16-005

#define QUOTE(arg) #arg
#define QUOTE2(arg) QUOTE(arg)

namespace Rivet {

  class CMS_2016_I1610623: public Analysis {
  public:
    int ievt = 0;
    enum histIds { kNjets_exc, kNjets_inc,
		   kPtj1, kPtj2, kPtj3, kPtj4,
		   kYj1, kYj2, kYj3, kYj4,
		   kHt1, kHt2, kHt3, kHt4,
		   kDphi1, kDphi2, kDphi3, kDphi4,
		   kDr, nHists};

    /// Constructor
    CMS_2016_I1610623()
      : Analysis("CMS_2016_I1610623"), ievt(0), _h(nHists)
    {  }


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs; ///< @todo No cuts?
      VisibleFinalState visfs(fs);

      VetoedFinalState fs_notaudeday(fs);
      fs_notaudeday.addDecayProductsVeto(PID::TAU);
      fs_notaudeday.addDecayProductsVeto(-PID::TAU);
      
      WFinder wFinder(fs_notaudeday, /*lepton cuts*/ Cuts::abseta < 2.4 && Cuts::pT > 25*GeV, PID::MUON, 
		      /*mT window*/ 50*GeV, sqrtS(),
		      /*missET cut*/ 0, 
		      /*dR dressing*/0.1, /*lepton dressing mode*/ WFinder::CLUSTERALL,
		      /*photon tracking*/ WFinder::NOTRACK, 
		      /*mass cut type*/ WFinder::TRANSMASS);

      addProjection(wFinder, "wmunuFinder");

      //      VetoedFinalState jetConstits(visfs);
      //      jetConstits.addVetoOnThisFinalState(wFinder);

      //      FastJets jets(jetConstits, FastJets::ANTIKT, 0.4);
      FastJets jets(visfs, FastJets::ANTIKT, 0.4);
      addProjection(jets, "AntiKt04Jets");

      int i = 1;
      for(auto& h: _h) h = bookHisto1D(i++, 1, 1);
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
      if(dressedLeptons.size() != 1){
	throw Exception("bug found  in file " __FILE__ ", line " QUOTE(__LINE__));
      }

      const Particle mu = dressedLeptons[0];

      // Cluster jets
      // NB. Veto has already been applied on leptons and photons used for dressing
      const FastJets& fj = applyProjection<FastJets>(event, "AntiKt04Jets");
      const Jets& jets = fj.jetsByPt(Cuts::absrap < 2.4 && Cuts::pT > 30*GeV);

      // Perform lepton-jet overlap rejection
      Jets goodjets;
      foreach (const Jet& j, jets) {
        // Decide if this jet is "good", i.e. isolated from the leptons
        /// @todo Nice use-case for any() and a C++11 lambda
        //bool overlap = false;
        //foreach (const Particle& l, dressedLeptons) {
        //  if (Rivet::deltaR(j, l) < 0.4) {
        //    overlap = true;
        //    break;
        //  }
        //}
	if(deltaR(j, mu) < 0.4) continue;

        // Fill good-jets collection
        goodjets.push_back(j);
      }


      // Weight to be used for histo filling
      const double w = event.weight();

      // Fill jet number integral histograms
      _h[kNjets_exc]->fill(goodjets.size(), w);
      /// @todo Could be better computed by toIntegral transform on exclusive histo
      for (size_t iJet = 0; iJet <= goodjets.size(); iJet++ )
        _h[kNjets_inc]->fill(iJet, w);


      if (goodjets.size() < 1) return;
      
      // Compute HT and min(\DeltaR(mu,j))
      double ht = 0;
      double drmin = 99.;
      for(auto j: goodjets){
	ht += j.pT();
	double dr = deltaR(mu, j);
	if(j.pT() > 100*GeV && dr < drmin) drmin = dr;
      }

      // Fill jet histograms
      for(unsigned ij = 0; ij < 4 && ij < goodjets.size(); ++ij){
	const auto& j = goodjets[ij];
	_h[kPtj1  + ij]->fill(j.pT() / GeV, w);
	_h[kYj1   + ij]->fill(j.absrapidity(), w);
	_h[kHt1   + ij]->fill(ht, w);
	double dphi = deltaPhi(mu, j);
	_h[kDphi1 + ij]->fill(dphi, w);
      }

      if(goodjets.size() > 0 && goodjets[0].pT() > 300*GeV) _h[kDr]->fill(drmin, w);
    }
    

    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = (sumOfWeights() != 0) ? crossSection()/sumOfWeights() : 1.0;

      MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
      MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
      MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());
      MSG_INFO("Norm factor   = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(6) << norm);

      MSG_INFO("My event count = " << std::setfill(' ') << std::setw(14) << std::fixed << ievt);

      for(auto h: _h) scale(h, norm);
    }


  private:

    /// @name Histograms
    
    std::vector<Histo1DPtr> _h;
  };


  DECLARE_RIVET_PLUGIN(CMS_2016_I1610623);


}
