#include "Pruner.h"
#include <iostream>
#include <vector>
#include "TLorentzVector.h"
#include <cstdlib>
using namespace std;

//Includes code generated by TTree:MakeClass():
//The fields of EventTree.h class will be used
//to access the EventTree contents. EventTree::Init(TTree*)
//should be called by the init(Tree*) method to link these
//fields to the EventTree contents.
#define EventTree_cxx
#include "EventTree.h"
void EventTree::Loop(){} //To make the linker happy.

/** Selection parameters **/
enum {SingleMu,SingleE, Ntrig}; //List of triggers used for our analysis
//The lines below correspond to the triggerBits used in our analysis. These are all defined in the trigger2016.h file of the ntuple producer.
static int trigSingleMu[2] = {11,16};
static int trigSingleE[1] = {9};
static float minLepPt = 10;
static const unsigned elIdMask = (1 <<10); //that is 1024;
class TagProbePruner: public Pruner, EventTree{

protected:
  void declareSubSelections();
  bool init(TChain* tree);
  bool filterEvent();
  bool filterEl(int iEl);
  bool filterMu(int iMu);
  bool isTightEle(std::vector<float> *ElPt, std::vector<float> *ElEta, std::vector<float> *ElPhi, std::vector<float> *ElE, std::vector<unsigned int> *ElId, std::vector<float> *ElEtaSc, std::vector<float> *ElPfIsoRho);
  bool isTightMuon(std::vector<float> *MuPt, std::vector<float> *MuEta, std::vector<float> *MuPhi, std::vector<float> *MuE, std::vector<unsigned int> *MuId, std::vector<unsigned int> *MuIdTight, std::vector<float> *MuPfIso);
  void makeFilterMask(bool (TagProbePruner::*filter)(int), std::vector<bool>& mask);
  void skimCollections();
  bool eventSelection();
  bool passTrigger(int trig);
  enum { Mu, Ele,NSubSels};
  static const int kEl = 11;
  static const int kMu = 13;
};

DECLARE_PRUNER(TagProbePruner, "Pruner of Tag & Probe analysis")

void TagProbePruner::declareSubSelections(){
  subSelections_.resize(NSubSels);
  subSelections_[Mu]		= SubSelection("Mu","muon selection for Tag & Probe");
  subSelections_[Ele]		= SubSelection("Ele","electron selection for Tag & Probe");	
  }

bool TagProbePruner::init(TChain* tree){
  if(tree == 0) return false;
  EventTree::Init(tree);
  return true;
  }

void TagProbePruner::makeFilterMask(bool (TagProbePruner::*filter)(int), std::vector<bool>& mask){
  for(unsigned i = 0; i < mask.size(); ++i){
    mask[i] = (this->*filter)(i);
    }
  }

bool TagProbePruner::filterEvent(){
  skimCollections();
  return eventSelection() && Pruner::filterEvent();
  }

void TagProbePruner::skimCollections(){
  std::vector<bool> mask;
  //Reco muon collections:
  mask.resize(MuPt->size());
  makeFilterMask(&TagProbePruner::filterMu, mask);
  filter(MuPt, mask);
  filter(MuEta, mask);
  filter(MuPhi, mask);
  filter(MuE, mask);
  filter(MuId, mask);
  filter(MuIdTight, mask);
  filter(MuCh, mask);
  filter(MuVtxZ, mask);
  filter(MuDxy, mask);
  filter(MuIsoRho, mask);
  filter(MuPfIso, mask);
  filter(MuType, mask);
  filter(MuIsoTkIsoAbs, mask);
  filter(MuIsoTkIsoRel, mask);
  filter(MuIsoCalAbs, mask);
  filter(MuIsoCombRel, mask);
  filter(MuTkNormChi2, mask);
  filter(MuTkHitCnt, mask);
  filter(MuMatchedStationCnt, mask);
  filter(MuDz, mask);
  filter(MuPixelHitCnt, mask);
  filter(MuTkLayerCnt, mask);
  filter(MuPfIsoChHad, mask);
  filter(MuPfIsoNeutralHad, mask);
  filter(MuPfIsoRawRel, mask);
  filter(MuHltMatch, mask);

  //Reco electron collections:
  mask.resize(ElPt->size());
  makeFilterMask(&TagProbePruner::filterEl, mask);
  filter(ElPt, mask);
  filter(ElEta, mask);
  filter(ElEtaSc, mask);
  filter(ElPhi, mask);
  filter(ElE, mask);
  filter(ElId, mask);
  filter(ElCh, mask);
  filter(ElMvaTrig, mask);
  filter(ElMvaNonTrig, mask);
  filter(ElMvaPresel, mask);
  filter(ElDEtaTkScAtVtx, mask);
  filter(ElDPhiTkScAtVtx, mask);
  filter(ElHoE, mask);
  filter(ElSigmaIetaIeta, mask);
  filter(ElSigmaIetaIetaFull5x5, mask);
  filter(ElEinvMinusPinv, mask);
  filter(ElD0, mask);
  filter(ElDz, mask);
  filter(ElExpectedMissingInnerHitCnt, mask);
  filter(ElPassConvVeto, mask);
  filter(ElHltMatch, mask);
  filter(ElPfIsoChHad, mask);
  filter(ElPfIsoNeutralHad, mask);
  filter(ElPfIsoIso, mask);
  filter(ElPfIsoPuChHad, mask);
  filter(ElPfIsoRaw, mask);
  filter(ElPfIsoDbeta, mask);
  filter(ElPfIsoRho, mask);
  filter(ElAEff, mask);
  }



bool TagProbePruner::isTightEle(std::vector<float> *ElPt, std::vector<float> *ElEta, std::vector<float> *ElPhi, std::vector<float> *ElE, std::vector<unsigned int> *ElId, std::vector<float> *ElEtaSc, std::vector<float> *ElPfIsoRho)
{
  bool Eallpass= false;
  for(unsigned int i = 0 ; i<ElPt->size() ; i++){
    bool passEta = false, passIso = false, passId = false, passPt = false;
    TLorentzVector currentLepton; currentLepton.SetPtEtaPhiE(ElPt->at(i),ElEta->at(i),ElPhi->at(i),ElE->at(i));
    passId = ElId->at(i) & (1<<17);
    int eta = fabs(ElEtaSc->at(i));//I took the supercluster eta since it's really the geometry which is taken here.
    passEta = (eta<=2.5 && (eta>=1.5660 || eta<=1.4442));
    if(eta>=1.5660 && ElPfIsoRho->at(i)<0.0646) passIso = true;
    if(eta<=1.4442 && ElPfIsoRho->at(i)<0.0354) passIso = true; //Numbers are taken from llvv_fwk and have not been checked.
    passPt = (currentLepton.Pt() >=30);
    Eallpass =Eallpass||(passEta && passIso && passId && passPt);
    }
  return Eallpass;
  }
bool TagProbePruner::isTightMuon(std::vector<float> *MuPt, std::vector<float> *MuEta, std::vector<float> *MuPhi, std::vector<float> *MuE, std::vector<unsigned int> *MuId, std::vector<unsigned int> *MuIdTight, std::vector<float> *MuPfIso){
  bool Mallpass = false;
  for(unsigned int i = 0 ; i<MuPt->size() ; i++){
    bool passEta = false, passIso = false, passId = false, passPt = false;
    TLorentzVector currentLepton; currentLepton.SetPtEtaPhiE(MuPt->at(i),MuEta->at(i),MuPhi->at(i),MuE->at(i));
    passId = MuIdTight->at(i) & (1<<0); //Look at the first vertex, hence the bit 0.
    int eta = fabs(MuEta->at(i));
    passEta = (eta<=2.4);
    passIso = (MuPfIso->at(i)<0.15); //Numbers are taken from llvv_fwk and have not been checked.
    passPt = (currentLepton.Pt() >=25);
    Mallpass= Mallpass||(passEta && passIso && passId && passPt);
    }
  return Mallpass;
  }  

bool TagProbePruner::passTrigger(int trig){
  vector<vector<int> > trigList(Ntrig);
  trigList[SingleMu].insert(trigList[SingleMu].end(),trigSingleMu,trigSingleMu+(sizeof(trigSingleMu)/sizeof(trigSingleMu[0])));
  trigList[SingleE].insert(trigList[SingleE].end(),trigSingleE,trigSingleE+(sizeof(trigSingleE)/sizeof(trigSingleE[0])));
  
  switch(trig){
  case SingleMu:
    for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return true;
    break;
  case SingleE:
    for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1<<trigList.at(SingleE).at(i))) return true;
    break;
  default:
    return false;
  }
return false; //If nothing found.
}


bool TagProbePruner::eventSelection(){

  bool passTightEle;
  bool passTightMu;
  switch(iSubSelection_){
  case Mu:
    passTightMu = isTightMuon(MuPt,MuEta,MuPhi,MuE,MuId,MuIdTight,MuPfIso);
  return MuPt->size() > 1 && passTightMu && passTrigger(SingleMu);
  case Ele:
    passTightEle = isTightEle(ElPt,ElEta,ElPhi,ElE,ElId,ElEtaSc,ElPfIsoRho);
  return ElPt->size() > 1 && passTightEle && passTrigger(SingleE);
  default:
  return false;
  }
};

bool TagProbePruner::filterEl(int iEl) {
	
  return ElPt->size() >0&&(*ElPt)[iEl] > minLepPt;
  }

bool TagProbePruner::filterMu(int iMu){  

  return MuPt->size() > 0 && ((*MuPt)[iMu] > minLepPt) ;
  }







