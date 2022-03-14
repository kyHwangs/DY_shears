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

static float minLepPt = 20;


enum {DoubleMu, SingleMu, DoubleE, HighPtE, SingleE, EMu, SinglePhoton, Ntrig}; //List of triggers used for our analysis
//The lines below correspond to the triggerBits used in our analysis. These are all defined in the trigger2016.h file of the ntuple producer.
static int trigDoubleMu[4] = {8,9,10,11};
static int trigSingleMu[4] = {10,11,15,16};
static int trigDoubleE[3] = {1,12,13};
static int trigHighPtE[1] = {17};//Located in DoubleElectron
static int trigSingleE[2] = {11,12};
static int trigEMu[2] = {0,3}; //all DZ paths are still missing, as well as Mu12Ele23
static int trigSinglePhoton[12] = {0,1,20,21,25,26,27,28,29,30,31,32};


/**************************/


class ZZ2l2vPruner: public Pruner, EventTree{
protected:
  void declareSubSelections();
  bool init(TChain* tree);
  bool filterEvent();
  
  bool filterEl(int iEl);
  bool filterMu(int iMu);

  bool passLepPt(vector<float> lepPt);
  
  void makeFilterMask(bool (ZZ2l2vPruner::*filter)(int), std::vector<bool>& mask);
  void skimCollections();
  bool eventSelection();
  bool passTrigger(int trig);
  enum {MC_DLep, DataDoubleMuDLep, DataSingleMuDLep, DataDoubleElDLep, DataSingleElDLep, DataElMuDLep, MC_DMu, MC_DE, DataDoubleMuDMu, DataDoubleMuDE, DataSingleMuDMu, DataSingleMuDE, DataDoubleElDMu, DataDoubleElDE, DataSingleElDMu, DataSingleElDE, DataElMuDMu, DataElMuDE, DataSinglePhoton, NSubSels};
};

DECLARE_PRUNER(ZZ2l2vPruner, "Pruner of ZZ2l2v analysis")

void ZZ2l2vPruner::declareSubSelections(){
  subSelections_.resize(NSubSels);
  subSelections_[MC_DLep]    = SubSelection("MC_DLep", "Dilepton selection for ZZ2l2v analysis, for MC samples.");
  subSelections_[DataDoubleMuDLep]     = SubSelection("DataDoubleMuDLep","Dilepton selection for ZZ2l2v analysis, for DoubleMu data");
  subSelections_[DataSingleMuDLep]     = SubSelection("DataSingleMuDLep","Dilepton selection for ZZ2l2v analysis, for SingleMu data");
  subSelections_[DataDoubleElDLep]     = SubSelection("DataDoubleElDLep","Dilepton selection for ZZ2l2v analysis, for DoubleEl data");
  subSelections_[DataSingleElDLep]     = SubSelection("DataSingleElDLep","Dilepton selection for ZZ2l2v analysis, for SingleEl data");
  subSelections_[DataElMuDLep]     = SubSelection("DataElMuDMu","Dilepton selection for ZZ2l2v analysis, for ElMu data");
  subSelections_[MC_DMu]    = SubSelection("MC_DMu", "Dimuon selection for ZZ2l2v analysis, for MC samples");
  subSelections_[MC_DE]     = SubSelection("MC_DE","Dielectron selection for ZZ2l2v analysis, for MC samples");
  subSelections_[DataDoubleMuDMu]     = SubSelection("DataDoubleMuDMu","Dimuon selection for ZZ2l2v analysis, for DoubleMu data");
  subSelections_[DataDoubleMuDE]     = SubSelection("DataDoubleMuDE","Dielectron selection for ZZ2l2v analysis, for DoubleMu data");
  subSelections_[DataSingleMuDMu]     = SubSelection("DataSingleMuDMu","Dimuon selection for ZZ2l2v analysis, for SingleMu data");
  subSelections_[DataSingleMuDE]     = SubSelection("DataSingleMuDE","Dielectron selection for ZZ2l2v analysis, for SingleMu data");
  subSelections_[DataDoubleElDMu]     = SubSelection("DataDoubleElDMu","Dimuon selection for ZZ2l2v analysis, for DoubleEl data");
  subSelections_[DataDoubleElDE]     = SubSelection("DataDoubleElDE","Dielectron selection for ZZ2l2v analysis, for DoubleEl data");
  subSelections_[DataSingleElDMu]     = SubSelection("DataSingleElDMu","Dimuon selection for ZZ2l2v analysis, for SingleEl data");
  subSelections_[DataSingleElDE]     = SubSelection("DataSingleElDE","Dielectron selection for ZZ2l2v analysis, for SingleEl data");
  subSelections_[DataElMuDMu]     = SubSelection("DataElMuDMu","Dimuon selection for ZZ2l2v analysis, for ElMu data");
  subSelections_[DataElMuDE]     = SubSelection("DataElMuDE","Dielectron selection for ZZ2l2v analysis, for ElMu data");
  subSelections_[DataSinglePhoton]     = SubSelection("DataSinglePhoton","Data passing only the SinglePhoton trigger");
}

bool ZZ2l2vPruner::init(TChain* tree){
  if(tree == 0) return false;
  EventTree::Init(tree);
  return true;
}

void ZZ2l2vPruner::makeFilterMask(bool (ZZ2l2vPruner::*filter)(int), std::vector<bool>& mask){
  for(unsigned i = 0; i < mask.size(); ++i){
    mask[i] = (this->*filter)(i);
  }
}

bool ZZ2l2vPruner::filterEvent(){
  //skimCollections(); //We don't want to run skimCollections because we need to keep all "bad" leptons in order to be able to apply a 3rd-lepton veto at the level of the analysis. However, if needed, this option can be reactivated simply by uncommenting this line.
  return eventSelection() && Pruner::filterEvent();
}
  

void ZZ2l2vPruner::skimCollections(){
  std::vector<bool> mask;
  
  //Reco muon collections:
  mask.resize(MuPt->size());
  makeFilterMask(&ZZ2l2vPruner::filterMu, mask);
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
  makeFilterMask(&ZZ2l2vPruner::filterEl, mask);
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


//to be run after skimCollections
bool ZZ2l2vPruner::eventSelection(){
  switch(iSubSelection_){
  case MC_DLep:
    return (passLepPt(*MuPt) || passLepPt(*ElPt)) && passTrigger(Ntrig); //For MC, requires any trigger to pass.
  case DataDoubleMuDLep:
    return (passLepPt(*MuPt) || passLepPt(*ElPt)) && passTrigger(DoubleMu);
  case DataSingleMuDLep:
    return (passLepPt(*MuPt) || passLepPt(*ElPt)) && passTrigger(SingleMu);
  case DataDoubleElDLep:
    return (passLepPt(*MuPt) || passLepPt(*ElPt)) && passTrigger(DoubleE);
  case DataSingleElDLep:
    return (passLepPt(*MuPt) || passLepPt(*ElPt)) && passTrigger(SingleE);
  case DataElMuDLep:
    return (passLepPt(*MuPt) || passLepPt(*ElPt)) && passTrigger(EMu);
  case MC_DMu:
    return passLepPt(*MuPt) && passTrigger(Ntrig);
  case DataDoubleMuDMu:
    return passLepPt(*MuPt) && passTrigger(DoubleMu);
  case DataSingleMuDMu:
    return passLepPt(*MuPt) && passTrigger(SingleMu);
  case DataDoubleElDMu:
    return passLepPt(*MuPt) && passTrigger(DoubleE);
  case DataSingleElDMu:
    return passLepPt(*MuPt) && passTrigger(SingleE);
  case DataElMuDMu:
    return passLepPt(*MuPt) && passTrigger(EMu);
  case MC_DE:
    return passLepPt(*ElPt) && !passLepPt(*MuPt) && passTrigger(Ntrig);
  case DataDoubleMuDE:
    return passLepPt(*ElPt) && !passLepPt(*MuPt) && passTrigger(DoubleMu);
  case DataSingleMuDE:
    return passLepPt(*ElPt) && !passLepPt(*MuPt) && passTrigger(SingleMu);
  case DataDoubleElDE:
    return passLepPt(*ElPt) && !passLepPt(*MuPt) && passTrigger(DoubleE);
  case DataSingleElDE:
    return passLepPt(*ElPt) && !passLepPt(*MuPt) && passTrigger(SingleE);
  case DataElMuDE:
    return passLepPt(*ElPt) && !passLepPt(*MuPt) && passTrigger(EMu);
  case DataSinglePhoton:
    return passTrigger(SinglePhoton);
  case NSubSels:
  default:
    return false;
  }
};

bool ZZ2l2vPruner::filterMu(int iMu){
  return ((*MuPt)[iMu] > minLepPt);
}

bool ZZ2l2vPruner::filterEl(int iEl){
  return ((*ElPt)[iEl] > minLepPt);
}

bool ZZ2l2vPruner::passLepPt(vector<float> lepPt){
  int goodLeptons = 0;
  for(unsigned int iLep = 0 ; iLep < lepPt.size() ; iLep++){
    if(lepPt[iLep] > minLepPt) goodLeptons ++;
  }
  return goodLeptons > 1;
}

bool ZZ2l2vPruner::passTrigger(int trig){
  vector<vector<int> > trigList(Ntrig);
  trigList[DoubleMu].insert(trigList[DoubleMu].end(),trigDoubleMu,trigDoubleMu+(sizeof(trigDoubleMu)/sizeof(trigDoubleMu[0])));
  trigList[SingleMu].insert(trigList[SingleMu].end(),trigSingleMu,trigSingleMu+(sizeof(trigSingleMu)/sizeof(trigSingleMu[0])));
  trigList[DoubleE].insert(trigList[DoubleE].end(),trigDoubleE,trigDoubleE+(sizeof(trigDoubleE)/sizeof(trigDoubleE[0])));
  trigList[HighPtE].insert(trigList[HighPtE].end(),trigHighPtE,trigHighPtE+(sizeof(trigHighPtE)/sizeof(trigHighPtE[0])));
  trigList[SingleE].insert(trigList[SingleE].end(),trigSingleE,trigSingleE+(sizeof(trigSingleE)/sizeof(trigSingleE[0])));
  trigList[EMu].insert(trigList[EMu].end(),trigEMu,trigEMu+(sizeof(trigEMu)/sizeof(trigEMu[0])));
  trigList[SinglePhoton].insert(trigList[SinglePhoton].end(),trigSinglePhoton,trigSinglePhoton+(sizeof(trigSinglePhoton)/sizeof(trigSinglePhoton[0])));
  switch(trig){
  case DoubleMu:
    for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return true;
    break;
  case SingleMu:
    for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return true;
    break;
  case DoubleE: //Includes also HighPtE
    for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return true;
    for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return true; //Accepted also for HighPtE
    break;
  case SingleE:
    for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1<<trigList.at(SingleE).at(i))) return true;
    break;
  case EMu:
    for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1<<trigList.at(SingleE).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[EMu].size() ; i++)  if(TrigHltElMu & (1<<trigList.at(EMu).at(i))) return true;
    break;
  case SinglePhoton:
    for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1<<trigList.at(SingleE).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[EMu].size() ; i++)  if(TrigHltElMu & (1<<trigList.at(EMu).at(i))) return false;
    for(unsigned int i = 0 ; i < trigList[SinglePhoton].size() ; i++)  if(TrigHltPhot & (1<<trigList.at(SinglePhoton).at(i))) return true;
    break;
  case Ntrig://In this case (used for MC), take if any trigger passed
    //for(unsigned int i = 0 ; i < trigList[SinglePhoton].size() ; i++)  if(TrigHltPhot & (1<<trigList.at(SinglePhoton).at(i))) return true; //It was decided not to take any photon triggers.
    for(unsigned int i = 0 ; i < trigList[EMu].size() ; i++)  if(TrigHltElMu & (1<<trigList.at(EMu).at(i))) return true;
    for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1<<trigList.at(SingleE).at(i))) return true;
    for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return true;
    for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return true;
    for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return true;
    for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return true;
    break;
  default:
    return false;
  }
  return false; //If nothing found.
}
