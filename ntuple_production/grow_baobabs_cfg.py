import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import re

process = cms.Process("GrowBaobabs")

# setup 'analysis'  options
opt = VarParsing.VarParsing ('analysis')
# Addition options.
# Note: if you add an option, update the code which write the values 
# in the configuration dump, you can find at the end of this file.
opt.register('minRun', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, 'Gives an indication on the minimum run number included in the input samples.')
opt.register('maxRun', 999999, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, 'Gives an indication on the maximum run number include in the input samples.')
opt.register('prodEra', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, 'Production run era. Label used to identify a run period whose data are processed together.')
opt.register('recoTag', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, 'Tag of the recontruction.')
opt.register('dataTier', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, 'Data tier of input dataset, typically AOD, AODSIM, MINIAOD or MINIAODSIM')
opt.register('isMC',    1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, 'Flag indicating if the input samples are from MC (1) or from the detector (0).')
opt.register('makeEdm', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, 'Switch for EDM output production. Use 0 (default) to disable it, 1 to enable it.')

#input files. Can be changed on the command line with the option inputFiles=...
opt.inputFiles = [
#'/store/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/00993A51-DF90-E611-A4EE-7845C4FC3650.root'
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver1-v1/110000/02B27BFE-1BEB-E611-8D50-001EC9B20ECB.root'
#"/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v1/70000/02477A4E-C586-E611-BC6F-02163E013D1C.root"
#'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/680BED0F-D919-E611-85E6-02163E01424F.root'
#'/store/mc/RunIIFall15MiniAODv2/TTbarDMJets_pseudoscalar_Mchi-1_Mphi-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/0A4E9031-7CB9-E511-8ABE-02163E00EA21.root'
#"file:/tmp/hbrun/theDYfile.root",
#"file:/tmp/hbrun/thePhotonData.root"
#'/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0C765598-8BD1-E511-BF63-20CF3027A566.root'
#'/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/0C6D4AB0-6F6C-E511-8A64-02163E0133CD.root'
#'/store/data/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/00039A2E-D7A7-E511-98EE-3417EBE64696.root'
]

#max number of events. #input files. Can be changed on the command line with the option maxEvents=...
opt.maxEvents = 1000

opt.parseArguments()


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.options.allowUnscheduled = cms.untracked.bool(True)

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.source = cms.Source("PoolSource",
                            fileNames =  cms.untracked.vstring(opt.inputFiles))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(opt.maxEvents))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('ntuple.root' )
)

if opt.isMC < 0 and len(process.source.fileNames) > 0:
  if re.match(r'.*/(MINI)?AODSIM/.*', process.source.fileNames[0]):
    print "MC dataset detected."
    opt.isMC = 1
  elif re.match(r'.*/(MINI)?AOD/.*', process.source.fileNames[0]):
    print "Real data dataset detected."
    opt.isMC = 0
  #endif
#endif

if opt.isMC < 0:
  raise Exception("Failed to detect data type. Data type need to be specify with the isMC cmsRun command line option")
#endif

if opt.prodEra in [ "13TeV_25ns", "13TeV_25ns_silver", "13TeV_25ns_silver"]: 
#for 76x:
    dataGlobalTag = '76X_dataRun2_16Dec2015_v0'
    mcGlobalTag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'
    triggerMenu = '2015'
    reapply_jec = True
    jec_file = False
    eg_corr = True   #photon and electron correction
    eg_corr_phot_file = "EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015"
    eg_corr_el_file   = "EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015"
else:
    #2016 data

    if opt.maxRun <= 280385: 
        #run B-G period
        dataGlobalTag = '80X_dataRun2_2016SeptRepro_v7'
    else:
        #run H period
        dataGlobalTag ='80X_dataRun2_Prompt_v16'
    #endif

    #mcGlobalTag = '80X_mcRun2_asymptotic_2016_v3'
    mcGlobalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
    triggerMenu = '2016'
    reapply_jec = False
    eg_corr = True
    eg_corr_phot_file = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_23Jan_ele"  
    eg_corr_el_file   = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_23Jan_ele"
    #eg_corr_phot_file = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/80X_ichepV2_2016_pho"
    #eg_corr_el_file   = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/80X_ichepV1_2016_ele"
    reRun_METfilter = True # reRun bad muons and bad charge hadrons filters

#endif

include_ak08 = True #switch to include anti-kt R=0.8 jets. ak(a) fatjet



#------------------------------------
#Condition DB tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if opt.isMC == 1:
  process.GlobalTag = GlobalTag(process.GlobalTag, mcGlobalTag, '')
else:
  process.GlobalTag = GlobalTag(process.GlobalTag, dataGlobalTag, '')
#fi

#--------------------------------------
#JEC
#
if reapply_jec:
  if jec_file:
    process.load("CondCore.DBCommon.CondDBCommon_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jec = cms.ESSource("PoolDBESSource",
          DBParameters = cms.PSet(
            messageLevel = cms.untracked.int32(0)
            ),
          timetype = cms.string('runnumber'),
          toGet = cms.VPSet(
          cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string(jec_file_tag),
                label  = cms.untracked.string('AK4PFchs')
                ),
          #Modified by Clement Leloup
          cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string(jec_file_tag),
                label  = cms.untracked.string('AK8PFchs')
                ),
          ), 
          connect = cms.string('sqlite:' + jec_file)
    )
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
  #endif jec_file

  jec_levels = ['L1FastJet', 'L2Relative', 'L3Absolute'] 
  
  if not opt.isMC:
    jec_levels.append('L2L3Residual')
  #endif


  #Modified by Clement Leloup
  from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
  process.load("RecoJets.JetProducers.PileupJetID_cfi")
  process.pileupJetIdUpdated = process.pileupJetId.clone(
    jets=cms.InputTag("slimmedJets"),
    inputIsCorrected=True,
    applyJec=True,
    vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
    )

  updateJetCollection(
    process,
    jetSource =cms.InputTag('slimmedJets'),
    labelName ='UpdatedJEC',
    jetCorrections = ('AK4PFchs', jec_levels, 'None')
  )

  jetSrc = 'updatedPatJetsUpdatedJEC'

  if include_ak08:
    updateJetCollection(
      process,
      jetSource = cms.InputTag('slimmedJetsAK8'),
      labelName = 'AK8UpdatedJEC',
      jetCorrections = ('AK8PFchs', jec_levels, 'None')
    )
    fatJetSrc = 'updatedPatJetsAK8UpdatedJEC'
    fatJetSw = 'on'
  else:
    fatJetSrc = ''
    fatJetSw = 'off'
  #endif //include_ak08
  puMvaName = 'pileupJetIdUpdated:fullDiscriminant'
  process.updatedPatJetsUpdatedJEC.userData.userFloats.src += [puMvaName]


  ### ---------------------------------------------------------------------------
  ### Removing the HF from the MET computation
  ### ---------------------------------------------------------------------------
  #process.noHFCands = cms.EDFilter("CandPtrSelector",
  #                                 src=cms.InputTag("packedPFCandidates"),
  #                                 cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
  #                                 )

  #jets are rebuilt from those candidates by the tools, no need to do anything else
  ### =================================================================================

  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
  process.load('Configuration.StandardSequences.MagneticField_38T_cff')

  #default configuration for miniAOD reprocessing, change the isData flag to run on data
  #for a full met computation, remove the pfCandColl input
  # runMetCorAndUncFromMiniAOD(process,
  #                           isData = not opt.isMC,
  #                           )

  # runMetCorAndUncFromMiniAOD(process,
  #                            isData=not opt.isMC,
  #                            pfCandColl=cms.InputTag("noHFCands"),
  #                            reclusterJets=True,   #needed for NoHF
  #                            recoMetFromPFCs=True, #needed for NoHF
  #                            postfix="NoHF"
  #                            )
else:
  jetSrc = "slimmedJets"

  #Modified by Clement Leloup
  if include_ak08:
    fatJetSrc = "slimmedJetsAK8"
    fatJetSw = "on"
  else:
    fatJetSrc = ""
    fatJetSw = "off"
  #endif include_ak08

  puMvaName = 'pileupJetId:fullDiscriminant'
#endif reapply_jec
#
#--------------------------------------------

# Photon and electron correction
#

if eg_corr:
  from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
  process = regressionWeights(process)
  process.load("EgammaAnalysis.ElectronTools.regressionApplication_cff")

  process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                     calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                         engineName = cms.untracked.string('TRandom3'),
                                                                                         ),
                                                     calibratedPatPhotons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                       engineName = cms.untracked.string('TRandom3'),
                                                                                       ),
                                                     )

  process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
  process.calibratedPatElectrons.isMC = cms.bool(opt.isMC != 0)
  process.calibratedPatElectrons.correctionFile = cms.string(eg_corr_el_file)

  process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
  process.calibratedPatPhotons.isMC = cms.bool(opt.isMC != 0)
  process.calibratedPatPhotons.correctionFile = cms.string(eg_corr_phot_file)


  #Old Version
  # copied from  EgammaAnalysis/ElectronTools/python/calibratedPhotonsRun2_cfi.py:
  #process.calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducerRun2",
  #                                              # input collections
  #                                              photons = cms.InputTag('slimmedPhotons'),                                                                                     # data or MC corrections
  #                                              # if isMC is false, data corrections are applied
  #                                              isMC = cms.bool(opt.isMC != 0),
  #                                              # set to True to get special "fake" smearing for synchronization. Use JUST in case of synchronization
  #                                              isSynchronization = cms.bool(False),
  #                                              correctionFile = cms.string(eg_corr_phot_file)
  #                                              )
  
  #Old Version
  #copied from  EgammaAnalysis/ElectronTools/python/calibratedElectronsRun2_cfi.py')
  #process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2", 
  #                                                # input collections
  #                                                electrons = cms.InputTag('slimmedElectrons'),
  #                                                gbrForestName = cms.string("gedelectron_p4combination_25ns"),
  #                                                # data or MC corrections
  #                                                # if isMC is false, data corrections are applied
  #                                                isMC = cms.bool(opt.isMC != 0),
  #                                                # set to True to get special "fake" smearing for synchronization. Use JUST in case of synchronization
  #                                                isSynchronization = cms.bool(False),
  #                                                correctionFile = cms.string(eg_corr_el_file)
  #                                                )
  electronSrc = "calibratedPatElectrons"
  photonSrc   = "calibratedPatPhotons"
else:
  electronSrc = "slimmedElectrons"
  photonSrc   = "slimmedPhotons"
#--------------------------------------------
#     MET filters (the 2 MET filters that need to be ruRun on top of the miniAOD for both data and MC)
#                 (cf https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Moriond_2017)
if reRun_METfilter:
    process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
    process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
    process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
    process.BadPFMuonFilter.taggingMode = cms.bool(True)

    process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
    process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
    process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
    process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
                                 ThePhotons = cms.InputTag("slimmedPhotons"),
                                 TheJets = cms.InputTag("slimmedJets"),
                                 L1Maps = cms.string("L1PrefiringMaps_new.root"), # update this line with the location of this file
                                 DataEra = cms.string("2016BtoH"), #Use 2016BtoH for 2016
                                 UseJetEMPt = cms.bool(False), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
                                 PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
                                 )


#--------------------------------------------

# Photon and Electron VID
#
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)


# define which IDs we want to produce
my_id_modules = [
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff'
                 ]

my_id_modulesPhotons = [
                 'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
                 ]

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

for idmod in my_id_modulesPhotons:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

#process.selectedPhotons = cms.EDFilter('PATPhotonSelector',
#                                       src = cms.InputTag(photonSrc),
#                                       cut = cms.string('pt>5 && abs(eta)')
#                                       )

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag(electronSrc) #we want to apply the selection on top of the calibrated photons and electrons
process.egmPhotonIDs.physicsObjectSrc = cms.InputTag(photonSrc)
process.egmPhotonIsolation.srcToIsolate = cms.InputTag(photonSrc)
process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag(photonSrc)
process.photonRegressionValueMapProducer.srcMiniAOD = cms.InputTag(photonSrc)
process.photonMVAValueMapProducer.srcMiniAOD = cms.InputTag(photonSrc)


#--------------------------------------------

  
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlineSlimmedPrimaryVertices')
    )

process.tupel = cms.EDAnalyzer("Tupel",
  triggerEvent = cms.InputTag( "patTriggerEvent" ),
  candidateSw  = cms.untracked.string("off"), #on, off, or withTrack
  candidateSrc = cms.untracked.InputTag("packedPFCandidates"),
  photonSw     = cms.untracked.string("on"), #on or off
  photonSrc    = cms.untracked.InputTag(photonSrc),
  electronSrc  = cms.untracked.InputTag(electronSrc),
  muonSrc      = cms.untracked.InputTag("slimmedMuons"),
  tauSrc       = cms.untracked.InputTag("slimmedTaus"),
  jetSrc       = cms.untracked.InputTag(jetSrc),
  fatJetSw     = cms.untracked.string(fatJetSw), #on or off
  fatJetSrc    = cms.untracked.InputTag(fatJetSrc),
  genSrc       = cms.untracked.InputTag("prunedGenParticles"),
  gjetSrc      = cms.untracked.InputTag('slimmedGenJets'),
  gfatJetSrc   = cms.untracked.InputTag('slimmedGenJetsAK8'),
  muonMatch    = cms.string( 'muonTriggerMatchHLTMuons' ),
  muonMatch2   = cms.string( 'muonTriggerMatchHLTMuons2' ),
  elecMatch    = cms.string( 'elecTriggerMatchHLTElecs' ),
  mSrcRho      = cms.untracked.InputTag('fixedGridRhoFastjetAll'),#arbitrary rho now
  CalojetLabel = cms.untracked.InputTag('slimmedJets'), #same collection now BB 
  metSrcs      = cms.VInputTag("slimmedMETs","slimmedMETsNoHF","slimmedMETsPuppi"),
  lheSrc       = cms.untracked.InputTag('externalLHEProducer'),
  puSrc        = cms.untracked.InputTag('slimmedAddPileupInfo'),
  puMvaName    = cms.untracked.string(puMvaName),
  puJetIdSrc   = cms.untracked.InputTag("pileupJetIdUpdated"),
  pvSrc        = cms.untracked.InputTag('goodOfflinePrimaryVertices'),
  reducedBarrelRecHitCollection = cms.InputTag("reducedEgamma","reducedEBRecHits"),
  reducedEndcapRecHitCollection = cms.InputTag("reducedEgamma","reducedEERecHits"),
  reducedPreshowerRecHitCollection = cms.InputTag("reducedEgamma","reducedESRecHits"),
  elecIDsMap = cms.VInputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1","egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose","egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium","egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight","egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
  phoIDsMap = cms.VInputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose","egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium","egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight"),
  checkOnFlyMETfilters = cms.untracked.string("off"),
  triggerObjectTag = cms.untracked.InputTag("selectedPatTrigger"),
  triggerMenu = cms.untracked.string(triggerMenu),
  effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt")
)

process.p = cms.Path()

if reapply_jec:
  process.p += cms.Sequence( process.pileupJetIdUpdated + process.patJetCorrFactorsUpdatedJEC + process.updatedPatJetsUpdatedJEC )
  if include_ak08:
    process.p += cms.Sequence( process.patJetCorrFactorsAK8UpdatedJEC + process.updatedPatJetsAK8UpdatedJEC )
  #endif include_ak08
#endif reapply_jec

if eg_corr:
  process.p += process.regressionApplication
  process.p += process.calibratedPatElectrons 
  process.p += process.calibratedPatPhotons

if reRun_METfilter:
  process.tupel.checkOnFlyMETfilters = cms.untracked.string("on")
  process.p += process.BadPFMuonFilter
  process.p += process.BadChargedCandidateFilter

process.p += process.goodOfflinePrimaryVertices
process.p += process.egmGsfElectronIDSequence
process.p += process.egmPhotonIDSequence
process.p += process.prefiringweight
process.p += process.tupel

if opt.makeEdm:
  process.out = cms.OutputModule("PoolOutputModule",
                                 fileName = cms.untracked.string('edm_out.root'),
                                 outputCommands = cms.untracked.vstring('keep *')
                                 )

  process.outpath = cms.EndPath(process.out)
#endif makeEdm

iFileName = "configDump_cfg.py"
file = open(iFileName,'w')
file.write(str(process.dumpPython()))
file.write('''
opt.minRun   = %d
opt.maxRun   = %d
opt.prodEra  = %s
opt.recoTag  = %s
opt.isMC     = %d
''' % (opt.minRun, opt.maxRun, opt.prodEra, opt.recoTag, opt.isMC))
file.close()

