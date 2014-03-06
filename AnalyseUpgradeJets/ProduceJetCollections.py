import FWCore.ParameterSet.Config as cms


process = cms.Process("L1UpgradeJet")

# Initialize MessageLogger and output report    
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000



process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
process.load("RecoParticleFlow.PFProducer.particleFlow_cff") 
process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")
process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")



#process.GlobalTag.globaltag = 'POSTLS161_V12::All'
#process.GlobalTag.globaltag = 'POSTLS261_V3::All'

# 6_2_2
process.GlobalTag.globaltag = 'START62_V1::All'

process.maxEvents = cms.untracked.PSet(
    # restrict number of events
#    input = cms.untracked.int32(100)
    # run over all events
    input = cms.untracked.int32(-1)
)





process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring(


#    '/store/data/Run2012C/SingleMu/RECO/PromptReco-v2/000/199/021/1864501B-13D1-E111-B851-0025901D5DF4.root'   
#    'file:/home/hep/mb1512/JetAnalyser/CMSSW_6_0_1_PostLS1v2_patch4/src/SingleMu.root'


    'file:/home/hep/mb1512/CMSSW/CMSSW_6_0_1_PostLS1v2_patch4/src/CopyData/ZeroBias4.root'



#'/store/data/Run2012C/SingleMu/RECO/PromptReco-v2/000/198/941/16D463F7-6BCF-E111-982E-001D09F25460.root'
#    'root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/data/Run2012C/ZeroBias/RAW/v1/000/198/609/F0C83922-74CA-E111-909B-003048673374.root',
#    '/store/data/Run2012C/ZeroBias/RAW/v1/000/198/588/24EC96E4-1BCA-E111-98AD-003048F118C2.root'


    # Local 140 PU Neutrino gun sample
 #   'file:/vols/cms04/jm1103/stuff.root'
#'/store/data/Run2012C/ZeroBias4/RAW-RECO/25Feb2013-v1/10000/001186BD-DD7F-E211-BDBE-0026189438E4.root'    
    #run locally at IC
#    'root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/data/Run2012C/ZeroBias/RECO/PromptReco-v1/000/198/609/0E478DDC-26CD-E111-90D2-001D09F2B30B.root',

#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_1000_1_goh.root',
#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_1001_1_J7Q.root',
#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_100_1_7he.root',
#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_1002_1_kzY.root',
#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_1003_1_hX7.root'
    
   ),

   # RAW files
   secondaryFileNames = cms.untracked.vstring(

#    '/store/data/Run2012C/SingleMu/RECO/PromptReco-v2/000/199/021/1864501B-13D1-E111-B851-0025901D5DF4.root'   
#    'file:/home/hep/mb1512/JetAnalyser/CMSSW_6_0_1_PostLS1v2_patch4/src/SingleMu.root'

    'file:/home/hep/mb1512/CMSSW/CMSSW_6_0_1_PostLS1v2_patch4/src/CopyData/ZeroBias4.root'

#'/store/data/Run2012C/SingleMu/RECO/PromptReco-v2/000/198/941/16D463F7-6BCF-E111-982E-001D09F25460.root'

#    'root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/data/Run2012C/ZeroBias/RAW/v1/000/198/609/F0C83922-74CA-E111-909B-003048673374.root',
#    '/store/data/Run2012C/ZeroBias/RAW/v1/000/198/588/24EC96E4-1BCA-E111-98AD-003048F118C2.root'
    
#   '/store/data/Run2012C/ZeroBias4/RAW-RECO/25Feb2013-v1/10000/001186BD-DD7F-E211-BDBE-0026189438E4.root'    
   	#run locally at IC
#    'root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/data/Run2012C/ZeroBias/RAW/v1/000/198/609/F0C83922-74CA-E111-909B-003048673374.root',
#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_1000_1_goh.root',
#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_1001_1_J7Q.root',
#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_100_1_7he.root',
#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_1002_1_kzY.root',
#     'root://gfe02.grid.hep.ph.ic.ac.uk:1094/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sives/SingleMu2012C_21Oct001_ALL/TT_calotowers_jets_dump_1003_1_hX7.root'


    			# 'root://eoscms//eos/cms///	
   )
)



# Make calotowers
process.towerMaker.hbheInput = cms.InputTag("hbheprereco")
process.towerMakerWithHO.hbheInput = cms.InputTag("hbheprereco")


# Get lumi
process.load("RecoLuminosity.LumiProducer.lumiProducer_cff")
from RecoLuminosity.LumiProducer.lumiProducer_cff import * 



# ***************************************************************************
# *                               Running processes                         *
# ***************************************************************************

                                                                                   

process.o1 = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('JetCollections.root'),
    #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('trigger_')),
    outputCommands = cms.untracked.vstring('drop *_*_*_*',

                                           # L1 objects
                                           # --------------------------------------------------
                                           
                                           # TT inputs
                                           'keep *_L1TestPatternCaloTowerProducer_*_*',
                                           'keep *_L1CaloTowerProducer_*_*',
                                           'keep *_L1RingSubtractionProducer_*_*',

                                           # GCT digis
                                           'keep *_gctDigis_*_*',
                                           # Current L1 emulation
#                                           'keep L1CaloRegions_l1GctHwDigis_*_*',
                                           
                                          # 'keep *_l1GctHwDigis_*_*',
                                           'keep *_valGctDigis_*_*',

                                           # L1 Upgrade jets
                                           'keep *_SLHCL1ExtraParticles_*_*',
                                           'keep *_L1TowerJetProducer_*_*',
                                           'keep *_L1TowerJetCentralityFilter_*_*',
                                           'keep *_L1TowerJetFilter1D_*_*',
                                           'keep *_L1TowerJetPUEstimator_*_*',
                                           'keep *_L1TowerJetPUSubtractedProducer_*_*',
                                           'keep *_L1CalibFilterTowerJetProducer_*_*',
                                           
 					   'keep *_l1extraParticles_*_*',

                                           # RECO objects
                                           # --------------------------------------------------

                                           # Calotowers
                                           'keep *_towerMaker_*_*',

                                           # RECO jet collections

                                           # kt6 calojets
                                           'keep *_kt6CaloJets_*_RECO',
                                           'keep *_ak5CaloJets_rho_L1UpgradeJet',
                                           # ak5 calojets
                                           'keep *_ak5CaloJets_*_RECO',
                                           'drop *_ak5CaloJets_rho_RECO',
                                           'keep *_ak5CaloJets_rho_L1UpgradeJet',
                                           'keep *_PUsubAK5CaloJetProducer_*_*',
                                           'keep *_PrePUsubAK5CaloJetProducer_*_*',
                                           'keep *_PUsubAK5RawCaloJetProducer_*_*',

                                           # ak5 PF jets
                                           'keep *_ak5PFJets_*_RECO',
                                           'drop *_ak5PFJets_rho_RECO',
                                           'keep *_ak5PFJets_rho_L1UpgradeJet',


                                           'keep *_offlinePrimaryVertices_*_*',
                                           'keep *_ak5JetID_*_*',
                                           'keep recoVertexs_offlinePrimaryVertices_*_RECO',
                                           'keep *_lumiProducer_*_*',

                                           'keep *_L1CaloTowerFastjetProducer_*_*',

                                           'keep *_PUsubAK5FastCaloJetProducer_*_*',
                                           'keep *_PrePUsubAK5FastCaloJetProducer_*_*',
                                           'keep *_PUsubAK5RawFastCaloJetProducer_*_*',

#                                           'keep *_PFcands_*_*',
                                           
),
)





# *******************************************************
# *               Offline jet collections               *
# *******************************************************

# #create collection of PrePU corrected Pre Calibrated AK5 Calo Jets
# process.PrePUsubAK5RawCaloJetProducer = cms.EDProducer('CalCaloProducer',
#                                     CaloJets = cms.InputTag("ak5CaloJets"),
# #                                    JetCorrector = cms.string('ak5CaloL1Fastjet'),          # Offset jet corrections
# )

#create collection of PU corrected Pre Calibrated AK5 Calo Jets
process.PUsubAK5RawCaloJetProducer = cms.EDProducer('CalCaloProducer',
                                   CaloJets = cms.InputTag("ak5CaloJets"),
                                   JetCorrector = cms.string('ak5CaloL1Fastjet'),          # Offset jet corrections
)

#create collection of PrePU corrected AK5 Calo Jets
process.PrePUsubAK5CaloJetProducer = cms.EDProducer('CalCaloProducer',
                                   CaloJets = cms.InputTag("ak5CaloJets"),
                                   JetCorrector = cms.string('ak5CaloL2L3Residual'),       # L2, L3 and residual jet corrections
)
#create collection of PU corrected AK5 Calo Jets
process.PUsubAK5CaloJetProducer = cms.EDProducer('CalCaloProducer',
                                   CaloJets = cms.InputTag("ak5CaloJets"),
                                   JetCorrector = cms.string('ak5CaloL1FastL2L3Residual'), # Offset, L2, L3 and residual jet corrections
)


#
#                              FastJet reconstructed with abs(eta) < 3
#

#create collection of PU corrected Pre Calibrated AK5 Calo Jets
process.PUsubAK5RawFastCaloJetProducer = cms.EDProducer('CalCaloProducer',
                                   CaloJets = cms.InputTag("L1CaloTowerFastjetProducer:CaloAk5CaloJet"),
                                   JetCorrector = cms.string('ak5CaloL1Fastjet'),          # Offset jet corrections
)

#create collection of PrePU corrected AK5 Calo Jets
process.PrePUsubAK5FastCaloJetProducer = cms.EDProducer('CalCaloProducer',
                                   CaloJets = cms.InputTag("L1CaloTowerFastjetProducer:CaloAk5CaloJet"),
                                   JetCorrector = cms.string('ak5CaloL2L3Residual'),       # L2, L3 and residual jet corrections
)
#create collection of PU corrected AK5 Calo Jets
process.PUsubAK5FastCaloJetProducer = cms.EDProducer('CalCaloProducer',
                                   CaloJets = cms.InputTag("L1CaloTowerFastjetProducer:CaloAk5CaloJet"),
                                   JetCorrector = cms.string('ak5CaloL1FastL2L3Residual'), # Offset, L2, L3 and residual jet corrections
)




from RecoJets.JetProducers.JetIDParams_cfi import *

process.ak5JetID = cms.EDProducer('JetIDProducer', JetIDParams,
        src = cms.InputTag('ak5CaloJets')

)

process.ak5PFJets.doRhoFastjet = cms.bool(True)


# ak5
process.ak5CaloJets.doPVCorrection = cms.bool( False )
process.ak5CaloJets.doAreaFastjet  = cms.bool( True )
process.ak5CaloJets.doRhoFastjet   = cms.bool( True )
process.ak5CaloJets.doPUOffsetCorr = cms.bool( False )
# NOTE: CURRENTLY RESTRICTING RHO DETERMINATION TO |eta| < 2.5
process.ak5CaloJets.Rho_EtaMax     = cms.double( 2.5 )
process.ak5CaloJets.jetPtMin       = cms.double( 0.00001 )

# kt6 (For rho estimation)
process.kt6CaloJets.doPVCorrection = cms.bool( False )
process.kt6CaloJets.rParam         = cms.double(0.6)
process.kt6CaloJets.doRhoFastjet   = cms.bool( True )
process.kt6CaloJets.doAreaFastjet  = cms.bool( True )
process.kt6CaloJets.doPUOffsetCorr = cms.bool( False )
# NOTE: CURRENTLY RESTRICTING RHO DETERMINATION TO |eta| < 2.4
process.kt6CaloJets.Rho_EtaMax     = cms.double( 2.4 )






process.PFcands = cms.EDProducer("L1CaloCandidate",
                                         ECALDigis      = cms.InputTag("ecalDigis:EcalTriggerPrimitives"),
                                         HCALDigis      = cms.InputTag("hcalDigis"),
                                         UseUpgradeHCAL = cms.bool(False),

                                     )






# # unpacker
# process.load('EventFilter.GctRawToDigi.l1GctHwDigis_cfi')
# process.l1GctHwDigis.unpackerVersion = cms.uint32(3)
# #process.l1GctHwDigis.unpackSharedRegions = cms.bool ( True )
# process.l1GctHwDigis.inputLabel = cms.InputTag( "source" )
# process.l1GctHwDigis.numberOfGctSamplesToUnpack = cms.uint32(5)
# process.l1GctHwDigis.hltMode = cms.bool( False )
# process.l1GctHwDigis.verbose = cms.untracked.bool ( False )
# process.l1GctHwDigis.unpackFibres = cms.untracked.bool ( True )
# process.l1GctHwDigis.unpackInternEm = cms.untracked.bool ( True )
# process.l1GctHwDigis.unpackInternJets = cms.untracked.bool ( True )

# GCT emulator
process.load('L1Trigger.Configuration.L1StartupConfig_cff')
import L1Trigger.GlobalCaloTrigger.gctDigis_cfi
process.valGctDigis = L1Trigger.GlobalCaloTrigger.gctDigis_cfi.gctDigis.clone()

#process.valGctDigis.inputLabel = cms.InputTag( "l1GctHwDigis" )
process.valGctDigis.inputLabel = cms.InputTag( 'gctDigis' )
process.valGctDigis.writeInternalData = cms.bool(True)
process.valGctDigis.useImprovedTauAlgorithm = cms.bool(False)
process.valGctDigis.preSamples = cms.uint32(0)
process.valGctDigis.postSamples = cms.uint32(0)


# ## GCT jet seed thresholds
# # -------------------------
# # Set cenjet threshold. Tau jet seed copies cenjet threshold
# JetFinderCentralJetSeed = 5.0
# # Set forjet threshold
# JetFinderForwardJetSeed = 5.0

# process.load('L1TriggerConfig.GctConfigProducers.l1GctConfig_cfi')
# process.L1GctConfigProducers.JetFinderCentralJetSeed = cms.double(JetFinderCentralJetSeed)
# process.L1GctConfigProducers.JetFinderForwardJetSeed = cms.double(JetFinderForwardJetSeed)


# # Remove unnecessary unpacking

process.RawToDigi.remove(process.csctfDigis)
process.RawToDigi.remove(process.dttfDigis)
#process.RawToDigi.remove(process.gctDigis)
#process.RawToDigi.remove(process.gtDigis)
#process.RawToDigi.remove(process.gtEvmDigis)
process.RawToDigi.remove(process.siPixelDigis)
process.RawToDigi.remove(process.siStripDigis)
# process.RawToDigi.remove(process.ecalDigis)
# process.RawToDigi.remove(process.ecalPreshowerDigis)
# process.RawToDigi.remove(process.hcalDigis)
process.RawToDigi.remove(process.muonCSCDigis)
process.RawToDigi.remove(process.muonDTDigis)
process.RawToDigi.remove(process.muonRPCDigis)

process.p1 = cms.Path(
                       process.RawToDigi
                       +process.lumiProducer

 #                      +process.L1GctConfigProducers
                       # Current emulation
#                       +process.l1GctHwDigis
                       +process.valGctDigis

                       # L1 Upgrade jets
                      +process.SLHCCaloTrigger

                       # RECO jets
                       +process.ak5CaloJets
                       +process.ak5PFJets
                       +process.kt6CaloJets

                      # Corrected jet collections
                       +process.PUsubAK5CaloJetProducer
                       +process.PrePUsubAK5CaloJetProducer
                       +process.PUsubAK5RawCaloJetProducer

                       # Corrected jet collections ( abs(eta) < 3 )
                       +process.PUsubAK5FastCaloJetProducer
                       +process.PrePUsubAK5FastCaloJetProducer
                       +process.PUsubAK5RawFastCaloJetProducer

                      +process.ak5JetID
#                       +process.PFcands
)                                       
process.outpath = cms.EndPath(process.o1)

process.options = cms.untracked.PSet(

    SkipEvent = cms.untracked.vstring('ProductNotFound')
    # Make the job crash in case of missing product
    #Rethrow = cms.untracked.vstring('ProductNotFound')
    
)

    
