import FWCore.ParameterSet.Config as cms


process = cms.Process("L1UpgradeJet")

# Initialize MessageLogger and output report    
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#from SimCalorimetry.Configuration.SimCalorimetry_cff import *

process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration/StandardSequences/L1HwVal_cff')
#process.load("Configuration.StandardSequences.RawToDigi_Data_cff") ###check this for MC!
#process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff") ###check this for MC!
#process.load("RecoParticleFlow.PFProducer.particleFlow_cff") 
#process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")
process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")



#process.GlobalTag.globaltag = 'POSTLS161_V12::All'
#process.GlobalTag.globaltag = 'POSTLS261_V3::All'

# 6_2_2
#process.GlobalTag.globaltag = 'START62_V1::All'
#process.GlobalTag.globaltag = 'START62_V1::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')

process.maxEvents = cms.untracked.PSet(
    # restrict number of events
#    input = cms.untracked.int32(100)
    # run over all events
    input = cms.untracked.int32(10000)
)

process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring(
#'file:FAB2AF55-B2D4-E211-B254-0025B3E05DDA.root'
#'/store/mc/UpgFall13d/HToTauTau_125_14TeV_powheg_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FAF4B2D5-0539-E311-B5C3-002618FDA279.root'
'/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FCF11D5B-0739-E311-9B6D-003048FFD796.root'
   ),
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
                                           
                                           'keep *_l1GctHwDigis_*_*',
                                           'keep *_ecalDigis_*_*',
                                           'keep *_hcalDigis_*_*',
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
                                           # Calotowers
                                           'keep *_towerMaker_*_*',

					   # Keep generator Jets 
					   'keep *_ak5GenJets_*_*'
),
)
process.PFcands = cms.EDProducer("L1CaloCandidate",
                                         ECALDigis      = cms.InputTag("ecalDigis:EcalTriggerPrimitives"),
                                         HCALDigis      = cms.InputTag("hcalDigis"),
                                         UseUpgradeHCAL = cms.bool(False),

                                     )
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
		       #process.PFcands
                       process.RawToDigi
                       #+process.lumiProducer
                       +process.valGctDigis

                       # L1 Upgrade jets
                       +process.SLHCCaloTrigger

)                                       
process.p1.insert(1, process.valHcalTriggerPrimitiveDigis)
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
process.valRctDigis.hcalDigis             = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis =  cms.InputTag("valHcalTriggerPrimitiveDigis")

process.outpath = cms.EndPath(process.o1)

process.options = cms.untracked.PSet(

    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

    
