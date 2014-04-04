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
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring(
#'file:FAB2AF55-B2D4-E211-B254-0025B3E05DDA.root'
#'/store/mc/UpgFall13d/HToTauTau_125_14TeV_powheg_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FAF4B2D5-0539-E311-B5C3-002618FDA279.root'

   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/00508A94-9E3A-E311-B0D8-002590593878.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0217D1CC-0239-E311-B29A-0025905964B4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0270AA59-0739-E311-9A91-00259059642A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/06685DF1-0639-E311-AA38-003048678F84.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0A374969-0239-E311-A37A-0026189438A0.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0A9CA0CF-0239-E311-9A71-003048FFD71A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0AF4B91D-0739-E311-B969-0026189438AB.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0C01A221-0739-E311-AA52-0025905938AA.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0C110B6C-0239-E311-86D2-00259059642E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0C3F21B4-0739-E311-A59B-003048FFD75C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0C4A0958-0739-E311-B7C1-002590593920.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0EC19659-0439-E311-920C-0025905938A4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0ED9B8CB-0239-E311-B157-003048678B38.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/104CC2F3-0639-E311-9BC1-002590593872.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/12204F22-0739-E311-AE07-0026189438B0.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/16B8DB58-0739-E311-B084-00259059391E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/16CC3126-0739-E311-9120-0025905964CC.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1814A01D-0739-E311-8110-00304867916E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/18DA0556-0739-E311-9347-0026189438A2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1A6106F0-0639-E311-B100-003048678BE8.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1AB32F1F-0739-E311-B2EB-003048678FA6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1C3D14B3-0739-E311-AE01-003048678E8A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1E34ACE5-0739-E311-A68F-0025905964C4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1E855BC7-0239-E311-889B-003048D3FC94.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/20ABBF1D-0739-E311-9336-003048678F9C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/20EF91CF-0239-E311-8F39-003048FFD71A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/22505A1D-0739-E311-9764-003048678BE6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/22AEFB1D-0739-E311-8FD7-0026189438CF.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/22E40C09-0239-E311-8964-003048FFCB9E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/240C0221-0739-E311-84B8-0025905938AA.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2441FD59-0739-E311-9FB0-0025905964C4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/24849826-0739-E311-87D0-003048FFD756.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/249F1CB7-9E3A-E311-AD0A-00261894386C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2635261D-0739-E311-90A4-00261894389E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2667FEE4-0739-E311-9F56-003048678B36.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/26EBB14A-0B39-E311-A992-002618943915.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/280787C8-4D3A-E311-99DA-002618943843.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/285108F0-0639-E311-A1C5-003048679296.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2A08B9B7-0739-E311-8655-0025905964BA.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2A20CEEF-0639-E311-A319-002618943905.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2ABDB21D-0739-E311-A32C-00304867BFB2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2C7C131F-0739-E311-9704-0026189437FE.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2C8A5D58-0739-E311-990B-002590596498.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2C8EDF09-0239-E311-B59E-003048FFCC2C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2C91EF65-0B39-E311-8A5C-003048678C62.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2E787A69-0239-E311-9081-00304867BFB2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2EB649C8-9E3A-E311-ACDE-0026189437FD.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3008DAC7-0239-E311-A5CD-003048678B72.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3009B3C8-9E3A-E311-8424-00304867918E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/30221AC8-0239-E311-B658-0026189437F9.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/30B6CE6B-0239-E311-8EA2-003048678B8E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/32543A06-0239-E311-A755-00304867BFAA.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3255F9F0-0639-E311-8492-003048678BE8.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/32CD2A54-0439-E311-B3CC-0026189438E1.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/32D55BE5-0739-E311-ACF2-0025905964C4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3470F755-0439-E311-A33B-003048679266.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3478A8B2-0739-E311-8273-0025905964C2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/34CA2CF1-0639-E311-9557-0025905964B6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/365ED357-0439-E311-8CE0-00304867904E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3AC5C5C9-0239-E311-BFCB-002618943902.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3AE03CC9-0239-E311-9DF1-0026189438CB.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3C467559-0439-E311-9554-003048FFD752.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3C661370-0239-E311-9FD6-002590593902.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3EB7A7F2-9B3A-E311-BA37-002618943978.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3EC31E1F-0739-E311-AD87-002618943951.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/40384522-0739-E311-B4B7-0026189437EC.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/409C2ECA-0239-E311-BA73-003048678FA6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/42A6BD5E-D53B-E311-BBCF-0025905964B6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/447FE702-0239-E311-A25D-0030486791DC.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/44A93CC9-0239-E311-85B4-0026189438CB.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/44ED4B45-503A-E311-AD18-002618943922.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/46DCA5B7-0739-E311-80E4-0025905964BA.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/482FE5C8-0239-E311-95DE-002618943951.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/48CF2ACA-0239-E311-957A-002590596484.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/4C529356-0439-E311-8C56-003048678B00.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/4E651026-0739-E311-ACDB-0025905822B6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/4EAAA959-0439-E311-9A5B-002618943934.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/4EBD6FB4-0739-E311-B95D-003048FFD7A2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/526FE71D-0739-E311-9963-0026189438EB.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5286201E-0739-E311-8885-003048678B34.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/528C1A1E-0739-E311-9D7D-002354EF3BE1.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/52A61958-0439-E311-B928-0025905964C4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/54320D1F-0739-E311-997B-002618943838.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/54381525-0739-E311-B2A2-00259059391E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/548E5DE5-0739-E311-9A77-002590593876.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/54B285B7-0739-E311-9E4B-0025905964BA.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/54E0ECCC-0239-E311-A297-003048FFD796.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/56F4D61D-0739-E311-AEE4-0026189438CE.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/581F65B2-0739-E311-962F-003048678DD6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5821F36D-0239-E311-891E-002618943870.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/58FCE46E-0239-E311-AC9E-003048FFD732.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5A1F71C7-0239-E311-9A74-002618943882.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5C8DCF71-0239-E311-8A21-003048FFCB8C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5CA00E3B-0039-E311-BD48-00261894391F.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5CA26BCA-0239-E311-B156-00261894394B.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5ED2856F-0239-E311-A95D-003048FFD7C2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/609690F4-0639-E311-8C2A-002590593878.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/620602F3-0639-E311-8C54-00261894390B.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/668426CC-9E3A-E311-878D-0025905938A4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/68AE06E5-0739-E311-90BB-00261894385D.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/6A79141F-0739-E311-B758-003048678F84.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/6AAAA01B-0739-E311-83EC-002618943910.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/6ADACA25-0739-E311-AC08-003048FFD7A2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/6C0E9F70-0239-E311-BE9B-003048FFCB96.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/6EF435B5-0739-E311-95E3-002590596486.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7033A2EF-0639-E311-970F-002618FDA26D.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7072D527-0739-E311-B606-0025905938B4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/74822ACD-0239-E311-A651-003048FFD75C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/74864ACC-0239-E311-92D4-003048FF9AC6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7810CA6E-0239-E311-9C75-003048FF9AC6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/78B73758-0439-E311-85B9-003048FFCB96.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/78CD386B-0239-E311-8FDD-00261894397E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7ABC2B5A-0739-E311-A5EF-003048FFD732.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7AC110CB-0639-E311-8158-00261894387E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7AC272B2-0739-E311-8B29-0025905964B2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7C3975F1-0639-E311-AD04-00304867BFAE.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7C9363F4-0639-E311-AD67-003048FFD754.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7EA7BEC8-0239-E311-BCE0-0026189438AF.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7ECDBB53-0439-E311-BCBC-002618FDA277.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7EEDEA24-0739-E311-BD5C-003048FF9AA6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/80D9611F-0739-E311-99B7-002618943864.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/823DD4AF-0739-E311-9B50-003048678ED2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/824799B3-0739-E311-8E7B-003048678A78.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/825F611F-0739-E311-8BAE-00261894382A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/848BB1C8-0239-E311-BB50-00261894386C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/84BC8F53-0439-E311-81B5-002618FDA207.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/86AFDAB8-9E3A-E311-A0AB-00259059391E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/8879AC6D-0239-E311-97F1-002590596490.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/8A9B47FA-0639-E311-B864-002590593920.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/8C04A4CD-0239-E311-9034-00304867926C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/8EC9BF57-0439-E311-8DE6-003048FF9AC6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/902FE91D-0739-E311-BE47-00248C0BE018.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9049FEF3-0639-E311-BD97-003048FFD7BE.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/923F2F55-0439-E311-95B9-00248C65A3EC.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9443256F-0239-E311-BB7B-002590596468.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/946346B5-0739-E311-9D59-003048FFD7BE.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/94A3DFF0-0639-E311-8F4D-002618FDA279.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/94D72625-0739-E311-AFC7-003048FF9AA6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/94E41F06-0239-E311-9CDE-002618943901.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9619586E-0239-E311-91EE-003048FFD71A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/96A7061E-0739-E311-B196-003048678F8E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/981DC7EF-0639-E311-BC55-002354EF3BDF.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9A3B8468-0B39-E311-BEC1-00261894398D.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9CD004F6-0639-E311-9DE5-003048FFCBA4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9E36D66C-0B39-E311-9B11-003048FFD732.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9E79351F-0739-E311-94A9-00261894386D.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A001D9B4-0739-E311-A107-003048FFCB96.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A020DE09-0239-E311-8865-003048FFCB96.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A243ADC9-0239-E311-9580-003048678FD6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A24FB638-4B3A-E311-8F95-00248C55CC62.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A2A17C24-0739-E311-BF8A-003048FFD744.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A65E2655-0739-E311-9EFD-00261894386D.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A85B9D25-0739-E311-AEAE-003048FFCC18.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/AC1936F3-0639-E311-A60D-003048679006.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/ACFFEB63-0439-E311-9EB7-00248C55CC3C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/AE45AFCD-0239-E311-96A8-00304867BFF2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/AEAC57B4-0739-E311-8FF9-002590593920.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/AEE947B4-0739-E311-83FC-003048678DA2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B003D0D8-973A-E311-B22F-0026189438BA.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B046AD6C-483A-E311-9A52-00261894396F.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B0BA07CC-0239-E311-AB73-002590596468.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B2C4B125-0739-E311-B388-0025905938B4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B2E269E5-0739-E311-AC1C-0025905964A6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B4819A69-0239-E311-BB34-002618B27F8A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B4E8BFB4-0739-E311-A882-003048FFCC2C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B4ECA520-0739-E311-8069-003048678C06.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B6B7F602-0239-E311-8C57-00304866C398.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B801E01D-0739-E311-B3B1-002354EF3BE6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B886311F-0739-E311-A425-0026189438E6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BA90ACBD-9E3A-E311-8781-0025905964BA.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BA91421F-0739-E311-8B30-0026189437F0.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BCA6056D-0239-E311-801B-002618FDA208.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BE5AD2B2-0739-E311-8665-00259059391E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BEFD5122-0739-E311-B3E1-0026189438B0.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C00E70F1-9B3A-E311-B9E4-002618943973.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C01F1D24-0739-E311-9A47-002590596490.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C226BB1D-0739-E311-846A-0030486791F2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C28859CA-9E3A-E311-B049-003048678B34.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C4FF02C8-0239-E311-B761-0030486792A8.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C6E26008-0239-E311-BF17-003048FFCC1E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C83CE1EF-0639-E311-ACDA-002618FDA248.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C88B4104-0239-E311-87C0-0026189438DB.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C8F0C6C6-9E3A-E311-A6EB-002590596468.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/CAEB3823-0739-E311-82AB-0025905964B6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/CC3A8F6C-0239-E311-A67B-002590596490.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/CE438B4C-973A-E311-AC73-003048FFD754.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/CE6AB825-0739-E311-94B5-0025905938B4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D002FB20-0739-E311-8F36-002618943934.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D016D76C-0239-E311-89ED-003048FFCB9E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D0386D6E-0239-E311-8238-0025905964C4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D0B79FEF-0639-E311-81E8-00261894382A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D2498421-0739-E311-B99A-002590593920.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D2537853-0439-E311-A7E0-002618943834.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D451916E-0239-E311-81C6-003048FFCC0A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D45596F2-0639-E311-B156-0025905964BC.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D4B5EC22-0739-E311-A629-00259059642E.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D4EAEAEF-0639-E311-8313-0025905964BE.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D6028863-0439-E311-9D5E-003048679228.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D80DFA72-0239-E311-B89A-003048FFD7D4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D89937F0-0639-E311-AA2C-002618943901.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E0077A1E-0739-E311-9765-00304867BFB2.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E0FD0D6D-0B39-E311-8B57-003048FF9AA6.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E633585C-493A-E311-867A-002590593878.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E88D1DF7-0639-E311-81DB-0025905964C0.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E8BA71AA-DB39-E311-83A6-00259059649C.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/EA710826-0739-E311-94C8-0025905964CC.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/EAE493CC-0239-E311-AFFC-0025905938D4.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/ECF277E8-0639-E311-9EC2-0026189438EB.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F0DF3421-0739-E311-8539-002590593902.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F2D61EC8-9E3A-E311-B06E-00261894384A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F4142B1F-0739-E311-A659-003048D3C010.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F8725D26-0739-E311-A655-003048FFCB74.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F891356B-0239-E311-8594-00261894390A.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FA490F6E-0239-E311-9665-003048D15E02.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FCF11D5B-0739-E311-9B6D-003048FFD796.root',
   'root://eoscms.cern.ch//eos/cms/store/mc/UpgFall13d/QCD_Pt_300to470_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FECA99E5-0739-E311-8CD8-0025905938B4.root'
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
    fileName = cms.untracked.string('JetCollectionsQcd.root'),
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


