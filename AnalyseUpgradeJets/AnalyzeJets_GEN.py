import FWCore.ParameterSet.Config as cms

from lutables_cfi import *

gFileName = cms.string('OutputJets.root')
#cms.options.SkipEvent = cms.untracked.vstring('ProductNotFound')
# **************************************************
# *                  Thresholds                    *
# **************************************************
# TT energy thresholds (E = ECAL, H = HCAL)
gTTEThreshold      = cms.vint32( 1, 2, 3, 4, 5 )
gTTHThreshold      = cms.vint32( 1, 2, 3, 4, 5 )
gTTEplusHThreshold = cms.vint32( 1, 2, 3, 4, 5 )
# **************************************************
# *                     Local                      *
# **************************************************
# Width of iEta rings to sample
gTTRingWidth       = cms.int32( 4 )
# PrePUS TT threshold fix - 2014-02-27_SingleMu_26Feb_11x11


process = cms.Process("Analyze")
process.load("FWCore.MessageService.MessageLogger_cfi")

# Stop stupid output after every event
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( 
  SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring( 'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/SingleMu_20Feb_11x11/JetCollections_224_1_l1V.root'),
                            fileNames = cms.untracked.vstring('file:JetCollections.root'),
   skipEvents = cms.untracked.uint32(0)
)
process.TFileService = cms.Service("TFileService",
fileName = gFileName,
# Dramatic decrease the time it takes to close the ROOT file
closeFileFast = cms.untracked.bool(True)                          
)

process.EventProducer = cms.EDProducer('EventProducer',

                                       # iEta range over which to analyse TTs
                                       TTiEtaRange     = cms.int32(28),
                                       # **************************************************
                                       # *               Calorimeter Towers                *
                                       # **************************************************
                                       # TEST PATTERN
                                       # CalorimeterTowers   = cms.InputTag("L1TestPatternCaloTowerProducer:"), #,"TEST","TEST_PATTERN"),
                                       # REAL TOWERS
                                       CalorimeterTowers   = cms.InputTag("L1CaloTowerProducer:","EVENT_DATA"),   
                                       # CaloTowers

                                       # **************************************************
                                       # *                  Thresholds                    *
                                       # **************************************************
                                       # TT energy thresholds (E = ECAL, H = HCAL)
             #                           TTEThreshold      = cms.vint32( 1, 2, 3, 4, 5 ),
#                                        TTHThreshold      = cms.vint32( 1, 2, 3, 4, 5 ),
#                                        TTEplusHThreshold = cms.vint32( 1, 2, 3, 4, 5 ),
                                       TTEThreshold      = cms.vint32( gTTEThreshold ),
                                       TTHThreshold      = cms.vint32( gTTHThreshold ),
                                       TTEplusHThreshold = cms.vint32( gTTEplusHThreshold ),
                                       
                                       # **************************************************
                                       # *                     Local                      *
                                       # **************************************************
                                       # Width of iEta rings to sample
                                       TTRingWidth       = gTTRingWidth,

)

process.JetProducer = cms.EDProducer('JetProducer',

                                     # **************************************************
                                     # *                      Jets                      *
                                     # **************************************************
                                     # ******************************
                                     # *       Cleaning cuts        *
                                     # ******************************
                                     
                                     minL1JetPt          = cms.double(1),
                                     maxL1JetEta         = cms.double(1.6),
                                     minRECOJetPt        = cms.double(1),
                                     maxRECOJetEta       = cms.double(1.6),


                                     # ******************************
                                     # *    TowerJet collections    *
                                     # ******************************
                                     ProtoTowerJet              = cms.InputTag("L1TowerJetProducer:"),
                                     #                                  FilteredCentralityJet= cms.InputTag("L1TowerJetCentralityFilter:FilteredTowerJets"),
                                     FilteredCentralityTowerJet = cms.InputTag("L1TowerJetCentralityFilter:"),
                                     Filtered1DTowerJet         = cms.InputTag("L1TowerJetFilter1D:"),

				     # L1 Jets to use ###############################################################
                                     PrePUSubTowerJet           = cms.InputTag("L1TowerJetPUSubtractedProducer:PrePUSubCenJets:L1UpgradeJet"),
				     ################################################################################

                                     # Current L1 Jets
                                     # ********************

                                     RCTEtaRegions = cms.vdouble( -3.0, -2.172, -1.74, -1.392, -1.044, -0.695, -0.348, 0.0,
                                                                  0.348, 0.695, 1.044, 1.392, 1.74, 2.172, 3.0),
                                     

                                     # Uncalibrated
                                     CurrentL1GCTUncalibCentralJet = cms.InputTag("valGctDigis:cenJets"),
                                     CurrentL1GCTUncalibTauJet     = cms.InputTag("valGctDigis:tauJets"),
                                     #extrajet = cms.VInputTag("l1extraParticles:Central","l1extraParticles:Tau",),

                                     # RECO jets
                                     # ********************

                                     # ak5 # Hijack one of these and make it read the ak5GenJets !
                                     PrePUSPreCalibCaloJets    = cms.InputTag("ak5GenJets"),
                                     #PUSPreCalibCaloJets       = cms.InputTag("PUsubAK5RawCaloJetProducer"),
                                     #PrePUSCaloJets            = cms.InputTag("PrePUsubAK5CaloJetProducer"),
                                     #PUsubCaloJets             = cms.InputTag("PUsubAK5CaloJetProducer"),
                                     # kt6
                                     #kt6Jets                   = cms.InputTag("kt6CaloJets"),

                                     #RecoVertices = cms.InputTag("offlinePrimaryVertices"),

)

# Jet histograms
process.JetHist = cms.EDAnalyzer('JetHist',

                                 # **************************************************
                                 # *                      Jets                      *
                                 # **************************************************


                                 # ******************************
                                 # *     FastJet collections    *
                                 # ******************************

                                 TTAk5L1Jet = cms.InputTag("L1FastjetProducer:TTAk5L1Jet"),
                                 
                                 
                                 # ******************************
                                 # *    TowerJet collections    *
                                 # ******************************
                                 ProtoTowerJet              = cms.InputTag("L1TowerJetProducer:"),
                                 #                                  FilteredCentralityJet= cms.InputTag("L1TowerJetCentralityFilter:FilteredTowerJets"),
                                 FilteredCentralityTowerJet = cms.InputTag("L1TowerJetCentralityFilter:"),
                                 Filtered1DTowerJet         = cms.InputTag("L1TowerJetFilter1D:"),
                                 PrePUSubTowerJet           = cms.InputTag("L1TowerJetPUSubtractedProducer:PrePUSubCenJets"),
                                 PUSubTowerJet              = cms.InputTag("L1CalibFilterTowerJetProducer:CenJets"),
                                 LocalPUSubTowerJet         = cms.InputTag("L1TowerJetPUSubtractedProducer:LocalPUSubCenJets"),
                                 
                                 # Calibrated jets
#                                 CalibratedPrePUSak5PrePUSTowerJet = cms.InputTag("JetCalibProducer:CalibratedTowerJetPrePUSak5PrePUSTower"),
                                 CalibratedPrePUSak5PUSTowerJet    = cms.InputTag("JetCalibProducer:CalibratedTowerJetPrePUSak5PUSTower"),
#                                  CalibratedPUSak5PUSTowerJet       = cms.InputTag("JetCalibProducer:CalibratedTowerJetPUSak5PUSTower"),
#                                  CalibratedLPUSak5PUSTowerJet      = cms.InputTag("JetCalibProducer:CalibratedTowerJetLPUSak5PUSTower"), 


                                 # Current L1 Jets
                                 # ********************
#                                 extrajet = cms.VInputTag("l1extraParticles:Central","l1extraParticles:Tau",),
                                     
                                 
                                 # RECO jets
                                 # ********************

                                 # Hand-calibrated ak5PUS 
                                 Calibratedak5PUSRawak5PUSL1Jet = cms.InputTag("JetCalibProducer:Calibratedak5PUSRawak5PUSL1Jet"),
                                 
                                 PUsubCaloJets = cms.InputTag("PUsubAK5CaloJetProducer"),
                                 jetCollection = cms.InputTag("ak5CaloJets"),                                  
                                 kt6Jets       = cms.InputTag("kt6CaloJets"),
                                 
                                 
                                 RecoVertices = cms.InputTag("offlinePrimaryVertices"),

				# Input from Previous step 
                                 PrePUSTowerJetL1Jet     = cms.InputTag("JetProducer:PrePUSTowerJetL1Jet"), # L1 Jets,
                                 PUSTowerJetL1Jet        = cms.InputTag("JetProducer:PUSTowerJetL1Jet"),
                                 LPUSTowerJetL1Jet       = cms.InputTag("JetProducer:LPUSTowerJetL1Jet"),

                                 # Calibrated jets
                                 # **********

                                 # PrePUS
#                                 CalibratedPrePUSak5PrePUSTowerJetL1Jet = cms.InputTag("JetCalibProducer:CalibratedTowerJetPrePUSak5PrePUSL1Jet"),
                                 CalibratedPrePUSak5PUSTowerJetL1Jet    = cms.InputTag("JetCalibProducer:CalibratedTowerJetPrePUSak5PUSL1Jet"),
                                 RecalibratedPrePUSak5PUSTowerJetL1Jet    = cms.InputTag("JetCalibProducer:RecalibratedTowerJetPrePUSak5PUSL1Jet"),
                                 NVTXRecalibratedPrePUSak5PUSTowerJetL1Jet    = cms.InputTag("JetCalibProducer:NVTXRecalibratedTowerJetPrePUSak5PUSL1Jet"),
                                 CalibratedPrePUSak5PUSTowerJetLt3L1Jet = cms.InputTag("JetCalibProducer:CalibratedTowerJetPrePUSak5PUSLt3L1Jet"),
#                                  CalibratedPrePUSak5PUSTowerJetNVTXLt15L1Jet = cms.InputTag("JetCalibProducer:CalibratedTowerJetPrePUSak5PUSNVTXLt15L1Jet"),
#                                  CalibratedPrePUSak5PUSTowerJetNVTXLt25L1Jet = cms.InputTag("JetCalibProducer:CalibratedTowerJetPrePUSak5PUSNVTXLt25L1Jet"),
#                                  CalibratedPrePUSak5PUSTowerJetNVTXLt50L1Jet = cms.InputTag("JetCalibProducer:CalibratedTowerJetPrePUSak5PUSNVTXLt50L1Jet"),
                                 # PUS
                                 # ----------
                                 CalibratedPUSak5PUSTowerJetL1Jet       = cms.InputTag("JetCalibProducer:CalibratedTowerJetPUSak5PUSL1Jet"),
                                 RecalibratedPUSak5PUSTowerJetL1Jet     = cms.InputTag("JetCalibProducer:RecalibratedTowerJetPUSak5PUSL1Jet"),

                                 # LPUS
                                 # ----------
                                 CalibratedLPUSak5PUSTowerJetL1Jet      = cms.InputTag("JetCalibProducer:CalibratedTowerJetLPUSak5PUSL1Jet"),
                                 RecalibratedLPUSak5PUSTowerJetL1Jet    = cms.InputTag("JetCalibProducer:RecalibratedTowerJetLPUSak5PUSL1Jet"), 


                                 
                                 CurrentUncalibJetL1Jet  = cms.InputTag("JetProducer:CurrentUncalibJetL1Jet"),
                                 CurrentJetL1Jet         = cms.InputTag("JetProducer:CurrentJetL1Jet"),
#                                 Ak5CaloJetL1Jet         = cms.InputTag("JetProducer:Ak5CaloJetL1Jet"),
#                                 Ak5CaloUncorrJetL1Jet   = cms.InputTag("JetProducer:Ak5CaloUncorrJetL1Jet"),


                                 PrePUSRawAk5CaloJetL1Jet  = cms.InputTag("JetProducer:PrePUSRawAk5CaloJetL1Jet"), 
                                 PUSRawAk5CaloJetL1Jet     = cms.InputTag("JetProducer:PUSRawAk5CaloJetL1Jet"),
                                 PrePUSAk5CaloJetL1Jet     = cms.InputTag("JetProducer:PrePUSAk5CaloJetL1Jet"),
                                # Ak5CaloJetL1Jet           = cms.InputTag("JetProducer:Ak5CaloJetL1Jet"),  


				# This is the collection which is used to caibrate to it seems, so put the genJets there, which are now called JetProducer:PrePUSRawAk5CaloJetL1Jet
                                 Ak5CaloJetL1Jet           = cms.InputTag("JetProducer:PrePUSRawAk5CaloJetL1Jet"),  
                                ################################################## 

                                 NumPrimaryVertices      = cms.InputTag("JetProducer:NVTX"),


                                 # Specify the boundaries of the eta regions in which to apply local PU subtraction
                                 LocalRhoEtaDivisions = cms.vdouble( -3.0, -1.3, 0.0, 1.3, 3.0 ),
                                 # Specify the boundaries of PT in which to apply local PT binning
                                 PTSlice              = cms.vdouble( 20., 30., 40., 50., 75., 100., 150. ),
                                 # Specify the boundaries of PU in which to apply local PU binning
#                                 PUSlice              = cms.vdouble( 0., 10., 20., 30., 40., 50., 75., 100., 150. ),
                                 #L1TriggerUpgradePerfStudies recommendations
                                 PUSlice              = cms.vdouble( 0., 15., 25., 150. ),                                
                                 # Specify the PT thresholds
                                 PTThreshold          = cms.vdouble( 0., 10., 20., 30., 40., 50., 60. ),
                                 # Eta region-level segmentation
#                                  EtaRegionSlice       = cms.vdouble( -3.0, -2.172, -1.74, -1.392, -1.044, -0.695, -0.348, 0.0,
#                                                                      0.348, 0.695, 1.044, 1.392, 1.74, 2.172, 3.0),

                                 # SUB-REGION LEVEL !!!!
#                                  EtaRegionSlice       = cms.vdouble( -3.0, -2.5, -2.172, -1.93, -1.74, -1.566, -1.392, -1.218,
#                                                                      -1.044, -0.87, -0.695, -0.522, -0.348, -0.174,
#                                                                      0.0,
#                                                                      0.174, 0.348, 0.522, 0.695, 0.87, 1.044,
#                                                                      1.218, 1.392, 1.566, 1.74, 1.93, 2.172, 2.5, 3.0),

                                 # SUB-REGION LEVEL !!!!
                                 EtaRegionSlice       = cms.vdouble(-3.0,-2.65,-2.5,-2.322,-2.172,-2.043,-1.93,-1.83,-1.74,-1.653,
                                                                    -1.566,-1.4790,-1.3920,-1.3050,-1.2180,-1.1310,-1.0440,-0.9570,
                                                                    -0.8700,-0.7830,-0.6950,-0.6090,-0.5220,-0.4350,-0.3480,-0.2610,
                                                                    -0.1740,-0.0870,
                                                                    0.,
                                                                    0.0870,0.1740,
                                                                    0.2610,0.3480,0.4350,0.5220,0.6090,0.6950,0.7830,0.8700,
                                                                    0.9570,1.0440,1.1310,1.2180,1.3050,1.3920,1.4790,1.566,
                                                                    1.653,1.74,1.83,1.93,2.043,2.172,2.322,2.5,2.65,3.0),
      
                                 # *************************************
                                 # *    Online-offline jet matching    *
                                 # *************************************

                                 # Matching parameters
                                 # -------------------
#                                 MaxDeltaR = cms.double(0.5),                                     
                                 MaxDeltaR = cms.double(0.7),                                     
                                 
                                 # iEta bin width to sample to obtain calibration factors
                                 iEtaCalibrationBinWidth = cms.double( 1 ),


                                 # *************************************
                                 # *        Turn on parameters         *
                                 # *************************************


                                 # Ak5 cleaning cuts (Loose)
                                 # ------------------------------------------------------------
                                 MaxEta         = cms.double(1.6),
#                                 MinPt          = cms.double(30),
                                 MinPt          = cms.double(1),

                                 # TODO: APPLY LOOSE JET SELECTIONS
                                 
                                 #
                                 # Jet selection here:
                                 #
                                 # https://cmssdt.cern.ch/SDT/doxygen/CMSSW_6_1_0/doc/html/d3/dda/SelectorUtils_2interface_2JetIDSelectionFunctor_8h_source.html


                                 onlineTurnOnSingleJetPT = cms.vint32( 30, 40, 50, 75, 100, 150 ),


)
# unpacker

process.p = cms.Path(
                      process.JetProducer
                      +process.JetHist
                     )

