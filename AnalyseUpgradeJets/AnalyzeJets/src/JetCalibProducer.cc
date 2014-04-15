//
//          MODIFIED BY A ELWOOD FOR THE PU140 MC SAMPLE
//          
//          INVERTED THE FIT FUNCTION, L1060 and L1379 correction -> 1.0/correction
//
//          ADDED OPTION TO INVERT THE FUNCTION
//
// -*- C++ -*-
//
// Package:    JetCalibProducer
// Class:      JetCalibProducer
// 
/**\class JetCalibProducer JetCalibProducer.cc AnalyseUpgradeJets/src/JetCalibProducer.cc

 Description: Produces jets and event variables for the JetAnalyzer validation framework.

 List of quantities produced: 

     - NVTX
     - 
     -
     -
     -
     -

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Baber
//         Created:  Tue Oct  1 10:53:36 BST 2013
// $Id$
//
//



// Print debugging messages
//#define VERBOSE





// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMacro.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"




#include "DataFormats/JetReco/interface/JetID.h" 
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "SimDataFormats/SLHC/interface/EtaPhiContainer.h"
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/TriggerTowerGeometry.h"





#include "AnalyseUpgradeJets/AnalyzeJets/interface/printing.h"
// #include "AnalyseUpgradeJets/AnalyzeJets/interface/JetMatch.h"
// #include "AnalyseUpgradeJets/AnalyzeJets/src/helperFunctions.cc"

#include <algorithm>  // for sorting

#include <fstream>

// Ranking function for sort
bool towerJetRankDescending ( l1slhc::L1TowerJet jet1, l1slhc::L1TowerJet jet2 ){      return ( jet1.p4().Pt() > jet2.p4().Pt() ); }
bool L1JetRankDescending ( l1extra::L1JetParticle jet1, l1extra::L1JetParticle jet2 ){ return ( jet1.p4().Pt() > jet2.p4().Pt() ); }


// Maximum jet pT to consider when dumping calibration values
const int MAX_L1_PT = 300;


  const Int_t TTetaBin = 56;

  // TT eta binning
  const double TTetaBins[] = {-3.0,-2.65,-2.5,-2.322,-2.172,-2.043,-1.93,-1.83,-1.74,-1.653,
				-1.566,-1.4790,-1.3920,-1.3050,-1.2180,-1.1310,-1.0440,-0.9570,
				-0.8700,-0.7830,-0.6950,-0.6090,-0.5220,-0.4350,-0.3480,-0.2610,
				-0.1740,-0.0870,
				0.,
				0.0870,0.1740,
				0.2610,0.3480,0.4350,0.5220,0.6090,0.6950,0.7830,0.8700,
				0.9570,1.0440,1.1310,1.2180,1.3050,1.3920,1.4790,1.566,
				1.653,1.74,1.83,1.93,2.043,2.172,2.322,2.5,2.65,3.0};
  


struct LUT{

  std::vector <double> LUTArr;
  int columns;

  // Default constructor
  LUT(): columns(0){}
  
  // Get a specific element. Note: numbering starts from 0
  double getElement( int row, int column ){
    int index = column + row*columns;
    return LUTArr[ index ];
  }

};

//
// class declaration
//

class JetCalibProducer : public edm::EDProducer {
   public:
      explicit JetCalibProducer(const edm::ParameterSet&);
      ~JetCalibProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  
      // ----------member functions ---------------------------
      
//   virtual void calibrateJet( edm::Handle<l1slhc::L1TowerJetCollection> const& UncalibJet_Tower, std::vector< double > const& iEtaPtOffset, 
// 			     std::vector< double > const& iEtaPtCorrection, 
// 			     std::auto_ptr<l1slhc::L1TowerJetCollection>     &outputCalibTowerJet_Tower,
// 			     std::auto_ptr<l1extra::L1JetParticleCollection> &outputCalibTowerJet_L1Jet );

  virtual void calibrateJet( edm::Handle<l1slhc::L1TowerJetCollection> const& UncalibJet_Tower,
			     std::vector< double > const& iEtaPtCorrections,
			     const int&  fitOrder,
			     std::auto_ptr<l1slhc::L1TowerJetCollection>     &outputCalibTowerJet_Tower,
			     std::auto_ptr<l1extra::L1JetParticleCollection> &outputCalibTowerJet_L1Jet,
			     bool newMethod );
  virtual void calibrateJet( edm::Handle<l1extra::L1JetParticleCollection> UncalibJet_L1Jet,
			     std::vector< double > const& iEtaPtCorrections,
			     const int&  fitOrder, 
			     std::auto_ptr<l1extra::L1JetParticleCollection> &outputCalibJet_L1Jet,
			     std::auto_ptr<l1extra::L1JetParticleCollection> &outputRecalibJet_L1Jet,
			     std::auto_ptr<l1extra::L1JetParticleCollection> &outputOffsetRecalibJet_L1Jet,
			     bool newMethod,
			     LUT lowPtLUT,
			     LUT highPtLUT,
			     LUT offsetPtLUT,
			     double transitionPt);
			     //int NVTX);


  virtual void loadLUT( TString filename, std::vector< double > &LUTArr );

  virtual double getRecalibFactor(  double correctedPt, int row, LUT recalibLUT );

  virtual void printCalibrationFactors( TString LUT, std::vector< double > const& iEtaPtCorrections, const int& fitOrder, int maxL1Pt, bool newMethod );



      // ----------member data ---------------------------
      edm::ParameterSet conf_;
      edm::Service<TFileService> fs;

      //  L1 pT calibration threshold, minimum L1 jet pT to apply correction  
      double pTCalibrationThreshold;
      // Width of iEta bins to sample to obtain pT corrections
      double iEtaCalibBinWidth;
      // Jet calibration iEta-binned pT corrections
      std::vector <double> iEtaPtCorrectionPrePUSak5PrePUS, iEtaPtOffsetPrePUSak5PrePUS;
      std::vector <double> iEtaPtCorrectionPrePUSak5PUS,    iEtaPtOffsetPrePUSak5PUS;
      std::vector <double> iEtaPtCorrectionPUSak5PUS,       iEtaPtOffsetPUSak5PUS;

      int fitOrder;
      std::vector <double> pTCalibration_PrePUS_ak5PUS;
      std::vector <double> pTCalibration_PrePUS_ak5PUSLt3;
      std::vector <double> pTCalibration_PrePUS_NVTXLt15;
      std::vector <double> pTCalibration_PrePUS_NVTXLt25;
      std::vector <double> pTCalibration_PrePUS_NVTXLt50;
      std::vector <double> pTCalibration_PUS_ak5PUS;
      std::vector <double> pTCalibration_LPUS_ak5PUS;
      std::vector <double> pTCalibration_ak5PUSRaw_ak5PUS;
      std::vector <double> pTCalibration_UncalibGCT_ak5PUS;

      // Histogram containers
      std::map< TString, TH2*> hist2D;

      // Eta region-level segmentation
      std::vector <double> etaRegionSlice;


  // Recalibration
  // --------------------------------------------------
  
  // PrePUS
  LUT recalibLowPt_PrePUS_ak5PUS;
  LUT recalibHighPt_PrePUS_ak5PUS;
  LUT recalibNVTX_PrePUS_ak5PUS;
  // PUS
  LUT recalibLowPt_PUS_ak5PUS;
  LUT recalibHighPt_PUS_ak5PUS;
  LUT recalibNVTX_PUS_ak5PUS;
  // LPUS
  LUT recalibLowPt_LPUS_ak5PUS;
  LUT recalibHighPt_LPUS_ak5PUS;
  LUT recalibNVTX_LPUS_ak5PUS;


  double LUTPtTransition;
  // Use when no recalibration is required
  LUT emptyLUT;

  // Number of reconstructed ak5 primary vertices
  int NVTX;

  bool invertFunction;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
JetCalibProducer::JetCalibProducer(const edm::ParameterSet& iConfig): conf_(iConfig)
{


  // Eta region-level segmentation 
  etaRegionSlice        = iConfig.getParameter< std::vector<double> >("EtaRegionSlice");
  sort (etaRegionSlice.begin(), etaRegionSlice.end());  // ensure the bins are in ascending order  


//   // PrePUSak5PrePUS
//   produces <l1slhc::L1TowerJetCollection>(     "CalibratedTowerJetPrePUSak5PrePUSTower" );
//   produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPrePUSak5PrePUSL1Jet" );

  // PrePUSak5PUS
  produces <l1slhc::L1TowerJetCollection>(     "CalibratedTowerJetPrePUSak5PUSTower" );
  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPrePUSak5PUSL1Jet" );
  produces <l1extra::L1JetParticleCollection>( "RecalibratedTowerJetPrePUSak5PUSL1Jet" );
  produces <l1extra::L1JetParticleCollection>( "NVTXRecalibratedTowerJetPrePUSak5PUSL1Jet" );

  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPrePUSak5PUSLt3L1Jet" );

  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPrePUSak5PUSNVTXLt15L1Jet" );
  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPrePUSak5PUSNVTXLt25L1Jet" );
  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPrePUSak5PUSNVTXLt50L1Jet" );
  // PUSak5US
  //  produces <l1slhc::L1TowerJetCollection>(     "CalibratedTowerJetPUSak5PUSTower" );
  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPUSak5PUSL1Jet" );
  produces <l1extra::L1JetParticleCollection>( "RecalibratedTowerJetPUSak5PUSL1Jet" );
  produces <l1extra::L1JetParticleCollection>( "NVTXRecalibratedTowerJetPUSak5PUSL1Jet" );

  // LPUSak5PUS
  //  produces <l1slhc::L1TowerJetCollection>(     "CalibratedTowerJetLPUSak5PUSTower" );
  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetLPUSak5PUSL1Jet" );
  produces <l1extra::L1JetParticleCollection>( "RecalibratedTowerJetLPUSak5PUSL1Jet" );
  produces <l1extra::L1JetParticleCollection>( "NVTXRecalibratedTowerJetLPUSak5PUSL1Jet" );


  // Closure test - ak5PUSRaw -> ak5PUS
  produces <l1extra::L1JetParticleCollection>( "Calibratedak5PUSRawak5PUSL1Jet" );


  produces <l1extra::L1JetParticleCollection>( "CalibratedGCTL1Jet" );

  // ****************************************************************************************************
  // *                                 Load configuration parameters                                    *
  // ****************************************************************************************************


  // ************************************************************
  // *                       Recalibration                      *
  // ************************************************************

  int LowPtParameters    = 10;
  int HighPtParameters   = 3;
  int NVTXParameters     = 2;

  // PrePUS
  recalibLowPt_PrePUS_ak5PUS.columns    = LowPtParameters;
  recalibHighPt_PrePUS_ak5PUS.columns   = HighPtParameters;
  recalibNVTX_PrePUS_ak5PUS.columns     = NVTXParameters;
  // PUS
  recalibLowPt_PUS_ak5PUS.columns       = LowPtParameters;
  recalibHighPt_PUS_ak5PUS.columns      = HighPtParameters;
  recalibNVTX_PUS_ak5PUS.columns        = NVTXParameters;
  // LPUS
  recalibLowPt_LPUS_ak5PUS.columns      = LowPtParameters;
  recalibHighPt_LPUS_ak5PUS.columns     = HighPtParameters;
  recalibNVTX_LPUS_ak5PUS.columns       = NVTXParameters;

  // Recalibration - 05/02/14

  // L1 Pt threshold for high Pt LUT
  LUTPtTransition = 0;

  // PrePUS
//  loadLUT( "/vols/cms04/mb1512/Batch/2014-02-05_SingleMu_12Dec_11x11_rev_5/Calibration_CalibPrePUS_ak5PUSLowPt.LUT",  recalibLowPt_PrePUS_ak5PUS.LUTArr );
//  loadLUT( "/vols/cms04/mb1512/Batch/2014-02-05_SingleMu_12Dec_11x11_rev_5/Calibration_CalibPrePUS_ak5PUSHighPt.LUT", recalibHighPt_PrePUS_ak5PUS.LUTArr );
//  loadLUT( "/vols/cms04/mb1512/Batch/2014-02-05_SingleMu_12Dec_11x11_rev_5/Calibration_CalibPrePUS_ak5PUSNVTXResponse.LUT", recalibNVTX_PrePUS_ak5PUS.LUTArr ); 

  // PUS
//  loadLUT( "/vols/cms04/mb1512/Batch/2014-02-07_SingleMu_12Dec_11x11/Calibration_CalibPUS_ak5PUSLowPt.LUT",        recalibLowPt_PUS_ak5PUS.LUTArr );
//  loadLUT( "/vols/cms04/mb1512/Batch/2014-02-07_SingleMu_12Dec_11x11/Calibration_CalibPUS_ak5PUSHighPt.LUT",       recalibHighPt_PUS_ak5PUS.LUTArr );
//  loadLUT( "/vols/cms04/mb1512/Batch/2014-02-07_SingleMu_12Dec_11x11/Calibration_CalibPUS_ak5PUSNVTXResponse.LUT", recalibNVTX_PUS_ak5PUS.LUTArr ); 

  // LPUS
//  loadLUT( "/vols/cms04/mb1512/Batch/2014-02-07_SingleMu_12Dec_11x11/Calibration_CalibLPUS_ak5PUSLowPt.LUT",        recalibLowPt_LPUS_ak5PUS.LUTArr );
//  loadLUT( "/vols/cms04/mb1512/Batch/2014-02-07_SingleMu_12Dec_11x11/Calibration_CalibLPUS_ak5PUSHighPt.LUT",       recalibHighPt_LPUS_ak5PUS.LUTArr );
//  loadLUT( "/vols/cms04/mb1512/Batch/2014-02-07_SingleMu_12Dec_11x11/Calibration_CalibLPUS_ak5PUSNVTXResponse.LUT", recalibNVTX_LPUS_ak5PUS.LUTArr ); 

  // ************************************************************
  // *                        Calibration                       *
  // ************************************************************
  SUBPRINT("Calibration")
    
  pTCalibrationThreshold = iConfig.getParameter< double >("pTCalibrationThreshold");
  iEtaCalibBinWidth = iConfig.getParameter< double >("iEtaCalibrationBinWidth");
    
  
  // Fit order used
  fitOrder = iConfig.getParameter< int >("nParams");                                                                                                       

  // PrePUSak5PUS
  pTCalibration_PrePUS_ak5PUS = iConfig.getParameter< std::vector< double > >("pTCalibration_PrePUS_ak5PUS");            
  //pTCalibration_PrePUS_ak5PUSLt3 = iConfig.getParameter< std::vector< double > >("pTCalibration_PrePUS_ak5PUSLt3");            
  //pTCalibration_PrePUS_NVTXLt15  = iConfig.getParameter< std::vector< double > >("pTCalibration_PrePUS_ak5PUSNVTXLt15");            
  //pTCalibration_PrePUS_NVTXLt25  = iConfig.getParameter< std::vector< double > >("pTCalibration_PrePUS_ak5PUSNVTXLt25");            
  //pTCalibration_PrePUS_NVTXLt50  = iConfig.getParameter< std::vector< double > >("pTCalibration_PrePUS_ak5PUSNVTXLt50");            
  // PUSak5PrePUS 
  pTCalibration_PUS_ak5PUS    = iConfig.getParameter< std::vector< double > >("pTCalibration_PUS_ak5PUS");            
  // PUSak5PUS
  //pTCalibration_LPUS_ak5PUS   = iConfig.getParameter< std::vector< double > >("pTCalibration_LPUS_ak5PUS");            
  // ak5PUSRaw
  //pTCalibration_ak5PUSRaw_ak5PUS = iConfig.getParameter< std::vector< double > >("pTCalibration_ak5PUSRaw_ak5PUS");

  //pTCalibration_UncalibGCT_ak5PUS = iConfig.getParameter< std::vector< double > >("pTCalibration_UncalibGCT_ak5PUS");

  //Invert the fit function?
  invertFunction = iConfig.getParameter< bool > ("invertFunction");

//   iEtaPtCorrectionPrePUSak5PrePUS  = iConfig.getParameter< std::vector< double > >("pTCorrectionPrePUSak5PrePUS");
//   iEtaPtOffsetPrePUSak5PrePUS      = iConfig.getParameter< std::vector< double > >("pTOffsetPrePUSak5PrePUS");
//   // PrePUSak5PUS
//   iEtaPtCorrectionPrePUSak5PUS     = iConfig.getParameter< std::vector< double > >("pTCorrectionPrePUSak5PUS");
//   iEtaPtOffsetPrePUSak5PUS         = iConfig.getParameter< std::vector< double > >("pTOffsetPrePUSak5PUS");
//   // PUSak5PUS
//   iEtaPtCorrectionPUSak5PUS        = iConfig.getParameter< std::vector< double > >("pTCorrectionPUSak5PUS");
//   iEtaPtOffsetPUSak5PUS            = iConfig.getParameter< std::vector< double > >("pTOffsetPUSak5PUS");




  TFileDirectory CalibDir = fs->mkdir( "JetCalibProducer" );


  // Make eta-pT calibration plots

  // PrePUS
  hist2D["pTCalibration_PrePUS_ak5PUS_eta"]    = CalibDir.make<TH2D>("pTCalibration_PrePUS_ak5PUS_eta",
								     "L1 PrePUS jet correction #eta vs L1 p_{T};L1 p_{T} (GeV);#eta", 
								     MAX_L1_PT,0,MAX_L1_PT, TTetaBin,TTetaBins);
  hist2D["pTCalibration_PrePUS_ak5PUS_row"]    = CalibDir.make<TH2D>("pTCalibration_PrePUS_ak5PUS_row",
								     "L1 PrePUS jet correction row vs L1 p_{T};L1 p_{T} (GeV);row", 
								     MAX_L1_PT,0,MAX_L1_PT, 57,-0.5,56.5);
  // PUS
  hist2D["pTCalibration_PUS_ak5PUS_eta"]       = CalibDir.make<TH2D>("pTCalibration_PUS_ak5PUS_eta",
								     "L1 PUS jet correction #eta vs L1 p_{T};L1 p_{T} (GeV);#eta", 
								     MAX_L1_PT,0,MAX_L1_PT, TTetaBin,TTetaBins);
  hist2D["pTCalibration_PUS_ak5PUS_row"]       = CalibDir.make<TH2D>("pTCalibration_PUS_ak5PUS_row",
								     "L1 PUS jet correction row vs L1 p_{T};L1 p_{T} (GeV);row", 
								     MAX_L1_PT,0,MAX_L1_PT, 57,-0.5,56.5);
  // LPUS
  hist2D["pTCalibration_LPUS_ak5PUS_eta"]      = CalibDir.make<TH2D>("pTCalibration_LPUS_ak5PUS_eta",
								     "L1 LPUS jet correction #eta vs L1 p_{T};L1 p_{T} (GeV);#eta", 
								     MAX_L1_PT,0,MAX_L1_PT, TTetaBin,TTetaBins);
  hist2D["pTCalibration_LPUS_ak5PUS_row"]      = CalibDir.make<TH2D>("pTCalibration_LPUS_ak5PUS_row",
								     "L1 LPUS jet correction row vs L1 p_{T};L1 p_{T} (GeV);row", 
								     MAX_L1_PT,0,MAX_L1_PT, 57,-0.5,56.5);
  

  hist2D["pTCalibration_ak5PUSRaw_ak5PUS_eta"] = CalibDir.make<TH2D>("pTCalibration_ak5PUSRaw_ak5PUS_eta",
								     "Ak5Pus Raw correction #eta vs ak5 p_{T};Ak5 p_{T} (GeV);#eta",
								     MAX_L1_PT,0,MAX_L1_PT, TTetaBin,TTetaBins);
  hist2D["pTCalibration_ak5PUSRaw_ak5PUS_row"] = CalibDir.make<TH2D>("pTCalibration_ak5PUSRaw_ak5PUS_row",
								     "Ak5Pus Raw correction row vs L1 p_{T};Ak5 p_{T} (GeV);row", 
								     MAX_L1_PT,0,MAX_L1_PT, 57,-0.5,56.5);


  hist2D["pTCalibration_UncalibGCT_ak5PUS_eta"] = CalibDir.make<TH2D>("pTCalibration_UncalibGCT_ak5PUS_eta",
								     "Ak5Pus Raw correction #eta vs ak5 p_{T};Ak5 p_{T} (GeV);#eta",
								     MAX_L1_PT,0,MAX_L1_PT, TTetaBin,TTetaBins);
  hist2D["pTCalibration_UncalibGCT_ak5PUS_row"] = CalibDir.make<TH2D>("pTCalibration_UncalibGCT_ak5PUS_row",
								     "Uncalibrated GCT jet correction row vs L1 p_{T};L1 p_{T} (GeV);row", 
								     MAX_L1_PT,0,MAX_L1_PT, 57,-0.5,56.5);


  // Make iEta-pT calibration plots
//   hist2D["pTCalibration_PrePUS_ak5PUS"]   = CalibDir.make<TH2D>("pTCalibration_PrePUS_ak5PUS","L1 jet correction i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PrePUS_ak5PUSLt3"]= CalibDir.make<TH2D>("pTCalibration_PrePUS_ak5PUSLt3","L1 jet correction i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PrePUS_NVTXLt15"] = CalibDir.make<TH2D>("pTCalibration_PrePUS_NVTXLt15","L1 jet correction i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PrePUS_NVTXLt25"] = CalibDir.make<TH2D>("pTCalibration_PrePUS_NVTXLt25","L1 jet correction i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PrePUS_NVTXLt50"] = CalibDir.make<TH2D>("pTCalibration_PrePUS_NVTXLt50","L1 jet correction i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PUS_ak5PUS"]      = CalibDir.make<TH2D>("pTCalibration_PUS_ak5PUS","L1 jet correction i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_LPUS_ak5PUS"]     = CalibDir.make<TH2D>("pTCalibration_LPUS_ak5PUS","L1 jet correction i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);

//    hist2D["pTCalibration_PrePUS_ak5PUS_PT"]      = CalibDir.make<TH2D>("pTCalibration_PrePUS_ak5PUS_PT","L1 corrected jet p_{T} i#eta vs L1 p_{T};Corrected L1 p_{T} (GeV);i#eta", 
// 								    MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PrePUS_ak5PUSLt3_PT"] = CalibDir.make<TH2D>("pTCalibration_PrePUS_ak5PUSLt3_PT","L1 corrected jet p_{T} i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PrePUS_NVTXLt15_PT"]   = CalibDir.make<TH2D>("pTCalibration_PrePUS_NVTXLt15_PT","L1 corrected jet p_{T} i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PrePUS_NVTXLt25_PT"]   = CalibDir.make<TH2D>("pTCalibration_PrePUS_NVTXLt25_PT","L1 corrected jet p_{T} i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PrePUS_NVTXLt50_PT"]    = CalibDir.make<TH2D>("pTCalibration_PrePUS_NVTXLt50_PT","L1 corrected jet p_{T} i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_PUS_ak5PUS_PT"]               = CalibDir.make<TH2D>("pTCalibration_PUS_ak5PUS_PT","L1 corrected jet p_{T} i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);
//   hist2D["pTCalibration_LPUS_ak5PUS_PT"]             = CalibDir.make<TH2D>("pTCalibration_LPUS_ak5PUS_PT","L1 corrected jet p_{T} i#eta vs L1 p_{T};L1 p_{T} (GeV);i#eta", 
// 								   MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);


// hist2D["pTCalibration_ak5PUSRaw_ak5PUS"]             = CalibDir.make<TH2D>("pTCalibration_ak5PUSRaw_ak5PUS",
// "Ak5Pus Raw corrected jet p_{T} i#eta vs L1 p_{T};L1 p_{T}(GeV);i#eta",
//									 MAX_L1_PT,0,MAX_L1_PT, 57,-28.5,28.5);


//     printCalibrationFactors( "pTCalibration_ak5PUSRaw_ak5PUS", pTCalibration_ak5PUSRaw_ak5PUS, fitOrder, MAX_L1_PT, true );

     printCalibrationFactors( "pTCalibration_PrePUS_ak5PUS",    pTCalibration_PrePUS_ak5PUS,    fitOrder, MAX_L1_PT, true );
//      // |eta| < 3
//      printCalibrationFactors( "pTCalibration_PrePUS_ak5PUSLt3", pTCalibration_PrePUS_ak5PUSLt3, fitOrder, MAX_L1_PT, true );
//      // NVTX < 15 
//      printCalibrationFactors( "pTCalibration_PrePUS_NVTXLt15",  pTCalibration_PrePUS_NVTXLt15,  fitOrder, MAX_L1_PT, true );
//      // 15 < NVTX < 25 
//      printCalibrationFactors( "pTCalibration_PrePUS_NVTXLt25",  pTCalibration_PrePUS_NVTXLt25,  fitOrder, MAX_L1_PT, true );
//      // 25 < NVTX < 50 
//      printCalibrationFactors( "pTCalibration_PrePUS_NVTXLt50",  pTCalibration_PrePUS_NVTXLt50,  fitOrder, MAX_L1_PT, true );
     // PUS 
//     printCalibrationFactors( "pTCalibration_PUS_ak5PUS",       pTCalibration_PUS_ak5PUS,       fitOrder, MAX_L1_PT, true );
     // LPUS 
//     printCalibrationFactors( "pTCalibration_LPUS_ak5PUS",      pTCalibration_LPUS_ak5PUS,      fitOrder, MAX_L1_PT, true );



     // Uncalib GCT
//     printCalibrationFactors( "pTCalibration_UncalibGCT_ak5PUS",      pTCalibration_UncalibGCT_ak5PUS,      fitOrder, MAX_L1_PT, true );

  
}


JetCalibProducer::~JetCalibProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
JetCalibProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;


   bool evValid = true;


//    std::auto_ptr<l1slhc::L1TowerJetCollection>     outputCalibTowerJet_Tower( new l1slhc::L1TowerJetCollection() );
//    std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJet_L1Jet( new l1extra::L1JetParticleCollection() );

//    // PrePUSak5PrePUS
//    std::auto_ptr<l1slhc::L1TowerJetCollection>     outputCalibTowerJetPrePUSak5PrePUS_Tower( new l1slhc::L1TowerJetCollection() );
//    std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPrePUSak5PrePUS_L1Jet( new l1extra::L1JetParticleCollection() );
   // PrePUSak5PUS   
   std::auto_ptr<l1slhc::L1TowerJetCollection>     outputCalibTowerJetPrePUSak5PUS_Tower( new l1slhc::L1TowerJetCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPrePUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputRecalibTowerJetPrePUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputNVTXRecalibTowerJetPrePUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );

   //   std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPrePUSak5PUSLt3_L1Jet( new l1extra::L1JetParticleCollection() );
//    std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPrePUSak5PUSNVTXLt15_L1Jet( new l1extra::L1JetParticleCollection() );
//    std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPrePUSak5PUSNVTXLt25_L1Jet( new l1extra::L1JetParticleCollection() );
//    std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPrePUSak5PUSNVTXLt50_L1Jet( new l1extra::L1JetParticleCollection() );
   // PUSak5PUS 
   std::auto_ptr<l1slhc::L1TowerJetCollection>     outputCalibTowerJetPUSak5PUS_Tower( new l1slhc::L1TowerJetCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputRecalibTowerJetPUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputNVTXRecalibTowerJetPUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );


   // LPUSak5PUS 
   std::auto_ptr<l1slhc::L1TowerJetCollection>     outputCalibTowerJetLPUSak5PUS_Tower( new l1slhc::L1TowerJetCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetLPUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputRecalibTowerJetLPUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputNVTXRecalibTowerJetLPUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );

   // ak5PUSRaw
   std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibAk5PUSRawak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );

   // uncaibGCT
   std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibGCT_L1Jet( new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputRecalibGCT_L1Jet( new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputNVTXRecalibGCT_L1Jet( new l1extra::L1JetParticleCollection() );


   // ****************************************************************************************************
   // *                                             Handles                                              *
   // ****************************************************************************************************
 
   //Need this for information about PU
   edm::Handle<int> vertices;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("NumPrimaryVertices"), vertices);
   if(!vertices.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("NumPrimaryVertices") << std::endl;
     evValid = false;
   }


   // ************************************************************ 
   // *                        RECO Jets                         * 
   // ************************************************************ 

   //PU subtracted AK5 calo jets-must be in root file read in                                                                                                  
//    edm::Handle<reco::CaloJetCollection> PUSRawAk5Jets;
//    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUSPreCalibCaloJets"), PUSRawAk5Jets);
//    if(!PUSRawAk5Jets.isValid()){
//      edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUSPreCalibCaloJets") << std::endl;
//      evValid = false;
//    }

//   edm::Handle<l1extra::L1JetParticleCollection> Ak5PUSRaw_L1Jet;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUSRawAk5CaloJetL1Jet"), Ak5PUSRaw_L1Jet);
//   if(!Ak5PUSRaw_L1Jet.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUSRawAk5CaloJetL1Jet") << std::endl;
//     evValid = false;
//   }


   // ************************************************************
   // *                          L1Jets                          *
   // ************************************************************


   edm::Handle<l1extra::L1JetParticleCollection> UncalibJetPrePUS_L1Jet;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UncalibratedPrePUSL1Jet"), UncalibJetPrePUS_L1Jet);
   if(!UncalibJetPrePUS_L1Jet.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("UncalibratedPrePUSL1Jet") << std::endl;
     evValid = false;
   }

   edm::Handle<l1extra::L1JetParticleCollection> UncalibJetPUS_L1Jet;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UncalibratedPUSL1Jet"), UncalibJetPUS_L1Jet);
   if(!UncalibJetPUS_L1Jet.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("UncalibratedPUSL1Jet") << std::endl;
     evValid = false;
   }

//   edm::Handle<l1extra::L1JetParticleCollection> UncalibJetLPUS_L1Jet;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UncalibratedLPUSL1Jet"), UncalibJetLPUS_L1Jet);
//   if(!UncalibJetLPUS_L1Jet.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("UncalibratedLPUSL1Jet") << std::endl;
//     evValid = false;
//   }




//   edm::Handle<l1extra::L1JetParticleCollection> UncalibJetGCT_L1Jet;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UncalibratedGCTL1Jet"), UncalibJetGCT_L1Jet);
//   if(!UncalibJetGCT_L1Jet.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("UncalibratedGCTL1Jet") << std::endl;
//     evValid = false;
//   }


   // ************************************************************
   // *                        Tower Jets                        *
   // ************************************************************
   SUBPRINT("Tower Jets")

//    // Uncalibrated TowerJets
//    edm::Handle<l1slhc::L1TowerJetCollection> UncalibJetPrePUS_Tower;
//    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UncalibratedPrePUSTowerJet"), UncalibJetPrePUS_Tower);
//    if(!UncalibJetPrePUS_Tower.isValid()){
//      evValid = false;
//      edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("UncalibratedPrePUSTowerJet") << std::endl;
//    }

//    edm::Handle<l1slhc::L1TowerJetCollection> UncalibJetPUS_Tower;
//    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UncalibratedPUSTowerJet"), UncalibJetPUS_Tower);
//    if(!UncalibJetPUS_Tower.isValid()){
//      evValid = false;
//      edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("UncalibratedPUSTowerJet") << std::endl;
//    }


   //   fitOrder = conf_.getParameter< int >("nParams");


  // PrePUSak5PrePUS
   //   pTCalibration_PUS_ak5PUS = conf_.getParameter< std::vector< double > >("pTCalibration_PUS_ak5PUS");

//   iEtaPtCorrectionPrePUSak5PrePUS  = conf_.getParameter< std::vector< double > >("pTCorrectionPrePUSak5PrePUS");
//   iEtaPtOffsetPrePUSak5PrePUS      = conf_.getParameter< std::vector< double > >("pTOffsetPrePUSak5PrePUS");
//   // PrePUSak5PUS
//   iEtaPtCorrectionPrePUSak5PUS     = conf_.getParameter< std::vector< double > >("pTCorrectionPrePUSak5PUS");
//   iEtaPtOffsetPrePUSak5PUS         = conf_.getParameter< std::vector< double > >("pTOffsetPrePUSak5PUS");
//   // PUSak5PUS
//   iEtaPtCorrectionPUSak5PUS        = conf_.getParameter< std::vector< double > >("pTCorrectionPUSak5PUS");
//   iEtaPtOffsetPUSak5PUS            = conf_.getParameter< std::vector< double > >("pTOffsetPUSak5PUS");


   if( evValid ){

     NVTX       = *vertices;

  
     // ****************************************************************************************************
     // *                                         Calibrate Jets                                           *
     // ****************************************************************************************************
     // Calibrate jets with LUT

//      for (l1slhc::L1TowerJetCollection::const_iterator Uncalib_It = UncalibJet_Tower->begin(); Uncalib_It != UncalibJet_Tower->end(); ++Uncalib_It ){

//        l1slhc::L1TowerJet calibJet = (*Uncalib_It);
//        int iPhi                    = calibJet.iPhi();
//        int iEta                    = calibJet.iEta();


//        // Jet pT threshold for calibration, only calibrate and store jets above threshold
//        if ( calibJet.Pt() < pTCalibrationThreshold ){
// 	 continue;
//        }
       

//        // Extract iEta dependent correction factors
//        // **************************************************
//        int iEtaIndex       = iEta + 28;
//        if (iEta > 0)  // Correct for missing iEta = 0
// 	 iEtaIndex--;
	  
//        double unCorrectedPt = calibJet.Pt(); 
//        double pTOffset      = iEtaPtOffset[iEtaIndex];
//        double pTCorrection  = iEtaPtCorrection[iEtaIndex];
//        double correctedPt   = unCorrectedPt*pTCorrection + pTOffset;//(unCorrectedPt - pTOffset)/pTCorrection;

//        // Create calibrated TowerJet
//        calibJet.setPt( correctedPt );
//        outputCalibTowerJet_Tower->insert( iEta, iPhi, calibJet );

//        // Create calibrated L1Jet
//        math::PtEtaPhiMLorentzVector tempJet;
//        tempJet.SetCoordinates( calibJet.p4().Pt(), calibJet.p4().eta(), calibJet.p4().phi(), calibJet.p4().M() );
//        outputCalibTowerJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );

//      }




     
//      // PrePUSak5PrePUS
//      calibrateJet( UncalibJetPrePUS_Tower, iEtaPtOffsetPrePUSak5PrePUS, iEtaPtCorrectionPrePUSak5PrePUS,
// 		   outputCalibTowerJetPrePUSak5PrePUS_Tower,
// 		   outputCalibTowerJetPrePUSak5PrePUS_L1Jet);


//      // PrePUSak5PUS
//      calibrateJet( UncalibJetPrePUS_Tower, iEtaPtOffsetPrePUSak5PUS, iEtaPtCorrectionPrePUSak5PUS,
// 		   outputCalibTowerJetPrePUSak5PUS_Tower,
// 		   outputCalibTowerJetPrePUSak5PUS_L1Jet);

     // ****************************************
     // *                 ak5PUS               *
     // ****************************************
//      calibrateJet( UncalibJetPUS_Tower, iEtaPtOffsetPUSak5PUS, iEtaPtCorrectionPUSak5PUS,
// 		   outputCalibTowerJetPUSak5PUS_Tower,
// 		   outputCalibTowerJetPUSak5PUS_L1Jet);

     // PrePUS
     calibrateJet( UncalibJetPrePUS_L1Jet, pTCalibration_PrePUS_ak5PUS, fitOrder,
		   outputCalibTowerJetPrePUSak5PUS_L1Jet, 
  	   outputRecalibTowerJetPrePUSak5PUS_L1Jet, 
		   outputNVTXRecalibTowerJetPrePUSak5PUS_L1Jet, 
		   true,
		   emptyLUT,//recalibLowPt_PrePUS_ak5PUS,
 		   emptyLUT,//recalibHighPt_PrePUS_ak5PUS,
 		   emptyLUT,//recalibNVTX_PrePUS_ak5PUS,
		   LUTPtTransition	   );


     // PUS
     calibrateJet( UncalibJetPUS_L1Jet, pTCalibration_PUS_ak5PUS, fitOrder,
		   outputCalibTowerJetPUSak5PUS_L1Jet, 
		   outputRecalibTowerJetPUSak5PUS_L1Jet, 
		   outputNVTXRecalibTowerJetPUSak5PUS_L1Jet, 
		   true,
		   emptyLUT,//recalibLowPt_PUS_ak5PUS,
		   emptyLUT,//recalibHighPt_PUS_ak5PUS,
		   emptyLUT,
		   LUTPtTransition	   );


     // LPUS
//     calibrateJet( UncalibJetLPUS_L1Jet, pTCalibration_LPUS_ak5PUS, fitOrder,
//		   outputCalibTowerJetLPUSak5PUS_L1Jet, 
//		   outputRecalibTowerJetLPUSak5PUS_L1Jet, 
//		   outputNVTXRecalibTowerJetLPUSak5PUS_L1Jet, 
//		   true,
//		   recalibLowPt_LPUS_ak5PUS,
//		   recalibHighPt_LPUS_ak5PUS,
//		   emptyLUT,
//		   LUTPtTransition	   );



     // Uncalib GCT
//     calibrateJet( UncalibJetGCT_L1Jet, pTCalibration_UncalibGCT_ak5PUS, fitOrder,
//                   outputCalibGCT_L1Jet,
//                   outputRecalibGCT_L1Jet,
//                   outputNVTXRecalibGCT_L1Jet,
//                   true,
//                   emptyLUT,
//                   emptyLUT,
//                   emptyLUT,
//                   LUTPtTransition );





//      calibrateJet( UncalibJetPrePUS_L1Jet, pTCalibration_PrePUS_ak5PUS, fitOrder,
// 		   outputCalibTowerJetPrePUSak5PUS_Tower,
//                    outputCalibTowerJetPrePUSak5PUS_L1Jet, true);




     
//      calibrateJet( UncalibJetPrePUS_Tower, pTCalibration_PrePUS_ak5PUSLt3, fitOrder,
// 		   outputCalibTowerJetPrePUSak5PUS_Tower,
//                    outputCalibTowerJetPrePUSak5PUSLt3_L1Jet, true);

     // NVTX binned
     // ********************************************************************************


//      // NVTX < 15
//      if (NVTX < 15){
//        calibrateJet( UncalibJetPUS_Tower, pTCalibration_PrePUS_NVTXLt15, fitOrder,
// 		     outputCalibTowerJetPrePUSak5PUS_Tower,
// 		     outputCalibTowerJetPrePUSak5PUSNVTXLt15_L1Jet, true);
//      }
//      // 15 < NVTX < 25
//      else if (NVTX < 25){
//        calibrateJet( UncalibJetPUS_Tower, pTCalibration_PrePUS_NVTXLt25, fitOrder,
// 		     outputCalibTowerJetPrePUSak5PUS_Tower,
// 		     outputCalibTowerJetPrePUSak5PUSNVTXLt25_L1Jet, true);
//      }
//      // 25 < NVTX < 50
// //      else if (NVTX < 50){
// //        calibrateJet( UncalibJetPUS_Tower, pTCalibration_PrePUS_NVTXLt50, fitOrder,
// // 		     outputCalibTowerJetPrePUSak5PUS_Tower,
// // 		     outputCalibTowerJetPrePUSak5PUSNVTXLt50_L1Jet, true);
// //      }

//      // PUS
//      calibrateJet( UncalibJetPUS_Tower, pTCalibration_PUS_ak5PUS, fitOrder,
// 		   outputCalibTowerJetPUSak5PUS_Tower,
//                    outputCalibTowerJetPUSak5PUS_L1Jet, true);
//      // LPUS
//      calibrateJet( UncalibJetPUS_Tower, pTCalibration_LPUS_ak5PUS, fitOrder,
// 		   outputCalibTowerJetLPUSak5PUS_Tower,
//                    outputCalibTowerJetLPUSak5PUS_L1Jet, true);



     // Create dummy collection to pass to the calibration function


//       calibrateJet( Ak5PUSRaw_L1Jet, pTCalibration_ak5PUSRaw_ak5PUS, fitOrder,
// 		    outputCalibAk5PUSRawak5PUS_L1Jet, true,
// 	 	    emptyLUT,
// 		    emptyLUT,
// 		    LUTPtTransition );
		    


     



     // ak5PUSRaw
//      calibrateJet( ak5PUSRaw, pTCalibration_ak5PUSRaw_ak5PUS, fitOrder,
//                    outputCalibTowerJetLPUSak5PUS_L1Jet, true);




     // Store Jets
//      iEvent.put( outputCalibTowerJet_Tower, "CalibratedTowerJetTower" );
//      iEvent.put( outputCalibTowerJet_L1Jet, "CalibratedTowerJetL1Jet");

     // PrePUSak5PrePUS
//      iEvent.put( outputCalibTowerJetPrePUSak5PrePUS_Tower, "CalibratedTowerJetPrePUSak5PrePUSTower" );
//      iEvent.put( outputCalibTowerJetPrePUSak5PrePUS_L1Jet, "CalibratedTowerJetPrePUSak5PrePUSL1Jet");
     // PrePUSak5PUS
     iEvent.put( outputCalibTowerJetPrePUSak5PUS_Tower, "CalibratedTowerJetPrePUSak5PUSTower" );
     iEvent.put( outputCalibTowerJetPrePUSak5PUS_L1Jet, "CalibratedTowerJetPrePUSak5PUSL1Jet");
     iEvent.put( outputRecalibTowerJetPrePUSak5PUS_L1Jet, "RecalibratedTowerJetPrePUSak5PUSL1Jet");
     iEvent.put( outputNVTXRecalibTowerJetPrePUSak5PUS_L1Jet, "NVTXRecalibratedTowerJetPrePUSak5PUSL1Jet");
     // PrePUSak5PUSLt3
     //     iEvent.put( outputCalibTowerJetPrePUSak5PUSLt3_L1Jet, "CalibratedTowerJetPrePUSak5PUSLt3L1Jet");
//      iEvent.put( outputCalibTowerJetPrePUSak5PUSNVTXLt15_L1Jet, "CalibratedTowerJetPrePUSak5PUSNVTXLt15L1Jet");
//      iEvent.put( outputCalibTowerJetPrePUSak5PUSNVTXLt25_L1Jet, "CalibratedTowerJetPrePUSak5PUSNVTXLt25L1Jet");
//      iEvent.put( outputCalibTowerJetPrePUSak5PUSNVTXLt50_L1Jet, "CalibratedTowerJetPrePUSak5PUSNVTXLt50L1Jet");
     // PUSak5PUS
     //     iEvent.put( outputCalibTowerJetPUSak5PUS_Tower, "CalibratedTowerJetPUSak5PUSTower" );
     iEvent.put( outputCalibTowerJetPUSak5PUS_L1Jet,       "CalibratedTowerJetPUSak5PUSL1Jet");
     iEvent.put( outputRecalibTowerJetPUSak5PUS_L1Jet,     "RecalibratedTowerJetPUSak5PUSL1Jet");
     iEvent.put( outputNVTXRecalibTowerJetPUSak5PUS_L1Jet, "NVTXRecalibratedTowerJetPUSak5PUSL1Jet");

     // LPUSak5PUS
     //     iEvent.put( outputCalibTowerJetLPUSak5PUS_Tower, "CalibratedTowerJetLPUSak5PUSTower" );
     iEvent.put( outputCalibTowerJetLPUSak5PUS_L1Jet,       "CalibratedTowerJetLPUSak5PUSL1Jet");
     iEvent.put( outputRecalibTowerJetLPUSak5PUS_L1Jet,     "RecalibratedTowerJetLPUSak5PUSL1Jet");
     iEvent.put( outputNVTXRecalibTowerJetLPUSak5PUS_L1Jet, "NVTXRecalibratedTowerJetLPUSak5PUSL1Jet");

     // ak5PUSRaw
     iEvent.put( outputCalibAk5PUSRawak5PUS_L1Jet, "Calibratedak5PUSRawak5PUSL1Jet" );

     // Calib GCT
     iEvent.put( outputCalibGCT_L1Jet, "CalibratedGCTL1Jet" );


   }

}

// ------------ method called once each job just before starting event loop  ------------
void 
JetCalibProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetCalibProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
JetCalibProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetCalibProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetCalibProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetCalibProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetCalibProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}




// Jet calibration
// void 
// JetCalibProducer::calibrateJet( edm::Handle<l1slhc::L1TowerJetCollection> const& UncalibJet_Tower, std::vector< double > const& iEtaPtOffset, 
// 				std::vector< double > const& iEtaPtCorrection, 
// 				std::auto_ptr<l1slhc::L1TowerJetCollection>     &outputCalibTowerJet_Tower,
// 				std::auto_ptr<l1extra::L1JetParticleCollection> &outputCalibTowerJet_L1Jet ){
void 
JetCalibProducer::calibrateJet( edm::Handle<l1slhc::L1TowerJetCollection> const& UncalibJet_Tower, 
				std::vector< double > const& iEtaPtCorrections,
				const int&  fitOrder, 
				std::auto_ptr<l1slhc::L1TowerJetCollection>     &outputCalibTowerJet_Tower,
				std::auto_ptr<l1extra::L1JetParticleCollection> &outputCalibTowerJet_L1Jet,
				bool newMethod ){



     // ****************************************************************************************************
     // *                                         Calibrate Jets                                           *
     // ****************************************************************************************************


  //     l1slhc::L1TowerJetCollection unsortedCalibratedTowerJets;
  std::vector <l1slhc::L1TowerJet> unsortedCalibratedTowerJets;

     // Calibrate jets with LUT
     for (l1slhc::L1TowerJetCollection::const_iterator Uncalib_It = UncalibJet_Tower->begin(); Uncalib_It != UncalibJet_Tower->end(); ++Uncalib_It ){

       l1slhc::L1TowerJet calibJet = (*Uncalib_It);
       //       int iPhi                    = calibJet.iPhi();
       int iEta                    = calibJet.iEta();

 
       // Jet pT threshold for calibration, only calibrate above threshold
       if ( calibJet.Pt() < pTCalibrationThreshold ){

	 // Create un-calibrated TowerJet 
	 //unsortedCalibratedTowerJets.push_back( (*Uncalib_It)  );

	 continue;
       }
       

       // Extract iEta dependent correction factors
       // **************************************************
       int iEtaIndex       = iEta + 28;
       if (iEta > 0)  // Correct for missing iEta = 0
	 iEtaIndex--;
       
       // Index in the LUT to search for the correction parameter
       int correctionIndex  = fitOrder*iEtaIndex;




       // Perform jet calibration
       // ------------------------------------------------------------
       double unCorrectedPt  = calibJet.Pt(); 
       double correctedPt    = 0;

    
       // Use the old method
       if (newMethod == false ){
	 // pT to the power of the required order
	 double pTOrder = unCorrectedPt;

	 //       std::cout << "iEta = " << iEta << "\nInitial pT = " << unCorrectedPt << "\n";
	 for (int iOrder = 0; iOrder < fitOrder; ++iOrder){
	 
	   uint index = correctionIndex + iOrder;

	   // Get correction for current order of pT
	   double correctionFactor = iEtaPtCorrections[ index ];
	   double correction       = 0;

	   if (iOrder == 0){
	     correction += correctionFactor;
	   }
	   else{
	     correction += correctionFactor*pTOrder;
	     // Generate pT to the power of the next order
	     pTOrder *= unCorrectedPt;
	   }
	 
	   // Calculate the corrected jet pT
	   correctedPt += correction;

	   //  	 std::cout << "CorrectionFactor = " << correctionFactor << "\tCorrection = " << correction 
	   // 		   << "\tCorrectedPt = " << correctedPt << "\n";

	 }
	 //       std::cout << "\n";
       }
       else{ // Use the new method

	 if ( fitOrder != 6 ){
	   edm::LogWarning("Calibration error") << "ERROR: Cannot perform the new fit method, require 6 parameters only " << fitOrder << " specfied.\n";
	 }
	 
	 // Get parameters
	 double p0 = iEtaPtCorrections[ correctionIndex ];
	 double p1 = iEtaPtCorrections[ correctionIndex + 1 ];
	 double p2 = iEtaPtCorrections[ correctionIndex + 2 ];
	 double p3 = iEtaPtCorrections[ correctionIndex + 3 ];
	 double p4 = iEtaPtCorrections[ correctionIndex + 4 ];
	 double p5 = iEtaPtCorrections[ correctionIndex + 5 ];

	 // JETMET uses log base 10
	 //	 double logPt = log10( unCorrectedPt );
	 	 double logPt = log( unCorrectedPt );

	 double term1 = p1 / ( logPt * logPt + p2 );
	 double term2 = p3 * exp( -p4*((logPt - p5)*(logPt - p5)) );

	 // Calculate the corrected Pt 
	 double correction = (p0 + term1 + term2);
   if(invertFunction) correction = 1.0/correction;
	 correctedPt    = correction*unCorrectedPt;
	 

// 	 std::cout << "\tp0 = " << p0 << "\n"
// 		   << "\tp1 = " << p1 << "\n"
// 		   << "\tp2 = " << p2 << "\n"
// 		   << "\tp3 = " << p3 << "\n"
// 		   << "\tp4 = " << p4 << "\n"
// 		   << "\tp5 = " << p5 << "\n\n";

// 	 std::cout << "iEta = " << iEta << "\tiPhi = " << iPhi 
// 		   << "\nl1Pt = " << unCorrectedPt << "\tCorrection = " << correction 
// 		   << "\tTerm 1 = " << term1 << "\tTerm 2 = " << term2 << "\tCorrectedPt = " << correctedPt << "\n\n\n";



       }

//        double pTOffset      = iEtaPtOffset[iEtaIndex];
//        double pTCorrection  = iEtaPtCorrection[iEtaIndex];
//         double correctedPt   = unCorrectedPt*pTCorrection + pTOffset;//(unCorrectedPt - pTOffset)/pTCorrection;

       // Create and store unsorted calibrated TowerJet
       calibJet.setPt( correctedPt );
       unsortedCalibratedTowerJets.push_back( calibJet );
       //       unsortedCalibratedTowerJets->insert( iEta, iPhilTowerEtaPhi.second, lJet );

     }

     // Sort calibrated jets
     std::sort( unsortedCalibratedTowerJets.begin(), unsortedCalibratedTowerJets.end(), towerJetRankDescending );

     // Store sorted and calibrated TowerJets
     for (l1slhc::L1TowerJetCollection::const_iterator Calib_It = unsortedCalibratedTowerJets.begin(); Calib_It != unsortedCalibratedTowerJets.end(); ++Calib_It ){

       l1slhc::L1TowerJet calibJet = (*Calib_It);
       int iPhi                    = calibJet.iPhi();
       int iEta                    = calibJet.iEta();
      
       outputCalibTowerJet_Tower->insert( iEta, iPhi, calibJet );

       // Create calibrated L1Jet
       math::PtEtaPhiMLorentzVector tempJet;
       tempJet.SetCoordinates( calibJet.p4().Pt(), calibJet.p4().eta(), calibJet.p4().phi(), calibJet.p4().M() );
       outputCalibTowerJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );

     }


}

// Dump the eta, l1 pT corrections
void 
JetCalibProducer::printCalibrationFactors( TString LUT, std::vector< double > const& iEtaPtCorrections, const int& fitOrder, int maxL1Pt, bool newMethod ){


     // ****************************************************************************************************
     // *                                         Calibrate Jets                                           *
     // ****************************************************************************************************

 
  int rows = iEtaPtCorrections.size()/fitOrder;
  // Calculate the TT grouping used in the calibration
  int groupedTTs = 56/rows;

  std::cout << LUT << "\t" << rows << "\t" << groupedTTs << "\n";

  for (int iRow = 0; iRow < rows; ++iRow){


      //  for ( int iEtaIndex = 0; iEtaIndex < 56; ++iEtaIndex ){

       // Extract iEta dependent correction factors
       // **************************************************
       
       // Index in the LUT to search for the correction parameter
       int correctionIndex  = fitOrder*iRow;


       for (int unCorrectedPt = 0; unCorrectedPt < maxL1Pt; ++unCorrectedPt){ 
   
	 // Use the old method
	 if (newMethod == false ){
// 	   // pT to the power of the required order
// 	   double pTOrder = unCorrectedPt;

// 	   //       std::cout << "iEta = " << iEta << "\nInitial pT = " << unCorrectedPt << "\n";
// 	   for (int iOrder = 0; iOrder < fitOrder; ++iOrder){
	     
// 	     uint index = correctionIndex + iOrder;
	     
// 	     // Get correction for current order of pT
// 	     double correctionFactor = iEtaPtCorrections[ index ];
// 	     double correction       = 0;
	     
// 	     if (iOrder == 0){
// 	       correction += correctionFactor;
// 	     }
// 	     else{
// 	       correction += correctionFactor*pTOrder;
// 	       // Generate pT to the power of the next order
// 	       pTOrder *= unCorrectedPt;
// 	     }
	     
// 	     // Calculate the corrected jet pT
// 	     correctedPt += correction;
	     
// 	     //  	 std::cout << "CorrectionFactor = " << correctionFactor << "\tCorrection = " << correction 
// 	     // 		   << "\tCorrectedPt = " << correctedPt << "\n";

// 	   }

	 }
	 else{ // Use the new method
	   
	   if ( fitOrder != 6 ){
	     edm::LogWarning("Calibration error") << "ERROR: Cannot perform the new fit method, require 6 parameters only " << fitOrder << " specfied.\n";
	   }
	 
	   // Get parameters
// 	   double p0 = PtCorrectionLUT.getElement( iRow, 0 );
//  	   double p1 = PtCorrectionLUT.getElement( iRow, 1 );
//  	   double p2 = PtCorrectionLUT.getElement( iRow, 2 );
//  	   double p3 = PtCorrectionLUT.getElement( iRow, 3 );
//  	   double p4 = PtCorrectionLUT.getElement( iRow, 4 );
//  	   double p5 = PtCorrectionLUT.getElement( iRow, 5 );
	   double p0 = iEtaPtCorrections[ correctionIndex ];
	   double p1 = iEtaPtCorrections[ correctionIndex + 1 ];
	   double p2 = iEtaPtCorrections[ correctionIndex + 2 ];
	   double p3 = iEtaPtCorrections[ correctionIndex + 3 ];
	   double p4 = iEtaPtCorrections[ correctionIndex + 4 ];
	   double p5 = iEtaPtCorrections[ correctionIndex + 5 ];
	   
	   
	   double logPt = log( unCorrectedPt );
	   
	   double term1 = p1 / ( logPt * logPt + p2 );
	   double term2 = p3 * exp( -p4*((logPt - p5)*(logPt - p5)) );

	   // Calculate the corrected Pt 
	   double correction = (p0 + term1 + term2);


	   // Print out LUT table row each calibration factor corresponds to for easy debugging
	   hist2D[ LUT + "_row" ]        ->Fill( unCorrectedPt, iRow, correction );

	   // Calculate the respective eta bins
	   // ----------------------------------------

	   // Fill grouped eta bins with the same data
	   for (int iTT = 0; iTT < groupedTTs; ++iTT){
	     // Determine which eta the bin corresponds to, add a little extra to ensure it is inside the bin
	     double eta = TTetaBins[ iRow*groupedTTs + iTT ] + 0.001;
	     // Fill correction plots
	     hist2D[ LUT + "_eta" ]        ->Fill( unCorrectedPt, eta, correction );


	     //	     std::cout << eta << "\t" << unCorrectedPt << "\t" << correction << "\n";
	   }

// 	   int iEta = iEtaIndex -28;
// 	   if (iEta >= 0){ // Correct for no iEta = 0
// 	     iEta++;
// 	   }

	   //	   std::cout << "iEta = " << iEta << "\tCorrection = " << correction << "\n";
	   
	   // Fill correction plots
	   //	   hist2D[ LUT ]        ->Fill( unCorrectedPt, iEta, correction );
	   //	   hist2D[ LUT + "_PT" ]->Fill( unCorrectedPt, iEta, correction*unCorrectedPt );
	   
	 }
	 
       }
       
       

  }
}

















void 
JetCalibProducer::calibrateJet( edm::Handle<l1extra::L1JetParticleCollection> UncalibJet_L1Jet,
				std::vector< double > const& iEtaPtCorrections,
				const int&  fitOrder, 
				std::auto_ptr<l1extra::L1JetParticleCollection> &outputCalibJet_L1Jet,
				std::auto_ptr<l1extra::L1JetParticleCollection> &outputRecalibJet_L1Jet,
				std::auto_ptr<l1extra::L1JetParticleCollection> &outputNVTXRecalibJet_L1Jet,
				bool newMethod,
				LUT lowPtLUT,
				LUT highPtLUT,
				LUT NVTXLUT,
				double transitionPt){
  //				int NVTX){


     // ****************************************************************************************************
     // *                                         Calibrate Jets                                           *
     // ****************************************************************************************************

  std::vector <l1extra::L1JetParticle> unsortedCalibratedL1Jets, unsortedRecalibratedL1Jets, unsortedNVTXRecalibratedL1Jets;
     
     // Calibrate jets with LUT
     for (l1extra::L1JetParticleCollection::const_iterator Uncalib_It = UncalibJet_L1Jet->begin(); Uncalib_It != UncalibJet_L1Jet->end(); ++Uncalib_It ){
    
       l1extra::L1JetParticle uncalibJet = (*Uncalib_It);
 
       // Jet pT threshold for calibration, only calibrate above threshold
       if ( uncalibJet.pt() < pTCalibrationThreshold ){
	 // Store un-calibrated L1Jet 
	 	 unsortedCalibratedL1Jets.push_back( uncalibJet );
	 continue;
       }
       
       // ***************************************************************** 
       // *                           Eta Binned                          * 
       // ***************************************************************** 
       // Find which in which Eta bin the jet resides 
       uint pEta;
       for (pEta = 1; pEta < etaRegionSlice.size(); ++pEta){

	 double eta = uncalibJet.eta();

	 // Get Eta bin lower and upper bounds 
	 double EtaLow  = etaRegionSlice[ pEta - 1];
   double EtaHigh = etaRegionSlice[ pEta ];
   // Eta resides within current boundary
   if ( (eta >= EtaLow) && (eta < EtaHigh) ){
     //	   std::cout <<  pEta << "\tEta = " << eta << "\t[" << EtaLow << ", " << EtaHigh << "]\n";
     break; // found the correct eta bin, break
   }

       }

       // Extract eta dependent correction factors
       // **************************************************
       int etaIndex       = pEta - 1;// Correct for array starting at zero

       // Index in the LUT to search for the correction parameter
       int correctionIndex  = fitOrder*etaIndex;


       // Perform jet calibration
       // ------------------------------------------------------------
       double unCorrectedPt  = uncalibJet.pt(); 
       double correctedPt              = 0;
       double recalibCorrectedPt       = 0;
       double recalibNVTXCorrectedPt = 0;

       // Use the old method
       if (newMethod == false ){
         // pT to the power of the required order
         double pTOrder = unCorrectedPt;

         //       std::cout << "iEta = " << iEta << "\nInitial pT = " << unCorrectedPt << "\n";
         for (int iOrder = 0; iOrder < fitOrder; ++iOrder){

           uint index = correctionIndex + iOrder;

           // Get correction for current order of pT
           double correctionFactor = iEtaPtCorrections[ index ];
           double correction       = 0;

           if (iOrder == 0){
             correction += correctionFactor;
           }
           else{
             correction += correctionFactor*pTOrder;
             // Generate pT to the power of the next order
             pTOrder *= unCorrectedPt;
           }

           // Calculate the corrected jet pT
           correctedPt += correction;

           //    	 std::cout << "CorrectionFactor = " << correctionFactor << "\tCorrection = " << correction 
           // 	   		   << "\tCorrectedPt = " << correctedPt << "\n";

         }
         //	        std::cout << "\n";
       }
       else{ // Use the new method

         if ( fitOrder != 6 ){
           edm::LogWarning("Calibration error") << "ERROR: Cannot perform the new fit method, require 6 parameters only " << fitOrder << " specfied.\n";
         }

         // Get parameters
         double p0 = iEtaPtCorrections[ correctionIndex ];
         double p1 = iEtaPtCorrections[ correctionIndex + 1 ];
         double p2 = iEtaPtCorrections[ correctionIndex + 2 ];
         double p3 = iEtaPtCorrections[ correctionIndex + 3 ];
         double p4 = iEtaPtCorrections[ correctionIndex + 4 ];
         double p5 = iEtaPtCorrections[ correctionIndex + 5 ];



         double logPt = log( unCorrectedPt );

         double term1 = p1 / ( logPt * logPt + p2 );
         double term2 = p3 * exp( -p4*((logPt - p5)*(logPt - p5)) );

         // Calculate the corrected Pt 
         double correction = (p0 + term1 + term2);
         if(invertFunction) correction = 1.0/correction;
         correctedPt    = correction*unCorrectedPt;

       }

       // ----------------------------------------
       // If applicable perform second calibration
       if ( lowPtLUT.columns != 0){

         double recalibFactor;

         if ( correctedPt < transitionPt ){ // Jet Pt is in low-Pt regime

           recalibFactor = getRecalibFactor( correctedPt, etaIndex, lowPtLUT );

         }
         else{ // Pt is in high-Pt regime

           recalibFactor = getRecalibFactor( correctedPt, etaIndex, highPtLUT );

         }

         // 	 std::cout << "CorrectionFactor = " << recalibFactor << "\tPt = " << correctedPt 
         // 		   << "\tCorrectedPt = " << recalibFactor*correctedPt << "\n";


         // Apply recalibration correction to jet Pt
         recalibCorrectedPt = recalibFactor*correctedPt;

       }
       // End secondary calibration
       // ----------------------------------------

       // ----------------------------------------
       // If applicable perform third (NVTX) calibration
       if ( NVTXLUT.columns != 0){


         double gradient = NVTXLUT.getElement( etaIndex, 1 );
         double constant = NVTXLUT.getElement( etaIndex, 0 );

         // Calculate the mean jet response for a given NVTX
         double response = gradient*NVTX + constant;

         // Calculate NVTX-corrected Pt
         recalibNVTXCorrectedPt = recalibCorrectedPt/response;


         // 	 std::cout << "DeltaPt = " << deltaPt << "\tPt = " << recalibCorrectedPt 
         // 		   << "\tOffsetCorrectedPt = " << recalibOffsetCorrectedPt << "\n";


       }
       // End tertiary calibration
       // ----------------------------------------





       // Create and store unsorted calibrated L1Jet

       // Create un-sorted calibrated L1Jet 
       math::PtEtaPhiMLorentzVector tempJet;
       tempJet.SetCoordinates( correctedPt, uncalibJet.p4().eta(), uncalibJet.p4().phi(), uncalibJet.p4().M() );
       unsortedCalibratedL1Jets.push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );

       // Create un-sorted recalibrated L1Jet 
       math::PtEtaPhiMLorentzVector tempJetRecalib;
       tempJetRecalib.SetCoordinates( recalibCorrectedPt, uncalibJet.p4().eta(), uncalibJet.p4().phi(), uncalibJet.p4().M() );
       unsortedRecalibratedL1Jets.push_back( l1extra::L1JetParticle( tempJetRecalib, l1extra::L1JetParticle::JetType::kCentral, 0 ) );

       // Create un-sorted recalibrated, NVTX-corrected L1Jet 
       math::PtEtaPhiMLorentzVector tempJetNVTXRecalib;
       tempJetNVTXRecalib.SetCoordinates( recalibNVTXCorrectedPt, uncalibJet.p4().eta(), uncalibJet.p4().phi(), uncalibJet.p4().M() );
       unsortedNVTXRecalibratedL1Jets.push_back( l1extra::L1JetParticle( tempJetNVTXRecalib, l1extra::L1JetParticle::JetType::kCentral, 0 ) );

     }

     // Sort calibrated jets
     std::sort( unsortedCalibratedL1Jets.begin(),         unsortedCalibratedL1Jets.end(),         L1JetRankDescending );
     std::sort( unsortedRecalibratedL1Jets.begin(),       unsortedRecalibratedL1Jets.end(),       L1JetRankDescending );
     std::sort( unsortedNVTXRecalibratedL1Jets.begin(), unsortedNVTXRecalibratedL1Jets.end(), L1JetRankDescending );

     // Store sorted and calibrated TowerJets
     for (l1extra::L1JetParticleCollection::const_iterator Calib_It = unsortedCalibratedL1Jets.begin(); Calib_It != unsortedCalibratedL1Jets.end(); ++Calib_It ){

       // Store calibrated L1Jet
       outputCalibJet_L1Jet->push_back( (*Calib_It) );

     }
     for (l1extra::L1JetParticleCollection::const_iterator Recalib_It = unsortedRecalibratedL1Jets.begin(); Recalib_It != unsortedRecalibratedL1Jets.end(); ++Recalib_It ){

       // Store calibrated L1Jet
       outputRecalibJet_L1Jet->push_back( (*Recalib_It) );

     }
     for (l1extra::L1JetParticleCollection::const_iterator NVTXRecalib_It = unsortedNVTXRecalibratedL1Jets.begin(); 
         NVTXRecalib_It != unsortedNVTXRecalibratedL1Jets.end(); ++NVTXRecalib_It ){

       // Store calibrated L1Jet
       outputNVTXRecalibJet_L1Jet->push_back( (*NVTXRecalib_It) );

     }


}

// Read a LUT stored in a textfile and store in a flat vector
void 
JetCalibProducer::loadLUT( TString filename, std::vector< double > &LUTArr ){

  std::ifstream iReadLUT( filename );

  if(!iReadLUT){
    throw cms::Exception("MissingLUT")  << "A LUT could not be loaded as the file: '" << filename << "' does not exist, check the path of all LUTs.\n";
  } 

  // Read and store the LUT in a flat vector
  while (!iReadLUT.eof()){
    TString element;
    iReadLUT >> element;
    if (element != ""){
      LUTArr.push_back( element.Atof() );
    }
  }

}


// Determine the recalibration Pt correction factor
double 
JetCalibProducer::getRecalibFactor( double correctedPt, int row, LUT recalibLUT ){

  double jetResponse = 0;
  double powerPt = 1; // Pt to the current power

  for (int power = 0; power < recalibLUT.columns; ++power){

    jetResponse += recalibLUT.getElement( row, power ) * powerPt;

    // Increase to next power of Pt
    powerPt *= correctedPt;
  }

  return 1/jetResponse;
}





//define this as a plug-in
DEFINE_FWK_MODULE(JetCalibProducer);
