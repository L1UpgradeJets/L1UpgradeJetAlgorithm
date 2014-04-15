// -*- C++ -*-
//
// Package:    JetHist
// Class:      JetHist
// 
/**\class JetHist JetHist.cc AnalyseUpgradeJets/src/JetHist.cc

 Description: Produces histograms of jets


 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Baber
//         Created:  Tue Oct  3 15:21:18 BST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TEfficiency.h"
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






#include <algorithm>  // for sorting

// Ranking function for sort
bool TLorentzVectorRankDescending( TLorentzVector jet1, TLorentzVector jet2 ){ return (jet1.Pt() > jet2.Pt() ); }




// Print debugging messages
//#define VERBOSE
#include "AnalyseUpgradeJets/AnalyzeJets/interface/printing.h"
#include "AnalyseUpgradeJets/AnalyzeJets/interface/histBooker.h"
#include "AnalyseUpgradeJets/AnalyzeJets/interface/JetMatch.h"
// #include "AnalyseUpgradeJets/AnalyzeJets/src/helperFunctions.cc"



// **************************************************
// Histogram switches
// **************************************************

// Switch to just produce the required calibration plots. Take care turning this off again as the number of histograms booked
// will greatly increase and, with the current pT binning, may exceed the memory limit!!!
#define CALIBRATION
#define L1_JET_PLOTS
#define BENCHMARKING

//#define FASTJET

// **************************************************



//
// class declaration
//

//class JetHist : public edm::EDProducer {
class JetHist : public edm::EDAnalyzer {
   public:
      explicit JetHist(const edm::ParameterSet&);
      ~JetHist();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
  //      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  
      // ----------member functions ---------------------------

      virtual void reverseCumulative( TH1* histogram, TH1* rCumulHist ); 
 

      virtual void fillL1Histograms( edm::Handle<l1extra::L1JetParticleCollection> const& L1Jets, TString prefix );
  //      virtual void calibrateL1RECO( edm::Handle<l1slhc::L1TowerJetCollection> const& L1Jets, edm::Handle<l1extra::L1JetParticleCollection> const& ak5Jets,
      virtual void calibrateL1RECO( edm::Handle<l1extra::L1JetParticleCollection> const& L1Jets, edm::Handle<l1extra::L1JetParticleCollection> const& ak5Jets,
				    TString prefix, double ak5PTthreshold, int maxMatchedJets, double ak5EtaRange);
  //				    TString prefix, double ak5PTthreshold, double L1PTthreshold, int maxMatchedJets, double ak5EtaRange );

//    virtual void matchOnlineOfflineJets( edm::Handle<l1slhc::L1TowerJetCollection> const& L1Jets, edm::Handle<l1extra::L1JetParticleCollection> const& ak5Jets,
// 					TString prefix, double ak5PTthreshold, double L1PTthreshold  );
      virtual void benchmarkJets( edm::Handle<l1extra::L1JetParticleCollection> const& L1Jets, edm::Handle<l1extra::L1JetParticleCollection> const& ak5Jets,
				  TString prefix, double ak5PTthreshold, double L1PTthreshold, int maxMatchedJets, double ak5EtaRange, double ak5MaxPT );


      virtual void turnOnL1Loop( std::vector<TLorentzVector> const &L1Vec, std::vector<TLorentzVector> const &ak5Vec, TString matchPre );
      virtual void turnOnAk5MatchedLoop( std::vector<TLorentzVector> const &L1Vec, std::vector<TLorentzVector> const &ak5Vec, TString matchPre );
      virtual void turnOnAk5Loop( std::vector<TLorentzVector> const &L1Vec, std::vector<TLorentzVector> const &ak5Vec, TString matchPre );

      virtual void studyMatchedJets( std::vector<TLorentzVector> const &L1Vec, std::vector<TLorentzVector> const &ak5Vec, TString matchPre, double maxDeltaR );

      // ----------member data ---------------------------
      edm::ParameterSet conf_;
      edm::Service<TFileService> fs;


      // Histogram containers
      std::map< TString, TH1*> hist1D;
      std::map< TString, TH1*> hist2D;
      std::map< TString, TH3*> hist3D;
      std::map< TString, TEfficiency*> histEff;

      // List of currently booked histograms
//       std::list<TString> histogram1DList;
//       std::list<TString> histogram2DList;

      // Number of reconstructed ak5 primary vertices 
      int NVTX;
      // Event number
      int Eventnr;

      // Jet collections utilised
      std::vector <TString> jetList;
  
      // Label for the 'ith' jet starting from zero!!! => LeadJet = 0
      std::vector <TString> jetIndexLabel;

      // Boundaries of the eta regions in which to apply local PU subtraction
      std::vector <double> localRhoEtaDivisions;
      // Boundaries of PT in which to apply local PT binning                                                             
      std::vector <double> ptSlice;
      // Specify the boundaries of PU in which to apply local PU binning                                                             
      std::vector <double> puSlice;
      // PT thresholds                                                                                                   
      std::vector <double> ptThreshold;
      // Eta region-level segmentation
      std::vector <double> etaRegionSlice;


      // ****************************************
      // *      Online-offline jet matching     *
      // ****************************************
    
      // Jet matching maximum delta R
      double maxDeltaR;

  
      // Online-Offline jet matching
  std::vector <TString> matchPrefix;
  std::vector <TString> matchLabel;

  std::vector <TString> calibPrefix;
  std::vector <TString> calibLabel;

      // Width of iEta bins to sample to obtain pT corrections
      double iEtaCalibBinWidth;
      
  int jetRankComparison;

  // NVTX binning
  std::vector<int> nvtxBins;

  // ****************************************
  // *         Turn on configuration        *
  // ****************************************

  // Online PT trigger PT thresholds
  std::vector < int > onlineTOSingleJetPT;



};
JetHist::JetHist(const edm::ParameterSet& iConfig): conf_(iConfig)
{


  // Associate the histogram container with the histogram booker, a class for handling histogram booking
  histBooker booker( &hist1D, &hist2D );



  // ****************************************************************************************************
  // *                                 Load configuration parameters                                    *
  // ****************************************************************************************************

  // Boundaries of the eta regions in which to apply local PU subtraction                                            
  localRhoEtaDivisions = iConfig.getParameter< std::vector<double> > ("LocalRhoEtaDivisions");
  sort (localRhoEtaDivisions.begin(), localRhoEtaDivisions.end());   // ensure the bins are in ascending order 

  // Boundaries of PT in which to apply local PT binning                                                             
  ptSlice              = iConfig.getParameter< std::vector<double> >("PTSlice");
  sort (ptSlice.begin(), ptSlice.end()); // ensure the bins are in ascending order

  // Specify the boundaries of PU in which to apply local PU binning                                                             
  puSlice              = iConfig.getParameter< std::vector<double> >("PUSlice");
  sort (puSlice.begin(), puSlice.end()); // ensure the bins are in ascending order  

  // PT thresholds                                                                                                   
  ptThreshold          = iConfig.getParameter< std::vector<double> >("PTThreshold");
  sort (ptThreshold.begin(), ptThreshold.end()); // ensure the bins are in ascending order

  // Eta region-level segmentation
  etaRegionSlice        = iConfig.getParameter< std::vector<double> >("EtaRegionSlice");
  sort (etaRegionSlice.begin(), etaRegionSlice.end());  // ensure the bins are in ascending order                                                                                                                    


  // NVTX binning
    nvtxBins.push_back( 0 );
    nvtxBins.push_back( 10 );
    nvtxBins.push_back( 15 );
    nvtxBins.push_back( 20 );
    nvtxBins.push_back( 25 );
    nvtxBins.push_back( 30 );
    std::sort(  nvtxBins.begin(), nvtxBins.end() );



  // ************************************************************
  // *                  Online-offline matching                 *
  // ************************************************************
  SUBPRINT("Jet matching")

  // Maximum delta R in online-offline jet matching
  maxDeltaR = iConfig.getParameter< double >("MaxDeltaR");

  iEtaCalibBinWidth = iConfig.getParameter< double >("iEtaCalibrationBinWidth");

  //Rank of jets to compare
  jetRankComparison = iConfig.getParameter< int >("JetRankComparison");

  // ************************************************************
  // *                   Turn on configuration                  *
  // ************************************************************
  SUBPRINT("Turn ons")

  // Online PT trigger PT thresholds
    onlineTOSingleJetPT = iConfig.getParameter< std::vector<int> >("onlineTurnOnSingleJetPT");
  // Offline turn on PT bins
  //  offlineTOSingleJetPTBin;
  // **********************************************************************************
  // *               Book histograms for current and upgrade algorithms               *
  // **********************************************************************************
  
  // Region-level eta binning
  const Int_t regionEtaBin = 14;   
  const double regionEtaBins[] = {  -3.0, -2.172, -1.74, -1.392, -1.044, -0.695, -0.348, 0.0,
    0.348, 0.695, 1.044, 1.392, 1.74, 2.172, 3.0};

  jetList.push_back("PrePUS");      // PrePUS upgrade Jets
  jetList.push_back("PUS");      	// PUS upgrade Jets
  jetList.push_back("LPUS");      	// PUS upgrade Jets
  jetList.push_back("CalibPrePUS");      // PrePUS upgrade Jets
  jetList.push_back("CalibPUS");      // PrePUS upgrade Jets
  jetList.push_back("CalibLPUS");      // PrePUS upgrade Jets
  jetList.push_back("Ak5");             // ak5 PUS calibrated Jets
  // ith Jet lists
  // ************************************************************

  jetIndexLabel.push_back("LeadJet");
  jetIndexLabel.push_back("SecondJet");
  jetIndexLabel.push_back("ThirdJet");
  jetIndexLabel.push_back("FourthJet");




  // ----------------------------------------
  //  Jet collections to be matched
  // ----------------------------------------

  // NOTE: WE ARE ON THE MEMORY BOUNDARY, ADDING ANY MORE
  //       COMPARISONS WILL LIKELY EXCEED THE MEMORY LIMIT : (
  //

#ifdef FASTJET
  calibPrefix.push_back("Calibration_TTAk5PrePUSRaw_ak5PrePUSRaw");
  calibLabel.push_back("TT Ak5 PrePUS Raw, Ak5 PrePUS Raw");
#endif

#ifdef CALIBRATION



  // PrePUS-ak5PUS
  calibPrefix.push_back("Calibration_PrePUS_ak5PUS");
  calibLabel.push_back("TowerJet PrePUS, Ak5 PUS");
  calibPrefix.push_back("Calibration_PrePUS_ak5PUS_3Jets");
  calibLabel.push_back("TowerJet PrePUS leading three jets, Ak5 PUS");
  calibPrefix.push_back("Calibration_PrePUS_ak5PUS_AllJets");
  calibLabel.push_back("TowerJet PrePUS, Ak5 PUS_AllJets");

  calibPrefix.push_back("Calibration_CalibPrePUS_ak5PUS");
  calibLabel.push_back("TowerJet CalibPrePUS, Ak5 PUS");
  calibPrefix.push_back("Calibration_CalibPrePUS_ak5PUS_3Jets");
  calibLabel.push_back("TowerJet CalibPrePUS leading three jets, Ak5 PUS");
  calibPrefix.push_back("Calibration_CalibPrePUS_ak5PUS_AllJets");
  calibLabel.push_back("TowerJet CalibPrePUS, Ak5 PUS_AllJets");

  calibPrefix.push_back("Calibration_PUS_ak5PUS");
  calibLabel.push_back("TowerJet PUS, Ak5 PUS");
  calibPrefix.push_back("Calibration_PUS_ak5PUS_3Jets");
  calibLabel.push_back("TowerJet PUS leading three jets, Ak5 PUS");
  calibPrefix.push_back("Calibration_PUS_ak5PUS_AllJets");
  calibLabel.push_back("TowerJet PUS, Ak5 PUS_AllJets");

  calibPrefix.push_back("Calibration_LPUS_ak5PUS");
  calibLabel.push_back("TowerJet LPUS, Ak5 PUS");
  calibPrefix.push_back("Calibration_LPUS_ak5PUS_3Jets");
  calibLabel.push_back("TowerJet LPUS leading three jets, Ak5 PUS");
  calibPrefix.push_back("Calibration_LPUS_ak5PUS_AllJets");
  calibLabel.push_back("TowerJet LPUS, Ak5 PUS_AllJets");

  calibPrefix.push_back("Calibration_CalibPUS_ak5PUS");
  calibLabel.push_back("TowerJet Calib PUS, Ak5 PUS");
  calibPrefix.push_back("Calibration_CalibPUS_ak5PUS_3Jets");
  calibLabel.push_back("TowerJet Calib PUS leading three jets, Ak5 PUS");
  calibPrefix.push_back("Calibration_CalibPUS_ak5PUS_AllJets");
  calibLabel.push_back("TowerJet Calib PUS, Ak5 PUS_AllJets");

  calibPrefix.push_back("Calibration_CalibLPUS_ak5PUS");
  calibLabel.push_back("TowerJet CalibLPUS, Ak5 PUS");
  calibPrefix.push_back("Calibration_CalibLPUS_ak5PUS_3Jets");
  calibLabel.push_back("TowerJet CalibLPUS leading three jets, Ak5 PUS");
  calibPrefix.push_back("Calibration_CalibLPUS_ak5PUS_AllJets");
  calibLabel.push_back("TowerJet CalibLPUS, Ak5 PUS_AllJets");

  //calibPrefix.push_back("Calibration_CalibPrePUS_ak5PUS");
  //calibLabel.push_back("TowerJet CalibPrePUS, Ak5 PUS");
  //calibPrefix.push_back("Calibration_Curr_ak5PUS");
  //calibLabel.push_back("Current, AK5 PUS");

  //calibPrefix.push_back("Calibration_UncalibCurr_ak5PUS");
  //calibLabel.push_back("Current uncalibrated, AK5 PUS");



#endif


  // ****************************************************************************************************
  // *                                           JET MATCHING                                           *
  // ****************************************************************************************************

#ifdef BENCHMARKING

  // Hand-calibrated ak5PUS
//   matchPrefix.push_back("CalibAk5PUSRawAk5PUS_ak5PUS");
//   matchLabel.push_back("CalibAk5PUSRaw, Ak5 PUS");



  // Current L1 GCT
  //matchPrefix.push_back("Curr_ak5PUS");
  //matchLabel.push_back("Current jet, Ak5 PUS");
  //matchPrefix.push_back("CurrUncalib_ak5PUS");
  //matchLabel.push_back("Current uncalibrated jet, Ak5 PUS");

  // Calibrated jets
  // --------------------

  // PrePUS-ak5PUS
  matchPrefix.push_back("PrePUS_ak5PUS");
  matchLabel.push_back("TowerJet PrePUS, Ak5 PUS");

  matchPrefix.push_back("PUS_ak5PUS");
  matchLabel.push_back("TowerJet PUS, Ak5 PUS");

  matchPrefix.push_back("CalibPrePUS_ak5PUS");
  matchLabel.push_back("TowerJet CalibPrePUS, Ak5 PUS");

  matchPrefix.push_back("CalibPUS_ak5PUS");
  matchLabel.push_back("TowerJet CalibPUS, Ak5 PUS");



#endif

  
  // ************************************************************
  // *                      L1 Jet plots                        *
  // ************************************************************

#ifdef L1_JET_PLOTS
  
  TFileDirectory etaBinnedDir = fs->mkdir( "etaBinned" );
  TFileDirectory jetDir = fs->mkdir( "Jets" );


  PRINT("Histogram jet booking")
     
      
    for ( uint iJet = 0; iJet < jetList.size(); iJet++ ){
      
      // Extract jet type prefix
      TString prefix           = jetList[iJet];
      TFileDirectory jetSubDir = jetDir.mkdir( prefix.Data() );
      
      // Print current prefix
      SUBPRINT(prefix)
	  
      // Human-readable jet labels
      TString jetLab;
      if      (prefix == "Curr")        jetLab = "Current";  
      else if (prefix == "PrePUS")      jetLab = "Upgrade PrePUS";
      else if (prefix == "PrePUSCalib") jetLab = "Upgrade PrePUS Calibrated";
      else if (prefix == "PrePUSRecalib") jetLab = "Upgrade PrePUS Recalibrated";
      else if (prefix == "PrePUSNVTXRecalib") jetLab = "Upgrade PrePUS NVTX-corrected, Recalibrated";
      else if (prefix == "PrePUSCalibLt3") jetLab = "Upgrade PrePUS Calibrated |#eta| < 3";
      else if (prefix == "PUS")         jetLab = "Upgrade Global PUS";
      else if (prefix == "PUSCalib")    jetLab = "Upgrade Global PUS Calibrated";
      else if (prefix == "LPUS")        jetLab = "Upgrade Local PUS";
      else if (prefix == "LPUSCalib")   jetLab = "Upgrade Local PUS Calibrated";
      else if (prefix == "Ak5PrePUS")   jetLab = "Ak5 PrePUS";
      else if (prefix == "Ak5PUS")      jetLab = "Ak5 PUS";
      


      // ******************************
      // General jet distributions
      // ******************************

      // Jet multiplicity
      booker.book1D( prefix + "_NJets", jetSubDir, "Number of Reconstructed Jets - " + jetLab + ";Jet Multiplicity;Events", pJetMultiplicity );

      
      // ******************************
      // 'ith' jet distributions
      // ******************************

      // PT
      booker.book1D( prefix + "_Jet_1_PT", jetSubDir, "Leading jet p_{T} - " + jetLab + ";p_{T} (GeV);Entries/(1 GeV)", pOffPTLarge );
      booker.book1D( prefix + "_Jet_2_PT", jetSubDir, "Second jet p_{T} - "  + jetLab + ";p_{T} (GeV);Entries/(1 GeV)", pOffPTLarge );
      booker.book1D( prefix + "_Jet_3_PT", jetSubDir, "Third jet p_{T} - "   + jetLab + ";p_{T} (GeV);Entries/(1 GeV)", pOffPTLarge );
      booker.book1D( prefix + "_Jet_4_PT", jetSubDir, "Fourth jet p_{T} - "  + jetLab + ";p_{T} (GeV);Entries/(1 GeV)", pOffPTLarge );

      // Eta
      booker.book1D( prefix + "_Jet_1_Eta", jetSubDir, "Leading jet #eta - " + jetLab + ";#eta;Entries", pOffEta );
      booker.book1D( prefix + "_Jet_2_Eta", jetSubDir, "Second jet #eta - "  + jetLab + ";#eta;Entries", pOffEta );
      booker.book1D( prefix + "_Jet_3_Eta", jetSubDir, "Third jet #eta - "   + jetLab + ";#eta;Entries", pOffEta );
      booker.book1D( prefix + "_Jet_4_Eta", jetSubDir, "Fourth jet #eta - "  + jetLab + ";#eta;Entries", pOffEta );

      // Phi
      booker.book1D( prefix + "_Jet_1_Phi", jetSubDir, "Leading jet #phi - " + jetLab + ";#phi;Entries", pOffPhi );
      booker.book1D( prefix + "_Jet_2_Phi", jetSubDir, "Second jet #phi - "  + jetLab + ";#phi;Entries", pOffPhi );
      booker.book1D( prefix + "_Jet_3_Phi", jetSubDir, "Third jet #phi - "   + jetLab + ";#phi;Entries", pOffPhi );
      booker.book1D( prefix + "_Jet_4_Phi", jetSubDir, "Fourth jet #phi - "  + jetLab + ";#phi;Entries", pOffPhi );
      
      // Phi vs Eta
      booker.book2D( prefix + "_Jet_1_Phi_vs_Eta", jetSubDir, "Number of reconstructed jets #phi vs #eta (Leading) - " + jetLab + ";#eta;#phi", pOffEta, pOffPhi );
      booker.book2D( prefix + "_Jet_2_Phi_vs_Eta", jetSubDir, "Number of reconstructed jets #phi vs #eta (Second) - "  + jetLab + ";#eta;#phi", pOffEta, pOffPhi );
      booker.book2D( prefix + "_Jet_3_Phi_vs_Eta", jetSubDir, "Number of reconstructed jets #phi vs #eta (Third) - "   + jetLab + ";#eta;#phi", pOffEta, pOffPhi );
      booker.book2D( prefix + "_Jet_4_Phi_vs_Eta", jetSubDir, "Number of reconstructed jets #phi vs #eta (Fourth) - "  + jetLab + ";#eta;#phi", pOffEta, pOffPhi );


      // ******************************
      // All jet distributions
      // ******************************

      // PT, eta, phi
      booker.book1D( prefix + "_PT",  jetSubDir, "All jet p_{T} - " + jetLab + ";p_{T} (GeV);Entries", pOffPTLarge );
      booker.book1D( prefix + "_Eta", jetSubDir, "All jet #eta - "  + jetLab + ";#eta;Entries", pOffEta );
      booker.book1D( prefix + "_Phi", jetSubDir, "All jet #phi - "  + jetLab + ";#phi;Entries", pOffPhi );

      // Phi vs Eta
      booker.book2D( prefix + "_Phi_vs_Eta", jetSubDir, "Number of reconstructed jets (Et Weighted) #phi vs #eta - " + jetLab + ";#eta;#phi", pOffEta, pOffPhi );
      booker.book2D( prefix + "_Eta_vs_Pt", jetSubDir, "Number of reconstructed jets #eta vs p_{T} - " + jetLab + ";p_{T} (GeV);#eta", 
		     pOffPTLarge, pOffEta );
      booker.book2D( prefix + "_Phi_vs_Pt", jetSubDir, "Number of reconstructed jets #phi vs p_{T} - " + jetLab + ";p_{T} (GeV);#phi", 
		     pOffPTLarge, pOffPhi );


	
    }

#endif





  // ********************************************************************************
  // *                                 Jet Matching                                 *
  // ********************************************************************************
  // This is primarily for calibration....
  PRINT("Jet matching")

    // Jet calibration code
  for (unsigned int iCalib = 0; iCalib < calibPrefix.size(); iCalib++){
    
    
    TString matchPre = calibPrefix[iCalib];
    TString matchLab = calibLabel[iCalib];
    TString matchDirName = "Calibration";

    TFileDirectory jetMatchingDir = fs->mkdir( matchDirName.Data() );
    // Create a calibration subsubdirectory
    TFileDirectory matchSubDir = jetMatchingDir.mkdir( matchPre.Data() );


    // pT Correlations
    //    hist3D[matchPre + "_L1PT_vs_OffPT_vs_NVTX"] = matchSubDir.make<TH3D>( matchPre + "_L1PT_vs_OffPT_vs_NVTX", "L1 P_{T} vs Offline P_{T} vs N_{VTX};Offline p_{T} (GeV);L1 P_{T} (GeV);N_{VTX}", pOffPT.bins, pOffPT.low, pOffPT.high, pOffPT.bins, pOffPT.low, pOffPT.high, pPU.bins, pPU.low, pPU.high );


    booker.book2DTProf( matchPre + "_L1PT_vs_OffPT" ,     matchSubDir, "L1 P_{T} vs Offline P_{T};Offline p_{T} (GeV);L1 P_{T} (GeV)", pOffPTLarge, pL1PTLarge);
    booker.book2DTProf( matchPre + "_OffPT_vs_L1PT" ,     matchSubDir, "Offline P_{T} vs L1 P_{T};L1 P_{T} (GeV);Offline p_{T} (GeV)", pOffPTLarge, pL1PTLarge);
    booker.book2DTProf( matchPre + "_DeltaPTRel_vs_OffPT"  , matchSubDir, "#DeltaP_{T}^{Rel} vs Offline P_{T};Offline p_{T} (GeV);#DeltaP_{T}^{Rel}", pOffPT, pDeltaPTRel); 
    booker.book2DTProf( matchPre + "_DeltaPTRel_vs_DeltaPT", matchSubDir, "#DeltaP_{T}^{Rel} vs #DeltaP_{T};#DeltaP_{T} (GeV);#DeltaP_{T}^{Rel}", pDeltaPT, pDeltaPTRel); 

    booker.book2DTProf( matchPre + "_DeltaPT_vs_OffPT"  , matchSubDir, "#DeltaP_{T} vs Offline P_{T};Offline p_{T} (GeV);#DeltaP_{T} (GeV)", pOffPT, pDeltaPT); 

    booker.book1D( matchPre + "_JetResponse",  matchSubDir, "Jet response;P_{T}^{L1}/P_{T}^{RECO};Entries", pResponse);
    booker.book2DTProf( matchPre + "_JetResponse_vs_OffPT",  matchSubDir, "Jet response vs Offline p_{T};Offline p_{T} (GeV);P_{T}^{L1}/P_{T}^{RECO}", pOffPTLarge, pResponse);
    booker.book2DTProf( matchPre + "_JetResponse_vs_L1PT",  matchSubDir, "Jet response vs L1 p_{T};L1 p_{T} (GeV);P_{T}^{L1}/P_{T}^{RECO}", pOffPT, pResponse);

    //    booker.book2DTProf( matchPre + "_JetResponse_vs_iEta",  matchSubDir, "Jet response vs i#eta;i#eta;P_{T}^{L1}/P_{T}^{RECO}", piEta, pResponse);
    booker.book2DTProf( matchPre + "_JetResponse_vs_L1Eta",  matchSubDir, "Jet response vs L1 #eta;i#eta;P_{T}^{L1}/P_{T}^{RECO}", pOffEta, pResponse);
    booker.book2DTProf( matchPre + "_JetResponse_vs_NVTX_L1Ptge12",  matchSubDir,
			"Jet response vs N_{VTX} (L1 P_{T} #ge 12GeV);N_{VTX};P_{T}^{L1}/P_{T}^{RECO}", pPU, pResponse);
    booker.book2DTProf( matchPre + "_JetResponse_vs_NVTX_L1Ptge30",  matchSubDir,
			"Jet response vs N_{VTX} (L1 P_{T} #ge 30GeV);N_{VTX};P_{T}^{L1}/P_{T}^{RECO}", pPU, pResponse);



    // ************************************************************
    // *                  Eta binned distributions                *
    // ************************************************************
    // Bin the calibration distributions in eta

    SUBPRINT("Eta binned")
      TFileDirectory EtaDir = matchSubDir.mkdir( "EtaBinned" );

    for (uint pEta = 1; pEta < etaRegionSlice.size(); ++pEta){

      // Get Eta bin lower and upper bounds
      double EtaLow  = etaRegionSlice[ pEta - 1];
      double EtaHigh = etaRegionSlice[ pEta ];
      TString EtaLowStr  = Form("%1.3f", EtaLow );
      TString EtaHighStr = Form("%1.3f", EtaHigh );

      TString EtaStr = "Eta_" + EtaLowStr + "_to_" + EtaHighStr; 
      TFileDirectory EtaSubDir = EtaDir.mkdir( EtaStr.Data() );

	booker.book2DTProf( matchPre + "_L1PT_vs_OffPT_"     + EtaStr, EtaSubDir, "L1 P_{T} vs Offline P_{T} (" + EtaStr + ");Offline p_{T} (GeV);L1 P_{T} (GeV)", 										 pOffPTLarge, pL1PTLarge);
	booker.book2DTProf( matchPre + "_OffPT_vs_L1PT_" + EtaStr, EtaSubDir, "Offline P_{T} vs L1 P_{T} (" + EtaStr + 
			    ");L1 P_{T} (GeV);Offline p_{T} (GeV)", pL1PTLarge, pOffPTLarge);

	booker.book2DTProf( matchPre + "_DeltaPTRel_vs_OffPT_"  + EtaStr, EtaSubDir, "#DeltaP_{T}^{Rel} vs Offline P_{T} (" + EtaStr + ");Offline p_{T} (GeV);#DeltaP_{T}^{Rel}", 
			    pOffPT, pDeltaPTRel); 

	booker.book2DTProf( matchPre + "_DeltaPTRel_vs_DeltaPT_"  + EtaStr, EtaSubDir, "#DeltaP_{T}^{Rel} vs #DeltaP_{T} (" + EtaStr + ");#DeltaP_{T} (GeV);#DeltaP_{T}^{Rel}", 
			    pDeltaPT, pDeltaPTRel); 

	booker.book2DTProf( matchPre + "_DeltaPT_vs_OffPT_"  + EtaStr, EtaSubDir, "#DeltaP_{T} vs Offline P_{T} (" + EtaStr + ");Offline p_{T} (GeV);#DeltaP_{T} (GeV)", 
			    pOffPT, pDeltaPT); 
	booker.book1D( matchPre + "_JetResponse_" + EtaStr,  EtaSubDir, "Jet response;P_{T}^{L1}/P_{T}^{RECO};Entries", pResponse);
	booker.book2DTProf( matchPre + "_JetResponse_vs_OffPT_" + EtaStr,  EtaSubDir, "Jet response vs Offline P_{T};Offline p_{T} (GeV);P_{T}^{L1}/P_{T}^{RECO}", pOffPTLarge, pResponse);


	// Recalib plots
	TFileDirectory RecalibEtaSubDir = EtaSubDir.mkdir( "Recalib" );

	booker.book2DTProf( matchPre + "_JetResponse_vs_L1PT_" + EtaStr,  RecalibEtaSubDir, 
			    "Jet response vs L1 P_{T};L1 p_{T} (GeV);P_{T}^{L1}/P_{T}^{RECO}", pL1PT, pResponse);
	booker.book2DTProf( matchPre + "_DeltaPT_vs_L1PT_"  + EtaStr, RecalibEtaSubDir, 
			    "#DeltaP_{T} vs L1 P_{T} (" + EtaStr + ");L1 p_{T} (GeV);#DeltaP_{T} (GeV)", pL1PT, pDeltaPT); 

	booker.book2DTProf( matchPre + "_JetResponse_vs_NVTX_" + EtaStr,  RecalibEtaSubDir, 
			    "Jet response vs N_{VTX};N_{VTX};P_{T}^{L1}/P_{T}^{RECO}", pPU, pResponse);
	booker.book2DTProf( matchPre + "_JetResponse_vs_NVTX_L1Ptge12_" + EtaStr,  RecalibEtaSubDir,
                            "Jet response vs N_{VTX} (L1 P_{T} #ge 12GeV);N_{VTX};P_{T}^{L1}/P_{T}^{RECO}", pPU, pResponse);
	booker.book2DTProf( matchPre + "_JetResponse_vs_NVTX_L1Ptge30_" + EtaStr,  RecalibEtaSubDir,
                            "Jet response vs N_{VTX} (L1 P_{T} #ge 30GeV);N_{VTX};P_{T}^{L1}/P_{T}^{RECO}", pPU, pResponse);


	booker.book2DTProf( matchPre + "_DeltaPT_vs_NVTX_"  + EtaStr, RecalibEtaSubDir, 
			    "#DeltaP_{T} vs N_{VTX} (" + EtaStr + ");N_{VTX};#DeltaP_{T} (GeV)", pPU, pDeltaPT); 


      
      } // End Eta loop




  } // End calibration loop





  for (unsigned int iMatch = 0; iMatch < matchPrefix.size(); iMatch++){
    
    
    TString matchPre = matchPrefix[iMatch];
    TString matchLab = matchLabel[iMatch];


    // Make directory to store plots depending on their contents
    TString matchDirName = "L1RECOMatching";

    TFileDirectory jetMatchingDir = fs->mkdir( matchDirName.Data() );

    // Create a match subsubdirectory
    TFileDirectory matchSubDir = jetMatchingDir.mkdir( matchPre.Data() );



    // Resolution distributions
    booker.book1D( matchPre + "_DeltaEta",   matchSubDir, "Lead four jets #Delta#eta - "  + matchLab + ";#Delta#eta;Entries", pDeltaEta);
    booker.book1D( matchPre + "_DeltaPhi",   matchSubDir, "Lead four jets #Delta#phi - "  + matchLab + ";#Delta#phi;Entries", pDeltaPhi);
    booker.book1D( matchPre + "_DeltaPTRel", matchSubDir, "Lead four jets #Deltap_{T} - " + matchLab + ";#Deltap_{T};Entries", pDeltaPTRel);

    
    // Specific jet distributions
    for (uint iJet = 0; iJet < jetIndexLabel.size(); ++iJet ){

      // jet labeling starting at 0. LeadJet = 0, ...                                                                  
      TString jetLabel = "_" + jetIndexLabel[ iJet ];
      TString jetText  = jetIndexLabel[ iJet ];

      // 'ith' jet distributions
      booker.book1D( matchPre + jetLabel + "_DeltaEta",   matchSubDir, jetText + " #Delta#eta - "  + matchLab + ";#Delta#eta;Entries", pDeltaEta);
      booker.book1D( matchPre + jetLabel + "_DeltaPhi",   matchSubDir, jetText + " #Delta#phi - "  + matchLab + ";#Delta#phi;Entries", pDeltaPhi);
      booker.book1D( matchPre + jetLabel + "_DeltaPTRel", matchSubDir, jetText + " #Deltap_{T} - " + matchLab + ";#Deltap_{T};Entries", pDeltaPTRel);
      booker.book2DTProf( matchPre + jetLabel + "_DeltaPTRel_vs_OffPT"  , matchSubDir, jetText + " #DeltaP_{T}^{Rel} vs Offline P_{T};Offline p_{T} (GeV);#DeltaP_{T}^{Rel}", pOffPT, pDeltaPTRel); 
      booker.book2DTProf( matchPre + jetLabel + "_DeltaPTRel_vs_L1PT"  , matchSubDir, jetText + " #DeltaP_{T}^{Rel} vs L1 P_{T};L1 p_{T} (GeV);#DeltaP_{T}^{Rel}", pOffPT, pDeltaPTRel); 
      booker.book2DTProf( matchPre + jetLabel + "_DeltaPTRel_vs_NVTX"  , matchSubDir, jetText + " #DeltaP_{T}^{Rel} vs N_{VTX};N_{VTX};#DeltaP_{T}^{Rel}", 
			  pPU, pDeltaPTRel);
      booker.book2DTProf( matchPre + jetLabel + "_DeltaPTRel_vs_Eta"  , matchSubDir, jetText + " #DeltaP_{T}^{Rel} vs #eta;L1 #eta;#DeltaP_{T}^{Rel}", 
			  pOffEta, pDeltaPTRel);

      booker.book2DTProf( matchPre + jetLabel + "_JetResponse_vs_OffPT"  , matchSubDir, jetText + " Jet response vs Offline P_{T};Offline p_{T} (GeV);P_{T}^{L1}/P_{T}^{RECO}", pOffPT, pResponse); 
      booker.book2DTProf( matchPre + jetLabel + "_JetResponse_vs_NVTX"  , matchSubDir, jetText + " Jet response vs N_{VTX};N_{VTX};P_{T}^{L1}/P_{T}^{RECO}", 
			  pPU, pResponse);
      booker.book2DTProf( matchPre + jetLabel + "_JetResponse_vs_Eta"  , matchSubDir, jetText + " Jet response vs #eta;L1 #eta;P_{T}^{L1}/P_{T}^{RECO}", 
			  pOffEta, pResponse);

      booker.book2D( matchPre + jetLabel + "_deltaPTRel_vs_deltaR", matchSubDir, "L1 " + jetText + " #Deltap_{T}^{rel} vs #DeltaR"
		     + matchLab + ";#DeltaR;#Deltap_{T}^{rel}", pDeltaR, pDeltaPTRel );



    }


    // Lead four jetss

    booker.book1D( matchPre + "_Eta",      matchSubDir, "Lead four jets #eta - "  + matchLab + ";#eta;Entries", pOffEta );
    booker.book1D( matchPre + "_Phi",      matchSubDir, "Lead four jets #phi - "  + matchLab + ";#phi;Entries", pOffPhi );
    booker.book1D( matchPre + "_PT",       matchSubDir, "Lead four jets p_{T} - " + matchLab + ";p_{T} (GeV);Entries", pOffPT );



    booker.book1D( matchPre + "_Ak5Eta",   matchSubDir, "Ak5 lead jet #eta - "  + matchLab + ";#eta;Entries", pOffEta );
    booker.book1D( matchPre + "_Ak5Phi",   matchSubDir, "Ak5 lead jet #phi - "  + matchLab + ";#phi;Entries", pOffPhi );
    booker.book1D( matchPre + "_Ak5PT",    matchSubDir, "Ak5 lead jet p_{T} - " + matchLab + ";p_{T} (GeV);Entries", pOffPT );



    // Energy sums
    // ------------------------------------------------------------
    TFileDirectory energySumsSubDir = matchSubDir.mkdir( "EnergySums" );

    booker.book1D( matchPre + "_HT",       energySumsSubDir, "H_{T} - " + matchLab + ";H_{T} (GeV);Entries", pL1HT );
    booker.book1D( matchPre + "_MHT",      energySumsSubDir, "#slash{H}_{T} - " + matchLab + ";#slash{H}_{T} (GeV);Entries", pL1MHT );
    booker.book1D( matchPre + "_Ak5HT",    energySumsSubDir, "H_{T} - " + matchLab + ";H_{T} (GeV);Entries", pL1HT );
    booker.book1D( matchPre + "_Ak5MHT",   energySumsSubDir, "#slash{H}_{T} - " + matchLab + ";#slash{H}_{T} (GeV);Entries", pL1MHT );

    booker.book2DTProf( matchPre + "_L1MHT_vs_L1HT",    energySumsSubDir, "L1 #slash{H}_{T} vs L1 H_{T};L1 H_{T} (GeV);L1 #slash{H}_{T} (GeV)", pL1HT, pL1MHT);
    booker.book2DTProf( matchPre + "_Ak5MHT_vs_Ak5HT",  energySumsSubDir, "Ak5 #slash{H}_{T} vs Ak5 H_{T};Ak5 H_{T} (GeV);Ak5 #slash{H}_{T} (GeV)", pOffHT, pOffMHT);
    booker.book2DTProf( matchPre + "_Ak5HT_vs_L1HT" ,   energySumsSubDir, "Ak5 H_{T} vs L1 H_{T};L1 H_{T} (GeV);Ak5 H_{T} (GeV)", pL1HT, pOffHT);
    booker.book2DTProf( matchPre + "_Ak5MHT_vs_L1MHT" , energySumsSubDir, "Ak5 #slash{H}_{T} vs L1 #slash{H}_{T};L1 #slash{H}_{T} (GeV);Ak5 #slash{H}_{T} (GeV)", pL1MHT, pOffMHT);


    // HT turnons
    histEff[ matchPre + "_TurnOn_HTgt75GeV" ]  = energySumsSubDir.make<TEfficiency>( matchPre + "_TurnOn_HTgt75GeV", 
										     "H_{T} > 75 GeV;Offline H_{T} (GeV);Efficiency",
										     131, -2.5, 652.5 );
    histEff[ matchPre + "_TurnOn_HTgt100GeV" ] = energySumsSubDir.make<TEfficiency>( matchPre + "_TurnOn_HTgt100GeV", 
										     "H_{T} > 100 GeV;Offline H_{T} (GeV);Efficiency",
										     131, -2.5, 652.5 );
    histEff[ matchPre + "_TurnOn_HTgt150GeV" ] = energySumsSubDir.make<TEfficiency>( matchPre + "_TurnOn_HTgt150GeV", 
										     "H_{T} > 150 GeV;Offline H_{T} (GeV);Efficiency",
										     131, -2.5, 652.5 );
    histEff[ matchPre + "_TurnOn_HTgt175GeV" ] = energySumsSubDir.make<TEfficiency>( matchPre + "_TurnOn_HTgt175GeV", 
										     "H_{T} > 175 GeV;Offline H_{T} (GeV);Efficiency",
										     131, -2.5, 652.5 );

    // MHT turnons
    histEff[ matchPre + "_TurnOn_MHTgt50GeV" ] = energySumsSubDir.make<TEfficiency>( matchPre + "_TurnOn_MHTgt50GeV", 
										     "#slash_{H}_{T} > 50 GeV;Offline #slash_{H}_{T} (GeV);Efficiency",
										     61, -2.5, 302.5 );
    histEff[ matchPre + "_TurnOn_MHTgt75GeV" ]  = energySumsSubDir.make<TEfficiency>( matchPre + "_TurnOn_MHTgt75GeV", 
										     "#slash_{H}_{T} > 75 GeV;Offline #slash_{H}_{T} (GeV);Efficiency",
										     61, -2.5, 302.5 );
    histEff[ matchPre + "_TurnOn_MHTgt100GeV" ] = energySumsSubDir.make<TEfficiency>( matchPre + "_TurnOn_MHTgt100GeV", 
										     "#slash_{H}_{T} > 100 GeV;Offline #slash_{H}_{T} (GeV);Efficiency",
										     61, -2.5, 302.5 );
    histEff[ matchPre + "_TurnOn_MHTgt150GeV" ] = energySumsSubDir.make<TEfficiency>( matchPre + "_TurnOn_MHTgt150GeV", 
										     "#slash_{H}_{T} > 150 GeV;Offline #slash_{H}_{T} (GeV);Efficiency",
										     61, -2.5, 302.5 );

    // ------------------------------------------------------------




    booker.book2DTProf( matchPre + "_DeltaPTRel_vs_L1PT" , matchSubDir, "#Deltap_{T} vs L1 p_{T};L1 p_{T} (GeV);#Deltap_{T}", pL1PTLarge, pDeltaPTRel);
//     booker.book2DTProf( matchPre + "_DeltaPhi_vs_DeltaEta" , matchSubDir, "#Delta#phi vs #Delta#eta;#Delta#eta;#Delta#phi", pDeltaEta, pDeltaPhi);
//     // Online-offline quantity correlations
//     booker.book2DTProf( matchPre + "_L1Eta_vs_OffEta" , matchSubDir, "L1 #eta vs Offline #eta;Offline #eta;L1 #eta", pOffEta, pOffEta);
//     booker.book2DTProf( matchPre + "_L1Phi_vs_OffPhi" , matchSubDir, "L1 #phi vs Offline #phi;Offline #phi;L1 #phi", pOffPhi, pOffPhi);

    
    // pT Correlations
    booker.book2DTProf( matchPre + "_L1PT_vs_OffPT" ,     matchSubDir, "L1 P_{T} vs Offline P_{T};Offline p_{T} (GeV);L1 P_{T} (GeV)", pOffPT, pOffPT);
    booker.book2DTProf( matchPre + "_OffPT_vs_L1PT" ,     matchSubDir, "Offline P_{T} vs L1 P_{T};L1 P_{T} (GeV);Offline p_{T} (GeV)", pOffPT, pOffPT);
    booker.book2DTProf( matchPre + "_DeltaPTRel_vs_OffPT"  , matchSubDir, "#DeltaP_{T}^{Rel} vs Offline P_{T};Offline p_{T} (GeV);#DeltaP_{T}^{Rel}", pOffPT, pDeltaPTRel);     
    booker.book2DTProf( matchPre + "_DeltaEta_vs_OffPT" , matchSubDir, "#Delta#eta vs Offline P_{T};Offline p_{T} (GeV);#Delta#eta",   pOffPT, pDeltaEta);
    booker.book2DTProf( matchPre + "_DeltaPhi_vs_OffPT" , matchSubDir, "#Delta#phi vs Offline P_{T};Offline p_{T} (GeV);#Delta#phi",   pOffPT, pDeltaPhi);
    // Eta Correlations
    booker.book2DTProf( matchPre + "_L1PT_vs_Eta"     , matchSubDir, "L1 P_{T} vs #eta;L1 #eta;L1 P_{T} (GeV)", pOffEta, pL1PTLarge);
    booker.book2DTProf( matchPre + "_OffPT_vs_Eta"     , matchSubDir, "Offline P_{T} vs #eta;Offline #eta;Offline P_{T} (GeV)", pOffEta, pOffPT);
    booker.book2DTProf( matchPre + "_DeltaPTRel_vs_Eta"  , matchSubDir, "#DeltaP_{T}^{Rel} vs #eta;L1 #eta;#DeltaP_{T}^{Rel}", pOffEta, pDeltaPTRel);
    booker.book2DTProf( matchPre + "_DeltaEta_vs_Eta" , matchSubDir, "#Delta#eta vs #eta;L1 #eta;#Delta#eta",   pOffEta, pDeltaEta);
    booker.book2DTProf( matchPre + "_DeltaPhi_vs_Eta" , matchSubDir, "#Delta#phi vs #eta;L1 #eta;#Delta#phi",   pOffEta, pDeltaPhi);
    // Phi Correlations
    booker.book2DTProf( matchPre + "_L1PT_vs_Phi"     , matchSubDir, "L1 P_{T} vs #phi;L1 #phi;L1 P_{T} (GeV)", pOffPhi, pL1PTLarge);
    booker.book2DTProf( matchPre + "_DeltaPTRel_vs_Phi"  , matchSubDir, "#DeltaP_{T}^{Rel} vs #phi;L1 #phi;#DeltaP_{T}^{Rel}", pOffPhi, pDeltaPTRel);
    booker.book2DTProf( matchPre + "_DeltaEta_vs_Phi" , matchSubDir, "#Delta#eta vs #phi;L1 #phi;#Delta#eta",   pOffPhi, pDeltaEta);
    booker.book2DTProf( matchPre + "_DeltaPhi_vs_Phi" , matchSubDir, "#Delta#phi vs #phi;L1 #phi;#Delta#phi",   pOffPhi, pDeltaPhi);
    // NVTX Correlations
    booker.book2DTProf( matchPre + "_L1PT_vs_NVTX"     , matchSubDir, "L1 P_{T} vs N_{VTX};N_{VTX};L1 P_{T} (GeV)", pPU, pL1PTLarge);
    booker.book2DTProf( matchPre + "_DeltaPTRel_vs_NVTX"  , matchSubDir, "#DeltaP_{T}^{Rel} vs N_{VTX};N_{VTX};#DeltaP_{T}^{Rel}", pPU, pDeltaPTRel);
    booker.book2DTProf( matchPre + "_DeltaEta_vs_NVTX" , matchSubDir, "#Delta#eta vs N_{VTX};N_{VTX};#Delta#eta",   pPU, pDeltaEta);
    booker.book2DTProf( matchPre + "_DeltaPhi_vs_NVTX" , matchSubDir, "#Delta#phi vs N_{VTX};N_{VTX};#Delta#phi",   pPU, pDeltaPhi);


    // Matching efficiency
    histEff[ matchPre + "_matchingEffAllJet" ] = matchSubDir.make<TEfficiency>( matchPre + "_matchingEffAllJet", "Matching efficiency all jets;Offline p_{T} (GeV);Matching efficiency", 
									  40, 0, 200 );
    histEff[ matchPre + "_matchingEffLeadJet" ] = matchSubDir.make<TEfficiency>( matchPre + "_matchingEffLeadJet", "Matching efficiency lead jet;Offline p_{T} (GeV);Matching efficiency", 
									  40, 0, 200 );



    TFileDirectory turnOnSubDir = matchSubDir.mkdir( "TurnOns" );
    TFileDirectory rateSubDir   = matchSubDir.mkdir( "Rates" );
    TFileDirectory nomatchSubDir   = matchSubDir.mkdir( "NoMatch" );
    TFileDirectory rateStudySubDir   = matchSubDir.mkdir( "RateStudy" );


    // Turn-on curves
    for ( uint iOnPT = 0; iOnPT < onlineTOSingleJetPT.size(); ++iOnPT ){
      // Get the current offline PT cut 
      int onlinePT     = onlineTOSingleJetPT[ iOnPT ];
      TString onlinePTStr = Form("%d", onlinePT );
      
      // Turn-on
      histEff[ matchPre + "_PT_TurnOn_OnPT_gt_" + onlinePTStr ] = turnOnSubDir.make<TEfficiency>( matchPre + "_PT_TurnOn_OnPT_gt_" + onlinePTStr, 
												 "Differential p_{T} turn-on all jets p_{T} " + onlinePTStr + ";Offline p_{T} (GeV);Efficiency", 40, 0, 200 );
      histEff[ matchPre + "_NoMatch_PT_TurnOn_OnPT_gt_" + onlinePTStr ] = turnOnSubDir.make<TEfficiency>( matchPre + "_NoMatch_PT_TurnOn_OnPT_gt_" + 
													 onlinePTStr, "Differential p_{T} turn-on all jets p_{T} " + onlinePTStr + ";Offline p_{T} (GeV);Efficiency", 40, 0, 200 );






      // Specific jet distributions
      for (uint iJet = 0; iJet < jetIndexLabel.size(); ++iJet ){
	
	// jet labeling starting at 0. LeadJet = 0, ...                                                                  
	TString jetLabel = "_" + jetIndexLabel[ iJet ];
	TString jetText  = jetIndexLabel[ iJet ];

	// L1 to RECO matched - L1 jets are defined as lead, second, ...
	histEff[ matchPre + jetLabel + "_PT_TurnOn_OnPT_gt_"+ onlinePTStr ] = turnOnSubDir.make<TEfficiency>( matchPre + jetLabel + "_PT_TurnOn_OnPT_gt_" +
	       onlinePTStr, "Differential p_{T} turn-on " + jetText + " p_{T} " + onlinePTStr + ";Offline " + jetText + 
              " p_{T} (GeV);Efficiency", 40, 0, 200 );

	// L1 Unmatched jets - L1 jets are defined as lead, second, ...
	histEff[ matchPre + "_NoMatch" + jetLabel + "_PT_TurnOn_OnPT_gt_"+ onlinePTStr ] = turnOnSubDir.make<TEfficiency>( matchPre + 
           "_NoMatch_" + jetLabel + "_PT_TurnOn_OnPT_gt_" + onlinePTStr, "Differential p_{T} turn-on " + jetText + " p_{T} " + 
           onlinePTStr + ";Offline " + jetText + " p_{T} (GeV);Efficiency", 40, 0, 200 );
      

	// Only book once....
	if (iOnPT == 0){



	  // Unmatched jets :(
	  booker.book1D( matchPre + jetLabel + "_NoMatch_PT", nomatchSubDir, "L1 " + jetText + " not matched p_{T} - "  + matchLab + 
			 ";L1 p_{T};Entries", pOffPT );

	  booker.book1D( matchPre + jetLabel + "_NoMatch_Eta", nomatchSubDir, "L1 " + jetText + " not matched #eta - "  + matchLab + 
			 ";L1 #eta;Entries", pOffEta );

	  booker.book1D( matchPre + jetLabel + "_NoMatch_Phi", nomatchSubDir, "L1 " + jetText + " not matched #phi - "  + matchLab + 
			 ";L1 #phi;Entries", pOffPhi );
	  
	  booker.book1D( matchPre + jetLabel + "_NoMatch_deltaPTRel", nomatchSubDir, "L1 " + jetText + " not matched #Deltap_{T}^{rel} with " 
			 + jetText + " ak5 - "  + matchLab + ";#Deltap_{T}^{rel};Entries", pDeltaPTRel );

	  booker.book1D( matchPre + jetLabel + "_NoMatch_deltaEta", nomatchSubDir, "L1 " + jetText + " not matched #Delta#eta with " 
			 + jetText + " ak5 - " + matchLab + ";#Delta#eta;Entries", pDeltaEta );

	  booker.book1D( matchPre + jetLabel + "_NoMatch_deltaPhi", nomatchSubDir, "L1 " + jetText + " not matched #Delta#phi with " 
			 + jetText + " ak5 - " + matchLab + ";#Delta#phi;Entries", pDeltaPhi );

	  booker.book1D( matchPre + jetLabel + "_NoMatch_deltaR", nomatchSubDir, "L1 " + jetText + " not matched #DeltaR  with " 
			 + jetText + " ak5 - " + matchLab + ";#DeltaR;Entries", pDeltaR );

	  booker.book2D( matchPre + jetLabel + "_NoMatch_deltaPTRel_vs_deltaR", nomatchSubDir, "L1 " + jetText + 
			 " not matched #Deltap_{T}^{rel} vs #DeltaR with " + jetText + " ak5 - " 
			 + matchLab + ";#DeltaR;#Deltap_{T}^{rel}", pDeltaR, pDeltaPTRel );



	  
	  booker.book1D( matchPre + "_NoMatch" + jetLabel + "_L1PT", rateSubDir, "L1 " + jetText + " p_{T} (No match) - "  + matchLab + ";L1 p_{T} threshold;Entries", pL1PT );
	  booker.book1D( matchPre + "_NoMatch" + jetLabel + "_L1PT_Cumul", rateSubDir, "Events triggered as function of L1 " + jetText + " p_{T} threshold (No match) - "  + matchLab + ";L1 p_{T} threshold;Triggered", pL1PT );

	  // New rate studies
	  // --------------------------------------------------

	  for (uint iMatch = 0; iMatch < 2; ++iMatch){

	    TString matchType;
	    if (iMatch == 0){ matchType = "_NoMatch"; }
	    else            { matchType = "_Matched"; }

	  histEff[ matchPre + matchType + "_ak5ge30_Eff_vs_L1Thresh" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + 
														       "_ak5ge30_Eff_vs_L1Thresh"
		   + jetLabel, "Trigger efficiency (ak5 p_{T} > 30 GeV) vs " + jetText + " L1 p_{T} threshold;L1 p_{T} threshold(GeV);Efficiency", 
														   41, -2.5, 202.5 );
	  histEff[ matchPre +  matchType + "_ak5ge40_Eff_vs_L1Thresh" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + 
															"_ak5ge40_Eff_vs_L1Thresh"
		   + jetLabel, "Trigger efficiency (ak5 p_{T} > 40 GeV) vs " + jetText + " L1 p_{T} threshold;L1 p_{T} threshold(GeV);Efficiency", 
														   41, -2.5, 202.5 );
	  histEff[ matchPre +  matchType + "_ak5ge50_Eff_vs_L1Thresh" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + 
															"_ak5ge50_Eff_vs_L1Thresh"
		   + jetLabel, "Trigger efficiency (ak5 p_{T} > 50 GeV) vs " + jetText + " L1 p_{T} threshold;L1 p_{T} threshold(GeV);Efficiency", 
														   41, -2.5, 202.5 );

	  histEff[ matchPre +  matchType + "_ak5ge30_Eff_vs_OffEta" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + "_ak5ge30_Eff_vs_OffEta"
		   + jetLabel, "Trigger efficiency (L1 p_{T} > 12 GeV, given ak5 p_{T} > 30 GeV) vs " + jetText + " Offline #eta;Offline #eta;Efficiency", 64, -3.2, 3.2 );
	  histEff[ matchPre +  matchType + "_ak5ge40_Eff_vs_OffEta" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + "_ak5ge40_Eff_vs_OffEta"
		   + jetLabel, "Trigger efficiency (L1 p_{T} > 12 GeV, given ak5 p_{T} > 40 GeV) vs " + jetText + " Offline #eta;Offline #eta;Efficiency", 64, -3.2, 3.2 );
	  histEff[ matchPre +  matchType + "_ak5ge50_Eff_vs_OffEta" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + "_ak5ge50_Eff_vs_OffEta"
		   + jetLabel, "Trigger efficiency (L1 p_{T} > 12 GeV, given ak5 p_{T} > 50 GeV) vs " + jetText + " Offline #eta;Offline #eta;Efficiency", 64, -3.2, 3.2 );
	  histEff[ matchPre +  matchType + "_ak5ge30_Eff_vs_OffEtaRegion" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + "_ak5ge30_Eff_vs_OffEtaRegion"
		   + jetLabel, "Trigger efficiency (L1 p_{T} > 12 GeV, given ak5 p_{T} > 30 GeV) vs " + jetText + " Offline #eta;Offline #eta;Efficiency", regionEtaBin, regionEtaBins );
	  histEff[ matchPre +  matchType + "_ak5ge40_Eff_vs_OffEtaRegion" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + "_ak5ge40_Eff_vs_OffEtaRegion"
		   + jetLabel, "Trigger efficiency (L1 p_{T} > 12 GeV, given ak5 p_{T} > 40 GeV) vs " + jetText + " Offline #eta;Offline #eta;Efficiency", regionEtaBin, regionEtaBins );
	  histEff[ matchPre +  matchType + "_ak5ge50_Eff_vs_OffEtaRegion" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + "_ak5ge50_Eff_vs_OffEtaRegion"
		   + jetLabel, "Trigger efficiency (L1 p_{T} > 12 GeV, given ak5 p_{T} > 50 GeV) vs " + jetText + " Offline #eta;Offline #eta;Efficiency", regionEtaBin, regionEtaBins );

	  histEff[ matchPre +  matchType + "_ak5ge30_Eff_vs_NVTX" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + "_ak5ge30_Eff_vs_NVTX"
		   + jetLabel, "Trigger efficiency (L1 p_{T} > 12 GeV, given ak5 p_{T} > 30 GeV) vs " + jetText + " N_{VTX};N_{VTX};Efficiency", 40, 0, 80 );
	  histEff[ matchPre +  matchType + "_ak5ge40_Eff_vs_NVTX" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + "_ak5ge40_Eff_vs_NVTX"
		   + jetLabel, "Trigger efficiency (L1 p_{T} > 12 GeV, given ak5 p_{T} > 40 GeV) vs " + jetText + " N_{VTX};N_{VTX};Efficiency", 40, 0, 80 );
	  histEff[ matchPre +  matchType + "_ak5ge50_Eff_vs_NVTX" + jetLabel ] = rateStudySubDir.make<TEfficiency>( matchPre +  matchType + "_ak5ge50_Eff_vs_NVTX"
		   + jetLabel, "Trigger efficiency (L1 p_{T} > 12 GeV, given ak5 p_{T} > 50 GeV) vs " + jetText + " N_{VTX};N_{VTX};Efficiency", 40, 0, 80 );



	  hist1D[ matchPre +  matchType + "" + jetLabel + "_TriggeredProto" ] = rateStudySubDir.make<TH1D>(  matchPre +  matchType + "" + jetLabel + "_TriggeredProto",
													jetText + " p_{T} (GeV) that triggered event (No match) - "
													+ matchLab + ";L1 p_{T};Entries", 41, -2.5, 202.5 );
	  
	  hist1D[ matchPre +  matchType + "" + jetLabel + "_Triggered" ] = rateStudySubDir.make<TH1D>(  matchPre +  matchType + "" + jetLabel + "_Triggered",
												   "Events triggered as a function of " + jetText + 
												   " p_{T} threshold (No match) - "  + matchLab + 
												   ";L1 p_{T} threshold;Events triggered", 41, -2.5, 202.5 );

	  } // End matchtype
	} // End book once	
      }
      
    }
  } // End matched jets
}


JetHist::~JetHist()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// **********************************************************************
// *                              analyze()                             *
// *                                                                    *
// *                                                                    *  
// **********************************************************************        
//JetHist::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
void
JetHist::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   bool evValid = true;

   // Get the event number
   Eventnr = iEvent.id().event();

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



   // ********************************************************************************
   // *                                    Jets                                      *
   // ********************************************************************************

  // ************************************************************
  // *                          FastJet                         *
  // ************************************************************


#ifdef FASTJET

  edm::Handle<l1extra::L1JetParticleCollection> TTAk5PrePUSRaw_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("TTAk5L1Jet"), TTAk5PrePUSRaw_L1Jet); 
  if(!TTAk5PrePUSRaw_L1Jet.isValid()){
    evValid = false;   
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("TTAk5L1Jet") << std::endl;
  }

#endif

 
  // ************************************************************
  // *                       Current Jets                       *
  // ************************************************************
  SUBPRINT("Current Jets")
    
    // // Get L1 jets to be used, currently Cental and Tau jets
    //    std::vector <InputTag> l1extraparticles = conf_.getParameter< std::vector < InputTag > >("extrajet");
   //  edm::Handle<l1extra::L1JetParticleCollection> Current_L1Jet;
//      // Loop through Central and tau jets
//     for (uint i = 0; i < l1extraparticles.size(); ++i){
    
//       iEvent.getByLabel(l1extraparticles[i], Current_L1Jet );
//       if(!Current_L1Jet.isValid()){
// 	evValid = false;
// 	edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("extrajet") << std::endl;
//       }
//    }

  // ************************************************************
  // *                        Tower Jets                        *
  // ************************************************************
  SUBPRINT("Tower Jets")
  edm::Handle<l1extra::L1JetParticleCollection> PrePUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSTowerJetL1Jet"), PrePUS_L1Jet); 
  if(!PrePUS_L1Jet.isValid()){
    evValid = false;   
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSTowerJetL1Jet") << std::endl;
  }
  edm::Handle<l1extra::L1JetParticleCollection> PUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUSTowerJetL1Jet"), PUS_L1Jet); 
  if(!PUS_L1Jet.isValid()){
    evValid = false;   
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUSTowerJetL1Jet") << std::endl;
  }
  edm::Handle<l1extra::L1JetParticleCollection> LPUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("LocalPUSTowerJetL1Jet"), LPUS_L1Jet); 
  if(!LPUS_L1Jet.isValid()){
    evValid = false;   
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("LocalPUSTowerJetL1Jet") << std::endl;
  }

  edm::Handle<l1extra::L1JetParticleCollection> CalibPrePUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibratedPrePUSak5PUSTowerJetL1Jet"), CalibPrePUS_L1Jet); 
  if(!CalibPrePUS_L1Jet.isValid()){
    evValid = false;   
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibratedPrePUSak5PUSTowerJetL1Jet") << std::endl;
  }
  edm::Handle<l1extra::L1JetParticleCollection> CalibPUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibratedPUSak5PUSTowerJetL1Jet"), CalibPUS_L1Jet); 
  if(!CalibPUS_L1Jet.isValid()){
    evValid = false;   
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibratedPUSak5PUSTowerJetL1Jet") << std::endl;
  }
  edm::Handle<l1extra::L1JetParticleCollection> CalibLPUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibratedLPUSak5PUSTowerJetL1Jet"), CalibLPUS_L1Jet); 
  if(!CalibLPUS_L1Jet.isValid()){
    evValid = false;   
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibratedLPUSak5PUSTowerJetL1Jet") << std::endl;
  }

  // ************************************************************
  // *                        RECO(Gen) Jets                         *
  // ************************************************************
  SUBPRINT("RECO Jets")
  edm::Handle<l1extra::L1JetParticleCollection> Ak5PUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("Ak5CaloJetL1Jet"), Ak5PUS_L1Jet);
  if(!Ak5PUS_L1Jet.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("Ak5CaloJetL1Jet") << std::endl;
    evValid = false;
  }

   evValid=true;
//  if (evValid){
   if(evValid){

   // ****************************************************************************************************
   // *                                        Primary vertices                                          *
   // ****************************************************************************************************
   TString prefix = "Vertex";

   
   // Extract number of reconstructed ak5 primary vertices
   NVTX       = *vertices;


    // ****************************************************************************************************
    // *                                   L1-Offline lead jet matching                                   *
    // ****************************************************************************************************

//#ifdef L1_JET_PLOTS
    PRINT("L1 plots")
 

       fillL1Histograms( PrePUS_L1Jet,   "PrePUS" );
       fillL1Histograms( PUS_L1Jet,      "PUS" );
       fillL1Histograms( LPUS_L1Jet,      "LPUS" );
       fillL1Histograms( CalibPrePUS_L1Jet,   "CalibPrePUS" );
       fillL1Histograms( CalibPUS_L1Jet,      "CalibPUS" );
       fillL1Histograms( CalibLPUS_L1Jet,      "CalibLPUS" );
       fillL1Histograms( Ak5PUS_L1Jet, "Ak5" );
//#endif


#ifdef FASTJET
      fillL1Histograms( TTAk5PrePUSRaw_L1Jet, "TTAk5PrePUSRaw");
#endif

      // ----------------------------------------------------------------------------------------------------                                                                    




      // ****************************************************************************************************
      // *                                            Calibration                                           *
      // ****************************************************************************************************

#ifdef FASTJET
      calibrateL1RECO(  TTAk5PrePUSRaw_L1Jet, Ak5PrePUSRaw_L1Jet, "Calibration_TTAk5PrePUSRaw_ak5PrePUSRaw", 1, 1, 2.5);
#endif

      // Closure test - hand-calibrated ak5PUS jets
      //      calibrateL1RECO(  CalibAk5PUSRawAk5PUS_L1Jet, Ak5PUS_L1Jet, "Calibration_CalibAk5PUSRawAk5PUS_ak5PUS", 1, 1, 2.5);



      // Closure test - use these known collections to test and optimise the calibration method
      //      calibrateL1RECO(  Ak5PUSRaw_L1Jet, Ak5PUS_L1Jet, "Calibration_ak5PUSRaw_ak5PUS", 1, 1, 2.5);

      //      calibrateL1RECO(  Ak5PrePUSRaw_L1Jet, Ak5PUSRaw_L1Jet, "Calibration_ak5PrePUSRaw_ak5PUSRaw", 1, 1, 2.5);
      //      calibrateL1RECO(  Ak5PrePUSRaw_L1Jet, Ak5PUS_L1Jet,    "Calibration_ak5PrePUSRaw_ak5PUS", 1, 1, 2.5);
      //      calibrateL1RECO(  Ak5PrePUSRaw_L1Jet, Ak5PrePUS_L1Jet, "Calibration_ak5PrePUSRaw_ak5PrePUS", 1, 1, 2.5);





      // |eta| < 3
      // PrePUS-ak5PUS
      //      calibrateL1RECO(  PrePUS_L1Jet, Ak5PUS_L1Jet, "Calibration_PrePUS_ak5PUS_EtaLt3", 1, 1, 3.0);

      // ************************************************************
      // *                     NVTX binned LUTs                     *
      // ************************************************************
//       // Really inefficient!!!
//       if (NVTX < 15){
// 	calibrateL1RECO(  PrePUS_L1Jet, Ak5PUS_L1Jet, "Calibration_PrePUS_ak5PUS_NVTXLt15", 1, 1, 2.5);
//       }
//       else if (NVTX < 25){
// 	calibrateL1RECO(  PrePUS_L1Jet, Ak5PUS_L1Jet, "Calibration_PrePUS_ak5PUS_NVTXLt25", 1, 1, 2.5);
//       }
//       else if (NVTX < 50){
// 	calibrateL1RECO(  PrePUS_L1Jet, Ak5PUS_L1Jet, "Calibration_PrePUS_ak5PUS_NVTXLt50", 1, 1, 2.5);
//       }
      
     



//#ifdef CALIBRATION


      
      // L1GCT-ak5PUS
      //calibrateL1RECO( Current_L1Jet,           Ak5PUS_L1Jet, "Calibration_Curr_ak5PUS",        1, 1, 2.5);
      // L1UncalibGCT-ak5PUS
      //calibrateL1RECO( CurrentUncalib_L1Jet,    Ak5PUS_L1Jet, "Calibration_UncalibCurr_ak5PUS", 1, 1, 2.5);
      // PrePUS
      // --------------------------------------------------

      // PrePUS-ak5PUS
      calibrateL1RECO(  PrePUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_PrePUS_ak5PUS",            1, jetRankComparison, 2.5);
      calibrateL1RECO(  PrePUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_PrePUS_ak5PUS_3Jets",            1, 3, 2.5);
      calibrateL1RECO(  PrePUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_PrePUS_ak5PUS_AllJets",            1, 999, 2.5);

      calibrateL1RECO(  CalibPrePUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_CalibPrePUS_ak5PUS",            1, jetRankComparison, 2.5);
      calibrateL1RECO(  CalibPrePUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_CalibPrePUS_ak5PUS_3Jets",            1, 3, 2.5);
      calibrateL1RECO(  CalibPrePUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_CalibPrePUS_ak5PUS_AllJets",            1, 999, 2.5);

      calibrateL1RECO(  PUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_PUS_ak5PUS",            1, jetRankComparison, 2.5);
      calibrateL1RECO(  PUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_PUS_ak5PUS_3Jets",            1, 3, 2.5);
      calibrateL1RECO(  PUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_PUS_ak5PUS_AllJets",            1, 999, 2.5);

      calibrateL1RECO(  LPUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_LPUS_ak5PUS",            1, jetRankComparison, 2.5);
      calibrateL1RECO(  LPUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_LPUS_ak5PUS_3Jets",            1, 3, 2.5);
      calibrateL1RECO(  LPUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_LPUS_ak5PUS_AllJets",            1, 999, 2.5);

      calibrateL1RECO(  CalibPUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_CalibPUS_ak5PUS",            1, jetRankComparison, 2.5);
      calibrateL1RECO(  CalibPUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_CalibPUS_ak5PUS_3Jets",            1, 3, 2.5);
      calibrateL1RECO(  CalibPUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_CalibPUS_ak5PUS_AllJets",            1, 999, 2.5);
      
      calibrateL1RECO(  CalibLPUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_CalibLPUS_ak5PUS",            1, jetRankComparison, 2.5);
      calibrateL1RECO(  CalibLPUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_CalibLPUS_ak5PUS_3Jets",            1, 3, 2.5);
      calibrateL1RECO(  CalibLPUS_L1Jet,                  Ak5PUS_L1Jet, "Calibration_CalibLPUS_ak5PUS_AllJets",            1, 999, 2.5);

//#endif


      // ----------------------------------------------------------------------------------------------------                                                                    


      // ****************************************************************************************************
      // *                                            Benchmarking                                          *
      // ****************************************************************************************************

//#ifdef BENCHMARKING

      // Closure test - hand-calibrated ak5PUS jets 
      //      benchmarkJets(  CalibAk5PUSRawAk5PUS_L1Jet, Ak5PUS_L1Jet, "CalibAk5PUSRawAk5PUS_ak5PUS", 1, 1, 4, 2.5, 200);




      // Current
      // --------------------------------------------------
      
      // Curr-ak5PUS
      //benchmarkJets( Current_L1Jet, Ak5PUS_L1Jet, "Curr_ak5PUS",   1, 12, 4, 2.5, 200 );
      // CurrUncalib-ak5PUS
      //benchmarkJets( CurrentUncalib_L1Jet, Ak5PUS_L1Jet, "CurrUncalib_ak5PUS",   1, 12, 4, 2.5, 200 );


      

      // PrePUS
      // --------------------------------------------------

      // PrePUS-ak5PUS
      benchmarkJets(  PrePUS_L1Jet,                 Ak5PUS_L1Jet, "PrePUS_ak5PUS",            1, 1, 999, 2.5, 1000);
      benchmarkJets( CalibPrePUS_L1Jet, Ak5PUS_L1Jet, "CalibPrePUS_ak5PUS", 1,1,999,2.5,1000);
      benchmarkJets(  PUS_L1Jet,                 Ak5PUS_L1Jet, "PUS_ak5PUS",            1, 1, 999, 2.5, 1000);
      benchmarkJets( CalibPUS_L1Jet, Ak5PUS_L1Jet, "CalibPUS_ak5PUS", 1,1,999,2.5,1000);
      benchmarkJets(  LPUS_L1Jet,                 Ak5PUS_L1Jet, "PUS_ak5PUS",            1, 1, 999, 2.5, 1000);
      benchmarkJets( CalibLPUS_L1Jet, Ak5PUS_L1Jet, "CalibPUS_ak5PUS", 1,1,999,2.5,1000);


//#endif

  }

   QUIT("Successfully completed analysis of single event")
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetHist::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetHist::endJob() {



  // Calculate the cumulative trigger rate for specific jet pT thresholds

  for (unsigned int iMatch = 0; iMatch < matchPrefix.size(); iMatch++){ 
 
    TString matchPre = matchPrefix[iMatch]; 
  
    
    // Specific jet distributions 
    for (uint iJet = 0; iJet < jetIndexLabel.size(); ++iJet ){
      
      // jet labeling starting at 0. LeadJet = 0, ... 
      TString jetLabel = "_" + jetIndexLabel[ iJet ];
      
      
       reverseCumulative( hist1D[ matchPre + "_NoMatch" + jetLabel + "_L1PT" ], hist1D[ matchPre + "_NoMatch" + jetLabel + "_L1PT_Cumul" ] );

      reverseCumulative( hist1D[ matchPre + "_NoMatch" + jetLabel + "_TriggeredProto" ], hist1D[ matchPre + "_NoMatch" + jetLabel + "_Triggered" ] );


//       // Calculate cumlative triggered events
//       TH1* h = (TH1*)hist1D[ matchPre + "_NoMatch" + jetLabel + "_L1PT" ]->Clone();
      
//       int nBins    = h->GetNbinsX();
//       int cumulBin = h->GetBinContent( nBins + 1 ); // Initialise with overflow bin
      
//       for ( int iBin = nBins; iBin > -1; --iBin ){
      
// 	cumulBin += h->GetBinContent( iBin );
// 	hist1D[ matchPre + "_NoMatch" + jetLabel + "_L1PT_Cumul" ]->SetBinContent( iBin, cumulBin );
	
//       }
      
    }

  }    
  
}

// ------------ method called when starting to processes a run  ------------
void 
JetHist::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetHist::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetHist::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetHist::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetHist::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}




// Calculate the reverse cumulative distribution
void
JetHist::reverseCumulative( TH1* histogram, TH1* rCumulHist ){

      
      // Calculate cumlative triggered events
  //      rCumulHist = (TH1*)histogram->Clone();
      
      int nBins    = histogram->GetNbinsX();
      int cumulBin = histogram->GetBinContent( nBins + 1 ); // Initialise with overflow bin
      
      for ( int iBin = nBins; iBin > -1; --iBin ){
      
	cumulBin += histogram->GetBinContent( iBin );
	rCumulHist->SetBinContent( iBin, cumulBin );
	
      }

}    




// ********************************************************************************
// *                  Fill all jet distributions for online jets                  *
// ******************************************************************************** 
void 
JetHist::fillL1Histograms( edm::Handle<l1extra::L1JetParticleCollection> const& L1Jets, TString prefix )
{
  // Event quantities
  int NJets = std::distance(L1Jets->begin(), L1Jets->end());
  
  // ***********************************************************************
  // *                  General jet distributions                          *
  // ***********************************************************************  
  
  hist1D[prefix + "_NJets"]->Fill( NJets );


  int jetIndex = 0;
  
  for ( l1extra::L1JetParticleCollection::const_iterator L1_It = L1Jets->begin(); L1_It != L1Jets->end(); ++L1_It ){
    
    
    double l1Pt         = L1_It->p4().Pt();
    double l1Eta        = L1_It->p4().Eta();
    double l1Phi        = L1_It->p4().Phi();

    

    // ***********************************************************************
    // *                 'ith' jet distributions                             *
    // ***********************************************************************
          
     // Jet triggers: First, second, third and fourth jet
     if( jetIndex < 4 ){
       TString jetNum = Form("%d", jetIndex + 1);

       // PT, eta, phi    
       hist1D[prefix + "_Jet_" + jetNum + "_PT"] ->Fill( l1Pt );
       hist1D[prefix + "_Jet_" + jetNum + "_Eta"]->Fill( l1Eta );
       hist1D[prefix + "_Jet_" + jetNum + "_Phi"]->Fill( l1Phi );

       // Phi vs Eta
       hist2D[prefix + "_Jet_" + jetNum + "_Phi_vs_Eta"]->Fill( l1Eta, l1Phi ); 
       
     }    
     

     // ***********************************************************************
     // *                  All jet distributions                              *
     // ***********************************************************************
     
     // PT, eta, phi    
     hist1D[prefix + "_PT"]        ->Fill( l1Pt  );
     hist1D[prefix + "_Eta"]       ->Fill( l1Eta );
     hist1D[prefix + "_Phi"]       ->Fill( l1Phi );
     

     // TEMPORARY UNTIL WE SORT OUT TH2 MAP SO WE CAN USE l1Pt WEIGHTING
     for (int iPt = 0;iPt < l1Pt;++iPt){
       // Phi vs Eta
       hist2D[prefix + "_Phi_vs_Eta"]->Fill( l1Eta, l1Phi ); 
       //       hist2D[prefix + "_Phi_vs_Eta"]->Fill( l1Eta, l1Phi, l1Pt ); 
     }

     hist2D[prefix + "_Eta_vs_Pt"]->Fill( l1Pt, l1Eta );
     hist2D[prefix + "_Phi_vs_Pt"]->Fill( l1Pt, l1Phi );



     
     // Increment the jet counter
     jetIndex++;
  } // End iterator loop


}









// ****************************************************************************************************
// *                                    Match Online to Offline Jets                                  *
// **************************************************************************************************** 
//
// QUICK AND DIRTY: REWRITE THIS CODE WHEN YOU HAVE MORE TIME...
//
//
void 
JetHist::calibrateL1RECO( edm::Handle<l1extra::L1JetParticleCollection> const& L1Jets, edm::Handle<l1extra::L1JetParticleCollection> const& ak5Jets,
			  TString prefix, double ak5PTthreshold, int maxMatchedJets, double ak5EtaRange){



// void 
// JetHist::calibrateL1RECO( edm::Handle<l1slhc::L1TowerJetCollection> const& L1Jets, edm::Handle<l1extra::L1JetParticleCollection> const& ak5Jets,
// 			  TString prefix, double ak5PTthreshold, int maxMatchedJets, double ak5EtaRange ){
			  //			  TString prefix, double ak5PTthreshold, double L1PTthreshold, int maxMatchedJets, double ak5EtaRange ){

  PRINT("Online-Offline Jet matching")

  // Find which histogram prefix to use
  int preIndex = -1;
  for (unsigned int iMatch = 0; iMatch < calibPrefix.size(); iMatch++){
    if (prefix == calibPrefix[iMatch]){
      preIndex = iMatch;
      break;
    }
  }
  if (preIndex == -1){
    std::cout << "ERROR - Prefix '" << prefix << "' was not found.\n";
    return;
  } 
  TString matchPre = calibPrefix[preIndex];
  TString matchLab = calibLabel[preIndex];

  SUBPRINT("Extracted prefix")
  //      ********************************************************************************
  //      *                          Online-Offline jet matching                         *
  //      ********************************************************************************
  //
  std::vector<TLorentzVector> ak5Vec, L1Vec;
  //  std::vector<int>            L1iEtaVec;

  //           *********************************************
  //           *          Fill ak5 TLorentzVector          *
  //           *********************************************

  //  for ( l1extra::L1JetParticleCollection::const_iterator ak5_It = ak5Jets->begin(); ak5_It != ak5Jets->end(); ++ak5_It ){
  for ( l1extra::L1JetParticleCollection::const_iterator ak5_It = ak5Jets->begin(); ak5_It != ak5Jets->end(); ++ak5_It ){

    double ak5Pt         = ak5_It->p4().Pt();
    double ak5Eta        = ak5_It->p4().Eta();
    double ak5Phi        = ak5_It->p4().Phi();
    double ak5M          = ak5_It->p4().M();

    // ========================================
    // Jet quality cuts
    // ========================================

    if( fabs( ak5Eta ) > ak5EtaRange ) continue;
    if( ak5Pt < ak5PTthreshold )       continue;

    // ========================================

    TLorentzVector ak5Jet;
    ak5Jet.SetPtEtaPhiM( ak5Pt, ak5Eta, ak5Phi, ak5M );
    ak5Vec.push_back( ak5Jet );

  }

  // <TEMPORARY> - Calibrated ak5 jets were previously not sorted, fixed only in new ntuples. Remove when transition to new ntuples is complete.
     // Sort calibrated jets
  //     std::sort( ak5Vec.begin(), ak5Vec.end(), TLorentzVectorRankDescending );
  // </TEMPORARY>


  //           *********************************************
  //           *          Fill L1 TLorentzVector           *
  //           *********************************************

  //  for (l1slhc::L1TowerJetCollection::const_iterator L1_It = L1Jets->begin(); L1_It != L1Jets->end(); ++L1_It ){
  for ( l1extra::L1JetParticleCollection::const_iterator L1_It = L1Jets->begin(); L1_It != L1Jets->end(); ++L1_It ){

    double L1Pt         = L1_It->p4().Pt();
    double L1Eta        = L1_It->p4().Eta();
    double L1Phi        = L1_It->p4().Phi();
    double L1M          = L1_It->p4().M();


    TLorentzVector L1Jet;
    L1Jet.SetPtEtaPhiM( L1Pt, L1Eta, L1Phi, L1M );
    L1Vec.push_back( L1Jet );

  }



  // ------------------------------ Perform jet matching ------------------------------


  SUBPRINT("Extracted jets")

    // THIS IS REALLY INEFFICIENT. WHEN YOU HAVE TIME WRITE A CLASS THAT RETURNS std::vector< std::pair< TLorentz*, TLorentz*> >
    // WITH SORT ALGORITHMS AS AN OPTION

    // Perform deltaR matching
    
    std::vector<pair_info> pairs = make_pairs( ak5Vec, L1Vec );
    std::sort(pairs.begin(), pairs.end(), sortDR);
    std::vector<int> L1MatchedIndex = analyse_pairs(pairs, L1Vec.size(), maxDeltaR);


  // ____________________
  // ALL JET MATCHING
  // ____________________
  SUBPRINT("All jet matching")
  int L1JetsMatched    = 0;
  int L1JetsNotMatched = 0;

  for(unsigned int i = 0; i < L1MatchedIndex.size(); i++) {

    int L1Index  = i;
    int ak5Index = L1MatchedIndex[i];

    // CHANGING DEFINITION - MAX JETS IN RANKED LIST i.e. 1 = leadJet only!!!!
    L1JetsMatched++;
    // Restrict to specified number of jets
    if ( L1JetsMatched > maxMatchedJets ){
      break;
    }


    if ( ak5Index == -1 ){
	//      std::cout << "No matched jet!!!\n";
      L1JetsNotMatched++;
      continue;
    }












    
//     // ----------------------------------------------------------------------------------------------------
//     // NVTX stuff
//     // ****************************************
//     TString NVTXStr;

//     if ( NVTX > nvtxBins[nvtxBins.size() - 1] ){ // NVTX exceeds last bin
//       int nvtxLow = nvtxBins[nvtxBins.size() - 1];
//       NVTXStr = "NVTX_" + TString(Form("%d", nvtxLow )) + "_to_inf";
//     }
//     else{ // NVTX is contained within one of the bins
//       for ( uint iNVTX = 0; iNVTX < nvtxBins.size() - 1; iNVTX++ ){
	
// 	int nvtxLow  = nvtxBins[ iNVTX ];
// 	int nvtxHigh = nvtxBins[ iNVTX + 1 ];
	
// 	if ( (NVTX >= nvtxLow) && (NVTX < nvtxHigh) ){
	  
// 	  // find which bin the nvtx value resides
// 	  NVTXStr = "NVTX_" + TString(Form("%d", nvtxLow )) +
// 	    "_to_" + TString(Form("%d", nvtxHigh ));
	  
// 	  break;
// 	}
//       }
//     }
//     // ----------------------------------------------------------------------------------------------------




//     // ----------------------------------------------------------------------------------------------------
//     // iEta stuff
//     // ************************************************************

//     // Get iEta
//     int iEta = L1iEtaVec[L1Index];
//     TString iEtaStr;
   
//     // iEta binned distributions    
//     // ****************************************
//     if (iEtaCalibBinWidth == 1){ // Individual iEta bins
//       iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEta)));
//     }
//     else{ // iEta bins of width greater than 1 

//       int iEtaBinAbs = (abs(iEta) - 1)/iEtaCalibBinWidth;
//       int sign = 1;
//       if (iEta < 0)
// 	sign = -1;
      
//       double iEtaBinLowAbs  = ( (iEtaBinAbs)*iEtaCalibBinWidth + 1 );
//       double iEtaBinHighAbs = ( (iEtaBinAbs + 1)*iEtaCalibBinWidth );

//       // Fix bin limits
//       if (iEtaBinLowAbs > 28)
// 	iEtaBinLowAbs = 28;
//       if (iEtaBinHighAbs > 28)
// 	iEtaBinHighAbs = 28;

//       double iEtaBinLow  = sign*iEtaBinLowAbs;
//       double iEtaBinHigh = sign*iEtaBinHighAbs;
      

//       if (iEtaBinLow > iEtaBinHigh)
// 	std::swap(iEtaBinLow, iEtaBinHigh);
      
//       iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEtaBinLow))) + "to" + TString(Form("%d",Int_t(iEtaBinHigh)));
//       //    std::cout << "\n\niEta = " << iEta << " is in the bin range:   " << iEtaBinLow << "to" << iEtaBinHigh << "\n";
//     }

    // ----------------------------------------------------------------------------------------------------


    //    TString iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEta)));
    
    double ak5Pt    = ak5Vec[ ak5Index ].Pt();
    //    double ak5Eta   = ak5Vec[ ak5Index ].Eta();
    //    double ak5Phi   = ak5Vec[ ak5Index ].Phi();

    double L1Pt     = L1Vec[ L1Index ].Pt();
    double L1Eta    = L1Vec[ L1Index ].Eta();
    //    double L1Phi    = L1Vec[ L1Index ].Phi();






//     // Restrict L1 PT
//     if ( L1Pt < L1PTthreshold){
//       continue;
//     }

    double deltaPt     = (L1Pt - ak5Pt);
    double deltaPtRel  = (L1Pt - ak5Pt)/ak5Pt; 
    double jetResponse = L1Pt/ak5Pt;
//     double deltaEta    = ak5Eta - L1Eta;
//     double deltaPhi    = ak5Vec[ ak5Index ].DeltaPhi( L1Vec[ L1Index ] );
    //    std::cout << "deltaPT = " << deltaPtRel << "\tdeltaEta = " << deltaEta << "\tdeltaPhi = " << deltaPhi << "\n";


    // ************************************************************
    // *                         Unbinned                         *
    // ************************************************************

//     // Resolution distributions
//     hist1D[matchPre + "_DeltaEta"]->Fill(deltaEta);
//     hist1D[matchPre + "_DeltaPhi"]->Fill(deltaPhi);
//     hist1D[matchPre + "_DeltaPTRel"] ->Fill(deltaPtRel); 

//     hist1D[matchPre + "_Eta"]     ->Fill(L1Eta);
//     hist1D[matchPre + "_Phi"]     ->Fill(L1Phi);
//     hist1D[matchPre + "_PT"]      ->Fill(L1Pt);

//     hist1D[matchPre + "_Ak5Eta"]  ->Fill(ak5Eta);
//     hist1D[matchPre + "_Ak5Phi"]  ->Fill(ak5Phi);
//     hist1D[matchPre + "_Ak5PT"]   ->Fill(ak5Pt);

    
//     // Online-offline quantity correlations
// //     hist2D[matchPre + "_L1Eta_vs_OffEta" ]         ->Fill(ak5Eta, L1Eta);       
// //     hist1D[matchPre + "_L1Eta_vs_OffEta_prof"]     ->Fill(ak5Eta, L1Eta);       
// //     hist2D[matchPre + "_L1Phi_vs_OffPhi" ]         ->Fill(ak5Phi, L1Phi);    
// //     hist1D[matchPre + "_L1Phi_vs_OffPhi_prof"]     ->Fill(ak5Phi, L1Phi);    
//     hist2D[matchPre + "_DeltaPhi_vs_DeltaEta" ]    ->Fill(deltaEta, deltaPhi);    
//     hist1D[matchPre + "_DeltaPhi_vs_DeltaEta_prof"]->Fill(deltaEta, deltaPhi);    
//     hist2D[matchPre + "_DeltaPTRel_vs_L1PT" ]         ->Fill(L1Pt, deltaPtRel);    
//     hist1D[matchPre + "_DeltaPTRel_vs_L1PT_prof"]     ->Fill(L1Pt, deltaPtRel);    

    // pT Correlations
    //    hist3D[matchPre + "_L1PT_vs_OffPT_vs_NVTX"]    ->Fill(ak5Pt, L1Pt, NVTX);
    hist2D[matchPre + "_L1PT_vs_OffPT"     ]       ->Fill(ak5Pt, L1Pt);
    hist1D[matchPre + "_L1PT_vs_OffPT_prof"]       ->Fill(ak5Pt, L1Pt);
    hist2D[matchPre + "_OffPT_vs_L1PT"     ]       ->Fill(L1Pt, ak5Pt);
    hist1D[matchPre + "_OffPT_vs_L1PT_prof"]       ->Fill(L1Pt, ak5Pt);
    hist2D[matchPre + "_DeltaPTRel_vs_OffPT"  ]    ->Fill(ak5Pt, deltaPtRel);
    hist1D[matchPre + "_DeltaPTRel_vs_OffPT_prof"] ->Fill(ak5Pt, deltaPtRel);

    hist2D[matchPre + "_DeltaPTRel_vs_DeltaPT"  ]    ->Fill(deltaPt, deltaPtRel);
    hist1D[matchPre + "_DeltaPTRel_vs_DeltaPT_prof"] ->Fill(deltaPt, deltaPtRel);

    hist2D[matchPre + "_DeltaPT_vs_OffPT"  ]       ->Fill(ak5Pt, deltaPt);
    hist1D[matchPre + "_DeltaPT_vs_OffPT_prof"]    ->Fill(ak5Pt, deltaPt);
    hist1D[matchPre + "_JetResponse"]              ->Fill(jetResponse);
    hist2D[matchPre + "_JetResponse_vs_OffPT"]     ->Fill(ak5Pt, jetResponse);
    hist1D[matchPre + "_JetResponse_vs_OffPT_prof"]->Fill(ak5Pt, jetResponse); 
    hist2D[matchPre + "_JetResponse_vs_L1PT"]     ->Fill(L1Pt, jetResponse);
    hist1D[matchPre + "_JetResponse_vs_L1PT_prof"]->Fill(L1Pt, jetResponse); 
    if ( L1Pt >= 12 ){
       hist2D[matchPre + "_JetResponse_vs_NVTX_L1Ptge12"]     ->Fill(NVTX, jetResponse);
       hist1D[matchPre + "_JetResponse_vs_NVTX_L1Ptge12_prof"]->Fill(NVTX, jetResponse);

       if ( L1Pt >= 30 ){
	 hist2D[matchPre + "_JetResponse_vs_NVTX_L1Ptge30"]     ->Fill(NVTX, jetResponse);
	 hist1D[matchPre + "_JetResponse_vs_NVTX_L1Ptge30_prof"]->Fill(NVTX, jetResponse);
       }
    }


//     hist2D[matchPre + "_DeltaEta_vs_OffPT" ]    ->Fill(ak5Pt, deltaEta);
//     hist1D[matchPre + "_DeltaEta_vs_OffPT_prof"]->Fill(ak5Pt, deltaEta);
//     hist2D[matchPre + "_DeltaPhi_vs_OffPT" ]    ->Fill(ak5Pt, deltaPhi);
//     hist1D[matchPre + "_DeltaPhi_vs_OffPT_prof"]->Fill(ak5Pt, deltaPhi);

//     // NVTX Correlations
//     hist2D[matchPre + "_L1PT_vs_NVTX"     ]    ->Fill(NVTX, L1Pt);
//     hist1D[matchPre + "_L1PT_vs_NVTX_prof"]    ->Fill(NVTX, L1Pt);
//     hist2D[matchPre + "_DeltaPTRel_vs_NVTX"  ]    ->Fill(NVTX, deltaPtRel);
//     hist1D[matchPre + "_DeltaPTRel_vs_NVTX_prof"] ->Fill(NVTX, deltaPtRel);
//     hist2D[matchPre + "_DeltaEta_vs_NVTX" ]    ->Fill(NVTX, deltaEta);
//     hist1D[matchPre + "_DeltaEta_vs_NVTX_prof"]->Fill(NVTX, deltaEta);
//     hist2D[matchPre + "_DeltaPhi_vs_NVTX" ]    ->Fill(NVTX, deltaPhi);
//     hist1D[matchPre + "_DeltaPhi_vs_NVTX_prof"]->Fill(NVTX, deltaPhi);

//       // Phi Correlations
//     hist2D[matchPre + "_L1PT_vs_Phi"     ]    ->Fill(L1Phi, L1Pt);    
//     hist1D[matchPre + "_L1PT_vs_Phi_prof"]    ->Fill(L1Phi, L1Pt);    
//     hist2D[matchPre + "_DeltaPTRel_vs_Phi"  ]    ->Fill(L1Phi, deltaPtRel);         
//     hist1D[matchPre + "_DeltaPTRel_vs_Phi_prof"] ->Fill(L1Phi, deltaPtRel);         
//     hist2D[matchPre + "_DeltaEta_vs_Phi" ]    ->Fill(L1Phi, deltaEta);        
//     hist1D[matchPre + "_DeltaEta_vs_Phi_prof"]->Fill(L1Phi, deltaEta);        
//     hist2D[matchPre + "_DeltaPhi_vs_Phi" ]    ->Fill(L1Phi, deltaPhi);        
//     hist1D[matchPre + "_DeltaPhi_vs_Phi_prof"]->Fill(L1Phi, deltaPhi);        
//       // Eta Correlations
//     hist2D[matchPre + "_L1PT_vs_Eta"     ]    ->Fill(L1Eta, L1Pt);          
//     hist1D[matchPre + "_L1PT_vs_Eta_prof"]    ->Fill(L1Eta, L1Pt);          
//     hist2D[matchPre + "_DeltaPTRel_vs_Eta"  ]    ->Fill(L1Eta, deltaPtRel);         
//     hist1D[matchPre + "_DeltaPTRel_vs_Eta_prof"] ->Fill(L1Eta, deltaPtRel);         
//     hist2D[matchPre + "_DeltaEta_vs_Eta" ]    ->Fill(L1Eta, deltaEta);        
//     hist1D[matchPre + "_DeltaEta_vs_Eta_prof"]->Fill(L1Eta, deltaEta);        
//     hist2D[matchPre + "_DeltaPhi_vs_Eta" ]    ->Fill(L1Eta, deltaPhi);        
//     hist1D[matchPre + "_DeltaPhi_vs_Eta_prof"]->Fill(L1Eta, deltaPhi);        


    // *****************************************************************
    // *                           Eta Binned                          *
    // *****************************************************************


    // ----------------------------------------------------------------------------------------------------
    // Eta stuff
    // ************************************************************
    // Find which in which Eta bin the L1 jet resides

    TString EtaStr;
    for (uint pEta = 1; pEta < etaRegionSlice.size(); ++pEta){

      // Get Eta bin lower and upper bounds
      double EtaLow  = etaRegionSlice[ pEta - 1];
      double EtaHigh = etaRegionSlice[ pEta ];
      
      if ( (L1Eta >= EtaLow) && (L1Eta < EtaHigh) ){
	TString EtaLowStr  = Form("%1.3f", EtaLow );
	TString EtaHighStr = Form("%1.3f", EtaHigh );
	EtaStr = "Eta_" + EtaLowStr + "_to_" + EtaHighStr; 
	break;
      }

    }



//     // Resolution distributions
//     hist1D[matchPre + "_DeltaEta_" + iEtaStr]->Fill(deltaEta);
//     hist1D[matchPre + "_DeltaPhi_" + iEtaStr]->Fill(deltaPhi);
//     hist1D[matchPre + "_DeltaPTRel_" + iEtaStr] ->Fill(deltaPtRel); 

//     hist1D[matchPre + "_Eta_" + iEtaStr]     ->Fill(L1Eta);
//     hist1D[matchPre + "_Phi_" + iEtaStr]     ->Fill(L1Phi);
//     hist1D[matchPre + "_PT_" + iEtaStr]      ->Fill(L1Pt);

//     hist1D[matchPre + "_Ak5Eta_" + iEtaStr]  ->Fill(ak5Eta);
//     hist1D[matchPre + "_Ak5Phi_" + iEtaStr]  ->Fill(ak5Phi);
//     hist1D[matchPre + "_Ak5PT_" + iEtaStr]   ->Fill(ak5Pt);

    // Online-offline quantity correlations
//     hist2D[matchPre + "_L1Eta_vs_OffEta_" + iEtaStr]               ->Fill(ak5Eta, L1Eta);       
//     hist1D[matchPre + "_L1Eta_vs_OffEta_" + iEtaStr + "_prof"]     ->Fill(ak5Eta, L1Eta);       
//     hist2D[matchPre + "_L1Phi_vs_OffPhi_" + iEtaStr]               ->Fill(ak5Phi, L1Phi);    
//     hist1D[matchPre + "_L1Phi_vs_OffPhi_" + iEtaStr + "_prof"]     ->Fill(ak5Phi, L1Phi);    
//     hist2D[matchPre + "_DeltaPhi_vs_DeltaEta_" + iEtaStr]          ->Fill(deltaEta, deltaPhi);    
//     hist1D[matchPre + "_DeltaPhi_vs_DeltaEta_" + iEtaStr + "_prof"]->Fill(deltaEta, deltaPhi);    
//     hist2D[matchPre + "_DeltaPTRel_vs_L1PT_" + iEtaStr]               ->Fill(L1Pt, deltaPtRel);    
//     hist1D[matchPre + "_DeltaPTRel_vs_L1PT_" + iEtaStr + "_prof"]     ->Fill(L1Pt, deltaPtRel);    

    // pT Correlations
    hist2D[matchPre + "_L1PT_vs_OffPT_"     + EtaStr]             ->Fill(ak5Pt, L1Pt);
    hist1D[matchPre + "_L1PT_vs_OffPT_"     + EtaStr + "_prof"]   ->Fill(ak5Pt, L1Pt);
    hist2D[matchPre + "_OffPT_vs_L1PT_"     + EtaStr]             ->Fill(L1Pt, ak5Pt);
    hist1D[matchPre + "_OffPT_vs_L1PT_"     + EtaStr + "_prof"]   ->Fill(L1Pt, ak5Pt);
    hist2D[matchPre + "_DeltaPTRel_vs_OffPT_"  + EtaStr]          ->Fill(ak5Pt, deltaPtRel);
    hist1D[matchPre + "_DeltaPTRel_vs_OffPT_"  + EtaStr + "_prof"]->Fill(ak5Pt, deltaPtRel);
    hist2D[matchPre + "_DeltaPTRel_vs_DeltaPT_"  + EtaStr]          ->Fill(deltaPt, deltaPtRel);
    hist1D[matchPre + "_DeltaPTRel_vs_DeltaPT_"  + EtaStr + "_prof"]->Fill(deltaPt, deltaPtRel);

    hist2D[matchPre + "_DeltaPT_vs_OffPT_"  + EtaStr]             ->Fill(ak5Pt, deltaPt);
    hist1D[matchPre + "_DeltaPT_vs_OffPT_"  + EtaStr + "_prof"]   ->Fill(ak5Pt, deltaPt);
    hist1D[matchPre + "_JetResponse_" + EtaStr ]                  ->Fill(jetResponse);
    hist2D[matchPre + "_JetResponse_vs_L1Eta"]                    ->Fill(L1Eta, jetResponse);
    hist1D[matchPre + "_JetResponse_vs_L1Eta_prof"]               ->Fill(L1Eta, jetResponse);

    hist2D[matchPre + "_JetResponse_vs_OffPT_" + EtaStr]          ->Fill(ak5Pt, jetResponse);
    hist1D[matchPre + "_JetResponse_vs_OffPT_" + EtaStr + "_prof"]->Fill(ak5Pt, jetResponse);


    hist2D[matchPre + "_DeltaPT_vs_L1PT_"  + EtaStr]             ->Fill(L1Pt, deltaPt);
    hist1D[matchPre + "_DeltaPT_vs_L1PT_"  + EtaStr + "_prof"]   ->Fill(L1Pt, deltaPt);

    hist2D[matchPre + "_JetResponse_vs_L1PT_" + EtaStr]          ->Fill(L1Pt, jetResponse);
    hist1D[matchPre + "_JetResponse_vs_L1PT_" + EtaStr + "_prof"]->Fill(L1Pt, jetResponse);

    hist2D[matchPre + "_JetResponse_vs_NVTX_" + EtaStr]          ->Fill(NVTX, jetResponse);
    hist1D[matchPre + "_JetResponse_vs_NVTX_" + EtaStr + "_prof"]->Fill(NVTX, jetResponse);

    if ( L1Pt >= 12 ){
       hist2D[matchPre + "_JetResponse_vs_NVTX_L1Ptge12_" + EtaStr]          ->Fill(NVTX, jetResponse);
       hist1D[matchPre + "_JetResponse_vs_NVTX_L1Ptge12_" + EtaStr + "_prof"]->Fill(NVTX, jetResponse);

       if ( L1Pt >= 30 ){
         hist2D[matchPre + "_JetResponse_vs_NVTX_L1Ptge30_" + EtaStr]          ->Fill(NVTX, jetResponse);
         hist1D[matchPre + "_JetResponse_vs_NVTX_L1Ptge30_" + EtaStr + "_prof"]->Fill(NVTX, jetResponse);
       }
    }
    hist2D[matchPre + "_DeltaPT_vs_NVTX_"  + EtaStr]             ->Fill(NVTX, deltaPt);
    hist1D[matchPre + "_DeltaPT_vs_NVTX_"  + EtaStr + "_prof"]   ->Fill(NVTX, deltaPt);

//     hist2D[matchPre + "_JetResponse_vs_Phi_" + EtaStr]           ->Fill(L1Phi, jetResponse);
//     hist1D[matchPre + "_JetResponse_vs_Phi_" + EtaStr + "_prof"] ->Fill(L1Phi, jetResponse);
//     hist2D[matchPre + "_DeltaPT_vs_Phi_"  + EtaStr]              ->Fill(L1Phi, deltaPt);
//     hist1D[matchPre + "_DeltaPT_vs_Phi_"  + EtaStr + "_prof"]    ->Fill(L1Phi, deltaPt);



//     hist2D[matchPre + "_DeltaEta_vs_OffPT_" + iEtaStr]          ->Fill(ak5Pt, deltaEta);
//     hist1D[matchPre + "_DeltaEta_vs_OffPT_" + iEtaStr + "_prof"]->Fill(ak5Pt, deltaEta);
//     hist2D[matchPre + "_DeltaPhi_vs_OffPT_" + iEtaStr]          ->Fill(ak5Pt, deltaPhi);
//     hist1D[matchPre + "_DeltaPhi_vs_OffPT_" + iEtaStr + "_prof"]->Fill(ak5Pt, deltaPhi);

//     // NVTX Correlations
//     hist2D[matchPre + "_L1PT_vs_NVTX_"     + iEtaStr]          ->Fill(NVTX, L1Pt);
//     hist1D[matchPre + "_L1PT_vs_NVTX_"     + iEtaStr + "_prof"]->Fill(NVTX, L1Pt);
//     hist2D[matchPre + "_DeltaPTRel_vs_NVTX_"  + iEtaStr]          ->Fill(NVTX, deltaPtRel);
//     hist1D[matchPre + "_DeltaPTRel_vs_NVTX_"  + iEtaStr + "_prof"]->Fill(NVTX, deltaPtRel);
//     hist2D[matchPre + "_DeltaEta_vs_NVTX_" + iEtaStr]          ->Fill(NVTX, deltaEta);
//     hist1D[matchPre + "_DeltaEta_vs_NVTX_" + iEtaStr + "_prof"]->Fill(NVTX, deltaEta);
//     hist2D[matchPre + "_DeltaPhi_vs_NVTX_" + iEtaStr]          ->Fill(NVTX, deltaPhi);
//     hist1D[matchPre + "_DeltaPhi_vs_NVTX_" + iEtaStr + "_prof"]->Fill(NVTX, deltaPhi);

//       // Phi Correlations
//     hist2D[matchPre + "_L1PT_vs_Phi_"     + iEtaStr]          ->Fill(L1Phi, L1Pt);    
//     hist1D[matchPre + "_L1PT_vs_Phi_"     + iEtaStr + "_prof"]->Fill(L1Phi, L1Pt);    
//     hist2D[matchPre + "_DeltaPTRel_vs_Phi_"  + iEtaStr]          ->Fill(L1Phi, deltaPtRel);         
//     hist1D[matchPre + "_DeltaPTRel_vs_Phi_"  + iEtaStr + "_prof"]->Fill(L1Phi, deltaPtRel);         
//     hist2D[matchPre + "_DeltaEta_vs_Phi_" + iEtaStr]          ->Fill(L1Phi, deltaEta);        
//     hist1D[matchPre + "_DeltaEta_vs_Phi_" + iEtaStr + "_prof"]->Fill(L1Phi, deltaEta);        
//     hist2D[matchPre + "_DeltaPhi_vs_Phi_" + iEtaStr]          ->Fill(L1Phi, deltaPhi);        
//     hist1D[matchPre + "_DeltaPhi_vs_Phi_" + iEtaStr + "_prof"]->Fill(L1Phi, deltaPhi);        
//       // Eta Correlations
//     hist2D[matchPre + "_L1PT_vs_Eta_"     + iEtaStr]          ->Fill(L1Eta, L1Pt);          
//     hist1D[matchPre + "_L1PT_vs_Eta_"     + iEtaStr + "_prof"]->Fill(L1Eta, L1Pt);          
//     hist2D[matchPre + "_DeltaPTRel_vs_Eta_"  + iEtaStr]          ->Fill(L1Eta, deltaPtRel);         
//     hist1D[matchPre + "_DeltaPTRel_vs_Eta_"  + iEtaStr + "_prof"]->Fill(L1Eta, deltaPtRel);         
//     hist2D[matchPre + "_DeltaEta_vs_Eta_" + iEtaStr]          ->Fill(L1Eta, deltaEta);        
//     hist1D[matchPre + "_DeltaEta_vs_Eta_" + iEtaStr + "_prof"]->Fill(L1Eta, deltaEta);        
//     hist2D[matchPre + "_DeltaPhi_vs_Eta_" + iEtaStr]          ->Fill(L1Eta, deltaPhi);        
//     hist1D[matchPre + "_DeltaPhi_vs_Eta_" + iEtaStr + "_prof"]->Fill(L1Eta, deltaPhi);        
    




    // *****************************************************************
    // *                          NVTX Binned                          *
    // *****************************************************************


//     // pT Correlations
//     hist2D[matchPre + "_L1PT_vs_OffPT_"     + NVTXStr]          ->Fill(ak5Pt, L1Pt);
//     hist1D[matchPre + "_L1PT_vs_OffPT_"     + NVTXStr + "_prof"]->Fill(ak5Pt, L1Pt);
//     hist2D[matchPre + "_OffPT_vs_L1PT_"     + NVTXStr]          ->Fill(L1Pt, ak5Pt);
//     hist1D[matchPre + "_OffPT_vs_L1PT_"     + NVTXStr + "_prof"]->Fill(L1Pt, ak5Pt);








    // Calibration Histograms
    // --------------------------------------------------------------------------------

//     hist1D[matchPre + "_DeltaPTRel_vs_GlobalRho_prof_" + iEtaStr]  ->Fill(L1Pt, ak5Pt);
//     hist1D[matchPre + "_DeltaPTRel_vs_EorHFired_prof_" + iEtaStr]  ->Fill(L1Pt, ak5Pt);
//     hist1D[matchPre + "_DeltaPTRel_vs_EplusHtotal_prof_" + iEtaStr]->Fill(L1Pt, ak5Pt);




  } // End matched jets




}






//
// Produce some benchmark plots for the online jet by matching to offline jets
//
void 
JetHist::benchmarkJets( edm::Handle<l1extra::L1JetParticleCollection> const& L1Jets, edm::Handle<l1extra::L1JetParticleCollection> const& ak5Jets,
			TString prefix, double ak5PTthreshold, double L1PTthreshold, int maxMatchedJets, double ak5EtaRange, double ak5MaxPT ){


  PRINT("Benchmark Jets");
  // Find which histogram prefix to use
  int preIndex = -1;
  for (unsigned int iMatch = 0; iMatch < matchPrefix.size(); iMatch++){
    if (prefix == matchPrefix[iMatch]){
      preIndex = iMatch;
      break;
    }
  }
  if (preIndex == -1){
    std::cout << "ERROR - Prefix '" << prefix << "' was not found.\n";
    return;
  } 
  TString matchPre = matchPrefix[preIndex];
  TString matchLab = matchLabel[preIndex];

  SUBPRINT("Extracted prefix")


  //      ********************************************************************************
  //      *                      PrePUS Online-Offline jet matching                      *
  //      ********************************************************************************
  //


  std::vector<TLorentzVector> ak5Vec, L1Vec;

  // Hardonic energy sums
  double L1HT(0), ak5HT(0), L1MHT(0), ak5MHT(0);
  TLorentzVector L1MHTVector; 
  TLorentzVector ak5MHTVector;

  //           *********************************************
  //           *          Fill ak5 TLorentzVector          *
  //           *********************************************

  for ( l1extra::L1JetParticleCollection::const_iterator ak5_It = ak5Jets->begin(); ak5_It != ak5Jets->end(); ++ak5_It ){

    double ak5Pt         = ak5_It->p4().Pt();
    double ak5Eta        = ak5_It->p4().Eta();
    double ak5Phi        = ak5_It->p4().Phi();
    double ak5M          = ak5_It->p4().M();

    TLorentzVector ak5Jet;
    ak5Jet.SetPtEtaPhiM( ak5Pt, ak5Eta, ak5Phi, ak5M );

    // Calculate energy sums
    ak5HT        += ak5Pt;
    ak5MHTVector += ak5Jet;

    // ========================================
    // Jet quality cuts
    // ========================================

    if( fabs( ak5Eta ) > ak5EtaRange ) continue;
    if( ak5Pt < ak5PTthreshold )       continue;

    // ========================================

    ak5Vec.push_back( ak5Jet );
  }
  
  // Calculate scalar projection of the MHT vector
  ak5MHT = ak5MHTVector.Pt();

  // <TEMPORARY> - Calibrated ak5 jets were previously not sorted, fixed only in new ntuples. Remove when transition to new ntuples is complete.
     // Sort calibrated jets
  //     std::sort( ak5Vec.begin(), ak5Vec.end(), TLorentzVectorRankDescending );
  // </TEMPORARY>


  //           *********************************************
  //           *          Fill L1 TLorentzVector           *
  //           *********************************************

  for ( l1extra::L1JetParticleCollection::const_iterator L1_It = L1Jets->begin(); L1_It != L1Jets->end(); ++L1_It ){

    double L1Pt         = L1_It->p4().Pt();
    double L1Eta        = L1_It->p4().Eta();
    double L1Phi        = L1_It->p4().Phi();
    double L1M          = L1_It->p4().M();

    TLorentzVector L1Jet;
    L1Jet.SetPtEtaPhiM( L1Pt, L1Eta, L1Phi, L1M );

    // Calculate energy sums
    L1HT        += L1Pt;
    L1MHTVector += L1Jet;

    // Determine whether to analyse jet
    if( L1Pt < L1PTthreshold) continue;
    if( fabs( L1Eta ) > ak5EtaRange ) continue;
    L1Vec.push_back( L1Jet );

  }

  // Calculate scalar projection of the MHT vector
  L1MHT = L1MHTVector.Pt(); // THIS IS WRONG, USE THE GCT HT CALCULATION


  // Energy sums
  // ------------------------------
    hist1D[matchPre + "_HT"]      ->Fill(L1HT);
    hist1D[matchPre + "_MHT"]     ->Fill(L1MHT);
    hist1D[matchPre + "_Ak5HT"]   ->Fill(ak5HT);
    hist1D[matchPre + "_Ak5MHT"]  ->Fill(ak5MHT);
    hist2D[ matchPre + "_L1MHT_vs_L1HT"]       ->Fill(L1HT,  L1MHT);
    hist1D[ matchPre + "_L1MHT_vs_L1HT_prof"]  ->Fill(L1HT,  L1MHT);
    hist2D[ matchPre + "_Ak5MHT_vs_Ak5HT"]     ->Fill(ak5HT, ak5MHT);
    hist1D[ matchPre + "_Ak5MHT_vs_Ak5HT_prof"]->Fill(ak5HT, ak5MHT);
    hist2D[ matchPre + "_Ak5HT_vs_L1HT"]       ->Fill(L1HT,  ak5HT);
    hist1D[ matchPre + "_Ak5HT_vs_L1HT_prof"]  ->Fill(L1HT,  ak5HT);
    hist2D[ matchPre + "_Ak5MHT_vs_L1MHT"]     ->Fill(L1MHT, ak5MHT);
    hist1D[ matchPre + "_Ak5MHT_vs_L1MHT_prof"]->Fill(L1MHT, ak5MHT);


    // HT turnons 
    // --------------------------------------------------
    bool passHT75(false), passHT100(false), passHT150(false), passHT175(false);
    if ( L1HT >= 75 ){
      passHT75 = true;
      if ( L1HT >= 100 ){
	passHT100 = true;
	if ( L1HT >= 150 ){
	  passHT150 = true;
	  if ( L1HT >= 175 ){
	    passHT175 = true;
	  }
	}
      }
    }
    histEff[ matchPre + "_TurnOn_HTgt75GeV" ] ->Fill( passHT75,  ak5HT );
    histEff[ matchPre + "_TurnOn_HTgt100GeV" ]->Fill( passHT100, ak5HT );
    histEff[ matchPre + "_TurnOn_HTgt150GeV" ]->Fill( passHT150, ak5HT );
    histEff[ matchPre + "_TurnOn_HTgt175GeV" ]->Fill( passHT175, ak5HT );

                                               
    // MHT turnons 
    // --------------------------------------------------
    bool passMHT50(false), passMHT75(false), passMHT100(false), passMHT150(false);
    if ( L1MHT >= 50 ){
      passMHT50 = true;
      if ( L1MHT >= 75 ){
	passMHT75 = true;
	if ( L1MHT >= 100 ){
	  passMHT100 = true;
	  if ( L1MHT >= 150 ){
	    passMHT150 = true;
	  }
	}
      }
    }
    histEff[ matchPre + "_TurnOn_MHTgt50GeV" ] ->Fill( passMHT50,  ak5MHT );
    histEff[ matchPre + "_TurnOn_MHTgt75GeV" ] ->Fill( passMHT75,  ak5MHT );
    histEff[ matchPre + "_TurnOn_MHTgt100GeV" ]->Fill( passMHT100, ak5MHT );
    histEff[ matchPre + "_TurnOn_MHTgt150GeV" ]->Fill( passMHT150, ak5MHT );










  // ------------------------------ Perform jet matching ------------------------------



  SUBPRINT("Extracted jets")

    // THIS IS REALLY INEFFICIENT. WHEN YOU HAVE TIME WRITE A CLASS THAT RETURNS std::vector< std::pair< TLorentz*, TLorentz*> >
    // WITH SORT ALGORITHMS AS AN OPTION


    // Perform deltaR matching, for calibration ONLY

  std::vector<pair_info> pairs = make_pairs( ak5Vec, L1Vec );
  std::sort(pairs.begin(), pairs.end(), sortDR);
  std::vector<int> L1MatchedIndex = analyse_pairs(pairs, L1Vec.size(), maxDeltaR);
  

  int L1JetsMatched    = 0;
  int L1JetsNotMatched = 0;
  
  // ****************************************************************************************************
  // *                                      Loop through L1 jets (Matched)                              *
  // ****************************************************************************************************
  for(unsigned int iRank = 0; iRank < L1MatchedIndex.size(); ++iRank) {

    int L1Index  = iRank;
    // Index of the matched ak5 jet
    int ak5Index = L1MatchedIndex[ L1Index ];

    double L1Pt     = L1Vec[ L1Index ].Pt();
    double L1Eta    = L1Vec[ L1Index ].Eta();
    double L1Phi    = L1Vec[ L1Index ].Phi();


    L1JetsMatched++;
    // Restrict to specified number of jets
    if ( L1JetsMatched > maxMatchedJets ){
      break;
    }

    
    if ( ak5Index == -1 ){
      
      // Fill mismatched ditributions
      if ( L1Index < int(jetIndexLabel.size()) ){
	
	// Get the current jet label, account for jet labeling starting at 0. LeadJet = 0, ...
	int jetIndex     = L1Index;
	TString jetLabel = "_" + jetIndexLabel[ jetIndex ];
	
	hist1D[matchPre + jetLabel + "_NoMatch_PT" ]  ->Fill( L1Pt  );	
	hist1D[matchPre + jetLabel + "_NoMatch_Eta" ] ->Fill( L1Eta );
	hist1D[matchPre + jetLabel + "_NoMatch_Phi" ] ->Fill( L1Phi );


	// Pick the corresponding ranked offline jet, if it exists
	if ( iRank < ak5Vec.size() ){

	  double ak5Pt    = ak5Vec[ iRank ].Pt();
	  double ak5Eta   = ak5Vec[ iRank ].Eta();
	  //    double ak5Phi   = ak5Vec[ ak5Index ].Phi();
	  double deltaPt     = (L1Pt - ak5Pt);
	  double deltaPtRel  = deltaPt/ak5Pt; 
	  //    double jetResponse = L1Pt/ak5Pt;
	  double deltaEta = ak5Eta - L1Eta;
	  double deltaPhi = ak5Vec[ iRank ].DeltaPhi( L1Vec[ L1Index ] );
	  
	  hist1D[matchPre + jetLabel + "_NoMatch_deltaPTRel" ]->Fill( deltaPtRel );	
	  hist1D[matchPre + jetLabel + "_NoMatch_deltaEta" ]  ->Fill( deltaEta   );
	  hist1D[matchPre + jetLabel + "_NoMatch_deltaPhi" ]  ->Fill( deltaPhi   );
	  
	  double deltaR = sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta );
	  hist1D[matchPre + jetLabel + "_NoMatch_deltaR" ]               ->Fill( deltaR );
	  hist2D[matchPre + jetLabel + "_NoMatch_deltaPTRel_vs_deltaR" ] ->Fill( deltaR, deltaPtRel );

	}
      }


      L1JetsNotMatched++;
      continue; // Jets cannot be matched, continue
    }

    double ak5Pt    = ak5Vec[ ak5Index ].Pt();
    double ak5Eta   = ak5Vec[ ak5Index ].Eta();
    double ak5Phi   = ak5Vec[ ak5Index ].Phi();
    double deltaPt     = (L1Pt - ak5Pt);
    double deltaPtRel  = deltaPt/ak5Pt; 
    double jetResponse = L1Pt/ak5Pt;
    double deltaEta = ak5Eta - L1Eta;
    double deltaPhi = ak5Vec[ ak5Index ].DeltaPhi( L1Vec[ L1Index ] );
    double deltaR = sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta );
//     //    std::cout << "deltaPT = " << deltaPtRel << "\tdeltaEta = " << deltaEta << "\tdeltaPhi = " << deltaPhi << "\n";


//     if ( ( prefix == "RecalibPrePUS_ak5PUS" ) && (i == 0) && (deltaPtRel > 0.3) ){
//       std::cout << "\tEvent: " << Eventnr << "\tdeltaPt = " << deltaPtRel << "\tpT = " << L1Pt << "\tEta = " << L1Eta  << "\n";
//      }
    


    // ****************************************************************************************************
    // *                                             Unbinned                                             *
    // ****************************************************************************************************


    // ************************************************************
    // *              Distributions for specific jets             *
    // ************************************************************

    if ( L1Index < int(jetIndexLabel.size()) ){

      // Get the current jet label, account for jet labeling starting at 0. LeadJet = 0, ...
      int jetIndex     = L1Index;
      TString jetLabel = "_" + jetIndexLabel[ jetIndex ];

      hist1D[matchPre + jetLabel + "_DeltaEta" ]  ->Fill( deltaEta );
      hist1D[matchPre + jetLabel + "_DeltaPhi" ]  ->Fill( deltaPhi );
      hist1D[matchPre + jetLabel + "_DeltaPTRel" ]->Fill( deltaPtRel );

      hist2D[matchPre + jetLabel + "_DeltaPTRel_vs_OffPT" ]     ->Fill(ak5Pt, deltaPtRel);    
      hist1D[matchPre + jetLabel + "_DeltaPTRel_vs_OffPT_prof"] ->Fill(ak5Pt, deltaPtRel);    
      hist2D[matchPre + jetLabel + "_DeltaPTRel_vs_L1PT" ]      ->Fill(L1Pt,  deltaPtRel);    
      hist1D[matchPre + jetLabel + "_DeltaPTRel_vs_L1PT_prof"]  ->Fill(L1Pt,  deltaPtRel);    

      hist2D[matchPre + jetLabel + "_DeltaPTRel_vs_NVTX"  ]     ->Fill(NVTX,  deltaPtRel);
      hist1D[matchPre + jetLabel + "_DeltaPTRel_vs_NVTX_prof"]  ->Fill(NVTX,  deltaPtRel);
      hist2D[matchPre + jetLabel + "_DeltaPTRel_vs_Eta"  ]      ->Fill(L1Eta, deltaPtRel);         
      hist1D[matchPre + jetLabel + "_DeltaPTRel_vs_Eta_prof"]   ->Fill(L1Eta, deltaPtRel);         

      hist2D[matchPre + jetLabel + "_JetResponse_vs_OffPT" ]    ->Fill(ak5Pt, jetResponse);    
      hist1D[matchPre + jetLabel + "_JetResponse_vs_OffPT_prof"]->Fill(ak5Pt, jetResponse);    
      hist2D[matchPre + jetLabel + "_JetResponse_vs_NVTX"  ]    ->Fill(NVTX,  jetResponse);
      hist1D[matchPre + jetLabel + "_JetResponse_vs_NVTX_prof"] ->Fill(NVTX,  jetResponse);
      hist2D[matchPre + jetLabel + "_JetResponse_vs_Eta"  ]     ->Fill(L1Eta, jetResponse);         
      hist1D[matchPre + jetLabel + "_JetResponse_vs_Eta_prof"]  ->Fill(L1Eta, jetResponse);         

      hist2D[matchPre + jetLabel + "_deltaPTRel_vs_deltaR" ]    ->Fill( deltaR, deltaPtRel );

     }

    // ************************************************************
    // *                         All jets                         *
    // ************************************************************
    hist1D[matchPre + "_DeltaEta"]->Fill(deltaEta);
    hist1D[matchPre + "_DeltaPhi"]->Fill(deltaPhi);
    hist1D[matchPre + "_DeltaPTRel"] ->Fill(deltaPtRel); 

    hist1D[matchPre + "_Eta"]     ->Fill(L1Eta);
    hist1D[matchPre + "_Phi"]     ->Fill(L1Phi);
    hist1D[matchPre + "_PT"]      ->Fill(L1Pt);

    hist1D[matchPre + "_Ak5Eta"]  ->Fill(ak5Eta);
    hist1D[matchPre + "_Ak5Phi"]  ->Fill(ak5Phi);
    hist1D[matchPre + "_Ak5PT"]   ->Fill(ak5Pt);

    // Online-offline quantity correlations
//     hist2D[matchPre + "_L1Eta_vs_OffEta" ]         ->Fill(ak5Eta, L1Eta);       
//     hist1D[matchPre + "_L1Eta_vs_OffEta_prof"]     ->Fill(ak5Eta, L1Eta);       
//     hist2D[matchPre + "_L1Phi_vs_OffPhi" ]         ->Fill(ak5Phi, L1Phi);    
//     hist1D[matchPre + "_L1Phi_vs_OffPhi_prof"]     ->Fill(ak5Phi, L1Phi);    
//     hist2D[matchPre + "_DeltaPhi_vs_DeltaEta" ]    ->Fill(deltaEta, deltaPhi);    
//     hist1D[matchPre + "_DeltaPhi_vs_DeltaEta_prof"]->Fill(deltaEta, deltaPhi);    

    // pT Correlations
    hist2D[matchPre + "_L1PT_vs_OffPT"     ]    ->Fill(ak5Pt, L1Pt);
    hist1D[matchPre + "_L1PT_vs_OffPT_prof"]    ->Fill(ak5Pt, L1Pt);
    hist2D[matchPre + "_OffPT_vs_L1PT"     ]    ->Fill(L1Pt, ak5Pt);
    hist1D[matchPre + "_OffPT_vs_L1PT_prof"]    ->Fill(L1Pt, ak5Pt);
    hist2D[matchPre + "_DeltaPTRel_vs_OffPT"  ]    ->Fill(ak5Pt, deltaPtRel);
    hist1D[matchPre + "_DeltaPTRel_vs_OffPT_prof"] ->Fill(ak5Pt, deltaPtRel);
    hist2D[matchPre + "_DeltaPTRel_vs_L1PT"  ]    ->Fill(L1Pt, deltaPtRel);
    hist1D[matchPre + "_DeltaPTRel_vs_L1PT_prof"] ->Fill(L1Pt, deltaPtRel);
    hist2D[matchPre + "_DeltaEta_vs_OffPT" ]    ->Fill(ak5Pt, deltaEta);
    hist1D[matchPre + "_DeltaEta_vs_OffPT_prof"]->Fill(ak5Pt, deltaEta);
    hist2D[matchPre + "_DeltaPhi_vs_OffPT" ]    ->Fill(ak5Pt, deltaPhi);
    hist1D[matchPre + "_DeltaPhi_vs_OffPT_prof"]->Fill(ak5Pt, deltaPhi);

    // NVTX Correlations
    hist2D[matchPre + "_L1PT_vs_NVTX"     ]    ->Fill(NVTX, L1Pt);
    hist1D[matchPre + "_L1PT_vs_NVTX_prof"]    ->Fill(NVTX, L1Pt);
    hist2D[matchPre + "_DeltaPTRel_vs_NVTX"  ]    ->Fill(NVTX, deltaPtRel);
    hist1D[matchPre + "_DeltaPTRel_vs_NVTX_prof"] ->Fill(NVTX, deltaPtRel);
    hist2D[matchPre + "_DeltaEta_vs_NVTX" ]    ->Fill(NVTX, deltaEta);
    hist1D[matchPre + "_DeltaEta_vs_NVTX_prof"]->Fill(NVTX, deltaEta);
    hist2D[matchPre + "_DeltaPhi_vs_NVTX" ]    ->Fill(NVTX, deltaPhi);
    hist1D[matchPre + "_DeltaPhi_vs_NVTX_prof"]->Fill(NVTX, deltaPhi);
    
    // Phi Correlations
    hist2D[matchPre + "_L1PT_vs_Phi"     ]    ->Fill(L1Phi, L1Pt);    
    hist1D[matchPre + "_L1PT_vs_Phi_prof"]    ->Fill(L1Phi, L1Pt);    
    hist2D[matchPre + "_DeltaPTRel_vs_Phi"  ]    ->Fill(L1Phi, deltaPtRel);         
    hist1D[matchPre + "_DeltaPTRel_vs_Phi_prof"] ->Fill(L1Phi, deltaPtRel);         
    hist2D[matchPre + "_DeltaEta_vs_Phi" ]    ->Fill(L1Phi, deltaEta);        
    hist1D[matchPre + "_DeltaEta_vs_Phi_prof"]->Fill(L1Phi, deltaEta);        
    hist2D[matchPre + "_DeltaPhi_vs_Phi" ]    ->Fill(L1Phi, deltaPhi);        
    hist1D[matchPre + "_DeltaPhi_vs_Phi_prof"]->Fill(L1Phi, deltaPhi);        
    // Eta Correlations
    hist2D[matchPre + "_L1PT_vs_Eta"     ]    ->Fill(L1Eta, L1Pt);          
    hist1D[matchPre + "_L1PT_vs_Eta_prof"]    ->Fill(L1Eta, L1Pt);          
    hist2D[matchPre + "_OffPT_vs_Eta"     ]    ->Fill(ak5Eta, ak5Pt);          
    hist1D[matchPre + "_OffPT_vs_Eta_prof"]    ->Fill(ak5Eta, ak5Pt);          
    hist2D[matchPre + "_DeltaPTRel_vs_Eta"  ]    ->Fill(L1Eta, deltaPtRel);         
    hist1D[matchPre + "_DeltaPTRel_vs_Eta_prof"] ->Fill(L1Eta, deltaPtRel);         
    hist2D[matchPre + "_DeltaEta_vs_Eta" ]    ->Fill(L1Eta, deltaEta);        
    hist1D[matchPre + "_DeltaEta_vs_Eta_prof"]->Fill(L1Eta, deltaEta);        
    hist2D[matchPre + "_DeltaPhi_vs_Eta" ]    ->Fill(L1Eta, deltaPhi);        
    hist1D[matchPre + "_DeltaPhi_vs_Eta_prof"]->Fill(L1Eta, deltaPhi);        


   } // End L1 matched loop





  turnOnL1Loop( L1Vec, ak5Vec, matchPre );
  turnOnAk5Loop( L1Vec, ak5Vec, matchPre );
   
  // APPEARS TO CAUSE A MEMORY ISSUE
  //        turnOnAk5MatchedLoop( L1Vec, ak5Vec, matchPre );



  //  histEff[ matchPre + "_matchingEff1Jet" ]->Fill( matchedJet, ak5Pt );
//       histEff[ matchPre + "_matchingEff2Jet" ]->Fill( matchedJet, ak5Pt );
//       histEff[ matchPre + "_matchingEff3Jet" ]->Fill( matchedJet, ak5Pt );
//       histEff[ matchPre + "_matchingEff4Jet" ]->Fill( matchedJet, ak5Pt );




  

}













// ********************************************************************************************************************************************
// *                                                                 Turn ons                                                                 *
// ********************************************************************************************************************************************

void 
JetHist::turnOnL1Loop( std::vector<TLorentzVector> const &L1Vec, std::vector<TLorentzVector> const &ak5Vec, TString matchPre ){



  // Produce turn ons requiring no matching between L1 and RECO jets


//   // ********************************************************************************
//   // Lead-second jet turn ons
//   // ********************************************************************************
//   //  if (( L1Vec.size() > 1) && ( ak5Vec.size() > 1)){
      
//  //    double leadL1Pt    =  L1Vec[ 0 ].Pt();
// //     double secondL1Pt  =  L1Vec[ 1 ].Pt();
// //     double leadL1Eta   =  L1Vec[ 0 ].Eta();
// //     double secondL1Eta =  L1Vec[ 1 ].Eta();
// //     double leadL1Phi   =  L1Vec[ 0 ].Phi();
// //     double secondL1Phi =  L1Vec[ 1 ].Phi();

// //     double leadAk5Pt    =  ak5Vec[ 0 ].Pt();
// //     double secondAk5Pt  =  ak5Vec[ 1 ].Pt();
// //     double leadAk5Eta   =  ak5Vec[ 0 ].Eta();
// //     double secondAk5Eta =  ak5Vec[ 1 ].Eta();
// //     double leadAk5Phi   =  ak5Vec[ 0 ].Phi();
// //     double secondAk5Phi =  ak5Vec[ 1 ].Phi();

// //     // Position quantities
// //     double deltaL1Eta  = fabs(leadL1Eta - secondL1Eta);
// //     double deltaL1Phi  = fabs(L1Vec[ 0 ].DeltaPhi(L1Vec[ 1 ]));
// //     double deltaAk5Eta = fabs(leadAk5Eta - secondAk5Eta);
// //     double deltaAk5Phi = fabs(ak5Vec[ 0 ].DeltaPhi(ak5Vec[ 1 ]));
    
// //     // Invariant mass
// //     double invMassL1   = 2*leadL1Pt * secondL1Pt*( cosh( deltaL1Eta)  - cos(deltaL1Phi) );
// //     double invMassAk5  = 2*leadAk5Pt*secondAk5Pt*( cosh( deltaAk5Eta) - cos(deltaAk5Phi) );
// //   }





    // Loop through both jet collections, starting with the leading jets of both
  for(unsigned int iRank = 0; iRank < L1Vec.size(); ++iRank) {

    int L1Index  = iRank;
    int ak5Index = -1; // Initially set to no jet
    // Check the corresponding 'ith' jet exists offline
    if (iRank < ak5Vec.size()){
      ak5Index = iRank; // ak5 jet exists
    }


    if ( L1Index < int(jetIndexLabel.size()) ){
      
      // Get the current jet label, account for jet labeling starting at 0. LeadJet = 0, ... 
      int jetIndex = L1Index;
      TString jetLabel = "_" + jetIndexLabel[ jetIndex ];
      
      double L1Pt     = L1Vec[ L1Index ].Pt();
      hist1D[ matchPre + "_NoMatch" + jetLabel + "_L1PT" ]->Fill( L1Pt );

      // Fill pT threshold below which event fired 
      hist1D[ matchPre + "_NoMatch" + jetLabel + "_TriggeredProto" ]->Fill( L1Pt );


      double ak5MatchedPt = -1;
      if ( ak5Index != -1 ){
	ak5MatchedPt = ak5Vec[ ak5Index ].Pt();
      }

	// New rate study 
	bool passMatchedAk5ge30(false), passMatchedAk5ge40(false), passMatchedAk5ge50(false);
	if ( ak5MatchedPt >= 30){
	  passMatchedAk5ge30 = true;

	  if ( ak5MatchedPt >= 40){
	    passMatchedAk5ge40 = true;

	    if ( ak5MatchedPt >= 50){
	      passMatchedAk5ge50 = true;

	    }
	  }
	}
	
	histEff[ matchPre + "_Matched_ak5ge30_Eff_vs_L1Thresh" + jetLabel ]->Fill( passMatchedAk5ge30, L1Pt );
	histEff[ matchPre + "_Matched_ak5ge40_Eff_vs_L1Thresh" + jetLabel ]->Fill( passMatchedAk5ge40, L1Pt );
	histEff[ matchPre + "_Matched_ak5ge50_Eff_vs_L1Thresh" + jetLabel ]->Fill( passMatchedAk5ge50, L1Pt );
  
    }
	
    // Require an 'ith' offline jet
    if (ak5Index != -1){
    
      double ak5Pt    = ak5Vec[ ak5Index ].Pt();
      //      double ak5Eta   = ak5Vec[ ak5Index ].Eta();
      //      double ak5Phi   = ak5Vec[ ak5Index ].Phi();

      double L1Pt     = L1Vec[ L1Index ].Pt();
      //      double L1Eta    = L1Vec[ L1Index ].Eta();
      //      double L1Phi    = L1Vec[ L1Index ].Phi();

      // ********************************************************************************
      // *                         Turn on curves (No matching)                         *
      // ********************************************************************************
      
      if ( L1Index < int(jetIndexLabel.size()) ){

        // Get the current jet label, account for jet labeling starting at 0. LeadJet = 0, ...                                                                
        int jetIndex = L1Index;
        TString jetLabel = "_" + jetIndexLabel[ jetIndex ];

	// New rate study 
	bool passAk5ge30(false), passAk5ge40(false), passAk5ge50(false);
	if ( ak5Pt >= 30){
	  passAk5ge30 = true;

	  if ( ak5Pt >= 40){
	    passAk5ge40 = true;

	    if ( ak5Pt >= 50){
	      passAk5ge50 = true;

	    }
	  }
	}
	
	histEff[ matchPre + "_NoMatch_ak5ge30_Eff_vs_L1Thresh" + jetLabel ]->Fill( passAk5ge30, L1Pt );
	histEff[ matchPre + "_NoMatch_ak5ge40_Eff_vs_L1Thresh" + jetLabel ]->Fill( passAk5ge40, L1Pt );
	histEff[ matchPre + "_NoMatch_ak5ge50_Eff_vs_L1Thresh" + jetLabel ]->Fill( passAk5ge50, L1Pt );

      }
           
      
      // REALLY SUPER INEFFICIENT...
      for ( uint iPT = 0; iPT < onlineTOSingleJetPT.size(); ++iPT ){
	
	// Get the current offline PT cut 
	int onlinePT        = onlineTOSingleJetPT[ iPT ];
	TString onlinePTStr = Form("%d", onlinePT );
	bool passesPTTrig   = L1Pt > onlinePT;

	// TODO: ACCOUNT FOR MISMATCHED RECO JETS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	histEff[ matchPre + "_NoMatch_PT_TurnOn_OnPT_gt_" + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );

	if ( L1Index < int(jetIndexLabel.size()) ){

	  // Get the current jet label, account for jet labeling starting at 0. LeadJet = 0, ... 
	  int jetIndex = L1Index;
	  TString jetLabel = "_" + jetIndexLabel[ jetIndex ];

	  histEff[ matchPre + "_NoMatch" + jetLabel + "_PT_TurnOn_OnPT_gt_"   + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );

	}
	
	
      }
    } // End ak5 jet requirement

  } // End L1 matched loop









}


void 
JetHist::turnOnAk5Loop( std::vector<TLorentzVector> const &L1Vec, std::vector<TLorentzVector> const &ak5Vec, TString matchPre ){




  // *********************************************************************************************
  // *                                   Iterate over ak5 jets                                   *
  // *********************************************************************************************

  // Efficiency relative to all offline jets

  // Iterate over ak5 jets
  std::vector<pair_info> pairsAk = make_pairs( L1Vec, ak5Vec );
  std::sort(pairsAk.begin(), pairsAk.end(), sortDR);
  std::vector<int> ak5MatchedIndex = analyse_pairs(pairsAk, ak5Vec.size(), maxDeltaR);


  SUBPRINT("Ak5 jet matching")

    // Loop over all ak5 jets
    for(unsigned int iRank = 0; iRank < ak5MatchedIndex.size(); iRank++) {

      int ak5Index = iRank;
      int L1Index  = ak5MatchedIndex[iRank];   
   
      double ak5Pt    = ak5Vec[ ak5Index ].Pt();
      double ak5Eta   = ak5Vec[ ak5Index ].Eta();
      //     double ak5Phi   = ak5Vec[ ak5Index ].Phi();


      int unmatchedL1Index = -1; // Initially set to no jet
      // Check the corresponding 'ith' jet exists online
      if (iRank < L1Vec.size()){
	unmatchedL1Index = iRank; // L1 jet exists
      }
      
      double unmatchedL1Pt = -1;
      if ( unmatchedL1Index != -1){
	unmatchedL1Pt     = L1Vec[ L1Index ].Pt();
      }


//       double L1Pt = -1;
//       if ( L1Index != -1){
// 	L1Pt     = L1Vec[ L1Index ].Pt();
//       }


      // ****************************************************************************************************
      // *                                            Unmatched jets                                        *
      // ****************************************************************************************************
      
      // Extracted the corresponding ranked jet online if it exists. No matching requirement is imposed just that the jets
      // online and offline are of the same rank. Avoids problems associated with mismatching but potentially suffers from
      // the problem of assessing the performance associated with the wrong jet. Online however we only care that if an
      // event fires it has a jet offline that it claimed to have triggered.
      
   
      if ( iRank < jetIndexLabel.size() ){
        // Get the current jet label, account for jet labeling starting at 0. LeadJet = 0, ...                                                                
        int jetIndex = iRank;
        TString jetLabel = "_" + jetIndexLabel[ jetIndex ];

	
	// ********************************************************************************
	// *                         Turn on curves (No matching)                         *
	// ********************************************************************************


	// Calculate efficiency to trigger an (unmatched) 12 GeV L1 jet given an offline jet is reconstructed
	bool passUnmatchedL1ge12(false);
	if (unmatchedL1Pt >= 12){
	  passUnmatchedL1ge12 = true;
	}
	if ( ak5Pt >= 30){
	  histEff[ matchPre + "_NoMatch_ak5ge30_Eff_vs_OffEta"           + jetLabel ]->Fill( passUnmatchedL1ge12, ak5Eta );
	  histEff[ matchPre + "_NoMatch_ak5ge30_Eff_vs_OffEtaRegion"     + jetLabel ]->Fill( passUnmatchedL1ge12, ak5Eta );
	  histEff[ matchPre + "_NoMatch_ak5ge30_Eff_vs_NVTX"             + jetLabel ]->Fill( passUnmatchedL1ge12, NVTX );
	  if ( ak5Pt >= 40){
	    histEff[ matchPre + "_NoMatch_ak5ge40_Eff_vs_OffEta"         + jetLabel ]->Fill( passUnmatchedL1ge12, ak5Eta );
	    histEff[ matchPre + "_NoMatch_ak5ge40_Eff_vs_OffEtaRegion"   + jetLabel ]->Fill( passUnmatchedL1ge12, ak5Eta );
	    histEff[ matchPre + "_NoMatch_ak5ge40_Eff_vs_NVTX"           + jetLabel ]->Fill( passUnmatchedL1ge12, NVTX );
	    if ( ak5Pt >= 50){
	      histEff[ matchPre + "_NoMatch_ak5ge50_Eff_vs_OffEta"       + jetLabel ]->Fill( passUnmatchedL1ge12, ak5Eta );
	      histEff[ matchPre + "_NoMatch_ak5ge50_Eff_vs_OffEtaRegion" + jetLabel ]->Fill( passUnmatchedL1ge12, ak5Eta );
	      histEff[ matchPre + "_NoMatch_ak5ge50_Eff_vs_NVTX"         + jetLabel ]->Fill( passUnmatchedL1ge12, NVTX );
	    }
	  }
	}
      
      }


//       // ****************************************************************************************************
//       // *                                             Matched jets                                         *
//       // ****************************************************************************************************
//       //
//       // Require a matching between online and offline jets. Rank matched jets by corresponding offline jet rank.
      
      
//       if ( L1Index < int(jetIndexLabel.size()) ){
//         // Get the current jet label, account for jet labeling starting at 0. LeadJet = 0, ...                                                                
//         int jetIndex = iRank;
//         TString jetLabel = "_" + jetIndexLabel[ jetIndex ];
		

// 	// ********************************************************************************
// 	// *                                Turn on curves (Matched)                      *
// 	// ********************************************************************************
      
// 	// Calculate efficiency to trigger a (matched) 12 GeV L1 jet given an offline jet is reconstructed
// 	bool passMatchedL1ge12(false);
// 	if (L1Pt >= 12){
// 	  passMatchedL1ge12 = true;
// 	}
// 	if ( ak5Pt >= 30){
// 	  histEff[ matchPre + "_Matched_ak5ge30_Eff_vs_OffEta" + jetLabel ]      ->Fill( passMatchedL1ge12, ak5Eta );
// 	  histEff[ matchPre + "_Matched_ak5ge30_Eff_vs_OffEtaRegion" + jetLabel ]->Fill( passMatchedL1ge12, ak5Eta );
// 	  histEff[ matchPre + "_Matched_ak5ge30_Eff_vs_NVTX" + jetLabel ]        ->Fill( passMatchedL1ge12, NVTX );
// 	  if ( ak5Pt >= 40){
// 	    histEff[ matchPre + "_Matched_ak5ge40_Eff_vs_OffEta" + jetLabel ]      ->Fill( passMatchedL1ge12, ak5Eta );
// 	    histEff[ matchPre + "_Matched_ak5ge40_Eff_vs_OffEtaRegion" + jetLabel ]->Fill( passMatchedL1ge12, ak5Eta );
// 	    histEff[ matchPre + "_Matched_ak5ge40_Eff_vs_NVTX" + jetLabel ]        ->Fill( passMatchedL1ge12, NVTX );
// 	    if ( ak5Pt >= 50){
// 	      histEff[ matchPre + "_Matched_ak5ge50_Eff_vs_OffEta" + jetLabel ]      ->Fill( passMatchedL1ge12, ak5Eta );
// 	      histEff[ matchPre + "_Matched_ak5ge50_Eff_vs_OffEtaRegion" + jetLabel ]->Fill( passMatchedL1ge12, ak5Eta );
// 	      histEff[ matchPre + "_Matched_ak5ge50_Eff_vs_NVTX" + jetLabel ]        ->Fill( passMatchedL1ge12, NVTX );
// 	    }
// 	  }
// 	}
	
      
//       } // End jet binning
    
      
//       // Turn on curves
//       for ( uint iPT = 0; iPT < onlineTOSingleJetPT.size(); ++iPT ){
	
// 	// Get the current offline PT cut                                                                                                                        
// 	int onlinePT        = onlineTOSingleJetPT[ iPT ];
// 	TString onlinePTStr = Form("%d", onlinePT );
// 	bool passesPTTrig   = false;
	
// 	// Check whether L1 jet (if matched) passed trigger threshold 
// 	if ( matchedJet == true ){
	  
// 	  double L1Pt     = L1Vec[ L1Index ].Pt();
	  
// 	  // Check whether the L1 jet passed the trigger
// 	  if (L1Pt > onlinePT){
// 	    passesPTTrig = true;
// 	  }
	  
// 	}
	
// 	// All jets
// 	histEff[ matchPre + "_PT_TurnOn_OnPT_gt_"             + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
	
// 	if (L1Index == 0){ // Lead jet
// 	  histEff[ matchPre + "_LeadJet_PT_TurnOn_OnPT_gt_"   + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
// 	}
// 	else if (L1Index == 1) { // Second jet
// 	  histEff[ matchPre + "_SecondJet_PT_TurnOn_OnPT_gt_" + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
// 	}
// 	else if (L1Index == 2) { // Third jet
// 	  histEff[ matchPre + "_ThirdJet_PT_TurnOn_OnPT_gt_"  + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
// 	}
// 	else if (L1Index == 3) { // Fourth jet
// 	  histEff[ matchPre + "_FourthJet_PT_TurnOn_OnPT_gt_" + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
// 	}
	
//       }




      
//     for ( uint iDeltaEta = 0; iDeltaEta < onlineTOLeadSecondJetDeltaEta.size(); ++iDeltaEta ){
      
//       // Get the current offline PT cut                                                                                                                        
//       int onlineDeltaEta        = onlineTOLeadSecondJetPT[ iDeltaEta ];
//       TString onlineDeltaEtaStr = Form("%d", onlineDeltaEta );
//       bool passesDeltaEtaTrig   = false;
     
//       // Check whether L1 jet (if matched) passed trigger threshold 
//       if ( matchedJet == true ){

// 	double L1Pt     = L1Vec[ L1Index ].Pt();
// 	double L1Eta    = L1Vec[ L1Index ].Eta();
// //       double L1Phi    = L1Vec[ L1Index ].Phi();

// 	// Check whether the L1 jet passed the trigger
// 	if (L1Pt > onlinePT){
// 	  passesPTTrig = true;
// 	}
	
//       }

//       histEff[ matchPre + "_PT_TurnOn_OnPT_gt_" + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
      
//     }



      
    } // End ak5matched loop
  



}

void 
JetHist::turnOnAk5MatchedLoop( std::vector<TLorentzVector> const &L1Vec, std::vector<TLorentzVector> const &ak5Vec, TString matchPre ){




  // *********************************************************************************************
  // *                                   Iterate over ak5 jets                                   *
  // *********************************************************************************************

  // Efficiency relative to all offline jets

  // Iterate over ak5 jets
  std::vector<pair_info> pairsAk = make_pairs( L1Vec, ak5Vec );
  std::sort(pairsAk.begin(), pairsAk.end(), sortDR);
  std::vector<int> ak5MatchedIndex = analyse_pairs(pairsAk, ak5Vec.size(), maxDeltaR);


  SUBPRINT("Ak5 jet matching")

    // Loop over all ak5 jets
    for(unsigned int iRank = 0; iRank < ak5MatchedIndex.size(); iRank++) {

      int ak5Index = iRank;
      int L1Index  = ak5MatchedIndex[iRank];   
      // Check whether L1 jet was matched
      bool matchedJet = false;
      if ( L1Index != -1 ){
	matchedJet = true;
      }
   
      double ak5Pt    = ak5Vec[ ak5Index ].Pt();
      double ak5Eta   = ak5Vec[ ak5Index ].Eta();
      //     double ak5Phi   = ak5Vec[ ak5Index ].Phi();


      // ********************************************************************************
      // *                        Matching efficiency (all jets)                        *
      // ********************************************************************************
      
            
      histEff[ matchPre + "_matchingEffAllJet" ]->Fill( matchedJet, ak5Pt );

      // Leadjet
      if (iRank == 0){
	histEff[ matchPre + "_matchingEffLeadJet" ]->Fill( matchedJet, ak5Pt );
      }


//       int unmatchedL1Index = -1; // Initially set to no jet
//       // Check the corresponding 'ith' jet exists online
//       if (iRank < L1Vec.size()){
// 	unmatchedL1Index = iRank; // L1 jet exists
//       }
      
//       double unmatchedL1Pt = -1;
//       if ( unmatchedL1Index != -1){
// 	unmatchedL1Pt     = L1Vec[ L1Index ].Pt();
//       }


      double L1Pt = -1;
      if ( L1Index != -1){
	L1Pt     = L1Vec[ L1Index ].Pt();
      }



      // ****************************************************************************************************
      // *                                             Matched jets                                         *
      // ****************************************************************************************************
      //
      // Require a matching between online and offline jets. Rank matched jets by corresponding offline jet rank.
      
      
      if ( L1Index < int(jetIndexLabel.size()) ){
        // Get the current jet label, account for jet labeling starting at 0. LeadJet = 0, ...                                                                
        int jetIndex = iRank;
        TString jetLabel = "_" + jetIndexLabel[ jetIndex ];
		

	// ********************************************************************************
	// *                                Turn on curves (Matched)                      *
	// ********************************************************************************
      
	// Calculate efficiency to trigger a (matched) 12 GeV L1 jet given an offline jet is reconstructed
	bool passMatchedL1ge12(false);
	if (L1Pt >= 12){
	  passMatchedL1ge12 = true;
	}
	if ( ak5Pt >= 30){
	  histEff[ matchPre + "_Matched_ak5ge30_Eff_vs_OffEta" + jetLabel ]      ->Fill( passMatchedL1ge12, ak5Eta );
	  histEff[ matchPre + "_Matched_ak5ge30_Eff_vs_OffEtaRegion" + jetLabel ]->Fill( passMatchedL1ge12, ak5Eta );
	  histEff[ matchPre + "_Matched_ak5ge30_Eff_vs_NVTX" + jetLabel ]        ->Fill( passMatchedL1ge12, NVTX );
	  if ( ak5Pt >= 40){
	    histEff[ matchPre + "_Matched_ak5ge40_Eff_vs_OffEta" + jetLabel ]      ->Fill( passMatchedL1ge12, ak5Eta );
	    histEff[ matchPre + "_Matched_ak5ge40_Eff_vs_OffEtaRegion" + jetLabel ]->Fill( passMatchedL1ge12, ak5Eta );
	    histEff[ matchPre + "_Matched_ak5ge40_Eff_vs_NVTX" + jetLabel ]        ->Fill( passMatchedL1ge12, NVTX );
	    if ( ak5Pt >= 50){
	      histEff[ matchPre + "_Matched_ak5ge50_Eff_vs_OffEta" + jetLabel ]      ->Fill( passMatchedL1ge12, ak5Eta );
	      histEff[ matchPre + "_Matched_ak5ge50_Eff_vs_OffEtaRegion" + jetLabel ]->Fill( passMatchedL1ge12, ak5Eta );
	      histEff[ matchPre + "_Matched_ak5ge50_Eff_vs_NVTX" + jetLabel ]        ->Fill( passMatchedL1ge12, NVTX );
	    }
	  }
	}
	
      
      } // End jet binning
    
      
      // Turn on curves
//       for ( uint iPT = 0; iPT < onlineTOSingleJetPT.size(); ++iPT ){
	
// 	// Get the current offline PT cut                                                                                                                        
// 	int onlinePT        = onlineTOSingleJetPT[ iPT ];
// 	TString onlinePTStr = Form("%d", onlinePT );
// 	bool passesPTTrig   = false;
	
// 	// Check whether L1 jet (if matched) passed trigger threshold 
// 	if ( matchedJet == true ){
	  
// 	  double L1Pt     = L1Vec[ L1Index ].Pt();
	  
// 	  // Check whether the L1 jet passed the trigger
// 	  if (L1Pt > onlinePT){
// 	    passesPTTrig = true;
// 	  }
	  
// 	}
	
// 	// All jets
// 	histEff[ matchPre + "_PT_TurnOn_OnPT_gt_"             + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
	
// 	if (L1Index == 0){ // Lead jet
// 	  histEff[ matchPre + "_LeadJet_PT_TurnOn_OnPT_gt_"   + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
// 	}
// 	else if (L1Index == 1) { // Second jet
// 	  histEff[ matchPre + "_SecondJet_PT_TurnOn_OnPT_gt_" + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
// 	}
// 	else if (L1Index == 2) { // Third jet
// 	  histEff[ matchPre + "_ThirdJet_PT_TurnOn_OnPT_gt_"  + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
// 	}
// 	else if (L1Index == 3) { // Fourth jet
// 	  histEff[ matchPre + "_FourthJet_PT_TurnOn_OnPT_gt_" + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
// 	}
	
//       }





      
//     for ( uint iDeltaEta = 0; iDeltaEta < onlineTOLeadSecondJetDeltaEta.size(); ++iDeltaEta ){
      
//       // Get the current offline PT cut                                                                                                                        
//       int onlineDeltaEta        = onlineTOLeadSecondJetPT[ iDeltaEta ];
//       TString onlineDeltaEtaStr = Form("%d", onlineDeltaEta );
//       bool passesDeltaEtaTrig   = false;
     
//       // Check whether L1 jet (if matched) passed trigger threshold 
//       if ( matchedJet == true ){

// 	double L1Pt     = L1Vec[ L1Index ].Pt();
// 	double L1Eta    = L1Vec[ L1Index ].Eta();
// //       double L1Phi    = L1Vec[ L1Index ].Phi();

// 	// Check whether the L1 jet passed the trigger
// 	if (L1Pt > onlinePT){
// 	  passesPTTrig = true;
// 	}
	
//       }

//       histEff[ matchPre + "_PT_TurnOn_OnPT_gt_" + onlinePTStr ]->Fill( passesPTTrig, ak5Pt );
      
//     }



      
    } // End ak5matched loop
  



}


// Look at deltaR distribution of matched jets
void JetHist::studyMatchedJets( std::vector<TLorentzVector> const &L1Vec, std::vector<TLorentzVector> const &ak5Vec, TString matchPre, double maxDeltaR ){



//   // Perform deltaR matching                                                                                                                                
//   std::vector<pair_info> pairs = make_pairs( ak5Vec, L1Vec );
//   std::sort(pairs.begin(), pairs.end(), sortDR);
//   // No restriction on max deltaR
//   std::vector<int> L1MatchedIndex = analyse_pairs(pairs, L1Vec.size(), maxDeltaR );


//   for(unsigned int i = 0; i < L1MatchedIndex.size(); i++) {

//     int L1Index  = i;
//     int ak5Index = L1MatchedIndex[i];

//     if ( ak5Index == -1 ){
//       continue;
//     }





//   }


}




//define this as a plug-in
DEFINE_FWK_MODULE(JetHist);
