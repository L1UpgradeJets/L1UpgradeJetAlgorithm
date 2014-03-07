// -*- C++ -*-
//
// Package:    JetCalibHist
// Class:      JetCalibHist
// 
/**\class JetCalibHist JetCalibHist.cc AnalyseUpgradeJets/src/JetCalibHist.cc

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
#include "AnalyseUpgradeJets/AnalyzeJets/interface/histBooker.h"
//#include "AnalyseUpgradeJets/AnalyzeJets/interface/JetMatch.h"
// #include "AnalyseUpgradeJets/AnalyzeJets/src/helperFunctions.cc"

#include <algorithm>  // for sorting



// Print debugging messages
//#define VERBOSE




//
// class declaration
//

//class JetCalibHist : public edm::EDProducer {
class JetCalibHist : public edm::EDAnalyzer {
   public:
      explicit JetCalibHist(const edm::ParameterSet&);
      ~JetCalibHist();

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
        
  
  virtual void matchOnlineOfflineJets( edm::Handle<l1slhc::L1TowerJetCollection> const& L1Jets, edm::Handle<l1extra::L1JetParticleCollection> const& ak5Jets,
				       TString prefix, int maxJetsMatched, double ak5PTthreshold );


      // ----------member data ---------------------------
      edm::ParameterSet conf_;
      edm::Service<TFileService> fs;


      // Histogram containers
      std::map< TString, TH1*> hist1D;
      std::map< TString, TH1*> hist2D;

      // Number of reconstructed ak5 primary vertices 
      int NVTX;


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

      // Width of iEta bins to sample to obtain pT corrections
      double iEtaCalibBinWidth;





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
JetCalibHist::JetCalibHist(const edm::ParameterSet& iConfig): conf_(iConfig)
{


  // Associate the histogram container with the histogram booker, a class for handling histogram booking
  histBooker booker( &hist1D, &hist2D );



  // ****************************************************************************************************
  // *                                 Load configuration parameters                                    *
  // ****************************************************************************************************

//   // L1 jet cleaning parameters
//   minL1JetPt  = iConfig.getParameter<double> ("minL1JetPt");
//   maxL1JetEta = iConfig.getParameter<double> ("maxL1JetEta");
//   // RECO jet cleaning parameters
//   minRECOJetPt  = iConfig.getParameter<double> ("minRECOJetPt");
//   maxRECOJetEta = iConfig.getParameter<double> ("maxRECOJetEta");


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



  // ************************************************************
  // *                  Online-offline matching                 *
  // ************************************************************
  SUBPRINT("Jet matching")

  // Maximum delta R in online-offline jet matching
  maxDeltaR = iConfig.getParameter< double >("MaxDeltaR");

  iEtaCalibBinWidth = iConfig.getParameter< double >("iEtaCalibrationBinWidth");




  // ***************************************************
  // *            Collections to be calibrated         *
  // ***************************************************

  // PrePUS-ak5PUS
  matchPrefix.push_back("PrePUS_ak5PUS");
  matchLabel.push_back("TowerJet PrePUS, Ak5 PUS");
  // PUS-ak5PUS
  matchPrefix.push_back("PUS_ak5PUS");
  matchLabel.push_back("TowerJet PUS, Ak5 PUS");
  // LPUS-ak5PUS
  matchPrefix.push_back("LPUS_ak5PUS");
  matchLabel.push_back("TowerJet PUS, Ak5 PUS");




    // ********************************************************************************
    // *                                 Jet Matching                                 *
    // ********************************************************************************
    PRINT("Jet matching")
    for (unsigned int iMatch = 0; iMatch < matchPrefix.size(); iMatch++){
    
      TString matchPre = matchPrefix[iMatch];
      TString matchLab = matchLabel[iMatch];
      
      // Create a match subsubdirectory
      TFileDirectory matchSubDir = fs->mkdir( matchPre.Data() );

      
    
      // pT Correlations
      booker.book2DTProf( matchPre + "_L1PT_vs_OffPT" ,     matchSubDir, "L1 P_{T} vs Offline P_{T};Offline p_{T} (GeV);L1 P_{T} (GeV)", pOffPT, pOffPT);
      booker.book2DTProf( matchPre + "_OffPT_vs_L1PT" ,     matchSubDir, "Offline P_{T} vs L1 P_{T};L1 P_{T} (GeV);Offline p_{T} (GeV)", pOffPT, pOffPT);

    
      // ************************************************************
      // *                 iEta binned distributions                *
      // ************************************************************

      SUBPRINT("iEta binned")
      TFileDirectory iEtaDir = matchSubDir.mkdir( "iEtaBinned" );

      TString previousiEtaBin = "";

      // SHOULD CUT THIS TO iETA <= 28 - (jetSize - 1)
      for (int iEta = -28;iEta <= 28;iEta++){
	
	if (iEta == 0)
	  continue;

	// Get current iEta bin
      // ****************************************
	TString iEtaStr;
      
	// iEta binned distributions    
	if (iEtaCalibBinWidth == 1){ // Individual iEta bins
	  iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEta)));
	}
	else{ // iEta bins of width greater than 1 
	  
	  int iEtaBinAbs = (abs(iEta) - 1)/iEtaCalibBinWidth;
	  int sign = 1;
	  if (iEta < 0)
	    sign = -1;
	  
	  double iEtaBinLowAbs  = ( (iEtaBinAbs)*iEtaCalibBinWidth + 1 );
	  double iEtaBinHighAbs = ( (iEtaBinAbs + 1)*iEtaCalibBinWidth );
	  
	  // Fix bin limits
	  if (iEtaBinLowAbs > 28)
	    iEtaBinLowAbs = 28;
	  if (iEtaBinHighAbs > 28)
	    iEtaBinHighAbs = 28;

	  double iEtaBinLow  = sign*iEtaBinLowAbs;
	  double iEtaBinHigh = sign*iEtaBinHighAbs;
	  
	
	  if (iEtaBinLow > iEtaBinHigh)
	    std::swap(iEtaBinLow, iEtaBinHigh);
	  
	  iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEtaBinLow))) + "to" + TString(Form("%d",Int_t(iEtaBinHigh)));

	  if (iEtaStr == previousiEtaBin){
	    continue;   // Avoid booking the same histogram again
	  }
	  else{
	    previousiEtaBin = iEtaStr;
	  }
	  
	}

	
	TFileDirectory iEtaSubDir = iEtaDir.mkdir( iEtaStr.Data() );



	// pT Correlations
	booker.book2DTProf( matchPre + "_L1PT_vs_OffPT_"     + iEtaStr, iEtaSubDir, "L1 P_{T} vs Offline P_{T} (" + iEtaStr + ");Offline p_{T} (GeV);L1 P_{T} (GeV)", 										 pOffPT, pOffPT);
	booker.book2DTProf( matchPre + "_OffPT_vs_L1PT_" + iEtaStr, iEtaSubDir, "Offline P_{T} vs L1 P_{T} (" + iEtaStr + 
			  ");L1 P_{T} (GeV);Offline p_{T} (GeV)", pOffPT, pOffPT);


      } // End iEta loop
      
    } // End matched jets
    
}



JetCalibHist::~JetCalibHist()
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
//JetCalibHist::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
void
JetCalibHist::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   bool evValid = true;




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
  // *                       Current Jets                       *
  // ************************************************************
//   SUBPRINT("Current Jets")
    
//     // Get L1 jets to be used, currently Cental and Tau jets
//     std::vector <InputTag> l1extraparticles = conf_.getParameter< std::vector < InputTag > >("extrajet");
  
  
//     Handle<l1extra::L1JetParticleCollection> Current_L1Jet;
//     // Loop through Central and tau jets
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

  // Pre-PU subtracted jets
  edm::Handle<l1slhc::L1TowerJetCollection> PrePUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSubTowerJet"), PrePUS_Tower);
  if(!PrePUS_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSubTowerJet") << std::endl;
  }
  // Global PU subtracted tower jets 
  edm::Handle<l1slhc::L1TowerJetCollection> PUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUSubTowerJet"), PUS_Tower);
  if(!PUS_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUSubTowerJet") << std::endl;
  }
  // Local PU subtracted tower jets 
  edm::Handle<l1slhc::L1TowerJetCollection> LPUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("LocalPUSubTowerJet"), LPUS_Tower);
  if(!LPUS_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("LocalPUSubTowerJet") << std::endl;
  }





  // ************************************************************
  // *                        RECO Jets                         *
  // ************************************************************
  SUBPRINT("RECO Jets")

  //PU subtracted AK5 calo jets-must be in root file read in
//   edm::Handle<reco::CaloJetCollection> Ak5PUS;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUsubCaloJets"), Ak5PUS);
//   if(!Ak5PUS.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUsubCaloJets") << std::endl;
//     evValid = false;
//   }
  edm::Handle<l1extra::L1JetParticleCollection> Ak5PUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("Ak5CaloJetL1Jet"), Ak5PUS_L1Jet);
  if(!Ak5PUS_L1Jet.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("Ak5CaloJetL1Jet") << std::endl;
    evValid = false;
  }

//   // uncorrected jets!
//   edm::Handle<edm::View< reco::CaloJet > > Ak5PrePUS;
//   //edm::Handle<reco::CaloJetCollection> PrePUSAk5Jets;
//   iEvent.getByLabel("ak5CaloJets", Ak5PrePUS );
//   if(!Ak5PrePUS.isValid()){
//     edm::LogWarning("MissingProduct") << "ak5CaloJets"<< std::endl;
//     evValid = false;
//   }
//   edm::Handle<l1extra::L1JetParticleCollection> Ak5PrePUS_L1Jet;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("Ak5CaloUncorrJetL1Jet"), Ak5PrePUS_L1Jet);
//   if(!Ak5PrePUS_L1Jet.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("Ak5CaloUncorrJetL1Jet") << std::endl;
//     evValid = false;
//   }



  if (evValid){
  
   
    
    // ****************************************************************************************************
    // *                                   L1-Offline lead jet matching                                   *
    // ****************************************************************************************************


      // ******************************** 
      // *            ak5 PUS           * 
      // ******************************** 

      // PrePUS-ak5PUS
      matchOnlineOfflineJets(  PrePUS_Tower, Ak5PUS_L1Jet, "PrePUS_ak5PUS", 1, 30 );
      // PUS-ak5PUS
      matchOnlineOfflineJets(  PUS_Tower,    Ak5PUS_L1Jet, "PUS_ak5PUS", 1, 30 );
      // LPUS-ak5PUS
      matchOnlineOfflineJets(  LPUS_Tower,   Ak5PUS_L1Jet, "LPUS_ak5PUS", 1, 30 );



  }









   QUIT("Successfully completed analysis of single event")
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetCalibHist::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetCalibHist::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
JetCalibHist::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetCalibHist::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetCalibHist::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetCalibHist::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetCalibHist::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}






// ****************************************************************************************************
// *                                    Match Online to Offline Jets                                  *
// **************************************************************************************************** 
void 
JetCalibHist::matchOnlineOfflineJets( edm::Handle<l1slhc::L1TowerJetCollection> const& L1Jets, edm::Handle<l1extra::L1JetParticleCollection> const& ak5Jets,
				      TString prefix, int maxJetsMatched, double ak5PTthreshold ){

//   PRINT("Online-Offline Jet matching")

//   // Find which histogram prefix to use
//   int preIndex = -1;
//   for (unsigned int iMatch = 0; iMatch < matchPrefix.size(); iMatch++){
//     if (prefix == matchPrefix[iMatch]){
//       preIndex = iMatch;
//       break;
//     }
//   }
//   if (preIndex == -1){
//     std::cout << "ERROR - Prefix '" << prefix << "' was not found.\n";
//     return;
//   } 
//   TString matchPre = matchPrefix[preIndex];
//   TString matchLab = matchLabel[preIndex];

//   SUBPRINT("Extracted prefix")

//   //      ********************************************************************************
//   //      *                      PrePUS Online-Offline jet matching                      *
//   //      ********************************************************************************
//   //


//   std::vector<TLorentzVector> ak5Vec, L1Vec;
//   std::vector<int>            L1iEtaVec;

//   //           *********************************************
//   //           *          Fill ak5 TLorentzVector          *
//   //           *********************************************

//   for ( l1extra::L1JetParticleCollection::const_iterator ak5_It = ak5Jets->begin(); ak5_It != ak5Jets->end(); ++ak5_It ){


//     double ak5Pt         = ak5_It->p4().Pt();
//     double ak5Eta        = ak5_It->p4().Eta();
//     double ak5Phi        = ak5_It->p4().Phi();
//     double ak5M          = ak5_It->p4().M();

//     // ========================================
//     // Jet quality cuts
//     // ========================================

//     //    if( ak5Pt < MIN_OFF_PT )   continue;
//     if( fabs( ak5Eta ) > 3.0 ) continue;
    
//     if( ak5Pt < ak5PTthreshold ) continue;

//     // ========================================

//     TLorentzVector ak5Jet;
//     ak5Jet.SetPtEtaPhiM( ak5Pt, ak5Eta, ak5Phi, ak5M );
//     ak5Vec.push_back( ak5Jet );

//   }


//   //           *********************************************
//   //           *          Fill L1 TLorentzVector           *
//   //           *********************************************


//   for (l1slhc::L1TowerJetCollection::const_iterator L1_It = L1Jets->begin(); L1_It != L1Jets->end(); ++L1_It ){

//     double L1Pt         = L1_It->p4().Pt();
//     double L1Eta        = L1_It->p4().Eta();
//     double L1Phi        = L1_It->p4().Phi();
//     double L1M          = L1_It->p4().M();

//     TLorentzVector L1Jet;
//     L1Jet.SetPtEtaPhiM( L1Pt, L1Eta, L1Phi, L1M );
//     L1Vec.push_back( L1Jet );

//     L1iEtaVec.push_back(L1_It->iEta());

//   }



//   // ------------------------------ Perform jet matching ------------------------------


//   SUBPRINT("Extracted jets")

//     // THIS IS REALLY INEFFICIENT. WHEN YOU HAVE TIME WRITE A CLASS THAT RETURNS std::vector< std::pair< TLorentz*, TLorentz*> >
//     // WITH SORT ALGORITHMS AS AN OPTION

//     std::vector<int> L1MatchedIndex;

//     // Perform deltaR matching, for calibration ONLY

//   std::vector<pair_info> pairs = make_pairs( ak5Vec, L1Vec );
//   std::sort(pairs.begin(), pairs.end(), sortDR);
//   L1MatchedIndex = analyse_pairs(pairs, L1Vec.size(), maxDeltaR);
  

//   // ____________________
//   // ALL JET MATCHING
//   // ____________________
//   SUBPRINT("All jet matching")
//   int jetsMatched = 0;

//   for(unsigned int i = 0; i < L1MatchedIndex.size(); i++) {

//     int L1Index  = i;
//     int ak5Index = L1MatchedIndex[i];

//     if ( ak5Index == -1 ){
//       //      std::cout << "No matched jet!!!\n";
//       continue;
//     }
//     jetsMatched++;
//     // Restrict to lead jet
//     if (jetsMatched > maxJetsMatched){
//       break;
//     }
     

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

    
//     double ak5Pt    = ak5Vec[ ak5Index ].Pt();
// //     double ak5Eta   = ak5Vec[ ak5Index ].Eta();
// //     double ak5Phi   = ak5Vec[ ak5Index ].Phi();

//     double L1Pt     = L1Vec[ L1Index ].Pt();
// //     double L1Eta    = L1Vec[ L1Index ].Eta();
// //     double L1Phi    = L1Vec[ L1Index ].Phi();


//     // ************************************************************
//     // *                         Unbinned                         *
//     // ************************************************************

//     // pT Correlations
//     hist2D[matchPre + "_L1PT_vs_OffPT"     ]    ->Fill(ak5Pt, L1Pt);
//     hist1D[matchPre + "_L1PT_vs_OffPT_prof"]    ->Fill(ak5Pt, L1Pt);

//     // ************************************************************
//     // *                          Binned                          *
//     // ************************************************************

//     // pT Correlations
//     hist2D[matchPre + "_L1PT_vs_OffPT_"     + iEtaStr]          ->Fill(ak5Pt, L1Pt);
//     hist1D[matchPre + "_L1PT_vs_OffPT_"     + iEtaStr + "_prof"]->Fill(ak5Pt, L1Pt);
    
//   }
  
  
}
 




//define this as a plug-in
DEFINE_FWK_MODULE(JetCalibHist);
