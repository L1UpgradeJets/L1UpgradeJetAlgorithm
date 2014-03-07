// -*- C++ -*-
//
// Package:    JetL1RECOMatchProducer
// Class:      JetL1RECOMatchProducer
// 
/**\class JetL1RECOMatchProducer JetL1RECOMatchProducer.cc AnalyseUpgradeJets/src/JetL1RECOMatchProducer.cc

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





//
// class declaration
//

class JetL1RECOMatchProducer : public edm::EDProducer {
   public:
      explicit JetL1RECOMatchProducer(const edm::ParameterSet&);
      ~JetL1RECOMatchProducer();

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


  virtual inline Double_t DeltaR(const l1extra::L1JetParticle & p1, const l1extra::L1JetParticle & p2);
  virtual inline bool sortDeltaR (const std::pair< l1extra::L1JetParticle, l1extra::L1JetParticle> jPair1,
				  const std::pair< l1extra::L1JetParticle, l1extra::L1JetParticle> jPair2 );


      // ----------member data ---------------------------
      edm::ParameterSet conf_;
      edm::Service<TFileService> fs;



      // ****************************************
      // *      Online-offline jet matching     *
      // ****************************************
    
      // Jet matching maximum delta R
      double maxDeltaR;


//       //  L1 pT calibration threshold, minimum L1 jet pT to apply correction  
//       double pTCalibrationThreshold;
//       // Width of iEta bins to sample to obtain pT corrections
//       double iEtaCalibBinWidth;
//       // Jet calibration iEta-binned pT corrections
//       std::vector <double> iEtaPtCorrectionPrePUSak5PrePUS, iEtaPtOffsetPrePUSak5PrePUS;
//       std::vector <double> iEtaPtCorrectionPrePUSak5PUS,    iEtaPtOffsetPrePUSak5PUS;
//       std::vector <double> iEtaPtCorrectionPUSak5PUS,       iEtaPtOffsetPUSak5PUS;


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
JetL1RECOMatchProducer::JetL1RECOMatchProducer(const edm::ParameterSet& iConfig): conf_(iConfig)
{

  // PrePUSak5PrePUS
  produces <l1slhc::L1TowerJetCollection>(     "CalibratedTowerJetPrePUSak5PrePUSTower" );
  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPrePUSak5PrePUSL1Jet" );
  // PrePUSak5PUS
  produces <l1slhc::L1TowerJetCollection>(     "CalibratedTowerJetPrePUSak5PUSTower" );
  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPrePUSak5PUSL1Jet" );
  // PUSak5PrePUS
  produces <l1slhc::L1TowerJetCollection>(     "CalibratedTowerJetPUSak5PUSTower" );
  produces <l1extra::L1JetParticleCollection>( "CalibratedTowerJetPUSak5PUSL1Jet" );


//   // ****************************************************************************************************
//   // *                                 Load configuration parameters                                    *
//   // ****************************************************************************************************

//   // ************************************************************
//   // *                        Calibration                       *
//   // ************************************************************
//   SUBPRINT("Calibration")
    
//   pTCalibrationThreshold = iConfig.getParameter< double >("pTCalibrationThreshold");
//   iEtaCalibBinWidth = iConfig.getParameter< double >("iEtaCalibrationBinWidth");
    
//   // PrePUSak5PrePUS
//   iEtaPtCorrectionPrePUSak5PrePUS  = iConfig.getParameter< std::vector< double > >("pTCorrectionPrePUSak5PrePUS");
//   iEtaPtOffsetPrePUSak5PrePUS      = iConfig.getParameter< std::vector< double > >("pTOffsetPrePUSak5PrePUS");
//   // PrePUSak5PUS
//   iEtaPtCorrectionPrePUSak5PUS     = iConfig.getParameter< std::vector< double > >("pTCorrectionPrePUSak5PUS");
//   iEtaPtOffsetPrePUSak5PUS         = iConfig.getParameter< std::vector< double > >("pTOffsetPrePUSak5PUS");
//   // PUSak5PUS
//   iEtaPtCorrectionPUSak5PUS        = iConfig.getParameter< std::vector< double > >("pTCorrectionPUSak5PUS");
//   iEtaPtOffsetPUSak5PUS            = iConfig.getParameter< std::vector< double > >("pTOffsetPUSak5PUS");


  // ************************************************************
  // *                  Online-offline matching                 *
  // ************************************************************
  SUBPRINT("Jet matching")

  // Maximum delta R in online-offline jet matching
  maxDeltaR = iConfig.getParameter< double >("MaxDeltaR");


  
}


JetL1RECOMatchProducer::~JetL1RECOMatchProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
JetL1RECOMatchProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   bool evValid = true;


//    std::auto_ptr<l1slhc::L1TowerJetCollection>     outputCalibTowerJet_Tower( new l1slhc::L1TowerJetCollection() );
//    std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJet_L1Jet( new l1extra::L1JetParticleCollection() );

   // PrePUSak5PrePUS
//    std::auto_ptr<l1slhc::L1TowerJetCollection>     outputCalibTowerJetPrePUSak5PrePUS_Tower( new l1slhc::L1TowerJetCollection() );
//    std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPrePUSak5PrePUS_L1Jet( new l1extra::L1JetParticleCollection() );
//    // PrePUSak5PUS   
//    std::auto_ptr<l1slhc::L1TowerJetCollection>     outputCalibTowerJetPrePUSak5PUS_Tower( new l1slhc::L1TowerJetCollection() );
//    std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPrePUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );
//    // PUSak5PUS 
//    std::auto_ptr<l1slhc::L1TowerJetCollection>     outputCalibTowerJetPUSak5PUS_Tower( new l1slhc::L1TowerJetCollection() );
//    std::auto_ptr<l1extra::L1JetParticleCollection> outputCalibTowerJetPUSak5PUS_L1Jet( new l1extra::L1JetParticleCollection() );


   // ****************************************************************************************************
   // *                                             Handles                                              *
   // ****************************************************************************************************
 
//    // ************************************************************
//    // *                        Tower Jets                        *
//    // ************************************************************
//    SUBPRINT("Tower Jets")

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

  edm::Handle<l1extra::L1JetParticleCollection> PrePUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSTowerJetL1Jet"), PrePUS_L1Jet); 
  if(!PrePUS_L1Jet.isValid()){
    evValid = false;   
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSTowerJetL1Jet") << std::endl;
  }
  edm::Handle<l1extra::L1JetParticleCollection> Ak5PUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("Ak5CaloJetL1Jet"), Ak5PUS_L1Jet);
  if(!Ak5PUS_L1Jet.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("Ak5CaloJetL1Jet") << std::endl;
    evValid = false;
  }






   if( evValid ){
  





  std::vector<TLorentzVector> ak5Vec, L1Vec;
  //  std::vector<int>            L1iEtaVec;

  //           *********************************************
  //           *          Fill ak5 TLorentzVector          *
  //           *********************************************

  //  for (edm::View< reco::CaloJet >::const_iterator ak5_It = ak5Jets->begin(); ak5_It != ak5Jets->end(); ++ak5_It ){
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




//      std::pair< std::vector<int>, std::vector<int> > l1RecoMatchedIndexPair;


//      std::vector<pair_info> pairs = make_pairs( ak5Vec, L1Vec );
//      //     std::sort(pairs.begin(), pairs.end(), sortDR);


//      // Matched ak5 index for every L1 jet in the event, -1 == no match within deltaR range.
//      std::vector<int> L1MatchedIndex = analyse_pairs(pairs, L1Vec.size(), maxDeltaR);



     










     std::vector< std::pair< l1extra::L1JetParticle, l1extra::L1JetParticle> > l1RecoJetPair;


     //     edm::Handle<l1extra::L1JetParticleCollection> const& L1Jets;

     for ( l1extra::L1JetParticleCollection::const_iterator L1_It = PrePUS_L1Jet->begin(); L1_It != PrePUS_L1Jet->end(); ++L1_It ){
       for ( l1extra::L1JetParticleCollection::const_iterator Ak5_It = Ak5PUS_L1Jet->begin(); Ak5_It != Ak5PUS_L1Jet->end(); ++Ak5_It ){


//        double L1Pt         = L1_It->p4().Pt();
//        double L1Eta        = L1_It->p4().Eta();
//        double L1Phi        = L1_It->p4().Phi();
//        double L1M          = L1_It->p4().M();

//        double Ak5Pt        = Ak5_It->p4().Pt();
//        double Ak5Eta       = Ak5_It->p4().Eta();
//        double Ak5Phi       = Ak5_It->p4().Phi();
//        double Ak5M         = Ak5_It->p4().M();

//        std::cout << "pT = " << L1Pt  << "\nEta = " << L1Eta  << "\tPhi = " << L1Phi  << "\n";
//        std::cout << "pT = " << Ak5Pt << "\nEta = " << Ak5Eta << "\tPhi = " << Ak5Phi << "\n";
//        std::cout << "dR = " << DeltaR( *L1_It, *Ak5_It ) << "\n\n";


	 l1RecoJetPair.push_back( std::make_pair(*L1_It, *Ak5_It) );
       }
     }



     for ( uint iL1RECO = 0; iL1RECO < l1RecoJetPair.size(); iL1RECO++ ){

       std::cout << DeltaR( l1RecoJetPair[ iL1RECO ].first, l1RecoJetPair[ iL1RECO ].second ) << "\n";

//        TLorentzVector L1Jet;
//        L1Jet.SetPtEtaPhiM( L1Pt, L1Eta, L1Phi, L1M );
//        L1Vec.push_back( L1Jet );

       
     }  
     //     std::sort( l1RecoJetPair.begin(), l1RecoJetPair.end(), sortDeltaR);
     std::cout << "\n";



//      std::vector< std::pair< *L1Jet, *L1Jet >


//   for(unsigned int i=0; i<gen_vec.size(); i++) {
//     for(unsigned int j=0; j<reco_vec.size(); j++) {
//       pinfo.push_back(  pair_info( i, j, gen_vec[i].DeltaR(reco_vec[j]) )  );
//     }
//   }

//   return pinfo;






//     std::vector<int> L1MatchedIndex;

//     // Build all possible pairs of RECO and L1 jets, sort in ascending order of dR
//     std::vector<pair_info> pairs = make_pairs( ak5Vec, L1Vec );
//     std::sort(pairs.begin(), pairs.end(), sortDR);
//     L1MatchedIndex = analyse_pairs(pairs, L1Vec.size(), maxDeltaR);
    
    

//     for(unsigned int i = 0; i < L1MatchedIndex.size(); i++) {
      
//       int L1Index  = i;
//       int ak5Index = L1MatchedIndex[i];
      
//       if ( ak5Index == -1 ){
// 	//      std::cout << "No matched jet!!!\n";
// 	continue;
//       }
//     }



 
   }

}

// ------------ method called once each job just before starting event loop  ------------
void 
JetL1RECOMatchProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetL1RECOMatchProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
JetL1RECOMatchProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetL1RECOMatchProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetL1RECOMatchProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetL1RECOMatchProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetL1RECOMatchProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}




// Jet calibration
// void 
// JetL1RECOMatchProducer::calibrateJet( edm::Handle<l1slhc::L1TowerJetCollection> const& UncalibJet_Tower, std::vector< double > const& iEtaPtOffset, 
// 				std::vector< double > const& iEtaPtCorrection, 
// 				std::auto_ptr<l1slhc::L1TowerJetCollection>     &outputCalibTowerJet_Tower,
// 				std::auto_ptr<l1extra::L1JetParticleCollection> &outputCalibTowerJet_L1Jet ){




//      // ****************************************************************************************************
//      // *                                         Calibrate Jets                                           *
//      // ****************************************************************************************************
//      // Calibrate jets with LUT

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



// }


inline Double_t JetL1RECOMatchProducer::DeltaR(const l1extra::L1JetParticle & p1, const l1extra::L1JetParticle & p2){

//inline Double_t LorentzVector::DeltaR(const LorentzVector & v) const {
  Double_t deta = p1.p4().Eta() - p2.p4().Eta();
  Double_t dphi = TVector2::Phi_mpi_pi( p1.p4().Phi() - p2.p4().Phi() );
  return TMath::Sqrt( deta*deta+dphi*dphi );
}


// Return whether jetPair1 < jetPair2
inline bool JetL1RECOMatchProducer::sortDeltaR(const std::pair< l1extra::L1JetParticle, l1extra::L1JetParticle> jPair1, 
			       const std::pair< l1extra::L1JetParticle, l1extra::L1JetParticle> jPair2 ){

  double dR1 = DeltaR( jPair1.first, jPair1.second );
  double dR2 = DeltaR( jPair2.first, jPair2.second );

  return dR1 < dR2;

}





//define this as a plug-in
DEFINE_FWK_MODULE(JetL1RECOMatchProducer);
