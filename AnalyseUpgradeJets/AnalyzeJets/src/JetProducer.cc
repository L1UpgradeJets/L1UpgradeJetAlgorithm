// -*- C++ -*-
//
// Package:    JetProducer
// Class:      JetProducer
// 
/**\class JetProducer JetProducer.cc AnalyseUpgradeJets/src/JetProducer.cc

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
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"



#include "DataFormats/JetReco/interface/JetID.h" 
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "SimDataFormats/SLHC/interface/EtaPhiContainer.h"
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/TriggerTowerGeometry.h"





#include "AnalyseUpgradeJets/AnalyzeJets/interface/printing.h"
// #include "AnalyseUpgradeJets/AnalyzeJets/src/JetMatch.h"
// #include "AnalyseUpgradeJets/AnalyzeJets/src/helperFunctions.cc"

#include <algorithm>  // for sorting

// Ranking function for sort
//bool caloJetRankDescending ( const reco::CaloJet jet1, const reco::CaloJet jet2 ){ return ( 1  ); }

bool CurrentL1JetRankDescending ( l1extra::L1JetParticle jet1, l1extra::L1JetParticle jet2 ){ return ( jet1.p4().Pt() > jet2.p4().Pt() ); }

const double PI = 3.141592654;


//
// class declaration
//

class JetProducer : public edm::EDProducer {
   public:
      explicit JetProducer(const edm::ParameterSet&);
      ~JetProducer();

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
      
      // Jet cleaning
      virtual bool cleanL1Jet(math::PtEtaPhiMLorentzVector jet);
      virtual bool cleanRECOJet(math::PtEtaPhiMLorentzVector jet);


      // ----------member data ---------------------------
      edm::ParameterSet conf_;
      edm::Service<TFileService> fs;


      // Histogram containers
      std::map< TString, TH1*> hist1D;
      std::map< TString, TH1*> hist2D;

      // L1 jet cleaning parameters
      double minL1JetPt;
      double maxL1JetEta;
      // RECO jet cleaning parameters
      double minRECOJetPt;
      double maxRECOJetEta;

      bool foldEta;
  
  // Eta region-level segmentation                                                                                                                        
  std::vector <double> RCTEtaRegions;


//       // Online jet collections
//       std::vector <TLorentzVector> currentJets;
      

//       std::vector <l1extra::L1JetParticle> prePUSTowerJets;
//   //      std::vector <TLorentzVector> prePUSTowerJets;
//       std::vector <TLorentzVector> PUSTowerJets;
//       std::vector <TLorentzVector> LPUSTowerJets;

//       // Offline jet collections
//       std::vector <TLorentzVector> ak5CaloJets;
//       std::vector <TLorentzVector> ak5CaloUncorrJets;
//       //std::vector <TLorentzVector> kt6Jets;



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
JetProducer::JetProducer(const edm::ParameterSet& iConfig): conf_(iConfig)
{


   // Eta region-level segmentation 
   RCTEtaRegions        = iConfig.getParameter< std::vector<double> >("RCTEtaRegions");
   sort (RCTEtaRegions.begin(), RCTEtaRegions.end());  // ensure the bins are in ascending order        



  produces<int>("NVTX");
//   produces<double>("prePusHT_Tower");
//    produces<double>("HT_Tower");
//    produces<double>("MHT_Tower");
//    produces<double>("PUS_Tower");
//  produces<int>("JetTTSize");


//  produces < std::vector<TLorentzVector> >("PrePUSTowerJetTLorentz");

  produces <l1extra::L1JetParticleCollection>("PrePUSTowerJetL1Jet");
  produces <l1extra::L1JetParticleCollection>("PUSTowerJetL1Jet");
  produces <l1extra::L1JetParticleCollection>("LPUSTowerJetL1Jet");
  produces <l1extra::L1JetParticleCollection>("CurrentUncalibJetL1Jet");
  produces <l1extra::L1JetParticleCollection>("CurrentJetL1Jet");

  produces <l1extra::L1JetParticleCollection>("PrePUSRawAk5CaloJetL1Jet");
  //produces <l1extra::L1JetParticleCollection>("PUSRawAk5CaloJetL1Jet");
  //produces <l1extra::L1JetParticleCollection>("PrePUSAk5CaloJetL1Jet");
  produces <l1extra::L1JetParticleCollection>("Ak5CaloJetL1Jet");
  //  produces <l1extra::L1JetParticleCollection>("Ak5CaloUncorrJetL1Jet");

//    produces < std::vector<TLorentzVector*> >("PrePUSTowerJetTLorentz");
//    produces < std::vector<TLorentzVector*> >("PUSTowerJetTLorentz");
//    produces < std::vector<TLorentzVector*> >("LPUSTowerJetTLorentz");
//    produces < std::vector<TLorentzVector*> >("CurrentJetTLorentz");
//    produces < std::vector<TLorentzVector*> >("Ak5CaloJetTLorentz");
//   produces < std::vector<TLorentzVector*> >("Ak5CaloUncorrJetTLorentz");






  // ****************************************************************************************************
  // *                                 Load configuration parameters                                    *
  // ****************************************************************************************************

  // L1 jet cleaning parameters
  minL1JetPt  = iConfig.getParameter<double> ("minL1JetPt");
  maxL1JetEta = iConfig.getParameter<double> ("maxL1JetEta");
  // RECO jet cleaning parameters
  minRECOJetPt  = iConfig.getParameter<double> ("minRECOJetPt");
  maxRECOJetEta = iConfig.getParameter<double> ("maxRECOJetEta");

  foldEta = iConfig.getParameter<bool> ("FoldEta");








  // ****************************************************************************************************
  // *                                      Vertex Distributions                                        *
  // ****************************************************************************************************

  PRINT("Vertex")

  TFileDirectory vtxDir = fs->mkdir( "Vertex" );
  

  TString prefix = "Vertex";

  hist1D[prefix + "_NVTX"]           = vtxDir.make<TH1F>(prefix + "_NVTX","Number of reconstructed ak5 primary vertices;N_{VTX};Entries", 300, 0., 300.);
  hist1D[prefix + "_Chi2"]           = vtxDir.make<TH1F>(prefix + "_Chi2","Vertex #chi^{2};#chi^{2};Entries",  201,-0.5,200.5);
  hist1D[prefix + "_Ndof"]           = vtxDir.make<TH1F>(prefix + "_Ndof","Vertex N_{DOF};N_{DOF};Entries",    201,-0.5,200.5);
  hist1D[prefix + "_Chi2_over_Ndof"] = vtxDir.make<TH1F>(prefix + "_Chi2_over_Ndof","Vertex #chi^{2}/N_{DOF};#chi^{2}/N_{DOF};Entries", 101,-0.5,10.5);
  hist2D[prefix + "_Ndof_vs_Chi2"]   = vtxDir.make<TH2D>(prefix + "_Ndof_vs_Chi2","Vertex Ndof vs #chi^{2} ;NDof;#chi^{2}", 201,-0.5,200.5, 201,-0.5,200.5);


  PRINT("Debug")

  TFileDirectory debugDir = fs->mkdir( "Debug" );
  hist1D["ak5_leadJet_minDeltaR"] = debugDir.make<TH1F>("ak5_leadJet_minDeltaR","Lead jet - jet minimum #DeltaR;Minimum #DeltaR;Entries", 301, 0.05, 3.05);

  hist2D["ak5_leadJet_minDeltaRPT"] = debugDir.make<TH2D>("ak5_leadJet_minDeltaRPT","Lead jet - Jet p_{T} vs jet minimum #DeltaR;Minimum #DeltaR;Jet p_{T}", 
							  301, 0.05, 3.05, 201, -0.5, 200.5);

  hist2D["ak5_leadJet_minDeltaRRelPT"] = debugDir.make<TH2D>("ak5_leadJet_minDeltaRRelPT","Lead jet - Relative Jet p_{T} vs jet minimum #DeltaR;Minimum #DeltaR;#Deltap_{T}/p_{T}", 
							  301, 0.05, 3.05, 601, -1.505, 4.505);


  hist1D["ak5_leadJet_JetDeltaR"] = debugDir.make<TH1F>("ak5_leadJet_JetDeltaR","Lead jet - jets #DeltaR;#DeltaR;Entries", 301, 0.05, 3.05);
  hist1D["ak5_leadJet_JetRelPT"]  = debugDir.make<TH1F>("ak5_leadJet_JetRelPT","Lead jet - Relative p_{T} of sum of jet p_{T} #DeltaR < 0.5;#Deltap_{T}/p_{T};Entries", 601, -3.005, 3.005);

  
}


JetProducer::~JetProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
JetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   bool evValid = true;





   std::auto_ptr<int>                      outputNVTX( new int() );
//   auto_ptr<double>                   outputprePusHT_Tower( new double() );
//   auto_ptr<double>                   outputHT_Tower(  new double() );
//   auto_ptr<double>                   outputMHT_Tower( new double() );

   std::auto_ptr<l1extra::L1JetParticleCollection> outputPrePUSTowerJet_L1Jet(   new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputPUSTowerJet_L1Jet(      new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputLPUSTowerJet_L1Jet(     new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputUncalibCurrentJet_L1Jet(new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputCurrentJet_L1Jet(       new l1extra::L1JetParticleCollection() );

   // RECO
   std::auto_ptr<l1extra::L1JetParticleCollection> outputPrePUSRawAk5CaloJet_L1Jet( new l1extra::L1JetParticleCollection() );
   //std::auto_ptr<l1extra::L1JetParticleCollection> outputPUSRawAk5CaloJet_L1Jet(    new l1extra::L1JetParticleCollection() );
   //std::auto_ptr<l1extra::L1JetParticleCollection> outputPrePUSAk5CaloJet_L1Jet(    new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputAk5CaloJet_L1Jet(          new l1extra::L1JetParticleCollection() );



//    std::auto_ptr< std::vector<TLorentzVector*> > outputPrePUSTowerJet_TLorentz(    new std::vector<TLorentzVector*>() );
//    std::auto_ptr< std::vector<TLorentzVector*> > outputPUSTowerJet_TLorentz(       new std::vector<TLorentzVector*>() );
//    std::auto_ptr< std::vector<TLorentzVector*> > outputLPUSTowerJet_TLorentz(      new std::vector<TLorentzVector*>() );
//    std::auto_ptr< std::vector<TLorentzVector*> > outputCurrentJet_TLorentz(        new std::vector<TLorentzVector*>() );
//    std::auto_ptr< std::vector<TLorentzVector*> > outputAk5CaloJet_TLorentz(        new std::vector<TLorentzVector*>() );
//    std::auto_ptr< std::vector<TLorentzVector*> > outputAk5CaloUncorrJet_TLorentz(  new std::vector<TLorentzVector*>() );



//    // Matched jets - TODO: MOVE TO A SEPERATE PRODUCER
//    std::auto_ptr< std::vector< std::pair<TLorentzVector*, TLorentzVector*> > > oPrePUS_Ak5PUSMatch( new std::vector< std::pair<TLorentzVector*, TLorentzVector*> >() );
//    std::auto_ptr< std::vector< std::pair<TLorentzVector*, TLorentzVector*> > > oPUS_Ak5PUSMatch(    new std::vector< std::pair<TLorentzVector*, TLorentzVector*> >() );
//    std::auto_ptr< std::vector< std::pair<TLorentzVector*, TLorentzVector*> > > oLPUS_Ak5PUSMatch(   new std::vector< std:
//												    :pair<TLorentzVector*, TLorentzVector*> >() );


   // ****************************************************************************************************
   // *                                             Handles                                              *
   // ****************************************************************************************************
     
   /*
   //Need this for information about PU
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RecoVertices"), vertices); 
   if(!vertices.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RecoVertices") << std::endl;
     evValid = false;
   }
   */


   // ********************************************************************************
   // *                                    Jets                                      *
   // ********************************************************************************



  // ************************************************************
  // *                       Current Jets                       *
  // ************************************************************
  SUBPRINT("Current Jets")

/*
    // Uncalibrated central jets
   edm::Handle<L1GctJetCandCollection>  L1GCTUncalibCenJet;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CurrentL1GCTUncalibCentralJet"), L1GCTUncalibCenJet );
   if(!L1GCTUncalibCenJet.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CurrentL1GCTUncalibCentralJet") << std::endl;
   }
   // Uncalibrated tau jets
   edm::Handle<L1GctJetCandCollection>  L1GCTUncalibTauJet;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CurrentL1GCTUncalibTauJet"), L1GCTUncalibTauJet );
   if(!L1GCTUncalibTauJet.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CurrentL1GCTUncalibTauJet") << std::endl;
   }
*/


    
   // Get L1 jets to be used, currently Cental and Tau jets
   // std::vector <InputTag> l1extraparticles = conf_.getParameter< std::vector < InputTag > >("extrajet");
   // Loop through Central and tau jets                                                                                                                         
   //for (uint i = 0; i < l1extraparticles.size(); ++i){
   //  Handle<l1extra::L1JetParticleCollection> currl1Jets;
   //  iEvent.getByLabel(l1extraparticles[i], currl1Jets );
   //  if(!currl1Jets.isValid()){
   //    evValid = false;
   //    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("extrajet") << std::endl;
   //  }
  // }


  // ************************************************************
  // *                        Tower Jets                        *
  // ************************************************************
  SUBPRINT("Tower Jets")

//   // Proto tower jets
//   edm::Handle<l1slhc::L1TowerJetCollection> Proto_Tower;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("ProtoTowerJet"), Proto_Tower);
//   if(!Proto_Tower.isValid()){
//     evValid = false;
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("ProtoTowerJet") << std::endl;
//   }
  
//   // Centrality filtered tower jets
//   edm::Handle<l1slhc::L1TowerJetCollection> FilteredCentrality_Tower;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCentralityTowerJet"), FilteredCentrality_Tower);
//   if(!FilteredCentrality_Tower.isValid()){
//     evValid = false;
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("FilteredCentralityTowerJet") << std::endl;
//   }

//   // 1D filtered jets
//   edm::Handle<l1slhc::L1TowerJetCollection> Filtered1D_Tower;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("Filtered1DTowerJet"), Filtered1D_Tower);
//   if(!Filtered1D_Tower.isValid()){
//     evValid = false;
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("Filtered1DTowerJet") << std::endl;
//   }

  // Pre-PU subtracted jets
  edm::Handle<l1slhc::L1TowerJetCollection> PrePUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSubTowerJet"), PrePUS_Tower);
  if(!PrePUS_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSubTowerJet") << std::endl;
  }  
  edm::Handle<l1slhc::L1TowerJetCollection> PUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUSubTowerJet"), PUS_Tower);
  if(!PUS_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUSubTowerJet") << std::endl;
  }  

  // ************************************************************
  // *                        RECO Jets                         *
  // ************************************************************
  SUBPRINT("RECO Jets")


  //PU subtracted AK5 calo jets-must be in root file read in 
  // Actually we hijacked this collection for the 
  // Generator jets
  edm::Handle<reco::GenJetCollection> PrePUSRawAk5Jets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSPreCalibCaloJets"), PrePUSRawAk5Jets);
  if(!PrePUSRawAk5Jets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSPreCalibCaloJets") << std::endl;
    evValid = false;
  }

/*
  //PU subtracted AK5 calo jets-must be in root file read in
  edm::Handle<reco::CaloJetCollection> PUSRawAk5Jets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUSPreCalibCaloJets"), PUSRawAk5Jets);
  if(!PUSRawAk5Jets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUSPreCalibCaloJets") << std::endl;
    evValid = false;
  }

  //PU subtracted AK5 calo jets-must be in root file read in
  edm::Handle<reco::CaloJetCollection> PrePUSAk5Jets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSCaloJets"), PrePUSAk5Jets);
  if(!PrePUSAk5Jets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSCaloJets") << std::endl;
    evValid = false;
  }

  //PU subtracted AK5 calo jets
  edm::Handle<reco::CaloJetCollection> PUSAk5Jets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUsubCaloJets"), PUSAk5Jets);
  if(!PUSAk5Jets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUsubCaloJets") << std::endl;
    evValid = false;
  }
*/
//   // uncorrected jets!
//   edm::Handle<edm::View< reco::CaloJet > > PrePUSAk5Jets;
//   //edm::Handle<reco::CaloJetCollection> PrePUSAk5Jets;
//   iEvent.getByLabel("ak5CaloJets", PrePUSAk5Jets );
//   if(!PrePUSAk5Jets.isValid()){
//     edm::LogWarning("MissingProduct") << "ak5CaloJets"<< std::endl;
//     evValid = false;
//   }




  if (evValid == false){

    edm::LogWarning("EventInvalid") << "Error, a collection is missing from the event.\n";
    
  }




  

   // ****************************************************************************************************
   // *                                        Primary vertices                                          *
   // ****************************************************************************************************
   TString prefix = "Vertex";

   
   // Extract number of reconstructed ak5 primary vertices
   //int NVTX    = vertices->size();
   
   //   std::cout << "NVTX = " << NVTX << "\n";
   
   
   //hist1D[prefix + "_NVTX"]->Fill(NVTX);

   /* 
   for ( reco::VertexCollection::const_iterator vertex = vertices->begin(); vertex != vertices->end(); ++vertex ) {
     

     double vtxChi2         = vertex->chi2();
     int    vtxNdof         = vertex->ndof();
     //      int    vtxTracks = vertex->tracks()->size();
     //      std::cout << "Chi2 = " << vtxChi2 << "\tNdof = " << vtxNdof << "\tChi2/ndof" << vtxChi2/vtxNdof << "\n";
     
     hist1D[prefix + "_Chi2"]          ->Fill(vtxChi2);
     hist1D[prefix + "_Ndof"]          ->Fill(vtxNdof);
     if (vtxNdof != 0) // Avoid division by zero 
       hist1D[prefix + "_Chi2_over_Ndof"]->Fill(vtxChi2/vtxNdof);
     hist2D[prefix + "_Ndof_vs_Chi2"]  ->Fill(vtxChi2, vtxNdof);
     
   }
   */




   // ****************************************************************************************************
   // *                                              Jets                                                *
   // ****************************************************************************************************



   // ********************************************************************************
   // *                                  Tower jets                                  *
   // ********************************************************************************

   // Pre PU subtraction
   //
   for (l1slhc::L1TowerJetCollection::const_iterator Tower_It = PrePUS_Tower->begin(); Tower_It != PrePUS_Tower->end(); ++Tower_It ){
     
     math::PtEtaPhiMLorentzVector tempJet;
     //<TEMP>
     //     double scaleFactor = 0.5;
       // </TEMP>

     double tempEta = Tower_It->p4().eta();
     if(foldEta) tempEta = fabs(tempEta);
     tempJet.SetCoordinates( Tower_It->p4().Pt(), tempEta, Tower_It->p4().phi(), Tower_It->p4().M() );



     // Only retain jet if it passes jet cleaning
     if ( !(cleanL1Jet(tempJet)) ) continue;
     outputPrePUSTowerJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );
     
//      // TODO APPLY CLEANING
//     TLorentzVector *tempJet2;
//     tempJet2->SetPtEtaPhiM( Tower_It->p4().Pt(), Tower_It->p4().eta(), Tower_It->p4().phi(), Tower_It->p4().M() );
//     outputPrePUSTowerJet_TLorentz->push_back( tempJet2 );


   }
   for (l1slhc::L1TowerJetCollection::const_iterator Tower_It = PUS_Tower->begin(); Tower_It != PUS_Tower->end(); ++Tower_It ){
     
     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( Tower_It->p4().Pt(), Tower_It->p4().eta(), Tower_It->p4().phi(), Tower_It->p4().M() );
     
     // Only retain jet if it passes jet cleaning
     if ( !(cleanL1Jet(tempJet)) ) continue;
     outputPUSTowerJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );

//      // TODO APPLY CLEANING
//      TLorentzVector *tempJet2;
//     tempJet2->SetPtEtaPhiM( Tower_It->p4().Pt(), Tower_It->p4().eta(), Tower_It->p4().phi(), Tower_It->p4().M() );
//     outputPUSTowerJet_TLorentz->push_back( tempJet2 );
     
   }
 

   // ********************************************************************************
   // *                                Current jets                                  *
   // ********************************************************************************


/*

   // Uncalibrated GCT jets
   std::vector <l1extra::L1JetParticle> unsortedUncalibCurrentL1Jets;



   // Tau jets
   for(unsigned int iTauJet = 0; iTauJet < L1GCTUncalibTauJet->size(); ++iTauJet ) {

     // GCT always returns four candidates, regardless of energy. Store only positive energy jets.
     if ( L1GCTUncalibTauJet->at(iTauJet).empty() ){ continue; }

     // Correct for 0.25 LSB
     double PtLSBScale = 4;
     double Pt   = PtLSBScale * L1GCTUncalibTauJet->at(iTauJet).rank();
     int RCTiEta = L1GCTUncalibTauJet->at(iTauJet).regionId().ieta();
     int RCTiPhi = L1GCTUncalibTauJet->at(iTauJet).regionId().iphi();
     
     // Calculate phi
     double phiDeg = RCTiPhi*20;
     if ( phiDeg > 160){ phiDeg -= 360; } // Get phi in degrees
     double phi = ( phiDeg*PI )/180.;     // Get phi in radians

     // Calculate eta
     int RCTEtaIndex = RCTiEta - 4;
     double eta =( RCTEtaRegions[ RCTEtaIndex ] + RCTEtaRegions[ RCTEtaIndex + 1 ] )/2;

     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( Pt, eta, phi, 0 );

     // Only retain jet if it passes jet cleaning                                                                                                           
     if ( !(cleanL1Jet(tempJet)) ) continue;
     unsortedUncalibCurrentL1Jets.push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );

     //          std::cout << "Uncalib TAU: Pt = " << Pt << "\tEta = " << eta << "\tPhi = " << phiDeg << "\t" << phi << "\n";

   }
   // Current uncalibrated L1 central jets
   for(unsigned int iCenJet = 0; iCenJet < L1GCTUncalibCenJet->size(); ++iCenJet ) {

     // GCT always returns four candidates, regardless of energy. Store only positive energy jets.
     if ( L1GCTUncalibCenJet->at(iCenJet).empty() ){ continue; }
     
     // Correct for 0.25 LSB
     double PtLSBScale = 4;
     double Pt   = PtLSBScale*L1GCTUncalibCenJet->at(iCenJet).rank();
     int RCTiEta = L1GCTUncalibCenJet->at(iCenJet).regionId().ieta();
     int RCTiPhi = L1GCTUncalibCenJet->at(iCenJet).regionId().iphi();
     
     // Calculate phi
     double phiDeg = RCTiPhi*20;
     if ( phiDeg > 160){ phiDeg -= 360; } // Get phi in degrees
     double phi = ( phiDeg*PI )/180.;     // Get phi in radians

     // Calculate eta
     int RCTEtaIndex = RCTiEta - 4;
     double eta =( RCTEtaRegions[ RCTEtaIndex ] + RCTEtaRegions[ RCTEtaIndex + 1 ] )/2;
     

     // Store jet
     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( Pt, eta, phi, 0 );

     // Only retain jet if it passes jet cleaning                                                                                                           
     if ( !(cleanL1Jet(tempJet)) ) continue;
     unsortedUncalibCurrentL1Jets.push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );

     //          std::cout << "Uncalib CEN: Pt = " << Pt << "\tEta = " << eta << "\tPhi = " << phiDeg << "\t" << phi << "\n";

   }
*/
  


   // Loop through central jet and tau jet collections
   //
   std::vector <l1extra::L1JetParticle> unsortedCurrentL1Jets;


   /*
   for (uint i = 0; i < l1extraparticles.size(); ++i){
     
     Handle<l1extra::L1JetParticleCollection> currl1Jets;
     iEvent.getByLabel(l1extraparticles[i],currl1Jets); 
     
     for(l1extra::L1JetParticleCollection::const_iterator itr = currl1Jets->begin();  itr != currl1Jets->end(); ++itr ) {
       
       math::PtEtaPhiMLorentzVector tempJet;
       tempJet.SetCoordinates( itr->p4().Pt(), itr->p4().eta(), itr->p4().phi(), itr->p4().M() );

       //        double phiInDeg = 180/PI* (itr->p4().phi()) ;
// 	      if ( fabs(itr->p4().eta() ) < 3 )
// 			std::cout << "L1 Jet Pt = " << itr->p4().Pt() << "\tEta = " << itr->p4().eta() << "\tPhi = " << phiInDeg  << "\n";

       // Only retain jet if it passes jet cleaning
       if ( !(cleanL1Jet(tempJet)) ) continue;
       //       outputCurrentJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );
       unsortedCurrentL1Jets.push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );


// 	// TODO APPLY CLEANING
// 	TLorentzVector *tempJet2;
// 	tempJet2->SetPtEtaPhiM( itr->p4().Pt(), itr->p4().eta(), itr->p4().phi(), itr->p4().M() );
// 	outputCurrentJet_TLorentz->push_back( tempJet2 );


     }
     
   }
   */
   
   //      std::cout << "--------------------------------------------------------------------------------\n";
  
   
   // Sort current jets
//   if ( unsortedCurrentL1Jets.size() > 1){
//     std::sort( unsortedCurrentL1Jets.begin(), unsortedCurrentL1Jets.end(),  CurrentL1JetRankDescending );
//   }
   // Sort uncalibrated current jets
//   if ( unsortedUncalibCurrentL1Jets.size() > 1){
//     std::sort( unsortedUncalibCurrentL1Jets.begin(), unsortedUncalibCurrentL1Jets.end(),  CurrentL1JetRankDescending );
//   }

   //int matchedUncalib(0);
   /*
   // Store sorted Current jets
   for (uint iJet = 0; iJet < unsortedCurrentL1Jets.size(); ++iJet){
       
       // Store now sorted L1Jets
       outputCurrentJet_L1Jet->push_back( unsortedCurrentL1Jets[ iJet ] );



       
       // 5 GeV seed threshold problem - Could not get the same number of jets to be retained when emulating the seed threshold
       //                                as this does not effect the clustering emulate the threshold by only storing uncalibrated
       //                                jets that can be matched to calibrated jets with the threshold applied

       // Store sorted uncalibrated Current jets
       for (uint iUncalibJet = 0; iUncalibJet < unsortedUncalibCurrentL1Jets.size(); ++iUncalibJet){
	 
	 // Check there is a corresponding calibrated jet
	 
 	double uncalibEta = unsortedUncalibCurrentL1Jets[ iUncalibJet ].p4().eta();
 	double calibEta   = unsortedCurrentL1Jets[ iJet ].p4().eta();
	double deltaEta   = uncalibEta - calibEta;

	 // Compare to 1 dp precision to avoid comparisons incorrectly not being matched
	if ( fabs(deltaEta) < 0.1 ){

	  double uncalibPhi = unsortedUncalibCurrentL1Jets[ iUncalibJet ].p4().phi();
	  double calibPhi = unsortedCurrentL1Jets[ iJet ].p4().phi();
	  double deltaPhi = uncalibPhi - calibPhi;

	  // Compare to 1 dp precision to avoid comparisons incorrectly not being matched. No need to worry about Phi-wrap around as we are only looking for
	  // identical phis
	  if ( fabs(deltaPhi) < 0.1 ){

// 	    std::cout << "MATCHED: " << "\t" << calibEta   << "\t" << uncalibEta << "\t" 
// 		      << uncalibPhi  << "\t" << uncalibPhi << "\n";
	    
	    // Store now sorted L1Jets
	    outputUncalibCurrentJet_L1Jet->push_back( unsortedUncalibCurrentL1Jets[ iUncalibJet ] );
	    matchedUncalib++;
	    break;
	  }
	}
       }

       
   }
   */

//    int calibTotal         = unsortedCurrentL1Jets.size();
//    double percentageMatch = 100.* double(matchedUncalib)/calibTotal;

//    if (percentageMatch < 100.){

//      std::cout << "Matched : " << matchedUncalib << "/" << calibTotal << "\t" 
// 	       << percentageMatch << "%\n";


//      for (uint iJet = 0; iJet < unsortedCurrentL1Jets.size(); ++iJet){
      

//        double calibEta   = unsortedCurrentL1Jets[ iJet ].p4().eta();
//        double calibPhi = unsortedCurrentL1Jets[ iJet ].p4().phi();
//        std::cout << "\t" << calibEta   <<  "\t" << calibPhi << "\n";
	     
//      }
//      std::cout << std::endl;

//      for (uint iUncalibJet = 0; iUncalibJet < unsortedUncalibCurrentL1Jets.size(); ++iUncalibJet){

//        double uncalibEta = unsortedUncalibCurrentL1Jets[ iUncalibJet ].p4().eta();
//        double uncalibPhi = unsortedUncalibCurrentL1Jets[ iUncalibJet ].p4().phi();
//        std::cout << "\t" << uncalibEta   <<  "\t" << uncalibPhi << "\n";

//      }

//    }

//    std::cout << "--------------------------------------------------------------------------------\n";
   
   // ********************************************************************************
   // *                                   RECO jets                                  *
   // ********************************************************************************

   // TEMPORARY
   TLorentzVector leadJet;
   double minDeltaR = 9999;
   double minDeltaRPT = 9999;
   double minDeltaRRelPT = 9999;
   //int jetIndex = 0;
   // TEMPORARY

   //   std::cout << "--------------------------------------------------------------------------------\n";






   // PrePUS raw jet collection
   //
   for (reco::GenJetCollection::const_iterator PrePUSRawak5_It = PrePUSRawAk5Jets->begin(); PrePUSRawak5_It != PrePUSRawAk5Jets->end(); ++PrePUSRawak5_It ){

     math::PtEtaPhiMLorentzVector tempJet;
     double tempEta=PrePUSRawak5_It->p4().eta();
     if(foldEta) tempEta=fabs(tempEta);
     tempJet.SetCoordinates( PrePUSRawak5_It->p4().Pt(), tempEta, PrePUSRawak5_It->p4().phi(), PrePUSRawak5_It->p4().M() );
     
     // Only retain jet if it passes jet cleaning
     if ( !(cleanRECOJet(tempJet)) ) continue;
     outputPrePUSRawAk5CaloJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );
      
   }

   // PUS raw jet collection
   //
   /*
   for (reco::CaloJetCollection::const_iterator PUSRawak5_It = PUSRawAk5Jets->begin(); PUSRawak5_It != PUSRawAk5Jets->end(); ++PUSRawak5_It ){

     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( PUSRawak5_It->p4().Pt(), PUSRawak5_It->p4().eta(), PUSRawak5_It->p4().phi(), PUSRawak5_It->p4().M() );
     
     // Only retain jet if it passes jet cleaning
     if ( !(cleanRECOJet(tempJet)) ) continue;
     outputPUSRawAk5CaloJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );
      
   }
*/
   //L2L3 PrePUS jet collection
   //



   //   reco::CaloJetCollection  unsortedPrePUSCaloJets;

/*
   for (reco::CaloJetCollection::const_iterator PrePUSak5_It = PrePUSAk5Jets->begin(); PrePUSak5_It != PrePUSAk5Jets->end(); ++PrePUSak5_It ){

     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( PrePUSak5_It->p4().Pt(), PrePUSak5_It->p4().eta(), PrePUSak5_It->p4().phi(), PrePUSak5_It->p4().M() );
     
     // Only retain jet if it passes jet cleaning
     if ( !(cleanRECOJet(tempJet)) ) continue;




          outputPrePUSAk5CaloJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );
	  //     // TEMPORARY (Sort jets) - Calibrated jets previously weren't sorted after calibration. Use this until new ntuples are produced.
	  //unsortedPrePUSCaloJets.push_back( *PrePUSak5_It );
      
   }
*/
//    // <TEMP>
//    std::sort( unsortedPrePUSCaloJets.begin(), unsortedPrePUSCaloJets.end(), caloJetRankDescending );
//    // Store the now sorted jet collection
//    for(reco::CaloJetCollection::const_iterator PrePUSak5_It = unsortedPrePUSCaloJets.begin(); PrePUSak5_It != unsortedPrePUSCaloJets.end() ; ++PrePUSak5_It) {

//      math::PtEtaPhiMLorentzVector tempJet;
//      tempJet.SetCoordinates( PrePUSak5_It->p4().Pt(), PrePUSak5_It->p4().eta(), PrePUSak5_It->p4().phi(), PrePUSak5_It->p4().M() );

//      outputPrePUSAk5CaloJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );
//    }
//    // </TEMP>


 //     outputPrePUSAk5CaloJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );


   //L2L3 PU corrected PUS jet collection
   //

   // Sort jets - Calibrated jets previously weren't sorted after calibration. Use this until new ntuples are produced.
   //   std::sort( PUSAk5Jets.begin(), PUSAk5Jets.end(), caloJetRankDescending );

//    double oldMax = 99999;
//    int count = 0;
/*
   for (reco::CaloJetCollection::const_iterator PUSak5_It = PUSAk5Jets->begin(); PUSak5_It != PUSAk5Jets->end(); ++PUSak5_It ){



//      if ( oldMax < PUSak5_It->p4().Pt()){
//        std::cout << count << " " << oldMax << " "<< PUSak5_It->p4().Pt() << " " << std::endl;
//      }
//      oldMax = PUSak5_It->p4().Pt();
//      count++;
     

//      double energyInHF = PUSak5_It->hadEnergyInHF() + PUSak5_It->emEnergyInHF();

//      std::cout << energyInHF << "\t" << PUSak5_It->p4().Pt() << "\t" << PUSak5_It->p4().eta() << "\n";


     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( PUSak5_It->p4().Pt(), PUSak5_It->p4().eta(), PUSak5_It->p4().phi(), PUSak5_It->p4().M() );
     
     //TEMPORARY
     if (jetIndex == 0){
       leadJet.SetPtEtaPhiM( PUSak5_It->p4().Pt(), PUSak5_It->p4().eta(), PUSak5_It->p4().phi(), PUSak5_It->p4().M() );
     }
     else{
       TLorentzVector  temporaryJet;
       temporaryJet.SetPtEtaPhiM( PUSak5_It->p4().Pt(), PUSak5_It->p4().eta(), PUSak5_It->p4().phi(), PUSak5_It->p4().M() );

       double deltaR =  leadJet.DeltaR( temporaryJet );
       if (deltaR < minDeltaR){
	 minDeltaR      = deltaR;
	 minDeltaRPT    = temporaryJet.Pt();
	 minDeltaRRelPT = ( temporaryJet.Pt() - leadJet.Pt() )/leadJet.Pt();



       }


       hist1D["ak5_leadJet_JetDeltaR"]->Fill( deltaR );
       
     }
     jetIndex++;
     //TEMPORARY


     // Only retain jet if it passes jet cleaning
     if ( !(cleanRECOJet(tempJet)) ) continue;
     outputAk5CaloJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );


//      // TODO APPLY CLEANING
//      TLorentzVector *tempJet2;
//      tempJet2->SetPtEtaPhiM( PUSak5_It->p4().Pt(), PUSak5_It->p4().eta(), PUSak5_It->p4().phi(), PUSak5_It->p4().M() );
//      outputAk5CaloJet_TLorentz->push_back( tempJet2 );
      
   }
*/

   //TEMPORARY
   if (minDeltaR != 9999){
     hist1D["ak5_leadJet_minDeltaR"]->Fill(minDeltaR);
     hist2D["ak5_leadJet_minDeltaRPT"]->Fill( minDeltaR, minDeltaRPT);
     hist2D["ak5_leadJet_minDeltaRRelPT"]->Fill( minDeltaR, minDeltaRRelPT );

   }
   //TEMPORARY


   // PrePUS ak5
   //    for (reco::CaloJetCollection::const_iterator PrePUSak5_It = PrePUSAk5Jets->begin(); PrePUSak5_It != PrePUSAk5Jets->end(); ++PrePUSak5_It ){
//    for ( edm::View< reco::CaloJet >::const_iterator PrePUSak5_It = PrePUSAk5Jets->begin(); PrePUSak5_It != PrePUSAk5Jets->end(); ++PrePUSak5_It ){
     
//      math::PtEtaPhiMLorentzVector tempJet;
//      tempJet.SetCoordinates( PrePUSak5_It->p4().Pt(), PrePUSak5_It->p4().eta(), PrePUSak5_It->p4().phi(), PrePUSak5_It->p4().M() );
     
//      // Only retain jet if it passes jet cleaning
//      if ( !(cleanRECOJet(tempJet)) ) continue;
//      outputAk5CaloUncorrJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );
     
//    }






   iEvent.put( outputNVTX, "NVTX" );
//    iEvent.put( outputprePusHT_Tower , "prePusHT_Tower" );
//    iEvent.put( outputHT_Tower,  "HT_Tower" );
//    iEvent.put( outputMHT_Tower, "MHT_Tower" );
//    iEvent.put( outputPUS_Tower_TLorentz, "PUS_Tower" );


   // Jets
   iEvent.put( outputPrePUSTowerJet_L1Jet,     "PrePUSTowerJetL1Jet"); // -> THESE will be out L1 Jets
   iEvent.put( outputPUSTowerJet_L1Jet,     "PUSTowerJetL1Jet"); // -> THESE will be out L1 Jets PU
   iEvent.put( outputUncalibCurrentJet_L1Jet,  "CurrentUncalibJetL1Jet");
   iEvent.put( outputCurrentJet_L1Jet,         "CurrentJetL1Jet");
   iEvent.put( outputPrePUSRawAk5CaloJet_L1Jet,"PrePUSRawAk5CaloJetL1Jet"); // -> These are the GenJet Collection
   //iEvent.put( outputPUSRawAk5CaloJet_L1Jet,   "PUSRawAk5CaloJetL1Jet");
   //iEvent.put( outputPrePUSAk5CaloJet_L1Jet,   "PrePUSAk5CaloJetL1Jet");
   //   iEvent.put( outputAk5CaloUncorrJet_L1Jet, "Ak5CaloUncorrJetL1Jet");
    

}

// ------------ method called once each job just before starting event loop  ------------
void 
JetProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
JetProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}





// L1 Jet cleaning
bool 
JetProducer::cleanL1Jet(math::PtEtaPhiMLorentzVector jet){
  
  // ****************************************
  // *           Jet quality cuts           *
  // ****************************************
  if( jet.Pt() < minL1JetPt ) return false;
  if( fabs(jet.Eta()) > maxL1JetEta ) return false;
    
  return true;

}

// RECO Jet cleaning
bool 
JetProducer::cleanRECOJet(math::PtEtaPhiMLorentzVector jet){
  
  // ****************************************
  // *           Jet quality cuts           *
  // ****************************************
  if( jet.Pt() < minRECOJetPt ) return false;
  if( fabs(jet.Eta()) > maxRECOJetEta ) return false;
    







     // LOOSE JET SELECTION 
     // -------------------

     // N90 hits > 1
//      if (!(PUSak5_It->n90 > 1)) continue;
//      if ( !(PUSak5_It->EnergyFractionHadronic < 0.99) ) continue;

//        if ( type == typeid(reco::PFJet) ) {
// 	 const reco::PFJet pfjet = static_cast<const reco::PFJet &>((*jets)[iJet]);
// 	 ThisIsClean=true;
// 	 //apply following only if |eta|<2.4: CHF>0, CEMF<0.99, chargedMultiplicity>0   
// 	 if(( pfjet.chargedHadronEnergy()/ pfjet.energy())<= 0.0  
// 	    && fabs(pfjet.eta())<2.4) ThisIsClean=false; 
// 	 if( (pfjet.chargedEmEnergy()/pfjet.energy())>= 0.99 
// 	     && fabs(pfjet.eta())<2.4 ) ThisIsClean=false;
// 	 if( pfjet.chargedMultiplicity()<=0 && fabs(pfjet.eta())<2.4 ) 
// 	   ThisIsClean=false;


// 	 // always require #Constituents > 1
// 	 if( pfjet.nConstituents() <=1 ) ThisIsClean=false;


// 	 if(ThisIsClean && 
// 	    (pfjet.neutralHadronEnergy()/pfjet.energy())< 0.99 
// 	    && (pfjet.neutralEmEnergy()/pfjet.energy())<0.99) 
// 	   passedId=true;
//        }







  return true;

}



//define this as a plug-in
DEFINE_FWK_MODULE(JetProducer);
