// -*- C++ -*-
//
// Package:    CaloJetFootprint
// Class:      CaloJetFootprint
// 
/**\class CaloJetFootprint CaloJetFootprint.cc AnalyseUpgradeJets/src/CaloJetFootprint.cc

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
//#include "RecoJets/CaloJetFootprints/interface/JetIDHelper.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "SimDataFormats/SLHC/interface/EtaPhiContainer.h"
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/TriggerTowerGeometry.h"





#include "AnalyseUpgradeJets/AnalyzeJets/interface/printing.h"
// #include "AnalyseUpgradeJets/AnalyzeJets/interface/JetMatch.h"
// #include "AnalyseUpgradeJets/AnalyzeJets/interface/helperFunctions.cc"

#include <algorithm>  // for sorting




//
// class declaration
//

class CaloJetFootprint : public edm::EDProducer {
   public:
      explicit CaloJetFootprint(const edm::ParameterSet&);
      ~CaloJetFootprint();

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
CaloJetFootprint::CaloJetFootprint(const edm::ParameterSet& iConfig): conf_(iConfig)
{








  // ****************************************************************************************************
  // *                                 Load configuration parameters                                    *
  // ****************************************************************************************************

  // L1 jet cleaning parameters
  minL1JetPt  = iConfig.getParameter<double> ("minL1JetPt");
  maxL1JetEta = iConfig.getParameter<double> ("maxL1JetEta");
  // RECO jet cleaning parameters
  minRECOJetPt  = iConfig.getParameter<double> ("minRECOJetPt");
  maxRECOJetEta = iConfig.getParameter<double> ("maxRECOJetEta");









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


CaloJetFootprint::~CaloJetFootprint()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CaloJetFootprint::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   bool evValid = true;




   std::auto_ptr<int>                      outputNVTX( new int() );

   // RECO
   std::auto_ptr<l1extra::L1JetParticleCollection> outputPrePUSRawAk5CaloJet_L1Jet( new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputPUSRawAk5CaloJet_L1Jet(    new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputPrePUSAk5CaloJet_L1Jet(    new l1extra::L1JetParticleCollection() );
   std::auto_ptr<l1extra::L1JetParticleCollection> outputAk5CaloJet_L1Jet(          new l1extra::L1JetParticleCollection() );


   // ****************************************************************************************************
   // *                                             Handles                                              *
   // ****************************************************************************************************
      
   //Need this for information about PU
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RecoVertices"), vertices); 
   if(!vertices.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RecoVertices") << std::endl;
     evValid = false;
   }



   // ********************************************************************************
   // *                                    Jets                                      *
   // ********************************************************************************





  // ************************************************************
  // *                        RECO Jets                         *
  // ************************************************************
  SUBPRINT("RECO Jets")


  //PU subtracted AK5 calo jets-must be in root file read in
  edm::Handle<reco::CaloJetCollection> PrePUSRawAk5Jets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSPreCalibCaloJets"), PrePUSRawAk5Jets);
  if(!PrePUSRawAk5Jets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSPreCalibCaloJets") << std::endl;
    evValid = false;
  }

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



  if (evValid == false){

    edm::LogWarning("EventInvalid") << "Error, a collection is missing from the event.\n";
    
  }




  

   // ****************************************************************************************************
   // *                                        Primary vertices                                          *
   // ****************************************************************************************************
   TString prefix = "Vertex";

   
   // Extract number of reconstructed ak5 primary vertices
   int NVTX    = vertices->size();
   
   //   std::cout << "NVTX = " << NVTX << "\n";
   
   
   hist1D[prefix + "_NVTX"]->Fill(NVTX);

   
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



   // ********************************************************************************
   // *                                Cleaning cuts                                 *
   // ********************************************************************************
   //
   //
   //TODO: Add cleaning cuts
   //
   //
   *outputNVTX = NVTX;




   
   // ********************************************************************************
   // *                                   RECO jets                                  *
   // ********************************************************************************



   for (reco::CaloJetCollection::const_iterator PUSak5_It = PUSAk5Jets->begin(); PUSak5_It != PUSAk5Jets->end(); ++PUSak5_It ){


     reco::CaloJet jet = (*PUSak5_It);

     std::cout << jet.eta() << "\t" << jet.getCaloConstituents().size() <<  "\n";
     for(unsigned int i = 0; i < jet.getCaloConstituents().size(); ++i){
     std::cout << i << "\n";

       std::cout << jet.getCaloConstituent(i)->eta()-jet.eta() << "\n";
       //       histContainer_["dEta"]->Fill( jet.getCaloConstituent(i)->eta()-jet.eta() );

     }



//      math::PtEtaPhiMLorentzVector tempJet;
//      tempJet.SetCoordinates( PUSak5_It->p4().Pt(), PUSak5_It->p4().eta(), PUSak5_It->p4().phi(), PUSak5_It->p4().M() );
     

//      // Only retain jet if it passes jet cleaning
//      if ( !(cleanRECOJet(tempJet)) ) continue;
//      outputAk5CaloJet_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );
      
   }



}

// ------------ method called once each job just before starting event loop  ------------
void 
CaloJetFootprint::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CaloJetFootprint::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
CaloJetFootprint::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CaloJetFootprint::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CaloJetFootprint::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CaloJetFootprint::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CaloJetFootprint::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}





// L1 Jet cleaning
bool 
CaloJetFootprint::cleanL1Jet(math::PtEtaPhiMLorentzVector jet){
  
  // ****************************************
  // *           Jet quality cuts           *
  // ****************************************
  if( jet.Pt() < minL1JetPt ) return false;
  if( fabs(jet.Eta()) > maxL1JetEta ) return false;
    
  return true;

}

// RECO Jet cleaning
bool 
CaloJetFootprint::cleanRECOJet(math::PtEtaPhiMLorentzVector jet){
  
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
DEFINE_FWK_MODULE(CaloJetFootprint);
