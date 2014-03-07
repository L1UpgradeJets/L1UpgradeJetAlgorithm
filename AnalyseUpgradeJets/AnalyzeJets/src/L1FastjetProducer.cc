// Original Author:  Mark Baber Imperial College, London

// Runs the fastjet algorithm on trigger towers to obtain L1-resolution jets using offline algorithms

#include "TH1F.h"
#include "TH2F.h"

#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"


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


// system include files                                                                                                                                       
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include <iostream>
#include <fstream>


#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"


#include "FWCore/ParameterSet/interface/FileInPath.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"


#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/TriggerTowerGeometry.h"
#include "SimDataFormats/SLHC/interface/EtaPhiContainer.h"


#include <algorithm> 
//#include <string> 
#include <vector>

//<fastjet>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <fastjet/tools/Filter.hh>


// #include "AnalyseUpgradeJets/fastjetInstall/include/fastjet/JetDefinition.hh"
// #include "AnalyseUpgradeJets/fastjetInstall/include/fastjet/PseudoJet.hh"
// #include "AnalyseUpgradeJets/fastjetInstall/include/fastjet/ClusterSequence.hh"
// #include "AnalyseUpgradeJets/fastjetInstall/include/fastjet/GhostedAreaSpec.hh"
// #include "AnalyseUpgradeJets/fastjetInstall/include/fastjet/ClusterSequenceArea.hh"
//</fastjet>

const double PI = 3.141592654;

using namespace fastjet;
 

class L1FastjetProducer : public edm::EDProducer {
  public:
        explicit L1FastjetProducer(const edm::ParameterSet&);
        ~L1FastjetProducer();
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

        // ----------member data ---------------------------
        edm::ParameterSet conf_;
        edm::Service<TFileService> fs;
  
        // Tower geometry converter 
        TriggerTowerGeometry mTowerGeo;

  // Histogram containers                                                                                                                                 
  std::map< TString, TH1*> hist1D;
  std::map< TString, TH2*> hist2D;

  int eventsProcessed;
  int maxEvents;


};

L1FastjetProducer::L1FastjetProducer( const edm::ParameterSet& iConfig):  conf_(iConfig){

//   produces<reco::CaloJetCollection>( "TTAk5CaloJet" );
//   produces<reco::CaloJetCollection>( "TTKt6CaloJet" );


  produces <l1extra::L1JetParticleCollection>("TTAk5L1Jet");
  produces <l1extra::L1JetParticleCollection>("TTKt6L1Jet");
  produces <l1extra::L1JetParticleCollection>("CaloAk5L1Jet");
  produces <l1extra::L1JetParticleCollection>("CaloKt6L1Jet");



//   TFileDirectory fastjetDir = fs->mkdir( "Fastjet" );


//   hist2D["Ak5TTJets_deltaPhi_vs_deltaEta"]  = fastjetDir.make<TH2D>("Ak5TTJets_deltaPhi_vs_deltaEta", 
// 								    "Constituent muliticity #Delta#phi vs #Delta#eta;#Delta#eta;#Delta#phi",
// 								    96, -0.59375, 0.59375, 96, -0.59375, 0.59375);

//   hist2D["Ak5TTJets_deltaPhi_vs_deltaEta_PtWeighted"]  = fastjetDir.make<TH2D>("Ak5TTJets_deltaPhi_vs_deltaEta_PtWeighted", 
// 								    "Constituent energy #Delta#phi vs #Delta#eta;#Delta#eta;#Delta#phi",
// 								    96, -0.59375, 0.59375, 96, -0.59375, 0.59375);



//   eventsProcessed = 0;
//   maxEvents = 100;

}

L1FastjetProducer::~L1FastjetProducer(  )
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
L1FastjetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){




  std::auto_ptr<l1extra::L1JetParticleCollection> outputTTAk5_L1Jet(   new l1extra::L1JetParticleCollection() );
  std::auto_ptr<l1extra::L1JetParticleCollection> outputTTKt6_L1Jet(   new l1extra::L1JetParticleCollection() );
  std::auto_ptr<l1extra::L1JetParticleCollection> outputCaloAk5_L1Jet( new l1extra::L1JetParticleCollection() );
  std::auto_ptr<l1extra::L1JetParticleCollection> outputCaloKt6_L1Jet( new l1extra::L1JetParticleCollection() );



  // ********************************************************************************
  // *                                   Towers                                     *
  // ********************************************************************************

  // TT collection
  //  SUBPRINT("Trigger towers")
    edm::Handle<l1slhc::L1CaloTowerCollection> caloTowers;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalorimeterTowers"), caloTowers);
  if(!caloTowers.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalorimeterTowers") << std::endl;
    //    evValid = false;
  }


  //  bool caloTowersPresent = true;
  // CaloTowers
  edm::Handle<CaloTowerCollection> caloFineTowers;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalorimeterFineTowers"), caloFineTowers);
  if(!caloFineTowers.isValid()){
    //    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalorimeterFineTowers") << std::endl;
    //    caloTowersPresent = false;
  }


  std::vector<PseudoJet> particlesTT;

  for( l1slhc::L1CaloTowerCollection::const_iterator lTT_It = caloTowers->begin() ;
       lTT_It != caloTowers->end() ; ++lTT_It ){

    // ****************************************
    // *    Load the calorimeter tower data   *
    // ****************************************
    int E      = lTT_It->E();
    int H      = lTT_It->H();
    int EplusH = E + H;
    int iEta   = lTT_It->iEta();
    int iPhi   = lTT_It->iPhi();
     double Eta = mTowerGeo.eta(iEta);
     double Phi = mTowerGeo.phi(iPhi);

    // Restrict to central TTs
    if (abs(iEta) > 28)
      continue;

    //<BUG?????>
    if (Phi > PI)
      Phi -= 2*PI;
    //<BUG?????>

  double Px = EplusH*cos(Phi);
  double Py = EplusH*sin(Phi);
  double Pz = EplusH*sinh(Eta);
  double Eproj  = EplusH*cosh(Eta); // Projection of energy along direction vector

  // Store TT as a particle: px py pz E
  particlesTT.push_back( PseudoJet( Px, Py, Pz, Eproj) );



  //    std::cout << iEta << "\t" << iPhi << "\t" << EplusH << "\n";

  }


  std::vector<PseudoJet> particlesCalo;

  for (CaloTowerCollection::const_iterator iCalo = caloFineTowers->begin(); iCalo != caloFineTowers->end(); ++iCalo) {

    // Restrict to central TTs                                                                                                                                
    if ( iCalo->ietaAbs() > 28 ){
      continue;
    }

    // **************************************** 
    // *    Load the calorimeter tower data   * 
    // ****************************************                                                                                                               
//     double E      = iCalo->emEt();
//     double H      = iCalo->hadEt();
    double EplusH = iCalo->et();
    int iEta   = iCalo->ieta();
    int iPhi   = iCalo->iphi();
    double Eta = mTowerGeo.eta(iEta);
    double Phi = mTowerGeo.phi(iPhi);

    if (Phi > PI)
      Phi -= 2*PI;

    double Px = EplusH*cos(Phi);
    double Py = EplusH*sin(Phi);
    double Pz = EplusH*sinh(Eta);
    double Eproj  = EplusH*cosh(Eta); // Projection of energy along direction vector
    
    // Store caloTower as a particle: px py pz E
    particlesCalo.push_back( PseudoJet( Px, Py, Pz, Eproj) );


    //    std::cout << iEta << "\t" << iPhi << "\t" << Eta << "\t" << Phi << "\t" << EplusH << "\n";

  }



  // *****************************************************
  // *               Define jet algorithms               *
  // *****************************************************

  // choose jet definitions
  double Ak5R = 0.5;
  JetDefinition Ak5(antikt_algorithm, Ak5R);
  double Kt6R = 0.6;
  JetDefinition Kt6(kt_algorithm, Kt6R);

  
  // run the clustering, extract the jets
  
  // ******************************
  // *             TT             *
  // ******************************
  ClusterSequence csAk5TT(particlesTT, Ak5 );
  ClusterSequence csKt6TT(particlesTT, Kt6 );
  std::vector<PseudoJet> Ak5TTJets = sorted_by_pt(csAk5TT.inclusive_jets());
  std::vector<PseudoJet> Kt6TTJets = sorted_by_pt(csKt6TT.inclusive_jets());

  // ******************************
  // *          CaloTower         *
  // ******************************
  ClusterSequence csAk5Calo(particlesCalo, Ak5 );
  ClusterSequence csKt6Calo(particlesCalo, Kt6 );
  std::vector<PseudoJet> Ak5CaloJets = sorted_by_pt(csAk5Calo.inclusive_jets());
  std::vector<PseudoJet> Kt6CaloJets = sorted_by_pt(csKt6Calo.inclusive_jets());




  // ********************************************************************
  // *                  Book Single event distributions                 *
  // ********************************************************************



//   int Eventnr   = iEvent.id().event();
//   TString eventStr = Form("%d", Eventnr ); 
  
//   // Create an event directory 
//   TFileDirectory fastjetDir = fs->mkdir( "Fastjet" );
//   TFileDirectory eventDir    = fastjetDir.mkdir( "Event" ); 

  
//   if (eventsProcessed < maxEvents){

//     TFileDirectory eventSubDir = eventDir.mkdir( eventStr.Data() ); 
//     hist2D[eventStr + "_Ak5TTJets_deltaPhi_vs_deltaEta"]             = eventSubDir.make<TH2D>("Event " + eventStr + "_Ak5TTJets_deltaPhi_vs_deltaEta", 
// 											      "Constituent muliticity #Delta#phi vs #Delta#eta;#Delta#eta;#Delta#phi",
// 											      96, -0.59375, 0.59375, 96, -0.59375, 0.59375);
    
//     hist2D[eventStr + "_Ak5TTJets_deltaPhi_vs_deltaEta_PtWeighted"]  = eventSubDir.make<TH2D>("Event " + eventStr + "_Ak5TTJets_deltaPhi_vs_deltaEta_PtWeighted", 
// 											      "Constituent energy #Delta#phi vs #Delta#eta;#Delta#eta;#Delta#phi",
// 											      96, -0.59375, 0.59375, 96, -0.59375, 0.59375);
//   }



//   // Analyse jet footprint
//   for (unsigned i = 0; i < Ak5TTJets.size(); i++) {


//     // RESTRICT TO LEAD JET FOR NOW
//     if (i != 0 ){ continue; }

//     //     double Pt  = Ak5TTJets[i].perp();
//      double Eta = Ak5TTJets[i].rap();
//      double Phi = Ak5TTJets[i].phi();
// //      double M   = 0;

//      std::vector<PseudoJet> constituents = Ak5TTJets[i].constituents();   

//      for (unsigned j = 0; j < constituents.size(); j++) {


//        double constEta = constituents[j].rap();
//        double constPhi = constituents[j].phi();
//        double constPt  = constituents[j].perp();
       
//        double deltaEta = constEta - Eta;
//        double deltaPhi = constPhi - Phi;


//        hist2D["Ak5TTJets_deltaPhi_vs_deltaEta"]           ->Fill( deltaEta, deltaPhi );
//        hist2D["Ak5TTJets_deltaPhi_vs_deltaEta_PtWeighted"]->Fill( deltaEta, deltaPhi, constPt );

//        // Fill single event plots
//        if ( eventsProcessed < maxEvents){
// 	 hist2D[eventStr + "_Ak5TTJets_deltaPhi_vs_deltaEta"]           ->Fill( deltaEta, deltaPhi );
// 	 hist2D[eventStr + "_Ak5TTJets_deltaPhi_vs_deltaEta_PtWeighted"]->Fill( deltaEta, deltaPhi, constPt );
// 	 eventsProcessed++;
//        }

//      }


//    }





  // ****************************************************************************************************
  // *                                            Store jets                                            *
  // ****************************************************************************************************


  // Ak5TT
   for (unsigned i = 0; i < Ak5TTJets.size(); i++) {

     double Pt  = Ak5TTJets[i].perp();
     double Eta = Ak5TTJets[i].rap();
     double Phi = Ak5TTJets[i].phi();
     double M   = 0;

     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( Pt, Eta, Phi, M );
     outputTTAk5_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );


   }

  // Kt6TT
   for (unsigned i = 0; i < Kt6TTJets.size(); i++) {

     double Pt  = Kt6TTJets[i].perp();
     double Eta = Kt6TTJets[i].rap();
     double Phi = Kt6TTJets[i].phi();
     double M   = 0;

     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( Pt, Eta, Phi, M );
     outputTTKt6_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );


   }


  // Ak5Calo
   for (unsigned i = 0; i < Ak5CaloJets.size(); i++) {

     double Pt  = Ak5CaloJets[i].perp();
     double Eta = Ak5CaloJets[i].rap();
     double Phi = Ak5CaloJets[i].phi();
     double M   = 0;

     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( Pt, Eta, Phi, M );
     outputCaloAk5_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );


   }

  // Kt6Calo
   for (unsigned i = 0; i < Kt6CaloJets.size(); i++) {

     double Pt  = Kt6CaloJets[i].perp();
     double Eta = Kt6CaloJets[i].rap();
     double Phi = Kt6CaloJets[i].phi();
     double M   = 0;

     math::PtEtaPhiMLorentzVector tempJet;
     tempJet.SetCoordinates( Pt, Eta, Phi, M );
     outputCaloKt6_L1Jet->push_back( l1extra::L1JetParticle( tempJet, l1extra::L1JetParticle::JetType::kCentral, 0 ) );


   }



//     for ( int lEta = lEtaMin; lEta <= lEtaMax; ++lEta ){
//       if (lEta == 0){ continue; }
      
//       for ( int lPhi = lPhiMin; lPhi <= lPhiMax; ++lPhi ){
	

// 	std::cout << lEta << "  " << lPhi << std::endl;


//       }
//     }



//   int lTowerIndex                     = mCaloTriggerSetup->getBin( aEta, aPhi );
//   std::pair < int, int > lTowerEtaPhi = mCaloTriggerSetup->getTowerEtaPhi( lTowerIndex );
  

//   // Construct a TowerJet object at the current iEta, iPhi position with given jet shape and size
//   l1slhc::L1TowerJet lJet( mJetDiameter, mJetShape , mJetShapeMap , lTowerEtaPhi.first , lTowerEtaPhi.second  );
//   // Calculate the geometric center of the jet
//   lJet.calculateJetCenter();




//     // Store jet in the output collection
//     mOutputCollection->insert( lTowerEtaPhi.first, lTowerEtaPhi.second, lJet );






  // Jets
  iEvent.put( outputTTAk5_L1Jet,   "TTAk5L1Jet");
  iEvent.put( outputTTKt6_L1Jet,   "TTKt6L1Jet");
  iEvent.put( outputCaloAk5_L1Jet, "CaloAk5L1Jet");
  iEvent.put( outputCaloKt6_L1Jet, "CaloKt6L1Jet");



}


// ------------ method called once each job just before starting event loop  ------------
void
L1FastjetProducer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1FastjetProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
L1FastjetProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
L1FastjetProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
L1FastjetProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
L1FastjetProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1FastjetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_EDM_PLUGIN( edm::MakerPluginFactory, edm::WorkerMaker < L1FastjetProducer >, "L1FastjetProducer" );
DEFINE_FWK_PSET_DESC_FILLER( L1FastjetProducer );
