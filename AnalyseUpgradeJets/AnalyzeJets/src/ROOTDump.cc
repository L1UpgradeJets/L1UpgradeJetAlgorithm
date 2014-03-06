#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
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
#include "TLorentzVector.h"
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

using namespace l1slhc;
using namespace edm;
using namespace std;
using namespace reco;

const Double_t PI = 3.141592654;


// Switch for TT-level info
#define TT_INFO


// ========================================
// class declaration
// ========================================

class ROOTDump : public edm::EDAnalyzer {

public:
  explicit ROOTDump(const edm::ParameterSet&);
  ~ROOTDump();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);




  // ----------member data ---------------------------
  ParameterSet conf_;
  edm::InputTag m_jetCollection;
  edm::InputTag m_jetIDMap;


  // TT variables
  int TT_Et, TT_Eta, TT_Phi;
  std::vector <int> TT_Et_Arr, TT_Eta_Arr, TT_Phi_Arr;
  int TT_index; // TT index in the 1D vector (start top left tower (eta,phi) = (-28,1) scanning lines of constant phi left to right
  //                            to bottom right TT   (eta,phi) = (28,72)

  // Should reserve space in the CONSTRUCTOR!!!!
  std::vector<int> TT_ArrayNew,muonPx;


  // ************************************************************************************
  // *                                       ROOT                                       *
  // ************************************************************************************


  // Create TTrees
  // ******************************

  // TowerJets
  TTree *ProtoJetTree;
  TTree *Filtered1DJetTree;
  TTree *PrePUSJetTree;
  TTree *PUSJetTree;
  TTree *LPUSJetTree;
  
  // TT
  TTree *ttTree;

  // Event
  TTree *evTree;
  

  // Jet tree branches
  // ****************************************
  int NJets;
  std::vector<double> jetsEt;
  //  std::vector<double> jetsEtPUSub;
  std::vector<int>    jetsIPhi;
  std::vector<int>    jetsIEta;
  std::vector<double> jetsPhiCentre;
  std::vector<double> jetsEtaCentre; 
  std::vector<double> jetsPhiWeighted;
  std::vector<double> jetsEtaWeighted;
  std::vector<int>    jetsAsymPhi;
  std::vector<int>    jetsAsymEta;
  std::vector<double> jetsRealArea;



  // TT tree branches
  // ****************************************

  std::vector<int>      ttE;
  std::vector<int>      ttH;
  std::vector<int>      ttiEta;
  std::vector<int>      ttiPhi;
  std::vector<int>      ttEcalFG;
  std::vector<int>      ttHcalFG;


  // Event tree branches
  // ****************************************

  int NVertices;
  int Runnr;
  int LumiBlock;
  int Eventnr;




  // Store towerjet collection
  void fillTowerJetROOTTree(edm::Handle<l1slhc::L1TowerJetCollection> const& TowerJet, TTree *t);

};





ROOTDump::ROOTDump(const edm::ParameterSet& iConfig): conf_(iConfig){




  // Allocate enough space to store all the TT in the barrel and endcaps
  TT_ArrayNew.resize(4032);


  //now do what ever initialization is needed
  edm::Service<TFileService> tfs;



  // Tower Jet trees
  // ****************************************


  ProtoJetTree = tfs->make<TTree>("ProtoJet","Jet");

  ProtoJetTree->Branch( "NJets",       &NJets);//, "Jet multiplicity/I" );
  ProtoJetTree->Branch( "Pt",          &jetsEt);//,"Jet p_{T}/D" );
  ProtoJetTree->Branch( "iPhi",        &jetsIPhi);//,"Jet i#Phi/I" );
  ProtoJetTree->Branch( "iEta",        &jetsIEta);//,"Jet i#eta/I" );
  ProtoJetTree->Branch( "PhiCentre",   &jetsPhiCentre);//,"Jet geometric centre #phi/D" );
  ProtoJetTree->Branch( "EtaCentre",   &jetsEtaCentre);//,"Jet geometric centre #eta/D" );
  ProtoJetTree->Branch( "PhiWeighted", &jetsPhiWeighted);//,"Jet weighted #phi/D" );
  ProtoJetTree->Branch( "EtaWeighted", &jetsEtaWeighted);//,"Jet weighted #eta/D" );
  ProtoJetTree->Branch( "AsymPhi",     &jetsAsymPhi);//,    "Jet A_{#phi}/I");
  ProtoJetTree->Branch( "AsymEta",     &jetsAsymEta);//,    "Jet A_{#eta}/I");
  ProtoJetTree->Branch( "JetRealArea", &jetsRealArea);//,    "Jet Real Area/D");

  
  // Copy the tree structure
  Filtered1DJetTree = ProtoJetTree->CloneTree();
  Filtered1DJetTree->SetName("1DFiltJet");
  PrePUSJetTree = ProtoJetTree->CloneTree();
  PrePUSJetTree->SetName("PrePUSJet");
  PUSJetTree = ProtoJetTree->CloneTree();
  PUSJetTree->SetName("PUSJet");
  LPUSJetTree = ProtoJetTree->CloneTree();
  LPUSJetTree->SetName("LPUSJet");
  

  // TT tree
  // ****************************************

  ttTree = tfs->make<TTree>("TT","TT");

  ttTree->Branch( "E",      &ttE);
  ttTree->Branch( "H",      &ttH);
  ttTree->Branch( "iEta",   &ttiEta);
  ttTree->Branch( "iPhi",   &ttiPhi);
  ttTree->Branch( "EcalFG", &ttEcalFG);
  ttTree->Branch( "HcalFG", &ttHcalFG);

  // Event tree
  // ****************************************

  evTree = tfs->make<TTree>("event","event");

  evTree->Branch( "NVertices", &NVertices, "Number of reconstructed AK5 Vertices/I" );
  evTree->Branch( "Runnr",         &Runnr, "Run Number/I" );
  evTree->Branch( "LumiBlock", &LumiBlock, "Luminosity Block/I" );
  evTree->Branch( "Eventnr",     &Eventnr, "Event Number/I" );

}



ROOTDump::~ROOTDump()
{
}



// ------------ method called once each job just before starting event loop  ------------
void ROOTDump::beginJob(){
}








// **********************************************************************
// *                              analyze()                             *
// *                                                                    *
// *                                                                    *
// **********************************************************************


void ROOTDump::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool evValid = true;

  /*

  //This product is missing in the latest ntuple. Commented out to allow the analyser to run.

  edm::Handle<L1TowerJetCollection> CalibJets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibCircle8"), CalibJets);
  if(!CalibJets.isValid()){
  evValid = false;
  edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibCircle8") << std::endl;
  }	
  */
  //calibrated tower MHT/HT collection: L1EtMissParticleCollection
  edm::Handle<l1extra::L1EtMissParticleCollection> L1MHT_up;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibCircle8MHT"),L1MHT_up);
  if(!L1MHT_up.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibCircle8MHT") << std::endl;
  }
  //calibrated tower jet collection: l1extra:: L1JetParticleCollection
  edm::Handle<l1extra::L1JetParticleCollection> CalibJets_L1extra;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibCircle8l1extra"), CalibJets_L1extra);
  if(!CalibJets_L1extra.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibCircle8l1extra") << std::endl;
  }


  //Current L1 extra jets
  vector<InputTag>  l1extraparticles= conf_.getParameter< vector < InputTag > >("extrajet"); 

  for (uint i = 0; i < l1extraparticles.size(); ++i){
    Handle<l1extra::L1JetParticleCollection> currl1Jets;
    iEvent.getByLabel(l1extraparticles[i], currl1Jets );
    if(!currl1Jets.isValid()){evValid = false;}
  }

  //Current HT&MHT
  edm::Handle<l1extra::L1EtMissParticleCollection> L1MHT_curr;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("L1extraMHT"),L1MHT_curr);
  if(!L1MHT_curr.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("L1extraMHT") << std::endl;

  }

  /*  // Rho information for the tower jets
  edm::Handle<double> CaloRho_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CaloRho"), CaloRho_Tower);
  if(!CaloRho_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CaloRho") << std::endl;
  }
  */

  // ************************************************************
  // *                        Tower Jets                        *
  // ************************************************************
  
  
  // Proto jets (TowerJetProducer)
  edm::Handle<l1slhc::L1TowerJetCollection> ProtoJets_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("ProtoJets"), ProtoJets_Tower);
  if(!ProtoJets_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("ProtoJets_Tower") << std::endl;
  }

  // 1DFilteredJets
  edm::Handle<l1slhc::L1TowerJetCollection> Filtered1DJets_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("Filtered1DJets"), Filtered1DJets_Tower);
  if(!Filtered1DJets_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("Filtered1DJets_Tower") << std::endl;
  }
    
  // Upgrade tower jets prior to PU subtraction
  edm::Handle<l1slhc::L1TowerJetCollection> PrePUSubUpgrCenJet_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSubUpgrCenJet"), PrePUSubUpgrCenJet_Tower);
    if(!PrePUSubUpgrCenJet_Tower.isValid()){
      evValid = false;
      edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSubUpgrCenJet") << std::endl;
    }

  // Upgrade tower jets after PU subtraction
  edm::Handle<l1slhc::L1TowerJetCollection> UpgrCenJet_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UpgrCenJet"), UpgrCenJet_Tower);
  if(!UpgrCenJet_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("UpgrCenJet") << std::endl;
  }
  

  edm::Handle<l1slhc::L1TowerJetCollection> LocalPUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("LocalPUSubUpgrCenJet"), LocalPUS_Tower);
  if(!LocalPUS_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("LocalPUSubUpgrCenJet") << std::endl;
  }

  /*
  // Towerjet information for the tower jets
  edm::Handle<double> TowerHT_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("TowerHT"), TowerHT_Tower);
  if(!TowerHT_Tower.isValid()){
  evValid = false;
  edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("TowerHT") << std::endl;
  }
  */
  
  // TT collection
  edm::Handle<l1slhc::L1CaloTowerCollection> caloTowers;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalorimeterTowers"), caloTowers);
  if(!caloTowers.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalorimeterTowers") << std::endl;
    evValid = false;
  }
  
  //////////////////////////////////
  // Handles to offline information
  //////////////////////////////////

  //Need this for information about PU
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RecoVertices"), vtx); 
  if(!vtx.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RecoVertices") << std::endl;
    evValid = false;
  }

  //PU subtracted AK5 calo jets-must be in root file read in
  edm::Handle<reco::CaloJetCollection> Calojets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUsubCaloJets"), Calojets);
  if(!Calojets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUsubCaloJets") << std::endl;
    evValid = false;
  }

  //information about offline calo rho-must be in root file read in
  edm::Handle<double> rhoCALOjets;  
  //  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RhoCaloJets"), rhoCALOjets);
  iEvent.getByLabel("ak5CaloJets","rho", rhoCALOjets);
  if(!rhoCALOjets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RhoCaloJets") << std::endl;
    evValid = false;
  }


  // handle to the jet ID variables
  edm::Handle<reco::JetIDValueMap> hJetIDMap;
  iEvent.getByLabel( conf_.getParameter<edm::InputTag>("jetIDHelperConfig") , hJetIDMap );
  if(!hJetIDMap.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("jetIDHelperConfig")<< std::endl;
    evValid = false;
  }

  edm::Handle<edm::View< reco::CaloJet > > hJets; // uncorrected jets!
  iEvent.getByLabel("ak5CaloJets", hJets );
  if(!hJets.isValid()){
    edm::LogWarning("MissingProduct") << "ak5CaloJets"<< std::endl;
    evValid = false;
  }





  // ------------------------------------------------------------------------------------------------------------------------------------------------------
  // -                                                              Fill Histograms                                                                       -
  // ------------------------------------------------------------------------------------------------------------------------------------------------------


 
  if(!evValid){
    std::cout << "Invalid event\n";	
  }
  else{   // valid event


    // Store event information
    // **************************************************

    // Extract the number of reconstructed AK5 vertices
    //    Int_t PVSIZE = vtx->size();
    NVertices = vtx->size();
    // Get the Luminosity block info from the event
    Runnr     = iEvent.id().run();
    LumiBlock = iEvent.id().luminosityBlock();
    Eventnr   = iEvent.id().event();


    // Store the event data
    evTree->Fill();




    // *****************************************************************************
    // *                     ProtoJet L1 Towerjet collection                       *
    // *****************************************************************************

    fillTowerJetROOTTree( ProtoJets_Tower, ProtoJetTree);


    // *****************************************************************************
    // *                  1D Eta Filtered L1 Towerjet collection                   *
    // *****************************************************************************
  
    fillTowerJetROOTTree( Filtered1DJets_Tower, Filtered1DJetTree);


    // *****************************************************************************
    // *                 Pre Pu-subtracted L1 Towerjet collection                  *
    // *****************************************************************************

    fillTowerJetROOTTree(PrePUSubUpgrCenJet_Tower, PrePUSJetTree);

    
    // *****************************************************************************
    // *                Global Pu-subtracted L1 Towerjet collection                *
    // *****************************************************************************

    fillTowerJetROOTTree(UpgrCenJet_Tower, PUSJetTree);

    // *****************************************************************************
    // *                 Local Pu-subtracted L1 Towerjet collection                *
    // *****************************************************************************

    fillTowerJetROOTTree(LocalPUS_Tower, LPUSJetTree);
    





    
    // ****************************************************************************************************
    // *                                        TT distributions                                          *
    // ****************************************************************************************************

    
#ifdef TT_INFO


    ttE.clear();
    ttH.clear();
    ttiEta.clear();
    ttiPhi.clear();
    ttEcalFG.clear();
    ttHcalFG.clear();
    
    



    for( l1slhc::L1CaloTowerCollection::const_iterator lTT_It = caloTowers->begin() ; 
	 lTT_It != caloTowers->end() ; ++lTT_It ){

      // ****************************************
      // *    Load the calorimeter tower data   *
      // ****************************************
      int E      = lTT_It->E();
      int H      = lTT_It->H();
      //      int EplusH = E + H;
      int iEta   = lTT_It->iEta();
      int iPhi   = lTT_It->iPhi();
      int EcalFG = lTT_It->EcalFG();
      int HcalFG = lTT_It->HcalFG();
    

      // Restrict to central TTs
      if (abs(iEta) > 28)
	continue;


      ttE.push_back(E);
      ttH.push_back(H);
      ttiEta.push_back(iEta);
      ttiPhi.push_back(iPhi);
      ttEcalFG.push_back(EcalFG);
      ttHcalFG.push_back(HcalFG);


    }


    ttTree->Fill();

#endif



    /*


  

    //Fill a vector of TLorentzVectors with offline PU corrected AK5 calo jets
    vector <TLorentzVector> offlineJets;
    TLorentzVector off;

    //    unsigned int idx;

    //L2L3 PU corrected jet collection
    reco::CaloJetCollection::const_iterator it = Calojets->begin() ;

    //uncorrected jet collection for jet ID reference
    for ( edm::View<reco::CaloJet>::const_iterator ibegin = hJets->begin(),
	    iend = hJets->end(), ijet = ibegin; ijet != iend; ++ijet , ++it) {

      //idx = ijet - ibegin; 

      //      edm::RefToBase<reco::CaloJet> jetRef = hJets->refAt(idx);
      //      reco::JetID const & jetId = (*hJetIDMap)[ jetRef ];
      
      //set the TLorentz vector with corrected energy
      off.SetPtEtaPhiM(it->p4().Pt(),it->p4().eta(),it->p4().phi(),it->p4().M());

      // COMMENT: Magic number
      //loose jet selections
      if (it->emEnergyFraction() < 0.01) continue;
      //if (jetId.fHPD>0.98  ) continue;
      //if (jetId.n90Hits<=1) continue;
      
      // pT and eta cuts
      if( fabs(off.Eta()) > 3.0 ) continue;

      // Fill vector with loose ID's calojets
      offlineJets.push_back(off);
    }

    // ****************************************************************************************************
    // *                             Fill (calibrated) jet vectors. Awkward...                            *
    // ****************************************************************************************************


    TLorentzVector up_l1, curr_l1; 
    vector<TLorentzVector> upgradeJets, currentJets;

    // ********************************************************************************
    // *  Iterate through jets, store the Pt, Eta, Phi and M of jets that pass the Pt *
    // *  threshold and Barrel + Endcap acceptance criteria                           *
    // ********************************************************************************

    // ******************************
    // *          Upgrade           *
    // ******************************

    for (l1extra::L1JetParticleCollection::const_iterator il1 = CalibJets_L1extra->begin(); il1!= CalibJets_L1extra->end() ; ++il1 ){

      up_l1.SetPtEtaPhiM(il1->p4().Pt(),il1->p4().Eta(),il1->p4().Phi(),il1->p4().M());

      if( fabs(up_l1.Eta()) > 3.0 ) continue;

      upgradeJets.push_back(up_l1);

    } // End upgrade jets

    // ******************************
    // *          Current           *
    // ******************************

    for (uint i = 0; i < l1extraparticles.size(); ++i){

      Handle<l1extra::L1JetParticleCollection> currl1Jets;
      iEvent.getByLabel(l1extraparticles[i],currl1Jets); 

      for(l1extra::L1JetParticleCollection::const_iterator itr = currl1Jets->begin();  itr != currl1Jets->end(); ++itr ) {

        curr_l1.SetPtEtaPhiM(itr->p4().Pt(),itr->p4().eta(),itr->p4().phi(),itr->p4().M());

        if( fabs(curr_l1.Eta()) > 3.0 ) continue;

	currentJets.push_back(curr_l1);

      }
    } // End current jets

    
*/








  } // end valid event

} // end analyser











  void ROOTDump::fillTowerJetROOTTree(edm::Handle<l1slhc::L1TowerJetCollection> const& TowerJet, TTree *t){


    // Clear previous data
    NJets = 0;
    jetsEt.clear();
    //    jetsEtPUSub.clear();
    jetsIPhi.clear();
    jetsIEta.clear();
    jetsPhiCentre.clear();
    jetsEtaCentre.clear(); 
    jetsPhiWeighted.clear();
    jetsEtaWeighted.clear();
    jetsAsymPhi.clear();
    jetsAsymEta.clear();
    jetsRealArea.clear();


    // Loop over jet collection
    for (l1slhc::L1TowerJetCollection::const_iterator Tower_It = TowerJet->begin(); Tower_It != TowerJet->end(); ++Tower_It ){
	
      // ****************************************
      // *    Load the TowerJet data            *
      // ****************************************
    
      int iEta           = Tower_It->iEta();
      int iPhi           = Tower_It->iPhi();
      //      int E              = Tower_It->E();
      //	bool central       = Tower_It->central();
      int AsymEta        = Tower_It->AsymEta();
      int AsymPhi        = Tower_It->AsymPhi();
      double Pt          = Tower_It->Pt();
      // Unweighted quantities (Geometric jet center)
      double Eta         = Tower_It->Eta();
      double Phi         = Tower_It->Phi();

      //<BUG?????>
      if (Phi > PI)
	Phi -= 2*PI;
      //<BUG?????>

      // Weighted quantities (Energy-weighted jet center)
      double WeightedEta = Tower_It->WeightedEta();
      double WeightedPhi = Tower_It->WeightedPhi();
      //	int JetSize        = Tower_It->JetSize();
      //	int JetArea        = Tower_It->JetArea();
      double JetRealArea = Tower_It->JetRealArea();
      //	double EcalMAD     = Tower_It->EcalMAD();
      //	double HcalMAD     = Tower_It->HcalMAD();
      //	double EnergyMAD   = Tower_It->EnergyMAD();


      // Store data in branches
      // ****************************************
      NJets++;
      jetsEt.push_back(Pt);
      //      jetsEtPUSub.push_back();
      jetsIPhi.push_back(iPhi);
      jetsIEta.push_back(iEta);
      jetsPhiCentre.push_back(Phi);
      jetsEtaCentre.push_back(Eta); 
      jetsPhiWeighted.push_back(WeightedPhi);
      jetsEtaWeighted.push_back(WeightedEta);
      jetsAsymPhi.push_back(AsymPhi);
      jetsAsymEta.push_back(AsymEta);
      jetsRealArea.push_back(JetRealArea);



    } // End TowerJet collection


    // Fill TTree
    // ****************************************

    t->Fill();

}





/*

// ********************************************************************************
// *                  Fill all jet distributions for online jets                  *
// ******************************************************************************** 

void ROOTDump::fillL1Histograms( vector<TLorentzVector> L1Jets, TString prefix, Int_t PVSIZE ){


  for (unsigned int iJet = 0; iJet < L1Jets.size(); iJet++){


    Double_t L1Pt   = L1Jets[iJet].Pt();
    Double_t L1Eta  = L1Jets[iJet].Eta();
    Double_t L1Phi  = L1Jets[iJet].Phi();


  }


}

*/
































// ------------ method called once each job just after ending the event loop  ------------
void ROOTDump::endJob(){
}





// ------------ method called when starting to processes a run  ------------
void 
ROOTDump::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ROOTDump::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ROOTDump::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ROOTDump::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ROOTDump::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}









// ------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------------------------------------




//define this as a plug-in
DEFINE_FWK_MODULE(ROOTDump);
