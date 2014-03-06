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
#include "TEfficiency.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TLorentzVector.h"
#include "TMacro.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"


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
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/TriggerTowerGeometry.h"
//#include "AnalyseUpgradeJets/AnalyzeJets/src/JetMatch.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include <algorithm>  // for sorting

using namespace l1slhc;
using namespace edm;
using namespace std;
using namespace reco;



//#defined VERBOSE



#include "AnalyseUpgradeJets/AnalyzeJets/interface/printing.h"
#include "AnalyseUpgradeJets/AnalyzeJets/interface/histBooker.h"

const double PI = 3.141592654;



// Data structure for jet visualisation
struct CanvasJet{
  
  TEllipse *jetCone;
  double pt;
  double area;

  // Constructor
  CanvasJet( TEllipse *aJetCone, double aPt, double aArea ): jetCone( aJetCone ), pt( aPt ), area( aArea ){}
  
//   // Comparison operators
//   bool operator <(CanvasJet const& rhs) {
//     return pt < rhs.pt;
//   }
//   bool operator >(CanvasJet const& rhs) {
//     return pt > rhs.pt;
//   }


};

bool jetAscending( CanvasJet* jet1, CanvasJet* jet2 ){ return ( jet1->pt < jet2->pt ); }


CanvasJet* makeJetCone( double eta, double phi, double pT, double R, double area, int color, bool transparent);

typedef std::vector<CanvasJet*> CanvasJetCollection;

// for (CaloTowerCollection::const_iterator iCalo = caloFineTowers->begin(); iCalo != caloFineTowers->end(); ++iCalo) {



void drawCanvasJetCollection( CanvasJetCollection canvasJet ){

  for (CanvasJetCollection::const_iterator iJet = canvasJet.begin(); iJet != canvasJet.end(); ++iJet) {                       
    (*iJet)->jetCone->Draw();
  }

}


std::vector <int> colourArr;

// Sort in asscending pT and colour jets
void rankColourCanvasJetCollection( CanvasJetCollection canvasJet ){

  // Order jet collections in reverse pt - Ensures the lead jet is printed on top (i.e. last)                                                                       
  std::sort ( canvasJet.begin(), canvasJet.end(), jetAscending );

  // Colour jets
  int colour;
  for (uint iJet = 0; iJet < canvasJet.size(); ++iJet){
    uint jetIndex =  canvasJet.size() - iJet - 1;
    // Change the colour if jet is a leading jet                                                                                                                    
    if ( jetIndex < colourArr.size() ){
      colour = colourArr[ jetIndex ];
    }
    else{
      colour = kMagenta-7;
    }
    canvasJet[ iJet ]->jetCone->SetLineColor( colour );
  }
}



// Fill a canvasJetCollection with TowerJets. Also make a phi vs eta pT-weighted TH2
CanvasJetCollection fillCanvasJetCollection( edm::Handle<l1slhc::L1TowerJetCollection> L1TowerJet, double etaCut, TH2* histogram ){

  CanvasJetCollection canvasJets;
  
  for (l1slhc::L1TowerJetCollection::const_iterator Tower_It = L1TowerJet->begin(); Tower_It != L1TowerJet->end(); ++Tower_It ){

	double eta = Tower_It->p4().eta();

	if ( fabs(eta) > etaCut ){
	  continue;
	}
	
	double phi  = Tower_It->p4().phi();
	double Pt   = Tower_It->p4().Pt();
 	double area = Tower_It->JetRealArea();
 	double R    = sqrt( area/PI );

	std::cout << "Making L1 cone: " << eta << "\t" << phi << "\t" << Pt << "\tArea = " << area <<  "\n";	
	canvasJets.push_back( makeJetCone( eta, phi, Pt, R, area, kRed, true) );

	// Fill a jet eta vs phi pt-weighted TH2 also...
	int PtRound  = int( Tower_It->p4().Pt() + 0.5 );
        histogram->Fill( eta, phi, PtRound );

      }
  
  // Make pT labels more visible
  histogram->SetMarkerSize( 1.5 );

  return canvasJets;

}


// Fill a canvasJetCollection with CaloJets. Also make a phi vs eta pT-weighted TH2
CanvasJetCollection fillCanvasJetCollection( edm::Handle<reco::CaloJetCollection> CaloJets, double etaCut, TH2* histogram ){

  CanvasJetCollection canvasJets;

  for (reco::CaloJetCollection::const_iterator CaloJet_It = CaloJets->begin(); CaloJet_It != CaloJets->end(); ++CaloJet_It ){ 
    
    double eta = CaloJet_It->p4().eta();
    
    if ( fabs(eta) > etaCut ){
      continue;
    }
    
    double phi = CaloJet_It->p4().phi();
    double Pt  = CaloJet_It->p4().Pt();
    double area = CaloJet_It->towersArea();
    double R    = sqrt( area/PI );
    
    std::cout << "Making ak5 cone: " << eta << "\t" << phi << "\t" << Pt << "\tArea = " << area <<  "\n";	
    canvasJets.push_back( makeJetCone( eta, phi, Pt, R, area, kRed, true) );
    
    // Fill a jet eta vs phi pt-weighted TH2 also...
    int PtRound  = int( CaloJet_It->p4().Pt() + 0.5 );
    histogram->Fill( eta, phi, PtRound );
    
  }

  // Make pT labels more visible
  histogram->SetMarkerSize( 1.5 );
  
  return canvasJets;
  
}

// Fill a canvasJetCollection with L1Jets. Also make a phi vs eta pT-weighted TH2
CanvasJetCollection fillCanvasJetCollection( edm::Handle<l1extra::L1JetParticleCollection> L1Jets, double etaCut, TH2* histogram ){

  CanvasJetCollection canvasJets;

  for (l1extra::L1JetParticleCollection::const_iterator L1Jet_It = L1Jets->begin(); L1Jet_It != L1Jets->end(); ++L1Jet_It ){ 
    
    double eta = L1Jet_It->p4().eta();
    
    if ( fabs(eta) > etaCut ){
      continue;
    }
    
    double phi = L1Jet_It->p4().phi();
    double Pt  = L1Jet_It->p4().Pt();
    double estR =  0.0807 * 11./2.; // No information about radius, estimate it
    double area = PI*( estR * estR );
    double R    = sqrt( area/PI );
    
    std::cout << "Making ak5 cone: " << eta << "\t" << phi << "\t" << Pt << "\tArea = " << area <<  "\n";	
    canvasJets.push_back( makeJetCone( eta, phi, Pt, R, area, kRed, true) );
    
    // Fill a jet eta vs phi pt-weighted TH2 also...
    int PtRound  = int( L1Jet_It->p4().Pt() + 0.5 );
    histogram->Fill( eta, phi, PtRound );
    
  }

  // Make pT labels more visible
  histogram->SetMarkerSize( 1.5 );
  
  return canvasJets;
  
}





// Fill a canvasJetCollection with L1extraparticles. Also make a phi vs eta pT-weighted TH2
CanvasJetCollection fillCanvasJetCollection( std::vector <InputTag> l1extraparticles, const edm::Event& iEvent, double etaCut, TH2* histogram ){

  CanvasJetCollection canvasJets;

  for (uint i = 0; i < l1extraparticles.size(); ++i){

    Handle<l1extra::L1JetParticleCollection> currentL1Jets;
    iEvent.getByLabel(l1extraparticles[i], currentL1Jets);

    for(l1extra::L1JetParticleCollection::const_iterator CurrL1_It = currentL1Jets->begin();  CurrL1_It != currentL1Jets->end(); ++CurrL1_It ) {

      
      double eta = CurrL1_It->p4().eta();
      
      if ( fabs(eta) > etaCut ){
	continue;
      }
      
      double phi = CurrL1_It->p4().phi();
      double Pt  = CurrL1_It->p4().Pt();
      double R   = 6*0.087;
      double area = PI*R*R;
      
      std::cout << "Making current GCT L1 cone: " << eta << "\t" << phi << "\t" << Pt << "\tArea = " << area <<  "\n";	
      canvasJets.push_back( makeJetCone( eta, phi, Pt, R, area, kRed, true) );
      
      // Fill a jet eta vs phi pt-weighted TH2 also...
      int PtRound  = int( CurrL1_It->p4().Pt() + 0.5 );
      histogram->Fill( eta, phi, PtRound );
      
      
      
    }
    
  }


  // Make pT labels more visible
  histogram->SetMarkerSize( 1.5 );
  
  return canvasJets;
  
}






// ========================================
// class declaration
// ========================================

class TTHist : public edm::EDAnalyzer {

public:
  explicit TTHist(const edm::ParameterSet&);
  ~TTHist();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  //------------member functions---------------------
  //  std::pair< double, int> Leading_Match(vector<TLorentzVector> offlineJets,vector<TLorentzVector> L1Jets );

  // ----------member data ---------------------------
  ParameterSet conf_;


  // Histogram containers
  std::map< TString, TH1*> hist1D;
  std::map< TString, TH2*> hist2D;
  std::map< TString, TEfficiency*> histEff;

  // Tower geometry converter
  TriggerTowerGeometry mTowerGeo;

  int NVTX;

  // Maximum and current number of single events to analyse
  int maxEvents;
  int eventCount;
  // Interesting events to analyse
  std::vector<int> eventsOfInterest;

  Service<TFileService> fs;

};





TTHist::TTHist(const edm::ParameterSet& iConfig): conf_(iConfig){

  // ****************************************
  // *              Start job               *
  // ****************************************
  PRINT("Constructor")


    // Associate the histogram container with the histogram booker, a class for handling histogram booking
    //    histBooker booker( &hist1D, &hist2D );


  maxEvents  = 0;
  eventCount = 0;



  //Interesting events
  //

//   // Lead jets
//   eventsOfInterest.push_back( 497615812 ); // Event: 497615812        deltaPt = 0.515791
//   eventsOfInterest.push_back( 498096141 ); // Event: 498096141        deltaPt = 0.504242
//   eventsOfInterest.push_back( 497206795 ); // Event: 497206795        deltaPt = 0.53912
//   eventsOfInterest.push_back( 497516143 ); // Event: 497516143        deltaPt = 0.538958

//   // Any jet
//   eventsOfInterest.push_back( 498505648 ); // Event: 498505648        deltaPt = 0.717097
//   eventsOfInterest.push_back( 498005346 ); // Event: 498005346        deltaPt = 0.936052
//   eventsOfInterest.push_back( 498306497 ); // Event: 498306497        deltaPt = 0.90324


  eventsOfInterest.push_back( 497321201 ); // Event: 497321201        deltaPt = 5.75415       pT = 39.1297    Eta = -0.8091
  eventsOfInterest.push_back( 497695401 ); // Event: 497695401        deltaPt = 0.580491      pT = 43.5172    Eta = -1.85379
  eventsOfInterest.push_back( 497654342 ); // Event: 497654342        deltaPt = 0.356738      pT = 55.7795    Eta = 0.476794

  eventsOfInterest.push_back( 497529475 ); // Event: 497529475        deltaPt = 1.28939       pT = 26.1262    Eta = -0.450818
  eventsOfInterest.push_back( 497497272 ); // Event: 497497272        deltaPt = 1.2716        pT = 29.8352    Eta = -0.837848

  eventsOfInterest.push_back( 497363803 ); // Event: 497363803        deltaPt = 1.12666       pT = 8.46527    Eta = 1.0875


  eventsOfInterest.push_back( 498114341 ); // Event: 498114341        deltaPt = 1.43579       pT = 10.3504    Eta = -0.0145
  eventsOfInterest.push_back( 497240815 ); // Event: 497240815        deltaPt = 1.73436       pT = 12.7855    Eta = -0.7047
  eventsOfInterest.push_back( 497584359 ); // Event: 497584359        deltaPt = 1.68263       pT = 10.7387    Eta = -1.88





  // Visualisation lead jet colours
      colourArr.push_back( kRed );
      colourArr.push_back( kOrange+9 );
      colourArr.push_back( kOrange+7 );
      colourArr.push_back( kOrange-3 );
      colourArr.push_back( kOrange );
      colourArr.push_back( kYellow );
      colourArr.push_back( kSpring+10 );
      colourArr.push_back( kSpring-3 );
      colourArr.push_back( kTeal+2 );




  // **************************
  // *  TT histogram binning  *
  // **************************
  const Int_t TTetaBin = 56;
  const Int_t TTphiBin = 72;

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

  // TT phi binning
  double phiDiv = 2*PI/(TTphiBin);
  double TTphiBins[TTphiBin + 1];
  for (Int_t iPhi = 0; iPhi < TTphiBin + 1; iPhi++) {
    double phi = iPhi*phiDiv - PI;
    TTphiBins[iPhi] = phi;
  }

  TFileDirectory TTDir = fs->mkdir( "TTHist" );




    // ************************************************************
    // *                  Correlation plots                       *
    // ************************************************************


  // CaloTower-TT calibration
  TFileDirectory TTCalib = TTDir.mkdir( "TTCalib" );
  hist2D["TT_vs_Calo_E"]       = TTCalib.make<TH2D>("TT_vs_Calo_E",      "TT vs CaloTower E;Calo E;TT E", 201,-0.05,20.05, 21,-0.25,20.25 );
  hist2D["TT_vs_Calo_H"]       = TTCalib.make<TH2D>("TT_vs_Calo_H",      "TT vs CaloTower H;Calo H;TT H", 201,-0.05,20.05, 21,-0.25,20.25 );
  hist2D["TT_vs_Calo_EplusH"]  = TTCalib.make<TH2D>("TT_vs_Calo_EplusH", "TT vs CaloTower E+H;Calo E+H;TT E+H", 401,-0.05,40.05, 41,-0.25,40.25 );


  hist2D["TT_over_Calo_E_vs_ieta"] = TTCalib.make<TH2D>("TT_over_Calo_E_vs_ieta", "TT/CaloTower E vs i#eta;i#eta;TT/Calo E", 57,-28.5,28.5, 41, -0.5,4.5 );
  hist2D["TT_over_Calo_H_vs_ieta"] = TTCalib.make<TH2D>("TT_over_Calo_H_vs_ieta", "TT/CaloTower H vs i#eta;i#eta;TT/Calo H", 57,-28.5,28.5, 41, -0.5,4.5 );

  hist1D["TT_over_Calo_E_vs_ieta_prof"] = TTCalib.make<TProfile>("TT_over_Calo_E_vs_ieta_prof", "TT/CaloTower E vs i#eta;i#eta;TT/Calo E", 57,-28.5,28.5, 0,4.5 );
  hist1D["TT_over_Calo_H_vs_ieta_prof"] = TTCalib.make<TProfile>("TT_over_Calo_H_vs_ieta_prof", "TT/CaloTower H vs i#eta;i#eta;TT/Calo H", 57,-28.5,28.5, 0,4.5 );


  // Thresholds
  hist2D["TT_over_Calo_E_vs_ieta_lt2GeV"] = TTCalib.make<TH2D>("TT_over_Calo_E_vs_ieta_lt2GeV", "TT/CaloTower E vs i#eta, E < 2 GeV;i#eta;TT/Calo E", 57,-28.5,28.5, 41, -0.5,4.5 );
  hist2D["TT_over_Calo_H_vs_ieta_lt2GeV"] = TTCalib.make<TH2D>("TT_over_Calo_H_vs_ieta_lt2GeV", "TT/CaloTower H vs i#eta, H < 2 GeV ;i#eta;TT/Calo H", 57,-28.5,28.5, 41, -0.5,4.5 );

  hist2D["TT_over_Calo_E_vs_ieta_ge2GeV"] = TTCalib.make<TH2D>("TT_over_Calo_E_vs_ieta_ge2GeV", "TT/CaloTower E vs i#eta, E #ge 2 GeV;i#eta;TT/Calo E", 57,-28.5,28.5, 41, -0.5,4.5 );
  hist2D["TT_over_Calo_H_vs_ieta_ge2GeV"] = TTCalib.make<TH2D>("TT_over_Calo_H_vs_ieta_ge2GeV", "TT/CaloTower H vs i#eta, H #ge 2 GeV ;i#eta;TT/Calo H", 57,-28.5,28.5, 41, -0.5,4.5 );





  // CaloTower-TT correlations    
  TFileDirectory TTCorr = TTDir.mkdir( "TTCorr" );

    TFileDirectory corrSubDir = TTCorr.mkdir( "Correlations" );

    // RCT
//     hist2D["RCT_E_vs_NVTX"]     = corrSubDir.make<TH2D>( "RCT_E_vs_NVTX",     "RCT region energy vs N_{VTX};N_{VTX};Region E (GeV)",       
// 							 51,-0.5,50.5, 51,-0.5,50.5 );
    hist2D["RCT_Etot_vs_NVTX"]  = corrSubDir.make<TH2D>( "RCT_Etot_vs_NVTX",  "RCT total region energy vs N_{VTX};N_{VTX};Region E (GeV)", 
							 51,-0.5,50.5, 101,-2,402 );
    hist2D["RCT_Fired_vs_NVTX"] = corrSubDir.make<TH2D>( "RCT_Fired_vs_NVTX", "RCT regions fired vs N_{VTX};N_{VTX};Regions fired",   
							 51,-0.5,50.5,     101,-1,201 );
    // TT
//     hist2D["TT_E_vs_NVTX"]     = corrSubDir.make<TH2D>( "TT_E_vs_NVTX",     "TT region energy vs N_{VTX};N_{VTX};Region E (GeV)",       
// 							51,-0.5,50.5, 51,-0.5,50.5 );
    hist2D["TT_Etot_vs_NVTX"]  = corrSubDir.make<TH2D>( "TT_Etot_vs_NVTX",  "TT total region energy vs N_{VTX};N_{VTX};Region E (GeV)", 
							51,-0.5,50.5, 101,-2,402 );
    hist2D["TT_Fired_vs_NVTX"] = corrSubDir.make<TH2D>( "TT_Fired_vs_NVTX", "TT regions fired vs N_{VTX};N_{VTX};Regions fired",       
							51,-0.5,50.5, 101,-1,201 );
    // Calo
//     hist2D["Calo_E_vs_NVTX"]     = corrSubDir.make<TH2D>( "Calo_E_vs_NVTX",     "Calo region energy vs N_{VTX};N_{VTX};Region E (GeV)",       
// 							  51,-0.5,50.5, 51,-0.5,50.5 );
    hist2D["Calo_Etot_vs_NVTX"]  = corrSubDir.make<TH2D>( "Calo_Etot_vs_NVTX",  "Calo total region energy vs N_{VTX};N_{VTX};Region E (GeV)",
							  51,-0.5,50.5, 101,-5,1005 );
    hist2D["Calo_Fired_vs_NVTX"] = corrSubDir.make<TH2D>( "Calo_Fired_vs_NVTX", "Calo regions fired vs N_{VTX};N_{VTX};Regions fired",      
							  51,-0.5,50.5, 151,-5,1505 );




  // TT-RCT region energy comparisons
  TFileDirectory RCTRegions = TTCorr.mkdir( "RCTRegions" );
  for (int RCTiEta = 4; RCTiEta < 18; ++RCTiEta){
    TString RCTiEtaStr = Form("%d", RCTiEta );

    // Region phi strips
    hist2D["RCT_vs_TT_ET_RCTiEta_" + RCTiEtaStr]   = RCTRegions.make<TH2D>("RCT_vs_TT_ET_RCTiEta_" + RCTiEtaStr,"RCT vs TT E_{T} RCT i#eta = " + RCTiEtaStr
								     + ";TT E;RCT E",  51,-0.5,50.5, 51,-0.5,50.5 );
    hist2D["Calo_vs_TT_ET_RCTiEta_" + RCTiEtaStr]  = RCTRegions.make<TH2D>("Calo_vs_TT_ET_RCTiEta_" + RCTiEtaStr,"Calo vs TT E_{T} RCT i#eta = " + RCTiEtaStr
								     + ";TT E;Calo E", 51,-0.5,50.5, 51,-0.5,50.5 );
    hist2D["Calo_vs_RCT_ET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH2D>("Calo_vs_RCT_ET_RCTiEta_" + RCTiEtaStr,"Calo vs RCT E_{T} RCT i#eta = " + RCTiEtaStr
								     + ";RCT E;Calo E", 51,-0.5,50.5, 51,-0.5,50.5 );
    
    hist1D["RCT_ET_Response_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("RCT_ET_Response_RCTiEta_" + RCTiEtaStr,"RCT E response RCT i#eta = " + RCTiEtaStr
									+ ";E^{RCT}/E^{Calo};Entries", 51,-0.05,5.05 );
    hist1D["TT_ET_Response_RCTiEta_"  + RCTiEtaStr] = RCTRegions.make<TH1D>("TT_ET_Response_RCTiEta_" + RCTiEtaStr,"TT E Response RCT i#eta = " + RCTiEtaStr
									+ ";E^{TT}/E^{Calo};Entries",  51,-0.05,5.05 );

//     hist1D["RCT_ET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("RCT_ET_RCTiEta_" + RCTiEtaStr,"RCT E RCT i#eta = " + RCTiEtaStr
// 									+ ";E^{RCT};Entries", 51,-0.5,50.5 );
//     hist1D["TT_ET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("TT_ET_RCTiEta_" + RCTiEtaStr,"TT E RCT i#eta = " + RCTiEtaStr
// 									+ ";E^{TT};Entries", 51,-0.5,50.5 );
//     hist1D["Calo_ET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("Calo_ET_RCTiEta_" + RCTiEtaStr,"Calo E RCT i#eta = " + RCTiEtaStr
// 									+ ";E^{Calo};Entries", 51,-0.5,50.5 );

    hist1D["RCT_RegionStripET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("RCT_RegionStripET_RCTiEta_" + RCTiEtaStr,"RCT region strip E RCT i#eta = " + RCTiEtaStr
									+ ";E^{RCT};Entries", 51,-0.5,50.5 );
    hist1D["TT_RegionStripET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("TT_RegionStripET_RCTiEta_" + RCTiEtaStr,"TT region strip E RCT i#eta = " + RCTiEtaStr
									+ ";E^{TT};Entries", 51,-0.5,50.5 );
    hist1D["Calo_RegionStripET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("Calo_RegionStripET_RCTiEta_" + RCTiEtaStr,"Calo region strip E RCT i#eta = " + RCTiEtaStr
									+ ";E^{Calo};Entries", 51,-0.5,50.5 );


    // Regions
    hist2D["RCT_vs_TT_RegionET_RCTiEta_" + RCTiEtaStr]   = RCTRegions.make<TH2D>("RCT_vs_TT_RegionET_RCTiEta_" + RCTiEtaStr,"RCT vs TT region E_{T} RCT i#eta = " 
									     + RCTiEtaStr
								     + ";TT E;RCT E",  51,-0.5,50.5, 51,-0.5,50.5 );
    hist2D["Calo_vs_TT_RegionET_RCTiEta_" + RCTiEtaStr]  = RCTRegions.make<TH2D>("Calo_vs_TT_RegionET_RCTiEta_" + RCTiEtaStr,"Calo vs TT region E_{T} RCT i#eta = "
									     + RCTiEtaStr
								     + ";TT E;Calo E", 51,-0.5,50.5, 51,-0.5,50.5 );
    hist2D["Calo_vs_RCT_RegionET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH2D>("Calo_vs_RCT_RegionETRCTiEta_" + RCTiEtaStr,"Calo vs RCT region E_{T} RCT i#eta = " 
									     + RCTiEtaStr
								     + ";RCT E;Calo E", 51,-0.5,50.5, 51,-0.5,50.5 );

    hist1D["RCT_RegionET_Response_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("RCT_RegionET_Response_RCTiEta_" + RCTiEtaStr,"RCT Region E response RCT i#eta = "
									      + RCTiEtaStr
									      + ";E^{RCT}/E^{Calo};Entries", 51,-0.05,5.05 );
    hist1D["TT_RegionET_Response_RCTiEta_"  + RCTiEtaStr] = RCTRegions.make<TH1D>("TT_RegionET_Response_RCTiEta_" + RCTiEtaStr,"TT Region E Response RCT i#eta = " 
									      + RCTiEtaStr
									      + ";E^{TT}/E^{Calo};Entries",  51,-0.05,5.05 );

    hist1D["RCT_RegionET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("RCT_RegionET_RCTiEta_" + RCTiEtaStr,"RCT region E RCT i#eta = " + RCTiEtaStr
									+ ";E^{RCT};Entries", 51,-0.5,50.5 );
    hist1D["TT_RegionET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("TT_RegionET_RCTiEta_" + RCTiEtaStr,"TT region E RCT i#eta = " + RCTiEtaStr
									+ ";E^{TT};Entries", 51,-0.5,50.5 );
    hist1D["Calo_RegionET_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH1D>("Calo_RegionET_RCTiEta_" + RCTiEtaStr,"Calo region E RCT i#eta = " + RCTiEtaStr
									+ ";E^{Calo};Entries", 51,-0.5,50.5 );



    hist2D["RCT_E_vs_NVTX_RCTiEta_" + RCTiEtaStr]     = corrSubDir.make<TH2D>( "RCT_E_vs_NVTX_RCTiEta_" + RCTiEtaStr, 
									       "RCT strip region energy vs N_{VTX} RCT i#eta = " + 
									       RCTiEtaStr + ";N_{VTX};Region E (GeV)",      51,-0.5,50.5,  51,-0.5,50.5 );

    hist2D["TT_E_vs_NVTX_RCTiEta_" + RCTiEtaStr]     = corrSubDir.make<TH2D>( "TT_E_vs_NVTX_RCTiEta_" + RCTiEtaStr, 
									       "TT strip region energy vs N_{VTX} RCT i#eta = " + 
									       RCTiEtaStr + ";N_{VTX};Region E (GeV)",      51,-0.5,50.5,  51,-0.5,50.5 );

    hist2D["Calo_E_vs_NVTX_RCTiEta_" + RCTiEtaStr]     = corrSubDir.make<TH2D>( "Calo_E_vs_NVTX_RCTiEta_" + RCTiEtaStr, 
									       "Calo strip region energy vs N_{VTX} RCT i#eta = " + 
									       RCTiEtaStr + ";N_{VTX};Region E (GeV)",      51,-0.5,50.5,  51,-0.5,50.5 );




    // deltaE RegionStrip

    hist2D["deltaRegionStripET_TTRCT_vs_TTFired_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH2D>("deltaRegionStripET_TTRCT_vs_TTFired_RCTiEta_" + RCTiEtaStr,
												"#DeltaRegionStripE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} vs TT fired;TT fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
												34,-0.5,33.5, 101,-0.25,50.25);
    hist1D["deltaRegionStripET_TTRCT_vs_TTFired_RCTiEta_" + RCTiEtaStr + "_prof"] = RCTRegions.make<TProfile>("deltaRegionStripET_TTRCT_vs_TTFired_RCTiEta_" 
													      + RCTiEtaStr + "_prof",
													      "#DeltaRegionStripE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) vs TT fired;TT fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
													      34,-0.5,33.5, -0.25,50.25);
    
 hist2D["deltaRegionStripET_TTRCT_vs_Calo0p5Fired_RCTiEta_" + RCTiEtaStr] = RCTRegions.make<TH2D>("deltaRegionStripET_TTRCT_vs_Calo0p5Fired_RCTiEta_" 
												  + RCTiEtaStr,
												  "#DeltaRegionStripE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} vs Calo > 0.5 GeV fired;Calo fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
												  71,-0.5,70.5, 101,-0.25,50.25);
 hist1D["deltaRegionStripET_TTRCT_vs_Calo0p5Fired_RCTiEta_" + RCTiEtaStr + "_prof"] = RCTRegions.make<TProfile>("deltaRegionStripET_TTRCT_vs_Calo0p5FiredRCTiEta_" + RCTiEtaStr + "_prof",
														"#DeltaRegionStripE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) vs Calo > 0.5 GeV fired;Calo fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
														71,-0.5,70.5, -0.25,50.25);
 
 


  }

    hist1D["RCT_RegionET"]  = TTCorr.make<TH1D>("RCT_RegionET","RCT region E;E^{RCT};Entries", 51,-0.5,50.5 );
    hist1D["TT_RegionET"]   = TTCorr.make<TH1D>("TT_RegionET","TT region E;E^{TT};Entries", 51,-0.5,50.5 );
    hist1D["Calo_RegionET"] = TTCorr.make<TH1D>("Calo_RegionET","Calo region E;E^{Calo};Entries", 51,-0.5,50.5 );

    hist2D["RCT_vs_TT_RegionET"]   = TTCorr.make<TH2D>("RCT_vs_TT_RegionET","RCT vs TT region E_{T};TT E;RCT E",  51,-0.5,50.5, 51,-0.5,50.5 );
    hist2D["Calo_vs_TT_RegionET"]  = TTCorr.make<TH2D>("Calo_vs_TT_RegionET","Calo vs TT region E_{T};TT E;Calo E", 51,-0.5,50.5, 51,-0.5,50.5 );
    hist2D["Calo_vs_RCT_RegionET"] = TTCorr.make<TH2D>("Calo_vs_RCT_ET_RegionET","Calo vs RCT region E_{T};RCT E;Calo E", 51,-0.5,50.5, 51,-0.5,50.5 );









    std::vector<TString> prefixArr;
    prefixArr.push_back("TT");
    prefixArr.push_back("Calo");


    for (uint iPre = 0; iPre < prefixArr.size(); ++iPre){

      TString TTPre = prefixArr[ iPre ];
      TString TTLab = TTPre;

      
      
      // Phi-Eta space
      hist2D[TTPre + "_E-Phi_vs_Eta"]          = TTDir.make<TH2D>(TTPre + "_E-Phi_vs_Eta","Ecal energy #phi vs #eta " +
								  TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      hist2D[TTPre + "_H-Phi_vs_Eta"]          = TTDir.make<TH2D>(TTPre + "_H-Phi_vs_Eta","Hcal energy #phi vs #eta " +
								  TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      hist2D[TTPre + "_EplusH-Phi_vs_Eta"]     = TTDir.make<TH2D>(TTPre + "_EplusH-Phi_vs_Eta","Ecal + Hcal energy #phi vs #eta " +
								  TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      //   hist2D[TTPre + "_EHratio-Phi_vs_Eta"]    = TTDir.make<TH2D>(TTPre + "_EHratio-Phi_vs_Eta","Ecal-Hcal energy ratio #phi vs #eta " +
      //                                                       TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      hist2D[TTPre + "_EFired-Phi_vs_Eta"]     = TTDir.make<TH2D>(TTPre + "_EFired-Phi_vs_Eta","Ecal fired #phi vs #eta " +
								  TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      hist2D[TTPre + "_HFired-Phi_vs_Eta"]     = TTDir.make<TH2D>(TTPre + "_HFired-Phi_vs_Eta","Hcal fired #phi vs #eta " +
								  TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      hist2D[TTPre + "_EorHFired-Phi_vs_Eta"]  = TTDir.make<TH2D>(TTPre + "_EorHFired-Phi_vs_Eta","Ecal or Hcal fired #phi vs #eta " +
								  TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      hist2D[TTPre + "_EandHFired-Phi_vs_Eta"] = TTDir.make<TH2D>(TTPre + "_EandHFired-Phi_vs_Eta","Ecal and Hcal fired #phi vs #eta " +
								  TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      // Energy sums
      hist1D[TTPre + "_Etotal"]        = TTDir.make<TH1D>(TTPre + "_Etotal","Total Ecal energy deposited " +                                                
							  TTLab + ";Total Ecal energy (GeV);Events", 126,-20,5020);  
      hist1D[TTPre + "_Htotal"]        = TTDir.make<TH1D>(TTPre + "_Htotal","Total Hcal energy deposited " +                                                
							  TTLab + ";Total Hcal energy (GeV);Events", 126,-20,5020);   
      hist1D[TTPre + "_EplusHtotal"]   = TTDir.make<TH1D>(TTPre + "_EplusHtotal","Total Ecal + Hcal energy deposited " +                                   
							  TTLab + ";Total Ecal + Hcal energy (GeV);Events", 101,-50,10050);

    }

    // Energy sum correlations
    hist2D["TT_vs_Calo_ET"]            = TTDir.make<TH2D>("TT_vs_Calo_ET","TT E_{T} vs Calo E_{T};Calo E_{T} (GeV);TT E_{T} (GeV)", 
							  151,-5,1505, 151,-5,1505);
    hist1D["TT_vs_Calo_ET_prof"]       = TTDir.make<TProfile>("TT_vs_Calo_ET_prof","TT E_{T} vs Calo E_{T};Calo E_{T} (GeV);TT E_{T} (GeV)", 
							  151,-5,1505, -5,1505);

    hist2D["RCT_vs_Calo_ET"]            = TTDir.make<TH2D>("RCT_vs_Calo_ET","RCT E_{T} vs Calo E_{T};Calo E_{T} (GeV);RT E_{T} (GeV)", 
							  151,-5,1505, 151,-5,1505);
    hist1D["RCT_vs_Calo_ET_prof"]       = TTDir.make<TProfile>("RCT_vs_Calo_ET_prof","RCT E_{T} vs Calo E_{T};Calo E_{T} (GeV);RT E_{T} (GeV)", 
							  151,-5,1505, -5,1505);

    hist2D["TT_vs_RCT_ET"]            = TTDir.make<TH2D>("TT_vs_RCT_ET","TT E_{T} vs RCT E_{T};RCT E_{T} (GeV);TT E_{T} (GeV)", 
							  151,-5,1505, 151,-5,1505);
    hist1D["TT_vs_RCT_ET_prof"]       = TTDir.make<TProfile>("TT_vs_RCT_ET_prof","TT E_{T} vs RCT E_{T};RCT E_{T} (GeV);TT E_{T} (GeV)", 
							  151,-5,1505, -5,1505);


    hist2D["deltaET_TTRCT_vs_TTFired"] = TTDir.make<TH2D>("deltaET_TTRCT_vs_TTFired",
							  "#DeltaE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} vs TT fired;TT fired;#DeltaE_{T}^{(RCT,TT)} (GeV)", 
							  151,-0.5,150.5, 151,-0.5,150.5);
    hist1D["deltaET_TTRCT_vs_TTFired_prof"] = TTDir.make<TProfile>("deltaET_TTRCT_vs_TTFired_prof",
								   "#DeltaE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) vs TT fired;TT fired;#DeltaE_{T}^{(RCT,TT)} (GeV)", 
								   151,-0.5,150.5, -0.5,150.5);

    hist2D["deltaET_TTRCT_vs_Calo0p5Fired"] = TTDir.make<TH2D>("deltaET_TTRCT_vs_Calo0p5Fired",
                                                          "#DeltaE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} vs Calo > 0.5 GeV fired;Calo fired;#DeltaE_{T}^{(RCT,TT)} (GeV)",
							       301,-0.5,300.5, 151,-0.5,150.5);
    hist1D["deltaET_TTRCT_vs_Calo0p5Fired_prof"] = TTDir.make<TProfile>("deltaET_TTRCT_vs_Calo0p5Fired_prof",
                                                                   "#DeltaE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) vs Calo > 0.5 GeV fired;Calo fired;#DeltaE_{T}^{(RCT,TT)} (GeV)",
                                                                   301,-0.5,300.5, -0.5,150.5);


    // Corrected for RCT energy loss from thresholding
    hist2D["deltaET_TTRCTCorrected_vs_TTFired"] = TTDir.make<TH2D>("deltaET_TTRCTCorrected_vs_TTFired",
                                                          "#DeltaE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} RCT corrected vs TT fired;TT fired;#DeltaE_{T}^{(RCT,TT)} (GeV)",
                                                          151,-0.5,150.5, 101,-0.5,100.5);
    hist1D["deltaET_TTRCTCorrected_vs_TTFired_prof"] = TTDir.make<TProfile>("deltaET_TTRCTCorrected_vs_TTFired_prof",
                                                                   "#DeltaE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) RCT corrected vs TT fired;TT fired;#DeltaE_{T}^{(RCT,TT)} (GeV)",
                                                                   151,-0.5,150.5, -0.5,100.5);

    hist2D["deltaET_TTRCTCorrected_vs_Calo0p5Fired"] = TTDir.make<TH2D>("deltaET_TTRCTCorrected_vs_Calo0p5Fired",
                                                          "#DeltaE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} RCT Corrected vs Calo > 0.5 GeV fired;Calo fired;#DeltaE_{T}^{(RCT,TT)} (GeV)",
							       301,-0.5,300.5, 151,-0.5,150.5);
    hist1D["deltaET_TTRCTCorrected_vs_Calo0p5Fired_prof"] = TTDir.make<TProfile>("deltaET_TTRCTCorrected_vs_Calo0p5Fired_prof",
                                                                   "#DeltaE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) RCT Corrected vs Calo > 0.5 GeV fired;Calo fired;#DeltaE_{T}^{(RCT,TT)} (GeV)",
                                                                   301,-0.5,300.5, -0.5,150.5);



    // DeltaRegionET
    hist2D["deltaRegionET_TTRCT_vs_TTFired"] = TTDir.make<TH2D>("deltaRegionET_TTRCT_vs_TTFired",
								"#DeltaRegionE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} vs TT fired;TT fired in region;#DeltaRegionE_{T}^{(RCT,TT)} (GeV)",
								17,-0.5,16.5, 121,-10.25,50.25);
    hist1D["deltaRegionET_TTRCT_vs_TTFired_prof"] = TTDir.make<TProfile>("deltaRegionET_TTRCT_vs_TTFired_prof",
									 "#DeltaRegionE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) vs TT fired;TT fired in region;#DeltaRegionE_{T}^{(RCT,TT)} (GeV)",
									 17,-0.5,16.5, -10.25,50.25);

    hist2D["deltaRegionET_TTRCT_vs_Calo0p5Fired"] = TTDir.make<TH2D>("deltaRegionET_TTRCT_vs_Calo0p5Fired",
								"#DeltaRegionE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} vs Calo > 0.5 GeV fired;Calo fired in region;#DeltaRegionE_{T}^{(RCT,TT)} (GeV)",
								17,-0.5,16.5, 121,-10.25,50.25);
    hist1D["deltaRegionET_TTRCT_vs_Calo0p5Fired_prof"] = TTDir.make<TProfile>("deltaRegionET_TTRCT_vs_Calo0p5Fired_prof",
									 "#DeltaRegionE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) vs Calo > 0.5 GeV fired;Calo fired in region;#DeltaRegionE_{T}^{(RCT,TT)} (GeV)",
									 17,-0.5,16.5, -10.25,50.25);



    // DeltaRegionStripET
    hist2D["deltaRegionStripET_TTRCT_vs_TTFired"] = TTDir.make<TH2D>("deltaRegionStripET_TTRCT_vs_TTFired",
								"#DeltaRegionStripE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} vs TT fired;TT fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
								     34,-0.5,33.5, 101,-0.25,50.25);
    hist1D["deltaRegionStripET_TTRCT_vs_TTFired_prof"] = TTDir.make<TProfile>("deltaRegionStripET_TTRCT_vs_TTFired_prof",
									 "#DeltaRegionStripE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) vs TT fired;TT fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
									      34,-0.5,33.5, -0.25,50.25);

    hist2D["deltaRegionStripET_TTRCT_vs_Calo0p5Fired"] = TTDir.make<TH2D>("deltaRegionStripET_TTRCT_vs_Calo0p5Fired",
								"#DeltaRegionStripE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} vs Calo > 0.5 GeV fired;Calo fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
									 71,-0.5,70.5, 101,-0.25,50.25);
    hist1D["deltaRegionStripET_TTRCT_vs_Calo0p5Fired_prof"] = TTDir.make<TProfile>("deltaRegionStripET_TTRCT_vs_Calo0p5Fired_prof",
									 "#DeltaRegionStripE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) vs Calo > 0.5 GeV fired;Calo fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
									 71,-0.5,70.5, -0.25,50.25);


    hist2D["deltaRegionStripETCorrected_TTRCT_vs_TTFired"] = TTDir.make<TH2D>("deltaRegionStripETCorrected_TTRCT_vs_TTFired",
								"#DeltaRegionStripE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} corrected vs TT fired;TT fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
								     34,-0.5,33.5, 101,-0.25,50.25);
    hist1D["deltaRegionStripETCorrected_TTRCT_vs_TTFired_prof"] = TTDir.make<TProfile>("deltaRegionStripETCorrected_TTRCT_vs_TTFired_prof",
									 "#DeltaRegionStripE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) corrected vs TT fired;TT fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
									      34,-0.5,33.5, -0.25,50.25);

    hist2D["deltaRegionStripETCorrected_TTRCT_vs_Calo0p5Fired"] = TTDir.make<TH2D>("deltaRegionStripETCorrected_TTRCT_vs_Calo0p5Fired",
								"#DeltaRegionStripE_{T}^{(RCT,TT)} RCT E_{T} - TT E_{T} corrected vs Calo > 0.5 GeV fired;Calo fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
									 71,-0.5,70.5, 101,-0.25,50.25);
    hist1D["deltaRegionStripETCorrected_TTRCT_vs_Calo0p5Fired_prof"] = TTDir.make<TProfile>("deltaRegionStripETCorrected_TTRCT_vs_Calo0p5Fired_prof",
									 "#DeltaRegionStripE_{T}^{(RCT,TT)} (RCT E_{T} - TT E_{T}) corrected vs Calo > 0.5 GeV fired;Calo fired in region;#DeltaRegionStripE_{T}^{(RCT,TT)} (GeV)",
									 71,-0.5,70.5, -0.25,50.25);



    // L1 Unmatched jets - L1 jets are defined as lead, second, ...
    histEff[ "TT_MET_TurnOn_gt30" ] = TTDir.make<TEfficiency>( "TT_MET_TurnOn_gt30", "TT MET > 30 GeV Differential turn-on;Calo MET (GeV);Efficiency", 20, 0, 200 );
    histEff[ "TT_MET_TurnOn_gt40" ] = TTDir.make<TEfficiency>( "TT_MET_TurnOn_gt40", "TT MET > 40 GeV Differential turn-on;Calo MET (GeV);Efficiency", 20, 0, 200 );
    histEff[ "TT_MET_TurnOn_gt50" ] = TTDir.make<TEfficiency>( "TT_MET_TurnOn_gt50", "TT MET > 50 GeV Differential turn-on;Calo MET (GeV);Efficiency", 20, 0, 200 );
    histEff[ "CurrL1_MET_TurnOn_gt30" ] = TTDir.make<TEfficiency>( "CurrL1_MET_TurnOn_gt30", "Current L1 MET > 30 GeV Differential turn-on;Calo MET (GeV);Efficiency", 20, 0, 200 );
    histEff[ "CurrL1_MET_TurnOn_gt40" ] = TTDir.make<TEfficiency>( "CurrL1_MET_TurnOn_gt40", "Current L1 MET > 40 GeV Differential turn-on;Calo MET (GeV);Efficiency", 20, 0, 200 );
    histEff[ "CurrL1_MET_TurnOn_gt50" ] = TTDir.make<TEfficiency>( "CurrL1_MET_TurnOn_gt50", "Current L1 MET > 50 GeV Differential turn-on;Calo MET (GeV);Efficiency", 20, 0, 200 );


    histEff[ "TT_ET_TurnOn_gt50" ]  = TTDir.make<TEfficiency>( "TT_ET_TurnOn_gt30", "TT E_{T} > 50 GeV Differential turn-on;Calo E_{T} (GeV);Efficiency", 50, 0, 500 );
    histEff[ "TT_ET_TurnOn_gt100" ] = TTDir.make<TEfficiency>( "TT_ET_TurnOn_gt40", "TT E_{T} > 150 GeV Differential turn-on;Calo E_{T} (GeV);Efficiency", 50, 0, 500 );
    histEff[ "TT_ET_TurnOn_gt150" ] = TTDir.make<TEfficiency>( "TT_ET_TurnOn_gt50", "TT E_{T} > 150 GeV Differential turn-on;Calo E_{T} (GeV);Efficiency", 50, 0, 500 );


//   TFileDirectory TTthresholdSubDir = kt6Dir.mkdir( "TTthreshold" );

//   for (unsigned int ttI = 0; ttI < ttThreshold.size(); ttI++){
      
//     // get TT threshold
//     int ttThresh   = ttThreshold[ttI];
//     TString ttThreshStr = Form("%d",int(ttThresh));
    
//     hist2D["TTFired_vs_kt6Rho_TT_gt_" + ttThreshStr]           = TTthresholdSubDir.make<TH2D>("TTFired_vs_kt6Rho_TT_gt_" + ttThreshStr, "TTs fired (E_{TT} > " + ttThreshStr + 
// 											      ") vs Offline kt6 #rho;kt6 #rho;TTs fired",
// 											      161,-0.125,40.125, 201,-0.5,200.5);
//     hist1D["TTFired_vs_kt6Rho_TT_gt_" + ttThreshStr + "_prof"] = TTthresholdSubDir.make<TProfile>("TTFired_vs_kt6Rho_TT_gt_" + ttThreshStr + "_prof",
// 												  "TTs fired (E_{TT} > " + ttThreshStr + 
// 												  ") vs Offline kt6 #rho profile;kt6 #rho;TTs fired",
// 												  161,-0.125,40.125, -0.5,200.5);
    
    
//   }



}



TTHist::~TTHist()
{
}



// ------------ method called once each job just before starting event loop  ------------
void TTHist::beginJob(){





}








// **********************************************************************
// *                              analyze()                             *
// *                                                                    *
// *                                                                    *
// **********************************************************************


void TTHist::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){


  //  bool evValid = true;

  //Need this for information about PU
  edm::Handle<int> vertices;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("NumPrimaryVertices"), vertices);
  if(!vertices.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("NumPrimaryVertices") << std::endl;
  }
  NVTX       = *vertices;


  // RCT regions 
  edm::Handle<L1CaloRegionCollection> caloRegions;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CaloRegions"), caloRegions );
   if(!caloRegions.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CaloRegions") << std::endl;
     //     RCTData = false; 

   } 








    // Current MET and ET
    SUBPRINT("Current L1 MET")
      edm::Handle< l1extra::L1EtMissParticleCollection > l1CurrMET;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CurrentL1MET"), l1CurrMET );
  if(!l1CurrMET.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CurrentL1MET") << std::endl;
    //    evValid = false;
  }


    // Current MHT and HT
    SUBPRINT("Current L1 MHT")
    edm::Handle< std::vector<l1extra::L1EtMissParticle> > l1CurrMHT;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CurrentL1MHT"), l1CurrMHT );
  if(!l1CurrMHT.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CurrentL1MHT") << std::endl;
    //    evValid = false;
  }


   // ********************************************************************************
   // *                                   Towers                                     *
   // ********************************************************************************

  // TT collection
  SUBPRINT("Trigger towers")
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




   // ********************************************************************************
   // *                                    Jets                                      *
   // ********************************************************************************



  // ************************************************************
  // *                       Current Jets                       *
  // ************************************************************
  SUBPRINT("Current Jets")
    
   // Get L1 jets to be used, currently Cental and Tau jets
    std::vector <InputTag> l1extraparticles = conf_.getParameter< std::vector < InputTag > >("extrajet");
   // Loop through Central and tau jets                                                                                                                         
   for (uint i = 0; i < l1extraparticles.size(); ++i){
     Handle<l1extra::L1JetParticleCollection> currl1Jets;
     iEvent.getByLabel(l1extraparticles[i], currl1Jets );
     if(!currl1Jets.isValid()){
       //evValid = false;
       edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("extrajet") << std::endl;
     }
   }


  // ************************************************************
  // *                        Tower Jets                        *
  // ************************************************************
  SUBPRINT("Tower Jets")

//   // Proto tower jets
//   edm::Handle<l1slhc::L1TowerJetCollection> Proto_Tower;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("ProtoTowerJet"), Proto_Tower);
//   if(!Proto_Tower.isValid()){
//     //evValid = false;
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("ProtoTowerJet") << std::endl;
//   }
  
//   // Centrality filtered tower jets
//   edm::Handle<l1slhc::L1TowerJetCollection> FilteredCentrality_Tower;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCentralityTowerJet"), FilteredCentrality_Tower);
//   if(!FilteredCentrality_Tower.isValid()){
//     //evValid = false;
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("FilteredCentralityTowerJet") << std::endl;
//   }

//   // 1D filtered jets
//   edm::Handle<l1slhc::L1TowerJetCollection> Filtered1D_Tower;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("Filtered1DTowerJet"), Filtered1D_Tower);
//   if(!Filtered1D_Tower.isValid()){
//     //evValid = false;
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("Filtered1DTowerJet") << std::endl;
//   }

  // Pre-PU subtracted jets
  edm::Handle<l1slhc::L1TowerJetCollection> PrePUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSubTowerJet"), PrePUS_Tower);
  if(!PrePUS_Tower.isValid()){
    //evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSubTowerJet") << std::endl;
  }
  
  // Global PU subtracted tower jets 
  edm::Handle<l1slhc::L1TowerJetCollection> PUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUSubTowerJet"), PUS_Tower);
  if(!PUS_Tower.isValid()){
    //evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUSubTowerJet") << std::endl;
  }

  // Local PU subtracted tower jets 
  edm::Handle<l1slhc::L1TowerJetCollection> LPUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("LocalPUSubTowerJet"), LPUS_Tower);
  if(!LPUS_Tower.isValid()){
    //evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("LocalPUSubTowerJet") << std::endl;
  }


  // Calibrated PrePUSak5PUS 
  edm::Handle<l1slhc::L1TowerJetCollection> CalibPrePUSak5PUS_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibratedPrePUSak5PUSTowerJet"), CalibPrePUSak5PUS_Tower);
  if(!CalibPrePUSak5PUS_Tower.isValid()){
    //    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibratedPrePUSak5PUSTowerJet") << std::endl;
  }

  edm::Handle<l1extra::L1JetParticleCollection> RecalibPrePUSak5PUS_L1Jet;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RecalibratedPrePUSak5PUSTowerJetL1Jet"), RecalibPrePUSak5PUS_L1Jet);
  if(!RecalibPrePUSak5PUS_L1Jet.isValid()){
    //    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RecalibratedPrePUSak5PUSTowerJetL1Jet") << std::endl;
  }



  // ************************************************************
  // *                        RECO Jets                         *
  // ************************************************************
  SUBPRINT("RECO Jets")

  // PU subtracted calibrated AK5 calo jets-must be in root file read in
  edm::Handle<reco::CaloJetCollection> PUSAk5Jets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUsubCaloJets"), PUSAk5Jets);
  if(!PUSAk5Jets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUsubCaloJets") << std::endl;
    //evValid = false;
  }

  // PrePU subtracted calibrated AK5 calo jets-must be in root file read in
  edm::Handle<reco::CaloJetCollection> PrePUSAk5Jets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSCaloJets"), PrePUSAk5Jets);
  if(!PrePUSAk5Jets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSCaloJets") << std::endl;
    //evValid = false;
  }






    // ******************************************************************** 
    // *                  Book Single event distributions                 * 
    // ******************************************************************** 
 

    int Eventnr   = iEvent.id().event(); 
    eventCount++; 


    // Check if the event is interesting
    bool interestingEvent = false;
    for (uint iInt = 0; iInt < eventsOfInterest.size(); ++iInt ){
    
      if ( Eventnr == eventsOfInterest[ iInt ] ){ 
	interestingEvent = true; 

	//Stop checking for this event again
	eventsOfInterest.erase ( eventsOfInterest.begin() + iInt);
	break;

      }
    


    }

 
//     // Create an event directory 
//     eventDir    = fs->mkdir( "Event" ); 

    if ( (eventCount <= maxEvents) || (interestingEvent) ){ // Process the first interesting events
 
      // **************************
      // *  TT histogram binning  *
      // **************************
      const Int_t TTetaBin = 56;
      const Int_t TTphiBin = 72;
      
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
      
      // TT phi binning
      double phiDiv = 2*PI/(TTphiBin);
      double TTphiBins[TTphiBin + 1];
      for (Int_t iPhi = 0; iPhi < TTphiBin + 1; iPhi++) {
	double phi = iPhi*phiDiv - PI;
	TTphiBins[iPhi] = phi;
      }




 
 
      TString eventStr = Form("%d", Eventnr ); 
      
      // Create an event directory 
      TFileDirectory eventDir    = fs->mkdir( "Event" ); 
      TFileDirectory eventSubDir = eventDir.mkdir( eventStr.Data() ); 

      // Tower inputs
      hist2D[eventStr + "_TT_EplusH-Phi_vs_Eta"]   = eventSubDir.make<TH2D>(eventStr + "_TT_EplusH-Phi_vs_Eta","TT Ecal + Hcal energy #phi vs #eta - Event " +
									    eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      hist2D[eventStr + "_Calo_EplusH-Phi_vs_Eta"] = eventSubDir.make<TH2D>(eventStr + "_Calo_EplusH-Phi_vs_Eta","TT Ecal + Hcal energy #phi vs #eta - Event " +
									    eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);

      // Reconstructed jets
      hist2D[eventStr + "_PrePUSAk5-Phi_vs_Eta"] = eventSubDir.make<TH2D>(eventStr + "_PrePUSAk5-Phi_vs_Eta","PrePUS Ak5 calojets #phi vs #eta - Event " +
									  eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
      hist2D[eventStr + "_PUSAk5-Phi_vs_Eta"]    = eventSubDir.make<TH2D>(eventStr + "_PUSAk5-Phi_vs_Eta","PUS Ak5 calojets #phi vs #eta - Event " +
									  eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);

      hist2D[eventStr + "_PrePUSL1-Phi_vs_Eta"] = eventSubDir.make<TH2D>(eventStr + "_PrePUSL1-Phi_vs_Eta","PrePUS L1 Towerjets #phi vs #eta - Event " +
									 eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);

      hist2D[eventStr + "_PUSL1-Phi_vs_Eta"] = eventSubDir.make<TH2D>(eventStr + "_PUSL1-Phi_vs_Eta","Globally PUS L1 Towerjets #phi vs #eta - Event " +
									 eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);

      hist2D[eventStr + "_LPUSL1-Phi_vs_Eta"] = eventSubDir.make<TH2D>(eventStr + "_LPUSL1-Phi_vs_Eta","Locally PUS L1 Towerjets #phi vs #eta - Event " +
									 eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);

      hist2D[eventStr + "_CurrL1-Phi_vs_Eta"]   = eventSubDir.make<TH2D>(eventStr + "_CurrL1-Phi_vs_Eta","Current L1 #phi vs #eta - Event " +
                                                                         eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);

      hist2D[eventStr + "_CalibPrePUSL1-Phi_vs_Eta"] = eventSubDir.make<TH2D>(eventStr + "_CalibPrePUSL1-Phi_vs_Eta","PrePUS L1 Towerjets #phi vs #eta - Event " +
									 eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);

      hist2D[eventStr + "_RecalibPrePUSL1-Phi_vs_Eta"] = eventSubDir.make<TH2D>(eventStr + "_RealibPrePUSL1-Phi_vs_Eta","Recalibrated PrePUS L1 Towerjets #phi vs #eta - Event " +
									 eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);


    }



    // L1 GCT ET
    //    double L1ET = l1CurrMET->at(0).etTotal();
    

    // ****************************************************************************************************
    // *                                       Calo distributions                                         *
    // ****************************************************************************************************








    // TT-RCT energy comparisons - Region strips
    std::vector<double> TTETRegion  (14, 0);
    std::vector<double> RCTETRegion (14, 0);
    std::vector<double> CaloETRegion (14, 0);

    // Multiplicity of the diffent calo objects in a region strip
    std::vector<int> TTRegionStripMultiplicity (14, 0);
    std::vector<int> RCTRegionStripMultiplicity (14, 0);
    std::vector<int> CaloRegionStripCaloAbove0p5Multiplicity (14, 0);


    // Region container [RCTiEta][RCTiPhi]
    std::vector< std::vector<double> > TTRegionEnergy   (14, std::vector<double> ( 18, 0 ) );
    std::vector< std::vector<double> > RCTRegionEnergy  (14, std::vector<double> ( 18, 0 ) );
    std::vector< std::vector<double> > CaloRegionEnergy (14, std::vector<double> ( 18, 0 ) );


    std::vector< std::vector<int> > TTRegionTTMultiplicity   (14, std::vector<int> ( 18, 0 ) );
  
    std::vector< std::vector<int> > CaloRegionCaloAbove0p5Multiplicity   (14, std::vector<int> ( 18, 0 ) );




    //    double RCTRegionETTotal(0);//, TTRegionETTotal(0), CaloRegionETTotal(0);
    int    RCTRegionFired(0),   TTRegionFired(0),   CaloRegionFired(0);


    //      // Calculate phi
    //      double phiDeg = RCTiPhi*20;
    //      if ( phiDeg > 160){ phiDeg -= 360; } // Get phi in degrees
    //      double phi = ( phiDeg*PI )/180.;     // Get phi in radians




    double RCTEt(0);
    
    for(unsigned int iRCT=0;iRCT < caloRegions->size(); ++iRCT ) {
      
      double RCTRegionET  = 0.5*caloRegions->at(iRCT).et();
      double RCTiEta      = caloRegions->at(iRCT).gctEta();
      double RCTiPhi      = caloRegions->at(iRCT).gctPhi();
      
      // Restrict to |eta| < 3
      if ( (RCTiEta < 4) || (RCTiEta > 17) ){ continue; }
      
      //	TString RCTiEtaStr = Form("%d", RCTiEta );
      int RCTindex = RCTiEta - 4;
      
      if ( RCTRegionET > 0 ){
	
	// Region
	RCTETRegion[ RCTindex ] += RCTRegionET;
	RCTRegionFired++;
	
	// Region strip
	RCTRegionEnergy[ RCTindex ][ RCTiPhi ] += RCTRegionET;
	RCTRegionStripMultiplicity[ RCTindex ]++; 
	
	// Total Et measured by RCT
	RCTEt += RCTRegionET;
	
      }
    }
    

    double TTEt(0);
    int TTsFired(0);
    
    for( l1slhc::L1CaloTowerCollection::const_iterator lTT_It = caloTowers->begin() ;
	 lTT_It != caloTowers->end() ; ++lTT_It ){
      
      // ****************************************
      // *    Load the calorimeter tower data   *
      // ****************************************
      double E      = 0.5*lTT_It->E();
      double H      = 0.5*lTT_It->H();
      int iEta      = lTT_It->iEta();
      int iPhi      = lTT_It->iPhi();
      
      // Restrict to central TTs
      if (abs(iEta) > 28)
	  continue;
     

      if ( E + H > 0 ){
	
	// Convert to RCT iEta
	int RCTiEtaIndex = iEta + 28;
	if ( iEta > 0 ){
	  RCTiEtaIndex--;
	}
	
	RCTiEtaIndex /= 4;
	
	// Convert to RCT iPhi
	int RCTiPhiIndex = (iPhi + 1)/4;
	if ( RCTiPhiIndex > 17 ){ RCTiPhiIndex = 0; }
	
	
	// Store TT region energy
	TTRegionEnergy[ RCTiEtaIndex ][ RCTiPhiIndex ] += E + H;
	TTRegionTTMultiplicity[ RCTiEtaIndex ][ RCTiPhiIndex ]++;
	
	// Region strip
	TTETRegion[ RCTiEtaIndex ] += E + H;
	TTRegionStripMultiplicity[RCTiEtaIndex]++;
	  
	  // Et according to TTs
	TTEt += E + H;
	TTsFired++;
	// 	  TTRegionETTotal += E + H;
	// 	  TTRegionFired++;
		 
	
      }
    }


      double caloEt(0);
      int caloAbove0p5Fired(0);
      
      for (CaloTowerCollection::const_iterator iCalo = caloFineTowers->begin(); iCalo != caloFineTowers->end(); ++iCalo) {
	
	// Restrict to central TTs
	if ( iCalo->ietaAbs() > 28 ){
	  continue;
	}
	
	// ****************************************
	// *    Load the calorimeter tower data   *
	// ****************************************
	double E   = iCalo->emEt();
	double H   = iCalo->hadEt();
	int iEta   = iCalo->ieta();
	int iPhi   = iCalo->iphi();
	
	if ( E + H > 0 ){
	  
	  // Convert to RCT iEta
	  int RCTiEtaIndex = iEta + 28;
	  if ( iEta > 0 ){
	    RCTiEtaIndex--;
	  }
	  
	  RCTiEtaIndex = RCTiEtaIndex/4;
	  
	  

	  // Convert to RCT iPhi
	  int RCTiPhiIndex = (iPhi + 1)/4;
	  if ( RCTiPhiIndex > 17 ){ RCTiPhiIndex = 0; }
	  //	std::cout << RCTiPhiIndex << " = " << iPhi << std::endl;
	  
	  
	  // Store calo region energy
	  CaloRegionEnergy[ RCTiEtaIndex ][ RCTiPhiIndex ] += E + H;
	  
	  //	  CaloRegionETTotal += E + H;
	  //	  CaloRegionFired++;
	  
	  caloEt += E + H;
	  
	  // Store calo region energy
	  CaloETRegion[ RCTiEtaIndex ] += E + H;
	  
	  
	  
	  if ( E + H > 0.5 ){
	    caloAbove0p5Fired++;
	    // Number of calotowers fired in region
	    CaloRegionCaloAbove0p5Multiplicity[ RCTiEtaIndex ][ RCTiPhiIndex ]++;
	    // Number of calotowers fired in region strip
	    CaloRegionStripCaloAbove0p5Multiplicity[ RCTiEtaIndex ]++;
	  }
	}
	
    }
      


      // Energy sum correlations
      hist2D["TT_vs_Calo_ET"]     ->Fill( caloEt, TTEt );
      hist1D["TT_vs_Calo_ET_prof"]->Fill( caloEt, TTEt );

      hist2D["RCT_vs_Calo_ET"]     ->Fill( caloEt, RCTEt );
      hist1D["RCT_vs_Calo_ET_prof"]->Fill( caloEt, RCTEt );

      hist2D["TT_vs_RCT_ET"]     ->Fill( RCTEt, TTEt );
      hist1D["TT_vs_RCT_ET_prof"]->Fill( RCTEt, TTEt );



      double deltaRCTTTEt       = RCTEt - TTEt;
      double deltaRCTTTEtplusRCTThresh = deltaRCTTTEt + 0.25*RCTRegionFired; // True RCT ET is higher by the energy loss through thresholding
//       std::cout << TTsFired <<  "   " << deltaRCTTTEt << "  " << deltaRCTTTEtplusRCTThresh << std::endl;
      
      hist2D["deltaET_TTRCT_vs_TTFired"]     ->Fill( TTsFired, deltaRCTTTEt );
      hist1D["deltaET_TTRCT_vs_TTFired_prof"]->Fill( TTsFired, deltaRCTTTEt );

      hist2D["deltaET_TTRCT_vs_Calo0p5Fired"]     ->Fill( caloAbove0p5Fired, deltaRCTTTEt );
      hist1D["deltaET_TTRCT_vs_Calo0p5Fired_prof"]->Fill( caloAbove0p5Fired, deltaRCTTTEt );


      hist2D["deltaET_TTRCTCorrected_vs_TTFired"]     ->Fill( TTsFired, deltaRCTTTEtplusRCTThresh );
      hist1D["deltaET_TTRCTCorrected_vs_TTFired_prof"]->Fill( TTsFired, deltaRCTTTEtplusRCTThresh );

      hist2D["deltaET_TTRCTCorrected_vs_Calo0p5Fired"]     ->Fill( caloAbove0p5Fired, deltaRCTTTEtplusRCTThresh );
      hist1D["deltaET_TTRCTCorrected_vs_Calo0p5Fired_prof"]->Fill( caloAbove0p5Fired, deltaRCTTTEtplusRCTThresh );




    for (int RCTiEta = 4; RCTiEta < 18; ++RCTiEta){

      int RCTindex = RCTiEta - 4;
      TString RCTiEtaStr = Form("%d", RCTiEta );


      // Region phi strips
      double RCTRegionStripET   = RCTETRegion[ RCTindex ];
      double TTRegionStripET    = TTETRegion[ RCTindex ];
      double CaloRegionStripET  = CaloETRegion[ RCTindex ];
      
      hist1D["RCT_RegionStripET_RCTiEta_"  + RCTiEtaStr]->Fill( RCTRegionStripET );
      hist1D["TT_RegionStripET_RCTiEta_"   + RCTiEtaStr]->Fill( TTRegionStripET );
      hist1D["Calo_RegionStripET_RCTiEta_" + RCTiEtaStr]->Fill( CaloRegionStripET );

      
      hist2D["RCT_vs_TT_ET_RCTiEta_" + RCTiEtaStr]  ->Fill( TTRegionStripET, RCTRegionStripET );

      hist2D["Calo_vs_TT_ET_RCTiEta_" + RCTiEtaStr] ->Fill( TTRegionStripET, CaloRegionStripET );
							
      hist2D["Calo_vs_RCT_ET_RCTiEta_" + RCTiEtaStr]->Fill( RCTRegionStripET, CaloRegionStripET );

      
      int TTsFiredInRegionStrip          = TTRegionStripMultiplicity[ RCTindex ];
      int RCTsFiredInRegionStrip         = RCTRegionStripMultiplicity[ RCTindex ];
      int caloAbove0p5FiredInRegionStrip = CaloRegionStripCaloAbove0p5Multiplicity[ RCTindex ];
      
      double deltaRCTTTRegionStripEt              = RCTRegionStripET - TTRegionStripET;
      double deltaRCTTTRegionStripEtplusRCTThresh = deltaRCTTTRegionStripEt + 0.25*RCTsFiredInRegionStrip; // True RCT ET is higher by the energy loss through thresholding   

      // deltaEt
      hist2D["deltaRegionStripET_TTRCT_vs_TTFired"]          ->Fill( TTsFiredInRegionStrip, deltaRCTTTRegionStripEt );
      hist1D["deltaRegionStripET_TTRCT_vs_TTFired_prof"]     ->Fill( TTsFiredInRegionStrip, deltaRCTTTRegionStripEt );

      hist2D["deltaRegionStripET_TTRCT_vs_Calo0p5Fired"]     ->Fill( caloAbove0p5FiredInRegionStrip, deltaRCTTTRegionStripEt );
      hist1D["deltaRegionStripET_TTRCT_vs_Calo0p5Fired_prof"]->Fill( caloAbove0p5FiredInRegionStrip, deltaRCTTTRegionStripEt );

      // Corrected deltaEt
      hist2D["deltaRegionStripETCorrected_TTRCT_vs_TTFired"]          ->Fill( TTsFiredInRegionStrip, deltaRCTTTRegionStripEtplusRCTThresh );
      hist1D["deltaRegionStripETCorrected_TTRCT_vs_TTFired_prof"]     ->Fill( TTsFiredInRegionStrip, deltaRCTTTRegionStripEtplusRCTThresh );

      hist2D["deltaRegionStripETCorrected_TTRCT_vs_Calo0p5Fired"]     ->Fill( caloAbove0p5FiredInRegionStrip, deltaRCTTTRegionStripEtplusRCTThresh );
      hist1D["deltaRegionStripETCorrected_TTRCT_vs_Calo0p5Fired_prof"]->Fill( caloAbove0p5FiredInRegionStrip, deltaRCTTTRegionStripEtplusRCTThresh );

      

//       // deltaEt - iEta binned
       hist2D["deltaRegionStripET_TTRCT_vs_TTFired_RCTiEta_" + RCTiEtaStr]          ->Fill( TTsFiredInRegionStrip, deltaRCTTTRegionStripEt );
       hist1D["deltaRegionStripET_TTRCT_vs_TTFired_RCTiEta_" + RCTiEtaStr + "_prof"]->Fill( TTsFiredInRegionStrip, deltaRCTTTRegionStripEt );

       hist2D["deltaRegionStripET_TTRCT_vs_Calo0p5Fired_RCTiEta_" + RCTiEtaStr]     ->Fill( caloAbove0p5FiredInRegionStrip, deltaRCTTTRegionStripEtplusRCTThresh );
       hist1D["deltaRegionStripET_TTRCT_vs_Calo0p5Fired_RCTiEta_" + RCTiEtaStr + "_prof"]->Fill( caloAbove0p5FiredInRegionStrip, deltaRCTTTRegionStripEtplusRCTThresh );
      




      if ( CaloRegionStripET > 0 ){
      
	double RCTETResponse = RCTRegionStripET/CaloRegionStripET;
	double TTETResponse  = TTRegionStripET/CaloRegionStripET;
	
	hist1D["RCT_ET_Response_RCTiEta_" + RCTiEtaStr]->Fill( RCTETResponse );
	hist1D["TT_ET_Response_RCTiEta_"  + RCTiEtaStr]->Fill( TTETResponse ); 

      }                 

      // NVTX correlations
      hist2D["RCT_E_vs_NVTX_RCTiEta_" + RCTiEtaStr] ->Fill( NVTX, RCTRegionStripET   );
      hist2D["TT_E_vs_NVTX_RCTiEta_" + RCTiEtaStr]  ->Fill( NVTX, TTRegionStripET    );
      hist2D["Calo_E_vs_NVTX_RCTiEta_" + RCTiEtaStr]->Fill( NVTX, CaloRegionStripET  );

      //      std::cout << RCTindex << "  " << RCTRegionStripET << "   " << TTRegionStripET << "  " << CaloRegionStripET << std::endl;

      for (int RCTiPhi = 0; RCTiPhi < 17; ++RCTiPhi){


	// Regions
	double RCTRegionET   = RCTRegionEnergy[ RCTindex ][ RCTiPhi ];
	double TTRegionET    = TTRegionEnergy[ RCTindex ][ RCTiPhi ];
	double CaloRegionET  = CaloRegionEnergy[ RCTindex ][ RCTiPhi ];

	int TTsFiredInRegion          = TTRegionTTMultiplicity[ RCTindex ][ RCTiPhi ];
	int caloAbove0p5FiredInRegion = CaloRegionCaloAbove0p5Multiplicity[ RCTindex ][ RCTiPhi ];


	double deltaRCTTTRegionEt       = RCTRegionET - TTRegionET;

	// Pseudo-region multiplicity
	if ( TTRegionET > 0){
	    TTRegionFired++;
	}
	if ( CaloRegionET > 0){
	    CaloRegionFired++;
	}


	//	double deltaRCTTTEtplusRCTThresh = deltaRCTTTEt - 0.25*TTRegionFired;
	//       std::cout << TTsFired <<  "   " << deltaRCTTTEt << "  " << deltaRCTTTEtplusRCTThresh << std::endl;

	hist2D["deltaRegionET_TTRCT_vs_TTFired"]     ->Fill( TTsFiredInRegion, deltaRCTTTRegionEt );
	hist1D["deltaRegionET_TTRCT_vs_TTFired_prof"]->Fill( TTsFiredInRegion, deltaRCTTTRegionEt );

	hist2D["deltaRegionET_TTRCT_vs_Calo0p5Fired"]     ->Fill( caloAbove0p5FiredInRegion, deltaRCTTTRegionEt );
	hist1D["deltaRegionET_TTRCT_vs_Calo0p5Fired_prof"]->Fill( caloAbove0p5FiredInRegion, deltaRCTTTRegionEt );



	hist1D["RCT_RegionET_RCTiEta_"  + RCTiEtaStr]->Fill( RCTRegionET );
	hist1D["TT_RegionET_RCTiEta_"   + RCTiEtaStr]->Fill( TTRegionET );
	hist1D["Calo_RegionET_RCTiEta_" + RCTiEtaStr]->Fill( CaloRegionET );

	hist1D["RCT_RegionET"]  ->Fill( RCTRegionET );
	hist1D["TT_RegionET"]   ->Fill( TTRegionET );
	hist1D["Calo_RegionET"] ->Fill( CaloRegionET );

	
	hist2D["RCT_vs_TT_RegionET_RCTiEta_" + RCTiEtaStr]  ->Fill( TTRegionET,  RCTRegionET );
	hist2D["Calo_vs_TT_RegionET_RCTiEta_" + RCTiEtaStr] ->Fill( TTRegionET,  CaloRegionET );
	hist2D["Calo_vs_RCT_RegionET_RCTiEta_" + RCTiEtaStr]->Fill( RCTRegionET, CaloRegionET );


	hist2D["RCT_vs_TT_RegionET"]  ->Fill( TTRegionET,  RCTRegionET );
	hist2D["Calo_vs_TT_RegionET"] ->Fill( TTRegionET,  CaloRegionET );
	hist2D["Calo_vs_RCT_RegionET"]->Fill( RCTRegionET, CaloRegionET );


	if ( CaloRegionStripET > 0 ){

	  double RCTRegionETResponse = RCTRegionET/CaloRegionET;
	  double TTRegionETResponse  = TTRegionET/CaloRegionET;

	  hist1D["RCT_RegionET_Response_RCTiEta_" + RCTiEtaStr]->Fill( RCTRegionETResponse );
	  hist1D["TT_RegionET_Response_RCTiEta_"  + RCTiEtaStr]->Fill( TTRegionETResponse ); 


	}


      }


    }



    hist2D["RCT_Etot_vs_NVTX"] ->Fill( NVTX, RCTEt );
    hist2D["RCT_Fired_vs_NVTX"]->Fill( NVTX, RCTRegionFired );
    // TT
    hist2D["TT_Etot_vs_NVTX"]  ->Fill( NVTX, TTEt );
    hist2D["TT_Fired_vs_NVTX"] ->Fill( NVTX, TTRegionFired );
    // Calo
    hist2D["Calo_Etot_vs_NVTX"] ->Fill( NVTX, caloEt );
    hist2D["Calo_Fired_vs_NVTX"]->Fill( NVTX, CaloRegionFired );






    PRINT("Calo distributions")
      
      TString TTPre = "Calo";
    
    int EFired(0), HFired(0), EorHFired(0), EandHFired(0);
    double maxE(0), maxH(0), totalE(0), totalH(0);
    
    // Energy sums
    TLorentzVector caloMETVec;
    //    double         caloET(0);
    
  
  for (CaloTowerCollection::const_iterator iCalo = caloFineTowers->begin(); iCalo != caloFineTowers->end(); ++iCalo) {

    // Restrict to central TTs                                                                                                                              
    if ( iCalo->ietaAbs() > 28 ){
      continue;
    }
    
    // ****************************************                                                                                                             
    // *    Load the calorimeter tower data   *                                                                                                             
    // ****************************************                                                                                                             
    double E      = iCalo->emEt();
    double H      = iCalo->hadEt();
    totalE += E;
    totalH += H;
    double EplusH = iCalo->et();
    int iEta   = iCalo->ieta();
    int iPhi   = iCalo->iphi();
    double Eta = mTowerGeo.eta(iEta);
    double Phi = mTowerGeo.phi(iPhi);

    if (Phi > PI)
      Phi -= 2*PI;

    // Energy sums
    TLorentzVector caloEtVec;
    caloEtVec.SetPtEtaPhiM( EplusH, Eta, Phi, 0 );

    // Calculate energy sums
    //    caloET     += EplusH;
    caloMETVec += caloEtVec;



    int EcalFired(0), HcalFired(0), EcalOrHcalFired(0), EcalAndHcalFired(0);                                                                            
    if (E > 0){
      EcalFired = 1;
      EFired++;
    }
    if (H > 0){
      HcalFired = 1;
      HFired++;
    }
    if ( (EcalFired + HcalFired) > 0 ){
      EcalOrHcalFired = 1;                                                                                                                          
      //      EorHFired++;
      if ( (EcalFired + HcalFired) == 2 ){
	EcalAndHcalFired = 1;                                                                                                                       
	//	EandHFired++;
      }
    }


    // Fill calo distributions
    hist2D[TTPre + "_EFired-Phi_vs_Eta"]    ->Fill( Eta, Phi, EcalFired );
    hist2D[TTPre + "_HFired-Phi_vs_Eta"]    ->Fill( Eta, Phi, HcalFired );
    hist2D[TTPre + "_EorHFired-Phi_vs_Eta"] ->Fill( Eta, Phi, EcalOrHcalFired );
    hist2D[TTPre + "_EandHFired-Phi_vs_Eta"]->Fill( Eta, Phi, EcalAndHcalFired );

    hist2D[TTPre + "_E-Phi_vs_Eta"]         ->Fill( Eta, Phi, E );
    hist2D[TTPre + "_H-Phi_vs_Eta"]         ->Fill( Eta, Phi, H );
    hist2D[TTPre + "_EplusH-Phi_vs_Eta"]    ->Fill( Eta, Phi, EplusH );




    // ******************************************************************** 
    // *                    Single event distributions                    * 
    // ******************************************************************** 

    // Make event distributions
    if ( (eventCount <= maxEvents) || (interestingEvent) ){

      TString eventStr = Form("%d", Eventnr );
      hist2D[eventStr + "_Calo_EplusH-Phi_vs_Eta"]->Fill( Eta, Phi, EplusH );
  
    }
    

  }

  hist1D[TTPre + "_Etotal"]     ->Fill( totalE );
  hist1D[TTPre + "_Htotal"]     ->Fill( totalH );
  hist1D[TTPre + "_EplusHtotal"]->Fill( totalE + totalH );

  //  double caloEt = totalE + totalH;

  // Calculate scalar projection of the MET vector
  double caloMET = caloMETVec.Pt();


    // ****************************************************************************************************
    // *                                        TT distributions                                          *
    // ****************************************************************************************************

    

    PRINT("TT distributions")
    TTPre = "TT";

    // Maximum event TT energies
    maxE = 0; maxH = 0; totalE = 0; totalH = 0;
    // Number of TTs that fired
    EFired = 0; HFired = 0; EorHFired = 0; EandHFired  = 0;

    int TTsBelow5GeV(0);

    // // Clear the array from the previous run
//     ttThresholdFired.clear();
//     ttThresholdFired.resize( ttThreshold.size() - 1 );





    // DEBUGGING - Correlation betwen TT and CaloTowers

    for (CaloTowerCollection::const_iterator iCalo = caloFineTowers->begin(); iCalo != caloFineTowers->end(); ++iCalo) {

      // Restrict to central TTs
      if ( iCalo->ietaAbs() > 28 ){
	continue;
      }

      // ****************************************
      // *    Load the calorimeter tower data   *
      // ****************************************
      double E      = iCalo->emEt();
      double H      = iCalo->hadEt();
      int iEta   = iCalo->ieta();
      int iPhi   = iCalo->iphi();

      for( l1slhc::L1CaloTowerCollection::const_iterator lTT_It = caloTowers->begin() ;
	   lTT_It != caloTowers->end() ; ++lTT_It ){

	// ****************************************
	// *    Load the calorimeter tower data   *
	// ****************************************
	double E2      = 0.5*lTT_It->E();
	double H2      = 0.5*lTT_It->H();
	int iEta2   = lTT_It->iEta();
	int iPhi2   = lTT_It->iPhi();

	double Eresponse = E2/E;
	double Hresponse = H2/H;

	// Restrict to central TTs
	if (abs(iEta2) > 28)
	  continue;

	if (( iEta == iEta2) && ( iPhi == iPhi2 ) ){
	  hist2D["TT_vs_Calo_E"]     ->Fill( E, E2 );
	  hist2D["TT_vs_Calo_H"]     ->Fill( H, H2 );
	  hist2D["TT_vs_Calo_EplusH"]->Fill( E + H, E2 + H2 );
	  hist2D["TT_over_Calo_E_vs_ieta"]->Fill( iEta, Eresponse );
	  hist2D["TT_over_Calo_H_vs_ieta"]->Fill( iEta, Hresponse );
	  hist1D["TT_over_Calo_E_vs_ieta_prof"]->Fill( iEta, Eresponse );
	  hist1D["TT_over_Calo_H_vs_ieta_prof"]->Fill( iEta, Hresponse );


	  // pT-binning
	  if ( E < 2){
	    hist2D["TT_over_Calo_E_vs_ieta_lt2GeV"]->Fill( iEta, Eresponse );
	  }
	  else{
	    hist2D["TT_over_Calo_E_vs_ieta_ge2GeV"]->Fill( iEta, Eresponse );
	  }
	  if ( H < 2){
	    hist2D["TT_over_Calo_H_vs_ieta_lt2GeV"]->Fill( iEta, Hresponse );
	  }
	  else{
	    hist2D["TT_over_Calo_H_vs_ieta_ge2GeV"]->Fill( iEta, Hresponse );
	  }


	  //	  std::cout << iEta << "\t" << iPhi << "\tEs: " << E << "\t" << E2 << "\tHs" << H << "\t" << H2 << "\n";
	  break;
	}

      }
      //      std::cout << "--------------------------------------------------------------------------------\n";

    }
    // DEBUGGING


    // Energy sums
    TLorentzVector L1METVec;
    double         L1TTET(0), L1TTET2GeV(0);


    for( l1slhc::L1CaloTowerCollection::const_iterator lTT_It = caloTowers->begin() ; 
	 lTT_It != caloTowers->end() ; ++lTT_It ){

      // ****************************************
      // *    Load the calorimeter tower data   *
      // ****************************************
      int E      = lTT_It->E();
      int H      = lTT_It->H();
      totalE += E;
      totalH += H;
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


      // Energy sums
      TLorentzVector ttEtVec;
      ttEtVec.SetPtEtaPhiM( double(EplusH)/2.0, Eta, Phi, 0 );
      
      // Calculate energy sums
      L1TTET2GeV += EplusH;
      L1METVec   += ttEtVec;

      
      // Store the maximum TT E and H
      if (E > maxE)
	maxE = E;
      if (H > maxH)
	maxH = H;

//       double EoverH(0);
//       if (H != 0){
// 	EoverH = double(E)/double(H);
//       }

//       if ( abs(iEta) > 20)
// 	std::cout << iPhi << "\t" << iEta << "\t" << EplusH << "\n";



      
      int EcalFired(0), HcalFired(0), EcalOrHcalFired(0), EcalAndHcalFired(0);
      if (E > 0){
	EcalFired = 1;
	EFired++;
      }
      if (H > 0){
	HcalFired = 1;
	HFired++;
      }
      if ( (EcalFired + HcalFired) > 0 ){
	EcalOrHcalFired = 1;
	EorHFired++;
	if ( (EcalFired + HcalFired) == 2 ){
	  EcalAndHcalFired = 1;
	  EandHFired++;
	}
      }


      if ( EplusH < 5)
	TTsBelow5GeV++;

//       hist2D[TTPre + "_Phi_vs_iPhi"]->Fill( iPhi, Phi );
//       hist2D[TTPre + "_Eta_vs_iEta"]->Fill( iEta, Eta ); 

      hist2D[TTPre + "_E-Phi_vs_Eta"]         ->Fill( Eta, Phi, E );
      hist2D[TTPre + "_H-Phi_vs_Eta"]         ->Fill( Eta, Phi, H );
      hist2D[TTPre + "_EplusH-Phi_vs_Eta"]    ->Fill( Eta, Phi, EplusH );

      //      hist2D[TTPre + "_EHratio-Phi_vs_Eta"]   ->Fill( Eta, Phi, EoverH );
//       hist2D[TTPre + "_EFG-Phi_vs_Eta"]       ->Fill( Eta, Phi, EcalFG );
//       hist2D[TTPre + "_HFG-Phi_vs_Eta"]       ->Fill( Eta, Phi, HcalFG );
      hist2D[TTPre + "_EFired-Phi_vs_Eta"]    ->Fill( Eta, Phi, EcalFired ); 
      hist2D[TTPre + "_HFired-Phi_vs_Eta"]    ->Fill( Eta, Phi, HcalFired ); 
      hist2D[TTPre + "_EorHFired-Phi_vs_Eta"] ->Fill( Eta, Phi, EcalOrHcalFired ); 
      hist2D[TTPre + "_EandHFired-Phi_vs_Eta"]->Fill( Eta, Phi, EcalAndHcalFired ); 


    // ******************************************************************** 
    // *                    Single event distributions                    * 
    // ******************************************************************** 

    // Make event distributions
      if ( (eventCount <= maxEvents) || (interestingEvent) ){

      TString eventStr = Form("%d", Eventnr );
      hist2D[eventStr + "_TT_EplusH-Phi_vs_Eta"]->Fill( Eta, Phi, EplusH );
  
    }
    



#ifdef TT_DIST

      hist2D[TTPre + "_E-iPhi_vs_iEta"]         ->Fill( iEta, iPhi, E );
      hist2D[TTPre + "_H-iPhi_vs_iEta"]         ->Fill( iEta, iPhi, H );
      hist2D[TTPre + "_EplusH-iPhi_vs_iEta"]    ->Fill( iEta, iPhi, EplusH );

      // Profile histograms
      hist1D[TTPre + "_E-iEta_prof"]         ->Fill( iEta, E );
      hist1D[TTPre + "_H-iEta_prof"]         ->Fill( iEta, H );
      hist1D[TTPre + "_EplusH-iEta_prof"]    ->Fill( iEta, EplusH );
      hist1D[TTPre + "_E-iPhi_prof"]         ->Fill( iPhi, E );
      hist1D[TTPre + "_H-iPhi_prof"]         ->Fill( iPhi, H );
      hist1D[TTPre + "_EplusH-iPhi_prof"]    ->Fill( iPhi, EplusH );

      //      hist2D[TTPre + "_EHratio-iPhi_vs_iEta"]   ->Fill( iEta, iPhi, EoverH ); 
      hist2D[TTPre + "_EFired-iPhi_vs_iEta"]    ->Fill( iEta, iPhi, EcalFired );
      hist2D[TTPre + "_HFired-iPhi_vs_iEta"]    ->Fill( iEta, iPhi, HcalFired ); 
      hist2D[TTPre + "_EorHFired-iPhi_vs_iEta"] ->Fill( iEta, iPhi, EcalOrHcalFired ); 
      hist2D[TTPre + "_EandHFired-iPhi_vs_iEta"]->Fill( iEta, iPhi, EcalAndHcalFired );  


      // iEta binned distributions
      TString iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEta)));
      hist1D[TTPre + "_E_" + iEtaStr]     ->Fill(E);
      hist1D[TTPre + "_H_" + iEtaStr]     ->Fill(H);
      hist1D[TTPre + "_EplusH_" + iEtaStr]->Fill(EplusH);


      //
      // CURRENTLY BROKEN
      //    
      bool ttBinned = false;
      // TT threshold plots
      for (unsigned int ttI = 0; ttI < ttThreshold.size(); ttI++){
	
	// get TT threshold
	int ttThresh   = ttThreshold[ttI];
	
	// TT energy is exceeded by the threshold
	if ( (EplusH <= ttThresh) && (ttI != 0) ){
	  
	  //	  std::cout << "E+H = " << EplusH << "\tThreshold = " << ttThresh << "\tCurrent = " << ttThresholdFired[ttI - 1] << "\n";
	  ttThresholdFired[ttI - 1]++;
	  ttBinned = true;
	  break;
  	}
	
      }
      if (ttBinned == false){
	ttThresholdFired[ttThresholdFired.size() - 1]++;
      }

#endif

    } // End TT loop


    // Correct for 0.5 LSB
    L1TTET = L1TTET2GeV/2.;


    hist1D[TTPre + "_Etotal"]     ->Fill( totalE );
    hist1D[TTPre + "_Htotal"]     ->Fill( totalH );
    hist1D[TTPre + "_EplusHtotal"]->Fill( totalE + totalH );
    
    //    double TTEt = totalE + totalH;
    
    // Calculate scalar projection of the MET vector
    double L1MET = L1METVec.Pt();





//     // Add turn ons!!!
//     // 30, 40, 50

//    std::cout << L1ET << "\t" << caloET << "\t\t\t" << L1MET << "\t" << caloMET << "\n";
    bool passesMetTrig = L1MET > 30;
    histEff[ "TT_MET_TurnOn_gt30" ]->Fill( passesMetTrig, caloMET );
    passesMetTrig = L1MET > 40;
    histEff[ "TT_MET_TurnOn_gt40" ]->Fill( passesMetTrig, caloMET );
    passesMetTrig = L1MET > 50;
    histEff[ "TT_MET_TurnOn_gt50" ]->Fill( passesMetTrig, caloMET );
    
    // Current L1 MET
    for( l1extra::L1EtMissParticleCollection::const_iterator lCurrMET_It = l1CurrMET->begin(); lCurrMET_It != l1CurrMET->end() ; ++lCurrMET_It ){

      double currL1MET = lCurrMET_It->etMiss();

      passesMetTrig = currL1MET > 30;
      histEff[ "CurrL1_MET_TurnOn_gt30" ]->Fill( passesMetTrig, caloMET );
      passesMetTrig = currL1MET > 40;
      histEff[ "CurrL1_MET_TurnOn_gt40" ]->Fill( passesMetTrig, caloMET );
      passesMetTrig = currL1MET > 50;
      histEff[ "CurrL1_MET_TurnOn_gt50" ]->Fill( passesMetTrig, caloMET );

    }


    bool passesEtTrig = L1TTET > 50;
    histEff[ "TT_ET_TurnOn_gt50" ] ->Fill( passesEtTrig, caloEt );
    passesEtTrig = L1TTET > 100;
    histEff[ "TT_ET_TurnOn_gt100" ]->Fill( passesEtTrig, caloEt );
    passesEtTrig = L1TTET > 150;
    histEff[ "TT_ET_TurnOn_gt150" ]->Fill( passesEtTrig, caloEt );










    // ******************************************************************** 
    // *                    Single event distributions                    * 
    // ******************************************************************** 

    // Make event distributions
    if ( (eventCount <= maxEvents) || (interestingEvent) ){

      TString eventStr = Form("%d", Eventnr );

      // Create an event directory 
      TFileDirectory eventDir    = fs->mkdir( "Event" );
      TFileDirectory eventSubDir = eventDir.mkdir( eventStr.Data() );
      TFileDirectory eventJetDir  = eventSubDir.mkdir( "Jet" );
      TFileDirectory eventTTDir   = eventSubDir.mkdir( "TT" );
      TFileDirectory eventCaloDir = eventSubDir.mkdir( "Calo" );






      // ****************************************
      // *             Visualisation            *
      // ****************************************                                                                                                                     











      // ******************************************************************************** 
      // *                                 Canvas jets                                  * 
      // ******************************************************************************** 

      // Eta cut on CanvasJets
      double etaCut = 3.0;

      //                Ak5
      // ****************************************

      // L2L3 PrePUS corrected jet collection 
      CanvasJetCollection PrePUSak5Cone = fillCanvasJetCollection( PrePUSAk5Jets, etaCut, hist2D[eventStr + "_PrePUSAk5-Phi_vs_Eta"] );
      // L2L3 PUS corrected jet collection 
      CanvasJetCollection PUSak5Cone    = fillCanvasJetCollection( PUSAk5Jets,    etaCut, hist2D[eventStr + "_PUSAk5-Phi_vs_Eta"] );

      //                L1
      // ****************************************
   
      // Pre PU subtraction 
      CanvasJetCollection PrePUSL1Cone = fillCanvasJetCollection( PrePUS_Tower, etaCut, hist2D[eventStr + "_PrePUSL1-Phi_vs_Eta"] );

      // Global PU subtraction
      CanvasJetCollection PUSL1Cone = fillCanvasJetCollection( PUS_Tower, etaCut, hist2D[eventStr + "_PUSL1-Phi_vs_Eta"] );

      // Local PU subtraction
      CanvasJetCollection LPUSL1Cone = fillCanvasJetCollection( LPUS_Tower, etaCut, hist2D[eventStr + "_LPUSL1-Phi_vs_Eta"] );
      
      // Current L1 GCT jets
      CanvasJetCollection CurrL1Cone   = fillCanvasJetCollection( l1extraparticles, iEvent, etaCut, hist2D[eventStr + "_CurrL1-Phi_vs_Eta"] );



      // Calibrated Pre PU subtraction 
      CanvasJetCollection CalibPrePUSL1Cone = fillCanvasJetCollection( CalibPrePUSak5PUS_Tower, etaCut, hist2D[eventStr + "_CalibPrePUSL1-Phi_vs_Eta"] );


      CanvasJetCollection RecalibPrePUSL1Cone = fillCanvasJetCollection( RecalibPrePUSak5PUS_L1Jet, etaCut, hist2D[eventStr + "_RecalibPrePUSL1-Phi_vs_Eta"] );

      // ************************************************************
      // *                   Rank and colour jets                   *
      // ************************************************************

      rankColourCanvasJetCollection( PrePUSak5Cone );
      rankColourCanvasJetCollection( PUSak5Cone    );
      rankColourCanvasJetCollection( PrePUSL1Cone  );
      rankColourCanvasJetCollection( CalibPrePUSL1Cone  );
      rankColourCanvasJetCollection( RecalibPrePUSL1Cone  );
      rankColourCanvasJetCollection( PUSL1Cone   );
      rankColourCanvasJetCollection( LPUSL1Cone  );
      rankColourCanvasJetCollection( CurrL1Cone  );

      // ************************************************************
      // *                           Jets                           *
      // ************************************************************


      // L1
      // ****************************************

      // PrePUS
      TCanvas* PrePUSL1Canv = eventJetDir.make<TCanvas>(eventStr + "_PrePUSL1", "", 800,600);
      hist2D[eventStr + "_PrePUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( PrePUSL1Cone );
      PrePUSL1Canv->Update();

      // CalibPrePUS
      TCanvas* CalibPrePUSL1Canv = eventJetDir.make<TCanvas>(eventStr + "_CalibPrePUSL1", "", 800,600);
      hist2D[eventStr + "_CalibPrePUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( CalibPrePUSL1Cone );
      CalibPrePUSL1Canv->Update();

      // RecalibPrePUS
      TCanvas* RecalibPrePUSL1Canv = eventJetDir.make<TCanvas>(eventStr + "_RecalibPrePUSL1", "", 800,600);
      hist2D[eventStr + "_RecalibPrePUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( RecalibPrePUSL1Cone );
      RecalibPrePUSL1Canv->Update();

      // PUS
      TCanvas* PUSL1Canv = eventJetDir.make<TCanvas>(eventStr + "_PUSL1", "", 800,600);
      hist2D[eventStr + "_PUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( PUSL1Cone );
      PUSL1Canv->Update();

//       // LPUS
//       TCanvas* LPUSL1Canv = eventJetDir.make<TCanvas>(eventStr + "_LPUSL1", "", 800,600);
//       hist2D[eventStr + "_LPUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
//       drawCanvasJetCollection( LPUSL1Cone );
//       LPUSL1Canv->Update();

      // Current L1
      TCanvas* CurrL1Canv = eventJetDir.make<TCanvas>(eventStr + "_CurrL1", "", 800,600);
      hist2D[eventStr + "_CurrL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( CurrL1Cone );
      CurrL1Canv->Update();

//       // L1 Compared to ak5PrePUS
//       TCanvas* PrePUSL1PrePUSAk5Canv = eventJetDir.make<TCanvas>(eventStr + "_PrePUSL1PrePUSAk5", "", 800,600);
//       hist2D[eventStr + "_PrePUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
//       drawCanvasJetCollection( PrePUSak5Cone );
//       PrePUSL1PrePUSAk5Canv->Update();

      // L1 Compared to ak5PUS
      TCanvas* PrePUSL1PUSAk5Canv = eventJetDir.make<TCanvas>(eventStr + "_PrePUSL1PUSAk5", "", 800,600);
      hist2D[eventStr + "_PrePUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( PUSak5Cone );
      PrePUSL1PUSAk5Canv->Update();

      TCanvas* CalibPrePUSL1PUSAk5Canv = eventJetDir.make<TCanvas>(eventStr + "_CalibPrePUSL1PUSAk5", "", 800,600);
      hist2D[eventStr + "_CalibPrePUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( PUSak5Cone );
      CalibPrePUSL1PUSAk5Canv->Update();

      TCanvas* PUSL1PUSAk5Canv = eventJetDir.make<TCanvas>(eventStr + "_PUSL1PUSAk5", "", 800,600);
      hist2D[eventStr + "_PUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( PUSak5Cone );
      PUSL1PUSAk5Canv->Update();

//       TCanvas* LPUSL1PUSAk5Canv = eventJetDir.make<TCanvas>(eventStr + "_LPUSL1PUSAk5", "", 800,600);
//       hist2D[eventStr + "_LPUSL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
//       drawCanvasJetCollection( PUSak5Cone );
//       LPUSL1PUSAk5Canv->Update();

      TCanvas* CurrL1PUSAk5Canv = eventJetDir.make<TCanvas>(eventStr + "_CurrL1PUSAk5", "", 800,600);
      hist2D[eventStr + "_CurrL1-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( PUSak5Cone );
      CurrL1PUSAk5Canv->Update();

      // RECO
      // ****************************************

//       // PrePUS
//       TCanvas* PrePUSAk5Canv = eventJetDir.make<TCanvas>(eventStr + "_PrePUSAk5", "", 800,600);
//       hist2D[eventStr + "_PrePUSAk5-Phi_vs_Eta"]->Draw("COLZ TEXT");
//       drawCanvasJetCollection( PrePUSak5Cone );
//       PrePUSAk5Canv->Update();

      TCanvas* PrePUSAk5PrePUSL1Canv = eventJetDir.make<TCanvas>(eventStr + "_PrePUSAk5PrePUSL1", "", 800,600);
      hist2D[eventStr + "_PrePUSAk5-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( PrePUSL1Cone );
      PrePUSAk5PrePUSL1Canv->Update();

      // PUS
      TCanvas* PUSAk5Canv = eventJetDir.make<TCanvas>(eventStr + "_PUSAk5", "", 800,600);
      hist2D[eventStr + "_PUSAk5-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( PUSak5Cone );
      PUSAk5Canv->Update();

      TCanvas* PUSAk5PrePUSL1Canv = eventJetDir.make<TCanvas>(eventStr + "_PUSAk5PrePUSL1", "", 800,600);
      hist2D[eventStr + "_PUSAk5-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( PrePUSL1Cone );
      PUSAk5PrePUSL1Canv->Update();

      TCanvas* PUSAk5RecalibPrePUSL1Canv = eventJetDir.make<TCanvas>(eventStr + "_PUSAk5RecalibPrePUSL1", "", 800,600);
      hist2D[eventStr + "_PUSAk5-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( RecalibPrePUSL1Cone );
      PUSAk5RecalibPrePUSL1Canv->Update();


      TCanvas* PUSAk5CurrL1Canv = eventJetDir.make<TCanvas>(eventStr + "_PUSAk5CurrL1", "", 800,600);
      hist2D[eventStr + "_PUSAk5-Phi_vs_Eta"]->Draw("COLZ TEXT");
      drawCanvasJetCollection( CurrL1Cone );
      PUSAk5CurrL1Canv->Update();

      // ************************************************************
      // *                      Trigger towers                      *
      // ************************************************************

      TCanvas* TTCanv   = eventTTDir.make<TCanvas>(eventStr + "_TT", "", 800,600);    
      hist2D[eventStr + "_TT_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      TTCanv->Update();

//       TCanvas* TTPrePUSL1Canv   = eventTTDir.make<TCanvas>(eventStr + "_TTL1Jets", "", 800,600);    
//       hist2D[eventStr + "_TT_EplusH-Phi_vs_Eta"]->Draw("COLZ");
//       drawCanvasJetCollection( PrePUSL1Cone );
//       TTPrePUSL1Canv->Update();

      TCanvas* TTPrePUSL1PTCanv   = eventTTDir.make<TCanvas>(eventStr + "_TTPrePUSL1PT", "", 800,600);    
      hist2D[eventStr + "_TT_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      drawCanvasJetCollection( PrePUSL1Cone );
      hist2D[eventStr + "_PrePUSL1-Phi_vs_Eta"]->Draw("SAME TEXT");
      TTPrePUSL1PTCanv->Update();

      TCanvas* TTRecalibPrePUSL1PTCanv   = eventTTDir.make<TCanvas>(eventStr + "_TTRecalibPrePUSL1PT", "", 800,600);    
      hist2D[eventStr + "_TT_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      drawCanvasJetCollection( RecalibPrePUSL1Cone );
      hist2D[eventStr + "_RecalibPrePUSL1-Phi_vs_Eta"]->Draw("SAME TEXT");
      TTRecalibPrePUSL1PTCanv->Update();

      TCanvas* TTCurrL1PTCanv   = eventTTDir.make<TCanvas>(eventStr + "_TTCurrL1PT", "", 800,600);    
      hist2D[eventStr + "_TT_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      drawCanvasJetCollection( CurrL1Cone );
      hist2D[eventStr + "_CurrL1-Phi_vs_Eta"]->Draw("SAME TEXT");
      TTCurrL1PTCanv->Update();

      TCanvas* TTPrePUSAk5PTCanv   = eventTTDir.make<TCanvas>(eventStr + "_TTPrePUSAk5PT", "", 800,600);    
      hist2D[eventStr + "_TT_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      drawCanvasJetCollection( PrePUSak5Cone );
      hist2D[eventStr + "_PrePUSAk5-Phi_vs_Eta"]->Draw("SAME TEXT");
      TTPrePUSAk5PTCanv->Update();

      TCanvas* TTPUSAk5PTCanv   = eventTTDir.make<TCanvas>(eventStr + "_TTPUSAk5PT", "", 800,600);    
      hist2D[eventStr + "_TT_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      drawCanvasJetCollection( PUSak5Cone );
      hist2D[eventStr + "_PUSAk5-Phi_vs_Eta"]->Draw("SAME TEXT");
      TTPUSAk5PTCanv->Update();


      // ************************************************************
      // *                        Calo towers                       *
      // ************************************************************
      TCanvas* CaloCanv   = eventCaloDir.make<TCanvas>(eventStr + "_Calo", "", 800,600);    
      hist2D[eventStr + "_Calo_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      CaloCanv->Update();


//       TCanvas* CaloPrePUSAk5Canv   = eventCaloDir.make<TCanvas>(eventStr + "_CaloPrePUSAk5", "", 800,600);    
//       hist2D[eventStr + "_Calo_EplusH-Phi_vs_Eta"]->Draw("COLZ");
//       drawCanvasJetCollection( PrePUSak5Cone );
//       CaloPrePUSAk5Canv->Update();

      TCanvas* CaloPrePUSAk5PTCanv   = eventCaloDir.make<TCanvas>(eventStr + "_CaloPrePUSAk5PT", "", 800,600);    
      hist2D[eventStr + "_Calo_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      drawCanvasJetCollection( PrePUSak5Cone );
      hist2D[eventStr + "_PrePUSAk5-Phi_vs_Eta"]->Draw("SAME TEXT");
      CaloPrePUSAk5PTCanv->Update();

      TCanvas* CaloPUSAk5PTCanv   = eventCaloDir.make<TCanvas>(eventStr + "_CaloPUSAk5PT", "", 800,600);    
      hist2D[eventStr + "_Calo_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      drawCanvasJetCollection( PUSak5Cone );
      hist2D[eventStr + "_PUSAk5-Phi_vs_Eta"]->Draw("SAME TEXT");
      CaloPUSAk5PTCanv->Update();

      TCanvas* CaloPrePUSL1PTCanv   = eventCaloDir.make<TCanvas>(eventStr + "_CaloPrePUSL1PT", "", 800,600);    
      hist2D[eventStr + "_Calo_EplusH-Phi_vs_Eta"]->Draw("COLZ");
      drawCanvasJetCollection( PrePUSL1Cone );
      hist2D[eventStr + "_PrePUSL1-Phi_vs_Eta"]->Draw("SAME TEXT");
      CaloPrePUSL1PTCanv->Update();
  


      if (eventCount == maxEvents){
	std::cout << "\n\nMade all event distributions\n\n";
      }

    }








	
#ifdef TT_DIST

    //    std::cout << "\nNEW EVENT\n";

    // Fill the maximum TT energies in the event
    hist1D[TTPre + "_MaxE"]->Fill(maxE);
    hist1D[TTPre + "_MaxH"]->Fill(maxH);
    // Fill TT multiplicities
    hist1D[TTPre + "_ETowersFired"]    ->Fill( EFired );
    hist1D[TTPre + "_HTowersFired"]    ->Fill( HFired );
    hist1D[TTPre + "_EorHTowersFired"] ->Fill( EorHFired );
    hist1D[TTPre + "_EandHTowersFired"]->Fill( EandHFired );

    
    // ************************************************************
    // *                       RECO kt6 rho                       *
    // ************************************************************
    hist2D["GlobalRho_vs_kt6Rho"]  ->Fill( *kt6CaloRho, *GlobalRho_Tower );
    hist2D["EFired_vs_kt6Rho"]     ->Fill( *kt6CaloRho, EFired );
    hist2D["HFired_vs_kt6Rho"]     ->Fill( *kt6CaloRho, HFired );
    hist2D["TTFired_vs_kt6Rho"]    ->Fill( *kt6CaloRho, EorHFired );
    hist2D["Etotal_vs_kt6Rho"]     ->Fill( *kt6CaloRho, totalE );
    hist2D["Htotal_vs_kt6Rho"]     ->Fill( *kt6CaloRho, totalH );
    hist2D["EplusHtotal_vs_kt6Rho"]->Fill( *kt6CaloRho, totalE + totalH );

    // Inter-quantity correlations
    hist2D["EplusHtotal_vs_TTFired"]->Fill( EorHFired, totalE + totalH );
    hist2D["EFired_vs_HFired"]      ->Fill( HFired, EFired );
    hist2D["ETotal_vs_EFired"]      ->Fill( EFired, totalE );
    hist2D["HTotal_vs_HFired"]      ->Fill( HFired, totalH );
//     hist2D["NJets_vs_EplusHtotal"]  ->Fill( totalE + totalH, upgradeJets.size() );
//     hist2D["NJets_vs_TTFired"]      ->Fill( EorHFired, upgradeJets.size() );
//  hist1D["NJets_vs_NVTX_prof"]

    hist2D["TTFired_vs_NVTX"]         ->Fill( NVTX, EorHFired );
    hist2D["TTFiredBelow5GeV_vs_NVTX"]->Fill( NVTX, TTsBelow5GeV );

    // Profiles
    hist1D["GlobalRho_vs_kt6Rho_prof"]  ->Fill( *kt6CaloRho, *GlobalRho_Tower );
    hist1D["EFired_vs_kt6Rho_prof"]     ->Fill( *kt6CaloRho, EFired );
    hist1D["HFired_vs_kt6Rho_prof"]     ->Fill( *kt6CaloRho, HFired );
    hist1D["TTFired_vs_kt6Rho_prof"]    ->Fill( *kt6CaloRho, EorHFired );
    hist1D["Etotal_vs_kt6Rho_prof"]     ->Fill( *kt6CaloRho, totalE );
    hist1D["Htotal_vs_kt6Rho_prof"]     ->Fill( *kt6CaloRho, totalH );
    hist1D["EplusHtotal_vs_kt6Rho_prof"]->Fill( *kt6CaloRho, totalE + totalH );
    hist1D["NVTX_vs_kt6Rho_prof"]       ->Fill( *kt6CaloRho, NVTX );   

    hist1D["GlobalRho_vs_NVTX_prof"]    ->Fill( NVTX, *GlobalRho_Tower );   
    hist1D["EFired_vs_NVTX_prof"]       ->Fill( NVTX, EFired );      
    hist1D["HFired_vs_NVTX_prof"]       ->Fill( NVTX, HFired );      
    hist1D["TTFired_vs_NVTX_prof"]      ->Fill( NVTX, EorHFired );
    hist1D["Etotal_vs_NVTX_prof"]       ->Fill( NVTX, totalE );      
    hist1D["Htotal_vs_NVTX_prof"]       ->Fill( NVTX, totalH );      
    hist1D["EplusHtotal_vs_NVTX_prof"]  ->Fill( NVTX, totalE + totalH ); 


    // Inter-quantity correlations
    hist1D["EplusHtotal_vs_TTFired_prof"]->Fill( EorHFired, totalE + totalH );
    hist1D["EFired_vs_HFired_prof"]      ->Fill( HFired, EFired );
    hist1D["ETotal_vs_EFired_prof"]      ->Fill( EFired, totalE );
    hist1D["HTotal_vs_HFired_prof"]      ->Fill( HFired, totalH );
//     hist1D["NJets_vs_EplusHtotal_prof"]  ->Fill( totalE + totalH, upgradeJets.size() );
//     hist1D["NJets_vs_TTFired_prof"]      ->Fill( EorHFired, upgradeJets.size() );

 
    // Fill TT threshold plots
    for (unsigned int ttI = 0; ttI < ttThreshold.size(); ttI++){
      
      // get TT threshold
      int ttThresh   = ttThreshold[ttI];  
      TString ttThreshStr = Form("%d",Int_t(ttThresh));
      int ttMultiplicity = ttThresholdFired[ttI];
      hist2D["TTFired_vs_kt6Rho_TT_gt_" + ttThreshStr]          ->Fill(*kt6CaloRho, ttMultiplicity );
      hist1D["TTFired_vs_kt6Rho_TT_gt_" + ttThreshStr + "_prof"]->Fill(*kt6CaloRho, ttMultiplicity );

    }

#endif










  QUIT("Successfully completed analysis of single event")

} // end analyser




// ------------ method called once each job just after ending the event loop  ------------
void TTHist::endJob(){



}

// ------------ method called when starting to processes a run  ------------
void
TTHist::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
TTHist::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TTHist::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TTHist::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTHist::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}







// ------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------------------------------------



CanvasJet* makeJetCone( double eta, double phi, double pT, double R, double area, int color, bool transparent){

  CanvasJet* jet = new CanvasJet( new TEllipse( eta, phi, R), pT, area );

  // Set fill attributes
  jet->jetCone->SetFillColor( color );
  if (transparent){
    jet->jetCone->SetFillStyle(4000);
  }
  else{
    jet->jetCone->SetFillStyle(3002);
  }
  
  // Set line attributes
  jet->jetCone->SetLineColor( color );
  jet->jetCone->SetLineWidth( 2 );

  return jet;

}







//define this as a plug-in
DEFINE_FWK_MODULE(TTHist);

//  LocalWords:  RCTiEtaStr
