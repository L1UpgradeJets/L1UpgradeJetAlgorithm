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





// ========================================
// class declaration
// ========================================

class TTGCTHist : public edm::EDAnalyzer {

public:
  explicit TTGCTHist(const edm::ParameterSet&);
  ~TTGCTHist();
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

  Service<TFileService> fs;

};





TTGCTHist::TTGCTHist(const edm::ParameterSet& iConfig): conf_(iConfig){

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

  TFileDirectory TTDir = fs->mkdir( "TTGCTHist" );




    // ************************************************************
    // *                  Correlation plots                       *
    // ************************************************************


  // CaloTower-TT correlations    
  TFileDirectory TTCorr = TTDir.mkdir( "TTCorr" );
  hist2D["TT_vs_Calo_E"]       = TTCorr.make<TH2D>("TT_vs_Calo_E",      "TT vs CaloTower E;Calo E;TT E", 201,-0.05,20.05, 21,-0.5,20.5 );
  hist2D["TT_vs_Calo_H"]       = TTCorr.make<TH2D>("TT_vs_Calo_H",      "TT vs CaloTower H;Calo H;TT H", 201,-0.05,20.05, 21,-0.5,20.5 );
  hist2D["TT_vs_Calo_EplusH"]  = TTCorr.make<TH2D>("TT_vs_Calo_EplusH", "TT vs CaloTower E+H;Calo E+H;TT E+H", 401,-0.05,40.05, 41,-0.5,40.5 );


  hist1D["TT_over_Calo_E_vs_ieta"] = TTCorr.make<TProfile>("TT_over_Calo_E_vs_ieta", "TT/CaloTower E vs i#eta;i#eta;TT/Calo E", 57,-28.5,28.5, 0,4 );
  hist1D["TT_over_Calo_H_vs_ieta"] = TTCorr.make<TProfile>("TT_over_Calo_H_vs_ieta", "TT/CaloTower H vs i#eta;i#eta;TT/Calo H", 57,-28.5,28.5, 0,4 );


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



TTGCTHist::~TTGCTHist()
{
}



// ------------ method called once each job just before starting event loop  ------------
void TTGCTHist::beginJob(){





}








// **********************************************************************
// *                              analyze()                             *
// *                                                                    *
// *                                                                    *
// **********************************************************************


void TTGCTHist::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){


  //  bool evValid = true;

//   //Need this for information about PU
//   edm::Handle<int> vertices;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("NumPrimaryVertices"), vertices);
//   if(!vertices.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("NumPrimaryVertices") << std::endl;
//   }
//   NVTX       = *vertices;


  // RCT regions 
  edm::Handle<L1CaloRegionCollection> caloRegions;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CaloRegions"), caloRegions );
   if(!caloRegions.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CaloRegions") << std::endl;
     //     RCTData = false; 

   } 








//     // Current MET and ET
//     SUBPRINT("Current L1 MET")
//       edm::Handle< l1extra::L1EtMissParticleCollection > l1CurrMET;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CurrentL1MET"), l1CurrMET );
//   if(!l1CurrMET.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CurrentL1MET") << std::endl;
//     //    evValid = false;
//   }


//     // Current MHT and HT
//     SUBPRINT("Current L1 MHT")
//     edm::Handle< std::vector<l1extra::L1EtMissParticle> > l1CurrMHT;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CurrentL1MHT"), l1CurrMHT );
//   if(!l1CurrMHT.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CurrentL1MHT") << std::endl;
//     //    evValid = false;
//   }


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





    // ******************************************************************** 
    // *                  Book Single event distributions                 * 
    // ******************************************************************** 
 

  //    int Eventnr   = iEvent.id().event(); 
    eventCount++; 

    

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

//       // NVTX correlations
//       hist2D["RCT_E_vs_NVTX_RCTiEta_" + RCTiEtaStr] ->Fill( NVTX, RCTRegionStripET   );
//       hist2D["TT_E_vs_NVTX_RCTiEta_" + RCTiEtaStr]  ->Fill( NVTX, TTRegionStripET    );
//       hist2D["Calo_E_vs_NVTX_RCTiEta_" + RCTiEtaStr]->Fill( NVTX, CaloRegionStripET  );

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



//     hist2D["RCT_Etot_vs_NVTX"] ->Fill( NVTX, RCTEt );
//     hist2D["RCT_Fired_vs_NVTX"]->Fill( NVTX, RCTRegionFired );
//     // TT
//     hist2D["TT_Etot_vs_NVTX"]  ->Fill( NVTX, TTEt );
//     hist2D["TT_Fired_vs_NVTX"] ->Fill( NVTX, TTRegionFired );
//     // Calo
//     hist2D["Calo_Etot_vs_NVTX"] ->Fill( NVTX, caloEt );
//     hist2D["Calo_Fired_vs_NVTX"]->Fill( NVTX, CaloRegionFired );








} // end analyser




// ------------ method called once each job just after ending the event loop  ------------
void TTGCTHist::endJob(){



}

// ------------ method called when starting to processes a run  ------------
void
TTGCTHist::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
TTGCTHist::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TTGCTHist::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TTGCTHist::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTGCTHist::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}







// ------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------------------------------------





//define this as a plug-in
DEFINE_FWK_MODULE(TTGCTHist);

//  LocalWords:  RCTiEtaStr
