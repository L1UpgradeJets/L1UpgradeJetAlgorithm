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
//#include "AnalyseUpgradeJets/AnalyzeJets/src/JetMatch.h"
#include "AnalyseUpgradeJets/AnalyzeJets/src/helperFunctions.cc"


#include <algorithm>  // for sorting

using namespace l1slhc;
using namespace edm;
using namespace std;
using namespace reco;




// ************************************************************ 
// *    Histogram switches - Enable histograms with defines   *
// ************************************************************

// Make upgrade jet distributions
#define UPGRADE
// Make current jet distributions
#define CURRENT

// ****************************************
//             Validation
// ****************************************

// Make vertex distributions
#define VERTEX
// // Make TT-level distributions
#define TT_DIST
// // Make all candidate jet distributions
// #define PROTO_JET
// #define FILT_1D
// #define FILT_CENT

// ****************************************

// Print debugging messages
//#define VERBOSE


#ifdef VERBOSE
#    define PRINT(outputStr) std::cout << "\n************************************************************\n" << "Making: " << (outputStr) << "\n" << "******\
******************************************************\n\n";
#    define SUBPRINT(outputStr) std::cout << "\t" << outputStr << "\n";
#    define QUIT(outputStr) std::cout << "\n\n\n" << outputStr << "\n\n\n"; exit(0);
#else
#    define PRINT(outputStr)
#    define SUBPRINT(outputStr)
#    define QUIT(outputStr)
#endif









// ********************************************************************************
// *                           Configuration parameters                           *
// ********************************************************************************

const double PI = 3.141592654;
// eta boundary of the barrel region
const double CEN_BOUNDARY = 1.3050;


// jet cone size
const double deltaR = 0.7;






// compares two TLorentzVectors, returns true if the first vector has a larger pT than the second
// Used for storting jets by Pt
bool sortTLorentz (TLorentzVector i,TLorentzVector j) { return ( i.Pt() > j.Pt() ); }


// return the azimuthal angle, constrained to the range (-PI,PI)
double getAzimuth(double angle){return angle - PI*Int_t(angle/PI);}

// Extract the median value of an array
//extern double Median( vector<double> aVec);













// ========================================
// class declaration
// ========================================

class AnalyzeJets : public edm::EDAnalyzer {

public:
  explicit AnalyzeJets(const edm::ParameterSet&);
  ~AnalyzeJets();
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
  std::pair< double, int> Leading_Match(vector<TLorentzVector> offlineJets,vector<TLorentzVector> L1Jets );



  // Fill offline-matched lead L1 jet histograms
  void fillL1OfflineHistograms( vector<TLorentzVector> L1Jets, vector<TLorentzVector> OfflineJets, TString prefix, Int_t NVTX );
  // Fill L1 all-jet histograms
  void fillL1Histograms( vector<TLorentzVector> L1Jets, TString prefix, Int_t NVTX );

  // Fill Tower jet correlation and dead jet histograms
  void fillL1CorrelationHistograms( edm::Handle<l1slhc::L1TowerJetCollection> const& JetColl1, edm::Handle<l1slhc::L1TowerJetCollection> const& JetColl2, 
				    TString prefix, Int_t NVTX );


  // Fill TowerJet histograms
  void fillTowerJetHistograms(edm::Handle<l1slhc::L1TowerJetCollection> const&, TString jetPrefix);

//   void matchOnlineOfflineJets( edm::Handle<l1slhc::L1TowerJetCollection> const&, edm::Handle<reco::CaloJetCollection> const&,
// 			       TString prefix, Int_t NVTX );

//   void matchOnlineOfflineJets( edm::Handle<l1slhc::L1TowerJetCollection> const&, edm::Handle<edm::View< reco::CaloJet > > const&, 
// 			       TString prefix, Int_t NVTX );

//   void matchOnlineOfflineJets( edm::Handle<l1slhc::L1TowerJetCollection> const&, edm::Handle<edm::View< reco::CaloJet > > const&, 
// 			       TString prefix, Int_t NVTX );
//   void matchOnlineOfflineJets( l1slhc::L1TowerJetCollection const&, edm::Handle<edm::View< reco::CaloJet > > const&, 
// 			       TString prefix, Int_t NVTX );



  // Benchmark offline to online jets
  //  void benchmark( edm::Handle<l1slhc::L1TowerJetCollection> const& L1Jets, edm::View< reco::CaloJet > const& offJets ){





  // ----------member data ---------------------------
  ParameterSet conf_;
  edm::InputTag m_jetCollection;
  edm::InputTag m_jetIDMap;





  // Version of algorithm utilised (current/upgrade); eta region: All, central, forward; jet list: PUSub, LPUSub, PrePUSub
  vector <TString> versionList, etaList, jetList;

  // Histogram containers
  std::map< TString, TH1*> hist1D;
  std::map< TString, TH2*> hist2D;
  std::map< TString, TH3*> hist3D;


  // Containers for pT and eta slices
  vector <double> ptSlice, etaSlice, PUSlice, ptThreshold, ttThreshold, ttThresholdFired;

  // Container for PU etaSlice
  vector <double> PUetaSlice;

  // Switch between analysing old (Pre Mk1) and new ntuples
  bool mUseOldNtuple;

  // Tower geometry converter
  TriggerTowerGeometry mTowerGeo;


  // Jet thresholds
  double MIN_L1_PT; 
  double MIN_OFF_PT;


  // Current event count
  Int_t eventCount;


  Service<TFileService> fs;


  // Online-Offline jet matching
  vector <TString> matchPrefix, matchLabel;
  //  vector < std::pair < TString, TString > > matchAxis;

  // Width of iEta bins to sample to obtain pT corrections
  double iEtaCalibBinWidth;
  // Jet calibration iEta-binned pT corrections
  vector <double> iEtaPtCorrection, iEtaPtOffset;

};





AnalyzeJets::AnalyzeJets(const edm::ParameterSet& iConfig): conf_(iConfig){

  // ****************************************
  // *              Start job               *
  // ****************************************
 PRINT("Constructor")






//    produces<int>("NVTX");
//    produces<double>("prePusHT_Tower");
//    produces<double>("HT_Tower");
//    produces<double>("MHT_Tower");
//    produces<double>("PUS_Tower");






    // ****************************************
    // Jet Quality selections
    // ****************************************
    
    // Set to zero...
    MIN_L1_PT  = 0; 
    MIN_OFF_PT = 0;
    
    // ****************************************


  // ************************************************************
  // *                        Calibration                       *
  // ************************************************************
  SUBPRINT("Calibration")

    iEtaCalibBinWidth = iConfig.getParameter< double >("iEtaCalibrationBinWidth");
    
    iEtaPtCorrection  = iConfig.getParameter< vector< double > >("pTCorrection");
    iEtaPtOffset      = iConfig.getParameter< vector< double > >("pTOffset");




    //  eventCount = 0;

  // Using pre-Mk1 ntuples?
  mUseOldNtuple =  iConfig.getParameter < bool > ( "UseOldNtuple" );

  if (mUseOldNtuple){
    std::cout << "WARNING: Running on old ntuple, the Pre-PU subtraction branch is not utilised.\n\n\n";
  }


  // Store the algorithm versions that are utilised
#ifdef CURRENT
  versionList.push_back("curr");
#endif
#ifdef UPGRADE
  versionList.push_back("upgrPrePUS");  // PrePUS
  versionList.push_back("upgr");        // Global PUS
  versionList.push_back("upgrLocal");   // Local PUS
  versionList.push_back("ak5PrePUS");   // ak5 pre PUS and L2L3
  versionList.push_back("PrePUSCalib"); // PrePUSCalib
#endif



  // Store the eta binning - All, central (|eta| < 1.74), forward ( 1.74 < |eta| < 3)
  etaList.push_back("");
  etaList.push_back("cen");
  etaList.push_back("for");

  // Store jet types
#ifdef PROTO_JET
  jetList.push_back("Proto");    // All jet canidates
#endif
#ifdef FILT_CENT
  jetList.push_back("FiltCent"); // Centrality filtered jets
#endif
#ifdef FILT_1D
  jetList.push_back("1DFilt");   // Post 1D filtered jets
#endif
  jetList.push_back("PrePUSub"); // Pre PU-subtracted jets
  jetList.push_back("PUSub");    // Global PU-subtracted jets
  jetList.push_back("LPUSub");   // Local  PU-subtracted jets




  // set the upper and lower limits of pT slices
  // the lower limit is set by the first element and the upper limit the last
  //ptSlice.push_back(10.);  
  ptSlice.push_back(30.);
  ptSlice.push_back(20.);  ptSlice.push_back(40.);  
  ptSlice.push_back(50.);  ptSlice.push_back(75.);
  ptSlice.push_back(100.); ptSlice.push_back(150.);
  // ensure the bins are in ascending order
  sort (ptSlice.begin(), ptSlice.end());


  // set pT threshold
  ptThreshold.push_back(0.);
  ptThreshold.push_back(10.);  ptThreshold.push_back(20.);  
  ptThreshold.push_back(30.);  ptThreshold.push_back(40.);
  ptThreshold.push_back(50.);  ptThreshold.push_back(60.);
  // ensure the bins are in ascending order
  sort (ptThreshold.begin(), ptThreshold.end());


//   for (int ttE = 1; ttE < 11; ttE++)
//     ttThreshold.push_back(ttE);
//   sort (ttThreshold.begin(), ttThreshold.end());
//   // Resize the TT multiplicity array to fit all the threhold values
//   ttThresholdFired.resize(ttThreshold.size());


  // set the upper and lower limits of eta slices
  // L1 eta: region-level segmentation
  etaSlice.push_back(-3.0);
  etaSlice.push_back(-2.172);
  etaSlice.push_back(-1.74);
  etaSlice.push_back(-1.392);
  etaSlice.push_back(-1.044);
  etaSlice.push_back(-0.695);
  etaSlice.push_back(-0.348);
  etaSlice.push_back(0.0);
  etaSlice.push_back(0.348);
  etaSlice.push_back(0.695);
  etaSlice.push_back(1.044);
  etaSlice.push_back(1.392);
  etaSlice.push_back(1.74);
  etaSlice.push_back(2.172);
  etaSlice.push_back(3.0);
  // ensure the bins are in ascending order
  sort (etaSlice.begin(), etaSlice.end());


  // Barrel-endcap
  PUetaSlice.push_back( -3. );
  PUetaSlice.push_back( -1.3 );
  PUetaSlice.push_back( 0. );
  PUetaSlice.push_back( 1.3 );
  PUetaSlice.push_back( 3. );
  // ensure the bins are in ascending order
  sort (PUetaSlice.begin(), PUetaSlice.end());


  // set the upper and lower limits of PU slices
  PUSlice.push_back(0);
  PUSlice.push_back(10);
  PUSlice.push_back(20);
  PUSlice.push_back(30);
  PUSlice.push_back(40);
  PUSlice.push_back(50);
  PUSlice.push_back(60);
  PUSlice.push_back(70);
  PUSlice.push_back(80);
  // ensure the bins are in ascending order
  sort (PUSlice.begin(), PUSlice.end());



}



AnalyzeJets::~AnalyzeJets()
{
}



// ------------ method called once each job just before starting event loop  ------------
void AnalyzeJets::beginJob(){








  // **********************************************************************************
  // *               Book histograms for current and upgrade algorithms               *
  // **********************************************************************************

  // to ensure bins are centred on the 'correct' value:
  // number of bins = [|Max_corr| + |Min_corr|]/binWidth, where Max_corr = Real max + 0.5*binWidth

  //  const double delPTbin(85),   delPTlow(-0.6125), delPThigh(1.5125);
  const double delPTbin(201),  delPTlow(-2.5125), delPThigh(2.5125);
  const double delEtaBin(121), delEtaLow(-0.605), delEtaHigh(0.605);
  const double delPhiBin(121), delPhiLow(-0.605), delPhiHigh(0.605);
  
  const double offPTbin(161),  offPTlow(-0.5),   offPThigh(160.5);
  const double L1PTbin(161),   L1PTlow(-0.5),    L1PThigh(160.5);
  const double offEtaBin(301), offEtaLow(-3.01), offEtaHigh(3.01);
  const double offPhiBin(315),  offPhiLow(-3.15), offPhiHigh(3.15);
  
  const double PUbin(81),     PUlow(-0.5),       PUhigh(80.5);
  const double JetBin(50),    JetLow(-0.5),      JetHigh(49.5);


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



  // ****************************************
  // *          Histogram booking           *
  // ****************************************
  PRINT("Histogram version booking")


  
  for (uint iVer = 0;iVer < versionList.size();iVer++){
    TString prefix = versionList[iVer]; // extract version prefix

    // Print current prefix
    SUBPRINT(prefix)


    TFileDirectory versionDir = fs->mkdir( prefix.Data() );


    // Human-readable version and eta labels
    TString version, etaLab;
    if (prefix == "curr")             version = "Current";  
    else if (prefix == "upgrPrePUS")  version = "Upgrade PrePUS";
    else if (prefix == "upgr")        version = "Upgrade Global PUS";
    else if (prefix == "upgrLocal")   version = "Upgrade Local PUS";
    else if (prefix == "ak5PrePUS")   version = "ak5 PrePUS";
    else if (prefix == "PrePUSCalib") version = "Upgrade PrePUS Calibrated";

    hist1D[prefix + "_HT"]  = versionDir.make <TH1F>(prefix + "_HT" , "H_{T} distribution - " + version + 
					      ";H_{T} (GeV);Entries/ 1 GeV", 300, 0., 300.);
    hist1D[prefix + "_MHT"] = versionDir.make <TH1F>(prefix + "_MHT", "#slash{H}_{T} distribution - " + version + 
					      ";#slash{H}_{T} (GeV);Entries/ 1 GeV", 300, 0., 300.);

    // lead jet distributions
    hist1D[prefix + "_Eta"] = versionDir.make<TH1F>(prefix + "_Eta", "Lead jet #eta - " + version + ";#eta;Entries", offEtaBin,offEtaLow,offEtaHigh);
    hist1D[prefix + "_Phi"] = versionDir.make<TH1F>(prefix + "_Phi", "Lead jet #phi - " + version + ";#phi;Entries", offPhiBin,offPhiLow,offPhiHigh);
    hist1D[prefix + "_PT"]  = versionDir.make<TH1F>(prefix + "_PT",  "Lead jet p_{T} - " + version + ";p_{T} (GeV);Entries", offPTbin, offPTlow, offPThigh);

    // jet pT distributions
    hist1D[prefix + "_Jet_1_PT"] = versionDir.make <TH1F>(prefix + "_Jet_1_PT", "Leading jet p_{T} - " + version + 
							  ";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);
    hist1D[prefix + "_Jet_2_PT"] = versionDir.make <TH1F>(prefix + "_Jet_2_PT", "Second jet p_{T} - " + version + 
							  ";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);
    hist1D[prefix + "_Jet_3_PT"] = versionDir.make <TH1F>(prefix + "_Jet_3_PT", "Third jet p_{T} - " + version + 
							  ";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);
    hist1D[prefix + "_Jet_4_PT"] = versionDir.make <TH1F>(prefix + "_Jet_4_PT", "Fourth jet p_{T} - " + version + 
							  ";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);
    // Resolutions
    hist1D[prefix + "_DeltaEta"] = versionDir.make<TH1F>(prefix + "_DeltaEta", "#Delta#eta - " + version + ";#Delta#eta;Entries", delEtaBin,delEtaLow,delEtaHigh);
    hist1D[prefix + "_DeltaPhi"] = versionDir.make<TH1F>(prefix + "_DeltaPhi", "#Delta#phi - " + version + ";#Delta#phi;Entries", delPhiBin,delPhiLow,delPhiHigh);
    hist1D[prefix + "_DeltaPT"]  = versionDir.make<TH1F>(prefix + "_DeltaPT",  "#Deltap_{T} - " + version + ";#Deltap_{T};Entries", delPTbin, delPTlow, delPThigh);
    

    

    // ************************************************************
    // *                Correlation distributions                 *
    // ************************************************************


    // Number of reconstructed Jets as a function of calorimeter position
    //    hist2D[prefix + "_NJETS_Phi_vs_Eta"]  = versionDir.make<TH2F>(prefix + "_NJETS_Phi_vs_Eta", "Number of reconstructed jets #phi vs #eta - " + version +
    //							   ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);


    hist2D[prefix + "_DeltaPT_vs_L1PT"]      = versionDir.make<TH2D>(prefix + "_DeltaPt_vs_L1PT", "#Deltap_{T} vs L1 p_{T} - " + version + 
								     ";L1 p_{T} (GeV);#Deltap_{T}", L1PTbin,L1PTlow,L1PThigh, delPTbin,delPTlow,delPThigh);
    hist2D[prefix + "_DeltaPhi_vs_DeltaEta"] = versionDir.make<TH2D>(prefix + "_DeltaPhi_vs_DeltaEta", "#Delta#phi vs #Delta#eta - " + version + 
								     ";#Delta#eta;#Delta#phi", delEtaBin,delEtaLow,delEtaHigh, delPhiBin,delPhiLow,delPhiHigh);
    // Online-offline quantity correlations
    hist2D[prefix + "_L1PT_vs_OffPT"]     = versionDir.make<TH2D>(prefix + "_L1PT_vs_OffPT", "L1 p_{T} vs Offline p_{T} - " + version + 
								  ";Offline p_{T} (GeV);L1 p_{T} (GeV)", offPTbin,offPTlow,offPThigh, L1PTbin,L1PTlow,L1PThigh);
    hist2D[prefix + "_L1Eta_vs_OffEta"]   = versionDir.make<TH2D>(prefix + "_L1Eta_vs_OffEta", "L1 #eta vs Offline #eta - " + version + 
								  ";Offline #eta;L1 #eta", offEtaBin,offEtaLow,offEtaHigh, offEtaBin,offEtaLow,offEtaHigh);
    hist2D[prefix + "_L1Phi_vs_OffPhi"]   = versionDir.make<TH2D>(prefix + "_L1Phi_vs_OffPhi", "L1 #phi vs Offline #phi - " + version + 
								  ";Offline #phi;L1 #phi", offPhiBin,offPhiLow,offPhiHigh, offPhiBin,offPhiLow,offPhiHigh);
    // pT Correlations
    hist2D[prefix + "_DeltaPT_vs_OffPT"]   = versionDir.make<TH2D>(prefix + "_DeltaPT_vs_OffPT","#DeltaP_{T} vs Offline P_{T} - " + version +
								   ";Offline p_{T} (GeV);#DeltaP_{T}", offPTbin,offPTlow,offPThigh, delPTbin,delPTlow,delPThigh);
    hist2D[prefix + "_DeltaEta_vs_OffPT"]  = versionDir.make<TH2D>(prefix + "_DeltaEta_vs_OffPT","#Delta#eta vs Offline P_{T} - " + version +
								   ";Offline p_{T} (GeV);#Delta#eta", offPTbin,offPTlow,offPThigh, delEtaBin,delEtaLow,delEtaHigh);
    hist2D[prefix + "_DeltaPhi_vs_OffPT"]  = versionDir.make<TH2D>(prefix + "_DeltaPhi_vs_OffPT","#Delta#phi vs Offline P_{T} - " + version +
								   ";Offline p_{T} (GeV);#Delta#phi", offPTbin,offPTlow,offPThigh, delPhiBin,delPhiLow,delPhiHigh);
    // NVTX Correlations
    hist2D[prefix + "_DeltaPT_vs_NVTX"]    = versionDir.make<TH2D>(prefix + "_DeltaPT_vs_NVTX","#DeltaP_{T} vs N_{VTX} - " + version +
								   ";N_{VTX};#DeltaP_{T}", PUbin,PUlow,PUhigh, delPTbin,delPTlow,delPThigh);
    hist2D[prefix + "_DeltaEta_vs_NVTX"]   = versionDir.make<TH2D>(prefix + "_DeltaEta_vs_NVTX","#Delta#eta vs N_{VTX} - " + version +
								   ";N_{VTX};#Delta#eta", PUbin,PUlow,PUhigh, delEtaBin,delEtaLow,delEtaHigh);
    hist2D[prefix + "_DeltaPhi_vs_NVTX"]   = versionDir.make<TH2D>(prefix + "_DeltaPhi_vs_NVTX","#Delta#phi vs N_{VTX} - " + version +
								   ";N_{VTX};#Delta#phi", PUbin,PUlow,PUhigh, delPhiBin,delPhiLow,delPhiHigh);
    // Phi Correlations
    hist2D[prefix + "_DeltaPT_vs_Phi"]     = versionDir.make<TH2D>(prefix + "_DeltaPT_vs_Phi","#DeltaP_{T} vs #phi - " + version +
								   ";#phi;#DeltaP_{T}", offPhiBin,offPhiLow,offPhiHigh, delPTbin,delPTlow,delPThigh);
    hist2D[prefix + "_DeltaEta_vs_Phi"]    = versionDir.make<TH2D>(prefix + "_DeltaEta_vs_Phi","#Delta#eta vs #phi - " + version +
								   ";#phi;#Delta#eta", offPhiBin,offPhiLow,offPhiHigh, delEtaBin,delEtaLow,delEtaHigh);
    hist2D[prefix + "_DeltaPhi_vs_Phi"]    = versionDir.make<TH2D>(prefix + "_DeltaPhi_vs_Phi","#Delta#phi vs #phi - " + version +
								   ";#phi;#Delta#phi", offPhiBin,offPhiLow,offPhiHigh, delPhiBin,delPhiLow,delPhiHigh);
    // Eta Correlations
    hist2D[prefix + "_DeltaPT_vs_Eta"]     = versionDir.make<TH2D>(prefix + "_DeltaPT_vs_Eta","#DeltaP_{T} vs #eta - " + version +
								   ";#eta;#DeltaP_{T}", offEtaBin,offEtaLow,offEtaHigh, delPTbin,delPTlow,delPThigh);
    hist2D[prefix + "_DeltaEta_vs_Eta"]    = versionDir.make<TH2D>(prefix + "_DeltaEta_vs_Eta","#Delta#eta vs #eta - " + version +
								   ";#eta;#Delta#eta", offEtaBin,offEtaLow,offEtaHigh, delEtaBin,delEtaLow,delEtaHigh);
    hist2D[prefix + "_DeltaPhi_vs_Eta"]    = versionDir.make<TH2D>(prefix + "_DeltaPhi_vs_Eta","#Delta#phi vs #eta - " + version +
								   ";#eta;#Delta#phi", offEtaBin,offEtaLow,offEtaHigh, delPhiBin,delPhiLow,delPhiHigh);


    
    // ************************************************************
    // Profile histograms
    // ************************************************************

    hist1D[prefix + "_DeltaPT_vs_L1PT_prof"]      = versionDir.make<TProfile>(prefix + "_DeltaPt_vs_L1PT_prof", "#Deltap_{T} vs L1 p_{T} profile - " + version + 
								     ";L1 p_{T} (GeV);#Deltap_{T}", L1PTbin,L1PTlow,L1PThigh, delPTlow,delPThigh);
    hist1D[prefix + "_DeltaPhi_vs_DeltaEta_prof"] = versionDir.make<TProfile>(prefix + "_DeltaPhi_vs_DeltaEta_prof", "#Delta#phi vs #Delta#eta profile - " + version + 
								     ";#Delta#eta;#Delta#phi", delEtaBin,delEtaLow,delEtaHigh, delPhiLow,delPhiHigh);
    // Online-offline quantity correlations
    hist1D[prefix + "_L1PT_vs_OffPT_prof"]     = versionDir.make<TProfile>(prefix + "_L1PT_vs_OffPT_prof", "L1 p_{T} vs Offline p_{T} profile - " + version + 
								  ";Offline p_{T} (GeV);L1 p_{T} (GeV)", offPTbin,offPTlow,offPThigh, L1PTlow,L1PThigh);
    hist1D[prefix + "_L1Eta_vs_OffEta_prof"]   = versionDir.make<TProfile>(prefix + "_L1Eta_vs_OffEta_prof", "L1 #eta vs Offline #eta profile - " + version + 
								  ";Offline #eta;L1 #eta", offEtaBin,offEtaLow,offEtaHigh, offEtaLow,offEtaHigh);
    hist1D[prefix + "_L1Phi_vs_OffPhi_prof"]   = versionDir.make<TProfile>(prefix + "_L1Phi_vs_OffPhi_prof", "L1 #phi vs Offline #phi profile - " + version + 
								  ";Offline #phi;L1 #phi", offPhiBin,offPhiLow,offPhiHigh, offPhiLow,offPhiHigh);
    // pT Correlations
    hist1D[prefix + "_DeltaPT_vs_OffPT_prof"]   = versionDir.make<TProfile>(prefix + "_DeltaPT_vs_OffPT_prof","#DeltaP_{T} vs Offline P_{T} profile - " + version +
								   ";Offline p_{T} (GeV);#DeltaP_{T}", offPTbin,offPTlow,offPThigh, delPTlow,delPThigh);
    hist1D[prefix + "_DeltaEta_vs_OffPT_prof"]  = versionDir.make<TProfile>(prefix + "_DeltaEta_vs_OffPT_prof","#Delta#eta vs Offline P_{T} profile - " + version +
								   ";Offline p_{T} (GeV);#Delta#eta", offPTbin,offPTlow,offPThigh, delEtaLow,delEtaHigh);
    hist1D[prefix + "_DeltaPhi_vs_OffPT_prof"]  = versionDir.make<TProfile>(prefix + "_DeltaPhi_vs_OffPT_prof","#Delta#phi vs Offline P_{T} profile - " + version +
								   ";Offline p_{T} (GeV);#Delta#phi", offPTbin,offPTlow,offPThigh, delPhiLow,delPhiHigh);
    // NVTX Correlations
    hist1D[prefix + "_DeltaPT_vs_NVTX_prof"]    = versionDir.make<TProfile>(prefix + "_DeltaPT_vs_NVTX_prof","#DeltaP_{T} vs N_{VTX} profile - " + version +
								   ";N_{VTX};#DeltaP_{T}", PUbin,PUlow,PUhigh, delPTlow,delPThigh);
    hist1D[prefix + "_DeltaEta_vs_NVTX_prof"]   = versionDir.make<TProfile>(prefix + "_DeltaEta_vs_NVTX_prof","#Delta#eta vs N_{VTX} profile - " + version +
								   ";N_{VTX};#Delta#eta", PUbin,PUlow,PUhigh, delEtaLow,delEtaHigh);
    hist1D[prefix + "_DeltaPhi_vs_NVTX_prof"]   = versionDir.make<TProfile>(prefix + "_DeltaPhi_vs_NVTX_prof","#Delta#phi vs N_{VTX} profile - " + version +
								   ";N_{VTX};#Delta#phi", PUbin,PUlow,PUhigh, delPhiLow,delPhiHigh);
    // Phi Correlations
    hist1D[prefix + "_DeltaPT_vs_Phi_prof"]     = versionDir.make<TProfile>(prefix + "_DeltaPT_vs_Phi_prof","#DeltaP_{T} vs #phi profile - " + version +
								   ";#phi;#DeltaP_{T}", offPhiBin,offPhiLow,offPhiHigh, delPTlow,delPThigh);
    hist1D[prefix + "_DeltaEta_vs_Phi_prof"]    = versionDir.make<TProfile>(prefix + "_DeltaEta_vs_Phi_prof","#Delta#eta vs #phi profile - " + version +
								   ";#phi;#Delta#eta", offPhiBin,offPhiLow,offPhiHigh, delEtaLow,delEtaHigh);
    hist1D[prefix + "_DeltaPhi_vs_Phi_prof"]    = versionDir.make<TProfile>(prefix + "_DeltaPhi_vs_Phi_prof","#Delta#phi vs #phi profile - " + version +
								   ";#phi;#Delta#phi", offPhiBin,offPhiLow,offPhiHigh, delPhiLow,delPhiHigh);
    // Eta Correlations
    hist1D[prefix + "_DeltaPT_vs_Eta_prof"]     = versionDir.make<TProfile>(prefix + "_DeltaPT_vs_Eta_prof","#DeltaP_{T} vs #eta profile - " + version +
								   ";#eta;#DeltaP_{T}", offEtaBin,offEtaLow,offEtaHigh, delPTlow,delPThigh);
    hist1D[prefix + "_DeltaEta_vs_Eta_prof"]    = versionDir.make<TProfile>(prefix + "_DeltaEta_vs_Eta_prof","#Delta#eta vs #eta profile - " + version +
								   ";#eta;#Delta#eta", offEtaBin,offEtaLow,offEtaHigh, delEtaLow,delEtaHigh);
    hist1D[prefix + "_DeltaPhi_vs_Eta_prof"]    = versionDir.make<TProfile>(prefix + "_DeltaPhi_vs_Eta_prof","#Delta#phi vs #eta profile - " + version +
									    ";#eta;#Delta#phi", offEtaBin,offEtaLow,offEtaHigh, delPhiLow,delPhiHigh);
    
    

    // ************************************************************
    // *                     Binned PT plots                      *
    // ************************************************************

    // create a pT binned subdirectory
    TFileDirectory pTbinnedDir = versionDir.mkdir( "pTBinned" );

    for (unsigned int iSlice = 0; iSlice < ptSlice.size() - 1; iSlice++){

      // get pT slice range
      double ptLow  = ptSlice[iSlice];
      double ptHigh = ptSlice[iSlice + 1];
      TString ptLowStr  = Form("%d",Int_t(ptLow));
      TString ptHighStr = Form("%d",Int_t(ptHigh));

      TString ptLabel = "pT(" + ptLowStr + "-" + ptHighStr + ")";
      TString ptTitle = "    (" + ptLowStr + " < p_{T}^{off} < " + ptHighStr + ")";


      // create a pT binned subsubdirectory
      TFileDirectory pTbinnedSubDir = pTbinnedDir.mkdir( ptLabel.Data() );


      hist1D[prefix + "_DeltaEta_" + ptLabel] = pTbinnedSubDir.make<TH1F>(prefix + "_DeltaEta_" + ptLabel, "#Delta#eta " + ptTitle + 
									  " - " + version + ";#Delta#eta;Entries", delEtaBin,delEtaLow,delEtaHigh);
      hist1D[prefix + "_DeltaPhi_" + ptLabel] = pTbinnedSubDir.make<TH1F>(prefix + "_DeltaPhi_" + ptLabel, "#Delta#phi" + ptTitle + 
									  " - " + version + ";#Delta#phi;Entries", delPhiBin,delPhiLow,delPhiHigh);
      hist1D[prefix + "_DeltaPT_"  + ptLabel] = pTbinnedSubDir.make<TH1F>(prefix + "_DeltaPT_"  + ptLabel, "#Deltap_{T}" + ptTitle + 
									  " - " + version + ";#Deltap_{T};Entries", delPTbin,delPTlow,delPThigh);


      hist1D[prefix + "_offEta_"   + ptLabel] = pTbinnedSubDir.make<TH1F>(prefix + "_offEta_" + ptLabel, "Offline #eta " + ptTitle + 
									  " - " + version + ";#eta;Entries", offEtaBin,offEtaLow,offEtaHigh);
      hist1D[prefix + "_offPhi_"   + ptLabel] = pTbinnedSubDir.make<TH1F>(prefix + "_offPhi_" + ptLabel, "Offline #phi" + ptTitle + 
									  " - " + version + ";#phi;Entries", offPhiBin,offPhiLow,offPhiHigh);

      hist1D[prefix + "_PU_"       + ptLabel] = pTbinnedSubDir.make<TH1F>(prefix + "_PU_"       + ptLabel, "PU" + ptTitle + 
									  " - " + version + ";PU;Entries", PUbin,PUlow,PUhigh);

      
    }

    // ************************************************************
    // *                    Binned eta plots                      *
    // ************************************************************

    // create an eta binned subdirectory
    //    TFileDirectory etaBinnedDir = fs->mkdir( "etaBinned" );
    TFileDirectory etaBinnedDir = versionDir.mkdir( "etaBinned" );


    for (unsigned int iSlice = 0;iSlice < etaSlice.size() - 1; iSlice++){

      // get pT slice range
      double etaLow  = etaSlice[iSlice];
      double etaHigh = etaSlice[iSlice + 1];
      TString etaLowStr  = Form("%1.3f", etaLow);
      TString etaHighStr = Form("%1.3f", etaHigh);

      TString etaLabel = "eta(" + etaLowStr + "-" + etaHighStr + ")";
      TString etaTitle = "    (" + etaLowStr + " < #eta^{off} < " + etaHighStr + ")";


      // create an eta binned subsubdirectory
      TFileDirectory etaBinnedSubDir = etaBinnedDir.mkdir( etaLabel.Data() );


      hist2D[prefix + "_Off_vs_On_PT_" + etaLabel] = etaBinnedSubDir.make<TH2D>(prefix + "Off_vs_On_PT_" + etaLabel, "Offline vs online p_{T} " + etaTitle + 
										" - " + version + ";Online p_{T} (GeV);Offline p_{T} (GeV)", 
										offPTbin,offPTlow,offPThigh, offPTbin,offPTlow,offPThigh);
      hist1D[prefix + "_Off_vs_On_PT_prof_" + etaLabel] = etaBinnedSubDir.make<TProfile>(prefix + "Off_vs_On_PT_prof_" + etaLabel, "Offline vs online p_{T} profile" 
											 + etaTitle + " - " + version + ";Online p_{T} (GeV);Offline p_{T} (GeV)", 
											 offPTbin,offPTlow,offPThigh, offPTlow,offPThigh);

      hist1D[prefix + "_DeltaEta_" + etaLabel] = etaBinnedSubDir.make<TH1F>(prefix + "_DeltaEta_" + etaLabel, "#Delta#eta " + etaTitle + 
									    " - " + version + ";#Delta#eta;Entries", delEtaBin,delEtaLow,delEtaHigh);
      hist1D[prefix + "_DeltaPhi_" + etaLabel] = etaBinnedSubDir.make<TH1F>(prefix + "_DeltaPhi_" + etaLabel, "#Delta#phi" + etaTitle + 
									    " - " + version + ";#Delta#phi;Entries", delPhiBin,delPhiLow,delPhiLow);
      hist1D[prefix + "_DeltaPT_"  + etaLabel] = etaBinnedSubDir.make<TH1F>(prefix + "_DeltaPT_"  + etaLabel, "#Deltap_{T}" + etaTitle + 
									    " - " + version + ";#Deltap_{T};Entries", delPTbin,delPTlow,delPThigh);

      hist1D[prefix + "_offEta_"   + etaLabel] = etaBinnedSubDir.make<TH1F>(prefix + "_offEta_" + etaLabel, "Offline #eta " + etaTitle + 
									    " - " + version + ";#eta;Entries", offEtaBin,offEtaLow,offEtaHigh);
      hist1D[prefix + "_offPhi_"   + etaLabel] = etaBinnedSubDir.make<TH1F>(prefix + "_offPhi_" + etaLabel, "Offline #phi" + etaTitle + 
									    " - " + version + ";#phi;Entries", offPhiBin,offPhiLow,offPhiHigh);
      hist1D[prefix + "_PU_"       + etaLabel] = etaBinnedSubDir.make<TH1F>(prefix + "_PU_"       + etaLabel, "PU" + etaTitle + 
									    " - " + version + ";PU;Entries", PUbin,PUlow,PUhigh);

    }


    // ************************************************************
    // *                    Eta region plots                      *
    // ************************************************************



    for (uint iEta = 0;iEta < etaList.size();iEta++){
	
      TString etaPre = etaList[iEta]; // extract eta region

      TString etaDir = etaPre;
      if (etaDir == ""){
	etaDir = "all" ;
      } 

      TFileDirectory etaRegionSubDir = etaBinnedDir.mkdir( etaDir.Data() );

      // 'central' jets
      if (etaPre == "cen"){ 
	etaLab = " (|#eta| < 1.3050)"; 
	etaPre = "_" + etaPre;
      }
      // 'forward' jets
      else if (etaPre == "for"){
	etaLab = " (1.3050 < |#eta| < 3)";
	etaPre = "_" + etaPre;
      }
      else {
	etaLab = "";
      }
	

      if ( etaLab == "" ){
	//Jet multiplicity
	hist1D[prefix + etaPre + "_NJETS"] = etaRegionSubDir.make<TH1F>(prefix + etaPre + "_NJETS","Number of Reconstructed Jets - " + version + etaLab + 
							    ";Jet Multiplicity;Events", JetBin,JetLow,JetHigh);
      }



      // reconstructed jet distribution
      hist2D[prefix + etaPre + "_L1Phi_vs_L1Eta"]     = etaRegionSubDir.make<TH2D>(prefix + etaPre + "_L1Phi_vs_L1Eta", 
								       "Number of reconstructed jets L1 #phi vs L1 #eta - " + version + etaLab + 
								       ";L1 #eta;L1 #phi", offEtaBin,offEtaLow,offEtaHigh, offPhiBin,offPhiLow,offPhiHigh);

      hist1D[prefix + etaPre + "_AllJet_L1Eta"]       = etaRegionSubDir.make<TH1F>(prefix + etaPre + "_AllJet_L1Eta", "All jets L1 #eta - " + version + etaLab + 
								       ";#eta;Entries", offEtaBin,offEtaLow,offEtaHigh);
	
      hist1D[prefix + etaPre + "_AllJet_L1Phi"]       = etaRegionSubDir.make<TH1F>(prefix + etaPre + "_AllJet_L1Phi", "All jets L1 #phi - " + version + etaLab + 
								       ";#phi;Entries", offPhiBin,offPhiLow,offPhiHigh);
	
      hist1D[prefix + etaPre + "_AllJet_L1PT"]        = etaRegionSubDir.make <TH1F>(prefix + etaPre + "_AllJet_L1PT", "All jets L1 p_{T} - " + version + etaLab + 
									";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);
	
      hist1D[prefix + etaPre + "_AllJetExcLead_L1PT"] = etaRegionSubDir.make <TH1F>(prefix + etaPre + "_AllJetExcLead_L1PT", 
									"All jets excluding lead jet, L1 p_{T} - " + version + etaLab + 
									";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);

      hist1D[prefix + etaPre + "_13Jet_L1PT"]         = etaRegionSubDir.make <TH1F>(prefix + etaPre + "_13Jet_L1PT", "13 jets L1 p_{T} - " + version + etaLab + 
									";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);
	
      hist1D[prefix + etaPre + "_13JetExcLead_L1PT"]  = etaRegionSubDir.make <TH1F>(prefix + etaPre +"_13JetExcLead_L1PT", 
									"13 jets excluding lead jet, L1 p_{T} - " + version + etaLab + 
									";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);


      
      
    } // end eta loop  

  } // end version loop
  



  // ************************************************************
  // *                       Offline rho                        *
  // ************************************************************


  TFileDirectory kt6Dir = fs->mkdir( "kt6" );


  hist2D["kt6RhoEst_vs_kt6Rho"]     = kt6Dir.make<TH2D>("kt6RhoEst_vs_kt6Rho", "Offline kt6 #rho estimate vs Offline kt6 #rho;kt6 #rho;Estimated kt6 #rho",
							161,-0.125,40.125, 161,-0.125,40.125);
  hist2D["GlobalRho_vs_kt6Rho"]      = kt6Dir.make<TH2D>("GlobalRho_vs_kt6Rho", "Global #rho vs Offline kt6 #rho;kt6 #rho;L1 #rho",
							 161,-0.125,40.125, 161,-0.125,40.125);
  hist2D["TTFired_vs_kt6Rho"]        = kt6Dir.make<TH2D>("TTFired_vs_kt6Rho", "TTs fired vs Offline kt6 #rho;kt6 #rho;TTs fired",
							 161,-0.125,40.125, 201,-0.5,200.5);
  hist2D["EFired_vs_kt6Rho"]         = kt6Dir.make<TH2D>("EFired_vs_kt6Rho", "Ecal cells fired vs Offline kt6 #rho;kt6 #rho;Ecal cells fired",
							 161,-0.125,40.125, 101,-0.5,100.5);
  hist2D["HFired_vs_kt6Rho"]         = kt6Dir.make<TH2D>("HFired_vs_kt6Rho", "Hcal cells fired vs Offline kt6 #rho;kt6 #rho;Hcal cells fired",
							 161,-0.125,40.125, 101,-0.5,100.5);
  hist2D["Etotal_vs_kt6Rho"]         = kt6Dir.make<TH2D>("Etotal_vs_kt6Rho", "Total Ecal energy vs Offline kt6 #rho;kt6 #rho;Total Ecal energy (GeV)",
							 161,-0.125,40.125, 1501,-0.5,1500.5);
  hist2D["Htotal_vs_kt6Rho"]         = kt6Dir.make<TH2D>("Htotal_vs_kt6Rho", "Total Hcal energy vs Offline kt6 #rho;kt6 #rho;Total Hcal energy (GeV)",
							 161,-0.125,40.125, 1501,-0.5,1500.5);
  hist2D["EplusHtotal_vs_kt6Rho"]    = kt6Dir.make<TH2D>("EplusHtotal_vs_kt6Rho", "Total TT energy vs Offline kt6 #rho;kt6 #rho;Total TT energy (GeV)",
							 161,-0.125,40.125, 1501,-0.5,1500.5);
  hist2D["HtPrePUS_vs_kt6Rho"]       = kt6Dir.make<TH2D>("HtPrePUS_vs_kt6Rho", "PrePUS H_{T} vs Offline kt6 #rho;kt6 #rho;PrePUS H_{T} (GeV)",
							 161,-0.125,40.125, 301,-0.5,300.5);
  hist2D["Ht_vs_HtPrePUS"]           = kt6Dir.make<TH2D>("Ht_vs_HtPrePUS", "H_{T} vs PrePUS H_{T};PrePUS H_{T};H_{T} (GeV)",
							 301,-0.5,300.5, 301,-0.5,300.5);
  hist2D["Ht_vs_kt6Rho"]             = kt6Dir.make<TH2D>("Ht_vs_kt6Rho", "H_{T} vs Offline kt6 #rho;kt6 #rho;H_{T} (GeV)",
							 161,-0.125,40.125, 301,-0.5,300.5);
  hist2D["MHt_vs_kt6Rho"]            = kt6Dir.make<TH2D>("MHt_vs_kt6Rho", "#slash{H}_{T} vs Offline kt6 #rho;kt6 #rho;#slash{H}_{T} (GeV)",
							 161,-0.125,40.125, 301,-0.5,300.5);
  hist2D["NJets_vs_kt6Rho"]          = kt6Dir.make<TH2D>("NJets_vs_kt6Rho", "Number of Reconstructed Jets vs Offline kt6 #rho;kt6 #rho;Jet Multiplicity",
							 161,-0.125,40.125, JetBin,JetLow,JetHigh);

  hist2D["EplusHtotal_vs_TTFired"]   = kt6Dir.make<TH2D>("EplusHtotal_vs_TTFired", "Total TT energy vs TTs fired;TTs fired;Total TT energy (GeV)",
							 201,-0.5,200.5, 1501,-0.5,1500.5);
  hist2D["EFired_vs_HFired"]         = kt6Dir.make<TH2D>("EFired_vs_HFired", "Ecal cells fired vs Hcal cells fired;Hcal cells fired;Ecal cells fired",
							 101,-0.5,100.5, 101,-0.5,100.5);
  hist2D["ETotal_vs_EFired"]         = kt6Dir.make<TH2D>("Etotal_vs_EFired", "Total Ecal energy vs Ecal cells fired;Ecal cells fired;Total Ecal energy (GeV)",
                                                         101,-0.5,100.5, 1501,-0.5,1500.5);
  hist2D["HTotal_vs_HFired"]         = kt6Dir.make<TH2D>("Htotal_vs_HFired", "Total Hcal energy vs Hcal cells fired;Hcal cells fired;Total Hcal energy (GeV)",
                                                         101,-0.5,100.5, 1501,-0.5,1500.5);
  hist2D["NJets_vs_EplusHtotal"]     = kt6Dir.make<TH2D>("NJets_vs_EplusHtotal", "Number of Reconstructed Jets vs Total TT energy;Jet Multiplicity;Total TT energy (GeV)",
							 1501,-0.5,1500.5, JetBin,JetLow,JetHigh);
  hist2D["NJets_vs_TTFired"]         = kt6Dir.make<TH2D>("NJets_vs_TTFired", "Number of Reconstructed Jets vs TTs fired;Jet Multiplicity;TTs fired",
							 201,-0.5,200.5, JetBin,JetLow,JetHigh);
    
  // Profiles 
  hist1D["kt6RhoEst_vs_kt6Rho_prof"]   = kt6Dir.make<TProfile>("kt6RhoEst_vs_kt6Rho_prof", "Offline kt6 #rho estimate vs Offline kt6 #rho profile;kt6 #rho;Estimated kt6 #rho",
							       161,-0.125,40.125, -0.125,40.125);
  hist1D["GlobalRho_vs_kt6Rho_prof"]   = kt6Dir.make<TProfile>("GlobalRho_vs_kt6Rho_prof", "Global #rho vs Offline kt6 #rho profile;kt6 #rho;L1 #rho",
							       161,-0.125,40.125, -0.125,40.125);
  hist1D["TTFired_vs_kt6Rho_prof"]     = kt6Dir.make<TProfile>("TTFired_vs_kt6Rho_prof", "TTs fired vs Offline kt6 #rho profile;kt6 #rho;TTs fired",
							       161,-0.125,40.125, -0.5,200.5);
  hist1D["EFired_vs_kt6Rho_prof"]      = kt6Dir.make<TProfile>("EFired_vs_kt6Rho_prof", "Ecal cells fired vs Offline kt6 #rho profile;kt6 #rho;Ecal cells fired",
							       161,-0.125,40.125, -0.5,100.5);
  hist1D["HFired_vs_kt6Rho_prof"]      = kt6Dir.make<TProfile>("HFired_vs_kt6Rho_prof", "Hcal cells fired vs Offline kt6 #rho profile;kt6 #rho;Hcal cells fired",
							       161,-0.125,40.125, -0.5,100.5);
  hist1D["Etotal_vs_kt6Rho_prof"]      = kt6Dir.make<TProfile>("Etotal_vs_kt6Rho_prof", "Total Ecal energy vs Offline kt6 #rho profile;kt6 #rho;Total Ecal energy (GeV)",
							       161,-0.125,40.125, -0.5,1500.5);
  hist1D["Htotal_vs_kt6Rho_prof"]      = kt6Dir.make<TProfile>("Htotal_vs_kt6Rho_prof", "Total Hcal energy vs Offline kt6 #rho profile;kt6 #rho;Total Hcal energy (GeV)",
							       161,-0.125,40.125, -0.5,1500.5);
  hist1D["EplusHtotal_vs_kt6Rho_prof"] = kt6Dir.make<TProfile>("EplusHtotal_vs_kt6Rho_prof", "Total TT energy vs Offline kt6 #rho profile;kt6 #rho;Total TT energy (GeV)",
							       161,-0.125,40.125, -0.5,1500.5);
  hist1D["HtPrePUS_vs_kt6Rho_prof"]    = kt6Dir.make<TProfile>("HtPrePUS_vs_kt6Rho_prof", "PrePUS H_{T} vs Offline kt6 #rho profile;kt6 #rho;PrePUS H_{T} (GeV)",
							       161,-0.125,40.125, -0.5,300.5);
  hist1D["Ht_vs_HtPrePUS_prof"]        = kt6Dir.make<TProfile>("Ht_vs_HtPrePUS_prof", "H_{T} vs PrePUS H_{T} profile;PrePUS H_{T};H_{T} (GeV)",
							       301,-0.5,300.5, -0.5,200.5);
  hist1D["Ht_vs_kt6Rho_prof"]          = kt6Dir.make<TProfile>("Ht_vs_kt6Rho_prof", "H_{T} vs Offline kt6 #rho profile;kt6 #rho;H_{T} (GeV)",
							       161,-0.125,40.125, -0.5,300.5);
  hist1D["MHt_vs_kt6Rho_prof"]         = kt6Dir.make<TProfile>("MHt_vs_kt6Rho_prof", "#slash{H}_{T} vs Offline kt6 #rho profile;kt6 #rho;#slash{H}_{T} (GeV)",
							       161,-0.125,40.125, -0.5,300.5);
  hist1D["NJets_vs_kt6Rho_prof"]       = kt6Dir.make<TProfile>("NJets_vs_kt6Rho_prof", "Number of Reconstructed Jets vs Offline kt6 #rho profile;kt6 #rho;Jet Multiplicity",
							       161,-0.125,40.125, JetLow,JetHigh);
  hist1D["NVTX_vs_kt6Rho_prof"]        = kt6Dir.make<TProfile>("NVTX_vs_kt6Rho_prof", "Number of Reconstructed Vertices vs Offline kt6 #rho profile;kt6 #rho;Vertex Multiplicity",
							       161,-0.125,40.125, PUlow,PUhigh);

  hist2D["TTFired_vs_NVTX"]     = kt6Dir.make<TH2D>("TTFired_vs_NVTX", "TTs fired vs Offline N_{VTX};N_{VTX};TTs fired",
						    PUbin,PUlow,PUhigh, 201,-0.5,200.5);
  hist2D["TTFiredBelow5GeV_vs_NVTX"]     = kt6Dir.make<TH2D>("TTFiredBelow5GeV_vs_NVTX", "TTs fired < 5 GeV vs Offline N_{VTX};N_{VTX};TTs fired",
						    PUbin,PUlow,PUhigh, 201,-0.5,200.5);


  //
  // NVTX correlations
  //

  hist1D["GlobalRho_vs_NVTX_prof"]   = kt6Dir.make<TProfile>("GlobalRho_vs_NVTX_prof", "Global #rho vs Offline N_{VTX} profile;N_{VTX};L1 #rho",
							     PUbin,PUlow,PUhigh, -0.125,40.125);
  hist1D["TTFired_vs_NVTX_prof"]     = kt6Dir.make<TProfile>("TTFired_vs_NVTX_prof", "TTs fired vs Offline N_{VTX} profile;N_{VTX};TTs fired",
							     PUbin,PUlow,PUhigh, -0.5,200.5);
  hist1D["EFired_vs_NVTX_prof"]      = kt6Dir.make<TProfile>("EFired_vs_NVTX_prof", "Ecal cells fired vs Offline N_{VTX} profile;N_{VTX};Ecal cells fired",
							     PUbin,PUlow,PUhigh, -0.5,100.5);
  hist1D["HFired_vs_NVTX_prof"]      = kt6Dir.make<TProfile>("HFired_vs_NVTX_prof", "Hcal cells fired vs Offline N_{VTX} profile;N_{VTX};Hcal cells fired",
							     PUbin,PUlow,PUhigh, -0.5,100.5);
  hist1D["Etotal_vs_NVTX_prof"]      = kt6Dir.make<TProfile>("Etotal_vs_NVTX_prof", "Total Ecal energy vs Offline N_{VTX} profile;N_{VTX};Total Ecal energy (GeV)",
							     PUbin,PUlow,PUhigh, -0.5,1500.5);
  hist1D["Htotal_vs_NVTX_prof"]      = kt6Dir.make<TProfile>("Htotal_vs_NVTX_prof", "Total Hcal energy vs Offline N_{VTX} profile;N_{VTX};Total Hcal energy (GeV)",
							     PUbin,PUlow,PUhigh, -0.5,1500.5);
  hist1D["EplusHtotal_vs_NVTX_prof"] = kt6Dir.make<TProfile>("EplusHtotal_vs_NVTX_prof", "Total TT energy vs Offline N_{VTX} profile;N_{VTX};Total TT energy (GeV)",
							     PUbin,PUlow,PUhigh, -0.5,1500.5);
  hist1D["HtPrePUS_vs_NVTX_prof"]    = kt6Dir.make<TProfile>("HtPrePUS_vs_NVTX_prof", "PrePUS H_{T} vs Offline N_{VTX} profile;N_{VTX};PrePUS H_{T} (GeV)",
							     PUbin,PUlow,PUhigh, -0.5,300.5);
  hist1D["Ht_vs_NVTX_prof"]          = kt6Dir.make<TProfile>("Ht_vs_NVTX_prof", "H_{T} vs Offline N_{VTX} profile;N_{VTX};H_{T} (GeV)",
							     PUbin,PUlow,PUhigh, -0.5,300.5);
  hist1D["MHt_vs_NVTX_prof"]         = kt6Dir.make<TProfile>("MHt_vs_NVTX_prof", "#slash{H}_{T} vs Offline N_{VTX} profile;N_{VTX};#slash{H}_{T} (GeV)",
							     PUbin,PUlow,PUhigh, -0.5,300.5);
  hist1D["NJets_vs_NVTX_prof"]       = kt6Dir.make<TProfile>("NJets_vs_NVTX_prof", "Number of Reconstructed Jets vs Offline N_{VTX} profile;N_{VTX};Jet Multiplicity",
							     PUbin,PUlow,PUhigh, JetLow,JetHigh);


  //
  // Other correlations
  //

  hist1D["EplusHtotal_vs_TTFired_prof"]   = kt6Dir.make<TProfile>("EplusHtotal_vs_TTFired_prof", "Total TT energy vs TTs fired;TTs fired profile;Total TT energy (GeV)",
								  201,-0.5,200.5, -0.5,1500.5);
  hist1D["EFired_vs_HFired_prof"]         = kt6Dir.make<TProfile>("EFired_vs_HFired_prof", "Ecal cells fired vs Hcal cells fired;Hcal cells fired;Ecal cells fired profile",
							 101,-0.5,100.5, -0.5,100.5);
  hist1D["ETotal_vs_EFired_prof"]         = kt6Dir.make<TProfile>("Etotal_vs_EFired_prof", "Total Ecal energy vs Ecal cells fired;Ecal cells fired;Total Ecal energy (GeV) profile",
								  101,-0.5,100.5, -0.5,1500.5);     
  hist1D["HTotal_vs_HFired_prof"]         = kt6Dir.make<TProfile>("Htotal_vs_HFired_prof", "Total Hcal energy vs Hcal cells fired;Hcal cells fired;Total Hcal energy (GeV) profile",
								  101,-0.5,100.5, -0.5,1500.5);     
  hist1D["NJets_vs_EplusHtotal_prof"]     = kt6Dir.make<TProfile>("NJets_vs_EplusHtotal_prof", "Number of Reconstructed Jets vs Total TT energy profile;Jet Multiplicity;Total TT energy (GeV)",
								  1501,-0.5,1500.5, JetLow,JetHigh);
  hist1D["NJets_vs_TTFired_prof"]         = kt6Dir.make<TProfile>("NJets_vs_TTFired_prof", "Number of Reconstructed Jets vs TTs fired profile;Jet Multiplicity;TTs fired",
								  201,-0.5,200.5,   JetLow,JetHigh);


  // *********************************************************
  // *               Eta binned distributions                *
  // *********************************************************

  for (unsigned int iSlice = 0;iSlice < PUetaSlice.size() - 1; iSlice++){
      
    TString jetLab = "kt6";
    TString jetPre = "kt6";
    
       // get eta slice range
       double PUetaLow  = PUetaSlice[iSlice];
       double PUetaHigh = PUetaSlice[iSlice + 1];
       TString PUetaLowStr  = Form("%1.3f", PUetaLow);
       TString PUetaHighStr = Form("%1.3f", PUetaHigh);
      
       TString PUetaLabel = "eta(" + PUetaLowStr + "-" + PUetaHighStr + ")";
       TString PUetaTitle = jetLab + "    (" + PUetaLowStr + " < #eta^{off} < " + PUetaHighStr + ")";

       // create an eta binned subsubdirectory
       TFileDirectory PUetaBinnedSubDir = kt6Dir.mkdir( PUetaLabel.Data() );

       hist1D[jetPre + "_Eta_" + PUetaLabel]    = PUetaBinnedSubDir.make<TH1F>(jetPre + "_AllJet_L1Eta_" + PUetaLabel, "All jets L1 #eta - " + PUetaTitle +
									       ";#eta;Entries", TTetaBin,TTetaBins);
       hist1D[jetPre + "_Phi_" + PUetaLabel]    = PUetaBinnedSubDir.make<TH1F>(jetPre + "_AllJet_L1Phi_" + PUetaLabel, "All jets L1 #phi - " + PUetaTitle +
									       ";#phi;Entries", TTphiBin,TTphiBins);
       hist1D[jetPre + "_PT_" + PUetaLabel]     = PUetaBinnedSubDir.make<TH1F>(jetPre + "_AllJet_L1PT_" + PUetaLabel, "All jets L1 p_{T} - " + PUetaTitle +
									       ";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);
       hist1D[jetPre + "_JetRho_" + PUetaLabel] = PUetaBinnedSubDir.make<TH1F>(jetPre + "_JetRho_" + PUetaLabel,"Jet energy density, #rho_{Jet} " + 
									       PUetaTitle + ";#rho_{Jet};Entries", 161, -0.125,40.125);
       hist1D[jetPre + "_Rho_" + PUetaLabel]    = PUetaBinnedSubDir.make<TH1F>(jetPre + "_Rho_" + PUetaLabel,"Median jet energy density, #rho " + 
									       PUetaTitle + ";#rho;Entries", 161, -0.125,40.125);
       hist1D[jetPre + "_NJETS_" + PUetaLabel]  = PUetaBinnedSubDir.make<TH1F>(jetPre + "_NJETS_" + PUetaLabel,"Number of Reconstructed Jets - " + 
									       PUetaTitle + ";Jet Multiplicity;Events", JetBin,JetLow,JetHigh);



       // Profile plots
       hist1D["LocalRho_vs_kt6LRho_prof_" + PUetaLabel]    = PUetaBinnedSubDir.make<TProfile>("LocalRho_vs_kt6LRho_prof_" + PUetaLabel,
											    "L1 Local #rho vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};L1 #rho_{Local}",
											    161,-0.125,40.125, -0.125,40.125);
       hist1D["kt6Rho_vs_kt6LRho_prof_" + PUetaLabel]      = PUetaBinnedSubDir.make<TProfile>("GlobalRho_vs_kt6LRho_prof_" + PUetaLabel,
											    "Global #rho vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};L1 #rho",
											    161,-0.125,40.125, -0.125,40.125);
       hist1D["TTFired_vs_kt6LRho_prof_" + PUetaLabel]     = PUetaBinnedSubDir.make<TProfile>("TTFired_vs_kt6LRho_prof_" + PUetaLabel,
											    "TTs fired vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};TTs fired",
											    161,-0.125,40.125, -0.5,200.5);
       hist1D["EFired_vs_kt6LRho_prof_" + PUetaLabel]      = PUetaBinnedSubDir.make<TProfile>("EFired_vs_kt6LRho_prof_" + PUetaLabel,
											    "Ecal cells fired vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};Ecal cells fired",
											    161,-0.125,40.125, -0.5,100.5);
       hist1D["HFired_vs_kt6LRho_prof_" + PUetaLabel]      = PUetaBinnedSubDir.make<TProfile>("HFired_vs_kt6LRho_prof_" + PUetaLabel,
											    "Hcal cells fired vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};Hcal cells fired",
											    161,-0.125,40.125, -0.5,100.5);
       hist1D["Etotal_vs_kt6LRho_prof_" + PUetaLabel]      = PUetaBinnedSubDir.make<TProfile>("Etotal_vs_kt6LRho_prof_" + PUetaLabel,
											    "Total Ecal energy vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};Total Ecal energy (GeV)",
											    161,-0.125,40.125, -0.5,1500.5);
       hist1D["Htotal_vs_kt6LRho_prof_" + PUetaLabel]      = PUetaBinnedSubDir.make<TProfile>("Htotal_vs_kt6LRho_prof_" + PUetaLabel,
											    "Total Hcal energy vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};Total Hcal energy (GeV)",
											    161,-0.125,40.125, -0.5,1500.5);
       hist1D["EplusHtotal_vs_kt6LRho_prof_" + PUetaLabel] = PUetaBinnedSubDir.make<TProfile>("EplusHtotal_vs_kt6LRho_prof_" + PUetaLabel,
											    "Total TT energy vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};Total TT energy (GeV)",
											    161,-0.125,40.125, -0.5,1500.5);
       hist1D["HtPrePUS_vs_kt6LRho_prof_" + PUetaLabel]    = PUetaBinnedSubDir.make<TProfile>("HtPrePUS_vs_kt6LRho_prof_" + PUetaLabel,
											    "PrePUS H_{T} vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};PrePUS H_{T} (GeV)",
											    161,-0.125,40.125, -0.5,300.5);
       hist1D["Ht_vs_kt6LRho_prof_" + PUetaLabel]          = PUetaBinnedSubDir.make<TProfile>("Ht_vs_kt6LRho_prof_" + PUetaLabel,
											    "H_{T} vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};H_{T} (GeV)",
											    161,-0.125,40.125, -0.5,300.5);
       hist1D["MHt_vs_kt6LRho_prof_" + PUetaLabel]         = PUetaBinnedSubDir.make<TProfile>("MHt_vs_kt6LRho_prof_" + PUetaLabel,
											    "#slash{H}_{T} vs Offline kt6 #rho_{local} profile - " + PUetaTitle + 
											    ";kt6 #rho_{local};#slash{H}_{T} (GeV)",
											    161,-0.125,40.125, -0.5,300.5);
       hist1D["LNJets_vs_kt6LRho_prof_" + PUetaLabel]       = PUetaBinnedSubDir.make<TProfile>("LNJets_vs_kt6LRho_prof_" + PUetaLabel,
											    "Number of locally Reconstructed Jets vs Offline kt6 #rho_{local} profile - " + 
											    PUetaTitle + ";kt6 #rho_{local};Local Jet Multiplicity",
											    161,-0.125,40.125, JetLow,JetHigh);


  } // End eta binning
    






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

 
  

  // ************************************************************
  // *                        Jet plots                         *
  // ************************************************************
  PRINT("Jet histograms")



  TFileDirectory jetDir = fs->mkdir( "TowerJet" );

  for (uint iJet = 0;iJet < jetList.size();iJet++){

    TString jetPre = jetList[iJet]; // extract jet type
    TString jetLab = jetPre;
    SUBPRINT(jetPre)

    TFileDirectory jetSubDir = jetDir.mkdir( jetPre.Data() );

	
    // Jet distributions
    hist1D[jetPre + "_AllJet_L1Eta"] = jetSubDir.make<TH1F>(jetPre + "_AllJet_L1Eta", "All jets L1 #eta - " + jetLab +  
						      ";#eta;Entries", TTetaBin,TTetaBins);
    //						    ";#eta;Entries", 601,-3.005,3.005);
    hist1D[jetPre + "_AllJet_L1Phi"] = jetSubDir.make<TH1F>(jetPre + "_AllJet_L1Phi", "All jets L1 #phi - " + jetLab +  
						      ";#phi;Entries", TTphiBin,TTphiBins);
    //						    ";#phi;Entries", 641,-3.205,3.205);
    hist1D[jetPre + "_AllJet_L1PT"]  = jetSubDir.make<TH1F>(jetPre + "_AllJet_L1PT", "All jets L1 p_{T} - " + jetLab +  
						      ";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);

    hist1D[jetPre + "_AllJet_L1Eta_PtWeight"] = jetSubDir.make<TH1F>(jetPre + "_AllJet_L1Eta_PtWeight", "All jets L1 p_{T}-weighted #eta - " + jetLab +  
						      ";#eta;p_{T}-weighted Entries", TTetaBin,TTetaBins);
    //						    ";#eta;Entries", 601,-3.005,3.005);
    hist1D[jetPre + "_AllJet_L1Phi_PtWeight"] = jetSubDir.make<TH1F>(jetPre + "_AllJet_L1Phi_PtWeight", "All jets L1 p_{T}-weighted #phi - " + jetLab +  
						      ";#phi;p_{T}-weighted Entries", TTphiBin,TTphiBins);


    hist2D[jetPre + "_AllJet_PhivsEta"]   = jetSubDir.make<TH2D>(jetPre + "_AllJet_PhivsEta", "All jets L1 #phi vs #eta - " + jetLab + 
								 ";#eta;#phi;", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
    hist2D[jetPre + "_AllJet_iPhivsiEta"] = jetSubDir.make<TH2D>(jetPre + "_AllJet_iPhivsiEta", "All jets L1 i#phi vs i#eta - " + jetLab + 
								 ";i#eta;i#phi;", 57,-28.5,28.5, 72,0.5,72.5);


    hist2D[jetPre + "_AllJet_PtvsEta"]      = jetSubDir.make<TH2D>(jetPre + "_AllJet_PtvsEta", "All jets L1 p_{T} vs #eta - " + jetLab + 
								   ";#eta;p_{T};", TTetaBin,TTetaBins, 201,-0.5,200.5);
    hist2D[jetPre + "_AllJet_PtvsPhi"]      = jetSubDir.make<TH2D>(jetPre + "_AllJet_PtvsPhi", "All jets L1 p_{T} vs #phi - " + jetLab + 
								   ";#phi;p_{T};", TTphiBin,TTphiBins, 201,-0.5,200.5);
    hist1D[jetPre + "_AllJet_PtvsEta_prof"] = jetSubDir.make<TProfile>(jetPre + "_AllJet_PtvsEta_prof","Profile of p_{T} vs #eta - " + jetLab + 
								       ";#eta;p_{T} (GeV)", TTetaBin,TTetaBins, 0,200);
    hist1D[jetPre + "_AllJet_PtvsPhi_prof"] = jetSubDir.make<TProfile>(jetPre + "_AllJet_PtvsPhi_prof","Profile of p_{T} vs #phi - " + jetLab + 
								       ";#phi;p_{T} (GeV)", TTphiBin,TTphiBins, 0,200);




    hist1D[jetPre + "_iEta"] = jetSubDir.make<TH1F>(jetPre + "_iEta","Jet i#eta " + jetLab + ";i#eta", 57,-28.5, 28.5);
    hist1D[jetPre + "_iPhi"] = jetSubDir.make<TH1F>(jetPre + "_iPhi","Jet i#phi " + jetLab + ";i#phi", 72,0.5,72.5);
    hist1D[jetPre + "_E"]    = jetSubDir.make<TH1F>(jetPre + "_E","Jet TT energy total " + jetLab + ";E (GeV)", 220,-0.5, 219.5);
    hist1D[jetPre + "_AsymEta"]    = jetSubDir.make<TH1F>(jetPre + "_AsymEta","Jet Eta Asymmetry parameter, A_{#eta} " + jetLab + ";A_{#eta}", 301,-150.5, 150.5);
    hist1D[jetPre + "_AsymPhi"]    = jetSubDir.make<TH1F>(jetPre + "_AsymPhi","Jet Phi Asymmetry parameter, A_{#phi} " + jetLab + ";A_{#phi}", 241,-120.5, 120.5);


    //    hist1D[jetPre + "_EcalMAD"]    = jetSubDir.make<TH1F>(jetPre + "_EcalMAD","Jet Ecal Median Absolute Deviation" + jetLab + ";MAD (GeV)", 301,-150.5, 150.5);



    hist1D[jetPre + "_UnweightedEta"] = jetSubDir.make<TH1F>(jetPre + "_UnweightedEta","Jet unweighted #eta " + jetLab + ";Unweighted #eta",  offEtaBin,offEtaLow,offEtaHigh);
    hist1D[jetPre + "_UnweightedPhi"] = jetSubDir.make<TH1F>(jetPre + "_UnweightedPhi","Jet unweighted #phi " + jetLab + ";Unweighted #phi",  offPhiBin,offPhiLow,offPhiHigh);



//     hist2D[jetPre + "_JetRealArea_Eta"] = jetSubDir.make<TH2D>(jetPre + "_JetRealArea_Eta","Variation of real #Delta#phi#Delta#eta area with weighted #eta " + 
// 							 jetLab + ";#eta;Jet area", TTetaBin,TTetaBins, 1500, 0., 1.5);
//     hist2D[jetPre + "_JetRealArea_Phi"] = jetSubDir.make<TH2D>(jetPre + "_JetRealArea_Phi","Variation of real #Delta#phi#Delta#eta area with weighted #phi " + 
// 							 jetLab + ";#phi;Jet area", TTphiBin,TTphiBins, 1500, 0., 1.5);
//     hist2D[jetPre + "_JetRealArea_iEta"] = jetSubDir.make<TH2D>(jetPre + "_JetRealArea_iEta","Variation of real #Delta#phi#Delta#eta area with i#eta " + 
// 							  jetLab + ";i#eta;Jet area", 57,-28.5, 28.5, 1500,0.,1.5);
//     hist2D[jetPre + "_JetRealArea_iPhi"] = jetSubDir.make<TH2D>(jetPre + "_JetRealArea_iPhi","Variation of real #Delta#phi#Delta#eta area with i#phi " + 
// 							  jetLab + ";i#phi;Jet area", 72,0.5,72.5,    1500,0.,1.5);
//     hist2D[jetPre + "_JetRealArea_Pt"]  = jetSubDir.make<TH2D>(jetPre + "_JetRealArea_Pt","Variation of p_{T} with real #Delta#phi#Delta#eta area " + 
// 							 jetLab + ";Jet area;p_{T} (GeV)", 350, .35,0.7, 202,-0.5,100.5);

  

    hist1D[jetPre + "_Rho"]             = jetSubDir.make<TH1F>(jetPre + "_Rho","Median jet energy density, #rho " + 
							       jetLab + ";#rho;Entries", 121, -0.25,60.25);

    hist2D[jetPre + "_Pt-AreaRatio_Eta"] = jetSubDir.make<TH2D>(jetPre + "_Pt-AreaRatio_Eta","Pt-JetArea Ratio vs #eta " + 
								jetLab + ";#eta;p_{T}/#Delta#ph#Delta#eta (GeV)", TTetaBin,TTetaBins, 51,-2,202);
    hist2D[jetPre + "_Pt-AreaRatio_Phi"] = jetSubDir.make<TH2D>(jetPre + "_Pt-AreaRatio_Phi","Pt-JetArea Ratio vs #phi " + 
								jetLab + ";#phi;p_{T}/#Delta#ph#Delta#eta (GeV)", TTphiBin,TTphiBins, 51,-2,202);
    hist2D[jetPre + "_Pt-AreaRatio_iEta"] = jetSubDir.make<TH2D>(jetPre + "_Pt-AreaRatio_iEta","Pt-JetArea Ratio vs i#eta " + 
								 jetLab + ";i#eta;p_{T}/#Delta#phi#Delta#eta (GeV)", 57,-28.5,28.5,   51,-2,202);
    hist2D[jetPre + "_Pt-AreaRatio_iPhi"] = jetSubDir.make<TH2D>(jetPre + "_Pt-AreaRatio_iPhi","Pt-JetArea Ratio vs i#phi " + 
								 jetLab + ";i#phi;p_{T}/#Delta#phi#Delta#eta (GeV)", 72,0.5,72.5,     51,-2,202);




    hist1D[jetPre + "_NJETS"] = jetSubDir.make<TH1F>(jetPre + "_NJETS","Number of Reconstructed Jets - " + 
					       jetLab + ";Jet Multiplicity;Events", JetBin,JetLow,JetHigh);
  

    




    // ************************************************************
    // *               Tower Jet Binned eta plots                 *
    // ************************************************************

    // create a PU eta binned subdirectory
    TFileDirectory PUetaBinnedDir =  jetSubDir.mkdir( "PUetaBinned" );


    hist2D[jetPre + "_RhoECvsRhoB-"]      = PUetaBinnedDir.make<TH2D>(jetPre + "_RhoECvsRhoB-",
								   "Negative #eta, Endcap #rho vs Barrel #rho;#rho_{Barrel};#rho_{Endcap}", 
								   161,-0.25,80.25, 161,-0.25,80.25);

    hist2D[jetPre + "_RhoECvsRhoB+"]      = PUetaBinnedDir.make<TH2D>(jetPre + "_RhoECvsRhoB+",
								   "Positive #eta, Endcap #rho vs Barrel #rho;#rho_{Barrel};#rho_{Endcap}", 
								   161,-0.25,80.25, 161,-0.25,80.25);

    hist2D[jetPre + "_RhoGlobalvsRhoB-"]  = PUetaBinnedDir.make<TH2D>(jetPre + "_RhoGlobalvsRhoB-",
								   "Negative #eta, Global #rho vs Barrel #rho;#rho_{Barrel};#rho_{Global}", 
								   161,-0.25,80.25, 161,-0.25,80.25);

    hist2D[jetPre + "_RhoGlobalvsRhoB+"]  = PUetaBinnedDir.make<TH2D>(jetPre + "_RhoGlobalvsRhoB+",
								   "Positive #eta, Global #rho vs Barrel #rho;#rho_{Barrel};#rho_{Global}", 
								   161,-0.25,80.25, 161,-0.25,80.25);

    hist2D[jetPre + "_RhoGlobalvsRhoEC-"] = PUetaBinnedDir.make<TH2D>(jetPre + "_RhoGlobalvsRhoEC-",
								   "Negative #eta, Global #rho vs Endcap #rho;#rho_{Endcap};#rho_{Global}", 
								   161,-0.25,80.25, 161,-0.25,80.25);

    hist2D[jetPre + "_RhoGlobalvsRhoEC+"] = PUetaBinnedDir.make<TH2D>(jetPre + "_RhoGlobalvsRhoEC+",
								   "Positive #eta, Global #rho vs Endcap #rho;#rho_{Endcap};#rho_{Global}", 
								   161,-0.25,80.25, 161,-0.25,80.25);

    hist1D[jetPre + "_deltaRhoBEC-"] = PUetaBinnedDir.make<TH1F>(jetPre + "_deltaRhoBEC-",
						       "Negative #eta, barrel #rho - endcap #rho, #Delta#rho_{B,EC}^{-};#Delta#rho_{B,EC}^{-};Entries", 
						       161, -40.25,40.25);
    
    hist1D[jetPre + "_deltaRhoBEC+"] = PUetaBinnedDir.make<TH1F>(jetPre + "_deltaRhoBEC+",
						       "Positive #eta, barrel #rho - endcap #rho, #Delta#rho_{B,EC}^{+};#Delta#rho_{B,EC}^{+};Entries", 
						       161, -40.25,40.25);




    // ************************************************************
    // *                    PU Eta slice plots                    *
    // ************************************************************
    

    for (unsigned int iSlice = 0;iSlice < PUetaSlice.size() - 1; iSlice++){
      
      // get eta slice range
      double PUetaLow  = PUetaSlice[iSlice];
      double PUetaHigh = PUetaSlice[iSlice + 1];
      TString PUetaLowStr  = Form("%1.3f", PUetaLow);
      TString PUetaHighStr = Form("%1.3f", PUetaHigh);
      
      TString PUetaLabel = "eta(" + PUetaLowStr + "-" + PUetaHighStr + ")";
      TString PUetaTitle = jetLab + "    (" + PUetaLowStr + " < #eta^{off} < " + PUetaHighStr + ")";


      // create an eta binned subsubdirectory
      TFileDirectory PUetaBinnedSubDir = PUetaBinnedDir.mkdir( PUetaLabel.Data() );

      hist1D[jetPre + "_Eta_" + PUetaLabel]    = PUetaBinnedSubDir.make<TH1F>(jetPre + "_AllJet_L1Eta_" + PUetaLabel, "All jets L1 #eta - " + PUetaTitle +
									      ";#eta;Entries", TTetaBin,TTetaBins);
      hist1D[jetPre + "_Phi_" + PUetaLabel]    = PUetaBinnedSubDir.make<TH1F>(jetPre + "_AllJet_L1Phi_" + PUetaLabel, "All jets L1 #phi - " + PUetaTitle +
									      ";#phi;Entries", TTphiBin,TTphiBins);
      hist1D[jetPre + "_PT_" + PUetaLabel]     = PUetaBinnedSubDir.make<TH1F>(jetPre + "_AllJet_L1PT_" + PUetaLabel, "All jets L1 p_{T} - " + PUetaTitle +
									      ";p_{T} (GeV);Entries/(1 GeV)", 201,-0.5,200.5);
      hist1D[jetPre + "_JetRho_" + PUetaLabel] = PUetaBinnedSubDir.make<TH1F>(jetPre + "_JetRho_" + PUetaLabel,"Jet energy density, #rho_{Jet} " + 
									      PUetaTitle + ";#rho_{Jet};Entries", 161, -0.25,80.25);
      hist1D[jetPre + "_Rho_" + PUetaLabel]    = PUetaBinnedSubDir.make<TH1F>(jetPre + "_Rho_" + PUetaLabel,"Median jet energy density, #rho " + 
									      PUetaTitle + ";#rho;Entries", 161, -0.25,80.25);

      hist1D[jetPre + "_RhoMin_" + PUetaLabel] = PUetaBinnedSubDir.make<TH1F>(jetPre + "_RhoMin_" + PUetaLabel,
									      "Minimum jet energy density, #rho_{Min} " +
									      PUetaTitle + ";#Delta#rho_{Min};Entries", 161, -0.25,80.25);
      
      hist1D[jetPre + "_RhoMax_" + PUetaLabel] = PUetaBinnedSubDir.make<TH1F>(jetPre + "_RhoMax_" + PUetaLabel,
									      "Maximum jet energy density, #Delta#rho_{Max} " +
									      PUetaTitle + ";#Delta#rho_{Max};Entries", 161, -0.25,80.25);
      
      hist1D[jetPre + "_deltaRhoMin_" + PUetaLabel] = PUetaBinnedSubDir.make<TH1F>(jetPre + "_deltaRhoMin_" + PUetaLabel,
										   "Minimum jet energy density - median, #Delta#rho_{Min,Median} " +
										   PUetaTitle + ";#Delta#rho_{Min,Median};Entries", 161, -40.25,40.25);
      
      hist1D[jetPre + "_deltaRhoMax_" + PUetaLabel] = PUetaBinnedSubDir.make<TH1F>(jetPre + "_deltaRhoMax_" + PUetaLabel,
										   "Maximum jet energy density - median, #Delta#rho_{Max,Median} " +
										   PUetaTitle + ";#Delta#rho_{Max,Median};Entries", 161, -40.25,40.25);
      


      hist1D[jetPre + "_NJETS_" + PUetaLabel]  = PUetaBinnedSubDir.make<TH1F>(jetPre + "_NJETS_" + PUetaLabel,"Number of Reconstructed Jets - " + 
									      PUetaTitle + ";Jet Multiplicity;Events", JetBin,JetLow,JetHigh);


    } // End eta binning
    


    // ************************************************************
    // *                    PT Threshold plots                    *
    // ************************************************************

    // create a pT binned subdirectory
    TFileDirectory pTThreshDir = jetSubDir.mkdir( "pTThreshold" );

    for (unsigned int iSlice = 0; iSlice < ptThreshold.size(); iSlice++){

      // get pT thresholdslice range
      double ptThresh   = ptThreshold[iSlice];
      TString ptThreshStr = Form("%d",Int_t(ptThresh));

      TString ptLabel = "pT>" + ptThreshStr;
      TString ptTitle = "    (" + ptLabel + ")";

      // create a pT binned subsubdirectory
      TFileDirectory pTThreshSubDir = pTThreshDir.mkdir( ptLabel.Data() );


      hist1D[jetPre + "_Eta_" + ptLabel] = pTThreshDir.make<TH1F>(jetPre + "_Eta_" + ptLabel, "#eta " + ptTitle + 
								  " - " + jetLab + ";#eta;Entries", TTetaBin,TTetaBins);
      hist1D[jetPre + "_Phi_" + ptLabel] = pTThreshDir.make<TH1F>(jetPre + "_Phi_" + ptLabel, "#phi" + ptTitle + 
								  " - " + jetLab + ";#phi;Entries", TTphiBin,TTphiBins);
      hist1D[jetPre + "_PT_"  + ptLabel] = pTThreshDir.make<TH1F>(jetPre + "_PT_"  + ptLabel, "p_{T}" + ptTitle + 
								  " - " + jetLab + ";p_{T} (GeV);Entries", 201,-0.5,200.5);
    }
    



  } // End jet plots

  
  



  // *************************************************************************************
  // *                             Correlations distributions                            *
  // *************************************************************************************
  PRINT("Correlation histograms")

  // Create a correlation subdirectory
  TFileDirectory correlationDir = fs->mkdir( "Correlations" );


  vector <TString> corrPrefix, corrLabel;
  vector < std::pair < TString, TString > > axisPrefix;

#ifdef PROTO_JET
  // Jet candidates vs 1D eta filtered jets
  corrPrefix.push_back("Proto_1DFilt");
  corrLabel.push_back("1D Filtered vs ProtoJet");
  axisPrefix.push_back( std::make_pair("ProtoJet","1D Filtered") );

  // Jet candidates vs 1D eta filtered jets
  corrPrefix.push_back("Proto_PrePUS");
  corrLabel.push_back("Pre PU subtracted vs ProtoJet");
  axisPrefix.push_back( std::make_pair("ProtoJet","Pre PU Subtraction") );

  // 1D eta filtered jets vs PrePUS
  corrPrefix.push_back("Proto_FiltCent");
  corrLabel.push_back("Centrality Filtered vs ProtoJet");
  axisPrefix.push_back( std::make_pair("ProtoJet","Centality Filtered") );

#endif
#ifdef FILT_CENT
  // Centrality filtered jets vs PrePUS
  corrPrefix.push_back("FiltCent_PrePUS");
  corrLabel.push_back("Pre PU subtracted vs Centrality Filtered");
  axisPrefix.push_back( std::make_pair("Centrality Filtered","Pre PU Subtraction") );
#endif
#ifdef FILT_1D
  // 1D eta filtered jets vs PrePUS
  corrPrefix.push_back("1DFilt_PrePUS");
  corrLabel.push_back("Pre PU subtracted vs 1D Filtered");
  axisPrefix.push_back( std::make_pair("1D Filtered","Pre PU Subtraction") );
#endif
  // PrePUS vs Global PUS
  corrPrefix.push_back("PrePUS_PUS");
  corrLabel.push_back("Post PU subtracted vs Pre PU subtracted");
  axisPrefix.push_back( std::make_pair("Pre PU Subtraction","Post PU Subtraction") );

  // PrePUS vs Local PUS
  corrPrefix.push_back("PrePUS_LPUS");
  corrLabel.push_back("Post local PU subtracted vs Pre PU subtracted");
  axisPrefix.push_back( std::make_pair("Pre PU Subtraction","Post local PU Subtraction") );

  // PUS vs Local PUS
  corrPrefix.push_back("PUS_LPUS");
  corrLabel.push_back("Post local PU subtracted vs Post global PU subtracted");
  axisPrefix.push_back( std::make_pair("Post global PU Subtraction","Post local PU Subtraction") );


  
  for (unsigned int iCorr = 0; iCorr < corrPrefix.size(); iCorr++){


    TString corrPre = corrPrefix[iCorr];
    TString corrLab = corrLabel[iCorr];
    TString axisXPre = axisPrefix[iCorr].first;
    TString axisYPre = axisPrefix[iCorr].second;

    // Create a correlation subsubdirectory
    TFileDirectory correlationSubDir = correlationDir.mkdir( corrPre.Data() );




    hist2D["Corr_" + corrPre + "_NJETS"] = correlationSubDir.make<TH2D>("Corr_NJETS","Number of Reconstructed Jets - " + corrLab + 
					  ";" + axisXPre + " N_{Jets};" + axisYPre + " N_{Jets}", JetBin,JetLow,JetHigh, JetBin,JetLow,JetHigh);
    

	
    // Jet distributions
    hist2D["Corr_" + corrPre + "_AllJet_L1Eta"] = correlationSubDir.make<TH2D>("Corr_" + corrPre + "_AllJet_L1Eta", "All jets L1 #eta - " + corrLab +  
							       ";" + axisXPre + " #eta;" + axisYPre + " #eta", TTetaBin,TTetaBins, TTetaBin,TTetaBins);
    hist2D["Corr_" + corrPre + "_AllJet_L1Phi"] = correlationSubDir.make<TH2D>("Corr_" + corrPre + "_AllJet_L1Phi", "All jets L1 #phi - " + corrLab +  
							       ";" + axisXPre + " #phi;" + axisYPre + " #phi", TTphiBin,TTphiBins, TTphiBin,TTphiBins);
    hist2D["Corr_" + corrPre + "_AllJet_L1PT"]  = correlationSubDir.make<TH2D>("Corr_" + corrPre + "_AllJet_L1PT", "All jets L1 p_{T} - " + corrLab +  
							       ";" + axisXPre + " p_{T} (GeV);" + axisYPre + " p_{T} (GeV)", 201,-0.5,200.5, 201,-0.5,200.5);

    // Dead-jet distributions
    TString deadLab = axisYPre + " Dead jets";
    hist1D["Dead_" + corrPre + "_AllJet_L1Eta"] = correlationSubDir.make<TH1F>("Dead_" + corrPre + "_AllJet_L1Eta", "All jets L1 #eta - " + deadLab +
							       ";#eta;Entries", TTetaBin,TTetaBins);
    hist1D["Dead_" + corrPre + "_AllJet_L1Phi"] = correlationSubDir.make<TH1F>("Dead_" + corrPre + "_AllJet_L1Phi", "All jets L1 #phi - " + deadLab +
							       ";#phi;Entries", TTphiBin,TTphiBins);
    hist1D["Dead_" + corrPre + "_AllJet_L1PT"]  = correlationSubDir.make<TH1F>("Dead_" + corrPre + "_AllJet_L1Pt", "All jets L1 p_{T} - " + deadLab +
							       ";p_{T} (GeV);Entries", 201,-0.05,20.05);

    hist1D["Dead_" + corrPre + "_AllJet_L1PTvsL1Eta"]  = correlationSubDir.make<TH2D>("Dead_" + corrPre + "_AllJet_L1PtvsL1Eta", "All jets L1 p_{T} vs #eta - " + deadLab +
								      ";#eta;p_{T} (GeV);", TTetaBin,TTetaBins, 201,-0.05,20.05);
    hist1D["Dead_" + corrPre + "_AllJet_L1PhivsL1Eta"] = correlationSubDir.make<TH2D>("Dead_" + corrPre + "_AllJet_L1PhivsL1Eta", "All jets L1 #phi vs #eta - " + deadLab +
								      ";#eta;#phi;", TTetaBin,TTetaBins, TTphiBin,TTphiBins);


  }





  // ****************************************************************************************************
  // *                                      Event Distributions                                         *
  // ****************************************************************************************************
  PRINT("Event histograms")

  TFileDirectory eventDir = fs->mkdir( "Event" );

  TString prefix = "Event";

  const double evBin(10001), evLow(-50000.0), evHigh(10050000);
  const double lumiBin(201), lumiLow(-0.5), lumiHigh(200.5);
  const double runnrBin(24), runnrLow(198586.5), runnrHigh(198610.5);
  
  
  hist1D["NVTX"]                = eventDir.make<TH1F>("NVTX","Number of ak5 reconstructed primary vertices;PU;Events", PUbin,PUlow,PUhigh);  
  hist1D[prefix + "_Runnr"]     = eventDir.make<TH1F>(prefix + "_Runnr","Run Number;Runnr;Events", runnrBin,runnrLow,runnrHigh);
  hist1D[prefix + "_LumiBlock"] = eventDir.make<TH1F>(prefix + "_LumiBlock","Luminosity block;LumiBlock;Events", lumiBin,lumiLow,lumiHigh);
  hist1D[prefix + "_Eventnr"]   = eventDir.make<TH1F>(prefix + "_Eventnr", "Event number;Eventnr;Events", evBin,evLow,evHigh);
  
  hist2D[prefix + "_NJETS"] = eventDir.make<TH2D>(prefix + "_NJETS","NJETS vs Eventnr;Eventnr;Jet Multiplicity", 
						  evBin,evLow,evHigh, JetBin,JetLow,JetHigh);
  hist2D[prefix + "_NVTX"]  = eventDir.make<TH2D>(prefix + "_NVTX","AK5 reconstructed vertices vs Eventnr;Eventnr;AK5 Vertices", 
						  evBin,evLow,evHigh, JetBin,JetLow,JetHigh);
  hist2D[prefix + "_Rho"]   = eventDir.make<TH2D>(prefix + "_Rho", "Global #rho vs Eventnr;Eventnr;#rho", 
						  evBin,evLow,evHigh, 161,-0.25,80.25);


  const double htBin(301), htLow(-0.5), htHigh(300.5);
  hist2D[prefix + "_HT"]    = eventDir.make<TH2D>(prefix + "_HT", "HT vs Eventnr;Eventnr;H_{T}",
                                                  evBin,evLow,evHigh, htBin,htLow,htHigh);
  hist2D[prefix + "_MHT"]   = eventDir.make<TH2D>(prefix + "_MHT", "MHT vs Eventnr;Eventnr;#slash{H}_{T}",
                                                  evBin,evLow,evHigh, htBin,htLow,htHigh);




//   // ****************************************************************************************************
//   // *                                      Vertex Distributions                                        *
//   // ****************************************************************************************************
// #ifdef VERTEX
//   PRINT("Vertex")

//   TFileDirectory vtxDir = fs->mkdir( "Vertex" );
  

//   prefix = "Vertex";


//   hist1D[prefix + "_Chi2"]           = vtxDir.make<TH1F>(prefix + "_Chi2","Vertex #chi^{2};#chi^{2};Entries",  201,-0.5,200.5);
//   hist1D[prefix + "_Ndof"]           = vtxDir.make<TH1F>(prefix + "_Ndof","Vertex N_{DOF};N_{DOF};Entries",    201,-0.5,200.5);
//   hist1D[prefix + "_Chi2_over_Ndof"] = vtxDir.make<TH1F>(prefix + "_Chi2_over_Ndof","Vertex #chi^{2}/N_{DOF};#chi^{2}/N_{DOF};Entries", 101,-0.5,10.5);
//   hist2D[prefix + "_Ndof_vs_Chi2"]   = vtxDir.make<TH2D>(prefix + "_Ndof_vs_Chi2","Vertex Ndof vs #chi^{2} ;NDof;#chi^{2}", 201,-0.5,200.5, 201,-0.5,200.5);


// #endif


//   // ****************************************************************************************************
//   // *                                        TT Distributions                                          *
//   // ****************************************************************************************************

// #ifdef TT_DIST
//   PRINT("TT histograms")

//   TString TTPre = "TT";
//   TString TTLab = "TT";


//   TFileDirectory TTDir = fs->mkdir( "TT" );

//   // 1D distributions
//   hist1D[TTPre + "_ETowersFired"]     = TTDir.make<TH1F>(TTPre + "_ETowersFired","Total number of Ecal TTs fired " +
// 						       TTLab + ";TTs Fired;Events", 501,-0.5,500.5);
//   hist1D[TTPre + "_HTowersFired"]     = TTDir.make<TH1F>(TTPre + "_HTowersFired","Total number of Hcal TTs fired " +
// 						       TTLab + ";TTs Fired;Events", 501,-0.5,500.5);
//   hist1D[TTPre + "_EorHTowersFired"]  = TTDir.make<TH1F>(TTPre + "_EorHTowersFired","Total number of Ecal or Hcal TTs fired " +
// 						       TTLab + ";TTs Fired;Events", 1001,-0.5,1000.5);
//   hist1D[TTPre + "_EandHTowersFired"] = TTDir.make<TH1F>(TTPre + "_EandHTowersFired","Total number of Ecal and Hcal TTs fired " +
// 						       TTLab + ";TTs Fired;Events", 501,-0.5,500.5);
//   hist1D[TTPre + "_MaxE"]             = TTDir.make<TH1F>(TTPre + "_MaxE","Largest Ecal energy deposit " +
// 						       TTLab + ";Ecal energy (GeV);Events", 201,-0.5,200.5);
//   hist1D[TTPre + "_MaxH"]             = TTDir.make<TH1F>(TTPre + "_MaxH","Largest Hcal energy deposit " +
// 						       TTLab + ";Hcal energy (GeV);Events", 201,-0.5,200.5);
//   hist1D[TTPre + "_Etotal"]           = TTDir.make<TH1F>(TTPre + "_Etotal","Total Ecal energy deposited " +
// 							 TTLab + ";Total Ecal energy (GeV);Events", 1501,-0.5,1500.5);
//   hist1D[TTPre + "_Htotal"]           = TTDir.make<TH1F>(TTPre + "_Htotal","Total Hcal energy deposited " +
// 							 TTLab + ";Total Hcal energy (GeV);Events", 1501,-0.5,1500.5);
//   hist1D[TTPre + "_EplusHtotal"]      = TTDir.make<TH1F>(TTPre + "_EplusHtotal","Total Ecal + Hcal energy deposited " +
// 							 TTLab + ";Total Ecal + Hcal energy (GeV);Events", 1501,-0.5,1500.5);



//   // 2D distributions

// //   hist2D[TTPre + "_Phi_vs_iPhi"]           = TTDir.make<TH2D>(TTPre + "_Phi_vs_iPhi", "_" + 
// // 							    TTLab + ";;", 72,0.5,72.5, TTphiBin,TTphiBins);
// //   hist2D[TTPre + "_Eta_vs_iEta"]           = TTDir.make<TH2D>(TTPre + "_Eta_vs_iEta","_ " + 
// // 							    TTLab + ";;", 57,-28.5,28.5, TTetaBin,TTetaBins);
  
// // Phi-Eta space
//   hist2D[TTPre + "_E-Phi_vs_Eta"]          = TTDir.make<TH2D>(TTPre + "_E-Phi_vs_Eta","Ecal energy #phi vs #eta " + 
// 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//   hist2D[TTPre + "_H-Phi_vs_Eta"]          = TTDir.make<TH2D>(TTPre + "_H-Phi_vs_Eta","Hcal energy #phi vs #eta " + 
// 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//   hist2D[TTPre + "_EplusH-Phi_vs_Eta"]     = TTDir.make<TH2D>(TTPre + "_EplusH-Phi_vs_Eta","Ecal + Hcal energy #phi vs #eta " + 
// 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
// //   hist2D[TTPre + "_EHratio-Phi_vs_Eta"]    = TTDir.make<TH2D>(TTPre + "_EHratio-Phi_vs_Eta","Ecal-Hcal energy ratio #phi vs #eta " + 
// // 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//   hist2D[TTPre + "_EFired-Phi_vs_Eta"]     = TTDir.make<TH2D>(TTPre + "_EFired-Phi_vs_Eta","Ecal fired #phi vs #eta " + 
// 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//   hist2D[TTPre + "_HFired-Phi_vs_Eta"]     = TTDir.make<TH2D>(TTPre + "_HFired-Phi_vs_Eta","Hcal fired #phi vs #eta " + 
// 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//   hist2D[TTPre + "_EorHFired-Phi_vs_Eta"]  = TTDir.make<TH2D>(TTPre + "_EorHFired-Phi_vs_Eta","Ecal or Hcal fired #phi vs #eta " + 
// 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//   hist2D[TTPre + "_EandHFired-Phi_vs_Eta"] = TTDir.make<TH2D>(TTPre + "_EandHFired-Phi_vs_Eta","Ecal and Hcal fired #phi vs #eta " + 
// 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
// // iPhi-iEta space
//   hist2D[TTPre + "_E-iPhi_vs_iEta"]          = TTDir.make<TH2D>(TTPre + "_E-iPhi_vs_iEta","Ecal i#phi vs i#eta " + 
// 							    TTLab + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
//   hist2D[TTPre + "_H-iPhi_vs_iEta"]          = TTDir.make<TH2D>(TTPre + "_H-iPhi_vs_iEta","Hcal i#phi vs i#eta " + 
// 							    TTLab + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
//   hist2D[TTPre + "_EplusH-iPhi_vs_iEta"]     = TTDir.make<TH2D>(TTPre + "_EplusH-iPhi_vs_iEta","Ecal + Hcal energy i#phi vs i#eta " + 
// 							    TTLab + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
// //   hist2D[TTPre + "_EHratio-iPhi_vs_iEta"]    = TTDir.make<TH2D>(TTPre + "_EHratio-iPhi_vs_iEta","Ecal-Hcal energy ratio i#phi vs i#eta " + 
// // 							    TTLab + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
//   hist2D[TTPre + "_EFired-iPhi_vs_iEta"]     = TTDir.make<TH2D>(TTPre + "_EFired-iPhi_vs_iEta","Ecal fired i#phi vs i#eta " + 
// 							    TTLab + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
//   hist2D[TTPre + "_HFired-iPhi_vs_iEta"]     = TTDir.make<TH2D>(TTPre + "_HFired-iPhi_vs_iEta","Hcal fired i#phi vs i#eta " + 
// 							    TTLab + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
//   hist2D[TTPre + "_EorHFired-iPhi_vs_iEta"]  = TTDir.make<TH2D>(TTPre + "_EorHFired-iPhi_vs_iEta","Ecal or Hcal fired i#phi vs i#eta " + 
// 							    TTLab + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
//   hist2D[TTPre + "_EandHFired-iPhi_vs_iEta"] = TTDir.make<TH2D>(TTPre + "_EandHFired-iPhi_vs_iEta","Ecal and Hcal fired i#phi vs i#eta " + 
// 							    TTLab + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);

//   // Profile histograms
//   hist1D[TTPre + "_E-iEta_prof"]             = TTDir.make<TProfile>(TTPre + "_E-iEta_prof","Profile of E vs i#eta - " + TTLab +
// 								    ";i#eta;E (L1 GeV)",   57,-28.5,28.5, 0, 100);
//   hist1D[TTPre + "_H-iEta_prof"]             = TTDir.make<TProfile>(TTPre + "_H-iEta_prof","Profile of H vs i#eta - " + TTLab +
// 								    ";i#eta;H (L1 GeV)",   57,-28.5,28.5, 0, 100);
//   hist1D[TTPre + "_EplusH-iEta_prof"]        = TTDir.make<TProfile>(TTPre + "_EplusH-iEta_prof","Profile of E+H vs i#eta - " + TTLab +
// 								    ";i#eta;E+H (L1 GeV)", 57,-28.5,28.5, 0, 100);
//   hist1D[TTPre + "_E-iPhi_prof"]             = TTDir.make<TProfile>(TTPre + "_E-iPhi_prof","Profile of E vs i#phi - " + TTLab +
// 								    ";i#phi;E (L1 GeV)",   72,0.5,72.5, 0, 100);
//   hist1D[TTPre + "_H-iPhi_prof"]             = TTDir.make<TProfile>(TTPre + "_H-iPhi_prof","Profile of H vs i#phi - " + TTLab +
// 								    ";i#phi;H (L1 GeV)",   72,0.5,72.5, 0, 100);
//   hist1D[TTPre + "_EplusH-iPhi_prof"]        = TTDir.make<TProfile>(TTPre + "_EplusH-iPhi_prof","Profile of E+H vs i#phi - " + TTLab +
// 								    ";i#phi;E+H (L1 GeV)", 72,0.5,72.5, 0, 100);


//   // ************************************************************
//   // *                 iEta binned distributions                *
//   // ************************************************************

//   SUBPRINT("iEta binned")
//  TFileDirectory ttiEtaDir = TTDir.mkdir( "iEtaBinned" );
//     for (int iEta = -28;iEta <= 28;iEta++){

//       if (iEta == 0)
// 	continue;

//       // Get current iEta
//       TString iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEta)));
//       TFileDirectory ttiEtaSubDir = ttiEtaDir.mkdir( iEtaStr.Data() );


//       hist1D[TTPre + "_E_" + iEtaStr]      = ttiEtaSubDir.make<TH1D>(TTPre + "_E_" + iEtaStr,"E (" + iEtaStr + ") - " + TTLab +
// 								     ";E (L1 GeV);Entries",   101,-0.5,100);
//       hist1D[TTPre + "_H_" + iEtaStr]      = ttiEtaSubDir.make<TH1D>(TTPre + "_H_" + iEtaStr,"H (" + iEtaStr + ") - " + TTLab +
// 								     ";H (L1 GeV);Entries",   101,-0.5,100);
//       hist1D[TTPre + "_EplusH_" + iEtaStr] = ttiEtaSubDir.make<TH1D>(TTPre + "_EplusH_" + iEtaStr,"EplusH (" + iEtaStr + ") - " + TTLab +
// 								     ";E+H (L1 GeV);Entries", 101,-0.5,100);

//     }



// //   hist2D[TTPre + "_EFG-Phi_vs_Eta"]        = TTDir.make<TH2D>(TTPre + "_EFG-Phi_vs_Eta","Ecal FG #phi vs #eta " + 
// // 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
// //   hist2D[TTPre + "_HFG-Phi_vs_Eta"]        = TTDir.make<TH2D>(TTPre + "_HFG-Phi_vs_Eta","Hcal FG #phi vs #eta " + 
// // 							    TTLab + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);

// #endif











// ********************************************************************************
// *                                 BENCHMARKING                                 *
// ******************************************************************************** 

  PRINT("Benchmark")

  // Create an event directory
  TFileDirectory benchDir       = fs->mkdir( "Benchmark" );



  TFileDirectory benchAk5SubDir = benchDir.mkdir( "Ak5" );
  
  
  hist1D["Ak5_JetArea"]         = benchAk5SubDir.make<TH1F>("Ak5_JetArea", "Ak5 Jet area;Jet area;Entries", 151,-0.05,1.55);
  hist2D["Ak5_pT_vs_JetArea"]   = benchAk5SubDir.make<TH2D>("Ak5_pT_vs_JetArea", "Ak5 p_{T} vs Jet area;Jet area;Jet p_{T} (GeV)", 
							    151,-0.05,1.55, offPTbin,offPTlow,offPThigh);
  hist1D["Ak5_TowersArea"]      = benchAk5SubDir.make<TH1F>("Ak5_TowersArea", "Ak5 Tower area;Tower area;Entries", 151,-0.05,1.55);
  hist2D["Ak5_pT_vs_TowersArea"]= benchAk5SubDir.make<TH2D>("Ak5_pT_vs_TowersArea", "Ak5 p_{T} vs Tower area;Tower area;Jet p_{T} (GeV)", 
							    151,-0.05,1.55, offPTbin,offPTlow,offPThigh);
  


//   hist2D["DeltaPt_Lead_vs_All"]= benchAk5SubDir.make<TH2D>("DeltaPt_Lead_vs_All", "Match comparison Lead vs All;Lead Match #Deltap_{T};All Match #Deltap_{T}",
// 							   delPTbin,delPTlow,delPThigh, delPTbin,delPTlow,delPThigh);
//   hist2D["DeltaR_Lead_vs_All"]= benchAk5SubDir.make<TH2D>("DeltaR_Lead_vs_All", "Match comparison Lead vs All;Lead Match #DeltaR;All Match #DeltaR",
// 							  51,-0.025,0.525, 51,-0.025,0.525);





  // ----------------------------------------
  //  Jet collections to be matched
  // ----------------------------------------

  // NOTE: WE ARE ON THE MEMORY BOUNDARY, ADDING ANY MORE
  //       COMPARISONS WILL LIKELY EXCEED THE MEMORY LIMIT : (
  //
  
  // PrePUS-ak5PrePUS
   matchPrefix.push_back("PrePUS_ak5PrePUS");
   matchLabel.push_back("TowerJet PrePUS, ak5 PrePUS");
  // PrePUS-ak5PUS
//   matchPrefix.push_back("PrePUS_ak5PUS");
//   matchLabel.push_back("TowerJet PrePUS, ak5 PUS");
  // PUS-ak5PUS
//   matchPrefix.push_back("PUS_ak5PUS");
//   matchLabel.push_back("TowerJet PUS, ak5 PUS");
  // PrePUSCalib-ak5PrePUS
//   matchPrefix.push_back("PrePUSCalib_ak5PrePUS");
//   matchLabel.push_back("TowerJet PrePUS Calib, ak5 PrePUS");
  // PrePUSCalib+PUS-ak5PUS
  // PrePUSCalib-ak5PUS
  // PUSCalib-ak5PUS



  //  matchAxis.push_back( std::make_pair("Pre PU Subtraction","Post local PU Subtraction") );



  for (unsigned int iMatch = 0; iMatch < matchPrefix.size(); iMatch++){
    
    
    TString matchPre = matchPrefix[iMatch];
    TString matchLab = matchLabel[iMatch];

    // Create a match subsubdirectory
    TFileDirectory matchSubDir = benchDir.mkdir( matchPre.Data() );

    hist2D[matchPre + "_L1PT_vs_OffPT"]        = matchSubDir.make<TH2D>(matchPre + "_L1PT_vs_OffPT",
									"L1 P_{T} vs Offline P_{T};Offline p_{T} (GeV);L1 P_{T} (GeV)", 
									offPTbin,offPTlow,offPThigh, offPTbin,offPTlow,offPThigh);
    hist1D[matchPre + "_OffPT_vs_L1PT_prof"]   = matchSubDir.make<TProfile>(matchPre + "_OffPT_vs_L1PT_prof",
									    "Offline P_{T} vs L1 P_{T};L1 P_{T} (GeV);Offline p_{T} (GeV)",
									    offPTbin,offPTlow,offPThigh, offPTlow,offPThigh);
    hist1D[matchPre + "_DeltaPT_vs_NVTX_prof"] = matchSubDir.make<TProfile>(matchPre + "_DeltaPT_vs_NVTX_prof","#DeltaP_{T} vs N_{VTX};N_{VTX};#DeltaP_{T}", 
									    PUbin,PUlow,PUhigh, delPTlow,delPThigh);
    hist1D[matchPre + "_DeltaPT_vs_Phi_prof"]  = matchSubDir.make<TProfile>(matchPre + "_DeltaPT_vs_Phi_prof","#DeltaP_{T} vs #phi;L1 #phi;#DeltaP_{T}", 
									    offPhiBin,offPhiLow,offPhiHigh, 
									    delPTlow,delPThigh);
    hist1D[matchPre + "_DeltaPT_vs_Eta_prof"]  = matchSubDir.make<TProfile>(matchPre + "_DeltaPT_vs_Eta_prof","#DeltaP_{T} vs #eta;L1 #eta;#DeltaP_{T}", 
									    offEtaBin,offEtaLow,offEtaHigh, 
									    delPTlow,delPThigh);
    



	




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

	//    std::cout << "\n\niEta = " << iEta << " is in the bin range:   " << iEtaBinLow << "to" << iEtaBinHigh << "\n";
      }






    // TString iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEta)));
      TFileDirectory iEtaSubDir = iEtaDir.mkdir( iEtaStr.Data() );



      hist2D[matchPre + "_DeltaPT_vs_L1PT_" + iEtaStr]      = iEtaSubDir.make<TH2D>(matchPre + "DeltaPt_vs_L1PT_" + iEtaStr, 
										    "#Deltap_{T} vs L1 p_{T} (" + iEtaStr +
										    ");L1 p_{T} (GeV);#Deltap_{T}", L1PTbin,L1PTlow,L1PThigh, 
										    delPTbin,delPTlow,delPThigh);
      hist2D[matchPre + "_DeltaPhi_vs_DeltaEta_" + iEtaStr] = iEtaSubDir.make<TH2D>(matchPre + "DeltaPhi_vs_DeltaEta_" + iEtaStr, 
										    "#Delta#phi vs #Delta#eta (" + iEtaStr +
										    ");#Delta#eta;#Delta#phi", delEtaBin,delEtaLow,delEtaHigh, 
										    delPhiBin,delPhiLow,delPhiHigh);
      // Online-offline quantity correlations
      hist2D[matchPre + "_L1Eta_vs_OffEta_" + iEtaStr]   = iEtaSubDir.make<TH2D>(matchPre + "L1Eta_vs_OffEta_" + iEtaStr, 
										 "L1 #eta vs Offline #eta (" + iEtaStr +
										 ");Offline #eta;L1 #eta", offEtaBin,offEtaLow,offEtaHigh, 
										 offEtaBin,offEtaLow,offEtaHigh);
      hist2D[matchPre + "_L1Phi_vs_OffPhi_" + iEtaStr]   = iEtaSubDir.make<TH2D>(matchPre + "L1Phi_vs_OffPhi_" + iEtaStr, 
										 "L1 #phi vs Offline #phi (" + iEtaStr +
										 ");Offline #phi;L1 #phi", offPhiBin,offPhiLow,offPhiHigh, 
										 offPhiBin,offPhiLow,offPhiHigh);


      // pT Correlations
      hist2D[matchPre + "_L1PT_vs_OffPT_" + iEtaStr]     = iEtaSubDir.make<TH2D>(matchPre + "L1PT_vs_OffPT_" + iEtaStr,"L1 P_{T} vs Offline P_{T} (" + 
										 iEtaStr + ");Offline p_{T} (GeV);L1 P_{T} (GeV)", 
										 offPTbin,offPTlow,offPThigh, offPTbin,offPTlow,offPThigh);
      hist2D[matchPre + "_DeltaPT_vs_OffPT_" + iEtaStr]  = iEtaSubDir.make<TH2D>(matchPre + "DeltaPT_vs_OffPT_" + iEtaStr,"#DeltaP_{T} vs Offline P_{T} (" + 
										 iEtaStr + ");Offline p_{T} (GeV);#DeltaP_{T}", 
										 offPTbin,offPTlow,offPThigh, delPTbin,delPTlow,delPThigh); 
      hist2D[matchPre + "_DeltaEta_vs_OffPT_" + iEtaStr] = iEtaSubDir.make<TH2D>(matchPre + "DeltaEta_vs_OffPT_" + iEtaStr,"#Delta#eta vs Offline P_{T} (" + 
										 iEtaStr + ");Offline p_{T} (GeV);#Delta#eta", 
										 offPTbin,offPTlow,offPThigh, delEtaBin,delEtaLow,delEtaHigh);
      hist2D[matchPre + "_DeltaPhi_vs_OffPT_" + iEtaStr] = iEtaSubDir.make<TH2D>(matchPre + "DeltaPhi_vs_OffPT_" + iEtaStr,"#Delta#phi vs Offline P_{T} (" + 
										 iEtaStr + ");Offline p_{T} (GeV);#Delta#phi", 
										 offPTbin,offPTlow,offPThigh, delPhiBin,delPhiLow,delPhiHigh);
      // NVTX Correlations
      hist2D[matchPre + "_L1PT_vs_NVTX_" + iEtaStr]      = iEtaSubDir.make<TH2D>(matchPre + "L1PT_vs_NVTX_"     + iEtaStr, "L1 P_{T} vs N_{VTX} (" + iEtaStr + 
										 ");N_{VTX};L1 P_{T} (GeV)", PUbin,PUlow,PUhigh, offPTbin,offPTlow,offPThigh);
      hist2D[matchPre + "_DeltaPT_vs_NVTX_" + iEtaStr]   = iEtaSubDir.make<TH2D>(matchPre + "DeltaPT_vs_NVTX_"  + iEtaStr, "#DeltaP_{T} vs N_{VTX} (" + 
										 iEtaStr + ");N_{VTX};#DeltaP_{T}", 
										 PUbin,PUlow,PUhigh, delPTbin,delPTlow,delPThigh);
      hist2D[matchPre + "_DeltaEta_vs_NVTX_" + iEtaStr]  = iEtaSubDir.make<TH2D>(matchPre + "DeltaEta_vs_NVTX_" + iEtaStr, "#Delta#eta vs N_{VTX} (" + 
										 iEtaStr + ");N_{VTX};#Delta#eta", 
										 PUbin,PUlow,PUhigh, delEtaBin,delEtaLow,delEtaHigh);
      hist2D[matchPre + "_DeltaPhi_vs_NVTX_" + iEtaStr]  = iEtaSubDir.make<TH2D>(matchPre + "DeltaPhi_vs_NVTX_" + iEtaStr, "#Delta#phi vs N_{VTX} (" + 
										 iEtaStr + ");N_{VTX};#Delta#phi", 
										 PUbin,PUlow,PUhigh, delPhiBin,delPhiLow,delPhiHigh);
      // Phi Correlations
      hist2D[matchPre + "_L1PT_vs_Phi_" + iEtaStr]        = iEtaSubDir.make<TH2D>(matchPre + "L1PT_vs_Phi_"    + iEtaStr, "L1 P_{T} vs #phi (" + iEtaStr + 
								     ");L1 #phi;L1 P_{T} (GeV)", offPhiBin,offPhiLow,offPhiHigh, offPTbin,offPTlow,offPThigh);
      hist2D[matchPre + "_DeltaPT_vs_Phi_" + iEtaStr]     = iEtaSubDir.make<TH2D>(matchPre + "DeltaPT_vs_Phi"  + iEtaStr, "#DeltaP_{T} vs #phi (" + iEtaStr +
								     ");L1 #phi;#DeltaP_{T}", offPhiBin,offPhiLow,offPhiHigh, delPTbin,delPTlow,delPThigh);
      hist2D[matchPre + "_DeltaEta_vs_Phi_" + iEtaStr]    = iEtaSubDir.make<TH2D>(matchPre + "DeltaEta_vs_Phi" + iEtaStr, "#Delta#eta vs #phi (" + iEtaStr +
								     ");L1 #phi;#Delta#eta", offPhiBin,offPhiLow,offPhiHigh, delEtaBin,delEtaLow,delEtaHigh);
      hist2D[matchPre + "_DeltaPhi_vs_Phi_" + iEtaStr]    = iEtaSubDir.make<TH2D>(matchPre + "DeltaPhi_vs_Phi" + iEtaStr, "#Delta#phi vs #phi (" + iEtaStr +
								     ");L1 #phi;#Delta#phi", offPhiBin,offPhiLow,offPhiHigh, delPhiBin,delPhiLow,delPhiHigh);
      // Eta Correlations
      hist2D[matchPre + "_L1PT_vs_Eta_" + iEtaStr]        = iEtaSubDir.make<TH2D>(matchPre + "L1PT_vs_Eta_"    + iEtaStr, "L1 P_{T} vs #eta (" + iEtaStr + 
								     ");L1 #eta;L1 P_{T} (GeV)", offEtaBin,offEtaLow,offEtaHigh, offPTbin,offPTlow,offPThigh);
      hist2D[matchPre + "_DeltaPT_vs_Eta_" + iEtaStr]     = iEtaSubDir.make<TH2D>(matchPre + "DeltaPT_vs_Eta"  + iEtaStr, "#DeltaP_{T} vs #eta (" + iEtaStr +
								     ");L1 #eta;#DeltaP_{T}", offEtaBin,offEtaLow,offEtaHigh, delPTbin,delPTlow,delPThigh);
      hist2D[matchPre + "_DeltaEta_vs_Eta_" + iEtaStr]    = iEtaSubDir.make<TH2D>(matchPre + "DeltaEta_vs_Eta" + iEtaStr, "#Delta#eta vs #eta (" + iEtaStr +
								     ");L1 #eta;#Delta#eta", offEtaBin,offEtaLow,offEtaHigh, delEtaBin,delEtaLow,delEtaHigh);
      hist2D[matchPre + "_DeltaPhi_vs_Eta_" + iEtaStr]    = iEtaSubDir.make<TH2D>(matchPre + "DeltaPhi_vs_Eta" + iEtaStr, "#Delta#phi vs #eta (" + iEtaStr +
 								     ");L1 #eta;#Delta#phi", offEtaBin,offEtaLow,offEtaHigh, delPhiBin,delPhiLow,delPhiHigh);



      // Calibration TProfiles
      // --------------------------------------------------------------------------------
      hist1D[matchPre + "_OffPT_vs_L1PT_prof_" + iEtaStr]   = iEtaSubDir.make<TProfile>(matchPre + "OffPT_vs_L1PT_prof_" + iEtaStr,
											"Offline P_{T} vs L1 P_{T} (" +
											iEtaStr + ");L1 P_{T} (GeV);Offline p_{T} (GeV)", 
											offPTbin,offPTlow,offPThigh, offPTlow,offPThigh);
      hist1D[matchPre + "_DeltaPT_vs_NVTX_prof_" + iEtaStr] = iEtaSubDir.make<TProfile>(matchPre + "DeltaPT_vs_NVTX_prof_" + iEtaStr,
											"#DeltaP_{T} vs N_{VTX} (" + 
											iEtaStr + ");N_{VTX};#DeltaP_{T}", PUbin,PUlow,PUhigh, 
											delPTlow,delPThigh);
      hist1D[matchPre + "_DeltaPT_vs_Phi_prof_" + iEtaStr]  = iEtaSubDir.make<TProfile>(matchPre + "DeltaPT_vs_Phi_prof_" + iEtaStr,"#DeltaP_{T} vs #phi (" + 
											iEtaStr + ");L1 #phi;#DeltaP_{T}", offPhiBin,offPhiLow,offPhiHigh, 
											delPTlow,delPThigh);


//       hist1D[matchPre + "_DeltaPT_vs_GlobalRho_prof_" + iEtaStr]  = iEtaSubDir.make<TProfile>(matchPre + "DeltaPT_vs_GlobalRho_prof_" + iEtaStr,
// 											     "#DeltaP_{T} vs Global #rho (" + iEtaStr +
// 											     ");Global #rho (GeV);#DeltaP_{T}", 161,-0.25,80.25, 
// 											     delPTlow,delPThigh);
//       hist1D[matchPre + "_DeltaPT_vs_EorHFired_prof_" + iEtaStr]  = iEtaSubDir.make<TProfile>(matchPre + "DeltaPT_vs_EorHFired_prof_" + iEtaStr,
// 											     "#DeltaP_{T} vs Total TTs fired (" + iEtaStr +
// 											     ");TTs fired;#DeltaP_{T}", 1001,-0.5,1000.5, 
// 											     delPTlow,delPThigh);
//       hist1D[matchPre + "_DeltaPT_vs_EplusHtotal_prof_" + iEtaStr] = iEtaSubDir.make<TProfile>(matchPre + "DeltaPT_vs_EplusHtotal_prof_" + iEtaStr,
// 											      "#DeltaP_{T} vs Total TT Energy (" 
// 											      + iEtaStr + ");Total TT Energy (GeV);#DeltaP_{T}", 
// 											      1501,-0.5,1500.5, delPTlow,delPThigh);


     } // End iEta loop


    // ********************************************************************************
    // *                             Calibrated jet plots                             *
    // ********************************************************************************

    TFileDirectory calibDir       = fs->mkdir( "Calibration" );
    hist1D["PrePUS_DeltaPT"] = calibDir.make<TH1F>("PrePUS_DeltaPT", "#Deltap_{T} - PrePUS;#Deltap_{T} (Calibrated - Uncalibrated)/Uncalibrated;Entries", 
						   delPTbin, delPTlow, delPThigh);
    hist2D["PrePUS_DeltaPT_vs_UncorrPt"] = calibDir.make<TH2D>("PrePUS_DeltaPT_vs_UncorrPt", 
							       "#Deltap_{T} vs Uncorrected p_{T} - PrePUS;Uncorrected p_{T} (GeV);#Deltap_{T} (Calibrated - Uncalibrated)/Uncalibrated;Entries", 
							       offPTbin,offPTlow,offPThigh, delPTbin,delPTlow,delPThigh);
    
    hist2D["PrePUS_DeltaPT_vs_iEta"] = calibDir.make<TH2D>("PrePUS_DeltaPT_vs_iEta", 
							   "#Deltap_{T} vs i#eta - PrePUS;i#eta;#Deltap_{T} (Calibrated - Uncalibrated)/Uncalibrated;Entries", 
							   57,-28.5,28.5, delPTbin,delPTlow,delPThigh);
    hist2D["PrePUS_DeltaPT_vs_NVTX"] = calibDir.make<TH2D>("PrePUS_DeltaPT_vs_NVTX",
                                                           "#Deltap_{T} vs N_{VTX} - PrePUS;N_{VTX};#Deltap_{T} (Calibrated - Uncalibrated)/Uncalibrated;Entries",
                                                           PUbin,PUlow,PUhigh, delPTbin,delPTlow,delPThigh);




  } // End matched jets


}









// **********************************************************************
// *                              analyze()                             *
// *                                                                    *
// *                                                                    *
// **********************************************************************


void AnalyzeJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){


  bool evValid = true;


//   auto_ptr<int>                      outputNVTX( new int() );
//   auto_ptr<double>                   outputprePusHT_Tower( new double() );
//   auto_ptr<double>                   outputHT_Tower(  new double() );
//   auto_ptr<double>                   outputMHT_Tower( new double() );
//   auto_ptr< vector<TLorentzVector> > outputPUS_Tower_TLorentz( new vector<TLorentzVector>() );




 

  PRINT("Handles")

  // ************************************************
  // *                     HT/MHT                   *
  // ************************************************
#ifdef CURRENT
  //Current HT&MHT
  edm::Handle<l1extra::L1EtMissParticleCollection> L1MHT_curr;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("L1extraMHT"),L1MHT_curr);
  if(!L1MHT_curr.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("L1extraMHT") << std::endl;

  }
#endif 
  //calibrated tower MHT/HT collection: L1EtMissParticleCollection
  edm::Handle<l1extra::L1EtMissParticleCollection> L1MHT_up;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("TowerJetMHT"),L1MHT_up);
  if(!L1MHT_up.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("TowerJetMHT") << std::endl;
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



  //calibrated tower jet collection: l1extra:: L1JetParticleCollection
  edm::Handle<l1extra::L1JetParticleCollection> CalibJets_L1extra;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibCircle8l1extra"), CalibJets_L1extra);
  if(!CalibJets_L1extra.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibCircle8l1extra") << std::endl;
  }

  #ifdef CURRENT
  //Current L1 extra jets

  // Get L1 jets to be used, currently Cental and Tau jets
  vector <InputTag> l1extraparticles = conf_.getParameter< vector < InputTag > >("extrajet"); 
  // Loop through Central and tau jets
  for (uint i = 0; i < l1extraparticles.size(); ++i){
    Handle<l1extra::L1JetParticleCollection> currl1Jets;
    iEvent.getByLabel(l1extraparticles[i], currl1Jets );
    if(!currl1Jets.isValid()){
      evValid = false;
      edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("extrajet") << std::endl;
    }
  }
  #endif



// #ifdef TT_DIST
//   // Rho information for the tower jets
//   edm::Handle<double> CaloRho_Tower;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CaloRho"), CaloRho_Tower);
//   if(!CaloRho_Tower.isValid()){
//     evValid = false;
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CaloRho") << std::endl;
//   }
// #endif


  // ************************************************************
  // *                        Tower Jets                        *
  // ************************************************************
  SUBPRINT("Tower Jets")

  // Upgrade tower jets prior to PU subtraction
  edm::Handle<l1slhc::L1TowerJetCollection> ProtoJets_Tower;
  edm::Handle<l1slhc::L1TowerJetCollection> FilteredCentralityJets_Tower;
  edm::Handle<l1slhc::L1TowerJetCollection> Filtered1DJets_Tower;
  edm::Handle<l1slhc::L1TowerJetCollection> PrePUSubUpgrCenJet_Tower;
  
  // Don't include handles that don't exist in pre Mk1 ntuples
  if (!mUseOldNtuple){
    //#ifdef PROTO_JET
    // Proto jets (TowerJetProducer)
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("ProtoJets"), ProtoJets_Tower);
    if(!ProtoJets_Tower.isValid()){
      evValid = false;
      edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("ProtoJets_Tower") << std::endl;
    }
    //#endif

    // FilteredCentralityJets
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCentralityJet"), FilteredCentralityJets_Tower);
    if(!FilteredCentralityJets_Tower.isValid()){
      evValid = false;
      edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("FilteredCentralityJets_Tower") << std::endl;
    }

    // Filtered1DJets
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("Filtered1DJets"), Filtered1DJets_Tower);
    if(!Filtered1DJets_Tower.isValid()){
      evValid = false;
      edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("Filtered1DJets_Tower") << std::endl;
    }

    // Pre-PU subtracted jets
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PrePUSubUpgrCenJet"), PrePUSubUpgrCenJet_Tower);
    if(!PrePUSubUpgrCenJet_Tower.isValid()){
      evValid = false;
      edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PrePUSubUpgrCenJet") << std::endl;
    }
    
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




  edm::Handle<double> GlobalRho_Tower;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("GlobalRho"), GlobalRho_Tower);
  if(!GlobalRho_Tower.isValid()){
    evValid = false;
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("GlobalRho") << std::endl;
  }



  
//   // TT collection
//   SUBPRINT("Trigger towers")
//   edm::Handle<l1slhc::L1CaloTowerCollection> caloTowers;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalorimeterTowers"), caloTowers);

//   //  std::cout << conf_.getParameter<edm::InputTag>("CalorimeterTowers").process() << "\n\n\n";

//   if(!caloTowers.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalorimeterTowers") << std::endl;
//     evValid = false;
//   }
  


  //////////////////////////////////
  // Handles to offline information
  //////////////////////////////////
  SUBPRINT("RECO")

  //Need this for information about PU
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RecoVertices"), vertices); 
  if(!vertices.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RecoVertices") << std::endl;
    evValid = false;
  }

  //PU subtracted AK5 calo jets-must be in root file read in
  edm::Handle<reco::CaloJetCollection> PUSAk5Jets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUsubCaloJets"), PUSAk5Jets);
  if(!PUSAk5Jets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUsubCaloJets") << std::endl;
    evValid = false;
  }

  edm::Handle<edm::View< reco::CaloJet > > PrePUSAk5Jets; // uncorrected jets!
  //edm::Handle<reco::CaloJetCollection> PrePUSAk5Jets;
  iEvent.getByLabel("ak5CaloJets", PrePUSAk5Jets );
  if(!PrePUSAk5Jets.isValid()){
    edm::LogWarning("MissingProduct") << "ak5CaloJets"<< std::endl;
    evValid = false;
  }



//   // kt6 jets
//   edm::Handle<reco::CaloJetCollection> kt6Jets;
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("kt6Jets"), kt6Jets);
//   if(!kt6Jets.isValid()){
//     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("kt6Jets") << std::endl;
//     evValid = false;
//   }

  //information about offline calo rho-must be in root file read in
  edm::Handle<double> rhoCALOjets;  
  //  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RhoCaloJets"), rhoCALOjets);
  iEvent.getByLabel("ak5CaloJets","rho", rhoCALOjets);
  if(!rhoCALOjets.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RhoCaloJets") << std::endl;
    evValid = false;
  }


  

  // kt6 CaloJet rho
  edm::Handle<double> kt6CaloRho;
  iEvent.getByLabel("kt6CaloJets","rho", kt6CaloRho);
  //iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RhoKt6CaloJets"), kt6CaloRho);
  if(!kt6CaloRho.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RhoKt6CaloJets") << std::endl;
    evValid = false;
  }

  


  // handle to the jet ID variables
  edm::Handle<reco::JetIDValueMap> hJetIDMap;
  iEvent.getByLabel( conf_.getParameter<edm::InputTag>("jetIDHelperConfig") , hJetIDMap );
  if(!hJetIDMap.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("jetIDHelperConfig")<< std::endl;
    evValid = false;
  }









  // ------------------------------------------------------------------------------------------------------------------------------------------------------
  // -                                                              Fill Histograms                                                                       -
  // ------------------------------------------------------------------------------------------------------------------------------------------------------


 
  if(!evValid){
    std::cout << "Invalid event\n";	
  }
  else{   // valid event
    

    // Extract the number of reconstructed AK5 vertices
    //
    //
    // TODO: CHECK THESE ARE PROPERLY CLEANED
    //
    //
    int    NVTX      = vertices->size();
    //    *outputNVTX = NVTX;


    // TODO: ADD THESE PARAMETERS
//     double prePusHT_Tower, HT_Tower, MHT_Tower;
//     *outputprePusHT_Tower = prePusHT_Tower;
//     *outputHT_Tower       = HT_Tower;
//     *outputMHT_Tower      = MHT_Tower;






//     // ****************************************************************************************************
//     // *                                      Vertex distributions                                        *
//     // ****************************************************************************************************

//     #ifdef VERTEX
    

//     TString vtxPre = "Vertex";

//     for ( reco::VertexCollection::const_iterator vertex = vertices->begin(); vertex != vertices->end(); ++vertex ) {


//       double vtxChi2         = vertex->chi2();
//       int    vtxNdof         = vertex->ndof();
//       //      int    vtxTracks = vertex->tracks()->size();
//       //      std::cout << "Chi2 = " << vtxChi2 << "\tNdof = " << vtxNdof << "\tChi2/ndof" << vtxChi2/vtxNdof << "\n";

//       hist1D[vtxPre + "_Chi2"]          ->Fill(vtxChi2);
//       hist1D[vtxPre + "_Ndof"]          ->Fill(vtxNdof);
//       if (vtxNdof != 0) // Avoid division by zero 
// 	hist1D[vtxPre + "_Chi2_over_Ndof"]->Fill(vtxChi2/vtxNdof);
//       hist2D[vtxPre + "_Ndof_vs_Chi2"]  ->Fill(vtxChi2, vtxNdof);

//     }

//     #endif






    // Ht pre-PUS
    double prePusHT = 0;
    for (l1slhc::L1TowerJetCollection::const_iterator prePus_It = PrePUSubUpgrCenJet_Tower->begin(); prePus_It != PrePUSubUpgrCenJet_Tower->end(); ++prePus_It ){
      prePusHT += prePus_It->Pt();
    }

    //    *outputprePusHT_Tower = prePusHT;






    
    
//     // ****************************************************************************************************
//     // *                                        TT distributions                                          *
//     // ****************************************************************************************************

    
// #ifdef TT_DIST
//     PRINT("TT distributions")
    
//       //      std::cout << "\n\niPhi\tiEta\tE+H\n";


//     //    l1slhc::L1CaloTowerCollection> caloTowers;
//     TString TTPre = "TT";

//     // Maximum event TT energies
//     int maxE(0), maxH(0), totalE(0), totalH(0);
//     // Number of TTs that fired
//     int EFired(0), HFired(0), EorHFired(0), EandHFired(0);

//     int TTsBelow5GeV(0);

//     // Clear the array from the previous run
//     ttThresholdFired.clear();
//     ttThresholdFired.resize( ttThreshold.size() - 1 );

//     for( l1slhc::L1CaloTowerCollection::const_iterator lTT_It = caloTowers->begin() ; 
// 	 lTT_It != caloTowers->end() ; ++lTT_It ){

//       // ****************************************
//       // *    Load the calorimeter tower data   *
//       // ****************************************
//       int E      = lTT_It->E();
//       int H      = lTT_It->H();
//       totalE += E;
//       totalH += H;
//       int EplusH = E + H;
//       int iEta   = lTT_It->iEta();
//       int iPhi   = lTT_It->iPhi();
//       double Eta = mTowerGeo.eta(iEta);
//       double Phi = mTowerGeo.phi(iPhi);
// //       int EcalFG = lTT_It->EcalFG();
// //       int HcalFG = lTT_It->HcalFG();

//       // Restrict to central TTs
//       if (abs(iEta) > 28)
// 	continue;
      
//       // Store the maximum TT E and H
//       if (E > maxE)
// 	maxE = E;
//       if (H > maxH)
// 	maxH = H;

// //       double EoverH(0);
// //       if (H != 0){
// // 	EoverH = double(E)/double(H);
// //       }

// //       if ( abs(iEta) > 20)
// // 	std::cout << iPhi << "\t" << iEta << "\t" << EplusH << "\n";


//       //<BUG?????>
//       if (Phi > PI)
// 	Phi -= 2*PI;
//       //<BUG?????>

      
//       int EcalFired(0), HcalFired(0); //, EcalOrHcalFired(0), EcalAndHcalFired(0);
//       if (E > 0){
// 	EcalFired = 1;
// 	EFired++;
//       }
//       if (H > 0){
// 	HcalFired = 1;
// 	HFired++;
//       }
//       if ( (EcalFired + HcalFired) > 0 ){
// 	//	EcalOrHcalFired = 1;
// 	EorHFired++;
// 	if ( (EcalFired + HcalFired) == 2 ){
// 	  //	  EcalAndHcalFired = 1;
// 	  EandHFired++;
// 	}
//       }


//       if ( EplusH < 5)
// 	TTsBelow5GeV++;

// //       hist2D[TTPre + "_Phi_vs_iPhi"]->Fill( iPhi, Phi );
// //       hist2D[TTPre + "_Eta_vs_iEta"]->Fill( iEta, Eta ); 

//       hist2D[TTPre + "_E-Phi_vs_Eta"]         ->Fill( Eta, Phi, E );
//       hist2D[TTPre + "_H-Phi_vs_Eta"]         ->Fill( Eta, Phi, H );
//       hist2D[TTPre + "_EplusH-Phi_vs_Eta"]    ->Fill( Eta, Phi, EplusH );

//       //      hist2D[TTPre + "_EHratio-Phi_vs_Eta"]   ->Fill( Eta, Phi, EoverH );
// //       hist2D[TTPre + "_EFG-Phi_vs_Eta"]       ->Fill( Eta, Phi, EcalFG );
// //       hist2D[TTPre + "_HFG-Phi_vs_Eta"]       ->Fill( Eta, Phi, HcalFG );
//       hist2D[TTPre + "_EFired-Phi_vs_Eta"]    ->Fill( Eta, Phi, EFired ); 
//       hist2D[TTPre + "_HFired-Phi_vs_Eta"]    ->Fill( Eta, Phi, HFired ); 
//       hist2D[TTPre + "_EorHFired-Phi_vs_Eta"] ->Fill( Eta, Phi, EorHFired ); 
//       hist2D[TTPre + "_EandHFired-Phi_vs_Eta"]->Fill( Eta, Phi, EandHFired ); 


//       hist2D[TTPre + "_E-iPhi_vs_iEta"]         ->Fill( iEta, iPhi, E );
//       hist2D[TTPre + "_H-iPhi_vs_iEta"]         ->Fill( iEta, iPhi, H );
//       hist2D[TTPre + "_EplusH-iPhi_vs_iEta"]    ->Fill( iEta, iPhi, EplusH );

//       // Profile histograms
//       hist1D[TTPre + "_E-iEta_prof"]         ->Fill( iEta, E );
//       hist1D[TTPre + "_H-iEta_prof"]         ->Fill( iEta, H );
//       hist1D[TTPre + "_EplusH-iEta_prof"]    ->Fill( iEta, EplusH );
//       hist1D[TTPre + "_E-iPhi_prof"]         ->Fill( iPhi, E );
//       hist1D[TTPre + "_H-iPhi_prof"]         ->Fill( iPhi, H );
//       hist1D[TTPre + "_EplusH-iPhi_prof"]    ->Fill( iPhi, EplusH );

//       //      hist2D[TTPre + "_EHratio-iPhi_vs_iEta"]   ->Fill( iEta, iPhi, EoverH ); 
//       hist2D[TTPre + "_EFired-iPhi_vs_iEta"]    ->Fill( iEta, iPhi, EFired );
//       hist2D[TTPre + "_HFired-iPhi_vs_iEta"]    ->Fill( iEta, iPhi, HFired ); 
//       hist2D[TTPre + "_EorHFired-iPhi_vs_iEta"] ->Fill( iEta, iPhi, EorHFired ); 
//       hist2D[TTPre + "_EandHFired-iPhi_vs_iEta"]->Fill( iEta, iPhi, EandHFired );  


//       // iEta binned distributions
//       TString iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEta)));
//       hist1D[TTPre + "_E_" + iEtaStr]     ->Fill(E);
//       hist1D[TTPre + "_H_" + iEtaStr]     ->Fill(H);
//       hist1D[TTPre + "_EplusH_" + iEtaStr]->Fill(EplusH);


//       //
//       // CURRENTLY BROKEN
//       //    
//       bool ttBinned = false;
//       // TT threshold plots
//       for (unsigned int ttI = 0; ttI < ttThreshold.size(); ttI++){
	
// 	// get TT threshold
// 	int ttThresh   = ttThreshold[ttI];
	
// 	// TT energy is exceeded by the threshold
// 	if ( (EplusH <= ttThresh) && (ttI != 0) ){
	  
// 	  //	  std::cout << "E+H = " << EplusH << "\tThreshold = " << ttThresh << "\tCurrent = " << ttThresholdFired[ttI - 1] << "\n";
// 	  ttThresholdFired[ttI - 1]++;
// 	  ttBinned = true;
// 	  break;
//   	}
	
//       }
//       if (ttBinned == false){
// 	ttThresholdFired[ttThresholdFired.size() - 1]++;
//       }

//     } // End TT loop
	
//     //    std::cout << "\nNEW EVENT\n";

//     // Fill the maximum TT energies in the event
//     hist1D[TTPre + "_MaxE"]->Fill(maxE);
//     hist1D[TTPre + "_MaxH"]->Fill(maxH);
//     // Fill TT multiplicities
//     hist1D[TTPre + "_ETowersFired"]    ->Fill( EFired );
//     hist1D[TTPre + "_HTowersFired"]    ->Fill( HFired );
//     hist1D[TTPre + "_EorHTowersFired"] ->Fill( EorHFired );
//     hist1D[TTPre + "_EandHTowersFired"]->Fill( EandHFired );


//     hist1D[TTPre + "_Etotal"]     ->Fill( totalE );
//     hist1D[TTPre + "_Htotal"]     ->Fill( totalH );
//     hist1D[TTPre + "_EplusHtotal"]->Fill( totalE + totalH );
    
//     // ************************************************************
//     // *                       RECO kt6 rho                       *
//     // ************************************************************
//     hist2D["GlobalRho_vs_kt6Rho"]  ->Fill( *kt6CaloRho, *GlobalRho_Tower );
//     hist2D["EFired_vs_kt6Rho"]     ->Fill( *kt6CaloRho, EFired );
//     hist2D["HFired_vs_kt6Rho"]     ->Fill( *kt6CaloRho, HFired );
//     hist2D["TTFired_vs_kt6Rho"]    ->Fill( *kt6CaloRho, EorHFired );
//     hist2D["Etotal_vs_kt6Rho"]     ->Fill( *kt6CaloRho, totalE );
//     hist2D["Htotal_vs_kt6Rho"]     ->Fill( *kt6CaloRho, totalH );
//     hist2D["EplusHtotal_vs_kt6Rho"]->Fill( *kt6CaloRho, totalE + totalH );

//     // Inter-quantity correlations
//     hist2D["EplusHtotal_vs_TTFired"]->Fill( EorHFired, totalE + totalH );
//     hist2D["EFired_vs_HFired"]      ->Fill( HFired, EFired );
//     hist2D["ETotal_vs_EFired"]      ->Fill( EFired, totalE );
//     hist2D["HTotal_vs_HFired"]      ->Fill( HFired, totalH );
// //     hist2D["NJets_vs_EplusHtotal"]  ->Fill( totalE + totalH, upgradeJets.size() );
// //     hist2D["NJets_vs_TTFired"]      ->Fill( EorHFired, upgradeJets.size() );
// //  hist1D["NJets_vs_NVTX_prof"]

//     hist2D["TTFired_vs_NVTX"]         ->Fill( NVTX, EorHFired );
//     hist2D["TTFiredBelow5GeV_vs_NVTX"]->Fill( NVTX, TTsBelow5GeV );

//     // Profiles
//     hist1D["GlobalRho_vs_kt6Rho_prof"]  ->Fill( *kt6CaloRho, *GlobalRho_Tower );
//     hist1D["EFired_vs_kt6Rho_prof"]     ->Fill( *kt6CaloRho, EFired );
//     hist1D["HFired_vs_kt6Rho_prof"]     ->Fill( *kt6CaloRho, HFired );
//     hist1D["TTFired_vs_kt6Rho_prof"]    ->Fill( *kt6CaloRho, EorHFired );
//     hist1D["Etotal_vs_kt6Rho_prof"]     ->Fill( *kt6CaloRho, totalE );
//     hist1D["Htotal_vs_kt6Rho_prof"]     ->Fill( *kt6CaloRho, totalH );
//     hist1D["EplusHtotal_vs_kt6Rho_prof"]->Fill( *kt6CaloRho, totalE + totalH );
//     hist1D["NVTX_vs_kt6Rho_prof"]       ->Fill( *kt6CaloRho, NVTX );   

//     hist1D["GlobalRho_vs_NVTX_prof"]    ->Fill( NVTX, *GlobalRho_Tower );   
//     hist1D["EFired_vs_NVTX_prof"]       ->Fill( NVTX, EFired );      
//     hist1D["HFired_vs_NVTX_prof"]       ->Fill( NVTX, HFired );      
//     hist1D["TTFired_vs_NVTX_prof"]      ->Fill( NVTX, EorHFired );
//     hist1D["Etotal_vs_NVTX_prof"]       ->Fill( NVTX, totalE );      
//     hist1D["Htotal_vs_NVTX_prof"]       ->Fill( NVTX, totalH );      
//     hist1D["EplusHtotal_vs_NVTX_prof"]  ->Fill( NVTX, totalE + totalH ); 


//     // Inter-quantity correlations
//     hist1D["EplusHtotal_vs_TTFired_prof"]->Fill( EorHFired, totalE + totalH );
//     hist1D["EFired_vs_HFired_prof"]      ->Fill( HFired, EFired );
//     hist1D["ETotal_vs_EFired_prof"]      ->Fill( EFired, totalE );
//     hist1D["HTotal_vs_HFired_prof"]      ->Fill( HFired, totalH );
// //     hist1D["NJets_vs_EplusHtotal_prof"]  ->Fill( totalE + totalH, upgradeJets.size() );
// //     hist1D["NJets_vs_TTFired_prof"]      ->Fill( EorHFired, upgradeJets.size() );

 
//     // Fill TT threshold plots
//     for (unsigned int ttI = 0; ttI < ttThreshold.size(); ttI++){
      
//       // get TT threshold
//       int ttThresh   = ttThreshold[ttI];  
//       TString ttThreshStr = Form("%d",Int_t(ttThresh));
//       int ttMultiplicity = ttThresholdFired[ttI];
//       hist2D["TTFired_vs_kt6Rho_TT_gt_" + ttThreshStr]          ->Fill(*kt6CaloRho, ttMultiplicity );
//       hist1D["TTFired_vs_kt6Rho_TT_gt_" + ttThreshStr + "_prof"]->Fill(*kt6CaloRho, ttMultiplicity );

//     }

// #endif






  
    
    PRINT("Tower Jet distributions")


      // *****************************************************************************
      // *                 ProtoJets - All jet candidates (0.1 GeV threshold)        *
      // *****************************************************************************


#ifdef PROTO_JET
	SUBPRINT("Proto")
            fillTowerJetHistograms(ProtoJets_Tower, "Proto");
#endif
    
      // *****************************************************************************
      // *                     Filtered centrality Jets                              *
      // *****************************************************************************
      #ifdef FILT_CENT
      SUBPRINT("FiltCent")
	fillTowerJetHistograms(FilteredCentralityJets_Tower, "FiltCent");
      #endif

      // *****************************************************************************
      // *                 1D Filtered Jets - Eta filtered                           *
      // *****************************************************************************
      #ifdef FILT_1D
	SUBPRINT("Filt1D")
            fillTowerJetHistograms(Filtered1DJets_Tower, "1DFilt");
      #endif 

      // *****************************************************************************
      // *                 Pre Pu-subtracted L1 Towerjet collection                  *
      // *****************************************************************************
	SUBPRINT("PrePUS")	
      fillTowerJetHistograms(PrePUSubUpgrCenJet_Tower, "PrePUSub");


    

    // *****************************************************************************
    // *                Global Pu-subtracted L1 Towerjet collection                *
    // *****************************************************************************
    SUBPRINT("PUS")
    fillTowerJetHistograms(UpgrCenJet_Tower,  "PUSub");


    // *****************************************************************************
    // *                 Local Pu-subtracted L1 Towerjet collection                *
    // *****************************************************************************
    SUBPRINT("LPUS")
    fillTowerJetHistograms(LocalPUS_Tower, "LPUSub");


    // *****************************************************************************
    // *                          PUS-PrePUS correlations                          *
    // *****************************************************************************

    PRINT("Tower Jet Correlations")
#ifdef PROTO_JET
    fillL1CorrelationHistograms( ProtoJets_Tower, Filtered1DJets_Tower, "Proto_1DFilt", NVTX );

    fillL1CorrelationHistograms( ProtoJets_Tower, PrePUSubUpgrCenJet_Tower, "Proto_PrePUS", NVTX );

    fillL1CorrelationHistograms( ProtoJets_Tower, FilteredCentralityJets_Tower, "Proto_FiltCent", NVTX );
#endif
#ifdef FILT_CENT
    fillL1CorrelationHistograms( FilteredCentralityJets_Tower, PrePUSubUpgrCenJet_Tower, "FiltCent_PrePUS", NVTX );
#endif
#ifdef FILT_1D
    fillL1CorrelationHistograms( Filtered1DJets_Tower,     PrePUSubUpgrCenJet_Tower, "1DFilt_PrePUS", NVTX );
#endif
    fillL1CorrelationHistograms( PrePUSubUpgrCenJet_Tower, UpgrCenJet_Tower, "PrePUS_PUS", NVTX );

    fillL1CorrelationHistograms( PrePUSubUpgrCenJet_Tower, LocalPUS_Tower,  "PrePUS_LPUS", NVTX );

    fillL1CorrelationHistograms( UpgrCenJet_Tower,         LocalPUS_Tower,     "PUS_LPUS", NVTX );



    // ************************************************************************************************************************
    // *                                                  Calibrated jets                                                     *
    // ************************************************************************************************************************






    // *************************************************************************************
    // *                                        AK5                                        *
    // *************************************************************************************

    PRINT("RECO")
    hist1D["NVTX"]->Fill(NVTX);


    //Fill a vector of TLorentzVectors with offline PU corrected AK5 calo jets
        vector <TLorentzVector> offlineJets, offlineUncorrJets;//, kt6JetVec;
    //    TLorentzVector off;

    //    unsigned int idx;


//     // ****************************************
//     // *                    kt6
//     // ****************************************

//     vector<double> kt6JetPtAreaRatio;
//     vector< vector <double> > kt6LocalJetPtAreaRatio;
//     vector< int > kt6LocalNJets;
//     // Create a container for each eta slice to store the jet rho in order to calculate local rhos, don't need last bin
//     kt6LocalJetPtAreaRatio.resize(  PUetaSlice.size() - 1 );
//     kt6LocalNJets.resize(  PUetaSlice.size() - 1 );


//     // NOTE: THESE ARE NOT CORRECTED FOR PU
//     for (reco::CaloJetCollection::const_iterator kt6It = kt6Jets->begin();kt6It != kt6Jets->end(); ++kt6It ){

//       off.SetPtEtaPhiM(kt6It->p4().Pt(),kt6It->p4().eta(),kt6It->p4().phi(),kt6It->p4().M());
//       kt6JetVec.push_back(off);

//       // All availible quantities can be found here:
//       //    https://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_4_4/doc/html/d6/d7a/classreco_1_1Jet.html
//       double Pt          = kt6It->p4().Pt();
//       double JetRealArea = kt6It->jetArea(); 
//       double WeightedEta = kt6It->p4().Eta();
//       double WeightedPhi = kt6It->p4().Phi();
//       // ratio of the jet pt to its area
//       double ptAreaRatio = Pt/JetRealArea;
      
      

//       // Store jet energy densities for calculation of global rho
//       kt6JetPtAreaRatio.push_back( ptAreaRatio );

      
//       // ************************************************************
//       // *                         Local rho                        *
//       // ************************************************************

//       for (unsigned int iSlice = 0;iSlice < PUetaSlice.size() - 1; iSlice++){
	
// 	// get eta slice range
// 	double PUetaLow  = PUetaSlice[iSlice];
// 	double PUetaHigh = PUetaSlice[iSlice + 1];

	  
// 	// fill the respective eta bin
// 	if ( (WeightedEta >= PUetaLow) && (WeightedEta < PUetaHigh) ){
	    
// 	  // get eta slice range
// 	  double PUetaLow  = PUetaSlice[iSlice];
// 	  double PUetaHigh = PUetaSlice[iSlice + 1];
	  
// 	  TString PUetaLowStr  = Form("%1.3f", PUetaLow);
// 	  TString PUetaHighStr = Form("%1.3f", PUetaHigh);
// 	  TString PUetaLabel = "eta(" + PUetaLowStr + "-" + PUetaHighStr + ")";
	  
// 	  TString jetPrefix = "kt6";

// 	  // Fill jet distributions
// 	  hist1D[jetPrefix + "_Eta_"    + PUetaLabel]->Fill(WeightedEta);
// 	  hist1D[jetPrefix + "_Phi_"    + PUetaLabel]->Fill(WeightedPhi);
// 	  hist1D[jetPrefix + "_PT_"     + PUetaLabel]->Fill(Pt);
// 	  hist1D[jetPrefix + "_JetRho_" + PUetaLabel]->Fill(ptAreaRatio);
	  

// 	  // Fill the jet in respective vector for calculating rho (lower edge)
// 	  kt6LocalJetPtAreaRatio[iSlice].push_back(ptAreaRatio);
// 	  kt6LocalNJets[iSlice]++;

// 	}
//       }
      

//     }

//     // ************************************************************
//     // *                       Calculate rho                      *
//     // ************************************************************

//     // ********************
//     // *      Global      *
//     // ********************
//     double kt6GlobalRho = Median(kt6JetPtAreaRatio);

//     hist2D["kt6RhoEst_vs_kt6Rho"]     ->Fill(*kt6CaloRho, kt6GlobalRho );
//     hist1D["kt6RhoEst_vs_kt6Rho_prof"]->Fill(*kt6CaloRho, kt6GlobalRho );

//     //    std::cout << "kt6Rho = " << *kt6CaloRho << "\tEstimated Rho = " << kt6GlobalRho << "\n\n";


//     // ********************
//     // *      Local       *
//     // ********************

//     // Calculate rho for bins of eta
//     for (unsigned int iSlice = 0; iSlice < kt6LocalJetPtAreaRatio.size(); iSlice++){
      
//       double etaBinnedRho = Median(kt6LocalJetPtAreaRatio[iSlice]);
//       //      double deltaRhoMin, deltaRhoMax;
      
//       Int_t etaBinnednJets = kt6LocalNJets[iSlice];
      
//       // get eta slice range
//       double PUetaLow  = PUetaSlice[iSlice];
//       double PUetaHigh = PUetaSlice[iSlice + 1];
//       TString PUetaLowStr  = Form("%1.3f", PUetaLow);
//       TString PUetaHighStr = Form("%1.3f", PUetaHigh);
//       TString PUetaLabel = "eta(" + PUetaLowStr + "-" + PUetaHighStr + ")"; 

//       TString jetPrefix = "kt6";
  
//       // Fill the local rho
//       hist1D[jetPrefix + "_Rho_" + PUetaLabel]->Fill(etaBinnedRho);  
//       hist1D[jetPrefix + "_NJETS_"  + PUetaLabel]->Fill(etaBinnednJets);

//       // Profile
//       //      hist1D["LocalRho_vs_kt6LRho_prof_" + PUetaLabel]                 ;
//       hist1D["kt6Rho_vs_kt6LRho_prof_" + PUetaLabel]     ->Fill( etaBinnedRho, *kt6CaloRho );
//       hist1D["TTFired_vs_kt6LRho_prof_" + PUetaLabel]    ->Fill( etaBinnedRho, EorHFired );
//       hist1D["EFired_vs_kt6LRho_prof_" + PUetaLabel]     ->Fill( etaBinnedRho, EFired );
//       hist1D["HFired_vs_kt6LRho_prof_" + PUetaLabel]     ->Fill( etaBinnedRho, HFired );
//       hist1D["Etotal_vs_kt6LRho_prof_" + PUetaLabel]     ->Fill( etaBinnedRho, totalE );
//       hist1D["Htotal_vs_kt6LRho_prof_" + PUetaLabel]     ->Fill( etaBinnedRho, totalH );
//       hist1D["EplusHtotal_vs_kt6LRho_prof_" + PUetaLabel]->Fill( etaBinnedRho, totalE + totalH );
//       hist1D["HtPrePUS_vs_kt6LRho_prof_" + PUetaLabel]   ->Fill( etaBinnedRho, prePusHT );
// //       hist1D["Ht_vs_kt6LRho_prof_" + PUetaLabel]         ->Fill( etaBinnedRho, UP_HT);
// //       hist1D["MHt_vs_kt6LRho_prof_" + PUetaLabel]        ->Fill( etaBinnedRho, UP_MHT);
//       hist1D["LNJets_vs_kt6LRho_prof_" + PUetaLabel]     ->Fill( etaBinnedRho, etaBinnednJets );

      
//     }	
    
    





    // ****************************************
    // *                    ak5
    // ****************************************

    //L2L3 PU corrected jet collection
    //    reco::CaloJetCollection::const_iterator it = PUSAk5Jets->begin();








//     //L2L3 PU corrected PUS jet collection
//     for (reco::CaloJetCollection::const_iterator PUSak5_It = PUSAk5Jets->begin(); PUSak5_It != PUSAk5Jets->end(); ++PUSak5_It ){
      
//       //set the TLorentz vector with uncorrected energy
//       off.SetPtEtaPhiM(PUSak5_It->p4().Pt(),PUSak5_It->p4().eta(),PUSak5_It->p4().phi(),PUSak5_It->p4().M());

//       // Jet quality cuts
//       if( off.Pt() < MIN_OFF_PT ) continue;
//       if( fabs(off.Eta()) > 3.0 ) continue;
     
//       if ( ( fabs(off.Eta()) < 3.0 ) ){
// 	offlineJets.push_back(off);
//       }
     
 
//     }
    
//     // PrePUS ak5
//     //    for (reco::CaloJetCollection::const_iterator PrePUSak5_It = PrePUSAk5Jets->begin(); PrePUSak5_It != PrePUSAk5Jets->end(); ++PrePUSak5_It ){
//     for ( edm::View< reco::CaloJet >::const_iterator PrePUSak5_It = PrePUSAk5Jets->begin(); PrePUSak5_It != PrePUSAk5Jets->end(); ++PrePUSak5_It ){



//       //set the TLorentz vector with uncorrected energy
//       off.SetPtEtaPhiM(PrePUSak5_It->p4().Pt(),PrePUSak5_It->p4().eta(),PrePUSak5_It->p4().phi(),PrePUSak5_It->p4().M());

//       // Jet quality cuts
//       if( off.Pt() < MIN_OFF_PT ) continue;
//       if( fabs(off.Eta()) > 3.0 ) continue;
     
//       if ( ( fabs(off.Eta()) < 3.0 ) ){
// 	offlineUncorrJets.push_back(off);
//       }
    
//     }


//     //uncorrected jet collection for jet ID reference
//     for ( edm::View<reco::CaloJet>::const_iterator ibegin = hJets->begin(),
// 	    iend = hJets->end(), ijet = ibegin; ijet != iend; ++ijet , ++it) {


//       //set the TLorentz vector with uncorrected energy
//       off.SetPtEtaPhiM(it->p4().Pt(),it->p4().eta(),it->p4().phi(),it->p4().M());
      
//       if ( ( fabs(off.Eta()) < 3.0 ) ){
// 	offlineUncorrJets.push_back(off);
//       }
      
//       //set the TLorentz vector with corrected energy
//       off.SetPtEtaPhiM(it->p4().Pt(),it->p4().eta(),it->p4().phi(),it->p4().M());

//       // COMMENT: Magic number
//       //loose jet selections
//       //      if (it->emEnergyFraction() < 0.01) continue;
//       //if (jetId.fHPD>0.98  ) continue;
//       //if (jetId.n90Hits<=1) continue;
      
//       // pT and eta cuts
//       if( off.Pt() < MIN_OFF_PT ) continue;
//       if( fabs(off.Eta()) > 3.0 ) continue;

//       // Fill vector with loose ID's calojets
//       offlineJets.push_back(off);
//     }




    // ****************************************************************************************************
    // *                             Fill (calibrated) jet vectors. Awkward...                            *
    // ****************************************************************************************************


    TLorentzVector up_l1, curr_l1; 
    vector<TLorentzVector> upgradeJetsPrePUS, upgradeJets, upgradeJetsLocal, currentJets;

    // ********************************************************************************
    // *  Iterate through jets, store the Pt, Eta, Phi and M of jets that pass the Pt *
    // *  threshold and Barrel + Endcap acceptance criteria                           *
    // ********************************************************************************

//     // ******************************
//     // *          Upgrade           *
//     // ******************************

//     // Pre PU subtraction

//     //    for (l1extra::L1JetParticleCollection::const_iterator il1 = CalibJets_L1extra->begin(); il1!= CalibJets_L1extra->end() ; ++il1 ){
//     for (l1slhc::L1TowerJetCollection::const_iterator Tower_It = PrePUSubUpgrCenJet_Tower->begin(); Tower_It != PrePUSubUpgrCenJet_Tower->end(); ++Tower_It ){

//       up_l1.SetPtEtaPhiM(Tower_It->p4().Pt(),Tower_It->p4().Eta(),Tower_It->p4().Phi(),Tower_It->p4().M());

//       if( up_l1.Pt() < MIN_L1_PT ) continue;
//       if( fabs(up_l1.Eta()) > 3.0 ) continue;

//       upgradeJetsPrePUS.push_back(up_l1);

//     } // End upgrade jets


//     // Global subtraction

//     //    for (l1extra::L1JetParticleCollection::const_iterator il1 = CalibJets_L1extra->begin(); il1!= CalibJets_L1extra->end() ; ++il1 ){
//     for (l1slhc::L1TowerJetCollection::const_iterator Tower_It = UpgrCenJet_Tower->begin(); Tower_It != UpgrCenJet_Tower->end(); ++Tower_It ){

//       up_l1.SetPtEtaPhiM(Tower_It->p4().Pt(),Tower_It->p4().Eta(),Tower_It->p4().Phi(),Tower_It->p4().M());

//       if( up_l1.Pt() < MIN_L1_PT ) continue;
//       if( fabs(up_l1.Eta()) > 3.0 ) continue;

//       upgradeJets.push_back(up_l1);

//     } // End upgrade jets

//     //    *outputPUS_Tower_TLorentz = upgradeJets;



//     // Local subtraction

//     for (l1slhc::L1TowerJetCollection::const_iterator Tower_It = LocalPUS_Tower->begin(); Tower_It != LocalPUS_Tower->end(); ++Tower_It ){
//       //    for (l1extra::L1JetParticleCollection::const_iterator il1 = CalibJets_L1extra->begin(); il1!= CalibJets_L1extra->end() ; ++il1 ){

//       up_l1.SetPtEtaPhiM(Tower_It->p4().Pt(),Tower_It->p4().Eta(),Tower_It->p4().Phi(),Tower_It->p4().M());

//       if( up_l1.Pt() < MIN_L1_PT ) continue;
//       if( fabs(up_l1.Eta()) > 3.0 ) continue;

//       upgradeJetsLocal.push_back(up_l1);

//     } // End upgrade jets




//     // ******************************
//     // *          Current           *
//     // ******************************

// #ifdef CURRENT
//     for (uint i = 0; i < l1extraparticles.size(); ++i){

//       Handle<l1extra::L1JetParticleCollection> currl1Jets;
//       iEvent.getByLabel(l1extraparticles[i],currl1Jets); 

//       for(l1extra::L1JetParticleCollection::const_iterator itr = currl1Jets->begin();  itr != currl1Jets->end(); ++itr ) {

//         curr_l1.SetPtEtaPhiM(itr->p4().Pt(),itr->p4().eta(),itr->p4().phi(),itr->p4().M());

// 	if( curr_l1.Pt() < MIN_L1_PT ) continue;
//         if( fabs(curr_l1.Eta()) > 3.0 ) continue;

// 	currentJets.push_back(curr_l1);

//       }
//     } // End current jets
// #endif

//     // ************************************************
//     // *               Order jets by pT               *
//     // ************************************************

//     // ak5
//     if( offlineJets.size() != 0 )
//       sort(offlineJets.begin(), offlineJets.end(), sortTLorentz);
//     if( offlineUncorrJets.size() != 0 )
//       sort(offlineUncorrJets.begin(), offlineUncorrJets.end(), sortTLorentz);
//     //upgrade jets
//     if(upgradeJetsPrePUS.size() != 0)
//       sort(upgradeJetsPrePUS.begin(), upgradeJetsPrePUS.end(),sortTLorentz);
//     if(upgradeJets.size() != 0)
//       sort(upgradeJets.begin(), upgradeJets.end(),sortTLorentz);
//     if(upgradeJetsLocal.size() != 0)
//       sort(upgradeJetsLocal.begin(), upgradeJetsLocal.end(),sortTLorentz);

// #ifdef CURRENT
//     //Current L1 Jets
//     if(currentJets.size() != 0)
//       sort(currentJets.begin(), currentJets.end(),sortTLorentz);
// #endif



//     // ****************************************************************************************************
//     // *                                   L1-Offline lead jet matching                                   *
//     // ****************************************************************************************************
//     PRINT("L1 plots")

    

// #ifdef UPGRADE
//       SUBPRINT("UPGRADE")

//     double Up_HT    = L1MHT_up->begin()->etTotal();
//     double Up_MHT   = L1MHT_up->begin()->pt();
//     hist1D["upgr_HT"]   ->Fill(Up_HT);
//     hist1D["upgr_MHT"]  ->Fill(Up_MHT);


//     hist2D["HtPrePUS_vs_kt6Rho"]       ->Fill(*kt6CaloRho, prePusHT);
//     hist1D["HtPrePUS_vs_kt6Rho_prof"]  ->Fill(*kt6CaloRho, prePusHT);
//     hist1D["HtPrePUS_vs_NVTX_prof"]    ->Fill(NVTX, prePusHT);

//     hist2D["Ht_vs_HtPrePUS"]       ->Fill(prePusHT, Up_HT);
//     hist1D["Ht_vs_HtPrePUS_prof"]  ->Fill(prePusHT, Up_HT);


//     hist2D["Ht_vs_kt6Rho"]       ->Fill(*kt6CaloRho, Up_HT);
//     hist2D["MHt_vs_kt6Rho"]      ->Fill(*kt6CaloRho, Up_MHT);
//     hist1D["Ht_vs_kt6Rho_prof"]  ->Fill(*kt6CaloRho, Up_HT);
//     hist1D["MHt_vs_kt6Rho_prof"] ->Fill(*kt6CaloRho, Up_MHT);
//     hist1D["Ht_vs_NVTX_prof"]    ->Fill(NVTX, Up_HT);      
//     hist1D["MHt_vs_NVTX_prof"]   ->Fill(NVTX, Up_MHT);         



//     // Compare upgrade algorithm - Pre PU subtraction
//     if (upgradeJetsPrePUS.size() != 0){
      
//       if ( offlineUncorrJets.size() != 0 ){ 
// 	fillL1OfflineHistograms( upgradeJetsPrePUS, offlineUncorrJets, "upgrPrePUS", NVTX );
//       }
//       fillL1Histograms( upgradeJetsLocal, "upgrPrePUS", NVTX );
      
//     }

//     // Compare upgrade algorithm - Pre PU subtraction
//     if (offlineUncorrJets.size() != 0){
      
//       if ( offlineJets.size() != 0 ){ 
// 	fillL1OfflineHistograms( offlineUncorrJets, offlineJets, "ak5PrePUS", NVTX );
//       }
//       fillL1Histograms( upgradeJetsLocal, "ak5PrePUS", NVTX );
      
//     }




//     // Compare upgrade algorithm - Global subtraction
//     if (upgradeJets.size() != 0){
      
//       if ( offlineJets.size() != 0 ){ 
// 	fillL1OfflineHistograms( upgradeJets, offlineJets, "upgr", NVTX );
//       }
//       fillL1Histograms( upgradeJets, "upgr", NVTX );
      
//     }

//     // Compare upgrade algorithm - Local subtraction
//     if (upgradeJetsLocal.size() != 0){
      
//       if ( offlineJets.size() != 0 ){ 
// 	fillL1OfflineHistograms( upgradeJetsLocal, offlineJets, "upgrLocal", NVTX );
//       }
//       fillL1Histograms( upgradeJetsLocal, "upgrLocal", NVTX );
      
//     }


// #endif
    
// #ifdef CURRENT
//       SUBPRINT("CURRENT")

//     double Curr_HT  = L1MHT_curr->begin()->etTotal();
//     double Curr_MHT = L1MHT_curr->begin()->pt();
//     hist1D["curr_HT"]   ->Fill(Curr_HT);
//     hist1D["curr_MHT"]  ->Fill(Curr_MHT);

//     // Compare current algorithm
//     if (currentJets.size() != 0){
      
//       if ( offlineUncorrJets.size() != 0 ){ 
// 	fillL1OfflineHistograms( currentJets, offlineUncorrJets, "curr", NVTX );
//       }
//       fillL1Histograms( currentJets, "curr", NVTX );

//     }
// #endif

    

 









//     // ********************************************************************
//     // *                       Event distributions                        *
//     // ********************************************************************
//     PRINT("EVENT")

//     // Get the Luminosity block info from the event
//     Int_t Runnr     = iEvent.id().run();
//     Int_t LumiBlock = iEvent.id().luminosityBlock();
//     Int_t Eventnr   = iEvent.id().event();

//     hist1D["Event_Runnr"]    ->Fill(Runnr);
//     hist1D["Event_LumiBlock"]->Fill(LumiBlock);
//     hist1D["Event_Eventnr"]  ->Fill(Eventnr);
    
//     //
//     // Lead jet pt, eta, phi, Ht, mHt, et, met
//     //
//     hist2D["Event_NJETS"]          ->Fill( Eventnr,     upgradeJets.size() );
//     hist2D["NJets_vs_kt6Rho"]      ->Fill( *kt6CaloRho, upgradeJets.size() );
//     hist1D["NJets_vs_kt6Rho_prof"] ->Fill( *kt6CaloRho, upgradeJets.size() );
//     hist1D["NJets_vs_NVTX_prof"]   ->Fill( NVTX,        upgradeJets.size() );     
//     hist2D["Event_NVTX"]           ->Fill( Eventnr, NVTX );

//     //double test =  *GlobalRho_Tower;
//     hist2D["Event_Rho"]     ->Fill( Eventnr, *GlobalRho_Tower);


    
//     hist2D["Event_HT"]   ->Fill( Eventnr, Up_HT );
//     hist2D["Event_MHT"]  ->Fill( Eventnr, Up_MHT );





//     // ********************************************************************
//     // *                    Single event distributions                    *
//     // ********************************************************************

//     eventCount++;    

//     if (eventCount < 11){ // Process the first ten events

       

//       TString eventStr = Form("%d", Eventnr );

//       // Create an event directory
//       TFileDirectory eventDir    = fs->mkdir( "Event" );
//       TFileDirectory eventSubDir = eventDir.mkdir( eventStr.Data() );



//       hist2D["Proto_Pt_Phi_vs_Eta"]       = eventSubDir.make<TH2D>("Proto_Pt_Phi_vs_Eta", "Jet p_{T}, #phi vs #eta - Event " 
// 								   + eventStr + " Proto;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);
//       hist2D["Proto_Pt_iPhi_vs_iEta"]     = eventSubDir.make<TH2D>("Proto_Pt_iPhi_vs_iEta", "Jet p_{T}, i#phi vs i#eta - Event " 
// 								   + eventStr + " Proto;L1 i#eta;L1 i#phi", 
// 								   57,-28.5,28.5, 72,0.5,72.5);     

//       hist2D["FiltCent_Pt_Phi_vs_Eta"]    = eventSubDir.make<TH2D>("FiltCent_Pt_Phi_vs_Eta", "Jet p_{T}, #phi vs #eta - Event " 
// 								   + eventStr + " FiltCent;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);
//       hist2D["FiltCent_Pt_iPhi_vs_iEta"]  = eventSubDir.make<TH2D>("FiltCent_Pt_iPhi_vs_iEta", "Jet p_{T}, i#phi vs i#eta - Event " 
// 								   + eventStr + " FiltCent;L1 i#eta;L1 i#phi", 
// 								   57,-28.5,28.5, 72,0.5,72.5);     

//       hist2D["Filt1D_Pt_Phi_vs_Eta"]      = eventSubDir.make<TH2D>("Filt1D_Pt_Phi_vs_Eta", "Jet p_{T}, #phi vs #eta - Event " 
// 								   + eventStr + " Filt1D;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);
//       hist2D["Filt1D_Pt_iPhi_vs_iEta"]    = eventSubDir.make<TH2D>("Filt1D_Pt_iPhi_vs_iEta", "Jet p_{T}, i#phi vs i#eta - Event " 
// 								   + eventStr + " Filt1D;L1 i#eta;L1 i#phi", 
// 								   57,-28.5,28.5, 72,0.5,72.5);     

//       hist2D["PrePUS_Pt_Phi_vs_Eta"]      = eventSubDir.make<TH2D>("PrePUS_Pt_Phi_vs_Eta", "Jet p_{T}, #phi vs #eta - Event " 
// 								   + eventStr + " PrePUS;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);
//       hist2D["PrePUS_Pt_iPhi_vs_iEta"]    = eventSubDir.make<TH2D>("PrePUS_Pt_iPhi_vs_iEta", "Jet p_{T}, i#phi vs i#eta - Event " 
// 								   + eventStr + " PrePUS;L1 i#eta;L1 i#phi", 
// 								   57,-28.5,28.5, 72,0.5,72.5);     
//       hist2D["PrePUS_AsymPhi_Phi_vs_Eta"] = eventSubDir.make<TH2D>("PrePUS_AsymPhi_Phi_vs_Eta", "Jet A_{#phi}, #phi vs #eta - Event " 
// 								   + eventStr + " PrePUS;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);
//       hist2D["PrePUS_AsymEta_Phi_vs_Eta"] = eventSubDir.make<TH2D>("PrePUS_AsymEta_Phi_vs_Eta", "Jet A_{#phi}, #phi vs #eta - Event " 
// 								   + eventStr + " PrePUS;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);
      
//       hist2D["PUS_Pt_Phi_vs_Eta"]         = eventSubDir.make<TH2D>("PUS_Pt_Phi_vs_Eta", "Jet p_{T}, #phi vs #eta - Event " 
// 								   + eventStr + " PUS;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);
//       hist2D["PUS_Pt_iPhi_vs_iEta"]       = eventSubDir.make<TH2D>("PUS_Pt_iPhi_vs_iEta", "Jet p_{T}, i#phi vs i#eta - Event " 
// 								   + eventStr + " PUS;L1 i#eta;L1 i#phi", 
// 								   57,-28.5,28.5, 72,0.5,72.5);     

//       hist2D["LPUS_Pt_Phi_vs_Eta"]        = eventSubDir.make<TH2D>("LPUS_Pt_Phi_vs_Eta", "Jet p_{T}, #phi vs #eta - Event " 
// 								   + eventStr + " LPUS;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);
//       hist2D["LPUS_Pt_iPhi_vs_iEta"]      = eventSubDir.make<TH2D>("LPUS_Pt_iPhi_vs_iEta", "Jet p_{T}, i#phi vs i#eta - Event " 
// 								   + eventStr + " LPUS;L1 i#eta;L1 i#phi", 
// 								   57,-28.5,28.5, 72,0.5,72.5);     

//       hist2D["Ak5_Pt_Phi_vs_Eta"]         = eventSubDir.make<TH2D>("Ak5_Pt_Phi_vs_Eta", "Jet p_{T}, #phi vs #eta - Event " 
// 								   + eventStr + " Ak5;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);
      
//       hist2D["Curr_Pt_Phi_vs_Eta"]        = eventSubDir.make<TH2D>("Curr_Pt_Phi_vs_Eta", "Jet p_{T}, #phi vs #eta - Event " 
// 								   + eventStr + " Current;L1 #eta;L1 #phi", 
// 								   601,-3.005,3.005, 631,-3.155,3.155);

//       // **************************
//       // *  TT histogram binning  *
//       // **************************
//       const Int_t TTetaBin = 56;
//       const Int_t TTphiBin = 72;

//       // TT eta binning
//       const double TTetaBins[] = {-3.0,-2.65,-2.5,-2.322,-2.172,-2.043,-1.93,-1.83,-1.74,-1.653,
// 				    -1.566,-1.4790,-1.3920,-1.3050,-1.2180,-1.1310,-1.0440,-0.9570,
// 				    -0.8700,-0.7830,-0.6950,-0.6090,-0.5220,-0.4350,-0.3480,-0.2610,
// 				    -0.1740,-0.0870,
// 				    0.,
// 				    0.0870,0.1740,
// 				    0.2610,0.3480,0.4350,0.5220,0.6090,0.6950,0.7830,0.8700,
// 				    0.9570,1.0440,1.1310,1.2180,1.3050,1.3920,1.4790,1.566,
// 				    1.653,1.74,1.83,1.93,2.043,2.172,2.322,2.5,2.65,3.0};

//       // TT phi binning
//       double phiDiv = 2*PI/(TTphiBin);
//       double TTphiBins[TTphiBin + 1];
//       for (Int_t iPhi = 0; iPhi < TTphiBin + 1; iPhi++) {
// 	double phi = iPhi*phiDiv - PI;
// 	TTphiBins[iPhi] = phi;
//       }



//       //Eta-phi space
//       hist2D["Ev_TT_E-Phi_vs_Eta"]       = eventSubDir.make<TH2D>("Ev_TT_E-Phi_vs_Eta","Ecal #phi vs #eta - Event " +
// 							       eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//       hist2D["Ev_TT_H-Phi_vs_Eta"]       = eventSubDir.make<TH2D>("Ev_TT_H-Phi_vs_Eta","Hcal #phi vs #eta - Event " +
// 							       eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//       hist2D["Ev_TT_E+H-Phi_vs_Eta"]     = eventSubDir.make<TH2D>("Ev_TT_E+H-Phi_vs_Eta","Ecal + Hcal #phi vs #eta - Event " +
// 							       eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//       hist2D["Ev_TT_EHratio-Phi_vs_Eta"] = eventSubDir.make<TH2D>("Ev_TT_EHratio-Phi_vs_Eta","Ecal-Hcal energy ratio #phi vs #eta " +
// 							       eventStr + ";#eta;#phi", TTetaBin,TTetaBins, TTphiBin,TTphiBins);
//       //iEta-iPhi space
//       hist2D["Ev_TT_E-iPhi_vs_iEta"]       = eventSubDir.make<TH2D>("Ev_TT_E-iPhi_vs_iEta","Ecal i#phi vs i#eta - Event " +
// 								 eventStr + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
//       hist2D["Ev_TT_H-iPhi_vs_iEta"]       = eventSubDir.make<TH2D>("Ev_TT_H-iPhi_vs_iEta","Hcal i#phi vs i#eta - Event " +
// 								 eventStr + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
//       hist2D["Ev_TT_E+H-iPhi_vs_iEta"]     = eventSubDir.make<TH2D>("Ev_TT_E+H-iPhi_vs_iEta","Ecal + Hcal i#phi vs i#eta - Event " +
// 								 eventStr + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);
//       hist2D["Ev_TT_EHratio-iPhi_vs_iEta"] = eventSubDir.make<TH2D>("Ev_TT_EHratio-iPhi_vs_iEta","Ecal-Hcal energy ratio i#phi vs i#eta " +
// 								 eventStr + ";i#eta;i#phi", 57,-28.5,28.5, 72,0.5,72.5);




//       // ******************************************************
//       // *                  Trigger Towers                    *
//       // ******************************************************
//       SUBPRINT("TT")
      

// 	// Store TT pattern
// 	TMacro * mTT = new TMacro(eventStr + "_TTs", "Event TTs");


// 	for( l1slhc::L1CaloTowerCollection::const_iterator lTT_It = caloTowers->begin() ;
// 	   lTT_It != caloTowers->end() ; ++lTT_It ){

// 	// ****************************************
// 	// *    Load the calorimeter tower data   *
// 	// ****************************************
// 	int E      = lTT_It->E();
// 	int H      = lTT_It->H();
// 	int EplusH = E + H;
// 	int iEta   = lTT_It->iEta();
// 	int iPhi   = lTT_It->iPhi();
// 	double Eta = mTowerGeo.eta(iEta);
// 	double Phi = mTowerGeo.phi(iPhi);
// 	double EoverH(0);
// 	if ( (H != 0) && (E != 0) ){
// 	  EoverH = double(E)/double(H);
// 	}

// 	//	int EcalFG = lTT_It->EcalFG();
// 	//	int HcalFG = lTT_It->HcalFG();

	
// 	// Restrict to central TTs
// 	if (abs(iEta) > 28)
// 	  continue;
	
// 	//<BUG?????>
// 	if (Phi > PI)
// 	  Phi -= 2*PI;
// 	//<BUG?????>


// 	// Store the event TTs in the test pattern format:
// 	// iEta, iPhi, E, H, FG   -  Currently ignoring FG
// 	mTT->AddLine( TString(Form("%d", iEta)) + "\t" + TString(Form("%d", iPhi)) + "\t" + TString(Form("%d", E)) + "\t" + TString(Form("%d", H)) + "\t0");


// 	hist2D["Ev_TT_E-Phi_vs_Eta"]        ->Fill(Eta,Phi,E);
// 	hist2D["Ev_TT_H-Phi_vs_Eta"]        ->Fill(Eta,Phi,H);
// 	hist2D["Ev_TT_E+H-Phi_vs_Eta"]      ->Fill(Eta,Phi,EplusH);
// 	if (EoverH != 0)
// 	  hist2D["Ev_TT_EHratio-Phi_vs_Eta"]->Fill(Eta,Phi,EoverH);
// 	hist2D["Ev_TT_E-iPhi_vs_iEta"]      ->Fill(iEta,iPhi,E);
// 	hist2D["Ev_TT_H-iPhi_vs_iEta"]      ->Fill(iEta,iPhi,H);
// 	hist2D["Ev_TT_E+H-iPhi_vs_iEta"]    ->Fill(iEta,iPhi,EplusH);
// 	hist2D["Ev_TT_EHratio-iPhi_vs_iEta"]->Fill(iEta,iPhi,EoverH);

//       }

// 	// Save the event TTs to file
// 	mTT->Write();

//       // To be implemented - Event-by-event distributions for events meeting certain jet multiplicity criteria

// //       hist2D[jetPrefix + "_Pt-AreaRatio_iEta"] = fs->make<TH2D>("Pt-AreaRatio_iEta","Pt-JetArea Ratio vs i#eta;i#eta PU-Subtracted jets;p_{T}/#Delta#phi#Delta#eta (GeV)",
// //       57,-28.5,28.5, 402,-1,200);
	



// 	// Visualisation
// 	TCanvas* TTcanv = eventSubDir.make<TCanvas>("TTcanv" + eventStr, "", 800,600);
// 	TLatex* latex = new TLatex();
	
// 	hist2D["Ev_TT_E+H-iPhi_vs_iEta"]->Draw("COLZ");
// 	latex->DrawLatex(0.9, 0.1, "TEST TEST TEST");
	


//     }



      if (!mUseOldNtuple){
	


	// ******************************************************
        // *                     Proto jets                     *
        // ******************************************************
        SUBPRINT("Proto")
	for (l1slhc::L1TowerJetCollection::const_iterator Proto_It = ProtoJets_Tower->begin(); Proto_It != ProtoJets_Tower->end();
	     ++Proto_It ){

	  hist2D["Proto_Pt_Phi_vs_Eta"]    ->Fill( Proto_It->WeightedEta(), Proto_It->WeightedPhi(), Proto_It->Pt() );
	  hist2D["Proto_Pt_iPhi_vs_iEta"]  ->Fill( Proto_It->iEta(), Proto_It->iPhi(), Proto_It->Pt() );

	}

	// ******************************************************
        // *              Filtered centrality jets              *
        // ******************************************************
        SUBPRINT("FiltCent")
	for (l1slhc::L1TowerJetCollection::const_iterator FiltCent_It = FilteredCentralityJets_Tower->begin(); 
	     FiltCent_It != FilteredCentralityJets_Tower->end(); ++FiltCent_It ){

	  hist2D["FiltCent_Pt_Phi_vs_Eta"]    ->Fill( FiltCent_It->WeightedEta(), FiltCent_It->WeightedPhi(), FiltCent_It->Pt() );
	  hist2D["FiltCent_Pt_iPhi_vs_iEta"]  ->Fill( FiltCent_It->iEta(), FiltCent_It->iPhi(), FiltCent_It->Pt() );

	}



	// ******************************************************
	// *                  Filtered 1D jets                  *
	// ******************************************************
	SUBPRINT("Filt1D")
        for (l1slhc::L1TowerJetCollection::const_iterator Filt1D_It = Filtered1DJets_Tower->begin(); Filt1D_It != Filtered1DJets_Tower->end();
             ++Filt1D_It ){

	  hist2D["Filt1D_Pt_Phi_vs_Eta"]   ->Fill( Filt1D_It->WeightedEta(), Filt1D_It->WeightedPhi(), Filt1D_It->Pt() );
	  hist2D["Filt1D_Pt_iPhi_vs_iEta"] ->Fill( Filt1D_It->iEta(), Filt1D_It->iPhi(), Filt1D_It->Pt() );

	}


	// ******************************************************
	// *               Pre-PU subtracted jets               *
	// ******************************************************
	SUBPRINT("PrePUS")
	for (l1slhc::L1TowerJetCollection::const_iterator PrePUS_It = PrePUSubUpgrCenJet_Tower->begin(); PrePUS_It != PrePUSubUpgrCenJet_Tower->end(); 
	     ++PrePUS_It ){
	  
	  double weightedEta = PrePUS_It->WeightedEta();
	  double weightedPhi = PrePUS_It->WeightedPhi();
	  
	  hist2D["PrePUS_Pt_Phi_vs_Eta"]     ->Fill( weightedEta, weightedPhi, PrePUS_It->Pt() );
	  hist2D["PrePUS_Pt_iPhi_vs_iEta"]   ->Fill( PrePUS_It->iEta(), PrePUS_It->iPhi(), PrePUS_It->Pt() );
	  hist2D["PrePUS_AsymEta_Phi_vs_Eta"]->Fill( weightedEta, weightedPhi, PrePUS_It->AsymEta() );
	  hist2D["PrePUS_AsymPhi_Phi_vs_Eta"]->Fill( weightedEta, weightedPhi, PrePUS_It->AsymPhi() );

	}
      }



      
      // ******************************************************
      // *                 PU subtracted jets                 *
      // ******************************************************
      SUBPRINT("PUS")
      for (l1slhc::L1TowerJetCollection::const_iterator PUS_It = UpgrCenJet_Tower->begin(); PUS_It != UpgrCenJet_Tower->end(); ++PUS_It ){
	
	hist2D["PUS_Pt_Phi_vs_Eta"]->Fill( PUS_It->WeightedEta(), PUS_It->WeightedPhi(), PUS_It->Pt() );
	hist2D["PUS_Pt_iPhi_vs_iEta"]->Fill( PUS_It->iEta(), PUS_It->iPhi(), PUS_It->Pt() );

      }

      // ******************************************************
      // *              Local PU subtracted jets              *
      // ******************************************************
      SUBPRINT("LPUS")
      for (l1slhc::L1TowerJetCollection::const_iterator LPUS_It = LocalPUS_Tower->begin(); LPUS_It != LocalPUS_Tower->end(); ++LPUS_It ){
	
	hist2D["LPUS_Pt_Phi_vs_Eta"]->Fill( LPUS_It->WeightedEta(), LPUS_It->WeightedPhi(), LPUS_It->Pt() );
	hist2D["LPUS_Pt_iPhi_vs_iEta"]->Fill( LPUS_It->iEta(), LPUS_It->iPhi(), LPUS_It->Pt() );

      }




//       // ******************************************************
//       // *                     Ak5 Jets                       *
//       // ******************************************************
//       SUBPRINT("Ak5")
//       for (unsigned int iOff = 0; iOff < offlineJets.size(); iOff++){

// 	hist2D["Ak5_Pt_Phi_vs_Eta"]->Fill( offlineJets[iOff].Eta(), offlineJets[iOff].Phi(), offlineJets[iOff].Pt() );
	
//       }

//       // ******************************************************
//       // *                   Current Jets                     *
//       // ******************************************************

// #ifdef CURRENT
// 	SUBPRINT("Current")
//       for (unsigned int iCur = 0; iCur < currentJets.size(); iCur++){

// 	hist2D["Curr_Pt_Phi_vs_Eta"]->Fill( currentJets[iCur].Eta(), currentJets[iCur].Phi(), currentJets[iCur].Pt() );
	
//       }
// #endif












// ********************************************************************************
// *                                                                              *
// *                                 BENCHMARKING                                 *
// *                                                                              *
// ******************************************************************************** 


	// Calibrate PrePUS jets with LUT
	l1slhc::L1TowerJetCollection PrePUSCalib_Tower;
	l1slhc::L1TowerJetCollection PrePUS;
	
	for (l1slhc::L1TowerJetCollection::const_iterator PrePUS_It = PrePUSubUpgrCenJet_Tower->begin(); PrePUS_It != PrePUSubUpgrCenJet_Tower->end();
             ++PrePUS_It ){
	  


	  L1TowerJet calibJet = (*PrePUS_It);
	  int iPhi            = calibJet.iPhi();
	  int iEta            = calibJet.iEta();

	  // HACK TO GET THE FUNCTION TO ACCEPT MY JETS :(
	  PrePUS.insert( iEta, iPhi, (*PrePUS_It) );

	  // Extract iEta dependent correction factors
	  // **************************************************
	  int iEtaIndex       = iEta + 28;
	  if (iEta > 0)  // Correct for missing iEta = 0
	    iEtaIndex--;
	  
	  double unCorrectedPt = calibJet.Pt(); 
	  double pTOffset      = iEtaPtOffset[iEtaIndex];
	  double pTCorrection  = iEtaPtCorrection[iEtaIndex];
	  
	  double correctedPt   = unCorrectedPt*pTCorrection + pTOffset;

	  calibJet.setPt( correctedPt );

// 	  std::cout << "iEta = " << iEta << "\tpTUncorr = " << unCorrectedPt << "\tpTCorr = " << pTCorrection 
// 		    << "\tpTOffset = " << pTOffset << "\tCalibrated pT = " << correctedPt << "\tpTdiff = " 
// 		    << correctedPt - unCorrectedPt << "\tpTFracDiff = " << (correctedPt - unCorrectedPt)/unCorrectedPt 
// 		    << "\n";

	  // Store the calibrated jet
	  PrePUSCalib_Tower.insert( iEta, iPhi, calibJet );


	  // Calibration plots
	  // **************************************************
	  double deltaPt = (correctedPt - unCorrectedPt)/ unCorrectedPt;
	  hist1D["PrePUS_DeltaPT"]            ->Fill(deltaPt);
	  hist2D["PrePUS_DeltaPT_vs_UncorrPt"]->Fill(unCorrectedPt, deltaPt);
	  hist2D["PrePUS_DeltaPT_vs_iEta"]    ->Fill(iEta, deltaPt);
	  hist2D["PrePUS_DeltaPT_vs_NVTX"]    ->Fill(NVTX, deltaPt);

	}


	// --------------------------------------------------------------------------------
	// Make TLorentzVector of Calibrated PrePUS jets
	
	//	l1slhc::L1TowerJetCollection PrePUS = (*PrePUSubUpgrCenJet_Tower);

	// 
	vector<TLorentzVector> PrePUSCalib;

	for (l1slhc::L1TowerJetCollection::const_iterator Tower_It = PrePUSCalib_Tower.begin(); Tower_It != PrePUSCalib_Tower.end(); ++Tower_It ){
	  //    for (l1extra::L1JetParticleCollection::const_iterator il1 = CalibJets_L1extra->begin(); il1!= CalibJets_L1extra->end() ; ++il1 ){

	  up_l1.SetPtEtaPhiM(Tower_It->p4().Pt(),Tower_It->p4().Eta(),Tower_It->p4().Phi(),Tower_It->p4().M());

	  if( up_l1.Pt() < MIN_L1_PT ) continue;
	  if( fabs(up_l1.Eta()) > 3.0 ) continue;

	  PrePUSCalib.push_back(up_l1);

	} // End PrePUSCalib jets


	if(PrePUSCalib.size() != 0)
	  sort(PrePUSCalib.begin(), PrePUSCalib.end(),sortTLorentz);

	// Compare calibration
	if (PrePUSCalib.size() != 0){
	  

	  if ( offlineUncorrJets.size() != 0 ){

	    fillL1OfflineHistograms( PrePUSCalib, offlineUncorrJets, "PrePUSCalib", NVTX );
	  }
	  
	  //	  fillL1Histograms( PrePUSCalib, "PrePUSCalib", NVTX );
	}


	// --------------------------------------------------------------------------------

// 	// PrePUS-ak5PrePUS
// 	//		matchOnlineOfflineJets( PrePUSubUpgrCenJet_Tower, PrePUSAk5Jets, "PrePUS_ak5PrePUS", NVTX );
// 	matchOnlineOfflineJets( PrePUS, PrePUSAk5Jets, "PrePUS_ak5PrePUS", NVTX );
// // 	// PrePUS-ak5PUS
// // 	matchOnlineOfflineJets( PrePUSubUpgrCenJet_Tower, PUSAk5Jets, "PrePUS_ak5PUS", NVTX );
// // 	// PUS-ak5PUS
// // 	matchOnlineOfflineJets( UpgrCenJet_Tower, PUSAk5Jets, "PUS_ak5PUS", NVTX );


	// TO BE ADDED WITH CALIBRATED JETS:
	// PrePUSCalib-ak5PrePUS
	//	matchOnlineOfflineJets( PrePUSCalib_Tower, PrePUSAk5Jets, "PrePUSCalib_ak5PrePUS", NVTX );
	// PrePUSCalib+PUS-ak5PUS
	// PrePUSCalib-ak5PUS
	// PUSCalib-ak5PUS


















//    iEvent.put( outputNVTX, "NVTX" );
//    iEvent.put( outputprePusHT_Tower , "prePusHT_Tower" );
//    iEvent.put( outputHT_Tower,  "HT_Tower" );
//    iEvent.put( outputMHT_Tower, "MHT_Tower" );
//    iEvent.put( outputPUS_Tower_TLorentz, "PUS_Tower" );



    PRINT("End event")
  } // end valid event

  QUIT("Successfully completed analysis of single event")

} // end analyser




// ------------ method called once each job just after ending the event loop  ------------
void AnalyzeJets::endJob(){

}







//////////////  member functions  ///////////////

// Matched jet: Distance between jets is dist_ref_up. First and ref of matched L1 jet is dist_ref_up.second

std::pair< double, int> AnalyzeJets::Leading_Match( vector<TLorentzVector> offlineJets, vector<TLorentzVector> L1Jets ){

  std::pair< double, int> dist_ref;
  dist_ref.first  = 999; //this is the distance between jets
  dist_ref.second = 999; //dummy initialization of jet reference

  for(unsigned int i = 0; i < L1Jets.size(); ++i){

    //    if ( (offlineJets[0].Pt() <= 0) || (L1Jets[i].Pt() <= 0) ) continue;

    // Check whether deltaR of jets is smaller than a previous match
    if( L1Jets[i].DeltaR(offlineJets[0]) < dist_ref.first ) {
      dist_ref.first  = L1Jets[i].DeltaR(offlineJets[0]);
      dist_ref.second = i;
    }


  }

  return dist_ref;
}


// ------------ method called when starting to processes a run  ------------
void 
AnalyzeJets::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
AnalyzeJets::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
AnalyzeJets::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
AnalyzeJets::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnalyzeJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// // Median calculating code used in L1TowerJetPUEstimator
// double Median( vector<double> aVec){

//   // Order vector collection
//   sort( aVec.begin(), aVec.end() );

//   double median(0);
//   int size = aVec.size();
//   int halfSize = size/2;
//   if( size == 0 ){
//     median = 0;
//   }
//   else if( size == 1 ){
//     median = aVec[0];
//   }
//   else if( size%2 == 0 ){
//     // Even number of entries, take average of the values around center
//     median = ( aVec[ halfSize - 1 ] + aVec[ halfSize ] ) * 0.5;
//   }
//   else{
//     // Odd number of entries, halfSize is central element
//     median = aVec[ halfSize ];
//   }

//   return median;
// }








// ------------------------------------------------------------------------------------------------------------------------------------------------------
// -                                                                Histogram functions                                                                 -
// ------------------------------------------------------------------------------------------------------------------------------------------------------


// ****************************************************************************************************
// *                                       TowerJet Histograms                                        *
// ****************************************************************************************************

void AnalyzeJets::fillTowerJetHistograms(edm::Handle<l1slhc::L1TowerJetCollection> const& TowerJet, TString jetPrefix){


  // <TEMPORARY>
  unsigned int jetIndex(0);
  vector<double> jetPtAreaRatio2;
  unsigned int skipJetsIndex(0); // Skip the lead jet
  // </TEMPORARY>

  int nJets(0);
  

  // Container for local rho determination
  std::vector< std::vector<double> > PUrhoSlice;
  std::vector<double> jetRhoMin, jetRhoMax;
  // Container for local jet multiplicities
  std::vector< Int_t > PUnjetSlice;
  // Create a container for each eta slice to store the jet rho in order to calculate local rhos, don't need last bin
  PUrhoSlice.resize(  PUetaSlice.size() - 1 );
  jetRhoMin.resize(   PUetaSlice.size() - 1, 999 ); // Fill with dummy data
  jetRhoMax.resize(   PUetaSlice.size() - 1 );
  PUnjetSlice.resize( PUetaSlice.size() - 1 );
 






		

  for (l1slhc::L1TowerJetCollection::const_iterator Tower_It = TowerJet->begin(); Tower_It != TowerJet->end(); ++Tower_It ){
	

    // ****************************************
    // *    Load the TowerJet data            *
    // ****************************************
    
    int iEta           = Tower_It->iEta();
    int iPhi           = Tower_It->iPhi();
    int E              = Tower_It->E();
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

    // ratio of the jet pt to its area
    double ptAreaRatio = Pt/JetRealArea;

    
    // <TEMPORARY>

    // The following code is temporary in order to validate the calculation of rho. Remove ASAP.
    // No rho calibration is performed and only the leading jet is skipped.
    // THIS SHOULD NOT BE PERFORMED ON PU SUBTRACTED JETS

    // Skip the specified number of jets in the calculation of the median for rho
    // numSkipJets = 1    =>   Skip leading jet only
    if( jetIndex > skipJetsIndex ){
      jetPtAreaRatio2.push_back( ptAreaRatio );
    }
    jetIndex++;

    // </TEMPORARY>

    SUBPRINT("\tiDistributions")
    // iEta, iPhi distributions
    hist1D[jetPrefix + "_iPhi"]->Fill(iPhi);
    hist1D[jetPrefix + "_iEta"]->Fill(iEta);

    hist1D[jetPrefix + "_E"]      ->Fill(E);
    hist1D[jetPrefix + "_AsymEta"]->Fill(AsymEta);
    hist1D[jetPrefix + "_AsymPhi"]->Fill(AsymPhi);

    SUBPRINT("\tArea distributions")
    // Fill the jet area distributions
//     hist2D[jetPrefix + "_JetRealArea_Phi"] ->Fill(WeightedPhi, JetRealArea );
//     hist2D[jetPrefix + "_JetRealArea_Eta"] ->Fill(WeightedEta, JetRealArea );
//     hist2D[jetPrefix + "_JetRealArea_iPhi"]->Fill(iPhi,        JetRealArea );
//     hist2D[jetPrefix + "_JetRealArea_iEta"]->Fill(iEta,        JetRealArea );
//     hist2D[jetPrefix + "_JetRealArea_Pt"]  ->Fill(JetRealArea, Pt );


    hist2D[jetPrefix + "_Pt-AreaRatio_Phi"]    ->Fill(Phi, ptAreaRatio );
    hist2D[jetPrefix + "_Pt-AreaRatio_Eta"]    ->Fill(Eta, ptAreaRatio );
    hist2D[jetPrefix + "_Pt-AreaRatio_iPhi"]->Fill(iPhi, ptAreaRatio );
    hist2D[jetPrefix + "_Pt-AreaRatio_iEta"]->Fill(iEta, ptAreaRatio );

	
    // jet quantities
    SUBPRINT("\tJet quantities")
    hist1D[jetPrefix + "_UnweightedEta"]      ->Fill(Eta);
    hist1D[jetPrefix + "_UnweightedPhi"]      ->Fill(Phi);
    hist1D[jetPrefix + "_AllJet_L1Eta"]       ->Fill(WeightedEta);
    hist1D[jetPrefix + "_AllJet_L1Phi"]       ->Fill(WeightedPhi);
    hist1D[jetPrefix + "_AllJet_L1PT"]        ->Fill(Pt);
    hist2D[jetPrefix + "_AllJet_PhivsEta"]    ->Fill( WeightedEta, WeightedPhi );
    hist2D[jetPrefix + "_AllJet_iPhivsiEta"]  ->Fill( iEta, iPhi );
    hist2D[jetPrefix + "_AllJet_PtvsEta"]     ->Fill( WeightedEta, Pt );
    hist2D[jetPrefix + "_AllJet_PtvsPhi"]     ->Fill( WeightedPhi, Pt );
    hist1D[jetPrefix + "_AllJet_PtvsEta_prof"]->Fill( WeightedEta, Pt );
    hist1D[jetPrefix + "_AllJet_PtvsPhi_prof"]->Fill( WeightedPhi, Pt );
						
    hist1D[jetPrefix + "_AllJet_L1Eta_PtWeight"]->Fill(WeightedEta, Pt);
    hist1D[jetPrefix + "_AllJet_L1Phi_PtWeight"]->Fill(WeightedPhi, Pt);


    // ************************************************************
    // *                    PT Threshold plots                    *
    // ************************************************************
    SUBPRINT("\tpT threshold distributions")
    for (unsigned int iSlice = 0; iSlice < ptThreshold.size(); iSlice++){


      // get pT threshold
      double ptThresh   = ptThreshold[iSlice];

      // fill the respective pT bin
      if ( Pt > ptThresh ){

	// get pT threshold
	TString ptThreshStr = Form("%d",Int_t(ptThresh));

	TString ptLabel = "pT>" + ptThreshStr;
	TString ptTitle = "    (" + ptLabel + ")";

	hist1D[jetPrefix + "_Eta_" + ptLabel]->Fill(WeightedEta);
	hist1D[jetPrefix + "_Phi_" + ptLabel]->Fill(WeightedPhi);
	hist1D[jetPrefix + "_PT_"  + ptLabel]->Fill(Pt);

      }
    }




    // ************************************************************
    // *                   Eta binned histograms                  *
    // ************************************************************
	
    // Bin TowerJets in eta
    SUBPRINT("\tEta binned distributions")	
    for (unsigned int iSlice = 0;iSlice < PUetaSlice.size() - 1; iSlice++){
	  
      // get eta slice range
      double PUetaLow  = PUetaSlice[iSlice];
      double PUetaHigh = PUetaSlice[iSlice + 1];

	  
      // fill the respective eta bin
      if ( (WeightedEta >= PUetaLow) && (WeightedEta < PUetaHigh) ){
	    
	// get eta slice range
	double PUetaLow  = PUetaSlice[iSlice];
	double PUetaHigh = PUetaSlice[iSlice + 1];
	
	TString PUetaLowStr  = Form("%1.3f", PUetaLow);
	TString PUetaHighStr = Form("%1.3f", PUetaHigh);
	TString PUetaLabel = "eta(" + PUetaLowStr + "-" + PUetaHighStr + ")";
	    

	// Fill the jet in respective vector for calculating rho (lower edge)
	PUrhoSlice[iSlice].push_back(ptAreaRatio);
	PUnjetSlice[iSlice]++;
	
	if (ptAreaRatio < jetRhoMin[iSlice])
	  jetRhoMin[iSlice] = ptAreaRatio;

	if (ptAreaRatio > jetRhoMax[iSlice])
	  jetRhoMax[iSlice] = ptAreaRatio;
	     
	// Fill jet distributions
	hist1D[jetPrefix + "_Eta_"    + PUetaLabel]->Fill(WeightedEta);  
	hist1D[jetPrefix + "_Phi_"    + PUetaLabel]->Fill(WeightedPhi);  
	hist1D[jetPrefix + "_PT_"     + PUetaLabel]->Fill(Pt);   
	hist1D[jetPrefix + "_JetRho_" + PUetaLabel]->Fill(ptAreaRatio);
	    
	break;
      }
	  
    } // end eta binning
	


    // Increment the jet count
    nJets++;

	

	
  } // End TowerJet collection


    SUBPRINT("\tRho")
    // <TEMPORARY>
    // Calculate rho, the median of the jet Pt-Area ratio
  double globalRho = Median(jetPtAreaRatio2);
  hist1D[jetPrefix + "_Rho"]->Fill(globalRho);
  // </TEMPORARY>


  // Calculate rho for bins of eta
  for (unsigned int iSlice = 0; iSlice < PUrhoSlice.size(); iSlice++){
      
    double etaBinnedRho = Median(PUrhoSlice[iSlice]);
    double deltaRhoMin, deltaRhoMax;

    Int_t etaBinnednJets = PUnjetSlice[iSlice];

    // get eta slice range
    double PUetaLow  = PUetaSlice[iSlice];
    double PUetaHigh = PUetaSlice[iSlice + 1];
    TString PUetaLowStr  = Form("%1.3f", PUetaLow);
    TString PUetaHighStr = Form("%1.3f", PUetaHigh);
    TString PUetaLabel = "eta(" + PUetaLowStr + "-" + PUetaHighStr + ")"; 


    // Fill the local rho
    hist1D[jetPrefix + "_Rho_"         + PUetaLabel]->Fill(etaBinnedRho);  
    
    // Fill only if a minimum or maximum exist
    if ( etaBinnedRho != 0 ){

      double rhoMin = jetRhoMin[iSlice];
      double rhoMax = jetRhoMax[iSlice];
      deltaRhoMin = rhoMin - etaBinnedRho; 
      deltaRhoMax = rhoMax - etaBinnedRho; 

      hist1D[jetPrefix + "_RhoMin_" + PUetaLabel]->Fill(rhoMin);
      hist1D[jetPrefix + "_RhoMax_" + PUetaLabel]->Fill(rhoMax);

      hist1D[jetPrefix + "_deltaRhoMin_" + PUetaLabel]->Fill(deltaRhoMin);
      hist1D[jetPrefix + "_deltaRhoMax_" + PUetaLabel]->Fill(deltaRhoMax);  
    }


    hist1D[jetPrefix + "_NJETS_"  + PUetaLabel]->Fill(etaBinnednJets);

  }	


  // Incredibly inefficient...
  if (PUrhoSlice.size() == 4){

    // Delta rho, barrel minus endcap
    double negRhoEC, negRhoB, posRhoB, posRhoEC;
    negRhoEC = Median(PUrhoSlice[0]);
    negRhoB  = Median(PUrhoSlice[1]);
    posRhoB  = Median(PUrhoSlice[2]);      
    posRhoEC = Median(PUrhoSlice[3]);


    hist2D[jetPrefix + "_RhoECvsRhoB-"]     ->Fill( negRhoB,  negRhoEC );
    hist2D[jetPrefix + "_RhoECvsRhoB+"]     ->Fill( posRhoB,  posRhoEC );
    hist2D[jetPrefix + "_RhoGlobalvsRhoB-"] ->Fill( negRhoB,  globalRho );
    hist2D[jetPrefix + "_RhoGlobalvsRhoB+"] ->Fill( posRhoB,  globalRho );
    hist2D[jetPrefix + "_RhoGlobalvsRhoEC-"]->Fill( negRhoEC, globalRho );
    hist2D[jetPrefix + "_RhoGlobalvsRhoEC+"]->Fill( posRhoEC, globalRho );

    hist1D[jetPrefix + "_deltaRhoBEC-"]->Fill(negRhoB - negRhoEC);
    hist1D[jetPrefix + "_deltaRhoBEC+"]->Fill(posRhoB - posRhoEC);


  }
  else{
    //    std::cout << "Disable this code, it assumes 4 eta bins in total (2 barrel, 2 endcap)\n";
  }

  // Multiplicity
  hist1D[jetPrefix + "_NJETS"]->Fill(nJets);


  SUBPRINT("\tEnd")
}


// Produces correlation plots between two towerjet collections using iEta and iPhi matching. Jets found in JetColl1 but not in JetColl2 are stored as 
// dead jets. Hence JetColl1 should be a superset of JetColl2.

void AnalyzeJets::fillL1CorrelationHistograms( edm::Handle<l1slhc::L1TowerJetCollection> const& JetColl1, 
					       edm::Handle<l1slhc::L1TowerJetCollection> const& JetColl2, 
					       TString prefix, Int_t NVTX ){



    // Iterate through the first jet collection first as this is a superset of the second
    for (l1slhc::L1TowerJetCollection::const_iterator JetColl1_It = JetColl1->begin(); 
	 JetColl1_It != JetColl1->end(); ++JetColl1_It ){

      bool matchedJet = false;

      // Use iEta and iPhi as identifiers. The overlap filtering in the algorithm ensures these should be unique
      // and these quantities will not change after PU subtraction.
      int iEta1 = JetColl1_It->iEta();
      int iPhi1 = JetColl1_It->iPhi();

      // Try and match to a jet in the second collection
      for (l1slhc::L1TowerJetCollection::const_iterator JetColl2_It = JetColl2->begin(); 
	   JetColl2_It != JetColl2->end(); ++JetColl2_It ){

	// Get the second collection identifiers
	int iEta2 = JetColl2_It->iEta();
	int iPhi2 = JetColl2_It->iPhi();


	if ( (iEta1 == iEta2) && (iPhi1 == iPhi2) ){ // Jets match

	  // ************************************************************
	  // *          Produce the correlation distributions           *
	  // ************************************************************

	  // Jet distributions
	  hist2D["Corr_" + prefix + "_AllJet_L1Eta"]->Fill(JetColl1_It->WeightedEta(), JetColl2_It->WeightedEta());
	  hist2D["Corr_" + prefix + "_AllJet_L1Phi"]->Fill(JetColl1_It->WeightedPhi(), JetColl2_It->WeightedPhi());
	  hist2D["Corr_" + prefix + "_AllJet_L1PT"] ->Fill(JetColl1_It->Pt(),          JetColl2_It->Pt());

	  matchedJet = true;
	  break;
	}
	
      } // End JetColl2

      if ( matchedJet == false ){ // Dead jet
	
	// These jets were killed in the second collection, RIP
	hist1D["Dead_" + prefix + "_AllJet_L1Eta"]       ->Fill(JetColl1_It->WeightedEta());
	hist1D["Dead_" + prefix + "_AllJet_L1Phi"]       ->Fill(JetColl1_It->WeightedPhi());
	hist1D["Dead_" + prefix + "_AllJet_L1PT"]        ->Fill(JetColl1_It->Pt());
	hist1D["Dead_" + prefix + "_AllJet_L1PhivsL1Eta"]->Fill(JetColl1_It->WeightedEta(), JetColl1_It->WeightedPhi());
	hist1D["Dead_" + prefix + "_AllJet_L1PTvsL1Eta"] ->Fill(JetColl1_It->WeightedEta(), JetColl1_It->Pt());
      
      }

    } // End JetColl1


    //    hist2D["Corr_NJETS"]->Fill(nJetsPrePUS, nJetsPostPUS);

}



// ********************************************************************************
// *          Fill lead jet distributions for online-offline matched jets         *
// ******************************************************************************** 

void AnalyzeJets::fillL1OfflineHistograms( vector<TLorentzVector> L1Jets, vector<TLorentzVector> OfflineJets, TString prefix, Int_t NVTX ){


  // deltaR distance between jets and reference of matched L1 jet
  std::pair< double, int> dist_ref_up = Leading_Match(OfflineJets, L1Jets);


  if (dist_ref_up.first < deltaR){


    // Index of the L1 jet matched to the offline jet
    Int_t L1Match = dist_ref_up.second;

    // Extract offline and online quantities for the leading jet
    double OffPt  = OfflineJets[0].Pt();
    double OffEta = OfflineJets[0].Eta();
    double OffPhi = OfflineJets[0].Phi();
    double l1Pt   = L1Jets[L1Match].Pt();
    double l1Eta  = L1Jets[L1Match].Eta();
    double l1Phi  = L1Jets[L1Match].Phi();
    double deltaPt  = (l1Pt - OffPt)/OffPt;
    double deltaEta = l1Eta - OffEta;
    double deltaPhi = getAzimuth(l1Phi - OffPhi);


    // Resolution distributions
    hist1D[prefix + "_DeltaEta"]->Fill(deltaEta);
    hist1D[prefix + "_DeltaPhi"]->Fill(deltaPhi);
    hist1D[prefix + "_DeltaPT"] ->Fill(deltaPt); 

    hist1D[prefix + "_Eta"]     ->Fill(l1Eta);
    hist1D[prefix + "_Phi"]     ->Fill(l1Phi);
    hist1D[prefix + "_PT"]      ->Fill(l1Pt);


    // correlation distributions
    hist2D[prefix + "_DeltaPT_vs_L1PT"]          ->Fill(l1Pt,   deltaPt);
    hist2D[prefix + "_DeltaPhi_vs_DeltaEta"]     ->Fill(deltaPhi, deltaEta);
    hist1D[prefix + "_DeltaPT_vs_L1PT_prof"]     ->Fill(l1Pt,   deltaPt);
    hist1D[prefix + "_DeltaPhi_vs_DeltaEta_prof"]->Fill(deltaPhi, deltaEta);

    // Online-offline quantity correlations
    hist2D[prefix + "_L1PT_vs_OffPT"]         ->Fill(OffPt,  l1Pt);
    hist2D[prefix + "_L1Eta_vs_OffEta"]       ->Fill(OffEta, l1Eta);
    hist2D[prefix + "_L1Phi_vs_OffPhi"]       ->Fill(OffPhi, l1Phi);
    hist1D[prefix + "_L1PT_vs_OffPT_prof"]    ->Fill(OffPt,  l1Pt);
    hist1D[prefix + "_L1Eta_vs_OffEta_prof"]  ->Fill(OffEta, l1Eta);
    hist1D[prefix + "_L1Phi_vs_OffPhi_prof"]  ->Fill(OffPhi, l1Phi);
    // pT Correlations
    hist2D[prefix + "_DeltaPT_vs_OffPT"]      ->Fill(OffPt,  deltaPt);
    hist2D[prefix + "_DeltaEta_vs_OffPT"]     ->Fill(OffPt,  deltaEta);
    hist2D[prefix + "_DeltaPhi_vs_OffPT"]     ->Fill(OffPt,  deltaPhi);
    hist1D[prefix + "_DeltaPT_vs_OffPT_prof"] ->Fill(OffPt,  deltaPt);
    hist1D[prefix + "_DeltaEta_vs_OffPT_prof"]->Fill(OffPt,  deltaEta);
    hist1D[prefix + "_DeltaPhi_vs_OffPT_prof"]->Fill(OffPt,  deltaPhi);
    // NVTX Correlations
    hist2D[prefix + "_DeltaPT_vs_NVTX"]       ->Fill(NVTX,  deltaPt);
    hist2D[prefix + "_DeltaEta_vs_NVTX"]      ->Fill(NVTX,  deltaEta);
    hist2D[prefix + "_DeltaPhi_vs_NVTX"]      ->Fill(NVTX,  deltaPhi);
    hist1D[prefix + "_DeltaPT_vs_NVTX_prof"]  ->Fill(NVTX,  deltaPt);
    hist1D[prefix + "_DeltaEta_vs_NVTX_prof"] ->Fill(NVTX,  deltaEta);
    hist1D[prefix + "_DeltaPhi_vs_NVTX_prof"] ->Fill(NVTX,  deltaPhi);
    // Phi Correlations
    hist2D[prefix + "_DeltaPT_vs_Phi"]      ->Fill(OffPhi,  deltaPt);
    hist2D[prefix + "_DeltaEta_vs_Phi"]     ->Fill(OffPhi,  deltaEta);
    hist2D[prefix + "_DeltaPhi_vs_Phi"]     ->Fill(OffPhi,  deltaPhi);
    hist1D[prefix + "_DeltaPT_vs_Phi_prof"] ->Fill(OffPhi,  deltaPt);
    hist1D[prefix + "_DeltaEta_vs_Phi_prof"]->Fill(OffPhi,  deltaEta);
    hist1D[prefix + "_DeltaPhi_vs_Phi_prof"]->Fill(OffPhi,  deltaPhi);
    // Eta Correlations
    hist2D[prefix + "_DeltaPT_vs_Eta"]      ->Fill(OffEta,  deltaPt);
    hist2D[prefix + "_DeltaEta_vs_Eta"]     ->Fill(OffEta,  deltaEta);
    hist2D[prefix + "_DeltaPhi_vs_Eta"]     ->Fill(OffEta,  deltaPhi);
    hist1D[prefix + "_DeltaPT_vs_Eta_prof"] ->Fill(OffEta,  deltaPt);
    hist1D[prefix + "_DeltaEta_vs_Eta_prof"]->Fill(OffEta,  deltaEta);
    hist1D[prefix + "_DeltaPhi_vs_Eta_prof"]->Fill(OffEta,  deltaPhi);

//     // profile histograms
//     hist1D[prefix + "_Delta_PT_vs_OffPT_prof"]->Fill(OffPt,  deltaPt);
//     hist1D[prefix + "_Delta_PT_vs_PU_prof"]   ->Fill(NVTX, deltaPt);


    // ************************************************************
    // *                    pT binned histograms                  *
    // ************************************************************

    for (unsigned int iSlice = 0;iSlice < ptSlice.size() - 1; iSlice++){

      // get pT slice range
      double ptLow  = ptSlice[iSlice];
      double ptHigh = ptSlice[iSlice + 1];

      // fill the respective pT bin
      if ( (OffPt >= ptLow) && (OffPt < ptHigh) ){

	TString ptLowStr  = Form("%d",Int_t(ptLow));
	TString ptHighStr = Form("%d",Int_t(ptHigh));
	TString ptLabel = "pT(" + ptLowStr + "-" + ptHighStr + ")";

	hist1D[prefix + "_DeltaEta_" + ptLabel]->Fill(deltaEta);
	hist1D[prefix + "_DeltaPhi_" + ptLabel]->Fill(deltaPhi);
	hist1D[prefix + "_DeltaPT_"  + ptLabel]->Fill(deltaPt);
	hist1D[prefix + "_offEta_"   + ptLabel]->Fill(OffEta); 
	hist1D[prefix + "_offPhi_"   + ptLabel]->Fill(OffPhi);
	hist1D[prefix + "_PU_"       + ptLabel]->Fill(NVTX);

	break;
      }

    } // end pT binning


    // ************************************************************
    // *                   eta binned histograms                  *
    // ************************************************************

    for (unsigned int iSlice = 0;iSlice < etaSlice.size() - 1; iSlice++){

      // get eta slice range
      double etaLow  = etaSlice[iSlice];
      double etaHigh = etaSlice[iSlice + 1];

      // fill the respective eta bin
      if ( (OffEta >= etaLow) && (OffEta < etaHigh) ){

	TString etaLowStr  = Form("%1.3f", etaLow);
	TString etaHighStr = Form("%1.3f", etaHigh);
	TString etaLabel = "eta(" + etaLowStr + "-" + etaHighStr + ")";

	hist2D[prefix + "_Off_vs_On_PT_" + etaLabel]     ->Fill(l1Pt, OffPt);
	hist1D[prefix + "_Off_vs_On_PT_prof_" + etaLabel]->Fill(l1Pt, OffPt);

	hist1D[prefix + "_DeltaEta_" + etaLabel]->Fill(deltaEta);
	hist1D[prefix + "_DeltaPhi_" + etaLabel]->Fill(deltaPhi);
	hist1D[prefix + "_DeltaPT_"  + etaLabel]->Fill(deltaPt);
	hist1D[prefix + "_offEta_"   + etaLabel]->Fill(OffEta);
	hist1D[prefix + "_offPhi_"   + etaLabel]->Fill(OffPhi);
	hist1D[prefix + "_PU_"       + etaLabel]->Fill(NVTX);

	break;
      }

    } // end eta binning

    // Jet multiplicities
    //hist1D[prefix + "_cen_NJETS"]->Fill(cenJets);	  
    //hist1D[prefix + "_for_NJETS"]->Fill(forJets);	  

  } // Inside deltaR cone


} // End L1-Offline jet comparisons





// // ********************************************************************************
// // *                  Fill all jet distributions for online jets                  *
// // ******************************************************************************** 

// void AnalyzeJets::fillL1Histograms( vector<TLorentzVector> L1Jets, TString prefix, Int_t NVTX ){

//   // Event quantities

//   hist1D[prefix + "_NJETS"]->Fill(L1Jets.size());





//   TString etaPre;
//   // Number of forward and central jets
//   Int_t cenJets(0), forJets(0);
//   unsigned int skippedJets(0);


//   for (unsigned int iJet = 0; iJet < L1Jets.size(); iJet++){


//     double L1Pt   = L1Jets[iJet].Pt();
//     double L1Eta  = L1Jets[iJet].Eta();
//     double L1Phi  = L1Jets[iJet].Phi();


//     if ( fabs(L1Eta) > 3.0 ){
//       skippedJets++; // account for any jets that are outside the barrel + endcap acceptance
//       continue;
//     }
	  
//     // determine whether the jet is in barrel
//     if ( fabs(L1Eta) < CEN_BOUNDARY ){
//       etaPre = "_cen";
//       cenJets++;
//     }
//     else{
//       etaPre = "_for";
//       forJets++;
//     }

//     // Jet triggers: First, second, third and fourth jet
//     if(iJet < 4){
//       TString jetNum = Form("%d",iJet + 1);
//       hist1D[prefix + "_Jet_" + jetNum + "_PT"]->Fill( L1Pt );
//     }    

//     // ***********************************************************************
//     // *                    Fill central and forward jets                    *
//     // ***********************************************************************

//     hist2D[prefix + etaPre + "_L1Phi_vs_L1Eta"]->Fill( L1Eta, L1Phi ); 
//     hist1D[prefix + etaPre + "_AllJet_L1Eta"]  ->Fill( L1Eta );
//     hist1D[prefix + etaPre + "_AllJet_L1Phi"]  ->Fill( L1Phi );
//     hist1D[prefix + etaPre + "_AllJet_L1PT"]   ->Fill( L1Pt );

//     if (iJet < 13 + skippedJets) // first 13 jets
//       hist1D[prefix + etaPre + "_13Jet_L1PT"]  ->Fill( L1Pt );	 

//     // Excluding lead jet
//     if (iJet != 0){
//       hist1D[prefix + etaPre + "_AllJetExcLead_L1PT"] ->Fill( L1Pt );
//       if (iJet < 13 + skippedJets) // first 13 jets
// 	hist1D[prefix + etaPre + "_13JetExcLead_L1PT"]->Fill( L1Pt );
//     }

//     // ***********************************************************************
//     // *                    Fill entire range                                *
//     // ***********************************************************************

//     hist2D[prefix + "_L1Phi_vs_L1Eta"]->Fill( L1Eta, L1Phi ); 
//     hist1D[prefix + "_AllJet_L1Eta"]  ->Fill( L1Eta );
//     hist1D[prefix + "_AllJet_L1Phi"]  ->Fill( L1Phi );
//     hist1D[prefix + "_AllJet_L1PT"]   ->Fill( L1Pt );

//     hist1D[prefix + "_AllJet_L1PT"]   ->Fill( L1Pt );
//     if (iJet < 13) // first 13 jets
//       hist1D[prefix + "_13Jet_L1PT"]  ->Fill( L1Pt );	  

//     // Excluding lead jet
//     if (iJet != 0){
//       hist1D[prefix + "_AllJetExcLead_L1PT"] ->Fill( L1Pt );
//       if (iJet < 13) // first 13 jets
// 	hist1D[prefix + "_13JetExcLead_L1PT"]->Fill( L1Pt );
//     }


//   }


// }






// ********************************************************************************
// *                                                                              *
// *                                 BENCHMARKING                                 *
// *                                                                              *
// ******************************************************************************** 



//void AnalyzeJets::matchOnlineOfflineJets( edm::Handle<l1slhc::L1TowerJetCollection> const& L1Jets, edm::Handle<reco::CaloJetCollection> const& ak5Jets,
//					  TString prefix, Int_t NVTX ){
// void AnalyzeJets::matchOnlineOfflineJets( edm::Handle<l1slhc::L1TowerJetCollection> const& L1Jets, edm::Handle<edm::View< reco::CaloJet > > const& ak5Jets, 
// 					  TString prefix, Int_t NVTX ){

// void AnalyzeJets::matchOnlineOfflineJets( edm::Handle<l1slhc::L1TowerJetCollection> const& L1Jets, edm::Handle<edm::View< reco::CaloJet > > const& ak5Jets, 
// 					  TString prefix, Int_t NVTX ){
// void AnalyzeJets::matchOnlineOfflineJets( l1slhc::L1TowerJetCollection const& L1Jets, edm::Handle<edm::View< reco::CaloJet > > const& ak5Jets, 
// 					  TString prefix, Int_t NVTX ){

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


//   //  int ak5TightJetCount = 0;
//   double maxDeltaR     = 0.5; // Maximum allowed deltaR between matched online and offline jets


//   //      ********************************************************************************
//   //      *                      PrePUS Online-Offline jet matching                      *
//   //      ********************************************************************************
//   //


//   std::vector<TLorentzVector> ak5Vec, L1Vec;
//   std::vector<int>            L1iEtaVec;

//   //           *********************************************
//   //           *          Fill ak5 TLorentzVector          *
//   //           *********************************************

//   for (edm::View< reco::CaloJet >::const_iterator ak5_It = ak5Jets->begin(); ak5_It != ak5Jets->end(); ++ak5_It ){

//     double ak5Pt         = ak5_It->p4().Pt();
//     double ak5Eta        = ak5_It->p4().Eta();
//     double ak5Phi        = ak5_It->p4().Phi();
//     double ak5M          = ak5_It->p4().M();

//     // ========================================
//     // Jet quality cuts
//     // ========================================

//     if( ak5Pt < MIN_OFF_PT )   continue;
//     if( fabs( ak5Eta ) > 3.0 ) continue;
    
//     // ========================================

//     TLorentzVector ak5Jet( ak5Pt, ak5Eta, ak5Phi, ak5M );
//     ak5Vec.push_back( ak5Jet );

//   }


//   //           *********************************************
//   //           *          Fill L1 TLorentzVector           *
//   //           *********************************************


//   for (l1slhc::L1TowerJetCollection::const_iterator L1_It = L1Jets.begin(); L1_It != L1Jets.end(); ++L1_It ){

//     double L1Pt         = L1_It->p4().Pt();
//     double L1Eta        = L1_It->p4().Eta();
//     double L1Phi        = L1_It->p4().Phi();
//     double L1M          = L1_It->p4().M();

//     TLorentzVector L1Jet( L1Pt, L1Eta, L1Phi, L1M );
//     L1Vec.push_back( L1Jet );

//     L1iEtaVec.push_back(L1_It->iEta());

//   }



//   // ------------------------------ Perform jet matching ------------------------------


//   SUBPRINT("Extracted jets")

//   std::vector<pair_info> pairs = make_pairs( ak5Vec, L1Vec );

//   std::sort(pairs.begin(), pairs.end(), sortDR);

//   std::vector<int> L1MatchedIndex = analyse_pairs(pairs, L1Vec.size(), maxDeltaR);

// //   // ____________________
// //   // LEAD JET MATCHING
// //   // ____________________

// //   int leadAk5Index = -1;
// //   double minDeltaR = 9999;

// //   for(unsigned int iAk5 = 0; iAk5 < ak5Vec.size(); iAk5++) {
// //     double curDeltaR = L1Vec[0].DeltaR(ak5Vec[iAk5]);
// //     if (curDeltaR < minDeltaR){
// //       minDeltaR = curDeltaR;
// //       leadAk5Index = iAk5;
// //     }
// //   }


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
//     if (jetsMatched > 1){
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


//     //    TString iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEta)));


//     double deltaEta = L1Vec[ L1Index ].Eta() - ak5Vec[ ak5Index ].Eta();
//     double deltaPhi = ak5Vec[ ak5Index ].DeltaPhi( L1Vec[ L1Index ] );
    
//     double ak5Pt    = ak5Vec[ ak5Index ].Pt();
//     double ak5Eta   = ak5Vec[ ak5Index ].Eta();
//     double ak5Phi   = ak5Vec[ ak5Index ].Phi();

//     double L1Pt     = L1Vec[ L1Index ].Pt();
//     double L1Eta    = L1Vec[ L1Index ].Eta();
//     double L1Phi    = L1Vec[ L1Index ].Phi();

//     double deltaPt  = (L1Pt - ak5Pt)/ak5Pt; 


//     // ************************************************************
//     // *                         Unbinned                         *
//     // ************************************************************

//    hist2D[matchPre + "_L1PT_vs_OffPT"]       ->Fill(ak5Pt, L1Pt);
//    hist1D[matchPre + "_OffPT_vs_L1PT_prof"]  ->Fill(L1Pt, ak5Pt);
//    hist1D[matchPre + "_DeltaPT_vs_NVTX_prof"]->Fill(NVTX, deltaPt);
//    hist1D[matchPre + "_DeltaPT_vs_Phi_prof"] ->Fill(L1Phi, deltaPt);
//    hist1D[matchPre + "_DeltaPT_vs_Eta_prof"] ->Fill(L1Eta, deltaPt);
    
  

//     // ************************************************************
//     // *                          Binned                          *
//     // ************************************************************


//     // Online-offline quantity correlations
    
//     hist2D[matchPre + "_L1Eta_vs_OffEta_" + iEtaStr]     ->Fill(ak5Eta, L1Eta);       
//     hist2D[matchPre + "_L1Phi_vs_OffPhi_" + iEtaStr]     ->Fill(ak5Phi, L1Phi);    
//     hist2D[matchPre + "_DeltaPhi_vs_DeltaEta_" + iEtaStr]->Fill(deltaEta, deltaPhi);    
//     hist2D[matchPre + "_DeltaPT_vs_L1PT_" + iEtaStr]     ->Fill(L1Pt, deltaPt);    

//     // pT Correlations
//     hist2D[matchPre + "_L1PT_vs_OffPT_" + iEtaStr]     ->Fill(ak5Pt, L1Pt);
//     hist2D[matchPre + "_DeltaPT_vs_OffPT_" + iEtaStr]  ->Fill(ak5Pt, deltaPt);
//     hist2D[matchPre + "_DeltaEta_vs_OffPT_" + iEtaStr] ->Fill(ak5Pt, deltaEta);
//     hist2D[matchPre + "_DeltaPhi_vs_OffPT_" + iEtaStr] ->Fill(ak5Pt, deltaPhi);

//     // NVTX Correlations
//     hist2D[matchPre + "_L1PT_vs_NVTX_" + iEtaStr]    ->Fill(NVTX, L1Pt);
//     hist2D[matchPre + "_DeltaPT_vs_NVTX_" + iEtaStr] ->Fill(NVTX, deltaPt);
//     hist2D[matchPre + "_DeltaEta_vs_NVTX_" + iEtaStr]->Fill(NVTX, deltaEta);
//     hist2D[matchPre + "_DeltaPhi_vs_NVTX_" + iEtaStr]->Fill(NVTX, deltaPhi);

//       // Phi Correlations
//     hist2D[matchPre + "_L1PT_vs_Phi_" + iEtaStr]     ->Fill(L1Phi, L1Pt);    
//     hist2D[matchPre + "_DeltaPT_vs_Phi_" + iEtaStr]  ->Fill(L1Phi, deltaPt);         
//     hist2D[matchPre + "_DeltaEta_vs_Phi_" + iEtaStr] ->Fill(L1Phi, deltaEta);        
//     hist2D[matchPre + "_DeltaPhi_vs_Phi_" + iEtaStr] ->Fill(L1Phi, deltaPhi);        
//       // Eta Correlations
//     hist2D[matchPre + "_L1PT_vs_Eta_" + iEtaStr]     ->Fill(L1Eta, L1Pt);          
//     hist2D[matchPre + "_DeltaPT_vs_Eta_" + iEtaStr]  ->Fill(L1Eta, deltaPt);         
//     hist2D[matchPre + "_DeltaEta_vs_Eta_" + iEtaStr] ->Fill(L1Eta, deltaEta);        
//     hist2D[matchPre + "_DeltaPhi_vs_Eta_" + iEtaStr] ->Fill(L1Eta, deltaPhi);        
    


//     // Calibration TProfiles
//     // --------------------------------------------------------------------------------
//     hist1D[matchPre + "_OffPT_vs_L1PT_prof_"+ iEtaStr]          ->Fill(L1Pt, ak5Pt);
//     hist1D[matchPre + "_DeltaPT_vs_NVTX_prof_" + iEtaStr]       ->Fill(NVTX, deltaPt);
//     hist1D[matchPre + "_DeltaPT_vs_Phi_prof_" + iEtaStr]        ->Fill(L1Phi,deltaPt);
// //     hist1D[matchPre + "_DeltaPT_vs_GlobalRho_prof_" + iEtaStr]  ->Fill(L1Pt, ak5Pt);
// //     hist1D[matchPre + "_DeltaPT_vs_EorHFired_prof_" + iEtaStr]  ->Fill(L1Pt, ak5Pt);
// //     hist1D[matchPre + "_DeltaPT_vs_EplusHtotal_prof_" + iEtaStr]->Fill(L1Pt, ak5Pt);
      



//   }
  
//   // --------------------------------------------------------------------------------





//   // Ak5 jets are ordered by pT => 1st in list is lead jet

// //   for (edm::View< reco::CaloJet >::const_iterator PrePUSak5_It = PrePUSAk5Jets->begin(); PrePUSak5_It != PrePUSAk5Jets->end(); ++PrePUSak5_It ){

    
// //     double ak5Pt         = PrePUSak5_It->p4().Pt();
// //     double ak5Eta        = PrePUSak5_It->p4().Eta();
// //     double ak5Phi        = PrePUSak5_It->p4().Phi();
// //     double ak5M          = PrePUSak5_It->p4().M();
// //     //    TLorentzVector ak5Vec = 

// //     double ak5JetArea    = PrePUSak5_It->jetArea();
// //     double ak5TowersArea = PrePUSak5_It->towersArea();




// //     // ========================================
// //     // Jet quality cuts
// //     // ========================================

// //     if( ak5Pt < MIN_OFF_PT )  continue;
// //     if( fabs( ak5Eta ) > 3.0 ) continue;
    
// //      hist1D["Ak5_JetArea"]         ->Fill(ak5JetArea);
// //      hist2D["Ak5_pT_vs_JetArea"]   ->Fill(ak5JetArea, ak5Pt);			       
// //      hist1D["Ak5_TowersArea"]      ->Fill(ak5TowersArea);
// //      hist2D["Ak5_pT_vs_TowersArea"]->Fill(ak5TowersArea, ak5Pt);			       
// //   }


// }
















// ------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------------------------------------




//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzeJets);

//  LocalWords:  iEtaSubDir
