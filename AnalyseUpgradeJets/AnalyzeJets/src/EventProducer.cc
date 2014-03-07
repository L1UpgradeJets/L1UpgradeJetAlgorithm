// -*- C++ -*-
//
// Package:    EventProducer
// Class:      EventProducer
// 
/**\class EventProducer EventProducer.cc AnalyseUpgradeJets/src/EventProducer.cc

 Description: Produces event variables for the JetAnalyzer validation framework.

 List of quantities produced: 

     - 
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

// REMOVE THIS..................
//const double PI = 3.141592654;




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
#include "AnalyseUpgradeJets/AnalyzeJets/interface/histBooker.h"
// #include "AnalyseUpgradeJets/AnalyzeJets/interface/JetMatch.h"
#include "AnalyseUpgradeJets/AnalyzeJets/src/helperFunctions.cc"

#include <algorithm>  // for sorting





//
// class declaration
//

class EventProducer : public edm::EDProducer {
   public:
      explicit EventProducer(const edm::ParameterSet&);
      ~EventProducer();

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
    
      // Fills the respective aboveThreshold bin 
      virtual void FillThreshold( int energy, std::vector<int> ttThreshold, std::vector<int> &aboveThreshold );
      // Makes the threshold bins cumulative
      virtual void CumulateThreshold( std::vector<int> &aboveThreshold );


      // ----------member data ---------------------------
      edm::ParameterSet conf_;
      edm::Service<TFileService> fs;

      // TT iEta range
      int ttiEtaRange;
      // TT thresholds
      std::vector <int> ttEThreshold;
      std::vector <int> ttHThreshold;
      std::vector <int> ttEplusHThreshold;
      // iEta ring width
      int ttRingWidth;
      int lBins;
      // Histogram containers
      std::map< TString, TH1*> hist1D;
      std::map< TString, TH1*> hist2D;

  


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
EventProducer::EventProducer(const edm::ParameterSet& iConfig): conf_(iConfig)
{

   produces <int>("TTEFired");
   produces <int>("TTHFired");
   produces <int>("TTEorHFired");
   produces <int>("TTEandHFired");

   // Thresholds
   produces < std::vector<int> >("TTEFiredAboveThreshold");
   produces < std::vector<int> >("TTHFiredAboveThreshold");
//    produces < std::vector<int> >("TTEorHFiredAboveThreshold");
//    produces < std::vector<int> >("TTEandHFiredAboveThreshold");
   produces < std::vector<int> >("TTEplusHFiredAboveThreshold");

   // ************************************************************
   // *   `                       TT energy                      *
   // ************************************************************

   produces <int>("TTETotal");
   produces <int>("TTHTotal");
//   produces <int>("TTEoverH");
   produces <int>("TTEMax");
   produces <int>("TTHMax");

  // ********************************************************************************
  // *                          TT information local                                *
  // ********************************************************************************

  produces < std::vector<int> >("TTLETotal");
  produces < std::vector<int> >("TTLHTotal");
  produces < std::vector<int> >("TTLEMax");
  produces < std::vector<int> >("TTLHMax");

  produces < std::vector<int> >("TTLEMedian");
  produces < std::vector<int> >("TTLHMedian");
  produces < std::vector<int> >("TTLEplusHMedian");




  // ****************************************************************************************************
  // *                                 Load configuration parameters                                    *
  // ****************************************************************************************************
  PRINT("TT thresholds")

  // TT |iEta| range to analyse
  ttiEtaRange       = iConfig.getParameter< int > ("TTiEtaRange");

  // TT thresholds
  ttEThreshold      = iConfig.getParameter< std::vector<int> > ("TTEThreshold");
  sort (ttEThreshold.begin(), ttEThreshold.end());   // ensure the bins are in ascending order 
  ttHThreshold      = iConfig.getParameter< std::vector<int> > ("TTHThreshold");
  sort (ttHThreshold.begin(), ttHThreshold.end());   // ensure the bins are in ascending order 
  ttEplusHThreshold = iConfig.getParameter< std::vector<int> > ("TTEplusHThreshold");
  sort (ttHThreshold.begin(), ttHThreshold.end());   // ensure the bins are in ascending order 

  // iEta ring width
  ttRingWidth       = iConfig.getParameter<int>("TTRingWidth");

  if ( 56 % ttRingWidth != 0){
    edm::LogWarning("TT ring width error") << "Error, iEta ring width of " << ttRingWidth 
					   << " does yield integer ring widths. " 
					   << "Specify a value of 'TTRingWidth; that is a factor of 56"
					   << " (1, 2, 4, 7, 8, 14, 28, 56).\n";
  }

  // Number of local bins
  lBins = 56/ttRingWidth;

  

  // Associate the histogram container with the histogram booker, a class for handling histogram booking
  histBooker booker( &hist1D, &hist2D );


  // create a TT subdirectory
  TFileDirectory TTSubDir = fs->mkdir( "TT" );


  booker.book1D( "TTEFired", TTSubDir, "Ecal cells fired", pEFired );
//      booker.book1D("TTEFired",
//      booker.book1D("TTEFired",
//      booker.book1D("TTEFired",
//      booker.book1D("TTEFired",


  

}


EventProducer::~EventProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
EventProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   bool evValid = true;

   
  // ********************************************************************************
  // *                          TT information global                               *
  // ********************************************************************************
   
   // ************************************************************
   // *   `                   TT multiplicity                    *
   // ************************************************************

   std::auto_ptr<int> outputTTEFired(     new int() );
   std::auto_ptr<int> outputTTHFired(     new int() );
   std::auto_ptr<int> outputTTEorHFired(  new int() );
   std::auto_ptr<int> outputTTEandHFired( new int() );

   // Thresholds
   std::auto_ptr< std::vector<int> > outputTTEFiredAboveThreshold(      new std::vector<int>() );
   std::auto_ptr< std::vector<int> > outputTTHFiredAboveThreshold(      new std::vector<int>() );
//    std::auto_ptr< std::vector<int> > outputTTEorHFiredAboveThreshold(   new std::vector<int>() );
//    std::auto_ptr< std::vector<int> > outputTTEandHFiredAboveThreshold(  new std::vector<int>() );
   std::auto_ptr< std::vector<int> > outputTTEplusHFiredAboveThreshold( new std::vector<int>() );


   // ************************************************************
   // *   `                       TT energy                      *
   // ************************************************************

   std::auto_ptr<int> outputTTETotal(      new int() );
   std::auto_ptr<int> outputTTHTotal(      new int() );
//   std::auto_ptr<int> outputTTEoverH(      new int() );
   std::auto_ptr<int> outputTTEMax(        new int() );
   std::auto_ptr<int> outputTTHMax(        new int() );

  // ********************************************************************************
  // *                          TT information local                                *
  // ********************************************************************************

  std::auto_ptr< std::vector<int> > outputTTLETotal(       new std::vector<int>() );
  std::auto_ptr< std::vector<int> > outputTTLHTotal(       new std::vector<int>() );
  std::auto_ptr< std::vector<int> > outputTTLEMax(         new std::vector<int>() );
  std::auto_ptr< std::vector<int> > outputTTLHMax(         new std::vector<int>() );

  std::auto_ptr< std::vector<int> > outputTTLEMedian(      new std::vector<int>() );
  std::auto_ptr< std::vector<int> > outputTTLHMedian(      new std::vector<int>() );
  std::auto_ptr< std::vector<int> > outputTTLEplusHMedian( new std::vector<int>() );


  // ****************************************************************************************************
  // *                                             Handles                                              *
  // ****************************************************************************************************

   // TT collection                                                                                                                                          
   SUBPRINT("Trigger towers")                                                                                                                                
   edm::Handle<l1slhc::L1CaloTowerCollection> TriggerTowers;
   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalorimeterTowers"), TriggerTowers);
   if(!TriggerTowers.isValid()){
     edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalorimeterTowers") << std::endl;
     evValid = false;
   }   



   if( evValid ){


     // ****************************************************************************************************
     // *                                         Trigger Towers                                           *
     // ****************************************************************************************************
     
     int totalE(0), totalH(0);
     int maxE(0), maxH(0);
     int EFiredTotal(0), HFiredTotal(0), EorHFiredTotal(0), EandHFiredTotal(0);



     // TTs fired above the threshold specified in the CMSSW config
     std::vector<int> EAboveThreshold( ttEThreshold.size() );
     std::vector<int> HAboveThreshold( ttEThreshold.size() );
//      std::vector<int> EorHAboveThreshold( ttEThreshold.size() );
//      std::vector<int> EandHAboveThreshold( ttEThreshold.size() );
     std::vector<int> EplusHAboveThreshold( ttEThreshold.size() );



     // ****************************************
     // *             Local binning            *
     // ****************************************
     std::vector<int> LETotal( lBins );
     std::vector<int> LHTotal( lBins );
     std::vector<int> LEMax( lBins );
     std::vector<int> LHMax( lBins );

     std::vector<int> LEMedian( lBins );
     std::vector<int> LHMedian( lBins );
     std::vector<int> LEplusHMedian( lBins );

     std::vector< std::vector<int> > LEMedianVec( lBins );
     std::vector< std::vector<int> > LHMedianVec( lBins );
     std::vector< std::vector<int> > LEplusHMedianVec( lBins );


     for( l1slhc::L1CaloTowerCollection::const_iterator TT_It = TriggerTowers->begin(); TT_It != TriggerTowers->end(); ++TT_It ){

       int iEta   = TT_It->iEta();
       // Restrict to iEta specified TTs
       if (abs(iEta) > ttiEtaRange)
	 continue;

       //       int iPhi   = TT_It->iPhi();
       //        double Eta = mTowerGeo.eta(iEta);
       //        double Phi = mTowerGeo.phi(iPhi);
       //       int EcalFG = lTT_It->EcalFG();
       //       int HcalFG = lTT_It->HcalFG();

       // Extract the iEta bin index for local binning
       int iEtaBinAbs = (abs(iEta) - 1)/ttRingWidth + 1;
       int iEtaBin;
       if (iEta > 0){ 
	 iEtaBin = lBins/2 + iEtaBinAbs - 1; 
       }
       else{ 
	 iEtaBin = lBins/2 - iEtaBinAbs; 
       }
       //      std::cout << iEtaBin << "\t" << iEta << "\n";


       // ****************************************
       // *    Load the calorimeter tower data   *
       // ****************************************
       int E      = TT_It->E();
       int H      = TT_It->H();
       totalE           += E;
       totalH           += H;
       LETotal[iEtaBin] += E;
       LHTotal[iEtaBin] += H;

       // Median - ZERO SUPPRESSED
       if ( E > 0 ){
	 LEMedianVec[iEtaBin].push_back( E );
       }
       if (H > 0){
	 LHMedianVec[iEtaBin].push_back( H );
       }
       int EplusH = E + H;
       if ( EplusH > 0 ){
	 LEplusHMedianVec[iEtaBin].push_back( EplusH );
       }
           
       // Store the event maximum TT E and H
       if (E > maxE)
	 maxE = E;
       if (H > maxH)
	 maxH = H;
       // local
       if (E > LEMax[iEtaBin])
         LEMax[iEtaBin] = E;
       if (H > LHMax[iEtaBin])
         LHMax[iEtaBin] = H;


//       double EoverH(0);
//       if (H != 0){
// 	EoverH = double(E)/double(H);
//       }

//       // Restrict phi to range [-pi,pi]
//       if (Phi > PI)
// 	Phi -= 2*PI;



//   // ************************************************************
//   // *                 iEta binned distributions                *
//   // ************************************************************

// //   SUBPRINT("iEta binned")
// //     TFileDirectory iEtaDir = matchSubDir.mkdir( "iEtaBinned" );

// //   TString previousiEtaBin = "";

// //   // SHOULD CUT THIS TO iETA <= 28 - (jetSize - 1)
// //      for (int iEta = -28;iEta <= 28;iEta++){

// //        if (iEta == 0)
// //  	continue;

//       // Get current iEta bin
//       // ****************************************
//       TString iEtaStr;
      
//       // iEta binned distributions    
//       if (ttRingWidth == 1){ // Individual iEta bins
// 	iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEta)));
//       }
//       else{ // iEta bins of width greater than 1 
	
// 	int iEtaBinAbs = (abs(iEta) - 1)/ttRingWidth;
// 	int sign = 1;
// 	if (iEta < 0)
// 	  sign = -1;
	
// 	double iEtaBinLowAbs  = ( (iEtaBinAbs)*ttRingWidth + 1 );
// 	double iEtaBinHighAbs = ( (iEtaBinAbs + 1)*ttRingWidth );

// // 	// Fix bin limits
// // 	if (iEtaBinLowAbs > 28)
// // 	  iEtaBinLowAbs = 28;
// // 	if (iEtaBinHighAbs > 28)
// // 	  iEtaBinHighAbs = 28;

// 	double iEtaBinLow  = sign*iEtaBinLowAbs;
// 	double iEtaBinHigh = sign*iEtaBinHighAbs;
	
// 	if (iEtaBinLow > iEtaBinHigh)
// 	  std::swap(iEtaBinLow, iEtaBinHigh);
	
// 	iEtaStr = "iEta_" + TString(Form("%d",Int_t(iEtaBinLow))) + "to" + TString(Form("%d",Int_t(iEtaBinHigh)));


// 	//std::cout << "\n\niEta = " << iEta << " is in the bin range:   " << iEtaBinLow << "to" << iEtaBinHigh << "\n";
//       }

//       //     }







      
      bool EFired(false), HFired(false); //, EcalOrHcalFired(0), EcalAndHcalFired(0);
      if (E > 0){
	EFired = true;;
	EFiredTotal++;

      }
      if (H > 0){
	HFired = true;
	HFiredTotal++;
      }
      if ( (EFired || HFired) ){
	EorHFiredTotal++;

	if ( (EFired && HFired) ){
	  EandHFiredTotal++;
	}
      }


      // ****************************************
      // *              Thresholds              *
      // ****************************************

      // Fill TT threshold threshold vectors
      FillThreshold( E,      ttEThreshold,      EAboveThreshold );
      FillThreshold( H,      ttHThreshold,      HAboveThreshold );
      FillThreshold( EplusH, ttEplusHThreshold, EplusHAboveThreshold );



      // ****************************************
      // *               Local TT               *
      // ****************************************

      // Calculate the median of quantities for each iEtaBin
      for ( int iEtaBin = 0; iEtaBin < lBins; ++iEtaBin ){
	
	LEMedian[iEtaBin]      = Median( LEMedianVec[iEtaBin] );
	LHMedian[iEtaBin]      = Median( LHMedianVec[iEtaBin] );
	LEplusHMedian[iEtaBin] = Median( LEplusHMedianVec[iEtaBin] );

				   
      }

//   // ********************************************************************************
//   // *                          TT information local                                *
//   // ********************************************************************************

//   std::auto_ptr< std::vector<int> > outputLTTETotal(       new std::vector<int>() );
//   std::auto_ptr< std::vector<int> > outputLTTHTotal(       new std::vector<int>() );

//   std::auto_ptr< std::vector<int> > outputLTTEMax(         new std::vector<int>() );
//   std::auto_ptr< std::vector<int> > outputLTTHMax(         new std::vector<int>() );

//   std::auto_ptr< std::vector<int> > outputLTTEMedian(      new std::vector<int>() );
//   std::auto_ptr< std::vector<int> > outputLTTHMedian(      new std::vector<int>() );
//   std::auto_ptr< std::vector<int> > outputLTTEplusHMedian( new std::vector<int>() );


//   // iEta ring width
//   ttRingWidth       = iConfig.getParameter<int>("TTRingWidth");




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

     } // End TT loop


     // ****************************************
     // *   Make threshold vectors cumulative  *
     // ****************************************

     CumulateThreshold( EAboveThreshold );
     CumulateThreshold( HAboveThreshold );
     CumulateThreshold( EplusHAboveThreshold );















     // Fill histograms
     hist1D["TTEFired"]->Fill( EFiredTotal );



     // Store event quantities
     *outputTTEFired     = EFiredTotal;
     *outputTTHFired     = HFiredTotal;
     *outputTTEorHFired  = EorHFiredTotal;
     *outputTTEandHFired = EandHFiredTotal;

     *outputTTETotal     = totalE;
     *outputTTHTotal     = totalH;
     *outputTTEMax       = maxE;
     *outputTTHMax       = maxH;




     // Thresholds
     *outputTTEFiredAboveThreshold      = EAboveThreshold;
     *outputTTHFiredAboveThreshold      = HAboveThreshold;
//      *outputTTEorHFiredAboveThreshold   =
//      *outputTTEandHFiredAboveThreshold  =
     *outputTTEplusHFiredAboveThreshold = EplusHAboveThreshold;

     // Local
     *outputTTLETotal      = LETotal;
     *outputTTLHTotal      = LHTotal;
     *outputTTLEMax        = LEMax;
     *outputTTLHMax        = LHMax;
     *outputTTLEMedian      = LEMedian;
     *outputTTLHMedian      = LHMedian;
     *outputTTLEplusHMedian = LEplusHMedian;

	
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




     // TODO: ADD
     //            - LUMI?
     //            - RUN ?
     //            


 
     // ************************************************************
     // *   `                Event hisograms                       *
     // ************************************************************
     //     booker.book1D("TTEFired",
//      booker.book1D("TTEFired",
//      booker.book1D("TTEFired",
//      booker.book1D("TTEFired",
//      booker.book1D("TTEFired",




     iEvent.put( outputTTEFired,     "TTEFired" );     
     iEvent.put( outputTTHFired,     "TTHFired" );     
     iEvent.put( outputTTEorHFired,  "TTEorHFired" ) ;
     iEvent.put( outputTTEandHFired, "TTEandHFired" );

     // Thresholds
     iEvent.put( outputTTEFiredAboveThreshold,      "TTEFiredAboveThreshold");
     iEvent.put( outputTTHFiredAboveThreshold,      "TTHFiredAboveThreshold");
//      iEvent.put( outputTTEorHFiredAboveThreshold,   "TTEorHFiredAboveThreshold");
//      iEvent.put( outputTTEandHFiredAboveThreshold,  "TTEandHFiredAboveThreshold");
     iEvent.put( outputTTEplusHFiredAboveThreshold, "TTEplusHFiredAboveThreshold");

     
     // ************************************************************
     // *   `                       TT energy                      *
     // ************************************************************

     iEvent.put( outputTTETotal, "TTETotal");
     iEvent.put( outputTTHTotal, "TTHTotal");    
     //   iEvent.put( outputTTEoverH, "TTEoverHTotal");
     iEvent.put( outputTTEMax, "TTEMax");
     iEvent.put( outputTTHMax, "TTHMax");

     // ********************************************************************************
     // *                          TT information local                                *
     // ********************************************************************************
     
     iEvent.put( outputTTLETotal,      "TTLETotal");      
     iEvent.put( outputTTLHTotal,      "TTLHTotal");      
     iEvent.put( outputTTLEMax,        "TTLEMax");
     iEvent.put( outputTTLHMax,        "TTLHMax");
     
     iEvent.put( outputTTLEMedian,      "TTLEMedian");
     iEvent.put( outputTTLHMedian,      "TTLHMedian");
     iEvent.put( outputTTLEplusHMedian, "TTLEplusHMedian");

       
   }


}

// ------------ method called once each job just before starting event loop  ------------
void 
EventProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
EventProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
EventProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
EventProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
EventProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EventProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void
EventProducer::FillThreshold( int energy, std::vector<int> ttThreshold, std::vector<int> &aboveThreshold ) {


      // Store the first bin which the TT energy exceeds
      bool exceedsThreshold = false;
      if ( energy > ttThreshold[ 0 ] ){
	
	for( uint iE = 0; iE < ttThreshold.size() - 1; ++iE){
	  
	  int thresholdLow  = ttThreshold[ iE ];
	  int thresholdHigh = ttThreshold[ iE + 1 ];
	  
	  if ( (energy > thresholdLow) && (energy <= thresholdHigh) ){
	    
	    //	    std::cout << "E = " << energy << ",\t thresholdLow = " << thresholdLow << "\n";
	    aboveThreshold[ iE ]++;
	    exceedsThreshold = true;
	    break;
	  }
	  
	}
	if ( exceedsThreshold == false ){ // Energy exceeds all the thresholds
	  aboveThreshold[ aboveThreshold.size() - 1 ]++;
	  //	  std::cout << "E = " << energy << ",\t thresholdLow = " << ttEThreshold[ aboveThreshold.size() - 1 ] << "\n";
	}
      }
}



void
EventProducer::CumulateThreshold( std::vector<int> &aboveThreshold ) {

     int cumulative(0);
     for( int i = aboveThreshold.size() - 1; i > -1; --i){
	 
       // Make each bin a cumulative threshold
       aboveThreshold[ i ] += cumulative;
       cumulative           = aboveThreshold[ i ];

     }

}






//define this as a plug-in
DEFINE_FWK_MODULE(EventProducer);
