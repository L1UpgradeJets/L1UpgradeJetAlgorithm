/* L1TestPatternCaloTowerProducer 



   ~ Preliminary ~
   Extracts calorimeter tower information from a test pattern file and outputs
   the corresponding TT to an output collection. Enables a fine level debugging
   of the upgrade jet algorithm.


   Mark Baber Imperial College, London 

*/
/*
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

// Includes for the Calo Scales
#include "CondFormats/DataRecord/interface/L1CaloEcalScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "FWCore/Framework/interface/EventSetup.h"
*/
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
/*
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
*/
#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "SimDataFormats/SLHC/interface/L1CaloTowerFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloTriggerSetup.h"
#include "SimDataFormats/SLHC/interface/L1CaloTriggerSetupRcd.h"

#include <map>
#include <deque>
#include <fstream>
#include <string>
#include <iostream>

class L1TestPatternCaloTowerProducer:public edm::EDProducer
{
  public:
	explicit L1TestPatternCaloTowerProducer( const edm::ParameterSet & );
	 ~L1TestPatternCaloTowerProducer(  );

  private:

	virtual void produce( edm::Event &, const edm::EventSetup & );

        // Collection to store test pattern calotowers 
        std::auto_ptr < l1slhc::L1CaloTowerCollection > mCaloTowers;

        // Test pattern reader and filename information 
        std::ifstream   testPatternReader;
        edm::FileInPath testPatternDirectory;



};




L1TestPatternCaloTowerProducer::L1TestPatternCaloTowerProducer( const edm::ParameterSet & aConfig ):
  mCaloTowers( NULL )
{

  std::cout << "\n\n----------------------------------------\nBegin: TestPatternCaloTowerProdcuer\n----------------------------------------\n\n";

  // Extract the test pattern filename from the config file
  testPatternDirectory = aConfig.getParameter<edm::FileInPath> ("TestPatternFile");
  // Open the test pattern file
  testPatternReader.open( testPatternDirectory.fullPath().c_str() );


  // Register Product
  //  produces < l1slhc::L1CaloTowerCollection > ( "TestPatternCaloTowers" );
  produces < l1slhc::L1CaloTowerCollection > ( );



}


L1TestPatternCaloTowerProducer::~L1TestPatternCaloTowerProducer(  )
{
}



void L1TestPatternCaloTowerProducer::produce( edm::Event & aEvent, const edm::EventSetup & aSetup )
{


	// Create a new l1slhc::L1CaloTowerCollection (auto_ptr should handle deletion of the last one correctly)
	mCaloTowers = std::auto_ptr < l1slhc::L1CaloTowerCollection > ( new l1slhc::L1CaloTowerCollection );


	if(!testPatternReader){
	  throw cms::Exception("Cannot open test pattern file")
	    << "ERROR: Test pattern file '" << testPatternDirectory.fullPath().c_str()
	    << "' could not be opened, check the input in the configuration file is valid.\n";
	}

	
	// Print header
	//	std::cout << "\niEta\tiPhi\tE\tH\tFG\n";

	int iEta(0), iPhi(0), E(0), H(0), FG(0);
	bool reading(false);

	// Read in the test pattern data
	//	while ( !(testPatternReader.eof()) ) { 
	
	std::string line;

	while( std::getline(testPatternReader, line)){
         
	  //	  std::cout << line << "\n";
	  

	  // The '=' Character is used for comments and to split each testpattern
	  //	  if (testPatternReader.peek() == '='){
	  if ( line[0] == '='){
	   
	    // Pattern is already being read, '=' character signifies the end of a testpattern
	    if (reading == true){
	      reading = false;
	      break;
	    }
	    // Character is a comment, keep reading to find test pattern data
	    else{

	      continue;
	      //	      testPatternReader.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	    }

	  }
	  else{
	    reading = true;
	    sscanf(line.c_str(), "%d%d%d%d%d", &iEta, &iPhi, &E, &H, &FG);
	  }


// 	  testPatternReader >> iEta >> iPhi 
// 			    >> E    >> H    >> FG;

// 	  std::cout << iEta << "\t" <<  iPhi << "\t" 
// 		    <<  E   << "\t" <<  H    << "\t" 
// 		    <<  FG  << "\n";


	  l1slhc::L1CaloTower lCaloTower( iEta, iPhi );

	  // Store the E, H and FG information in the TT
	  lCaloTower.setEcal( E, FG );
	  lCaloTower.setHcal( H, FG );

	  // Store in the collection
	  mCaloTowers->insert( iEta , iPhi , lCaloTower );

	}
	testPatternReader.close();

	std::cout << "\nTest patterns loaded successfully.\n";


	aEvent.put( mCaloTowers);//, "TestPatternCaloTowers" );

}



DEFINE_EDM_PLUGIN( edm::MakerPluginFactory, edm::WorkerMaker < L1TestPatternCaloTowerProducer >,
		   "L1TestPatternCaloTowerProducer" );
DEFINE_FWK_PSET_DESC_FILLER( L1TestPatternCaloTowerProducer );
