//
// Producer to calibrate calojets given a specified calibration scheme
// Fix: Jets are sorted after calibration to give correct pT ordering
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include <algorithm> 

using namespace edm;
using namespace std;
using namespace reco;

// Ranking function for sort
bool caloJetRankDescending ( reco::CaloJet jet1, reco::CaloJet jet2 ){ return ( jet1.pt() > jet2.pt() ); }



class CalCaloProducer : public edm::EDProducer {
   public:
      explicit CalCaloProducer(const edm::ParameterSet&);
      ~CalCaloProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
     ParameterSet conf_;


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
CalCaloProducer::CalCaloProducer(const edm::ParameterSet& iConfig):
conf_(iConfig)
{
   //register your products
produces<CaloJetCollection>();
  
}


CalCaloProducer::~CalCaloProducer()
{
 

}


//
// member functions
//


// ------------ method called to produce the data  ------------
void
CalCaloProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  const JetCorrector* jetCorrector = JetCorrector::getJetCorrector(conf_.getParameter<string>("JetCorrector"),iSetup);
  edm::Handle<reco::CaloJetCollection> Calojets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CaloJets"), Calojets);

auto_ptr< CaloJetCollection > outputColl(new CaloJetCollection());


 reco::CaloJetCollection  unsortedCaloJets;

 // Perform JEC on uncalibrated jets. Note these may no longer be pT ordered! 
  for(reco::CaloJetCollection::const_iterator it = Calojets->begin(); it != Calojets->end() ; ++it) {
    
    //     double scale = jetCorrector->correction(*it,iEvent,iSetup);
    //     CaloJet* oldCorrectedJet = new CaloJet(scale*it->p4(),it->vertex(),it->getSpecific ());
    //      outputColl->push_back(CaloJet(scale*it->p4(),it->vertex(),it->getSpecific ()));

    // Recommended method:  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetAnalysis#CorrOnTheFly
    CaloJet correctedJet = *it; 
    double jec = jetCorrector->correction(*it,iEvent,iSetup); 
    correctedJet.scaleEnergy(jec);   

    unsortedCaloJets.push_back( correctedJet );

//     if ( (oldCorrectedJet->pt() != correctedJet.pt() ) ){
 //       std::cout << oldCorrectedJet->pt() << "\t" << correctedJet.pt() << "\t" << (oldCorrectedJet->pt() == correctedJet.pt() ) << "\n";
//     }

   }

  // Sort jets
  std::sort( unsortedCaloJets.begin(), unsortedCaloJets.end(), caloJetRankDescending );

  // Store the now sorted jet collection
  for(reco::CaloJetCollection::const_iterator it = unsortedCaloJets.begin(); it != unsortedCaloJets.end() ; ++it) {
    outputColl->push_back( *it );
  }


 iEvent.put(outputColl);

 
}

// ------------ method called once each job just before starting event loop  ------------
void 
CalCaloProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CalCaloProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
CalCaloProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CalCaloProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CalCaloProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CalCaloProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CalCaloProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CalCaloProducer);
