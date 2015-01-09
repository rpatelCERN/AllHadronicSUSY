// -*- C++ -*-
//
// Package:    RA2Selection
// Class:      RA2Selection
// 
/**\class RA2Selection RA2Selection.cc RA2Classic/RA2Selection/src/RA2Selection.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arne-Rasmus Draeger,,,uni-hamburg
//         Created:  Mon Oct 22 11:46:48 CEST 2012
// $Id: RA2Selection.cc,v 1.3 2012/12/06 15:28:37 adraeger Exp $
//
//


// system include files
#include "AllHadronicSUSY/Utils/interface/RA2Selection.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RA2Selection::RA2Selection(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   HtJetsTag_ = iConfig.getParameter<edm::InputTag>("HTJets");
   MhtTag_ = iConfig.getParameter<edm::InputTag>("MhtTag");
   MhtJetTag_ = iConfig.getParameter<edm::InputTag>("MhtJetTag");
   HTMin_ = iConfig.getParameter <double> ("HTMin");
   MHTMin_ = iConfig.getParameter <double> ("MHTMin");
   deltaPhi1_ = iConfig.getParameter <double> ("deltaPhi1");
   deltaPhi2_ = iConfig.getParameter <double> ("deltaPhi2");
   deltaPhi3_ = iConfig.getParameter <double> ("deltaPhi3");
   nJets_ = iConfig.getParameter <unsigned int> ("nJets");
}


RA2Selection::~RA2Selection()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
RA2Selection::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  double deltaPhi1=0;
  double deltaPhi2=0;
  double deltaPhi3=99;
  double ht=0;
  edm::Handle< edm::View<reco::Candidate> > htJets;
  iEvent.getByLabel(HtJetsTag_,htJets);
  if(htJets->size() < nJets_) return false;
   for( edm::View <reco::Candidate>::const_iterator HTJetsCan = htJets->begin(); HTJetsCan!=htJets->end();HTJetsCan++)
  	{
   		ht+= HTJetsCan->pt();
 	}

  if(ht < HTMin_) return false;
  edm::Handle< edm::View<reco::Candidate> > mht;
  iEvent.getByLabel(MhtTag_,mht);
  if( mht->at(0).pt() < MHTMin_ ) return false;
  edm::Handle< edm::View<reco::Candidate> > mhtJets;
  iEvent.getByLabel(MhtJetTag_,mhtJets);


  if( mhtJets.isValid() && mht.isValid() ) 
  {
    	if( mhtJets->size() > 0 ) 
	{
      		deltaPhi1 = std::abs(reco::deltaPhi(mhtJets->at(0).phi(),mht->at(0).phi()));
		if(deltaPhi1 < deltaPhi1_) return false;
      		if( mhtJets->size() > 1 ) 
		{	

			deltaPhi2 = std::abs(reco::deltaPhi(mhtJets->at(1).phi(),mht->at(0).phi()));
			if(deltaPhi2 < deltaPhi2_) return false;
			if( mhtJets->size() > 2 ) 
			{

	  			deltaPhi3 = std::abs(reco::deltaPhi(mhtJets->at(2).phi(),mht->at(0).phi()));
				if (deltaPhi3 < deltaPhi3_) return false;
			}
      		}
    	}
  }


   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
RA2Selection::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RA2Selection::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
RA2Selection::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
RA2Selection::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
RA2Selection::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
RA2Selection::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RA2Selection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(RA2Selection);
