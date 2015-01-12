// -*- C++ -*-
//
// Package:    GenLeptonRecoCand
// Class:      GenLeptonRecoCand
// 
/**\class GenLeptonRecoCand GenLeptonRecoCand.cc RA2Classic/GenLeptonRecoCand/src/GenLeptonRecoCand.cc
 * 
 * Description: [one line class summary]
 * 
 * Implementation:
 *     [Notes on implementation]
 */
//
// Original Author:  Arne-Rasmus Draeger,68/111,4719,
//         Created:  Fri Apr 11 16:35:33 CEST 2014
// $Id$
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
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

//
// class declaration
//

class GenInfo : public edm::EDProducer {
public:
	explicit GenInfo(const edm::ParameterSet&);
	~GenInfo();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	edm::InputTag PrunedGenParticleTag_;
	//int pdgID_;
	
	//const reco::GenParticle* TopFound(const reco::GenParticle * particle);
	//const reco::GenParticle* TauFound(const reco::GenParticle * particle);
	
	
	// ----------member data ---------------------------
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
GenInfo::GenInfo(const edm::ParameterSet& iConfig)
{
	//register your produc
	PrunedGenParticleTag_ 				= 	iConfig.getParameter<edm::InputTag >("PrunedGenParticleTag");

	//now do what ever other initialization is needed
    produces<double>("topweightOfficial");
    produces<double>("topweight");

	
}


GenInfo::~GenInfo()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GenInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	Handle<edm::View<reco::GenParticle> > pruned;
	iEvent.getByLabel(PrunedGenParticleTag_,pruned);
    float topPt=-1;
    float topBarPt=-1;

	for(size_t i=0; i<pruned->size();i++)
	{
        if(abs((*pruned)[i].pdgId() ) != 6)continue;
           if((*pruned)[i].pdgId() == 6)topPt=(*pruned)[i].pt();
           if((*pruned)[i].pdgId() == -6)topBarPt=(*pruned)[i].pt();
          if (topPt>=0 && topBarPt >=0) break;
	}
           const double a = 0.156; //combined 8 TeV values
           const double b =  -0.00137 ;
           //overflow
           if (topPt >400) topPt=400;
           if (topBarPt >400) topBarPt=400;

           double SFt = exp(a + b*topPt);
           double SFtbar = exp(a + b*topBarPt);
           
           double topweightOfficial = sqrt( SFt * SFtbar );
           
           if (topPt<0 || topBarPt <0)topweightOfficial=1.0;
           
           const  double p0 = 1.18246e+00;
           const  double p1 = 4.63312e+02;
           const  double p2 = 2.10061e-06;
           
           double x=topPt;
           if ( x>p1 ) x = p1; //use flat factor above 463 GeV
           double topweight = p0 + p2 * x * ( x - 2 * p1 );
           std::auto_ptr<double> htp1(new double(topweightOfficial));
           iEvent.put(htp1, "topweightOfficial");
    
    std::auto_ptr<double> htp2(new double(topweight));
    iEvent.put(htp2, "topweight");
    
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenInfo::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenInfo::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
GenInfo::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenInfo::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenInfo::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenInfo::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenInfo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenInfo);
