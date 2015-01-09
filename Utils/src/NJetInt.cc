// -*- C++ -*-
//
// Package:    NJetInt
// Class:      NJetInt
// 
/**\class NJetInt NJetInt.cc RA2Classic/NJetInt/src/NJetInt.cc
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

//
// class declaration
//

class NJetInt : public edm::EDProducer {
public:
	explicit NJetInt(const edm::ParameterSet&);
	~NJetInt();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    bool GoodJets(unsigned int i, edm::Handle< edm::View<pat::Jet> > Jets );

	edm::InputTag JetTag_;
	
	
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
NJetInt::NJetInt(const edm::ParameterSet& iConfig)
{
	//register your produc
	JetTag_ = iConfig.getParameter<edm::InputTag>("JetTag");
	
	produces<int>("NJets");
	produces<int>("NJets50");
	produces<int>("NJets30");
	produces<int>("NJets20");	
	/* Examples
	 *   produces<ExampleData2>();
	 * 
	 *   //if do put with a label
	 *   produces<ExampleData2>("label");
	 * 
	 *   //if you want to put into the Run
	 *   produces<ExampleData2,InRun>();
	 */
	//now do what ever other initialization is needed
	
}


NJetInt::~NJetInt()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
NJetInt::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	int NJets=0;
	int NJets50=0;
	int NJets30=0;
	int NJets20=0;
	edm::Handle< edm::View<pat::Jet> > Jets;
	iEvent.getByLabel(JetTag_,Jets);
	if( Jets.isValid() ) {
		
		NJets=Jets->size();
		for(int j=0; j<NJets; ++j){
            bool isGood=GoodJets(j, Jets);
            if(!isGood)continue;
		if(Jets->at(j).pt()<20 )continue;
		++NJets20;
		if(Jets->at(j).pt()<30)continue;
		++NJets30;
		if(Jets->at(j).pt()<50)continue;
		++NJets50;
		}

	}
	else std::cout<<"NJetInt::Invlide Tag: "<<JetTag_.label()<<std::endl;
	std::auto_ptr<int> htp(new int(NJets));
	iEvent.put(htp,"NJets");
	
        std::auto_ptr<int> htp1(new int(NJets20));
        iEvent.put(htp1,"NJets20");
	
        std::auto_ptr<int> htp2(new int(NJets30));
        iEvent.put(htp2,"NJets30");	

	std::auto_ptr<int> htp3(new int(NJets50));
        iEvent.put(htp3,"NJets50");
}

bool NJetInt::GoodJets(unsigned int i, edm::Handle< edm::View<pat::Jet> > Jets ){
    bool isGood=false;
    if( Jets.isValid() ) {
        //   float pt=Jets->at(i).pt();
        float eta=Jets->at(i).eta();
        float neufrac=Jets->at(i).neutralHadronEnergyFraction();//gives raw energy in the denominator
        float phofrac=Jets->at(i).neutralEmEnergyFraction();//gives raw energy in the denominator
        float chgfrac=Jets->at(i).chargedHadronEnergyFraction();
        float chgEMfrac=Jets->at(i).chargedEmEnergyFraction();
        
        // int nconstit=Jets->at(i).getPFConstituents().size();
        int chgmulti=Jets->at(i).chargedHadronMultiplicity();
        if( fabs(eta)<2.4 && neufrac<0.99 && phofrac<0.99 &&chgmulti>0 && chgfrac>0 && chgEMfrac<0.99)isGood=true;
        
    }
    return isGood;
}

// ------------ method called once each job just before starting event loop  ------------
void 
NJetInt::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NJetInt::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
NJetInt::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
NJetInt::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
NJetInt::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
NJetInt::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NJetInt::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NJetInt);
