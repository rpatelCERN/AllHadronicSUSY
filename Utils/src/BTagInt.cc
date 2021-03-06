// -*- C++ -*-
//
// Package:    BTagInt
// Class:      BTagInt
// 
/**\class BTagInt BTagInt.cc RA2Classic/BTagInt/src/BTagInt.cc
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
#include "DataFormats/Candidate/interface/Candidate.h"

//
// class declaration
//

class BTagInt : public edm::EDProducer {
public:
	explicit BTagInt(const edm::ParameterSet&);
	~BTagInt();
	
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
	std::string   btagname_;
	double        btagvalue_;
	
	
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
BTagInt::BTagInt(const edm::ParameterSet& iConfig)
{
	//register your produc
	JetTag_ = iConfig.getParameter<edm::InputTag>("JetTag");
	btagname_ = iConfig.getParameter<std::string>  ("BTagInputTag");
	btagvalue_   = iConfig.getParameter<double>       ("BTagCutValue");
	
	produces<int>("BTags");
	produces<int>("BTags20");
	produces<int>("BTags30");
	produces<int>("BTags50");


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


BTagInt::~BTagInt()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
BTagInt::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;
	using namespace pat;
	int BTags=0;
	int BTags20=0;
	      int BTags30=0;
	      int BTags50=0;
	std::auto_ptr< std::vector<double> > bDiscrim(new std::vector<double>);	
	edm::Handle< edm::View<pat::Jet> > Jets;
	iEvent.getByLabel(JetTag_,Jets);
	if( Jets.isValid() ) {
		for(unsigned int i=0; i<Jets->size();i++)
		{
            bool isGood=GoodJets(i, Jets);
            if(!isGood)continue;
		  if(Jets->at(i).bDiscriminator(btagname_) >btagvalue_)BTags++;
		  else continue;
		  if(Jets->at(i).pt()<20)continue;
		   ++BTags20;
		  if(Jets->at(i).pt()<30)continue;
		  ++BTags30;
		  if(Jets->at(i).pt()<50)continue;
                  ++BTags50;
		}
	}
	else std::cout<<"BTagInt::Invlide Tag: "<<JetTag_.label()<<std::endl;
	std::auto_ptr<int> htp(new int(BTags));
	iEvent.put(htp,"BTags");

	        std::auto_ptr<int> htp1(new int(BTags20));
        iEvent.put(htp1,"BTags20");
	
	                std::auto_ptr<int> htp2(new int(BTags30));
        iEvent.put(htp2,"BTags30");

	                std::auto_ptr<int> htp3(new int(BTags50));
        iEvent.put(htp3,"BTags50");
	
}

bool BTagInt::GoodJets(unsigned int i, edm::Handle< edm::View<pat::Jet> > Jets ){
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
BTagInt::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BTagInt::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
BTagInt::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
BTagInt::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
BTagInt::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BTagInt::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BTagInt::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BTagInt);
