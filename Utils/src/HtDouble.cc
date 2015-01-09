// -*- C++ -*-
//
// Package:    HTDouble
// Class:      HTDouble
// 
/**\class HTDouble HTDouble.cc RA2Classic/HTDouble/src/HTDouble.cc
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

//
// class declaration
//

class HTDouble : public edm::EDProducer {
public:
	explicit HTDouble(const edm::ParameterSet&);
	~HTDouble();
	
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
HTDouble::HTDouble(const edm::ParameterSet& iConfig)
{
	//register your produc
	JetTag_ = iConfig.getParameter<edm::InputTag>("JetTag");
	
	produces<double>("HT");
	produces<double>("genHT");
	produces<int>("genJetMatch");
	produces<double>("HT30");
	produces<double>("HT50");
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


HTDouble::~HTDouble()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HTDouble::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	double ht_=0; double genht_=0;
	double ht30_=0; double ht50_=0;
	int genMatches=0;
//	edm::Handle< edm::View<reco::Candidate> > Jets;
//	iEvent.getByLabel(JetTag_,Jets);
        edm::Handle< edm::View<pat::Jet> > Jets;
        iEvent.getByLabel(JetTag_,Jets);
	if( Jets.isValid() ) {
		for(unsigned int i=0; i<Jets->size();i++)
		{
            bool isGood=GoodJets(i, Jets);
            if(!isGood)continue;
			if(Jets->at(i).genJet()!=0)
			{
			++genMatches;
			genht_+=Jets->at(i).genJet()->pt();
			}
			ht_ +=Jets->at(i).pt();

			if(Jets->at(i).pt()<30)continue;

			ht30_+=Jets->at(i).pt();
			if(Jets->at(i).pt()<50)continue;
			ht50_+=Jets->at(i).pt();
			
		}
	}
	else std::cout<<"HTDouble::Invlide Tag: "<<JetTag_.label()<<std::endl;
	std::auto_ptr<double> htp(new double(ht_));
	iEvent.put(htp, "HT");

	 std::auto_ptr<double> htp1(new double(genht_));
        iEvent.put(htp1, "genHT");	

	
	         std::auto_ptr<double> htp2(new double(ht30_));
        iEvent.put(htp2, "HT30");

	                 std::auto_ptr<double> htp3(new double(ht50_));
        iEvent.put(htp3, "HT50");
	
	         std::auto_ptr<int> htp4(new int(genMatches));
        iEvent.put(htp4, "genJetMatch");

}

bool HTDouble::GoodJets(unsigned int i, edm::Handle< edm::View<pat::Jet> > Jets ){
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
HTDouble::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HTDouble::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
HTDouble::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HTDouble::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HTDouble::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HTDouble::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HTDouble::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HTDouble);
