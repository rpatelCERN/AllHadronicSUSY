// -*- C++ -*-
//
// Package:    METDouble
// Class:      METDouble
// 
/**\class METDouble METDouble.cc RA2Classic/METDouble/src/METDouble.cc
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

class METDouble : public edm::EDProducer {
public:
	explicit METDouble(const edm::ParameterSet&);
	~METDouble();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	edm::InputTag metTag_;
	
	
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
METDouble::METDouble(const edm::ParameterSet& iConfig)
{
	//register your produc
	metTag_ = iConfig.getParameter<edm::InputTag> ("METTag");
	
	produces<double>("MET");
        
	produces<double>("metEnergy");
	
	produces<double>("metPhi");
        produces<double>("mEtSig");	
	produces<double>("metSumEt");	
	
	produces<double>("genMET");
        produces<double>("genmetEt");
        produces<double>("genmetPhi");
	produces<double>("genmetE");	
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


METDouble::~METDouble()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
METDouble::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	double met_=0;
	double mEnergy_=0;
	double metPhi_=0;double mEtSig_=0;
	double metSumEt=0;
	double genMet=0; 
	double genmetEt=0;
	double genmetE=0;double genmetPhi=0;
	edm::Handle< edm::View<pat::MET> > MET;
	iEvent.getByLabel(metTag_,MET); 
	if(MET.isValid() ){
	met_=MET->at(0).pt();

	mEnergy_=MET->at(0).energy();
	metPhi_=MET->at(0).phi();
	mEtSig_=MET->at(0).mEtSig();
	metSumEt=MET->at(0).sumEt();
	genMet=MET->at(0).genMET()->pt();

	genmetEt=MET->at(0).genMET()->et();
	genmetPhi=MET->at(0).genMET()->phi();
	genmetE=MET->at(0).genMET()->energy();

	//std::cout<<MET->at(0).mEtSig()<<std::endl;

	}
	else std::cout<<"METDouble::Invlide Tag: "<<metTag_.label()<<std::endl;
	std::auto_ptr<double> htp(new double(met_));
	iEvent.put(htp,"MET");

        std::auto_ptr<double> htp1(new double(mEnergy_));
        iEvent.put(htp1, "metEnergy");

        std::auto_ptr<double> htp2(new double(metPhi_));
        iEvent.put(htp2, "metPhi");	

	std::auto_ptr<double> htp3(new double(mEtSig_));
        iEvent.put(htp3, "mEtSig");

	std::auto_ptr<double> htp4(new double(metSumEt));
	iEvent.put(htp4, "metSumEt");

        std::auto_ptr<double> htp5(new double(genMet));
        iEvent.put(htp5, "genMET");

	std::auto_ptr<double> htp6(new double(genmetEt));
        iEvent.put(htp6, "genmetEt");

                std::auto_ptr<double> htp7(new double(genmetE));
        iEvent.put(htp7, "genmetE");

                std::auto_ptr<double> htp8(new double(genmetPhi));
        iEvent.put(htp8, "genmetPhi");

}

// ------------ method called once each job just before starting event loop  ------------
void 
METDouble::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METDouble::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
METDouble::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
METDouble::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
METDouble::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
METDouble::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
METDouble::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(METDouble);
