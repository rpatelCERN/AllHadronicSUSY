// -*- C++ -*-
//
// Package:    DeltaPhiDouble
// Class:      DeltaPhiDouble
// 
/**\class DeltaPhiDouble DeltaPhiDouble.cc RA2Classic/DeltaPhiDouble/src/DeltaPhiDouble.cc
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
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
//
// class declaration
//

class DeltaPhiDouble : public edm::EDProducer {
public:
	explicit DeltaPhiDouble(const edm::ParameterSet&);
	~DeltaPhiDouble();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	edm::InputTag MHTJetTag_, DeltaPhiJetTag_, metTag_;
	
	
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
DeltaPhiDouble::DeltaPhiDouble(const edm::ParameterSet& iConfig)
{
	//register your produce
	metTag_ = iConfig.getParameter<edm::InputTag> ("METTag");
	MHTJetTag_ = iConfig.getParameter<edm::InputTag> ("MHTJets");
	DeltaPhiJetTag_ = iConfig.getParameter<edm::InputTag> ("DeltaPhiJets");
	
	produces<double>("DeltaPhi1");
	produces<double>("DeltaPhi2");
	produces<double>("DeltaPhi3");
	        produces<double>("DeltaPhiMET1");

        produces<double>("DeltaPhiMET2");
        produces<double>("DeltaPhiMET3");
	produces<double>("DeltaPhiMinMet");
	produces<double>("DeltaPhiMinMHt");
//edm::Handle<JetCorrectorParametersCollection> JetCorParColl;
//const JetCorrector* corrector = 0;	

	/* Examples
 *
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


DeltaPhiDouble::~DeltaPhiDouble()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DeltaPhiDouble::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	double deltaphi1=10;
	double deltaphi2=10;
	double deltaphi3=10;
	edm::Handle< edm::View<reco::Candidate> > MHTJets;
	iEvent.getByLabel(MHTJetTag_,MHTJets);
	reco::MET::LorentzVector mhtLorentz(0,0,0,0);
	if( MHTJets.isValid() ) {
		for(unsigned int i=0; i<MHTJets->size();i++)
		{
			mhtLorentz -=MHTJets->at(i).p4();
		}
	}
	else std::cout<<"DeltaPhiDouble::Invlide MHT Jet Tag: "<<MHTJetTag_.label()<<std::endl;
	reco::MET::LorentzVector metLorentz(0,0,0,0);
        edm::Handle< edm::View<reco::MET> > MET;
        iEvent.getByLabel(metTag_,MET);
        if(MET.isValid() )metLorentz=MET->at(0).p4();
	edm::Handle< edm::View<reco::Candidate> > DeltaPhiJets;
	iEvent.getByLabel(DeltaPhiJetTag_,DeltaPhiJets);
	float minDelPhiMET=99;
	float minDelPhiMHT=99;
	        double deltaphimet1=10;
        double deltaphimet2=10;
        double deltaphimet3=10;
//	float sum=0;
//	float jetres=0.1;
	if( DeltaPhiJets.isValid() ) {
		int count=0;
		for(unsigned int i=0; i<DeltaPhiJets->size();i++)
		{
			
			if(count==3)break;
			if(count==0){ 	
			deltaphi1 = std::abs(reco::deltaPhi(DeltaPhiJets->at(i).phi(),mhtLorentz.phi()));
			 deltaphimet1 = std::abs(reco::deltaPhi(DeltaPhiJets->at(i).phi(),metLorentz.phi()));}
			if(count==1){
		 	deltaphi2 = std::abs(reco::deltaPhi(DeltaPhiJets->at(i).phi(),mhtLorentz.phi()));
                        deltaphimet2 = std::abs(reco::deltaPhi(DeltaPhiJets->at(i).phi(),metLorentz.phi()));
}
			if(count==2){
		 	deltaphi3 = std::abs(reco::deltaPhi(DeltaPhiJets->at(i).phi(),mhtLorentz.phi()));
			deltaphimet3 = std::abs(reco::deltaPhi(DeltaPhiJets->at(i).phi(),metLorentz.phi()));
		}
			count++;
	         	float dpmht=std::abs(reco::deltaPhi(DeltaPhiJets->at(i).phi(),mhtLorentz.phi()));	
		     float dpmet=std::abs(reco::deltaPhi(DeltaPhiJets->at(i).phi(),metLorentz.phi()));
		     if(minDelPhiMHT>dpmht)minDelPhiMHT=dpmht;
		     if(minDelPhiMET>dpmet)minDelPhiMET=dpmet;
	
		    // sum += pow( jetres*(DeltaPhiJets->at(i).px()*DeltaPhiJets->at(i).py() - ), 2);
			//if(count==3) break;
		}
	}
	else std::cout<<"DeltaPhiDouble::Invlide DeltaPhiJets Jet Tag: "<<DeltaPhiJetTag_.label()<<std::endl;
	
	std::auto_ptr<double> htp1(new double(deltaphi1));
	iEvent.put(htp1,"DeltaPhi1");
	std::auto_ptr<double> htp2(new double(deltaphi2));
	iEvent.put(htp2,"DeltaPhi2");
	std::auto_ptr<double> htp3(new double(deltaphi3));
	iEvent.put(htp3,"DeltaPhi3");
	std::auto_ptr<double> metp1(new double(deltaphimet1));
        iEvent.put(metp1,"DeltaPhiMET1");
        std::auto_ptr<double> metp2(new double(deltaphimet2));
        iEvent.put(metp2,"DeltaPhiMET2");
        std::auto_ptr<double> metp3(new double(deltaphimet3));
        iEvent.put(metp3,"DeltaPhiMET3");

	std::auto_ptr<double> metdp(new double(minDelPhiMET));
        iEvent.put(metdp,"DeltaPhiMinMet");	
	std::auto_ptr<double> mhtdp(new double(minDelPhiMHT));
	iEvent.put(mhtdp,"DeltaPhiMinMHt");
}

// ------------ method called once each job just before starting event loop  ------------
void 
DeltaPhiDouble::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DeltaPhiDouble::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
DeltaPhiDouble::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DeltaPhiDouble::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DeltaPhiDouble::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DeltaPhiDouble::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DeltaPhiDouble::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DeltaPhiDouble);
