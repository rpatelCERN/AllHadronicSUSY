#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
// additional headers
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TString.h"
#include "TTree.h"
#include <TFile.h>
#include <vector>
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//
// class declaration
//

class TreeMaker : public edm::EDProducer {
public:
	explicit TreeMaker(const edm::ParameterSet&);
	~TreeMaker();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
	virtual void beginJob() override;
	virtual void produce(edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;
	// ----------member data ---------------------------
	const unsigned int nMaxCandidates_;
	void setBranchVariablesToDefault();
	TString treeName_;
	TTree* tree_;	
	bool debug_;
	// generell event information
	UInt_t runNum_;      
	UInt_t lumiBlockNum_;
	UInt_t evtNum_;
	// any float precision varialbes
	std::vector<edm::InputTag> varsDoubleTags_;
	std::vector<std::string> varsDoubleNames_;
	std::vector<Float_t> varsDouble_;
	// any int precision varialbes
	std::vector<edm::InputTag> varsIntTags_;
	std::vector<std::string> varsIntNames_;
	std::vector<Int_t> varsInt_;
	// any bool precision varialbes
	std::vector<edm::InputTag> varsBoolTags_;
	std::vector<std::string> varsBoolNames_;
	std::vector<UChar_t> varsBool_;
	// any TLorentzVector precision varialbes
	std::vector<edm::InputTag> varsTLorentzVectorTags_;
	std::vector<std::string> varsTLorentzVectorNames_;
	std::vector<TLorentzVector> varsTLorentzVector_;
	// any TLorentzVector precision varialbes
	std::vector<edm::InputTag> vectorTLorentzVectorTags_;
	std::vector<std::string> vectorTLorentzVectorNames_;
	std::vector<std::vector< TLorentzVector> > vectorTLorentzVector_;
	// any reco candidate plus addiation doubles
	std::vector<edm::InputTag> varsRecoCandTags_;
	std::vector<std::string> varsRecoCandNames_;
	std::vector<UShort_t> RecoCandN_;
	std::vector<Float_t*> RecoCandPt_;
	std::vector<Float_t*> RecoCandEta_;
	std::vector<Float_t*> RecoCandPhi_;
	std::vector<Float_t*> RecoCandE_;
	std::vector<TLorentzVector*> RecoCandLorentzVector_;
	std::vector<UShort_t> RecoCandAdditionalBoolVariablesN_,RecoCandAdditionalIntVariablesN_,RecoCandAdditionalFloatVariablesN_;
	std::vector< std::vector<edm::InputTag> > RecoCandAdditionalBoolVariablesTags_;
	std::vector< std::vector<UChar_t*> > RecoCandAdditionalBoolVariables_;
	std::vector< std::vector<edm::InputTag> > RecoCandAdditionalIntVariablesTags_;
	std::vector< std::vector<Int_t*> > RecoCandAdditionalIntVariables_;
	std::vector< std::vector<edm::InputTag> > RecoCandAdditionalFloatVariablesTags_;
	std::vector< std::vector<Float_t*> > RecoCandAdditionalFloatVariables_;
};
