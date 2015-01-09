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
#include <DataFormats/Math/interface/deltaR.h>
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//
// class declaration
//

class RA2TreeMaker : public edm::EDProducer {
public:
	explicit RA2TreeMaker(const edm::ParameterSet&);
	~RA2TreeMaker();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle);
	const reco::Candidate* promtParticle(const reco::Candidate * particle);
	const reco::Candidate* TauFound(const reco::Candidate * particle);
	const reco::Candidate* WFound(const reco::Candidate * particle);
	
private:
	virtual void beginJob() override;
	virtual void produce(edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;
	
	//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	
	// ----------member data ---------------------------
	const unsigned int nMaxCandidates_;
	void setBranchVariablesToDefault();
	
	TString treeName_;
	TTree* tree_;
	bool MC_, QCD_, LostLepton_, UseAll_;
	
	// CMSSW 7
	UShort_t debug_;

	// event information
	UInt_t runNum_;      
	UInt_t lumiBlockNum_;
	UInt_t evtNum_;
	edm::InputTag vertexCollectionTag_;
	UShort_t nVtx_;
	// Any double-precision variables
	std::vector<edm::InputTag> varsDoubleTags_;
	std::vector<std::string> varsDoubleNamesInTree_;
	std::vector<Float_t> varsDouble_;
	// Any boolean variables
	std::vector<edm::InputTag> filterDecisionTags_;
	std::vector<UChar_t> filterDecisions_;
	// main search variables
	edm::InputTag ra2JetsTag_;
	edm::InputTag metTag_;
	std::string   btagname_;
	double        btagvalue_;
	Float_t ht_, mht_, mhtEta_, mhtPhi_;
	Float_t minHT_, minMHT_, maxEtaHTJets_, maxEtaMHTJets_, minPTHTJets_, minPTMHTJets_;
	int minNJets_;
	Float_t metPt_, metEta_, metPhi_;
	
	UShort_t nJets_, BTags_;
	UShort_t nIsoLeptons_, nIsoTracks_;
	Float_t jet1Pt_, jet2Pt_, jet3Pt_;
	Float_t jet1Eta_, jet2Eta_, jet3Eta_;
	Float_t deltaPhi1_, deltaPhi2_, deltaPhi3_;
	// ra2jets vector
	std::vector<Float_t> ra2Jetsht_, ra2Jetsmht_;
	std::vector<UShort_t> ra2JetsnJets_, ra2JetsBTags_;
	std::vector<Float_t> ra2Jetsjet1Pt_, ra2Jetsjet2Pt_, ra2Jetsjet3Pt_;
	std::vector<Float_t> ra2Jetsjet1Eta_, ra2Jetsjet2Eta_, ra2Jetsjet3Eta_;
	std::vector<Float_t> ra2JetsdeltaPhi1_, ra2JetsdeltaPhi2_, ra2JetsdeltaPhi3_;
	std::vector<edm::InputTag> ra2JetsCollectionInputTag_;
	std::vector<std::string>   ra2JetsCollectionNameInTree_;
	std::vector<UShort_t> ra2JetsN_;
	std::vector<Float_t*> ra2JetsPt_;
	std::vector<Float_t*> ra2JetsEta_;
	std::vector<Float_t*> ra2JetsPhi_;
	std::vector<Float_t*> ra2JetsE_;
	std::vector<std::string> ra2JetsBTagInputTag_;
	std::vector<double> ra2JetsBTagValueInput_;
	std::vector<Float_t*> ra2JetsBTagValue_;
	std::vector<UShort_t*> ra2JetsBTag_;
	// lepton stuff
	std::vector<edm::InputTag> leptonTag_;
	std::vector<std::string>   leptonTagName_;
	std::vector<UShort_t> leptonN_;
	std::vector<Float_t*> leptonPt_;
	std::vector<Float_t*> leptonEta_;
	std::vector<Float_t*> leptonPhi_;
	std::vector<Float_t*> leptonE_;
	
	// isotracks stuff
	std::vector<edm::InputTag> IsoTrackTag_;
	std::vector<std::string>   IsoTrackTagName_;
	std::vector<UShort_t> isoTrackN_;
	std::vector<Float_t*> isoTrackPt_;
	std::vector<Float_t*> isoTrackEta_;
	std::vector<Float_t*> isoTrackPhi_;
	std::vector<Float_t*> isoTrackE_;
	// MC gen information CMSSW 7
	edm::InputTag prunedGenToken_;
	edm::InputTag packedGenToken_;
	UShort_t GenWNum_;
	Float_t GenWPt_[25];
	Float_t GenWPhi_[25];
	Float_t GenWEta_[25];
	UShort_t GenMuNum_;
	UShort_t GenMuFromTau_[25];
	Float_t GenMuPt_[25];
	Float_t GenMuPhi_[25];
	Float_t GenMuEta_[25];
	UShort_t GenElecNum_;
	UShort_t GenElecFromTau_[25];
	Float_t GenElecPt_[25];
	Float_t GenElecPhi_[25];
	Float_t GenElecEta_[25];
	UShort_t GenTauNum_;
	Float_t GenTauPt_[25];
	Float_t GenTauPhi_[25];
	Float_t GenTauEta_[25];
	edm::InputTag genra2JetsTag_;
	Float_t genht_, genmht_;
	UShort_t gennJets_;
	Float_t genjet1Pt_, genjet2Pt_, genjet3Pt_;
	Float_t genjet1Eta_, genjet2Eta_, genjet3Eta_;
	Float_t gendeltaPhi1_, gendeltaPhi2_, gendeltaPhi3_;
	UShort_t GenJetNum_;
	Float_t GenJetPt_[200];
	Float_t GenJetPhi_[200];
	Float_t GenJetEta_[200];
	Float_t GenJetE_[200];

	
};
