// -*- C++ -*-
//
// Package:    AllHadronicSUSY/RA2TreeMaker
// Class:      RA2TreeMaker
// 
/**\class RA2TreeMaker RA2TreeMaker.cc AllHadronicSUSY/RA2TreeMaker/plugins/RA2TreeMaker.cc
 * 
 * Description: [one line class summary]
 * 
 * Implementation:
 *     [Notes on implementation]
 */
//
// Original Author:  Arne-Rasmus Draeger
//         Created:  Fri, 19 Sep 2014 13:48:35 GMT
//
//
#include "AllHadronicSUSY/RA2TreeMaker/interface/RA2TreeMaker.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <memory>

// system include files

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
RA2TreeMaker::RA2TreeMaker(const edm::ParameterSet& iConfig)
: nMaxCandidates_(200), tree_(0)
{
  // generell
  treeName_ = iConfig.getParameter<std::string>("TreeName");
  vertexCollectionTag_ = iConfig.getParameter<edm::InputTag>("VertexCollection");
  varsDoubleTags_ = iConfig.getParameter< std::vector<edm::InputTag> >("VarsDouble");
  varsDoubleNamesInTree_= iConfig.getParameter< std::vector<std::string> >  ("VarsDoubleNamesInTree");
  filterDecisionTags_ = iConfig.getParameter< std::vector<edm::InputTag> >("Filters");
  MC_ = iConfig.getParameter<bool> ("MC");
  QCD_ = iConfig.getParameter<bool> ("QCD");
  LostLepton_ = iConfig.getParameter<bool> ("LostLepton");
  if(!QCD_ && !LostLepton_) UseAll_=true;
  else UseAll_=false;
  debug_ = iConfig.getParameter<bool> ("debug");
  if(MC_)std::cout<<"Running on MC. Lepton gen information will be stored"<<std::endl;
  // cmssw 7 tags
  prunedGenToken_ = iConfig.getParameter<edm::InputTag> ("prunedGenParticles");
  packedGenToken_ = iConfig.getParameter<edm::InputTag>("packedGenParticles");
  ra2JetsTag_ = iConfig.getParameter<edm::InputTag>("RA2DefaultJetsTag");
  btagname_    = iConfig.getParameter<std::string>  ("bTagName");
  btagvalue_   = iConfig.getParameter<double>       ("bTagValue");
  metTag_ = iConfig.getParameter<edm::InputTag> ("METTag");
  // selection criteria
  minHT_   = iConfig.getParameter<double>       ("MinHT");
  minMHT_   = iConfig.getParameter<double>       ("MinMHT");
  //minNJets_   = iConfig.getParameter<int>       ("MinNJets");
  maxEtaHTJets_ = iConfig.getParameter<double>       ("MaxEtaHTJets");
  minPTHTJets_ = iConfig.getParameter<double>       ("MinPTHTJets");
  maxEtaMHTJets_ = iConfig.getParameter<double>       ("MaxEtaMHTJets");
  minPTMHTJets_ = iConfig.getParameter<double>       ("MinPTMHTJets");
  std::cout<<"RA2TreeMaker::Basic CUTs: "<<std::endl;
  std::cout<<"minHT_"<<minHT_<<std::endl;
  std::cout<<"minMHT_"<<minMHT_<<std::endl;
  std::cout<<"maxEtaHTJets_"<<maxEtaHTJets_<<std::endl;
  std::cout<<"minPTHTJets_"<<minPTHTJets_<<std::endl;
  std::cout<<"maxEtaMHTJets_"<<maxEtaMHTJets_<<std::endl;
  std::cout<<"minPTMHTJets_"<<minPTMHTJets_<<std::endl;
  ra2JetsCollectionInputTag_          	= 	iConfig.getParameter< std::vector<edm::InputTag> >("ra2JetsCollectionInputTag");
  ra2JetsCollectionNameInTree_        	= 	iConfig.getParameter< std::vector<std::string> >  ("ra2JetsCollectionNameInTree");
  ra2JetsBTagInputTag_        		= 	iConfig.getParameter< std::vector<std::string> >  ("ra2JetsBTagInputTag");
  ra2JetsBTagValueInput_        		= 	iConfig.getParameter< std::vector<double> >  ("ra2JetsBTagValueInput_");
  leptonTag_ 				= 	iConfig.getParameter< std::vector<edm::InputTag> >("LeptonTag");
  leptonTagName_        		= 	iConfig.getParameter< std::vector<std::string> >  ("LeptonTagName");
  IsoTrackTag_ 				= 	iConfig.getParameter< std::vector<edm::InputTag> >("IsoTrackTag");
  IsoTrackTagName_        			= 	iConfig.getParameter< std::vector<std::string> >  ("IsoTrackTagName");
  genra2JetsTag_				=	iConfig.getParameter<edm::InputTag>  ("GenJetTag");
  varsDouble_ = std::vector<Float_t>(varsDoubleTags_.size(),1.);
  filterDecisions_ = std::vector<UChar_t>(filterDecisionTags_.size(),0);
  for(unsigned int i = 0; i < ra2JetsCollectionInputTag_.size(); ++i) {
    ra2JetsN_.push_back(0);
    ra2JetsPt_.push_back (new Float_t[nMaxCandidates_]);
    ra2JetsEta_.push_back(new Float_t[nMaxCandidates_]);
    ra2JetsPhi_.push_back(new Float_t[nMaxCandidates_]);
    ra2JetsE_.push_back  (new Float_t[nMaxCandidates_]);
    ra2JetsBTagValue_.push_back (new Float_t[nMaxCandidates_]);
    ra2JetsBTag_.push_back(new UShort_t[nMaxCandidates_]);
    ra2Jetsht_.push_back(0);
    ra2Jetsmht_.push_back(0);
    ra2JetsnJets_.push_back(0);
    ra2JetsBTags_.push_back(0);
    ra2Jetsjet1Pt_.push_back(0); ra2Jetsjet2Pt_.push_back(0); ra2Jetsjet3Pt_.push_back(0);
    ra2Jetsjet1Eta_.push_back(0); ra2Jetsjet2Eta_.push_back(0); ra2Jetsjet3Eta_.push_back(0);
    ra2JetsdeltaPhi1_.push_back(0); ra2JetsdeltaPhi2_.push_back(0); ra2JetsdeltaPhi3_.push_back(0);
  }
  // ra2Leptons
  for(unsigned int i = 0; i < leptonTag_.size(); ++i) {
    leptonN_.push_back(0);
    leptonPt_.push_back (new Float_t[25]);
    leptonEta_.push_back(new Float_t[25]);
    leptonPhi_.push_back(new Float_t[25]);
    leptonE_.push_back  (new Float_t[25]);
  }
  // iso tracks
  for(unsigned int i = 0; i < leptonTag_.size(); ++i) {
    isoTrackN_.push_back(0);
    isoTrackPt_.push_back (new Float_t[25]);
    isoTrackEta_.push_back(new Float_t[25]);
    isoTrackPhi_.push_back(new Float_t[25]);
    isoTrackE_.push_back  (new Float_t[25]);
  }
}


RA2TreeMaker::~RA2TreeMaker()
{
  for (unsigned int i=0; i < ra2JetsCollectionInputTag_.size();i++)
  {
    delete []   ra2JetsPt_.at(i);
    delete []   ra2JetsEta_.at(i);
    delete []   ra2JetsPhi_.at(i);
    delete []   ra2JetsE_.at(i);
    delete []   ra2JetsBTag_.at(i);
    delete []   ra2JetsBTagValue_.at(i);
  }
  for (unsigned int i=0; i < leptonTag_.size();i++)
  {
    delete []   leptonPt_.at(i);
    delete []   leptonEta_.at(i);
    delete []   leptonPhi_.at(i);
    delete []   leptonE_.at(i);
  }
  for (unsigned int i=0; i < IsoTrackTag_.size();i++)
  {
    delete []   isoTrackPt_.at(i);
    delete []   isoTrackEta_.at(i);
    delete []   isoTrackPhi_.at(i);
    delete []   isoTrackE_.at(i);
  }
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
RA2TreeMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace pat;
  setBranchVariablesToDefault();
  
  // Event information
  edm::EventAuxiliary aux = iEvent.eventAuxiliary();
  runNum_       = aux.run();
  lumiBlockNum_ = aux.luminosityBlock();
  evtNum_       = aux.event();
  
  // Number of vertices
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexCollectionTag_,vertices);
  if( vertices.isValid() ) {
    nVtx_ = vertices->size();
  }
  // Double-precision variables
  for(unsigned int i = 0; i < varsDoubleTags_.size(); ++i) {
    edm::Handle<double> var;
    iEvent.getByLabel(varsDoubleTags_.at(i),var);
    if( var.isValid() ) {
      varsDouble_.at(i) = *var;
    }
  }
  for(unsigned int i = 0; i < filterDecisionTags_.size(); ++i) {
    edm::Handle<bool> dec;
    iEvent.getByLabel(filterDecisionTags_.at(i),dec);
    if( dec.isValid() ) {
      if( *dec ) filterDecisions_.at(i) = 1;
      else filterDecisions_.at(i) = 0;
    }
  }
  edm::Handle< edm::View<reco::MET> > MET;
  iEvent.getByLabel(metTag_,MET);  
  metPt_=MET->at(0).pt();
  metEta_=MET->at(0).eta();
  metPhi_=MET->at(0).phi();
  
  edm::Handle< edm::View<pat::Jet> > ra2JetsCands;
  iEvent.getByLabel(ra2JetsTag_,ra2JetsCands);
  reco::MET::LorentzVector mhtLorentz(0,0,0,0);
  for (unsigned int i = 0; i < ra2JetsCands->size();i++)
  {
    
    // select HT Jets
    if(ra2JetsCands->at(i).pt() > minPTHTJets_ && abs(ra2JetsCands->at(i).eta() ) <maxEtaHTJets_ )/// check values!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
      ht_+=ra2JetsCands->at(i).pt();
      nJets_++;
      if(ra2JetsCands->at(i).bDiscriminator(btagname_) >btagvalue_) 
      {
//	      const reco::TrackRefVector *tracks = &ra2JetsCands->at(i).associatedTracks();
//	      std::cout<<"BTag associated tracks size:"<<tracks->size()<<std::endl;
	      BTags_++;
      }
    }
    if(ra2JetsCands->at(i).pt() > minPTMHTJets_ && abs(ra2JetsCands->at(i).eta() ) <maxEtaMHTJets_ ) /// check values!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
      mhtLorentz -=ra2JetsCands->at(i).p4();
    }
  }
  mht_=mhtLorentz.pt();
  mhtEta_=mhtLorentz.eta();
  mhtPhi_=mhtLorentz.phi();
  int count=0;
  for (unsigned int i = 0; i < ra2JetsCands->size();i++)
  {
    if( ra2JetsCands->at(i).pt() > minPTMHTJets_ &&abs(ra2JetsCands->at(i).eta() ) <maxEtaMHTJets_ )
    {
      if(count==0)
      {
	deltaPhi1_ = std::abs(reco::deltaPhi(ra2JetsCands->at(i).phi(),mhtLorentz.phi()));
	jet1Pt_ = ra2JetsCands->at(i).pt();
	jet1Eta_ = ra2JetsCands->at(i).eta();
	
      }
      if(count==1)
      {
	deltaPhi2_ = std::abs(reco::deltaPhi(ra2JetsCands->at(i).phi(),mhtLorentz.phi()));
	jet2Pt_ = ra2JetsCands->at(i).pt();
	jet2Eta_ = ra2JetsCands->at(i).eta();
      }
      if(count==2)
      {
	deltaPhi3_ = std::abs(reco::deltaPhi(ra2JetsCands->at(i).phi(),mhtLorentz.phi()));
	jet3Pt_ = ra2JetsCands->at(i).pt();
	jet3Eta_ = ra2JetsCands->at(i).eta();
      }
      count++;
      if(count==3) break;
    }
  }
  count=0;
  for(unsigned int i = 0; i < ra2JetsCollectionInputTag_.size(); ++i) {
    edm::Handle< edm::View<pat::Jet> > ra2JetsCands;
    iEvent.getByLabel(ra2JetsCollectionInputTag_.at(i),ra2JetsCands);
    if( ra2JetsCands.isValid() ) 
    {
      reco::MET::LorentzVector mhtLorentz2(0,0,0,0);
      ra2JetsN_[i]=ra2JetsCands->size();
      
      for(unsigned int j = 0; j < ra2JetsN_[i]; ++j) 
      {
	if(j<nMaxCandidates_)
	{
	  ra2JetsPt_.at(i)[j] = ra2JetsCands->at(j).pt();
	  ra2JetsEta_.at(i)[j] = ra2JetsCands->at(j).eta();
	  ra2JetsPhi_.at(i)[j] = ra2JetsCands->at(j).phi();
	  ra2JetsE_.at(i)[j] = ra2JetsCands->at(j).energy();
	  ra2JetsBTagValue_.at(i)[j]=ra2JetsCands->at(j).bDiscriminator(ra2JetsBTagInputTag_.at(i).c_str());
	  if(ra2JetsCands->at(j).bDiscriminator(ra2JetsBTagInputTag_.at(i).c_str()) >ra2JetsBTagValueInput_.at(i)) ra2JetsBTag_.at(i)[j]=1;
	  else ra2JetsBTag_.at(i)[j]=0;
	  // compute ht mht njets and deltaphi
	  if(ra2JetsCands->at(i).pt() > minPTHTJets_ && abs(ra2JetsCands->at(i).eta() ) <maxEtaHTJets_ ) // ht definitions
	  {
	    ra2Jetsht_.at(i) +=ra2JetsCands->at(j).pt();
	    ra2JetsnJets_.at(i)++;
	    if(ra2JetsCands->at(j).bDiscriminator(ra2JetsBTagInputTag_.at(i).c_str()) >ra2JetsBTagValueInput_.at(i)) ra2JetsBTags_.at(i)++;
	  }
	  if( ra2JetsCands->at(i).pt() > minPTMHTJets_ &&abs(ra2JetsCands->at(i).eta() ) <maxEtaMHTJets_ ) /// check values!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  {
	    mhtLorentz2 -=ra2JetsCands->at(j).p4();
	  }
	}
	
      }
      
      int count=0;
      for(unsigned int j = 0; j<ra2JetsN_[i];j++)
      {
	if(count==0)
	{
	  ra2JetsdeltaPhi1_.at(i) = std::abs(reco::deltaPhi(ra2JetsCands->at(j).phi(),mhtLorentz2.phi()));
	  ra2Jetsjet1Pt_.at(i) = ra2JetsCands->at(j).pt();
	  ra2Jetsjet1Eta_.at(i) = ra2JetsCands->at(j).eta();
	}
	if(count==1)
	{
	  ra2JetsdeltaPhi2_.at(i) = std::abs(reco::deltaPhi(ra2JetsCands->at(j).phi(),mhtLorentz2.phi()));
	  ra2Jetsjet2Pt_.at(i) = ra2JetsCands->at(j).pt();
	  ra2Jetsjet2Eta_.at(i) = ra2JetsCands->at(j).eta();
	}
	if(count==2)
	{
	  ra2JetsdeltaPhi3_.at(i) = std::abs(reco::deltaPhi(ra2JetsCands->at(j).phi(),mhtLorentz2.phi()));
	  ra2Jetsjet3Pt_.at(i) = ra2JetsCands->at(j).pt();
	  ra2Jetsjet3Eta_.at(i) = ra2JetsCands->at(j).eta();
	}
	count++;
	if(count==3) break;
	
      }
      ra2Jetsmht_.at(i)=mhtLorentz2.pt();
    }
  }
  // ra2Leptons
  for(unsigned int i = 0; i < leptonTag_.size(); ++i) {
    edm::Handle< edm::View<reco::Candidate> > cands;
    iEvent.getByLabel(leptonTag_.at(i),cands);
    
    if( cands.isValid() ) {
      std::string name = leptonTag_.at(i).label();
      if(leptonTagName_.size()== leptonTag_.size() )
      {
	name = leptonTagName_.at(i);
      }
      if(name.find("RecoIso") !=std::string::npos)
      {
	nIsoLeptons_+=cands->size();
      }
      leptonN_[i] = cands->size();
      for(unsigned int j = 0; j < cands->size(); ++j) {
	leptonPt_.at(i)[j] = cands->at(j).pt();
	leptonEta_.at(i)[j] = cands->at(j).eta();
	leptonPhi_.at(i)[j] = cands->at(j).phi();
	leptonE_.at(i)[j] = cands->at(j).energy();
      }
    }
  }
  // Isolated Tracks
  for(unsigned int i = 0; i < IsoTrackTag_.size(); ++i) {
    edm::Handle< edm::View<pat::PackedCandidate> > cands;
    iEvent.getByLabel(IsoTrackTag_.at(i),cands);
    
    if( cands.isValid() ) {
      std::string name = IsoTrackTag_.at(i).label();
      if(IsoTrackTagName_.size()== IsoTrackTag_.size() )
      {
	name = IsoTrackTagName_.at(i);
      }
	
      if(name.find("SelectedIsoTracks") !=std::string::npos)
      {
	nIsoTracks_+=cands->size();
      }
      isoTrackN_[i] = cands->size();
      for(unsigned int j = 0; j < cands->size(); ++j) {
	isoTrackPt_.at(i)[j] = cands->at(j).pt();
	isoTrackEta_.at(i)[j] = cands->at(j).eta();
	isoTrackPhi_.at(i)[j] = cands->at(j).phi();
	isoTrackE_.at(i)[j] = cands->at(j).energy();
      }
    }
  }
  // extract gen information on leptons
  if(MC_)
  {
    // store gen jets:
    edm::Handle< edm::View<reco::GenJet> > ra2GenJetsCands;
    iEvent.getByLabel(genra2JetsTag_,ra2GenJetsCands);
    reco::MET::LorentzVector genmhtLorentz(0,0,0,0);
    GenJetNum_=ra2GenJetsCands->size();
    for (unsigned int i = 0; i < ra2GenJetsCands->size();i++)
    {
      GenJetPt_[i]=ra2GenJetsCands->at(i).pt();
      GenJetEta_[i]=ra2GenJetsCands->at(i).eta();
      GenJetPhi_[i]=ra2GenJetsCands->at(i).phi();
      GenJetE_[i]=ra2GenJetsCands->at(i).energy();
      std::vector <const GenParticle*> Constituents = ra2GenJetsCands->at(i).getGenConstituents();
      //std::cout<<"GenJetConstituentsSize: "<<Constituents.size()<<std::endl;
      for(unsigned int ii=0; ii < Constituents.size();ii++)
      {
	if(true) std::cout<<"GenParticlesofJet["<<i<<"] constituent["<<ii<<"] pdgId:"<<Constituents[ii]->pdgId()<<std::endl;
      }
      // select HT Jets
      if(ra2GenJetsCands->at(i).pt() > minPTHTJets_ && abs(ra2GenJetsCands->at(i).eta() ) <maxEtaHTJets_ )/// check values!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      {
	genht_+=ra2GenJetsCands->at(i).pt();
	gennJets_++;
      }
      if(ra2GenJetsCands->at(i).pt() > minPTMHTJets_ && abs(ra2GenJetsCands->at(i).eta() ) <maxEtaMHTJets_ ) /// check values!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      {
	genmhtLorentz -=ra2GenJetsCands->at(i).p4();
      }
    }
    genmht_=genmhtLorentz.pt();
    int count=0;
    for (unsigned int i = 0; i < ra2GenJetsCands->size();i++)
    {
      if( ra2GenJetsCands->at(i).pt() > minPTMHTJets_ &&abs(ra2GenJetsCands->at(i).eta() ) <maxEtaMHTJets_ )
      {
	if(count==0)
	{
	  gendeltaPhi1_ = std::abs(reco::deltaPhi(ra2GenJetsCands->at(i).phi(),genmhtLorentz.phi()));
	  genjet1Pt_ = ra2GenJetsCands->at(i).pt();
	  genjet1Eta_ = ra2GenJetsCands->at(i).eta();
	  
	}
	if(count==1)
	{
	  gendeltaPhi2_ = std::abs(reco::deltaPhi(ra2GenJetsCands->at(i).phi(),genmhtLorentz.phi()));
	  genjet2Pt_ = ra2GenJetsCands->at(i).pt();
	  genjet2Eta_ = ra2GenJetsCands->at(i).eta();
	}
	if(count==2)
	{
	  gendeltaPhi3_ = std::abs(reco::deltaPhi(ra2GenJetsCands->at(i).phi(),genmhtLorentz.phi()));
	  genjet3Pt_ = ra2GenJetsCands->at(i).pt();
	  genjet3Eta_ = ra2GenJetsCands->at(i).eta();
	}
	count++;
	if(count==3) break;
      }
    }


    // Packed particles are all the status 1, so usable to remake jets
    // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
    Handle<edm::View<reco::GenParticle> > pruned;
    iEvent.getByLabel(prunedGenToken_,pruned);
    for(size_t i=0; i<pruned->size();i++)
    {
      if( abs((*pruned)[i].pdgId() ) == 24 && (*pruned)[i].status()==22) // W from the hard interaction 
      {
	//std::cout<<"W found with daughters N: "<<(*pruned)[i].numberOfDaughters()<<std::endl;
	const Candidate * FinalW = WFound(&(*pruned)[i]);
	size_t wDaugthers = FinalW->numberOfDaughters();
	for(size_t ii=0;ii< wDaugthers; ii++)
	{
	  if(abs( FinalW->daughter(0)->pdgId() )==24 )// true if the w has another daugther
	  {
	    std::cout<<"NEVERTRUE";
	    FinalW = FinalW->daughter(0);
	  }
	}
	//std::cout<<"W of W with daughters N: "<<FinalW->numberOfDaughters()<<std::endl;
	for(size_t ii=0;ii< FinalW->numberOfDaughters(); ii++)
	{
	  if(abs( FinalW->daughter(ii)->pdgId() )==24 )// true if the w has another daugther
	  {
	    FinalW = (*pruned)[i].daughter(ii);
	    std::cout<<"NEVERTRUE2";
	  }
	}
	
	// if(FinalW == nullptr) FinalW = &(*pruned)[i];
	wDaugthers =FinalW->numberOfDaughters();
	GenWNum_++;
	GenWPt_[GenWNum_-1]=FinalW->pt();
	GenWEta_[GenWNum_-1]=FinalW->eta();
	GenWPhi_[GenWNum_-1]=FinalW->phi();
	for(size_t ii=0;ii< wDaugthers; ii++)
	{
	  //std::cout<<"New method daugther["<<ii<<"] id:"<<FinalW->daughter(ii)->pdgId()<<" ";
	  if(abs(FinalW->daughter(ii)->pdgId())== 11) // electron
	  {
	    GenElecNum_++;
	    GenElecFromTau_[GenElecNum_-1]=0;
	    GenElecPt_[GenElecNum_-1]=FinalW->daughter(ii)->pt();
	    GenElecEta_[GenElecNum_-1]=FinalW->daughter(ii)->eta();
	    GenElecPhi_[GenElecNum_-1]=FinalW->daughter(ii)->phi();
	  }
	  if(abs(FinalW->daughter(ii)->pdgId())== 13) // muon
	  {
	    GenMuNum_++;
	    GenMuFromTau_[GenMuNum_-1]=0;
	    GenMuPt_[GenMuNum_-1]=FinalW->daughter(ii)->pt();
	    GenMuEta_[GenMuNum_-1]=FinalW->daughter(ii)->eta();
	    GenMuPhi_[GenMuNum_-1]=FinalW->daughter(ii)->phi();
	  }
	  if(abs(FinalW->daughter(ii)->pdgId())== 15) // tau
	  {
	    GenTauNum_++;
	    GenTauPt_[GenTauNum_-1]=FinalW->daughter(ii)->pt();
	    GenTauEta_[GenTauNum_-1]=FinalW->daughter(ii)->eta();
	    GenTauPhi_[GenTauNum_-1]=FinalW->daughter(ii)->phi();
	    const Candidate * FinalTauDecay = TauFound(FinalW->daughter(ii));
	    // std::cout<<"FinalTauID:"<<FinalTauDecay->pdgId()<<" ";
	    for(size_t iii=0; iii<FinalTauDecay->numberOfDaughters();iii++)
	    {
	      //std::cout<<"FinalTauDaugther["<<iii<<"] pdgId:"<<FinalTauDecay->daughter(iii)->pdgId()<<", ";
	      if(abs(FinalTauDecay->daughter(iii)->pdgId())== 11) // electron
	      {
		GenElecNum_++;
		GenElecFromTau_[GenElecNum_-1]=1;
		GenElecPt_[GenElecNum_-1]=FinalTauDecay->daughter(iii)->pt();
		GenElecEta_[GenElecNum_-1]=FinalTauDecay->daughter(iii)->eta();
		GenElecPhi_[GenElecNum_-1]=FinalTauDecay->daughter(iii)->phi();
	      }
	      if(abs(FinalTauDecay->daughter(iii)->pdgId())== 13) // muon
	      {
		GenMuNum_++;
		GenMuFromTau_[GenMuNum_-1]=1;
		GenMuPt_[GenMuNum_-1]=FinalTauDecay->daughter(iii)->pt();
		GenMuEta_[GenMuNum_-1]=FinalTauDecay->daughter(iii)->eta();
		GenMuPhi_[GenMuNum_-1]=FinalTauDecay->daughter(iii)->phi();
	      }
	    }
	    // std::cout<<std::endl;
	  }
	  
	}
	// std::cout<<std::endl;
      }
    }
  //  std::cout<<"___________________________"<<std::endl;
  }
  // Fill variables into tree
  tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
RA2TreeMaker::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ) {
    throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
  }
  tree_ = fs->make<TTree>(treeName_,treeName_);  
  tree_->SetAutoSave(10000000000);
  tree_->SetAutoFlush(1000000);
  tree_->Branch("RunNum",&runNum_,"RunNum/i");
  tree_->Branch("LumiBlockNum",&lumiBlockNum_,"LumiBlockNum/i");
  tree_->Branch("EvtNum",&evtNum_,"EvtNum/i");
  tree_->Branch("NVtx",&nVtx_,"NVtx/s");
  for(unsigned int i = 0; i < varsDouble_.size(); ++i) {
    TString name = varsDoubleTags_.at(i).label();
    if( varsDoubleNamesInTree_.size() == varsDoubleTags_.size() ) {
      name = varsDoubleNamesInTree_.at(i);
    }
    tree_->Branch(name,&(varsDouble_.at(i)),name+"/F");
  }
  if(LostLepton_ || UseAll_)
  {
  tree_->Branch("HT",&ht_,"HT/F");
  tree_->Branch("MHT",&mht_,"MHT/F");
  tree_->Branch("MHTEta",&mhtEta_,"MHTEta/F");
  tree_->Branch("MHTPhi",&mhtPhi_,"MHTPhi/F");
  tree_->Branch("NJets",&nJets_,"NJets/s");
  tree_->Branch("BTags",&BTags_,"BTags/s");
  }
  tree_->Branch("Leptons",&nIsoLeptons_,"Leptons/s");
  tree_->Branch("IsolatedTracks",&nIsoTracks_,"IsolatedTracks/s");
  tree_->Branch("METPt",&metPt_,"METPt/F");
  tree_->Branch("METEta",&metEta_,"METEta/F");
  tree_->Branch("METPhi",&metPhi_,"METPhi/F");
  if(LostLepton_ || UseAll_)
  {
  tree_->Branch("Jet1Pt",&jet1Pt_,"Jet1Pt/F");
  tree_->Branch("Jet2Pt",&jet2Pt_,"Jet2Pt/F");
  tree_->Branch("Jet3Pt",&jet3Pt_,"Jet3Pt/F");
  tree_->Branch("Jet1Eta",&jet1Eta_,"Jet1Eta/F");
  tree_->Branch("Jet2Eta",&jet2Eta_,"Jet2Eta/F");
  tree_->Branch("Jet3Eta",&jet3Eta_,"Jet3Eta/F");
  tree_->Branch("DeltaPhi1",&deltaPhi1_,"DeltaPhi1/F");
  tree_->Branch("DeltaPhi2",&deltaPhi2_,"DeltaPhi2/F");
  tree_->Branch("DeltaPhi3",&deltaPhi3_,"DeltaPhi3/F");
  }
  for(unsigned int i = 0; i < ra2JetsCollectionInputTag_.size(); ++i) {
    std::string name = ra2JetsCollectionInputTag_.at(i).label();
    if( i < ra2JetsCollectionNameInTree_.size() ) {
      name = ra2JetsCollectionNameInTree_.at(i);
    }
    if(LostLepton_ || UseAll_)
    {
    tree_->Branch((name+"HT").c_str(),&(ra2Jetsht_.at(i)),(name+"HT/F").c_str());
    tree_->Branch((name+"MHT").c_str(),&(ra2Jetsmht_.at(i)),(name+"MHT/F").c_str());
    tree_->Branch((name+"NJets").c_str(),&(ra2JetsnJets_.at(i)),(name+"NJets/s").c_str());
    tree_->Branch((name+"BTags").c_str(),&(ra2JetsBTags_.at(i)),(name+"BTags/s").c_str());
    tree_->Branch((name+"Jet1Pt").c_str(),&(ra2Jetsjet1Pt_.at(i)),(name+"Jet1Pt/F").c_str());
    tree_->Branch((name+"Jet2Pt").c_str(),&(ra2Jetsjet2Pt_.at(i)),(name+"Jet2Pt/F").c_str());
    tree_->Branch((name+"Jet3Pt").c_str(),&(ra2Jetsjet3Pt_.at(i)),(name+"Jet3Pt/F").c_str());
    tree_->Branch((name+"Jet1Eta").c_str(),&(ra2Jetsjet1Eta_.at(i)),(name+"Jet1Eta/F").c_str());
    tree_->Branch((name+"Jet2Eta").c_str(),&(ra2Jetsjet2Eta_.at(i)),(name+"Jet2Eta/F").c_str());
    tree_->Branch((name+"Jet3Eta").c_str(),&(ra2Jetsjet3Eta_.at(i)),(name+"Jet3Eta/F").c_str());
    tree_->Branch((name+"DeltaPhi1").c_str(),&(ra2JetsdeltaPhi1_.at(i)),(name+"DeltaPhi1/F").c_str());
    tree_->Branch((name+"DeltaPhi2").c_str(),&(ra2JetsdeltaPhi2_.at(i)),(name+"DeltaPhi2/F").c_str());
    tree_->Branch((name+"DeltaPhi3").c_str(),&(ra2JetsdeltaPhi3_.at(i)),(name+"DeltaPhi3/F").c_str());
    }
    // jet information
    tree_->Branch((name+"Num").c_str(),&(ra2JetsN_.at(i)),(name+"Num/s").c_str());
    tree_->Branch((name+"Pt").c_str(), ra2JetsPt_.at(i), (name+"Pt["+name+"Num]/F").c_str());
    tree_->Branch((name+"Eta").c_str(),ra2JetsEta_.at(i),(name+"Eta["+name+"Num]/F").c_str());
    tree_->Branch((name+"Phi").c_str(),ra2JetsPhi_.at(i),(name+"Phi["+name+"Num]/F").c_str());
    tree_->Branch((name+"E").c_str(),  ra2JetsE_.at(i),  (name+"E["+name+"Num]/F").c_str());
    tree_->Branch((name+"BTagValue").c_str(),  ra2JetsBTagValue_.at(i),  (name+"BTagValue["+name+"Num]/F").c_str());
    tree_->Branch((name+"BTag").c_str(),  ra2JetsBTag_.at(i),  (name+"BTag["+name+"Num]/s").c_str());
  }
  // ra2Leptons
  if(LostLepton_ || UseAll_)
  {
  for(unsigned int i = 0; i < leptonTag_.size(); ++i) {
    std::string name = leptonTag_.at(i).label();
    if( i < leptonTagName_.size() ) {
      name = leptonTagName_.at(i);
    }
    tree_->Branch((name+"Num").c_str(),&(leptonN_[i]),(name+"Num/s").c_str());
    tree_->Branch((name+"Pt").c_str(), leptonPt_.at(i), (name+"Pt["+name+"Num]/F").c_str());
    tree_->Branch((name+"Eta").c_str(),leptonEta_.at(i),(name+"Eta["+name+"Num]/F").c_str());
    tree_->Branch((name+"Phi").c_str(),leptonPhi_.at(i),(name+"Phi["+name+"Num]/F").c_str());
    tree_->Branch((name+"E").c_str(),  leptonE_.at(i),  (name+"E["+name+"Num]/F").c_str());
  }
  for(unsigned int i = 0; i < IsoTrackTag_.size(); ++i) {
    std::string name = IsoTrackTag_.at(i).label();
    if( i < IsoTrackTagName_.size() ) {
      name = IsoTrackTagName_.at(i);
    }
    tree_->Branch((name+"Num").c_str(),&(isoTrackN_[i]),(name+"Num/s").c_str());
    tree_->Branch((name+"Pt").c_str(), isoTrackPt_.at(i), (name+"Pt["+name+"Num]/F").c_str());
    tree_->Branch((name+"Eta").c_str(),isoTrackEta_.at(i),(name+"Eta["+name+"Num]/F").c_str());
    tree_->Branch((name+"Phi").c_str(),isoTrackPhi_.at(i),(name+"Phi["+name+"Num]/F").c_str());
    tree_->Branch((name+"E").c_str(),  isoTrackE_.at(i),  (name+"E["+name+"Num]/F").c_str());
  }
  }
  if(MC_)
  {
    if(LostLepton_ || UseAll_)
    {
    tree_->Branch("GenWNum",&GenWNum_,"GenWNum/s");
    tree_->Branch("GenWPt", GenWPt_,"GenWPt[GenWNum]/F");
    tree_->Branch("GenWEta", GenWEta_,"GenWEta[GenWNum]/F");
    tree_->Branch("GenWPhi", GenWPhi_,"GenWPhi[GenWNum]/F");
    tree_->Branch("GenMuNum",&GenMuNum_,"GenMuNum/s");
    tree_->Branch("GenMuFromTau",GenMuFromTau_,"GenMuFromTau[GenMuNum]/s");
    tree_->Branch("GenMuPt", GenMuPt_,"GenMuPt[GenMuNum]/F");
    tree_->Branch("GenMuEta", GenMuEta_,"GenMuEta[GenMuNum]/F");
    tree_->Branch("GenMuPhi", GenMuPhi_,"GenMuPhi[GenMuNum]/F");
    tree_->Branch("GenElecNum",&GenElecNum_,"GenElecNum/s");
    tree_->Branch("GenElecFromTau",GenElecFromTau_,"GenElecFromTau[GenElecNum]/s");
    tree_->Branch("GenElecPt", GenElecPt_,"GenElecPt[GenElecNum]/F");
    tree_->Branch("GenElecEta", GenElecEta_,"GenElecEta[GenElecNum]/F");
    tree_->Branch("GenElecPhi", GenElecPhi_,"GenElecPhi[GenElecNum]/F");
    tree_->Branch("GenTauNum",&GenTauNum_,"GenTauNum/s");
    tree_->Branch("GenTauPt", GenTauPt_,"GenTauPt[GenTauNum]/F");
    tree_->Branch("GenTauEta", GenTauEta_,"GenTauEta[GenTauNum]/F");
    tree_->Branch("GenTauPhi", GenTauPhi_,"GenTauPhi[GenTauNum]/F");
    
    tree_->Branch("genHT",&genht_,"genHT/F");
    tree_->Branch("genMHT",&genmht_,"genMHT/F");
    tree_->Branch("genNJets",&gennJets_,"genNJets/s");
    //tree_->Branch("BTags",&genBTags_,"BTags/s");
    tree_->Branch("genJet1Pt",&genjet1Pt_,"genJet1Pt/F");
    tree_->Branch("genJet2Pt",&genjet2Pt_,"genJet2Pt/F");
    tree_->Branch("genJet3Pt",&genjet3Pt_,"genJet3Pt/F");
    tree_->Branch("genJet1Eta",&genjet1Eta_,"genJet1Eta/F");
    tree_->Branch("genJet2Eta",&genjet2Eta_,"genJet2Eta/F");
    tree_->Branch("genJet3Eta",&genjet3Eta_,"genJet3Eta/F");
    tree_->Branch("genDeltaPhi1",&gendeltaPhi1_,"genDeltaPhi1/F");
    tree_->Branch("genDeltaPhi2",&gendeltaPhi2_,"genDeltaPhi2/F");
    tree_->Branch("genDeltaPhi3",&gendeltaPhi3_,"genDeltaPhi3/F");
    }
    if(QCD_ || UseAll_)
    {
    std::string name = genra2JetsTag_.label();
    tree_->Branch((name+"Num").c_str(),&(GenJetNum_),(name+"Num/s").c_str());
    tree_->Branch((name+"Pt").c_str(), GenJetPt_, (name+"Pt["+name+"Num]/F").c_str());
    tree_->Branch((name+"Eta").c_str(),GenJetEta_,(name+"Eta["+name+"Num]/F").c_str());
    tree_->Branch((name+"Phi").c_str(),GenJetPhi_,(name+"Phi["+name+"Num]/F").c_str());
    tree_->Branch((name+"E").c_str(),  GenJetE_,  (name+"E["+name+"Num]/F").c_str());
    }
  }
  for(unsigned int i = 0; i < filterDecisionTags_.size(); ++i) {
    if(debug_) std::cout<<"Filter"<<i<<" with name"<<filterDecisionTags_.at(i).label()<<" has been selected"<<std::endl;
    TString name = "Filter_";
    name += filterDecisionTags_.at(i).label();
    tree_->Branch(name,&(filterDecisions_.at(i)),name+"/b");
  }
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RA2TreeMaker::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
 * void
 * RA2TreeMaker::beginRun(edm::Run const&, edm::EventSetup const&)
 * {
 * }
 */

// ------------ method called when ending the processing of a run  ------------
/*
 * void
 * RA2TreeMaker::endRun(edm::Run const&, edm::EventSetup const&)
 * {
 * }
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
 * void
 * RA2TreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 * {
 * }
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
 * void
 * RA2TreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 * {
 * }
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RA2TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void 
RA2TreeMaker::setBranchVariablesToDefault() 
{
  // event information
  runNum_=0;
  lumiBlockNum_=0;
  evtNum_=0;
  nVtx_=0;
  for(unsigned int i = 0; i < varsDouble_.size(); ++i) {
    varsDouble_.at(i) = 9999.;
  }
  // main search variables
  ht_=0.;
  mht_=0.;
  mhtEta_=0.;
  mhtPhi_=0.;
  metPt_=0.;
  metEta_=0.;
  metPhi_=0.;
  nJets_=0;
  BTags_=0;
  BTags_=0;
  nIsoLeptons_=0;
  nIsoTracks_=0;
  jet1Pt_=-10;
  jet2Pt_=-10;
  jet3Pt_=-10;
  jet1Eta_=-10;
  jet2Eta_=-10;
  jet3Eta_=-10;
  deltaPhi1_=-10;
  deltaPhi2_=-10;
  deltaPhi3_=-10;
  // ra2Jets
  for(unsigned int i = 0; i < ra2JetsCollectionInputTag_.size(); ++i) {
    ra2JetsN_.at(i) = 0;
    ra2Jetsht_.at(i)=0;
    ra2Jetsmht_.at(i)=0;
    ra2JetsnJets_.at(i)=0;
    ra2JetsBTags_.at(i)=0;
    ra2Jetsjet1Pt_.at(i)=-9999.;
    ra2Jetsjet2Pt_.at(i)=-9999.;
    ra2Jetsjet3Pt_.at(i)=-9999.;
    ra2Jetsjet1Eta_.at(i)=-9999.;
    ra2Jetsjet2Eta_.at(i)=-9999.;
    ra2Jetsjet3Eta_.at(i)=-9999.;
    ra2JetsdeltaPhi1_.at(i)=-1.;
    ra2JetsdeltaPhi2_.at(i)=-1.;
    ra2JetsdeltaPhi3_.at(i)=-1.;
    for(unsigned int j = 0; j < nMaxCandidates_; ++j) {
      
      ra2JetsPt_.at(i)[j]  = -9999.;
      ra2JetsEta_.at(i)[j] = -9999.;
      ra2JetsPhi_.at(i)[j] = -9999.;
      ra2JetsE_.at(i)[j]   = -9999.;
      ra2JetsBTagValue_.at(i)[j]   = -9999.;
      ra2JetsBTag_.at(i)[j]   = 1000.;
    }
  }
  for(unsigned int i = 0; i < filterDecisions_.size(); ++i) {
    filterDecisions_.at(i) = 0;
  }
  for(unsigned int i = 0; i < varsDouble_.size(); ++i) {
    varsDouble_.at(i) = 9999.;
  }
  for(unsigned int i = 0; i < leptonTag_.size(); ++i) {
    leptonN_[i] = 0;
    for(unsigned int j = 0; j < 25; ++j) {
      leptonPt_.at(i)[j]  = -9999.;
      leptonEta_.at(i)[j] = -9999.;
      leptonPhi_.at(i)[j] = -9999.;
      leptonE_.at(i)[j]   = -9999.;
    }
  }
  for(unsigned int i = 0; i < IsoTrackTag_.size(); ++i) {
    isoTrackN_[i] = 0;
    for(unsigned int j = 0; j < 25; ++j) {
      isoTrackPt_.at(i)[j]  = -9999.;
      isoTrackEta_.at(i)[j] = -9999.;
      isoTrackPhi_.at(i)[j] = -9999.;
      isoTrackE_.at(i)[j]   = -9999.;
    }
  }
  if(MC_)
  {
    GenWNum_=0;
    GenMuNum_=0;
    GenElecNum_=0;
    GenTauNum_=0;
    for(unsigned int j = 0; j < 25; ++j) {
      GenWPt_[j]  = -9999.;
      GenWPhi_[j]  = -9999.;
      GenWEta_[j]  = -9999.;
      GenMuFromTau_[j]  = 1000;
      GenMuPt_[j]  = -9999.;
      GenMuPhi_[j]  = -9999.;
      GenMuEta_[j]  = -9999.;
      GenElecFromTau_[j]  = 1000;
      GenElecPt_[j]  = -9999.;
      GenElecPhi_[j]  = -9999.;
      GenElecEta_[j]  = -9999.;
      GenTauPt_[j]  = -9999.;
      GenTauPhi_[j]  = -9999.;
      GenTauEta_[j]  = -9999.;
      
    }
    genmht_=0;
    genht_=0;
    gennJets_=0;
    GenJetNum_=0;
    for(unsigned int j =0; j<200;j++)
    {
      GenJetPt_[j]=-9999;
      GenJetEta_[j]=-9999;
      GenJetPhi_[j]=-9999;
      GenJetE_[j]=-9999;
    }
  }
  
}
const reco::Candidate* RA2TreeMaker::promtParticle(const reco::Candidate * particle)
{
  int particleId = particle->pdgId();
  std::cout<<"promtParticle::Particle Id:"<<particleId<<std::endl;
  for(size_t i=0;i< particle->numberOfMothers();i++)
  {
    if(particle->mother(i)->pdgId()==particleId) return promtParticle( particle->mother(i) );
    if(abs(particle->mother(i)->pdgId())==24 || abs(particle->mother(i)->pdgId())==15 ) return particle->mother(i);
  }
  return 0;
  
}
bool RA2TreeMaker::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;
  
  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++)
  {
    if(isAncestor(ancestor,particle->mother(i))) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}
const reco::Candidate* RA2TreeMaker::TauFound(const reco::Candidate * particle)
{
  for(size_t i=0;i< particle->numberOfDaughters();i++)
  {
    if(abs(particle->daughter(i)->pdgId() )== 24 || abs(particle->daughter(i)->pdgId() )== 15) return TauFound(particle->daughter(i));
  }
  return particle;
  
}
const reco::Candidate* RA2TreeMaker::WFound(const reco::Candidate * particle)
{
  for(size_t i=0;i< particle->numberOfDaughters();i++)
  {
    if(abs(particle->daughter(i)->pdgId() )== 24) return WFound(particle->daughter(i));
  }
  return particle;
  
}
//define this as a plug-in
DEFINE_FWK_MODULE(RA2TreeMaker);
