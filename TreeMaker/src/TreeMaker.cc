// -*- C++ -*-
//
// Package:    AllHadronicSUSY/TreeMaker
// Class:      TreeMaker
// 
/**\class TreeMaker TreeMaker.cc AllHadronicSUSY/TreeMaker/plugins/TreeMaker.cc
 * 
 * Description: [one line class summary]
 * 
 * Implementation:
 *     [Notes on implementation]
 */
//
// Original Author:  Arne-Rasmus Draeger
//         Created:  Fri, 03 Dec 2014 13:48:35 GMT
//
//
#include "AllHadronicSUSY/TreeMaker/interface/TreeMaker.h"
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
TreeMaker::TreeMaker(const edm::ParameterSet& iConfig)
: nMaxCandidates_(200), tree_(0)
{
  // generell
  treeName_ = iConfig.getParameter<std::string>("TreeName");
	debug_ = iConfig.getParameter<bool> ("debug");
	varsDoubleTags_ = iConfig.getParameter< std::vector<edm::InputTag> >("VarsDouble");
	varsDoubleNames_= iConfig.getParameter< std::vector<std::string> >  ("VarsDoubleNamesInTree");
	varsIntTags_ = iConfig.getParameter< std::vector<edm::InputTag> >("VarsInt");
	varsIntNames_= iConfig.getParameter< std::vector<std::string> >  ("VarsIntNamesInTree");
	varsBoolTags_ = iConfig.getParameter< std::vector<edm::InputTag> >("VarsBool");
	varsBoolNames_= iConfig.getParameter< std::vector<std::string> >  ("VarsBoolNamesInTree");
	varsTLorentzVectorTags_ = iConfig.getParameter< std::vector<edm::InputTag> >("VarsTLorentzVector");
	varsTLorentzVectorNames_= iConfig.getParameter< std::vector<std::string> >  ("VarsTLorentzVectorNamesInTree");
	vectorTLorentzVectorTags_ = iConfig.getParameter< std::vector<edm::InputTag> >("VectorTLorentzVector");
	vectorTLorentzVectorNames_= iConfig.getParameter< std::vector<std::string> >  ("VectorTLorentzVectorNamesInTree");
  // input tags for float variables eg HT MHT MET or what not
  varsRecoCandNames_= iConfig.getParameter< std::vector<std::string> >  ("VarsRecoCand");
}


TreeMaker::~TreeMaker()
{
  
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TreeMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
	//Float variables
	for(unsigned int i = 0; i < varsDoubleTags_.size(); ++i) {
		edm::Handle<double> var;
		iEvent.getByLabel(varsDoubleTags_.at(i),var);
		if( var.isValid() ) {
			varsDouble_.at(i) = *var;
		}
	}
	
	for(unsigned int i = 0; i < varsIntTags_.size(); ++i) {
		edm::Handle<int> var;
		iEvent.getByLabel(varsIntTags_.at(i),var);
		if( var.isValid() ) {
			varsInt_.at(i) = *var;
		}
	}
	for(unsigned int i = 0; i < varsBoolTags_.size(); ++i) {
		edm::Handle<bool> dec;
		iEvent.getByLabel(varsBoolTags_.at(i),dec);
		if( dec.isValid() ) {
			if( *dec ) varsBool_.at(i) = 1;
			else varsBool_.at(i) = 0;
		}
	}
	
	for(unsigned int i = 0; i < varsTLorentzVectorTags_.size(); ++i) {
		edm::Handle<TLorentzVector> var;
		iEvent.getByLabel(varsTLorentzVectorTags_.at(i),var);
		if( var.isValid() ) {
			varsTLorentzVector_.at(i) = *var;
		}
	}
	for(unsigned int i = 0; i < varsTLorentzVectorTags_.size(); ++i) {
		edm::Handle<std::vector<TLorentzVector> > var;
		iEvent.getByLabel(varsTLorentzVectorTags_.at(i),var);
		if( var.isValid() ) 
		{
			for(unsigned int ii=0; ii< var->size();ii++)
			{
				vectorTLorentzVector_.at(i).push_back(var->at(ii));
			}
		}
	}
	for(unsigned int i = 0; i < varsRecoCandTags_.size(); ++i) 
	{
		edm::Handle< edm::View<reco::Candidate> > cands;
		iEvent.getByLabel(varsRecoCandTags_.at(i),cands);
		if( cands.isValid() ) 
		{
			RecoCandN_[i] = (unsigned int)(cands->size());
			for(unsigned int ii=0; ii < cands->size();ii++)
			{
				RecoCandPt_.at(i)[ii] = cands->at(ii).pt();
				RecoCandEta_.at(i)[ii] = cands->at(ii).eta();
				RecoCandPhi_.at(i)[ii] = cands->at(ii).phi();
				RecoCandE_.at(i)[ii] = cands->at(ii).energy();
				RecoCandLorentzVector_.at(i)[ii].SetXYZT(cands->at(ii).p4().X(),cands->at(ii).p4().Y(),cands->at(ii).p4().Z(),cands->at(ii).p4().T());
			}
			for(unsigned int ii=0; ii< RecoCandAdditionalFloatVariablesTags_[i].size();ii++)
			{
				edm::Handle<std::vector<double> > FloatVar;
				iEvent.getByLabel(RecoCandAdditionalFloatVariablesTags_[i].at(ii),FloatVar);
				if( !FloatVar.isValid() )
				{
					if(debug_)std::cout<<"Warning Float variable with lable: "<<RecoCandAdditionalFloatVariablesTags_[i].at(ii).label()<<" not found!!!"<<std::endl;
					break;
				}
				for(unsigned int iii=0; iii<cands->size();iii++)
				{
					RecoCandAdditionalFloatVariables_.at(i)[ii][iii] = FloatVar->at(iii);
				}
			}
			// loop over int variables
			for(unsigned int ii=0; ii< RecoCandAdditionalIntVariablesTags_[i].size();ii++)
			{
				edm::Handle<std::vector<int> > IntVar;
				iEvent.getByLabel(RecoCandAdditionalIntVariablesTags_[i].at(ii),IntVar);
				if( !IntVar.isValid() )
				{
					if(debug_)std::cout<<"Warning Int variable with lable: "<<RecoCandAdditionalIntVariablesTags_[i].at(ii).label()<<" not found!!!"<<std::endl;
					break;
				}
				for(unsigned int iii=0; iii<cands->size();iii++)
				{
					RecoCandAdditionalIntVariables_.at(i)[ii][iii] = IntVar->at(iii);
				}
			}
			// loop over bool variables
			for(unsigned int ii=0; ii< RecoCandAdditionalBoolVariablesTags_[i].size();ii++)
			{
				edm::Handle<std::vector<bool> > boolVar;
				iEvent.getByLabel(RecoCandAdditionalBoolVariablesTags_[i].at(ii),boolVar);
				if( !boolVar.isValid() )
				{
					if(debug_)std::cout<<"Warning bool variable with lable: "<<RecoCandAdditionalBoolVariablesTags_[i].at(ii).label()<<" not found!!!"<<std::endl;
					break;
				}
				for(unsigned int iii=0; iii<cands->size();iii++)
				{
					if( boolVar->at(iii) ) RecoCandAdditionalBoolVariables_.at(i)[ii][iii] = 1;
					else RecoCandAdditionalBoolVariables_.at(i)[ii][iii] = 0;
				}
			}
		}
		else if(debug_)std::cout<<"Warning recoCand with tag: "<<varsRecoCandTags_.at(i).label()<<" not found!"<<std::endl;
		
	}
	
// std::cout<<"Filling tree with:"<<std::endl;
// for(unsigned int i=0; i<RecoCandN_.size();i++)std::cout<<"RecoCand["<<i<<"] Numberofentries: "<<RecoCandN_.at(i)<<std::endl;
  tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
TreeMaker::beginJob()
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
	
	varsDouble_ = std::vector<Float_t>(varsDoubleTags_.size(),1.);
	for(unsigned int i = 0; i < varsDouble_.size(); ++i) {
		std::string  name = varsDoubleTags_.at(i).label();
		if( varsDoubleNames_.size() == varsDoubleTags_.size() ) {
			name = varsDoubleNames_.at(i);
		}
		name.erase (std::remove (name.begin(), name.end(), ':'), name.end());
		TString namet = name;
		tree_->Branch(namet,&(varsDouble_.at(i)),namet+"/F");
	}
	varsInt_ = std::vector<int>(varsIntTags_.size(),1);
	for(unsigned int i = 0; i < varsIntTags_.size(); ++i) {
		std::string  name = varsIntTags_.at(i).label();
		if( varsIntNames_.size() == varsIntTags_.size() ) {
			name = varsIntNames_.at(i);
		}
		name.erase (std::remove (name.begin(), name.end(), ':'), name.end());
		TString namet = name;
		tree_->Branch(namet,&(varsInt_.at(i)),namet+"/I");
	}
	varsBool_ = std::vector<UChar_t>(varsBoolTags_.size(),0);
	for(unsigned int i = 0; i < varsBoolTags_.size(); ++i) {
		std::string name = "";
		name += varsBoolTags_.at(i).label();
		if( varsBoolNames_.size() == varsDoubleTags_.size() ) {
			name = varsBoolNames_.at(i);
		}
		name.erase (std::remove (name.begin(), name.end(), ':'), name.end());
		TString namet = name;
		tree_->Branch(namet,&(varsBool_.at(i)),namet+"/b");
	}
  // loop over input varsFloat string to extract optional names in tree
	RecoCandN_ = std::vector<UShort_t>(varsRecoCandNames_.size(),0);
  for(unsigned int i=0; i<varsRecoCandNames_.size();i++)
  {
    RecoCandPt_.push_back (new Float_t[200]);
    RecoCandEta_.push_back(new Float_t[200]);
    RecoCandPhi_.push_back(new Float_t[200]);
    RecoCandE_.push_back  (new Float_t[200]);
		RecoCandLorentzVector_.push_back(new TLorentzVector[200]);

    std::string temp = varsRecoCandNames_[i];
		std::string nameInTree = "";
		std::string ttemp ="";
		

    std::cout<<"RecoCand Setup: item["<<i<<"] full string: "<<temp<<std::endl;
    if(temp.find('|')<temp.size() ) temp = temp.substr(0,temp.find("|") );
		if(temp.find('(')<temp.size() && temp.find(')')<temp.size() ) 
		{
			nameInTree = temp.substr(temp.find('(')+1, temp.find(')')-temp.find('(')-1);
			temp=temp.substr(0,temp.find('('));
		}
		else nameInTree = temp;
    std::cout<<"RecoCand Tag: "<<temp<<std::endl;
    varsRecoCandTags_.push_back(edm::InputTag(temp));
		std::cout<<"RecoCand stored name in tree: "<<nameInTree<<std::endl;
		temp=nameInTree;
		temp.erase (std::remove (temp.begin(), temp.end(), ':'), temp.end());
    tree_->Branch((temp+"Num").c_str(),&(RecoCandN_.at(i)),(temp+"Num/s").c_str());
    tree_->Branch((temp+"Pt").c_str(), RecoCandPt_.at(i), (temp+"Pt["+temp+"Num]/F").c_str());
    tree_->Branch((temp+"Eta").c_str(),RecoCandEta_.at(i),(temp+"Eta["+temp+"Num]/F").c_str());
    tree_->Branch((temp+"Phi").c_str(),RecoCandPhi_.at(i),(temp+"Phi["+temp+"Num]/F").c_str());
    tree_->Branch((temp+"E").c_str(),  RecoCandE_.at(i),  (temp+"E["+temp+"Num]/F").c_str());
		tree_->Branch((temp+"TLorentzVector").c_str(),  RecoCandLorentzVector_.at(i),  (temp+"TLorentzVector["+temp+"Num]/F").c_str());
		std::string mainNameInTree=temp;
    temp = varsRecoCandNames_[i];
		unsigned int countBool=0;
		unsigned int countInt=0;
		unsigned int countFloat=0;
		std::vector<UChar_t*> vecUChart;
		RecoCandAdditionalBoolVariables_.push_back(vecUChart);
		std::vector<Int_t*> vecInt;
		RecoCandAdditionalIntVariables_.push_back(vecInt);
    std::vector<Float_t*> vecFloat;
		RecoCandAdditionalFloatVariables_.push_back(vecFloat);
    std::vector<edm::InputTag> tagvec;
		RecoCandAdditionalBoolVariablesTags_.push_back(tagvec);
		RecoCandAdditionalIntVariablesTags_.push_back(tagvec);
		RecoCandAdditionalFloatVariablesTags_.push_back(tagvec);
    std::string temp2="";
		std::string tag="";
		int typ=-1;
    while (temp.find('|')<temp.size() ) // loop over the posible additonal variables
    {
			typ=-1;
      temp = temp.substr(temp.find("|")+1 );
      temp2="";
      if(temp.find('|')<temp.size() ) temp2 = temp.substr(0, temp.find('|') );
      else temp2=temp;
			// check for typ definition and for optional naming
			if(temp2.find('(')<temp2.size() && temp2.find(')')<temp2.size() )
			{
			//	std::cout<<"POINT1::temp2: "<<temp2<<std::endl;
				// check for optional naming
				if(temp2.find('_')<temp2.size())
				{
					nameInTree = temp2.substr(temp2.find('_')+1, temp2.find(')')-temp2.find('_')-1);
					tag = temp2.substr(0, temp2.find('('));
					if(temp2.find("b_")<temp2.size() )typ = 0;
					if(temp2.find("I_")<temp2.size() )typ = 1;
					if(temp2.find("F_")<temp2.size() )typ = 2;
				}
				else 
				{
					nameInTree = temp2.substr(0, temp2.find('('));
					tag=temp2.substr(0, temp2.find('('));
					if(temp2.find('b')<temp2.size() )typ = 0;
					if(temp2.find('I')<temp2.size() )typ = 1;
					if(temp2.find('F')<temp2.size() )typ = 2;

				}
			}
			else if(typ==-1)std::cout<<"Warning no typ selected for additonal input object: "<<temp2<<" of main varialbe: "<<nameInTree<<"Please use: tag(x_Name) with x=b,I,F (bool, int float) and optional Name for naming in the tree"<<std::endl;
     // std::cout<<"| TagName: "<<temp2<<std::endl;
			nameInTree.erase (std::remove (nameInTree.begin(), nameInTree.end(), ':'), nameInTree.end());
			std::cout<<"Sub Typ: Tag: "<<tag<<", typ: "<<typ<< ", nameIn Tree: "<<nameInTree<<std::endl;
			if(typ==0)
			{
				
				RecoCandAdditionalBoolVariablesTags_[i].push_back(edm::InputTag(tag ) );
				RecoCandAdditionalBoolVariables_[i].push_back(new UChar_t[200]);
				tree_->Branch((mainNameInTree+"_"+nameInTree).c_str(), RecoCandAdditionalBoolVariables_.at(i).at(countBool), (mainNameInTree+"_"+nameInTree+"["+mainNameInTree+"Num]/b").c_str());
// 				tree_->Branch((nameInTree).c_str(), RecoCandAdditionalBoolVariables_.at(i).at(countBool), (nameInTree+"["+mainNameInTree+"Num]/b").c_str());
				countBool++;
			}
			if(typ==1)
			{
				
				RecoCandAdditionalIntVariablesTags_[i].push_back(edm::InputTag(tag ) );
				RecoCandAdditionalIntVariables_[i].push_back(new Int_t[200]);
 				tree_->Branch((mainNameInTree+"_"+nameInTree).c_str(), RecoCandAdditionalIntVariables_.at(i).at(countInt), (mainNameInTree+"_"+nameInTree+"["+mainNameInTree+"Num]/I").c_str());
// 				tree_->Branch((nameInTree).c_str(), RecoCandAdditionalIntVariables_.at(i).at(countInt), (nameInTree+"["+mainNameInTree+"Num]/I").c_str());
				countInt++;
			}
			if(typ==2)
			{
				
			RecoCandAdditionalFloatVariablesTags_[i].push_back(edm::InputTag(tag ) );
			RecoCandAdditionalFloatVariables_[i].push_back(new Float_t[200]);
			tree_->Branch((mainNameInTree+"_"+nameInTree).c_str(), RecoCandAdditionalFloatVariables_.at(i).at(countFloat), (mainNameInTree+"_"+nameInTree+"["+mainNameInTree+"Num]/F").c_str());
      countFloat++;
			}
			if(typ>2)std::cout<<"Error typ is: "<<typ<<" which is not defined!!! check!"<<std::endl;
    }
    RecoCandAdditionalBoolVariablesN_.push_back(countBool);
		RecoCandAdditionalIntVariablesN_.push_back(countInt);
		RecoCandAdditionalFloatVariablesN_.push_back(countFloat);
  }
 	varsTLorentzVector_ = std::vector<TLorentzVector>(varsTLorentzVectorTags_.size(),TLorentzVector(0.,0.,0.,0.));
 	for(unsigned int i = 0; i < varsTLorentzVectorTags_.size(); ++i) {
		std::string  name = "TLorentzVector_";
		name += varsTLorentzVectorTags_.at(i).label();
		if( varsTLorentzVectorTags_.size() == varsTLorentzVectorNames_.size() ) {
			name = "TLorentzVector_"+varsTLorentzVectorNames_.at(i);
		}
		name.erase (std::remove (name.begin(), name.end(), ':'), name.end());
		TString namet = name;
		tree_->Branch(namet,namet,&(varsTLorentzVector_.at(i)));
	}
	for(unsigned int i=0; i< vectorTLorentzVectorTags_.size();i++)
	{
		std::vector<TLorentzVector> vector;
		vectorTLorentzVector_.push_back(vector);
		std::string  name = "VectorTLorentzVector_";
		name += vectorTLorentzVectorTags_.at(i).label();
		if(vectorTLorentzVectorNames_.size() == vectorTLorentzVectorTags_.size())
		{
			name = "VectorTLorentzVector_"+vectorTLorentzVectorNames_.at(i);
		}
		name.erase (std::remove (name.begin(), name.end(), ':'), name.end());
		TString namet = name;
		tree_->Branch(namet,namet,&(vectorTLorentzVector_.at(i)));
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeMaker::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void 
TreeMaker::setBranchVariablesToDefault() 
{
	// event information
	runNum_=0;
	lumiBlockNum_=0;
	evtNum_=0;
  for(unsigned int i = 0; i < varsRecoCandTags_.size(); ++i) 
	{
		RecoCandN_.at(i)=0;
		for(unsigned int ii=0; ii< nMaxCandidates_; ii++)
		{
			RecoCandPt_.at(i)[ii]=-999.;
			RecoCandEta_.at(i)[ii]=-999.;
			RecoCandPhi_.at(i)[ii]=-999.;
			RecoCandE_.at(i)[ii]=-999.;
			RecoCandLorentzVector_.at(i)[ii]=TLorentzVector (0.,0.,0.,0.);
		}
		// loop over the additonal variables
		for(unsigned int iii=0; iii< nMaxCandidates_; iii++)
		{
			for(unsigned int ii=0; ii < RecoCandAdditionalBoolVariablesTags_.at(i).size();ii++)
			{
				RecoCandAdditionalBoolVariables_.at(i)[ii][iii]=0;
			}
			for(unsigned int ii=0; ii < RecoCandAdditionalIntVariablesTags_.at(i).size();ii++)
			{
				RecoCandAdditionalIntVariables_.at(i)[ii][iii]=-9999;
			}
			for(unsigned int ii=0; ii < RecoCandAdditionalFloatVariablesTags_.at(i).size();ii++)
			{
				RecoCandAdditionalFloatVariables_.at(i)[ii][iii]=-9999.;
			}
		}
	}
	for(unsigned int i = 0; i < varsBool_.size(); ++i) {
		varsBool_.at(i) = 0;
	}
	for(unsigned int i = 0; i < varsInt_.size(); ++i) {
		varsInt_.at(i) = 9999;
	}
	for(unsigned int i = 0; i < varsDouble_.size(); ++i) {
		varsDouble_.at(i) = 9999.;
	}
	for(unsigned int i=0; i < vectorTLorentzVectorTags_.size();i++)
	{
		vectorTLorentzVector_.at(i).clear();
	}

  
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeMaker);
