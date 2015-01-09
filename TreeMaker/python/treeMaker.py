import FWCore.ParameterSet.Config as cms

TreeMaker = cms.EDProducer(
'TreeMaker',
# Name of the output tree
TreeName          = cms.string('RA2Tree'),
## might help if something isn working wilil produce couts
debug = cms.bool(False),
# List of InputTags for Float_t variables (Float_t) stored in
# optional names to store in the tree can be defined  not only will have an effect if the number of input tags for variable is exactly the same as number of optional names!!
VarsDouble        = cms.VInputTag(),
VarsDoubleNamesInTree = cms.vstring(),
# List of InputTags for Int variables (Int) stored in
# optional names to store in the tree can be defined  not only will have an effect if the number of input tags for variable is exactly the same as number of optional names!!
VarsInt        = cms.VInputTag(),
VarsIntNamesInTree = cms.vstring(),
# List of InputTags for bool variables (bool) stored in
# optional names to store in the tree can be defined  not only will have an effect if the number of input tags for variable is exactly the same as number of optional names!!
VarsBool        = cms.VInputTag(),
VarsBoolNamesInTree = cms.vstring(),

# List of InputTags for TLorentz variables (TLorentz) stored in
# optional names to store in the tree can be defined  not only will have an effect if the number of input tags for variable is exactly the same as number of optional names!!
VarsTLorentzVector        = cms.VInputTag(),
VarsTLorentzVectorNamesInTree = cms.vstring(),

# List of InputTags for vectors of TLorentz variables (std::vec<TLorentz>) stored in
# optional names to store in the tree can be defined  not only will have an effect if the number of input tags for variable is exactly the same as number of optional names!!
VectorTLorentzVector        = cms.VInputTag(),
VectorTLorentzVectorNamesInTree = cms.vstring(),
# list of reco candidate objects eg leptons: For each reco cand collection the pt eta phi e and Tlorentzvector will be stored in arrays. In addition optional float int or bool variables can be stored
# syntax example:  VarsRecoCand = cms.vstring('selectedIDIsoMuons(selectedIDIsoMuonsName)|Isolation(F_NameOfIsolationInTreeisolation)|Jup(b_nameofJup)|wwww(F)','selectedIDMuons','selectedRecoIsoElec','selectedRecoElec'),
# selectedIDIsoMuons would be the tag of the reco candiate 
# (selectedIDIsoMuonsName ) optoinal naming in the tree if not defined tag name will be used
# separated with | addition variables int float or bool defined by: TagName(typ) always needed optional TagName(typ_NameInTree) 
VarsRecoCand = cms.vstring(),
)
