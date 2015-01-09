import FWCore.ParameterSet.Config as cms

RA2TreeMaker = cms.EDProducer(
'RA2TreeMaker',
# Name of the output tree
TreeName          = cms.string('RA2Tree'),

# collection from which the tree variable "NumVtx" is determined
VertexCollection  = cms.InputTag('goodVertices'),
# List of InputTags for double-precision variables (double) stored in
# the event. (For space reason, they are stored as Float_t in the tree.)
VarsDouble        = cms.VInputTag(),
# Names of the double-precision variables as stored in the tree. If
# this vector is not specified, the generic names "<InputTag::label()>"
# are used.
VarsDoubleNamesInTree = cms.vstring(),
# list of filter decisions (bool) written from filters in tag mode
# will be stored as "Filter_..."
Filters           = cms.VInputTag(),

MC = cms.bool(False),
QCD = cms.bool(False),
LostLepton = cms.bool(False),
debug = cms.bool(False),
# CMSSW 7 varaibles miniAOD
prunedGenParticles = cms.InputTag("prunedGenParticles"),
packedGenParticles = cms.InputTag("packedGenParticles"),   
RA2DefaultJetsTag = cms.InputTag("patJetsAK5PFCHS"), 
METTag  = cms.InputTag("slimmedMETs"),
bTagName     = cms.string("combinedSecondaryVertexBJetTags"),
bTagValue    = cms.double(0.679),
# jet selection for HT MHT NJets and Btag
MinHT   = cms.double(0),
MinMHT   = cms.double(0),
#minNJets_   = cms.int(2),
MaxEtaHTJets = cms.double(2.5),
MinPTHTJets = cms.double(50),
MaxEtaMHTJets = cms.double(5),
MinPTMHTJets = cms.double(30),
# ra2 jet collections
ra2JetsCollectionInputTag = cms.VInputTag(),
ra2JetsCollectionNameInTree = cms.vstring(),
ra2JetsBTagInputTag = cms.vstring(),
ra2JetsBTagValueInput_ = cms.vdouble(),
# ra2 leptons
LeptonTag = cms.VInputTag(),
LeptonTagName = cms.vstring(),
# iso tracks
IsoTrackTag = cms.VInputTag(),
IsoTrackTagName = cms.vstring(),
GenJetTag = cms.InputTag("slimmedGenJets"),
)
