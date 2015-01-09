import FWCore.ParameterSet.Config as cms

deltaPhiProducer = cms.EDProducer('DeltaPhiProducer',
	MHT               = cms.InputTag('mhtPF'),
	MHTJets            = cms.InputTag('MHTJets'),  
	MinJetPt      = cms.double(30),
	MaxJetEta     = cms.double(5)
)
