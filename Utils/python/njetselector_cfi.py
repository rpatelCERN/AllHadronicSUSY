import FWCore.ParameterSet.Config as cms

nJetProducer = cms.EDProducer('NJetSelector',
  JetCollection = cms.InputTag("patJets"),
  MinJetPt      = cms.double(50),
  MaxJetEta     = cms.double(2.5)
)