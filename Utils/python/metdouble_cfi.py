import FWCore.ParameterSet.Config as cms

metdouble = cms.EDProducer('METDouble',
METTag  = cms.InputTag("slimmedMETs"),
)
