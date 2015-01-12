import FWCore.ParameterSet.Config as cms

GenInfo = cms.EDProducer('GenInfo',
  PrunedGenParticleTag  = cms.InputTag("prunedGenParticles"),
)
