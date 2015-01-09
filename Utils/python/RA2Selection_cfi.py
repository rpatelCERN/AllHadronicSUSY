
import FWCore.ParameterSet.Config as cms

RA2Selection = cms.EDFilter(
    'RA2Selection',
    	HTJets          = cms.InputTag('HTJets'),
	MhtJetTag	= cms.InputTag('MHTJets'),
	MhtTag		= cms.InputTag('mhtPFchs'),
	HTMin		= cms.double(500),
	MHTMin		= cms.double(200),
	deltaPhi1	= cms.double(0.5),
	deltaPhi2	= cms.double(0.5),
	deltaPhi3	= cms.double(0.3),
	nJets		= cms.uint32 (2)
)