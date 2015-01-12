# $Id: makeTreeFromPAT_cff.py,v 1.16 2013/01/24 15:42:53 mschrode Exp $
#

import FWCore.ParameterSet.Config as cms

def makeTreeTreeFromMiniADO(process,
outFileName,
NJetsMin=2,
HTMin=350.,
MHTMin=0.,
reportEveryEvt=10,
testFileName="",
Global_Tag="",
MC=False,
debug = False,
QCD=False,
LostLepton=False,
numProcessedEvt=1000,
skip=0
):

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = Global_Tag

    ## --- Log output ------------------------------------------------------
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.MessageLogger.cerr = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
        )
    process.MessageLogger.cout = cms.untracked.PSet(
        INFO = cms.untracked.PSet(reportEvery = cms.untracked.int32(reportEveryEvt))
        )
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        ) 


    ## --- Files to process ------------------------------------------------
    process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(numProcessedEvt)
        )
    process.source = cms.Source("PoolSource",

        fileNames = cms.untracked.vstring(
'/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/4836AA58-5A6B-E411-8865-20CF305B053E.root',
'/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/829D372D-7F6B-E411-81B1-0025907B5048.root',
'/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/C43D68C5-8D6B-E411-B030-0025907750A0.root'

),
	skipEvents=cms.untracked.uint32(skip),
 #       fileNames = cms.untracked.vstring(
 #		'file:/nfs/dust/cms/user/csander/LHE/workdir/simulation_test/T1qqqqHV/output_66.root'
#		)
        )
        
    #hltPath=['HLT_PFHT350_PFMET100_v*','HLT_PFNoPUHT350_PFMET100_v*']
    #process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
    #process.hltHighLevel.HLTPaths = cms.vstring(hltPath)
    #process.hltHighLevel.andOr = cms.bool(True)
    #process.hltHighLevel.throw = cms.bool(False)

    #process.HLTSelection = cms.Sequence(
     #   process.hltHighLevel
      #  )
   # if MC:
    #    print "Running on MC: removing HLT selection"
     #   process.HLTSelection.remove(process.hltHighLevel)
    #elif not hltPath:
     #   print "Empty list of HLT paths: removing HLT selection"
      #  process.HLTSelection.remove(process.hltHighLevel)
    ## --- Output file -----------------------------------------------------
    process.TFileService = cms.Service(
        "TFileService",
        fileName = cms.string(outFileName+".root")
        )
	    
    ## --- Selection sequences ---------------------------------------------
    # leptons
    process.load("PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi")
    process.load("PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi")
    process.selectedIDIsoMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
    (pfIsolationR04().sumChargedHadronPt+
    max(0.,pfIsolationR04().sumNeutralHadronEt+
    pfIsolationR04().sumPhotonEt-
    0.50*pfIsolationR04().sumPUPt))/pt < 0.20 &&
    (isPFMuon && (isGlobalMuon || isTrackerMuon))'''))
    process.selectedIDMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
    (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))
    process.selectedIDIsoElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
    gsfTrack.isAvailable() &&
    gsfTrack.hitPattern().numberOfLostHits('MISSING_INNER_HITS')<2 &&
    (pfIsolationVariables().sumChargedHadronPt+
    max(0.,pfIsolationVariables().sumNeutralHadronEt+
    pfIsolationVariables().sumPhotonEt-
    0.5*pfIsolationVariables().sumPUPt))/pt < 0.20'''))
    process.selectedIDElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
    gsfTrack.isAvailable() &&
    gsfTrack.hitPattern().numberOfLostHits('MISSING_INNER_HITS')<2'''))
    
    
       ## --- Setup of TreeMaker ----------------------------------------------
    FilterNames = cms.VInputTag()
 #   FilterNames.append(cms.InputTag("HBHENoiseFilterRA2","HBHENoiseFilterResult","PAT"))
 #   FilterNames.append(cms.InputTag("beamHaloFilter"))
 #   FilterNames.append(cms.InputTag("eeNoiseFilter"))
 #   FilterNames.append(cms.InputTag("trackingFailureFilter"))
  #  FilterNames.append(cms.InputTag("inconsistentMuons"))
 #   FilterNames.append(cms.InputTag("greedyMuons"))
 #   FilterNames.append(cms.InputTag("ra2EcalTPFilter"))
#    FilterNames.append(cms.InputTag("ra2EcalBEFilter"))
#    FilterNames.append(cms.InputTag("hcalLaserEventFilter"))
#    FilterNames.append(cms.InputTag("ecalLaserCorrFilter"))
#    FilterNames.append(cms.InputTag("eeBadScFilter"))
#    FilterNames.append(cms.InputTag("PBNRFilter"))
#    FilterNames.append(cms.InputTag("HCALLaserEvtFilterList2012"))
#    FilterNames.append(cms.InputTag("manystripclus53X"))
#    FilterNames.append(cms.InputTag("toomanystripclus53X"))
#    FilterNames.append(cms.InputTag("logErrorTooManyClusters"))
#    FilterNames.append(cms.InputTag("RA2HONoiseFilter"))
    
    
    ## --- Setup WeightProducer -------------------------------------------
    from AllHadronicSUSY.WeightProducer.getWeightProducer_cff import getWeightProducer
    process.WeightProducer = getWeightProducer(testFileName)
    process.WeightProducer.Lumi                       = cms.double(5000)
    process.WeightProducer.PU                         = cms.int32(0) # PU S10 3 for S10 2 for S7
    process.WeightProducer.FileNamePUDataDistribution = cms.string("NONE")
    print process.WeightProducer.PU

    from RecoBTag.Configuration.RecoBTag_cff import *
    from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import *
    process.slimmedJetsPFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
      j2tParametersVX,
      jets = cms.InputTag("iterativeCone5PFJets")
    )
    process.slimmedJetsPFJetTracksAssociatorAtVertex.jets = "slimmedJets"
    process.slimmedJetsPFJetTracksAssociatorAtVertex.tracks = "generalTracks"
    
    process.slimmedJetsPFImpactParameterTagInfos = impactParameterTagInfos.clone()
    process.slimmedJetsPFImpactParameterTagInfos.jetTracks = "slimmedJetsPFJetTracksAssociatorAtVertex"
    process.slimmedJetsPFSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
    process.slimmedJetsPFSecondaryVertexTagInfos.trackIPTagInfos = "slimmedJetsPFImpactParameterTagInfos"
    #slimmedJetsPFSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
    #slimmedJetsPFSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("slimmedJetsPFSecondaryVertexTagInfos") )
    process.slimmedJetsPFCombinedSecondaryVertexBJetTags = combinedInclusiveSecondaryVertexV2BJetTags.clone()
    process.slimmedJetsPFStandardCombinedSecondaryVertex = combinedSecondaryVertex.clone()
    process.slimmedJetsPFCombinedSecondaryVertexBJetTags.jetTagComputer = cms.string('slimmedJetsPFStandardCombinedSecondaryVertex')
    process. slimmedJetsPFCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("slimmedJetsPFImpactParameterTagInfos"), cms.InputTag("slimmedJetsPFSecondaryVertexTagInfos") )
    
    process.slimmedJetsPFJetBtaggingSV = cms.Sequence(
    	process.slimmedJetsPFImpactParameterTagInfos *
    process.slimmedJetsPFSecondaryVertexTagInfos *
    # slimmedJetsPFStandardCombinedSecondaryVertex *
    process.slimmedJetsPFCombinedSecondaryVertexBJetTags
    )
    process.slimmedJetsPFJetsBtag = cms.Sequence(
    process.slimmedJetsPFJetTracksAssociatorAtVertex *
    process.slimmedJetsPFJetBtaggingSV
    )
    
    ## isotrack producer
    from AllHadronicSUSY.Utils.trackIsolationMaker_cfi import trackIsolationFilter
    from AllHadronicSUSY.Utils.trackIsolationMaker_cfi import trackIsolationCounter
    ## default
    process.IsolatedTracks = trackIsolationFilter.clone(
      doTrkIsoVeto= False,
      vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidatesTag = cms.InputTag("packedPFCandidates"),
      dR_ConeSize         = cms.double(0.3),
      dz_CutValue         = cms.double(0.05),
      minPt_PFCandidate   = cms.double(15.0),
      isoCut              = cms.double(0.1),
      )
    #study
    process.IsolatedTracksPT10 = trackIsolationFilter.clone(
      doTrkIsoVeto= False,
      vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidatesTag = cms.InputTag("packedPFCandidates"),
      dR_ConeSize         = cms.double(0.3),
      dz_CutValue         = cms.double(0.05),
      minPt_PFCandidate   = cms.double(10.0),
      isoCut              = cms.double(0.1),
      )
    process.IsolatedTracksPT10IsoCut08 = trackIsolationFilter.clone(
      doTrkIsoVeto= False,
      vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidatesTag = cms.InputTag("packedPFCandidates"),
      dR_ConeSize         = cms.double(0.3),
      dz_CutValue         = cms.double(0.05),
      minPt_PFCandidate   = cms.double(10.0),
      isoCut              = cms.double(0.08),
      )
    process.IsolatedTracksPT10IsoCut12 = trackIsolationFilter.clone(
      doTrkIsoVeto= False,
      vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidatesTag = cms.InputTag("packedPFCandidates"),
      dR_ConeSize         = cms.double(0.3),
      dz_CutValue         = cms.double(0.05),
      minPt_PFCandidate   = cms.double(10.0),
      isoCut              = cms.double(0.12),
      )
    process.CountIsoTracks = trackIsolationCounter.clone(
      src = cms.InputTag("IsolatedTracks"),
      minNumber = 1,
      )

    # Producers
    from AllHadronicSUSY.Utils.subJetSelection_cfi import SubJetSelection
    process.HTJets = SubJetSelection.clone(
    JetTag  = cms.InputTag('slimmedJets'),
    MinPt								  = cms.double(0),
    MaxEta								  = cms.double(2.5),
    )
    from AllHadronicSUSY.Utils.htdouble_cfi import htdouble
    process.HT = htdouble.clone(
    JetTag  = cms.InputTag('HTJets'),
    )
    from AllHadronicSUSY.Utils.njetint_cfi import njetint
    process.NJets = njetint.clone(
    JetTag  = cms.InputTag('HTJets'),
    )
    from AllHadronicSUSY.Utils.btagint_cfi import btagint
    process.BTags = btagint.clone(
    JetTag  = cms.InputTag('HTJets'),
    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    BTagCutValue					= cms.double(0.679)
    )
    from AllHadronicSUSY.Utils.subJetSelection_cfi import SubJetSelection
    process.MHTJets = SubJetSelection.clone(
    JetTag  = cms.InputTag('slimmedJets'),
    MinPt								  = cms.double(30),
    MaxEta								  = cms.double(5.0),
    )
    from AllHadronicSUSY.Utils.jetproperties_cfi import jetproperties
    process.MHTJetsProperties = jetproperties.clone(
    JetTag  = cms.InputTag('MHTJets'),
    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    )
    from AllHadronicSUSY.Utils.mhtdouble_cfi import mhtdouble
    process.MHT = mhtdouble.clone(
    JetTag  = cms.InputTag('MHTJets'),
    )
    from AllHadronicSUSY.Utils.deltaphidouble_cfi import deltaphidouble
    process.DeltaPhi = deltaphidouble.clone(
    DeltaPhiJets  = cms.InputTag('HTJets'),
    MHTJets  = cms.InputTag("MHTJets"),
    )
    from AllHadronicSUSY.Utils.metdouble_cfi import metdouble
    process.MET = metdouble.clone(
    METTag  = cms.InputTag("slimmedMETs"),
    )
    from AllHadronicSUSY.Utils.leptonint_cfi import leptonint
    process.Leptons = leptonint.clone(
    LeptonTag = cms.VInputTag(cms.InputTag('selectedIDIsoMuons'),cms.InputTag('selectedIDIsoElectrons')),
    srcEle = cms.InputTag("slimmedElectrons"),
    srcMuon = cms.InputTag("slimmedMuons"),
    )
    from AllHadronicSUSY.Utils.primaryverticies_cfi import primaryverticies
    process.NVtx = primaryverticies.clone(
    VertexCollection  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    )
    from AllHadronicSUSY.Utils.genLeptonRecoCand_cfi import genLeptonRecoCand
    process.GenLeptons = genLeptonRecoCand.clone(
    PrunedGenParticleTag  = cms.InputTag("prunedGenParticles"),
    )
    RecoCandVector = cms.vstring()
    #RecoCandVector.extend(['selectedIDIsoMuons','selectedIDIsoElectrons','IsolatedTracks']) # basic muons electrons and isoalted tracks
    #RecoCandVector.extend(['selectedIDMuons','selectedIDElectrons']) # mu and e no isolation cuts
    #RecoCandVector.extend(['GenLeptons:Boson(GenBoson)|GenLeptons:BosonPDGId(I_GenBosonPDGId)','GenLeptons:Muon(GenMu)|GenLeptons:MuonTauDecay(I_GenMuFromTau)' ,'GenLeptons:Electron(GenElec)|GenLeptons:ElectronTauDecay(I_GenElecFromTau)','GenLeptons:Tau(GenTau)|GenLeptons:TauHadronic(I_GenTauHad)'] ) # gen information on leptons
    #RecoCandVector.extend(['MHTJetsProperties(MHTJets)|MHTJetsProperties:bDiscriminator(F_bDiscriminator)|MHTJetsProperties:chargedEmEnergyFraction(F_chargedEmEnergyFraction)|MHTJetsProperties:chargedHadronEnergyFraction(F_chargedHadronEnergyFraction)|MHTJetsProperties:chargedHadronMultiplicity(I_chargedHadronMultiplicity)|MHTJetsProperties:electronMultiplicity(I_electronMultiplicity)|MHTJetsProperties:jetArea(F_jetArea)|MHTJetsProperties:muonEnergyFraction(F_muonEnergyFraction)|MHTJetsProperties:muonMultiplicity(I_muonMultiplicity)|MHTJetsProperties:neutralEmEnergyFraction(F_neutralEmEnergyFraction)|MHTJetsProperties:neutralHadronMultiplicity(I_neutralHadronMultiplicity)|MHTJetsProperties:photonEnergyFraction(F_photonEnergyFraction)|MHTJetsProperties:photonMultiplicity(I_photonMultiplicity)|MHTJetsProperties:CorrFactor(F_CorrFactor)|MHTJetsProperties:genEmE(F_genEmE)|MHTJetsProperties:genHadE(F_genHadE)|MHTJetsProperties:genInvE(F_genInvE)|MHTJetsProperties:genMatch(I_genMatch)|MHTJetsProperties:flavour(I_flavour)|MHTJetsProperties:isGood(I_isGood)'] ) # jet information on various variables

    
    from AllHadronicSUSY.TreeMaker.treeMaker import TreeMaker
    process.TreeMaker2 = TreeMaker.clone(
    	TreeName          = cms.string("ReducedTree"),
    	VarsRecoCand = RecoCandVector,
    	#VarsRecoCand = cms.vstring('selectedIDIsoMuons','selectedIDIsoElectrons','IsolatedTracks','HTJets'),
    	VarsDouble        = cms.VInputTag(cms.InputTag('WeightProducer:weight'),cms.InputTag('MHT:MHT'), cms.InputTag('MHT:genMHT'),cms.InputTag('MET:MET'),cms.InputTag('MET:metEnergy'), cms.InputTag('MET:metPhi'),cms.InputTag('MET:mEtSig'),cms.InputTag("MET:metSumEt"),cms.InputTag("MET:genMET"),cms.InputTag("MET:genmetEt"),cms.InputTag("MET:genmetE"),cms.InputTag("MET:genmetPhi"),cms.InputTag('HT:HT'),cms.InputTag('HT:HT30'), cms.InputTag('HT:HT50'),cms.InputTag('DeltaPhi:DeltaPhi1'),cms.InputTag('DeltaPhi:DeltaPhi2'),cms.InputTag('DeltaPhi:DeltaPhi3'), cms.InputTag('DeltaPhi:DeltaPhiMinMet'),cms.InputTag('DeltaPhi:DeltaPhiMinMHt'),cms.InputTag('MHTJetsProperties:DeltaPhiN1'), cms.InputTag('MHTJetsProperties:DeltaPhiN2'), cms.InputTag('MHTJetsProperties:DeltaPhiN3'), cms.InputTag('MHTJetsProperties:minDeltaPhiN')),
    	VarsDoubleNamesInTree = cms.vstring('WeightProducer','MHT', 'genMHT','MET','metEnergy', 'metPhi', 'mEtSig','metSumEt','genMET','genmetEt', 'genmetE','genmetPhi','HT', 'HT30','HT50','DeltaPhi1','DeltaPhi2','DeltaPhi3','DeltaPhiMinMet', 'DeltaPhiMinMHt', 'DeltaPhiN1', 'DeltaPhiN2','DeltaPhiN3', 'minDeltaPhiN'),
    	VarsInt = cms.VInputTag(cms.InputTag('NJets:NJets'), cms.InputTag('NJets:NJets20'),cms.InputTag('NJets:NJets30'),cms.InputTag('NJets:NJets50'),cms.InputTag('BTags:BTags'), cms.InputTag('BTags:BTags20'), cms.InputTag('BTags:BTags30'),cms.InputTag('BTags:BTags50'),cms.InputTag('Leptons:Leptons'),cms.InputTag('Leptons:Electrons'),cms.InputTag('Leptons:Muons'),cms.InputTag('NVtx'), cms.InputTag('HT:genJetMatch')),
	VarsIntNamesInTree=cms.vstring('NJets', 'NJets20', 'NJets30', 'NJets50', 'BTags', 'BTags20', 'BTags30', 'BTags50','Leptons', 'Electrons', 'Muons','NVtx', 'genJetMatch',)	
    #	VarsDoubleNamesInTree = cms.vstring('WeightProducer'),
    	)

    ## --- Final paths ----------------------------------------------------

    process.dump = cms.EDAnalyzer("EventContentAnalyzer")
    process.WriteTree = cms.Path(
    	process.selectedIDIsoMuons *
    	process.selectedIDMuons *
    	process.selectedIDIsoElectrons *
    	process.selectedIDElectrons *
    	process.WeightProducer *
    	process.IsolatedTracks *
 #   	process.IsolatedTracksPT10 *
 #   	process.IsolatedTracksPT10IsoCut08 *
 #   	process.IsolatedTracksPT10IsoCut12 *
  #  	process.slimmedJetsPFCombinedSecondaryVertexBJetTags *
      process.HTJets *
      process.HT *
      process.NJets *
      process.BTags *
      process.MHTJets *
      process.MHTJetsProperties *
      process.MHT *
      process.Leptons *
      process.MET *
      process.DeltaPhi *
      process.NVtx *
      process.GenLeptons *
    	#process.dump *
 #   	process.CountIsoTracks *
 #   	process.PrintDecay *
    	process.TreeMaker2

        )
