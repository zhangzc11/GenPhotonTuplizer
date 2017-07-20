import FWCore.ParameterSet.Config as cms
import os, sys, imp, re


CMSSW_VERSION=os.getenv("CMSSW_VERSION")
process = cms.Process("GenPhotonTuplizer")


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'MCRUN2_71_V1::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.MessageLogger.cerr.FwkReport.reportEvery = 100

#load input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
	#'root://cms-xrd-global.cern.ch//store/mc/RunIIWinter15GenOnly/GluGluHToGG_M125_13TeV_powheg2_minloHJJ_pythia8/GEN/MCRUN2_71_V1-v2/110000/04485D22-A46A-E711-B25D-0025904E41E4.root'
	'root://cms-xrd-global.cern.ch//store/mc/RunIIWinter15GenOnly/GluGluHToGG_M125_13TeV_powheg2_minloHJJ_pythia8/GEN/MCRUN2_71_V1-v2/110000/1E9D9DAB-616A-E711-B44A-002590E50AFE.root'
    )
)

#process.Out = cms.OutputModule("PoolOutputModule",
#         fileName = cms.untracked.string ("MyOutputFile.root")
#)

#define output file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("genPhotonNtuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

#########################paratmeters for the tuplizer##############################
process.ntuples = cms.EDAnalyzer('GenPhotonTuplizer',
genParticles = cms.InputTag("genParticles"),
genJets = cms.InputTag("ak4GenJets"),
genInfo = cms.InputTag("generator", "", "GEN")
#lheInfo = cms.InputTag("externalLHEProducer", "", "LHE"),
)

#define path
process.p = cms.Path()
process.p *= process.ntuples

#process.end = cms.EndPath(process.Out)
