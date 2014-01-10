import FWCore.ParameterSet.Config as cms

import sys

startEv=0
if len(sys.argv)>2 : startEv=int(sys.argv[2])
print '[runMyEDMtoHEPMCAnalyzer_cfg.py] starting at event #%d'%startEv

process = cms.Process("Analysis")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring('file:pp_WZtoAnything_13TeV.root'),
                             skipEvents = cms.untracked.uint32(startEv)
                             )
process.HepMCconvert = cms.EDAnalyzer( "MyEDMtoHEPMCAnalyzer" )
process.p = cms.Path( process.HepMCconvert )
