#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-s', '--short-queue',    dest='squeue'             , help='short batch queue'            , default='1nh')
parser.add_option('-q', '--long-queue' ,    dest='lqueue'             , help='long batch queue'             , default='1nd')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='V00-00-00')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,      type=int)
parser.add_option('-R', '--nRuns'       ,    dest='nRuns'              , help='number of runs'               , default=0,      type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=3,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=3,      type=int)
parser.add_option('-a', '--alpha'       ,    dest='alpha'              , help='incidence angle in rad'       , default=0,      type=float)
parser.add_option('-p', '--phi'         ,    dest='phi'                , help='incidence phi angle in pi unit' , default=0.5,      type=float)
parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0,      type=float)
parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
parser.add_option('-f', '--datafile'    ,    dest='datafile'           , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to process' , default=0,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-E', '--eosin'       ,    dest='eosin'              , help='eos path to read input root file from EOS',  default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

redofit=1
label='200u'

workdir='/afs/cern.ch/work/a/amagnan/PFCalEEAna/'

#nPuVtxset=[0]#,140]
nPuVtxset=[0]
#[0,140,200]

#runlist=[0,1,2,3,4,5,6,7,8,9]
#runlist=[0,1,3,6,8,9]
#runlist=[0,1,5,6,7]
interCalibList=[3] #0,1,2,3,4,5,10,15,20,50]

for nPuVtx in nPuVtxset :
    for interCalib in interCalibList:

        if nPuVtx>0 :
            suffix='Pu%d_IC%d'%(nPuVtx,interCalib)
        else :
            suffix='IC%d'%(interCalib)

        nevents=opt.nevts
        myqueue=opt.lqueue
        bval="BOFF"
        if opt.Bfield>0 : bval="BON" 
    
        outDir='%s/git%s/version%d/model%d/%s/pu%s'%(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,nPuVtx)
        #outDir='%s/git%s/version%d/%s/run%s'%(opt.out,opt.gittag,opt.version,opt.datatype,nPuVtx,run)
        eosDir='%s/git%s/%s'%(opt.eos,opt.gittag,opt.datatype)
        eosDirIn='%s/git%s/%s'%(opt.eosin,opt.gittag,opt.datatype)

        outlog='hreso'
        g4log='hresojob.log'
        os.system('mkdir -p %s/%s'%(workdir,outDir))
    #clean up old batch outputs
        os.system('rm -f %s/%s/*.*.out'%(workdir,outDir))
    #wrapper
        scriptFile = open('%s/%s/runHResoJob.sh'%(workdir,outDir), 'w')
        scriptFile.write('#!/bin/bash\n')
        scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
    #scriptFile.write('cd %s\n'%(outDir))
        outTag='_version%d_model%d_%s'%(opt.version,opt.model,bval)
        #outTag='%s_run%d'%(outTag,run)
        if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)

        scriptFile.write('localdir=`pwd`\n')
        scriptFile.write('cp -r %s/data .\n'%os.getcwd())
        scriptFile.write('cp -r %s/scripts .\n'%os.getcwd())
        scriptFile.write('mkdir -p %s\n'%outDir)
        if (opt.nRuns==0) :
            scriptFile.write('%s/bin/higgsResoWithTruth -c scripts/DefaultConfigHiggs.cfg -n %d -i root://eoscms//eos/cms%s --digifilePath=root://eoscms//eos/cms%s  -s HGcal_%s.root -r Digi%s_%s%s.root -o %s.root --redoStep=%s | tee %s.log\n'%(os.getcwd(),opt.nevts,eosDirIn,eosDir,outTag,suffix,label,outTag,outDir,redofit,outlog))
        else:
            scriptFile.write('%s/bin/higgsResoWithTruth -c scripts/DefaultConfigHiggs.cfg -n %d --nRuns=%d -i root://eoscms//eos/cms%s --digifilePath=root://eoscms//eos/cms%s  -s HGcal_%s -r Digi%s_%s%s -o %s.root --redoStep=%s | tee %s.log\n'%(os.getcwd(),opt.nevts,opt.nRuns,eosDirIn,eosDir,outTag,suffix,label,outTag,outDir,redofit,outlog))

        scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
        scriptFile.write('ls * >> %s\n'%(g4log))
        scriptFile.write('echo "--deleting core files: too heavy!!" >> %s\n'%(g4log))
        scriptFile.write('rm -f core.* >> %s\n'%(g4log))
        scriptFile.write('cp * %s/%s/\n'%(workdir,outDir))
        scriptFile.write('cp %s/* %s/%s/\n'%(outDir,workdir,outDir))
        scriptFile.write('cp %s.root %s/%s.root\n'%(outDir,workdir,outDir))
        scriptFile.write('echo "All done"\n')
        scriptFile.close()
    
    #submit
        os.system('chmod u+rwx %s/%s/runHResoJob.sh'%(workdir,outDir))
        if opt.nosubmit : os.system('echo bsub -q %s %s/%s/runHResoJob.sh'%(myqueue,workdir,outDir)) 
        else: os.system("bsub -q %s \'%s/%s/runHResoJob.sh\'"%(myqueue,workdir,outDir))
    
