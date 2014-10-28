#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-s', '--short-queue',    dest='squeue'             , help='short batch queue'            , default='1nd')
parser.add_option('-q', '--long-queue' ,    dest='lqueue'             , help='long batch queue'             , default='2nw')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='V00-00-00')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,      type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=3,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=3,      type=int)
parser.add_option('-a', '--alpha'       ,    dest='alpha'              , help='incidence angle in rad'       , default=0,      type=float)
parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0,      type=float)
parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
parser.add_option('-f', '--datafile'    ,    dest='datafile'           , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to process' , default=0,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

enlist=[0]
if opt.dogun : enlist=[20,30,40,50,60,70,80,90,100,125,150,175,200]

#alphaset=[0.361,0.297,0.244,0.200,0.164,0.134,0.110]
alphaset=[0.361]
nPuVtx=140
#etaset=[17,19,21,23,25,27,29]
etaset=[17]

counter=0

for alpha in alphaset :
    eta=etaset[counter]
    counter=counter+1
    for et in enlist :
        
        nevents=opt.nevts
        myqueue=opt.lqueue
        bval="BOFF"
        if opt.Bfield>0 : bval="BON" 
        
        outDir='%s/%s/git%s/version%d/%s/200um/eta%s_et%s_pu%s'%(os.getcwd(),opt.out,opt.gittag,opt.version,opt.datatype,eta,et,nPuVtx)
        eosDir='%s/git%s/%s'%(opt.eos,opt.gittag,opt.datatype)
        
        outlog='%s/egreso_eta%s_et%s_pu%s.log'%(os.getcwd(),eta,et,nPuVtx)
        g4log='egresojob%s.log'%(nPuVtx)
        os.system('mkdir -p %s'%outDir)
   
        #wrapper
        scriptFile = open('%s/runEGResoJob.sh'%(outDir), 'w')
        scriptFile.write('#!/bin/bash\n')
        scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
        #scriptFile.write('cd %s\n'%(outDir))
        outTag='version%d_model%d_%s'%(opt.version,opt.model,bval)
        if et>0 : outTag='%s_et%d'%(outTag,et)
        if alpha>0 : outTag='%s_alpha%3.3f'%(outTag,alpha) 
        if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)
        scriptFile.write('localdir=`pwd`\n')
        scriptFile.write('cp -r %s/data .\n'%os.getcwd())
        if (nPuVtx==0) :
            scriptFile.write('%s/bin/egammaResolution %s root://eoscms//eos/cms%s/ HGcal_%s.root Digi_%s.root %s.root 2 | tee %s\n'%(os.getcwd(),opt.nevts,eosDir,outTag,outTag,outDir,outlog))
        else:
           scriptFile.write('%s/bin/egammaResolution %s root://eoscms//eos/cms%s/ HGcal_%s.root PuMix%s_%s.root %s.root 2 | tee %s\n'%(os.getcwd(),opt.nevts,eosDir,outTag,nPuVtx,outTag,outDir,outlog)) 

        scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
        scriptFile.write('ls * >> %s\n'%(g4log))
        scriptFile.write('cp * %s/\n'%(outDir))
        scriptFile.write('echo "All done"\n')
        scriptFile.close()
    
        #submit
        os.system('chmod u+rwx %s/runEGResoJob.sh'%(outDir))
        if opt.nosubmit : os.system('echo bsub -q %s %s/runEGResoJob.sh'%(myqueue,outDir)) 
        else: os.system("bsub -q %s \'%s/runEGResoJob.sh\'"%(myqueue,outDir))

