#!/usr/bin/env python

import os,sys
import optparse
import commands


usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                 , default='1nh')
parser.add_option('-R', '--R'          ,    dest='requirementtoBatch' , help='requirement for batch queue' , default='pool>30000')
parser.add_option('-i', '--i'          ,    dest='input'              , help='input directory'             , default='/castor/cern.ch/grid/ilc/prod/clic/1.4tev/h_nunu/ILD/DST/00001504/000/')
parser.add_option('-o', '--out'        ,    dest='output'             , help='output directory'            , default='/store/cmst3/user/psilva/PFcal/')
parser.add_option('-l', '--localout'   ,    dest='locoutput'          , help='local output directory'      , default='${PWD}/FARM/')
parser.add_option('-s', '--steer'      ,    dest='steer'              , help='steer file'                  , default='test/fastjet.xml')
(opt, args) = parser.parse_args()


#prepare output
os.system('cmsMkdir %s'%opt.output)
outDir=os.path.expandvars(opt.locoutput)
os.system('mkdir -p %s'%outDir)

isCastor=False
lsCommand='cmsLs %s | awk \'{print $5}\''%opt.input
if opt.input.find('castor/cern.ch')>=0 :
    isCastor=True
    lsCommand='rfdir %s | awk \'{print $9}\''%opt.input

j=0
fList=commands.getstatusoutput(lsCommand)[1].split('\n')
for f in fList:

    #get the base name of the file (normalize castor and EOS outputs)
    f=os.path.basename(f)
    if f.find('.slcio')<0 : continue

    j=j+1
    jobName='ilcsim%d'%j
    
    #create a standalone Marlin job
    scriptFile = open('%s/runJob_%d.sh'%(outDir,j), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source /afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/init_ilcsoft.sh\n')
    if isCastor : scriptFile.write('rfcp %s/%s /tmp/%s \n'%(opt.input,f,f))
    else        : scriptFile.write('cmsStage %s/%s /tmp/%s \n'%(opt.input,f,f))
    scriptFile.write('Marlin %s/../%s --MyLCIOOutputProcessor.LCIOOutputFile=/tmp/fj_%s --global.LCIOInputFiles=/tmp/%s \n'%(outDir,opt.steer,f,f))
    scriptFile.write('cmsStage /tmp/fj_%s %s \n'%(f,opt.output))
    scriptFile.write('rm /tmp/*.slcio \n')
    scriptFile.write('echo "All done for job %d" \n'%j)
    scriptFile.close()

    os.system('chmod u+rwx %s/runJob_%d.sh'%(outDir,j))
    os.system("bsub -q %s -R \"%s\" -J %s \'%s/runJob_%d.sh\'"%(opt.queue,opt.requirementtoBatch,jobName,outDir,j))
    
