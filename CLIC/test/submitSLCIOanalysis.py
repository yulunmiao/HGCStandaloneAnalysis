#!/usr/bin/env python

import os,sys
import optparse
import commands


usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                 , default='8nm')
parser.add_option('-R', '--R'          ,    dest='requirementtoBatch' , help='requirement for batch queue' , default='pool>30000')
parser.add_option('-i', '--i'          ,    dest='input'              , help='input directory'             , default='/store/cmst3/user/psilva/PFcal/h_nunu')
parser.add_option('-l', '--localout'   ,    dest='locoutput'          , help='local output directory'      , default='${PWD}/FARM/')
parser.add_option('-s', '--script'     ,    dest='script'             , help='analysis script'             , default='${PWD}/test/runSLCIOanalysis.py')
parser.add_option('-o', '--options'    ,    dest='options'            , help='options'                     , default='-j CA8')
(opt, args) = parser.parse_args()


#prepare output
outDir=os.path.expandvars(opt.locoutput)
os.system('mkdir -p %s/output'%outDir)


j=0
fList=commands.getstatusoutput('cmsLs %s | awk \'{print $5}\''%opt.input)[1].split('\n')

for f in fList:
    if len(f)==0 : continue

    j=j+1
    jobName='ilcsim%d'%j

    #create the steer for the job

    #create a standalone job
    scriptFile = open('%s/runJob_%d.sh'%(outDir,j), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source /afs/cern.ch/eng/clic/software/setupSLC6.sh \n')
    scriptFile.write('python %s -i %s -o %s/output/job_%d.root %s \n'%(opt.script,f,outDir,j,opt.options))
    scriptFile.write('echo "All done for job %d" \n'%j)
    scriptFile.close()

    os.system('chmod u+rwx %s/runJob_%d.sh'%(outDir,j))
    os.system("bsub -q %s -R \"%s\" -J %s \'%s/runJob_%d.sh\'"%(opt.queue,opt.requirementtoBatch,jobName,outDir,j))

