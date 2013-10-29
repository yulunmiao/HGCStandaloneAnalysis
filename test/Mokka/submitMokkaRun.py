#!/usr/bin/env python

import os,sys
import optparse
import commands


usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                 , default='1nw')
parser.add_option('-R', '--R'          ,    dest='requirementtoBatch' , help='requirement for batch queue' , default='pool>30000')
parser.add_option('-i', '--i'          ,    dest='input'              , help='input directory'             , default='/store/cmst3/user/psilva/PFcal/pp_WZ/HepEvt')
parser.add_option('-o', '--out'        ,    dest='output'             , help='output directory'            , default='/store/cmst3/user/psilva/PFcal/pp_WZ/')
parser.add_option('-l', '--localout'   ,    dest='locoutput'          , help='local output directory'      , default='${PWD}/FARM/')
(opt, args) = parser.parse_args()


#prepare output
os.system('cmsMkdir %s/DST'%opt.output)
os.system('cmsMkdir %s/REC'%opt.output)
outDir=os.path.expandvars(opt.locoutput)
os.system('mkdir -p %s'%outDir)


j=0
fList=commands.getstatusoutput('cmsLs %s | awk \'{print $5}\''%opt.input)[1].split('\n')
for f in fList:

    j=j+1
    jobName='ilcsim%d'%j

    basef=os.path.basename(f)
    basef=basef.replace('.hepevt','_DST.slcio')
    baserecf=os.path.basename(f)
    baserecf=baserecf.replace('.hepevt','.slcio')

    #create a standalone Mokka and Marlin job
    scriptFile = open('%s/runJob_%d.sh'%(outDir,j), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source /afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/init_ilcsoft.sh\n')
    scriptFile.write('mkdir /tmp/ilcsim%d\n'%(j))
    scriptFile.write('cd /tmp/ilcsim%d\n'%(j))
    scriptFile.write('gtar -xzvf /afs/cern.ch/user/p/psilva/CLIC/PFCal/test/Mokka/MokkaSim.tar.gz -C ./\n')
    scriptFile.write('cmsStage %s ./bbudsc_3evt.hepevt\n'%(f))
    scriptFile.write('echo "Running simulation"\n')
    scriptFile.write('Mokka -M ILD_o1_v05 bbudsc_3evt.steer\n')
    scriptFile.write('cmsStage bbudsc_3evt.slcio %s/REC/%s\n'%(opt.output,baserecf))
    scriptFile.write('echo "Running reconstruction"\n')
    scriptFile.write('Marlin bbudsc_3evt_stdreco.xml\n')
    scriptFile.write('cmsStage bbudsc_3evt_DST.slcio %s/DST/%s\n'%(opt.output,basef))
    scriptFile.write('cd -\n')
    scriptFile.write('rm -rf /tmp/ilcsim%d\n'%(j))
    scriptFile.write('echo "All done for job %d" \n'%j)
    scriptFile.close()

    os.system('chmod u+rwx %s/runJob_%d.sh'%(outDir,j))
    os.system("bsub -q %s -R \"%s\" -J %s \'%s/runJob_%d.sh\'"%(opt.queue,opt.requirementtoBatch,jobName,outDir,j))

