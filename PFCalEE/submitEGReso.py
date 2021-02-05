#!/usr/bin/env python
import os, sys
import argparse
from utils import SubmitBase, create_dir

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-t', '--git-tag'  , dest='gittag'    , help='git tag version', default='V00-00-00')
parser.add_argument(      '--nRuns'    , dest='nRuns'     , type=int,   help='number of run, 0-indexed', default=1)
parser.add_argument('-v', '--version'  , dest='version'   , type=int,   help='detector version', required=True)
parser.add_argument('-m', '--model'    , dest='model'     , type=int,   help='detector model', required=True)
parser.add_argument('-a', '--etas'     , dest='etas'      , type=float, help='incidence eta', nargs='+')
parser.add_argument('-p', '--phi'      , dest='phi'       , type=float, help='incidence phi angle in pi unit', default=0.5)
# 1 = hexagons, 2=diamonds, 3=triangles, 4=squares
parser.add_argument(      '--shape'    , dest='shape'     , type=int,  help='shape', default=1) 
parser.add_argument('-b', '--Bfield'   , dest='Bfield'    , type=float, help='B field value in Tesla', required=True)
parser.add_argument('-d', '--datatype' , dest='datatype'  , help='data type or particle to shoot', default='e-')
parser.add_argument('-f', '--datafile' , dest='datafile'  , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_argument('-n', '--nevts'    , dest='nevts'     , type=int, help='number of events to generate', default=1000)
parser.add_argument('-o', '--out'      , dest='out'       , help='output directory', default=os.getcwd())
parser.add_argument(      '--nPuVtx'   , dest='nPuVtxList', type=int, help='pileup scenarios (csv) ', nargs='+', default=[0])
parser.add_argument('-e', '--eosOut'   , dest='eosout'    , help='eos path to save root file to EOS', default='')
parser.add_argument('-E', '--eosIn'    , dest='eosin'     , help='eos path to read input root file from EOS (if empty, it is equal to `--eosOut`.', default='')
parser.add_argument('-g', '--gun'      , dest='dogun'     , help='use particle gun.', action='store_true')
parser.add_argument('-S', '--no-submit', dest='nosubmit'  , help='Do not submit batch job.', action='store_true')
parser.add_argument('--enList'         , dest='enList'    , type=int, help='E_T list to use with gun', nargs='+', default=[5,10,20,30,40,60,80,100,150,200])
parser.add_argument('--interCalib'     , dest='iCalibList', type=int, help='inter calibration list in percentage', nargs='+', default=[3]) #0,1,2,3,4,5,10,15,20,50]
parser.add_argument('--etamean'        , dest='etamean'   , help='mean value of eta ring to save', default=0,  type=float)
parser.add_argument('--deta'           , dest='deta'      , help='width of eta ring', default=0, type=float)
parser.add_argument('--inPathPU'       , dest='inPathPU'  , help='input path for PU files (overrides defaults)', action='store_true')
parser.add_argument('--workdir'        , dest='workdir'  , help='Workdir for outputs (avoid using eos for later plotting)', required=True) #/afs/cern.ch/work/a/amagnan/PFCalEEAna/
opt, _ = parser.parse_known_args()

###################################################################################################
###################################################################################################
###################################################################################################
class SubmitAnalysis(SubmitBase):
    def __init__(self, nSiLayers, label, redofit, **kwargs):
            super(SubmitAnalysis, self).__init__(**kwargs)
            
            self.condor_submit_name_ = 'condorSubmitEGReso.sub'
            self.jobName_ = 'runEGResoJob.sh'
            self.vtx_tag = '$(NPUVTX)'
            self.ic_tag = '$(IC)'
            
            self.nSiLayers = nSiLayers
            self.label = label
            self.redofit = redofit

            self.tags = (self.en_tag, self.eta_tag, self.run_tag, self.ic_tag, self.vtx_tag)
            self.labels = ('energy', 'eta', 'run', 'ic', 'npuvtx')
                        
    @property
    def condor_submit_name(self):
        return self.condor_submit_name_

    @property
    def jobName(self):
        return self.jobName_

    def write_shell_script_file(self):
        with open('{}/{}/{}'.format(self.p.workdir,self.outDir,self.jobName_), 'w') as s:
            s.write('#!/bin/bash\n')

            #input arguments: energy, eta, run, intercalibration and #pu vertices
            s.write('ARGS=`getopt -o "" -l ",energy:,eta:,run:,ic:,npuvtx:" -n "getopts_${0}" -- "$@"`\n')
            s.write('eval set -- "$ARGS"\n')
            s.write('while true; do\n')
            s.write('case "$1" in\n')
            for l,t in zip(self.labels, self.tags):
                s.write('--'+l+')\n')
                s.write('if [ -n "$2" ]; then\n')
                if l=='eta':
                    tmp = "$(echo ${2} | sed 's/\.//');"
                    s.write('{}='.format(self.clean_tag(t))+tmp+'\n')
                else:
                    s.write('{}="${{2}}";\n'.format(self.clean_tag(t)))
                s.write('echo "'+l+': {}";\n'.format(self.shellify_tag(t)))
                s.write('fi\n')
                s.write('shift 2;;\n')
            s.write('--)\n')
            s.write('shift\n')
            s.write('break;;\n')
            s.write('esac\n')
            s.write('done\n\n')
            
            s.write('cd {}\n'.format(os.getcwd()))
            s.write('source g4env.sh\n')
            s.write('cd -\n') 
            s.write('echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"\n')

            s.write('localdir=`pwd`\n')
            s.write('cp -r {}/analysis/data .\n'.format(os.getcwd()))
            s.write('cp -r {}/analysis/scripts .\n'.format(os.getcwd()))

            outTag='version{}_model{}_{}'.format(opt.version,opt.model,bval)
            inTag = outTag
            addToTag_ = '_en{}_eta{}'.format(self.shellify_tag(self.en_tag),self.shellify_tag(self.eta_tag))
            outTag += '_npuvtx{}_ic{}'.format(self.shellify_tag(self.vtx_tag),self.shellify_tag(self.ic_tag))
            inTag  += addToTag_
            outTag += addToTag_
            if self.p.phi!=0.5:
                inTag  += '_phi{n:.{r}f}pi'.format(n=self.p.phi,r=3)
                outTag += '_phi{n:.{r}f}pi'.format(n=self.p.phi,r=3)
            outTag += '_run{}'.format(self.shellify_tag(self.run_tag))
            inTag += '_run{}'.format(self.shellify_tag(self.run_tag))
            
            scriptOutDir = self.outDir + outTag
            create_dir(scriptOutDir)
            
            s.write('{}/analysis/bin/egammaResoWithTruth -c scripts/DefaultConfig.cfg -n {} -i root://eoscms//eos/cms{}/ --digifilePath=root://eoscms//eos/cms{}/ -s HGcal_{}.root -r Digi_{}{}.root -o {}.root --redoStep={}\n'.format(os.getcwd(),self.p.nevts,self.eosDirIn,self.eosDirOut,inTag,self.label,outTag,scriptOutDir,int(self.redofit)))
                        
            s.write('echo "--Local directory is " $localdir\n')
            s.write('ls *\n')
            s.write('echo "--deleting core files: too heavy!!"\n')
            s.write('rm -f core.*\n')
            s.write('echo "--deleting evt-by-evt files: too many!!"\n')
            s.write('cp {}/* {}/{}/\n'.format(scriptOutDir,self.p.workdir,scriptOutDir))
            _tmp = '' if self.redofit else '_nofit'
            s.write('cp {}.root {}/{}' + _tmp + '.root\n'.format(scriptOutDir,self.p.workdir,scriptOutDir))
            s.write('echo "All done"\n')
                

    def write_condor_submission_file(self):
        with open('{}/{}/{}'.format(self.p.workdir,self.outDir,self.condor_submit_name_), 'w') as s:
            s.write('universe = vanilla\n')
            s.write('Executable = {}/{}/{}\n'.format(self.p.workdir,self.outDir,self.jobName))
            s.write('Arguments = ')
            for l,t in zip(self.labels, self.tags):
                s.write('--'+l+' '+t+' ')
            s.write('\n')
            s.write('Requirements = (OpSysAndVer =?= "CentOS7")\n')

            kw = dict(pre='ana_', npuvtx=self.vtx_tag, ic=self.ic_tag, en=self.en_tag,
                      eta=self.eta_tag, run=self.run_tag)
            out_name = self._unique_name(ext='out', **kw)
            err_name = self._unique_name(ext='err', **kw)
            log_name = self._unique_name(ext='log', **kw)
            
            s.write('Output = {}/{}/{}\n'.format(self.p.workdir,self.outDir,out_name))
            s.write('Error = {}/{}/{}\n'.format(self.p.workdir,self.outDir,err_name))
            s.write('Log = {}/{}/{}\n'.format(self.p.workdir,self.outDir,log_name))
            s.write('RequestMemory = 10MB\n')
            s.write('+JobFlavour = "espresso"\n')
            s.write('Queue {nruns} {n}, {ic}, {en}, {eta} from (\n'.format( nruns=self.p.nRuns, n=self.clean_tag(self.vtx_tag), ic=self.clean_tag(self.ic_tag), en=self.clean_tag(self.en_tag), eta=self.clean_tag(self.eta_tag) ))

            for nvid in self.p.nPuVtxList:
                for icid in self.p.iCalibList:
                    for et in self.p.enList:
                        for eta in self.p.etas:
                            s.write('{}, {}, {}, {}\n'.format(nvid,icid,et,str(eta)))
            s.write(')')


###################################################################################################
###################################################################################################
###################################################################################################
bval = 'BON' if opt.Bfield>0 else 'BOFF'
lab = '200u'
redofit = True
nSiLayers = 2

odir = '{}/git{}/version{}/model{}/{}/{}/'.format(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,lab)
if opt.phi!=0.5:
    odir += '_phi{n:.{r}f}pi/'.format(n=self.p.phi,r=3)
edirout='{}/git{}/{}'.format(opt.eosout,opt.gittag,opt.datatype)
edirin='{}/git{}/{}'.format(opt.eosin,opt.gittag,opt.datatype) if opt.eosin != '' else edirout

create_dir('{}/{}'.format(opt.workdir,odir))
      
subana = SubmitAnalysis(nSiLayers=nSiLayers, label=lab, redofit=redofit,
                        outDir=odir, eosDirIn=edirin, eosDirOut=edirout, bfield=bval, params=opt)
subana.write_shell_script_file()
subana.write_condor_submission_file()

os.system('chmod u+rwx {}/{}/{}'.format(opt.workdir,odir,subana.jobName))
if opt.nosubmit:
    os.system('echo condor_submit {}/{}/{}'.format(opt.workdir,odir,subana.condor_submit_name)) 
else:
    os.system('condor_submit {}/{}/{}'.format(opt.workdir,odir,subana.condor_submit_name)) 
