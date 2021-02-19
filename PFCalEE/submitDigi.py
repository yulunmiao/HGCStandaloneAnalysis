#!/usr/bin/env python

import os, sys
import argparse
from utils import SubmitBase

git_tag=os.popen('git describe --tags --abbrev=0').read()

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-t', '--git-tag'  , dest='gittag'    , help='git tag version', default='V00-00-00')
parser.add_argument(      '--nRuns'    , dest='nRuns'     , type=int,   help='number of run, 0-indexed', default=-1)
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
opt, _ = parser.parse_known_args()

###################################################################################################
###################################################################################################
###################################################################################################
class SubmitDigi(SubmitBase):
    def __init__(self, nSiLayers, label, addNoise,
                 granularity, threshold, noise, pathPU, **kwargs):
            super(SubmitDigi, self).__init__(**kwargs)
            
            self.condor_submit_name_ = 'condorSubmitDigi.sub'
            self.jobName_ = 'runDigiJob.sh'
            self.vtx_tag = '$(NPUVTX)'
            self.ic_tag = '$(IC)'
            
            self.nSiLayers = nSiLayers
            self.addNoise = 'false' if addNoise==False else 'true'
            self.granularity = granularity
            self.threshold = threshold
            self.noise = noise
            self.label = label
            self.pathPU = pathPU

            self.tags = (self.en_tag, self.eta_tag, self.run_tag, self.ic_tag, self.vtx_tag)
            self.labels = ('energy', 'eta', 'run', 'ic', 'npuvtx')

    @property
    def condor_submit_name(self):
        return self.condor_submit_name_

    @property
    def jobName(self):
        return self.jobName_

    def write_shell_script_file(self):
        with open('{}/{}'.format(self.outDir,self.jobName_), 'w') as s:
            s.write('#!/usr/bin/env bash\n')
                    
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
                        
            s.write('localdir=`pwd`\n')
            s.write('export HOME={}\n'.format(os.environ['HOME']))
            s.write('cd {}\n'.format(os.getcwd()))
            s.write('source g4env.sh\n')
            s.write('echo $PATH\n')
            s.write('cd $localdir\n')
            
            outTag = 'version{}_model{}_{}'.format(self.p.version,self.p.model,bval)
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
            substr = '{cwd}/userlib/bin/digitizer -c {cwd}/userlib/DigiConfig.cfg -n {n} -i {i}/HGcal_{tag}.root -o $localdir/ --granulStr={g}  --noiseStr={noise} --threshStr={thresh} --interCalib={ic} --nSiLayers={nl} --nPU={npu} --puPath={path} '.format(cwd=os.getcwd(),n=self.p.nevts,i=self.eosDirIn,tag=inTag,g=self.granularity,noise=self.noise,thresh=self.threshold,ic=self.shellify_tag(self.ic_tag),nl=self.nSiLayers,npu=self.shellify_tag(self.vtx_tag),path=self.pathPU)
            if self.p.etamean>1.3:
                s.write(substr+'--etamean={em1:.{em2}f} --deta={p1:.{p2}f}'.format(em1=self.p.etamean,em2=2,p1=self.p.deta,p2=2))
            else:
                s.write(substr)
            s.write('-a {}\n'.format(self.addNoise))
            s.write('echo "--Local directory is " $localdir\n')
            s.write('echo home=$HOME\n')
            s.write('echo path=$PATH\n')
            s.write('echo ldlibpath=$LD_LIBRARY_PATH\n')
            s.write('ls *\n')
            if len(self.p.eosout)>0:
                s.write('eos mkdir -p {}\n'.format(self.eosDirOut))
                s.write('eos cp $localdir/DigiPFcal.root {}/Digi_{}{}.root\n'.format(self.eosDirOut,self.label,outTag))
                s.write('if (( "$?" != "0" )); then\n')
                s.write('echo " --- Problem with copy of file DigiPFcal.root to EOS. Keeping locally."\n')
                s.write('else\n')
                s.write("eossize=`eos ls -l {}/Digi_{}{}.root | awk \'{{print $5}}\'`\n".format(self.eosDirOut,self.label,outTag))
                s.write("localsize=`ls -l DigiPFcal.root | awk \'{print $5}\'`\n")
                s.write('if [ ${eossize} != ${localsize} ]; then\n')
                s.write('echo " --- Copy of digi file to eos failed. Localsize = ${localsize}, eossize = ${eossize}. Keeping locally..."\n')
                s.write('else\n')
                s.write('echo " --- Size check done: Localsize = ${localsize}, eossize = ${eossize}"\n')
                s.write('echo " --- File DigiPFcal.root successfully copied to EOS: {}/Digi_{}{}.root"\n'.format(self.eosDirOut,self.label,outTag))
                s.write('rm DigiPFcal.root\n')
                s.write('fi\n')
                s.write('fi\n')
            else:
                s.write('mv DigiPFcal.root Digi_{}{}.root\n'.format(self.label,outTag))
                s.write('echo "--deleting core files: too heavy!!"\n')
                s.write('rm core.*\n')
                s.write('cp * {}/\n'.format(self.outDir))
                s.write('echo "All done"\n')


    def write_condor_submission_file(self):
        with open('{}/{}'.format(self.outDir,self.condor_submit_name_), 'w') as s:
                s.write('universe = vanilla\n')
                s.write('Executable = {}/{}\n'.format(self.outDir,self.jobName_))
                s.write('Arguments = ')
                for l,t in zip(self.labels, self.tags):
                    s.write('--'+l+' '+t+' ')
                s.write('\n')
                s.write('Requirements = (OpSysAndVer =?= "CentOS7")\n')

                t = ( ('prefix', 'digi'), ('npuvtx', self.vtx_tag), ('ic', self.ic_tag),
                      ('en', self.en_tag), ('eta', self.eta_tag), ('run', self.run_tag) )
                out_name = self._unique_name( t + (('ext', 'out'),) )
                err_name = self._unique_name( t + (('ext', 'err'),) )
                log_name = self._unique_name( t + (('ext', 'log'),) )
                
                s.write('Output = {}/{}\n'.format(self.outDir,out_name))
                s.write('Error = {}/{}\n'.format(self.outDir,err_name))
                s.write('Log = {}/{}\n'.format(self.outDir,log_name))
                s.write('RequestMemory = 30MB\n')
                s.write('+JobFlavour = "microcentury"\n')
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
odir = '{}/git{}/version_{}/model_{}/{}/{}/{}'.format(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,bval,lab)
edirout = '/eos/cms{}/git{}/{}'.format(opt.eosout,opt.gittag,opt.datatype)
edirin = 'root://eoscms//eos/cms{}/git{}/{}'.format(opt.eosin,opt.gittag,opt.datatype) if opt.eosin != '' else edirout
nSiLayers = 2

nmult = ('0-27:0.27', '0-27:0.13', '0-27:0.07', '0-27:0.13')
n63 = ('0-51:0.27,53-68:0.15', '0-51:0.13,53-68:0.15', '0-51:0.07,53-68:0.15', '0-51:0.13,53-68:0.15')
n70 = ('0-25:0.27', '0-25:0.13', '0-25:0.07', '0-25:0.13')
n73 = ('0-46:0.27,48-61:0.15', '0-46:0.13,48-61:0.15', '0-46:0.07,48-61:0.15', '0-46:0.13,48-61:0.15')
def get_noise(noise, l):
    if l=='100u':   return noise[0]
    elif l=='200u': return noise[1]
    elif l=='300u': return noise[2]
    else:           return noise[3]

pudflt = 'root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V12/MinBias/'
vdict = {8:   dict(puFile='root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V12/MinBias/',
                   granularity='0-20:4,21-30:6', noise='0-30:0.14', threshold='0-30:2'),
         13:  dict(puFile='root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V13/MinBias/',
                   granularity='0-29:4,30-65:4', noise='0-65:0.15', threshold='0-65:5'),
         21:  dict(puFile=pudflt, granularity='0-23:6,24-33:8', noise='0-33:0.14', threshold='0-33:2'),
         22:  dict(puFile=pudflt, granularity='0-9:8', noise='0-9:0.14', threshold='0-9:2'),
         23:  dict(puFile=pudflt, granularity='0-53:12', noise='0-53:0.14', threshold='0-53:2'),
         24:  dict(puFile=pudflt, granularity='0-23:6,24-33:8', noise='0-33:0.14', threshold='0-33:2'),
         25:  dict(puFile='root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V25/MinBias/',
                   granularity='0-29:4,30-41:4,42-53:8', noise='0-41:0.14,42-53:0.2', threshold='0-53:5'),
         26:  dict(puFile=pudflt, granularity='0-29:4,30-41:4,42-53:8', noise='0-41:0.14,42-53:0.2', threshold='0-53:5'),
         27:  dict(puFile=pudflt, granularity='0-11:4,12-23:8', noise='0-11:0.14,12-23:0.2', threshold='0-23:5'),
         28:  dict(puFile=pudflt, granularity='0-11:8', noise='0-11:0.2', threshold='0-11:5'),
         30:  dict(puFile=pudflt, granularity='0-27:4', noise='0-27:0.14', threshold='0-27:5'),
         31:  dict(puFile=pudflt, granularity='0-11:4,12-23:8', noise='0-11:0.14,12-23:0.2', threshold='0-23:5'),
         32:  dict(puFile=pudflt, granularity='0-11:8', noise='0-11:0.2', threshold='0-11:5'),
         33:  dict(puFile='root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V33/MinBias/pile/gitV00-03-07/e-/',
                   granularity='0-27:1,28-39:1,40-51:1', noise='0-39:0.,40-51:0.', threshold='0-51:5'),
         34:  dict(puFile=pudflt, granularity='0-23:4', noise='0-23:0.14', threshold='0-23:5'),
         35:  dict(puFile=pudflt, granularity='0-17:4', noise='0-17:0.14', threshold='0-17:5'),
         36:  dict(puFile=pudflt, granularity='0-23:4,24-34:4,35-46:8', noise='0-34:0.14,35-46:0.2', threshold='0-46:5'),
         37:  dict(puFile=pudflt, granularity='0-17:4,18-26:4,27-38:8', noise='0-26:0.14,27-38:0.2', threshold='0-38:5'),
         38:  dict(puFile=pudflt, granularity='0-10:4,11-22:8', noise='0-10:0.14,11-22:0.2', threshold='0-22:5'),
         39:  dict(puFile=pudflt, granularity='0-8:4,9-20:8', noise='0-8:0.14,9-20:0.2', threshold='0-20:5'),
         60:  dict(puFile=( 'root://eoscms//eos/cms/store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR/gitV08-05-00/MinBiasLarge/' if lab=='' else
                            'root://eoscms//eos/cms/store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR/gitV08-05-00/MinBiasSmall/' ),
                   granularity='0-27:1', threshold='0-27:5', noise=get_noise(nmult,lab)),     
         61:  dict(puFile=pudflt, granularity='0-39:1', noise='0-23:0.12,24-39:0.15', threshold='0-39:5'),
         62:  dict(puFile=pudflt, granularity='0-15:1', noise='0-15:0.15', threshold='0-15:5'),
         63:  dict(puFile=('root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalTDR/gittestV8/MinBiasSmall/'
                           if lab==''
                           else 'root://eoscms//eos/cms/store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR/gitV08-01-00/MinBiasLarge/'),
                   granularity='0-68:1', threshold='0-68:5', noise=get_noise(n63,lab)),
         64:  dict(puFile=pudflt, granularity='0-27:1', threshold='0-27:5', noise=get_noise(nmult,lab)),
         65:  dict(puFile=pudflt, granularity='0-27:1', threshold='0-27:5', noise=get_noise(nmult,lab)),
         66:  dict(puFile=pudflt, granularity='0-23:1', noise='0-23:0.12', threshold='0-23:5'),
         67:  dict(puFile=( 'root://eoscms//eos/cms/store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR/gitV08-05-00/MinBiasLarge/' if lab=='' else
                            'root://eoscms//eos/cms/store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR/gitV08-05-00/MinBiasSmall/' ),
                   granularity='0-27:1', threshold='0-27:5', noise=get_noise(nmult,lab)),
         70:  dict(puFile=pudflt, granularity='0-25:1', threshold='0-25:5', noise=get_noise(n70,lab)),
         73:  dict(puFile=pudflt, granularity='0-61:1', threshold='0-61:5', noise=get_noise(n70,lab)),
         83:  dict(puFile=pudflt, granularity='0-61:1', threshold='0-61:5', noise=get_noise(n70,lab)),
         100: dict(puFile=pudflt, granularity='0-27:4', noise='0-27:0.14', threshold='0-27:5'),
         110: dict(puFile=pudflt, granularity='0-27:4', noise='0-27:0.14', threshold='0-27:5'),
}

gran = vdict[opt.version]['granularity']
thr  = vdict[opt.version]['threshold']
noi  = vdict[opt.version]['noise']
if opt.inPathPU:
    pathPU = 'root://eoscms//eos/cms/{}'.format(opt.inPathPU)
    print('Default value for `pathPU` overriden with {}.'.format(self.pathPU))
else:
    pathPU = vdict[opt.version]['puFile']

subdigi = SubmitDigi(nSiLayers=nSiLayers, label=lab, addNoise=False,
                     granularity=gran, threshold=thr, noise=noi, pathPU=pathPU,
                     outDir=odir, eosDirIn=edirin, eosDirOut=edirout, bfield=bval, params=opt)
subdigi.write_shell_script_file()
subdigi.write_condor_submission_file()

os.system('chmod u+rwx {}/{}'.format(odir, subdigi.jobName))
if opt.nosubmit:
        os.system('echo condor_submit {}/{}'.format(odir, subdigi.condor_submit_name)) 
else:
        os.system('condor_submit {}/{}'.format(odir, subdigi.condor_submit_name))
