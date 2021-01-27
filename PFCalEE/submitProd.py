#!/usr/bin/env python

import os, sys, errno
import argparse
import math
import random

git_tag=os.popen('git describe --tags --abbrev=0').read()

usage = 'usage: %prog [options]'
parser = argparse.ArgumentParser(usage)
parser.add_argument('-s', '--short-queue' , dest='squeue'     , help='short batch queue [%default]' , default='tomorrow')
parser.add_argument('-q', '--long-queue'  , dest='lqueue'     , help='long batch queue [%default]'  , default='nextweek')
parser.add_argument('-t', '--git-tag'     , dest='gittag'     , help='git tag version [%default]'   , default=git_tag)
parser.add_argument(      '--nRuns'       , dest='nRuns'      , type=int,   help='number of run, 0-indexed', default=-1)
parser.add_argument('-v', '--version'     , dest='version'    , type=int,   help='detector version', default=3)
parser.add_argument('-m', '--model'       , dest='model'      , type=int,   help='detector model', default=3)
parser.add_argument(      '--granularity' , dest='granularity', type=int,   help='lateral granularity (0=HD,1=LD) [%default]', default=1)
parser.add_argument('-a', '--etas'        , dest='etas'       , type=float, help='incidence eta', nargs='+')
parser.add_argument('-p', '--phi'         , dest='phi'        , type=float, help='incidence phi angle in pi unit' , default=0.5)
parser.add_argument(      '--shape'       , dest='shape'      , type=int,   help='shape', default=1) # 1 = hexagons, 2=diamonds, 3=triangles, 4=squares
parser.add_argument('-b', '--Bfield'      , dest='Bfield'     , type=float, help='B field value in Tesla'       , default=0)
parser.add_argument('-d', '--datatype'    , dest='datatype'   ,             help='data type or particle to shoot', default='e-')
parser.add_argument('-f', '--datafile'    , dest='datafile'   ,             help='full path to HepMC input file', default='') #data/example_MyPythia.dat
parser.add_argument('-F', '--datafileeos' , dest='datafileeos',             help='EOS path to HepMC input file', default='') #/eos/cms/store/cmst3/group/hgcal/HGCalMinbias/Pythia8/
parser.add_argument('-n', '--nevts'       , dest='nevts'      , type=int,   help='number of events to generate' , default=1000)
parser.add_argument('-o', '--out'         , dest='out'        ,             help='output directory'             , default=os.getcwd() )
parser.add_argument('-e', '--eos'         , dest='eos'        ,             help='eos path to save root file to EOS',         default='')
parser.add_argument('-g', '--gun'         , dest='dogun'      ,             help='use particle gun.', action="store_true")
parser.add_argument(      '--enList'      , dest='enList'     , type=int,   help='E_T list to use with gun [%default]', nargs='+', default=[5,10,20,30,40,60,80,100,150,200])
parser.add_argument('-S', '--no-submit'   , dest='nosubmit'   ,             help='Do not submit batch job.', action="store_true")
(opt, args) = parser.parse_known_args()


###################################################################################################
###################################################################################################
###################################################################################################
class SubmitProd:
    def __init__(self, outDir, eosDir, bfield, params):
        #variables
        self.outDir = outDir
        self.eosDir = eosDir
        self.p = params
        self.bfield = bfield
        if self.bfield not in ('BON', 'BOFF'):
            raise ValueError('[submitProd.py] The magnetic filed must be either ON or OFF.')
        self.en_tag = '$(ENERGY)'
        self.eta_tag = '$(ETA)'
        self.run_tag = '$(Process)'
        self.condor_submit_name = 'condorSubmitProd.sub'

        #lambda functions
        self.mac_name = lambda e,a,r: 'g4steer_en' + e + '_eta' + a + '_run' + r + '.mac'
        self.clean_tag = lambda t: t.strip('$').strip('(').strip(')')
        self.shellify_tag = lambda t: t.replace('(', '{').replace(')','}')
        self.remove_dot = lambda t: t.replace('.', '')

        #other operations
        self.create_dir(self.outDir)

    def create_dir(self, d):
        try:        
            os.makedirs(d)              
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def gen_uniform_int_random_seeds_(low, high, size):
        random.seed()
        r = random.uniform(low=low, high=high, size=size)
        return [int(x) for x in r]
        
    def write_shell_script_file(self):
        """
        Writes a general shell file which runs the neration step for particular set of run,
        energy and eta values. The file is stored under `d`/runJob.sh

        Args: -p:      Command line arguments
              -d:      Output directory
        """
        with open('{}/runJob.sh'.format(self.outDir), 'w') as s:
            s.write('#!/usr/bin/env bash\n')

            #input arguments: energy, eta and run
            s.write('ARGS=`getopt -o "" -l ",energy:,eta:,run:" -n "getopts_${0}" -- "$@"`\n')
            s.write('eval set -- "$ARGS"\n')
            s.write('while true; do\n')
            s.write('case "$1" in\n')
            s.write('--energy)\n')
            s.write('if [ -n "$2" ]; then\n')
            s.write('{}="${{2}}";\n'.format(self.clean_tag(self.en_tag)))
            s.write('echo "Energy: {}";\n'.format(self.shellify_tag(self.en_tag)))
            s.write('fi\n')
            s.write('shift 2;;\n')
            s.write('--eta)\n')
            s.write('if [ -n "$2" ]; then\n')
            s.write("{}=$(echo ${{2}} | sed 's/\.//');\n".format(self.clean_tag(self.eta_tag)))
            s.write("echo \"Eta: {}\";\n".format(self.shellify_tag(self.eta_tag)))
            s.write('fi\n')
            s.write('shift 2;;\n')
            s.write('--run)\n')
            s.write('if [ -n "$2" ]; then\n')
            s.write('{}="${{2}}";\n'.format(self.clean_tag(self.run_tag)))
            s.write('echo "Run: {}";\n'.format(self.shellify_tag(self.run_tag)))
            s.write('fi\n')
            s.write('shift 2;;\n')
            s.write('--)\n')
            s.write('shift\n')
            s.write('break;;\n')
            s.write('esac\n')
            s.write('done\n\n')
            
            s.write('localdir=`pwd`\n')
            s.write('export HOME={}\n'.format(os.environ['HOME']))
            s.write('cd {}/\n'.format(os.getcwd()))
            s.write('source g4env.sh\n')
            s.write('cd $localdir\n')

            if len(self.p.datafileeos)>0:
                s.write('eos cp %s/%s %s\n'%(self.p.datafileeos,self.p.datafile,self.p.datafile))

            mac_shell_name = self.mac_name(self.shellify_tag(self.en_tag), self.shellify_tag(self.eta_tag),
                                           self.shellify_tag(self.run_tag))
            s.write('cp {}/{} .\n'.format(self.outDir, mac_shell_name))
            cmd = ( 'PFCalEE {} --model {} --version {} --eta {} --shape {}'
                    .format(mac_shell_name, self.p.model, self.p.version,
                            self.shellify_tag(self.eta_tag), self.p.shape) )
            if not self.p.granularity: cmd += ' --fineGranularity'
            s.write(cmd + '\n')
            outTag = 'version{}_model{}_{}'.format(self.p.version, self.p.model, self.bfield)
            outTag += '_en{}_eta{}'.format(self.shellify_tag(self.en_tag),self.shellify_tag(self.eta_tag)) 
            if self.p.phi != 0.5: outTag += '_phi{n:.{r}f}pi'.format(n=self.p.phi,r=3)
            outTag += '_run{}'.format(self.shellify_tag(self.run_tag))
            s.write('mv PFcal.root HGcal_{}.root\n'.format(outTag))
            s.write('echo "--Local directory is " $localdir >> g4_{}.log\n'.format(outTag))
            s.write('echo home=$HOME >> g4_{}.log\n'.format(outTag))
            s.write('echo path=$PATH >> g4_{}.log\n'.format(outTag))
            s.write('echo ldlibpath=$LD_LIBRARY_PATH >> g4_{}.log\n'.format(outTag))
            s.write('ls -ltrh * >> g4_{}.log\n'.format(outTag))
            if len(self.p.eos)>0:
                #s.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setuself.p.sh | sed "s/alias /export my/" > eosenv.sh\n')
                #s.write('source eosenv.sh\n')
                s.write('eos mkdir -p %s\n'%eosDir)
                s.write('eos cp HGcal_%s.root %s/HGcal_%s.root\n'%(outTag,eosDir,outTag))
                s.write('if (( "$?" != "0" )); then\n')
                s.write('echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> g4{}.log\n'.format(outTag))
                s.write('else\n')
                s.write('eossize=`eos ls -l %s/HGcal_%s.root | awk \'{print $5}\'`\n'%(eosDir,outTag))
                s.write('localsize=`ls -l HGcal_%s.root | awk \'{print $5}\'`\n'%(outTag))
                s.write('if [ $eossize != $localsize ]; then\n')
                s.write('echo " --- Copy of sim file to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> g4_{}.log\n'.format(outTag))
                s.write('else\n')
                s.write('echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> g4_{}.log\n'.format(outTag))
                s.write('echo " --- File PFcal.root successfully copied to EOS: {ed}/HGcal_{ot}.root" >> g4_{ot}.log\n'.format(ed=eosDir,ot=outTag))
                s.write('rm HGcal_{}.root\n'.format(outTag))
                s.write('fi\n')
                s.write('fi\n')

            s.write('echo "--deleting core files and hepmc files: too heavy!!"\n')
            s.write('rm core.*\n')
            if len(self.p.datafileeos)>0:
                s.write('rm {}\n'.format(self.p.datafile))
            s.write('cp * {}/\n'.format(self.outDir))
            s.write('echo "All done"\n')

    def write_geant4_files(self):
        """
        Writes all required geant4 input files, one
        for each run (different seed) and energy.
        """
        niters = len(self.p.nRuns)*len(self.p.enList)*len(self.p.etas)
        gen_kwargs = dict(low=0, high=100000, size=niters)
        seeds1 = gen_uniform_int_random_seeds_(**gen_kwargs)
        seeds2 = gen_uniform_int_random_seeds_(**gen_kwargs)
        
        for run in range(self.p.nRuns):
            for iet,et in enumerate(self.p.enList):
                for ieta,eta in enumerate(self.p.etas):
                    gen_idx = ( ( run * len(self.p.enList) * len(self.p.etas) ) +
                                ( iet * len(self.p.etas) ) +
                                ( ieta ) )
                    assert(gen_idx < niters)
                    
                    with open('{}/{}'.format(self.outDir, self.mac_name(str(et), self.remove_dot(str(eta)), str(run))), 'w') as s:
                        s.write('/control/verbose 0\n')
                        s.write('/control/saveHistory\n')
                        s.write('/run/verbose 0\n')
                        s.write('/event/verbose 0\n')
                        s.write('/tracking/verbose 0\n')
                        s.write('/N03/det/setField {n:.{r}f} T\n'.format(n=self.p.Bfield,r=1))
                        s.write('/N03/det/setModel {}\n'.format(self.p.model))
                        s.write('/random/setSeeds {} {}\n'.format(seeds1[gen_idx], seeds2[gen_idx]) )
                        if self.p.dogun :
                            s.write('/generator/select particleGun\n')
                            s.write('/gun/particle {} \n'.format(self.p.datatype))
                            en = et*math.cosh(eta) if eta<5 else et
                            s.write('/gun/energy {n:.{r}f} GeV\n'.format(n=en, r=6))
                            if self.p.model != 2:
                                alpha = 2*math.atan(math.exp(-1.*eta));
                                s.write('/gun/direction {} {} {}\n'.format(math.cos(math.pi*self.p.phi)*math.sin(alpha),math.sin(math.pi*self.p.phi)*math.sin(alpha),math.cos(alpha)))
                        else :
                            s.write('/generator/select hepmcAscii\n')
                            s.write('/generator/hepmcAscii/open {}\n'.format(self.p.datafile))
                            s.write('/generator/hepmcAscii/verbose 0\n')
                        s.write('/run/beamOn {}\n'.format(self.p.nevts))


    def write_condor_submission_file(self):
        """
        Writes one single condor submission file, which is expanded to multiple
        jobs for different energies, etas and runs.
        """
        with open('{}/{}'.format(self.outDir,self.condor_submit_name), 'w') as s:
            s.write('universe = vanilla\n')
            s.write('Executable = {}/runJob.sh\n'.format(self.outDir))
            s.write('Arguments = --energy {} --eta {} --run {}\n'.format(self.en_tag, self.eta_tag, self.run_tag))
            s.write('Requirements = (OpSysAndVer =?= "CentOS7")\n')
            s.write('Output = {}/condorTree.out\n'.format(self.outDir))
            s.write('Error = {}/condorTree.err\n'.format(self.outDir))
            s.write('Log = {}/condorTree.log\n'.format(self.outDir))
            s.write('RequestMemory = 10MB\n')
            s.write('+JobFlavour = "nextweek"\n')
            s.write('Queue {nruns} {entag}, {etatag} from (\n'.format( nruns=self.p.nRuns, entag=self.clean_tag(self.en_tag),
                                                                       etatag=self.clean_tag(self.eta_tag) ))
            for et in self.p.enList:
                for eta in self.p.etas:
                    s.write('{}, {}\n'.format(et,str(eta),r=3))
            s.write(')')
###################################################################################################
###################################################################################################
###################################################################################################

bval = 'BON' if opt.Bfield>0 else 'BOFF'
outDir = '{}git_{}/version_{}/model_{}/{}/{}'.format(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,bval)
if opt.phi != 0.5: outDir='{out}/phi_{n:.{r}f}pi'.format(out=outDir,n=opt.phi,r=3)
eosDir = '/eos/cms{}/git{}/{}'.format(opt.eos,opt.gittag,opt.datatype)

subprod = SubmitProd(outDir=outDir, eosDir=eosDir, bfield=bval, params=opt)
subprod.write_shell_script_file()
subprod.write_geant4_files()
subprod.write_condor_submission_file()

os.system('chmod u+rwx {}/runJob.sh'.format(outDir))
if opt.nosubmit:
    os.system('echo condor_submit {}/{}'.format(outDir, subprod.condor_submit_name)) 
else:
    os.system('condor_submit {}/{}'.format(outDir, subprod.condor_submit_name))
