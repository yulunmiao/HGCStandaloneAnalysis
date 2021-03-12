export USERBASE=`pwd`
#slc6 setup
#ARCH=x86_64-slc6-gcc46-opt
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/${ARCH}/setup.sh 
#export QTHOME=/afs/cern.ch/sw/lcg/external/qt/4.8.4/${ARCH}/
#export G4BASE=/afs/cern.ch/sw/lcg/external/geant4
#export DAWNHOME=/afs/cern.ch/sw/lcg/external/dawn/3_88a/x86_64-slc5-gcc43-opt/
#export G4DAWNFILE_DEST_DIR=${USERBASE}/DawnFiles/
#export HEPMC_DIR=/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/${ARCH}/
#export FASTJET_INSTALL=/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/${ARCH}/

export BASEINSTALL=/cvmfs/sft.cern.ch/lcg/views/LCG_97/x86_64-centos7-gcc8-opt/

source ${BASEINSTALL}/setup.sh

#cd $G4INSTALL/share/Geant4-10.6.1/geant4make/
#source geant4make.sh
#cd - &> /dev/null

export G4DIR=$G4INSTALL/lib64
mkdir -p $USERBASE/g4build
export G4Build=$USERBASE/g4build

export HepMC_DIR=/cvmfs/sft.cern.ch/lcg/releases/HepMC/2.06.10-1a364/x86_64-centos7-gcc8-opt/


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${BASEINSTALL}/lib64:${HepMC_DIR}/lib:$USERBASE/userlib/lib:$USERBASE/analysis/lib
#slc6 setup
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh
#cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/${ARCH}/root/
#cd /cvmfs/sft.cern.ch/lcg/releases/ROOT/v6-20-00-patches-5b35b/${ARCH}/
#source bin/thisroot.sh
#cd - &> /dev/null
#export PATH=$DAWNHOME/bin:$PATH:$FASTJET_INSTALL/bin
export PATH=$PATH:${BASEINSTALL}/bin:${G4Build}
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:${USERBASE}/userlib/include

