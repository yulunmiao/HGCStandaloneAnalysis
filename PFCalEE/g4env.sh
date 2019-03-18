#!/bin/bash
export USERBASE=`pwd`
ARCH=x86_64-slc6-gcc46-opt
BASESW=/cvmfs/sft.cern.ch/lcg/external

#gcc
source $BASESW/gcc/4.6.3/x86_64-slc6/setup.sh

#hepmc,fastjet
export HEPMC_DIR=$BASESW/HepMC/2.06.08/${ARCH}
export FASTJET_INSTALL=${BASESW}/fastjet/3.0.3/${ARCH}


#Geant4
cd /cvmfs/geant4.cern.ch/geant4/9.6.p04/${ARCH}/share/Geant4-9.6.4/geant4make/
source geant4make.sh
cd - &> /dev/null

export BOOSTSYS=/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$XERCESCROOT/lib:$HEPMC_DIR/lib:$USERBASE/userlib/lib:$USERBASE/analysis/lib:$FASTJET_INSTALL/lib:$BOOSTSYS/lib

source $BASESW/ROOT/5.34.00/x86_64-slc5-gcc46-opt/root/bin/thisroot.sh

export PATH=$PATH:$FASTJET_INSTALL/bin


#visualization disable for the moment
#export LD_LIBRARY_PATH=/usr/local/geantExternalLibs:$LD_LIBRARY_PATH
#export QTHOME=${BASESW}qt/4.8.7-cms/ 
#export DAWNHOME=/afs/cern.ch/sw/lcg/external/dawn/3_88a/x86_64-slc5-gcc43-opt/
#export G4DAWNFILE_DEST_DIR=${USERBASE}/DawnFiles/
#export PATH=$DAWNHOME/bin:$PATH
