source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/setup.sh 
export QTHOME=/afs/cern.ch/sw/lcg/external/qt/4.8.4/x86_64-slc6-gcc46-opt/
export G4BASE=/afs/cern.ch/sw/lcg/external/geant4
source $G4BASE/9.6.p02/x86_64-slc6-gcc46-opt/share/Geant4-9.6.2/geant4make/geant4make.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$XERCESCROOT/lib
source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh