# PFCalEE

Geant4 simulation of a Si/Sci-base sampling calorimeter

Check https://twiki.cern.ch/twiki/bin/view/CMS/HGCalPerformanceStudiesWithG4

Geometry implementation is instantiated by detector versions in an enum - cf. src/DetectorConstruction.cc and src/SamplingSection.cc

A small ntuple is stored with the energy deposits - cf. src/EventAction.cc 

When changing the TTree content, adding homemade classes rather than simple objects, the classes should be added to userlib. 
The dictionary for root to understand the classes also need to be remade. 
Use "make dictionary" before make, inside of userlib/. 
Follow instructions in:
https://twiki.cern.ch/twiki/bin/view/Sandbox/AnnemarieMagnanSandbox
for what needs to be put in the homemade classes and makefile for root to
understand them (see also example in class userlib/include/HGCSSSimHit.hh).

## Setup the environment

source g4env.sh

## Compile

```
mkdir -p userlib/{lib,obj,bin} && cd userlib && make dictionary && make -j 5 && cd - && ./makeG4
```

NB: whenever you pull the code from git again you may get unadvertly the userlib/dict.cc updated.
It's better to make sure that the dictionary is built again by repeating the steps above.

## Jobs submission

You can submit in parallel the runs with `submitProd.py`. 
The `-h` option will print out all the possible arguments to use. Some special ones are

* `-S` to not submit automatically to batch queues
* `-g` to do particleGun (by opposition to hepmc file, see example below)

In case of options conflicts, you can replace python submitProd.py by ./submitProd.py...

### An example example with particle gun:

* edit submitProd.py to set the energy loop to the values wanted.
* For loop is to generate several samples with same stat in parallel.
```
for i in `seq 0 5`; do 
  python submitProd.py -s tomorrow -q nextweek -t V08-07-01 -g -r ${i} -v 3 -m 0 -e /store/cmst3/group/hgcal/Geant4 -o ~/work/ntuples -d e- -n 2500; 
done
```

### an example with hepmc file:

```
./submitProd.py -S -q 2nd -t V00-00-00 -f /afs/cern.ch/work/a/amagnan/CMSSW_6_2_0_SLHC8/src/UserCode/Gen2HepMC/test/VBFH_sel.dat  -v 20 -m 2 -e /store/cmst3/group/hgcal/Geant4 -o ~/work/ntuples -d VBFH -n 1000
```
