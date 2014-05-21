#!/bin/sh

INPATH=/afs/cern.ch/work/a/amagnan/public/HGCalEEGeant4/

for v in 20;
do
    for s in 0 ; #1 2 3 4 5 6;
    do
	for pu in SimpleSignal #PedroPU #SimplePU #SimpleSignal #PU GammaPU
	do
	    for eta in 25 #20 30 35
	    do
		mkdir -p PLOTS/version_$v/scenario_$s/$pu/eta$eta/
		mkdir -p PLOTS/version_$v/scenario_$s/$pu/PU/eta$eta/
		#mkdir -p PLOTS/version_$v/scenario_$s/$pu"_pipm/eta"$eta/
		#mkdir -p PLOTS/version_$v/scenario_$s/$pu"_pi0/eta"$eta/
		#mkdir -p PLOTS/version_$v/scenario_$s/$pu/eta$eta/xySimHits/
		./bin/fillHistos 0 /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi/version_$v/scenario_$s/$pu/eta$eta/ DigiPFcal.root
		#./bin/fillHistos 0 /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi/version_$v/scenario_$s/$pu/eta$eta/ DigiPFcal.root 0 0 1
		#mkdir -p PLOTS/version_$v/$pu/eta$eta/
		#mkdir -p PLOTS/version_$v/$pu/PU/eta$eta/
		#mkdir -p PLOTS/version_$v/$pu/eta$eta/xySimHits/
		#./bin/fillHistos 10 $INPATH/version_$v/$pu/eta$eta/ PFcal.root 0 1
		#./bin/fillHistos 20 $INPATH/version_$v/$pu/eta$eta/ PFcal.root 0 0 1
	    done
	done
    done
done
