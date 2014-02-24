#!/bin/sh


for eta in 20 25 30 35
do
    for type in PedroPU #SimpleSignal PedroPU SimplePU
    do
	for justpileup in `ls /afs/cern.ch/work/p/pdauncey/public/annemarie/$type/eta$eta/*.dat`;
	do
	    echo "$justpileup" > tmp1
	    sed "s|/afs/cern.ch/work/p/pdauncey/public/annemarie/${type}/eta${eta}/pileup0\{0,3\}||" tmp1 > tmp2
	    subdir=e_`sed "s|[A-Z]\.dat||" tmp2`

	    echo " -- processing file "$justpileup 
	    echo " -- subdir " $subdir
	
	    ./submitRunHEPMC.py -q 2nd -v 23 -f $justpileup -d $type/eta$eta/$subdir -n 10 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEGeant4/

	done
    done
done
