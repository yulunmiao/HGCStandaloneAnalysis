#!/bin/sh

#for eta in 20 25 30 35
#do
    for type in Signal #PedroPU
    do
	#for justpileup in `ls /afs/cern.ch/work/p/pdauncey/public/annemarie/$type/eta$eta/*.dat`;
	for justpileup in `ls /afs/cern.ch/work/p/pdauncey/public/annemarie/version_4/$type/*.dat`;
	do
	    echo "$justpileup" > tmp1
	    #sed "s|/afs/cern.ch/work/p/pdauncey/public/annemarie/${type}/eta${eta}/pileup0\{0,3\}||" tmp1 > tmp2
	    sed "s|/afs/cern.ch/work/p/pdauncey/public/annemarie/version_4/${type}/Photon0\{0,3\}||" tmp1 > tmp2
	    #subdir=eta$eta"_e"`sed "s|[A-Z]\.dat||" tmp2`
	    subdir=e`sed "s|[A-Za-z]\{0,3\}\.dat||" tmp2`

	    echo " -- processing file "$justpileup 
	    echo " -- subdir " $subdir
	
	    ./submitRunHEPMC.py -q 2nd -v 3 -f $justpileup -d $type -s $subdir -n 1000 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEGeant4/ -e /store/user/amagnan/HGCalEEGeant4

	done
    done
#done
