#!/bin/sh

INPATH=root://eoscms//eos/cms/store/user/amagnan/HGCalEEGeant4/
#INPATH=root://eoscms//eos/cms/store/user/msun/V12
GITVERSION=gitV00-02-09

etaset=(17 19 21 23 25 27 29)

counter=0

for alpha in 0.361 #0.297 0.244 0.200 0.164 0.134 0.110
do
    eta=${etaset[$counter]}
    echo "counter:$counter, eta=$eta"
    let counter=$counter+1
    for v in 12
    do
	for p in gamma
	do
	    t=2
	    MYDIR=PLOTS/$GITVERSION/version$v/$p/${t}00um/
	    mkdir -p $MYDIR
	    #for r in 0
	    #do
	    for et in 20 30 40 50 60 70 80 90 100 125 150 175 200
	    do
		echo "Processing v${v}_p${p}_et${et}"#_run${r}"
		#if (( "$et"!="80" )); then
		    ./bin/egammaResolution 0 $INPATH/$GITVERSION/$p/ HGcal_version${v}_model2_BOFF_et${et}_alpha${alpha}.root Digi_version${v}_model2_BOFF_et${et}_alpha${alpha}.root $MYDIR/eta${eta}_et${et}.root ${t} >& egreso_eta${eta}_et${et}.log &
		#else
		#    ./bin/egammaResolution 0 $INPATH/$GITVERSION/$p/ HGcal_version${v}_model2_BOFF_et${et}_alpha${alpha}.root Digi_version${v}_model2_BOFF_et${et}_alpha${alpha}.root $MYDIR/eta${eta}_et${et}.root ${t} >& egreso_eta${eta}_et${et}.log
		#fi
	    done
	done
    done
done
