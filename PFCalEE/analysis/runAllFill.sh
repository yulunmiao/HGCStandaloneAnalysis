#!/bin/sh

#INPATH=root://eoscms//eos/cms/store/user/amagnan/HGCalEEGeant4/gitV00-02-02/
INPATH=root://eoscms//eos/cms/store/user/msun/V12
GITVERSION=gitV00-02-07

for eta in 20 #25 30
do
    alpha=0.135
    if (( "$eta"=="25" )); then
	alpha=0.082
    elif (( "$eta"=="30" )); then
	alpha=0.050
    fi

    for v in 12
    do
	for p in e-
	do
	    for t in 2
	    do
		MYDIR=PLOTS/$GITVERSION/version$v/$p/${t}00um/
		mkdir -p $MYDIR
		for r in 0
		do
	    #for e in 10 15 18 20 25 30 35 40 45 50 60 80
	    #for e in 10 15 18 20 25 30 35 40 45 50 60 80
	    #for e in 50 60 80 100 150 200 300 400 500
		#for e in 5 10 15 20 25 30 40 50 60 100 150 200 300 400 500 1000 2000
		    for e in 10 50 100 200
		    do
		    echo "Processing v${v}_p${p}_e${e}_run${r}"
		    #./bin/validation 0 $INPATH/eta$eta/$GITVERSION/$p/ HGcal_version${v}_model1_BOFF_e${e}.root Digi_version${v}_model1_BOFF_e${e}.root $MYDIR/validation_${t}00um_e${e}.root ${t} >& val_si${t}_v${v}_p${p}_e${e}.log &
		    ./bin/egammaResolution 0 $INPATH/eta$eta/$GITVERSION/$p/ HGcal_version${v}_model2_BOFF_e${e}_alpha${alpha}.root Digi_version${v}_model2_BOFF_e${e}_alpha${alpha}.root $MYDIR/eta${eta}_e${e}.root ${t} >& egreso_eta${eta}_e${e}.log &
		done
	    done
	done
    done
done
done
