# =====
# PFCal
# =====
#
# Check https://twiki.cern.ch/twiki/bin/view/CMS/PFForwardCalorimeterStudies
#

# Run FastJet over a set of SLCIO files in castor/eos (SLC5)

python test/submitMarlinRun.py -i /castor/cern.ch/grid/ilc/prod/clic/1.4tev/h_nunu/ILD/DST/00002022/000 -o /store/cmst3/user/psilva/PFcal/hbb_nunu_ILD -s test/fastjet_newpfos.xml -l ${PWD}/FARMILD

python test/submitMarlinRun.py -i /store/cmst3/user/psilva/PFcal/pp_WZ/DST/  -o /store/cmst3/user/psilva/PFcal/pp_WZ/DST/FJ -s test/fastjet_pp.xml -l ${PWD}/FARMWZ

# Run simple analyzer over SLCIO files in EOS (SLC6)

python test/submitSLCIOanalysis.py -s `pwd`/test/runSLCIOanalysis.py -o "-j SelectedKt075PF" -i /store/cmst3/user/psilva/PFcal/hbb_nunu_ILD -l `pwd`/hbb_nunu_ILD -q 8nh

python test/submitSLCIOanalysis.py -s `pwd`/test/runSLCIOanalysis.py -o "-j Kt4PF"   -i /store/cmst3/user/psilva/PFcal/pp_WZ/DST/FJ -l `pwd`/pp_WZ_AK4 -q 8nh
python test/submitSLCIOanalysis.py -s `pwd`/test/runSLCIOanalysis.py -o "-j Kt4PF"   -i /store/cmst3/user/psilva/PFcal/pp_quarkGun/DST/FJ -l `pwd`/pp_quarkGun_AK4 -q 8nh

# Merge outputs and show plots

hadd hbb_nunu_ILD.root hbb_nunu_ILD/output/*.root
rm -rf hbb_nunu_ILD
python test/studyVqqResolution.py -i hbb_nunu_ILD.root