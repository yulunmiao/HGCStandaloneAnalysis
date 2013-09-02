# =====
# PFCal
# =====
#
# Check https://twiki.cern.ch/twiki/bin/view/CMS/PFForwardCalorimeterStudies
#

# Run FastJet over a set of SLCIO files in castor

python test/submitMarlinRun.py -i /castor/cern.ch/grid/ilc/prod/clic/1.4tev/h_nunu/SID/DST/00002466/000 -o /store/cmst3/user/psilva/PFcal/hbb_nunu_SID

python test/submitMarlinRun.py -i /castor/cern.ch/grid/ilc/prod/clic/1.4tev/h_nunu/ILD/DST/00002022/000 -o /store/cmst3/user/psilva/PFcal/hbb_nunu_ILD -s test/fastjet_newpfos.xml -l ${PWD}/FARMILD

# Run simple analyzer over SLCIO files in EOS

python test/submitSLCIOanalysis.py -s `pwd`/test/runSLCIOanalysis.py -o "-j SelectedKt075PF" -i /store/cmst3/user/psilva/PFcal/hbb_nunu_ILD -l `pwd`/hbb_nunu_ILD -q 8nh

# Merge outputs and show plots

hadd hbb_nunu_ILD.root hbb_nunu_ILD/output/*.root
rm -rf hbb_nunu_ILD
python test/studyVqqResolution.py -i hbb_nunu_ILD.root