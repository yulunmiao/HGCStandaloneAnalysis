# =====
# PFCal
# =====
#
# Check https://twiki.cern.ch/twiki/bin/view/CMS/PFForwardCalorimeterStudies
#

# Run FastJet over a set of SLCIO files in castor

python test/submitMarlinRun.py -i /castor/cern.ch/grid/ilc/prod/clic/1.4tev/h_nunu/ILD/DST/00001504/000 -o /store/cmst3/user/psilva/PFcal/h_nunu

# Run simple analyzer over SLCIO files in EOS

python test/submitSLCIOanalysis.py -s `pwd`/test/runSLCIOanalysis.py -o "-j CA8PF" -i /store/cmst3/user/psilva/PFcal/qqqq -l `pwd`/qqqq_CA8PF -q 8nh

# Merge outputs and show plots

hadd qqqq_CA8PF.root qqqq_CA8PF/output/*.root
rm -rf qqqq_CA8PF/output/*.root
python test/studyVqqResolution.py -i qqqq_CA8PF.root