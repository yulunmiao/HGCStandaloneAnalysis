# =====
# PFCal
# =====
#
# Check https://twiki.cern.ch/twiki/bin/view/CMS/PFForwardCalorimeterStudies
#

# Run FastJet over a set of SLCIO files in castor

python test/submitMarlinRun.py -i /castor/cern.ch/grid/ilc/prod/clic/1.4tev/h_nunu/ILD/DST/00001504/000 -o /store/cmst3/user/psilva/PFcal/h_nunu

# Run simple analyzer over SLCIO files in EOS

python test/submitSLCIOanalysis.py -s `pwd`/test/runSLCIOanalysis.py -o CA8 -i /store/cmst3/user/psilva/PFcal/h_nunu -l `pwd`/hnunu_CA8

# Merge outputs and show plots

hadd hnunu_CA8.root hnunu_CA8/*.root
rm -rf hnunu_CA8/*.root
python test/studyVqqResolution.py -i hnunu_CA8.root