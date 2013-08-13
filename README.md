PFCal
=====

PF Calorimeter studies

Check https://twiki.cern.ch/twiki/bin/view/CMS/PFForwardCalorimeterStudies


Run FastJet over a SLCIO file (will create AK5PF and CA8PF)

Marlin test/fastjet.xml --MyLCIOOutputProcessor.LCIOOutputFile=qqqq_dst_2163_1_ak5pf.slcio --global.LCIOInputFiles=qqqq_dst_2163_1.slcio

Run simple analyzer over SLCIO file

python test/studyVqqResolution.py -i qqqq_dst_2163_1_fj.slcio -j AK5PF