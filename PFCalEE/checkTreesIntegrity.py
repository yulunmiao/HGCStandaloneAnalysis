import ROOT
import os, sys
import glob

path=sys.argv[1]
step=int(sys.argv[2])

d = {1: ['HGcal_', 'HGCSSTree'],
     2: ['Digi_', 'RecoTree']}
listing = glob.glob(os.path.join(path, d[step][0]+'*'))
for f in listing:
    fIn=ROOT.TFile.Open(os.path.join(path,f))
    nhits=fIn.Get(d[step][1]).GetEntriesFast()
    if nhits==0:
        print f
    fIn.Close()
