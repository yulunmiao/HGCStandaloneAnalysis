#! /usr/bin/env python

from pyLCIO import IOIMPL
from ROOT import TFile, TH1F, TH2F
import optparse
import math
import re
import os

def findRecoMatch(genP,jetHandle,cone=0.5):

    toReturn=[genP,None]
    if jetHandle is None: return toReturn
    if genP.getLorentzVec().Pt()==0 : return toReturn

    minDR=9999.
    matchedJet=None
    for jet in jetHandle :
        if jet.getLorentzVec().Pt()==0 : continue
        dR=jet.getLorentzVec().DeltaR(genP.getLorentzVec())
        if dR>minDR : continue
        minDR=dR
        matchedJet=jet
    if minDR>cone : matchedJet=None
    toReturn[1]=matchedJet

    return toReturn


def analyze( dstFileName, genCollName, jetCollName, output, minPt=0.5 ):

    #parse the number
    jetCone = float( (re.findall(r'[0-9]+', jetCollName))[0] )/10
    print 'Using %s jet collection with matching R=%f'%(jetCollName,jetCone)

    #open root file to store the histograms
    rootFile=TFile(output,'recreate')
    histos={}
    cats=['j','mj','v']
    for c in cats:

        histos[c+'n']        = TH1F(c+'n',       ';Multiplicity;Events',                     5,0,5)
        histos[c+'pt']       = TH1F(c+'pt',      ';Transverse momentum [GeV];Jets/quarks',   100,0,500)
        histos[c+'eta']      = TH1F(c+'eta',     ';Pseudo-rapidity;Jets/quarks',             50,0,5)
        histos[c+'mass']     = TH1F(c+'mass',    ';Mass [GeV];Jets/quarks',                  200,0,500)
        histos[c+'ptvsmass'] = TH2F(c+'ptvsmass',';Transverse momentum [GeV];Mass [GeV]',    100,0,500,200,0,500)

        if c!='mj' : continue
        histos[c+'dptvspt']  = TH2F(c+'dptvspt', ';Generated transverse momentum [GeV];p_{T} response;Jets', 50,0,500,100,0,5)
        histos[c+'dmvspt']   = TH2F(c+'dmtvspt', ';Generated transverse momentum [GeV];Mass response;Jets',  50,0,500,100,0,5)
        histos[c+'dptvseta'] = TH2F(c+'dptvseta', ';Generated pseudo-rapidity;p_{T} response;Jets',          50,0,5,  100,0,5)        
        histos[c+'dmvseta']  = TH2F(c+'dmvseta', ';Generated pseudo-rapidity;Mass response;Jets',            50,0,5,  100,0,5)        

    for key in histos:
        histos[key].Sumw2()
        histos[key].SetMarkerStyle(20)
        histos[key].SetMarkerColor(1)
        histos[key].SetLineColor(1)
        if key[0]=='v':
            histos[key].SetFillStyle(3002)
            histos[key].SetFillColor(32)

    # create the LCReader and open the input file
    lcReader = IOIMPL.LCFactory.getInstance().createLCReader()
    lcReader.open(dstFileName)
    
    # filling the tree
    for event in lcReader:

        genHandle=None
        jetHandle=None
        for collectionName, collection in event:
            if collectionName == genCollName : genHandle=collection
            if collectionName == jetCollName : jetHandle=collection
 
        #search for prompt V->qq'
        vqq=[]
        for p in genHandle:
            if math.fabs(p.getPDG())==23 or math.fabs(p.getPDG())==24 :

                decStr=' -> %s -> '%p.getPDG()

                vDecaysHad=False
                for d in p.getDaughters():
                    decStr=decStr+' %d '%d.getPDG()
                    if not (math.fabs(d.getPDG())<6 or math.fabs(d.getPDG())>100): continue
                    vDecaysHad=True

                vIsPrompt=False
                for m in p.getParents() :
                    decStr = ' %d '%m.getPDG() + decStr
                    if not (math.fabs(d.getPDG())==11 or math.fabs(d.getPDG())==25): continue
                    vIsPrompt=True
                if len(p.getParents())==0 : vIsPrompt=True

                print '%s prompt=%s decaysHad=%s'%(decStr,vIsPrompt,vDecaysHad)

                if not (vIsPrompt and vDecaysHad): continue
                vqq.append(p)

        hasVqq=(len(vqq)!=0)


        #kinematics of the reconstructed jets
        selJets=[]
        for j in jetHandle:
            if j.getLorentzVec().Pt()<30 or math.fabs(j.getLorentzVec().Eta())>3: continue

            tags=['j']
            if hasVqq : tags.append('mj')
            for t in tags :
                histos[t+'pt'].Fill(j.getLorentzVec().Pt())
                histos[t+'eta'].Fill(math.fabs(j.getLorentzVec().Eta()))
                histos[t+'mass'].Fill(j.getLorentzVec().M())
                histos[t+'ptvsmass'].Fill(j.getLorentzVec().Pt(),j.getLorentzVec().M())

            selJets.append(j)
        histos['jn'].Fill(len(selJets))
        if hasVqq : histos['mjn'].Fill(len(selJets))

        #now consider only matched events with V->qq'
        if not hasVqq: continue

        #fill gen level plots
        histos['vn'].Fill(len(vqq))
        for v in vqq:
            histos['vpt'].Fill(v.getLorentzVec().Pt())            
            histos['veta'].Fill(math.fabs(v.getLorentzVec().Eta()))
            histos['vmass'].Fill(v.getLorentzVec().M())            
            histos['vptvsmass'].Fill(v.getLorentzVec().Pt(),v.getLorentzVec().M())            

        #resolution for matched jets
        matchedJets=[ findRecoMatch(v,jetHandle,jetCone) for v in vqq ]
        for m in matchedJets :
            v=m[0]
            j=m[1]
            if j is None or v is None : continue
            dpt=j.getLorentzVec().Pt()/v.getLorentzVec().Pt()
            dm=j.getLorentzVec().M()/v.getLorentzVec().M()
            histos['mjdptvspt'].Fill(v.getLorentzVec().Pt(),dpt)
            histos['mjdmvspt'].Fill(v.getLorentzVec().Pt(),dm)
            histos['mjdptvseta'].Fill(math.fabs(v.getLorentzVec().Eta()),dpt)
            histos['mjdmvseta'].Fill(math.fabs(v.getLorentzVec().Eta()),dm)
            

    # write and close the file
    rootFile.Write()
    rootFile.Close()    


if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',    dest='input'  , help='input'                , default='../Test/bbudsc_3evt_DST.slcio')
    parser.add_option('-o', '--output',   dest='output',  help='output'               , default='analysis.root')
    parser.add_option('-g', '--gen'  ,    dest='genColl', help='gen level collection' , default='MCParticlesSkimmed')
    parser.add_option('-j', '--jet'  ,    dest='jetColl', help='pf jet collection'    , default='AK5PF')
    (opt,args)=parser.parse_args()

    fToProcess=opt.input
    if(fToProcess.find('cmst3')>0):
        print 'Copying input locally to tmp'
        locF='/tmp/%s'%os.path.basename(fToProcess)
        os.system('cmsStage %s %s'%(fToProcess,locF))
        fToProcess=locF
        
    analyze( dstFileName=fToProcess, genCollName=opt.genColl, jetCollName=opt.jetColl, output=opt.output)
    print 'Analysis ended'

    if(fToProcess.find('/tmp/')==0): 
        print 'Removing local input'
        os.system('rm %s'%fToProcess)




