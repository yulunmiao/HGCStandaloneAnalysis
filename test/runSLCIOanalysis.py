#! /usr/bin/env python

from pyLCIO import IOIMPL
from ROOT import TFile, TH1F, TH2F, TLorentzVector
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

    print '%d (%3.1f,%3.1f,%3.1f) \t (%3.1f,%3.1f,%3.1f) \t\t %3.1f'%(
        genP.getPDG(),
        genP.getLorentzVec().Pt(),genP.getLorentzVec().Eta(),genP.getLorentzVec().Phi(),
        matchedJet.getLorentzVec().Pt(),matchedJet.getLorentzVec().Eta(),matchedJet.getLorentzVec().Phi(),
        minDR
        )

    if minDR>cone : matchedJet=None
    toReturn[1]=matchedJet

    return toReturn


def findAllRecoMatch(genPcoll,jetHandle,cone=0.5):

    toReturn=[]
    if jetHandle is None or genPcoll is None: return toReturn
 
    for genP in genPcoll:
        j=findRecoMatch(genP,jetHandle,cone)[1]
        if j is None: continue

        isAdded=False
        for mj in toReturn:
            dR=j.getLorentzVec().DeltaR(mj.getLorentzVec())
            if dR<cone: continue
            isAdded=True
        if isAdded : continue

        toReturn.append(j)

    return toReturn


def analyze( dstFileName, genCollName, jetCollName, output, minPt=0.5 ):

    #parse the number
    if jetCollName=="JetOut": jetCone=0.75
    else :
        jetCone = float( (re.findall(r'[0-9]+', jetCollName))[0] )/10
        if jetCone==7.5 : jetCone=0.75 

    print 'Using %s jet collection with matching R=%f'%(jetCollName,jetCone)

    #open root file to store the histograms
    rootFile=TFile(output,'recreate')
    histos={}
    histos['sel'] = TH1F('sel',  ';Selection;Events',  3,0,3)
    histos['sel'].GetXaxis().SetBinLabel(1,'RECO')
    histos['sel'].GetXaxis().SetBinLabel(2,'Gen X->qq')
    histos['sel'].GetXaxis().SetBinLabel(3,'Reco X->qq')
    cats=['j1','j2']
    for c in cats:
        histos[c+'dr']       = TH1F(c+'dr',      ';#Delta R(jet #1,b);Jets',   40,0,4)
        histos[c+'dpt']      = TH1F(c+'dpt',      ';#Delta p_{T}/p_{T};Jets',   50,-2,2)
    histos['nmatches'] = TH1F('nmatches',  ';Jet-b matches found;Events',  4,0,4)
    histos['mjj']      = TH1F('mjj',       ';Dijet mass [GeV];Events',       100,0,250)
    histos['dmjj']     = TH1F('dmjj',      ';#Delta m/m;Events',  50,-2,2)

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
    #print 'Generated \t\t Reconstructed \t\t\t Delta'
    ievent=0
    for event in lcReader:
        ievent=ievent+1
        print ievent
        genHandle=None
        jetHandle=None
        for collectionName, collection in event:
            #print collectionName
            if collectionName == genCollName : genHandle=collection
            if collectionName == jetCollName : jetHandle=collection
            
        if jetHandle is None : continue

        histos['sel'].Fill(0)

        #search for prompt X->qq'
        genX=[]
        matchedJets=[]
        for p in genHandle:
            if math.fabs(p.getPDG())==25  :
            #if math.fabs(p.getPDG())==25 or math.fabs(p.getPDG())==23 or math.fabs(p.getPDG())==24 :
            #if math.fabs(p.getPDG())==25 or math.fabs(p.getPDG())==24 :

                #vIsPrompt=False
                #for m in p.getParents() :
                #    if not (math.fabs(m.getPDG())==11 or math.fabs(m.getPDG())==25): continue
                #    vIsPrompt=True
                #if len(p.getParents())==0 : vIsPrompt=True
                vIsPrompt=True

                vDecaysHad=False
                iq=0
                for d in p.getDaughters():
                    if not (math.fabs(d.getPDG())<6 or math.fabs(d.getPDG())>100): continue
                    vDecaysHad=True
                    if vIsPrompt :
                        iq=iq+1
                        j=findRecoMatch(d,jetHandle,jetCone)[1]
                        if j is None: continue
                        if j not in matchedJets :
                            matchedJets.append(j)
                            c='j%d'%iq
                            histos[c+'dr'].Fill(j.getLorentzVec().DeltaR(d.getLorentzVec()))
                            histos[c+'dpt'].Fill(j.getLorentzVec().Pt()/d.getLorentzVec().Pt()-1)

                if vDecaysHad : genX.append(p)

        if len(genX)!=1 : continue
        histos['sel'].Fill(1)

        histos['nmatches'].Fill(len(matchedJets))
        if len(matchedJets)==2:
            histos['sel'].Fill(2)
            dijet=TLorentzVector(matchedJets[0].getLorentzVec())
            dijet+=matchedJets[1].getLorentzVec()
            
            print '%d (%3.1f,%3.1f,%3.1f) \t (%3.1f,%3.1f,%3.1f) \t deltaM=%f'%(
                genX[0].getPDG(),
                genX[0].getLorentzVec().Pt(),genX[0].getLorentzVec().Eta(),genX[0].getLorentzVec().Phi(),
                dijet.Pt(),dijet.Eta(),dijet.Phi(),
                dijet.M()-genX[0].getLorentzVec().M()
                )

            histos['mjj'].Fill(dijet.M())
            histos['dmjj'].Fill(dijet.M()/genX[0].getLorentzVec().M()-1)


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




