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

    #if not matchedJet is None:    
    #    print '%d (%3.1f,%3.1f,%3.1f) \t (%3.1f,%3.1f,%3.1f) \t\t %3.1f'%(
    #        genP.getPDG(),
    #        genP.getLorentzVec().Pt(),genP.getLorentzVec().Eta(),genP.getLorentzVec().Phi(),
    #        matchedJet.getLorentzVec().Pt(),matchedJet.getLorentzVec().Eta(),matchedJet.getLorentzVec().Phi(),
    #        minDR
    #        )

    if minDR>cone : matchedJet=None
    toReturn[1]=matchedJet

    return toReturn

"""
"""
def isSemiLepJet(genHandle,j,cone=0.5) :
    isSemiLep=False
    for p in genHandle:
        if p.getGeneratorStatus() != 1 : continue
        pid=math.fabs(p.getPDG())
        if pid<11 or pid>16 : continue
        dR=j.getLorentzVec().DeltaR(p.getLorentzVec())
        if dR>cone : continue
        isSemiLep=True
    return isSemiLep


def analyze( dstFileName, genCollName, jetCollName, output, minPt=0.5 ):

    jetCone=0.5
    #parse the number
    #if jetCollName=="JetOut": jetCone=0.75
    #else :
    #    jetCone = float( (re.findall(r'[0-9]+', jetCollName))[0] )/10
    #    if jetCone==7.5 : jetCone=0.75 

    print 'Using %s jet collection with matching R=%f'%(jetCollName,jetCone)

    #open root file to store the histograms
    rootFile=TFile(output,'recreate')
    histos={}
    histos['sel'] = TH1F('sel',  ';Selection;Events',  3,0,3)
    histos['sel'].GetXaxis().SetBinLabel(1,'RECO')
    histos['sel'].GetXaxis().SetBinLabel(2,'Gen X->qq')
    histos['sel'].GetXaxis().SetBinLabel(3,'Reco X->qq')
    histos['nmatches'] = TH1F('nmatches',  ';Jet-b matches found;Events',  4,0,4)

    cats=['j1','j2','j']
    kin=['','30to50','50to100','100toInf']
    reg=['','barrel','endcap']
    sel=['','inm']
    for c in cats:
        for k in kin:
            for r in reg:
                for s in sel:
                    histos[c+s+'dr'+k+r]       = TH1F(c+s+'dr'+k+r,      ';#Delta R(jet #1,b);Jets',   40,0,2)
                    histos[c+s+'dpt'+k+r]      = TH1F(c+s+'dpt'+k+r,      ';#Delta p_{T}/p_{T};Jets',   100,-2,2)
                    histos[c+s+'den'+k+r]      = TH1F(c+s+'den'+k+r,      ';#Delta E/E;Jets',   100,-2,2)
    
    kin=['','30','50','100']
    reg=['','bb','ee','eb']
    for k in kin:
        for r in reg:
            histos['mjj'+k+r]      = TH1F('mjj'+k+r,       ';Dijet mass [GeV];Events',       100,0,250)
            histos['dmjj'+k+r]     = TH1F('dmjj'+k+r,      ';#Delta m/m;Events',  50,-2,2)

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
        #print ievent
        
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
                vIsPrompt=True

                vDecaysHad=False
                iq=0
                for d in p.getDaughters():
                    #if not (math.fabs(d.getPDG())<6 or math.fabs(d.getPDG())>100): continue
                    if not (math.fabs(d.getPDG())<6): continue
                    vDecaysHad=True
                    if vIsPrompt :
                        iq=iq+1
                        j=findRecoMatch(d,jetHandle,jetCone)[1]
                        if j is None: continue

                        semiLepJet=isSemiLepJet(genHandle,j,jetCone)
                        if semiLepJet: continue

                        #save jet and compute its resolution
                        if j not in matchedJets :
                            matchedJets.append([j,d])
                            cats=['j','j%d'%iq]
                            reg=['']
                            if math.fabs(j.getLorentzVec().Eta())<1.5 : reg.append('barrel')
                            else : reg.append('endcap')
                            kin=['']
                            if j.getLorentzVec().Pt()>30 and j.getLorentzVec().Pt()<50 : kin.append('30to50')
                            if j.getLorentzVec().Pt()>50 and j.getLorentzVec().Pt()<100 : kin.append('50to100')
                            if j.getLorentzVec().Pt()>100 : kin.append('100toInf')
                            for c in cats:
                                for k in kin:
                                    for r in reg:
                                        histos[c+'dr'+k+r].Fill(j.getLorentzVec().DeltaR(d.getLorentzVec()))
                                        histos[c+'dpt'+k+r].Fill(j.getLorentzVec().Pt()/d.getLorentzVec().Pt()-1)
                                        histos[c+'den'+k+r].Fill(j.getLorentzVec().E()/d.getLorentzVec().E()-1)
                                    
                if vDecaysHad : genX.append(p)

        if len(genX)!=1 : continue
        histos['sel'].Fill(1)

        histos['nmatches'].Fill(len(matchedJets))
        if len(matchedJets)==2:

            #further selection
            j=[matchedJets[0][0],matchedJets[1][0]]
            q=[matchedJets[0][1],matchedJets[1][1]]
            if math.fabs(j[0].getLorentzVec().Eta())>3 or math.fabs(j[1].getLorentzVec().Eta())>3 : continue
            histos['sel'].Fill(2)

            dijet=TLorentzVector(j[0].getLorentzVec())
            dijet+=j[1].getLorentzVec()

            #mass resolution
            kin=['']
            if j[0].getLorentzVec().Pt()>30 and j[1].getLorentzVec().Pt()>30 : kin.append('30') 
            if j[0].getLorentzVec().Pt()>50 and j[1].getLorentzVec().Pt()>50 : kin.append('50') 
            if j[0].getLorentzVec().Pt()>100 and j[1].getLorentzVec().Pt()>100 : kin.append('100') 
            
            reg=['']
            regStr=''
            if   math.fabs(j[0].getLorentzVec().Eta())<1.5 : regStr='b'
            elif math.fabs(j[0].getLorentzVec().Eta())<3 :   regStr='e'
            if   math.fabs(j[1].getLorentzVec().Eta())<1.5 : regStr=regStr+'b'
            elif math.fabs(j[1].getLorentzVec().Eta())<3 :   regStr=regStr+'e'
            if regStr=='be' : regStr='eb'
            reg.append(regStr)

            for k in kin:
                for r in reg:
                    histos['mjj'+k+r].Fill(dijet.M())
                    histos['dmjj'+k+r].Fill(dijet.M()/genX[0].getLorentzVec().M()-1)


            #print '%d (%3.1f,%3.1f,%3.1f) \t (%3.1f,%3.1f,%3.1f) \t deltaM=%f'%(
            #    genX[0].getPDG(),
            #    genX[0].getLorentzVec().Pt(),genX[0].getLorentzVec().Eta(),genX[0].getLorentzVec().Phi(),
            #    dijet.Pt(),dijet.Eta(),dijet.Phi(),
            #    dijet.M()-genX[0].getLorentzVec().M()
            #    )

            #resolution in tight mass window
            if dijet.M()>110 and dijet.M()<140 :
                for iq in [0,1] :
                    cats=['j','j%d'%(iq+1)]
                    reg=['']
                    
                    if math.fabs(j[iq].getLorentzVec().Eta())<1.5 : reg.append('barrel')
                    else : reg.append('endcap')
                    kin=['']
                    if j[iq].getLorentzVec().Pt()>30 and j[iq].getLorentzVec().Pt()<50  : kin.append('30to50')
                    if j[iq].getLorentzVec().Pt()>50 and j[iq].getLorentzVec().Pt()<100 : kin.append('50to100')
                    if j[iq].getLorentzVec().Pt()>100                                   : kin.append('100toInf')
                    for c in cats:
                        for k in kin:
                            for r in reg:
                                histos[c+'inmdr'+k+r].Fill(j[iq].getLorentzVec().DeltaR(q[iq].getLorentzVec()))
                                histos[c+'inmdpt'+k+r].Fill(j[iq].getLorentzVec().Pt()/q[iq].getLorentzVec().Pt()-1)
                                histos[c+'inmden'+k+r].Fill(j[iq].getLorentzVec().E()/q[iq].getLorentzVec().E()-1)



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




