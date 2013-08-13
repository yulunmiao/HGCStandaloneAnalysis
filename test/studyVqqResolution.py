#! /usr/bin/env python

from pyLCIO import IOIMPL
from ROOT import TFile, TH1F, TH2F,TCanvas, TProfile, TPaveText, gStyle
import optparse
import math
import re

def findRecoMatch(q,jetHandle,cone=0.5):

    if jetHandle is None: return None
    if q.getLorentzVec().Pt()==0 : return None
    
    for jet in jetHandle :
        if jet.getLorentzVec().Pt()==0 : continue
        dR=jet.getLorentzVec().DeltaR(q.getLorentzVec())
        if dR<cone : return jet 

    return None


def analyze( dstFileName, pfCollName, genCollName, jetCollName, minPt=0.5 ):

    #parse the number
    jetCone = float( (re.findall(r'[0-9]+', jetCollName))[0] )/10

    #open root file to store the histograms
    rootFileName=dstFileName.replace('.slcio','.root')
    rootFile=TFile(rootFileName,'recreate')
    histos={}
    histos['njets']=TH1F('njets',';Jet multiplicity;Events',6,0,6)
    histos['nmergedjets']=histos['njets'].Clone('nmergedjets')
    histos['jpt']=TH1F('jpt',';Transverse momentum [GeV];Jets',50,0,500)
    histos['jeta']=TH1F('jeta',';Pseudo-rapidity;Jets',50,0,5)
    histos['dptvspt']=TProfile('dptvspt',';Generated transverse momentum [GeV];response;Jets',50,0,500,0,10)
    histos['dptvseta']=TProfile('dptvseta',';Generated pseudo-rapidity;Response;Jets',50,0,5,0,10)

#    histos['mjj']=TH1F('mjj',';M_{jj} [GeV];Jet pairs',100,0,200)
#    histos['wmjj']=histos['mjj'].Clone('wmjj')
#    histos['zmjj']=histos['mjj'].Clone('zmjj')



    # create the LCReader and open the input file
    lcReader = IOIMPL.LCFactory.getInstance().createLCReader()
    lcReader.open( dstFileName )
    
    # filling the tree
    for event in lcReader:

        pfHandle=None
        genHandle=None
        jetHandle=None
        for collectionName, collection in event:
            if collectionName == pfCollName : pfHandle=collection
            if collectionName == genCollName : genHandle=collection
            if collectionName == jetCollName : jetHandle=collection

        #gen particles
        matchedJets={}
        if genHandle is not None:
            for p in genHandle:

                #only quarks
                if math.fabs(p.getPDG())>6: continue

                #only quarks from hard process or W/Z decays
                isPrompt=False
                for m in p.getParents() :
                    if math.fabs(m.getPDG())==11 or math.fabs(m.getPDG())==23 or math.fabs(m.getPDG())==24: isPrompt=True
                if not isPrompt: continue
                
                #require minimum pT
                if p.getLorentzVec().Pt()<minPt : continue

                #find match and store
                jet=findRecoMatch(p,jetHandle,jetCone)
                if jet is None: continue
                if jet in matchedJets : matchedJets[jet].append([p])
                else : matchedJets[jet]=[p]
        
        njets=0
        nmergedJets=0
        for j in matchedJets :

            jet=j.getLorentzVec()
            if jet.Pt()<30 : continue

            nmatches=0
            for p in matchedJets[j]:
                gen=p.getLorentzVec()
                nmatches=nmatches+1
            if nmatches>1 : nmergedJets=nmergedJets+1

            njets=njets+1
            histos['jpt'].Fill(jet.Pt())            
            histos['jeta'].Fill(math.fabs(jet.Eta()))
            histos['dptvspt'].Fill(gen.Pt(),jet.Pt()/gen.Pt())
            histos['dptvseta'].Fill(math.fabs(gen.Eta()),jet.Pt()/gen.Pt())
        histos['njets'].Fill(njets)
        histos['nmergedjets'].Fill(nmergedJets)

#        #match possible V->qq' decays
#            #for p in genOfInterest:
#             #       if math.fabs(p.getPDG())!=math.fabs(mcp.getPDG()): continue
#                    if p.getPDG()==mcp.getPDG() : continue
#                    qqSystem=p.getLorentzVec()+mcp.getLorentzVec()
#
#                    deltaMZ=math.fabs(qqSystem.M()-91.1876)
#                    deltaMW=math.fabs(qqSystem.M()-80.385)
#                    if(deltaMZ<deltaMW) :
#                        if(deltaMZ<100) : 
#                            genZqq.append([p,mcp])
#                            jjSystem=findRecoMatch(p,mcp,jetHandle)
#                            if len(jjSystem)==1 : 
#                                recoZjj.append(jjSystem)
#                                print '%f %f'%(qqSystem.Pt(),jjSystem[0].getLorentzVec().Pt())
#                    else:
#                        if(deltaMW<100) : 
#                            genWqq.append([p,mcp])
#                            jjSystem=findRecoMatch(p,mcp,jetHandle)
#                            if len(jjSystem)==1 : 
#                                recoWjj.append(jjSystem)
#                                print '%f %f'%(qqSystem.Pt(),jjSystem[0].getLorentzVec().Pt())
#
#        for jj in recoZjj:
#            jjP4=jj[0].getLorentzVec()
#            histos['mjj'].Fill(jjP4.M())
#            histos['zmjj'].Fill(jjP4.M())
#        for jj in recoWjj:
#            jjP4=jj[0].getLorentzVec()
#            histos['mjj'].Fill(jjP4.M())
#            histos['wmjj'].Fill(jjP4.M())
#
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)

    c=TCanvas('c','c',1000,500)
    c.Divide(2,2)
    c.cd(1)
    histos['njets'].Draw('hist')
    #histos['nmergedjets'].SetFillStyle(1001)
    #histos['nmergedjets'].SetFillColor(34)
    #histos['nmergedjets'].Draw('histsame')
    pt=TPaveText(0.08,0.9,0.5,0.99,'brNDC')
    pt.SetBorderSize(0)
    pt.SetFillStyle(0)
    pt.AddText(jetCollName+ ' jets simulation')
    pt.Draw()
    p=c.cd(3)
    p.Divide(2,1)
    p.cd(1)
    histos['jpt'].Draw('hist')
    p.cd(2)
    histos['jeta'].Draw('hist')
    c.cd(2)
    histos['dptvspt'].Draw('e1')
    c.cd(4)
    histos['dptvseta'].Draw('e1')
    c.SaveAs('jetkin.png')

    # write and close the file
    rootFile.Write()
    rootFile.Close()    


if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',    dest='input'  , help='input'                , default='../Test/bbudsc_3evt_DST.slcio')
    parser.add_option('-g', '--gen'  ,    dest='genColl', help='gen level collection' , default='MCParticlesSkimmed')
    parser.add_option('-p', '--pf'   ,    dest='pfColl' , help='pf level collection'  , default='PandoraPFANewPFOs')
    parser.add_option('-j', '--jet'  ,    dest='jetColl', help='pf jet collection'    , default='AK5PF')
    (opt,args)=parser.parse_args()

    analyze( dstFileName=opt.input, pfCollName=opt.pfColl, genCollName=opt.genColl, jetCollName=opt.jetColl)




