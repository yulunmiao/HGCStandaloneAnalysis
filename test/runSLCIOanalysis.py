#! /usr/bin/env python

from pyLCIO import IOIMPL
from ROOT import TFile, TH1F, TH2F, TProfile
import optparse
import math
import re
import os

def findRecoMatch(q,jetHandle,cone=0.5):

    if jetHandle is None: return None
    if q.getLorentzVec().Pt()==0 : return None

    minDR=9999.
    matchedJet=-1
    ijet=-1
    for jet in jetHandle :
        ijet=ijet+1
        if jet.getLorentzVec().Pt()==0 : continue
        dR=jet.getLorentzVec().DeltaR(q.getLorentzVec())
        if dR>minDR : continue
        minDR=dR
        matchedJet=ijet
    if minDR>cone : matchedJet=-1
    return matchedJet


def analyze( dstFileName, genCollName, jetCollName, output, minPt=0.5 ):

    #parse the number
    jetCone = float( (re.findall(r'[0-9]+', jetCollName))[0] )/10

    #open root file to store the histograms
    rootFile=TFile(output,'recreate')
    histos={}
    cands=['j','q']
    for c in cands:
        histos['n'+c]=TH1F('n'+c,';Multiplicity;Events',6,0,6)
        histos[c+'pt']=TH1F(c+'pt',';Transverse momentum [GeV];Jets/quarks',50,0,500)
        histos[c+'mass']=TH1F(c+'mass',';Mass [GeV];Jets/quarks',50,0,250)
        histos[c+'eta']=TH1F(c+'eta',';Pseudo-rapidity;Jets/quarks',50,0,5)
        histos[c+'mindr']=TH1F(c+'mindr',';min #Delta R;Jets/quarks',20,0,6)
    histos['nq'].SetFillStyle(3002)
    histos['nq'].SetFillColor(32)
    histos['qpt'].SetFillStyle(3002)
    histos['qpt'].SetFillColor(32)
    histos['qeta'].SetFillStyle(3002)
    histos['qeta'].SetFillColor(32)
    histos['qmass'].SetFillStyle(3002)
    histos['qmass'].SetFillColor(32)
    histos['dptvspt']=TProfile('dptvspt',';Generated transverse momentum [GeV];Response;Jets',50,0,500,0,10)
    histos['dptvseta']=TProfile('dptvseta',';Generated pseudo-rapidity;Response;Jets',50,0,5,0,10)
    for key in histos:
        histos[key].SetMarkerStyle(20)
        histos[key].SetMarkerColor(1)
        histos[key].SetLineColor(1)

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


        #gen particles
        resolvedQuarks=[]
        matchedJets={}
        if genHandle is None:
            for j in jetHandle : matchedJets[j]=[j]
        else:
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
                resolvedQuarks.append(p)

                #find match and store
                ijet=findRecoMatch(p,jetHandle,jetCone)
                if ijet<0: continue
                try:
                    matchedJets[ijet].append(p)
                except:
                    matchedJets[ijet]=[p]

        #check selected quarks (if DR<jetCone consider the sum of the two)
        selQuarks=[]
        for iq in xrange(0,len(resolvedQuarks)):

            mergedQuark=resolvedQuarks[iq].getLorentzVec()
            for jq in xrange(iq+1,len(resolvedQuarks)):
                if iq>=jq : continue
                dR=resolvedQuarks[iq].getLorentzVec().DeltaR(resolvedQuarks[iq].getLorentzVec())
                if dR>jetCone : continue
                mergedQuark+=resolvedQuarks[jq].getLorentzVec()
            selQuarks.append( mergedQuark )

        nquarks=0
        for q in selQuarks:
            if q.Pt()<30 : continue
            nquarks=nquarks+1
            histos['qpt'].Fill(q.Pt())            
            histos['qmass'].Fill(q.M())            
            histos['qeta'].Fill(math.fabs(q.Eta()))
        if nquarks==0 : continue
        histos['nq'].Fill(nquarks)
                    
        #kinematics of the matched jets
        njets=0
        for j in matchedJets :

            jet=jetHandle.at(j).getLorentzVec()
            if jet.Pt()<30 : continue
            
            gen=TLorentzVector(0,0,0,0)
            for p in matchedJets[j]:
                gen+=p.getLorentzVec()
            #if gen.Pt()/jet.Pt()<0.5 : continue

            njets=njets+1
            histos['jpt'].Fill(jet.Pt())            
            histos['jmass'].Fill(jet.M())            
            histos['jeta'].Fill(math.fabs(jet.Eta()))
            histos['dptvspt'].Fill(gen.Pt(),jet.Pt()/gen.Pt())
            histos['dptvseta'].Fill(math.fabs(gen.Eta()),jet.Pt()/gen.Pt())

        histos['nj'].Fill(njets)

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




