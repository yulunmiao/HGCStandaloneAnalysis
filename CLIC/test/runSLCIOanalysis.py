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
    histos['sel'] = TH1F('sel',  ';Selection;Events',  4,0,4)
    histos['sel'].GetXaxis().SetBinLabel(1,'Total')
    histos['sel'].GetXaxis().SetBinLabel(2,'Gen V')
    histos['sel'].GetXaxis().SetBinLabel(3,'Gen V#rightarrowqq')
    histos['sel'].GetXaxis().SetBinLabel(4,'Reco V#rightarrowqq')

    histos['nvtx'] = TH1F('nvtx',  ';Vertices;Events',20,0,20)

    histos['pfsumpt'] = TH1F('pfsumpt',  ';Scalar transverse momentum sum [GeV];Events', 50,0,500)
    histos['pfmet'] = TH1F('pfmet',  ';Raw PF MET [GeV];Events', 50,0,100)
    histos['pfmetoversumpt'] = TH1F('pfmetoversumpt',  ';MET/#Sigma p_{T};Events', 50,0,0.5)

    kin=['','30to40','40to50','50to75','75to100','100toInf']
    reg=['','barrel','endcap']
    for k in kin:
        for r in reg:
            histos['jden'+k+r]  = TH1F('den'+k+r, ';#DeltaE/E;Jets',           100,-2,2)
            histos['jdpt'+k+r]  = TH1F('dpt'+k+r, ';#Deltap_{T}/p_{T};Jets',   100,-2,2)
            histos['jdphi'+k+r] = TH1F('dphi'+k+r, ';|#Delta#phi| [rad];Jets', 25,0,0.5)
            histos['jdeta'+k+r] = TH1F('deta'+k+r, ';|#Delta#eta|;Jets',       25,0,0.5)
            histos['jdr'+k+r]   = TH1F('dr'+k+r,  ';#DeltaR(jet,q);Jets',      25,0,0.5)
    
    bos=['','h','w','z']
    kin=['','50','100']
    reg=['','bb','ee','eb']
    for ib in bos:
        histos[ib+'mqq']      = TH1F(ib+'mqq',       ';Diquark mass [GeV];Events',       250,0,250)
        histos[ib+'qpt']      = TH1F(ib+'qpt',       ';Quark transverse momentum [GeV];Events',       250,0,250)
        histos[ib+'qeta']     = TH1F(ib+'qeta',     ';Quark pseudo-rapidity;Events',       50,0,3.0)
        for k in kin:
            for r in reg:
                histos[ib+'mjj'+k+r]      = TH1F(ib+'mjj'+k+r,       ';Dijet mass [GeV];Events', 400,0,200)
                histos[ib+'dmjj'+k+r]     = TH1F(ib+'dmjj'+k+r,      ';#Delta m/m;Events',  100,-2.02,1.98)

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
        print 'Event %d'%ievent
        
        genHandle=None
        jetHandle=None
        pfHandle=None
        pfClusterHandle=None
        pvHandle=None
        for collectionName, collection in event:
            if collectionName == genCollName : genHandle=collection
            if collectionName == jetCollName : jetHandle=collection
            if collectionName == "PandoraPFOs" : pfHandle=collection
            if collectionName == "PrimaryVertex" : pvHandle=collection
        if jetHandle is None : continue

        histos['nvtx'].Fill(len(pvHandle))
        
        pfsumpt=0.
        pfmet=TLorentzVector(0,0,0,0)
        for pf in pfHandle:
            pfmet = pfmet + pf.getLorentzVec()
            pfsumpt = pfsumpt+pf.getLorentzVec().Pt()
        histos['pfsumpt'].Fill(pfsumpt)
        histos['pfmet'].Fill(pfmet.Pt())
        histos['pfmetoversumpt'].Fill(pfmet.Pt()/pfsumpt)

        
        #Hard process: status 3 particles 
        #save only after at least one boson was found, otherwise incoming partons will be considered
        genBosons=[]
        promptQuarks=[]
        for p in genHandle:
            if p.getGeneratorStatus()!=3 : continue            
            if math.fabs(p.getPDG())==23 or math.fabs(p.getPDG())==24 or math.fabs(p.getPDG())==25:
                genBosons.append(p)
            if len(genBosons)==0: continue
            #if math.fabs(p.getPDG())>6 : continue
            if math.fabs(p.getPDG())>4 : continue
            promptQuarks.append( p )

        #if none found, check if from ILD official production (status is different)
        if len(genBosons)==0 :
            for p in genHandle:
                #if math.fabs(p.getPDG())!=23 and math.fabs(p.getPDG())!=24 and math.fabs(p.getPDG())!=25: continue
                if math.fabs(p.getPDG())!=25 : continue
                if p.getGeneratorStatus()!=2 : continue
                genBosons.append(p)
                for d in p.getDaughters() :
                    if math.fabs(d.getPDG())>6 : continue
                    promptQuarks.append(d)
        
        #if still not found may be the quark/gluon gun: take all quarks status 2
        if len(genBosons)==0 :
            for p in genHandle:
                if not(math.fabs(p.getPDG())<4 or math.fabs(p.getPDG())==21) : continue 
                if p.getGeneratorStatus()!=2 : continue
                histos['qpt'].Fill( p.getLorentzVec().Pt() )
                histos['qeta'].Fill( math.fabs( p.getLorentzVec().Eta() ) )
                if math.fabs( p.getLorentzVec().Eta() ) > 2.5 : continue
                jet=findRecoMatch(p,jetHandle,jetCone)[1]
                if jet is None: continue

                pt=p.getLorentzVec().Pt()
                eta=math.fabs(p.getLorentzVec().Eta())
                if pt<20 or eta>2.5 : continue
                
                jetKin=['']
                if pt<40    : jetKin.append('30to40')
                elif pt<50  : jetKin.append('40to50')
                elif pt<75  : jetKin.append('50to75')
                elif pt<100 : jetKin.append('75to100')
                else        : jetKin.append('100toInf')

                jetReg=['']
                if eta<1.5 : jetReg.append('barrel')
                else       : jetReg.append('endcap')
                
                for k in jetKin:
                    for r in jetReg:
                        histos['jden'+k+r].Fill(jet.getLorentzVec().E()/p.getLorentzVec().E()-1)
                        histos['jdpt'+k+r].Fill(jet.getLorentzVec().Pt()/p.getLorentzVec().Pt()-1)
                        histos['jdr'+k+r].Fill(p.getLorentzVec().DeltaR(jet.getLorentzVec()))
                        histos['jdeta'+k+r].Fill(math.fabs(p.getLorentzVec().Eta()-jet.getLorentzVec().Eta()))
                        histos['jdphi'+k+r].Fill(math.fabs(p.getLorentzVec().DeltaPhi(jet.getLorentzVec())))


        #match di-quarks to bosons, match- each quark to a reco jet (if available)
        bosonDecays=[]
        bosonRecoDecays=[]
        goodBosonRecoFound=False
        nPromptQuarks=len(promptQuarks)
        matchedPromptQuarksIdx=[]
        for b in genBosons:

            requireSameFlavor=( math.fabs(b.getPDG())==23 or math.fabs(b.getPDG())==25 )
            bosonCharge=b.getCharge()
            bosonMass=b.getLorentzVec().M()

            for i in xrange(0,nPromptQuarks):
                if i in matchedPromptQuarksIdx : continue
                for j in xrange(i+1,nPromptQuarks):
                    if j in matchedPromptQuarksIdx : continue

                    #di-quark
                    qqCharge=promptQuarks[i].getCharge()+promptQuarks[j].getCharge()
                    qq=promptQuarks[i].getLorentzVec()+promptQuarks[j].getLorentzVec()
                    drqq2b=qq.DeltaR(b.getLorentzVec())
                    qqMass=qq.M()
                    dm=math.fabs(qqMass-bosonMass)
                    qqFlavor=math.fabs(promptQuarks[i].getPDG()+promptQuarks[j].getPDG())

                    #total charge
                    if int(qqCharge) != int(bosonCharge) : continue

                    #check if both have the same flavor (Z/H do not mix different weak isospin components)
                    if requireSameFlavor :
                        if qqFlavor!=0 : continue

                    #match direction
                    if drqq2b>0.1 : continue

                    bosName=['']
                    if math.fabs(b.getPDG())==25 : bosName.append('h')
                    if math.fabs(b.getPDG())==23 : bosName.append('z')
                    if math.fabs(b.getPDG())==24 : bosName.append('w')
                    for ibname in bosName : histos[ibname+'mqq'].Fill(qq.M())

                    #match mass within 5 GeV
                    if dm>2.5 : continue

                    #acceptance cuts for individual quarks
                    for ibname in bosName : 
                        histos[ib+'qpt'].Fill( promptQuarks[i].getLorentzVec().Pt() )
                        histos[ib+'qpt'].Fill( promptQuarks[j].getLorentzVec().Pt() )
                        histos[ib+'qeta'].Fill( math.fabs( promptQuarks[i].getLorentzVec().Eta() ) )
                        histos[ib+'qeta'].Fill( math.fabs( promptQuarks[j].getLorentzVec().Eta() ) )
                    if promptQuarks[i].getLorentzVec().Pt()<20 : continue
                    if promptQuarks[j].getLorentzVec().Pt()<20 : continue
                    if math.fabs( promptQuarks[i].getLorentzVec().Eta() ) > 2.5 : continue
                    if math.fabs( promptQuarks[j].getLorentzVec().Eta() ) > 2.5 : continue

                    #ok to analyze this one
                    bosonDecays.append([b,promptQuarks[i],promptQuarks[j]])
                    j1=findRecoMatch(promptQuarks[i],jetHandle,jetCone)[1]
                    j2=findRecoMatch(promptQuarks[j],jetHandle,jetCone)[1]
                    if j1 is not None and j2 is not None: goodBosonRecoFound=True
                    bosonRecoDecays.append([b,j1,j2])
                    matchedPromptQuarksIdx.append(i)
                    matchedPromptQuarksIdx.append(j)
                    break

        #event selection
        histos['sel'].Fill(0)
        if len(genBosons)==0 : continue
        histos['sel'].Fill(1)
        if len(bosonDecays)==0: continue
        histos['sel'].Fill(2)
        if goodBosonRecoFound==True:  histos['sel'].Fill(3)

        #match the generated decays at reconstruction level
        for i in xrange(0,len(bosonDecays)):
            b=bosonDecays[i][0]
            q1=bosonDecays[i][1]
            q2=bosonDecays[i][2]
            j1=bosonRecoDecays[i][1]
            j2=bosonRecoDecays[i][2]
            if j1 is None: continue
            if j2 is None: continue

            dijet=j1.getLorentzVec()+j2.getLorentzVec()
            dRjj=j1.getLorentzVec().DeltaR( j2.getLorentzVec() )
            if dRjj < jetCone:
                print 'Merged jets are not considered...'
                continue


            inMwindow = (math.fabs(dijet.M()-b.getLorentzVec().M())<15)
            pt1=j1.getLorentzVec().Pt()
            eta1=math.fabs(j1.getLorentzVec().Eta())
            pt2=j2.getLorentzVec().Pt()
            eta2=math.fabs(j2.getLorentzVec().Eta())

            jjKin=['']
            if pt1>50 or pt2>50 : jjKin.append('50')
            if pt1>100 or pt2>100 : jjKin.append('100')

            jjReg='b'
            if eta1>1.5: jjReg='e'
            if eta2<=1.5: jjReg='b'+jjReg
            if eta2>1.5: jjReg='e'+jjReg
            if jjReg=='be' : jjReg='eb'
            jjRegs=['',jjReg]

            bosNames=['']
            if math.fabs(b.getPDG())==25 : bosNames.append('h')
            if math.fabs(b.getPDG())==23 : bosNames.append('z')
            if math.fabs(b.getPDG())==24 : bosNames.append('w')
                
            for ib in bosNames:
                for k in jjKin:
                    for r in jjRegs:
                        histos[ib+'mjj'+k+r].Fill(dijet.M())
                        histos[ib+'dmjj'+k+r].Fill(dijet.M()/b.getLorentzVec().M()-1)


            #individual jet resolutions
            for qj in [[q1,j1],[q2,j2]]:
                quark=qj[0]
                jet=qj[1]

                pt=quark.getLorentzVec().Pt()
                eta=math.fabs(quark.getLorentzVec().Eta())
                if pt<20 or eta>2.5 : continue
                
                jetKin=['']
                if pt<40    : jetKin.append('30to40')
                elif pt<50  : jetKin.append('40to50')
                elif pt<75  : jetKin.append('50to75')
                elif pt<100 : jetKin.append('75to100')
                else        : jetKin.append('100toInf')

                jetReg=['']
                if eta<1.5 : jetReg.append('barrel')
                else       : jetReg.append('endcap')
                

                for k in jetKin:
                    for r in jetReg:
                        histos['jden'+k+r].Fill(jet.getLorentzVec().E()/quark.getLorentzVec().E()-1)
                        histos['jdpt'+k+r].Fill(jet.getLorentzVec().Pt()/quark.getLorentzVec().Pt()-1)
                        histos['jdr'+k+r].Fill(quark.getLorentzVec().DeltaR(jet.getLorentzVec()))
                        histos['jdeta'+k+r].Fill(math.fabs(quark.getLorentzVec().Eta()-jet.getLorentzVec().Eta()))
                        histos['jdphi'+k+r].Fill(math.fabs(quark.getLorentzVec().DeltaPhi(jet.getLorentzVec())))

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
    if(fToProcess.find('/store/')==0):
        print 'Copying input locally to tmp'
        locF='/tmp/%s'%os.path.basename(fToProcess)
        os.system('cmsStage %s %s'%(fToProcess,locF))
        fToProcess=locF
        
    analyze( dstFileName=fToProcess, genCollName=opt.genColl, jetCollName=opt.jetColl, output=opt.output)
    print 'Analysis ended'

    if(fToProcess.find('/tmp/')==0): 
        print 'Removing local input'
        os.system('rm %s'%fToProcess)




