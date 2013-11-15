#! /usr/bin/env python

from ROOT import TFile, TF1, TH1F, TH2F,TCanvas, TProfile, TPaveText, TLorentzVector, TLegend,TGraphErrors, gStyle, gROOT
from ROOT import RooArgSet, RooRealVar,RooGaussian, RooDataHist, RooAbsData, RooArgList, RooCBShape, RooAddPdf, RooAbsPdf, RooFit
import optparse
import math

allHeaders=[]

def drawHeader() :
    allHeaders.append(TPaveText(0.08,0.92,0.6,0.98,'brNDC'))
    ih=len(allHeaders)-1
    allHeaders[ih].SetBorderSize(0)
    allHeaders[ih].SetFillStyle(0)
    allHeaders[ih].SetTextAlign(12)
    allHeaders[ih].AddText('ILD simulation, pp collisions @ 14 TeV')
    allHeaders[ih].Draw()



def studyVqqResolution(rootFile):

    #get all from file
    histos={}
    inF=TFile.Open(rootFile)
    keys=inF.GetListOfKeys()
    for k in keys:
        obj=inF.Get(k.GetName())
        obj.SetDirectory(0)
        histos[k.GetName()]=obj
    inF.Close()

    #plot
    gROOT.SetBatch()
    gROOT.SetStyle('Plain')
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(1111)
    gStyle.SetOptTitle(0)
    gStyle.SetStatFont(42)

    kin=['','30to40','40to50','50to75','75to100','100toInf']
    for k in kin:        
        c=TCanvas('c','c',600,600)
        c.cd()
        c.SetCanvasSize(1000,500)
        c.SetWindowSize(1000,500)
        c.Divide(2,1)
        c.cd(1)
        histos['deta'+k+'barrel'].SetLineWidth(2)
        histos['deta'+k+'barrel'].SetTitle('barrel')
        histos['deta'+k+'barrel'].Draw('hist')
        histos['deta'+k+'endcap'].SetLineWidth(2)
        histos['deta'+k+'endcap'].SetLineStyle(7)
        histos['deta'+k+'endcap'].SetTitle('endcap')
        histos['deta'+k+'endcap'].Draw('histsame')
        leg=TLegend(0.6,0.92,0.9,0.98)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.AddEntry(histos['deta'+k+'barrel'],'barrel','f')
        leg.AddEntry(histos['deta'+k+'endcap'],'endcap','f')
        leg.SetNColumns(2)
        leg.Draw()
        drawHeader()
        c.cd(2)
        histos['dphi'+k+'barrel'].SetLineWidth(2)
        histos['dphi'+k+'barrel'].SetTitle('barrel')
        histos['dphi'+k+'barrel'].Draw('hist')
        histos['dphi'+k+'endcap'].SetLineWidth(2)
        histos['dphi'+k+'endcap'].SetLineStyle(7)
        histos['dphi'+k+'endcap'].SetTitle('endcap')
        histos['dphi'+k+'endcap'].Draw('histsame')
        c.Modified()
        c.Update()
        c.SaveAs('dr_%s.png'%k)


    labels=[]
    responseVars=['dpt','den','dphi','deta','dr']
    for r in responseVars:
        barrelResponse=TGraphErrors()
        barrelResponse.SetName(r+'barrelresponse')
        barrelResponse.SetLineWidth(2)
        barrelResponse.SetFillStyle(0)
        barrelResponse.SetMarkerStyle(20)
        barrelCoreResponse=barrelResponse.Clone(r+'barrelcoreresponse')
        endcapResponse=TGraphErrors()
        endcapResponse.SetName(r+'endcapresponse')
        endcapResponse.SetLineWidth(2)
        endcapResponse.SetFillStyle(0)
        endcapResponse.SetMarkerStyle(24)
        endcapCoreResponse=endcapResponse.Clone(r+'endcapresponse')
        for k in kin: 
            c.cd()
            c.Clear()
            c.SetWindowSize(1000,500)
            c.Divide(2,1)
            for i in [1,2] :
                c.cd(i)
                reg='barrel'
                if i==2: reg='endcap' 

                h=histos[r+k+reg]
                x=RooRealVar("x", h.GetXaxis().GetTitle(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
                data=RooDataHist("data", "dataset with x", RooArgList(x), h)
                frame=x.frame()
                RooAbsData.plotOn( data, frame, RooFit.DataError(RooAbsData.SumW2) )

                mean1=RooRealVar("mean1","mean1",0,-0.5,0.5);
                sigma1=RooRealVar("sigma1","sigma1",0.1,0.01,1.0);
                gauss1=RooGaussian("g1","g",x,mean1,sigma1)
                
                if r=='dpt' or r=='den' :
                    mean2=RooRealVar("mean2","mean2",0,-0.5,0.5);
                    sigma2=RooRealVar("sigma2","sigma2",0.1,0.01,1.0);
                    alphacb=RooRealVar("alphacb","alphacb",1,0.1,3);
                    ncb=RooRealVar("ncb","ncb",4,1,100)
                    gauss2 = RooCBShape("cb2","cb",x,mean2,sigma2,alphacb,ncb);
                else:
                    mean1.setRange(0,0.5)
                    mean2=RooRealVar("mean2","mean",0,0,1);
                    sigma2=RooRealVar("sigma2","sigma",0.1,0.01,1.0);
                    gauss2=RooGaussian("g2","g",x,mean2,sigma2) ;

                frac = RooRealVar("frac","fraction",0.9,0.0,1.0)
                if data.sumEntries()<100 :
                    frac.setVal(1.0)
                    frac.setConstant(True)
                model = RooAddPdf("sum","g1+g2",RooArgList(gauss1,gauss2), RooArgList(frac))

                status=model.fitTo(data,RooFit.Save()).status()
                if status!=0 : continue

                model_cdf=model.createCdf(RooArgSet(x)) ;
                cl=0.90
                ul=0.5*(1.0+cl)
                closestToCL=1.0
                closestToUL=-1
                closestToMedianCL=1.0
                closestToMedian=-1
                for ibin in xrange(1,h.GetXaxis().GetNbins()*10):
                    xval=h.GetXaxis().GetXmin()+(ibin-1)*h.GetXaxis().GetBinWidth(ibin)/10.
                    x.setVal(xval)
                    cdfValToCL=math.fabs(model_cdf.getVal()-ul)
                    if cdfValToCL<closestToCL:
                        closestToCL=cdfValToCL
                        closestToUL=xval
                    cdfValToCL=math.fabs(model_cdf.getVal()-0.5)
                    if cdfValToCL<closestToMedianCL:
                        closestToMedianCL=cdfValToCL
                        closestToMedian=xval

                RooAbsPdf.plotOn(model,frame)
                frame.Draw()

                if i==1: drawHeader()
                labels.append( TPaveText(0.6,0.92,0.9,0.98,'brNDC') )
                ilab=len(labels)-1
                labels[ilab].SetName(r+k+'txt')
                labels[ilab].SetBorderSize(0)
                labels[ilab].SetFillStyle(0)
                labels[ilab].SetTextFont(42)
                labels[ilab].SetTextAlign(12)
                kinReg=k.replace('to','-')
                kinReg=kinReg.replace('Inf','#infty')
                labels[ilab].AddText('['+reg+'] '+kinReg)
                labels[ilab].Draw()
                
                resolutionVal=math.fabs(closestToUL-closestToMedian)
                responseGr=barrelResponse
                responseCoreGr=barrelCoreResponse
                coreResolutionVal=sigma1.getVal()
                coreResolutionErr=sigma1.getError()
                if frac.getVal()<0.7 and (sigma2.getVal()<sigma1.getVal()) :
                    coreResolutionVal=sigma2.getVal()
                    coreResolutionErr=sigma2.getError()


                if i==2 : 
                    responseGr=endcapResponse
                    responseCoreGr=endcapCoreResponse
                if k!='' :
                    nrespPts=responseGr.GetN()
                    kinAvg=150
                    kinWidth=50
                    if k=='30to40' : 
                        kinAvg=35
                        kinWidth=5
                    if k=='40to50' : 
                        kinAvg=45
                        kinWidth=5
                    if k=='50to75' : 
                        kinAvg=62.5
                        kinWidth=12.5
                    elif k=='75to100' :
                        kinAvg=87.5
                        kinWidth=12.5
                    responseGr.SetPoint(nrespPts,kinAvg,resolutionVal)
                    responseCoreGr.SetPoint(nrespPts,kinAvg,coreResolutionVal)
                    responseCoreGr.SetPointError(nrespPts,kinWidth,coreResolutionErr)

                labels.append( TPaveText(0.15,0.7,0.4,0.9,'brNDC') )
                ilab=len(labels)-1
                labels[ilab].SetName(r+k+'fitrestxt')
                labels[ilab].SetBorderSize(0)
                labels[ilab].SetFillStyle(0)
                labels[ilab].SetTextFont(42)
                labels[ilab].SetTextAlign(12)
                labels[ilab].AddText('Gaussian #1 (f=%3.3f)'%frac.getVal())
                labels[ilab].AddText('#mu=%3.3f#pm%3.3f'%(mean1.getVal(),mean1.getError()))
                labels[ilab].AddText('#sigma=%3.3f#pm%3.3f'%(sigma1.getVal(),sigma1.getError()))
                labels[ilab].AddText('Gaussian #2 (f=%3.3f)'%(1-frac.getVal()))
                labels[ilab].AddText('#mu=%3.3f#pm%3.3f'%(mean2.getVal(),mean2.getError()))
                labels[ilab].AddText('#sigma=%3.3f#pm%3.3f'%(sigma2.getVal(),sigma2.getError()))
                
                labels[ilab].Draw()

            c.Modified()
            c.Update()
            c.SaveAs(r+'res_'+k+'.png')
        
        frame=TGraphErrors()
        frame.SetPoint(0,0,0)
        frame.SetPoint(1,200,0.3)
        frame.SetMarkerStyle(1)
        frame.SetFillStyle(0)
        frame.SetName('frame')
        cresp=TCanvas('cresp','cresp',500,500)
        cresp.cd()
        frame.Draw('ap')
        barrelResponse.Draw('pl')
        endcapResponse.Draw('pl')
        frame.GetXaxis().SetTitle("Quark transverse momentum [GeV]") 
        frame.GetYaxis().SetTitle("Resolution %3.2f C.L."%cl )
        frame.GetYaxis().SetTitleOffset(1.4)
        frame.GetYaxis().SetNdivisions(10)
        drawHeader()
        leg=TLegend(0.6,0.92,0.9,0.98)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.AddEntry(barrelResponse,'barrel','fp')
        leg.AddEntry(endcapResponse,'endcap','fp')
        leg.SetNColumns(2)
        leg.Draw()
        cresp.Modified()
        cresp.Update()
        cresp.SaveAs(r+'res_evol.png')

        frameCore=frame.Clone('framecore')
        cresp.Clear()
        frameCore.Draw('ap')
        barrelCoreResponse.Draw('pl')
        endcapCoreResponse.Draw('pl')
        frameCore.GetXaxis().SetTitle("Quark transverse momentum [GeV]") 
        frameCore.GetYaxis().SetTitle("Core resolution")
        frameCore.GetYaxis().SetTitleOffset(1.4)
        frameCore.GetYaxis().SetNdivisions(10)
        frameCore.GetYaxis().SetRangeUser(0,0.2)
        drawHeader()
        leg=TLegend(0.6,0.92,0.9,0.98)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.AddEntry(barrelCoreResponse,'barrel','fp')
        leg.AddEntry(endcapCoreResponse,'endcap','fp')
        leg.SetNColumns(2)
        leg.Draw()
        cresp.Modified()
        cresp.Update()
        cresp.SaveAs(r+'rescore_evol.png')

    bosons=['h','z','w']
    kin=['','50','100']
    region=['','bb','eb','ee']
    for k in kin:        
        for r in region:

            c=TCanvas('c','c',600,600)
            c.cd()
            histos['mjj'+k+r].Rebin()
            histos['mjj'+k+r].Draw()
            ic=1
            leg=TLegend(0.6,0.92,0.9,0.98)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextFont(42)
            leg.AddEntry(histos['mjj'+k+r],'inclusive','f')
            for b in bosons:
                if histos[b+'mjj'+k+r].Integral()<=0 : continue 
                ic=ic+1
                histos[b+'mjj'+k+r].Rebin()
                histos[b+'mjj'+k+r].SetLineColor(ic)
                histos[b+'mjj'+k+r].SetLineWidth(2)
                histos[b+'mjj'+k+r].SetMarkerColor(ic)
                histos[b+'mjj'+k+r].SetMarkerStyle(1)
                histos[b+'mjj'+k+r].SetFillStyle(3000+ic)
                histos[b+'mjj'+k+r].SetFillColor(ic)
                histos[b+'mjj'+k+r].Draw('histsame')
                leg.AddEntry(histos[b+'mjj'+k+r],b,"f")
            leg.SetNColumns(ic)
            leg.Draw()
            drawHeader()
            labels.append( TPaveText(0.65,0.8,0.9,0.9,'brNDC') )
            ilab=len(labels)-1
            labels[ilab].SetName(k+r+'mjj')
            labels[ilab].SetBorderSize(0)
            labels[ilab].SetFillStyle(0)
            labels[ilab].SetTextFont(42)
            labels[ilab].SetTextAlign(12)
            regionTitle="inclusive"
            if r == 'bb' : regionTitle='barrel-barrel'
            if r == 'eb' : regionTitle='endcap-barrel'
            if r == 'ee' : regionTitle='endcap-endcap'
            labels[ilab].AddText(regionTitle)
            ptthreshold=30
            if k!='' : ptthreshold=float(k)
            labels[ilab].AddText('p_{T}>%3.0f GeV'%ptthreshold)
            labels[ilab].Draw()
            
            c.Modified()
            c.Update()
            c.SaveAs('mjj'+k+r+'.png')


    massResolutionGrs=[]
    for r in region:
        massResolution=TGraphErrors()
        massResolution.SetName(r+'dm')
        massResolution.SetLineWidth(2)
        massResolution.SetFillStyle(0)
        massResolution.SetMarkerStyle(20+len(massResolutionGrs))
        massResolution.SetMarkerColor(1+len(massResolutionGrs))
        massResolution.SetLineColor(1+len(massResolutionGrs))
        massResolution.SetFillColor(1+len(massResolutionGrs))
        massResolutionGrs.append(massResolution)
        
        for k in kin:        

            c=TCanvas('c','c',600,600)
            c.cd()
            h=histos['dmjj'+k+r]
            x=RooRealVar("x", h.GetXaxis().GetTitle(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
            data=RooDataHist("data", "dataset with x", RooArgList(x), h)
            frame=x.frame()
            RooAbsData.plotOn( data, frame, RooFit.DataError(RooAbsData.SumW2) )
            
            mean1=RooRealVar("mean1","mean1",0,-0.5,0.5);
            sigma1=RooRealVar("sigma1","sigma1",0.1,0.01,1.0);
            gauss1=RooGaussian("g1","g",x,mean1,sigma1)
            mean2=RooRealVar("mean2","mean2",0,-0.5,0.5);
            sigma2=RooRealVar("sigma2","sigma2",0.1,0.01,1.0);
            alphacb=RooRealVar("alphacb","alphacb",1,0.1,3);
            ncb=RooRealVar("ncb","ncb",4,1,100)
            gauss2 = RooCBShape("cb2","cb",x,mean2,sigma2,alphacb,ncb);
            frac = RooRealVar("frac","fraction",0.9,0.0,1.0)
            model = RooAddPdf("sum","g1+g2",RooArgList(gauss1,gauss2), RooArgList(frac))
            status=model.fitTo(data,RooFit.Save()).status()
            if status!=0 : continue
            RooAbsPdf.plotOn(model,frame)
            frame.Draw()

            labels.append( TPaveText(0.6,0.65,0.85,0.9,'brNDC') )
            ilab=len(labels)-1
            labels[ilab].SetName(r+k+'dmfitrestxt')
            labels[ilab].SetBorderSize(0)
            labels[ilab].SetFillStyle(0)
            labels[ilab].SetTextFont(42)
            labels[ilab].SetTextAlign(12)
            labels[ilab].AddText('Gaussian #1 (f=%3.3f)'%frac.getVal())
            labels[ilab].AddText('#mu=%3.3f#pm%3.3f'%(mean1.getVal(),mean1.getError()))
            labels[ilab].AddText('#sigma=%3.3f#pm%3.3f'%(sigma1.getVal(),sigma1.getError()))
            labels[ilab].AddText('Gaussian #2 (f=%3.3f)'%(1-frac.getVal()))
            labels[ilab].AddText('#mu=%3.3f#pm%3.3f'%(mean2.getVal(),mean2.getError()))
            labels[ilab].AddText('#sigma=%3.3f#pm%3.3f'%(sigma2.getVal(),sigma2.getError()))
            labels[ilab].Draw()

            drawHeader()
            labels.append( TPaveText(0.15,0.8,0.4,0.9,'brNDC') )
            ilab=len(labels)-1
            labels[ilab].SetName(k+r+'dmjj')
            labels[ilab].SetBorderSize(0)
            labels[ilab].SetFillStyle(0)
            labels[ilab].SetTextFont(42)
            labels[ilab].SetTextAlign(12)
            regionTitle="inclusive"
            if r == 'bb' : regionTitle='barrel-barrel'
            if r == 'eb' : regionTitle='endcap-barrel'
            if r == 'ee' : regionTitle='endcap-endcap'
            labels[ilab].AddText(regionTitle)
            ptthreshold=30
            if k!='' : ptthreshold=float(k)
            labels[ilab].AddText('p_{T}>%3.0f GeV'%ptthreshold)
            labels[ilab].Draw()

            c.Modified()
            c.Update()
            c.SaveAs('dmjj'+k+r+'.png')

            massResolution.SetTitle(regionTitle)
            ip=massResolution.GetN()
            x=40
            xerr=10
            if k=='50' :
                x=75
                xerr=25
            elif k=='100':
                x=150
                xerr=50
            y=sigma1.getVal()
            yerr=sigma1.getError()
            if frac.getVal()<0.8:
                if sigma2.getVal()<sigma1.getVal():
                    y=sigma2.getVal()
                    ey=sigma2.getError()
            massResolution.SetPoint(ip,x,y)
            massResolution.SetPointError(ip,xerr,yerr)
            

    frame=TGraphErrors()
    frame.SetPoint(0,0,0)
    frame.SetPoint(1,200,0.2)
    frame.SetMarkerStyle(1)
    frame.SetFillStyle(0)
    frame.SetName('dmframe')
    cdmevol=TCanvas('cdmevol','cdmevol',500,500)
    cdmevol.cd()
    frame.Draw('ap')
    leg=TLegend(0.6,0.92,0.9,0.98)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    for dmGr in massResolutionGrs :
        dmGr.Draw('pl')
        leg.AddEntry(dmGr,dmGr.GetTitle(),'fp')
    frame.GetXaxis().SetTitle("Leading quark transverse momentum [GeV]") 
    frame.GetYaxis().SetTitle("Core resolution")
    frame.GetYaxis().SetTitleOffset(1.4)
    frame.GetYaxis().SetNdivisions(10)
    drawHeader()
    leg.SetNColumns(2)
    leg.Draw()
    cdmevol.Modified()
    cdmevol.Update()
    cdmevol.SaveAs('dm_evol.png')


    c=TCanvas('c','c',600,600)
    c.cd()
    histos['sel'].Draw('histtext')
    drawHeader()
    c.Modified()
    c.Update()
    c.SaveAs('selection.png')


    return


if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',    dest='input'  , help='input'                , default='analysis.root')
    (opt,args)=parser.parse_args()

    studyVqqResolution(rootFile=opt.input)

