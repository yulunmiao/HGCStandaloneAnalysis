#! /usr/bin/env python

from ROOT import TFile, TF1, TH1F, TH2F,TCanvas, TProfile, TPaveText, TLorentzVector, TLegend,TGraphErrors, gStyle, gROOT
from ROOT import RooRealVar,RooGaussian, RooChebychev, RooPolynomial, RooDataHist, RooAbsData, RooArgList, RooCBShape, RooLandau,RooAddPdf, RooAbsPdf, RooFit
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

    c=TCanvas('c','c',600,600)
    #histos['sel'].Draw('histtext')
    #drawHeader()
    #c.Modified()
    #c.Update()
    #c.SaveAs('selection.png')
    #raw_input()

    kin=['','30to50','50to100','100toInf']
    for k in kin:        
        c.Clear()
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
        barrelResponse.SetMarkerStyle(20)
        endcapResponse=TGraphErrors()
        endcapResponse.SetName(r+'endcapresponse')
        endcapResponse.SetMarkerStyle(24)
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

                mean1=RooRealVar("mean1","cb_mean",0,-1,1);
                sigma1=RooRealVar("sigma1","cb_sigma",0.1,0.,1.0);
                gauss1=RooGaussian("g1","g",x,mean1,sigma1)

                mean2=RooRealVar("mean2","mean",0,-1,1);
                sigma2=RooRealVar("sigma2","sigma",0.1,0.,1.0);
                gauss2=RooGaussian("g2","g",x,mean2,sigma2) ;

                frac = RooRealVar("frac","fraction",0.9,0.0,1.0) 
                model = RooAddPdf("sum","g1+g2",RooArgList(gauss1,gauss2), RooArgList(frac))

                model.fitTo(data)
                
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
                
                responseGr=barrelResponse
                if i==2 : responseGr=endcapResponse
                if k!='' :
                    nrespPts=responseGr.GetN()
                    resolutionVal=sigma1.getVal()
                    resolutionErr=sigma1.getError()
                    if math.fabs(mean2.getVal())<math.fabs(mean1.getVal()) :
                        resolutionVal=sigma2.getVal()
                        resolutionErr=sigma2.getError()
                    kinAvg=150
                    kinWidth=50
                    if k=='30to50' : 
                        kinAvg=40
                        kinWidth=10
                    elif k=='50to100' :
                        kinAvg=75
                        kinWidth=25
                    responseGr.SetPoint(nrespPts,kinAvg,resolutionVal)
                    responseGr.SetPointError(nrespPts,kinWidth,resolutionErr)
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
        
        cresp=TCanvas('cresp','cresp',500,500)
        cresp.cd()
        barrelResponse.Draw('ap')
        endcapResponse.Draw('p')
        barrelResponse.GetXaxis().SetTitle("Quark transverse momentum [GeV]") 
        barrelResponse.GetYaxis().SetTitle("Resolution")
        drawHeader()
        leg=TLegend(0.6,0.92,0.9,0.98)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.AddEntry(barrelResponse,'barrel','p')
        leg.AddEntry(endcapResponse,'endcap','p')
        leg.SetNColumns(2)
        leg.Draw()
        cresp.Modified()
        cresp.Update()
        cresp.SaveAs(r+'res_evol.png')
      

    return
    c.Clear()
    histos['mjj'].Rebin(10)
    histos['mjj'].Draw()
    bosons=['h','z','w']
    ic=1
    leg=TLegend(0.15,0.75,0.5,0.9)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.AddEntry(histos['mjj'],'inclusive','f')
    for b in bosons:
        if histos[b+'mjj'].Integral()<=0 : continue 
        ic=ic+1
        histos[b+'mjj'].Rebin(10)
        histos[b+'mjj'].SetLineColor(ic)
        histos[b+'mjj'].SetLineWidth(2)
        histos[b+'mjj'].SetMarkerColor(ic)
        histos[b+'mjj'].SetMarkerStyle(1)
        histos[b+'mjj'].Draw('histsame')
        leg.AddEntry(histos[b+'mjj'],b,"f")
    leg.Draw()
    drawHeader()
    c.Modified()
    c.Update()
    c.SaveAs('mjj.png')
    raw_input()

    # resolutions
    resFunc=TF1('resfunc','[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])',-2,2)
    resFunc.SetLineColor(4)
    resFunc.SetParName(0,'N_{1}')
    resFunc.SetParName(1,'#mu_{1}')
    resFunc.SetParName(2,'#sigma_{1}')
    resFunc.SetParName(3,'N_{2}')
    resFunc.SetParName(4,'#mu_{2}')
    resFunc.SetParName(5,'#sigma_{2}')


    allVars=['dpt','den','inmdpt','inmden']
    kin=['30to50','50to100','100toInf']
    reg=['','barrel','endcap']
    for v in allVars:
        c.Clear()
        c.SetWindowSize(1500,1000)
        c.Divide(len(kin),len(reg))
        
        i=0
        for r in reg:
            for k in kin:
                i=i+1
                p=c.cd(i)
                h=histos[v+k+r]
                h.Draw()
                resFunc.SetParameter(0,h.GetMaximum())
                if k!='' : resFunc.SetParLimits(0,h.GetMaximum()*0.8,h.GetMaximum())
                else : resFunc.SetParLimits(0,h.GetMaximum()*0.2,h.GetMaximum()*1.2)
                offset=h.GetXaxis().GetBinCenter(h.GetMaximumBin())
                #if(math.fabs(h.GetMean())>0.2) : offset=h.GetMean()
                resFunc.SetParameter(1,offset)
                resFunc.SetParLimits(1,offset-2*h.GetMeanError(),offset+2*h.GetMeanError())
                resFunc.SetParLimits(2,0.0,0.15)
                resFunc.SetParameter(4,h.GetMaximumBin())
                resFunc.SetParLimits(4,-1,1)
                resFunc.SetParLimits(5,0.1,2)
                h.Fit(resFunc,'RQ','',-2,2)
                
                if i==1 : drawHeader()
                labels.append( TPaveText(0.1,0.9,0.5,0.8,'brNDC') )
                ilab=len(labels)-1
                labels[ilab].SetName(k+r+'txt')
                labels[ilab].SetBorderSize(0)
                labels[ilab].SetFillStyle(0)
                labels[ilab].SetTextFont(42)
                labels[ilab].SetTextAlign(12)
                if r=='' : labels[ilab].AddText('[inclusive]')           
                else     : labels[ilab].AddText('['+r+']')
                kinReg=k.replace('to','-')
                kinReg=kinReg.replace('Inf','#infty')
                labels[ilab].AddText('p_{T} ' + kinReg)
                labels[ilab].Draw()
                p.Modified()
                p.Update()

        c.Modified()
        c.Update()
        c.SaveAs(v+'res.png')
        raw_input()


if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',    dest='input'  , help='input'                , default='analysis.root')
    (opt,args)=parser.parse_args()

    studyVqqResolution(rootFile=opt.input)

