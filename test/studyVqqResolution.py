#! /usr/bin/env python

from ROOT import TFile, TF1, TH1F, TH2F,TCanvas, TProfile, TPaveText, TLorentzVector, TLegend, gStyle, gROOT
import optparse

allHeaders=[]

def drawHeader() :
    allHeaders.append(TPaveText(0.08,0.94,0.5,0.98,'brNDC'))
    ih=len(allHeaders)-1
    allHeaders[ih].SetBorderSize(0)
    allHeaders[ih].SetFillStyle(0)
    allHeaders[ih].SetTextAlign(12)
    allHeaders[ih].AddText('ILD simulation')
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
    gROOT.SetStyle('Plain')
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(1111)
    gStyle.SetOptTitle(0)
    gStyle.SetStatFont(42)

    c=TCanvas('c','c',600,600)
    histos['sel'].Draw('histtext')
    drawHeader()
    c.Modified()
    c.Update()
    c.SaveAs('selection.png')

    c.Clear()
    histos['j1dr'].SetLineWidth(2)
    histos['j1dr'].SetMarkerStyle(1)
    histos['j1dr'].SetTitle('Jet #1')
    histos['j1dr'].Draw('hist')
    histos['j2dr'].SetLineWidth(2)
    histos['j2dr'].SetMarkerStyle(1)
    histos['j2dr'].SetMarkerColor(2)
    histos['j2dr'].SetLineColor(2)
    histos['j2dr'].SetTitle('Jet #2')
    histos['j2dr'].Draw('histsame')
    leg=c.BuildLegend()
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    drawHeader()
    c.Modified()
    c.Update()
    c.SaveAs('dr.png')

    c.Clear()
    histos['mjj'].Draw()
    drawHeader()
    c.Modified()
    c.Update()
    c.SaveAs('mjj.png')


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
    #cats=['j1','j2','j']
    cats=['j']
    kin=['30to50','50to100','100toInf']
    reg=['','barrel','endcap']
    labels=[]
    for v in allVars:
        for cat in cats:
            c.Clear()
            c.SetWindowSize(1500,1000)
            c.Divide(len(kin),len(reg))

            i=0
            for r in reg:
                for k in kin:
                    i=i+1
                    p=c.cd(i)
                    h=histos[cat+v+k+r]
                    h.Draw()
                    resFunc.SetParameter(0,h.GetMaximum())
                    if k!='' : resFunc.SetParLimits(0,h.GetMaximum()*0.8,h.GetMaximum())
                    else : resFunc.SetParLimits(0,h.GetMaximum()*0.2,h.GetMaximum()*1.2)
                    resFunc.SetParameter(1,0.0)
                    resFunc.SetParLimits(2,0.0,0.15)
                    resFunc.SetParameter(4,0.0)
                    resFunc.SetParLimits(4,-1,1)
                    resFunc.SetParLimits(5,0.1,2)
                    h.Fit(resFunc,'RQ','',-2,2)
                
                    if i==1 : drawHeader()
                    labels.append( TPaveText(0.1,0.9,0.5,0.8,'brNDC') )
                    ilab=len(labels)-1
                    labels[ilab].SetName(cat+k+r+'txt')
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
            c.SaveAs(cat+v+'res.png')
            raw_input()


if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',    dest='input'  , help='input'                , default='analysis.root')
    (opt,args)=parser.parse_args()

    studyVqqResolution(rootFile=opt.input)

