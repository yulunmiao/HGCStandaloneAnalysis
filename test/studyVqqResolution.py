#! /usr/bin/env python

from ROOT import TFile, TF1, TH1F, TH2F,TCanvas, TProfile, TPaveText, TLorentzVector, TLegend, gStyle, gROOT
import optparse

def drawHeader() :
    pt=TPaveText(0.08,0.9,0.5,0.99,'brNDC')
    pt.SetBorderSize(0)
    pt.SetFillStyle(0)
    pt.AddText('ILD simulation')
    pt.Draw()



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

    c=TCanvas('c','c',600,600)
    histos['sel'].Draw('histtext')
    drawHeader()
    c.Modified()
    c.Update()
    c.SaveAs('selection.png')
    #raw_input()

    c.Clear()
    histos['j1dr'].SetLineWidth(2)
    histos['j1dr'].SetTitle('Jet #1')
    histos['j1dr'].Draw('hist')
    histos['j2dr'].SetLineWidth(2)
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
    #raw_input()


    resFunc=TF1('resfunc','[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])',-2,2)
    resFunc.SetLineColor(4)

    cats=['j1','j2','j']
    kin=['','30to50','50to100','100toInf']
    reg=['','barrel','endcap']
    for cat in cats:
        c.Clear()
        c.SetWindowSize(1500,1000)
        c.Divide(len(kin),len(reg))

        i=0
        for r in reg:
            for k in kin:
                i=i+1
                p=c.cd(i)
                h=histos[cat+'dpt'+k+r]
                h.Draw()
                #h.Fit('gaus','','',-0.25,0.25)
               # gaus=h.GetFunction('gaus')
                #gaus.SetLineColor(4)
                #gaus.SetLineWidth(1)
                #resFunc.SetParameter(0,gaus.GetParameter(0))
                resFunc.SetParameter(0,h.GetMaximum())
                resFunc.SetParLimits(0,h.GetMaximum()*0.6,h.GetMaximum())
                resFunc.SetParameter(1,0.0)
                resFunc.SetParLimits(2,0.0,0.15)
                resFunc.SetParameter(4,0.0)
                resFunc.SetParLimits(5,0.1,2)
                #h.Fit(resFunc,'QI','',-0.25,1)
                #h.Fit(resFunc,'QI','',-0.25,2)
                h.Fit(resFunc,'RQ','',-2,2)
                

                
                if i==1 : drawHeader()
                pt=TPaveText(0.1,0.9,0.5,0.8,'brNDC')
                pt.SetBorderSize(0)
                pt.SetFillStyle(0)
                pt.SetTextFont(42)
                pt.SetTextAlign(12)
                pt.AddText(k)
                pt.AddText(r)
                pt.Draw()
                p.Modified()
                p.Update()

        c.Modified()
        c.Update()
        c.SaveAs(cat+'res.png')
        raw_input()


    c.Clear()
    histos['mjj'].Draw()
    drawHeader()
    c.Modified()
    c.Update()
    c.SaveAs('mjj.png')
    raw_input()


if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',    dest='input'  , help='input'                , default='analysis.root')
    (opt,args)=parser.parse_args()

    studyVqqResolution(rootFile=opt.input)

