#! /usr/bin/env python

from ROOT import TFile, TH1F, TH2F,TCanvas, TProfile, TPaveText, TLorentzVector, TLegend, gStyle
import optparse

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
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)

    c=TCanvas('c','c',1200,600)
    c.Divide(2,2)
    c.cd(1)
    histos['nj'].Draw('histe1')
    histos['nq'].Draw('histsame')
    pt=TPaveText(0.08,0.9,0.5,0.99,'brNDC')
    pt.SetBorderSize(0)
    pt.SetFillStyle(0)
    pt.AddText('ILD simulation')
    pt.Draw()
    leg=TLegend(0.15,0.6,0.5,0.8,'p_{T}>30 GeV','brNDC')
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.AddEntry(histos['nj'],'jet','f')
    leg.AddEntry(histos['nq'],'quark','f')
    leg.Draw()
    p=c.cd(3)
    p.Divide(3,1)
    p.cd(1)
    histos['jpt'].Draw('histe1')
    histos['qpt'].Draw('histsame')
    p.cd(2)
    histos['jeta'].Draw('histe1')
    histos['qeta'].Draw('histsame')
    p.cd(3)
    histos['jmass'].Draw('histe1')
    histos['qmass'].Draw('histsame')
    c.cd(2)
    histos['dptvspt'].Draw('e1')
    c.cd(4)
    histos['dptvseta'].Draw('e1')
    c.SaveAs('jetkin.png')
    
    raw_input()


if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',    dest='input'  , help='input'                , default='analysis.root')
    (opt,args)=parser.parse_args()

    studyVqqResolution(rootFile=opt.input)

