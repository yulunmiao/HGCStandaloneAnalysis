#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"

#include "TDRStyle.h"

int plotSR(){

  SetTdrStyle();

  const unsigned nLayers = 30;

  TCanvas *myc = new TCanvas("myc","",1);
  gStyle->SetOptStat(0);
  myc->Divide(2,1);

  TH2F *h5 = new TH2F("h5","SR 5;width (pads);layer",
		     5,-2.5,2.5,nLayers,0,nLayers);
  TH2F *h6 = new TH2F("h6","SR 6;width (pads);layer",
		     5,-2.5,2.5,nLayers,0,nLayers);

  unsigned iSR = 5;

  for (unsigned iL(0);iL<nLayers;++iL){	      
    if (iL==0) h5->Fill(0.,iL);
    else if (iL<5) h5->Fill(0.,iL);
    else if (iL<10) {h5->Fill(0.,iL);h5->Fill(1.,iL);}
    else if (iL<15) {h5->Fill(-1.,iL);h5->Fill(0.,iL);h5->Fill(1.,iL);}
    else if (iL<20) {h5->Fill(-1.,iL);h5->Fill(0.,iL);h5->Fill(1.,iL);h5->Fill(2.,iL);}
    else {h5->Fill(-2.,iL);h5->Fill(-1.,iL);h5->Fill(0.,iL);h5->Fill(1.,iL);h5->Fill(2.,iL);}
  }

  myc->cd(1);
  h5->Draw("col");

  TLine *line = new TLine(0.2,0,0.2,nLayers);
  line->Draw();

  TLatex lat;
  lat.DrawLatex(0.3,15,"photon");
  lat.DrawLatexNDC(0.4,0.96,"SR 5");

  iSR = 6;
  for (unsigned iL(0);iL<nLayers;++iL){	      
    if (iL==0) h6->Fill(0.,iL);
    else if (iL<5) {h6->Fill(0.,iL);h6->Fill(1.,iL);}
    else if (iL<12) {h6->Fill(-1.,iL);h6->Fill(0.,iL);h6->Fill(1.,iL);}
    else {h6->Fill(-2.,iL);h6->Fill(-1.,iL);h6->Fill(0.,iL);h6->Fill(1.,iL);h6->Fill(2.,iL);}
  }

  myc->cd(2);
  h6->Draw("col");
  line->Draw();
  lat.DrawLatex(0.3,15,"photon");
  lat.DrawLatexNDC(0.4,0.96,"SR 6");


  return 0;
}//main
