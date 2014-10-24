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

int plotPositionReso(){//main

  std::string plotDir = "../PLOTS/gitV00-02-07/version12/e-/200um/";

  TFile *file = TFile::Open((plotDir+"eta24_e50.root").c_str());

  std::string suffix = "";

  file->cd();

  TCanvas *mycL = new TCanvas("mycL","mycL",1500,1000);
  mycL->Divide(2,2);
  TCanvas *mycD = new TCanvas("mycD","mycD",1500,1000);

  //TCanvas *myc = new TCanvas("myc","myc",1);

  TH2D *p_errorMatrix = (TH2D*)gDirectory->Get("p_errorMatrix");
  TH1F *p_chi2overNDF = (TH1F*)gDirectory->Get("p_chi2overNDF");
  TH1F *p_impactY = (TH1F*)gDirectory->Get("p_impactY");
  TH1F *p_impactY_truth = (TH1F*)gDirectory->Get("p_impactY_truth");
  TH1F *p_angleY = (TH1F*)gDirectory->Get("p_angleY");
  TH1F *p_angleY_truth = (TH1F*)gDirectory->Get("p_angleY_truth");
  
  TH1F *p_nLayersFit = (TH1F*)gDirectory->Get("p_nLayersFit");
  //TH2F *p_etavsphi_max = (TH2F*)gDirectory->Get("p_etavsphi_max");
  TH2F *p_recoXvsLayer = (TH2F*)gDirectory->Get("p_recoXvsLayer");
  TH2F *p_recoYvsLayer = (TH2F*)gDirectory->Get("p_recoYvsLayer");
  TH2F *p_recoZvsLayer = (TH2F*)gDirectory->Get("p_recoZvsLayer");
  TH2F *p_truthXvsLayer = (TH2F*)gDirectory->Get("p_truthXvsLayer");
  TH2F *p_truthYvsLayer = (TH2F*)gDirectory->Get("p_truthYvsLayer");
  TH2F *p_fitXvsLayer = (TH2F*)gDirectory->Get("p_fitXvsLayer");
  TH2F *p_fitYvsLayer = (TH2F*)gDirectory->Get("p_fitYvsLayer");


  mycL->cd(1);
  gPad->SetLogz(1);
  p_errorMatrix->SetStats(0);
  p_errorMatrix->SetMinimum(0.01);
  p_errorMatrix->Draw("colz");

  mycL->cd(2);
  gPad->SetLogy(1);
  gStyle->SetOptStat("eMRuo");
  //gStyle->SetStatW(0.4);
  //gStyle->SetStatH(0.3);
  p_chi2overNDF->Draw();

  mycL->cd(3);
  gPad->SetLogy(1);
  p_impactY->SetLineColor(1);
  p_impactY->SetMarkerColor(1);
  p_impactY->SetMarkerStyle(21);
  p_impactY->StatOverflows(0);
  p_impactY->GetXaxis()->SetRangeUser(580,670);
  p_impactY->SetMaximum(p_impactY_truth->GetMaximum()*1.1);
  p_impactY->Draw("PE");
  p_impactY_truth->SetLineColor(2);
  p_impactY_truth->Draw("same");

  mycL->cd(4);
  gPad->SetLogy(1);
  p_angleY->SetLineColor(1);
  p_angleY->SetMarkerColor(1);
  p_angleY->SetMarkerStyle(21);
  p_angleY->GetXaxis()->SetRangeUser(0,0.36);
  p_angleY->SetMaximum(p_angleY_truth->GetMaximum()*1.1);
  p_angleY->Draw("PE");
  p_angleY_truth->SetLineColor(2);
  p_angleY_truth->Draw("same");

  mycL->Update();
  mycL->Print((plotDir+"PositionFitSummary"+suffix+".pdf").c_str());

  mycD->Divide(3,2);
  mycD->cd(1);
  gPad->SetLogy(1);
  p_nLayersFit->Draw();
  //p_etavsphi_max->Draw("colz");
  //p_etavsphi_max->SetStats(0);
  mycD->cd(4);
  p_recoZvsLayer->Draw("colz");
  p_recoZvsLayer->SetStats(0);
  mycD->cd(2);
  gPad->SetLogz(1);
  p_recoXvsLayer->RebinY(2);
  p_truthXvsLayer->RebinY(2);
  p_recoXvsLayer->Draw("colz");
  p_recoXvsLayer->SetStats(0);
  p_truthXvsLayer->Draw("same");
  p_truthXvsLayer->SetStats(0);
  mycD->cd(5);
  gPad->SetLogz(1);
  p_recoYvsLayer->RebinY(4);
  p_truthYvsLayer->RebinY(4);
  p_recoYvsLayer->Draw("colz");
  p_recoYvsLayer->SetStats(0);
  p_truthYvsLayer->Draw("same");
  p_truthYvsLayer->SetStats(0);
  mycD->cd(3);
  gPad->SetLogz(1);
  //p_fitXvsLayer->RebinY(2);
  p_fitXvsLayer->Draw("colz");
  p_fitXvsLayer->SetStats(0);
  mycD->cd(6);
  gPad->SetLogz(1);
  //p_fitYvsLayer->RebinY(4);
  p_fitYvsLayer->Draw("colz");
  p_fitYvsLayer->SetStats(0);

  mycD->Update();
  mycD->Print((plotDir+"PositionFitDebug"+suffix+".pdf").c_str());



  return 0;

}//main
