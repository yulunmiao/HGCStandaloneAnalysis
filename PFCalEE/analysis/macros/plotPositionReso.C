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

  std::string plotDir = "../PLOTS/gitV00-02-09/version12/gamma/200um/";
  plotDir += "eta17_et30_pu140";

  TFile *file = TFile::Open((plotDir+".root").c_str());
  
  std::string suffix = "";
  
  file->cd();

  TCanvas *mycL = new TCanvas("mycL","mycL",1500,1000);
  TCanvas *mycD = new TCanvas("mycD","mycD",1500,1000);
  TCanvas *mycR = new TCanvas("mycR","mycR",1500,1000);

  //TCanvas *myc = new TCanvas("myc","myc",1);

  TH2D *p_errorMatrix = (TH2D*)gDirectory->Get("p_errorMatrix");
  TH1F *p_chi2overNDF = (TH1F*)gDirectory->Get("p_chi2overNDF");
  TH1F *p_impactX = (TH1F*)gDirectory->Get("p_impactX");
  TH1F *p_impactX_residual = (TH1F*)gDirectory->Get("p_impactX_residual");
  TH1F *p_impactX_truth = (TH1F*)gDirectory->Get("p_impactX_truth");
  TH1F *p_tanAngleX = (TH1F*)gDirectory->Get("p_tanAngleX");
  TH1F *p_tanAngleX_residual = (TH1F*)gDirectory->Get("p_tanAngleX_residual");
  TH1F *p_tanAngleX_truth = (TH1F*)gDirectory->Get("p_tanAngleX_truth");
  TH1F *p_impactY = (TH1F*)gDirectory->Get("p_impactY");
  TH1F *p_impactY_residual = (TH1F*)gDirectory->Get("p_impactY_residual");
  TH1F *p_impactY_truth = (TH1F*)gDirectory->Get("p_impactY_truth");
  TH1F *p_tanAngleY = (TH1F*)gDirectory->Get("p_tanAngleY");
  TH1F *p_tanAngleY_residual = (TH1F*)gDirectory->Get("p_tanAngleY_residual");
  TH1F *p_tanAngleY_truth = (TH1F*)gDirectory->Get("p_tanAngleY_truth");
  
  TH1F *p_nLayersFit = (TH1F*)gDirectory->Get("p_nLayersFit");
  //TH2F *p_etavsphi_max = (TH2F*)gDirectory->Get("p_etavsphi_max");
  TH2F *p_recoXvsLayer = (TH2F*)gDirectory->Get("p_recoXvsLayer");
  TH2F *p_recoYvsLayer = (TH2F*)gDirectory->Get("p_recoYvsLayer");
  TH2F *p_recoZvsLayer = (TH2F*)gDirectory->Get("p_recoZvsLayer");
  TH2F *p_truthXvsLayer = (TH2F*)gDirectory->Get("p_truthXvsLayer");
  TH2F *p_truthYvsLayer = (TH2F*)gDirectory->Get("p_truthYvsLayer");
  TH2F *p_fitXvsLayer = (TH2F*)gDirectory->Get("p_fitXvsLayer");
  TH2F *p_fitYvsLayer = (TH2F*)gDirectory->Get("p_fitYvsLayer");
  TH1F *p_positionReso = (TH1F*)gDirectory->Get("p_positionReso");
  TH1F *p_angularReso = (TH1F*)gDirectory->Get("p_angularReso");


  mycL->Divide(3,2);
  mycL->cd(1);
  gPad->SetLogz(1);
  p_errorMatrix->SetStats(0);
  p_errorMatrix->SetMinimum(0.01);
  p_errorMatrix->Draw("colz");

  mycL->cd(4);
  gPad->SetLogy(1);
  gStyle->SetOptStat("eMRuo");
  gStyle->SetStatW(0.4);
  //gStyle->SetStatH(0.3);
  p_chi2overNDF->GetXaxis()->SetRangeUser(0,20);
  p_chi2overNDF->Draw();

  mycL->cd(2);
  //gPad->SetLogy(1);
  p_impactX->SetLineColor(1);
  p_impactX->SetMarkerColor(1);
  p_impactX->SetMarkerStyle(21);
  p_impactX->StatOverflows(0);
  double minX = p_impactX_truth->GetMean()-20;
  double maxX = p_impactX_truth->GetMean()+20;
  p_impactX->GetXaxis()->SetRangeUser(minX,maxX);
  //p_impactX->SetMaximum(p_impactX_truth->GetMaximum()*1.1);
  p_impactX->Draw("PE");
  p_impactX_truth->SetLineColor(2);
  p_impactX_truth->Draw("same");

  mycL->cd(3);
  //gPad->SetLogy(1);
  p_impactY->SetLineColor(1);
  p_impactY->SetMarkerColor(1);
  p_impactY->SetMarkerStyle(21);
  p_impactY->StatOverflows(0);
  double minY = p_impactY_truth->GetMean()-20;
  double maxY = p_impactY_truth->GetMean()+60;
  p_impactY->GetXaxis()->SetRangeUser(minY,maxY);
  //p_impactY->SetMaximum(p_impactY_truth->GetMaximum()*1.1);
  p_impactY->Draw("PE");
  p_impactY_truth->SetLineColor(2);
  p_impactY_truth->Draw("same");

  mycL->cd(5);
  gPad->SetLogy(1);
  p_tanAngleX->SetLineColor(1);
  p_tanAngleX->SetMarkerColor(1);
  p_tanAngleX->SetMarkerStyle(21);
  p_tanAngleX->GetXaxis()->SetRangeUser(p_tanAngleX_truth->GetMean()-0.2,p_tanAngleX_truth->GetMean()+0.2);
  p_tanAngleX->SetMaximum(p_tanAngleX_truth->GetMaximum()*1.1);
  p_tanAngleX->Draw("PE");
  p_tanAngleX_truth->SetLineColor(2);
  p_tanAngleX_truth->Draw("same");

  mycL->cd(6);
  gPad->SetLogy(1);
  p_tanAngleY->SetLineColor(1);
  p_tanAngleY->SetMarkerColor(1);
  p_tanAngleY->SetMarkerStyle(21);
  p_tanAngleY->GetXaxis()->SetRangeUser(p_tanAngleY_truth->GetMean()-0.2,p_tanAngleY_truth->GetMean()+0.2);
  p_tanAngleY->SetMaximum(p_tanAngleY_truth->GetMaximum()*1.1);
  p_tanAngleY->Draw("PE");
  p_tanAngleY_truth->SetLineColor(2);
  p_tanAngleY_truth->Draw("same");

  mycL->Update();
  mycL->Print((plotDir+"/PositionFitSummary"+suffix+".pdf").c_str());

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
  p_recoXvsLayer->GetYaxis()->SetRangeUser(p_recoXvsLayer->GetMean(2)-100,p_recoXvsLayer->GetMean(2)+100);
  p_recoXvsLayer->Draw("colz");
  p_recoXvsLayer->SetStats(0);
  p_truthXvsLayer->Draw("same");
  p_truthXvsLayer->SetStats(0);
  mycD->cd(5);
  gPad->SetLogz(1);
  p_recoYvsLayer->RebinY(4);
  p_truthYvsLayer->RebinY(4);
  p_recoYvsLayer->GetYaxis()->SetRangeUser(p_recoYvsLayer->GetMean(2)-200,p_recoYvsLayer->GetMean(2)+200);
  p_recoYvsLayer->Draw("colz");
  p_recoYvsLayer->SetStats(0);
  p_truthYvsLayer->Draw("same");
  p_truthYvsLayer->SetStats(0);
  mycD->cd(3);
  gPad->SetLogz(1);
  //p_fitXvsLayer->RebinY(2);
  p_fitXvsLayer->GetYaxis()->SetRangeUser(p_fitXvsLayer->GetMean(2)-20,p_fitXvsLayer->GetMean(2)+20);
  p_fitXvsLayer->Draw("colz");
  p_fitXvsLayer->SetStats(0);
  mycD->cd(6);
  gPad->SetLogz(1);
  //p_fitYvsLayer->RebinY(4);
  p_fitYvsLayer->GetYaxis()->SetRangeUser(p_fitYvsLayer->GetMean(2)-100,p_fitYvsLayer->GetMean(2)+100);
  p_fitYvsLayer->Draw("colz");
  p_fitYvsLayer->SetStats(0);

  mycD->Update();
  mycD->Print((plotDir+"/PositionFitDebug"+suffix+".pdf").c_str());

  mycR->Divide(3,2);
  mycR->cd(1);
  p_positionReso->GetXaxis()->SetRangeUser(p_positionReso->GetMean()-5*p_positionReso->GetRMS(),
					    p_positionReso->GetMean()+5*p_positionReso->GetRMS());
  p_positionReso->Draw();
  mycR->cd(4);
  p_angularReso->GetXaxis()->SetRangeUser(p_angularReso->GetMean()-5*p_angularReso->GetRMS(),
					    p_angularReso->GetMean()+5*p_angularReso->GetRMS());
  p_angularReso->Draw();
  mycR->cd(2);
  p_impactX_residual->GetXaxis()->SetRangeUser(p_impactX_residual->GetMean()-5*p_impactX_residual->GetRMS(),
					    p_impactX_residual->GetMean()+5*p_impactX_residual->GetRMS());
  p_impactX_residual->Draw();
  mycR->cd(5);
  p_impactY_residual->GetXaxis()->SetRangeUser(p_impactY_residual->GetMean()-5*p_impactY_residual->GetRMS(),
					    p_impactY_residual->GetMean()+5*p_impactY_residual->GetRMS());
  p_impactY_residual->Draw();
  mycR->cd(3);
  p_tanAngleX_residual->GetXaxis()->SetRangeUser(p_tanAngleX_residual->GetMean()-5*p_tanAngleX_residual->GetRMS(),
					    p_tanAngleX_residual->GetMean()+5*p_tanAngleX_residual->GetRMS());
  p_tanAngleX_residual->Draw();
  mycR->cd(6);
  p_tanAngleY_residual->GetXaxis()->SetRangeUser(p_tanAngleY_residual->GetMean()-5*p_tanAngleY_residual->GetRMS(),
						 p_tanAngleY_residual->GetMean()+5*p_tanAngleY_residual->GetRMS());
  p_tanAngleY_residual->Draw();


  mycR->Update();
  mycR->Print((plotDir+"/PositionFitReso"+suffix+".pdf").c_str());


  return 0;

}//main
