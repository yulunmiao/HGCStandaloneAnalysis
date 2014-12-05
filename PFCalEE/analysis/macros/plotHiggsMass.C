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


int plotHiggsMass(){//main

  const unsigned nP = 1;
  unsigned puVal[nP] = {140};//,140};

  const unsigned nV = 7;

  std::string label[nV] = {
    "True E, Shower pos",
    "True E, Vtx smear",
    "True E, Shower angle",
    "Reco E, True pos",
    "Reco E, Shower pos",
    "Reco E, Vtx smear",
    "Reco E, Shower angle",
  };

  double val[nV];
  double valerr[nV];
  for (unsigned iV(0);iV<nV;++iV){
    val[iV] = iV+0.5;
    valerr[iV] = 0;
  }

  //fill arrays
  double mean[nP][nV];
  double meanerr[nP][nV];
  double sigma[nP][nV];
  double sigmaerr[nP][nV];
  
  TFile *f[nP];
  f[0] = TFile::Open("../PLOTS/gitV00-02-10/version12/HggSmall/200um/pu0.root");
  f[1] = TFile::Open("../PLOTS/gitV00-02-10/version12/HggSmall/200um/pu140.root");

  TCanvas *myc = new TCanvas("myc","myc",1);
  TCanvas *myc2 = new TCanvas("myc2","myc2",1500,1000);
  myc2->Divide(4,2);
  gStyle->SetOptStat("eMRuo");
  gStyle->SetOptFit(1111);
  //gStyle->SetStatH(0.2);
  //gStyle->SetStatW(0.4);
  for (unsigned ipu(0);ipu<nP;++ipu){//loop on pu
    if (!f[ipu]) {
      std::cout << " -- Input file for pu " << puVal[ipu] << " not found ! " << std::endl;
      return 1;
    }
    f[ipu]->cd("HiggsMass");
    TH1F *hMass[nV];
    hMass[0] = (TH1F*)gDirectory->Get("p_position_trueE");
    hMass[1] = (TH1F*)gDirectory->Get("p_position_vtxsmear_trueE");
    hMass[2] = (TH1F*)gDirectory->Get("p_angle_trueE");
    hMass[3] = (TH1F*)gDirectory->Get("p_trueDir_recoE");
    hMass[4] = (TH1F*)gDirectory->Get("p_position_recoE");
    hMass[5] = (TH1F*)gDirectory->Get("p_position_vtxsmear_recoE");
    hMass[6] = (TH1F*)gDirectory->Get("p_angle_recoE");
   
    TH1F *trueM = (TH1F*)gDirectory->Get("p_trueDir_trueE");
    myc2->cd(1);
    trueM->Draw();

    for (unsigned iV(0);iV<nV;++iV){//loop on masses
      myc2->cd(2+iV);
      hMass[iV]->Draw();
      hMass[iV]->Fit("gaus","+","same",110,160);
      TF1 *fit = (TF1*)hMass[iV]->GetFunction("gaus");
      mean[ipu][iV] = fit->GetParameter(1);
      meanerr[ipu][iV] = fit->GetParError(1);
      sigma[ipu][iV] = fit->GetParameter(2);
      sigmaerr[ipu][iV] = fit->GetParError(2);
      //mean[ipu][iV] = hMass[iV]->GetMean();
      //meanerr[ipu][iV] = hMass[iV]->GetMeanError();
      //sigma[ipu][iV] = hMass[iV]->GetRMS();
      //sigmaerr[ipu][iV] = hMass[iV]->GetRMSError();
    }
    std::ostringstream lsave;
    lsave << "PLOTS/HiggsMasses_pu" << puVal[ipu] << ".pdf";
    myc2->Print(lsave.str().c_str());

  }

  myc->cd();
  gPad->SetGridy(1);

  gStyle->SetOptStat(0);
  TGraphErrors *grMass[nP];
  TGraphErrors *grSigma[nP];
 
  for (unsigned iH(0); iH<2; ++iH){
    TGraphErrors *gr[nP];

    for (unsigned iP(0); iP<nP; ++iP){
      
      grMass[iP] = new TGraphErrors(nV,val,mean[iP],valerr,meanerr[iP]);
      grSigma[iP] = new TGraphErrors(nV,val,sigma[iP],valerr,sigmaerr[iP]);
      
      grMass[iP]->SetTitle(";;Mass (GeV)");
      grSigma[iP]->SetTitle(";;#sigma_{M} (GeV)");
      
      gr[iP] = iH==0? grMass[iP] : grSigma[iP];
      gr[iP]->SetMarkerStyle(21+iP);
      gr[iP]->SetMarkerColor(iP+1);
      if (iH==1) {
	gr[iP]->SetMinimum(0);
	gr[iP]->SetMaximum(4.5);
      }
      else {
	gr[iP]->SetMinimum(120);
	gr[iP]->SetMaximum(126);
      }
      TAxis *ax = gr[iP]->GetHistogram()->GetXaxis();
      Double_t x1 = ax->GetBinLowEdge(1);
      Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
      gr[iP]->GetHistogram()->GetXaxis()->Set(nV,x1,x2);
      
      for(Int_t k=0;k<nV;k++){
	gr[iP]->GetHistogram()->GetXaxis()->SetBinLabel(k+1,label[k].c_str());
      }
      
      if (iP==0) gr[iP]->Draw("AP");
      else gr[iP]->Draw("P");

      TLatex lat;
      lat.SetTextColor(iP+1);
      char buf[100];
      sprintf(buf,"PU = %d",puVal[iP]);
      lat.DrawLatexNDC(0.2+0.2*iP,0.2,buf);
    }//loop on PU

    TLatex lat;
    lat.DrawLatexNDC(0.1,0.92,"Pythia gg->Higgs");
    myc->Update();
    if (iH==0) myc->Print("PLOTS/SummaryHiggsMass.pdf");
    else myc->Print("PLOTS/SummaryHiggsReso.pdf");

  }//loop on histos


  return 0;


}//main
