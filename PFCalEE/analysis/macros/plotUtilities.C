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
#include "TProfile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"

struct FitResult{
  double chi2;
  unsigned ndf;
  double mean;
  double sigma;
  double meanerr;
  double sigmaerr;
};

double E(const unsigned pT, const unsigned eta){
  return pT*cosh(eta/10.);
};

double pT(const unsigned E, const unsigned eta){
  return E/cosh(eta/10.);
};

void drawChi2(TCanvas *myc,TH1F ** p_chi2ndf,const unsigned nSR){
  
  gStyle->SetOptStat("eMRuo");
  gStyle->SetStatH(0.4);
  gStyle->SetStatW(0.4);
  for (unsigned iSR(0); iSR<nSR;++iSR){
    myc->cd(iSR+1);
    if (p_chi2ndf[iSR]) p_chi2ndf[iSR]->Draw();
  }
  
  myc->Update();
  std::ostringstream lsave;
  lsave << "PLOTS/EnergyFitQuality.pdf";
  myc->Print(lsave.str().c_str());
}

TPad* plot_ratio(TPad *canv, bool up){
  canv->SetFillColor      (0);
  canv->SetBorderMode     (0);
  canv->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canv->SetLeftMargin     (0.18);
  canv->SetLeftMargin     (0.17);
  canv->SetRightMargin    (0.05);
  canv->SetTopMargin      (0.06);
  canv->SetBottomMargin   (0.18);
  // Setup a frame which makes sense
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);      

  canv->cd();
  TPad *pad = 0;
  if (up){
    pad = new TPad("upper","pad",0, 0.26 ,1 ,1);
    pad->SetBottomMargin(0.01);
    pad->SetTopMargin(0.06);
    pad->Draw();
    pad->cd();
    return pad;
  }
  else {
    pad = new TPad("lower","pad",0, 0   ,1 ,0.26);  
    pad->SetTopMargin(0.01);
    pad->SetBottomMargin(0.25);
    pad->Draw();
    return pad;
  }

};


void getTotalEnergyString(const unsigned nLayers,
			  const unsigned nBack,
			  std::string & totE,
			  std::string & backE,
			  const unsigned iSR){

  std::ostringstream lNameTot,lNameBack;
  lNameTot.str("");
  lNameBack.str("");

  if (iSR<6){
    for (unsigned iL(0);iL<nLayers;++iL){	      
      if (iL==0) lNameTot << "absweight_" << iL  << "*energy_" << iL << "_SR" << iSR ;
      else lNameTot << "+absweight_" << iL  << "*energy_" << iL << "_SR" << iSR ;
    }
  }
  else lNameTot << "wgtEtotal";///" << tanh(etaval[ieta]);

  for (unsigned iL(nLayers-nBack); iL<nLayers;++iL){
    if (iL==nLayers-nBack) lNameBack << "absweight_" << iL << "*energy_" << iL << "_SR" << iSR ;
    else lNameBack << "+absweight_" << iL << "*energy_" << iL << "_SR" << iSR ;
  }
  
  
  totE = lNameTot.str();
  backE = lNameBack.str();
  
};
