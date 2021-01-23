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

int makeBackLeakCor(const unsigned nLayers,
		    const unsigned nBack,
		    const unsigned iSR,
		    const unsigned pT,
		    const unsigned eta,
		    const unsigned pu,
		    const double offset,
		    const double calib,
		    double & backLeakCor,
		    TCanvas * mycE2D,
		    TTree *ltree,
		    TFile *outfile,
		    TGraphErrors * corrBackLeakFit,
		    const TString & plotDir
		    ){

  double Eval = E(pT,eta);

  
  mycE2D->cd();
  gPad->SetRightMargin(0.15);
  std::ostringstream lName;
  std::string lNameTot,lNameBack;
  getTotalEnergyString(nLayers,nBack,lNameTot,lNameBack,iSR);

  lName.str("");
  lName << std::setprecision(6);
  lName << "(";
  lName << lNameTot;
  lName << " - " << offset << ")/" << calib;
  lName << ":(" << lNameBack << ")/(" << lNameTot << ")";

  //std::cout << lName.str().c_str() << std::endl;

  ltree->Draw(lName.str().c_str(),"","colz");


  lName.str("");
  lName << "energy" << pT << "_vsBackFraction";
  TH2F * p_ErecovsEback = (TH2F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str()); // 2D
  
  if (!p_ErecovsEback){
    std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null." << std::endl;
    return 1;
  }

  p_ErecovsEback->Rebin(2);
  p_ErecovsEback->SetTitle(";E_{back}/E_{tot};E_{tot} (GeV)");
  lName << "_pfx";
  TProfile *tmpProf = p_ErecovsEback->ProfileX(lName.str().c_str());
  tmpProf->SetMarkerStyle(20);
  tmpProf->SetMarkerColor(1);
  tmpProf->SetLineColor(1);
  p_ErecovsEback->Draw("colz");
  tmpProf->Draw("PEsame");
  tmpProf->Fit("pol1","","same");
  
  //tmpProf->SetStats(1);
  //gStyle->SetOptFit(1111);
  TF1 *fitcor = (TF1*)tmpProf->GetFunction("pol1");
  if (!fitcor) {
    std::cout << " Fit failed for back leakage correction" << std::endl;
    mycE2D->Update();
    backLeakCor = 0;
    if (Eval < 50) {
      //for low energies: ignore the fit and return success....
      return 0;
    }
    return 1;
  }
  else {
    backLeakCor = fitcor->GetParameter(1);
    std::cout << " ---- back leakage correction factor: " << backLeakCor << " +/- " << fitcor->GetParError(1) << std::endl;
    char buf[500];
    TLatex lat;
    sprintf(buf,"E=%3.3f #times f_{back} + %3.3f",backLeakCor,fitcor->GetParameter(0));
    lat.DrawLatexNDC(0.2,0.85,buf);
		    
    Int_t np=corrBackLeakFit->GetN();
    //if (!dovsE) corrBackLeakFit[iSR]->SetPoint(np,genEn[iE],backLeakCor[oldIdx[iE]][iSR]);
    corrBackLeakFit->SetPoint(np,Eval,backLeakCor);
    corrBackLeakFit->SetPointError(np,0.0,fitcor->GetParError(1));
    
  }
  
  std::ostringstream saveName;
  saveName.str("");
  saveName << plotDir << "/ErecovsbackFraction_eta" << eta << "_pu" << pu;
  saveName << "_E" << pT << "_SR" << iSR;
  mycE2D->Update();
  mycE2D->Print((saveName.str()+".pdf").c_str());
  mycE2D->Print((saveName.str()+".C").c_str());
  
  
  outfile->cd();
  p_ErecovsEback->Write();
  tmpProf->Write();

  
  return 0;
};

int plotBackLeakFit(const TString & plotDir,
		    TGraphErrors *corrBackLeakFit,
		    const unsigned eta,
		    const unsigned pu){
  
  TCanvas *myc = new TCanvas("mycBF","mycBF",1);
  myc->cd();
  corrBackLeakFit->SetTitle(";E (GeV);back cor");
  corrBackLeakFit->Draw("APE");
  myc->Update();
  std::ostringstream lsave;
  lsave.str("");

  if (system(TString("mkdir -p ")+plotDir+TString("/BackCor"))) return 1;

  lsave << plotDir << "/BackCor/";
  lsave << "BackLeakCor";
  lsave << "_eta" << eta << "_pu" << pu ;
  lsave << "_vsE";
  myc->Print((lsave.str()+".pdf").c_str());
  myc->Print((lsave.str()+".C").c_str());

  return 0;
};
