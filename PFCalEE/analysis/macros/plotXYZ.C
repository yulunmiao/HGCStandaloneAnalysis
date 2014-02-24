#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"


int plotXYZ(){//main  

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 2;
  std::string scenario[nS] = {
    "quark_u/eta30/",
    "quark_u/PU/eta30/"
  };
  
  const unsigned nV = 1;
  TString version[nV] = {"23"};
  const double Emip = 0.0548;//in MeV

  const unsigned mipThresh = 10;

  const unsigned nEvts = 1000;
  unsigned event[nEvts];// = {190};//,6,12};
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
    event[ievt] = ievt;
  }

  const unsigned nLayers = 60;
  const unsigned nEcalLayers = 30;

  TCanvas *mycECAL = new TCanvas("mycECAL","mycECAL",1500,1000);
  TCanvas *mycHCAL = new TCanvas("mycHCAL","mycHCAL",1500,1000);
  TCanvas *mycAll = new TCanvas("mycAll","mycAll",1500,1000);
  const unsigned nPads = nEvts>10 ? 10 : nEvts;
  mycAll->Divide(static_cast<unsigned>(nPads/2.+0.5),nEvts/2<1 ? 1 : 2);


  const unsigned nCanvas = nS;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
  }
  
  std::ostringstream saveName;
  
  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
    
      TString plotDir = "../PLOTS/version_"+version[iV]+"/"+scenario[iS]+"/";
      

      bool isRECO = false;
      //if (scenario[iS].find("scenario_") != scenario[iS].npos) isRECO=true;

      if (isRECO) plotDir += "Reco/";

      for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
	std::ostringstream lName;
	lName << plotDir << "CalibHistos_E100_evt" << event[ievt] << ".root";
	TFile *inputFile = TFile::Open(lName.str().c_str());
	if (!inputFile) {
	  std::cout << " -- Error, input file " << lName.str() << " cannot be opened. Going to next..." << std::endl;
	  continue;
	  //return 1;
	}
	TH2F *p_xy[nLayers];
	TH3F *p_xyz = 0;
	if (!isRECO) p_xyz = (TH3F*)gDirectory->Get("p_xyz")->Clone();
	else p_xyz = (TH3F*)gDirectory->Get("p_recoxyz")->Clone();
      
	if (!p_xyz) {
	  std::cout << " -- ERROR, pointer for XYZ histogram is null. Exiting..." << std::endl;
	  return 1;
	}
	p_xyz->Sumw2();
	double EmaxEcal = 0;
	double EmaxHcal = 0;

	for (unsigned iL(0); iL<nLayers; ++iL){
	  lName.str("");
	  if (!isRECO) lName << "p_xy_" << iL;
	  else lName << "p_recoxy_" << iL;
	  p_xy[iL] = (TH2F*)gDirectory->Get(lName.str().c_str());
	  if (!p_xy[iL]) {
	    std::cout << " -- ERROR, pointer for histogram is null for layer: " << iL << ". Exiting..." << std::endl;
	    return 1;
	  }
	  p_xy[iL]->Scale(1./Emip);
	  if (!isRECO) {
	    p_xy[iL]->RebinX(4);
	    p_xy[iL]->RebinY(4);
	  }
	  double Etot = p_xy[iL]->GetMaximum();
	  if (Etot > EmaxEcal && iL<nEcalLayers) EmaxEcal = Etot;
	  if (Etot > EmaxHcal && iL>=nEcalLayers) EmaxHcal = Etot;
	}


	gStyle->SetOptStat(0);
	std::ostringstream ltitle;
	ltitle << "Event #" << event[ievt] ;

	myc[iS]->cd();
	p_xyz->Scale(1./Emip);
	if (!isRECO) {
	  p_xyz->RebinY(4);
	  p_xyz->RebinZ(4);
	}
	//p_xyz->SetMinimum(100);
	for (int xb(1); xb<p_xyz->GetNbinsX()+1;++xb){
	  for (int yb(1); yb<p_xyz->GetNbinsY()+1;++yb){
	    for (int zb(1); zb<p_xyz->GetNbinsZ()+1;++zb){
	      //std::cout << xb << " " << yb << " " << zb << " " << p_xyz->GetBinContent(xb,yb,zb) << std::endl;
	      if (p_xyz->GetBinContent(xb,yb,zb) < mipThresh) p_xyz->SetBinContent(xb,yb,zb,0);
	    }
	  }
	}


	p_xyz->GetXaxis()->SetLabelSize(0.05);
	p_xyz->GetYaxis()->SetLabelSize(0.05);
	p_xyz->GetZaxis()->SetLabelSize(0.05);
	p_xyz->GetXaxis()->SetTitleSize(0.05);
	p_xyz->GetYaxis()->SetTitleSize(0.05);
	p_xyz->GetZaxis()->SetTitleSize(0.05);
	p_xyz->SetTitle(ltitle.str().c_str());
	p_xyz->Draw("");
	
	myc[iS]->Update();
	saveName.str("");
	saveName << plotDir << "/xyzSimHits_evt" << event[ievt] << "_mipThresh" << mipThresh;
	myc[iS]->Print((saveName.str()+".png").c_str());
	myc[iS]->Print((saveName.str()+".pdf").c_str());
	
	//continue;

	mycAll->cd(ievt%nPads+1);
	gPad->SetLogz(1);
	TString lLabel = "yz_";
	lLabel += ievt;
	TH2D *proj = (TH2D*)p_xyz->Project3D(lLabel);
	proj->GetZaxis()->SetLabelSize(0.05);
	proj->GetZaxis()->SetTitleSize(0.05);
	proj->GetZaxis()->SetTitle("E(MIP)");
	proj->GetZaxis()->SetTitleOffset(-0.5);
	proj->SetTitle(ltitle.str().c_str());
	//proj->SetMinimum(0.1);
	//proj->Draw("LEGO2z");
	proj->Draw("colz");

	if (ievt%nPads == nPads-1){
	  mycAll->Update();
	  saveName.str("");
	  saveName << plotDir << "/xySimHits_proj_" << ievt/nPads << "_mipThresh" << mipThresh;
	  mycAll->Print((saveName.str()+".png").c_str());
	  mycAll->Print((saveName.str()+".pdf").c_str());
	}
	//continue;

	mycECAL->Clear();
	unsigned counter = 1;
	mycECAL->Divide(3,3);

	mycHCAL->Clear();
	mycHCAL->Divide(5,3);


	for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
	  //std::cout << " -- Processing layer " << iL << std::endl;
	  if (iL%3==2 && counter < 10 && iL<nEcalLayers) {
	    std::cout << " -- Ecal layer " << iL << " counter = " << counter << std::endl;
	    mycECAL->cd(counter);
	    counter++;
	  }
	  else if (counter<25 && iL>=nEcalLayers){
	    std::cout << " -- Hcal layer " << iL << " counter = " << counter << std::endl;
	    mycHCAL->cd(counter-9);
	    counter++;
	  }
	  else if (iL==nEcalLayers+15){
	    saveName.str("");
	    saveName << plotDir << "/xySimHits_HCAL_evt" << event[ievt] << "_subset1";
	    mycHCAL->Update();
	    mycHCAL->Print((saveName.str()+".png").c_str());
	    mycHCAL->Print((saveName.str()+".pdf").c_str());
	    counter = 10;
	    std::cout << " -- Hcal layer " << iL << " counter = " << counter << std::endl;
	    mycHCAL->cd(counter-9);
	    counter++;
	  } 
	  else continue;
	  gPad->SetLogz(1);
	  p_xy[iL]->SetMaximum(iL<nEcalLayers ? EmaxEcal : EmaxHcal);
	  p_xy[iL]->GetXaxis()->SetLabelSize(0.05);
	  p_xy[iL]->GetYaxis()->SetLabelSize(0.05);
	  p_xy[iL]->GetXaxis()->SetTitleSize(0.05);
	  p_xy[iL]->GetYaxis()->SetTitleSize(0.05);
	  p_xy[iL]->GetZaxis()->SetTitleOffset(-0.5);
	  p_xy[iL]->GetZaxis()->SetTitle("E(MIP)");

	  //p_xy[iL]->Draw("colz");
	  char buf[500];
	  sprintf(buf,"Layer %d",iL);
	  p_xy[iL]->SetTitle(buf);

	  //p_xy[iL]->Draw("LEGO2z");
	  p_xy[iL]->Draw("colz");
	}//loop on layers
	saveName.str("");
	saveName << plotDir << "/xySimHits_ECAL_evt" << event[ievt] << "_subset";
	mycECAL->Update();
	mycECAL->Print((saveName.str()+".png").c_str());
	mycECAL->Print((saveName.str()+".pdf").c_str());
	saveName.str("");
	saveName << plotDir << "/xySimHits_HCAL_evt" << event[ievt] << "_subset2";
	mycHCAL->Update();
	mycHCAL->Print((saveName.str()+".png").c_str());
	mycHCAL->Print((saveName.str()+".pdf").c_str());


      }//loop on events

 
    }//loop on scenarios

  }//loop on versions
  
  return 0;


}//main
