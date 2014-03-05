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


int plotPions(){//main  

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 1;
  std::string scenario[nS] = {
    //"pi-/twiceSampling/GeVCal/"
    "pi-/concept/GeVCal/"
  };
  
  const unsigned nV = 1;
  TString version[nV] = {"20"};

  std::string pDetector = "HCALvsECAL";

  //const double EmipSi = 0.0822;//in MeV

  const unsigned nLayers = 64;//64 //54
  const unsigned nEcalLayers = 31;

  //unsigned genEn[]={5,10,25,40,50,60,80,100,150,200,300,400,500};//,1000,2000};
  unsigned genEn[]={50,100,150,200};
  //unsigned genEn[]={10};
  unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);

  TCanvas *mycECAL = new TCanvas("mycECAL","mycECAL",1500,1000);
  TCanvas *myc = new TCanvas("myc","myc",1);
  mycECAL->Divide(5,3);

  std::ostringstream saveName;
  
  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
    
      TString plotDir = "../PLOTS/version"+version[iV]+"/"+scenario[iS]+"/";
      
      std::ostringstream lName;
      lName << plotDir << "CalibHistos.root";
      TFile *inputFile = TFile::Open(lName.str().c_str());
      if (!inputFile) {
	std::cout << " -- Error, input file " << lName.str() << " cannot be opened. Going to next..." << std::endl;
	continue;
	//return 1;
      }


      TH2F *p_HCALvsECAL[nGenEn];

      double slope[nGenEn];
      double slopeErr[nGenEn];
      double genEnErr[nGenEn];
      double grgenEn[nGenEn];

      for (unsigned iE(0); iE<nGenEn; ++iE){

	std::cout << " -- Processing energy " << genEn[iE] << std::endl;
	genEnErr[iE] = 0;
	grgenEn[iE] = genEn[iE];

	lName.str("");
	lName << "p_HCALvsECAL_" << genEn[iE];
	p_HCALvsECAL[iE] = (TH2F*)gDirectory->Get(lName.str().c_str());
	if (!p_HCALvsECAL[iE]) {
	  std::cout << " -- ERROR, pointer for histogram is null for energy: " << genEn[iE] << ". Exiting..." << std::endl;
	  return 1;
	}
	if (pDetector.find("ScivsSi") != pDetector.npos) {
	  p_HCALvsECAL[iE]->GetXaxis()->SetTitle("HCAL Si");
	  p_HCALvsECAL[iE]->GetYaxis()->SetTitle("HCAL Scint");
	}
	p_HCALvsECAL[iE]->RebinX(4);
	p_HCALvsECAL[iE]->RebinY(4);

	gStyle->SetOptStat(0);

	mycECAL->cd(iE+1);
	//gPad->SetLogz(1);
	p_HCALvsECAL[iE]->GetXaxis()->SetLabelSize(0.04);
	p_HCALvsECAL[iE]->GetYaxis()->SetLabelSize(0.04);
	p_HCALvsECAL[iE]->GetXaxis()->SetTitleSize(0.04);
	p_HCALvsECAL[iE]->GetYaxis()->SetTitleSize(0.04);

	double maxX = 0;
	double minX = 10000;
	double maxY = 0;
	int binIntercept = p_HCALvsECAL[iE]->GetMaximumBin();
	int binXint,binYint,binZint;
	p_HCALvsECAL[iE]->GetBinXYZ(binIntercept,binXint,binYint,binZint);
	double Xintercept = p_HCALvsECAL[iE]->GetXaxis()->GetBinLowEdge(binXint);

	for (unsigned iX(1); iX<p_HCALvsECAL[iE]->GetNbinsX()+1; ++iX){
	  double tmpX = p_HCALvsECAL[iE]->GetXaxis()->GetBinLowEdge(iX+1);
	  for (unsigned iY(1); iY<p_HCALvsECAL[iE]->GetNbinsY()+1; ++iY){
	    double tmpY = p_HCALvsECAL[iE]->GetYaxis()->GetBinLowEdge(iY+1);
	    double content = p_HCALvsECAL[iE]->GetBinContent(iX,iY);
	    if (content>0 && tmpX>maxX) {
	      maxX = tmpX;
	    }
	    if (content>0 && tmpX<minX) {
	      minX = tmpX;
	    }
	    if (content>0 && tmpY>maxY) {
	      maxY = tmpY;
	    }
	  }
	}

	std::cout << " -- minX = " << minX << " maxX = " << maxX << " "
		  << " maxY = " << maxY << " "
		  << " Xintercept = " << Xintercept << " "
		  << std::endl;


	maxX = maxX*1.3;
	maxY = maxY*1.3;
	p_HCALvsECAL[iE]->GetXaxis()->SetRangeUser(0,maxX);
	p_HCALvsECAL[iE]->GetYaxis()->SetRangeUser(0,maxY);
	p_HCALvsECAL[iE]->GetYaxis()->SetTitleOffset(1.3);
	//p_HCALvsECAL[iE]->GetZaxis()->SetTitle("Events");
	
	//p_HCALvsECAL[iE]->Draw("colz");
	char buf[500];
	sprintf(buf,"#pi^{+} %d GeV",genEn[iE]);
	p_HCALvsECAL[iE]->SetTitle(buf);

	//p_HCALvsECAL[iE]->Draw("LEGO2z");
	p_HCALvsECAL[iE]->Draw("colz");
	TProfile *prof = p_HCALvsECAL[iE]->ProfileX();
	prof->SetMarkerStyle(23);
	prof->SetMarkerColor(1);
	prof->Draw("PEsame");
	TF1 *mypol1 = new TF1("mypol1","[0]+[1]*x",0,2000);
	mypol1->SetParameters(maxY,-1.4);
	mypol1->SetParLimits(1,-2,-1);
	double tmpRange = Xintercept-minX;
	prof->Fit("mypol1","R+","same",minX+tmpRange/5.,Xintercept-tmpRange/5.);//,maxX/20,maxX);

	TLatex lat;
	sprintf(buf,"y = #alpha + #beta #times x");

	double range = maxX;
	double xpos = p_HCALvsECAL[iE]->GetXaxis()->GetXmin()+range/10.;
	double ypos = maxY;
	lat.DrawLatex(xpos,ypos*0.9,buf);
	sprintf(buf,"#alpha = %3.3f #pm %3.3f",
		mypol1->GetParameter(0),mypol1->GetParError(0));
	lat.DrawLatex(xpos,ypos*0.8,buf);
	sprintf(buf,"#beta = %3.3f #pm %3.3f",
		mypol1->GetParameter(1),mypol1->GetParError(1));
	lat.DrawLatex(xpos,ypos*0.7,buf);
	sprintf(buf,"Chi2/NDF = %3.1f/%d = %3.1f",mypol1->GetChisquare(),mypol1->GetNDF(),mypol1->GetChisquare()/mypol1->GetNDF());
	lat.DrawLatex(xpos,ypos*0.6,buf);

	slope[iE] = mypol1->GetParameter(1);
	slopeErr[iE] = mypol1->GetParError(1);

      }//loop on energies
 
      saveName.str("");
      saveName << plotDir << "/" << pDetector;
      mycECAL->Update();
      mycECAL->Print((saveName.str()+".png").c_str());
      mycECAL->Print((saveName.str()+".pdf").c_str());

      //TGraphErrors *gr = new TGraphErrors((unsigned)nGenEn,(double*)genEn,slope,genEnErr,slopeErr);
      TGraphErrors *gr = new TGraphErrors(nGenEn,grgenEn,slope,genEnErr,slopeErr);
      myc->cd();
      myc->SetLogx(1);
      //gr->SetMinimum(-2);
      //gr->SetMaximum(0);
      gr->GetXaxis()->SetTitle("Energy (GeV");
      gr->GetYaxis()->SetTitle("Slope HCAL/ECAL");
      if (pDetector.find("ScivsSi") != pDetector.npos) gr->GetYaxis()->SetTitle("Slope Scint/Si");
      gr->SetTitle("");
      gr->Draw("AP");

      gStyle->SetOptFit(1111);
      gr->Fit("pol0","R+","same",50,2000);
      TF1 *fitslope = gr->GetFunction("pol0");

      saveName.str("");
      saveName << plotDir << "/SlopevsE_" << pDetector;
      myc->Update();
      myc->Print((saveName.str()+".png").c_str());
      myc->Print((saveName.str()+".pdf").c_str());

    }//loop on scenarios
    
  }//loop on versions
  
  return 0;


}//main
