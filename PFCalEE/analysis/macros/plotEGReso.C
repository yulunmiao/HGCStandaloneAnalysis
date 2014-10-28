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

TPad* plot_ratio(TCanvas *canv, bool up){
  canv->SetFillColor      (0);
  canv->SetBorderMode     (0);
  canv->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canv->SetLeftMargin     (0.18);
  canv->SetLeftMargin     (0.17);
  canv->SetRightMargin    (0.05);
  canv->SetTopMargin      (0.12);
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
    pad->SetBottomMargin(0.05);
    pad->SetTopMargin(0.09);
    pad->Draw();
    pad->cd();
    return pad;
  }
  else {
    pad = new TPad("lower","pad",0, 0   ,1 ,0.26);  
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.24);
    pad->Draw();
    return pad;
  }

};

int plotEGReso(){//main

  const unsigned nPu = 2;
  unsigned pu[nPu] = {0,140};//,140};

  const unsigned nS = 1;
  std::string scenario[nS] = {
    "gamma/200um/"
  };

  const unsigned eta = 17;

  const unsigned nEvtMin = 500;

  TString pSuffix = "";

  unsigned rebinReco = 8;

  bool addNoiseTerm = false;
  
  const unsigned nV = 1;
  TString version[nV] = {"12"};//,"0"};
  
  const unsigned nLayers = 30;

  const unsigned nSR = 4;//reject SR0

  const bool doVsE = true;

  std::string unit = "MIPs";
  const char* unitStr = unit.c_str();

  const unsigned MAX = 4;
  TString type[MAX];
  Float_t sigmaStoch[nPu][nSR][MAX];
  Float_t sigmaStochErr[nPu][nSR][MAX];
  Float_t sigmaConst[nPu][nSR][MAX];
  Float_t sigmaConstErr[nPu][nSR][MAX];
  Float_t sigmaNoise[nPu][nSR][MAX];
  Float_t sigmaNoiseErr[nPu][nSR][MAX];
  
  std::ostringstream saveName;
  unsigned genEn[]={20,30,40,50,60,70,80,90,100,125,150,175,200};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  unsigned rebin[20] = {4,4,4,6,6,
			6,6,8,8,10,
			10,10,100,100,100,
			6,6,6,6,6};


  unsigned nx=0,ny=0;

  if (nGenEn>12) {nx=5;ny=3;}
  else if (nGenEn > 10)
    {nx=4;ny=3;}
  else if (nGenEn > 6)
    {nx=5;ny=2;}
  else if (nGenEn > 4)
    {nx=3;ny=2;}
  else if (nGenEn > 2)
    {nx=2;ny=2;}
  else 
    {nx=nGenEn;ny=1;}
  
  //canvas so they are created only once
  TCanvas *mycE[nSR];
  for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
    std::ostringstream lName;
    lName << "mycE" << iSR;
    mycE[iSR] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    mycE[iSR]->Divide(nx,ny);
  }

  const unsigned nCanvas = 4;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
    myc[iC]->Divide(2,2);
  }

  TH1F *p_chi2overNDF = new TH1F("p_chi2overNDF",";#chi^{2}/NDF",500,0,500);
  
  //TPad *upper = plot_ratio(myc[1], true);
  //TPad *lower = plot_ratio(myc[1], false);
  //if (!upper || !lower){
  //std::cout << " Pb..." << upper << " " << lower << std::endl;
  //return 1;
  //}

  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
      
      TString plotDir = "../PLOTS/gitV00-02-09/version"+version[iV]+"/"+scenario[iS]+"/";
      
      for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
	unsigned puOption = pu[ipu];
	
       
	TH1F *p_Ereco[nGenEn][nSR];

	for (unsigned iE(0); iE<nGenEn; ++iE){
	  
	  std::cout << "- Processing energy : " << genEn[iE] 
		    << std::endl;
	  
	  TFile *inputFile = 0;
	  std::ostringstream linputStr;
	  linputStr << plotDir << "eta" << eta << "_et" << genEn[iE] << "_pu" << pu[ipu] << ".root";
	  inputFile = TFile::Open(linputStr.str().c_str());
	  if (!inputFile) {
	    std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Exiting..." << std::endl;
	    return 1;
	  }
	  else std::cout << " -- File " << inputFile->GetName() << " sucessfully opened." << std::endl;
	  
	  for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	    std::cout << " --Processing signal region: " << iSR << std::endl;
	    std::ostringstream lName;
	    lName.str("");
	    lName << "p_wgtSubtractESR" << iSR+1 ;
	    p_Ereco[iE][iSR] = (TH1F*)gDirectory->Get(lName.str().c_str());
	    if (!p_Ereco[iE][iSR]){
	      std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null." << std::endl;
	      return 1;
	    }
	    std::cout << " --- Reco E = entries " << p_Ereco[iE][iSR]->GetEntries() 
		      << " mean " << p_Ereco[iE][iSR]->GetMean() 
		      << " rms " << p_Ereco[iE][iSR]->GetRMS() 
		      << " overflows " << p_Ereco[iE][iSR]->GetBinContent(p_Ereco[iE][iSR]->GetNbinsX()+1)
		      << std::endl;
	      
	    //p_Ereco[iE]->Rebin(rebin[iE]);
	    p_Ereco[iE][iSR]->Rebin(genEn[iE]<50?rebinReco:genEn[iE]<100?2*rebinReco:4*rebinReco);
	  }//loop on SR
	}//loop on energies

      
	gStyle->SetOptStat(0);

	TGraphErrors *calibReco[nSR];
	TGraphErrors *resoReco[nSR];
	TGraphErrors *calibRecoFit[nSR];
	TGraphErrors *deltaRecoFit[nSR];
	TGraphErrors *resoRecoFit[nSR];

	for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	  std::cout << " --Processing signal region: " << iSR << std::endl;
	  TString srStr = "";
	  srStr += iSR;
	  //draw calibration curves
	  calibReco[iSR] = new TGraphErrors();
	  calibReco[iSR]->SetName("calibReco"+srStr);
	  calibReco[iSR]->SetTitle("");
	  calibReco[iSR]->SetMarkerStyle(20+iSR);
	  calibReco[iSR]->SetMarkerColor(iSR!=4? iSR+1 : iSR+2);
	  calibReco[iSR]->SetLineColor(iSR!=4? iSR+1 : iSR+2);
	  resoReco[iSR] = (TGraphErrors *) calibReco[iSR]->Clone("resoReco"+srStr);
	  calibRecoFit[iSR] = (TGraphErrors *) calibReco[iSR]->Clone("calibRecoFit"+srStr);
	  deltaRecoFit[iSR] = (TGraphErrors *) calibReco[iSR]->Clone("deltaRecoFit"+srStr);
	  resoRecoFit[iSR] = (TGraphErrors *) calibReco[iSR]->Clone("resoRecoFit"+srStr);
	  
	  //simhits
	  for (unsigned iE(0); iE<nGenEn; ++iE){

	    //skip data with too little stat
	    if (p_Ereco[iE][iSR]->GetEntries()<nEvtMin) continue;

	    std::cout << "- Processing energy : " << genEn[iE] << std::endl;
	    mycE[iSR]->cd(iE+1);
	    double eMin = p_Ereco[iE][iSR]->GetMean()-10*p_Ereco[iE][iSR]->GetRMS();
	    double eMax = p_Ereco[iE][iSR]->GetMean()+10*p_Ereco[iE][iSR]->GetRMS();
	    p_Ereco[iE][iSR]->GetXaxis()->SetRangeUser(eMin,eMax);
	    p_Ereco[iE][iSR]->Draw("PE");
	    char buf[500];
	    sprintf(buf,"Photons, E_{T}=%d GeV",genEn[iE]);
	    p_Ereco[iE][iSR]->SetTitle(buf);
	    double nRMS = iSR<2? 3 : 2;
	    p_Ereco[iE][iSR]->Fit("gaus","LR0","",
				  p_Ereco[iE][iSR]->GetMean()-nRMS*p_Ereco[iE][iSR]->GetRMS(),
				  p_Ereco[iE][iSR]->GetMean()+nRMS*p_Ereco[iE][iSR]->GetRMS());
	    
	    TF1 *fitResult = p_Ereco[iE][iSR]->GetFunction("gaus");
	    nRMS = iSR<2? 1 : 2;
	    p_Ereco[iE][iSR]->Fit("gaus","LR+","same",
				  fitResult->GetParameter(1)-nRMS*fitResult->GetParameter(2),
				  fitResult->GetParameter(1)+2*fitResult->GetParameter(2));
	    fitResult = p_Ereco[iE][iSR]->GetFunction("gaus");

	    p_chi2overNDF->Fill(fitResult->GetChisquare()/fitResult->GetNDF());

	    TLatex lat;
	    double latx = std::max(0.,p_Ereco[iE][iSR]->GetMean()-5*p_Ereco[iE][iSR]->GetRMS());
	    double laty = p_Ereco[iE][iSR]->GetMaximum();
	    sprintf(buf,"<E_{T}> = %3.3f %s",p_Ereco[iE][iSR]->GetMean(),unitStr);
	    lat.DrawLatex(latx,laty*0.9,buf);
	    sprintf(buf,"RMS = %3.3f #pm %3.1f %s",p_Ereco[iE][iSR]->GetRMS(),p_Ereco[iE][iSR]->GetRMSError(),unitStr);
	    lat.DrawLatex(latx,laty*0.8,buf);
	    sprintf(buf,"RMS/mean = %3.3f",p_Ereco[iE][iSR]->GetRMS()/p_Ereco[iE][iSR]->GetMean());
	    lat.DrawLatex(latx,laty*0.7,buf);
	    sprintf(buf,"<E_{T}fit> = %3.3f +/- %3.3f %s",fitResult->GetParameter(1),fitResult->GetParError(1),unitStr);
	    lat.DrawLatex(latx,laty*0.6,buf);
	    sprintf(buf,"RMSfit = %3.3f +/- %3.3f %s",fitResult->GetParameter(2),fitResult->GetParError(2),unitStr);
	    lat.DrawLatex(latx,laty*0.5,buf);
	    sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
	    lat.DrawLatex(latx,laty*0.4,buf);
	    
	    sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
	    lat.DrawLatex(latx,laty*0.3,buf);
      
	    Int_t np=calibReco[iSR]->GetN();
	    calibReco[iSR]->SetPoint(np,genEn[iE],p_Ereco[iE][iSR]->GetMean());
	    calibReco[iSR]->SetPointError(np,0.0,p_Ereco[iE][iSR]->GetMeanError());
	    resoReco[iSR]->SetPoint(np,doVsE?genEn[iE] :1/sqrt(genEn[iE]),p_Ereco[iE][iSR]->GetRMS()/p_Ereco[iE][iSR]->GetMean());
	    resoReco[iSR]->SetPointError(np,0,p_Ereco[iE][iSR]->GetRMSError()/p_Ereco[iE][iSR]->GetMean());
	    calibRecoFit[iSR]->SetPoint(np,genEn[iE],fitResult->GetParameter(1));
	    calibRecoFit[iSR]->SetPointError(np,0.0,fitResult->GetParError(1));
	    deltaRecoFit[iSR]->SetPoint(np,genEn[iE],( ((fitResult->GetParameter(1)-0)/200)-genEn[iE])/genEn[iE]);
	    deltaRecoFit[iSR]->SetPointError(np,0.0,fitResult->GetParError(1)/200*1./genEn[iE]);
	    resoRecoFit[iSR]->SetPoint(np,doVsE?genEn[iE] :1/sqrt(genEn[iE]),fitResult->GetParameter(2)/fitResult->GetParameter(1));
	    double errFit = fitResult->GetParameter(2)/fitResult->GetParameter(1)*sqrt(pow(fitResult->GetParError(2)/fitResult->GetParameter(2),2)+pow(fitResult->GetParError(1)/fitResult->GetParameter(1),2));
	    resoRecoFit[iSR]->SetPointError(np,0,errFit);
      
	  }//loop on energies

	  saveName.str("");
	  saveName << plotDir << "/Ereco_eta" << eta << "_pu" << puOption << "_SR" << iSR;
	  mycE[iSR]->Update();
	  mycE[iSR]->Print((saveName.str()+".png").c_str());
	  mycE[iSR]->Print((saveName.str()+".pdf").c_str());
	  /*
	}//loop on signal regions

	myc[4]->cd();
	p_chi2overNDF->Draw();
	
	myc[4]->Update();
	return 1;

      }//loop on pu options

      
    }//loop on scenarios
    
  }//loop on versions
	  */
	
	  //return 1;

	//draw calib

	const unsigned imax = 4;

	for(unsigned i=0; i<imax ; i++)
	  {
	    //if (i==1) upper->cd();
	    //else 
	    myc[i]->cd(iSR+1);
	    type[i] = i==0 ? "calibReco" : i==1 ? "calibRecoFit" : i==2 ? "resoReco" : "resoRecoFit";
	    
	    if (nPu > 1){
	      type[i] += "_pu";
	      type[i] += puOption;
	    }

	    std::cout << "- Processing type : " << type[i] << std::endl;

	    TGraphErrors * gr =( i==0 ? calibReco[iSR] : i==1 ? calibRecoFit[iSR] : i==2 ? resoReco [iSR]: resoRecoFit[iSR]);

	    gr->GetXaxis()->SetLabelSize(0.06);
	    gr->GetXaxis()->SetTitleSize(0.06);
	    gr->GetYaxis()->SetLabelSize(0.06);
	    gr->GetYaxis()->SetTitleSize(0.06);
	    gr->GetXaxis()->SetTitleOffset(0.7);
	    gr->GetYaxis()->SetTitleOffset(0.8);

	    gr->Draw("ap");
	    //gr->GetYaxis()->SetRangeUser(0,i%4<2?100 : 0.2);
	    

	    if(i<2){
	      if (i!=1) gr->GetXaxis()->SetTitle("Beam E_{T} [GeV]");
	      else gr->GetXaxis()->SetTitle("");
	      gr->GetYaxis()->SetTitle("Average energy deposited ["+TString(unitStr)+"]"); 
	    }
	    else { 
	      gr->GetXaxis()->SetTitle(doVsE?"Beam E_{T} [GeV]" :"1/#sqrt{Beam E_{T}} [1/#sqrt{GeV}]"); 
	      gr->GetYaxis()->SetTitle("Relative energy resolution");     
	    }
	    char buf[500];
	    if(i<2){
	      TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
	      //TF1 *fitFunc=new TF1("calib","[0]+[1]*x",20,100);
	      fitFunc->SetLineColor(iSR!=4? iSR+1 : iSR+2);
	      gr->Fit(fitFunc,"RME");
	      TLatex lat;
	      lat.SetTextColor(iSR!=4? iSR+1 : iSR+2);
	      sprintf(buf,"<E> #propto a + b #times E_{T} ");
	      if (i<2) lat.DrawLatex(genEn[0],gr->GetYaxis()->GetXmax()*0.9,buf);
	      sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0),unitStr);
	      lat.DrawLatex(genEn[0]+i/4*50,gr->GetYaxis()->GetXmax()*(0.8-i/2*0.25),buf);
	      sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1),unitStr);
	      lat.DrawLatex(genEn[0]+i/4*50,gr->GetYaxis()->GetXmax()*(0.7-i/2*0.25),buf);
	      sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
	      lat.DrawLatex(genEn[0]+i/4*50,gr->GetYaxis()->GetXmax()*(0.6-i/2*0.25),buf);
	      //draw deltaE/E vs E
	      /*if (i==1){
		lower->cd();
		gPad->SetLogx(0);
		gPad->SetGridx(1);
		gPad->SetGridy(1);
		TGraphErrors * grDelta = deltaRecoFit[iSR];
		grDelta->SetTitle("");
		grDelta->SetMinimum(-0.1);
		grDelta->SetMaximum(0.1);
		grDelta->GetXaxis()->SetLabelSize(0.15);
		grDelta->GetXaxis()->SetTitleSize(0.15);
		grDelta->GetYaxis()->SetLabelSize(0.12);
		grDelta->GetYaxis()->SetTitleSize(0.15);
		grDelta->GetXaxis()->SetTitleOffset(0.5);
		grDelta->GetYaxis()->SetTitleOffset(0.3);

		grDelta->Draw(i==1? "ap" : "p");
		//grDelta->GetYaxis()->SetRangeUser(0,grDelta->GetYaxis()->GetXmax());
		grDelta->GetXaxis()->SetTitle("Beam energy [GeV]");
		grDelta->GetYaxis()->SetTitle("(#Delta E)/E");
		}*/

	    }
	    else
	      {
		//TF1 *fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1]+[2]*x*x*x*x)",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		TF1 *fitFunc2;
		if (doVsE){
		  fitFunc2 =new TF1("reso","sqrt([0]*[0]/x+[1]*[1])",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		  if (addNoiseTerm) fitFunc2 =new TF1("reso","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		}
		else {
		  fitFunc2 =new TF1("reso","sqrt([0]*[0]*x*x+[1]*[1])",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		  if (addNoiseTerm) fitFunc2 =new TF1("reso","sqrt([0]*[0]*x*x+[1]*[1]+[2]*[2]*x*x*x*x)",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		}
		//if (i<4) fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1])",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		//else fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1])",gr->GetXaxis()->GetXmin(),0.2);//gr->GetXaxis()->GetXmax());
		
		fitFunc2->SetParameter(0,0.215);
		fitFunc2->SetParLimits(0,0,1);
		fitFunc2->SetParameter(1,0.01);
		fitFunc2->SetParLimits(1,0,1);
		if (addNoiseTerm) {
		  fitFunc2->SetParLimits(2,0,2);
		  fitFunc2->FixParameter(2,0.06);
		}
		fitFunc2->SetLineColor(iSR!=4? iSR+1 : iSR+2);

		gr->Fit(fitFunc2,"RME");
		sigmaStoch[ipu][iSR][i] = (fitFunc2->GetParameter(0));
		sigmaStochErr[ipu][iSR][i] = fitFunc2->GetParError(0)/(2*sigmaStoch[ipu][iSR][i]);
		sigmaConst[ipu][iSR][i] = (fitFunc2->GetParameter(1));
		sigmaConstErr[ipu][iSR][i] = fitFunc2->GetParError(1)/(2*sigmaConst[ipu][iSR][i]);
		if (addNoiseTerm) {
		  sigmaNoise[ipu][iSR][i] = (fitFunc2->GetParameter(2));
		  sigmaNoiseErr[ipu][iSR][i] = fitFunc2->GetParError(2)/(2*sigmaNoise[ipu][iSR][i]);
		}
		TLatex lat;
		lat.SetTextColor(iSR!=4? iSR+1 : iSR+2);
		if (addNoiseTerm) sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E_{T}}} #oplus c #oplus #frac{n}{E_{T}}");
		else sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E_{T}}} #oplus c");
		double Emin = doVsE?40 : 1/sqrt(genEn[nGenEn-1]);
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax(),buf);
		sprintf(buf,"s=%3.3f #pm %3.3f",sigmaStoch[ipu][iSR][i],sigmaStochErr[ipu][iSR][i]);
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.9-i/6*0.25),buf);
		sprintf(buf,"c=%3.3f #pm %3.3f",sigmaConst[ipu][iSR][i],sigmaConstErr[ipu][iSR][i]);
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.82-i/6*0.25),buf);
		sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc2->GetChisquare(),fitFunc2->GetNDF(),fitFunc2->GetChisquare()/fitFunc2->GetNDF());
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.66-i/6*0.25),buf);
		//if (i>3){
		if (addNoiseTerm) {
		  sprintf(buf,"n=%3.3f #pm %3.3f",sigmaNoise[ipu][iSR][i],sigmaNoiseErr[ipu][iSR][i]);
		  lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.74-i/6*0.25),buf);
		}
	      }

	  }
	}//loop on SR

	for(unsigned i=0; i<4 ; i++){
	  myc[i]->Update();
	  if (i>1 && doVsE){
	    myc[i]->Print(plotDir+"/"+type[i]+"_vsET.pdf");
	    myc[i]->Print(plotDir+"/"+type[i]+"_vsET.png");
	  }
	  else {
	    myc[i]->Print(plotDir+"/"+type[i]+".pdf");
	    myc[i]->Print(plotDir+"/"+type[i]+".png");
	  }
	}
      }//loop on pu options
      
      
    }//loop on scenarios
    
  }//loop on versions
  
  std::cout << "scenario & type & sigmaStoch & sigmaConst";
  if (addNoiseTerm) std::cout << " & sigmaNoise";
  std::cout << " \\\\ \n" 
	    << "\\hline\n";
  for (unsigned iS(0); iS<nSR;++iS){
    for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
      std::cout << "pu " << pu[ipu] << std::endl;
      for (unsigned i(2); i<MAX;++i){
	std::cout << "SR" << iS << " & " << type[i] << " & " <<  std::setprecision(3)
		  << "$" << sigmaStoch[ipu][iS][i] << "\\pm" << sigmaStochErr[ipu][iS][i] << "$ & "
		  << "$" << sigmaConst[ipu][iS][i] << "\\pm" << sigmaConstErr[ipu][iS][i] << "$";
	if (addNoiseTerm) std::cout << " & $" << sigmaNoise[ipu][iS][i] << "\\pm" << sigmaNoiseErr[ipu][iS][i] << "$";
	std::cout << "\n";
      }
      std::cout <<"\\hline\n";
    }
  }


  return 0;
  
}//main
