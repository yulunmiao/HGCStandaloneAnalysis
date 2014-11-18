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

struct FitResult{
  double chi2;
  unsigned ndf;
  double mean;
  double sigma;
  double meanerr;
  double sigmaerr;
};

TPad* plot_ratio(TPad *canv, bool up){
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


unsigned fitEnergy(TH1F *hist,
		   TPad *pad,
		   std::string unitStr,
		   FitResult & lres,
		   unsigned isr){
  
  pad->cd();
  //double eMin = hist->GetMean()-5*hist->GetRMS();
  //double eMax = hist->GetMean()+5*hist->GetRMS();
  //hist->GetXaxis()->SetRangeUser(eMin,eMax);
  hist->Draw("PE");


  double nRMSm = isr<1? 1 : 2;
  double nRMSp = 2;
  
  TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  fitResult->SetParameters(hist->GetBinContent(hist->GetMaximumBin()),
			   hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin()),
			   hist->GetRMS());

  std::cout << " Initial params: "  << fitResult->GetParameter(0) << " "<< fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
	    << std::endl;


  int status = hist->Fit("fitResult","LR0+","",
			 fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
			 fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
  
  
  std::cout << " First fit: " << status << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
	    << std::endl;

  if (status != 0 || fitResult->GetChisquare()/fitResult->GetNDF()>5){
    std::cout << " -- Bad fit ! Try again..." << std::endl;
    status = hist->Fit("fitResult","LR0+","",
		       fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
		       fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
    
    std::cout << " Second fit: " << status << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
	      << std::endl;
  }
  
  std::cout << " Final fit: " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
	    << std::endl;
  
  fitResult->SetLineColor(2);
  fitResult->Draw("same");
  
  if (status != 0) {
    std::cout << " ERROR! Fit failed ! Please have a look. I stop here..." << std::endl;
    //totalE for pu140 is expected to be pathological :/
    if (isr!=7) return 1;
  }

  char buf[500];
  TLatex lat;
  double latx = hist->GetXaxis()->GetXmin()+(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin())/20.;
  double laty = hist->GetMaximum();
  sprintf(buf,"<E_{T}fit> = %3.3f +/- %3.3f %s",fitResult->GetParameter(1),fitResult->GetParError(1),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.9,buf);
  sprintf(buf,"RMSfit = %3.3f +/- %3.3f %s",fitResult->GetParameter(2),fitResult->GetParError(2),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.8,buf);
  sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
  lat.DrawLatex(latx,laty*0.7,buf);
  
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
  lat.DrawLatex(latx,laty*0.6,buf);
  
  lres.chi2 = fitResult->GetChisquare();
  lres.ndf = fitResult->GetNDF();
  lres.mean = fitResult->GetParameter(1);
  lres.meanerr = fitResult->GetParError(1);
  lres.sigma = fitResult->GetParameter(2);
  lres.sigmaerr = fitResult->GetParError(2);

  return 0;
};

TPad* plotCalibration(TGraphErrors *gr,TPad *pad,bool doRatio, TGraphErrors *grDelta,std::string unit, double & calib,double & calibErr, double & offset, double & offsetErr){

  TPad *upper = 0;
  TPad *lower = 0;

  if (!doRatio) pad->cd();
  else {
    pad->Clear();
    upper = plot_ratio(pad, true);
    lower = plot_ratio(pad, false);
    upper->cd();
  }
  gr->GetXaxis()->SetLabelSize(0.0);
  gr->GetXaxis()->SetTitleSize(0.0);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(7);
  gr->GetYaxis()->SetTitleOffset(0.9);
  
  gr->Draw("ap");

  if (!doRatio) gr->GetXaxis()->SetTitle("E_{T} (GeV)");
  else gr->GetXaxis()->SetTitle("");
  
  gr->GetYaxis()->SetTitle(("Average energy deposited ("+unit+")").c_str()); 
  char buf[500];
  TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
  fitFunc->SetLineColor(6);
  gr->Fit(fitFunc,"RME");
  TLatex lat;
  lat.SetTextColor(6);
  sprintf(buf,"<E> #propto a + b #times E_{T} ");
  lat.DrawLatex(30,gr->GetYaxis()->GetXmax()*0.9,buf);
  sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0),unit.c_str());
  lat.DrawLatex(30,gr->GetYaxis()->GetXmax()*0.8,buf);
  sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1),unit.c_str());
  lat.DrawLatex(30,gr->GetYaxis()->GetXmax()*0.7,buf);
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatex(30,gr->GetYaxis()->GetXmax()*0.6,buf);

  calib = fitFunc->GetParameter(1);
  offset = fitFunc->GetParameter(0);
  calibErr = fitFunc->GetParError(1);
  offsetErr = fitFunc->GetParError(0);

  if (doRatio){
    //draw deltaE/E vs E
    lower->cd();
    gPad->SetLogx(0);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    double loffset = fitFunc->GetParameter(0);
    double lslope = fitFunc->GetParameter(1);
    double range = 1;
    if (unit=="GeV") {
      loffset=0;
      lslope=1;
      range = 10;
    }

    //fill delta
    for (int ip(0);ip<gr->GetN();++ip){
      double x=0;
      double y=0;
      gr->GetPoint(ip,x,y);
      grDelta->SetPoint(ip,x,(y-loffset)/lslope-x);
      double err = gr->GetErrorY(ip)/lslope;
      grDelta->SetPointError(ip,0,err);
      std::cout << "Calib " << ip << " Egen=" << x << " Erec=" << y << " delta=" << (y-loffset)/lslope-x << std::endl;
    }
    grDelta->SetTitle("");
    grDelta->SetMinimum(-1.*range);
    grDelta->SetMaximum(range);
    grDelta->GetXaxis()->SetLabelSize(0.15);
    grDelta->GetXaxis()->SetTitleSize(0.15);
    grDelta->GetYaxis()->SetLabelSize(0.12);
    grDelta->GetYaxis()->SetTitleSize(0.15);
    grDelta->GetXaxis()->SetTitleOffset(0.5);
    grDelta->GetYaxis()->SetTitleOffset(0.3);
    
    grDelta->Draw("ap");
    //grDelta->GetYaxis()->SetRangeUser(0,grDelta->GetYaxis()->GetXmax());
    grDelta->GetXaxis()->SetTitle("E_{T} (GeV)");
    grDelta->GetYaxis()->SetTitle("(#Delta E)/E");

    TLine *line = new TLine(grDelta->GetXaxis()->GetXmin(),0,grDelta->GetXaxis()->GetXmax(),0);
    line->SetLineColor(2);//kYellow+4);
    line->Draw();
    

    lower->Update();

  }

  return upper;
};

void plotResolution(TGraphErrors *gr,TPad *pad,
		    const unsigned ipu,
		    const double & stoch0,
		    const double & const0,
		    const double & noise0,
		    double & stoch,double & stochErr, 
		    double & constant, double & constErr,
		    double & noise,double & noiseErr){

  pad->cd();
  gr->GetXaxis()->SetLabelSize(0.06);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(0.7);
  gr->GetYaxis()->SetTitleOffset(0.8);
  
  gr->Draw("ap");
  gr->GetXaxis()->SetTitle("E_{T} (GeV)");
  gr->GetYaxis()->SetTitle("#sigma/E");

  TF1 *fitFunc =new TF1("reso","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());

  fitFunc->SetParameter(0,stoch0);
  fitFunc->SetParLimits(0,0,1);
  fitFunc->SetParameter(1,const0);
  fitFunc->SetParLimits(1,0,1);
  fitFunc->SetParameter(2,noise0);
  
  if (ipu<2) fitFunc->FixParameter(2,0.0);
  if (ipu==2) fitFunc->FixParameter(1,const0);

  int status = gr->Fit(fitFunc,"RME0");

  if (fitFunc->GetChisquare()/fitFunc->GetNDF()>2){
    status = gr->Fit(fitFunc,"RME0","",20,100);
  }

  fitFunc->SetLineColor(6);
  fitFunc->Draw("same");

  stoch = fitFunc->GetParameter(0);
  stochErr = fitFunc->GetParError(0);
  constant = fitFunc->GetParameter(1);
  constErr = fitFunc->GetParError(1);
  noise = fitFunc->GetParameter(2);
  noiseErr = fitFunc->GetParError(2);

  char buf[500];
  TLatex lat;
  lat.SetTextColor(6);
  sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E_{T}}} #oplus c #oplus #frac{n}{E_{T}}");

  double Emin = 30;
  lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*0.9,buf);
  sprintf(buf,"s=%3.2f #pm %3.2f",stoch,stochErr);
  lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*0.8,buf);
  sprintf(buf,"c=%3.2f #pm %3.2f",constant,constErr);
  lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*0.7,buf);
  sprintf(buf,"n=%3.2f #pm %3.2f",noise,noiseErr);
  lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*0.6,buf);
  sprintf(buf,"status = %d, #chi^{2}/N = %3.1f/%d = %3.1f",status,fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*0.5,buf);
  
};

void retrievePuSigma(TTree *atree, TTree *atreePu, 
		     std::vector<double> & aSigmaVec, 
		     const std::string plotdir,
		     double * calib){

  unsigned nLayers = 30;
  unsigned nSR=5;
  aSigmaVec.resize(nSR+2,0.5);
  if (!atree || !atreePu) return;
  else std::cout << " Trees found." << std::endl;

  unsigned evtIndex = 0;
  unsigned evtIndexPu = 0;
  atree->SetBranchAddress("eventIndex",&evtIndex);
  atreePu->SetBranchAddress("eventIndex",&evtIndexPu);

  std::vector<std::vector<double> > energySR[2];
  //  std::vector<std::vector<double> > subtractedenergySR[2];
 
  std::vector<double> emptyvec;
  emptyvec.resize(nSR,0);
  energySR[0].resize(nLayers,emptyvec);
  energySR[1].resize(nLayers,emptyvec);
  //subtractedenergySR[0].resize(nLayers,emptyvec);

  std::ostringstream label;
  for (unsigned iL(0); iL<nLayers;++iL){
    for (unsigned iSR(0);iSR<nSR;++iSR){
      label.str("");
      label << "energy_" << iL << "_SR" << iSR;
      atree->SetBranchAddress(label.str().c_str(),&energySR[0][iL][iSR]);
      atreePu->SetBranchAddress(label.str().c_str(),&energySR[1][iL][iSR]);
      //label.str("");
      //label << "subtractedenergy_" << iL << "_SR" << iSR;
      //outtree_->SetBranchAddress(label.str().c_str(),&subtractedenergySR[iL][iSR]);
    }
  }
  int nEvtsPu = atreePu->GetEntries();
  int nEvts = atree->GetEntries();
  if (nEvts!=nEvtsPu) {
    std::cout << " -- pb, not all events found! nEvts = " << nEvts << " nEvtsPu=" << nEvtsPu << std::endl;
    nEvts = std::min(nEvts,nEvtsPu);
    return;
  }

  TH1F *p_sigma[nSR+2];
  for (unsigned iSR(0); iSR<nSR+2;++iSR){
    label.str("");
    label << "p_sigma_" << iSR;
    p_sigma[iSR] = new TH1F(label.str().c_str(),";PuE-E (GeV)",1000,0,200);
    p_sigma[iSR]->StatOverflows();
  }

  for (int ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    atree->GetEntry(ievt);
    atreePu->GetEntry(ievt);

    if (evtIndex != evtIndexPu) {
      //std::cout << " -- Different events found: " << evtIndex << " " << evtIndexPu << ". Skipping." << std::endl;
      continue;
    }

    double puE[nSR+2];
    double E[nSR+2];
    for (unsigned iSR(0); iSR<nSR+2;++iSR){//loop on signal region
      puE[iSR]=0;
      E[iSR]=0;
      for (unsigned iL(0);iL<nLayers;++iL){
	if (iSR<nSR){
	  E[iSR]   += energySR[0][iL][iSR];
	  puE[iSR] += energySR[1][iL][iSR];
	}
	else if (iSR==nSR) {
	  if (iL<5) {
	    E[iSR] += energySR[0][iL][0];
	    puE[iSR] += energySR[1][iL][0];
	  }
	  else if (iL<10) {
	    E[iSR] += energySR[0][iL][1];
	    puE[iSR] += energySR[1][iL][1];
	      }
	  else if (iL<15) {
	    E[iSR] += energySR[0][iL][2];
	    puE[iSR] += energySR[1][iL][2];
	  } 
	  else if (iL<20) {
	    E[iSR] += energySR[0][iL][3];
	    puE[iSR] += energySR[1][iL][3];
	  } 
	  else  {
	    E[iSR] += energySR[0][iL][4];
	    puE[iSR] += energySR[1][iL][4];
	  }
	}
	else if (iSR==nSR+1) {
	  if (iL<5) {
	    E[iSR] += energySR[0][iL][0];
	    puE[iSR] += energySR[1][iL][0];
	  }
	  else if (iL<12) {
	    E[iSR] += energySR[0][iL][2];
	    puE[iSR] += energySR[1][iL][2];
	  }
	  else  {
	    E[iSR] += energySR[0][iL][4];
	    puE[iSR] += energySR[1][iL][4];
	  } 
	}
      }//loop on layers
      p_sigma[iSR]->Fill((puE[iSR]-E[iSR])/calib[iSR]);

    }//loop on SR

  }//loop on events
  TCanvas *mycS = new TCanvas("mycS","sigma Pu",1500,1000);
  mycS->Divide(4,2);
  gStyle->SetOptStat("eMRuo");
  for (unsigned iSR(0); iSR<nSR+2;++iSR){//loop on signal region
    aSigmaVec[iSR] = p_sigma[iSR]->GetRMS();
    mycS->cd(iSR+1);
    p_sigma[iSR]->Draw();
  }
  
  mycS->Update();
  mycS->Print((plotdir+"_PuSubtractedE.pdf").c_str());

};

int plotEGReso(){//main

  SetTdrStyle();

  const unsigned nPu = 3;
  unsigned pu[nPu] = {0,0,140};

  const unsigned nS = 1;
  std::string scenario[nS] = {
    "gamma/200um/"
  };

  const unsigned neta = 2;//7;
  unsigned eta[neta]={17,25};//19,21,23,25,27,29};

  double etaval[neta];
  double etaerr[neta];

  const unsigned nEvtMin = 150;

  TString pSuffix = "";

  
  const unsigned nV = 1;
  TString version[nV] = {"12"};//,"0"};
  
  const unsigned nLayers = 30;

  const unsigned nSR = 8;
  double srval[nSR];
  double srerr[nSR];
  for (unsigned iSR(0); iSR<nSR;++iSR){
    srval[iSR] = iSR*1.;
    srerr[iSR] = 0.;
  }

  double calib[nPu][neta][nSR];
  double calibErr[nPu][neta][nSR];
  double offset[nPu][neta][nSR];
  double offsetErr[nPu][neta][nSR];

  double sigmaStoch[nPu][neta][nSR];
  double sigmaStochErr[nPu][neta][nSR];
  double sigmaConst[nPu][neta][nSR];
  double sigmaConstErr[nPu][neta][nSR];
  double sigmaNoise[nPu][neta][nSR];
  double sigmaNoiseErr[nPu][neta][nSR];
  
  std::ostringstream saveName;

  //unsigned genEnAll[]={3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};
  unsigned genEnAll[]={20};//,40,60,80,100};
  const unsigned nGenEnAll=sizeof(genEnAll)/sizeof(unsigned);

  //canvas so they are created only once
  TCanvas *mycE[nGenEnAll];
  for (unsigned iE(0); iE<nGenEnAll;++iE){
    std::ostringstream lName;
    lName << "mycE" << genEnAll[iE];
    mycE[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    mycE[iE]->Divide(4,2);
  }

  const unsigned nCanvas = 2;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    myc[iC]->Divide(4,2);
  }

  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
      
      TString plotDir = "../PLOTS/gitV00-02-12/version"+version[iV]+"/"+scenario[iS]+"/";
      TTree *ltree[neta][nPu][nGenEnAll];
      
      for (unsigned ieta(0);ieta<neta;++ieta){//loop on eta
	
	etaval[ieta] = eta[ieta]/10.;
	etaerr[ieta] = 0;
	
	for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
	  unsigned puOption = pu[ipu];

	  std::string unit = "MIPS";
	  if (ipu>0) unit = "GeV";

	  //identify valid energy values
	  bool skip[nGenEnAll];
	  unsigned nValid = 0;
	  for (unsigned iE(0); iE<nGenEnAll; ++iE){	  
	    skip[iE] = false;
	    TFile *inputFile = 0;
	    std::ostringstream linputStr;
	    linputStr << plotDir << "eta" << eta[ieta] << "_et" << genEnAll[iE] << "_pu" << pu[ipu] << pSuffix << ".root";
	    inputFile = TFile::Open(linputStr.str().c_str());
	    if (!inputFile) {
	      std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Skipping..." << std::endl;
	      //	    return 1;
	      skip[iE] = true;
	    }
	    else {
	      inputFile->cd("Energies");
	      ltree[ieta][ipu][iE] = (TTree*)gDirectory->Get("Ereso");
	      
	      if (!ltree[ieta][ipu][iE]){
		std::cout << " -- File " << inputFile->GetName() << " sucessfully opened but tree Ereso not found! Skipping." << std::endl;
		skip[iE] = true;
	      } else { 
	      std::cout << " -- File " << inputFile->GetName() << " sucessfully opened and tree found." << std::endl;
	      nValid++;
	      }
	    }
	  }
	  
	  unsigned newidx = 0;
	  const unsigned nGenEn = nValid;
	  unsigned genEn[nGenEn];
	  for (unsigned iE(0); iE<nGenEnAll; ++iE){	  
	    if (!skip[iE]) {
	      genEn[newidx]=genEnAll[iE];
	      newidx++;
	    }
	  }
	  
	  TH1F *p_Ereco[nGenEn][nSR];
	  TGraphErrors *calibRecoFit[nSR];
	  TGraphErrors *calibRecoDelta[nSR];
	  TGraphErrors *resoRecoFit[nSR];

	  //draw calibration curves
	  for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	    TString srStr = "";
	    srStr += iSR;
	    calibRecoFit[iSR] = new TGraphErrors();
	    calibRecoFit[iSR]->SetName("calibRecoFit"+srStr);
	    calibRecoFit[iSR]->SetTitle("");
	    calibRecoFit[iSR]->SetMarkerStyle(20);
	    calibRecoFit[iSR]->SetMarkerColor(1);
	    calibRecoFit[iSR]->SetLineColor(1);
	    calibRecoDelta[iSR] = (TGraphErrors *) calibRecoFit[iSR]->Clone("calibRecoDelta"+srStr);
	    resoRecoFit[iSR] = (TGraphErrors *) calibRecoFit[iSR]->Clone("resoRecoFit"+srStr);
	  }

	  //get calib and offset from 0 pu file for each SR.
	  gStyle->SetOptStat(0);
	  gStyle->SetOptFit(0);

	  for (unsigned iE(0); iE<nGenEn; ++iE){
	    
	    std::cout << "- Processing energy : " << genEn[iE] 
		      << std::endl;
	    
	    TFile *inputFile = 0;
	    std::ostringstream linputStr;
	    linputStr << plotDir << "eta" << eta[ieta] << "_et" << genEn[iE] << "_pu" << pu[ipu] << pSuffix << ".root";
	    inputFile = TFile::Open(linputStr.str().c_str());
	    inputFile->cd("Energies");
	    ltree[ieta][ipu][iE] = (TTree*)gDirectory->Get("Ereso");

	    std::cout << " -- Tree entries for eta=" << eta[ieta] << " pu=" << pu[ipu] << " : " << ltree[ieta][ipu][iE]->GetEntries() << std::endl;

	    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	      std::cout << " --Processing signal region: " << iSR << std::endl;

	      mycE[iE]->cd(iSR+1);

	      std::ostringstream lName;
	      lName.str("");
	      if (ipu>0) lName << "(";
	      if (iSR<5){
		for (unsigned iL(0);iL<nLayers;++iL){	      
		  if (iL==0) lName << "subtractedenergy_" << iL << "_SR" << iSR ;
		  else lName << "+subtractedenergy_" << iL << "_SR" << iSR ;
		}
	      }
	      else if (iSR==5) {
		for (unsigned iL(0);iL<nLayers;++iL){	      
		  if (iL==0) lName << "subtractedenergy_" << iL << "_SR0" ;
		  else if (iL<5) lName << "+subtractedenergy_" << iL << "_SR0";
		  else if (iL<10) lName << "+subtractedenergy_" << iL << "_SR1";
		  else if (iL<15) lName << "+subtractedenergy_" << iL << "_SR2";
		  else if (iL<20) lName << "+subtractedenergy_" << iL << "_SR3";
		  else lName << "+subtractedenergy_" << iL << "_SR4";
		}
	      }
	      else if (iSR==6) {
		for (unsigned iL(0);iL<nLayers;++iL){	      
		  if (iL==0) lName << "subtractedenergy_" << iL << "_SR0" ;
		  else if (iL<5) lName << "+subtractedenergy_" << iL << "_SR0";
		  else if (iL<12) lName << "+subtractedenergy_" << iL << "_SR2";
		  else lName << "+subtractedenergy_" << iL << "_SR4";
		}
	      }
	      else lName << "wgtEtotal";
	      if (ipu>0) lName << " - " << offset[0][ieta][iSR] << ")/" << calib[0][ieta][iSR];
	      ltree[ieta][ipu][iE]->Draw(lName.str().c_str(),"","");
	      lName.str("");
	      lName << "energy" << genEn[iE] << "_SR" << iSR ;
	      p_Ereco[iE][iSR] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str()); // 1D
	      if (!p_Ereco[iE][iSR]){
		std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null." << std::endl;
		return 1;
	      }
	      std::cout << " --- Reco E = entries " << p_Ereco[iE][iSR]->GetEntries() 
			<< " mean " << p_Ereco[iE][iSR]->GetMean() 
			<< " rms " << p_Ereco[iE][iSR]->GetRMS() 
			<< " overflows " << p_Ereco[iE][iSR]->GetBinContent(p_Ereco[iE][iSR]->GetNbinsX()+1)
			<< std::endl;
	      
	      
	      p_Ereco[iE][iSR]->SetTitle((";E ("+unit+");events").c_str());

	      //take min 20 bins
	      if(p_Ereco[iE][iSR]->GetNbinsX()>40) p_Ereco[iE][iSR]->Rebin(2);

	      //skip data with too little stat
	      if (p_Ereco[iE][iSR]->GetEntries()<nEvtMin) {
		gPad->Clear();
		continue;
	      }

	      TPad *lpad = (TPad*)(mycE[iE]->cd(iSR+1));
	      FitResult lres;
	      if (fitEnergy(p_Ereco[iE][iSR],lpad,unit,lres,iSR)!=0) return 1;
	      lpad->cd();
	      char buf[500];
	      sprintf(buf,"#gamma E_{T}=%d GeV + PU %d",genEn[iE],pu[ipu]);
	      TLatex lat;
	      lat.SetTextSize(0.05);
	      lat.DrawLatexNDC(0.25,0.965,buf);
	      sprintf(buf,"#eta=%3.1f, SR %d",etaval[ieta],iSR);
	      lat.SetTextSize(0.06);
	      lat.DrawLatexNDC(0.15,0.87,buf);
	      if (iSR==4) lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

	      Int_t np=calibRecoFit[iSR]->GetN();
	      calibRecoFit[iSR]->SetPoint(np,genEn[iE],lres.mean);
	      calibRecoFit[iSR]->SetPointError(np,0.0,lres.meanerr);

	      double reso = fabs(lres.sigma/lres.mean);
	      resoRecoFit[iSR]->SetPoint(np,genEn[iE],reso);
	      double errFit = reso*sqrt(pow(lres.sigmaerr/lres.sigma,2)+pow(lres.meanerr/lres.mean,2));
	      resoRecoFit[iSR]->SetPointError(np,0,errFit);


	    }//loop on SR

	    saveName.str("");
	    saveName << plotDir << "/Ereco_eta" << eta[ieta] << "_pu" << puOption;
	    if (ipu==0) saveName << "raw";
	    saveName << "_E" << genEn[iE] << pSuffix;
	    mycE[iE]->Update();
	    mycE[iE]->Print((saveName.str().c_str()+pSuffix)+".pdf");
	    

	    std::vector<double> sigmaVec;
	    if (ipu>1) {
	      std::ostringstream lsave;
	      lsave << plotDir << "/eta"<< eta[ieta] << "_pu" << puOption << pSuffix;
	      retrievePuSigma(ltree[ieta][1][iE], ltree[ieta][ipu][iE], 
			      sigmaVec, 
			      lsave.str(),
			      calib[0][ieta]);
	    }


	  }//loop on energies


	  //plot and fit calib
	  for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	    TPad *lpad = (TPad*)(myc[0]->cd(iSR+1));
	    TPad *upper = plotCalibration(calibRecoFit[iSR],lpad,
					  true,calibRecoDelta[iSR],
					  unit,
					  calib[ipu][ieta][iSR],
					  calibErr[ipu][ieta][iSR],
					  offset[ipu][ieta][iSR],
					  offsetErr[ipu][ieta][iSR]);
	    upper->cd();
	    char buf[500];
	    sprintf(buf,"#gamma #eta=%3.1f + PU %d",etaval[ieta],pu[ipu]);
	    TLatex lat;
	    lat.SetTextSize(0.07);
	    lat.DrawLatexNDC(0.4,0.15,buf);
	    sprintf(buf,"SR %d",iSR);
	    lat.DrawLatexNDC(0.6,0.22,buf);
	    if (iSR==0) lat.DrawLatexNDC(0.01,0.95,"HGCAL G4 standalone");
	    
	  }

	  myc[0]->Update();
	  std::ostringstream lsave;
	  lsave << plotDir << "/";
	  if (ipu==0) lsave << "CalibMipToGeV";
	  else lsave << "Calib";
	  lsave << "_eta" << eta[ieta] << "_pu" << puOption << pSuffix;
	  myc[0]->Print((lsave.str()+".pdf").c_str());

	  //plot reso
	  for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	    TPad *lpad = (TPad*)(myc[1]->cd(iSR+1));

	    double stoch0 = pu[ipu]==0? 0.14 : sigmaStoch[1][ieta][iSR];
	    double const0 = pu[ipu]==0? 0.01 : sigmaConst[1][ieta][iSR];
	    double noise0 = pu[ipu]==0? 0 : 0.5;//sigmaVec[iSR];

	    plotResolution(resoRecoFit[iSR],lpad,
			   ipu,
			   stoch0,const0,noise0,
			   sigmaStoch[ipu][ieta][iSR],
			   sigmaStochErr[ipu][ieta][iSR],
			   sigmaConst[ipu][ieta][iSR],
			   sigmaConstErr[ipu][ieta][iSR],
			   sigmaNoise[ipu][ieta][iSR],
			   sigmaNoiseErr[ipu][ieta][iSR]);
	    lpad->cd();
	    char buf[500];
	    sprintf(buf,"#gamma #eta=%3.1f + PU %d",etaval[ieta],pu[ipu]);
	    TLatex lat;
	    lat.SetTextSize(0.07);
	    lat.DrawLatexNDC(0.25,0.965,buf);
	    sprintf(buf,"SR %d",iSR);
	    lat.DrawLatexNDC(0.5,0.87,buf);
	    if (iSR==4) lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
	    
	  }

	  myc[1]->Update();
	  lsave.str("");
	  lsave << plotDir << "/";
	  if (ipu==0) lsave << "ResoRaw";
	  else lsave << "Reso";
	  lsave << "_eta" << eta[ieta] << "_pu" << puOption << pSuffix;
	  myc[1]->Print((lsave.str()+".pdf").c_str());



	}//loop on pu
      }//loop on eta

      if (nGenEnAll==1) continue;

      TCanvas *mycR = new TCanvas("mycR","Sampling",1500,1000);
      TCanvas *mycC = new TCanvas("mycC","Constant",1500,1000);
      TCanvas *mycN = new TCanvas("mycN","Noise",1500,1000);
      TLegend *leg = new TLegend(0.75,0.8,0.94,0.94);
      leg->SetFillColor(10);

      mycR->Divide(2,1);
      mycC->Divide(2,1);
      mycN->Divide(2,1);
      TPad *mypad[3][neta];
      TPad *left[3];
      TPad *right[3];
      left[0] = (TPad*)mycR->cd(1);
      right[0] = (TPad*)mycR->cd(2);
      left[1] = (TPad*)mycC->cd(1);
      right[1] = (TPad*)mycC->cd(2);
      left[2] = (TPad*)mycN->cd(1);
      right[2] = (TPad*)mycN->cd(2);
      for (unsigned iC(0);iC<3;++iC){
	left[iC]->Divide(1,4);
	right[iC]->Divide(1,3);
	for (unsigned ieta=0; ieta<neta;++ieta){//loop on pt values
	  if (ieta<4) mypad[iC][ieta] = (TPad*)left[iC]->GetPad(ieta+1);
	  else mypad[iC][ieta] = (TPad*)right[iC]->GetPad((ieta-4)+1);
	}
      }

      TGraphErrors *grStoch[nPu-1][neta];
      TGraphErrors *grConst[nPu-1][neta];
      TGraphErrors *grNoise[nPu-1][neta];
      for (unsigned ieta(0); ieta<neta;++ieta){
	
	for (unsigned ipu(0); ipu<( (nPu>1)?(nPu-1):nPu ); ++ipu){//loop on pu
	  
	  grStoch[ipu][ieta] = new TGraphErrors(nSR,srval,sigmaStoch[ipu+1][ieta],srerr,sigmaStochErr[ipu+1][ieta]);
	  grConst[ipu][ieta] = new TGraphErrors(nSR,srval,sigmaConst[ipu+1][ieta],srerr,sigmaConstErr[ipu+1][ieta]);
	  grNoise[ipu][ieta] = new TGraphErrors(nSR,srval,sigmaNoise[ipu+1][ieta],srerr,sigmaNoiseErr[ipu+1][ieta]);

	  TGraphErrors *gr=0;
	  for (unsigned iP(0);iP<3;++iP){
	    gr = (iP==0) ? grStoch[ipu][ieta] : (iP==1) ? grConst[ipu][ieta] : grNoise[ipu][ieta];
	    if (!gr) continue;
	    if (iP==0) mypad[0][ieta]->cd();
	    else if (iP==1) mypad[1][ieta]->cd();
	    else mypad[2][ieta]->cd();
	    gPad->SetGridy(1);
	    gr->SetLineColor(ipu+1);
	    gr->SetMarkerColor(ipu+1);
	    gr->SetMarkerStyle(ipu+21);
	    gr->SetMinimum(0);
	    if (iP<2) gr->SetMaximum(iP==0?0.3:0.1);
	    if (iP==0) gr->SetTitle(";SR;sampling term (GeV^{#frac{1}{2}})");
	    else if (iP==1) gr->SetTitle(";SR;constant term");
	    else gr->SetTitle(";SR;noise term (GeV)");
	    gr->Draw( (ipu==0) ? "APEL" : "PEL");
	    std::ostringstream label;
	    if (ieta==0 && iP==0){
	      label.str("");
	      label << "PU " << pu[ipu+1];
	      leg->AddEntry(gr,label.str().c_str(),"P");
	    }
	    TLatex lat;
	    if (ipu==nPu-2){
	      label.str("");
	      label << "#eta=" << etaval[ieta];
	      lat.SetTextSize(0.06);
	      lat.DrawLatexNDC(0.2,0.85,label.str().c_str());
	      leg->Draw("same");
	    }
	    if (ieta==3 && ipu==0) lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
	  }
	}
      }
      
      mycR->Update();
      std::ostringstream lsave;
      lsave << plotDir << "/SamplingTerm_vsSR" << pSuffix;
      mycR->Print((lsave.str()+".pdf").c_str());
      mycR->Print((lsave.str()+".png").c_str());

      mycC->Update();
      lsave.str("");
      lsave << plotDir << "/ConstantTerm_vsSR" << pSuffix;
      mycC->Print((lsave.str()+".pdf").c_str());
      mycC->Print((lsave.str()+".png").c_str());
      
      mycN->Update();
      lsave.str("");
      lsave << plotDir << "/NoiseTerm_vsSR" << pSuffix;
      mycN->Print((lsave.str()+".pdf").c_str());
      mycN->Print((lsave.str()+".png").c_str());
      
      
    }//loop on scenarios
    
  }//loop on versions

  return 0;
  
}//main
