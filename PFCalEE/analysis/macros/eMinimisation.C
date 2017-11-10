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
#include "TPaveStats.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "effSigmaMacro.C"

int eMinimisation(){

  bool doPi = true;
  bool doHgg = false;

  bool addCst = false;

  TFile *output = 0;
  std::ostringstream outname;
  outname << "Eminimisation";
  if (doPi) outname << "_pi";
  else if (doHgg) outname << "_Hgg";
  if (addCst) outname << "_extraCstForFit";
  outname << ".root";
  output = TFile::Open(outname.str().c_str(),"RECREATE");

  TCanvas *myc = new TCanvas("myc","myc",1500,1000);
  myc->Divide(2,2);
  TCanvas *mycR = new TCanvas("mycR","mycR",1500,1000);
  mycR->Divide(2,2);
  TCanvas *mycW = new TCanvas("mycW","mycW",1500,1000);

  TLatex lat;
  std::ostringstream label;

  //for (unsigned iL(0); iL<nL; ++iL){//loop on layers
  // std::ostringstream label;
  //label << "energy_" << iL << "_SR5";
    //factory->AddVariable(label.str().c_str(),label.str().c_str(),"GeV",'F');
  // }//loop on layers
  //factory->AddVariable( "(energy_0_SR5)", "E0", "Mips" );
  //factory->AddVariable( "(energy_2_SR5+energy_4_SR5+energy_6_SR5+energy_8_SR5)", "E1even", "Mips" );
  //factory->AddVariable( "(energy_1_SR5+energy_3_SR5+energy_5_SR5+energy_7_SR5+energy_9_SR5)", "E1odd", "Mips" );
  //factory->AddVariable( "(energy_10_SR5+energy_12_SR5+energy_14_SR5+energy_16_SR5+energy_18_SR5)", "E2even", "Mips" );
  //factory->AddVariable( "(energy_11_SR5+energy_13_SR5+energy_15_SR5+energy_17_SR5+energy_19_SR5)", "E2odd", "Mips" );
  //factory->AddVariable( "(energy_20_SR5+energy_22_SR5+energy_24_SR5+energy_26_SR5)", "E3even", "Mips" );
  //factory->AddVariable( "(energy_21_SR5+energy_23_SR5+energy_25_SR5+energy_27_SR5)", "E3odd", "Mips" );

  const unsigned nL = doPi ? 52 : 28;
  const unsigned nE = (doHgg || doPi)? 1 : 7;
  unsigned E[7] = {5,10,20,30,50,70,100};
  double max[7] = {3000,5000,10000,20000,30000,40000,50000};
  if (doPi) {
    E[0] = 50;
    max[0] = 30000;
  }
  const double eta = 2.0;

  const unsigned nFit = addCst ? nL+1 : nL;

  TFile *fin[nE];
  TTree *tree[nE];
  double eLayer[nE][nFit];
  double absw[nL];
  double trueE[nE];
  double wgtEtotal[nE];

  double calib[nE];
  double calibValeri[nE];

  output->cd();
  TH2D *p_Matrix;
  TH1F *p_wgts[nE];
  TH1F *p_wgtsValeri[nE];
  TH1F *p_wgtEtotal[nE];
  TH1F *p_wgtEtotalValeri[nE];
  TH2F *p_wgtEtotal2D[nE];
  TH2F *p_wgtEtotalValeri2D[nE];
  TH1F *p_wgtEtotalCaliboverEtrue[nE];
  TH1F *p_wgtEtotalValeriCaliboverEtrue[nE];
  TH1F *p_optEoverEtrue[nE];
  TH2F *p_optE2D[nE];

  TH2F *p_EperLayer[nE];

  TGraphErrors *grLin[3];
  TGraphErrors *grReso[3];
  TGraphErrors *grSigmaEff[3];
  for (unsigned iG(0); iG<3; ++iG){
    grLin[iG] = new TGraphErrors();
    grLin[iG]->SetTitle("; E_{true} (GeV) ; E (GeV)");
    grReso[iG] = new TGraphErrors();
    grReso[iG]->SetTitle("; E_{true} (GeV) ; #sigma/E");
    grSigmaEff[iG] = new TGraphErrors();
    grSigmaEff[iG]->SetTitle("; E_{true} (GeV) ; #sigma_{eff}");
  }
  label.str("");
  label << "p_Matrix";//_" << E[iE];
  p_Matrix = new TH2D(label.str().c_str(),";i;j",nFit,0,nFit,nFit,0,nFit);

  //add constant term
  TMatrixDSym matrix(nFit);
  TVectorD v(nFit);

  for (unsigned iE(0); iE<nE; ++iE){//loop on energies
    output->cd();
    label.str("");
    label << "p_wgtEtotal_" << E[iE];
    p_wgtEtotal[iE] = new TH1F(label.str().c_str(),";E_{absW} (Mips)",1000,0,doHgg?100000:max[iE]);
    label.str("");
    label << "p_wgtEtotalCaliboverEtrue_" << E[iE];
    p_wgtEtotalCaliboverEtrue[iE] = new TH1F(label.str().c_str(),";E_{absW} (GeV)",100,0,2);//0.8,1.2);
    label.str("");
    label << "p_wgtEtotalValeriCaliboverEtrue_" << E[iE];
    p_wgtEtotalValeriCaliboverEtrue[iE] = new TH1F(label.str().c_str(),";E_{absW}(Valeri) (GeV)",100,0,2);//0.8,1.2);
    label.str("");
    label << "p_optEoverEtrue_" << E[iE];
    p_optEoverEtrue[iE] = new TH1F(label.str().c_str(),";E_{opt}/E_{true}",100,0,2);//0.8,1.2);

    label.str("");
    label << "p_wgtEtotalValeri_" << E[iE];
    p_wgtEtotalValeri[iE] = new TH1F(label.str().c_str(),";E_{absW}(Valeri) (Mips)",1000,0,doHgg?100000:max[iE]);

    label.str("");
    label << "p_wgtEtotal2D_" << E[iE];
    p_wgtEtotal2D[iE] = new TH2F(label.str().c_str(),";E_{absW} (Mips); E_{true} (GeV)",1000,0,doHgg?100000:max[iE],100,0,1000);
   label.str("");
    label << "p_wgtEtotalValeri2D_" << E[iE];
    p_wgtEtotalValeri2D[iE] = new TH2F(label.str().c_str(),";E_{absW}(Valeri) (Mips); E_{true} (GeV)",1000,0,doHgg?100000:max[iE],100,0,1000);
    label.str("");
    label << "p_optE2D_" << E[iE];
    p_optE2D[iE] = new TH2F(label.str().c_str(),";E_{opt} (GeV); E_{true} (GeV)",100,0,1000,100,0,1000);
    label.str("");
    label << "p_EperLayer_" << E[iE];
    p_EperLayer[iE] = new TH2F(label.str().c_str(),";layer;E_{layer} (Mips)",nL,0,nL,1000,0,5000);

    label.str("");
    label << "p_wgts_" << E[iE];
    p_wgts[iE] = new TH1F(label.str().c_str(),";layer;weights/ref(dEdx)",nFit,0,nFit);
    label.str("");
    label << "p_wgtsValeri_" << E[iE];
    p_wgtsValeri[iE] = new TH1F(label.str().c_str(),";layer;weights/ref(ValeridEdx)",nFit,0,nFit);

    label.str("");
    if (doHgg) label  << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalHexa/githexaV03-01-01/version30/Hgg/pu0_gamma1.root";
    else if (doPi) label  << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalHexa/githexaV03-01-01/version33/model2/pi-/eta20_et" << E[iE] << "_pu0_IC3.root";
    else label << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalHexa/gittestHexa/version100/model3/gamma/eta20_et" << E[iE] << "_pu0_IC3_Si2.root";
    fin[iE]  = TFile::Open(label.str().c_str());
    fin[iE]->cd("Energies");  

    tree[iE] = (TTree*)gDirectory->Get("Ereso");
    for (unsigned iL(0); iL<nL; ++iL){//loop on layers
      label.str("");
      label << "energy_" << iL;
      if (doHgg) label << "_SR4";
      //else if (doPi) label << "_SR5";
      else label << "_SR5";
      tree[iE]->SetBranchAddress(label.str().c_str(),&eLayer[iE][iL]);
      label.str("");
      label << "absweight_" << iL ;
      tree[iE]->SetBranchAddress(label.str().c_str(),&absw[iL]);
    }
    eLayer[iE][nL] = 1;
    tree[iE]->SetBranchAddress("trueE",&trueE[iE]);
    tree[iE]->SetBranchAddress("wgtEtotal",&wgtEtotal[iE]);
 
    const unsigned nEvts = (doHgg?1:0.5)*tree[iE]->GetEntries();
    std::cout << "--- Processing: " << nEvts << " events" << std::endl;


    for (Long64_t ievt=0; ievt<nEvts;ievt++) {//loop on entries
      
      tree[iE]->GetEntry(ievt);
      if (ievt%100 == 0) {
	std::cout << "--- ... Processing event: " << ievt << std::endl;
      }

      double Etot = 0;
      double Evaleri = 0;
      for (unsigned iL(0); iL<nFit; ++iL){//loop on layers
	v[iL] += eLayer[iE][iL]*trueE[iE]/nEvts;
	if (iL<nL){
	  double valVal = iL<(nL-1) ? (absw[iL]+absw[iL+1])/2. : absw[iL];
	  Evaleri += eLayer[iE][iL]*valVal;
	  Etot += eLayer[iE][iL]*absw[iL];
	  p_EperLayer[iE]->Fill(iL,eLayer[iE][iL]);
	}
	for (unsigned jL(0); jL<nFit; ++jL){//loop on layers
	  matrix[iL][jL] += eLayer[iE][iL]*eLayer[iE][jL]*1./nEvts;
	}
      }
      if (doHgg){
	//p_wgtEtotal[iE]->Fill(Etot);
	p_wgtEtotal2D[iE]->Fill(Etot,trueE[iE]);
	p_wgtEtotalValeri2D[iE]->Fill(Evaleri,trueE[iE]);
      }
      if (!doHgg) {
	p_wgtEtotal[iE]->Fill(doPi?Etot : wgtEtotal[iE]);
	p_wgtEtotalValeri[iE]->Fill(Evaleri);
      }

    }//loop on entries

  }//loop on energies

  for (unsigned iL(0); iL<nFit; ++iL){//loop on layers
    for (unsigned jL(0); jL<nFit; ++jL){//loop on layers
      p_Matrix->Fill(iL,jL,matrix[iL][jL]);
    }
  }
  
  myc->cd(4);
  gPad->SetLogz(1);
  p_Matrix->SetStats(0);
  p_Matrix->Draw("colz");
  myc->Update();


  matrix.Invert();
  
  TVectorD weights(nFit);
  weights = matrix*v;

  
  for (unsigned iE(0); iE<nE; ++iE){//loop on energies
    
    gStyle->SetOptStat(0);
    if (!doHgg){
      myc->cd(1);
      p_wgtEtotal[iE]->Draw();
      p_wgtEtotal[iE]->Fit("gaus","","same");
      TF1 *fit = p_wgtEtotal[iE]->GetFunction("gaus");
      if (fit) calib[iE] = fit->GetParameter(1);
      myc->cd(2);
      p_wgtEtotalValeri[iE]->Draw();
      p_wgtEtotalValeri[iE]->Fit("gaus","","same");
      fit = p_wgtEtotalValeri[iE]->GetFunction("gaus");
      if (fit) calibValeri[iE] = fit->GetParameter(1);
      std::cout << " -- E = " << E[iE] << ", Gaussian mean = " << calib[iE] << " Valeri's method " << calibValeri[iE] << std::endl;
    }
    else {
      myc->cd(1);
      p_wgtEtotal2D[iE]->Draw("colz");
      p_wgtEtotal2D[iE]->Fit("pol1","","same");
      TF1 *fit = p_wgtEtotal2D[iE]->GetFunction("pol1");
      calib[iE] = fit->GetParameter(1);
      myc->cd(2);
      p_wgtEtotalValeri2D[iE]->Draw("colz");
      p_wgtEtotalValeri2D[iE]->Fit("pol1","","same");
      fit = p_wgtEtotalValeri2D[iE]->GetFunction("pol1");
      calibValeri[iE] = fit->GetParameter(1);
      std::cout << " -- E = " << E[iE] << ", Gaussian mean = " << calib[iE] << " Valeri's method " << calibValeri[iE] << std::endl;
      myc->Update();
    }

    double Etrue = doHgg?1:E[iE]*cosh(eta);
    std::cout << " --output weights: " << std::endl;
    for (unsigned iL(0); iL<nFit; ++iL){//loop on layers
      double ref = iL<nL?(doHgg?absw[iL]*calib[iE]:absw[iL]*Etrue/calib[iE]) : 1;
      double refVal = iL<(nL-2) ? (doHgg?(absw[iL]+absw[iL+1])/2.*calib[iE]:(absw[iL]+absw[iL+1])/2.*Etrue/calib[iE]) : iL<(nL-1) ? (doHgg?absw[iL]*calib[iE]:absw[iL]*Etrue/calib[iE]) : 1;
      //std::cout << iL << " & " << absw[iL] << " & " << ref << " & " << weights[iL] << "\\\\" << std::endl;
      p_wgts[iE]->Fill(iL,weights[iL]/ref);
      p_wgtsValeri[iE]->Fill(iL,weights[iL]/refVal);
    }
    
    if (doHgg) {
      label.str("");
      label  << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalHexa/githexaV03-01-01/version30/Hgg/pu0_gamma2.root";
      fin[iE]  = TFile::Open(label.str().c_str());
      fin[iE]->cd("Energies");  
      
      tree[iE] = (TTree*)gDirectory->Get("Ereso");
      for (unsigned iL(0); iL<nL; ++iL){//loop on layers
	label.str("");
	label << "energy_" << iL;
	if (doHgg) label << "_SR4";
	else label  << "_SR5";
	tree[iE]->SetBranchAddress(label.str().c_str(),&eLayer[iE][iL]);
	label.str("");
	label << "absweight_" << iL ;
	tree[iE]->SetBranchAddress(label.str().c_str(),&absw[iL]);
      }
      eLayer[iE][nL] = 1;
      tree[iE]->SetBranchAddress("trueE",&trueE[iE]);
      tree[iE]->SetBranchAddress("wgtEtotal",&wgtEtotal[iE]);
    }
    const unsigned nEvts = tree[iE]->GetEntries();

    for (Long64_t ievt=(doHgg?0:nEvts/2); ievt<nEvts;ievt++) {//loop on entries
      
      tree[iE]->GetEntry(ievt);
      if (ievt%100 == 0) {
	std::cout << "--- ... Processing event: " << ievt << std::endl;
      }
      double Evaleri = 0;
      double Eopt = 0;
      double Etot = 0;
      for (unsigned iL(0); iL<nFit; ++iL){//loop on layers
	Eopt += eLayer[iE][iL]*weights[iL];
	if (iL<nL){
	  double val = iL<(nL-1) ? (absw[iL]+absw[iL+1])/2. : absw[iL];
	  Evaleri += eLayer[iE][iL]*val;
	  Etot += eLayer[iE][iL]*absw[iL];
	  p_EperLayer[iE]->Fill(iL,eLayer[iE][iL]);
	}
      }
      p_wgtEtotalCaliboverEtrue[iE]->Fill(doHgg?Etot*calib[iE]/trueE[iE]:(doPi?Etot:wgtEtotal[iE])/calib[iE]);

      p_wgtEtotalValeriCaliboverEtrue[iE]->Fill(doHgg?Evaleri*calibValeri[iE]/trueE[iE]:Evaleri/calibValeri[iE]);
      p_optEoverEtrue[iE]->Fill(Eopt/trueE[iE]);
      p_optE2D[iE]->Fill(Eopt,trueE[iE]);
    }//loop on entries

    if (doHgg){
      myc->cd(3);
      p_optE2D[iE]->Draw("colz");
      p_optE2D[iE]->Fit("pol1","","same");
    }

    mycR->cd(1);
    p_optEoverEtrue[iE]->Fit("gaus");
    TF1 *fit = p_optEoverEtrue[iE]->GetFunction("gaus");
    if (fit){
      fit->SetLineColor(3);
      grReso[2]->SetPoint(iE,Etrue,fit->GetParameter(2)/fit->GetParameter(1));
      grReso[2]->SetPointError(iE,0,fit->GetParError(2)/fit->GetParameter(1));
      grLin[2]->SetPoint(iE,Etrue,fit->GetParameter(1));
      grLin[2]->SetPointError(iE,0,fit->GetParError(1));
    }
    p_wgtEtotalCaliboverEtrue[iE]->Fit("gaus");
    fit = p_wgtEtotalCaliboverEtrue[iE]->GetFunction("gaus");
    if (fit){
      fit->SetLineColor(1);
      grReso[0]->SetPoint(iE,Etrue,fit->GetParameter(2)/fit->GetParameter(1));
      grReso[0]->SetPointError(iE,0,fit->GetParError(2)/fit->GetParameter(1));
      grLin[0]->SetPoint(iE,Etrue,fit->GetParameter(1));
      grLin[0]->SetPointError(iE,0,fit->GetParError(1));
    }
    p_wgtEtotalValeriCaliboverEtrue[iE]->Fit("gaus");
    fit = p_wgtEtotalValeriCaliboverEtrue[iE]->GetFunction("gaus");
    if (fit){
      fit->SetLineColor(2);
      grReso[1]->SetPoint(iE,Etrue,fit->GetParameter(2)/fit->GetParameter(1));
      grReso[1]->SetPointError(iE,0,fit->GetParError(2)/fit->GetParameter(1));
      grLin[1]->SetPoint(iE,Etrue,fit->GetParameter(1));
      grLin[1]->SetPointError(iE,0,fit->GetParError(1));
    }

    double sigmaeff_ref = effSigmaMacro(p_wgtEtotalCaliboverEtrue[iE]);
    double sigmaeff_ref_valeri = effSigmaMacro(p_wgtEtotalValeriCaliboverEtrue[iE]);
    double sigmaeff_opt = effSigmaMacro(p_optEoverEtrue[iE]);

    std::cout << "E " << Etrue << " sigma_eff=" <<  sigmaeff_ref << " " << sigmaeff_ref_valeri << " " << sigmaeff_opt << std::endl;

    //if (doHgg){
    mycR->cd(4);
    p_optEoverEtrue[iE]->SetLineColor(3);
    p_optEoverEtrue[iE]->Draw();
    p_wgtEtotalCaliboverEtrue[iE]->SetLineColor(1);
      p_wgtEtotalCaliboverEtrue[iE]->Draw("same");
      p_wgtEtotalValeriCaliboverEtrue[iE]->SetLineColor(2);
      p_wgtEtotalValeriCaliboverEtrue[iE]->Draw("same");

      std::cout << "Opt: " <<p_optEoverEtrue[iE]->GetEntries()
		<< " " << p_optEoverEtrue[iE]->GetBinContent(0)
		<< " " << p_optEoverEtrue[iE]->GetBinContent(p_optEoverEtrue[iE]->GetNbinsX()+1)
		<< std::endl
		<< " absW " <<  p_wgtEtotalCaliboverEtrue[iE]->GetEntries()
		<< " " << p_wgtEtotalCaliboverEtrue[iE]->GetBinContent(0)
		<< " " << p_wgtEtotalCaliboverEtrue[iE]->GetBinContent(p_wgtEtotalCaliboverEtrue[iE]->GetNbinsX()+1)
		<< std::endl
		<< " Valeri " << p_wgtEtotalValeriCaliboverEtrue[iE]->GetEntries()
		<< " " << p_wgtEtotalValeriCaliboverEtrue[iE]->GetBinContent(0)
		<< " " << p_wgtEtotalValeriCaliboverEtrue[iE]->GetBinContent(p_wgtEtotalValeriCaliboverEtrue[iE]->GetNbinsX()+1)
		<< std::endl;

      //}


    grSigmaEff[0]->SetPoint(iE,Etrue,sigmaeff_ref);
    grSigmaEff[0]->SetPointError(iE,0,p_wgtEtotalCaliboverEtrue[iE]->GetRMS()/sqrt(2.*p_wgtEtotalCaliboverEtrue[iE]->GetEntries()));

    grSigmaEff[1]->SetPoint(iE,Etrue,sigmaeff_ref_valeri);
    grSigmaEff[1]->SetPointError(iE,0,p_wgtEtotalValeriCaliboverEtrue[iE]->GetRMS()/sqrt(2.*p_wgtEtotalValeriCaliboverEtrue[iE]->GetEntries()));

    grSigmaEff[2]->SetPoint(iE,Etrue,sigmaeff_opt);
    grSigmaEff[2]->SetPointError(iE,0,p_optEoverEtrue[iE]->GetRMS()/sqrt(2.*p_optEoverEtrue[iE]->GetEntries()));

  }//loop on energies
  if (!doHgg){
    myc->cd(1);
    p_wgtEtotal[nE-1]->Draw();
    for (unsigned iE(0); iE<nE; ++iE){//loop on energies
      p_wgtEtotal[iE]->Draw("same");
    }
    myc->cd(2);
    p_wgtEtotalValeri[nE-1]->Draw();
    for (unsigned iE(0); iE<nE; ++iE){//loop on energies
      p_wgtEtotalValeri[iE]->Draw("same");
    }
    myc->cd(3);
    gPad->SetLogz(1);
    p_EperLayer[0]->Draw("colz");
    p_EperLayer[0]->ProfileX();//"_pfx",1,-1,"d");
    TProfile *tmp = (TProfile*)gDirectory->Get("p_EperLayer_50_pfx");
    tmp->SetMarkerStyle(22);
    p_EperLayer[0]->GetYaxis()->SetRangeUser(0,350);
    tmp->Draw("PEsame");
    myc->Update();
  }

  myc->Update();
  std::ostringstream lsave;
  lsave << "Optimisation";
  if (doHgg) lsave << "_Hgg";
  else if (doPi) lsave << "_pi";
  if (addCst) lsave << "_extraCstForFit";
  lsave << ".pdf";

  myc->Print(lsave.str().c_str());

  TLegend *leg = new TLegend(0.6,0.6,0.99,0.99);
  leg->SetFillColor(0);
  leg->AddEntry(grReso[0],"Simple dEdx weights","P");
  leg->AddEntry(grReso[1],"Valeri's scheme","P");
  leg->AddEntry(grReso[2],"Optimised weights","P");

  mycR->cd(3);
  for (unsigned iG(0); iG<3; ++iG){
    grSigmaEff[iG]->SetLineColor(iG+1);
    grSigmaEff[iG]->SetMarkerColor(iG+1);
    grSigmaEff[iG]->SetMarkerStyle(20+iG);
    if (!doPi) grSigmaEff[iG]->GetYaxis()->SetRangeUser(0.018,0.022);
    else grSigmaEff[iG]->GetYaxis()->SetRangeUser(0.,0.1);
    grSigmaEff[iG]->Draw(iG==0?"APE":"PE");
  }
  leg->Draw("same");

  mycR->cd(1);
  gPad->SetGridy(1);
  //fitFunc->SetParameters();
  for (unsigned iG(0); iG<3; ++iG){
    grLin[iG]->SetLineColor(iG+1);
    grLin[iG]->SetMarkerColor(iG+1);
    grLin[iG]->SetMarkerStyle(20+iG);
    grLin[iG]->GetYaxis()->SetRangeUser(0.95,1.05);
    grLin[iG]->Draw(iG==0?"APE":"PE");
    if (!doHgg && !doPi) grLin[iG]->Fit("pol1","BIME0","same");
  }
  leg->Draw("same");

  mycR->cd(2);
  TF1 *fitFunc[3];
  //fitFunc->SetParameters();
  for (unsigned iG(0); iG<3; ++iG){
    fitFunc[iG] = new TF1("reso","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,1000);
    grReso[iG]->SetLineColor(iG+1);
    grReso[iG]->SetMarkerColor(iG+1);
    grReso[iG]->SetMarkerStyle(20+iG);
    if (!doPi) grReso[iG]->GetYaxis()->SetRangeUser(0.018,0.022);
    else grReso[iG]->GetYaxis()->SetRangeUser(0.,0.1);
    grReso[iG]->Draw(iG==0?"APE":"PE");
    if (!doHgg && !doPi){
      grReso[iG]->Fit(fitFunc[iG],"BIME0");
      fitFunc[iG]->SetLineColor(1+iG);
      fitFunc[iG]->Draw("same");
    }
  }
  leg->Draw("same");

  mycR->cd(4);
  gPad->SetGridx(1);
  leg->Draw("same");

  mycR->Update();
  lsave.str("");
  lsave << "OptimReso";
  if (doHgg) lsave << "_Hgg";
  else if (doPi) lsave << "_pi";
  if (addCst) lsave << "_extraCstForFit";
  lsave << ".pdf";
  mycR->Print(lsave.str().c_str());

  mycW->Divide(2,1);
  mycW->cd(1);
  gPad->SetGridy(1);
  for (unsigned iE(0); iE<nE; ++iE){//loop on energies
    p_wgts[iE]->SetLineColor(iE+1);
    p_wgts[iE]->SetMarkerColor(iE+1);
    p_wgts[iE]->SetMarkerStyle(20+iE);
    p_wgts[iE]->GetYaxis()->SetRangeUser(-1,doPi ? 3 : 10);
    p_wgts[iE]->Draw(iE==0?"PL":"PLsame");
  }
  mycW->cd(2);
  gPad->SetGridy(1);
  for (unsigned iE(0); iE<nE; ++iE){//loop on energies
    p_wgtsValeri[iE]->SetLineColor(iE+1);
    p_wgtsValeri[iE]->SetMarkerColor(iE+1);
    p_wgtsValeri[iE]->SetMarkerStyle(20+iE);
    p_wgtsValeri[iE]->GetYaxis()->SetRangeUser(-1,doPi ? 3 : 10);
    p_wgtsValeri[iE]->Draw(iE==0?"PL":"PLsame");
  }
  mycW->Update();

  lsave.str("");
  lsave << "OptimWeights";
  if (doHgg) lsave << "_Hgg";
  else if (doPi) lsave << "_pi";
  if (addCst) lsave << "_extraCstForFit";
  lsave << ".pdf";

  mycW->Print(lsave.str().c_str());

   // Save the output
  output->Write();
  //output->Close();



  return 0;
}//main
