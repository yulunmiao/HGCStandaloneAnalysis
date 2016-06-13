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
void testEntry(TTree *truth, unsigned & truthIdx,double & trueE,unsigned & skip){
  truth->GetEntry(truthIdx);
  if (trueE<0.5) {
    truthIdx++;
    skip++;
    testEntry(truth,truthIdx,trueE,skip);
  }
}

int eMinimisationHgg(){

  TFile *output = TFile::Open("Eminimisation_pu0.root","RECREATE");

  TCanvas *myc = new TCanvas("myc","myc",1500,1000);
  myc->Divide(2,2);
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

  const unsigned nL = 30;

  TFile *fin;
  TTree *tree;
  double eLayer[nL];
  double absw[nL];
  double trueE;
  double wgtEtotal;
  unsigned evtIdx;

  double calib;
  double calibValeri;

  output->cd();
  TH2D *p_Matrix;
  TH1F *p_wgts;
  TH1F *p_wgtsValeri;
  TH2F *p_wgtEopt;
  TH2F *p_wgtEtotal;
  TH2F *p_wgtEtotalValeri;
  TH1F *p_wgtEtotalCaliboverEtrue;
  TH1F *p_wgtEtotalValeriCaliboverEtrue;
  TH1F *p_optEoverEtrue;

  TH2F *grReso[3];
  TH2F *grSigmaEff[3];
  for (unsigned iG(0); iG<3; ++iG){
    label.str("");
    label << "grReso_" << iG;
    grReso[iG] = new TH2F(label.str().c_str(),
			  "; E_{true} (GeV) ; #sigma/E",
			  100,0,1000,100,0,0.3);
    label.str("");
    label << "grSigmaEff_" << iG;
    grSigmaEff[iG] = new TH2F(label.str().c_str(),
			      "; E_{true} (GeV) ; #sigma_{eff}",
			      100,0,1000,100,0,50);
  }
  label.str("");
  label << "p_Matrix";//_" << E[iE];
  p_Matrix = new TH2D(label.str().c_str(),";i;j",nL,0,nL,nL,0,nL);

  TMatrixDSym matrix(nL);
  TVectorD v(nL);

  output->cd();
  label.str("");
  label << "p_wgtEopt";
  p_wgtEopt = new TH2F(label.str().c_str(),";E_{opt} (GeV);E_{truth} (GeV)",100,0,1000,100,0,1000);
  label.str("");
  label << "p_wgtEtotal";
  p_wgtEtotal = new TH2F(label.str().c_str(),";E_{absW} (Mips);E_{truth} (GeV)",1000,0,100000,100,0,1000);
  label.str("");
  label << "p_wgtEtotalCaliboverEtrue";
  p_wgtEtotalCaliboverEtrue = new TH1F(label.str().c_str(),";E_{absW} (GeV)",100,0.8,1.2);
  label.str("");
  label << "p_wgtEtotalValeriCaliboverEtrue";
  p_wgtEtotalValeriCaliboverEtrue = new TH1F(label.str().c_str(),";E_{absW}(Valeri) (GeV)",100,0.8,1.2);
  label.str("");
  label << "p_optEoverEtrue";
  p_optEoverEtrue = new TH1F(label.str().c_str(),";E_{opt}/E_{true}",50,0.8,1.2);
  
  label.str("");
  label << "p_wgtEtotalValeri";
  p_wgtEtotalValeri = new TH2F(label.str().c_str(),";E_{absW}(Valeri) (Mips);E_{truth} (GeV)",1000,0,100000,100,0,1000);
  
  label.str("");
  label << "p_wgts";
  p_wgts = new TH1F(label.str().c_str(),";layer;weights/ref(dEdx)",nL,0,nL);
  label.str("");
  label << "p_wgtsValeri";
  p_wgtsValeri = new TH1F(label.str().c_str(),";layer;weights/ref(ValeridEdx)",nL,0,nL);
  
  label.str("");
  label << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-14/version12/Hgg/200um/pu0_save.root";
  fin  = TFile::Open(label.str().c_str());

  std::string dir[2] = {"Gamma1E","Gamma2E"};
  std::string var[2] = {"Egamma1","Egamma2"};

  fin->cd(dir[0].c_str());  
  tree = (TTree*)gDirectory->Get("Ereso");
  tree->SetBranchAddress("wgtEtotal",&wgtEtotal);
  tree->SetBranchAddress("eventIndex",&evtIdx);
  for (unsigned iL(0); iL<nL; ++iL){//loop on layers
    label.str("");
    label << "energy_" << iL << "_SR4";
    tree->SetBranchAddress(label.str().c_str(),&eLayer[iL]);
  }

  //unsigned run[10] = {1156,1180,1186,1137,1177,1189,1176,1237,1275,1122};

  TFile *ftruth = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-14/version12/Hgg/200um/HggTruthInfo.root");
  //TFile *ftruth = TFile::Open("HggTruthInfo_run0.root");
  ftruth->cd();
  TTree *truth = (TTree*)gDirectory->Get("Etruth");
  truth->SetBranchAddress(var[0].c_str(),&trueE);
  for (unsigned iL(0); iL<nL; ++iL){//loop on layers
    label.str("");
    label << "absweight" << iL ;
    truth->SetBranchAddress(label.str().c_str(),&absw[iL]);
  }
  
  unsigned nEvts = tree->GetEntries();
  std::cout << "--- Processing: " << nEvts << " events" << std::endl;
  
  //unsigned previd = 0;
  //unsigned iR = 0;
  //unsigned totevt = 0;
  for (Long64_t ievt=0; ievt<nEvts;ievt++) {//loop on entries
    
    //if (iR==6) {
      //totevt += run[iR];
      //iR++;
    //}
    tree->GetEntry(ievt);
    if (ievt%500 == 0) {
      std::cout << "--- ... Processing event: " << ievt << std::endl;
    }

    /*if (evtIdx < previd)  
      {
	totevt += run[iR];
	iR++;
	std::cout << "--- ... Processing event: " << ievt 
		  << " totevt " << totevt
		  << " evtidx " << evtIdx
		  << " previd " << previd
		  << " run " << iR-1
		  << " nEvt " << run[iR-1]
		  << std::endl;
      }
    previd=evtIdx;
    */
    unsigned truthIdx = evtIdx;//evtIdx+totevt;
    truth->GetEntry(truthIdx);


    if (ievt%100==0) std::cout << ievt << " " << truth->GetEntries() << " " << evtIdx << std::endl;


    if (ievt==0) {
      std::cout << " absorber weights: " << std::endl;
      for (unsigned iL(0); iL<nL; ++iL){//loop on layers
	std::cout << " Layer " << iL << " w=" << absw[iL] << std::endl;
      }
    }
    //get corresponding event in truth tree
    //testEntry(truth,truthIdx,trueE,skip);

    
    double Etot = 0;
    double Evaleri = 0;
    for (unsigned iL(0); iL<nL; ++iL){//loop on layers
      v[iL] += eLayer[iL]*trueE/nEvts;
      double valVal = iL<(nL-1) ? (absw[iL]+absw[iL+1])/2. : absw[iL];
      Evaleri += eLayer[iL]*valVal;
      Etot += eLayer[iL]*absw[iL];
      for (unsigned jL(0); jL<nL; ++jL){//loop on layers
	matrix[iL][jL] += eLayer[iL]*eLayer[jL]*1./nEvts;
      }
    }

    p_wgtEtotal->Fill(Etot,trueE);
    p_wgtEtotalValeri->Fill(Evaleri,trueE);

  }//loop on entries
  
  //return 1;
  
  for (unsigned iL(0); iL<nL; ++iL){//loop on layers
    for (unsigned jL(0); jL<nL; ++jL){//loop on layers
      p_Matrix->Fill(iL,jL,matrix[iL][jL]);
    }
  }
  
  myc->cd(4);
  p_Matrix->SetStats(0);
  p_Matrix->Draw("colz");
  myc->Update();

  
  matrix.Invert();
  
  TVectorD weights(nL);
  weights = matrix*v;
  
  
  gStyle->SetOptStat(0);
  myc->cd(1);
  p_wgtEtotal->Draw("colz");
  p_wgtEtotal->Fit("pol1","","same");
  TF1 *fit = p_wgtEtotal->GetFunction("pol1");
  calib = fit->GetParameter(1);
  myc->cd(2);
  p_wgtEtotalValeri->Draw("colz");
  p_wgtEtotalValeri->Fit("pol1","","same");
  fit = p_wgtEtotalValeri->GetFunction("pol1");
  calibValeri = fit->GetParameter(1);
  std::cout << "Calib = " << calib << " Valeri's method " << calibValeri << std::endl;
  myc->Update();

  std::cout << " --output weights: " << std::endl;
  for (unsigned iL(0); iL<nL; ++iL){//loop on layers
    double ref = absw[iL]*calib;
    double refVal = iL<(nL-1) ? (absw[iL]+absw[iL+1])/2.*calib : absw[iL]*calib;
    //std::cout << iL << " & " << absw[iL] << " & " << ref << " & " << weights[iL] << "\\\\" << std::endl;
    p_wgts->Fill(iL,weights[iL]/ref);
    p_wgtsValeri->Fill(iL,weights[iL]/refVal);
    }
    
  //apply to second photon
  ftruth->cd();
  truth->SetBranchAddress(var[1].c_str(),&trueE);

  fin->cd(dir[1].c_str());  
  tree = (TTree*)gDirectory->Get("Ereso");
  for (unsigned iL(0); iL<nL; ++iL){//loop on layers
    label.str("");
    label << "energy_" << iL << "_SR4";
    tree->SetBranchAddress(label.str().c_str(),&eLayer[iL]);
  }
  tree->SetBranchAddress("wgtEtotal",&wgtEtotal);
  tree->SetBranchAddress("eventIndex",&evtIdx);


  nEvts = tree->GetEntries();

  //previd = 0;
  //iR = 0;
  //totevt = 0;
  for (Long64_t ievt=0; ievt<nEvts;ievt++) {//loop on entries

    //if (iR==6) {
    //totevt += run[iR];
    //iR++;
    //}

    tree->GetEntry(ievt);
    if (ievt%500 == 0) {
      std::cout << "--- ... Processing event: " << ievt << std::endl;
    }

    //get corresponding event in truth tree

    /*if (evtIdx < previd)  
      {
	totevt += run[iR];
	iR++;
	std::cout << "--- ... Processing event: " << ievt 
		  << " totevt " << totevt
		  << " evtidx " << evtIdx
		  << " previd " << previd
		  << " run " << iR-1
		  << " nEvt " << run[iR-1]
		  << std::endl;
      }
    previd=evtIdx;
    */
    unsigned truthIdx = evtIdx;//evtIdx+totevt;
    truth->GetEntry(truthIdx);

    //if (ievt%100==0) std::cout << ievt << " " << totevt << " " << evtIdx << " " << truthIdx << std::endl;

    double Etot = 0;
    double Evaleri = 0;
    double Eopt = 0;
    for (unsigned iL(0); iL<nL; ++iL){//loop on layers
      Eopt += eLayer[iL]*weights[iL];
      double val = iL<(nL-1) ? (absw[iL]+absw[iL+1])/2. : absw[iL];
      Evaleri += eLayer[iL]*val;
      Etot += eLayer[iL]*absw[iL];
    }
    p_wgtEtotalCaliboverEtrue->Fill(Etot*calib/trueE);
    
    p_wgtEtotalValeriCaliboverEtrue->Fill(Evaleri*calibValeri/trueE);
    p_optEoverEtrue->Fill(Eopt/trueE);
    p_wgtEopt->Fill(Eopt,trueE);
  }//loop on entries


  double sigmaeff_ref = effSigmaMacro(p_wgtEtotalCaliboverEtrue);
  double sigmaeff_ref_valeri = effSigmaMacro(p_wgtEtotalValeriCaliboverEtrue);
  double sigmaeff_opt = effSigmaMacro(p_optEoverEtrue);

  std::cout << " sigma_eff=" <<  sigmaeff_ref << " " << sigmaeff_ref_valeri << " " << sigmaeff_opt << std::endl;

  myc->cd(2);
  p_wgtEopt->Draw("colz");
  myc->Update();
  
  myc->cd(3);
  p_optEoverEtrue->SetLineColor(3);
  p_optEoverEtrue->Draw();
  p_wgtEtotalCaliboverEtrue->Draw("same");
  p_wgtEtotalValeriCaliboverEtrue->SetLineColor(2);
  p_wgtEtotalValeriCaliboverEtrue->Draw("same");
  TLegend *leg = new TLegend(0.5,0.5,0.89,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(p_wgtEtotalCaliboverEtrue,"Simple dEdx weights","P");
  leg->AddEntry(p_wgtEtotalValeriCaliboverEtrue,"Valeri's scheme","P");
  leg->AddEntry(p_optEoverEtrue,"Optimised weights","P");
  leg->Draw("same");

  myc->Update();
  myc->Print("HggPu0_Optimisation.pdf");

  mycW->Divide(2,1);
  mycW->cd(1);
  p_wgts->SetLineColor(1);
  p_wgts->SetMarkerColor(1);
  p_wgts->SetMarkerStyle(20);
  p_wgts->GetYaxis()->SetRangeUser(0,3.7);
  p_wgts->Draw("PL");

  mycW->cd(2);
  p_wgtsValeri->SetLineColor(1);
  p_wgtsValeri->SetMarkerColor(1);
  p_wgtsValeri->SetMarkerStyle(20);
  p_wgtsValeri->GetYaxis()->SetRangeUser(0,3.7);
  p_wgtsValeri->Draw("PL");
  
  mycW->Update();
  mycW->Print("HggPu0_OptimWeights.pdf");

   // Save the output
  output->Write();
  //output->Close();



  return 0;
}//main
