#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"


int main(int argc, char** argv){//main


  TFile *inputFile = TFile::Open("root://eoscms//eos/cms/store/user/amagnan/HGCalEEGeant4/mu-/HGcal_version3_e50.root");
  if (!inputFile) {
    std::cout << " -- Error, input file cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
  if (!lTree){
    std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  
  TCanvas *myc = new TCanvas("myc","myc",500,500);
  
  TH2F *p_nHits = new TH2F("nHits","; layer; Number of hits; Events",30,0,30,20,0,20);
  TH1F *p_hitEnergy = new TH1F("hitEnergy",";E (MeV);SimHits",250,0,1);
  TH1F *p_hitEnergySel = new TH1F("hitEnergySel",";E (MeV);SimHits",100,0.01,0.5);

  std::vector<HGCSSSimHit> * simhitvec = 0;
  float volNb = 0;
  lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lTree->SetBranchAddress("volNb",&volNb);

  const unsigned nEvts = lTree->GetEntries();

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

    lTree->GetEntry(ievt);
    
    if (ievt/30%100==0) std::cout << std::endl << "entry " << ievt << " volNb = " << volNb << " : ";

    unsigned nHits = 0;
    double energySel = 0;
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];
      double energy = lHit.energy();
      if (energy>0) {
	p_hitEnergy->Fill(energy);
	nHits++;
	if (nHits==1) energySel = energy;
      }
    }//loop on hits
    p_nHits->Fill(volNb,nHits);
    if (nHits==1) p_hitEnergySel->Fill(energySel);

  }//loop on entries

  
  myc->cd();
  gPad->SetLogz(1);
  gStyle->SetOptStat(1111110);
  p_nHits->Draw("colz");

  myc->Update();
  myc->Print("PLOTS/version_3/mu-/mipHits.png");
  myc->Print("PLOTS/version_3/mu-/mipHits.pdf");
  myc->Print("PLOTS/version_3/mu-/mipHits.C");


  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergy->Draw();
  p_hitEnergy->Fit("landau","R+","",0.035,1);
  
  myc->Update();
  myc->Print("PLOTS/version_3/mu-/mipDepositAll.png");
  myc->Print("PLOTS/version_3/mu-/mipDepositAll.pdf");
  myc->Print("PLOTS/version_3/mu-/mipDepositAll.C");

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel->Draw();
  p_hitEnergySel->Fit("landau","LR+","",0.02,0.5);
  
  myc->Update();
  myc->Print("PLOTS/version_3/mu-/mipDepositSel.png");
  myc->Print("PLOTS/version_3/mu-/mipDepositSel.pdf");
  myc->Print("PLOTS/version_3/mu-/mipDepositSel.C");



  return 0;

}//main
