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
  
  TString genEn = "50";
  const double Emip = 0.0559;//in MeV


  TFile *inputFile = TFile::Open("DigiFiles/version_20/scenario_0/e-/e_"+genEn+"/DigiPFcal.root");
  if (!inputFile) {
    std::cout << " -- Error, input file cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  TTree *lTree = (TTree*)inputFile->Get("RecoTree");
  if (!lTree){
    std::cout << " -- Error, tree RecoTree cannot be opened either. Exiting..." << std::endl;
    return 1;
  }

  TCanvas *myc = new TCanvas("myc","myc",1000,500);

  TH1F *hitEnergy = new TH1F("hitEnergy",";E (MIPs)",25,0,50);
  TH1F *simhitEnergy = new TH1F("simhitEnergy",";E (MIPs)",25,0,50);

  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

  const unsigned nEvts = lTree->GetEntries();

  double maxEhit = 0;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (ievt%100==0) std::cout << "... Processing entry: " << ievt << std::endl;

    lTree->GetEntry(ievt);

    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];
      simhitEnergy->Fill(lHit.energy()/Emip);
    }//loop on hits


    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      hitEnergy->Fill(lHit.energy());
      if (lHit.energy()>maxEhit) maxEhit = lHit.energy();
    }//loop on rechits

  }//loop on entries

  std::cout << " -- max hit energy = " << maxEhit << std::endl;

  myc->cd();
  myc->Divide(2,1);
  myc->cd(1);
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  float yMax = simhitEnergy->GetMaximum();
  simhitEnergy->Draw();
  TLatex lat;
  lat.DrawLatex(0,yMax*2,"SIM e-, "+genEn+" GeV, 2.5#times2.5 mm^{2} cells");

  myc->cd(2);
  gPad->SetLogy(1);
  yMax = hitEnergy->GetMaximum();
  hitEnergy->Draw();
  lat.DrawLatex(0,yMax*2,"RECO e-, "+genEn+" GeV, 1#times1 cm^{2} cells");

  myc->Update();
  myc->Print("PLOTS/HitEnergy_1x1_"+genEn+"GeV.png");
  myc->Print("PLOTS/HitEnergy_1x1_"+genEn+"GeV.pdf");

  return 0;

}//main
