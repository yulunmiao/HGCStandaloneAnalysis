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
#include "TChain.h"
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


#include "utilities.h"

int main(int argc, char** argv){//main

  const unsigned nPar = static_cast<unsigned>(argc);
  if (nPar < 5) {
    std::cout << " Usage: "
              << argv[0] << " <nEvts to process (0=all)>"
	      << " <eta value: 1.600,2.000,2.850>"
	      << " <output dir>"
	      << " <BON ? 1>" << std::endl;
    return 1;
  }
  
  unsigned nEvts = atoi(argv[1]);
  std::string etaval = argv[2];
  std::string plotBase = argv[3];
  bool bfieldON = atoi(argv[4]);

  std::string inFilePath = "root://eoscms//eos/cms/store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR/gitV08-02-00/mu-/HGcal__version63_model2_";
  if (bfieldON) inFilePath += "BON";
  else inFilePath += "BOFF";
  //if (etaval=="1.600") inFilePath += "_et20_eta";
  //else inFilePath += "_et13_eta";
  inFilePath += "_et0_eta";

  double expMPV = etaval=="1.600" ? 0.09 : etaval=="2.000"? 0.055 : 0.03;

  const unsigned nSi = etaval=="1.600" ? 3 : etaval=="2.000"? 2 : 1;
  const unsigned firstLayerScint = 52;

  std::string outFilePath = plotBase+"/mipcalib_E0GeV_eta"+etaval;
  if (bfieldON) outFilePath += "BON";
  outFilePath += ".root";

  TChain *lTree = new TChain("HGCSSTree");
  TFile * mipFile = 0;
  for (unsigned i(0);i<10;++i){
    std::ostringstream lstrrec;
    lstrrec << inFilePath << etaval << "_run" << i << ".root";
    if (!testInputFile(lstrrec.str(),mipFile)) continue;
    lTree->AddFile(lstrrec.str().c_str());
  }
  
  
  if (!lTree){
    std::cout << " -- Error, tree MipTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  
  std::cout << " Trees added." << std::endl;
  
  
  TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();
  
  const unsigned nLayers = 69;//33;
  
  TTree *outtree = new TTree("MipFit","Tree for Mip calibration");
  std::vector<double> EperLayer;
  EperLayer.resize(nLayers,0);
  std::vector<double> nHitsperLayer;
  nHitsperLayer.resize(nLayers,0);
  for (unsigned iL(0); iL<nLayers; ++iL){
    EperLayer[iL] = 0;
    std::ostringstream label;
    label << "E_layer" << iL;
    outtree->Branch(label.str().c_str(),&EperLayer[iL]);
    nHitsperLayer[iL] = 0;
    label.str("");
    label << "nHits_layer" << iL;
    outtree->Branch(label.str().c_str(),&nHitsperLayer[iL]);    
  }

  TH2F *p_nHits = new TH2F("nHits","; layer; Number of hits; Events",nLayers,0,nLayers,20,0,20);
  TH1F *p_hitEnergy_si = new TH1F("hitEnergy_si",";E (MeV);SimHits",1000,0,0.2);
  TH1F *p_hitEnergySel_si = new TH1F("hitEnergySel_si",";E (MeV);SimHits",1000,0,0.2);
  TH1F *p_hitEnergySel_si_EE = new TH1F("hitEnergySel_si_EE",";E (MeV);SimHits",1000,0,0.2);
  TH1F *p_hitEnergy_scint = new TH1F("hitEnergy_scint",";E (MeV);SimHits",1000,0,10);
  TH1F *p_hitEnergySel_scint = new TH1F("hitEnergySel_scint",";E (MeV);SimHits",1000,0,10);

  std::vector<HGCSSSimHit> * simhitvec = 0;
  lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);

  if (nEvts==0 || nEvts > lTree->GetEntries()) nEvts = lTree->GetEntries();

  std::cout << " -- Processing " << nEvts << " entries." << std::endl;
  
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    
    lTree->GetEntry(ievt);
    
    if (ievt%1000==0) std::cout << "- Event " << ievt << std::endl;
    
    for (unsigned iL(0); iL<nLayers; ++iL){
      nHitsperLayer[iL] = 0;
      EperLayer[iL] = 0;
    }
    
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];
      double energy = lHit.energy();
      if (energy<=0) continue;
      unsigned layer = lHit.layer();
      if (lHit.silayer() >= nSi) continue;
      if (layer < firstLayerScint) p_hitEnergy_si->Fill(energy);
      else p_hitEnergy_scint->Fill(energy);
      nHitsperLayer[layer]++;
    }//loop on hits

    for (unsigned iL(0); iL<nLayers; ++iL){
      p_nHits->Fill(iL,nHitsperLayer[iL]);
    }
    
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];
      unsigned layer = lHit.layer();
      if (layer < firstLayerScint && nHitsperLayer[layer]!=nSi) continue;
      if (layer > firstLayerScint && nHitsperLayer[layer]!=1) continue;
      if (lHit.silayer() >= nSi) continue;
      EperLayer[layer] += lHit.energy();
    }

    for (unsigned iL(0); iL<nLayers; ++iL){
      if (EperLayer[iL] <0.001) continue;
      if (iL < firstLayerScint) p_hitEnergySel_si->Fill(EperLayer[iL]);
      if (iL < 28) p_hitEnergySel_si_EE->Fill(EperLayer[iL]);
      else  p_hitEnergySel_scint->Fill(EperLayer[iL]);
    }

    outtree->Fill();
  }//loop on entries
  

  TCanvas *myc = new TCanvas("myc","myc",1);
  
  myc->cd();
  gPad->SetLogz(1);
  gStyle->SetOptStat(1111110);
  p_nHits->Draw("colz");

  myc->Update();
  myc->Print((plotBase+"/mipHits.pdf").c_str());
  myc->Print((plotBase+"/mipHits.C").c_str());


  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergy_si->Draw();
  p_hitEnergy_si->Fit("landau","+","",expMPV*0.8,expMPV*2);

  myc->Update();
  myc->Print((plotBase+"/mipDepositAll_si.pdf").c_str());
  myc->Print((plotBase+"/mipDepositAll_si.C").c_str());

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel_si->Draw();
  p_hitEnergySel_si->Fit("landau","L+","",expMPV*0.8,expMPV*2);

  myc->Update();
  myc->Print((plotBase+"/mipDepositSel_si.pdf").c_str());
  myc->Print((plotBase+"/mipDepositSel_si.C").c_str());

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel_si_EE->Draw();
  p_hitEnergySel_si_EE->Fit("landau","L+","",expMPV*0.8,expMPV*2);

  myc->Update();
  myc->Print((plotBase+"/mipDepositSel_si_EE.pdf").c_str());
  myc->Print((plotBase+"/mipDepositSel_si_EE.C").c_str());


  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergy_scint->Draw();
  p_hitEnergy_scint->Fit("landau","+");
  
  myc->Update();
  myc->Print((plotBase+"/mipDepositAll_scint.pdf").c_str());
  myc->Print((plotBase+"/mipDepositAll_scint.C").c_str());

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel_scint->Draw();
  p_hitEnergySel_scint->Fit("landau","L+");
  
  myc->Update();
  myc->Print((plotBase+"/mipDepositSel_scint.pdf").c_str());
  myc->Print((plotBase+"/mipDepositSel_scint.C").c_str());

  outputFile->cd();
  outtree->Write();
  outputFile->Write();
  outputFile->Close();

  return 0;

}//main
