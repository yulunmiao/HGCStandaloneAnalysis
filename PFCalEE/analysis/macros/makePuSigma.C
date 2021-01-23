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

bool retrievePuSigma(const unsigned pu,
		     const TString version, 
		     const unsigned nLayers,
		     TTree *atree, TTree *atreePu, 
		     const std::vector<double> & calib,
		     const double & etaval,
		     std::vector<double> & noisePu){
  
  unsigned nSR=6;
  if (!atree || !atreePu) {
    std::cout << " -- Info, both trees were not found:" << atree << " " << atreePu << std::endl;
    return false;
  }
  else std::cout << " Trees found." << std::endl;

  std::ostringstream foutname;
  foutname << "PLOTS/PuSubtraction_v" << version << ".root";
  TFile *fout = TFile::Open(foutname.str().c_str(),"RECREATE");
  
  std::ostringstream lsave;
  lsave << "eta"<< etaval  << "_pu" << pu;
  fout->mkdir(lsave.str().c_str());
  fout->cd(lsave.str().c_str());
  gStyle->SetOptStat("eMRuo");
  TH1F *p_sigma[nSR];
  TH2F *p_EvspuE[nSR];
  std::ostringstream label;
  if (pu>0){
    for (unsigned iSR(0); iSR<nSR;++iSR){
      label.str("");
      label << "p_sigma_" << iSR << "_" << lsave.str();
      //if (iSR<(nSR-1) ) 
      p_sigma[iSR] = new TH1F(label.str().c_str(),";PuE-E (GeV)",500,-20,20);
      //else p_sigma[iSR] = new TH1F(label.str().c_str(),";PuE-E (GeV)",5000,0,10000);
      label.str("");
      label << "p_EvspuE_" << iSR << "_" << lsave.str();
      p_EvspuE[iSR] = new TH2F(label.str().c_str(),";PuE-E (GeV);PuE-PuEsubtr (GeV)",
			       500,-20,20,500,-20,20);
      //p_sigma[iSR]->StatOverflows();
    }
  }

  unsigned evtIndex = 0;
  unsigned evtIndexPu = 0;
  atree->SetBranchAddress("eventIndex",&evtIndex);
  atreePu->SetBranchAddress("eventIndex",&evtIndexPu);

  if (atree->GetBranchStatus("eventIndex") != 1 || 
      atreePu->GetBranchStatus("eventIndex") != 1) {
    std::cout << " -- Error! Branch eventIndex not found: " << atree->GetBranchStatus("eventIndex") << " " << atreePu->GetBranchStatus("eventIndex") << std::endl;
    return false;
  }
  double totalE = 0;
  double totalEpu = 0;
  atree->SetBranchAddress("wgtEtotal",&totalE);
  atreePu->SetBranchAddress("wgtEtotal",&totalEpu);

  std::vector<std::vector<double> > energySR[2];
  std::vector<std::vector<double> > subtractedenergySR;
  std::vector<double> absweight;

  std::vector<double> emptyvec;
  emptyvec.resize(nSR,0);
  energySR[0].resize(nLayers,emptyvec);
  energySR[1].resize(nLayers,emptyvec);
  subtractedenergySR.resize(nLayers,emptyvec);
  absweight.resize(nLayers,1);

  for (unsigned iL(0); iL<nLayers;++iL){
    absweight[iL] = 10.0166;
    //label.str("");
    //label << "absweight_" << iL;
    //atree->SetBranchAddress(label.str().c_str(),&absweight[iL]);
    //atreePu->SetBranchAddress(label.str().c_str(),&absweight[iL]);

    for (unsigned iSR(0);iSR<nSR;++iSR){
      label.str("");
      label << "energy_" << iL << "_SR" << iSR;
      atree->SetBranchAddress(label.str().c_str(),&energySR[0][iL][iSR]);
      atreePu->SetBranchAddress(label.str().c_str(),&energySR[1][iL][iSR]);
      label.str("");
      label << "subtractedenergy_" << iL << "_SR" << iSR;
      atreePu->SetBranchAddress(label.str().c_str(),&subtractedenergySR[iL][iSR]);
    }
  }
  absweight[0] = 20.3628;
  absweight[nLayers-1] = 13.0629;
  int nEvtsPu = atreePu->GetEntries();
  int nEvts = atree->GetEntries();
  if (nEvts!=nEvtsPu) {
    std::cout << " -- Warning, not all events found! nEvts = " << nEvts << " nEvtsPu=" << nEvtsPu << std::endl;
    //return;
  }

  int ievtpu(0);

  unsigned nSkipped = 0;

  for (int ievt(0); ievt<nEvts && ievtpu<nEvtsPu; ++ievt,++ievtpu){//loop on entries
    atree->GetEntry(ievt);
    atreePu->GetEntry(ievtpu);
    
    if (ievt%100 == 0) {
      //if (nSkipped < 10) 
      //std::cout << "... Processing entry: evt " << ievt << " idx=" << evtIndex << " / pu evt " << ievtpu << " idx=" << evtIndexPu << std::endl;
    }

    if (evtIndexPu!=evtIndex) {
      //std::cout << " -- Different events found: " << evtIndex << " " << evtIndexPu;
      if (evtIndexPu<evtIndex) {
	//std::cout << ". Skipping in 140pu tree." << std::endl;
	ievt--;
      }
      else {
	//std::cout << ". Skipping in 0 pu tree." << std::endl;
	ievtpu--;
      }
      nSkipped++;
      continue;
    }
    //if (energySR[1][14][2]<energySR[0][14][2] && evtIndex<10) {
    //std::cout << "nopu " << ievt << " idx=" << evtIndex << " pu " << ievtpu << " idx=" << evtIndexPu << std::endl
    //		<< "SR2 E[14] = " << energySR[0][14][2] << " " << energySR[1][14][2] << std::endl;
    //	}

    double puE[nSR];
    double E[nSR];
    double subtrE[nSR];
    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
      puE[iSR]=0;
      E[iSR]=0;
      subtrE[iSR]=0;
      for (unsigned iL(0);iL<nLayers;++iL){
	E[iSR]   += absweight[iL]*energySR[0][iL][iSR];
	subtrE[iSR]   += absweight[iL]*subtractedenergySR[iL][iSR];
	puE[iSR] += absweight[iL]*energySR[1][iL][iSR];
      }//loop on layers

      if (iSR==4 && E[iSR]>puE[iSR] && evtIndex<100){
	std::cout << "nopu " << ievt << " idx=" << evtIndex << " pu " << ievtpu << " idx=" << evtIndexPu << std::endl
		  << "noPU E = " << E[iSR]
		  << " PU E = " << puE[iSR]
		  << std::endl;
      }
      //cut outliers to have correct mean
      if (fabs(puE[iSR]-E[iSR])/calib[iSR]>20) continue;
      p_sigma[iSR]->Fill((puE[iSR]-E[iSR])/calib[iSR]);
      p_EvspuE[iSR]->Fill(puE[iSR]/calib[iSR],E[iSR]/calib[iSR]);
    }//loop on SR
    
    //p_sigma[nSR+2]->Fill((totalEpu-totalE)/tanh(etaval)*1./calib[nSR+2]);
    if (p_sigma[nSR]) p_sigma[nSR]->Fill((totalEpu-totalE)/calib[nSR]);
    
  }//loop on events

  for (unsigned iSR(0); iSR<nSR+1;++iSR){
    if (p_sigma[iSR]) {
      std::cout << " -- SR " << iSR << " - Found " << p_sigma[iSR]->GetEntries() << " events for pu subtraction sigma! nEvts = " << nEvts << " nEvtsPu=" << nEvtsPu << " nskipped = " << nSkipped << std::endl;
      noisePu.push_back(fabs(p_sigma[iSR]->GetRMS()));
    }
  }

  fout->Write();
  TCanvas *mycSig = new TCanvas("mycSig","puSigma",1);
  //mycSig->Divide(2,2);

  mycSig->cd();
  p_EvspuE[2]->Draw("colz");



  return true;

};
