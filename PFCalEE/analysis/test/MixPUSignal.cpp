#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TChain.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"

int main(int argc, char** argv){//main  

  if (argc < 5) {
    std::cout << " Usage: "
              << argv[0] << " <nEvts to process (0=all)>"
              << " <full path to signal file>"
              << " <signal file name (DigiPFcal.root)>"
              << " <full path to MinBias file>"
              << " <MinBias file name>"
              << std::endl;
    return 1;
  }


  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned nLayers = 1;//64;//Calice 54;//Scint 9;//HCAL 33;//All 64
  const unsigned nEcalLayers = 31;//31;

  double minX=-1700,maxX=1700;
  double minY=-1700,maxY=1700;
  double minZ=3170,maxZ=5070;
  
  unsigned nEvtsOut = 50;

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned pNevts = atoi(argv[1]);
  std::string pilePath = argv[2];
  std::string pileName = argv[3];
  std::string signalPath = argv[4];
  std::string signalName = argv[5];//"PFcal.root"; 
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;
  
  TRandom3 lRndm(0);

  std::cout << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----- Hardcoded configuration is : " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " -- N layers = " << nLayers << std::endl
	    << " -- N ECAL layers = " << nEcalLayers << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " -----------------------------------" << std::endl;

  bool isG4Tree = false;

  //***************** Get Signal Tree *******************************//  
  std::ostringstream signalInput;
  signalInput << signalPath << signalName ;
  TFile *signalFile = TFile::Open(signalInput.str().c_str());
  TTree *signalTree = (TTree*)signalFile->Get("RecoTree");

  //***************** Get MinBias Tree *************** **************//
  //If filename contains any wildcard, a bash script will be used to list the files.
  TChain *lTree = new TChain("RecoTree");
  if(pileName.find("*")!=pileName.npos){ 
     ofstream myscript;
     myscript.open("eosls.sh");
     myscript<<"#!/bin/bash" << std::endl;
     myscript<<"/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls " << pilePath << std::endl; 
     myscript.close();
     FILE *script = popen("bash eosls.sh", "r");
     char eoslsName[100];
     while(fgets(eoslsName, 100, script)) {
        std::ostringstream pileInput;
        std::string temp = std::string(eoslsName).substr(0,strlen(eoslsName)-1);
        if(temp.find("HG")!=temp.npos)continue;
        pileInput << pilePath << temp;
        lTree->AddFile(pileInput.str().c_str());
     }
     pclose(script);  
     system("rm ./eosls.sh");
   }
   else {
     std::ostringstream pileInput;
     pileInput << pilePath << pileName;
     lTree->AddFile(pileInput.str().c_str());
   }
 
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSRecoHit> * signalhitvec = 0;

  lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  signalTree->SetBranchAddress("HGCSSRecoHitVec",&signalhitvec);
  
  const unsigned nEvts = 
    isG4Tree ? 
    ((pNevts > lTree->GetEntries()/nLayers || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()/nLayers) : pNevts) :
    ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;
  
  std::cout << "- Processing = " << nEvts  << " events out of " ;
  if (isG4Tree) std::cout << lTree->GetEntries()/nLayers << std::endl;
  else std::cout << lTree->GetEntries() << std::endl;
  

  TFile *outputFile = TFile::Open("../140PU/PFcal_140PU_EE.root","CREATE");

  if (!outputFile) {
    std::cout << " -- Error, cannot open output file. Exiting..." << std::endl;
    return 1;
  }

  TTree *outputTree = new TTree("PUTree","140 PU tree");
  HGCSSRecoHitVec lRecoHits;
  outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);

  for (unsigned evtOut(0); evtOut<nEvtsOut;++evtOut){//loop on output events
    
    std::vector<unsigned> ipuevt;
    unsigned lVtx = 0;
    
    //get poisson <140>
    lVtx = lRndm.Poisson(10);
    ipuevt.resize(lVtx,1);

    //temp vector to save all hits
    HGCSSRecoHitVec tmpvec[nLayers];

    //Random signal
    unsigned iSig = 0;
    iSig = lRndm.Integer(1000);
    std::cout << " -- PU Random number for signal " << " is " << iSig << std::endl;
    signalTree->GetEntry(iSig);
    for(unsigned iH(0); iH<(*signalhitvec).size(); ++iH){
        HGCSSRecoHit lHit = (*signalhitvec)[iH];         
        tmpvec[0].push_back(lHit);
    }

    for (unsigned iV(0); iV<lVtx; ++iV){//loop on interactions
      ipuevt[iV] = lRndm.Integer(nEvts);
      //get random PU events among available;
      std::cout << " -- PU Random number #" << iV << " for event " << evtOut << " is " << ipuevt[iV] << std::endl;
    }
    
    for (unsigned iV(0); iV<lVtx; ++iV){//loop on interactions
      std::cout << " Processing interaction " << iV << std::endl;
      unsigned ievt = isG4Tree? nLayers*ipuevt[iV] : ipuevt[iV];	
      for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
	if (!isG4Tree && iL>0) continue;

	lTree->GetEntry(ievt+iL);
	//std::cout << iL << " " << (*simhitvec).size() << std::endl;
	for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	  HGCSSRecoHit lHit = (*rechitvec)[iH];
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  double posz = lHit.get_z();
	  

	  bool inFid = posx > minX && posx < maxX &&
	    posy > minY && posy < maxY &&
	    posz > minZ && posz < maxZ;
	  
	  if (inFid) tmpvec[iL].push_back(lHit);
	  
	}//loop on hits
      }//loop on layers
    }//loop on interactions
    for (unsigned iL(0); iL<nLayers; ++iL){
      lRecoHits.clear();
      //std::cout << " -- Layer " << iL << " has " << tmpvec[iL].size() << " hits." << std::endl;
      lRecoHits.reserve(tmpvec[iL].size());
      for (unsigned iH(0); iH<tmpvec[iL].size(); ++iH){
	lRecoHits.push_back(tmpvec[iL][iH]);
      }
      outputTree->Fill();
      tmpvec[iL].clear();
    }
  }//loop on out events

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();

  return 0;
  
  
}//main
