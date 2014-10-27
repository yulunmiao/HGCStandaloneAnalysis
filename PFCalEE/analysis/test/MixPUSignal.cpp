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
#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

int main(int argc, char** argv){//main  

  if (argc < 7) {
    std::cout << " Usage: "
              << argv[0] << " <nEvts to process (0=all)>"
              << " <path to MinBias file>"
              << " <MinBias file name>"
              << " <path to signal file>"
	      << " <signal file name>"
 	      << " <full path to output file>"
	      << " <Number of PU to add (140)>"
              << std::endl;
    return 1;
  }


  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  bool concept = true;

  double etamin = 1.4;
  double etamax = 3.0;
  //double minZ=3170,maxZ=5070;
  

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned pNevts = atoi(argv[1]);
  std::string pilePath = argv[2];
  std::string pileName = argv[3];
  std::string signalPath = argv[4];
  std::string signalName = argv[5];
  std::string outPath = argv[6];
  unsigned nPU = atoi(argv[7]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input minbias file path: " << pilePath << std::endl
	    << " -- minbias file name: " << pileName << std::endl
	    << " -- Input signal file path: " << signalPath << std::endl
	    << " -- signal file name: " << signalName << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Adding Poisson(" << nPU << ") interactions."  << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;


  //***************** Get Signal Tree *******************************//  
  std::ostringstream signalInput;
  signalInput << signalPath << signalName ;
  TFile *signalFile = TFile::Open(signalInput.str().c_str());
  if (!signalFile) {
    std::cout << " -- Error, input file " << signalInput.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else std::cout << " -- input file " << signalFile->GetName() << " successfully opened." << std::endl;
  
  TTree *signalTree = (TTree*)signalFile->Get("RecoTree");
  if (!signalTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  //***************** Get MinBias Tree *************** **************//
  //If filename contains any wildcard, a bash script will be used to list the files.
  TChain *puTree = new TChain("RecoTree");
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
        puTree->AddFile(pileInput.str().c_str());
     }
     pclose(script);  
     system("rm ./eosls.sh");
   }
   else {
     std::ostringstream pileInput;
     pileInput << pilePath << pileName;
     puTree->AddFile(pileInput.str().c_str());
   }
 
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSRecoHit> * signalhitvec = 0;

  puTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  signalTree->SetBranchAddress("HGCSSRecoHitVec",&signalhitvec);
  
  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  HGCSSInfo * info=(HGCSSInfo*)signalFile->Get("Info");
  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  //models 0,1 or 3.
  //bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  //extract input energy

  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();

  unsigned indices[7] = {0,0,0,0,0,0,0};
  //fill layer indices
  if (versionNumber==22) {
    indices[4] = 0;
    indices[5] = 10;
    indices[6] = 10;
  }
  else if (isCaliceHcal) {
    indices[3] = 0;
    indices[4] = 38;
    indices[5] = 47;
    indices[6] = 54;
  }
  else if (versionNumber==21) {
    indices[3] = 0;
    indices[4] = 24;
    indices[5] = 34;
    indices[6] = 34;
  }
  else if (versionNumber < 20){
    indices[0] = 0;
    indices[1] = versionNumber==8?11:10;
    indices[2] = versionNumber==8?21:20;
    indices[3] = versionNumber==8?31:30;
    indices[4] = indices[3];
    indices[5] = indices[3];
    indices[6] = indices[3];
  }
  else {
    indices[0] = 0;
    indices[1] = 11;
    indices[2] = 21;
    indices[3] = 31;
    indices[4] = 55;
    indices[5] = 65;
    indices[6] = 65;
  }


  myDetector.buildDetector(indices,concept,isCaliceHcal);

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;
  
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Event loop /////////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

 const unsigned nEvts = ((pNevts > signalTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(signalTree->GetEntries()) : pNevts) ;
   
 std::cout << "- Processing = " << nEvts  << " events out of " ;
 std::cout << signalTree->GetEntries() << std::endl;
  
 const unsigned nPuEvts = puTree->GetEntries();
 std::cout << "- Number of PU events available: " << nPuEvts  << std::endl;

 TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }

  outputFile->cd();

  TTree *outputTree = new TTree("PUTree","140 PU tree");
  
  HGCSSRecoHitVec lRecoHits;
  unsigned nPuVtx;
  outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);
  outputTree->Branch("nPuVtx",&nPuVtx);

  for (unsigned ievt(0); ievt<nEvts;++ievt){//loop on output events
    
    lRecoHits.clear();
    lRecoHits.reserve((*signalhitvec).size());
    //get signal event
    signalTree->GetEntry(ievt);
    for(unsigned iH(0); iH<(*signalhitvec).size(); ++iH){
      //copy to fill output vec
      HGCSSRecoHit lHit = (*signalhitvec)[iH];         
      //double posz = lHit.get_z();
      double eta = lHit.eta();
      bool inFid = fabs(eta) > etamin && fabs(eta) < etamax;
      // && posz>minZ && posz<maxZ;
      if (inFid) lRecoHits.push_back(lHit);
    }
    
    //get PU events
    std::vector<unsigned> ipuevt;
    
    //get poisson <140>
    nPuVtx = lRndm.Poisson(nPU);
    ipuevt.resize(nPuVtx,1);

    for (unsigned iV(0); iV<nPuVtx; ++iV){//loop on interactions
      ipuevt[iV] = lRndm.Integer(nPuEvts);
      //get random PU events among available;
      std::cout << " -- PU Random number #" << iV << " for event " << ievt << " is " << ipuevt[iV] << std::endl;

      puTree->GetEntry(ipuevt[iV]);
      lRecoHits.reserve(lRecoHits.size()+(*rechitvec).size());
      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	HGCSSRecoHit lHit = (*rechitvec)[iH];
	//double posz = lHit.get_z();
	double eta = lHit.eta();
	bool inFid = fabs(eta) > etamin && fabs(eta) < etamax;
	// && posz>minZ && posz<maxZ;
	if (inFid) lRecoHits.push_back(lHit);
	
      }//loop on hits
    }//loop on interactions

    //fill tree
    outputTree->Fill();
    
  }//loop on out events

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();

  return 0;
  
}//main
