#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"

int main(int argc, char** argv){//main  

  /////////////////////////////////////////////////////////////
  //parameters
  /////////////////////////////////////////////////////////////
  const unsigned nRequiredArgs = 4;
  if (static_cast<unsigned>(argc) < nRequiredArgs) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <full path to input file>"
	      << " <granularities \"layer_i-layer_j:factor,layer:factor,...\">"
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string granulStr = argv[3];
  unsigned debug = 0;
  if (static_cast<unsigned>(argc) >nRequiredArgs) {
    debug = atoi(argv[4]);
    std::cout << " -- DEBUG output is enabled " << debug << std::endl;
  }


  const unsigned nLayers = N_LAYERS;
  unsigned granularity[nLayers];
  for (unsigned iL(0); iL<nLayers; ++iL){
    granularity[iL] = 1;
  }

  std::vector<std::string> layVec;
  boost::split( layVec, granulStr, boost::is_any_of(","));

  for (unsigned iE(0); iE<layVec.size(); ++iE){//loop on elements
    std::vector<std::string> lPair;
    boost::split( lPair, layVec[iE], boost::is_any_of(":"));
    if (lPair.size() != 2) {
      std::cout << " -- Wrong string for granularities given as input:" << layVec[iE] << " Try again, expecting exactly one symbol \":\" between two \",\" ..." << std::endl;
      return 1;
    }
    std::vector<std::string> lLay;
    boost::split( lLay, lPair[0], boost::is_any_of("-"));
    if (lLay.size() > 2) {
      std::cout << " -- Wrong string for granularities given as input:" << lPair[0] << " Try again, expecting at most one symbol \"-\"." << std::endl;
      return 1;
    }
    unsigned beginIdx =  atoi(lLay[0].c_str());
    unsigned endIdx = lLay.size() == 1 ? beginIdx :  atoi(lLay[1].c_str());
    for (unsigned iL(beginIdx); iL<endIdx+1; ++iL){
      granularity[iL] = atoi(lPair[1].c_str());
    }
  }//loop on elements

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Processing " << pNevts << " events." << std::endl
	    << " -- Granularities are setup like this:" << std::endl;
  for (unsigned iL(0); iL<nLayers; ++iL){
    std::cout << "Layer " ;
    if (iL<10) std::cout << " ";
    std::cout << iL << " : " << granularity[iL] << " - ";
    if (iL%5==4) std::cout << std::endl;
  }
        

  /////////////////////////////////////////////////////////////
  //output
  /////////////////////////////////////////////////////////////

  std::string outputStr = filePath+"/DigiPFcal.root";
  TFile *outputFile = TFile::Open(outputStr.c_str(),"RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << outputStr << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  TTree *outputTree = new TTree("RecoTree","HGC Standalone simulation reco tree");
  HGCSSSimHitVec lSimHits;
  HGCSSRecoHitVec lRecoHits;
  unsigned maxSimHits = 0;
  unsigned maxRecoHits = 0;
  unsigned lEvent = 0;
  outputTree->Branch("event",&lEvent);
  outputTree->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&lSimHits);
  outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::string inputStr = filePath + "/PFcal.root";
  TFile *inputFile = TFile::Open(inputStr.c_str());
  
  if (!inputFile) {
    std::cout << " -- Error, input file " << inputStr << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  TTree *inputTree = (TTree*)inputFile->Get("HGCSSTree");
  if (!inputTree){
    std::cout << " -- Error, tree HGCSSTree  cannot be opened. Exiting..." << std::endl;
    return 1;
  }


  /////////////////////////////////////////////////////////////
  //input tree
  /////////////////////////////////////////////////////////////

  float event;
  float volNb;
  std::vector<HGCSSSimHit> * hitvec = 0;

  inputTree->SetBranchAddress("event",&event);
  inputTree->SetBranchAddress("volNb",&volNb);
  inputTree->SetBranchAddress("HGCSSSimHitVec",&hitvec);
    
  /////////////////////////////////////////////////////////////
  //Loop on events
  /////////////////////////////////////////////////////////////


  const unsigned nEvts = (pNevts > inputTree->GetEntries()/30. || pNevts==0) ? static_cast<unsigned>(inputTree->GetEntries()/30.) : pNevts;

  std::cout << "- Processing = " << nEvts  << " events out of " << inputTree->GetEntries()/30. << std::endl;

  //create map used to assemble hits per event.
  std::map<unsigned,HGCSSRecoHit> lHitMap;
  std::pair<std::map<unsigned,HGCSSRecoHit>::iterator,bool> isInserted;

  for (unsigned ientry(0); ientry<nEvts*30; ++ientry){//loop on entries

    unsigned ievt =  ientry/30;

    inputTree->GetEntry(ientry);
    lEvent = event;
    unsigned layer = volNb;
    
    if (debug) {
      std::cout << " **DEBUG** Processing evt " << ievt << " layer " << layer << std::endl;
    }
    else if (ievt%50 == 0 && layer==0) std::cout << "... Processing event: " << ievt << std::endl;
    

    for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*hitvec)[iH];
      lSimHits.push_back(lHit);
      HGCSSRecoHit lRecHit(lHit,granularity[layer]);
      //double posx = lHit.get_x();
      //double posy = lHit.get_y();

      isInserted = lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lRecHit.fullcellid(),lRecHit));
      if (!isInserted.second) isInserted.first->second.Add(lHit);

    }//loop on input simhits

    if (debug) {
      std::cout << " **DEBUG** simhits = " << (*hitvec).size() << " " << lSimHits.size() << " recohits = " << lHitMap.size() << std::endl;
    }

    if (layer == nLayers-1){
      //save event
      std::map<unsigned,HGCSSRecoHit>::iterator lIter = lHitMap.begin();
      lRecoHits.reserve(lHitMap.size());
      for (; lIter != lHitMap.end(); ++lIter){
	lRecoHits.push_back(lIter->second);
      }

      outputTree->Fill();
      //reserve necessary space and clear vectors.
      if (lSimHits.size() > maxSimHits) {
	maxSimHits = 2*lSimHits.size();
	std::cout << " -- INFO: event " << ievt << " maxSimHits updated to " << maxSimHits << std::endl;
      }
      if (lRecoHits.size() > maxRecoHits) {
	maxRecoHits = 2*lRecoHits.size();
	std::cout << " -- INFO: event " << ievt << " maxRecoHits updated to " << maxRecoHits << std::endl;
      }
      lSimHits.clear();
      lRecoHits.clear();
      lHitMap.clear();
      lSimHits.reserve(maxSimHits);
      lRecoHits.reserve(maxRecoHits);
    }
    
  }//loop on entries

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();

  return 0;

}//main
