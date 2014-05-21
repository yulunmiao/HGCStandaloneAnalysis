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
#include "TRandom3.h"
#include "Math/Vector4D.h"

#include "fastjet/ClusterSequence.hh"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"

using namespace fastjet;

template <class T>
void extractParameterFromStr(std::string aStr,T vec){ 
  if (aStr == "") return;
  std::vector<std::string> layVec;
  boost::split( layVec, aStr, boost::is_any_of(","));

  for (unsigned iE(0); iE<layVec.size(); ++iE){//loop on elements
    std::vector<std::string> lPair;
    boost::split( lPair, layVec[iE], boost::is_any_of(":"));
    if (lPair.size() != 2) {
      std::cout << " -- Wrong string for parameter given as input:" << layVec[iE] << " Try again, expecting exactly one symbol \":\" between two \",\" ..." << std::endl;
      exit(1);
    }
    std::vector<std::string> lLay;
    boost::split( lLay, lPair[0], boost::is_any_of("-"));
    if (lLay.size() > 2) {
      std::cout << " -- Wrong string for granularities given as input:" << lPair[0] << " Try again, expecting at most one symbol \"-\"." << std::endl;
      exit(1);
    }
    unsigned beginIdx =  atoi(lLay[0].c_str());
    unsigned endIdx = lLay.size() == 1 ? beginIdx :  atoi(lLay[1].c_str());
    for (unsigned iL(beginIdx); iL<endIdx+1; ++iL){
      std::istringstream(lPair[1])>>vec[iL];
    }
  }//loop on elements
}


void simToDigi(HGCSSRecoHit & aHit, 
	       const double & aNoise, 
	       const double & aMeVtoMip, 
	       const unsigned & aMipToADC){
  //convert to MIP
  double oldE = aHit.energy()*aMeVtoMip;
  double newE = oldE+aNoise;
  aHit.noiseFraction(newE > 0 ? 
		     (fabs(aNoise) < oldE ? aNoise/oldE : 1) : -1);
  aHit.energy(newE > 0 ? newE : 0);
  //round to the lower integer, not nearer...
  aHit.adcCounts(static_cast<unsigned>(aHit.energy()*aMipToADC));
}

bool digiToReco(HGCSSRecoHit & aRecHit,
		const unsigned & aMipToADC,
		const unsigned & aThreshold){
  aRecHit.energy(aRecHit.adcCounts()*1.0/aMipToADC);
  if (aRecHit.adcCounts() > aThreshold) 
    return true;
  return false;
}

int main(int argc, char** argv){//main  

  /////////////////////////////////////////////////////////////
  //parameters
  /////////////////////////////////////////////////////////////
  const unsigned nReqA = 8;
  const unsigned nPar = static_cast<unsigned>(argc);
  if (nPar < nReqA) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"<< std::endl
	      << "<full path to input file>"<< std::endl
	      << "<full path to output file>"<< std::endl
	      << "<granularities \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
	      << "<noise (in Mips) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
	      << "<threshold (in ADC counts) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
	      << "<mipToADC (in ADC counts/MIP)> " << std::endl
	      << std::endl
	      << "<optional: randomSeed (default=0)> "  << std::endl
	      << "<optional: debug (default=0)>" << std::endl
	      << "<optional: save sim hits (default=1)> " << std::endl
	      << "<optional: save digi hits (default=0)> " << std::endl
	      << "<optional: maxTime (default=2ns)>" << std::endl
	      << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  std::string inFilePath = argv[2];
  std::string outFilePath = argv[3];
  std::string granulStr = argv[4];
  std::string noiseStr = argv[5];
  std::string threshStr = argv[6];
  unsigned pMipToADC = 50;
  std::istringstream(argv[7])>>pMipToADC;

  std::cout << " ----------------------------------------" << std::endl
	    << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inFilePath << std::endl
	    << " -- Output file path: " << outFilePath << std::endl
	    << " -- Processing " << pNevts << " events." << std::endl
	    << " -- Granularities: " << granulStr << std::endl
	    << " -- noise: " << noiseStr << std::endl
	    << " -- thresholds: " << threshStr << std::endl
	    << " -- Mip to ADC conversion: " << pMipToADC << " ADC counts / MIP." << std::endl;
	    
  double pMaxTime = 2.0;//ns
  unsigned debug = 0;
  unsigned pSeed = 0;
  bool pSaveDigis = 0;
  bool pSaveSims = 1;
  if (nPar > nReqA) std::istringstream(argv[nReqA+1])>>pSeed;
  if (nPar > nReqA+1) {
    debug = atoi(argv[nReqA+1]);
    std::cout << " -- DEBUG output is set to " << debug << std::endl;
  }
  if (nPar > nReqA+2) std::istringstream(argv[nReqA+2])>>pSaveDigis;
  if (nPar > nReqA+3) std::istringstream(argv[nReqA+3])>>pSaveSims;
  if (nPar > nReqA+4) std::istringstream(argv[nReqA+4])>>pMaxTime;
  
  std::cout << " -- Random seed will be set to : " << pSeed << std::endl;
  if (pSaveDigis) std::cout << " -- DigiHits are saved." << std::endl;
  if (!pSaveSims) std::cout << " -- SimHits are not saved." << std::endl;
  std::cout << " -- Max time integrated: " << pMaxTime << " ns." << std::endl
	    << " ----------------------------------------" << std::endl;
  

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned nLayers = 54; //64;//Calice 54;//Scint 9;//HCAL 33;//All 64
  const double xWidth = 500; //200
  const unsigned nEcalLayers = 31;//31;
  const unsigned nHcalSiLayers = 47;//concept 24;//calice 47

  bool concept = false;

  bool makeJets = false;
  // choose a jet definition
  double R = 0.5;
  JetDefinition jet_def(antikt_algorithm, R);

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  unsigned granularity[nLayers];
  double pNoiseInMips[nLayers];
  unsigned pThreshInADC[nLayers];
  for (unsigned iL(0); iL<nLayers; ++iL){
    granularity[iL] = 1;
    pNoiseInMips[iL] = 0.1;
    pThreshInADC[iL] = static_cast<unsigned>(5*pNoiseInMips[iL]*pMipToADC+0.5);
  }

  extractParameterFromStr<unsigned[(unsigned)nLayers]>(granulStr,granularity);
  extractParameterFromStr<double[nLayers]>(noiseStr,pNoiseInMips);
  extractParameterFromStr<unsigned[(unsigned)nLayers]>(threshStr,pThreshInADC);

  //unsigned nbCells = 0;

  std::cout << " -- Granularities and noise are setup like this:" << std::endl;
  for (unsigned iL(0); iL<nLayers; ++iL){
    std::cout << "Layer " ;
    if (iL<10) std::cout << " ";
    std::cout << iL << " : " << granularity[iL] << ", " << pNoiseInMips[iL] << " mips, " << pThreshInADC[iL] << " adc - ";
    if (iL%5==4) std::cout << std::endl;
    
    //nbCells += N_CELLS_XY_MAX/(granularity[iL]*granularity[iL]);
  }
        
  //std::cout << " -- Total number of cells = " << nbCells << std::endl;

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(pSeed);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  //initialise calibration class

  HGCSSCalibration mycalib(inFilePath,concept);

  std::cout << " -- nEcalLayers = " << nEcalLayers 
	    << ", mip weights = " << mycalib.mipWeight(0) << " " << mycalib.mipWeight(1) << " " << mycalib.mipWeight(2)
	    << ", GeV weights = " << mycalib.gevWeight(0) << " offset " << mycalib.gevOffset(0)
	    << std::endl
	    << " -- nHcalSiLayers  = " << nHcalSiLayers 
	    << ", mip weights = " << mycalib.mipWeight(3) << " " << mycalib.mipWeight(4) << " " << mycalib.mipWeight(5) 
	    << ", GeV weights = " << mycalib.gevWeight(1) << " " << mycalib.gevWeight(2) << " offsets: " << mycalib.gevOffset(1) << " " <<   mycalib.gevOffset(2)
	    << std::endl
	    << " -- conversions: HcalToEcalConv = " <<mycalib.HcalToEcalConv() << " BHcalToFHcalConv = " << mycalib.BHcalToFHcalConv() << std::endl
	    << " -----------------------------------" << std::endl;

  /////////////////////////////////////////////////////////////
  //output
  /////////////////////////////////////////////////////////////

  std::ostringstream outputStr;
  outputStr << outFilePath << "/DigiPFcal" ;
  // outputStr << "_" << granulStr << "_";
  // if (noiseStr.size()>0) outputStr << noiseStr << "_";
  // else outputStr << pNoiseInMips[0] << "_";
  // if (threshStr.size()>0) outputStr << threshStr << "_";
  // else outputStr << pThreshInADC[0] << "_";
  // outputStr << pMipToADC << "ADCPerMip";
  // //outputStr << "_" << pMaxTime << "ns";
  if (pSaveDigis)  outputStr << "_withDigiHits";
  if (!pSaveSims)  outputStr << "_withoutSimHits";
  outputStr << ".root";
  
  TFile *outputFile = TFile::Open(outputStr.str().c_str(),"RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << outputStr.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- File will be saved as " << outputStr.str() << std::endl;
  }

  TTree *outputTree = new TTree("RecoTree","HGC Standalone simulation reco tree");
  HGCSSSimHitVec lSimHits;
  HGCSSRecoHitVec lDigiHits;
  HGCSSRecoHitVec lRecoHits;
  HGCSSRecoJetVec lCaloJets;
  unsigned maxSimHits = 0;
  unsigned maxRecHits = 0;
  unsigned maxRecJets = 0;
  unsigned lEvent = 0;
  float lVolX0 = 0;
  float lVolLambda = 0;
  outputTree->Branch("event",&lEvent);
  outputTree->Branch("volX0",&lVolX0);
  outputTree->Branch("volLambda",&lVolLambda);
  if (pSaveSims) outputTree->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&lSimHits);
  if (pSaveDigis) outputTree->Branch("HGCSSDigiHitVec","std::vector<HGCSSRecoHit>",&lDigiHits);
  outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);
  if (makeJets) outputTree->Branch("HGCSSRecoJetVec","std::vector<HGCSSRecoJet>",&lCaloJets);
  TH1F * p_noise = new TH1F("noiseCheck",";noise (MIPs)",100,-2,2);

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::string inputStr = inFilePath ;
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

  float event=0;
  float volNb=0;
  float volX0=0;
  //float volLambda = 0;
  std::vector<HGCSSSimHit> * hitvec = 0;

  inputTree->SetBranchAddress("event",&event);
  inputTree->SetBranchAddress("volNb",&volNb);
  inputTree->SetBranchAddress("volX0trans",&volX0);
  //inputTree->SetBranchAddress("volLambda",&volLambda);
  inputTree->SetBranchAddress("HGCSSSimHitVec",&hitvec);
    
  /////////////////////////////////////////////////////////////
  //Loop on events
  /////////////////////////////////////////////////////////////


  const unsigned nEvts = (pNevts > inputTree->GetEntries()/nLayers || pNevts==0) ? static_cast<unsigned>(inputTree->GetEntries()/nLayers) : pNevts;

  std::cout << "- Processing = " << nEvts  << " events out of " << inputTree->GetEntries()/nLayers << std::endl;

  //create map used to assemble hits per event.
  std::map<unsigned,HGCSSRecoHit> lHitMap;
  std::pair<std::map<unsigned,HGCSSRecoHit>::iterator,bool> isInserted;
  std::vector<PseudoJet> lParticles;

  for (unsigned ientry(0); ientry<nEvts*nLayers; ++ientry){//loop on entries

    unsigned ievt =  ientry/nLayers;

    inputTree->GetEntry(ientry);
    lEvent = event;
    lVolX0 = volX0;
    //lVolLambda = volLambda;
    unsigned layer = volNb;
    
    if (debug>1) {
      std::cout << " **DEBUG** Processing evt " << ievt << " layer " << layer << std::endl;
    }
    else if (debug == 1 && layer==0) std::cout << "... Processing event: " << ievt << std::endl;
    else if (ievt%50 == 0 && layer==0) std::cout << "... Processing event: " << ievt << std::endl;
    

    for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*hitvec)[iH];

      //do not save hits with 0 energy...
      if (lHit.energy()>0) lSimHits.push_back(lHit);

      //C-AMM: TO DO ?
      //discard simhits not in the right time window
      //but needs to take into account propagation of the particle from interaction point...
      //for the moment accept everything.
      //In CALICE, shaper time much larger means everything is integrated.
      //if (lHit.time() > pMaxTime) continue;

      HGCSSRecoHit lRecHit(lHit,granularity[layer]);
      //double posx = lHit.get_x();
      //double posy = lHit.get_y();

      isInserted = lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lRecHit.fullcellid(),lRecHit));
      if (!isInserted.second) isInserted.first->second.Add(lHit);

    }//loop on input simhits

    if (debug>1) {
      std::cout << " **DEBUG** simhits = " << (*hitvec).size() << " " << lSimHits.size() << " recohits = " << lHitMap.size() << std::endl;
    }

    //create digihits, everywhere to have also pure noise.
    for (unsigned iX(0); iX<xWidth/(CELL_SIZE_X*granularity[layer])/2;++iX){
      for (unsigned iY(0); iY<xWidth/(CELL_SIZE_Y*granularity[layer])/2;++iY){
	HGCSSRecoHit lNoiseHit;
	lNoiseHit.layer(layer);
	lNoiseHit.encodeCellId(true,true,iX,iY,granularity[layer]);
	lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lNoiseHit.fullcellid(),lNoiseHit));
	lNoiseHit.encodeCellId(true,false,iX,iY,granularity[layer]);
	lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lNoiseHit.fullcellid(),lNoiseHit));
	lNoiseHit.encodeCellId(false,true,iX,iY,granularity[layer]);
	lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lNoiseHit.fullcellid(),lNoiseHit));
	lNoiseHit.encodeCellId(false,false,iX,iY,granularity[layer]);
	lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lNoiseHit.fullcellid(),lNoiseHit));
      }
    }

    
    if (layer == nLayers-1){
      //save event
      std::map<unsigned,HGCSSRecoHit>::iterator lIter = lHitMap.begin();
      lDigiHits.reserve(lHitMap.size());

      for (; lIter != lHitMap.end(); ++lIter){
	HGCSSRecoHit & lHit = lIter->second;
	//convert to MIP and add noise.
	double lNoise = lRndm->Gaus(0,pNoiseInMips[layer]);
	simToDigi(lHit,lNoise,mycalib.MeVToMip(layer),pMipToADC);
	lDigiHits.push_back(lHit);
	p_noise->Fill(lNoise);
      }
      //apply threshold
      lRecoHits.reserve(lDigiHits.size());
      for (unsigned iH(0); iH<lDigiHits.size(); ++iH){
	//copy digihit to modify it
	HGCSSRecoHit lHit = lDigiHits[iH];
	if (digiToReco(lHit,pMipToADC,pThreshInADC[layer])){
	  lRecoHits.push_back(lHit);
	  //TOFIX: inverse y and z to have eta=0...
	  if (lHit.get_z()!=0) lParticles.push_back( PseudoJet(lHit.px(),lHit.pz(),lHit.py(),lHit.E()));
	  //ROOT::Math::XYZTVector lcheck(lHit.get_x(),lHit.get_y(),lHit.get_z(),0);
	  // std::cout << lHit.get_x() << " " 
	  // 	    << lHit.get_y() << " " 
	  // 	    << lHit.get_z() << " " 
	  // 	    << lHit.px() << " " 
	  // 	    << lHit.py() << " "
	  // 	    << lHit.pz() << " "
	  // 	    << lHit.E() << " " 
	  // 	    << lHit.eta() << " ("
	  //   //<< lcheck.eta() << ") "
	  // 	    << lHit.phi() //<< "("
	  //   //<< lcheck.phi() << ")"
	  // 	    << std::endl;
	}
      }

      if (debug) {
	std::cout << " **DEBUG** sim-digi-reco hits = " << lSimHits.size() << "-" << lDigiHits.size() << "-" << lRecoHits.size() << std::endl;
      }


      if (makeJets){//makeJets

	// run the clustering, extract the jets
	ClusterSequence cs(lParticles, jet_def);
	std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());


	// print the jets
	std::cout <<   "-- evt " << ievt << ": found " << jets.size() << " Jets." << std::endl;
	for (unsigned i = 0; i < jets.size(); i++) {
	  const PseudoJet & lFastJet = jets[i];
	  //TOFIX // inverted y and z...
	  HGCSSRecoJet ljet(lFastJet.px(),
			    lFastJet.pz(),
			    lFastJet.py(),
			    lFastJet.E());
	  if (lFastJet.has_constituents()) ljet.nConstituents(lFastJet.constituents().size());
	  if (lFastJet.has_area()){
	    ljet.area(lFastJet.area());
	    ljet.area_error(lFastJet.area_error());
	  }

	  lCaloJets.push_back(ljet);
	  std::cout << " -------- jet " << i << ": "
		    << lFastJet.E() << " " 
		    << lFastJet.perp() << " " 
		    << lFastJet.rap() << " " << lFastJet.phi() << " "
		    << lFastJet.constituents().size() << std::endl;
	  // std::vector<PseudoJet> constituents = lFastJet.constituents();
	  // for (unsigned j = 0; j < constituents.size(); j++) {
	  //   std::cout << "    constituent " << j << "'s pt: " << constituents[j].perp()
	  // 	      << std::endl;
	  // }
	}

      }//makeJets
      
      outputTree->Fill();
      //reserve necessary space and clear vectors.
      if (lSimHits.size() > maxSimHits) {
	maxSimHits = 2*lSimHits.size();
	std::cout << " -- INFO: event " << ievt << " maxSimHits updated to " << maxSimHits << std::endl;
      }
      if (lRecoHits.size() > maxRecHits) {
	maxRecHits = 2*lRecoHits.size();
	std::cout << " -- INFO: event " << ievt << " maxRecHits updated to " << maxRecHits << std::endl;
      }
      if (lCaloJets.size() > maxRecJets) {
	maxRecJets = 2*lCaloJets.size();
	std::cout << " -- INFO: event " << ievt << " maxRecJets updated to " << maxRecJets << std::endl;
      }
      lSimHits.clear();
      lDigiHits.clear();
      lRecoHits.clear();
      lCaloJets.clear();
      lHitMap.clear();
      lParticles.clear();
      lSimHits.reserve(maxSimHits);
      lParticles.reserve(maxRecHits);
      lCaloJets.reserve(maxRecJets);
    }
    
  }//loop on entries

  outputFile->cd();
  outputTree->Write();
  p_noise->Write();
  outputFile->Close();

  return 0;

}//main
