#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "Math/Vector4D.h"

#include "fastjet/ClusterSequence.hh"

#include "HGCSSEvent.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

using namespace fastjet;

template <class T>
void extractParameterFromStr(std::string aStr,T & vec){ 
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


int main(int argc, char** argv){//main  

  /////////////////////////////////////////////////////////////
  //parameters
  /////////////////////////////////////////////////////////////
  const unsigned nReqA = 7;
  const unsigned nPar = static_cast<unsigned>(argc);
  if (nPar < nReqA) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"<< std::endl
	      << "<full path to input file>"<< std::endl
	      << "<full path to output file>"<< std::endl
	      << "<granularities \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
	      << "<noise (in Mips) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
	      << "<threshold (in ADC counts) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
	      << std::endl
	      << "<optional: randomSeed (default=0)> "  << std::endl
	      << "<optional: debug (default=0)>" << std::endl
	      << "<optional: save sim hits (default=0)> " << std::endl
	      << "<optional: save digi hits (default=0)> " << std::endl
	      << "<optional: make jets (default=0)> " << std::endl
	      << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  std::string inFilePath = argv[2];
  std::string outFilePath = argv[3];
  std::string granulStr = argv[4];
  std::string noiseStr = argv[5];
  std::string threshStr = argv[6];

  std::cout << " ----------------------------------------" << std::endl
	    << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inFilePath << std::endl
	    << " -- Output file path: " << outFilePath << std::endl
	    << " -- Processing " << pNevts << " events." << std::endl
	    << " -- Granularities: " << granulStr << std::endl
	    << " -- noise: " << noiseStr << std::endl
	    << " -- thresholds: " << threshStr << std::endl
    ;
	    
  unsigned debug = 0;
  unsigned pSeed = 0;
  bool pSaveDigis = 0;
  bool pSaveSims = 0;
  bool pMakeJets = false;
  if (nPar > nReqA) std::istringstream(argv[nReqA])>>pSeed;
  if (nPar > nReqA+1) {
    debug = atoi(argv[nReqA+1]);
    std::cout << " -- DEBUG output is set to " << debug << std::endl;
  }
  if (nPar > nReqA+2) std::istringstream(argv[nReqA+2])>>pSaveDigis;
  if (nPar > nReqA+3) std::istringstream(argv[nReqA+3])>>pSaveSims;
  if (nPar > nReqA+4) std::istringstream(argv[nReqA+4])>>pMakeJets;
  
  std::cout << " -- Random seed will be set to : " << pSeed << std::endl;
  if (pSaveDigis) std::cout << " -- DigiHits are saved." << std::endl;
  if (pSaveSims) std::cout << " -- SimHits are saved." << std::endl;
  if (pMakeJets) std::cout << " -- Making jets." << std::endl;
  std::cout << " ----------------------------------------" << std::endl;
  
  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
  //for HGCAL, true means only 12 FHCAL layers considered (24 are simulated)
  bool concept = true;

  // choose a jet definition
  double R = 0.5;
  JetDefinition jet_def(antikt_algorithm, R);

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  //initialise detector
  HGCSSDetector & myDetector = theDetector();
  bool isScintOnly =  inFilePath.find("version22")!=inFilePath.npos;
  bool isHCALonly = inFilePath.find("version21")!=inFilePath.npos || isScintOnly;
  bool isCaliceHcal = inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  unsigned indices[7] = {0,0,0,0,0,0,0};
  //fill layer indices
  if (isScintOnly) {
    indices[4] = 0;
    indices[5] = 9;
    indices[6] = 9;
  }
  else if (isCaliceHcal) {
    indices[3] = 0;
    indices[4] = 38;
    indices[5] = 47;
    indices[6] = 54;
  }
  else if (isHCALonly) {
    indices[3] = 0;
    indices[4] = 24;
    indices[5] = 33;
    indices[6] = 33;
  }
  else {
    indices[0] = 0;
    indices[1] = 11;
    indices[2] = 21;
    indices[3] = 31;
    indices[4] = 55;
    indices[5] = 64;
    indices[6] = 64;
  }

  myDetector.buildDetector(indices,concept,isCaliceHcal);

  //initialise calibration class
  HGCSSCalibration mycalib(inFilePath);

  const unsigned nLayers = myDetector.nLayers();

  HGCSSGeometryConversion geomConv(inFilePath);
  //const double xWidth = geomConv.getXYwidth();

  HGCSSDigitisation myDigitiser;

  std::vector<unsigned> granularity;
  granularity.resize(nLayers,1);
  std::vector<double> pNoiseInMips;
  pNoiseInMips.resize(nLayers,0.12);
  std::vector<unsigned> pThreshInADC;
  pThreshInADC.resize(nLayers,25);

  extractParameterFromStr<std::vector<unsigned> >(granulStr,granularity);
  extractParameterFromStr<std::vector<double> >(noiseStr,pNoiseInMips);
  extractParameterFromStr<std::vector<unsigned> >(threshStr,pThreshInADC);

  //unsigned nbCells = 0;

  std::cout << " -- Granularities and noise are setup like this:" << std::endl;
  for (unsigned iL(0); iL<nLayers; ++iL){
    std::cout << "Layer " ;
    if (iL<10) std::cout << " ";
    std::cout << iL << " : " << granularity[iL] << ", " << pNoiseInMips[iL] << " mips, " << pThreshInADC[iL] << " adc - ";
    if (iL%5==4) std::cout << std::endl;
    myDigitiser.setNoise(iL,pNoiseInMips[iL]);
    //nbCells += N_CELLS_XY_MAX/(granularity[iL]*granularity[iL]);
  }
        
  //std::cout << " -- Total number of cells = " << nbCells << std::endl;

  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(pSeed);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  /////////////////////////////////////////////////////////////
  //output
  /////////////////////////////////////////////////////////////

  std::ostringstream outputStr;
  outputStr << outFilePath << "/DigiPFcal" ;
  if (pSaveDigis)  outputStr << "_withDigiHits";
  if (pSaveSims)  outputStr << "_withSimHits";
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
  HGCSSEvent lEvent;
  outputTree->Branch("HGCSSEvent",&lEvent);
  if (pSaveSims) outputTree->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&lSimHits);
  if (pSaveDigis) outputTree->Branch("HGCSSDigiHitVec","std::vector<HGCSSRecoHit>",&lDigiHits);
  outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);
  if (pMakeJets) outputTree->Branch("HGCSSRecoJetVec","std::vector<HGCSSRecoJet>",&lCaloJets);
  TH1F * p_noise = new TH1F("noiseCheck",";noise (MIPs)",100,-5,5);

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

  HGCSSEvent * event=0;
  std::vector<HGCSSSimHit> * hitvec = 0;

  inputTree->SetBranchAddress("HGCSSEvent",&event);
  inputTree->SetBranchAddress("HGCSSSimHitVec",&hitvec);
    
  /////////////////////////////////////////////////////////////
  //Loop on events
  /////////////////////////////////////////////////////////////

  const unsigned nEvts = (pNevts > inputTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(inputTree->GetEntries()) : pNevts;

  std::cout << "- Processing = " << nEvts  << " events out of " << inputTree->GetEntries() << std::endl;

  std::vector<PseudoJet> lParticles;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

    inputTree->GetEntry(ievt);
    lEvent.eventNumber(event->eventNumber());
    const double cellSize = event->cellSize();
    lEvent.cellSize(cellSize);
    
    //sanity check
    if (geomConv.cellSize() != cellSize){
      std::cerr << " -- Warning ! Cellsize is not as expected. Reinitialising 2-D histograms... " << std::endl;
      geomConv.cellSize(cellSize);
      geomConv.initialiseHistos(true);
    }

    //unsigned layer = volNb;
    
    if (debug>0) {
      std::cout << " **DEBUG** Processing evt " << ievt << std::endl;
    }
    else if (ievt%50 == 0) std::cout << "... Processing event: " << ievt << std::endl;
    

    unsigned prevLayer = 10000;
    DetectorEnum type = DetectorEnum::FECAL;
    unsigned subdetLayer=0;
    for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*hitvec)[iH];

      //do not save hits with 0 energy...
      if (lHit.energy()>0 && pSaveSims) lSimHits.push_back(lHit);
      unsigned layer = lHit.layer();
      if (layer != prevLayer){
	const HGCSSSubDetector & subdet = myDetector.subDetector(layer);
	type = subdet.type;
	subdetLayer = layer-subdet.layerIdMin;
	prevLayer = layer;
	if (debug > 1) std::cout << " - layer " << layer << " " << subdet.name << " " << subdetLayer << std::endl;
      }      
      double energy = lHit.energy()*mycalib.MeVToMip(layer);
      double posx = lHit.get_x(cellSize);
      double posy = lHit.get_y(cellSize);
      double posz = lHit.get_z();
      double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);
      if (energy>0 && 
	  lHit.silayer() < geomConv.getNumberOfSiLayers(type)//,lHit.eta()) 
	  )
	geomConv.fill(type,subdetLayer,energy,realtime,posx,posy,posz);

    }//loop on input simhits

    if (debug>1) {
      std::cout << " **DEBUG** simhits = " << (*hitvec).size() << " " << lSimHits.size() << std::endl;
    }

    //create hits, everywhere to have also pure noise
    //digitise
    //apply threshold
    //save
    unsigned nTotBins = 0;
    for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
      TH2D *histE = geomConv.get2DHist(iL,"E");
      TH2D *histTime = geomConv.get2DHist(iL,"Time");
      TH2D *histZ = geomConv.get2DHist(iL,"Z");
      DetectorEnum adet = myDetector.subDetector(iL).type;
      bool isScint = myDetector.subDetector(iL).isScint;
      bool isSi = myDetector.subDetector(iL).isSi;
      nTotBins += histE->GetNbinsX()*histE->GetNbinsY();
      if (pSaveDigis) lDigiHits.reserve(nTotBins);

      double meanZpos = geomConv.getAverageZ(iL);

      if (debug>1){
	std::cout << " -- Layer " << iL << " " << myDetector.subDetector(iL).name << " z=" << meanZpos
		  << " totbins = " << nTotBins << " histE entries = " << histE->GetEntries() << std::endl;
      }

      for (int iX(1); iX<histE->GetNbinsX()+1;++iX){
	for (int iY(1); iY<histE->GetNbinsY()+1;++iY){
	  double digiE = 0;
	  double simE = histE->GetBinContent(iX,iY);
	  double time = 0;
	  if (simE>0) time = histTime->GetBinContent(iX,iY)/simE;

	  bool passTime = myDigitiser.passTimeCut(adet,time);
	  if (!passTime) continue;

	  double posz = 0;
	  if (simE>0) posz = histZ->GetBinContent(iX,iY)/simE;
	  else posz = meanZpos;
	  double x = histE->GetXaxis()->GetBinCenter(iX);
	  double y = histE->GetYaxis()->GetBinCenter(iY);

	  //if (fabs(x) > 500 || fabs(y)>500) std::cout << " x=" << x << ", y=" << y << std::endl;

	  //correct for particle angle in conversion to MIP
	  double simEcor = myDigitiser.mipCor(simE,x,y,posz);
	  digiE = simEcor;

	  if (isScint && simEcor>0) {
	    digiE = myDigitiser.digiE(simEcor);
	  }
	  myDigitiser.addNoise(digiE,iL,p_noise);

	  double noiseFrac = 1.0;
	  if (simEcor>0) noiseFrac = (digiE-simEcor)/simEcor;

	  //for silicon-based Calo
	  unsigned adc = 0;
	  if (isSi){
	    adc = myDigitiser.adcConverter(digiE,adet);
	    digiE = myDigitiser.adcToMIP(adc,adet);
	  }
	  bool aboveThresh = 
	    (isSi && adc > pThreshInADC[iL]) ||
	    (isScint && digiE > pThreshInADC[iL]*myDigitiser.adcToMIP(1,adet));
	  //histE->SetBinContent(iX,iY,digiE);
	  if ((!pSaveDigis && aboveThresh) ||
	      pSaveDigis)
	    {//save hits

	      HGCSSRecoHit lRecHit;
	      lRecHit.layer(iL);
	      lRecHit.energy(digiE);
	      lRecHit.adcCounts(adc);
	      lRecHit.x(x);
	      lRecHit.y(y);
	      lRecHit.z(posz);
	      lRecHit.noiseFraction(noiseFrac);
	      //unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*granularity[iL]));
	      //unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*granularity[iL]));
	      //lRecHit.encodeCellId(x>0,y>0,x_cell,y_cell,granularity[iL]);

	      if (pSaveDigis) lDigiHits.push_back(lRecHit);

	      lRecoHits.push_back(lRecHit);

	      if (pMakeJets){
		if (posz>0) lParticles.push_back( PseudoJet(lRecHit.px(),lRecHit.py(),lRecHit.pz(),lRecHit.E()));
	      }

	    }//save hits
	}//loop on y
      }//loop on x
    }//loop on layers

    if (debug) {
      std::cout << " **DEBUG** sim-digi-reco hits = " << (*hitvec).size() << "-" << lDigiHits.size() << "-" << lRecoHits.size() << std::endl;
    }
    
    
    if (pMakeJets){//pMakeJets
      
      // run the clustering, extract the jets
      ClusterSequence cs(lParticles, jet_def);
      std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
      
      // print the jets
      std::cout <<   "-- evt " << ievt << ": found " << jets.size() << " Jets." << std::endl;
      for (unsigned i = 0; i < jets.size(); i++) {
	const PseudoJet & lFastJet = jets[i];
	//TOFIX // inverted y and z...
	HGCSSRecoJet ljet(lFastJet.px(),
			  lFastJet.py(),
			  lFastJet.pz(),
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
      
    }//pMakeJets
    
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
    geomConv.initialiseHistos();
    lParticles.clear();
    if (pSaveSims) lSimHits.reserve(maxSimHits);
    lRecoHits.reserve(maxRecHits);
    lParticles.reserve(maxRecHits);
    lCaloJets.reserve(maxRecJets);
    
  }//loop on entries

  outputFile->cd();
  outputTree->Write();
  p_noise->Write();
  outputFile->Close();

  return 0;

}//main
