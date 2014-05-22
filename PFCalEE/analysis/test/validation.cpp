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

#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

int main(int argc, char** argv){//main  

  if (argc < 5) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <path to input file>"
	      << " <suffix to input file>"
	      << " <full path to output file>" 
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
  //for HGCAL, true means only 12 FHCAL layers considered (24 are simulated)
  bool concept = true;

  double minX=-1700,maxX=1700;
  double minY=-1700,maxY=1700;
  double minZ=3170,maxZ=5070;
  //double minX=-510,maxX=510;
  //double minY=-510,maxY=510;
  //double minZ=-1000,maxZ=1000;

  unsigned nX=(maxX-minX)/10,nY=(maxY-minY)/10;
  unsigned nZ=maxZ-minZ;

  // choose a jet definition
  //double R = 0.5;
  //JetDefinition jet_def(antikt_algorithm, R);

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string fileName = argv[3];

  std::string inFilePath = filePath+fileName;

  std::string outPath = argv[4];
  unsigned debug = 0;
  if (argc >5) debug = atoi(argv[5]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

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
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;

  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  
  
  std::cout << " -- 2-D histograms: " << std::endl
	    << " -- X: " << nX << " " << minX << " " << maxX << std::endl
	    << " -- Y: " << nY << " " << minY << " " << maxY << std::endl
	    << " -- Z: " << nZ << " " << minZ << " " << maxZ << std::endl
    ;
  outputFile->cd();

  TH2F *p_EsimvsLayer = new TH2F("p_EsimvsLayer",";layer ; Esim (MIPs)",
				 nLayers,0,nLayers,
				 1000,0,50000);
  TH2F *p_ErecovsLayer = new TH2F("p_ErecovsLayer",";layer ; Ereco (MIPs)",
				  nLayers,0,nLayers,
				  1000,0,50000);
  TH1F *p_timeSim = new TH1F("p_timeSim",";G4 time (ns)",1000,0,1000);
  TH2F *p_HCALvsECAL = new TH2F("p_HCALvsECAL",";ECAL (GeV);HCAL (GeV)",
				500,0,500,
				500,0,500);
  TH2F *p_BHCALvsFHCAL = new TH2F("p_BHCALvsFHCAL",";FHCAL (GeV);BHCAL (GeV)",
				  500,0,500,
				  500,0,500);

  TH1F *p_nGenPart = new TH1F("p_nGenPart","n(genParticles)",200,0,200);
  TH1F *p_genPartId = new TH1F("p_genPartId","pdgid",12000,-6000,6000);

  TH1F *p_nSimHits = new TH1F("p_nSimHits","n(SimHits)",
			      1000,0,5000000);
  TH1F *p_nRecHits = new TH1F("p_nRecHits","n(RecHits)",
			      1000,0,50000);

  TH2F *p_xy[nLayers];
  TH2F *p_recoxy[nLayers];
  TH1F *p_EfracSim[nLayers];
  TH1F *p_EfracReco[nLayers];
  TH1F *p_Esim[nSections];
  TH1F *p_Ereco[nSections];

  std::ostringstream lName;
  for (unsigned iL(0); iL<nLayers; ++iL){
    lName.str("");
    lName << "p_xy_" << iL;
    p_xy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
			nX,minX,maxX,
			nY,minY,maxY);
    lName.str("");
    lName << "p_recoxy_" << iL;
    p_recoxy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
			    nX,minX,maxX,
			    nY,minY,maxY);
    lName.str("");
    lName << "p_EfracSim_" << iL;
    p_EfracSim[iL] = new TH1F(lName.str().c_str(),";integrated sim E_{layer}/E_{total}",101,0,1.01);
    lName.str("");
    lName << "p_EfracReco_" << iL;
    p_EfracReco[iL] = new TH1F(lName.str().c_str(),";integrated reco E_{layer}/E_{total}",101,0,1.01); 
  }
  
  for (unsigned iD(0); iD<nSections; ++iD){
    lName.str("");
    lName << "p_Esim_" << myDetector.detName(iD);
    if (myDetector.detType(iD)==DetectorEnum::BHCAL1 || myDetector.detType(iD)==DetectorEnum::BHCAL2) p_Esim[iD] = new TH1F(lName.str().c_str(),";Esim (MIPs)",2000,0,20000);
    else p_Esim[iD] = new TH1F(lName.str().c_str(),";Esim (MIPs)",2000,0,200000);
    lName.str("");
    lName << "p_Ereco_" << myDetector.detName(iD);
    if (myDetector.detType(iD)==DetectorEnum::BHCAL1 || myDetector.detType(iD)==DetectorEnum::BHCAL2)  p_Ereco[iD] = new TH1F(lName.str().c_str(),";Ereco (MIPs)",2000,0,20000);
    else p_Ereco[iD] = new TH1F(lName.str().c_str(),";Ereco (MIPs)",2000,0,200000);
  }


  std::ostringstream input;
  input << filePath << "HGcal_" << fileName;

  TFile *simFile = TFile::Open(input.str().c_str());

  if (!simFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  TTree *lSimTree = (TTree*)simFile->Get("HGCSSTree");
  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  input.str("");
  input << filePath << "Digi_" << fileName;
  
  TFile *recFile = TFile::Open(input.str().c_str());

  if (!recFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  TTree *lRecTree = (TTree*)recFile->Get("RecoTree");
  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }


  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  
  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;
  
  std::cout << "- Processing = " << nEvts  << " events out of " << lSimTree->GetEntries() << std::endl;
  
  //Initialise histos
  //necessary to have overflows ?
  gStyle->SetOptStat(1111111);
  double EtotSim[nLayers];
  double EtotRec[nLayers];
  
  for (unsigned iL(0);iL<nLayers;++iL){
    EtotSim[iL] = 0;
    EtotRec[iL] = 0;
  }
  double Esim[nSections];
  double Ereco[nSections];
  for (unsigned iD(0); iD<nSections; ++iD){
    Esim[iD] = 0;
    Ereco[iD] = 0;
  }

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%10 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);

    if (debug){
      std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
    }

    for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles
      p_genPartId->Fill((*genvec)[iP].pdgid());
    }//loop on gen particles

    p_nGenPart->Fill((*genvec).size());

    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];
      unsigned layer = lHit.layer();
      unsigned sec =  myDetector.getSection(layer);

      double posx = lHit.get_x();
      double posy = lHit.get_y();
      double posz = lHit.get_z();
      double radius = sqrt(posx*posx+posy*posy);
	
      double energy = lHit.energy()*mycalib.MeVToMip(layer);
      if (debug>1) {
	std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		  << " --  position x,y " << posx << "," << posy << std::endl;
	lHit.Print(std::cout);
      }

      p_xy[layer]->Fill(posx,posy,energy);
      //correct for time of flight
      double lRealTime = mycalib.correctTime(lHit.time(),posx,posy,posz);
      p_timeSim->Fill(lRealTime);
      
      EtotSim[layer] += energy;
      if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << EtotSim[layer];

      Esim[sec] += energy;
      
    }//loop on hits

    p_nSimHits->Fill((*simhitvec).size());

    if (debug)  std::cout << std::endl;

    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      if (debug>1) {
	std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		  << " --  position x,y " << lHit.get_x() << "," << lHit.get_y() << std::endl;
	lHit.Print(std::cout);
      }
      
      double posx = lHit.get_x();
      double posy = lHit.get_y();
      double posz = lHit.get_z();

      double energy = lHit.energy();//in MIP already...
      unsigned layer = lHit.layer();
      unsigned sec =  myDetector.getSection(layer);
      
      p_recoxy[layer]->Fill(posx,posy,energy);
      EtotRec[layer] += energy;
      if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << EtotRec[layer];

      Ereco[sec] += energy;
    }//loop on rechits
    
    p_nRecHits->Fill((*rechitvec).size());

    //fill histos
    double EtmpSim[nSections];
    double EtmpRec[nSections];  
    for (unsigned iD(0); iD<nSections; ++iD){
      EtmpSim[iD] = 0;
      EtmpRec[iD] = 0;
    }
    for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
      if (debug) std::cout << " -- Layer " << iL 
			   << " total sim E = " << EtotSim[iL] 
			   << " total rec E = " << EtotRec[iL] 
			   << std::endl;
      unsigned sec =  myDetector.getSection(iL);
      p_EsimvsLayer->Fill(iL,EtotSim[iL]);
      p_ErecovsLayer->Fill(iL,EtotRec[iL]);
      EtmpSim[sec] += EtotSim[iL];
      EtmpRec[sec] += EtotRec[iL];
      if (Esim[sec]>0) p_EfracSim[iL]->Fill(EtmpSim[sec]/Esim[sec]);
      else p_EfracSim[iL]->Fill(0);
      if (Ereco[sec]>0) p_EfracReco[iL]->Fill(EtmpRec[sec]/Ereco[sec]);
      else p_EfracReco[iL]->Fill(0);
      EtotSim[iL] = 0;
      EtotRec[iL] = 0;
	
    }//loop on layers

    double Eecal = 0;
    if (myDetector.getSection(DetectorEnum::FECAL)<nSections) Eecal += Ereco[myDetector.getSection(DetectorEnum::FECAL)];
    if (myDetector.getSection(DetectorEnum::MECAL)<nSections) Eecal += Ereco[myDetector.getSection(DetectorEnum::MECAL)];
    if (myDetector.getSection(DetectorEnum::BECAL)<nSections) Eecal += Ereco[myDetector.getSection(DetectorEnum::BECAL)];
    double Efhcal = 0;
    if (myDetector.getSection(DetectorEnum::FHCAL)<nSections) Efhcal += Ereco[myDetector.getSection(DetectorEnum::FHCAL)];
    double Ebhcal = 0;
    if (myDetector.getSection(DetectorEnum::BHCAL1)<nSections) Ebhcal += Ereco[myDetector.getSection(DetectorEnum::BHCAL1)];
    if (myDetector.getSection(DetectorEnum::BHCAL2)<nSections) Ebhcal += Ereco[myDetector.getSection(DetectorEnum::BHCAL2)];

    p_HCALvsECAL->Fill(Eecal,Efhcal+Ebhcal);
    p_BHCALvsFHCAL->Fill(Efhcal,Ebhcal);

    for (unsigned iD(0); iD<nSections; ++iD){
      p_Esim[iD]->Fill(Esim[iD]);
      p_Ereco[iD]->Fill(Ereco[iD]);
      Esim[iD]=0;
      Ereco[iD]=0;
    }
  }//loop on entries
  
  //write
  // for (unsigned iL(0); iL<nLayers; ++iL){
  //     outputFile->cd();
  //     p_xy[iL]->Write();
  //     p_recoxy[iL]->Write();
  //     p_EfracSim[iL]->Write();
  //     p_EfracReco[iL]->Write();
  //   }
  // p_timeSim->Write();
  // p_EsimvsLayer->Write();
  // p_ErecovsLayer->Write();
  // p_nGenPart->Write();
  // p_genPartId->Write();
  for (unsigned iD(0); iD<nSections; ++iD){
    //p_Esim[iD]->Write();
    //p_Ereco[iD]->Write();
    std::cout << " -- Summary of sim energies " 
	      << myDetector.detName(iD) << std::endl
	      <<  p_Esim[iD]->GetEntries() 
	      << " mean " << p_Esim[iD]->GetMean() 
	      << " rms " << p_Esim[iD]->GetRMS() 
	      << " rms/mean " << p_Esim[iD]->GetRMS()/p_Esim[iD]->GetMean()
	      << " underflows " << p_Esim[iD]->GetBinContent(0)
	      << " overflows " << p_Esim[iD]->GetBinContent(p_Esim[iD]->GetNbinsX()+1)
	      << std::endl;
    std::cout << " -- Summary of reco energies " 
	      << myDetector.detName(iD) << std::endl
	      <<  p_Ereco[iD]->GetEntries() 
	      << " mean " << p_Ereco[iD]->GetMean() 
	      << " rms " << p_Ereco[iD]->GetRMS() 
	      << " rms/mean " << p_Ereco[iD]->GetRMS()/p_Ereco[iD]->GetMean()
	      << " underflows " << p_Ereco[iD]->GetBinContent(0)
	      << " overflows " << p_Ereco[iD]->GetBinContent(p_Ereco[iD]->GetNbinsX()+1)
	      << std::endl;
  }
  //p_HCALvsECAL->Write();
  //p_BHCALvsFHCAL->Write();
  //p_nSimHits->Write();
  //p_nRecHits->Write();
  
  outputFile->Write();
  
  return 0;


}//main
