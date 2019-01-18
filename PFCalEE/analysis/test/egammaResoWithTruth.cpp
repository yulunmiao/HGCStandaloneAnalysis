#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSPUenergy.hh"

#include "PositionFit.hh"
#include "SignalRegion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

using boost::lexical_cast;
namespace po=boost::program_options;

struct BigHit{
  int index;
  double Emax;
  double Esum;
  BigHit(int idx,double Em,double Es){
    index=idx;
    Emax=Em;
    Esum=Es;
  };
};

bool testInputFile(std::string input, TFile* & file){
  file = TFile::Open(input.c_str());
  
  if (!file) {
    std::cout << " -- Error, input file " << input.c_str() << " cannot be opened. Skipping..." << std::endl;
    return false;
  }
  else std::cout << " -- input file " << file->GetName() << " successfully opened." << std::endl;
  return true;
};


int main(int argc, char** argv){//main  

  bool doPaul = false;

  //Input output and config options
  std::string cfg;
  bool concept;
  //size of signal region to perform Chi2 position fit.
  //in units of 2.5mm cells to accomodate different granularities
  unsigned nSR;
  //maximum value of residuals to use in error matrix: discard positions that are too far away 
  double residualMax;//mm
  unsigned pNevts;
  std::string filePath;
  std::string digifilePath;
  unsigned nRuns;
  std::string simFileName;
  std::string recoFileName;
  std::string outPath;
  unsigned nSiLayers;
  //0:do just the energies, 1:do fit+energies, 2: do zpos+fit+energies
  unsigned redoStep;
  unsigned debug;
  bool applyPuMixFix;
  unsigned g4trackID;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("concept",        po::value<bool>(&concept)->default_value(true))
    ("nSR",            po::value<unsigned>(&nSR)->default_value(3))
    ("residualMax",    po::value<double>(&residualMax)->default_value(25))
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("digifilePath", po::value<std::string>(&digifilePath)->default_value(""))
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("outPath,o",      po::value<std::string>(&outPath)->required())
    ("nSiLayers",      po::value<unsigned>(&nSiLayers)->default_value(2))
    ("redoStep",       po::value<unsigned>(&redoStep)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("g4trackID",        po::value<unsigned>(&g4trackID)->default_value(1))
    ("applyPuMixFix",  po::value<bool>(&applyPuMixFix)->default_value(false))
    ;

  // ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);



  std::string inFilePath = filePath+simFileName;

  size_t end=outPath.find_last_of(".");
  std::string outFolder = outPath.substr(0,end);


  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Digi Input file path: " << digifilePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Output folder: " << outFolder << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Number cells in signal region for fit: " << nSR << " cells" << std::endl
	    << " -- Residual max considered for filling matrix and fitting: " << residualMax << " mm" << std::endl
	    << " -- Apply PUMix fix? " << applyPuMixFix << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events per run." << std::endl;

  TRandom3 lRndm(1);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;
  std::ostringstream inputrec;
  if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
  else 
    inputrec << digifilePath << "/" << recoFileName;

  std::cout << inputsim.str() << " " << inputrec.str() << std::endl;

  HGCSSInfo * info;

  TChain *lSimTree = new TChain("HGCSSTree");
  TChain *lRecTree = 0;
  
  TFile * simFile = 0;
  TFile * recFile = 0;

  if (recoFileName.find("Digi") != recoFileName.npos) 
    lRecTree = new TChain("RecoTree");
  else lRecTree = new TChain("PUTree");

  if (nRuns == 0){
    if (!testInputFile(inputsim.str(),simFile)) return 1;
    lSimTree->AddFile(inputsim.str().c_str());
    if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
    else {
      std::cout << " -- Error in getting information from simfile!" << std::endl;
      return 1;
    }
    if (!testInputFile(inputrec.str(),recFile)) return 1;
    lRecTree->AddFile(inputrec.str().c_str());
  }
  else {
    for (unsigned i(0);i<nRuns;++i){
      std::ostringstream lstrsim;
      std::ostringstream lstrrec;
      lstrsim << inputsim.str() << "_run" << i << ".root";
      if (testInputFile(lstrsim.str(),simFile)){  
	if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
	else {
	  std::cout << " -- Error in getting information from simfile!" << std::endl;
	  return 1;
	}
      }
      else continue;
      lstrrec << inputrec.str() << "_run" << i << ".root";
      if (!testInputFile(lstrrec.str(),recFile)) continue;
      lSimTree->AddFile(lstrsim.str().c_str());
      lRecTree->AddFile(lstrrec.str().c_str());
    }
  }

  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }


  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  //double calorSizeXY = info->calorSizeXY();
  //double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
    //CAMM-dirty fix
  //const unsigned shape = inputsim.str().find("Diamond")!=inputsim.str().npos? 2: inputsim.str().find("Triangle")!=inputsim.str().npos?3 :1 ;//
  //const double cellSize = 6.496345;//info->cellSize();
  //const double calorSizeXY = 2800*2;//info->calorSizeXY();
  const unsigned shape = info->shape();
  const double cellSize = info->cellSize();
  const double calorSizeXY = info->calorSizeXY();
 
  //if (calorSizeXY<1 || calorSizeXY>6000) calorSizeXY=495;

  bool doHexa = fabs(cellSize-2.5)>0.01;
  //if (!doHexa) cellSize = 4*cellSize;

  //models 0,1 or 3.
  //bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  //extract input energy

  std::cout << " -- Calor size XY = " << calorSizeXY
	    << ", version number = " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();
 
  myDetector.buildDetector(versionNumber,model,true,false,false);
  //myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;


  HGCSSGeometryConversion geomConv(model,cellSize,false,3);
  //set granularity to get cellsize for PU subtraction
  std::vector<unsigned> granularity;
  granularity.resize(nLayers,1);
  geomConv.setGranularity(granularity);
  geomConv.setXYwidth(calorSizeXY);
  geomConv.setVersion(versionNumber);
  //if (doHexa) geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  //else geomConv.initialiseSquareMap(calorSizeXY,cellSize);
  if (shape==2) geomConv.initialiseDiamondMap(calorSizeXY,10.);
  else if (shape==3) geomConv.initialiseTriangleMap(calorSizeXY,10.*sqrt(2.));
  else if (shape==1) geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  else if (shape==4) geomConv.initialiseSquareMap(calorSizeXY,10.);
  //square map for BHCAL
  geomConv.initialiseSquareMap1(1.3,3.0,-1.*TMath::Pi(),TMath::Pi(),2./360.);//eta phi segmentation
  geomConv.initialiseSquareMap2(1.3,3.0,-1.*TMath::Pi(),TMath::Pi(),2./288.);//eta phi segmentation
  if (doPaul) geomConv.initialiseSquareMap(calorSizeXY,cellSize>5?83.:62);

  geomConv.initialiseHistos();


  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();


  ///initialise PU density object

  HGCSSPUenergy puDensity("data/EnergyDensity.dat");


    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////// positionFit /////////////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
  
  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;
  
  std::cout << " -- Processing " << nEvts << " events out of " << lSimTree->GetEntries() << " " << lRecTree->GetEntries() << std::endl;

  //perform first loop over simhits to find z positions of layers
  //PositionFit lChi2Fit(nSR,residualMax,nLayers,nSiLayers,applyPuMixFix,debug);
  //lChi2Fit.initialise(outputFile,"PositionFit",outFolder,geomConv,puDensity);
  //if (!lChi2Fit.getZpositions(versionNumber)) 
  //lChi2Fit.getZpositions(versionNumber,lSimTree,500);

  //std::cout << " -- positionfit initilisation done." << std::endl;
  
  //perform second loop over events to find positions to fit and get energies

  std::vector<double> zpos;
  zpos.resize(myDetector.nLayers(),0);
  for (unsigned iL(0); iL<myDetector.nLayers(); ++iL){//loop on layers
    zpos[iL] = myDetector.sensitiveZ(iL);
  }


  SignalRegion SignalEnergy(outFolder, nLayers, zpos, nEvts, geomConv, puDensity,applyPuMixFix,versionNumber,doHexa,g4trackID);
  SignalEnergy.initialise(outputFile,"Energies");

  std::cout << " -- sigenergy initialisation done." << std::endl;

  //loop on events
  HGCSSEvent * event = 0;
  HGCSSEvent * eventRec = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;

  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  lRecTree->SetBranchAddress("HGCSSEvent",&eventRec);
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lRecTree->GetBranch("nPuVtx")) lRecTree->SetBranchAddress("nPuVtx",&nPuVtx);

  unsigned ievtRec = 0;
  unsigned nSkipped = 0;
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (ievtRec>=lRecTree->GetEntries()) continue;
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievtRec);
    //std::cout << " Getting entries " << ievt << " " << ievtRec;
    if (nPuVtx>0 && eventRec->eventNumber()==0 && event->eventNumber()!=0) {
      std::cout << " skip !" << ievt << " " << ievtRec << std::endl;
      nSkipped++;
      continue;
    }

    //modify hits on the fly
    //CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (doPaul){
      std::map<int,BigHit> lBoxMap[nLayers];
      std::pair<std::map<int,BigHit>::iterator,bool> lInsert;
      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	HGCSSRecoHit & lHit = (*rechitvec)[iH];
	unsigned lay = lHit.layer();
	if (lay<nLayers) {
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  int id = geomConv.squareMap()->FindBin(posx,posy);
	  //if (id==7036 && lay==6) std::cout << " Hit " << iH << " layer " << lHit.layer() << " id " << id << " energy " << lHit.energy() << std::endl;
	  lInsert = lBoxMap[lay].insert(std::pair<int,BigHit>(id,BigHit(iH,lHit.energy(),lHit.energy())));
	  std::map<int,BigHit>::iterator lEle = lInsert.first;
	  if (!lInsert.second) {
	    int & maxidx = lEle->second.index;
	    double & maxE = lEle->second.Emax;
	    double & sumE = lEle->second.Esum;
	    //if (id==7036 && lay==6) std::cout << "already filled with: " << maxidx << " " << maxE << std::endl;
	    if (maxE < lHit.energy()){
	      //change element in map
	      //update hit with sum of energy
	      maxE = lHit.energy();
	      sumE += lHit.energy();
	      lHit.energy(sumE);
	      (*rechitvec)[maxidx].energy(0);
	      maxidx = iH;
	      //if (id==7036 && lay==6) std::cout << "Changed ele in map: " << maxidx << " " << (*rechitvec)[maxidx].energy() << " " << iH << " " << (*rechitvec)[iH].energy() << " maxE " << maxE << " sumE " << sumE << std::endl;
	    }
	    else {
	      //change hit, update element in map
	      sumE += lHit.energy();
	      (*rechitvec)[maxidx].energy(sumE);
	      (*rechitvec)[iH].energy(0);
	      //if (id==7036 && lay==6) std::cout << "Changed hit: " << maxidx << " " << (*rechitvec)[maxidx].energy() << " " << iH << " " << (*rechitvec)[iH].energy() << std::endl;
	    }
	    //(lInsert.first)->second.second += lHit.energy();
	  }
	}
      }
      /*      std::cout << " Afterwards: " << std::endl;
      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	HGCSSRecoHit & lHit = (*rechitvec)[iH];
	unsigned lay = lHit.layer();
	if (lay>=nLayers) continue;
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	int id = geomConv.squareMap()->FindBin(posx,posy);
	if (id==7036 && lay==6) std::cout << " Hit " << iH << " layer " << lHit.layer() << " id " << id << " energy " << lHit.energy() << std::endl;
	}*/

    }

    //std::cout  << std::endl;
    SignalEnergy.fillEnergies(ievt,(*event),(*genvec),(*ssvec),(*simhitvec),(*rechitvec),nPuVtx);
    ievtRec++;

  }//loop on entries

  //finalise

  SignalEnergy.finalise();

  outputFile->Write();
  //outputFile->Close();
  
  std::cout << " - Skipped " << nSkipped << " events." << std::endl;
  std::cout << " - End of egammaResoWithTruth program." << std::endl;

  return 0;
  

}//main
