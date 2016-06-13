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
#include "TVector3.h"

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

using boost::lexical_cast;
namespace po=boost::program_options;

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

  //Input output and config options
  std::string cfg;
   std::string filePath;
  std::string simFileName;
  unsigned pNevts;
  unsigned nRuns;
  std::string outPath;
  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("outPath,o",      po::value<std::string>(&outPath)->required())
    ;
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);
  std::string inFilePath = filePath+simFileName;
  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl;
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;

  HGCSSInfo * info;
  TChain *lSimTree = new TChain("HGCSSTree");
  TFile * simFile = 0;

  if (nRuns == 0){
    if (!testInputFile(inputsim.str(),simFile)) return 1;
    lSimTree->AddFile(inputsim.str().c_str());
    if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
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
      lSimTree->AddFile(lstrsim.str().c_str());
    }
  }

  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  //double calorSizeXY = info->calorSizeXY();
  double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  std::cout //<< " -- Calor size XY = " << calorSizeXY
	    << ", version number = " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;

  //initialise detector
  HGCSSDetector & myDetector = theDetector();
 
  myDetector.buildDetector(versionNumber,true,false);

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;

  //output file
  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();
  TTree *outtree = new TTree("Etruth","Tree to save truth info and abs weights");
  unsigned evtIdx=0;
  double Egamma1=0;
  double Egamma2 = 0;
  std::vector<double> absweight;
  absweight.resize(nLayers,0);
  outtree->Branch("eventIndex",&evtIdx);
  outtree->Branch("Egamma1",&Egamma1);
  outtree->Branch("Egamma2",&Egamma2);
  for (unsigned iL(0); iL<nLayers;++iL){
    std::ostringstream label; 
    label.str("");
    label << "absweight" << iL;
    outtree->Branch(label.str().c_str(),&absweight[iL]);
  }


  //loop on events
  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;

  std::cout << " -- Processing " << nEvts << " events out of " << lSimTree->GetEntries() << std::endl;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (ievt%100 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    lSimTree->GetEntry(ievt);
    
    for(unsigned iL(0); iL<(*ssvec).size(); iL++){
      absweight[iL] = (*ssvec)[iL].volX0trans()/(*ssvec)[1].volX0trans();
      //(*ssvec)[iL].voldEdx()/(*ssvec)[1].voldEdx();
    }
    evtIdx = ievt;
    Egamma1 = 0;
    Egamma2 = 0;

   for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
     if ((*genvec)[iP].trackID()==1) 
       Egamma1 = (*genvec)[iP].E()/1000.;
     if ((*genvec)[iP].trackID()==2) 
       Egamma2 = (*genvec)[iP].E()/1000.;
     
   }

   outtree->Fill();

  }//loop on entries

  outputFile->Write();

  return 0;

}//main
