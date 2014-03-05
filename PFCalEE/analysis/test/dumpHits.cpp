#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"

#include "HGCSSSimHit.hh"

int main(int argc, char** argv){//main  

  if (argc < 4) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <full path to input file: root://eoscms//eos/cms/store/user/amagnan/HGCalEEGeant4/Signal>"
	      << " <output path>"
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }
  const unsigned nLayers = 30;

  std::cout << " -- N layers = " << nLayers << std::endl;

  unsigned genEn[]={5,10,25,50,75,100,150,200,300,500,1000};
  //unsigned genEn[]={5,20,50,100,150,200};
  const unsigned nGenEn = sizeof(genEn)/sizeof(unsigned);
  
  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string outPath = argv[3];
  unsigned debug = 0;
  if (argc >4) debug = atoi(argv[4]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output files path: " << outPath << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  for (unsigned iE(0); iE<nGenEn; ++iE){

    std::cout << "- Processing energy : " << genEn[iE] << std::endl;

    std::ostringstream input;
    input << filePath ;
    bool pEOS = true;
    if (pEOS) {
      input << "_e" ;
      input << genEn[iE] << ".root";
    }
    else {
      input << "/e_" ;
      input << genEn[iE] << "/PFcal.root" ;
    }
    TFile *inputFile = TFile::Open(input.str().c_str());

    if (!inputFile) {
      std::cout << " -- Error, input file " << input.str() << " cannot be opened. Trying next one..." << std::endl;
      continue;
    }

    TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
    if (!lTree){
      std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting ..." << std::endl;
      return 1;
    }

    TString genEnStr;
    genEnStr += genEn[iE];
    std::ofstream hitsOut;
    hitsOut.open(outPath+"/hits_"+genEnStr+"GeV.dat");

    if (!hitsOut.is_open()){
      std::cout << " -- Cannot open output file for writting ! Please create directory: " << outPath << std::endl;
      return 1;
    }

    float volNb = 0;
    float volX0 = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;

    lTree->SetBranchAddress("volNb",&volNb);
    lTree->SetBranchAddress("volX0",&volX0);
    lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);

    const unsigned nEvts = 
      (pNevts > lTree->GetEntries()/nLayers || pNevts==0) ? 
      static_cast<unsigned>(lTree->GetEntries()/nLayers) : 
      pNevts;

    std::cout << "- Processing = " << nEvts  << " events out of " ;
    std::cout << lTree->GetEntries()/nLayers << std::endl;


    for (unsigned ievt(0); ievt<nEvts*nLayers; ++ievt){//loop on entries
      if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
      else if (ievt%1000 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
      lTree->GetEntry(ievt);

      if (debug){
	std::cout << "... Processing layer " << volNb << " with " << (*simhitvec).size() << " simhits " << std::endl;
      }
      
      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	HGCSSSimHit lHit = (*simhitvec)[iH];
	unsigned layer = lHit.layer();
	if (layer==volNb+1) {
	  if (ievt/nLayers==0) std::cout << " -- Warning, applying patch to layer number..." << std::endl; 
	  layer = volNb;
	}

	double posx = lHit.get_x();
	double posy = lHit.get_y();
	double posz = lHit.get_z();
	double energy = lHit.energy();
	if (debug>1) {
	  std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		    << " --  position x,y " << posx << "," << posy << std::endl;
	  lHit.Print(std::cout);
	}

	if (energy>0) {
	  hitsOut << static_cast<unsigned>(ievt/nLayers);
	  hitsOut << " " << layer << " " << posx << " " << posy << " " << posz << " " << energy << std::endl;
	}
	
      }//loop on hits

    }//loop on entries
  
    hitsOut.close();

  }//loop on energies

  return 0;

}//main
