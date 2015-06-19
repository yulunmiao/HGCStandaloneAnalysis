#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TChain.h"

#include "HGCSSInfo.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"

#include "HadEnergy.hh"

#include "TRandom3.h"

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

bool getValueFromString(std::vector<std::string>& outputString, std::string& sourceString){
     boost::split( outputString, sourceString, boost::is_any_of(","));
};

bool getValueFromString(std::vector<int>& outputNumber, std::string& sourceString){
     std::vector<std::string> outputString;
     boost::split( outputString, sourceString, boost::is_any_of(","));
     for(unsigned i(0); i < outputString.size(); i++) outputNumber.push_back(atoi(outputString[i].c_str()));
};

void getTotalEnergy(std::vector<std::vector<double>>& Etotal, std::ostringstream& inputsim, std::ostringstream& inputrec){

    TFile * simFile = 0;
    TFile * recFile = 0;
 
    TChain  *lSimTree = new TChain("HGCSSTree");
    TChain  *lRecTree = new TChain("RecoTree");

    if (!testInputFile(inputsim.str(),simFile)) return;
    lSimTree->AddFile(inputsim.str().c_str());
    if (!testInputFile(inputrec.str(),recFile)) return;
    lRecTree->AddFile(inputrec.str().c_str());
     
    HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
    const unsigned versionNumber = info->version();
    HGCSSDetector & myDetector = theDetector();
    myDetector.buildDetector(versionNumber,false,false);
    const unsigned nSec = myDetector.nSections();    
    double recSum[nSec];
    double sumEE(0), sumFH(0), sumBH(0);   
    std::vector<double> EE;
    std::vector<double> EFH;
    std::vector<double> EBH;
    for(unsigned iS(0); iS <nSec; iS++){
       recSum[iS] = 0;
       EE.clear();
       EFH.clear();
       EBH.clear();
    }        
 
    HGCSSEvent * event = 0;
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;
    
    lSimTree->SetBranchAddress("HGCSSEvent",&event);
    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  
    const unsigned nEvts = lSimTree->GetEntries(); 
    
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

      lSimTree->GetEntry(ievt);
      lRecTree->GetEntry(ievt);
      for(unsigned iS(0); iS <nSec; iS++){
        recSum[iS] = 0;
      }        
      sumEE = 0;
      sumFH = 0;
      sumBH = 0;

      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
        const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
        unsigned layer = lHit.layer();
        unsigned sec =  myDetector.getSection(layer);
        double absweight = (*ssvec)[layer].volX0trans()/(*ssvec)[0].volX0trans();
        double energy = lHit.energy();

        recSum[sec] += energy*absweight;
    }
    
    for(unsigned iS(0); iS < nSec; iS++){
      //get total energy for each sub-detector
      DetectorEnum type = myDetector.detType(iS); 
      if(type == DetectorEnum::FECAL || type == DetectorEnum::MECAL || type == DetectorEnum::BECAL)
         sumEE += recSum[iS];
      else if(type == DetectorEnum::FHCAL){
         sumFH += recSum[iS];
      }
      else if(type == DetectorEnum::BHCAL1 || type == DetectorEnum::BHCAL2)
         sumBH += recSum[iS];
    }
    EE.push_back(sumEE);
    EFH.push_back(sumFH);
    EBH.push_back(sumBH);
  }
  Etotal.push_back(EE);
  Etotal.push_back(EFH);
  Etotal.push_back(EBH);
}

void setCalibFactor(double calibSlope, double calibOffset, std::vector<std::pair<std::string, std::string>> fileName1, std::vector<std::pair<std::string, std::string>> fileName2, TCanvas *c){

    c->cd();
    TGraph * slopesummary = new TGraph();
    const unsigned nP = fileName1.size();
    if(fileName1.size() != fileName2.size()){std::cout << "Got different number of EE/FH/BH files when set calibration factor" << std::endl; return;}
    std::vector<std::vector<double>> Etotal;
    double Energy1[nP], Energy2[nP];
    TH1F * p_Energy1=new TH1F("","",1000,0,100000);
    TH1F * p_Energy2=new TH1F("","",1000,0,100000);

    std::ostringstream inputsim, inputrec;
    unsigned iP(0);
    for(std::vector<std::pair<std::string, std::string>>::iterator iF = fileName1.begin() ; iF != fileName1.end(); ++iF){
      
      inputsim.str("");
      inputrec.str("");
      inputsim << iF->first; 
      inputrec << iF->second;
     
      getTotalEnergy(Etotal, inputsim, inputrec);
      for(unsigned iE(0); iE < Etotal[0].size(); iE++)p_Energy1->Fill(Etotal[0][iE]);
      Energy1[iP] = p_Energy1->GetMean();
      std::cout << "Total energy  = " << Energy1[iP] << "(MIPS)" << std::endl;
      p_Energy1->Reset();
      Etotal.clear();
      iP++;
    } 

    iP = 0;
    for(std::vector<std::pair<std::string, std::string>>::iterator iF = fileName2.begin() ; iF != fileName2.end(); ++iF){
      
      inputsim.str("");
      inputrec.str("");
      inputsim << iF->first; 
      inputrec << iF->second;
     
      getTotalEnergy(Etotal, inputsim, inputrec);
      for(unsigned iE(0); iE < Etotal[1].size(); iE++)p_Energy2->Fill(Etotal[1][iE] + Etotal[2][iE]);
      Energy2[iP] = p_Energy2->GetMean();
      std::cout << "Total energy  = " << Energy2[iP] << "(MIPS)" << std::endl;
      p_Energy2->Reset();
      Etotal.clear();
      iP++;
    } 
  
    for(unsigned iP(0); iP < nP; iP++){
      slopesummary->SetPoint(iP,Energy1[iP],Energy2[iP]);
    }
    slopesummary->SetMarkerSize(1);
    slopesummary->SetMarkerStyle(20);
    slopesummary->SetMarkerColor(1);
    slopesummary->Draw("AP");
    slopesummary->Fit("pol1");
    TF1 *fit = (TF1*)slopesummary->GetFunction("pol1");
   // fit->Draw("same");
    calibOffset = fit->GetParameter(0);
    calibSlope = fit->GetParameter(1);
    return;
}

void setCalibFactor(double calibFactor, std::vector<std::pair<std::string, std::string>>& fileName,double FHtoEslope,double FHtoEoffset,double BHtoEslope,double BHtoEoffset,  TCanvas *c){

    c->cd();
    TGraph *slopesummary = new TGraph();
    const unsigned nP = fileName.size();
    TH2F *p_FHvsBH = new TH2F("","",1000,0,10000,1000,0,10000);
    std::vector<std::vector<double>> Etotal;

    std::ostringstream inputsim, inputrec;
    for(std::vector<std::pair<std::string, std::string>>::iterator iF = fileName.begin() ; iF != fileName.end(); ++iF){
      
      inputsim.str("");
      inputrec.str("");
      inputsim << iF->first; 
      inputrec << iF->second;
     
      getTotalEnergy(Etotal, inputsim, inputrec);
      unsigned nevt = Etotal[1].size();
      for(unsigned ievt(0); ievt < nevt; ievt++){
         p_FHvsBH->Fill(Etotal[0][ievt]+(Etotal[1][ievt]-FHtoEoffset)/FHtoEslope, (Etotal[2][ievt]-BHtoEoffset)/BHtoEslope);
      }

      TProfile *prof = p_FHvsBH->ProfileX();
      prof->Fit("pol1");
      TF1 *f1 = (TF1*)prof->GetFunction("pol1");
      double slope = f1->GetParameter(1);
  
      Int_t np=slopesummary->GetN();  
      slopesummary->SetPoint(np,np,slope);
    }
 
    slopesummary->Draw("AP"); 
    double slopeMean = fabs(slopesummary->GetMean(2));
    calibFactor = slopeMean;
} 
         

int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  bool concept;
  //size of signal region to perform Chi2 position fit.
  unsigned pNevts;
  bool doFHEMcalib;
  bool doBHEMcalib;
  bool doFHBHcalib;
  std::string eeSimFiles;
  std::string fhSimFiles;
  std::string bhSimFiles;
  std::string eeRecoFiles;
  std::string fhRecoFiles;
  std::string bhRecoFiles;
  std::string filePath;
  std::string digifilePath;
  unsigned nRuns;
  std::string genEnergy;
  std::string simFileName;
  std::string recoFileName;
  std::string outPath;
  unsigned nSiLayers;
  //0:do just the energies, 1:do fit+energies, 2: do zpos+fit+energies
  unsigned debug;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("concept",        po::value<bool>(&concept)->default_value(true))
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("doFHEMcalib", po::value<bool>(&doFHEMcalib)->default_value(false))
    ("doBHEMcalib", po::value<bool>(&doBHEMcalib)->default_value(false))
    ("doFHBHcalib",  po::value<bool>(&doFHBHcalib)->default_value(false))
    ("eeSimFiles",     po::value<std::string>(&eeSimFiles)->default_value(""))
    ("fhSimFiles",     po::value<std::string>(&fhSimFiles)->default_value(""))
    ("bhSimFiles",     po::value<std::string>(&bhSimFiles)->default_value(""))
    ("eeRecoFiles",     po::value<std::string>(&eeRecoFiles)->default_value(""))
    ("fhRecoFiles",     po::value<std::string>(&fhRecoFiles)->default_value(""))
    ("bhRecoFiles",     po::value<std::string>(&bhRecoFiles)->default_value(""))
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("genEnergy",    po::value<std::string>(&genEnergy)->required()) 
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("outPath,o",      po::value<std::string>(&outPath)->required())
    ("nSiLayers",      po::value<unsigned>(&nSiLayers)->default_value(2))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ;

  // ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
//  bool isCalibed = true;
  bool isCalibed = false;
  //////////////////////////////////////////////////////////
  //// Hardcoded factor ////////////////////////////////////
  //////////////////////////////////////////////////////////
  double FHtoEslope = 11.085; 
  double FHtoEoffset = -4.34;
  double BHtoEslope = 4.079;
  //double BHtoEoffset = 293.55;
  double BHtoEoffset = 0;
  double ECALslope = 116.897;
  double ECALoffset = -101.775;

  double FHtoBHslope  = 0.2; 
 
  double HcalPionOffset = 0.; 
  double HcalPionCalib = 0.92;
////////////////////////////////////
//
  if(doFHEMcalib){
    if(eeSimFiles == ""|| eeRecoFiles == "" || fhSimFiles == "" || fhRecoFiles == "") {std::cout << "EE/FH files cannot be empty for FHvsEE calibration" << std::endl;}
    else{
      std::vector<std::string> eeSimFileVec, eeRecoFileVec, fhSimFileVec, fhRecoFileVec;
      getValueFromString( eeSimFileVec, eeSimFiles);
      getValueFromString( eeRecoFileVec,eeRecoFiles);
      getValueFromString( fhSimFileVec,fhSimFiles);
      getValueFromString( fhRecoFileVec,fhRecoFiles);
      std::vector<std::pair<std::string, std::string>> eeFileName, fhFileName;
      for(unsigned iF(0); iF<eeSimFileVec.size(); iF++){
         eeFileName.push_back(make_pair(eeSimFileVec[iF], eeRecoFileVec[iF]));
         fhFileName.push_back(make_pair(fhSimFileVec[iF], fhRecoFileVec[iF]));
      }
      TCanvas *c1 = new TCanvas("FHcalibtoEM","FHcalibtoEM",800,600);
      setCalibFactor(FHtoEslope,FHtoEoffset,eeFileName, fhFileName, c1); 
      c1->SaveAs("FHcalibtoEM.png");
    }
 }
  if(doBHEMcalib){
    if(eeSimFiles == ""|| eeRecoFiles == "" || bhSimFiles == "" || bhRecoFiles == "") {std::cout << "EE/FH files cannot be empty for FHvsEE calibration" << std::endl;}
    else{
      std::vector<std::string> eeSimFileVec, eeRecoFileVec, bhSimFileVec, bhRecoFileVec;
      getValueFromString( eeSimFileVec,eeSimFiles);
      getValueFromString( eeRecoFileVec,eeRecoFiles);
      getValueFromString( bhSimFileVec,bhSimFiles);
      getValueFromString( bhRecoFileVec,bhRecoFiles);
      std::vector<std::pair<std::string, std::string>> eeFileName, bhFileName;
      for(unsigned iF(0); iF<eeSimFileVec.size(); iF++){
         eeFileName.push_back(make_pair(eeSimFileVec[iF], eeRecoFileVec[iF]));
         bhFileName.push_back(make_pair(bhSimFileVec[iF], bhRecoFileVec[iF]));
      }
      TCanvas *c2 = new TCanvas("BHcalibtoEM","BHcalibtoEM",800,600);
      setCalibFactor(BHtoEslope,BHtoEoffset,eeFileName, bhFileName, c2); 
      c2->SaveAs("BHcalibtoEM.png");
    }
 }
  if(doFHBHcalib){
    if(fhSimFiles == ""|| fhRecoFiles == "") {std::cout << "FH files cannot be empty for FHvsEE calibration" << std::endl;}
    else{
      std::vector<std::string> fhSimFileVec, fhRecoFileVec;
      getValueFromString( fhSimFileVec,fhSimFiles);
      getValueFromString( fhRecoFileVec,fhRecoFiles);
      std::vector<std::pair<std::string, std::string>> fhFileName;
      for(unsigned iF(0); iF<fhSimFileVec.size(); iF++){
         fhFileName.push_back(make_pair(fhSimFileVec[iF], fhRecoFileVec[iF]));
      }
      TCanvas *c3 = new TCanvas("BHcalibtoFH","BHcalibtoFH",800,600);
      setCalibFactor(FHtoBHslope, fhFileName, FHtoEslope,FHtoEoffset,BHtoEslope,BHtoEoffset,c3); 
    }
 }
      

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(1234);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  std::string inFilePath = filePath+simFileName;

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  std::vector<int> genEn;
  getValueFromString(genEn, genEnergy);
  const unsigned nGenEn=genEn.size();
  double eta = 2.00;

  /////////////////////////////////////////////////////////////
  //input file format
  ///////////////////////////////////////////////////////////////

  std::size_t begin = simFileName.find("_et")+3;
  std::size_t middle = simFileName.find("_eta");
  std::size_t end = simFileName.find_last_of(".root")+1;
  std::size_t run(0);
  if(nRuns>0)run = simFileName.find("_run");
  std::string simHeader = simFileName.substr(0,begin);
  std::string simAppend("");
  if(nRuns>0)simAppend = simFileName.substr(middle,run-middle);
  else simAppend = simFileName.substr(middle,end-middle);

  std::size_t begin_rec = recoFileName.find("_et")+3;
  std::size_t middle_rec = recoFileName.find("_eta");
  std::size_t end_rec = recoFileName.find_last_of(".root")+1;
  std::size_t run_rec(0);
  if(nRuns>0)run_rec = recoFileName.find("_run");
  std::string recHeader = recoFileName.substr(0,begin_rec);
  std::string recAppend;
  if(nRuns>0)recAppend = recoFileName.substr(middle_rec,run_rec-middle_rec);
  else recAppend = recoFileName.substr(middle_rec,end_rec-middle_rec);


  std::vector<double> EE[nGenEn];
  std::vector<double> EFHCAL[nGenEn];
  std::vector<double> EBHCAL[nGenEn];
  std::vector<double> GlobalC[nGenEn];

  for(unsigned iGen(0); iGen < nGenEn; iGen++){
 
     TFile * simFile = 0;
     TFile * recFile = 0;
  
     TChain  *lSimTree = new TChain("HGCSSTree");
     TChain  *lRecTree = new TChain("RecoTree");
   
     std::ostringstream inputsim, inputrec;
     if (nRuns == 0){
       inputsim << filePath << "/" << simHeader << genEn[iGen] << simAppend;
       inputrec << filePath << "/" << recHeader << genEn[iGen] << recAppend;;
       if (!testInputFile(inputsim.str(),simFile)) return 1;
       lSimTree->AddFile(inputsim.str().c_str());
       if (!testInputFile(inputrec.str(),recFile)) return 1;
       lRecTree->AddFile(inputrec.str().c_str());
     }
     else {
       for (unsigned i(0);i<nRuns;++i){
         inputsim.str("");
         inputsim << filePath << "/" << simHeader << genEn[iGen] << simAppend << "_run" << i << ".root";
         inputrec.str("");
         inputrec << filePath << "/" << recHeader << genEn[iGen] << recAppend << "_run" << i << ".root";
         if (!testInputFile(inputsim.str(),simFile)) return 1;
         lSimTree->AddFile(inputsim.str().c_str());
         if (!testInputFile(inputrec.str(),recFile)) return 1;
         lRecTree->AddFile(inputrec.str().c_str());
       }
     }


     /////////////////////////////////////////////////////////////
     //output
     ///////////////////////////////////////////////////////////////
     std::ostringstream outFile;
     outFile.str("");
     outFile << outPath << "pion_reso_et" << genEn[iGen] << ".root";
     TFile *outputFile = TFile::Open(outFile.str().c_str(),"RECREATE");

     if (!outputFile) {
       std::cout << " -- Error, output file " << outFile.str() << " cannot be opened. Please create output directory. Exiting..." << std::endl;
     return 1;
     }
     else {
       std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
     }
   
     ///////////////////////////////
     // Info    ////////////////////
     ///////////////////////////////
   
     HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
     const double cellSize = info->cellSize();
     const unsigned versionNumber = info->version();
     const unsigned model = info->model();

     bool isCaliceHcal = versionNumber==23;   
     std::cout << " -- Version number is : " << versionNumber 
   	    << ", model = " << model
   	    << ", cellSize = " << cellSize
   	    << std::endl;

     //initialise detector
     HGCSSDetector & myDetector = theDetector();
 
     myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

     const unsigned nLayers = myDetector.nLayers();
     const unsigned nSections = myDetector.nSections();

     std::cout << " -- N layers = " << nLayers << std::endl
      	       << " -- N sections = " << nSections << std::endl;
     
     HadEnergy myhadReso(myDetector, lSimTree, lRecTree, outputFile, pNevts);
     myhadReso.addLimMIP(5.5);
     myhadReso.setFHtoE(FHtoEslope, FHtoEoffset);
     myhadReso.setBHtoE(BHtoEslope, BHtoEoffset);
     myhadReso.setEEcalib(ECALslope, ECALoffset);
     
 
     myhadReso.fillEnergies(); 

     outputFile->Close();
   }// loop over GenEn 
  
  return 0;
  
  
}//main
