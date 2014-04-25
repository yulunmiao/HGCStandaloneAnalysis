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

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"

#include "TRandom3.h"

bool fillHistoConditions(const bool selectEarlyDecay,
			 const double E_FHCAL,
			 const double Etotal){
  if (selectEarlyDecay) {
    if (E_FHCAL < 0.95*Etotal) return false; 
  }
  return true;
}

void findMeanEnergy(TTree *lTree,
		    const bool isG4Tree,
		    std::vector<HGCSSSimHit>* & simhitvec,
		    float & volNb,
		    const unsigned nLayers,
		    const unsigned nEvts,
		    const bool selectEarlyDecay,
		    HGCSSCalibration & mycalib,
		    TH1F *p_Etotal){

  double Etotal[3] = {0,0,0};
  double Etot[nLayers];
  for (unsigned iL(0);iL<nLayers;++iL){
    Etot[iL] = 0;
  }

  for (unsigned ievt(0); ievt<(isG4Tree?nEvts*nLayers:nEvts); ++ievt){//loop on entries
    if (ievt%1000 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
    lTree->GetEntry(ievt);

    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];
      unsigned layer = lHit.layer();
      double energy = lHit.energy();
      
      double weightedE = energy*mycalib.MeVToMip(layer);
      
      Etot[layer] += energy;
      //if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << Etot[layer];
      mycalib.incrementEnergy(layer,weightedE,Etotal[0],Etotal[1],Etotal[2]);
    }//loop on hits

    if (!isG4Tree || (isG4Tree && (volNb == nLayers-1))){//fill histos
      //if (debug) std::cout << " -- Filling histograms..." << std::endl;

      double EReco = mycalib.recoEnergyUncor(Etotal[0],Etotal[1],Etotal[2]);
      double EReco_FHCAL = mycalib.recoEnergyUncor(Etotal[0],Etotal[1],Etotal[2],true);
      bool doFill = fillHistoConditions(selectEarlyDecay,EReco_FHCAL,EReco);

      if (doFill) {
	p_Etotal->Fill(EReco);
      }

      for (unsigned iD(0); iD<3; ++iD){
	Etotal[iD] = 0;
      }
      for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
	Etot[iL] = 0;
      }

    }//save histos

  }//loop on entries
  
}//findMeanEnergy


// void fillHitSpectrums(TTree *lTree,
// 		      const bool isG4Tree,
// 		      std::vector<HGCSSSimHit> * & simhitvec,
// 		      float & volNb,
// 		      const unsigned nLayers,
// 		      const unsigned nEvts,
// 		      const bool selectEarlyDecay,
// 		      HGCSSCalibration & mycalib,
// 		      const double & meanRecoE,
// 		      const double & rmsRecoE,
// 		      TH1F *p_hitSpectrum_lowTail,
// 		      TH1F *p_hitSpectrum_highTail){

//   double Etotal[3] = {0,0,0};
//   double Etot[nLayers];
//   for (unsigned iL(0);iL<nLayers;++iL){
//     Etot[iL] = 0;
//   }
//   TH2F * EmipHits = new TH2F("EmipHitsSpectra",";x(mm);y(mm)",17,-255,255,17,-255,255);

//   for (unsigned ievt(0); ievt<(isG4Tree?nEvts*nLayers:nEvts); ++ievt){//loop on entries
//     if (ievt%1000 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
//     lTree->GetEntry(ievt);

//     for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
//       HGCSSSimHit lHit = (*simhitvec)[iH];
//       unsigned layer = lHit.layer();
//       double energy = lHit.energy();

//       double posx = lHit.get_x();
//       double posy = lHit.get_y();

//       double weightedE = energy*mycalib.MeVToMip(layer);
      
//       Etot[layer] += energy;
//       //if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << Etot[layer];
//       mycalib.incrementEnergy(layer,weightedE,Etotal[0],Etotal[1],Etotal[2]);

//       //calculate mean Emip
//       EmipHits->Fill(posx,posy,weightedE);


//     }//loop on hits

//     if (!isG4Tree || (isG4Tree && (volNb == nLayers-1))){//fill histos
//       //if (debug) std::cout << " -- Filling histograms..." << std::endl;
//       double EReco = mycalib.recoEnergyUncor(Etotal[0],Etotal[1],Etotal[2]);
//       double EReco_FHCAL = mycalib.recoEnergyUncor(Etotal[0],Etotal[1],Etotal[2],true);
//       bool doFill = fillHistoConditions(selectEarlyDecay,EReco_FHCAL,EReco);

//       if (doFill) {
// 	bool lowTail = EReco < (meanRecoE-rmsRecoE);
// 	bool highTail = EReco > (meanRecoE+rmsRecoE);
// 	for (int ix(1); ix<EmipHits->GetNbinsX()+1; ++ix){
// 	  for (int iy(1); iy<EmipHits->GetNbinsY()+1; ++iy){
// 	    double eTmp = EmipHits->GetBinContent(ix,iy);
// 	    if (eTmp<0.5) continue;
// 	    if (lowTail) p_hitSpectrum_lowTail->Fill(eTmp);
// 	    if (highTail) p_hitSpectrum_highTail->Fill(eTmp);
// 	  }
// 	}
//       }

//       for (unsigned iD(0); iD<3; ++iD){
// 	Etotal[iD] = 0;
//       }
//       for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
// 	Etot[iL] = 0;
//       }
//       EmipHits->Reset();
      
//     }//save histos

//   }//loop on entries
  
// }//fillHitSpectrums


int main(int argc, char** argv){//main  

  if (argc < 4) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <full path to input file>"
	      << " <file name (PFcal.root, or DigiPFcal.root)>" 
	      << " <optional: list of energies: 5,10,25>"
	      << " <optional: debug (default=0)>"
	      << " <optional: doMeanRms (default=0)>"
      //<< " <optional: doHitSpectra (default=0)>"
	      << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  bool calibrate = true;//do MIP to GeV for total E.
  bool concept = false;

  bool selectEarlyDecay = true;//>95% of E in HCAL

  const unsigned nLimits = 15;
  const double pElim[nLimits] = {10,5,6,7,8,9,9.5,10.5,11,11.5,12,13,14,15,20};

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(1);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  //for basic control plots
  TCanvas *mycE = new TCanvas("mycE","mycE",1500,1000);
  gStyle->SetOptStat("eMRuoi");

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string fileName = argv[3];//"PFcal.root";
  unsigned debug = 0;
  if (argc >5) debug = atoi(argv[5]);
  unsigned doMeanRms = 0;
  //unsigned doHitSpectra = 0;
  if (argc >6) doMeanRms = atoi(argv[6]);
  //if (argc >7) doHitSpectra = atoi(argv[7]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Input file name: " << fileName << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;
  if (doMeanRms) std::cout << " -- Making data file with mean-RMS" << std::endl;
  //if (doHitSpectra) std::cout << " -- Making data file with hit spectra mean" << std::endl;


  unsigned tmpEn[]={30};
  unsigned tmpSize = sizeof(tmpEn)/sizeof(unsigned);
  std::vector<unsigned> genEn;

  std::vector<std::string> enVec;
  std::string genEnStr;
  if (argc >4) {
    genEnStr = argv[4];
    boost::split( enVec, genEnStr, boost::is_any_of(","));
    for (unsigned iE(0); iE<enVec.size();++iE){
      unsigned tmpE = 0;
      std::istringstream(enVec[iE])>>tmpE;
      genEn.push_back(tmpE);
    }
  } else {
    for (unsigned iE(0); iE<tmpSize;++iE){
      genEn.push_back(tmpEn[iE]);
    }
  }
  const unsigned nGenEn=genEn.size();
  
  std::cout << " -- Setup to run on " << nGenEn << " energy points." << std::endl;


  TString plotDir = "PLOTS/";
  std::vector<std::string> pathVec;
  bool pEOS = false;
  if (filePath.find("eos") != filePath.npos) {
    pEOS = true;
  }
   
  //if (filePath.find("version_23") != filePath.npos) pEOS = false;

  if (pEOS) boost::split( pathVec, filePath, boost::is_any_of("/_"));
  else  boost::split( pathVec, filePath, boost::is_any_of("/"));
  std::string pVersion("");
  std::string pScenario("");
  std::string pParticle("");
  std::string pEta("");

  for (unsigned i(0);i<pathVec.size();++i){
    std::string lTmp = pathVec[pEOS? pathVec.size()-i-1 : i];
    if (lTmp.find("version") != lTmp.npos) {
      pVersion = lTmp;//+"_"+pathVec[pathVec.size()-i];
    }
     else if (lTmp.find("scenario") != lTmp.npos) {
      pScenario = lTmp;
    }
    else if (pEOS && pVersion.size()!=0 && lTmp.find("HGcal")!=lTmp.npos){
      continue;
    }
    else if (pVersion.size()!=0 && pParticle.size()==0 && lTmp.find("Geant")==lTmp.npos)  {
      pParticle = lTmp;
      if (pEOS) break;
    }
    else if (lTmp.find("eta") != lTmp.npos) {
      pEta = lTmp;
      if (!pEOS) break;
    }
    else if (pVersion.size()!=0) {
      break;
    }
  }

  plotDir += pVersion;
  plotDir += "/";
  plotDir += pScenario;
  plotDir += "/";
  plotDir += pParticle;
  plotDir += "/";
  plotDir += pEta;
  plotDir += "/";

  if (concept) plotDir += "concept/";
  else plotDir += "twiceSampling/";
  if (calibrate) plotDir += "GeVCal/";
  //plotDir += "simpleWeights/";
  if (selectEarlyDecay) plotDir += "EarlyDecay/";


  std::cout << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----- Hardcoded configuration is : " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " -- calibrate = " <<calibrate << std::endl
	    << " -- concept = " <<concept << std::endl;

  //initialise calibration class

  HGCSSCalibration mycalib(filePath,concept,calibrate);

  const unsigned nLayers = mycalib.nLayers();//64;//Calice 54;//Scint 9;//HCAL 33;//All 64
  const unsigned nEcalLayers = mycalib.nEcalLayers();//31;
  const unsigned nFHcalLayers = mycalib.nFHcalLayers();//concept 24;//calice 47

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- nEcalLayers = " << nEcalLayers 
	    << ", mip weights = " << mycalib.mipWeight(0) << " " << mycalib.mipWeight(1) << " " << mycalib.mipWeight(2)
	    << ", GeV weights = " << mycalib.gevWeight(0) << " offset " << mycalib.gevOffset(0)
	    << std::endl
	    << " -- nFHcalLayers  = " << nFHcalLayers 
	    << ", mip weights = " << mycalib.mipWeight(3) << " " << mycalib.mipWeight(4) << " " << mycalib.mipWeight(5) 
	    << ", GeV weights = " << mycalib.gevWeight(1) << " " << mycalib.gevWeight(2) << " offsets: " << mycalib.gevOffset(1) << " " <<   mycalib.gevOffset(2)
	    << std::endl
	    << " -- conversions: HcalToEcalConv = " <<mycalib.HcalToEcalConv() << " BHcalToFHcalConv = " << mycalib.BHcalToFHcalConv() << std::endl
	    << " -- Elim for Cglobal = " << pElim[0] << "-" << pElim[nLimits-1] << " Mips" << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " -----------------------------------" << std::endl;

  std::cout << " -- Output file directory is : " << plotDir << std::endl;


  TFile *outputFile = 0;
  if (!doMeanRms){ 
    outputFile = TFile::Open(plotDir+"/CalibHcalHistos_test.root","RECREATE");
    
    if (!outputFile) {
      std::cout << " -- Error, output file " << outputFile->GetName() << " cannot be opened. Please create output directory : " << plotDir << "  Exiting..." << std::endl;
      return 1;
    }
    std::cout << " -- Created output file " << outputFile->GetName() << std::endl;
  }

  TH1F *p_Etotal[nGenEn];
  TH1F *p_EshowerCor[nGenEn];
  TH1F *p_hitSpectrum_lowTail[nGenEn];
  TH1F *p_hitSpectrum_highTail[nGenEn];
  TH1F *p_meanHitSpectrum_lowTail[nGenEn];
  TH1F *p_meanHitSpectrum_highTail[nGenEn];
  TH1F *p_Cglobal[nGenEn][nLimits];
  TH1F *p_Eshower[nGenEn][nLimits];
  TH2F *p_EvsCglobal[nGenEn][nLimits];

  std::ostringstream lName;

  bool isG4Tree = true;

  std::ofstream meanRmsOut;
  std::ifstream meanRmsIn;

  std::map<unsigned,double> meanRecoE;
  std::map<unsigned,double> rmsRecoE;
  //std::map<unsigned,double> EmipMean;

  if (doMeanRms) {
    meanRmsOut.open(plotDir+"/MeanRmsHcal.dat");
  }
  else {
    meanRmsIn.open(plotDir+"/MeanRmsHcal.dat");
    if(!meanRmsIn.is_open()){
      std::cerr << "Unable to open file for writting meanRMS for HCAL. " << std::endl;
      return 1;
    }
    while(1){
      unsigned genE = 0;
      double mean = 0;
      double rms = 0;
      meanRmsIn>>genE>>mean>>rms;
      meanRecoE[genE] = mean;
      rmsRecoE[genE] = rms;
      std::cout << genE << " " << mean << " " << rms << std::endl;
      if(meanRmsIn.eof()){
	break; 
      }
    }
    std::cout << " -- Number of energy points found: " << meanRecoE.size() << std::endl;
    meanRmsIn.close();
  }


  // std::ofstream hitSpectraOut;
  // std::ifstream hitSpectraIn;
  // if (doHitSpectra) hitSpectraOut.open(plotDir+"/HitSpectraHcal.dat");
  // else {
  //   hitSpectraIn.open(plotDir+"/HitSpectraHcal.dat");

  //   if(!hitSpectraIn.is_open()){
  //     std::cerr << "Unable to open file for writting hitSpectra for HCAL. " << std::endl;
  //     return 1;
  //   }
  //   while(1){
  //     unsigned genE = 0;
  //     double mean = 0;
  //     hitSpectraIn>>genE>>mean;
  //     EmipMean[genE] = mean;
  //     std::cout << genE << " " << mean << std::endl;
  //     if(hitSpectraIn.eof()){
  // 	break; 
  //     }
  //   }
  //   std::cout << " -- Number of energy points found: " << EmipMean.size() << std::endl;

  // }

  for (unsigned iE(0); iE<nGenEn; ++iE){

    std::cout << "- Processing energy : " << genEn[iE] << std::endl;

    std::ostringstream input;
    input << filePath ;
    //if (filePath.find("pi+")!=filePath.npos) input << "/pt_" ;
    //HGcal_version_8_e50.root
    if (pEOS) {
      input << "_e" ;
      input << genEn[iE] << fileName;
    }
    else {
      input << "/e_" ;
      input << genEn[iE];
      input << "/" << fileName;
    }
    TFile *inputFile = TFile::Open(input.str().c_str());

    if (!inputFile) {
      std::cout << " -- Error, input file " << input.str() << " cannot be opened. Trying next one..." << std::endl;
      continue;
    }

    TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
    if (!lTree){
      std::cout << " -- Error, tree HGCSSTree cannot be opened. Trying RecoTree instead." << std::endl;
      isG4Tree = false;
      lTree = (TTree*)inputFile->Get("RecoTree");
      if (!lTree){
	std::cout << " -- Error, tree RecoTree cannot be opened either. Exiting..." << std::endl;
	return 1;
      }
    }


    unsigned event = 0;
    float volNb = 0;
    float volX0 = 0;
    //float volX0ref = 0;
    //float volLambda = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;

    if (isG4Tree){
      //lTree->SetBranchAddress("event",&event);
      lTree->SetBranchAddress("volNb",&volNb);
      lTree->SetBranchAddress("volX0trans",&volX0);
      //lTree->SetBranchAddress("volLambda",&volLambda);
      lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    }
    else {
      lTree->SetBranchAddress("event",&event);
      lTree->SetBranchAddress("volX0",&volX0);
      //lTree->SetBranchAddress("volLambda",&volLambda);
      lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
      lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
    }

    const unsigned nEvts = 
      isG4Tree ? 
      ((pNevts > lTree->GetEntries()/nLayers || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()/nLayers) : pNevts) :
      ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;



    std::cout << "- Processing = " << nEvts  << " events out of " ;
    if (isG4Tree) std::cout << lTree->GetEntries()/nLayers << std::endl;
    else std::cout << lTree->GetEntries() << std::endl;
    
    lName.str("");
    lName << "p_Etotal_" << genEn[iE];
    p_Etotal[iE] = new TH1F(lName.str().c_str(),";G4 Etotal (GeV)",1000,0,2*genEn[iE]);
    p_Etotal[iE]->StatOverflows();


    bool firstEvent = true;

    //first training loop
    if (doMeanRms){
      findMeanEnergy(lTree,isG4Tree,
		     simhitvec,volNb,
		     nLayers,nEvts,
		     selectEarlyDecay,
		     mycalib,
		     p_Etotal[iE]);
      
      meanRecoE[genEn[iE]] = p_Etotal[iE]->GetMean();
      rmsRecoE[genEn[iE]] = p_Etotal[iE]->GetRMS();
      
      meanRmsOut <<  genEn[iE] << " " 
		 << meanRecoE[genEn[iE]] << " " 
		 << rmsRecoE[genEn[iE]] 
		 << std::endl;

      std::cout << " E=" << genEn[iE] << " Entries: " << p_Etotal[iE]->GetEntries() << ", <Ereco> " << meanRecoE[genEn[iE]] << " RMS " << rmsRecoE[genEn[iE]] << " GeV" << std::endl;

      mycE->cd();
      mycE->SetLogy(0);
      p_Etotal[iE]->GetXaxis()->SetRangeUser(meanRecoE[genEn[iE]]-5*rmsRecoE[genEn[iE]],
					     meanRecoE[genEn[iE]]+5*rmsRecoE[genEn[iE]]);
      std::ostringstream saveName;
      saveName.str("");
      saveName << plotDir << "/ControlERecoTotal_" << genEn[iE] << "GeV";
      mycE->Update();
      mycE->Print((saveName.str()+".png").c_str());
      mycE->Print((saveName.str()+".pdf").c_str());

      //if (!doHitSpectra) 
      continue;
    }
    else {
      std::cout << " From input data file: E=" << genEn[iE] << ": <Ereco> " << meanRecoE[genEn[iE]] << " RMS " << rmsRecoE[genEn[iE]] << " GeV" << std::endl;
    }


    // if (doHitSpectra){
    //   fillHitSpectrums(lTree,isG4Tree,
    // 		       simhitvec,volNb,
    // 		       nLayers,nEvts,
    // 		       selectEarlyDecay,
    // 		       mycalib,
    // 		       meanRecoE[genEn[iE]],rmsRecoE[genEn[iE]],
    // 		       p_hitSpectrum_lowTail[iE],
    // 		       p_hitSpectrum_highTail[iE]);
      

    //   p_hitSpectrum_lowTail[iE]->Scale(1./p_hitSpectrum_lowTail[iE]->GetEntries());
    //   p_hitSpectrum_highTail[iE]->Scale(1./p_hitSpectrum_highTail[iE]->GetEntries());

    //   //get mean only in the range 0-15...???
    //   p_hitSpectrum_lowTail[iE]->GetXaxis()->SetRangeUser(0,15);
    //   p_hitSpectrum_highTail[iE]->GetXaxis()->SetRangeUser(0,15);
    //   double meanHitLowTail = p_hitSpectrum_lowTail[iE]->GetMean();
    //   double meanHitHighTail = p_hitSpectrum_highTail[iE]->GetMean();
    //   EmipMean[genEn[iE]] = (meanHitLowTail+meanHitHighTail)/2.;
      
    //   hitSpectraOut << genEn[iE] << " " 
    // 		    << EmipMean[genEn[iE]]
    // 		    << std::endl;

    //   std::cout << " E=" << genEn[iE] << ": <Ehit> lowTail " << meanHitLowTail << " hightail " << meanHitHighTail << " avg " << EmipMean[genEn[iE]] << " MIPs." << std::endl;
      
    //   mycE->cd();
    //   mycE->SetLogy(1);
    //   p_hitSpectrum_lowTail[iE]->GetXaxis()->SetRangeUser(0,200);
    //   p_hitSpectrum_highTail[iE]->GetXaxis()->SetRangeUser(0,200);
    //   p_hitSpectrum_lowTail[iE]->SetLineColor(4);
    //   p_hitSpectrum_highTail[iE]->SetLineColor(2);
    //   p_hitSpectrum_lowTail[iE]->Draw();
    //   p_hitSpectrum_highTail[iE]->Draw("same");

    //   std::ostringstream saveName;
    //   saveName.str("");
    //   saveName << plotDir << "/ControlHitSpectra_log_" << genEn[iE] << "GeV";
    //   mycE->Update();
    //   mycE->Print((saveName.str()+".png").c_str());
    //   mycE->Print((saveName.str()+".pdf").c_str());

    //   mycE->SetLogy(0);
    //   //p_hitSpectrum_lowTail[iE]->Rebin(2);
    //   //p_hitSpectrum_highTail[iE]->Rebin(2);
    //   p_hitSpectrum_lowTail[iE]->GetXaxis()->SetRangeUser(0,15);
    //   p_hitSpectrum_highTail[iE]->GetXaxis()->SetRangeUser(0,15);
    //   p_hitSpectrum_lowTail[iE]->SetLineColor(4);
    //   p_hitSpectrum_highTail[iE]->SetLineColor(2);
    //   p_hitSpectrum_lowTail[iE]->Draw();
    //   p_hitSpectrum_highTail[iE]->Draw("same");

    //   saveName.str("");
    //   saveName << plotDir << "/ControlHitSpectra_" << genEn[iE] << "GeV";
    //   mycE->Update();
    //   mycE->Print((saveName.str()+".png").c_str());
    //   mycE->Print((saveName.str()+".pdf").c_str());

    //   TH1F *htmp = (TH1F*)p_hitSpectrum_highTail[iE]->Clone();
    //   htmp->Divide(p_hitSpectrum_highTail[iE],p_hitSpectrum_lowTail[iE]);
    //   htmp->GetYaxis()->SetTitle("Ratio highTail/LowTail");
    //   htmp->GetXaxis()->SetRangeUser(0,15);
    //   mycE->SetGridx(1);
    //   mycE->SetGridy(1);
    //   htmp->Draw("");
    //   saveName.str("");
    //   saveName << plotDir << "/ControlHitSpectraRatio_" << genEn[iE] << "GeV";
    //   mycE->Update();
    //   mycE->Print((saveName.str()+".png").c_str());
    //   mycE->Print((saveName.str()+".pdf").c_str());

    //   mycE->SetGridx(0);
    //   mycE->SetGridy(0);


    //   continue;
    // }
    // else {
    //   std::cout << " --Using from input data file: E=" << genEn[iE] << ": <Ehit> avg " << EmipMean[genEn[iE]] << " MIPs." << std::endl;
    // }

    //Initialise histos
    //necessary to have overflows ?
    gStyle->SetOptStat(1111111);
    double Etotal[3] = {0,0,0};
    double Etot[nLayers];
    for (unsigned iL(0);iL<nLayers;++iL){
      Etot[iL] = 0;
    }

    TH2F * EmipHits = new TH2F("EmipHits",";x(mm);y(mm)",17,-255,255,17,-255,255);
    //TH2F * EmipHits = new TH2F("EmipHits",";x(mm);y(mm)",16,-240,240,16,-240,240);
    unsigned nTotal = 0;
    unsigned nTotalSignal = 0;

    lName.str("");
    lName << "p_hitSpectrum_lowTail_" << genEn[iE];
    p_hitSpectrum_lowTail[iE] =  new TH1F(lName.str().c_str(),";E^{hit}_{HCAL} (MIPs)",1000,0,200);
    p_hitSpectrum_lowTail[iE]->Sumw2();

    lName.str("");
    lName << "p_hitSpectrum_highTail_" << genEn[iE];
    p_hitSpectrum_highTail[iE] =  new TH1F(lName.str().c_str(),";E^{hit}_{HCAL} (MIPs)",1000,0,200);
    p_hitSpectrum_highTail[iE]->Sumw2();

    lName.str("");
    lName << "p_meanHitSpectrum_lowTail_" << genEn[iE];
    p_meanHitSpectrum_lowTail[iE] =  new TH1F(lName.str().c_str(),";<E_{hit}^{HCAL}> (MIPs)",1000,0,50);
    p_meanHitSpectrum_lowTail[iE]->Sumw2();

    lName.str("");
    lName << "p_meanHitSpectrum_highTail_" << genEn[iE];
    p_meanHitSpectrum_highTail[iE] =  new TH1F(lName.str().c_str(),";<E_{hit}^{HCAL}> (MIPs)",1000,0,50);
    p_meanHitSpectrum_highTail[iE]->Sumw2();

    lName.str("");
    lName << "p_EshowerCor_" << genEn[iE];
    p_EshowerCor[iE] = new TH1F(lName.str().c_str(),";E_{shower} (GeV)",1000,0,2*genEn[iE]);

    for (unsigned iLim(0); iLim<nLimits;++iLim){
      lName.str("");
      lName << "p_Cglobal_" << genEn[iE] << "_" << iLim;
      p_Cglobal[iE][iLim] =  new TH1F(lName.str().c_str(),";C_{global}",1000,0,2);
      
      lName.str("");
      lName << "p_Eshower_" << genEn[iE] << "_" << iLim;
      p_Eshower[iE][iLim] = new TH1F(lName.str().c_str(),";E_{shower} (GeV)",1000,0,2*genEn[iE]);

      lName.str("");
      lName << "p_EvsCglobal_" << genEn[iE] << "_" << iLim;
      p_EvsCglobal[iE][iLim] =  new TH2F(lName.str().c_str(),";C_{global};Etot (GeV)",
					 1000,0,2,
					 1000,0,2*genEn[iE]);

    }

    for (unsigned ievt(0); ievt<(isG4Tree?nEvts*nLayers:nEvts); ++ievt){//loop on entries
      if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
      else if (ievt%1000 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
      //unsigned lEvt = isG4Tree? ievt/nLayers : ievt;

      lTree->GetEntry(ievt);

      //if (volNb==0) volX0ref = volX0;
      //else if (volNb == nEcalLayers) volX0ref = volX0;

      if (debug){
	if (isG4Tree) {
	  std::cout << "... Processing layer " << volNb << " with " << (*simhitvec).size() << " simhits, X0 = " << volX0 << " weight = " << mycalib.MeVToMip(static_cast<unsigned>(volNb)) << std::endl;
	}
	else {
	  std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
	}
      }
      

      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	HGCSSSimHit lHit = (*simhitvec)[iH];
	unsigned layer = lHit.layer();
	if (isG4Tree && layer==volNb+1) {
	  if (firstEvent) std::cout << " -- Warning, applying patch to layer number..." << std::endl; 
	  layer = volNb;
	}

	double posx = lHit.get_x();
	double posy = lHit.get_y();
	//double posz = lHit.get_z();

	double energy = lHit.energy();
	if (debug>1) {
	  std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		    << " --  position x,y " << posx << "," << posy << std::endl;
	  lHit.Print(std::cout);
	}

	double weightedE = energy*mycalib.MeVToMip(layer);

	Etot[layer] += energy;
	if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << Etot[layer];
	mycalib.incrementEnergy(layer,weightedE,Etotal[0],Etotal[1],Etotal[2]);
	//calculate mean Emip just for HCAL layers
	if (layer<nFHcalLayers) EmipHits->Fill(posx,posy,weightedE);

      }//loop on hits

      nTotal += (*simhitvec).size();
      nTotalSignal += (*simhitvec).size();

      if (debug)  std::cout << std::endl;

      if (!isG4Tree || (isG4Tree && (volNb == nLayers-1))){//fill histos
	//if (debug) std::cout << " -- Filling histograms..." << std::endl;

	double EReco = mycalib.recoEnergyUncor(Etotal[0],Etotal[1],Etotal[2]);
	double EReco_FHCAL = mycalib.recoEnergyUncor(Etotal[0],Etotal[1],Etotal[2],true);
	bool doFill = fillHistoConditions(selectEarlyDecay,EReco_FHCAL,EReco);

	if (doFill) {
	  bool lowTail = EReco < (meanRecoE[genEn[iE]]-rmsRecoE[genEn[iE]]);
	  bool highTail = EReco > (meanRecoE[genEn[iE]]+rmsRecoE[genEn[iE]]);
	  //fill Cglobal histos
	  unsigned nHitsCountNum[nLimits];
	  unsigned nHitsCountDen = 0;

	  //get mean hit energy
	  double EmipMean = 0;
	  unsigned nHits = 0;
	  for (int ix(1); ix<EmipHits->GetNbinsX()+1; ++ix){
	    for (int iy(1); iy<EmipHits->GetNbinsY()+1; ++iy){
	      double lNoise = lRndm->Gaus(0,0.12);//0.12mip noise
	      EmipHits->SetBinContent(ix,iy,EmipHits->GetBinContent(ix,iy)+lNoise);
	      double eTmp = EmipHits->GetBinContent(ix,iy);
	      if (eTmp < 0.5) continue;
	      EmipMean += eTmp;
	      nHits++;
	      if (lowTail) {
		p_hitSpectrum_lowTail[iE]->Fill(eTmp);
	      }
	      if (highTail) {
		p_hitSpectrum_highTail[iE]->Fill(eTmp);
	      }
	    }
	  }
	  EmipMean = EmipMean/nHits;
	  if (lowTail) p_meanHitSpectrum_lowTail[iE]->Fill(EmipMean);
	  if (highTail) p_meanHitSpectrum_highTail[iE]->Fill(EmipMean);

	  for (unsigned iLim(0); iLim<nLimits;++iLim){
	    nHitsCountNum[iLim] = 0;
	    for (int ix(1); ix<EmipHits->GetNbinsX()+1; ++ix){
	      for (int iy(1); iy<EmipHits->GetNbinsY()+1; ++iy){
		double eTmp = EmipHits->GetBinContent(ix,iy);
		if (eTmp < 0.5) continue;
		if (eTmp<pElim[iLim]) nHitsCountNum[iLim]++;
		if (iLim==0 && eTmp<EmipMean) nHitsCountDen++;
	      }
	    }
	    double Cglobal = 0;
	    if (nHitsCountDen>0) Cglobal = nHitsCountNum[iLim]*1.0/nHitsCountDen;
	    p_Cglobal[iE][iLim]->Fill(Cglobal);
	    p_EvsCglobal[iE][iLim]->Fill(Cglobal,EReco_FHCAL);
	    p_Eshower[iE][iLim]->Fill(mycalib.hcalShowerEnergy(Cglobal,Etotal[1],Etotal[2],false));
	    if (iLim==0) p_EshowerCor[iE]->Fill(mycalib.showerEnergy(Cglobal,Etotal[0],Etotal[1],Etotal[2]));
	  }
	  p_Etotal[iE]->Fill(EReco);
	}

	for (unsigned iD(0); iD<3; ++iD){
	  Etotal[iD] = 0;
	}
	for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
	  Etot[iL] = 0;
	}

	EmipHits->Reset();
	
	nTotal = 0;
	nTotalSignal = 0;
	firstEvent = false;
      }//save histos

    }//loop on entries
  
    outputFile->cd();

    p_hitSpectrum_lowTail[iE]->Write();
    p_hitSpectrum_highTail[iE]->Write();
    p_meanHitSpectrum_lowTail[iE]->Write();
    p_meanHitSpectrum_highTail[iE]->Write();
    for (unsigned iLim(0); iLim<nLimits;++iLim){
      p_Cglobal[iE][iLim]->Write();
      p_EvsCglobal[iE][iLim]->Write();
      p_Eshower[iE][iLim]->Write();
    }
    p_EshowerCor[iE]->Write();
    p_Etotal[iE]->Write();

    std::cout << " -- Summary of total Energy: " << std::endl
	      << " ---- entries " << p_Etotal[iE]->GetEntries() 
	      << " mean " << p_Etotal[iE]->GetMean() 
	      << " rms " << p_Etotal[iE]->GetRMS() 
	      << " underflows " << p_Etotal[iE]->GetBinContent(0)
	      << " overflows " << p_Etotal[iE]->GetBinContent(p_Etotal[iE]->GetNbinsX()+1)
	      << std::endl
	      << " -- Summary of total Energy corrected: " << std::endl
	      << " ---- entries " << p_EshowerCor[iE]->GetEntries() 
	      << " mean " << p_EshowerCor[iE]->GetMean() 
	      << " rms " << p_EshowerCor[iE]->GetRMS() 
	      << " underflows " << p_EshowerCor[iE]->GetBinContent(0)
	      << " overflows " << p_EshowerCor[iE]->GetBinContent(p_EshowerCor[iE]->GetNbinsX()+1)
	      << std::endl
	      << " -- Summary of hit spectra low tail: " << std::endl
	      << " ---- entries " << p_hitSpectrum_lowTail[iE]->GetEntries() 
	      << " mean " << p_hitSpectrum_lowTail[iE]->GetMean() 
	      << " rms " << p_hitSpectrum_lowTail[iE]->GetRMS() 
	      << " underflows " << p_hitSpectrum_lowTail[iE]->GetBinContent(0)
	      << " overflows " << p_hitSpectrum_lowTail[iE]->GetBinContent(p_hitSpectrum_lowTail[iE]->GetNbinsX()+1)
	      << std::endl
	      << " -- Summary of hit spectra high tail: " << std::endl
	      << " ---- entries " << p_hitSpectrum_highTail[iE]->GetEntries() 
	      << " mean " << p_hitSpectrum_highTail[iE]->GetMean() 
	      << " rms " << p_hitSpectrum_highTail[iE]->GetRMS() 
	      << " underflows " << p_hitSpectrum_highTail[iE]->GetBinContent(0)
	      << " overflows " << p_hitSpectrum_highTail[iE]->GetBinContent(p_hitSpectrum_highTail[iE]->GetNbinsX()+1)
	      << std::endl;


  }//loop on energies

  if (doMeanRms) meanRmsOut.close();
  //if (!doMeanRms && doHitSpectra) hitSpectraOut.close();

  //if (!doMeanRms && !doHitSpectra) {
  if (!doMeanRms) {
    outputFile->Write();
    std::cout << " Outputfile written: " << outputFile->GetName() << std::endl;
  }

  return 0;


}//main
