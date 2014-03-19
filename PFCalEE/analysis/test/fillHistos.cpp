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

// double zlayer(const unsigned layer) {
//   double w1b = 1.6+3+0.1;
//   double w2b = 3.3+3+0.1;
//   double w3b = 5.6+3+0.1;
//   double wa = 0.1+4;
//   double W1 = 10*(w1b+wa);
//   double W2 = 10*(w2b+wa);
//   //double W3 = 10*(w3b+wa);

//   double totalWidth = 10*(w1b+w2b+w3b)+30*wa;
//   unsigned il = layer%10+1;

//   double offset = layer/10==0 ? (il*w1b + (il-1)*wa) :
//     ( layer/10==1 ? W1 + il*w2b + (il-1)*wa :
//       W1 + W2 + il*w3b + (il-1)*wa );

//   return -totalWidth/2. + offset; //in mm
  
// }

double getBinWidth(const double etaMin,const double etaMax, std::string etaStr, const double zlaymm){
  double eta0 = 0;
  if (etaStr.find("eta20")!=etaStr.npos) eta0 = 2.0;
  else if (etaStr.find("eta25")!=etaStr.npos) eta0 = 2.5;
  else if (etaStr.find("eta30")!=etaStr.npos) eta0 = 3.0;
  else if (etaStr.find("eta35")!=etaStr.npos) eta0 = 3.5;
  else {
     std::cout << " ERROR: new eta point ! Please implement value to convert y position... Exiting..." << std::endl;
    exit(1);
  }

  double xDist = 310+zlaymm/10.;//distance in cm from beam spot to layer.
  double theta0 = 2*atan(exp(-eta0));
  double y0 = tan(theta0)*xDist;//cm

  double theta1 = 2*atan(exp(-etaMin));
  double ypos1 = xDist*tan(theta1)-y0;

  double theta2 = 2*atan(exp(-etaMax));
  double ypos2 = xDist*tan(theta2)-y0;

  double dy = ypos1-ypos2;//cm

  // std::cout << " - etaBin " << etaStr 
  // 	    << " etaMin = " << etaMin 
  // 	    << " etaMax = " << etaMax
  // 	    << " dy = " << dy
  // 	    << std::endl;

  return dy;

}

double getEta(double ypos, std::string etaStr, const double zlaymm){
  double eta0 = 0;
  if (etaStr.find("20")!=etaStr.npos) eta0 = 2.0;
  else if (etaStr.find("25")!=etaStr.npos) eta0 = 2.5;
  else if (etaStr.find("30")!=etaStr.npos) eta0 = 3.0;
  else if (etaStr.find("35")!=etaStr.npos) eta0 = 3.5;
  else {
    //std::cout << " ERROR: new eta point ! Please implement value to convert y position... Setting it to 10..." << std::endl;
     //exit(1);
     eta0=10;
  }
  //find height using center of detector
  double xDist0 = 3100;
  double theta0 = 2*atan(exp(-eta0));
  double y0 = tan(theta0)*xDist0;//mm
  //find theta at center of layer.
  double xDist = 3100+zlaymm;//distance from beam spot to layer in mm
  double tanTheta = (y0+ypos)/xDist;
  double theta = atan(tanTheta);
  double eta = -log(tan(theta/2.));
  //std::cout << " - etaBin " << etaStr << " y0 = " << y0 << " y = " << ypos << " eta = " << eta << std::endl;

  return eta;

}

double getWeight(unsigned layer,std::string aVersion, double volX0ref, double volX0, bool toGeV=false){
  double mipHCALSi = 0.0849;//MeV
  double mipHCALSci = 1.49;//MeV
  //double calibHCALSi = 1/3.144;//MeV/GeV
  double calibHCALSi = 1/1.496;//MeV/GeV
  double calibHCALSci = 1/14.778;//MeV/GeV
  if (aVersion.find("22")!=aVersion.npos){
    //just scint HCAL
    return toGeV? calibHCALSci : 1;
  }
  if (aVersion.find("21")!=aVersion.npos){
    if (layer < 24)
      return toGeV? calibHCALSi: 1;
    //return 1;
    else if (layer < 33)
      //return toGeV? calibHCAL*mipHCAL*volX0/volX0ref : volX0/volX0ref;
      return toGeV? calibHCALSci*92.196/65.235 : 92.196/65.235 ;
    //return calibHCAL*mipHCAL*volX0/volX0ref;
  }
  if (aVersion.find("23")!=aVersion.npos){
    calibHCALSci = 1/29.815;//MeV/GeV
    if (layer < 47)
      return toGeV? calibHCALSci: 1;
    //return 1;
    else if (layer < 54)
      //return toGeV? calibHCAL*mipHCAL*volX0/volX0ref : volX0/volX0ref;
      return toGeV? calibHCALSci*104./22. : 104./22. ;
    //return calibHCAL*mipHCAL*volX0/volX0ref;
  }

  if (layer<11) 
    return toGeV? 107.84 : 1;
  //if (aVersion.find("20")!= aVersion.npos || 
  //aVersion.find("21") != aVersion.npos ||
  //aVersion.find("22") != aVersion.npos ||
  //aVersion.find("23") != aVersion.npos){
  else if (layer < 21) 
    return toGeV? 136.31 : 1.264;
  //return 8.001/5.848;
  //return 0.8/0.5;
  else if (layer<31) 
    return toGeV? 218.79 : 2.029;
  //return 10.854/5.848;
  //return 1.2/0.5;
  else if (layer < 55)
    return toGeV? calibHCALSi: 1;
  //return 65.235/5.848;
  else if (layer < 64)
    return toGeV? calibHCALSci*92.196/65.235 : 92.196/65.235;
  //return 40.779/5.848;

  return 1;
  //}
  //else {
  //  if (layer < 20) return 0.8/0.4;
  //  else if (layer < 30) return 1.2/0.4;
  //  else return 1;
  //}
}


int main(int argc, char** argv){//main  

  if (argc < 4) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <full path to input file>"
	      << " <file name (PFcal.root, or DigiPFcal.root)>" 
	      << " <optional: list of energies: 5,10,25>"
	      << " <optional: debug (default=0)>"
	      << " <optional: overlayPU (default=0)>"
	      << " <optional: saveEventByEventThresh (default=0)>"
	      << " <optional: doChargedPionsOnly (default=0)>"
	      << " <optional: doNeutralPionsOnly (default=0)>"
	      << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  bool applyWeight = true;
  bool saveTxtOutput = false;
  bool doFiducialCuts = false;
  bool selectPiHCAL = false;
  bool calibrate = false;
  bool concept = true;

  bool selectEarlyDecay = false;

  const unsigned nLayers = 64;//64;//Calice 54;//Scint 9;//HCAL 33;//All 64
  const unsigned nEcalLayers = 31;//31;
  const unsigned nHcalSiLayers = 24;//concept 24;//calice 47
  const double signalRegionInX=20;
  //const double Emip = 0.0548;//MeV
  const double HcalToEcalConv = 1;//1/1.772;
  const double HcalSciToHcalSiConv = 1;//1/1.33;//1/5.12;

  //double minX=-250,maxX=250;
  double minX=-1000,maxX=1000;
  //double minY=150,maxY=800;
  //double minY=1100,maxY=1500;
  double minY=-1000,maxY=1000;
  double minZ=3170,maxZ=3500;
  unsigned nX=(maxX-minX)/2.5,nY=(maxY-minY)/2.5;
  if (nLayers != nEcalLayers) {
    maxZ = 5070;
  }
  unsigned nZ=maxZ-minZ;

  std::cout << " -- 2-D histograms: " << std::endl
	    << " -- X: " << nX << " " << minX << " " << maxX << std::endl
	    << " -- Y: " << nY << " " << minY << " " << maxY << std::endl
	    << " -- Z: " << nZ << " " << minZ << " " << maxZ << std::endl
    ;

  const unsigned nOcc = 1;//4;
  const unsigned occThreshold[nOcc] = {1};//,5,10,20};

  const unsigned nSmear = 1;
  const double smearFrac[nSmear] = {0};//,0.01,0.02,0.03,0.04,0.05,0.07,0.10,0.15,0.20};

  const unsigned nLat = 1;//20;
  const double latSize[nLat] = {100};//{10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};//+/- in mm

  bool saveLLR = false;
  double etaMinLLR = 3.8;
  double etaMaxLLR = 3.9;

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string fileName = argv[3];//"PFcal.root";
  unsigned debug = 0;
  if (argc >5) debug = atoi(argv[5]);
  bool doChargedPions = false;
  bool doNeutralPions = false;
  bool saveEventByEvent = false;
  bool overlayPU = false;
  if (argc>6) overlayPU = atoi(argv[6]);
  unsigned saveEvtThresh = 0;
  if (argc>7) {
    saveEvtThresh = atoi(argv[7]);
    if (saveEvtThresh>0) saveEventByEvent = true;
  }
  if (argc>8) doChargedPions = atoi(argv[8]);
  if (argc>9) doNeutralPions = atoi(argv[9]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Input file name: " << fileName << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  if (overlayPU) std::cout << " -- Overlaying PU events from Minbias bank" << std::endl;
  if (doChargedPions) std::cout << " -- Using hits from charged pions only" << std::endl;
  if (doNeutralPions) std::cout << " -- Using hits from neutral pions only" << std::endl;

  unsigned nMaxHits = saveEvtThresh;//10000;
  if (saveEventByEvent) std::cout << " -- saving event by event for events with nHits above " << nMaxHits <<std::endl;

  TRandom3 lRndm(0);

  //unsigned tmpEn[]={5,10,20,25,50,75,100,125,150,175,200,300,500};
  //unsigned tmpEn[]={5,10,25,40,50,60,80,100,200,300,400,500,1000,2000};
  //unsigned tmpEn[]={5,10,20,25,40,50,60,75,80,100,125,150,175,200,300,400,500,1000};//,1000,2000};
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
  
  bool isHCALonly = filePath.find("version21")!=filePath.npos;
  bool isCalice = filePath.find("version23")!=filePath.npos;

  const double totalWidth=(filePath.find("version20")==filePath.npos && 
			   filePath.find("version21")==filePath.npos) ? 200 : 500;


  std::cout << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----- Hardcoded configuration is : " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " -- applyWeight = " << applyWeight << std::endl
	    << " -- saveTxtOutput = " <<saveTxtOutput << std::endl
	    << " -- doFiducialCuts = " <<doFiducialCuts << std::endl
	    << " -- selectPiHCAL = " <<selectPiHCAL << std::endl
	    << " -- calibrate = " <<calibrate << std::endl
	    << " -- concept = " <<concept << std::endl
	    << " -- N layers = " << nLayers << std::endl
	    << " -- nEcalLayers = " << nEcalLayers << " weights " << getWeight(5,isHCALonly?"version3":filePath,1,1,calibrate) << " " << getWeight(15,isHCALonly?"version3":filePath,1,1,calibrate) << " " << getWeight(25,isHCALonly?"version3":filePath,1,1,calibrate) << std::endl
	    << " -- nHcalSiLayers  = " << nHcalSiLayers << " weight " << getWeight(isHCALonly?1:nEcalLayers+1,filePath,1,1,calibrate) << " Scint " << getWeight(isHCALonly?nHcalSiLayers+1:nEcalLayers+nHcalSiLayers+1,filePath,1,1,calibrate) << std::endl
	    << " -- conversions: HcalToEcalConv = " <<HcalToEcalConv << " HcalSciToHcalSiConv = " << HcalSciToHcalSiConv << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " -----------------------------------" << std::endl;

  TString plotDir = "PLOTS/";
  std::vector<std::string> pathVec;
  bool pEOS = false;
  if (filePath.find("eos") != filePath.npos) {
    pEOS = true;
  }
   
  if (filePath.find("version_23") != filePath.npos) pEOS = false;

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
  if (overlayPU) plotDir += "/PU";
  if (!doNeutralPions && doChargedPions) plotDir += "_pipm";
  else if (!doChargedPions && doNeutralPions) plotDir += "_pi0";
  plotDir += "/";
  plotDir += pEta;
  plotDir += "/";

  bool isMB = false;
  if (pParticle.find("PedroPU")!=pParticle.npos) isMB=true;

  if (!applyWeight) plotDir += "noWeights/";
  if (selectPiHCAL) plotDir += "selectPiHCAL/";
  if (concept) plotDir += "concept/";
  else if (isHCALonly) plotDir += "twiceSampling/";
  if (calibrate) plotDir += "GeVCal/";
  //plotDir += "simpleWeights/";
  if (selectEarlyDecay) plotDir += "EarlyDecay/";

  std::cout << " -- Output file directory is : " << plotDir << std::endl;

  TFile *outputFile = TFile::Open(plotDir+"/CalibHistos.root","RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << plotDir << "/CalibHistos.root cannot be opened. Please create output directory : " << plotDir << "  Exiting..." << std::endl;
    return 1;
  }


  TH2F *p_xy[nGenEn][nLayers];
  TH3F *p_xyz[nGenEn];
  TH2F *p_xy_evt[nLayers];
  TH2F *p_xz_evt[nLayers];
  TH2F *p_zy_evt[nLayers];
  TH3F *p_xyz_evt = 0;
  TH2F *p_recoxy[nGenEn][nLayers];
  TH3F *p_recoxyz[nGenEn];
  TH2F *p_recoxy_evt[nLayers];

  //std::cout << " -- Checking layer-eta conversion:" << std::endl;

  for (unsigned iL(0); iL<nLayers; ++iL){
    p_recoxy_evt[iL] = 0;
    p_xy_evt[iL] = 0;
    p_xz_evt[iL] = 0;
    p_zy_evt[iL] = 0;
    //std::cout << iL << " " << zlay[iL] << " " 
    //<< llrBinWidth[iL] << " y=0@eta3 " 
    //<< getEta(0,"eta30",zlay[iL]) << std::endl;
  }

  TH3F *p_recoxyz_evt = 0;
  TH2F *p_EvsLayer[nGenEn];
  TH1F *p_Etot[nGenEn][nLayers];
  TH1F *p_Efrac[nGenEn][nLayers];
  TH1F *p_time[nGenEn][nLayers];
  TH1F *p_Etotal[nGenEn][3];
  TH1F *p_EfracLateral[nGenEn][2][nLat];
  TH2F *p_HCALvsECAL[nGenEn];

  TH1F *p_Ereco[nGenEn][nSmear];
  TH1F *p_nSimHits[nGenEn];
  TH1F *p_nRecHits[nGenEn];

  TH1F *p_Edensity[nGenEn][nLayers];
  TH2F *p_occupancy[nOcc][nLayers];

  std::ostringstream lName;

  for (unsigned iO(0); iO<nOcc; ++iO){
    for (unsigned iL(0); iL<nLayers; ++iL){
      lName.str("");
      lName << "p_occupancy_" << occThreshold[iO] << "_" << iL;
      p_occupancy[iO][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",totalWidth/10.,-totalWidth/2.,totalWidth/2.,totalWidth/10.,-totalWidth/2.,totalWidth/2.);
      //p_occupancy[iO][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",100,getEta(totalWidth/2.,pEta),getEta(-totalWidth/2.,pEta),100,);
    }
  }


  bool isG4Tree = true;
  
  //for basic control plots
  TCanvas *mycE = new TCanvas("mycE","mycE",1500,1000);
  mycE->Divide(5,2);
  gStyle->SetOptStat("eMRuoi");

  for (unsigned iE(0); iE<nGenEn; ++iE){

    if (isMB) genEn[iE] = iE;

    std::cout << "- Processing energy : " << genEn[iE] << std::endl;
    TString genEnStr;
    genEnStr += genEn[iE];
    std::ofstream hitsOut;
    if (saveTxtOutput) hitsOut.open(plotDir+"/hits_"+genEnStr+"GeV.dat");
    std::ofstream optim;
    if (saveTxtOutput) optim.open(plotDir+"/optimisation_"+genEnStr+"GeV.dat");
    std::ofstream edensity;
    genEnStr = "";
    genEnStr += iE;
    if (saveLLR) edensity.open(plotDir+"/edensity_"+genEnStr+".dat");

    double E1=0;
    double E2=0;
    double E3=0;

    std::ostringstream input;
    input << filePath ;
    //if (filePath.find("pi+")!=filePath.npos) input << "/pt_" ;
    //HGcal_version_8_e50.root
    if (pEOS) {
      input << "_e" ;
      input << genEn[iE] << fileName;
    }
    else {
      //input << "/e_" ;
      //input << genEn[iE] << "/" << fileName ;
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
    float volX0ref = 0;
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
    
    TFile *puFile = 0;
    TTree *lPuTree = 0;
    std::vector<HGCSSSimHit> * pusimhitvec = 0;
    std::vector<HGCSSRecoHit> * purechitvec = 0;

    unsigned nPuEvts= 0;

    if (overlayPU){
      std::string inStr;// = "/afs/cern.ch/work/a/amagnan/public/";
      
      //if (isG4Tree) inStr += "HGCalEEGeant4/"+pVersion+"/"+pScenario+"/PedroPU/"+pEta+"/PFcal.root";
      //else  inStr += "HGCalEEDigi/"+pVersion+"/"+pScenario+"/PedroPU/"+pEta+"/DigiPFcal.root";
      if (isG4Tree) inStr += "140PU/PFcal_140PU_EEHE.root";
      puFile = TFile::Open(inStr.c_str());
      if (!puFile) {
	std::cout << " -- Error, input file " << inStr << " for PU cannot be opened. Exiting..." << std::endl;
	return 1;
      }
      else std::cout << " -- using file " << inStr << " for PU." << std::endl;
      lPuTree = (TTree*)puFile->Get("PUTree");
      if (!lPuTree){
	std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
	return 1;
      }

      if (isG4Tree){
	lPuTree->SetBranchAddress("HGCSSSimHitVec",&pusimhitvec);
      }
      else {
	lPuTree->SetBranchAddress("HGCSSSimHitVec",&pusimhitvec);
	lPuTree->SetBranchAddress("HGCSSRecoHitVec",&purechitvec);
      }
      nPuEvts =  isG4Tree ? 
	static_cast<unsigned>(lPuTree->GetEntries()/nLayers):
	static_cast<unsigned>(lPuTree->GetEntries());
      //nPuEvts =  static_cast<unsigned>(lPuTree->GetEntries());

      std::cout << "- For PU: " << nPuEvts  << " single-interaction events are available." << std::endl;
      
    }


    //Initialise histos
    //necessary to have overflows ?
    gStyle->SetOptStat(1111111);
    double Etot[nLayers];
    double eLLR[nLayers];
    double zlay[nLayers];

    for (unsigned iL(0);iL<nLayers;++iL){
      Etot[iL] = 0;
      eLLR[iL] = 0;
      zlay[iL] = 0;
    }
    double Etotal[3] = {0,0,0};
    double Elateral[3][20];
    for (unsigned iD(0); iD<3; ++iD){
      for (unsigned iR(0); iR<nLat; ++iR){
	Elateral[iD][iR] = 0;
      }
    }
    unsigned nTotal = 0;
    unsigned nTotalSignal = 0;
    double Ereco[nSmear];
    for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
      Ereco[iSmear]= 0;
    }
 

    lName.str("");
    lName << "p_xyz_" << genEn[iE];
    p_xyz[iE] = new TH3F(lName.str().c_str(),";z(mm);x(mm);y(mm)",2000,-1000,1000,200,-250,250,200,-250,250);
    if (!isG4Tree){
      lName.str("");
      lName << "p_recoxyz_" << genEn[iE];
      p_recoxyz[iE] = new TH3F(lName.str().c_str(),";z(mm);x(mm);y(mm)",2000,-1000,1000,50,-250,250,50,-250,250);
    }
    else p_recoxyz[iE] = 0;
    p_Etotal[iE][0] = 0;
    p_Etotal[iE][1] = 0;
    p_Etotal[iE][2] = 0;
    p_nSimHits[iE] = 0;
    p_EvsLayer[iE] = 0;
    p_nRecHits[iE] = 0;
    for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
      p_Ereco[iE][iSmear] = 0;
    }
    for (unsigned iR(0); iR<nLat; ++iR){
      lName.str("");
      lName << "p_EfracLateral_" << genEn[iE] << "_ECAL_" << iR;
      p_EfracLateral[iE][0][iR] =  new TH1F(lName.str().c_str(),";E_{total}^{x #times x}/E_{total}",101,0,1.01);
      lName.str("");
      lName << "p_EfracLateral_" << genEn[iE] << "_HCAL_" << iR;
      p_EfracLateral[iE][1][iR] =  new TH1F(lName.str().c_str(),";E_{total}^{x #times x}/E_{total}",101,0,1.01);
    }
    for (unsigned iL(0); iL<nLayers; ++iL){
      lName.str("");
      lName << "p_xy_" << genEn[iE] << "_" << iL;
      p_xy[iE][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",200,-250,250,200,-250,250);
      if (!isG4Tree){
	lName.str("");
	lName << "p_recoxy_" << genEn[iE] << "_" << iL;
	p_recoxy[iE][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",50,-250,250,50,-250,250);
      }
      else p_recoxy[iE][iL] = 0;
      Etot[iL] = 0;
      lName.str("");
      lName << "p_Efrac_" << genEn[iE] << "_" << iL;
      p_Efrac[iE][iL] = new TH1F(lName.str().c_str(),";integrated E_{layer}/E_{total}",101,0,1.01);
      lName.str("");
      lName << "p_time_" << genEn[iE] << "_" << iL;
      p_time[iE][iL] = new TH1F(lName.str().c_str(),";G4 time (ns)",500,0,20);
      lName.str("");
      lName << "p_Edensity_" << genEn[iE] << "_" << iL;
      p_Edensity[iE][iL] = new TH1F(lName.str().c_str(),";#eta",130,1.8,4.4);
    }



    bool firstEvent = true;
    unsigned ipuevt = 0;


    for (unsigned ievt(0); ievt<(isG4Tree?nEvts*nLayers:nEvts); ++ievt){//loop on entries
      if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
      else if (ievt%1000 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
      lTree->GetEntry(ievt);

      if (volNb==0) volX0ref = volX0;
      else if (volNb == nEcalLayers) volX0ref = volX0;

      if (debug){
	if (isG4Tree) {
	  std::cout << "... Processing layer " << volNb << " with " << (*simhitvec).size() << " simhits, X0 = " << volX0 << " weight = " << getWeight(volNb,pVersion,volX0ref,volX0,calibrate) << std::endl;
	}
	else {
	  std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
	}
      }
      
      double Emax = 0;

      unsigned lEvt = isG4Tree? ievt/nLayers : ievt;
 
      if (saveEventByEvent && 
	  (!isG4Tree || (isG4Tree && (volNb == 0)))
	  ){


	if (p_xyz_evt) p_xyz_evt->Delete();
	p_xyz_evt = new TH3F("p_xyz",";z(mm);x(mm);y(mm)",nZ/2,minZ,maxZ,nX/4,minX,maxX,nY/4,minY,maxY);
	if (!isG4Tree) {
	  if (p_recoxyz_evt) p_recoxyz_evt->Delete();
	  p_recoxyz_evt = new TH3F("p_recoxyz",";z(mm);x(mm);y(mm)",nZ,minZ,maxZ,nX/4,minX,maxX,nY/4,minY,maxY);
	}
	
	for (unsigned iL(0); iL<nLayers; ++iL){
	  lName.str("");
	  lName << "p_xy_" << iL;
	  if (p_xy_evt[iL]) p_xy_evt[iL]->Delete();
	  p_xy_evt[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",nX,minX,maxX,nY,minY,maxY);
	  lName.str("");
	  lName << "p_xz_" << iL;
	  if (p_xz_evt[iL]) p_xz_evt[iL]->Delete();
	  p_xz_evt[iL] = new TH2F(lName.str().c_str(),";x(mm);z(mm)",nX,minX,maxX,nZ,minZ,maxZ);
	  lName.str("");
	  lName << "p_zy_" << iL;
	  if (p_zy_evt[iL]) p_zy_evt[iL]->Delete();
	  p_zy_evt[iL] = new TH2F(lName.str().c_str(),";z(mm);y(mm)",nZ,minZ,maxZ,nY,minY,maxY);
	  if (!isG4Tree){
	    lName.str("");
	    lName << "p_recoxy_" << iL;
	    if (p_recoxy_evt[iL]) p_recoxy_evt[iL]->Delete();
	    p_recoxy_evt[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",nX/4,minX,maxX,nY/4,minY,maxY);

	  }
	}
      }

      //add PU hits
      if (overlayPU){//add PU hits
	if (isG4Tree? (ievt%nLayers==0) : true) {
	  ipuevt = lRndm.Integer(nPuEvts);
	  //get random PU events among available;
	  std::cout << " -- PU Random number for entry " << ievt << " is " << ipuevt << std::endl;
	}
	
	lPuTree->GetEntry(isG4Tree? nLayers*ipuevt+(ievt%nLayers) : ipuevt);
	for (unsigned iH(0); iH<(*pusimhitvec).size(); ++iH){//loop on hits
	  HGCSSSimHit lHit = (*pusimhitvec)[iH];
	  unsigned layer = lHit.layer();
	  if (isG4Tree && layer==volNb+1) {
	    if (firstEvent) std::cout << " -- Warning PU hit#" << iH << ", applying patch to layer number..." << layer << " " << volNb << std::endl; 
	    layer = volNb;
	  }
	    
	  //select "flavour"
	  if (doChargedPions){
	    //if ((lHit.nHadrons() == 0 && lHit.nMuons() == 0)) continue;
	    //HACK !!
	    if (lHit.nProtons() == 0) continue;
	  }
	  else if (doNeutralPions){
	    //if (lHit.nGammas()==0 && lHit.nElectrons()==0 ) continue;
	    //HACK !!
	    if (lHit.nNeutrons()==0 ) continue;
	  }
	    
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  double posz = lHit.get_z();
	  //double posx = lHit.get_y();
	  //double posy = lHit.get_x();
	  double energy = lHit.energy();
	  if (debug>1) {
	    std::cout << " --  SimHit " << iH << "/" << (*pusimhitvec).size() << " --" << std::endl
		      << " --  position x,y " << posx << "," << posy << std::endl;
	    lHit.Print(std::cout);
	  }
	  double weightedE = energy;
	  if (applyWeight) weightedE *= getWeight(layer,pVersion,volX0ref,volX0,calibrate);
	  //reject lower region for boundary effects
	  if (doFiducialCuts){
	    if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
		 (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	    //region left and right also for boundary effects
	    if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
		 (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	  }
	  p_xy[iE][layer]->Fill(posx,posy,energy);
	  p_xyz[iE]->Fill(posz,posx,posy,energy);
	    
	  if (saveEventByEvent){
	    bool pass = (layer < nEcalLayers) ||
	      (layer>= nEcalLayers && 
	       layer<nEcalLayers+nHcalSiLayers &&
	       (!concept || (concept && layer%2==1))
	       ) ||
	      layer>=nEcalLayers+nHcalSiLayers;
	    if (pass){
	      p_xy_evt[layer]->Fill(posx,posy,energy);
	      p_xz_evt[layer]->Fill(posx,posz,energy);
	      p_zy_evt[layer]->Fill(posz,posy,energy);
	      p_xyz_evt->Fill(posz,posx,posy,energy);
	    }
	  }
	  p_time[iE][layer]->Fill(lHit.time());
	  //restrict in y and x to have just the influence of the eta bin wanted
	  //(at high eta, surface of the cone in x is smaller than detector size)
	  if (fabs(posx)<signalRegionInX/2.){
	    double lEtaTmp = getEta(posy,pEta,posz);
	    p_Edensity[iE][layer]->Fill(lEtaTmp,energy);
	  }
	  Etot[layer] += energy;
	  if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << Etot[layer];
	    
	  if (isCalice){
	    if (layer<nHcalSiLayers) {
	      Etotal[0] += weightedE;
	      Etotal[2] += weightedE;
	    }
	    else {
	      Etotal[1] += weightedE;
	      Etotal[2] += weightedE*HcalSciToHcalSiConv;
	    }
	    for (unsigned iR(0); iR<nLat; ++iR){
	      if (fabs(posx)<latSize[iR] &&
		  fabs(posy)<latSize[iR]) Elateral[0][iR]+= weightedE;
	    }
	  }
	  else if (isHCALonly) {
	    if (layer<nHcalSiLayers) {
	      if (!concept || 
		  (concept && layer%2==1)){
		Etotal[0] += weightedE;
		Etotal[2] += weightedE;
		for (unsigned iR(0); iR<nLat; ++iR){
		  if (fabs(posx)<latSize[iR] &&
		      fabs(posy)<latSize[iR]) Elateral[0][iR]+= weightedE;
		}
		  
	      }
	    }
	    else {
	      //twice sampling: correct weight.
	      if (!concept) weightedE = weightedE*2.;
	      Etotal[1] += weightedE;
	      Etotal[2] += weightedE*HcalSciToHcalSiConv;
	      for (unsigned iR(0); iR<nLat; ++iR){
		if (fabs(posx)<latSize[iR] &&
		    fabs(posy)<latSize[iR]) Elateral[0][iR]+= weightedE;
	      }
	    }
	      
	  }//HCALonly
	  else {//ECALonly or ECAL+HCAL
	    if (layer<nEcalLayers) {
	      Etotal[0] += weightedE;
	      Etotal[2] += weightedE;
	      for (unsigned iR(0); iR<nLat; ++iR){
		if (fabs(posx)<latSize[iR] &&
		    fabs(posy)<latSize[iR]) Elateral[0][iR]+= weightedE;
	      }
	    }//ECAL
	    else {
	      if (layer<nEcalLayers+nHcalSiLayers) {
		if (!concept ||
		    (concept && layer%2==1)){
		  Etotal[1] += weightedE;
		  Etotal[2] += weightedE*HcalToEcalConv;
		  for (unsigned iR(0); iR<nLat; ++iR){
		    if (fabs(posx)<latSize[iR] &&
			fabs(posy)<latSize[iR]) Elateral[1][iR]+= weightedE;
		  }
		}
	      }
	      else {
		//twice sampling: correct weight.
		if (!concept) weightedE = weightedE*2.;
		Etotal[1] += weightedE*HcalSciToHcalSiConv;
		Etotal[2] += weightedE*HcalSciToHcalSiConv*HcalToEcalConv;
		for (unsigned iR(0); iR<nLat; ++iR){
		  if (fabs(posx)<latSize[iR] &&
		      fabs(posy)<latSize[iR]) Elateral[1][iR]+= weightedE;
		}
	      }
	    }//HCAL
	  }//ECALonly or ECAL+HCAL
	    
	}//loop on hits
	  
	  
	nTotal += (*pusimhitvec).size();
	if (debug)  std::cout << std::endl;
	if (!isG4Tree){
	  for (unsigned iSmear(0); iSmear<nSmear;++iSmear){	  
	    Ereco[iSmear] = 0;
	  }
	  for (unsigned iH(0); iH<(*purechitvec).size(); ++iH){//loop on rechits
	    HGCSSRecoHit lHit = (*purechitvec)[iH];
	    if (debug>1) {
	      std::cout << " --  RecoHit " << iH << "/" << (*purechitvec).size() << " --" << std::endl
			<< " --  position x,y " << lHit.get_x() << "," << lHit.get_y() << std::endl;
	      lHit.Print(std::cout);
	    }
	    double posx = lHit.get_x();
	    double posy = lHit.get_y();
	    double posz = lHit.get_z();
	    if (doFiducialCuts){
	      //reject lower region for boundary effects
	      if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
		   (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	      //region left and right also for boundary effects
	      if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
		   (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	    }
	      
	    double energy = lHit.energy();//in MIPs already
	    unsigned layer = lHit.layer();
	    p_recoxy[iE][layer]->Fill(posx,posy,energy);
	    p_recoxyz[iE]->Fill(posz,posx,posy,energy);
	    for (unsigned iO(0); iO<nOcc; ++iO){
	      if (energy > occThreshold[iO]) p_occupancy[iO][layer]->Fill(posx,posy);
	    }
	    double weightedE = energy;
	    if (applyWeight)  weightedE *= getWeight(layer,pVersion,volX0ref,volX0,calibrate);
	    for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	      double smearFactor = lRndm.Gaus(0,smearFrac[iSmear]);
	      Ereco[iSmear] += weightedE*(1+smearFactor);
	    }
	      
	    if (saveEventByEvent){
	      p_recoxy_evt[layer]->Fill(posx,posy,energy);
	      p_recoxyz_evt->Fill(posz,posx,posy,energy);
	    }
	      
	      
	  }
	}
      }//add PU hits

      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	HGCSSSimHit lHit = (*simhitvec)[iH];
	unsigned layer = lHit.layer();
	if (isG4Tree && layer==volNb+1) {
	  if (firstEvent) std::cout << " -- Warning, applying patch to layer number..." << std::endl; 
	  layer = volNb;
	}

	//select "flavour"
	if (doChargedPions){
	  //if ((lHit.nHadrons() == 0 && lHit.nMuons() == 0)) continue;
	  //HACK !!
	  if (lHit.nProtons()==0 ) continue;
	}
	else if (doNeutralPions){
	  //HACK !!
	  if (lHit.nNeutrons()==0 ) continue;
	  //if (lHit.nGammas()==0 && lHit.nElectrons()==0 ) continue;
	}

	double posx = lHit.get_x();
	double posy = lHit.get_y();
	double posz = lHit.get_z();
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!! WARNING, quick and dirty fix !!!!!!!!!!!!!!!!!!
	//double posx = lHit.get_y();
	//double posy = lHit.get_x();
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double energy = lHit.energy();
	if (debug>1) {
	  std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		    << " --  position x,y " << posx << "," << posy << std::endl;
	  lHit.Print(std::cout);
	}
	if (layer<10) E1 += energy;
	else if (layer<20) E2 += energy;
	else E3 += energy;

	if (energy>0 && saveTxtOutput) {
	  if (isG4Tree) hitsOut << static_cast<unsigned>(ievt/nLayers);
	  else hitsOut << ievt;
	  hitsOut << " " << layer << " " << posx << " " << posy << " " << energy << std::endl;
	}

	double weightedE = energy;
	if (applyWeight) weightedE *= getWeight(layer,pVersion,volX0ref,volX0,calibrate);

	if (doFiducialCuts){
	  //reject lower region for boundary effects
	  if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
	       (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	  //region left and right also for boundary effects
	  if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
	       (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	}

	p_xy[iE][layer]->Fill(posx,posy,energy);
	p_xyz[iE]->Fill(posz,posx,posy,energy);
	if (saveEventByEvent){
	  bool pass = (layer < nEcalLayers) ||
	    (layer>= nEcalLayers && 
	     layer<nEcalLayers+nHcalSiLayers &&
	     (!concept || (concept && layer%2==1))
	     ) ||
	    layer>=nEcalLayers+nHcalSiLayers;
	  if (pass){
	    p_xy_evt[layer]->Fill(posx,posy,energy);
	    p_xz_evt[layer]->Fill(posx,posz,energy);
	    p_zy_evt[layer]->Fill(posz,posy,energy);
	    p_xyz_evt->Fill(posz,posx,posy,energy);
	  }
	}
	p_time[iE][layer]->Fill(lHit.time());
	if (fabs(posx)<signalRegionInX/2.){
	  double lEtaTmp = getEta(posy,pEta,posz);
	  p_Edensity[iE][layer]->Fill(lEtaTmp,energy);
	  if (lEtaTmp>=etaMinLLR && lEtaTmp<etaMaxLLR) eLLR[layer]+=energy;
	}
	Etot[layer] += energy;
	if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << Etot[layer];

 
	if (isCalice){
	  if (layer<nHcalSiLayers) {
	    Etotal[0] += weightedE;
	    Etotal[2] += weightedE;
	  }
	  else {
	    Etotal[1] += weightedE;
	    Etotal[2] += weightedE*HcalSciToHcalSiConv;
	  }
	  for (unsigned iR(0); iR<nLat; ++iR){
	    if (fabs(posx)<latSize[iR] &&
		fabs(posy)<latSize[iR]) Elateral[0][iR]+= weightedE;
	  }
	}
	else if (isHCALonly) {
	  if (layer<nHcalSiLayers) {
	    if (!concept || 
		(concept && layer%2==1)){
	      Etotal[0] += weightedE;
	      Etotal[2] += weightedE;
	      for (unsigned iR(0); iR<nLat; ++iR){
		if (fabs(posx)<latSize[iR] &&
		    fabs(posy)<latSize[iR]) Elateral[0][iR]+= weightedE;
	      }
	      
	    }
	  }
	  else {
	    //twice sampling: correct weight.
	    if (!concept) weightedE = weightedE*2.;
	    Etotal[1] += weightedE;
	    Etotal[2] += weightedE*HcalSciToHcalSiConv;
	    for (unsigned iR(0); iR<nLat; ++iR){
	      if (fabs(posx)<latSize[iR] &&
		  fabs(posy)<latSize[iR]) Elateral[0][iR]+= weightedE;
	    }
	  }
	  
	}//HCALonly
	else {//ECALonly or ECAL+HCAL
	  if (layer<nEcalLayers) {
	    Etotal[0] += weightedE;
	    Etotal[2] += weightedE;
	    for (unsigned iR(0); iR<nLat; ++iR){
	      if (fabs(posx)<latSize[iR] &&
		  fabs(posy)<latSize[iR]) Elateral[0][iR]+= weightedE;
	    }
	  }//ECAL
	  else {
	    if (layer<nEcalLayers+nHcalSiLayers) {
	      if (!concept ||
		  (concept && layer%2==1)){
		Etotal[1] += weightedE;
		Etotal[2] += weightedE*HcalToEcalConv;
		for (unsigned iR(0); iR<nLat; ++iR){
		  if (fabs(posx)<latSize[iR] &&
		      fabs(posy)<latSize[iR]) Elateral[1][iR]+= weightedE;
		}
	      }
	    }
	    else {
	      //twice sampling: correct weight.
	      if (!concept) weightedE = weightedE*2.;
	      Etotal[1] += weightedE*HcalSciToHcalSiConv;
	      Etotal[2] += weightedE*HcalSciToHcalSiConv*HcalToEcalConv;
	      for (unsigned iR(0); iR<nLat; ++iR){
		if (fabs(posx)<latSize[iR] &&
		    fabs(posy)<latSize[iR]) Elateral[1][iR]+= weightedE;
	      }
	    }
	  }//HCAL
	}//ECALonly or ECAL+HCAL
       	
	zlay[layer] = posz;

      }//loop on hits

      nTotal += (*simhitvec).size();
      nTotalSignal += (*simhitvec).size();
      if (debug)  std::cout << std::endl;
      if (!isG4Tree){
	if (!overlayPU) {
	  for (unsigned iSmear(0); iSmear<nSmear;++iSmear){	  
	    Ereco[iSmear] = 0;
	  }
	}
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

	  if (doFiducialCuts){
	    //reject lower region for boundary effects
	    if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
		 (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	    //region left and right also for boundary effects
	    if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
		 (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	  }

	  double energy = lHit.energy();//in MIP already...
	  unsigned layer = lHit.layer();
	  for (unsigned iO(0); iO<nOcc; ++iO){
	    if (energy > occThreshold[iO]) p_occupancy[iO][layer]->Fill(posx,posy);
	  }
	  double weightedE = energy;
	  if (applyWeight)  weightedE *= getWeight(layer,pVersion,volX0ref,volX0,calibrate);
	  for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	    double smearFactor = lRndm.Gaus(0,smearFrac[iSmear]);
	    Ereco[iSmear] += weightedE*(1+smearFactor);
	  }

	  p_recoxy[iE][layer]->Fill(posx,posy,energy);
	  p_recoxyz[iE]->Fill(posz,posx,posy,energy);
	  if (saveEventByEvent){
	    p_recoxy_evt[layer]->Fill(posx,posy,energy);
	    p_recoxyz_evt->Fill(posz,posx,posy,energy);
	  }
	}
	if (firstEvent){
	  for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	    lName.str("");
	    lName << "p_Ereco_" << genEn[iE] << "_smear" << iSmear;
	    p_Ereco[iE][iSmear] = new TH1F(lName.str().c_str(),";Reco Etotal (MIPs)",1000,0.1*Ereco[iSmear],2*Ereco[iSmear]);
	  }
	  lName.str("");
	  lName << "p_nRecHits_" << genEn[iE];
	  unsigned nRecHits = (*rechitvec).size();
	  if (overlayPU) nRecHits += (*purechitvec).size();
	  p_nRecHits[iE] = new TH1F(lName.str().c_str(),"; nRecHits",1000,static_cast<unsigned>(nRecHits/10.),nRecHits*6);
	}
	for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	  p_Ereco[iE][iSmear]->Fill(Ereco[iSmear]);
	  Ereco[iSmear] = 0;
	}
	p_nRecHits[iE]->Fill((*rechitvec).size());
	if (debug) std::cout << "... recoE = " << Ereco << std::endl;
      }




      if (!isG4Tree || (isG4Tree && (volNb == nLayers-1))){//fill histos
	//if (debug) std::cout << " -- Filling histograms..." << std::endl;
	bool doFill = true;
	if (selectEarlyDecay) {
	  double tmp = 0;
	  double tmp2 = 0;
	  for (unsigned iL(0);iL<5;++iL){
            tmp += Etot[iL];
	    tmp2 += Etot[iL+5];
	  }
	  if (tmp2 < 10*tmp) doFill = false; 
	}
	if (selectPiHCAL && Etotal[0]>10) doFill = false;

	double Etmp = 0;

	if (saveTxtOutput) optim << E1 << " " << E2 << " " << E3 << std::endl;
	E1=0;
	E2=0;
	E3=0;

	if (firstEvent){
	  for (unsigned iL(0);iL<nLayers;++iL){
	    if (Etot[iL]>Emax) Emax = Etot[iL];
	  }
	  for (unsigned iL(0);iL<nLayers;++iL){
	    lName.str("");
	    lName << "p_Etot_" << genEn[iE] << "_" << iL;
	    p_Etot[iE][iL] = new TH1F(lName.str().c_str(),";G4 Etot (MeV)",1000,0,3*Emax);
	  }
	  lName.str("");
	  lName << "p_EvsLayer_" << genEn[iE];
	  p_EvsLayer[iE] = new TH2F(lName.str().c_str(),";layer;E (MeV)",nLayers,0,nLayers,1000,0,3*Emax);
	  lName.str("");
	  lName << "p_Etotal_" << genEn[iE] << "_ECAL";
	  p_Etotal[iE][0] = new TH1F(lName.str().c_str(),";G4 Etotal ECAL (MeV)",1000,0,10*Etotal[0]);
	  p_Etotal[iE][0]->StatOverflows();
	  lName.str("");
	  lName << "p_Etotal_" << genEn[iE] << "_HCAL";
	  if (isCalice) p_Etotal[iE][1] = new TH1F(lName.str().c_str(),";G4 Etotal HCAL (MeV)",1000,0,10*Etotal[0]/HcalSciToHcalSiConv);
	  else if (isHCALonly) p_Etotal[iE][1] = new TH1F(lName.str().c_str(),";G4 Etotal HCAL (MeV)",1000,0,10*Etotal[0]/HcalSciToHcalSiConv);
	  else p_Etotal[iE][1] = new TH1F(lName.str().c_str(),";G4 Etotal HCAL (MeV)",1000,0,10*Etotal[0]/HcalToEcalConv);
	  p_Etotal[iE][1]->StatOverflows();
	  lName.str("");
	  lName << "p_Etotal_" << genEn[iE] << "_ECALHCAL";
	  p_Etotal[iE][2] = new TH1F(lName.str().c_str(),";G4 Etotal ECAL+HCAL (MeV)",1000,0,10*Etotal[2]);
	  p_Etotal[iE][2]->StatOverflows();
	  lName.str("");
	  lName << "p_HCALvsECAL_" << genEn[iE];
	  if (isCalice) p_HCALvsECAL[iE] = new TH2F(lName.str().c_str(),";G4 Etotal HCAL-Sci (MeV); G4 Etotal HCAL-Sci-TCMT (MeV)",1000,0,5*Etotal[2],1000,0,10*Etotal[0]/HcalSciToHcalSiConv);
	  else if (isHCALonly) p_HCALvsECAL[iE] = new TH2F(lName.str().c_str(),";G4 Etotal HCAL-Si (MeV); G4 Etotal HCAL-Sci (MeV)",1000,0,5*Etotal[2],1000,0,10*Etotal[0]/HcalSciToHcalSiConv);
	  else p_HCALvsECAL[iE] = new TH2F(lName.str().c_str(),";G4 Etotal ECAL (MeV); G4 Etotal HCAL (MeV)",1000,0,3*Etotal[2],1000,0,10*Etotal[0]/HcalToEcalConv);
	  p_HCALvsECAL[iE]->StatOverflows();
	  lName.str("");
	  lName << "p_nSimHits_" << genEn[iE];
	  p_nSimHits[iE] = new TH1F(lName.str().c_str(),"; nSimHits",1000,static_cast<unsigned>(nTotal/10.),nTotal*10);
	}

	for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
	  if (debug) std::cout << " -- Layer " << iL << " total E = " << Etot[iL] << std::endl;
	  if (doFill) {
	    p_Etot[iE][iL]->Fill(Etot[iL]);
	    p_EvsLayer[iE]->Fill(iL,Etot[iL]);
	    Etmp += Etot[iL];
	    if (Etotal[0]+Etotal[1] > 0) p_Efrac[iE][iL]->Fill(Etmp/(Etotal[0]+Etotal[1]));
	    else p_Efrac[iE][iL]->Fill(0);
	  }
	  Etot[iL] = 0;

	  //////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////
	  //For LLR study
	  //get total density in one eta bin
	  if (saveLLR){// && iL>0 && iL<30){
	    double llrBinWidth = getBinWidth(etaMinLLR,etaMaxLLR,"eta35",zlay[iL]);
	    double density = eLLR[iL]/signalRegionInX*1/llrBinWidth;
	    if (isG4Tree) edensity << static_cast<unsigned>(ievt/nLayers);
	    else edensity << ievt;
	    edensity << " " << iL << " " << density << std::endl;
	    eLLR[iL] = 0;
	  }
	  //end-For LLR study
	  //////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////

	}//loop on layers
	if (doFill) {

	  p_HCALvsECAL[iE]->Fill(Etotal[0],Etotal[1]);

	  for (unsigned iD(0); iD<3; ++iD){
	    p_Etotal[iE][iD]->Fill(Etotal[iD]);
	    for (unsigned iR(0); iR<nLat; ++iR){
	      if (isCalice || isHCALonly) {
		if (iD==0 && (Etotal[0]+Etotal[1] > 0) ) p_EfracLateral[iE][0][iR]->Fill(Elateral[0][iR]/(Etotal[0]+Etotal[1]));
	      }
	      else if (iD<2) {
		if (Etotal[iD] > 0) p_EfracLateral[iE][iD][iR]->Fill(Elateral[iD][iR]/Etotal[iD]);
	      }
	    }
	  }

	  p_nSimHits[iE]->Fill(nTotal);

	  if (saveEventByEvent && nTotalSignal >= nMaxHits){
	    std::ostringstream evtName;
	    evtName << plotDir << "/CalibHistos_E" << genEn[iE] << "_evt" << lEvt << ".root";
	    TFile *outputEvt = TFile::Open(evtName.str().c_str(),"RECREATE");
	    
	    if (!outputEvt) {
	      std::cout << " -- Error, output file for evt " << lEvt << " cannot be opened. Please create output directory : " << plotDir << ". Exiting..." << std::endl;
	      return 1;
	    }
	    else {
	      std::cout << " -- Opening event file for evt " << lEvt << std::endl;
	    }
	    outputEvt->cd();
	    for (unsigned iL(0); iL<nLayers; ++iL){
	      p_xy_evt[iL]->Write();
	      p_xz_evt[iL]->Write();
	      p_zy_evt[iL]->Write();
	      if (!isG4Tree) p_recoxy_evt[iL]->Write();
	    }
	    if (p_xyz_evt) p_xyz_evt->Write();
	    if (!isG4Tree) p_recoxyz_evt->Write();
	    outputEvt->Close();
	  }
	}
	for (unsigned iD(0); iD<3; ++iD){
	  Etotal[iD] = 0;
	  for (unsigned iR(0); iR<nLat; ++iR){
	    Elateral[iD][iR] = 0;
	  }
	}
	nTotal = 0;
	nTotalSignal = 0;
	firstEvent = false;
      }//save histos

    }//loop on entries
  
    if (saveTxtOutput) {
      hitsOut.close();
      optim.close();
    }
    if (saveLLR) edensity.close();
    
    for (unsigned iL(0); iL<nLayers; ++iL){
      outputFile->cd();
      p_xy[iE][iL]->Write();
      if (!isG4Tree) p_recoxy[iE][iL]->Write();
      p_time[iE][iL]->Write();
      p_Edensity[iE][iL]->Write();
      p_Etot[iE][iL]->Write();
      p_Efrac[iE][iL]->Write();
    }
    p_xyz[iE]->Write();
    p_EvsLayer[iE]->Write();
    for (unsigned iD(0); iD<3; ++iD){
      p_Etotal[iE][iD]->Write();
      for (unsigned iR(0); iR<nLat; ++iR){
	if (iD<2) p_EfracLateral[iE][iD][iR]->Write();
      }
    }
    p_HCALvsECAL[iE]->Write();
    p_nSimHits[iE]->Write();
    if (!isG4Tree) {
      for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	p_Ereco[iE][iSmear]->Write();
      }
      p_nRecHits[iE]->Write();
      p_recoxyz[iE]->Write();
    }

    std::cout << " -- Summary of energies ECAL: " << std::endl
	      << " ---- SimHits: entries " << p_Etotal[iE][0]->GetEntries() 
	      << " mean " << p_Etotal[iE][0]->GetMean() 
	      << " rms " << p_Etotal[iE][0]->GetRMS() 
	      << " underflows " << p_Etotal[iE][0]->GetBinContent(0)
	      << " overflows " << p_Etotal[iE][0]->GetBinContent(p_Etotal[iE][0]->GetNbinsX()+1)
	      << std::endl;

    std::cout << " -- Summary of energies HCAL: " << std::endl
	      << " ---- SimHits: entries " << p_Etotal[iE][1]->GetEntries() 
	      << " mean " << p_Etotal[iE][1]->GetMean() 
	      << " rms " << p_Etotal[iE][1]->GetRMS() 
	      << " underflows " << p_Etotal[iE][1]->GetBinContent(0)
	      << " overflows " << p_Etotal[iE][1]->GetBinContent(p_Etotal[iE][1]->GetNbinsX()+1)
	      << std::endl;

    std::cout << " -- Summary of energies ECAL+HCAL: " << std::endl
	      << " ---- SimHits: entries " << p_Etotal[iE][2]->GetEntries() 
	      << " mean " << p_Etotal[iE][2]->GetMean() 
	      << " rms " << p_Etotal[iE][2]->GetRMS() 
	      << " underflows " << p_Etotal[iE][2]->GetBinContent(0)
	      << " overflows " << p_Etotal[iE][2]->GetBinContent(p_Etotal[iE][2]->GetNbinsX()+1)
	      << std::endl;
    if (!isG4Tree) {
      std::cout << " ---- RecHits: entries " << p_Ereco[iE][0]->GetEntries() 
		<< " mean " << p_Ereco[iE][0]->GetMean() 
		<< " rms " << p_Ereco[iE][0]->GetRMS() 
		<< " underflows " << p_Ereco[iE][0]->GetBinContent(0)
		<< " overflows " << p_Ereco[iE][0]->GetBinContent(p_Ereco[iE][0]->GetNbinsX()+1)
		<< std::endl;
    }
    std::cout << " -- Summary of hits: " << std::endl
	      << " ---- SimHits: entries " << p_nSimHits[iE]->GetEntries() 
	      << " mean " << p_nSimHits[iE]->GetMean() 
	      << " rms " << p_nSimHits[iE]->GetRMS() 
	      << " underflows " << p_nSimHits[iE]->GetBinContent(0)
	      << " overflows " << p_nSimHits[iE]->GetBinContent(p_nSimHits[iE]->GetNbinsX()+1)
	      << std::endl;
    if (!isG4Tree) {
      std::cout << " ---- RecHits: entries " << p_nRecHits[iE]->GetEntries() 
		<< " mean " << p_nRecHits[iE]->GetMean() 
		<< " rms " << p_nRecHits[iE]->GetRMS() 
		<< " underflows " << p_nRecHits[iE]->GetBinContent(0)
		<< " overflows " << p_nRecHits[iE]->GetBinContent(p_nRecHits[iE]->GetNbinsX()+1)
		<< std::endl;
    }

    // //create summary control plots
    // mycE->cd(iE+1);
    // p_Etotal[iE]->Rebin(10);
    // p_Etotal[iE]->Draw();
    // TLatex lat;
    // char buf[500];
    // sprintf(buf,"Egen = %d GeV",genEn[iE]);
    // lat.DrawLatex(p_Etotal[iE]->GetMean(),p_Etotal[iE]->GetMaximum()*1.1,buf);

  }//loop on energies

 if (!isG4Tree) {
   outputFile->cd();
   for (unsigned iO(0); iO<nOcc; ++iO){
     for (unsigned iL(0); iL<nLayers; ++iL){
       p_occupancy[iO][iL]->Write();
     }
   }
 }

  std::cout << "\\hline \n";
  std::cout << " Energy (GeV) & Eecal (RMS) & Ehcal (RMS) & Etot (RMS) "
    // << " & nSimHits (RMS) "
	    << "\\\\ \n"; 
  for (unsigned iE(0); iE<nGenEn; ++iE){//loop on energies
    if (!p_Etotal[iE][0] || !p_Etotal[iE][1] || !p_Etotal[iE][2]) continue;
    std::cout << genEn[iE] << " & "
	      << p_Etotal[iE][0]->GetMean() << " (" << p_Etotal[iE][0]->GetRMS() 
	      << ") & "
	      << p_Etotal[iE][1]->GetMean() << " (" << p_Etotal[iE][1]->GetRMS() 
	      << ") & "
	      << p_Etotal[iE][2]->GetMean() << " (" << p_Etotal[iE][2]->GetRMS() 
        //<< ") & "
      //<< p_nSimHits[iE]->GetMean() << " (" << p_nSimHits[iE]->GetRMS() 
	      << ") \\\\ \n ";
  }
  std::cout << "\\hline \n";



  // std::ostringstream saveName;
  // saveName.str("");
  // saveName << plotDir << "/ControlG4Etotal";
  // mycE->Update();
  // mycE->Print((saveName.str()+".png").c_str());
  // mycE->Print((saveName.str()+".pdf").c_str());

  // for (unsigned iE(0); iE<nGenEn; ++iE){//loop on energies
  //   p_nSimHits[iE]->Draw();
  //   p_nSimHits[iE]->Rebin(10);
  //   TLatex lat;
  //   char buf[500];
  //   sprintf(buf,"Egen = %d GeV",genEn[iE]);
  //   lat.DrawLatex(p_nSimHits[iE]->GetMean(),p_nSimHits[iE]->GetMaximum()*1.1,buf);
    
  // }//loop on energies

  // saveName.str("");
  // saveName << plotDir << "/ControlG4nSimHits";
  // mycE->Update();
  // mycE->Print((saveName.str()+".png").c_str());
  // mycE->Print((saveName.str()+".pdf").c_str());

  
  outputFile->Write();

  return 0;


}//main
