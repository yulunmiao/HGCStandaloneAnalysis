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

double xLayer(const unsigned layer) {
  double w1b = 1.6+3+0.1;
  double w2b = 3.3+3+0.1;
  double w3b = 5.6+3+0.1;
  double wa = 0.1+4;
  double W1 = 10*(w1b+wa);
  double W2 = 10*(w2b+wa);
  double W3 = 10*(w3b+wa);

  double totalWidth = 10*(w1b+w2b+w3b)+30*wa;
  unsigned il = layer%10+1;

  double offset = layer/10==0 ? (il*w1b + (il-1)*wa) :
    ( layer/10==1 ? W1 + il*w2b + (il-1)*wa :
      W1 + W2 + il*w3b + (il-1)*wa );

  return -totalWidth/2. + offset; //in mm
  
}

double getBinWidth(const double etaMin,const double etaMax, std::string etaStr, const double xLaymm){
  double eta0 = 0;
  if (etaStr.find("eta20")!=etaStr.npos) eta0 = 2.0;
  else if (etaStr.find("eta25")!=etaStr.npos) eta0 = 2.5;
  else if (etaStr.find("eta30")!=etaStr.npos) eta0 = 3.0;
  else if (etaStr.find("eta35")!=etaStr.npos) eta0 = 3.5;
  else {
     std::cout << " ERROR: new eta point ! Please implement value to convert y position... Exiting..." << std::endl;
    exit(1);
  }

  double xDist = 310+xLaymm/10.;//distance in cm from beam spot to layer.
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

double getEta(double ypos, std::string etaStr, const double xLaymm){
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
  double xDist = 3100+xLaymm;//distance from beam spot to layer in mm
  double tanTheta = (y0+ypos)/xDist;
  double theta = atan(tanTheta);
  double eta = -log(tan(theta/2.));
  //std::cout << " - etaBin " << etaStr << " y0 = " << y0 << " y = " << ypos << " eta = " << eta << std::endl;

  return eta;

}

double getWeight(unsigned layer,std::string aVersion){
  if (layer<10) return 1;
  if (aVersion.find("20")!= aVersion.npos || 
      aVersion.find("21") != aVersion.npos ||
      aVersion.find("22") != aVersion.npos ||
      aVersion.find("23") != aVersion.npos){
    if (layer < 20) return 1.43;
    //return 0.8/0.5;
    else if (layer<30) return 2.15;
    //return 1.2/0.5;
    else return 1;
  }
  else {
    if (layer < 20) return 0.8/0.4;
    else if (layer < 30) return 1.2/0.4;
    else return 1;
  }
}


int main(int argc, char** argv){//main  

  if (argc < 4) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <full path to input file>"
	      << " <file name (PFcal.root, or DigiPFcal.root)>" 
	      << " <optional: debug (default=0)>"
	      << " <optional: overlayPU (default=0)>"
	      << " <optional: saveEventByEvent (default=0)>"
	      << " <optional: doChargedPionsOnly (default=0)>"
	      << " <optional: doNeutralPionsOnly (default=0)>"
	      << std::endl;
    return 1;
  }

  bool applyWeight = true;
  bool saveTxtOutput = true;
  bool doFiducialCuts = false;

  const unsigned nLayers = 30;
  const unsigned nEcalLayers = 30;
  const double signalRegionInX=20;
  const double Emip = 0.0548;//MeV

  const unsigned nOcc = 4;
  const unsigned occThreshold[nOcc] = {1,5,10,20};

  const unsigned nSmear = 10;
  const double smearFrac[nSmear] = {0,0.01,0.02,0.03,0.04,0.05,0.07,0.10,0.15,0.20};

  bool saveEventByEvent = false;
  bool overlayPU = false;

  TRandom3 lRndm(0);

  bool saveLLR = false;
  double etaMinLLR = 3.8;
  double etaMaxLLR = 3.9;
  double llrBinWidth[nLayers];


  std::cout << " -- N layers = " << nLayers << std::endl;

  unsigned genEn[]={5,10,20,25,50,75,100,125,150,175,200,250,300,500};
  //unsigned genEn[]={100};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  
  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string fileName = argv[3];//"PFcal.root";
  unsigned debug = 0;
  if (argc >4) debug = atoi(argv[4]);
  bool doChargedPions = false;
  bool doNeutralPions = false;
  if (argc>5) overlayPU = atoi(argv[5]);
  if (argc>6) saveEventByEvent = atoi(argv[6]);
  if (argc>7) doChargedPions = atoi(argv[7]);
  if (argc>8) doNeutralPions = atoi(argv[8]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Input file name: " << fileName << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  if (overlayPU) std::cout << " -- Overlaying PU events from Minbias bank" << std::endl;
  if (doChargedPions) std::cout << " -- Using hits from charged pions only" << std::endl;
  if (doNeutralPions) std::cout << " -- Using hits from neutral pions only" << std::endl;

  unsigned nMaxHits = 10000;
  if (saveEventByEvent) std::cout << " -- saving event by event for events with nHits above " << nMaxHits <<std::endl;

  TString plotDir = "PLOTS/";
  std::vector<std::string> pathVec;
  boost::split( pathVec, filePath, boost::is_any_of("/_"));
  std::string pVersion("");
  std::string pScenario("");
  std::string pParticle("");
  std::string pEta("");

  for (unsigned i(0);i<pathVec.size();++i){
    std::string lTmp = pathVec[i];
    if (lTmp.find("version") != lTmp.npos) {
      pVersion = lTmp;
    }
    else if (lTmp.find("scenario") != lTmp.npos) {
      pScenario = lTmp;
    }
    else if (pVersion.size()!=0 && pParticle.size()==0)  {
      pParticle = lTmp;
    }
    else if (lTmp.find("eta") != lTmp.npos) {
      pEta = lTmp;
      break;
    }
    else if (pVersion.size()!=0) {
      break;
    }
  }

  plotDir += pVersion;
  plotDir += "";
  //plotDir += "/";
  plotDir += pScenario;
  //plotDir += "/";
  plotDir += "_";
  plotDir += pParticle;
  if (overlayPU) plotDir += "/PU";
  if (!doNeutralPions && doChargedPions) plotDir += "_pipm";
  else if (!doChargedPions && doNeutralPions) plotDir += "_pi0";
  plotDir += "/";
  plotDir += pEta;
  plotDir += "/";

  const double totalWidth=(pVersion.find("23")==pVersion.npos) ? 200 : 500;

  bool isMB = false;
  if (pParticle.find("Pedro")!=pParticle.npos) isMB=true;

  if (!applyWeight) plotDir += "noWeights/";
  //plotDir += "simpleWeights/";

  std::cout << " -- Output file directory is : " << plotDir << std::endl;

  TFile *outputFile = TFile::Open(plotDir+"/CalibHistos.root","RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << plotDir << "/CalibHistos.root cannot be opened. Please create output directory : " << plotDir << ". Exiting..." << std::endl;
    return 1;
  }


  double xLay[nLayers];

  TH2F *p_xy[nGenEn][nLayers];
  TH3F *p_xyz[nGenEn];
  TH2F *p_xy_evt[nLayers];
  TH3F *p_xyz_evt = 0;
  TH2F *p_recoxy[nGenEn][nLayers];
  TH3F *p_recoxyz[nGenEn];
  TH2F *p_recoxy_evt[nLayers];

  std::cout << " -- Checking layer-eta conversion:" << std::endl;

  for (unsigned iL(0); iL<nLayers; ++iL){
    p_recoxy_evt[iL] = 0;
    p_xy_evt[iL] = 0;
    xLay[iL] = xLayer(iL);
    llrBinWidth[iL] = getBinWidth(etaMinLLR,etaMaxLLR,"eta35",xLay[iL]);
    std::cout << iL << " " << xLay[iL] << " " 
	      << llrBinWidth[iL] << " y=0@eta3 " 
	      << getEta(0,"eta30",xLay[iL]) << std::endl;
  }

  TH3F *p_recoxyz_evt = 0;
  TH2F *p_EvsLayer[nGenEn];
  TH1F *p_Etot[nGenEn][nLayers];
  TH1F *p_Efrac[nGenEn][nLayers];
  TH1F *p_time[nGenEn][nLayers];
  TH1F *p_Etotal[nGenEn];

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
    //input << "/e_" ;
    //input << genEn[iE] << "/" << fileName ;
    input << "_e" ;
    input << genEn[iE] << fileName;
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
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;

    if (isG4Tree){
      //lTree->SetBranchAddress("event",&event);
      lTree->SetBranchAddress("volNb",&volNb);
      lTree->SetBranchAddress("volX0",&volX0);
      lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    }
    else {
      lTree->SetBranchAddress("event",&event);
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
      std::string inStr = "/afs/cern.ch/work/a/amagnan/public/";
      
      if (isG4Tree) inStr += "HGCalEEGeant4/"+pVersion+"/"+pScenario+"/PedroPU/"+pEta+"/PFcal.root";
      else  inStr += "HGCalEEDigi/"+pVersion+"/"+pScenario+"/PedroPU/"+pEta+"/DigiPFcal.root";
      puFile = TFile::Open(inStr.c_str());
      if (!puFile) {
	std::cout << " -- Error, input file " << inStr << " for PU cannot be opened. Exiting..." << std::endl;
	return 1;
      }
      else std::cout << " -- using file " << inStr << " for PU." << std::endl;
      lPuTree = isG4Tree ? (TTree*)puFile->Get("HGCSSTree") : (TTree*)puFile->Get("RecoTree");
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

      std::cout << "- For PU: " << nPuEvts  << " events are available." << std::endl;
      
    }


    //Initialise histos
    //necessary to have overflows ?
    gStyle->SetOptStat(1111111);
    double Etot[nLayers];
    double eLLR[nLayers];
    for (unsigned iL(0);iL<nLayers;++iL){
      Etot[iL] = 0;
      eLLR[iL] = 0;
    }
    double Etotal = 0;
    unsigned nTotal = 0;
    unsigned nTotalSignal = 0;
    double Ereco[nSmear];
    for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
      Ereco[iSmear]= 0;
    }
 

    lName.str("");
    lName << "p_xyz_" << genEn[iE];
    p_xyz[iE] = new TH3F(lName.str().c_str(),";layer;x(mm);y(mm)",nLayers,0,nLayers,200,-250,250,200,-250,250);
    if (!isG4Tree){
      lName.str("");
      lName << "p_recoxyz_" << genEn[iE];
      p_recoxyz[iE] = new TH3F(lName.str().c_str(),";layer;x(mm);y(mm)",nLayers,0,nLayers,50,-250,250,50,-250,250);
    }
    else p_recoxyz[iE] = 0;
    p_Etotal[iE] = 0;
    p_nSimHits[iE] = 0;
    p_EvsLayer[iE] = 0;
    p_nRecHits[iE] = 0;
    for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
      p_Ereco[iE][iSmear] = 0;
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

      if (debug){
	if (isG4Tree) {
	  std::cout << "... Processing layer " << volNb << " with " << (*simhitvec).size() << " simhits, weight = " << getWeight(volNb,pVersion) << std::endl;
	  Etot[static_cast<unsigned>(volNb)] = 0;
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
	p_xyz_evt = new TH3F("p_xyz",";layer;x(mm);y(mm)",nLayers,0,nLayers,200,-250,250,200,-250,250);
	if (!isG4Tree) {
	  if (p_recoxyz_evt) p_recoxyz_evt->Delete();
	  p_recoxyz_evt = new TH3F("p_recoxyz",";layer;x(mm);y(mm)",nLayers,0,nLayers,50,-250,250,50,-250,250);
	}
	
	for (unsigned iL(0); iL<nLayers; ++iL){
	  lName.str("");
	  lName << "p_xy_" << iL;
	  if (p_xy_evt[iL]) p_xy_evt[iL]->Delete();
	  p_xy_evt[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",200,-250,250,200,-250,250);
	  if (!isG4Tree){
	    lName.str("");
	    lName << "p_recoxy_" << iL;
	    if (p_recoxy_evt[iL]) p_recoxy_evt[iL]->Delete();
	    p_recoxy_evt[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",50,-250,250,50,-250,250);

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
	    if ((lHit.nHadrons() == 0 && lHit.nMuons() == 0)) continue;
	  }
	  else if (doNeutralPions){
	    if (lHit.nGammas()==0 && lHit.nElectrons()==0 ) continue;
	  }
	  
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  //double posx = lHit.get_y();
	  //double posy = lHit.get_x();
	  double energy = lHit.energy();
	  if (debug>1) {
	    std::cout << " --  SimHit " << iH << "/" << (*pusimhitvec).size() << " --" << std::endl
		      << " --  position x,y " << posx << "," << posy << std::endl;
	    lHit.Print(std::cout);
	  }
	  double weightedE = energy;
	  if (applyWeight) weightedE *= getWeight(layer,pVersion);
	  //reject lower region for boundary effects
	  if (doFiducialCuts){
	    if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
		 (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	    //region left and right also for boundary effects
	    if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
		 (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	  }
	  p_xy[iE][layer]->Fill(posx,posy,energy);
	  p_xyz[iE]->Fill(layer,posx,posy,energy);

	  if (saveEventByEvent){
	    p_xy_evt[layer]->Fill(posx,posy,energy);
	    p_xyz_evt->Fill(layer,posx,posy,energy);
	  }
	  p_time[iE][layer]->Fill(lHit.time());
	  //restrict in y and x to have just the influence of the eta bin wanted
	  //(at high eta, surface of the cone in x is smaller than detector size)
	  if (fabs(posx)<signalRegionInX/2.){
	    double lEtaTmp = getEta(posy,pEta,xLay[layer]);
	    p_Edensity[iE][layer]->Fill(lEtaTmp,energy);
	  }
	  Etot[layer] += weightedE;
	  if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << weightedE << " " << Etot[layer];
	  Etotal += weightedE;
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
	    p_recoxyz[iE]->Fill(layer,posx,posy,energy);
	    for (unsigned iO(0); iO<nOcc; ++iO){
	      if (energy > occThreshold[iO]) p_occupancy[iO][layer]->Fill(posx,posy);
	    }
	    double weightedE = energy;
	    if (applyWeight)  weightedE *= getWeight(layer,pVersion);
	    for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	      double smearFactor = lRndm.Gaus(0,smearFrac[iSmear]);
	      Ereco[iSmear] += weightedE*(1+smearFactor);
	    }
	      
	    if (saveEventByEvent){
	      p_recoxy_evt[layer]->Fill(posx,posy,energy);
	      p_recoxyz_evt->Fill(layer,posx,posy,energy);
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
	  if ((lHit.nHadrons() == 0 && lHit.nMuons() == 0)) continue;
	}
	else if (doNeutralPions){
	  if (lHit.nGammas()==0 && lHit.nElectrons()==0 ) continue;
	}

	double posx = lHit.get_x();
	double posy = lHit.get_y();
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
	if (applyWeight) weightedE *= getWeight(layer,pVersion);

	if (doFiducialCuts){
	  //reject lower region for boundary effects
	  if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
	       (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	  //region left and right also for boundary effects
	  if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
	       (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	}

	p_xy[iE][layer]->Fill(posx,posy,energy);
	p_xyz[iE]->Fill(layer,posx,posy,energy);
	if (saveEventByEvent){
	  p_xy_evt[layer]->Fill(posx,posy,energy);
	  p_xyz_evt->Fill(layer,posx,posy,energy);
	}
	p_time[iE][layer]->Fill(lHit.time());
	if (fabs(posx)<signalRegionInX/2.){
	  double lEtaTmp = getEta(posy,pEta,xLay[layer]);
	  p_Edensity[iE][layer]->Fill(lEtaTmp,energy);
	  if (lEtaTmp>=etaMinLLR && lEtaTmp<etaMaxLLR) eLLR[layer]+=energy;
	}
	Etot[layer] += weightedE;
	if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << weightedE << " " << Etot[layer];
	Etotal += weightedE;
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
	  if (applyWeight)  weightedE *= getWeight(layer,pVersion);
	  for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	    double smearFactor = lRndm.Gaus(0,smearFrac[iSmear]);
	    Ereco[iSmear] += weightedE*(1+smearFactor);
	  }

	  p_recoxy[iE][layer]->Fill(posx,posy,energy);
	  p_recoxyz[iE]->Fill(layer,posx,posy,energy);
	  if (saveEventByEvent){
	    p_recoxy_evt[layer]->Fill(posx,posy,energy);
	    p_recoxyz_evt->Fill(layer,posx,posy,energy);
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

	double Etmp = 0;

	if (saveTxtOutput) optim << E1 << " " << E2 << " " << E3 << std::endl;
	E1=0;
	E2=0;
	E3=0;

	if (firstEvent){
	  for (unsigned iL(0);iL<nLayers;++iL){
	    if (Etot[iL]>Emax) Emax = Etot[iL];
	    lName.str("");
	    lName << "p_Etot_" << genEn[iE] << "_" << iL;
	    p_Etot[iE][iL] = new TH1F(lName.str().c_str(),";G4 Etot (MeV)",1000,Etot[iL]*0.1,3*Etot[iL]);
	  }
	  lName.str("");
	  lName << "p_EvsLayer_" << genEn[iE];
	  p_EvsLayer[iE] = new TH2F(lName.str().c_str(),";layer;E (MeV)",nLayers,0,nLayers,1000,0,3*Emax);
	  lName.str("");
	  lName << "p_Etotal_" << genEn[iE];
	  p_Etotal[iE] = new TH1F(lName.str().c_str(),";G4 Etotal (MeV)",1000,0.1*Etotal,2*Etotal);
	  p_Etotal[iE]->StatOverflows();
	  lName.str("");
	  lName << "p_nSimHits_" << genEn[iE];
	  p_nSimHits[iE] = new TH1F(lName.str().c_str(),"; nSimHits",1000,static_cast<unsigned>(nTotal/10.),nTotal*10);
	}

	for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
	  if (debug) std::cout << " -- Layer " << iL << " total E = " << Etot[iL] << std::endl;
	  p_Etot[iE][iL]->Fill(Etot[iL]);
	  p_EvsLayer[iE]->Fill(iL,Etot[iL]);
	  Etmp += Etot[iL];
	  if (Etotal > 0) p_Efrac[iE][iL]->Fill(Etmp/Etotal);
	  else p_Efrac[iE][iL]->Fill(0);
	  Etot[iL] = 0;

	  //////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////
	  //For LLR study
	  //get total density in one eta bin
	  if (saveLLR){// && iL>0 && iL<30){
	    double density = eLLR[iL]/signalRegionInX*1/llrBinWidth[iL];
	    if (isG4Tree) edensity << static_cast<unsigned>(ievt/nLayers);
	    else edensity << ievt;
	    edensity << " " << iL << " " << density << std::endl;
	    eLLR[iL] = 0;
	  }
	  //end-For LLR study
	  //////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////

	}//loop on layers
	p_Etotal[iE]->Fill(Etotal);
	p_nSimHits[iE]->Fill(nTotal);

	if (saveEventByEvent && nTotalSignal > nMaxHits){
	  std::ostringstream evtName;
	  evtName << plotDir << "/CalibHistos_E" << genEn[iE] << "_evt" << lEvt << ".root";
	  TFile *outputEvt = TFile::Open(evtName.str().c_str(),"RECREATE");
	  
	  if (!outputEvt) {
	    std::cout << " -- Error, output file for evt " << lEvt << " cannot be opened. Please create output directory : " << plotDir << ". Exiting..." << std::endl;
	    return 1;
	  }
	  outputEvt->cd();
	  
	  for (unsigned iL(0); iL<nLayers; ++iL){
	    p_xy_evt[iL]->Write();
	    if (!isG4Tree) p_recoxy_evt[iL]->Write();
	  }
	  p_xyz_evt->Write();
	  if (!isG4Tree) p_recoxyz_evt->Write();
	  outputEvt->Close();
	}

	Etotal = 0;
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
    p_Etotal[iE]->Write();
    p_nSimHits[iE]->Write();
    if (!isG4Tree) {
      for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	p_Ereco[iE][iSmear]->Write();
      }
      p_nRecHits[iE]->Write();
      p_recoxyz[iE]->Write();
    }

    std::cout << " -- Summary of energies: " << std::endl
	      << " ---- SimHits: entries " << p_Etotal[iE]->GetEntries() 
	      << " mean " << p_Etotal[iE]->GetMean() 
	      << " rms " << p_Etotal[iE]->GetRMS() 
	      << " underflows " << p_Etotal[iE]->GetBinContent(0)
	      << " overflows " << p_Etotal[iE]->GetBinContent(p_Etotal[iE]->GetNbinsX()+1)
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
  std::cout << " Energy (GeV) & Esim (RMS) "
    // << " & nSimHits (RMS) "
	    << "\\\\ \n"; 
  for (unsigned iE(0); iE<nGenEn; ++iE){//loop on energies
    std::cout << genEn[iE] << " & "
	      << p_Etotal[iE]->GetMean() << " (" << p_Etotal[iE]->GetRMS() 
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
