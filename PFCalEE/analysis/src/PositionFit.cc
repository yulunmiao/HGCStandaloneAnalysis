#include <iomanip>

#include "PositionFit.hh"
#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSDetector.hh"

std::pair<unsigned, std::pair<double,double> > PositionFit::findMajorityValue(std::vector<std::pair<double,double> > & values) const
{
  
  unsigned lTot = values.size();
  if (!lTot) return std::pair<unsigned, std::pair<double,double> >(0,std::pair<double,double>(0,0));
  
  std::sort(values.begin(),values.end());

  unsigned lMajorityCounter = 0;
  std::pair<double,double> lMaj = std::pair<double,double>(0,0);
  
  std::vector<std::pair<double,double> >::iterator lIter = values.begin();
  for ( ; lIter != values.end(); ) {
    //std::cout << lIter->first << " " << lIter->second << std::endl;
    unsigned lCounter = std::count(lIter,values.end(),*lIter);
    if (lCounter > lMajorityCounter) {
      lMajorityCounter = lCounter;
      lMaj = *lIter;
    }
    lIter += lCounter;
  }
    
  //std::cout << " -- Found majority value " << lMaj.first << " " << lMaj.second << " for " << lMajorityCounter << " elements out of " << values.size() << "." << std::endl;

  return std::pair<unsigned, std::pair<double,double> >(lMajorityCounter,lMaj);
  
}

PositionFit::PositionFit(const unsigned nSR,
			 const double & residualMax, 
			 const unsigned nLayers, 
			 const unsigned nSiLayers,
			 const bool applyPuMixFix,
			 const unsigned debug
			 ){
  nSR_ = nSR;
  residualMax_ = residualMax;
  chi2ndfmax_ = 20;
  nLayers_ = nLayers;
  nSiLayers_ = nSiLayers;
  debug_ = debug;
  useMeanPU_ = true;
  fixForPuMixBug_ = applyPuMixFix;

  p_nGenParticles = 0;
  p_hitMeanPuContrib = 0;
  p_hitEventPuContrib = 0;
  p_diffPuContrib = 0;

  p_etavsphi = 0;
  p_etavsphi_max = 0;
  p_etavsphi_truth = 0;
  p_yvsx_max = 0;
  p_yvsx_truth = 0;

  p_residuals_x = 0;
  p_residuals_y = 0;
  p_errorMatrix = 0;
  p_corrMatrix = 0;
  p_chi2[0] = 0;
  p_chi2[1] = 0;

  p_chi2overNDF[0] = 0;
  p_impactX[0] = 0;
  p_impactY[0] = 0;
  p_impactX_residual = 0;
  p_impactY_residual = 0;
  p_tanAngleX[0] = 0;
  p_tanAngleY[0] = 0;
  p_tanAngleX_residual = 0;
  p_tanAngleY_residual = 0;
  p_positionReso = 0;
  p_angularReso = 0;
  p_chi2overNDF[1] = 0;
  p_impactX[1] = 0;
  p_impactY[1] = 0;
  p_tanAngleX[1] = 0;
  p_tanAngleY[1] = 0;

}

void PositionFit::initialise(TFile *outputFile,
			     const std::string outputDir,
			     const std::string outFolder, 
			     const HGCSSGeometryConversion & geomConv, 
			     const HGCSSPUenergy & puDensity){
  outputDir_ = outputDir;
  setOutputFile(outputFile);

  outFolder_ = outFolder;

  geomConv_ = geomConv;
  puDensity_ = puDensity;

  nL_mean_.resize(nLayers_,0);
  mean_[0].resize(nLayers_,0);
  mean_[1].resize(nLayers_,0);

  nL_sigma_.resize(nLayers_,nL_mean_);
  sigma_[0].resize(nLayers_,mean_[0]);
  sigma_[1].resize(nLayers_,mean_[0]);
  for (unsigned iL(0);iL<nLayers_;++iL){
    nL_sigma_[iL].resize(nLayers_,0);
    sigma_[0][iL].resize(nLayers_,0);
    sigma_[1][iL].resize(nLayers_,0);
  }

}

void PositionFit::initialisePositionHistograms(){
  //for yvsx plots
  unsigned nX=2*1700/10-2;
  double minX=-1.*nX*5-5,maxX=nX*5+5;
  double minY=minX,maxY=maxX;
  nX += 1;
  unsigned nY = nX;

  outputFile_->cd(outputDir_.c_str());

  p_nGenParticles = new TH1F("p_nGenParticles",";nGenParticles",10,0,10);

  p_genxy.resize(nLayers_,0);
  p_recoxy.resize(nLayers_,0);
  std::ostringstream lName;
  for (unsigned iL(0); iL<nLayers_; ++iL){
    lName.str("");
    lName << "p_genxy_" << iL;
    p_genxy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
			   nX*10,minX,maxX,
			   nY*10,minY,maxY);
    lName.str("");
    lName << "p_recoxy_" << iL;
    p_recoxy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
   			    nX,minX,maxX,
   			    nY,minY,maxY);
  }

  p_hitMeanPuContrib = new TH2F("p_hitMeanPuContrib",";layer;E_{PU} (MIPs) from mean;hits",nLayers_,0,nLayers_,1000,0,50);
  p_hitMeanPuContrib->StatOverflows();

  p_hitEventPuContrib = new TH2F("p_hitEventPuContrib",";layer;E_{PU} (MIPs) from RC;hits",nLayers_,0,nLayers_,1000,0,50);
  p_hitEventPuContrib->StatOverflows();

  p_diffPuContrib = new TH1F("p_diffPuContrib",";E_{PU}^{RC}-E_{PU}^{avg} (MIPs) in SR2;events",1000,-2000,2000);
  p_diffPuContrib->StatOverflows();

  p_etavsphi = new TH2F("p_etavsphi",";#phi_{hit};#eta_{hit};n_{events}",
			900,-3.1416,3.1416,
			250,1.4,3.0);

  p_etavsphi_max = new TH2F("p_etavsphi_max",";#phi_{max};#eta_{max};n_{events}",
			900,-3.1416,3.1416,
			250,1.4,3.0);
  p_etavsphi_truth = new TH2F("p_etavsphi_truth",";#phi_{gen};#eta_{gen};n_{events}",
			900,-3.1416,3.1416,
			250,1.4,3.0);

  p_yvsx_max = new TH2F("p_yvsx_max",";x_{max};y_{max};n_{events}",
			nX,minX,maxX,nY,minY,maxY);
  p_yvsx_truth = new TH2F("p_yvsx_truth",";x_{gen};y_{gen};n_{events}",
			  nX,minX,maxX,nY,minY,maxY);


  p_residuals_x = new TH1F("p_residuals_x",";xreco-xtruth (mm)",1000,-50,50);
  p_residuals_y = new TH1F("p_residuals_y",";yreco-ytruth (mm)",1000,-50,50);
  p_residuals_x->StatOverflows();
  p_residuals_y->StatOverflows();

}

void PositionFit::initialiseFitHistograms(){

  outputFile_->cd(outputDir_.c_str());

  //check if already defined
  if (!p_chi2[0]){

    p_recoXvsLayer = new TH2F("p_recoXvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,200,-100,100);
    p_recoYvsLayer = new TH2F("p_recoYvsLayer",";layer;weighted y (mm);n_{events}",nLayers_,0,nLayers_,1200,300,1500);
    p_recoZvsLayer = new TH2F("p_recoZvsLayer",";layer;avg z (mm);n_{events}",nLayers_,0,nLayers_,3000,3170,3470);
    p_truthXvsLayer = new TH2F("p_truthXvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,200,-100,100);
    p_truthYvsLayer = new TH2F("p_truthYvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,1200,300,1500);
    p_fitXvsLayer = new TH2F("p_fitXvsLayer",";layer;fit x (mm);n_{events}",nLayers_,0,nLayers_,200,-100,100);
    p_fitYvsLayer = new TH2F("p_fitYvsLayer",";layer;fit y (mm);n_{events}",nLayers_,0,nLayers_,1200,300,1500);
    p_nLayersFit = new TH1F("p_nLayersFit",";#layers in fit;n_{events}",31,-0.5,30.5);

    p_chi2[0] = new TH1F("p_chi2",";#chi^{2};n_{events}",1000,0,5000);
    p_chi2[1] = new TH1F("p_chi2_truth",";#chi^{2};n_{events}",1000,0,5000);
    p_chi2overNDF[0] = new TH1F("p_chi2overNDF",";#chi^{2}/NDF;n_{events}",1000,0,500);
    p_chi2overNDF[1] = new TH1F("p_chi2overNDF_truth",";#chi^{2}/NDF;n_{events}",1000,0,500);
    for (unsigned rt(0); rt<2;++rt){
      p_chi2[rt]->StatOverflows();
      p_chi2overNDF[rt]->StatOverflows();
    }
    
    p_impactX[0] = new TH1F("p_impactX",";x front face impact (mm);n_{events}",500,-100,100);
    p_impactX[1] = new TH1F("p_impactX_truth",";x front face impact (mm);n_{events}",500,-100,100);
    p_impactY[0] = new TH1F("p_impactY",";y front face impact (mm);n_{events}",1200,300,1500);
    p_impactY[1] = new TH1F("p_impactY_truth",";y front face impact (mm);n_{events}",1200,300,1500);
    p_tanAngleX[0] = new TH1F("p_tanAngleX",";x direction tanAngle (rad);n_{events}",500,-1,1);
    p_tanAngleX[1] = new TH1F("p_tanAngleX_truth",";x direction tanAngle (rad);n_{events}",500,-1,1);
    p_tanAngleY[0] = new TH1F("p_tanAngleY",";y direction tanAngle (rad);n_{events}",500,-1,1);
    p_tanAngleY[1] = new TH1F("p_tanAngleY_truth",";y direction tanAngle (rad);n_{events}",500,-1,1);

    p_impactX_residual = new TH1F("p_impactX_residual",";residual x front face impact (mm);n_{events}",200,-10,10);
    p_tanAngleX_residual = new TH1F("p_tanAngleX_residual",";residual x direction tanAngle (rad);n_{events}",200,-0.1,0.1);
    p_impactY_residual = new TH1F("p_impactY_residual",";residual y front face impact (mm);n_{events}",200,-10,10);
    p_tanAngleY_residual = new TH1F("p_tanAngleY_residual",";residual y direction tanAngle (rad);n_{events}",200,-0.1,0.1);

    
    p_positionReso = new TH1F("p_positionReso",";#sigma_{x,y} (mm);n_{events}",500,0,50);
    p_angularReso = new TH1F("p_angularReso",";#sigma_{#theta} (rad);n_{events}",500,0,1);
  }
}

bool PositionFit::getGlobalMaximum(const unsigned ievt, const unsigned nVtx, std::vector<HGCSSRecoHit> *rechitvec,double & aPhimax,double & aEtamax){

  bool oneresult = true;
  TH2F *letavsphi = new TH2F("letavsphi",";#phi;#eta;hits",
			     900,-3.1416,3.1416,
			     250,1.4,3.0);
  
  /*
  TH2F *letavsphi_layer[nLayers_];
  if (nVtx>0){
    std::ostringstream lName;
    for (unsigned iL(0); iL<nLayers_; ++iL){
      lName.str("");
      lName << "letavsphi_layer" << iL;
      letavsphi_layer[iL] = new TH2F(lName.str().c_str(),";#phi;#eta;hits",
				     900,-3.1416,3.1416,
				     250,1.4,3.0);
    }
  }
  */

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;
    double posz = lHit.get_z();
    //double radius = sqrt(posx*posx+posy*posy);
    double energy = lHit.energy();
    unsigned layer = lHit.layer();
    //if (energy>1) std::cout << "Hit " << layer << " " << posx << " " << posy << " " << posz << " " << energy << std::endl;
    
    if (debug_>1) {
      std::cout << " --  RecHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		<< " --  position x,y " << posx << "," << posy << std::endl;
      lHit.Print(std::cout);
    }
    
    ROOT::Math::XYZVector pos(posx,posy,posz);
    letavsphi->Fill(pos.phi(),pos.eta(),energy);

    //if (nVtx>0){
    //letavsphi_layer[layer]->Fill(pos.phi(),pos.eta(),energy);
    //}
    p_etavsphi->Fill(pos.phi(),pos.eta(),energy);
  }//loop on hits

  if (debug_)  std::cout << std::endl;

  //add histograms for all layers but iL
  //for (unsigned iL(0); iL<nLayers_; ++iL){
  //letavsphi_layer[iL]->Add(letavsphi,letavsphi_layer[iL],1,-1);
  //}
  
  //get position of maximum E tower
  int maxbin = letavsphi->GetMaximumBin();

  int binx,biny,binz;
  letavsphi->GetBinXYZ(maxbin,binx,biny,binz);
  
  aPhimax =letavsphi->GetXaxis()->GetBinCenter(binx); 
  aEtamax =letavsphi->GetYaxis()->GetBinCenter(biny); 

  /*
  if (nVtx>0){
    double binxsize =  letavsphi->GetXaxis()->GetBinWidth(binx);
    double binysize =  letavsphi->GetYaxis()->GetBinWidth(biny);
    
    //allow for +/- 1 bin in each direction
    double dRmax = 2*sqrt(pow(binxsize,2)+pow(binysize,2));
    
    //unsigned counter = 0;
    //std::vector<double> layervec;
    //std::vector<double> etavec;
    //std::vector<double> phivec;
    std::vector<std::pair<double,double> > etaphivec;

    for (unsigned iL(0); iL<nLayers_; ++iL){
      int maxbin_layer = letavsphi_layer[iL]->GetMaximumBin();
      int binxl,binyl,binzl;
      letavsphi_layer[iL]->GetBinXYZ(maxbin_layer,binxl,binyl,binzl);
      double leta = letavsphi_layer[iL]->GetYaxis()->GetBinCenter(binyl);
      double lphi = letavsphi_layer[iL]->GetXaxis()->GetBinCenter(binxl);
      etaphivec.push_back(std::pair<double,double>(leta,lphi));
      //double deta = leta-aEtamax;
      //double dphi = lphi-aPhimax;
      //if (dphi<-1.*TMath::Pi()) dphi += 2*TMath::Pi();
      //double dR = sqrt(pow(deta,2)+pow(dphi,2));
      //if (dR>dRmax) {
      //layervec.push_back(iL);
      //etavec.push_back(leta);
      //phivec.push_back(lphi);
      //counter++;
      //}
    }
    
    std::pair<unsigned, std::pair<double,double> > lmaj = findMajorityValue(etaphivec);
    double leta = lmaj.second.first;
    double lphi = lmaj.second.second;
    double deta = leta-aEtamax;
    double dphi = lphi-aPhimax;
    if (dphi<-1.*TMath::Pi()) dphi += 2.*TMath::Pi();
    double dR = sqrt(pow(deta,2)+pow(dphi,2));
    
    if (dR>dRmax) {
      std::cout << " -- Warning ! Event " << ievt << " with " << lmaj.first << " majority layers dRmax=" << dRmax << " away, probably from PU." << std::endl
		<< " All layers 0-29 eta-phi = " << aEtamax << " " << aPhimax << std::endl
		<< " Maj value = " << leta << " " << lphi << std::endl;
      oneresult = false;
    }
    // if (counter > 0){
    //std::cout << " -- Warning ! Event " << ievt << " with " << counter << " layers dRmax=" << dRmax << " away, probably from PU." << std::endl
    // 		<< " All layers 0-29 eta-phi = " << aEtamax << " " << aPhimax << std::endl;
    //   for (unsigned iC(0); iC<counter;++iC){
    // 	std::cout << " Removing layer " << layervec[iC] << " eta-phi = " << etavec[iC] << " " << phivec[iC] << std::endl;
    //   }
  //}
  }//if nvtx>0
  */

  if (debug_) std::cout << " MaxE cell eta,phi = " << aEtamax << " " << aPhimax << std::endl;

  letavsphi->Delete();
  /*if (nVtx>0) {
    for (unsigned iL(0); iL<nLayers_; ++iL){
      letavsphi_layer[iL]->Delete();
    }
    }*/

  return oneresult;
}

bool PositionFit::getTruthPosition(std::vector<HGCSSGenParticle> *genvec,std::vector<ROOT::Math::XYPoint> & truthPos){

  bool found = false;
  for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
    //if ((*genvec).size()!= 1) (*genvec)[iP].Print(std::cout);
    if ((*genvec)[iP].trackID()==1){
      found = true;
      double x0 = (*genvec)[iP].x();
      double y0 = (*genvec)[iP].y();
      double z0 = (*genvec)[iP].z();
      //double p = sqrt(pow((*genvec)[iP].px(),2)+pow((*genvec)[iP].py(),2)+pow((*genvec)[iP].pz(),2));
      //std::cout << "init : " << x0 << " " << y0 << " " << z0 << std::endl;
      //fill layers by propagating with momentum
      //ROOT::Math::XYZVector unit((*genvec)[iP].px()/p,(*genvec)[iP].py()/p,(*genvec)[iP].pz()/p);
      
      //std::cout << " Gen particle eta,phi = " << unit.eta() << " " << unit.phi() << std::endl;
      p_etavsphi_truth->Fill((*genvec)[iP].phi(),(*genvec)[iP].eta());
      //std::cout << " Truth pos eta-phi = " << (*genvec)[iP].eta() << " " << (*genvec)[iP].phi() << std::endl;
      
      for (unsigned iL(0); iL<nLayers_; ++iL){
	//double xy = (avgZ_[iL]-z0)/sinh(unit.eta());
	//double x = xy*cos(unit.phi())+x0;
	//double y = xy*sin(unit.phi())+y0;
	double x = x0+(avgZ_[iL]-z0)*(*genvec)[iP].px()/(*genvec)[iP].pz();
	double y = y0+(avgZ_[iL]-z0)*(*genvec)[iP].py()/(*genvec)[iP].pz();
	//std::cout << "Lay " << iL << ": " << x << " " << y << " " << avgZ_[iL] << std::endl;
	p_genxy[iL]->Fill(x,y,1);
	truthPos[iL] = ROOT::Math::XYPoint(x,y);
      }
      
    }
  }//loop on gen particles

  if (!found){
    std::cout << " - Info: no G4trackID=1 found, already converted..." << std::endl;
    //std::cout << " -- Number of genparticles: " << (*genvec).size() << std::endl;
    //for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles 
    //std::cout << " --- particle " << iP << std::endl;
    //(*genvec)[iP].Print(std::cout);
    //}
  }
  else {
    p_nGenParticles->Fill((*genvec).size());
  }
  return found;

}

bool PositionFit::getZpositions(){
  std::ifstream fin;
  std::ostringstream finname;
  finname << outFolder_ << "/zPositions.dat";
  fin.open(finname.str());
  if (!fin.is_open()){
    std::cout << " Cannot open input file " << finname.str() << "! Refilling now..." << std::endl;
    return false;
  }
  
  std::cout << " Reading z position per layer from input file " << finname.str() << std::endl;

  std::vector<unsigned> layerId;
  std::vector<double> posz;
  layerId.reserve(nLayers_);
  posz.reserve(nLayers_);

  while (!fin.eof()){
    unsigned l=nLayers_;
    double z=0;
    fin>>l>>z;
    if (l<nLayers_){
      avgZ_.push_back(z);
      std::cout << " Layer " << l << ", z = " << z << std::endl;
    }
  }
  
  if (avgZ_.size() != nLayers_) {
    std::cout << " -- Warning! Problem in extracting z positions, did not find one value per layer. Please check input file: " << finname.str() << std::endl
	      << " Proceeding to refilling them from simtree." << std::endl;
    return false;
  }

  fin.close();
  return true;
 
}

void PositionFit::getZpositions(TTree *aSimTree,
				const unsigned nEvts){
  

  std::vector<HGCSSSimHit> * simhitvec = 0;
  aSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);

  std::ofstream fout;
  std::ostringstream foutname;
  foutname << outFolder_ << "/zPositions.dat";
  fout.open(foutname.str());
  if (!fout.is_open()){
    std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }
  
  std::cout << "--- Filling z positions:" << std::endl
	    << "- Processing = " << nEvts  << " events out of " << aSimTree->GetEntries() << std::endl;

  avgZ_.resize(nLayers_,0);


  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    aSimTree->GetEntry(ievt);
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on rechits
      HGCSSSimHit lHit = (*simhitvec)[iH];
      unsigned layer = lHit.layer();
      if (layer >= nLayers_) {
	continue;
      }

      //discard some si layers...
      if (lHit.silayer() >= nSiLayers_) continue; 

      double posz = lHit.get_z();
      //get z position of hits
      if (avgZ_[layer]<posz) avgZ_[layer]=posz;
    }
  }

  std::cout << " --- Z positions of layers: " << std::endl;
  for (unsigned iL(0); iL<nLayers_;++iL){
    std::cout << " Layer " << iL << ", z = " << avgZ_[iL] << std::endl;
    fout << iL << " " << avgZ_[iL] << std::endl;
  }

  fout.close();
  
}

void PositionFit::getInitialPositions(TTree *aSimTree, 
				      TTree *aRecTree,
				      const unsigned nEvts){

  initialisePositionHistograms();
  //HGCSSDetector & myDetector = theDetector();

  //////////////////////////////////////////////////
  ///////// Event loop /////////////////////////////
  //////////////////////////////////////////////////

  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;

  aSimTree->SetBranchAddress("HGCSSEvent",&event);
  aSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  aSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  aSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  
  aRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (aRecTree->GetBranch("nPuVtx")) aRecTree->SetBranchAddress("nPuVtx",&nPuVtx);

  std::cout << "- Processing = " << nEvts  << " events out of " << aSimTree->GetEntries() << std::endl;

  //bool firstEvent = true;

  unsigned nConvertedPhotons = 0;
  unsigned nMultipleMax = 0;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    aSimTree->GetEntry(ievt);

    aRecTree->GetEntry(ievt);

    if (debug_) std::cout << " nPuVtx = " << nPuVtx << std::endl;

    if (debug_){
      std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< " gen " << (*genvec).size() << std::endl;
    }

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    //////// output files to save position for chi2 fit //////////
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    double phimax = 0;
    double etamax = 0;
    bool oneresult = getGlobalMaximum(ievt,nPuVtx,rechitvec,phimax,etamax);
    if (!oneresult) nMultipleMax++;
    p_etavsphi_max->Fill(phimax,etamax);

    std::vector<ROOT::Math::XYPoint> truthPos;
    truthPos.resize(nLayers_,ROOT::Math::XYPoint(0,0));
    bool found = getTruthPosition(genvec,truthPos);

    if (!found) {
      nConvertedPhotons++;
      continue;
    }
    
    p_yvsx_truth->Fill(truthPos[10].X(),truthPos[10].Y());

    std::vector<double> xmax;
    xmax.resize(nLayers_,0);
    std::vector<double> ymax;
    ymax.resize(nLayers_,0);
    getMaximumCell(rechitvec,phimax,etamax,xmax,ymax);
    p_yvsx_max->Fill(xmax[10],ymax[10]);

    //get PU contrib from elsewhere in the event
    //loop over phi with same etamax
    //take average per layer: not all 9 cells of 3*3 area have hits...
    std::vector<double> puE;
    puE.resize(nLayers_,0);
    if (nPuVtx>0){
      unsigned nRandomCones = 50;
      double phistep = TMath::Pi()/nRandomCones;
      if (debug_) std::cout << "--- etamax = " << etamax << " phimax=" << phimax << " phistep = " << phistep << std::endl;
      for (unsigned ipm(0);ipm<nRandomCones;++ipm){
	std::vector<double> xmaxrc;
	xmaxrc.resize(nLayers_,0);
	std::vector<double> ymaxrc;
	ymaxrc.resize(nLayers_,0);
	double phirc = phimax-TMath::Pi();
	if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
	if (ipm%2==0) phirc += ipm/2*phistep+phistep/2.;
	else  phirc = phirc - ipm/2*phistep-phistep/2.;
	if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
	//take from geom to not be biased by hit having PU, because
	//not from geom means find cell with a hit closest to maxpos...
	getMaximumCellFromGeom(phirc,etamax,xmaxrc,ymaxrc);
	if (debug_) std::cout << "rc #" << ipm << " phirc=" << phirc << " xmax[10]=" << xmaxrc[10] << " ymax[10]=" << ymaxrc[10] << " r=" << sqrt(pow(xmaxrc[10],2)+pow(ymaxrc[10],2)) << std::endl;
	getPuContribution(rechitvec,xmaxrc,ymaxrc,puE);
      }
      
      //normalise to one cell: must count cells with 0 hit !
      //use cell size at etamax...
      for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
	if (debug_) std::cout << "layer " << iL ;
	unsigned nCells = nRandomCones*nSR_*geomConv_.cellSize()/geomConv_.cellSize(iL,etamax);
	puE[iL] = puE[iL]/nCells;
	if (debug_) std::cout << " Epu=" << puE[iL] << std::endl;	
      }

    }//if PU

    std::vector<ROOT::Math::XYPoint> recoPos;
    recoPos.resize(nLayers_,ROOT::Math::XYPoint(0,0));
    std::vector<unsigned> nHits;
    nHits.resize(nLayers_,0);

    //get energy-weighted position around maximum
    getEnergyWeightedPosition(rechitvec,nPuVtx,xmax,ymax,recoPos,nHits,puE);

    std::ofstream fout;
    std::ostringstream foutname;
    foutname << outFolder_ << "/initialPos_evt" << ievt << ".dat";
    fout.open(foutname.str());
    if (!fout.is_open()){
      std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
      exit(1);
    }

    if (debug_) std::cout << " Summary of reco and truth positions:" << std::endl;
    for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
      if (debug_) std::cout << iL << " nHits=" << nHits[iL] << " Max=(" << xmax[iL] << "," << ymax[iL] << ")\t Reco=(" << recoPos[iL].X() << "," << recoPos[iL].Y() << ")\t Truth=(" << truthPos[iL].X() << "," << truthPos[iL].Y() << ")" << std::endl;
      fout << iL << " " << recoPos[iL].X() << " " << recoPos[iL].Y() << " " << truthPos[iL].X() << " " << truthPos[iL].Y() << std::endl;
    }

    fout.close();

    fillErrorMatrix(recoPos,truthPos,nHits);

    //firstEvent = false;
  }//loop on entries

  std::cout << " -- Number of converted photons: " << nConvertedPhotons << std::endl;
  std::cout << " -- Number of events with multiple choice for maximum eta-phi bin: " << nMultipleMax << std::endl;

}

void PositionFit::getMaximumCellFromGeom(const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax){

  for (unsigned iL(0); iL<nLayers_;++iL){
    double theta = 2*atan(exp(-etamax));
    double rho = avgZ_[iL]/cos(theta);
    xmax[iL] = rho*tan(theta)*cos(phimax);
    ymax[iL] = rho*tan(theta)*sin(phimax);
  }//loop on layers

}

void PositionFit::getMaximumCell(std::vector<HGCSSRecoHit> *rechitvec,const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax){
  
  std::vector<double> dRmin;
  dRmin.resize(nLayers_,10);
  
  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    
    unsigned layer = lHit.layer();
    if (layer >= nLayers_) {
      continue;
    }
    
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;
    double posz = lHit.get_z();
    if (debug_>1) {
      std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		<< " --  position x,y " << posx << "," << posy << std::endl;
      lHit.Print(std::cout);
    }


    ROOT::Math::XYZVector pos(posx,posy,posz);
    double deta = fabs(pos.eta()-etamax);
    double dphi = pos.phi()-phimax;
    if (dphi<-1.*TMath::Pi()) dphi += 2*TMath::Pi();
    double dR = sqrt(pow(deta,2)+pow(dphi,2));
    if (dR<dRmin[layer]) {
      dRmin[layer] = dR;
      xmax[layer] = posx;
      ymax[layer] = posy;
    }
    
    double energy = lHit.energy();
    p_recoxy[layer]->Fill(posx,posy,energy);
    
  }//loop on rechits
    
}

void PositionFit::getEnergyWeightedPosition(std::vector<HGCSSRecoHit> *rechitvec,const unsigned nPU, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<ROOT::Math::XYPoint> & recoPos,std::vector<unsigned> & nHits,std::vector<double> & puE,const bool puSubtracted){
  
  std::vector<double> eSum;
  double eSum_puMean = 0;
  double eSum_puEvt = 0;

  eSum.resize(nLayers_,0);
  double step = geomConv_.cellSize()*nSR_/2.+0.1;//+0.1 to accomodate double precision
  if (debug_) std::cout << "step = " << step << std::endl;

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    double energy = lHit.energy();//in MIP already...
    unsigned layer = lHit.layer();
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;


    if (fabs(posx-xmax[layer]) < step && 
	fabs(posy-ymax[layer]) < step){
      if (puSubtracted) {
	double leta = lHit.eta();
	if (debug_>1) std::cout << " -- Hit " << iH << ", eta=" << leta << ", energy before PU subtraction: " << energy << " after: " ;
	double lCorMean =  puDensity_.getDensity(leta,layer,geomConv_.cellSizeInCm(layer,leta),nPU);
	eSum_puMean += lCorMean;
	p_hitMeanPuContrib->Fill(layer,lCorMean);
	double lCorEvent = puE[layer];
	eSum_puEvt += lCorEvent;
	p_hitEventPuContrib->Fill(layer,lCorEvent);
	double lCor = 0;
	if (useMeanPU_) lCor = lCorMean;
	else lCor = lCorEvent;
	energy = std::max(0.,energy - lCor);
	if (debug_>1) std::cout << energy << std::endl;
      }
      recoPos[layer].SetX(recoPos[layer].X() + posx*energy);
      recoPos[layer].SetY(recoPos[layer].Y() + posy*energy);
      eSum[layer] += energy;
      if (energy>0) nHits[layer]++;
    }
    
  }//loop on rechits

  p_diffPuContrib->Fill(eSum_puEvt-eSum_puMean);

  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    if (nHits[iL]==0) continue;
    recoPos[iL].SetX(recoPos[iL].X()/eSum[iL]);
    recoPos[iL].SetY(recoPos[iL].Y()/eSum[iL]);
  }
  
}

void PositionFit::getPuContribution(std::vector<HGCSSRecoHit> *rechitvec, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<double> & puE){

  double step = geomConv_.cellSize()*nSR_/2.+0.1;

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    double energy = lHit.energy();//in MIP already...
    unsigned layer = lHit.layer();
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;

    //std::cout << "- iH " << iH << " x=" << posx << " xmax=" << xmax[layer] << " y=" << posy << " ymax=" << ymax[layer] << " step " << step << std::endl;
    if (fabs(posx-xmax[layer]) < step && 
	fabs(posy-ymax[layer]) < step){
      if (debug_>1) std::cout << "- iH " << iH 
			     << " x=" << posx 
			     << " xmax=" << xmax[layer] 
			     << " y=" << posy 
			     << " ymax=" << ymax[layer] 
			     << " step " << step 
			     << " --- Pass, layer" << layer 
			     << " E="<< energy << std::endl;
      puE[layer] += energy;
    }

  }//loop on rechits
  //exit(1);
}

void PositionFit::fillErrorMatrix(const std::vector<ROOT::Math::XYPoint> & recoPos,const std::vector<ROOT::Math::XYPoint> & truthPos,const std::vector<unsigned> & nHits){

  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    if (nHits[iL]==0) continue;
    double residual_xi = recoPos[iL].X()-truthPos[iL].X();
    double residual_yi = recoPos[iL].Y()-truthPos[iL].Y();
    p_residuals_x->Fill(residual_xi);
    p_residuals_y->Fill(residual_yi);
    if (fabs(residual_xi)>residualMax_ || fabs(residual_yi)>residualMax_) continue;
    mean_[0][iL] += residual_xi;
    mean_[1][iL] += residual_yi;
    ++nL_mean_[iL];
    for (unsigned jL(0);jL<nLayers_;++jL){//loop on layers
      if (nHits[jL]==0) continue;
      double residual_xj = recoPos[jL].X()-truthPos[jL].X();
      double residual_yj = recoPos[jL].Y()-truthPos[jL].Y();
      if (fabs(residual_xj)>residualMax_ || fabs(residual_yj)>residualMax_) continue;
      double sigma_x = residual_xi*residual_xj;
      double sigma_y = residual_yi*residual_yj;
      sigma_[0][iL][jL] += sigma_x;
      sigma_[1][iL][jL] += sigma_y;
      ++nL_sigma_[iL][jL];
    }//loop on layers
  }//loop on layers

}


void PositionFit::finaliseErrorMatrix(){
  //finalise error matrix

  std::ofstream fmatrix;
  std::ostringstream fmatrixname;
  fmatrixname << outFolder_ << "/errorMatrix.dat";
  fmatrix.open(fmatrixname.str());
  if (!fmatrix.is_open()){
    std::cout << " Cannot open outfile " << fmatrixname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }

  matrix_.ResizeTo(nLayers_,nLayers_);

  //set mean values first
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    mean_[0][iL] = mean_[0][iL]/nL_mean_[iL];
    mean_[1][iL] = mean_[1][iL]/nL_mean_[iL];
  }
  //set sigmas and fill matrix
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    for (unsigned jL(0);jL<nLayers_;++jL){//loop on layers
      sigma_[0][iL][jL] = sigma_[0][iL][jL]/nL_sigma_[iL][jL];
      sigma_[1][iL][jL] = sigma_[1][iL][jL]/nL_sigma_[iL][jL];
      //consider average of both x and y in one matrix
      matrix_[iL][jL] = 0.5*(sigma_[0][iL][jL]-mean_[0][iL]*mean_[0][jL]+
			     sigma_[1][iL][jL]-mean_[1][iL]*mean_[1][jL]);
      //matrix_[jL][iL] = matrix_[iL][jL];
      
      fmatrix << iL << " " << jL << " " << std::setprecision(15) << matrix_[iL][jL] << std::endl;

      //if (iL!=jL){
      //p_matrix->Fill(jL,iL,matrix_[iL][jL]);
      //}
    }
  }

  fmatrix.close();

}

void PositionFit::fillCorrelationMatrix(){
  std::cout << " -- Filling correlation matrix" << std::endl;
  outputFile_->cd(outputDir_.c_str());

  p_errorMatrix = new TH2D("p_errorMatrix",";i;j;M_{ij}",
			   nLayers_,0,nLayers_,
			   nLayers_,0,nLayers_);
  p_corrMatrix = new TH2D("p_corrMatrix",";i;j;M_{ij}",
			  nLayers_,0,nLayers_,
			  nLayers_,0,nLayers_);

  corrMatrix_.ResizeTo(nLayers_,nLayers_);
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    for (unsigned jL(0);jL<nLayers_;++jL){//loop on layers
      p_errorMatrix->Fill(iL,jL,matrix_[iL][jL]);
      if (matrix_[iL][iL]!=0 && matrix_[jL][jL]!= 0){
	corrMatrix_[iL][jL] =matrix_[iL][jL]/sqrt(matrix_[iL][iL]*matrix_[jL][jL]); 
	p_corrMatrix->Fill(iL,jL,corrMatrix_[iL][jL]);
      }
    }
  }
}


bool PositionFit::fillMatrixFromFile(){

  std::ifstream fmatrix;
  std::ostringstream fmatrixname;
  fmatrixname << outFolder_ << "/errorMatrix.dat";
  fmatrix.open(fmatrixname.str());
  if (!fmatrix.is_open()){
    std::cout << " -- Cannot open outfile " << fmatrixname.str() << "! Refilling the matrix..." << std::endl;
    return false;
  }

  matrix_.ResizeTo(nLayers_,nLayers_);
  if (debug_) std::cout << " -- Error matrix: " << std::endl;
  while (!fmatrix.eof()){
    unsigned iL=nLayers_;
    unsigned jL=nLayers_;
    double m=0;
    fmatrix>>iL>>jL>>m;
    if (iL<nLayers_ && jL<nLayers_){
      if (debug_) std::cout << std::setprecision(15) << iL << " " << jL << " " << m << std::endl;
      matrix_[iL][jL] = m;
    }
  }
  
  fmatrix.close();

  return true;
}

bool PositionFit::performLeastSquareFit(TTree *aSimTree, 
					const unsigned nEvts){
  
  initialiseFitHistograms();

  //try reading matrix from file, if fail
  //return false and refill the matrix.
  if (!fillMatrixFromFile()) return false;
  
  //fill matrices
  fillCorrelationMatrix();

  //get back data for each event and perform chi2 fit:
  std::cout << " -- Performing chi2 fit for each event" << std::endl;

  unsigned nInvalidFits=0;
  //unsigned nFailedFits=0;
  unsigned nFailedFitsAfterCut=0;

  //open new file to save accurate positions
  std::ofstream fout;
  std::ostringstream foutname;
  foutname << outFolder_ << "/accuratePos.dat";
  fout.open(foutname.str());
  if (!fout.is_open()){
    std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }

  //std::vector<HGCSSRecoHit> * rechitvec = 0;
  //aRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  //std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;

  aSimTree->SetBranchAddress("HGCSSEvent",&event);
  aSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  aSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  aSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    aSimTree->GetEntry(ievt);

    //cut outliers
    unsigned fit = fitEvent(ievt,nInvalidFits,fout,true);
    if (fit==1){
      // std::cout << " -- Number of genparticles: " << (*genvec).size() << std::endl;
      // for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles 
      // 	std::cout << " --- particle " << iP << std::endl;
      // 	(*genvec)[iP].Print(std::cout);
      // }
      std::cout << " -- Event " << ievt << " skipped." << std::endl;
    }
    if (fit>1) {
      //std::cout << " ---- Fit failed ! Reprocess event " << ievt << " cutting outliers" << std::endl;
      //nFailedFits++;
      //if (fitEvent(ievt,nInvalidFits,fout,true)>1)
      std::cout << " -- Event " << ievt << " failed fit." << std::endl;
      nFailedFitsAfterCut++;
    }

  }//loop on entries
  
  fout.close();    


  std::cout << " -- Number of invalid fits: " << nInvalidFits << std::endl;
  //std::cout << " -- Number of fits failed: " << nFailedFits << std::endl;
  std::cout << " -- Number of fits failed after cutting outliers: " << nFailedFitsAfterCut << std::endl;
  return true;

}

unsigned PositionFit::fitEvent(const unsigned ievt, unsigned & nInvalidFits, std::ofstream & fout, const bool cutOutliers){

  std::vector<unsigned> layerId;
  std::vector<double> posx;
  std::vector<double> posy;
  std::vector<double> posz;
  std::vector<double> posxtruth;
  std::vector<double> posytruth;
  layerId.reserve(nLayers_);
  posx.reserve(nLayers_);
  posy.reserve(nLayers_);
  posz.reserve(nLayers_);
  posxtruth.reserve(nLayers_);
  posytruth.reserve(nLayers_);
  
  if (!getPositionFromFile(ievt,
			   layerId,posx,posy,posz,
			   posxtruth,posytruth,
			   cutOutliers)){
    nInvalidFits++;
    return 1;
  }

  const unsigned nL = layerId.size();
  
  
  //fill some control histograms
  //only once
  if (cutOutliers){
    for (unsigned iL(0); iL<nL;++iL){
      //std::cout << layerId[iL] << " " << posx[iL] << " " << posy[iL] << " " << posz[iL] << std::endl;
      p_recoXvsLayer->Fill(layerId[iL],posx[iL]);
      p_recoYvsLayer->Fill(layerId[iL],posy[iL]);
      p_recoZvsLayer->Fill(layerId[iL],posz[iL]);
      p_truthXvsLayer->Fill(layerId[iL],posxtruth[iL]);
      p_truthYvsLayer->Fill(layerId[iL],posytruth[iL]);
    }
  }

  p_nLayersFit->Fill(nL);

  //if less than 3 valid layers: no point doing a fit !!
  if (nL<3){
    nInvalidFits++;
    return 1;
  }
  
  //number of points: x and y per layer minus number of parameters: 2 for x + 2 for y.
  double ndf = 2*nL-4;
  
  //Get error matrix removing lines with zero hits
  TMatrixDSym e(nL);
  TVectorD u(nL),z(nL),x(nL),y(nL);
  
  for(unsigned i(0);i<nL;++i) {
    u(i)=1.0;
    z(i)=posz[i];
    //std::cout << "fit() z(" << i << ") = " << z(i) << std::endl;
    
    for(unsigned j(i);j<nL;++j) {
      e(i,j)=matrix_(layerId[i],layerId[j]);
      e(j,i)=matrix_(layerId[j],layerId[i]);
    }
  }
  
  e.Invert();
  
  //do fit for reco and truth
  double positionFF[2][2];
  double TanAngle[2][2];
 
  for (unsigned rt(0); rt<2;++rt){
    if (debug_) {
      std::cout << "... Processing ";
      if (rt==0) std::cout << " fit to reco position.";
      else std::cout << " fit to truth position.";
      std::cout << std::endl;
    }
    double chiSq(0.0);
    double position[2];
    
    TMatrixD fitMatrix(4,4);
    for (unsigned ii(0);ii<4;++ii){
      for (unsigned ij(0);ij<4;++ij){
	fitMatrix[ii][ij]=0;
      }
    }

    //resolve equation for x and y separately
    for(unsigned xy(0);xy<2;xy++) {//loop on x or y
      if (debug_) {
	std::cout << "... Processing ";
	if (xy==0) std::cout << " fit to x position.";
	else std::cout << " fit to y position.";
	std::cout << std::endl;
      }
      for(unsigned i(0);i<nL;i++) {
	x(i)= rt==0 ? ((xy==0) ? posx[i] : posy[i]) : ((xy==0) ? posxtruth[i] : posytruth[i]);
	//std::cout << "fit() x(" << i << ") = " << x(i) << std::endl;
      }
      
      TMatrixD w(2,2);
      TVectorD v(2),p(2);
      
      w(0,0)=u*(e*u);
      w(0,1)=u*(e*z);
      w(1,0)=z*(e*u);
      w(1,1)=z*(e*z);
      
      v(0)=u*(e*x);
      v(1)=z*(e*x);
      
      w.Invert();
      
      p=w*v;
      if (debug_) {
	std::cout << "fit() w(0,0) = " << w(0,0) << std::endl;
	std::cout << "fit() w(0,1) = " << w(0,1) << std::endl;
	std::cout << "fit() w(1,0) = " << w(1,0) << std::endl;
	std::cout << "fit() w(1,1) = " << w(1,1) << std::endl;	
	std::cout << "fit() p(0) = " << p(0) << std::endl;
	std::cout << "fit() p(1) = " << p(1) << std::endl;
      }
      
      position[xy] = p(0);
      positionFF[rt][xy] = p(0)+p(1)*posz[0];
      TanAngle[rt][xy] = p(1);

      //sanity check for nan values
      if (w(0,0)==w(0,0)) fitMatrix[2*xy][2*xy]=w(0,0);
      if (w(0,1)==w(0,1)) fitMatrix[2*xy][2*xy+1]=w(0,1);
      if (w(1,0)==w(1,0)) fitMatrix[2*xy+1][2*xy]=w(1,0);
      if (w(1,1)==w(1,1)) fitMatrix[2*xy+1][2*xy+1]=w(1,1);
      
      
      TVectorD dp(nL);
      for(unsigned i(0);i<nL;i++) {
	dp(i)=x(i)-p(0)-p(1)*z(i);
      }
      
      chiSq+=dp*(e*dp);
    }//loop on x or y
    
    //chi2 test
    if (chiSq/ndf>chi2ndfmax_) {
      std::cout << " ---- Fit failed for event " << ievt << std::endl;
      std::cout << "Chi2/ndf = " << chiSq << "/" << ndf << "=" << chiSq/ndf << std::endl;
      //std::cout << "fitw(0,0) = " << fitMatrix[0][0] << std::endl;
      //std::cout << "fitw(1,1) = " << fitMatrix[1][1] << std::endl;
      //std::cout << "fitw(2,2) = " << fitMatrix[2][2] << std::endl;
      //std::cout << "fitw(3,3) = " << fitMatrix[3][3] << std::endl;	
      std::cout << "ecal frontface position = " << positionFF[0][0] << " " << positionFF[0][1] << std::endl;
      std::cout << "ecal frontface tanAngle = " << TanAngle[0][0] << " " << TanAngle[0][1] << std::endl;
      return 2;
    }

    p_chi2[rt]->Fill(chiSq);
    p_chi2overNDF[rt]->Fill(chiSq/ndf);
    p_impactX[rt]->Fill(positionFF[rt][0]);
    p_tanAngleX[rt]->Fill(TanAngle[rt][0]);
    p_impactY[rt]->Fill(positionFF[rt][1]);
    p_tanAngleY[rt]->Fill(TanAngle[rt][1]);

    if (rt==0) {
      p_positionReso->Fill(sqrt(fitMatrix[0][0]));
      p_angularReso->Fill(sqrt(fitMatrix[1][1]));
      
      fout << ievt << " " 
	   << position[0] << " " 
	   << sqrt(fitMatrix[0][0]) << " " 
	   << TanAngle[0][0] << " " 
	   << sqrt(fitMatrix[1][1]) << " "
	   << position[1] << " " 
	   << sqrt(fitMatrix[2][2]) << " "
	   << TanAngle[0][1] << " "
	   << sqrt(fitMatrix[3][3])
	   << std::endl;
    
      for (unsigned iL(0); iL<nL;++iL){
	p_fitXvsLayer->Fill(layerId[iL],position[0]+TanAngle[0][0]*posz[iL]);
	p_fitYvsLayer->Fill(layerId[iL],position[1]+TanAngle[0][1]*posz[iL]);
      }
    }

  }//reco or truth

  p_impactX_residual->Fill(positionFF[0][0]-positionFF[1][0]);
  p_tanAngleX_residual->Fill(TanAngle[0][0]-TanAngle[1][0]);
  p_impactY_residual->Fill(positionFF[0][1]-positionFF[1][1]);
  p_tanAngleY_residual->Fill(TanAngle[0][1]-TanAngle[1][1]);
      

  return 0;
}

bool PositionFit::getPositionFromFile(const unsigned ievt,
				      std::vector<unsigned> & layerId,
				      std::vector<double> & posx,
				      std::vector<double> & posy,
				      std::vector<double> & posz,
				      std::vector<double> & posxtruth,
				      std::vector<double> & posytruth,
				      bool cutOutliers){


  std::ifstream fin;
  std::ostringstream finname;
  finname << outFolder_ << "/initialPos_evt" << ievt << ".dat";
  fin.open(finname.str());
  if (!fin.is_open()){
      std::cout << " Cannot open input file " << finname.str() << "!" << std::endl;
      return false;
  }
  
  while (!fin.eof()){
    unsigned l=nLayers_;
    double xr=0,yr=0,xt=0,yt=0;
    fin>>l>>xr>>yr>>xt>>yt;
    if (l<nLayers_){
      bool pass=fabs(xr-xt)<residualMax_ && fabs(yr-yt)<residualMax_;
      if (!cutOutliers || (cutOutliers && pass)){
	layerId.push_back(l);
	posx.push_back(xr);
	posy.push_back(yr);
	posz.push_back(avgZ_[l]);
	posxtruth.push_back(xt);
	posytruth.push_back(yt);	
      }
    }
  }
  
  fin.close();
  /*
  //@TODO to use something else than truth info :/
  if (cutOutliers){
    //
    for (unsigned i(0);i<layerId.size();++i){
      double xt=?;
      double yt=?;
      double xr=posx[i];
      double yr=posy[i];
      bool pass=fabs(xr-xt)<residualMax_ && fabs(yr-yt)<residualMax_;
      std::cout << i << " l=" << layerId[i] << " pass=" << pass << " xr=" << xr << " yr=" << yr << std::endl;
      if (!pass) {
	std::cout << " ---- erase layer " << *(layerId.begin()+i) << std::endl;
	layerId.erase(layerId.begin()+i);
	posx.erase(posx.begin()+i);
	posy.erase(posy.begin()+i);
	posxtruth.erase(posxtruth.begin()+i);
	posytruth.erase(posytruth.begin()+i);
	i--;
      }
    }
  }
  */

  return true;
}




