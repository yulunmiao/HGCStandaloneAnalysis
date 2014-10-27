#include <iomanip>

#include "PositionFit.hh"
#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSDetector.hh"


PositionFit::PositionFit(const unsigned nSR,
			 const double & residualMax, 
			 const unsigned nLayers, 
			 const unsigned nSiLayers,
			 unsigned debug){
  nSR_ = nSR;
  residualMax_ = residualMax;
  chi2ndfmax_ = 20;
  nLayers_ = nLayers;
  nSiLayers_ = nSiLayers;
  debug_ = debug;

  p_residuals_x = 0;
  p_residuals_y = 0;
  p_errorMatrix = 0;
  p_corrMatrix = 0;
  p_chi2[0] = 0;
  p_chi2[1] = 0;

  p_chi2overNDF[0] = 0;
  p_impactX[0] = 0;
  p_impactY[0] = 0;
  p_angleX[0] = 0;
  p_angleY[0] = 0;
  p_positionReso[0] = 0;
  p_angularReso[0] = 0;
  p_chi2overNDF[1] = 0;
  p_impactX[1] = 0;
  p_impactY[1] = 0;
  p_angleX[1] = 0;
  p_angleY[1] = 0;
  p_positionReso[1] = 0;
  p_angularReso[1] = 0;

}

void PositionFit::initialise(TFile *outputFile, 
			     std::string outFolder, 
			     const HGCSSGeometryConversion & geomConv, 
			     const HGCSSPUenergy & puDensity){
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
  //for xvsy plots
  double minX=-1700,maxX=1700;
  double minY=-1700,maxY=1700;
  //double minZ=3170,maxZ=3370;
  //double minX=-510,maxX=510;
  //double minY=-510,maxY=510;
  //double minZ=-1000,maxZ=1000;
  unsigned nX=(maxX-minX)/10,nY=(maxY-minY)/10;
  //unsigned nZ=maxZ-minZ;

  outputFile_->cd();
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

  p_etavsphi_max = new TH2F("p_etavsphi_max",";#phi_{max};#eta_{max};n_{events}",100,-3.1416,3.1416,100,1.4,3.6);

  p_residuals_x = new TH1F("p_residuals_x",";xreco-xtruth (mm)",1000,-50,50);
  p_residuals_y = new TH1F("p_residuals_y",";yreco-ytruth (mm)",1000,-50,50);
  p_residuals_x->StatOverflows();
  p_residuals_y->StatOverflows();

}

void PositionFit::initialiseFitHistograms(){

  outputFile_->cd();

  //check if already defined
  if (!p_chi2[0]){

    p_recoXvsLayer = new TH2F("p_recoXvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,200,-100,100);
    p_recoYvsLayer = new TH2F("p_recoYvsLayer",";layer;weighted y (mm);n_{events}",nLayers_,0,nLayers_,700,300,1000);
    p_recoZvsLayer = new TH2F("p_recoZvsLayer",";layer;avg z (mm);n_{events}",nLayers_,0,nLayers_,3000,3170,3470);
    p_truthXvsLayer = new TH2F("p_truthXvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,200,-100,100);
    p_truthYvsLayer = new TH2F("p_truthYvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,700,300,1000);
    p_fitXvsLayer = new TH2F("p_fitXvsLayer",";layer;fit x (mm);n_{events}",nLayers_,0,nLayers_,200,-100,100);
    p_fitYvsLayer = new TH2F("p_fitYvsLayer",";layer;fit y (mm);n_{events}",nLayers_,0,nLayers_,700,300,1000);
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
    p_impactY[0] = new TH1F("p_impactY",";y front face impact (mm);n_{events}",700,300,1000);
    p_impactY[1] = new TH1F("p_impactY_truth",";y front face impact (mm);n_{events}",700,300,1000);
    p_angleX[0] = new TH1F("p_angleX",";x direction angle (rad);n_{events}",500,-1,1);
    p_angleX[1] = new TH1F("p_angleX_truth",";x direction angle (rad);n_{events}",500,-1,1);
    p_angleY[0] = new TH1F("p_angleY",";y direction angle (rad);n_{events}",500,-1,1);
    p_angleY[1] = new TH1F("p_angleY_truth",";y direction angle (rad);n_{events}",500,-1,1);
    
    p_positionReso[0] = new TH1F("p_positionResoX",";#sigma_{x,y} (mm);n_{events}",500,0,20);
    p_positionReso[1] = new TH1F("p_positionResoY",";#sigma_{x,y} (mm);n_{events}",500,0,20);
    p_angularReso[0] = new TH1F("p_angularResoX",";#sigma_{#theta} (rad);n_{events}",500,0,1);
    p_angularReso[1] = new TH1F("p_angularResoY",";#sigma_{#theta} (rad);n_{events}",500,0,1);
  }
}

void PositionFit::getGlobalMaximum(std::vector<HGCSSRecoHit> *rechitvec,double & aPhimax,double & aEtamax){

  TH2F *etavsphi = new TH2F("etavsphi",";#phi;#eta;hits",150,-3.1416,3.1416,160,1.4,3.0);
  

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    
    double posx = lHit.get_x();
    double posy = lHit.get_y();
    double posz = lHit.get_z();
    //double radius = sqrt(posx*posx+posy*posy);
    double energy = lHit.energy();
    
    //if (energy>1) std::cout << "Hit " << layer << " " << posx << " " << posy << " " << posz << " " << energy << std::endl;
    
    if (debug_>1) {
      std::cout << " --  RecHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		<< " --  position x,y " << posx << "," << posy << std::endl;
      lHit.Print(std::cout);
    }
    
    ROOT::Math::XYZVector pos(posx,posy,posz);
    etavsphi->Fill(pos.phi(),pos.eta(),energy);
    
  }//loop on hits

  if (debug_)  std::cout << std::endl;

  
  //get position of maximum E tower
  int maxbin = etavsphi->GetMaximumBin();
  int binx,biny,binz;
  etavsphi->GetBinXYZ(maxbin,binx,biny,binz);
  aPhimax =etavsphi->GetXaxis()->GetBinCenter(binx); 
  aEtamax =etavsphi->GetYaxis()->GetBinCenter(biny); 

  if (debug_) std::cout << " MaxE cell eta,phi = " << aEtamax << " " << aPhimax << std::endl;

  etavsphi->Delete();
}

void PositionFit::getTruthPosition(std::vector<HGCSSGenParticle> *genvec,std::vector<ROOT::Math::XYPoint> & truthPos){

  for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
    //if ((*genvec).size()!= 1) (*genvec)[iP].Print(std::cout);
    if ((*genvec)[iP].trackID()==1){
      double x0 = (*genvec)[iP].x();
      double y0 = (*genvec)[iP].y();
      double z0 = (*genvec)[iP].z();
      double p = sqrt(pow((*genvec)[iP].px(),2)+pow((*genvec)[iP].py(),2)+pow((*genvec)[iP].pz(),2));
      //std::cout << "init : " << x0 << " " << y0 << " " << z0 << std::endl;
      //fill layers by propagating with momentum
      ROOT::Math::XYZVector unit((*genvec)[iP].px()/p,(*genvec)[iP].py()/p,(*genvec)[iP].pz()/p);
      
      //std::cout << " Gen particle eta,phi = " << unit.eta() << " " << unit.phi() << std::endl;
      
      for (unsigned iL(0); iL<nLayers_; ++iL){
	double xy = (avgZ_[iL]-z0)/sinh(unit.eta());
	double x = xy*cos(unit.phi())+x0;
	double y = xy*sin(unit.phi())+y0;
	
	//std::cout << "Lay " << iL << ": " << x << " " << y << " " << avgZ_[iL] << std::endl;
	p_genxy[iL]->Fill(x,y,1);
	truthPos[iL] = ROOT::Math::XYPoint(x,y);
      }
      
    }
  }//loop on gen particles
}

bool PositionFit::getZpositions(){
  std::ifstream fin;
  std::ostringstream finname;
  finname << outFolder_ << "_zPositions.dat";
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

  std::ofstream fout;
  std::ostringstream foutname;
  foutname << outFolder_ << "_zPositions.dat";
  fout.open(foutname.str());
  if (!fout.is_open()){
    std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
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
  aRecTree->SetBranchAddress("nPuVtx",&nPuVtx);

  std::cout << "- Processing = " << nEvts  << " events out of " << aSimTree->GetEntries() << std::endl;

  //bool firstEvent = true;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    aSimTree->GetEntry(ievt);
    aRecTree->GetEntry(ievt);

    if (debug_){
      std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
    }

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    //////// output files to save position for chi2 fit //////////
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    double phimax = 0;
    double etamax = 0;
    getGlobalMaximum(rechitvec,phimax,etamax);
    p_etavsphi_max->Fill(phimax,etamax);

    std::vector<ROOT::Math::XYPoint> truthPos;
    truthPos.resize(nLayers_,ROOT::Math::XYPoint());
    getTruthPosition(genvec,truthPos);

    std::vector<double> xmax;
    xmax.resize(nLayers_,0);
    std::vector<double> ymax;
    ymax.resize(nLayers_,0);
    getMaximumCell(rechitvec,phimax,etamax,xmax,ymax);

    std::vector<ROOT::Math::XYPoint> recoPos;
    recoPos.resize(nLayers_,ROOT::Math::XYPoint(0,0));
    std::vector<unsigned> nHits;
    nHits.resize(nLayers_,0);

    //get energy-weighted position around maximum
    getEnergyWeightedPosition(rechitvec,nPuVtx,xmax,ymax,recoPos,nHits);

    std::ofstream fout;
    std::ostringstream foutname;
    foutname << outFolder_ << "_initialPos_evt" << ievt << ".dat";
    fout.open(foutname.str());
    if (!fout.is_open()){
      std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
      exit(1);
    }

    if (debug_) std::cout << " Summary of reco and truth positions:" << std::endl;
    for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
      if (debug_) std::cout << iL << " nHits=" << nHits[iL] << " Max=(" << xmax[iL] << "," << ymax[iL] << ")\t Reco=(" << recoPos[iL].X() << "," << recoPos[iL].Y() << ")\t Truth=(" << truthPos[iL].X() << "," << truthPos[iL].Y() << ")" << std::endl;
      fout << iL << " " << recoPos[iL].X() << " " << recoPos[iL].Y() << " " << avgZ_[iL] << " " << truthPos[iL].X() << " " << truthPos[iL].Y() << std::endl;
    }

    fout.close();

    fillErrorMatrix(recoPos,truthPos,nHits);

    //firstEvent = false;
  }//loop on entries


}

void PositionFit::getMaximumCell(std::vector<HGCSSRecoHit> *rechitvec,const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax){
  
  std::vector<double> dRmin;
  dRmin.resize(nLayers_,10);
  
  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    if (debug_>1) {
      std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		<< " --  position x,y " << lHit.get_x() << "," << lHit.get_y() << std::endl;
      lHit.Print(std::cout);
    }
    
    unsigned layer = lHit.layer();
    if (layer >= nLayers_) {
      continue;
    }
    
    double posx = lHit.get_x();
    double posy = lHit.get_y();
    double posz = lHit.get_z();
    ROOT::Math::XYZVector pos(posx,posy,posz);
    double deta = fabs(pos.eta()-etamax);
    double dphi = fabs(pos.phi()-phimax);
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

void PositionFit::getEnergyWeightedPosition(std::vector<HGCSSRecoHit> *rechitvec,const unsigned nPU, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<ROOT::Math::XYPoint> & recoPos,std::vector<unsigned> & nHits,const bool puSubtracted){
  
  std::vector<double> eSum;
  eSum.resize(nLayers_,0);
  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    double energy = lHit.energy();//in MIP already...
    unsigned layer = lHit.layer();
    double posx = lHit.get_x();
    double posy = lHit.get_y();

    double step = geomConv_.cellSize()*nSR_/2.+0.1;//+0.1 to accomodate double precision
    if (fabs(posx-xmax[layer]) < step && 
	fabs(posy-ymax[layer]) < step){
      if (puSubtracted) {
	double leta = lHit.eta();
	energy = std::max(0.,energy - puDensity_.getDensity(leta,layer,geomConv_.cellSize(layer,leta),nPU));
      }
      recoPos[layer].SetX(recoPos[layer].X() + posx*energy);
      recoPos[layer].SetY(recoPos[layer].Y() + posy*energy);
      eSum[layer] += energy;
      if (energy>0) nHits[layer]++;
    }
    
  }//loop on rechits
  
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    if (nHits[iL]==0) continue;
    recoPos[iL].SetX(recoPos[iL].X()/eSum[iL]);
    recoPos[iL].SetY(recoPos[iL].Y()/eSum[iL]);
  }
  
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
  fmatrixname << outFolder_ << "_errorMatrix.dat";
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
  outputFile_->cd();

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
  fmatrixname << outFolder_ << "_errorMatrix.dat";
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

bool PositionFit::performLeastSquareFit(TTree *aRecTree, 
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
  unsigned nFailedFits=0;
  unsigned nFailedFitsAfterCut=0;

  //open new file to save accurate positions
  std::ofstream fout;
  std::ostringstream foutname;
  foutname << outFolder_ << "_accuratePos.dat";
  fout.open(foutname.str());
  if (!fout.is_open()){
    std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }

  //std::vector<HGCSSRecoHit> * rechitvec = 0;
  //aRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    //aRecTree->GetEntry(ievt);
    
    if (!fitEvent(ievt,nInvalidFits,fout)) {
      std::cout << " ---- Fit failed ! Reprocess event " << ievt << " cutting outliers" << std::endl;
      nFailedFits++;
      if (!fitEvent(ievt,nInvalidFits,fout,true)) nFailedFitsAfterCut++;
    }

  }//loop on entries
  
  fout.close();    


  std::cout << " -- Number of invalid fits: " << nInvalidFits << std::endl;
  std::cout << " -- Number of fits failed: " << nFailedFits << std::endl;
  std::cout << " -- Number of fits failed after cutting outliers: " << nFailedFitsAfterCut << std::endl;
  return true;

}

bool PositionFit::fitEvent(const unsigned ievt, unsigned & nInvalidFits, std::ofstream & fout, const bool cutOutliers){

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
			   cutOutliers)) return false;
  
  const unsigned nL = layerId.size();
  
  
  //fill some control histograms
  //only once
  if (!cutOutliers){
    for (unsigned iL(0); iL<nL;++iL){
      //std::cout << layerId[iL] << " " << posx[iL] << " " << posy[iL] << " " << posz[iL] << std::endl;
      p_recoXvsLayer->Fill(layerId[iL],posx[iL]);
      p_recoYvsLayer->Fill(layerId[iL],posy[iL]);
      p_recoZvsLayer->Fill(layerId[iL],posz[iL]);
      p_truthXvsLayer->Fill(layerId[iL],posxtruth[iL]);
      p_truthYvsLayer->Fill(layerId[iL],posytruth[iL]);
    }
  }

  //if less than 3 valid layers: no point doing a fit !!
  if (nL<3){
    nInvalidFits++;
    return false;
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
  
  for (unsigned rt(0); rt<2;++rt){
    if (debug_) {
      std::cout << "... Processing ";
      if (rt==0) std::cout << " fit to reco position.";
      else std::cout << " fit to truth position.";
      std::cout << std::endl;
    }
    double chiSq(0.0);
    double position[2];
    double positionFF[2];
    double TanAngle[2];
    
    TMatrixD fitMatrix(4,4);
    
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
      positionFF[xy] = p(0)+p(1)*posz[0];
      TanAngle[xy] = p(1);
      
      fitMatrix[2*xy][2*xy]=w(0,0);
      fitMatrix[2*xy][2*xy+1]=w(0,1);
      fitMatrix[2*xy+1][2*xy]=w(1,0);
      fitMatrix[2*xy+1][2*xy+1]=w(1,1);
      
      
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
      std::cout << "fitw(0,0) = " << fitMatrix[0][0] << std::endl;
      std::cout << "fitw(1,1) = " << fitMatrix[1][1] << std::endl;
      std::cout << "fitw(2,2) = " << fitMatrix[2][2] << std::endl;
      std::cout << "fitw(3,3) = " << fitMatrix[3][3] << std::endl;	
      std::cout << "ecal frontface position = " << positionFF[0] << " " << positionFF[1] << std::endl;
      std::cout << "ecal frontface angle = " << atan(TanAngle[0]) << " " << atan(TanAngle[1]) << std::endl;
      return false;
    }

    p_chi2[rt]->Fill(chiSq);
    p_chi2overNDF[rt]->Fill(chiSq/ndf);
    p_impactX[rt]->Fill(positionFF[0]);
    p_angleX[rt]->Fill(atan(TanAngle[0]));
    p_impactY[rt]->Fill(positionFF[1]);
    p_angleY[rt]->Fill(atan(TanAngle[1]));
    
    if (rt==0) {
      p_positionReso[0]->Fill(sqrt(fitMatrix[0][0]));
      p_positionReso[1]->Fill(sqrt(fitMatrix[2][2]));
      p_angularReso[0]->Fill(sqrt(fitMatrix[1][1]));
      p_angularReso[1]->Fill(sqrt(fitMatrix[3][3]));
      
      
      fout << ievt << " " 
	   << position[0] << " " 
	   << sqrt(fitMatrix[0][0]) << " " 
	   << TanAngle[0] << " " 
	   << sqrt(fitMatrix[1][1]) << " "
	   << position[1] << " " 
	   << sqrt(fitMatrix[2][2]) << " "
	   << TanAngle[1] << " "
	   << sqrt(fitMatrix[3][3])
	   << std::endl;
    }
    
    p_nLayersFit->Fill(nL);
    for (unsigned iL(0); iL<nL;++iL){
      p_fitXvsLayer->Fill(layerId[iL],position[0]+TanAngle[0]*posz[iL]);
      p_fitYvsLayer->Fill(layerId[iL],position[1]+TanAngle[1]*posz[iL]);
    }

  }//reco or truth

  return true;
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
  finname << outFolder_ << "_initialPos_evt" << ievt << ".dat";
  fin.open(finname.str());
  if (!fin.is_open()){
      std::cout << " Cannot open input file " << finname.str() << "! Refilling..." << std::endl;
      return false;
  }
  
  while (!fin.eof()){
    unsigned l=nLayers_;
    double xr=0,yr=0,z=0,xt=0,yt=0;
    fin>>l>>xr>>yr>>z>>xt>>yt;
    if (l<nLayers_){
      bool pass = fabs(xr-xt)<residualMax_ && fabs(yr-yt)<residualMax_;
      if (!cutOutliers || (cutOutliers && pass)){
	layerId.push_back(l);
	posx.push_back(xr);
	posy.push_back(yr);
	posz.push_back(z);
	posxtruth.push_back(xt);
	posytruth.push_back(yt);	
      }
    }
  }
  
  fin.close();
  return true;
}




