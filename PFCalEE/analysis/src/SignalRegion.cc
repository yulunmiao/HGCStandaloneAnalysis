#include "SignalRegion.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"

SignalRegion:: SignalRegion(const std::string inputFolder,
                            const unsigned nLayers,
                            const unsigned nevt,
                            const HGCSSGeometryConversion & geomConv,
                            const HGCSSPUenergy & puDensity,
			    const bool applyPuMixFix){

  nSR_ = 5;
  nevt_ = nevt;
  nLayers_ = nLayers;
  geomConv_ = geomConv;
  puDensity_ = puDensity;
  fixForPuMixBug_ = applyPuMixFix;
  
  double zpos(0);    
  int layerIndex(0);
  
  std::ifstream fzpos;
  std::ostringstream finname;
  finname << inputFolder << "/zPositions.dat";
  fzpos.open(finname.str());
  if (!fzpos.is_open()){
    std::cout << " Cannot open input file " << finname.str() << "! Exiting..." << std::endl;
    exit(1);
  }
  for(unsigned iL(0);iL<nLayers_;iL++){
    fzpos >> layerIndex >> zpos;
    zPos_.push_back(zpos);
  }
  
  std::ifstream fxypos;
  finname.str("");
  finname << inputFolder << "/accuratePos.dat";
  fxypos.open(finname.str());
  if (!fxypos.is_open()){
    std::cout << " Cannot open input file " << finname.str() << "! Exiting..." << std::endl;
    exit(1);
  }

  //all events did not pass the chi2 fit: fill only those found.
  //keep failed ones to emptyvec so they are ignored afterwards
  accuratePos_.clear();
  std::vector<ROOT::Math::XYZVector> emptyvec;
  accuratePos_.resize(nevt_,emptyvec);

  while (!fxypos.eof()){
    unsigned eventIndex = nevt_;
    double xpos(0),ypos(0),xangle(0),yangle(0);
    double fitMatrix[4] = {0,0,0,0};
    fxypos >> eventIndex >> xpos >> fitMatrix[0] >> xangle >> fitMatrix[1] >> ypos >> fitMatrix[2] >> yangle >> fitMatrix[3];
    //testing for nan
    if ( eventIndex != eventIndex || xpos != xpos || fitMatrix[0]!=fitMatrix[0] || xangle!=xangle || fitMatrix[1]!=fitMatrix[1] || ypos!=ypos || fitMatrix[2]!=fitMatrix[2] || yangle!=yangle || fitMatrix[3]!=fitMatrix[3]){
      std::cout << " Found nan ! Fix code !" << std::endl;
      std::cout << eventIndex << " " << xpos << " " << fitMatrix[0] << " " << xangle << " " << fitMatrix[1] << " " << ypos << " " << fitMatrix[2] << " " << yangle << " " << fitMatrix[3]<< std::endl;
      exit(1);
    }
    if (eventIndex<nevt_) {
      std::vector<ROOT::Math::XYZVector> tmpXYZ;
      tmpXYZ.reserve(nLayers_);
      for(unsigned iL(0);iL<nLayers_;iL++){
	ROOT::Math::XYZVector position( xpos + xangle*zPos_[iL], ypos + yangle*zPos_[iL], zPos_[iL]);
	tmpXYZ.push_back(position);
      }
      accuratePos_[eventIndex] = tmpXYZ;
    }
    else break;
  }
  
}


SignalRegion::~SignalRegion(){
}

void SignalRegion::initialise(TTree *aSimTree, TTree *aRecTree, 
			      TFile *outputFile){

  //mycalib_ = mycalib;
  setOutputFile(outputFile);

  initialiseHistograms();
  
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

  std::cout << " -- Now filling signal region histograms..." << std::endl;
  unsigned nSkipped = 0;
    for (unsigned ievt(0); ievt< nevt_; ++ievt){//loop on entries
        if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
        aSimTree->GetEntry(ievt);
        aRecTree->GetEntry(ievt);
 
	//fill weights for first event only: same in all events
        if (ievt==0){
	  absweight_.clear();
	  absweight_.reserve(nLayers_);
	  std::cout << " -- Absorber weights used for total energy:" << std::endl;
	  for(unsigned iL(0); iL<nLayers_; iL++){
	    double w = (*ssvec)[iL].volX0trans()/(*ssvec)[1].volX0trans();
	    std::cout << " - Layer " << iL << " w=" << w << std::endl;
            absweight_.push_back(w);
	  }
	}
	if (absweight_.size()!=nLayers_) {
	  std::cout << " -- Error! Not all layers found! Fix code." << std::endl;
	  exit(1);
	}
        std::vector<ROOT::Math::XYZVector> eventPos = accuratePos_[ievt];

        if(eventPos.size()!=nLayers_) {
	  std::cout << " -- Event " << ievt << " skipped, accurate position size = " << eventPos.size() << std::endl;
	  nSkipped++;
	  continue;
	}

	//initialise values for current event
        totalE_ = 0;
	wgttotalE_ = 0;

	for (unsigned iL(0); iL<nLayers_;++iL){
	  for (unsigned iSR(0);iSR<nSR_;++iSR){
	    energySR_[iL][iSR] = 0;
	    subtractedenergySR_[iL][iSR] = 0;
	  }
	}


	//get event-by-event PU
	//get PU contrib from elsewhere in the event
	//loop over phi with same etamax
	//take average per layer: not all 9 cells of 3*3 area have hits...
	/*	std::vector<double> puE;
	puE.resize(nLayers_,0);
	if (nPuVtx>0){
	  unsigned nRandomCones = 50;
	  double phistep = TMath::Pi()/nRandomCones;
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
	    getPuContribution(rechitvec,xmaxrc,ymaxrc,puE);
	  }
      
	  //normalise to one cell: must count cells with 0 hit !
	  //use cell size at etamax...
	  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
	    unsigned nCells = nRandomCones*nSR_*geomConv_.cellSize()/geomConv_.cellSize(iL,etamax);
	    puE[iL] = puE[iL]/nCells;
	  }

	}//if PU
	*/


	// Define different signal region and sum over energy
	for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	    const HGCSSRecoHit & lHit = (*rechitvec)[iH];

	    unsigned layer = lHit.layer();
	    if (layer >= nLayers_) {
		continue;
	    }
	    double posx = lHit.get_x();
	    if (fixForPuMixBug_) posx-=1.25;
	    double posy = lHit.get_y();
	    if (fixForPuMixBug_) posy-=1.25;
	    double energy = lHit.energy();
            double leta = lHit.eta();

            totalE_ += energy;
            wgttotalE_ += energy*absweight_[layer];    

	    double puE = puDensity_.getDensity(leta,layer,geomConv_.cellSizeInCm(layer,leta),nPuVtx);
            double subtractedenergy = std::max(0.,energy - puE);
            double halfCell = 0.5*geomConv_.cellSize(layer,leta);

	    double dx = eventPos[layer].x()-posx;
	    double dy = eventPos[layer].y()-posy;
	    
	    //SR0-4
	    for (unsigned isr(0); isr<nSR_;++isr){
	      if ( (fabs(dx) <= ((isr+1)*halfCell)) && (fabs(dy) <= ((isr+1)*halfCell))){
		energySR_[layer][isr] += energy;
		subtractedenergySR_[layer][isr] += subtractedenergy;
	      }
	    }
	}//loop on hits
	
	fillHistograms();
	outtree_->Fill();

    }//loop on events

    std::cout << " -- Histograms for signal regions have been filled !" << std::endl;
    std::cout << " -- Number of skipped events: " << nSkipped << std::endl;
}


void SignalRegion::initialiseHistograms(){

    outputFile_->cd();

    outtree_ = new TTree("Ereso","Tree to save energies in signal regions");

    outtree_->Branch("rawEtotal",&totalE_);
    outtree_->Branch("wgtEtotal",&wgttotalE_);

    std::vector<double> emptyvec;
    emptyvec.resize(nSR_,0);
    energySR_.resize(nLayers_,emptyvec);
    subtractedenergySR_.resize(nLayers_,emptyvec);

    std::ostringstream label;
    for (unsigned iL(0); iL<nLayers_;++iL){
      for (unsigned iSR(0);iSR<nSR_;++iSR){
	label.str("");
	label << "energy_" << iL << "_SR" << iSR;
	outtree_->Branch(label.str().c_str(),&energySR_[iL][iSR]);
	label.str("");
	label << "subtractedenergy_" << iL << "_SR" << iSR;
	outtree_->Branch(label.str().c_str(),&subtractedenergySR_[iL][iSR]);
      }
    }

    p_rawEtotal = new TH1F("p_rawEtotal", "Total E (MIP)", 5000,0,200000);
    p_wgtEtotal = new TH1F("p_wgtEtotal", "Total weighted E (MIP)",5000, 0, 200000);

    p_rawESR.resize(nSR_,0);
    p_wgtESR.resize(nSR_,0);
    p_rawSubtractESR.resize(nSR_,0);
    p_wgtSubtractESR.resize(nSR_,0);

    for (unsigned iSR(0);iSR<nSR_;++iSR){
      label.str("");
      label << "rawESR" << iSR;
      p_rawESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR} (MIPs);events", 5000,0,200000);
      p_rawESR[iSR]->StatOverflows();

      label.str("");
      label << "wgtESR" << iSR;
      p_wgtESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR} (MIPs);events", 5000,0,200000);
      p_wgtESR[iSR]->StatOverflows();
      
      label.str("");
      label << "rawSubtractESR" << iSR;
      p_rawSubtractESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR}^{PUsubtr} (MIPs);events", 5000,0,200000);
      p_rawSubtractESR[iSR]->StatOverflows();
      label.str("");
      label << "wgtSubtractESR" << iSR;
      p_wgtSubtractESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR}^{PUsubtr} (MIPs);events", 5000,0,200000);
      p_wgtSubtractESR[iSR]->StatOverflows();
    }//loop on sr

}

void SignalRegion::fillHistograms(){
        
  p_rawEtotal->Fill(totalE_);
  p_wgtEtotal->Fill(wgttotalE_);

  for (unsigned iSR(0);iSR<nSR_;++iSR){
    //Fill energy without PU subtraction
    bool subtractPU = false;
    p_rawESR[iSR]->Fill( getEtotalSR(iSR, subtractPU));

    double wgtESR = 0;
    for(unsigned iL(0); iL < nLayers_;iL++){
      wgtESR += getSR(iSR, iL, subtractPU)*absweight(iL);

    }
    p_wgtESR[iSR]->Fill(wgtESR);

    //Fill energy after PU subtraction
    subtractPU = true;
    p_rawSubtractESR[iSR]->Fill( getEtotalSR(iSR, subtractPU));
    double wgtSubtractESR = 0;
    for(unsigned iL(0); iL < nLayers_;iL++){
      wgtSubtractESR += getSR(iSR, iL, subtractPU)*absweight(iL);
    }
    p_wgtSubtractESR[iSR]->Fill(wgtSubtractESR);

  }//loop on SR

}




