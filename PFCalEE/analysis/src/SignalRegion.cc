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

    nevt_ = nevt;
    nLayers_ = nLayers;
    geomConv_ = geomConv;
    puDensity_ = puDensity;
    fixForPuMixBug_ = applyPuMixFix;

    double xpos(0),ypos(0),zpos(0),xangle(0),yangle(0);    
    int layerIndex(0),eventIndex(0);

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
    for(unsigned ievt(0);ievt<nevt_;ievt++){ 
        std::string fitMatrix[4];
        fxypos >> eventIndex >> xpos >> fitMatrix[0] >> xangle >> fitMatrix[1] >> ypos >> fitMatrix[2] >> yangle >> fitMatrix[3];
        std::vector<ROOT::Math::XYZVector> tmpXYZ;
        for(unsigned iL(0);iL<nLayers_;iL++){
           ROOT::Math::XYZVector position( xpos + xangle*zPos_[iL], ypos + yangle*zPos_[iL], zPos_[iL]);
           tmpXYZ.push_back(position);
        }
        accuratePos_.push_back(tmpXYZ);
    }

}


SignalRegion::~SignalRegion(){
};



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

    for (unsigned ievt(0); ievt< nevt_; ++ievt){//loop on entries
        if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
        aSimTree->GetEntry(ievt);
        aRecTree->GetEntry(ievt);
 
	//fill weights for first event only: same in all events
        if (ievt==0){
	  absweight_.clear();
	  absweight_.reserve(nLayers_);
	  for(unsigned iL(0); iL<nLayers_; iL++){
            absweight_.push_back((*ssvec)[iL].volX0trans()/(*ssvec)[1].volX0trans() );
	  }
	}
        std::vector<ROOT::Math::XYZVector> eventPos = accuratePos_[ievt];
        if(eventPos.size()!=nLayers_) continue; 
	std::vector<double> accurateX;
	std::vector<double> accurateY;

        double totalE(0), wgttotalE(0);

	std::vector<double> signalSR0, signalSR1, signalSR2, signalSR3, signalSR4;
	std::vector<double> subtractSR0, subtractSR1, subtractSR2, subtractSR3, subtractSR4;
	signalSR0.resize(eventPos.size(),0);
	signalSR1.resize(eventPos.size(),0);
	signalSR2.resize(eventPos.size(),0);
	signalSR3.resize(eventPos.size(),0);
	signalSR4.resize(eventPos.size(),0);
	subtractSR0.resize(eventPos.size(),0);
	subtractSR1.resize(eventPos.size(),0);
	subtractSR2.resize(eventPos.size(),0);
	subtractSR3.resize(eventPos.size(),0);
	subtractSR4.resize(eventPos.size(),0);
	// set the accurate coordinate for each layer
	for(unsigned iL(0); iL< eventPos.size(); iL++){
	    accurateX.push_back(eventPos[iL].x());
	    accurateY.push_back(eventPos[iL].y());
	}
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

            totalE += energy;
            wgttotalE += energy*absweight_[layer];    

            double subtractedenergy = std::max(0.,energy - puDensity_.getDensity(leta,layer,geomConv_.cellSizeInCm(layer,leta),nPuVtx));
            double halfCell = 0.5*geomConv_.cellSize(layer,leta);

	    //SR0
	    if(fabs(posx + halfCell-accurateX[layer])< halfCell && fabs(posy+ halfCell-accurateY[layer])< halfCell){
                signalSR0[layer] += energy;
	        subtractSR0[layer] += subtractedenergy;
            }
            //SR1
	    if(fabs(posx + halfCell-accurateX[layer]) < 2*halfCell && fabs(posy+ halfCell-accurateY[layer]) < 2*halfCell){
	        signalSR1[layer] += energy; 
	        subtractSR1[layer] += subtractedenergy;
            }
	    //SR2
	    if(fabs(posx + halfCell-accurateX[layer]) < 3*halfCell && fabs(posy+ halfCell-accurateY[layer]) < 3*halfCell){
		signalSR2[layer] += energy;
	        subtractSR2[layer] += subtractedenergy;
            }
	    //SR3
	    if(fabs(posx + halfCell-accurateX[layer]) < 4*halfCell && fabs(posy+ halfCell-accurateY[layer]) < 4*halfCell){
		signalSR3[layer] += energy;
	        subtractSR3[layer] += subtractedenergy;
            }
	    //SR4
	    if(fabs(posx + halfCell-accurateX[layer]) < 5*halfCell && fabs(posy+ halfCell-accurateY[layer]) < 5*halfCell){
		signalSR4[layer] += energy;
	        subtractSR4[layer] += subtractedenergy;
            }

	 }

         totalE_.push_back(totalE);
         wgttotalE_.push_back(wgttotalE); 
	 energySR0_.push_back(signalSR0); 
	 energySR1_.push_back(signalSR1); 
	 energySR2_.push_back(signalSR2); 
	 energySR3_.push_back(signalSR3); 
	 energySR4_.push_back(signalSR4); 
	 subtractedenergySR0_.push_back(subtractSR0); 
	 subtractedenergySR1_.push_back(subtractSR1); 
	 subtractedenergySR2_.push_back(subtractSR2); 
	 subtractedenergySR3_.push_back(subtractSR3); 
	 subtractedenergySR4_.push_back(subtractSR4); 

     }
}


void SignalRegion::initialiseHistograms(){

    outputFile_->cd();

  //check if already defined
    p_rawEtotal = new TH1F("p_rawEtotal", "Total E (MIP)", 5000,0,200000);
    p_wgtEtotal = new TH1F("p_wgtEtotal", "Total weighted E (MIP)",5000, 0, 200000);

    p_rawESR0 = new TH1F("p_rawESR0"," E, SR0", 5000,0,200000);
    p_rawESR1 = new TH1F("p_rawESR1"," E, SR1", 5000,0,200000);
    p_rawESR2 = new TH1F("p_rawESR2"," E, SR2", 5000,0,200000);
    p_rawESR3 = new TH1F("p_rawESR3"," E, SR3", 5000,0,200000);
    p_rawESR4 = new TH1F("p_rawESR4"," E, SR4", 5000,0,200000);
    p_rawESR0->StatOverflows();
    p_rawESR1->StatOverflows();
    p_rawESR2->StatOverflows();
    p_rawESR3->StatOverflows();
    p_rawESR4->StatOverflows();
    
    p_wgtESR0 = new TH1F("p_wgtESR0"," E, SR0", 5000,0,200000);
    p_wgtESR1 = new TH1F("p_wgtESR1"," E, SR1", 5000,0,200000);
    p_wgtESR2 = new TH1F("p_wgtESR2"," E, SR2", 5000,0,200000);
    p_wgtESR3 = new TH1F("p_wgtESR3"," E, SR3", 5000,0,200000);
    p_wgtESR4 = new TH1F("p_wgtESR4"," E, SR4", 5000,0,200000);
    p_wgtESR0->StatOverflows();
    p_wgtESR1->StatOverflows();
    p_wgtESR2->StatOverflows();
    p_wgtESR3->StatOverflows();
    p_wgtESR4->StatOverflows();
   
    p_rawSubtractESR0 = new TH1F("p_rawSubtractESR0"," E (with PU subtraction), SR0", 5000,0,200000);
    p_rawSubtractESR1 = new TH1F("p_rawSubtractESR1"," E (with PU subtraction), SR1", 5000,0,200000);
    p_rawSubtractESR2 = new TH1F("p_rawSubtractESR2"," E (with PU subtraction), SR2", 5000,0,200000);
    p_rawSubtractESR3 = new TH1F("p_rawSubtractESR3"," E (with PU subtraction), SR3", 5000,0,200000);
    p_rawSubtractESR4 = new TH1F("p_rawSubtractESR4"," E (with PU subtraction), SR4", 5000,0,200000);
    p_rawSubtractESR0->StatOverflows();
    p_rawSubtractESR1->StatOverflows();
    p_rawSubtractESR2->StatOverflows();
    p_rawSubtractESR3->StatOverflows();
    p_rawSubtractESR4->StatOverflows();

    p_wgtSubtractESR0 = new TH1F("p_wgtSubtractESR0"," E (with PU subtraction), SR0", 5000,0,200000);
    p_wgtSubtractESR1 = new TH1F("p_wgtSubtractESR1"," E (with PU subtraction), SR1", 5000,0,200000);
    p_wgtSubtractESR2 = new TH1F("p_wgtSubtractESR2"," E (with PU subtraction), SR2", 5000,0,200000);
    p_wgtSubtractESR3 = new TH1F("p_wgtSubtractESR3"," E (with PU subtraction), SR3", 5000,0,200000);
    p_wgtSubtractESR4 = new TH1F("p_wgtSubtractESR4"," E (with PU subtraction), SR4", 5000,0,200000);
    p_wgtSubtractESR0->StatOverflows();
    p_wgtSubtractESR1->StatOverflows();
    p_wgtSubtractESR2->StatOverflows();
    p_wgtSubtractESR3->StatOverflows();
    p_wgtSubtractESR4->StatOverflows();
}

void SignalRegion::fillHistograms(){
        
    //Fill energy without PU subtraction
    bool subtractPU = false;
    for(unsigned ievt(0); ievt < nevt_; ievt++){

        p_rawEtotal->Fill(totalE_[ievt]);
        p_wgtEtotal->Fill(wgttotalE_[ievt]);

        p_rawESR0->Fill( getEtotalSR0(ievt, subtractPU));
        p_rawESR1->Fill( getEtotalSR1(ievt, subtractPU));
        p_rawESR2->Fill( getEtotalSR2(ievt, subtractPU));
        p_rawESR3->Fill( getEtotalSR3(ievt, subtractPU));
        p_rawESR4->Fill( getEtotalSR4(ievt, subtractPU));

        double wgtESR0(0), wgtESR1(0), wgtESR2(0), wgtESR3(0), wgtESR4(0);
        for(unsigned iL(0); iL < nLayers_;iL++){
           wgtESR0 += getSR0(ievt, iL, subtractPU)*absweight(iL);
           wgtESR1 += getSR1(ievt, iL, subtractPU)*absweight(iL);
           wgtESR2 += getSR2(ievt, iL, subtractPU)*absweight(iL);
           wgtESR3 += getSR3(ievt, iL, subtractPU)*absweight(iL);
           wgtESR4 += getSR4(ievt, iL, subtractPU)*absweight(iL);
        }
        p_wgtESR0->Fill(wgtESR0);
        p_wgtESR1->Fill(wgtESR1);
        p_wgtESR2->Fill(wgtESR2);
        p_wgtESR3->Fill(wgtESR3);
        p_wgtESR4->Fill(wgtESR4);
    }
 

   //Fill energy after PU subtraction
    subtractPU = true;
    for(unsigned ievt(0); ievt < nevt_; ievt++){
        p_rawSubtractESR0->Fill( getEtotalSR0(ievt, subtractPU));
        p_rawSubtractESR1->Fill( getEtotalSR1(ievt, subtractPU));
        p_rawSubtractESR2->Fill( getEtotalSR2(ievt, subtractPU));
        p_rawSubtractESR3->Fill( getEtotalSR3(ievt, subtractPU));
        p_rawSubtractESR4->Fill( getEtotalSR4(ievt, subtractPU));

        double wgtSubtractESR0(0), wgtSubtractESR1(0), wgtSubtractESR2(0), wgtSubtractESR3(0), wgtSubtractESR4(0);
        for(unsigned iL(0); iL < nLayers_;iL++){
           wgtSubtractESR0 += getSR0(ievt, iL, subtractPU)*absweight(iL);
           wgtSubtractESR1 += getSR1(ievt, iL, subtractPU)*absweight(iL);
           wgtSubtractESR2 += getSR2(ievt, iL, subtractPU)*absweight(iL);
           wgtSubtractESR3 += getSR3(ievt, iL, subtractPU)*absweight(iL);
           wgtSubtractESR4 += getSR4(ievt, iL, subtractPU)*absweight(iL);
        }
        p_wgtSubtractESR0->Fill(wgtSubtractESR0);
        p_wgtSubtractESR1->Fill(wgtSubtractESR1);
        p_wgtSubtractESR2->Fill(wgtSubtractESR2);
        p_wgtSubtractESR3->Fill(wgtSubtractESR3);
        p_wgtSubtractESR4->Fill(wgtSubtractESR4);
    }

    std::cout << " -- Histograms for signal regions have been filled !" << std::endl;

}




