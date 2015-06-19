#include "HadEnergy.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"

HadEnergy::HadEnergy(HGCSSDetector & myDetector, TChain* lSimTree, TChain* lRecTree, TFile *outputFile, unsigned pNevts): myDetector_(myDetector),lSimTree_(lSimTree),lRecTree_(lRecTree), pNevts_(pNevts)
{ 
  nLayers_ = myDetector_.nLayers();
  nSections_ = myDetector_.nSections();

  FHtoEslope_ = 1.; 
  FHtoEoffset_ = 0.;
  BHtoEslope_ = 1.;
  BHtoEoffset_ = 0.;
  ECALslope_ = 1.;
  ECALoffset_ = 0.;

  FHtoBHslope_ = 1.;
  
  outputFile_ = outputFile;
  debug_ = false;

  bookHist(outputFile);
}


HadEnergy::~HadEnergy(){
}


void HadEnergy::bookHist(TFile *outputFile){

    outputFile_->cd();

    outtree_ = new TTree("Ereso","Tree to save energies");

    outtree_->Branch("wgtEtotal",&wgttotalE_);
    outtree_->Branch("correctedEtotal",&correctedtotalE_);
    outtree_->Branch("globalC",&Cglobal_);

    energy_.resize(nLayers_);
    std::ostringstream label;
    for (unsigned iL(0); iL<nLayers_;++iL){
	label.str("");
	label << "energy_" << iL;
	outtree_->Branch(label.str().c_str(),&energy_[iL]);
    }
 
     std::ostringstream lName;
     p_spectrum = new TH1F("p_spectrum","p_spectrum",1000,0,250);
     p_spectrumByLayer = new TH2F("p_spectrumByLayer","p_spectrumByLayer",nLayers_,0,nLayers_-1,100,0,5);
     p_spectrum_hightail = new TH1F("p_spectrum_hightail","p_spectrum_hightail",1000,0,250);
     p_spectrum_lowtail = new TH1F("p_spectrum_lowtail","p_spectrum_lowtail",1000,0,250);

}


bool HadEnergy::fillEnergies(){

  p_spectrum->Reset();
  std::vector<TH1F*> spectrumCollection;
  std::vector<double> Etotal;
  TH1F* p_Etotal = new TH1F("p_Etotal","p_Etotal",1000,0,100000);
  int nSec = nSections_;
  double EE(0), EFHCAL(0), EBHCAL(0);
  double EmipMeanFH(0);
  int    nhitsFH(0); 
  double recSum[nSec];

  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
     
  lSimTree_->SetBranchAddress("HGCSSEvent",&event);
  lSimTree_->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree_->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
     
  lRecTree_->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
   
  const unsigned nEvts = ((pNevts_ > lRecTree_->GetEntries() || pNevts_==0) ? static_cast<unsigned>(lRecTree_->GetEntries()) : pNevts_) ;
     
  std::cout << "- Processing = " << nEvts  << " events out of " << lRecTree_->GetEntries() << std::endl;
  
  spectrumCollection.reserve(nEvts);
  for (unsigned ievt(0); ievt<nEvts; ++ievt){// loop on entries
    if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
    lSimTree_->GetEntry(ievt);
    lRecTree_->GetEntry(ievt);
   
     if (debug_){
         std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
     }

    EmipMeanFH = 0;
    EE = EFHCAL = EBHCAL = 0;
    nhitsFH = 0;
    spectrumCollection[ievt] = new TH1F("","",1000,0,250);
    energy_.clear();
    wgttotalE_ = 0;
    correctedtotalE_.clear();
    for(unsigned iS(0); iS < nSec; iS++){
      recSum[iS] = 0;
    }

    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
      const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
      unsigned layer = lHit.layer();
      if (layer >= nLayers_) {
         std::cout << " WARNING! RecoHits with layer " << layer << " outside of detector's definition range ! Please fix the digitiser or the detector definition used here. Ignoring..." << std::endl;
       continue;
      }
     
      double absweight = (*ssvec)[layer].volX0trans()/(*ssvec)[0].volX0trans();
      unsigned sec =  myDetector_.getSection(layer); 

      double energy = lHit.energy();

      p_spectrumByLayer->Fill(layer,energy);

      energy_[layer] += energy*absweight; 
      recSum[sec] += energy*absweight;

      DetectorEnum type = myDetector_.detType(sec);
      if(type == DetectorEnum::FHCAL){
         spectrumCollection[ievt]->Fill(energy);
         EmipMeanFH += energy;
         nhitsFH += 1;
      }
    }//loop on hits

    for(unsigned iS(0); iS < nSec; iS++){
    //get total energy for each sub-detector
      DetectorEnum type = myDetector_.detType(iS); 
      if(type == DetectorEnum::FECAL || type == DetectorEnum::MECAL || type == DetectorEnum::BECAL)
         EE += (recSum[iS]-ECALoffset_)/ECALslope_ ;
      else if(type == DetectorEnum::FHCAL){
         EFHCAL += (recSum[iS]-FHtoEoffset_)/FHtoEslope_;
      }
      else if(type == DetectorEnum::BHCAL1 || type == DetectorEnum::BHCAL2)
         EBHCAL += (recSum[iS]-BHtoEoffset_)/BHtoEslope_; 
      else {
         std::cout << "the subdetector type is not defined" << std::endl;
       }
    }

    wgttotalE_ = EE + (EFHCAL-FHtoEoffset_)/FHtoEslope_ + ((EBHCAL-BHtoEoffset_)/BHtoEslope_)/FHtoBHslope_;
    Etotal.push_back(wgttotalE_);
    p_Etotal->Fill(wgttotalE_);
    if(nhitsFH!=0)EmipMeanFH = EmipMeanFH/nhitsFH; 

  //double GenE = genEn_*cosh(eta_);
  //if(EFHCAL < 0.9*GenE/1.25)p_spectrum_lowtail->Add(p_spectrum);
  //else if(EFHCAL > 1.1*GenE/1.25)p_spectrum_hightail->Add(p_spectrum);

    for(unsigned iLim(0); iLim < LimMIP_.size(); iLim++){
       double globalC = calcGlobalC(LimMIP_[iLim], EmipMeanFH, spectrumCollection[ievt]);
       Cglobal_.push_back(globalC); 
       correctedtotalE_.push_back(wgttotalE_*globalC);
    }

    outtree_->Fill();
  }
   
  double EMean = p_Etotal->GetMean();
  double ERMS = p_Etotal->GetRMS();
  for (unsigned ievt(0); ievt<nEvts; ++ievt){
    p_spectrum->Add(spectrumCollection[ievt]);
    if(Etotal[ievt] < EMean-ERMS)p_spectrum_lowtail->Add(spectrumCollection[ievt]);
    else if(Etotal[ievt] > EMean+ERMS)p_spectrum_hightail->Add(spectrumCollection[ievt]);
  }
 
  spectrumCollection.clear();
  outtree_->Write();
  p_spectrum->Write(); 
  p_spectrum_lowtail->Write();
  p_spectrum_hightail->Write();  
}

double HadEnergy::calcGlobalC(double LimMIP, double EmipMean, TH1F* spectrum){
  Int_t binLim = spectrum->GetXaxis()->FindBin(LimMIP);
  Int_t binAve = spectrum->GetXaxis()->FindBin(EmipMean);
  float countLim(0);
  for(unsigned iB(0); iB < binLim; iB++){
     countLim += spectrum->GetBinContent(iB);
  }
  float countAve(0);
  for(unsigned iB(0); iB < binAve; iB++){
     countAve += spectrum->GetBinContent(iB);
  }
  double globalC(0);
  if(countAve!=0)globalC = countLim/countAve;
 
  return globalC;
}

