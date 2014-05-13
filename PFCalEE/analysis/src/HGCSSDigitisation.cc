#include "HGCSSDigitisation.hh"
#include <cmath>
#include <sstream>
#include <iostream>

unsigned HGCSSDigitisation::nRandomPhotoElec(const double & aMipE){
  double mean = aMipE*npe_;
  int result = rndm_.Poisson(mean);
  if (result<0){
    std::cout << "WARNING!! HGCSSDigitisation::nRandomPhotoElec Poisson return negative number!! " << aMipE << " " << mean << " " << result << std::endl;

  }
  return static_cast<unsigned>(result);
}

double HGCSSDigitisation::nPixels(const double & aMipE){
  unsigned npe = nRandomPhotoElec(aMipE);
  double x = exp(-1.*npe/nTotal_);
  double npix = nTotal_*(1-x)/(1-crossTalk_*x);
  double result = positiveRandomGaus(npix);
  if (result<0) std::cout << "WARNING!! HGCSSDigitisation::nPixels negative result!! " << npe << " " << x << " " << npix << " " << result << std::endl;
  return result;
}

double HGCSSDigitisation::positiveRandomGaus(const double & mean){
  double result = rndm_.Gaus(mean,sigmaPix_);
  if (result<0) result = positiveRandomGaus(mean);
  return result;
}

double HGCSSDigitisation::digiE(const double & aMipE){
  double npix = nPixels(aMipE);
  double result = nTotal_/npe_*log((nTotal_-crossTalk_*npix)/(nTotal_-npix));
  if (result<0) std::cout << "WARNING!! HGCSSDigitisation::digiE negative result!! " << npix << " " << result << std::endl;
  return result;
}

void HGCSSDigitisation::addNoise(double & aDigiE, const Detector & adet){
  bool print = false;
  //if (aDigiE>0) print = true;
  if (print) std::cout << "HGCSSDigitisation::addNoise " << aDigiE << " ";
  double lNoise = rndm_.Gaus(0,noise_[adet]);
  aDigiE += lNoise;
  if (aDigiE<0) aDigiE = 0;
  if (print) std::cout << lNoise << " " << aDigiE << std::endl;
}

void HGCSSDigitisation::digitiseECAL(std::vector<TH2D *> & aHistVec){
  //To be ported from test/digitizer.cpp!
  for (unsigned iL(0); iL<aHistVec.size();++iL){
    for (int ix(1); ix<aHistVec[iL]->GetNbinsX()+1; ++ix){
      for (int iy(1); iy<aHistVec[iL]->GetNbinsY()+1; ++iy){
	double eTmp = aHistVec[iL]->GetBinContent(ix,iy);
	addNoise(eTmp,Detector::ECAL);
	double eDigi = adcConverter(eTmp,Detector::ECAL);
	aHistVec[iL]->SetBinContent(ix,iy,eDigi);
      }
    }
  }
}

double HGCSSDigitisation::adcConverter(const double & eMIP, const Detector & adet){
  return static_cast<unsigned>(eMIP*mipToADC_[adet])*1.0/mipToADC_[adet];
}

void HGCSSDigitisation::digitiseFHCAL(std::vector<TH2D *> & aHistVec){
  for (unsigned iL(0); iL<aHistVec.size();++iL){
    for (int ix(1); ix<aHistVec[iL]->GetNbinsX()+1; ++ix){
      for (int iy(1); iy<aHistVec[iL]->GetNbinsY()+1; ++iy){
	double eTmp = aHistVec[iL]->GetBinContent(ix,iy);
	double eDigi = 0;
	if (eTmp>0) {
	  eDigi = digiE(eTmp);
	  //std::cout << iL << " " << ix << " " << iy << " " << eTmp << " " << eDigi << std::endl;
	}
	addNoise(eDigi,Detector::FHCAL);
	aHistVec[iL]->SetBinContent(ix,iy,eDigi);
      }
    }
  }
}

void HGCSSDigitisation::digitiseBHCAL(std::vector<TH2D *> & aHistVec){
  //To be implemented!
  digitiseFHCAL(aHistVec);
}

void HGCSSDigitisation::Print(std::ostream & aOs) const{
  aOs << "====================================" << std::endl
      << " = Random seed: " << rndm_.GetSeed() << std::endl
      << " = Nphoto-electrons: " << npe_ << std::endl
      << " = cross-talk: " << crossTalk_ << std::endl
      << " = Npixels total: " << nTotal_ << std::endl
      << " = sigmaPixel: " << sigmaPix_ << std::endl
      << " = sigmaNoise: " << noise_.find(Detector::ECAL)->second << " " << noise_.find(Detector::FHCAL)->second << " " << noise_.find(Detector::BHCAL)->second << std::endl
      << " = MIPtoADC conversions: " << mipToADC_.find(Detector::ECAL)->second << " " << mipToADC_.find(Detector::FHCAL)->second << std::endl
      << "====================================" << std::endl;
};
