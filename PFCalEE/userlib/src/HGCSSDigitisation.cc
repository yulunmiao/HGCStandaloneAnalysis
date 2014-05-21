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

double HGCSSDigitisation::mipCor(const double & aMipE,
				 const double & posx, 
				 const double & posy,
				 const double & posz){
  double costheta = fabs(posz)/sqrt(posz*posz+posx*posx+posy*posy);
  if (costheta>0) return aMipE*costheta;
  return aMipE;
}

double HGCSSDigitisation::digiE(const double & aMipE){
  double npix = nPixels(aMipE);
  double result = nTotal_/npe_*log((nTotal_-crossTalk_*npix)/(nTotal_-npix));
  if (result<0) std::cout << "WARNING!! HGCSSDigitisation::digiE negative result!! " << npix << " " << result << std::endl;
  return result;
}

void HGCSSDigitisation::addNoise(double & aDigiE, const unsigned & alay ,
				 TH1F * & hist){
  bool print = false;
  //if (aDigiE>0) print = true;
  if (print) std::cout << "HGCSSDigitisation::addNoise " << aDigiE << " ";
  double lNoise = rndm_.Gaus(0,noise_[alay]);
  if (hist) hist->Fill(lNoise);
  aDigiE += lNoise;
  if (aDigiE<0) aDigiE = 0;
  if (print) std::cout << lNoise << " " << aDigiE << std::endl;
}

unsigned HGCSSDigitisation::adcConverter(double eMIP, DetectorEnum adet){
  if (eMIP<0) eMIP=0;
  return static_cast<unsigned>(eMIP*mipToADC_[adet]);
}

double HGCSSDigitisation::adcToMIP(const unsigned adcCounts, DetectorEnum adet){
  return adcCounts*1.0/mipToADC_[adet];
}

double HGCSSDigitisation::sumBins(const std::vector<TH2D *> & aHistVec,
				  const double & aMipThresh)
{
  double energy = 0;
  for (unsigned iL(0); iL<aHistVec.size();++iL){
    for (int ix(1); ix<aHistVec[iL]->GetNbinsX()+1; ++ix){
      for (int iy(1); iy<aHistVec[iL]->GetNbinsY()+1; ++iy){
	double eTmp = aHistVec[iL]->GetBinContent(ix,iy);
	if (eTmp > aMipThresh) energy+=eTmp;
      }
    }
  }
  return energy;
}
void HGCSSDigitisation::Print(std::ostream & aOs) const{
  aOs << "====================================" << std::endl
      << " = Random seed: " << rndm_.GetSeed() << std::endl
      << " = Nphoto-electrons: " << npe_ << std::endl
      << " = cross-talk: " << crossTalk_ << std::endl
      << " = Npixels total: " << nTotal_ << std::endl
      << " = sigmaPixel: " << sigmaPix_ << std::endl
    //<< " = sigmaNoise: ECAL " << noise_.find(DetectorEnum::ECAL)->second << ", FHCAL " << noise_.find(DetectorEnum::FHCAL)->second << ", BHCAL " << noise_.find(DetectorEnum::BHCAL)->second << std::endl
      << " = MIPtoADC conversions: ECAL " << mipToADC_.find(DetectorEnum::FECAL)->second << ", FHCAL " << mipToADC_.find(DetectorEnum::FHCAL)->second << std::endl
      << " = Time cut: ECAL " << timeCut_.find(DetectorEnum::FECAL)->second << ", FHCAL " << timeCut_.find(DetectorEnum::FHCAL)->second << ", BHCAL " << timeCut_.find(DetectorEnum::BHCAL1)->second << std::endl
      << "====================================" << std::endl;
};
