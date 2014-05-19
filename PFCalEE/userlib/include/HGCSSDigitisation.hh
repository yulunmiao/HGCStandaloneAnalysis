#ifndef HGCSSDigitisation_h
#define HGCSSDigitisation_h


#include <string>
#include <vector>
#include <map>
#include "TRandom3.h"
#include "TH2D.h"
#include "HGCSSDetector.hh"

class HGCSSDigitisation {

public:

  HGCSSDigitisation():
    seed_(0),
    npe_(11),
    crossTalk_(0.25),
    nTotal_(1156),
    sigmaPix_(3)
  {
    rndm_.SetSeed(seed_);
    //noise_[DetectorEnum::ECAL] = 0.12;
    //noise_[DetectorEnum::FHCAL] = 0.12;
    //noise_[DetectorEnum::BHCAL] = 0.12;
    mipToADC_[DetectorEnum::FECAL] = 50;
    mipToADC_[DetectorEnum::MECAL] = 50;
    mipToADC_[DetectorEnum::BECAL] = 50;
    mipToADC_[DetectorEnum::FHCAL] = 50;
    mipToADC_[DetectorEnum::BHCAL1] = 50;
    mipToADC_[DetectorEnum::BHCAL2] = 50;
    timeCut_[DetectorEnum::FECAL] = 200;//ns
    timeCut_[DetectorEnum::MECAL] = 200;//ns
    timeCut_[DetectorEnum::BECAL] = 200;//ns
    timeCut_[DetectorEnum::FHCAL] = 200;//ns
    timeCut_[DetectorEnum::BHCAL1] = 200;//ns
    timeCut_[DetectorEnum::BHCAL2] = 200;//ns
  };

  ~HGCSSDigitisation(){};

  inline void setRandomSeed(const unsigned aSeed){
    seed_ = aSeed;
    rndm_.SetSeed(seed_);
  };

  inline void setNpe(const unsigned aNpe){
    npe_ = aNpe;
  };

  inline void setCrossTalk(const double & aCrossTalk){
    crossTalk_ = aCrossTalk;
  };

  inline void setNTotalPixels(const unsigned & aN){
    nTotal_ = aN;
  };

  inline void setSigmaPix(const unsigned aSigma){
    sigmaPix_ = aSigma;
  };

  inline void setNoise(const double & aNoise, const unsigned & alay){
    noise_[alay] = aNoise;
  };

  inline void setMipToADC(const double & aMipToADC, DetectorEnum adet){
    mipToADC_[adet] = aMipToADC;
  };

  inline void setTimeCut(const double & aTimeCut, DetectorEnum adet){
    timeCut_[adet] = aTimeCut;
  };

  inline bool passTimeCut(const double aTime, DetectorEnum adet){
    return (aTime < timeCut_[adet]);
  };

  unsigned nRandomPhotoElec(const double & aMipE);

  double nPixels(const double & aMipE);

  double positiveRandomGaus(const double & mean);

  double mipCor(const double & aMipE,
		const double & posx, 
		const double & posy,
		const double & posz);

  double digiE(const double & aMipE);

  void addNoise(double & aDigiE, const unsigned & alay, TH1F * & hist);
  
  unsigned adcConverter(double eMIP, DetectorEnum adet);

  double adcToMIP(const unsigned acdCounts, DetectorEnum adet);

  double sumBins(const std::vector<TH2D *> & aHistVec,
		 const double & aMipThresh);

  void Print(std::ostream & aOs) const;

private:
  TRandom3 rndm_;
  unsigned seed_;
  unsigned npe_;
  double crossTalk_;
  unsigned nTotal_;
  double sigmaPix_;
  std::map<DetectorEnum,double> timeCut_;
  std::map<unsigned,double> noise_;
  std::map<DetectorEnum,unsigned> mipToADC_;

};

#endif
