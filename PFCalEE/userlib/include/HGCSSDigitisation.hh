#ifndef HGCSSDigitisation_h
#define HGCSSDigitisation_h


#include <string>
#include <vector>
#include <map>
#include "TRandom3.h"
#include "TH2D.h"
#include "HGCSSSubDetector.hh"

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
    //noise_[HGCSSDetector::ECAL] = 0.12;
    //noise_[HGCSSDetector::FHCAL] = 0.12;
    //noise_[HGCSSDetector::BHCAL] = 0.12;
    mipToADC_[HGCSSDetector::ECAL] = 50;
    mipToADC_[HGCSSDetector::FHCAL] = 50;
    mipToADC_[HGCSSDetector::BHCAL] = 50;
    timeCut_[HGCSSDetector::ECAL] = 200;//ns
    timeCut_[HGCSSDetector::FHCAL] = 200;//ns
    timeCut_[HGCSSDetector::BHCAL] = 200;//ns
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

  inline void setMipToADC(const double & aMipToADC, const HGCSSDetector & adet){
    mipToADC_[adet] = aMipToADC;
  };

  inline void setTimeCut(const double & aTimeCut, const HGCSSDetector & adet){
    timeCut_[adet] = aTimeCut;
  };

  inline bool passTimeCut(const double aTime, const HGCSSDetector & adet){
    return (aTime < timeCut_[adet]);
  };

  unsigned nRandomPhotoElec(const double & aMipE);

  double nPixels(const double & aMipE);

  double positiveRandomGaus(const double & mean);

  double digiE(const double & aMipE);

  void addNoise(double & aDigiE, const HGCSSDetector & adet, TH1F * & hist);

  void digitiseSiCAL(std::vector<TH2D *> & aHistVec,const HGCSSDetector & adet);

  double adcConverter(const double & eMIP, const HGCSSDetector & adet);

  void digitiseSciCAL(std::vector<TH2D *> & aHistVec,const HGCSSDetector & adet);

  void Print(std::ostream & aOs) const;

private:
  TRandom3 rndm_;
  unsigned seed_;
  unsigned npe_;
  double crossTalk_;
  unsigned nTotal_;
  double sigmaPix_;
  std::map<HGCSSDetector,double> timeCut_;
  std::map<unsigned,double> noise_;
  std::map<HGCSSDetector,unsigned> mipToADC_;

};

#endif
