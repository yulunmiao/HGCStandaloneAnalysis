#ifndef HGCSSDigitisation_h
#define HGCSSDigitisation_h


#include <string>
#include <vector>
#include <map>
#include "TRandom3.h"
#include "TH2D.h"

class HGCSSDigitisation {

public:

  enum Detector {
    ECAL=0,
    FHCAL=1,
    BHCAL=2
  };

  HGCSSDigitisation():
    seed_(0),
    npe_(11),
    crossTalk_(0.25),
    nTotal_(1156),
    sigmaPix_(3)
  {
    rndm_.SetSeed(seed_);
    noise_[Detector::ECAL] = 0.12;
    noise_[Detector::FHCAL] = 0.12;
    noise_[Detector::BHCAL] = 0.12;
    mipToADC_[Detector::ECAL] = 50;
    mipToADC_[Detector::FHCAL] = 50;

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

  inline void setNoise(const double & aNoise, const Detector & adet){
    noise_[adet] = aNoise;
  };

  unsigned nRandomPhotoElec(const double & aMipE);

  double nPixels(const double & aMipE);

  double positiveRandomGaus(const double & mean);

  double digiE(const double & aMipE);

  void addNoise(double & aDigiE, const Detector & adet);

  void digitiseECAL(std::vector<TH2D *> & aHistVec);

  double adcConverter(const double & eMIP, const Detector & adet);

  void digitiseFHCAL(std::vector<TH2D *> & aHistVec);

  void digitiseBHCAL(std::vector<TH2D *> & aHistVec);

  void Print(std::ostream & aOs) const;

private:
  TRandom3 rndm_;
  unsigned seed_;
  unsigned npe_;
  double crossTalk_;
  unsigned nTotal_;
  double sigmaPix_;
  std::map<Detector,double> noise_;
  std::map<Detector,unsigned> mipToADC_;

};

#endif
