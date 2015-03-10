#ifndef HGCSSCalibration_h
#define HGCSSCalibration_h


#include <string>
#include <vector>
#include "TH2D.h"

#include "HGCSSDetector.hh"

class HGCSSCalibration {

public:
  HGCSSCalibration(std::string filePath);
  ~HGCSSCalibration();

  double correctTime(const double & aTime,
		     const double & posx,
		     const double & posy,
		     const double & posz);

  double MeVToMip(const unsigned layer,
		  const bool absWeight=false) const;

  //double MeVToMip(const unsigned layer, const double aEta,
  //const bool absWeight=false) const;

  double MeVToMip(const unsigned layer, const double aRadius,
		  const bool absWeight=false) const;

private:
  HGCSSCalibration(){};

  double vtx_x_;
  double vtx_y_;
  double vtx_z_;

  bool isHCALonly_;
  bool isCaliceHcal_;

};



#endif









