#ifndef HGCSSCalibration_h
#define HGCSSCalibration_h


#include <string>
#include <vector>
#include "TH2D.h"

class HGCSSCalibration {

public:
  HGCSSCalibration(std::string filePath,
		   const bool concept, 
		   const bool calibrate=true);
  ~HGCSSCalibration();

  double correctTime(const double & aTime,
		     const double & posx,
		     const double & posy,
		     const double & posz);

  double MeVToMip(const unsigned layer,
		  const bool absWeight=false) const;

  inline HGCSSDetector & detector() {
    return detector_;
  };

private:
  HGCSSCalibration(){};

  double vtx_x_;
  double vtx_y_;
  double vtx_z_;

  bool concept_;
  bool isHCALonly_;
  bool isCaliceHcal_;

  HGCSSDetector detector_;


};



#endif









