#ifndef HGCSSCalibration_h
#define HGCSSCalibration_h


#include <string>
#include <vector>
#include "TH2D.h"

#include "HGCSSDetector.hh"

class HGCSSCalibration {

public:
  HGCSSCalibration(){
    vtx_x_ = 0;
    vtx_y_ = 0;
    vtx_z_ = 0;
  };

  HGCSSCalibration(std::string filePath);
  ~HGCSSCalibration();

  double addTimeOfFlight(const double & aTime,
			 const double & posx,
			 const double & posy,
			 const double & posz,
			 const double & vtxx=0,
			 const double & vtxy=0,
			 const double & vtxz=0);

  double correctTime(const double & aTime,
		     const double & posx,
		     const double & posy,
		     const double & posz,
		     const double & vtxx=0,
		     const double & vtxy=0,
		     const double & vtxz=0);

  double MeVToMip(const unsigned layer,
		  const bool absWeight=false) const;

  //double MeVToMip(const unsigned layer, const double aEta,
  //const bool absWeight=false) const;

  double MeVToMip(const unsigned layer, const double aRadius,
		  const bool absWeight=false) const;

private:

  double vtx_x_;
  double vtx_y_;
  double vtx_z_;

  bool isHCALonly_;
  bool isCaliceHcal_;

};



#endif









