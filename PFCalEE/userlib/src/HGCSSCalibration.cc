#include "HGCSSCalibration.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSCalibration::HGCSSCalibration(std::string filePath,
				   const bool concept, 
				   const bool calibrate){

  vtx_x_ = 0;
  vtx_y_ = 0;
  vtx_z_ = 0;
  
  concept_ = concept;

  bool isScintOnly =  filePath.find("version22")!=filePath.npos;
  isHCALonly_ = filePath.find("version21")!=filePath.npos || isScintOnly;
  isCaliceHcal_ = filePath.find("version23")!=filePath.npos || filePath.find("version_23")!=filePath.npos;
  bool isECALHCAL = filePath.find("version20")!=filePath.npos;

  if (isCaliceHcal_) vtx_z_ = -1091.75;//mm
  else if (isScintOnly) vtx_z_ = -430.65;//mm
  else if (isHCALonly_) vtx_z_ = -824.01;//mm
  else if (isECALHCAL) vtx_z_ = -985.875;//mm
  else vtx_z_ = -161.865;//mm

  std::vector<unsigned> indices;
  indices.resize(7,0);

  //fill layer indices
  if (isScintOnly) {
    indices[4] = 0;
    indices[5] = 9;
    indices[6] = 9;
  }
  else if (isCaliceHcal_) {
    indices[3] = 0;
    indices[4] = 38;
    indices[5] = 47;
    indices[6] = 54;
  }
  else if (isHCALonly_) {
    indices[3] = 0;
    indices[4] = 24;
    indices[5] = 33;
    indices[6] = 33;
  }
  else {
    indices[0] = 0;
    indices[1] = 11;
    indices[2] = 21;
    indices[3] = 31;
    indices[4] = 55;
    indices[5] = 64;
    indices[6] = 64;
  }

  HGCSSSubDetector FECAL;
  FECAL.type = DetectorEnum::FECAL;
  FECAL.name = "FECAL";
  FECAL.layerIdMin = indices[0];
  FECAL.layerIdMax = indices[1];
  FECAL.mipWeight = 1./0.0548;//mip
  FECAL.absWeight = 1.;//ratio of abs dedx
  FECAL.gevWeight = 1.0;
  FECAL.gevOffset = 0.0;
  FECAL.isSi = true;
  detector_.addSubdetector(FECAL);

  HGCSSSubDetector MECAL;
  MECAL.type = DetectorEnum::MECAL;
  MECAL.name = "MECAL";
  MECAL.layerIdMin = indices[1];
  MECAL.layerIdMax = indices[2];
  MECAL.mipWeight = 1./0.0548;//mip
  MECAL.absWeight = 8.001/5.848;//ratio of abs dedx
  MECAL.gevWeight = 1.0;
  MECAL.gevOffset = 0.0;
  MECAL.isSi = true;
  detector_.addSubdetector(MECAL);

  HGCSSSubDetector BECAL;
  BECAL.type = DetectorEnum::BECAL;
  BECAL.name = "BECAL";
  BECAL.layerIdMin = indices[2];
  BECAL.layerIdMax = indices[3];
  BECAL.mipWeight = 1./0.0548;//mip
  BECAL.absWeight = 10.854/5.848;//ratio of abs dedx
  BECAL.gevWeight = 1.0;
  BECAL.gevOffset = 0.0;
  BECAL.isSi = true;
  detector_.addSubdetector(BECAL);

  HGCSSSubDetector FHCAL;
  FHCAL.type = DetectorEnum::FHCAL;
  FHCAL.name = "FHCAL";
  FHCAL.layerIdMin = indices[3];
  FHCAL.layerIdMax = indices[4];
  FHCAL.mipWeight = 1./0.0849;
  FHCAL.absWeight = 65.235/5.848;//ratio of abs dedx
  if (!concept) FHCAL.absWeights = 0.5*65.235/5.848;
  FHCAL.gevWeight = 1.;
  FHCAL.gevOffset = 0.;
  FHCAL.isSi = true;
  if (isCaliceHcal_) {
    FHCAL.mipWeight = 1./0.807;
    FHCAL.absWeights = 1.;
    FHCAL.gevWeight = 1./41.69;//MIPtoGeV
    FHCAL.gevOffset = -4.3/41.69;//offset in GeV
    FHCAL.isScint = true;
    FHCAL.isSi = false;
  }
  detector_.addSubdetector(FHCAL);

  HGCSSSubDetector BHCAL;
  BHCAL.type = DetectorEnum::BHCAL;
  BHCAL.name = "BHCAL";
  BHCAL.layerIdMin = indices[4];
  BHCAL.layerIdMax = indices[5];
  BHCAL.mipWeight = 1./1.49;
  BHCAL.absWeight = 92.196/5.848;
  BHCAL.gevWeight = 1.0;
  BHCAL.gevOffset = 0.0;
  BHCAL.isScint = true;

  if (isCaliceHcal_) {
    BHCAL.type = DetectorEnum::BHCAL1;
    BHCAL.name = "BHCAL1";
    BHCAL.mipWeight = 1./0.807;
    BHCAL.absWeights = 1.;
    BHCAL.gevWeight = 1./41.69;//MIPtoGeV
    BHCAL.gevOffset = 0.0;
  }

  HGCSSSubDetector BHCAL2;
  BHCAL2.type = DetectorEnum::BHCAL2;
  BHCAL2.name = "BHCAL2";
  BHCAL2.layerIdMin = indices[5];
  BHCAL2.layerIdMax = indices[6];
  BHCAL2.mipWeight = 1./0.807;
  BHCAL2.absWeight = 104./21.;
  BHCAL2.gevWeight = 1./41.69;//MIPtoGeV
  BHCAL2.gevOffset = 0.0;
  BHCAL2.isScint = true;

  if (isCaliceHcal_){
    detector_.addSubdetector(BHCAL);
    detector_.addSubdetector(BHCAL2,true);
    HcalToEcalConv_ = 1/0.914;//1/0.955;//e/pi
    HcalToEcalConv_offset_ = -1.04;//-0.59;//-0.905;//e/pi
    BHcalToFHcalConv_  = 1;
  }
  else detector_.addSubdetector(BHCAL,true);

}


HGCSSCalibration::~HGCSSCalibration(){
}

double HGCSSCalibration::correctTime(const double & aTime,
				     const double & posx,
				     const double & posy,
				     const double & posz){
  double distance = sqrt(pow(posx-vtx_x_,2)+
			 pow(posy-vtx_y_,2)+
			 pow(posz-vtx_z_,2)
			 );
  double c = 299.792458;//3.e8*1000./1.e9;//in mm / ns...
  double cor = distance/c;
  // if (aTime>0 && cor > aTime) std::cout << " -- Problem ! Time correction is too large ";
  // if (aTime>0) std::cout << " -- hit time,x,y,z,cor = " 
  // 			 << aTime << " " << posx << " " << posy << " " << posz << " " 
  // 			 << cor << std::endl;
  double result = aTime-cor;
  if (result<0) result = 0;
  return result;
}

double HGCSSCalibration::MeVToMip(const unsigned layer, const bool absWeight) const{
  if (layer < detector_.nLayers())
    return detector_.subDetector(layer).mipWeight
      *(absWeight?detector_.subDetector(layer).absWeight : 1.0);
  return 1;
}
