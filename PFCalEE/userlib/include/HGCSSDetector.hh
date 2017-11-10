#ifndef HGCSSDetector_h
#define HGCSSDetector_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TH2D.h"

enum DetectorEnum {
  FECAL,
  MECAL,
  BECAL,
  FHCAL,
  BHCAL1,
  BHCAL2
};

class HGCSSSubDetector {

public:
  HGCSSSubDetector():
    type(FECAL),
    name(""),
    layerIdMin(0),
    layerIdMax(0),
    mipWeight(1),
    absWeight(1),
    gevWeight(1),
    gevOffset(0),
    isSi(false),
    isScint(false),
    radiusLim(0)
  {};
  ~HGCSSSubDetector(){};

  DetectorEnum type;
  std::string name;
  unsigned layerIdMin;
  unsigned layerIdMax;
  double mipWeight;
  double absWeight;
  double gevWeight;
  double gevOffset;
  bool isSi;
  bool isScint;
  double radiusLim;

  inline unsigned nLayers() const{
    return (layerIdMax-layerIdMin);
  };

private:

};

class HGCSSDetector {

public:
  friend HGCSSDetector & theDetector();

  inline void initialiseIndices(const unsigned versionNumber){
    
    indices_.clear();
    indices_.resize(7,0);
    //fill layer indices
    if (versionNumber==22){
      indices_[4] = 0;
      indices_[5] = 10;
      indices_[6] = 10;
    }
    else if (versionNumber==28 || versionNumber==32) {
      indices_[4] = 0;
      indices_[5] = 12;
      indices_[6] = 12;
    }
    else if (versionNumber==23) {
      indices_[3] = 0;
      indices_[4] = 38;
      indices_[5] = 47;
      indices_[6] = 54;
    }
    else if (versionNumber==21) {
      indices_[3] = 0;
      indices_[4] = 24;
      indices_[5] = 34;
      indices_[6] = 34;
    }
    else if (versionNumber==27 || versionNumber==31) {
      indices_[3] = 0;
      indices_[4] = 12;
      indices_[5] = 24;
      indices_[6] = 24;
    }
    else if (versionNumber==38) {
      indices_[3] = 0;
      indices_[4] = 11;
      indices_[5] = 23;
      indices_[6] = 23;
    }
    else if (versionNumber==39) {
      indices_[3] = 0;
      indices_[4] = 9;
      indices_[5] = 21;
      indices_[6] = 21;
    }
    else if (versionNumber < 20){
      indices_[0] = 0;
      indices_[1] = versionNumber==8?11:10;
      indices_[2] = versionNumber==8?21:20;
      indices_[3] = versionNumber==8?31:30;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 30 || versionNumber == 60 || (versionNumber >= 100 && versionNumber < 104)){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 33){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = 40;
      indices_[5] = 52;
      indices_[6] = 52;
    }
    else if (versionNumber == 34){
      indices_[0] = 0;
      indices_[1] = 8;
      indices_[2] = 16;
      indices_[3] = 24;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 35){
      indices_[0] = 0;
      indices_[1] = 6;
      indices_[2] = 12;
      indices_[3] = 18;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 36){
      indices_[0] = 0;
      indices_[1] = 8;
      indices_[2] = 16;
      indices_[3] = 24;
      indices_[4] = 35;
      indices_[5] = 47;
      indices_[6] = 47;
    }
    else if (versionNumber == 37){
      indices_[0] = 0;
      indices_[1] = 6;
      indices_[2] = 12;
      indices_[3] = 18;
      indices_[4] = 27;
      indices_[5] = 39;
      indices_[6] = 39;
    }
    else if (versionNumber == 61){
      indices_[3] = 0;
      indices_[4] = 25;
      indices_[5] = 41;
      indices_[6] = 41;
    }
    else if (versionNumber == 62){
      indices_[4] = 0;
      indices_[5] = 16;
      indices_[6] = 16;
    }
    else if (versionNumber == 63){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = 53;
      indices_[5] = 57;
      indices_[6] = 69;
      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3198;
      sensitiveZ_[1] = 3207.1;
      sensitiveZ_[2] = 3222.4;
      sensitiveZ_[3] = 3231.5;
      sensitiveZ_[4] = 3246.8;
      sensitiveZ_[5] = 3255.9;
      sensitiveZ_[6] = 3271.2;
      sensitiveZ_[7] = 3280.3;
      sensitiveZ_[8] = 3295.6;
      sensitiveZ_[9] = 3304.7;
      sensitiveZ_[10] = 3320;
      sensitiveZ_[11] = 3329.1;
      sensitiveZ_[12] = 3344.4;
      sensitiveZ_[13] = 3353.5;
      sensitiveZ_[14] = 3368.8;
      sensitiveZ_[15] = 3377.9;
      sensitiveZ_[16] = 3393.2;
      sensitiveZ_[17] = 3402.3;
      sensitiveZ_[18] = 3417.6;
      sensitiveZ_[19] = 3426.7;
      sensitiveZ_[20] = 3442;
      sensitiveZ_[21] = 3451.1;
      sensitiveZ_[22] = 3466.4;
      sensitiveZ_[23] = 3475.5;
      sensitiveZ_[24] = 3490.8;
      sensitiveZ_[25] = 3499.9;
      sensitiveZ_[26] = 3515.2;
      sensitiveZ_[27] = 3524.3;
      sensitiveZ_[28] = 3577.4;
      sensitiveZ_[29] = 3626.4;
      sensitiveZ_[30] = 3675.4;
      sensitiveZ_[31] = 3724.4;
      sensitiveZ_[32] = 3773.4;
      sensitiveZ_[33] = 3822.4;
      sensitiveZ_[34] = 3871.4;
      sensitiveZ_[35] = 3920.4;
      sensitiveZ_[36] = 3969.4;
      sensitiveZ_[37] = 4020.3;
      sensitiveZ_[38] = 4071.2;
      sensitiveZ_[39] = 4122.1;
      sensitiveZ_[40] = 4206;
      sensitiveZ_[41] = 4289.9;
      sensitiveZ_[42] = 4373.8;
      sensitiveZ_[43] = 4457.7;
      sensitiveZ_[44] = 4541.6;
      sensitiveZ_[45] = 4625.5;
      sensitiveZ_[46] = 4709.4;
      sensitiveZ_[47] = 4793.3;
      sensitiveZ_[48] = 4877.2;
      sensitiveZ_[49] = 4961.1;
      sensitiveZ_[50] = 5045;
      sensitiveZ_[51] = 5128.9;
      sensitiveZ_[52] = 6.91046e-310;
      sensitiveZ_[53] = 3971.2;
      sensitiveZ_[54] = 4022.1;
      sensitiveZ_[55] = 4073;
      sensitiveZ_[56] = 4123.9;
      sensitiveZ_[57] = 4207.8;
      sensitiveZ_[58] = 4291.7;
      sensitiveZ_[59] = 4375.6;
      sensitiveZ_[60] = 4459.5;
      sensitiveZ_[61] = 4543.4;
      sensitiveZ_[62] = 4627.3;
      sensitiveZ_[63] = 4711.2;
      sensitiveZ_[64] = 4795.1;
      sensitiveZ_[65] = 4879;
      sensitiveZ_[66] = 4962.9;
      sensitiveZ_[67] = 5046.8;
      sensitiveZ_[68] = 5130.7;
    }
    else if (versionNumber == 110){
      indices_[0] = 0;
      indices_[1] = 4;
      indices_[2] = indices_[1];
      indices_[3] = indices_[1];
      indices_[4] = indices_[1];
      indices_[5] = indices_[1];
      indices_[6] = indices_[1];
    }
    else {
      indices_[0] = 0;
      indices_[1] = versionNumber==24?11:10;
      indices_[2] = versionNumber==24?21:20;
      indices_[3] = versionNumber==24?31:30;
      indices_[4] = versionNumber==24?55:42;
      indices_[5] = versionNumber==24?65:54;
      indices_[6] = versionNumber==24?65:54;
    }
    
  };

  void buildDetector(const unsigned versionNumber,
		     bool concept=true,
		     bool isCaliceHcal=false,
		     bool bypassR=false);

  const HGCSSSubDetector & subDetectorByLayer(const unsigned aLayer);

  unsigned getSection(const unsigned aLayer) const;
  inline unsigned section(const DetectorEnum adet){
    if (enumMap_.find(adet) != enumMap_.end())
      return enumMap_[adet];
    return nSections();
  };

  void addSubdetector(const HGCSSSubDetector & adet);
  
  void finishInitialisation();

  inline bool isMixedLayer(const unsigned versionNumber,const unsigned aLayer){
    if (versionNumber!=63) return false;
    if (aLayer<36) return false;
    else return true;
  };

  inline unsigned nLayers(const unsigned aSection) const{
    return subdets_[aSection].nLayers();
  };

  inline unsigned nLayers(DetectorEnum adet){
    return subdets_[enumMap_[adet]].nLayers();
  };

  const HGCSSSubDetector & subDetectorByEnum(DetectorEnum adet);
  inline const HGCSSSubDetector & subDetectorBySection(const unsigned aSection) const{
    return subdets_[aSection];
  };

  inline unsigned nLayers() const{
    return nLayers_;
  };

  inline unsigned nSections() const{
    return nSections_;
  };

  inline DetectorEnum detType(const unsigned aSection) const{
    return subdets_[aSection].type;
  };
  
  inline DetectorEnum detTypeLayer(const unsigned aLayer) const{
    return subdets_[getSection(aLayer)].type;
  };
  
  inline std::string detName(const unsigned aSection) const{
    return subdets_[aSection].name;
  };


  void reset();

  void printDetector(std::ostream & aOs) const ;

  inline double sensitiveZ(const unsigned layer){
    if (layer<sensitiveZ_.size()) return sensitiveZ_[layer];
    else {
      std::cout << " ERROR! Trying to access layer " << layer << " outside of range: " << sensitiveZ_.size() << " nLayers " << nLayers_ << std::endl;
      exit(1);
    }
  };

private:
  HGCSSDetector(){
    bypassRadius_ = false;
  };

  ~HGCSSDetector(){
    reset();
  };
  
  std::vector<HGCSSSubDetector> subdets_;
  std::vector<unsigned> indices_;
  std::vector<unsigned> section_;
  std::map<DetectorEnum,unsigned> enumMap_;

  unsigned nLayers_;
  unsigned nSections_;
  bool bypassRadius_;

  std::vector<double> sensitiveZ_;

};

HGCSSDetector & theDetector();




#endif
