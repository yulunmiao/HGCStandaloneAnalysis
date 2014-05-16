#ifndef HGCSSDetector_h
#define HGCSSDetector_h

#include <string>
#include <vector>
#include "TH2D.h"

enum DetectorEnum {
  ECAL,
  FECAL,
  MECAL,
  BECAL,
  FHCAL,
  BHCAL,
  BHCAL1,
  BHCAL2
};

class HGCSSDetector {

public:
  HGCSSDetector(){};

  ~HGCSSDetector(){
    subdets_.clear();
    enumMap_.clear();
    indices_.clear();
    section_.clear();
  };
  
  const HGCSSSubDetector & subDetector(const unsigned aLayer){
    unsigned section = getSection(aLayer);
    return subdets_[section];
  };

  inline unsigned getSection(const unsigned aLayer) const{
    if (aLayer>=nLayers_) {
      std::cerr << " -- Error ! Trying to access layer " << aLayer 
		<< " outside of range. nLayers = " << nLayers_
		<< std::endl;
      exit(1);
    }
    return section_[aLayer];
  };

  inline const HGCSSSubDetector & getSubDetector(const unsigned index){
    return subdets_[index];
  };

  inline void addSubdetector(const HGCSSSubDetector & adet, const bool is_last){
    subdets_.push_back(adet);
    enumMap_[adet.type]=subdets_.size()-1;
    indices_.push_back(adet.layerIdMin);
    if (is_last) {
      nSections_ = subdets_.size();
      indices_.push_back(adet.layerIdMax);
      nLayers_ = adet.layerIdMax;
      //initialise layer-section conversion
      unsigned lastEle = indices_.size()-1;
      section_.resize(nLayers_,0);
      for (unsigned iL(0); iL<nLayers_;++iL){
	for (unsigned i(0); i<lastEle;++i){
	  if (iL >= indices_[i] && iL < indices_[i+1]) section_[iL] = i;
	}
      }
    }
  };

  inline DetectorEnum detType(const unsigned aSection) const{
    return subdets_[aSection].type;
  };

  inline DetectorEnum detType(const unsigned aLayer) const{
    return subdets_[getSection(aLayer)].type;
  };

  inline std::string detName(const unsigned aSection) const{
    return subdets_[aSection].name;
  };

  inline unsigned nLayers(const unsigned aSection) const{
    return subdets_[aSection].nLayers();
  };

  inline unsigned nLayers(const DetectorEnum & adet) const{
    return subdets_[enumMap_[adet]].nLayers();
  };

  const HGCSSSubDetector & subDetector(const DetectorEnum & adet){
    if (enumMap_.find(adet) == enumMap_.end()){
      std::cerr << " -- Error ! Trying to access subdetector enum not present in this detector: "
		<< adet 
		<< std::endl;
      exit(1);
    } 
    return subdets_[enumMap_[adet]];
  };

  inline unsigned nLayers() const{
    return nLayers_;
  };

  inline unsigned nSections() const{
    return nSections_;
  };

private:
  std::vector<HGCSSSubDetector> subdets_;
  std::vector<unsigned> indices_;
  std::vector<unsigned> section_;
  std::map<DetectorEnum,unsigned> enumMap_;

  unsigned nLayers_;
  unsigned nSections_;
};


class HGCSSSubDetector {

public:
  HGCSSSubDetector():
    type(DetectorEnum::FECAL),
    name(""),
    layerIdMin(0),
    layerIdMax(0),
    mipWeight(1),
    absWeight(1),
    gevWeight(1),
    gevOffset(0),
    isSi(false),
    isScint(false)
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

  inline unsigned nLayers() const{
    return (layerIdMax-layerIdMin);
  };

private:

};



#endif
