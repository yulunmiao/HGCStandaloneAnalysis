#ifndef HGCSSGeometryConversion_h
#define HGCSSGeometryConversion_h


#include <string>
#include <vector>
#include <cmath>
#include "TH2D.h"

#include "HGCSSDetector.hh"

class HGCSSGeometryConversion{
  
public:
  HGCSSGeometryConversion(){};
  HGCSSGeometryConversion(std::string filePath,int model, double cellsize);

  ~HGCSSGeometryConversion();

  void setGranularity(const std::vector<unsigned> & granul);

  unsigned getGranularity(const unsigned aLayer, const HGCSSSubDetector & adet);

  inline double getXYwidth() const {
    return width_;
  };
  
  inline void setXYwidth(double width) {
    width_ = width;
  };
  
  inline double cellSize() const{
    return cellSize_;
  };

  //hardcode fine granularity at high eta ?
  inline double cellSize(const unsigned aLayer, const double aEta) const{
    if (fabs(aEta)<10) 
      return cellSize_*granularity_[aLayer];
    return cellSize_*3;
  };

  inline double cellSizeInCm(const unsigned aLayer, const double aEta) const{
    return cellSize(aLayer, aEta)/10.;
  };

  //unsigned getNumberOfSiLayers(const DetectorEnum type,
  //const double & eta=0) const;
  unsigned getNumberOfSiLayers(const DetectorEnum type,
			       const double & radius=10000) const;

  void initialiseHistos(const bool recreate=false,
			std::string uniqStr="",
			const bool print=true);

  void fill(const DetectorEnum type,
	    const unsigned newlayer,
	    const double & weightedE,
	    const double & aTime,
	    const double & posx,
	    const double & posy,
	    const double & posz);

  double getAverageZ(const unsigned layer);

  double sumBins(const std::vector<TH2D *> & aHistVec,
		 const double & aMipThresh=0.);

  void resetVector(std::vector<TH2D *> & aVec,
		   std::string aVar,
		   std::string aString,
		   const HGCSSSubDetector & aDet,
		   const unsigned nLayers,
		   bool recreate=false,
		   bool print=true);


  void deleteHistos(std::vector<TH2D *> & aVec);

  TH2D * get2DHist(const unsigned layer,std::string name);

  inline std::vector<TH2D *> & get2DEnergyVec(const DetectorEnum aDet){
    return HistMapE_[aDet];
  };

  inline std::vector<TH2D *> & get2DTimeVec(const DetectorEnum aDet){
    return HistMapTime_[aDet];
  };

  inline std::vector<TH2D *> & get2DZposVec(const DetectorEnum aDet){
    return HistMapZ_[aDet];
  };

private:

  double width_;
  double cellSize_;
  std::vector<unsigned> granularity_;
  unsigned model_;

  std::map<DetectorEnum,std::vector<TH2D *> > HistMapE_;
  std::map<DetectorEnum,std::vector<TH2D *> > HistMapTime_;
  std::map<DetectorEnum,std::vector<TH2D *> > HistMapZ_;
  std::map<DetectorEnum,std::vector<double> > avgMapZ_;
  std::map<DetectorEnum,std::vector<double> > avgMapE_;
};



#endif
