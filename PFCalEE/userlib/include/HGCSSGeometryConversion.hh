#ifndef HGCSSGeometryConversion_h
#define HGCSSGeometryConversion_h


#include <string>
#include <vector>
#include "TH2D.h"

#include "HGCSSDetector.hh"

class HGCSSGeometryConversion{
  
public:
  HGCSSGeometryConversion(std::string filePath,
			  const HGCSSDetector & aDet);

  ~HGCSSGeometryConversion();

  void setGranularity(const std::vector<unsigned> & granul);

  unsigned getGranularity(const unsigned aLayer, const DetectorEnum adet);

  inline double getXYwidth() const {
    return width_;
  };
  
  inline void cellSize(const double & asize){
    cellSize_ = asize;
  };

  inline double cellSize() const{
    return cellSize_;
  };

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
		   std::string aString,
		   const HGCSSSubDetector & aDet,
		   const unsigned nLayers,
		   bool recreate=false);


  void deleteHistos(std::vector<TH2D *> & aVec);

  TH2F * & get2DHist(const unsigned layer,std::string name);

  inline std::vector<TH2D *> & get2DEnergyVec(const DetectorEnum aDet){
    return 2DHistMapE_[aDet];
  };

  inline std::vector<TH2D *> & get2DTimeVec(const DetectorEnum aDet){
    return 2DHistMapTime_[aDet];
  };

  inline std::vector<TH2D *> & get2DZposVec(const DetectorEnum aDet){
    return 2DHistMapZ_[aDet];
  };

private:

  double width_;
  double cellSize_;
  HGCSSDetector detector_;

  std::map<DetectorEnum,std::vector<TH2D *> > 2DHistMapE_;
  std::map<DetectorEnum,std::vector<TH2D *> > 2DHistMapTime_;
  std::map<DetectorEnum,std::vector<TH2D *> > 2DHistMapZ_;
  std::map<DetectorEnum,std::vector<double> > avgMapZ_;
  std::map<DetectorEnum,std::vector<double> > avgMapE_;
};



#endif
