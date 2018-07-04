#ifndef HGCSSGeometryConversion_h
#define HGCSSGeometryConversion_h


#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "TH2D.h"
#include "TH2Poly.h"
#include "TMath.h"
#include "HGCSSDetector.hh"

struct MergeCells {
  double energy;
  double time;
  //double z;
};


class HGCSSGeometryConversion{
  
public:
  HGCSSGeometryConversion(){};
  HGCSSGeometryConversion(const unsigned model, const double cellsize, const bool bypassR=false, const unsigned nSiLayers=3);

  ~HGCSSGeometryConversion();

  static void convertFromEtaPhi(std::pair<double,double> & xy, const double & z);

  inline void setVersion(const unsigned aV){
    version_ = aV;
  };

  inline TH2Poly *hexagonMap(){
    static TH2Poly hc;
    return &hc;
  };

  inline TH2Poly *diamondMap(){
    static TH2Poly hc;
    return &hc;
  };

  inline TH2Poly *triangleMap(){
    static TH2Poly hc;
    return &hc;
  };

  inline TH2Poly *squareMap(){
    static TH2Poly hsq;
    return &hsq;
  };

  inline TH2Poly *squareMap1(){
    static TH2Poly hsq1;
    return &hsq1;
  };

  inline TH2Poly *squareMap2(){
    static TH2Poly hsq2;
    return &hsq2;
  };

  inline TH2Poly *hexagonMap(TH2Poly & hc){
    return &hc;
  };

  inline TH2Poly *diamondMap(TH2Poly & hc){
    return &hc;
  };

  inline TH2Poly *triangleMap(TH2Poly & hc){
    return &hc;
  };

  inline TH2Poly *squareMap(TH2Poly & hc){
    return &hc;
  };

  inline TH2Poly *squareMap1(TH2Poly & hc){
    return &hc;
  };

  inline TH2Poly *squareMap2(TH2Poly & hc){
    return &hc;
  };

  std::map<int,std::pair<double,double> > hexaGeom;
  std::map<int,std::pair<double,double> > diamGeom;
  std::map<int,std::pair<double,double> > triangleGeom;
  std::map<int,std::pair<double,double> > squareGeom;
  std::map<int,std::pair<double,double> > squareGeom1;
  std::map<int,std::pair<double,double> > squareGeom2;

  inline void copyhexaGeom(const std::map<int,std::pair<double,double> > & ageom) {
    hexaGeom = ageom;
  };

  inline void copydiamGeom(const std::map<int,std::pair<double,double> > &ageom){
    diamGeom = ageom;
  };

  inline void copytriangleGeom(const std::map<int,std::pair<double,double> > &ageom){
    triangleGeom = ageom;
  };

  inline void copysquareGeom(const std::map<int,std::pair<double,double> > &ageom){
    squareGeom = ageom;
  };

  inline void copysquareGeom1(const std::map<int,std::pair<double,double> > &ageom){
    squareGeom1 = ageom;
  };

  inline void copysquareGeom2(const std::map<int,std::pair<double,double> > &ageom){
    squareGeom2 = ageom;
  };

  void initialiseSquareMap(const double xymin, const double side);
  void initialiseSquareMap1(const double xmin, const double xmax, const double ymin, const double ymax, const double side);
  void initialiseSquareMap2(const double xmin, const double xmax, const double ymin, const double ymax, const double side);

  void initialiseSquareMap(TH2Poly *map, const double xymin, const double side, bool print);
  void initialiseSquareMap(TH2Poly *map, const double xmin, const double xmax, const double ymin, const double ymax, const double side, bool print);

  void initialiseDiamondMap(const double & xmin, const double & ymin, const double side);

  void initialiseDiamondMap(const double & xymin, const double side);

  void initialiseDiamondMap(TH2Poly *map, const double & xmin, const double & ymin, const double side, const unsigned nhexa, bool print);

  void initialiseTriangleMap(const double xymin, const double side);

  void initialiseTriangleMap(TH2Poly *map, const double xymin, const double side, bool print);

  void initialiseHoneyComb(const double xymin, const double side);
  void initialiseHoneyComb(const double xymin, const double side, double & xstart, double & ystart);

  void initialiseHoneyComb(TH2Poly *map, const double xymin, const double side, bool print);

  void initialiseHoneyComb(TH2Poly *map, const double xymin, const double side, bool print, double & xstart, double & ystart);

  void fillXY(TH2Poly* hist, std::map<int,std::pair<double,double> > & geom);

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

  inline void setCellSize(double size){
    cellSize_=size;
  };

  //hardcode fine granularity at high eta ?
  //  inline double cellSize(const unsigned aLayer, const double aEta) const{
  //  if (fabs(aEta)<10) 
  //    return cellSize_*granularity_[aLayer];
  //  return cellSize_*3;
  //};
  //  inline double cellSizeInCm(const unsigned aLayer, const double aEta) const{
  // return cellSize(aLayer, aEta)/10.;
  //};
  double cellSize(const unsigned aLayer, const double aR) const;

  double cellSizeInCm(const unsigned aLayer, const double aR) const;

  //unsigned getNumberOfSiLayers(const DetectorEnum type,
  //const double & eta=0) const;
  unsigned getNumberOfSiLayers(const DetectorEnum type,
			       const double & radius=10000) const;

  void initialiseHistos(const bool recreate=false,
			std::string uniqStr="",
			const bool print=true);

  //void fill(const DetectorEnum type,
  //const unsigned newlayer,
  //const double & weightedE,
  //const double & aTime,
  //const double & posx,
  //const double & posy,
  //const double & posz);

  void fill(const unsigned layer,
	    const double & weightedE,
	    const double & aTime,
	    const unsigned & cellid,
	    const double & posz);


  double getAverageZ(const unsigned layer);

  //with TH2Poly
  /*double sumBins(const std::vector<TH2Poly *> & aHistVec,
		 const double & aMipThresh=0.);

  void resetVector(std::vector<TH2Poly *> & aVec,
		   std::string aVar,
		   std::string aString,
		   const HGCSSSubDetector & aDet,
		   const unsigned nLayers,
		   bool recreate=false,
		   bool print=true);


  void deleteHistos(std::vector<TH2Poly *> & aVec);

  TH2Poly * get2DHist(const unsigned layer,std::string name);

  inline std::vector<TH2Poly *> & get2DEnergyVec(const DetectorEnum aDet){
    return HistMapE_[aDet];
  };

  inline std::vector<TH2Poly *> & get2DTimeVec(const DetectorEnum aDet){
    return HistMapTime_[aDet];
  };

  inline std::vector<TH2Poly *> & get2DZposVec(const DetectorEnum aDet){
    return HistMapZ_[aDet];
  };
  */

  //with maps of MergeCells struct
  double sumBins(const std::map<unsigned,MergeCells> & aHistVec,
		 const double & aMipThresh=0.);

  void resetVector(std::map<unsigned,MergeCells> & aVec);

  void deleteHistos(std::map<unsigned,MergeCells> & aVec);

  std::map<unsigned,MergeCells> & get2DHist(const unsigned layer);

private:

  void myHoneycomb(TH2Poly* map,
		   Double_t xstart,
		   Double_t ystart,
		   Double_t a,  // side length
		   Int_t k,     // # hexagons in a column
		   Int_t s);    // # columns
  
  bool dopatch_;
  double width_;
  double cellSize_;
  std::vector<unsigned> granularity_;
  unsigned model_;
  bool bypassRadius_;
  unsigned nSiLayers_;
  unsigned version_;
  //std::map<DetectorEnum,std::vector<TH2Poly *> > HistMapE_;
  //std::map<DetectorEnum,std::vector<TH2Poly *> > HistMapTime_;
  //std::map<DetectorEnum,std::vector<TH2Poly *> > HistMapZ_;
  //std::map<DetectorEnum,std::vector<double> > avgMapZ_;
  //std::map<DetectorEnum,std::vector<double> > avgMapE_;
  std::map<unsigned,std::map<unsigned,MergeCells> > HistMap_;
  std::map<unsigned,double> avgMapZ_;
  std::map<unsigned,double> avgMapE_;


};



#endif
