#ifndef HGCSSCalibration_h
#define HGCSSCalibration_h


#include <string>
#include <vector>
#include "TH2D.h"

class HGCSSCalibration {


public:
  HGCSSCalibration(std::string filePath,const bool concept, const bool calibrate=true);
  ~HGCSSCalibration();

  double correctTime(const double & aTime,
		     const double & posx,
		     const double & posy,
		     const double & posz);

  double MeVToMip(const unsigned layer) const;

  double mipWeight(const unsigned index) const;
  double gevWeight(const unsigned index) const;
  double gevOffset(const unsigned index) const;

  double hcalShowerEnergy(const double Cglobal,
			  const double E_1,
			  const double E_2,
			  const bool correctLinearity=true);
  
  double showerEnergy(const double Cglobal,
		      const double E_0,
		      const double E_1,
		      const double E_2,
		      const bool correctLinearity=true);

  double recoEnergyUncor(const double E_0,
			 const double E_1,
			 const double E_2,
			 bool FHCALonly=false);

  void incrementEnergy(const unsigned layer,
		       const double & weightedE,
		       const double & aTime,
		       const double & posx,
		       const double & posy);

  void incrementBlockEnergy(const unsigned layer,
			    const double & weightedE,
			    double & E_0,
			    double & E_1,
			    double & E_2);

  double sumBins(const std::vector<TH2D *> & aHistVec,
		 const double & aMipThresh=0.);

  double getECALenergy(const double & aMipThresh=0.);
  double getFHCALenergy(const double & aMipThresh=0.);
  double getBHCALenergy(const double & aMipThresh=0.);

  inline bool isHCALonly(){
    return isHCALonly_ || isCaliceHcal_;
  };

  inline double HcalToEcalConv() const{
    return HcalToEcalConv_;
  };

  inline double hcalTimeThreshold() const{
    return hcalTimeThresh_;
  };

  inline void setHcalTimeThreshold(const double aTime){
    hcalTimeThresh_ = aTime;
  };

  inline double BHcalToFHcalConv() const{
    return BHcalToFHcalConv_;
  };


  inline unsigned nEcalLayers() const {
    return indices_[3]-indices_[0];
  }
  
  inline unsigned nFHcalLayers() const {
    return indices_[4]-indices_[3];
  }
  
  inline unsigned nBHcalLayers() const {
    return indices_[6]-indices_[4];
  }
  
  inline unsigned nLayers() const {
    return indices_[6];
  }
  
  void reset2DHistos();

  void resetVector(std::vector<TH2D *> & aVec,
		   std::string aString,
		   const unsigned nLayers,
		   const unsigned nBins,
		   const double & min,
		   const double & max);

  void deleteHistos(std::vector<TH2D *> & aVec);


  inline std::vector<TH2D *> & getECAL2DHistoVec(){
    return E_ECAL_;
  };

  inline std::vector<TH2D *> & getFHCAL2DHistoVec(){
    return E_FHCAL_;
  };

  inline std::vector<TH2D *> & getBHCAL2DHistoVec(){
    return E_BHCAL_;
  };

private:
  HGCSSCalibration(){};

  double hcalTimeThresh_;

  double vtx_x_;
  double vtx_y_;
  double vtx_z_;

  double HcalToEcalConv_;
  double HcalToEcalConv_offset_;
  double BHcalToFHcalConv_;
  bool concept_;
  bool isHCALonly_;
  bool isCaliceHcal_;

  std::vector<unsigned> indices_;
  std::vector<unsigned> index_;
  std::vector<double> mipWeights_;
  std::vector<double> gevWeights_;
  std::vector<double> gevOffsets_;

  std::vector<TH2D *> E_ECAL_;
  std::vector<TH2D *> E_FHCAL_;
  std::vector<TH2D *> E_BHCAL_;


};



#endif









