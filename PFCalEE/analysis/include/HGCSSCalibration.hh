#ifndef HGCSSCalibration_h
#define HGCSSCalibration_h


#include <string>
#include <vector>

class HGCSSCalibration {


public:
  HGCSSCalibration(std::string filePath,const bool concept, const bool calibrate=true);
  ~HGCSSCalibration(){};

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
		       const double weightedE,
		       double & E_0,
		       double & E_1,
		       double & E_2);


  inline bool isHCALonly(){
    return isHCALonly_ || isCaliceHcal_;
  };

  inline double HcalToEcalConv() const{
    return HcalToEcalConv_;
  };

  inline double BHcalToFHcalConv() const{
    return BHcalToFHcalConv_;
  };


  inline const unsigned nEcalLayers() const {
    return indices_[3]-indices_[0];
  }
  
  inline const unsigned nFHcalLayers() const {
    return indices_[4]-indices_[3];
  }
  
  inline const unsigned nBHcalLayers() const {
    return indices_[6]-indices_[4];
  }
  
  inline const unsigned nLayers() const {
    return indices_[6];
  }
  
private:
  HGCSSCalibration(){};

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


};



#endif









