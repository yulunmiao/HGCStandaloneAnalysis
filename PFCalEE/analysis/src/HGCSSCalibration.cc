#include "HGCSSCalibration.hh"


HGCSSCalibration::HGCSSCalibration(std::string filePath,const bool concept, const bool calibrate){


  HcalToEcalConv_ = 1;//e/pi....1/1.772;
  HcalToEcalConv_offset_ = 0;//e/pi....1/1.772;
  BHcalToFHcalConv_  = 1;//1/1.33;//1/5.12;

  concept_ = concept;

  bool isScintOnly =  filePath.find("version22")!=filePath.npos;
  isHCALonly_ = filePath.find("version21")!=filePath.npos || isScintOnly;
  isCaliceHcal_ = filePath.find("version23")!=filePath.npos || filePath.find("version_23")!=filePath.npos;


  mipWeights_.resize(6,1);
  gevWeights_.resize(3,1);
  gevOffsets_.resize(3,0);

  //fill weights for : abs width + mip
  //ecal
  mipWeights_[0] = 1./0.0548;//mip
  mipWeights_[1] = 1./0.0548*8.001/5.848;//mip times ratio of abs dedx
  mipWeights_[2] = 1./0.0548*10.854/5.848;

  //HCAL Si or Calice Sci
  if (isCaliceHcal_) {
    mipWeights_[3] = 1./0.807;
    mipWeights_[4] = 1./0.807;
    mipWeights_[5] = 1./0.807*104/21.;//ratio of abs width
    if (calibrate) {
      gevWeights_[1] = 1./41.98;//MIPtoGeV
      gevWeights_[2] = 1./41.98;//MIPtoGeV
      gevOffsets_[1] = -1.7/41.98;//offset in GeV
    }
    HcalToEcalConv_ = 1./0.92;//e/pi
    HcalToEcalConv_offset_ = -1.06;//e/pi
    BHcalToFHcalConv_  = 1;
  }
  else {
    mipWeights_[3] = 1./0.0849*65.235/5.848;
    if (!concept) mipWeights_[3] = 1./0.0849*0.5*65.235/5.848;
    mipWeights_[4] = 1./1.49*92.196/5.848;
    mipWeights_[5] = 1;
  }

  indices_.resize(7,0);

  //fill layer indices
  if (isScintOnly) {
    indices_[4] = 0;
    indices_[5] = 9;
    indices_[6] = 9;
  }
  else if (isCaliceHcal_) {
    indices_[3] = 0;
    indices_[4] = 38;
    indices_[5] = 47;
    indices_[6] = 54;
  }
  else if (isHCALonly_) {
    indices_[3] = 0;
    indices_[4] = 24;
    indices_[5] = 33;
    indices_[6] = 33;
  }
  else {
    indices_[0] = 0;
    indices_[1] = 11;
    indices_[2] = 21;
    indices_[3] = 31;
    indices_[4] = 55;
    indices_[5] = 64;
    indices_[6] = 64;
  }

  //initialise layer-section conversion
  unsigned lastEle = indices_.size()-1;
  index_.resize(indices_[lastEle],0);
  for (unsigned iL(0); iL<indices_[lastEle];++iL){
    for (unsigned i(0); i<lastEle;++i){
      if (iL >= indices_[i] && iL < indices_[i+1]) index_[iL] = i;
    }
  }
  
}

double HGCSSCalibration::MeVToMip(const unsigned layer) const{
  if (layer < index_.size())
    return mipWeights_[index_[layer]];
  return 1;
}

double HGCSSCalibration::mipWeight(const unsigned index) const{
  if (index < mipWeights_.size()) 
    return mipWeights_[index];
  return 1;
}

double HGCSSCalibration::gevWeight(const unsigned section) const{
  if (section < gevWeights_.size()) 
    return gevWeights_[section];
  return 1;
}

double HGCSSCalibration::gevOffset(const unsigned section) const{
  if (section < gevOffsets_.size()) 
    return gevOffsets_[section];
  return 0;
}

double HGCSSCalibration::hcalShowerEnergy(const double Cglobal,
					  const double E_1,
					  const double E_2,
					  const bool correctLinearity){
  double a1 = 1.12;
  double a2 = 0.003;
  double a3 = 0;//1.2e-5;
  double Esh = Cglobal*(E_1*gevWeight(1)-gevOffset(1)) + BHcalToFHcalConv_*(E_2*gevWeight(2)-gevOffset(2));
  if (!correctLinearity) return Esh;
  return Esh*(a1+a2*Esh+a3*Esh*Esh);
}

double HGCSSCalibration::showerEnergy(const double Cglobal,
				      const double E_0,
				      const double E_1,
				      const double E_2,
				      const bool correctLinearity){
  return (E_0*gevWeight(0)-gevOffset(0)) + hcalShowerEnergy(Cglobal,E_1,E_2,correctLinearity);
}

double HGCSSCalibration::recoEnergyUncor(const double E_0,
					 const double E_1,
					 const double E_2,
					 bool FHCALonly){
  if (FHCALonly) return HcalToEcalConv_*(E_1*gevWeight(1)-gevOffset(1)-HcalToEcalConv_offset_);
  return (E_0*gevWeight(0)-gevOffset(0)) + HcalToEcalConv_*((E_1*gevWeight(1)-gevOffset(1)) + BHcalToFHcalConv_*(E_2*gevWeight(2)-gevOffset(2))-HcalToEcalConv_offset_);
}

void HGCSSCalibration::incrementEnergy(const unsigned layer,
				       const double weightedE,
				       double & E_0,
				       double & E_1,
				       double & E_2){
  if (layer<indices_[3]) {
    E_0 += weightedE;
  }//ECAL
  else {
    if (layer<indices_[4]) {
      if (isCaliceHcal_ ||
	  (!concept_ ||
	   (concept_ && layer%2==1))
	  ){
	  E_1 += weightedE;
      }
    }
    else {
      E_2 += weightedE;
    }
  }//HCAL
  
}


