#ifndef SignalRegion_h
#define SignalRegion_h

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "HGCSSRecoHit.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSPUenergy.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSCalibration.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"


class SignalRegion{

public:
    SignalRegion(const std::string inputFolder, 
                 const unsigned nLayers, 
                 const unsigned nevt,
                 const HGCSSGeometryConversion & geomConv,
                 const HGCSSPUenergy & puDensity,
		 const bool applyPuMixFix);

    ~SignalRegion();

    void initialise(TTree *aSimTree, TTree *aRecoTree, 
		    TFile *outputFile);
   
    void initialiseHistograms();

    void fillHistograms();

    inline void setOutputFile(TFile *outputFile){
      outputFile_ = outputFile;
    };

  double getEtotalSR(unsigned iSR, bool subtractPU){
    double Etotal(0);
    if (iSR>=nSR_) return 0;
    for(unsigned iL(0);iL<nLayers_;iL++){
      if(subtractPU) Etotal += subtractedenergySR_[iL][iSR];
      else Etotal += energySR_[iL][iSR];
    }
    return Etotal;
  };
  
  double getSR(unsigned iSR, unsigned layer, bool subtractPU){
    if(layer >= nLayers_) return 0;
    if (iSR>=nSR_) return 0;
    if(subtractPU) {return subtractedenergySR_[layer][iSR];}
    else {return energySR_[layer][iSR];}
  };

  double absweight(unsigned layer){
    if(layer >= nLayers_) return 0;
    return absweight_[layer];
  }
  
private:
  
  unsigned nSR_;
  unsigned nevt_;
  unsigned nLayers_;    
  
  TFile *outputFile_;
  TTree *outtree_;
  
  HGCSSGeometryConversion geomConv_;
  HGCSSPUenergy puDensity_;
  //HGCSSCalibration *mycalib_;
  
  bool fixForPuMixBug_;
  
  std::vector<double> zPos_;
  std::vector<std::vector<ROOT::Math::XYZVector> > accuratePos_;
  std::vector<double> absweight_;
  
  //for tree
  double totalE_;
  double wgttotalE_;
  std::vector<std::vector<double> > energySR_;
  std::vector<std::vector<double> > subtractedenergySR_;
  
  TH1F *p_rawEtotal;
  TH1F *p_wgtEtotal;
  std::vector<TH1F*> p_rawESR;
  std::vector<TH1F*> p_wgtESR;
  std::vector<TH1F*> p_rawSubtractESR;
  std::vector<TH1F*> p_wgtSubtractESR;

};

#endif
