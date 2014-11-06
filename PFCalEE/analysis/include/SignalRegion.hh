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

    double getEtotalSR0(unsigned ievt, bool subtractPU){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            if(subtractPU) Etotal += subtractedenergySR0_[ievt][iL];
            else Etotal += energySR0_[ievt][iL];
        }
        return Etotal;
    }

    double getEtotalSR1(unsigned ievt, bool subtractPU){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            if(subtractPU) Etotal += subtractedenergySR1_[ievt][iL];
            else Etotal += energySR1_[ievt][iL];
        }
        return Etotal;
    }

    double getEtotalSR2(unsigned ievt, bool subtractPU){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            if(subtractPU) Etotal += subtractedenergySR2_[ievt][iL];
            else Etotal += energySR2_[ievt][iL];
        }
        return Etotal;
    }

    double getEtotalSR3(unsigned ievt, bool subtractPU){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            if(subtractPU) Etotal += subtractedenergySR3_[ievt][iL];
            else Etotal += energySR3_[ievt][iL];
        }
        return Etotal;
    }

    double getEtotalSR4(unsigned ievt, bool subtractPU){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            if(subtractPU) Etotal += subtractedenergySR4_[ievt][iL];
            else Etotal += energySR4_[ievt][iL];
        }
        return Etotal;
    }

    double getSR0(unsigned ievt, unsigned layer, bool subtractPU){
        if(subtractPU) {return subtractedenergySR0_[ievt][layer];}
        else {return energySR0_[ievt][layer];}
    }

    double getSR1(unsigned ievt, unsigned layer, bool subtractPU){
        if(subtractPU){ return subtractedenergySR1_[ievt][layer];}
        else{ return energySR1_[ievt][layer];}
    }

    double getSR2(unsigned ievt, unsigned layer, bool subtractPU){
        if(subtractPU){ return subtractedenergySR2_[ievt][layer];}
        else{ return energySR2_[ievt][layer];}
    }

    double getSR3(unsigned ievt, unsigned layer, bool subtractPU){
        if(subtractPU){ return subtractedenergySR3_[ievt][layer];}
        else{ return energySR3_[ievt][layer];}
    }

    double getSR4(unsigned ievt, unsigned layer, bool subtractPU){
        if(subtractPU){ return subtractedenergySR4_[ievt][layer];}
        else{ return energySR4_[ievt][layer];}
    }

    double absweight(unsigned layer){
        if(layer >= nLayers_) return 0;
        return absweight_[layer];
    }

private:

    unsigned nevt_;
    unsigned nLayers_;    

    TFile *outputFile_;

    HGCSSGeometryConversion geomConv_;
    HGCSSPUenergy puDensity_;
  //HGCSSCalibration *mycalib_;

  bool fixForPuMixBug_;

    std::vector<double> zPos_;
    std::vector<std::vector<ROOT::Math::XYZVector> > accuratePos_;
    std::vector<double> absweight_;

    std::vector<double> totalE_;
    std::vector<double> wgttotalE_;

    std::vector<std::vector<double> > energySR0_;
    std::vector<std::vector<double> > energySR1_;
    std::vector<std::vector<double> > energySR2_;
    std::vector<std::vector<double> > energySR3_;
    std::vector<std::vector<double> > energySR4_;

    std::vector<std::vector<double> > subtractedenergySR0_;
    std::vector<std::vector<double> > subtractedenergySR1_;
    std::vector<std::vector<double> > subtractedenergySR2_;
    std::vector<std::vector<double> > subtractedenergySR3_;
    std::vector<std::vector<double> > subtractedenergySR4_;

    TH1F *p_rawEtotal;
    TH1F *p_wgtEtotal;

    TH1F *p_rawESR0;
    TH1F *p_rawESR1;
    TH1F *p_rawESR2;
    TH1F *p_rawESR3;
    TH1F *p_rawESR4;
    TH1F *p_wgtESR0;
    TH1F *p_wgtESR1;
    TH1F *p_wgtESR2;
    TH1F *p_wgtESR3;
    TH1F *p_wgtESR4;

    TH1F *p_rawSubtractESR0;
    TH1F *p_rawSubtractESR1;
    TH1F *p_rawSubtractESR2;
    TH1F *p_rawSubtractESR3;
    TH1F *p_rawSubtractESR4;
    TH1F *p_wgtSubtractESR0;
    TH1F *p_wgtSubtractESR1;
    TH1F *p_wgtSubtractESR2;
    TH1F *p_wgtSubtractESR3;
    TH1F *p_wgtSubtractESR4;
};

#endif
