#ifndef SignalRegion_h
#define SignalRegion_h

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"

#include "HGCSSRecoHit.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"


class SignalRegion{

public:
    SignalRegion(const std::string inputFolder, double cellSize, unsigned nLayers, unsigned nevt);
    ~SignalRegion();

    void initialise(std::vector<HGCSSRecoHit> *rechitvec, unsigned ievt);

    double getEtotalSR0(unsigned ievt){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            Etotal += energySR0_[ievt][iL];
        }
        return Etotal;
    }

    double getEtotalSR1(unsigned ievt){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            Etotal += energySR1_[ievt][iL];
        }
        return Etotal;
    }

    double getEtotalSR2(unsigned ievt){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            Etotal += energySR2_[ievt][iL];
        }
        return Etotal;
    }

    double getEtotalSR3(unsigned ievt){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            Etotal += energySR3_[ievt][iL];
        }
        return Etotal;
    }

    double getEtotalSR4(unsigned ievt){
        double Etotal(0);
        for(unsigned iL(0);iL<nLayers_;iL++){
            Etotal += energySR4_[ievt][iL];
        }
        return Etotal;
    }

    double getSR0(unsigned ievt, unsigned layer){
        return energySR0_[ievt][layer];
    }

    double getSR1(unsigned ievt, unsigned layer){
        return energySR1_[ievt][layer];
    }

    double getSR2(unsigned ievt, unsigned layer){
        return energySR2_[ievt][layer];
    }

    double getSR3(unsigned ievt, unsigned layer){
        return energySR3_[ievt][layer];
    }

    double getSR4(unsigned ievt, unsigned layer){
        return energySR4_[ievt][layer];
    }


private:

    unsigned nevt_;
    double cellSize_; 
    unsigned nLayers_;    
    std::vector<double> zPos_;
    std::vector<std::vector<ROOT::Math::XYZVector> > accuratePos_;
    std::vector<std::vector<double> > energySR0_;
    std::vector<std::vector<double> > energySR1_;
    std::vector<std::vector<double> > energySR2_;
    std::vector<std::vector<double> > energySR3_;
    std::vector<std::vector<double> > energySR4_;

};


#endif
