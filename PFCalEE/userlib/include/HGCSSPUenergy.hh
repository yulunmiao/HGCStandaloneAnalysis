#ifndef HGCSSPUenergy_h
#define HGCSSPUenergy_h

#include <string>
#include <fstream>
#include <vector>
#include "TF1.h"
#include "TMath.h"

class HGCSSPUenergy{

public:
    HGCSSPUenergy(std::string filePath);
    ~HGCSSPUenergy(); 
    double getDensity(double eta, int layer, double cellSize, int PU);

private:
    std::vector<double> p0_;
    std::vector<double> p1_;
};

#endif 
