#include "HGCSSPUenergy.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSPUenergy::HGCSSPUenergy(std::string filePath){
   
    std::string function("");
    std::string dataformat("");
    unsigned layer;
    double parameter0, parameter1; 

    std::ostringstream fileName;
    fileName << filePath; 
    std::ifstream EnergyDensity;
    EnergyDensity.open(fileName.str().c_str(),std::ios::in);
   
    std::getline(EnergyDensity,function);
    std::getline(EnergyDensity,dataformat);
    std::cout << "Energy Density function: " << function << std::endl;
    std::cout << dataformat << std::endl;

    EnergyDensity >> layer >> parameter0 >> parameter1;
    while(!EnergyDensity.eof( )){
       std::cout << "layer " << layer << ": " << parameter0 << " " <<parameter1 << std::endl;
       if(layer==p0_.size()){
           p0_.push_back(parameter0);
           p1_.push_back(parameter1);
       } else {
           std::cout << "error: try to add layer "<< layer <<" to existing " << p0_.size() << "layers" <<std::endl;
       }
       EnergyDensity >> layer >> parameter0 >> parameter1;
    }
}

HGCSSPUenergy::~HGCSSPUenergy(){
}

double HGCSSPUenergy::getDensity(double eta, int layer, double cellSize, int PU){
    
    return TMath::Exp(p0_[layer] + eta*p1_[layer])*cellSize*cellSize*PU;
}

 



 
