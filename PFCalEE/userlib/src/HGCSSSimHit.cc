#include "HGCSSSimHit.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSSimHit::HGCSSSimHit(const G4SiHit & aSiHit){
  energy_ = aSiHit.energy;
  time_ = aSiHit.time;
  zpos_ = aSiHit.hit_z;
  layer_ = aSiHit.layer;

  //coordinates in mm
  //double z = aSiHit.hit_x;
  double y = aSiHit.hit_y;
  double x = aSiHit.hit_x;
  //cellid encoding:
  bool x_side = x>0 ? true : false;
  bool y_side = y>0 ? true : false;
  unsigned x_cell = static_cast<unsigned>(fabs(x)/(CELL_SIZE_X*getGranularity()));
  unsigned y_cell = static_cast<unsigned>(fabs(y)/(CELL_SIZE_Y*getGranularity()));

  encodeCellId(x_side,y_side,x_cell,y_cell);

  nGammas_= 0;
  nElectrons_ = 0;
  nMuons_ = 0;
  nNeutrons_ = 0;
  nProtons_ = 0;
  nHadrons_ = 0;
  if(abs(aSiHit.pdgId)==22) nGammas_++;
  else if(abs(aSiHit.pdgId)==11) nElectrons_++;
  else if(abs(aSiHit.pdgId)==13) nMuons_++;
  else if(abs(aSiHit.pdgId)==2112) nNeutrons_++;
  else if(abs(aSiHit.pdgId)==2212) nProtons_++;
  else nHadrons_++;
}

void HGCSSSimHit::encodeCellId(const bool x_side,const bool y_side,const unsigned x_cell,const unsigned y_cell){
  cellid_ = 
    x_side | (x_cell<<1) |
    (y_side<<16) | (y_cell<<17);

  // std::cout << " Cross-check of encoding: cellid=" << cellid_ << std::endl
  // 	    << " x_side " << x_side << " " << get_x_side() << std::endl
  // 	    << " y_side " << y_side << " " << get_y_side() << std::endl
  // 	    << " x_cell " << x_cell << " " << get_x_cell() << std::endl
  // 	    << " y_cell " << y_cell << " " << get_y_cell() << std::endl
  //   ;
}

void HGCSSSimHit::Add(const G4SiHit & aSiHit){

  time_ = (time_*energy_ + aSiHit.time*aSiHit.energy)/(energy_+aSiHit.energy);

  if(abs(aSiHit.pdgId)==22) nGammas_++;
  else if(abs(aSiHit.pdgId)==11) nElectrons_++;
  else if(abs(aSiHit.pdgId)==13) nMuons_++;
  else if(abs(aSiHit.pdgId)==2112) nNeutrons_++;
  else if(abs(aSiHit.pdgId)==2212) nProtons_++;
  else nHadrons_++;

  energy_ += aSiHit.energy;
}

void HGCSSSimHit::Print(std::ostream & aOs) const{
  aOs << "====================================" << std::endl
      << " = Layer " << layer_ << " cellid " << cellid_ << std::endl
      << " = Energy " << energy_ << " time " << time_ << std::endl
      << " = g " << nGammas_ 
      << " e " << nElectrons_ 
      << " mu " << nMuons_ 
      << " neutron " << nNeutrons_ 
      << " proton " << nProtons_ 
      << " had " << nHadrons_ 
      << std::endl
      << "====================================" << std::endl;

}

void HGCSSSimHit::PrintGeometry(std::ostream & aOs) const{
  aOs << " =========================================================== " << std::endl
      << " ============= Printing of Transverse Geometry ============= " << std::endl
      << " =========================================================== " << std::endl
      << " = CELL_SIZE_X " << CELL_SIZE_X << " mm" << std::endl
      << " = CELL_SIZE_Y " << CELL_SIZE_Y << " mm" << std::endl
      << " =========================================================== " << std::endl
      << " =========================================================== " << std::endl;
}
