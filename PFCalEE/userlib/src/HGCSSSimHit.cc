#include "HGCSSSimHit.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSSimHit::HGCSSSimHit(const G4SiHit & aSiHit, const unsigned & asilayer, const float cellSize){
  energy_ = aSiHit.energy;
  //energy weighted time
  //PS: need to call calculateTime() after all hits 
  //have been added to have divided by totalE!!
  time_ = aSiHit.time*aSiHit.energy;
  zpos_ = aSiHit.hit_z;
  setLayer(aSiHit.layer,asilayer);

  //coordinates in mm
  //double z = aSiHit.hit_x;
  double y = aSiHit.hit_y;
  double x = aSiHit.hit_x;
  //cellid encoding:
  bool x_side = x>0 ? true : false;
  bool y_side = y>0 ? true : false;
  unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*getGranularity()));
  unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*getGranularity()));

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

  trackIDMainParent_ = aSiHit.parentId;
  energyMainParent_ = aSiHit.energy;

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

  time_ = time_ + aSiHit.time*aSiHit.energy;
  //PS: need to call calculateTime() after all hits 
  //have been added to have divided by totalE!!

  if(abs(aSiHit.pdgId)==22) nGammas_++;
  else if(abs(aSiHit.pdgId)==11) nElectrons_++;
  else if(abs(aSiHit.pdgId)==13) nMuons_++;
  else if(abs(aSiHit.pdgId)==2112) nNeutrons_++;
  else if(abs(aSiHit.pdgId)==2212) nProtons_++;
  else nHadrons_++;

  energy_ += aSiHit.energy;
  if (aSiHit.energy > energyMainParent_){
    trackIDMainParent_ = aSiHit.parentId;
    energyMainParent_ = aSiHit.energy;
  }

}

/*double HGCSSSimHit::eta() const {
  double x = get_x();
  double y = get_y();
  double theta = acos(fabs(zpos_)/sqrt(zpos_*zpos_+x*x+y*y));
  double leta = -log(tan(theta/2.));
  if (zpos_>0) return leta;
  else return -leta;
  }*/

double HGCSSSimHit::theta() const {
  return 2*atan(exp(-1.*eta()));
}

double HGCSSSimHit::eta() const {
  return position().eta();
}

double HGCSSSimHit::phi() const {
  return position().phi();
}

void HGCSSSimHit::Print(std::ostream & aOs) const{
  aOs << "====================================" << std::endl
      << " = Layer " << layer() << " siLayer " << silayer() << " cellid " << cellid_ << std::endl
      << " = Energy " << energy_ << " time " << time_ << std::endl
      << " = g " << nGammas_ 
      << " e " << nElectrons_ 
      << " mu " << nMuons_ 
      << " neutron " << nNeutrons_ 
      << " proton " << nProtons_ 
      << " had " << nHadrons_ 
      << std::endl
      << " = main parent: trackID " << trackIDMainParent_ << " efrac " << mainParentEfrac()
      << std::endl
      << "====================================" << std::endl;

}

