#include "HGCSSSimHit.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSSimHit::HGCSSSimHit(const G4SiHit & aSiHit,
			 const unsigned & cellid){
  energy_ = aSiHit.energy;
  time_ = aSiHit.time;
  layer_ = aSiHit.layer;
  cellid_ = cellid;
  nGammas_= 0;
  nElectrons_ = 0;
  nMuons_ = 0;
  nHadrons_ = 0;
  if(abs(aSiHit.pdgId)==22) nGammas_++;
  else if(abs(aSiHit.pdgId)==11) nElectrons_++;
  else if(abs(aSiHit.pdgId)==13) nMuons_++;
  else nHadrons_++;
}

void HGCSSSimHit::Add(const G4SiHit & aSiHit){

  time_ = (time_*energy_ + aSiHit.time*aSiHit.energy)/(energy_+aSiHit.energy);

  if(abs(aSiHit.pdgId)==22) nGammas_++;
  else if(abs(aSiHit.pdgId)==11) nElectrons_++;
  else if(abs(aSiHit.pdgId)==13) nMuons_++;
  else nHadrons_++;

  energy_ += aSiHit.energy;
}

void HGCSSSimHit::Print(std::ostream & aOs){
  aOs << "====================================" << std::endl
      << " = Layer " << layer_ << " cellid " << cellid_ << std::endl
      << " = Energy " << energy_ << " time " << time_ << std::endl
      << " = g " << nGammas_ << " e " << nElectrons_ << " mu " << nMuons_ << " had " << nHadrons_ << std::endl
      << "====================================" << std::endl;

}
