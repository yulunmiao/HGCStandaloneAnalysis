#ifndef _hgcsssimhit_hh_
#define _hgcsssimhit_hh_

#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <map>

#include "G4SiHit.hh"

static const float CELL_SIZE_X=2.5;//mm
static const float CELL_SIZE_Y=CELL_SIZE_X;

class HGCSSSimHit{

public:
  HGCSSSimHit():
    energy_(0),
    time_(0),
    zpos_(0),
    layer_(0),
    cellid_(0),
    nGammas_(0),
    nElectrons_(0),
    nMuons_(0),
    nNeutrons_(0),
    nProtons_(0),
    nHadrons_(0),
    trackIDMainParent_(0),
    energyMainParent_(0)
  {
    
  };
  HGCSSSimHit(const G4SiHit & aSiHit);

  ~HGCSSSimHit(){};

  inline double energy() const {
    return energy_;
  };

  inline double time() const {
    return time_;
  };

  inline void calculateTime() {
    if (energy_>0) time_ = time_/energy_;
  };

  inline unsigned layer() const {
    return layer_;
  };

  inline void layer(const unsigned & layer){
    layer_ = layer;
  };

  inline unsigned cellid() const {
    return cellid_;
  };

  inline unsigned nGammas() const {
    return nGammas_;
  };

  inline unsigned nElectrons() const {
    return nElectrons_;
  };

  inline unsigned nMuons() const {
    return nMuons_;
  };

  inline unsigned nNeutrons() const {
    return nNeutrons_;
  };

  inline unsigned nProtons() const {
    return nProtons_;
  };
  inline unsigned nHadrons() const {
    return nHadrons_;
  };
  inline unsigned numberOfParticles() const {
    return nGammas_+nElectrons_+nMuons_+nNeutrons_+nProtons_+nHadrons_;
  };

  inline double gFrac() const {
    return nGammas_/numberOfParticles();
  };

  inline double eFrac() const {
    return nElectrons_/numberOfParticles();
  };

  inline  double muFrac() const {
    return nMuons_/numberOfParticles();
  };

  inline double neutronFrac() const {
    return nNeutrons_/numberOfParticles();
  };

  inline double protonFrac() const {
    return nProtons_/numberOfParticles();
  };

  inline double hadFrac() const {
    return nHadrons_/numberOfParticles();
  };

  void Add(const G4SiHit & aSiHit);

  void encodeCellId(const bool x_side,const bool y_side,const unsigned x_cell,const unsigned y_cell);

  inline bool get_x_side() const{
    return cellid_ & 0x0001;
  };

  inline bool get_y_side() const {
    return (cellid_ & 0x00010000) >> 16;
  };

  inline unsigned get_x_cell() const {
    return (cellid_ & 0xFFFE) >> 1;
  };

  inline unsigned get_y_cell() const {
    return (cellid_ & 0xFFFE0000) >> 17;
  };

  inline double get_x() const {
    float sign = get_x_side() ? 1. : -1. ;
    if (sign > 0)
      return get_x_cell()*sign*CELL_SIZE_X*getGranularity()+CELL_SIZE_X*getGranularity()/2;
    else return get_x_cell()*sign*CELL_SIZE_X*getGranularity()-CELL_SIZE_X*getGranularity()/2;
  };

  inline double get_y() const {
    float sign = get_y_side() ? 1. : -1. ;
    if (sign > 0)
      return get_y_cell()*sign*CELL_SIZE_Y*getGranularity()+CELL_SIZE_Y*getGranularity()/2;
    else return get_y_cell()*sign*CELL_SIZE_Y*getGranularity()-CELL_SIZE_Y*getGranularity()/2;
  };

  inline bool get_x_side_old() const{
    return cellid_ & 0x0001;
  };

  inline bool get_y_side_old() const {
    return (cellid_ & 0x0100) >> 8;
  };

  inline unsigned get_x_cell_old() const {
    return (cellid_ & 0x00FE) >> 1;
  };

  inline unsigned get_y_cell_old() const {
    return (cellid_ & 0xFE00) >> 9;
  };

  inline double get_x_old() const {
    float sign = get_x_side_old() ? 1. : -1. ;
    if (sign > 0)
      return get_x_cell_old()*sign*CELL_SIZE_X*getGranularity()+CELL_SIZE_X*getGranularity()/2;
    else return get_x_cell_old()*sign*CELL_SIZE_X*getGranularity()-CELL_SIZE_X*getGranularity()/2;
  };

  inline double get_y_old() const {
    float sign = get_y_side_old() ? 1. : -1. ;
    if (sign > 0)
      return get_y_cell_old()*sign*CELL_SIZE_Y*getGranularity()+CELL_SIZE_Y*getGranularity()/2;
    else return get_y_cell_old()*sign*CELL_SIZE_Y*getGranularity()-CELL_SIZE_Y*getGranularity()/2;
  };

  inline double get_z() const {
    return zpos_;
  };

  inline unsigned getGranularity() const{
    return 1;
  };

  inline int mainParentTrackID() const{
    return trackIDMainParent_;
  }

  inline double mainParentEfrac() const {
    return energyMainParent_/energy_;
  }

  void PrintGeometry(std::ostream & aOs) const ;

  void Print(std::ostream & aOs) const ;

private:

  double energy_;
  double time_;
  double zpos_;
  unsigned layer_;
  unsigned cellid_;
  unsigned nGammas_;
  unsigned nElectrons_;
  unsigned nMuons_;
  unsigned nNeutrons_;
  unsigned nProtons_;
  unsigned nHadrons_;
  int trackIDMainParent_;
  double energyMainParent_;

  ClassDef(HGCSSSimHit,1);



};


typedef std::vector<HGCSSSimHit> HGCSSSimHitVec;



#endif
