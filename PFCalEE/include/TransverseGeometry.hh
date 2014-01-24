#ifndef _transversegeometry_hh_
#define _transversegeometry_hh_

#include <iomanip>
#include <vector>
#include <map>

#include "G4SiHit.hh"
#include "HGCSSSimHit.hh"

static const float SIZE_X=200;//mm
static const float SIZE_Y=SIZE_X;
static const float CELL_SIZE_X=2.5;//mm
static const float CELL_SIZE_Y=CELL_SIZE_X;
static const unsigned N_CELLS_XY_MAX=SIZE_X/CELL_SIZE_X*SIZE_Y/CELL_SIZE_Y;
static const unsigned N_LAYERS=30;
static const unsigned GRANULARITY[N_LAYERS]={
  1,1,1,1,1,
  1,1,1,1,1,
  1,1,1,1,1,
  1,1,1,1,1,
  1,1,1,1,1,
  1,1,1,1,1
};

class TransverseGeometry {

public:
  TransverseGeometry();
  ~TransverseGeometry(){};

  void SetHit(const G4SiHit & aSiHit);
  void SetHit(const HGCSSSimHit & aHit);

  inline unsigned cellid(){
    return cellid_;
  };
  
  inline unsigned layer(){
    return layer_;
  };

  inline double weight(){
    if (layer_ < 10) return 1;
    else if (layer_ < 20) return 2;
    else return 3;
  }

  void PrintGeometry(std::ostream & aOs);

  void encodeCellId(bool x_side,bool y_side,unsigned x_cell,unsigned y_cell);

  inline bool get_x_side(){
    return cellid_ & 0x0001;
  };

  inline bool get_y_side(){
    return (cellid_ & 0x0100) >> 8;
  };

  inline unsigned get_x_cell(){
    return (cellid_ & 0x00FE) >> 1;
  };

  inline unsigned get_y_cell(){
    return (cellid_ & 0xFE00) >> 9;
  };

  inline double get_x(){
    float sign = get_x_side() ? 1. : -1. ;
    if (sign > 0)
      return get_x_cell()*sign*CELL_SIZE_X*GRANULARITY[layer_]+CELL_SIZE_X*GRANULARITY[layer_]/2;
    else return get_x_cell()*sign*CELL_SIZE_X*GRANULARITY[layer_]-CELL_SIZE_X*GRANULARITY[layer_]/2;
  };

  inline double get_y(){
    float sign = get_y_side() ? 1. : -1. ;
    if (sign > 0)
      return get_y_cell()*sign*CELL_SIZE_Y*GRANULARITY[layer_]+CELL_SIZE_Y*GRANULARITY[layer_]/2;
    else return get_y_cell()*sign*CELL_SIZE_Y*GRANULARITY[layer_]-CELL_SIZE_Y*GRANULARITY[layer_]/2;
  };


private:
  unsigned cellid_;
  unsigned layer_;


};

#endif
