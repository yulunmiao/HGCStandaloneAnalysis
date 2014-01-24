#include "TransverseGeometry.hh"

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>

TransverseGeometry::TransverseGeometry(){
  cellid_ = 0;
  layer_=0;
}

void TransverseGeometry::SetHit(const G4SiHit & aSiHit){

  //coordinates in mm
  //double z = aSiHit.hit_x;
  double x = aSiHit.hit_y;
  double y = aSiHit.hit_z;

  //cellid encoding:
  bool x_side = x>0 ? true : false;
  bool y_side = y>0 ? true : false;

  layer_ = aSiHit.layer;
  
  unsigned x_cell = static_cast<unsigned>(fabs(x)/(CELL_SIZE_X*GRANULARITY[layer_]));
  unsigned y_cell = static_cast<unsigned>(fabs(y)/(CELL_SIZE_Y*GRANULARITY[layer_]));

  encodeCellId(x_side,y_side,x_cell,y_cell);

}

void TransverseGeometry::SetHit(const HGCSSSimHit & aHit){
  layer_ = aHit.layer();
  cellid_ = aHit.cellid();
}

void TransverseGeometry::PrintGeometry(std::ostream & aOs){
  aOs << " =========================================================== " << std::endl
      << " ============= Printing of Transverse Geometry ============= " << std::endl
      << " =========================================================== " << std::endl
      << " = SIZE_X " << SIZE_X << " mm" << std::endl
      << " = SIZE_Y " << SIZE_Y << " mm" << std::endl
      << " = N_LAYERS " << N_LAYERS << std::endl
      << " = CELL_SIZE_X " << CELL_SIZE_X << " mm" << std::endl
      << " = CELL_SIZE_Y " << CELL_SIZE_Y << " mm" << std::endl
      << " = N_CELLS_XY_MAX " << N_CELLS_XY_MAX << " per layer" << std::endl
    ;
  for (unsigned iL(0); iL<N_LAYERS; ++iL){
    aOs << " == LAYER " << iL << " GRANULARITY " << GRANULARITY[iL] << " N_CELLS " << SIZE_X/(CELL_SIZE_X*GRANULARITY[iL])*SIZE_Y/(CELL_SIZE_Y*GRANULARITY[iL]) << std::endl;
  }
  aOs << " =========================================================== " << std::endl
      << " =========================================================== " << std::endl;
}

void TransverseGeometry::encodeCellId(bool x_side,bool y_side,unsigned x_cell,unsigned y_cell){
  cellid_ = 
    x_side | (x_cell<<1) |
    (y_side<<8) | (y_cell<<9);

  // std::cout << " Cross-check of encoding: cellid=" << cellid_ << std::endl
  // 	    << " x_side " << x_side << " " << get_x_side() << std::endl
  // 	    << " y_side " << y_side << " " << get_y_side() << std::endl
  // 	    << " x_cell " << x_cell << " " << get_x_cell() << std::endl
  // 	    << " y_cell " << y_cell << " " << get_y_cell() << std::endl
  //   ;
}
