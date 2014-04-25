#ifndef _hgcssgenparticle_hh_
#define _hgcssgenparticle_hh_

#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <map>

class HGCSSGenParticle{

public:
  HGCSSGenParticle():
    time_(0),
    xpos_(0),
    ypos_(0),
    zpos_(0),
    mass_(0),
    px_(0),
    py_(0),
    pz_(0),
    pdgid_(0),
    charge_(0),
    trackID_(0)
  {
    
  };

  ~HGCSSGenParticle(){};

  inline double mass() const {
    return mass_;
  };
  inline double time() const {
    return time_;
  };
  inline double x() const {
    return xpos_;
  };
  inline double y() const {
    return ypos_;
  };
  inline double z() const {
    return zpos_;
  };
  inline double px() const {
    return px_;
  };
  inline double py() const {
    return py_;
  };
  inline double pz() const {
    return pz_;
  };
  inline int pdgid() const {
    return pdgid_;
  };
  inline double charge() const {
    return charge_;
  };
  inline int trackID() const {
    return trackID_;
  };

  inline void setPosition(const double & x, const double & y, const double & z) {
    xpos_ = x;
    ypos_ = y;
    zpos_ = z;
  };

  inline void setMomentum(const double & px, const double & py, const double & pz) {
    px_ = px;
    py_ = py;
    pz_ = pz;
  };


  inline void mass(const double & val) {
    mass_ = val;
  };
  inline void time(const double & val) {
    time_ = val;
  };
  inline void pdgid(const int & val) {
    pdgid_ = val;
  };
  inline void charge(const double & val) {
    charge_ = val;
  };
  inline void trackID(const int & val) {
    trackID_ = val;
  };

  inline bool isIncoming() const{
    return ((pdgid_ != 0 || mass_ > 0) && pz_ >= 0 && time_ < 10);
  };

  void Print(std::ostream & aOs) const ;

private:

  double time_;
  double xpos_;
  double ypos_;
  double zpos_;
  double mass_;
  double px_;
  double py_;
  double pz_;
  int pdgid_;
  double charge_;
  int trackID_;

  ClassDef(HGCSSGenParticle,1);



};


typedef std::vector<HGCSSGenParticle> HGCSSGenParticleVec;



#endif
