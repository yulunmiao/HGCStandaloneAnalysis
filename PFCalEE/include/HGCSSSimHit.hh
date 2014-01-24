#ifndef _hgcsssimhit_hh_
#define _hgcsssimhit_hh_

#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>

#include "G4SiHit.hh"

class HGCSSSimHit{

public:
  HGCSSSimHit():
    energy_(0),
    time_(0),
    layer_(0),
    cellid_(0),
    nGammas_(0),
    nElectrons_(0),
    nMuons_(0),
    nHadrons_(0)
  {};
  HGCSSSimHit(const G4SiHit & aSiHit,
	      const unsigned & cellid);

  ~HGCSSSimHit(){};

  inline double energy() const {
    return energy_;
  };
  inline double time() const {
    return time_;
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
  inline unsigned nHadrons() const {
    return nHadrons_;
  };

  inline unsigned numberOfParticles() const {
    return nGammas_+nElectrons_+nMuons_+nHadrons_;
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

  inline double hadFrac() const {
    return nHadrons_/numberOfParticles();
  };

  void Add(const G4SiHit & aSiHit);
  
  void Print(std::ostream & aOs);

private:

  double energy_;
  double time_;
  unsigned layer_;
  unsigned cellid_;
  unsigned nGammas_;
  unsigned nElectrons_;
  unsigned nMuons_;
  unsigned nHadrons_;

  ClassDef(HGCSSSimHit,1);



};


typedef std::vector<HGCSSSimHit> HGCSSSimHitVec;



#endif
