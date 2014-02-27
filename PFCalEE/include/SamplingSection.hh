#ifndef _samplingsection_hh_
#define _samplingsection_hh_

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <iomanip>
#include <vector>

#include "G4SiHit.hh"

class SamplingSection
{
public:
  //CTOR
  SamplingSection(G4double Pb, G4double Cu, G4double Si, G4double PCB, G4double Air)
  {
    Pb_thick=Pb;     Pb_X0=0; Pb_vol=0;
    Cu_thick=Cu;     Cu_X0=0; Cu_vol=0;
    Si_thick=Si;     Si_X0=0; Si_vol=0;
    PCB_thick=PCB;   PCB_X0=0; PCB_vol=0;
    Air_thick=Air;   Air_X0=0; Air_vol=0;
    Total_thick=Pb+Cu+Si+PCB+Air;
    Si_HitVec_size_max = 0;
    resetCounters();
  }

  //DTOR
  ~SamplingSection() { }
  
  //
  void add(G4double den, G4double dl, G4double globalTime,G4int pdgId,G4VPhysicalVolume* vol, 
	   const G4ThreeVector & position,
	   G4int trackID, G4int parentID,
	   G4int layerId);
  //void add(G4double den, G4double dl, G4double globalTime,G4int pdgId,G4VPhysicalVolume* vol, int iyiz);
  
  //reset
  inline void resetCounters()
  {
    Pb_den=0;  Pb_dl=0;
    Cu_den=0;  Cu_dl=0;
    Si_den=0;  Si_dl=0; Si_time=0;
    for(size_t i=0; i<81; i++) Si_dendydz[i]=0;
    Si_gFlux=0; Si_eFlux=0; Si_muFlux=0; Si_hadFlux=0; 
    PCB_den=0; PCB_dl=0;
    Air_den=0; Air_dl=0;
    //reserve some space based on first event....
    if (Si_HitVec.size() > Si_HitVec_size_max) {
      Si_HitVec_size_max = 2*Si_HitVec.size();
      G4cout << "-- SamplingSection::resetCounters(), space reserved for HitVec vector increased to " << Si_HitVec_size_max << G4endl;
    }
    Si_HitVec.clear();
    Si_HitVec.reserve(Si_HitVec_size_max);
  }
  
  //
  G4double getMeasuredEnergy(bool weighted=true);
  G4double getMeasuredEnergyInPos(int iyiz) { return Si_dendydz[iyiz]; }
  G4double getAbsorbedEnergy();
  G4double getTotalEnergy();
  G4double getAbsorberX0();  
  G4double getX0transversed();
  G4double getAbsorberLambda();
  G4double getHadronicFraction();
  G4double getNeutronFraction();
  G4double getMuonFraction();
  G4double getPhotonFraction();
  G4double getElectronFraction();
  G4double getAverageTime();

  const G4SiHitVec & getSiHitVec() const;
  void trackParticleHistory(const G4SiHitVec & incoming);

  //
  void report(bool header=false);

  //members
  G4double           Pb_thick, Cu_thick, Si_thick, PCB_thick, Air_thick;
  G4double           Pb_X0,    Cu_X0,    Si_X0,    PCB_X0,    Air_X0;
  G4double           Pb_L0,    Cu_L0;
  G4double           Pb_den,   Cu_den,   Si_den,   PCB_den,   Air_den, Si_dendydz[81];
  G4double           Pb_dl,    Cu_dl,    Si_dl,    PCB_dl,    Air_dl;
  G4VPhysicalVolume* Pb_vol,  *Cu_vol,  *Si_vol,  *PCB_vol,  *Air_vol;
  G4double           Si_gFlux, Si_eFlux, Si_muFlux, Si_neutronFlux, Si_hadFlux, Si_time;
  G4double Total_thick;
  G4SiHitVec Si_HitVec;
  unsigned Si_HitVec_size_max;


};

#endif
