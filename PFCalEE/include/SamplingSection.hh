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
  SamplingSection(G4double Pb, G4double Cu,
		  G4double Si0, G4double Si1, G4double Si2,
		  G4double PCB, G4double Air)
  {
    Pb_thick=Pb;     Pb_X0=0; Pb_vol=0;
    Cu_thick=Cu;     Cu_X0=0; Cu_vol=0;
    Si_thick[0]=Si0; Si_thick[1]=Si1; Si_thick[2]=Si2;
    Si_X0=0; 
    Si_vol[0]=0; Si_vol[1]=0; Si_vol[2]=0;
    PCB_thick=PCB;   PCB_X0=0; PCB_vol=0;
    Air_thick=Air;   Air_X0=0; Air_vol=0;
    Total_thick=Pb+Cu+Si0+Si1+Si2+PCB+Air;
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
    for (unsigned idx(0); idx<3; ++idx){
      Si_den[idx]=0;  Si_dl[idx]=0; Si_time[idx]=0;
      Si_gFlux[idx]=0; Si_eFlux[idx]=0; Si_muFlux[idx]=0; Si_hadFlux[idx]=0; 
      //reserve some space based on first event....
      if (Si_HitVec[idx].size() > Si_HitVec_size_max) {
	Si_HitVec_size_max = 2*Si_HitVec[idx].size();
	G4cout << "-- SamplingSection::resetCounters(), space reserved for HitVec vector increased to " << Si_HitVec_size_max << G4endl;
      }
      Si_HitVec[idx].clear();
      Si_HitVec[idx].reserve(Si_HitVec_size_max);
    }
    PCB_den=0; PCB_dl=0;
    Air_den=0; Air_dl=0;
  }
  
  //
  G4double getMeasuredEnergy(bool weighted=true);
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

  const G4SiHitVec & getSiHitVec(const unsigned & idx) const;
  void trackParticleHistory(const unsigned & idx,const G4SiHitVec & incoming);

  //
  void report(bool header=false);

  //members
  G4double           Pb_thick, Cu_thick, Si_thick[3], PCB_thick, Air_thick;
  G4double           Pb_X0,    Cu_X0,    Si_X0,    PCB_X0,    Air_X0;
  G4double           Pb_L0,    Cu_L0;
  G4double           Pb_den,   Cu_den,   Si_den[3],   PCB_den,   Air_den;
  G4double           Pb_dl,    Cu_dl,    Si_dl[3],    PCB_dl,    Air_dl;
  G4VPhysicalVolume* Pb_vol,  *Cu_vol,  *Si_vol[3], *PCB_vol,  *Air_vol;
  G4double           Si_gFlux[3], Si_eFlux[3], Si_muFlux[3], Si_neutronFlux[3], Si_hadFlux[3], Si_time[3];
  G4double Total_thick;
  G4SiHitVec Si_HitVec[3];
  unsigned Si_HitVec_size_max;


};

#endif
