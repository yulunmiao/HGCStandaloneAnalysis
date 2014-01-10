#ifndef _samplingsection_hh_
#define _samplingsection_hh_

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <iomanip>

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
    resetCounters();
  }

  //DTOR
  ~SamplingSection() { }
  
  //
  void add(G4double den, G4double dl, G4double globalTime,G4int pdgId,G4VPhysicalVolume* vol);
  
  //reset
  inline void resetCounters()
  {
    Pb_den=0;  Pb_dl=0;
    Cu_den=0;  Cu_dl=0;
    Si_den=0;  Si_dl=0; Si_time=0;
    Si_gFlux=0; Si_eFlux=0; Si_muFlux=0; Si_hadFlux=0; 
    PCB_den=0; PCB_dl=0;
    Air_den=0; Air_dl=0;
  }
  
  //
  G4double getMeasuredEnergy(bool weighted=true);
  G4double getAbsorbedEnergy();
  G4double getTotalEnergy();
  G4double getAbsorberX0();  
  G4double getX0transversed();
  G4double getHadronicFraction();
  G4double getMuonFraction();
  G4double getPhotonFraction();
  G4double getElectronFraction();
  G4double getAverageTime();

  //
  void report(bool header=false);

  //members
  G4double           Pb_thick, Cu_thick, Si_thick, PCB_thick, Air_thick;
  G4double           Pb_X0,    Cu_X0,    Si_X0,    PCB_X0,    Air_X0;
  G4double           Pb_den,   Cu_den,   Si_den,   PCB_den,   Air_den;
  G4double           Pb_dl,    Cu_dl,    Si_dl,    PCB_dl,    Air_dl;
  G4VPhysicalVolume* Pb_vol,  *Cu_vol,  *Si_vol,  *PCB_vol,  *Air_vol;
  G4double           Si_gFlux, Si_eFlux, Si_muFlux, Si_hadFlux, Si_time;
  G4double Total_thick;
};

#endif
