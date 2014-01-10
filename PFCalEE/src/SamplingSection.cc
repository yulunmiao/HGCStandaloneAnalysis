#include "G4VPhysicalVolume.hh"

#include "SamplingSection.hh"

//
void SamplingSection::add(G4double den, G4double dl, G4double globalTime, G4int pdgId, G4VPhysicalVolume* vol)
{
  if(Pb_vol && vol->GetName()==Pb_vol->GetName())        { Pb_den+=den;  Pb_dl+=dl;  }
  else if(Cu_vol && vol->GetName()==Cu_vol->GetName())   { Cu_den+=den;  Cu_dl+=dl;  }
  else if(Si_vol && vol->GetName()==Si_vol->GetName())   
    { 
      Si_den+=den;  Si_dl+=dl;  Si_time+=den*globalTime;
      
      //discriminate further by particle type
      if(abs(pdgId)==22)      Si_gFlux += den;
      else if(abs(pdgId)==11) Si_eFlux += den;
      else if(abs(pdgId)==13) Si_muFlux += den;
      else                    Si_hadFlux += den;
    }
  else if(PCB_vol && vol->GetName()==PCB_vol->GetName()) { PCB_den+=den; PCB_dl+=dl; }
  else if(Air_vol && vol->GetName()==Air_vol->GetName()) { Air_den+=den; Air_dl+=dl; }

}

//
void SamplingSection::report(bool header)
{
  if(header) G4cout << "E/[MeV]\t  Si\tAbsorber\tTotal\tSi g frac\tSi e frac\tSi mu frac\tSi had frac\tSi <t>" << G4endl;
  G4cout << std::setprecision(3) << "\t  " << getMeasuredEnergy(false) << "\t" << getAbsorbedEnergy() << "\t\t" << getTotalEnergy() << "\t"
	 << getPhotonFraction() << "\t" << getElectronFraction() << "\t" << getMuonFraction() << "\t" << getHadronicFraction() << "\t"
	 << getAverageTime()
	 << G4endl; 
}

G4double SamplingSection::getAverageTime()
{
  return Si_den>0 ? Si_time/Si_den : 0 ;
}

//
G4double SamplingSection::getPhotonFraction()
{
  return Si_den > 0 ? Si_gFlux/Si_den : 0;
}

//
G4double SamplingSection::getElectronFraction()
{
  return Si_den > 0 ? Si_eFlux/Si_den : 0;
}

//
G4double SamplingSection::getMuonFraction()
{
  return Si_den > 0 ? Si_muFlux/Si_den : 0;
}

//
G4double SamplingSection::getHadronicFraction()
{
  return Si_den > 0 ? Si_hadFlux/Si_den : 0;
}

//
G4double SamplingSection::getMeasuredEnergy(bool weighted)
{
  G4double weight=(weighted ? getAbsorberX0() : 1.0);
  return weight*Si_den;
}

//
G4double SamplingSection::getAbsorberX0()
{
  // G4cout << Pb_thick << " " << Pb_X0 << G4endl;
  return Pb_thick*Pb_X0+Cu_thick*Cu_X0;
}

//
G4double SamplingSection::getX0transversed()
{
  G4double x0Trans(0);
  if(Pb_X0>0) x0Trans += Pb_thick/Pb_X0;
  if(Cu_X0>0) x0Trans += Cu_thick/Cu_X0;
  return x0Trans;
}

//
G4double SamplingSection::getAbsorbedEnergy()
{
  return Pb_den+Cu_den;
}

//
G4double SamplingSection::getTotalEnergy()
{
  return Air_den+getMeasuredEnergy(false)+getAbsorbedEnergy();
}
