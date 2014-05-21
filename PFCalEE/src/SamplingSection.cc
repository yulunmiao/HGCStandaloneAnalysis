#include "G4VPhysicalVolume.hh"

#include "SamplingSection.hh"

//
void SamplingSection::add(G4double den, G4double dl, 
			  G4double globalTime, G4int pdgId, 
			  G4VPhysicalVolume* vol, 
			  const G4ThreeVector & position,
			  G4int trackID, G4int parentID,
			  G4int layerId)
//void SamplingSection::add(G4double den, G4double dl, G4double globalTime, G4int pdgId, G4VPhysicalVolume* vol, int iyiz)
{
  std::string lstr = vol->GetName();
  if(Pb_vol && lstr==Pb_vol->GetName())        { Pb_den+=den;  Pb_dl+=dl;  }
  else if(Cu_vol && lstr==Cu_vol->GetName())   { Cu_den+=den;  Cu_dl+=dl;  }
  else if(PCB_vol && lstr==PCB_vol->GetName()) { PCB_den+=den; PCB_dl+=dl; }
  else if(Air_vol && lstr==Air_vol->GetName()) { Air_den+=den; Air_dl+=dl; }
  else {
    unsigned idx = 3;
    if(Si_vol[0] && lstr==Si_vol[0]->GetName()) idx = 0;
    else if (Si_vol[1] && lstr==Si_vol[1]->GetName()) idx = 1;
    else if (Si_vol[2] && lstr==Si_vol[2]->GetName()) idx = 2;
    if (idx != 3){
      Si_den[idx]+=den;  Si_dl[idx]+=dl;  Si_time[idx]+=den*globalTime;
	
      //discriminate further by particle type
      if(abs(pdgId)==22)      Si_gFlux[idx] += den;
      else if(abs(pdgId)==11) Si_eFlux[idx] += den;
      else if(abs(pdgId)==13) Si_muFlux[idx] += den;
      else if (abs(pdgId)==2112) Si_neutronFlux[idx] += den;
      else {
	Si_hadFlux[idx] += den;
      }
	
      //add hit
      G4SiHit lHit;
      lHit.energy = den;
      lHit.time = globalTime;
      lHit.pdgId = pdgId;
      lHit.layer = layerId;
      lHit.hit_x = position.x();
      lHit.hit_y = position.y();
      lHit.hit_z = position.z();
      lHit.trackId = trackID;
      lHit.parentId = parentID;
      Si_HitVec[idx].push_back(lHit);
      
    }
  }

}

//
void SamplingSection::report(bool header)
{
  if(header) G4cout << "E/[MeV]\t  Si\tAbsorber\tTotal\tSi g frac\tSi e frac\tSi mu frac\tSi had frac\tSi <t> \t nG4SiHits" << G4endl;
  G4cout << std::setprecision(3) << "\t  " << getMeasuredEnergy(false) << "\t" << getAbsorbedEnergy() << "\t\t" << getTotalEnergy() << "\t"
	 << getPhotonFraction() << "\t" << getElectronFraction() << "\t" << getMuonFraction() << "\t" << getHadronicFraction() << "\t"
	 << getAverageTime() << "\t"
	 << Si_HitVec[0].size()+Si_HitVec[1].size()+Si_HitVec[2].size() << "\t"
	 << G4endl; 
}

G4double SamplingSection::getAverageTime()
{
  double etot = Si_den[0]+Si_den[1]+Si_den[2];
  return etot>0 ? (Si_time[0]+Si_time[1]+Si_time[2])/etot : 0 ;
}

//
G4double SamplingSection::getPhotonFraction()
{
  double etot = Si_den[0]+Si_den[1]+Si_den[2];
  return etot>0 ? (Si_gFlux[0]+Si_gFlux[1]+Si_gFlux[2])/etot : 0 ;
}

//
G4double SamplingSection::getElectronFraction()
{
  double etot = Si_den[0]+Si_den[1]+Si_den[2];
  return etot>0 ? (Si_eFlux[0]+Si_eFlux[1]+Si_eFlux[2])/etot : 0 ;
}

//
G4double SamplingSection::getMuonFraction()
{
  double etot = Si_den[0]+Si_den[1]+Si_den[2];
  return etot>0 ? (Si_muFlux[0]+Si_muFlux[1]+Si_muFlux[2])/etot : 0 ;
}

//
G4double SamplingSection::getNeutronFraction()
{
  double etot = Si_den[0]+Si_den[1]+Si_den[2];
  return etot>0 ? (Si_neutronFlux[0]+Si_neutronFlux[1]+Si_neutronFlux[2])/etot : 0 ;
}

//
G4double SamplingSection::getHadronicFraction()
{
  double etot = Si_den[0]+Si_den[1]+Si_den[2];
  return etot>0 ? (Si_hadFlux[0]+Si_hadFlux[1]+Si_hadFlux[2])/etot : 0 ;
}

//
G4double SamplingSection::getMeasuredEnergy(bool weighted)
{
  G4double weight=(weighted ? getAbsorberX0() : 1.0);
  double etot = Si_den[0]+Si_den[1]+Si_den[2];
  return weight*etot;
}

//
G4double SamplingSection::getAbsorberX0()
{
  // G4cout << Pb_thick << " " << Pb_X0 << G4endl;
  return Pb_thick*Pb_X0+Cu_thick*Cu_X0;
}

//
G4double SamplingSection::getAbsorberLambda()
{
  G4double lTrans(0);
  if(Pb_L0>0) lTrans += Pb_thick/Pb_L0;
  if(Cu_L0>0) lTrans += Cu_thick/Cu_L0;
  return lTrans;
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
  return Air_den+PCB_den+getMeasuredEnergy(false)+getAbsorbedEnergy();
}

const G4SiHitVec & SamplingSection::getSiHitVec(const unsigned & idx) const
{
  return Si_HitVec[idx];
}

void SamplingSection::trackParticleHistory(const unsigned & idx, const G4SiHitVec & incoming)
{
  for (unsigned iP(0); iP<Si_HitVec[idx].size(); ++iP){//loop on g4hits
    G4int parId = Si_HitVec[idx][iP].parentId;
    for (unsigned iI(0); iI<incoming.size(); ++iI){//loop on previous layer
      G4int trId = incoming[iI].trackId;
      if (trId == parId) Si_HitVec[idx][iP].parentId = incoming[iI].parentId;
    }//loop on previous layer
  }//loop on g4hits
  //  return Si_HitVec;
}
