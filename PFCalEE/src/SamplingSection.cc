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
  if(Pb_vol && vol->GetName()==Pb_vol->GetName())        { Pb_den+=den;  Pb_dl+=dl;  }
  else if(Cu_vol && vol->GetName()==Cu_vol->GetName())   { Cu_den+=den;  Cu_dl+=dl;  }
  else if(Si_vol && vol->GetName()==Si_vol->GetName())   
    { 
      Si_den+=den;  Si_dl+=dl;  Si_time+=den*globalTime;

      //Si_dendydz[iyiz] += den;

      //discriminate further by particle type
      if(abs(pdgId)==22)      Si_gFlux += den;
      else if(abs(pdgId)==11) Si_eFlux += den;
      else if(abs(pdgId)==13) Si_muFlux += den;
      else if (abs(pdgId)==2112) Si_neutronFlux += den;
      else {
	Si_hadFlux += den;
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
      Si_HitVec.push_back(lHit);
     
    }
  else if(PCB_vol && vol->GetName()==PCB_vol->GetName()) { PCB_den+=den; PCB_dl+=dl; }
  else if(Air_vol && vol->GetName()==Air_vol->GetName()) { Air_den+=den; Air_dl+=dl; }

}

//
void SamplingSection::report(bool header)
{
  if(header) G4cout << "E/[MeV]\t  Si\tAbsorber\tTotal\tSi g frac\tSi e frac\tSi mu frac\tSi had frac\tSi <t> \t nG4SiHits" << G4endl;
  G4cout << std::setprecision(3) << "\t  " << getMeasuredEnergy(false) << "\t" << getAbsorbedEnergy() << "\t\t" << getTotalEnergy() << "\t"
	 << getPhotonFraction() << "\t" << getElectronFraction() << "\t" << getMuonFraction() << "\t" << getHadronicFraction() << "\t"
	 << getAverageTime() << "\t"
	 << Si_HitVec.size() << "\t"
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
G4double SamplingSection::getNeutronFraction()
{
  return Si_den > 0 ? Si_neutronFlux/Si_den : 0;
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
  return Air_den+getMeasuredEnergy(false)+getAbsorbedEnergy();
}

const G4SiHitVec & SamplingSection::getSiHitVec() const
{
  return Si_HitVec;
}

void SamplingSection::trackParticleHistory(const G4SiHitVec & incoming)
{
  for (unsigned iP(0); iP<Si_HitVec.size(); ++iP){//loop on g4hits
    G4int parId = Si_HitVec[iP].parentId;
    for (unsigned iI(0); iI<incoming.size(); ++iI){//loop on previous layer
      G4int trId = incoming[iI].trackId;
      if (trId == parId) Si_HitVec[iP].parentId = incoming[iI].parentId;

    }//loop on previous layer


  }//loop on g4hits


  //  return Si_HitVec;
}
