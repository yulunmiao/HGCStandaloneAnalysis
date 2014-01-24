#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//
SteppingAction::SteppingAction()                                         
{
  eventAction_ = (EventAction*)G4RunManager::GetRunManager()->GetUserEventAction();               
  eventAction_->Add(  ((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getStructure() );
}

//
SteppingAction::~SteppingAction()
{ }

//
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) stepl = aStep->GetStepLength();
  const G4ThreeVector &pos=aStep->GetPreStepPoint()->GetPosition();

  Int_t iz=Int_t(pos.z()/eventAction_->GetCellSize())+4;
  if(iz<0) iz=0; if(iz>8) iz=8;
  Int_t iy=Int_t(pos.y()/eventAction_->GetCellSize())+4;
  if(iy<0) iy=0; if(iy>8) iy=8;
  Int_t iyiz(iz+iy*9);

  G4int pdgId=aStep->GetTrack()->GetDefinition()->GetPDGEncoding(); 
  G4double globalTime=aStep->GetTrack()->GetGlobalTime();
  eventAction_->Detect(edep,stepl,globalTime,pdgId,volume,iyiz);
}
