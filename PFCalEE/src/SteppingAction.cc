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
  // get PreStepPoint
  const G4StepPoint *thePreStepPoint = aStep->GetPreStepPoint();
  // get TouchableHandle
  G4TouchableHandle theTouchable = thePreStepPoint->GetTouchableHandle();


  G4VPhysicalVolume* volume = theTouchable->GetVolume();
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) stepl = aStep->GetStepLength();
  G4int pdgId=aStep->GetTrack()->GetDefinition()->GetPDGEncoding(); 
  G4double globalTime=aStep->GetTrack()->GetGlobalTime();

  G4ThreeVector position = thePreStepPoint->GetPosition();
  // get local hit position using touchable with theGlobalPos
  //G4ThreeVector LocalHitPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(position);
 
  //G4cout << "Pre   position " << volume->GetName() << " " << position << G4endl;
  //G4cout << "Post  position " << aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << " " << aStep->GetPostStepPoint()->GetPosition() << G4endl;
  //G4cout << "Local position " << LocalHitPos << G4endl;
  
  eventAction_->Detect(edep,stepl,globalTime,pdgId,volume,position);
}
