#include "EventAction.hh"

#include "RunAction.hh"
#include "EventActionMessenger.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//
EventAction::EventAction()
{
  runAct = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
  eventMessenger = new EventActionMessenger(this);
  printModulo = 5;
  outF_=TFile::Open("PFcal.root","RECREATE");
  ntuple_=new TNtuple("CaloStack","CaloStack","event:volNb:volX0:volX0trans:den:denWeight:denAbs:denTotal:gFrac:eFrac:muFrac:hadFrac:avgTime");
}

//
EventAction::~EventAction()
{
  outF_->cd();
  ntuple_->Write();
  outF_->Close();

  delete eventMessenger;
}

//
void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  evtNb_ = evt->GetEventID();
  if (evtNb_%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb_ << G4endl;
    CLHEP::HepRandom::showEngineStatus();
  }


}

//
void EventAction::Detect(G4double edep, G4double stepl,G4double globalTime, G4int pdgId, G4VPhysicalVolume *volume)
{
  for(size_t i=0; i<detector_->size(); i++) (*detector_)[i].add(edep,stepl,globalTime,pdgId,volume);
}

//
void EventAction::EndOfEventAction(const G4Event* evt)
{
  bool debug(evtNb_%printModulo == 0);

  for(size_t i=0; i<detector_->size(); i++) 
    {
      event_[0]=evtNb_;
      event_[1]=i;
      event_[2]=(*detector_)[i].getAbsorberX0();
      event_[3]=(*detector_)[i].getX0transversed();
      event_[4]=(*detector_)[i].getMeasuredEnergy(false);
      event_[5]=(*detector_)[i].getMeasuredEnergy(true);
      event_[6]=(*detector_)[i].getAbsorbedEnergy();
      event_[7]=(*detector_)[i].getTotalEnergy();
      event_[8]=(*detector_)[i].getPhotonFraction();
      event_[9]=(*detector_)[i].getElectronFraction();
      event_[10]=(*detector_)[i].getMuonFraction();
      event_[11]=(*detector_)[i].getHadronicFraction();
      event_[12]=(*detector_)[i].getAverageTime();
      ntuple_->Fill(event_);
      if(debug) (*detector_)[i].report( (i==0) );
      (*detector_)[i].resetCounters();
    }
}  
