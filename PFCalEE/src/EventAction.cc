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
  //ntuple_=new TNtuple("CaloStack","CaloStack","event:volNb:volX0:volX0trans:den:denWeight:denAbs:denTotal:gFrac:eFrac:muFrac:hadFrac:avgTime:nhits");
  tree_=new TTree("HGCSSTree","HGC Standalone simulation tree");
  tree_->Branch("event",&event_[0]);
  tree_->Branch("volNb",&event_[1]);
  tree_->Branch("volX0",&event_[2]);
  tree_->Branch("volX0trans",&event_[3]);
  tree_->Branch("den",&event_[4]);
  tree_->Branch("denWeight",&event_[5]);
  tree_->Branch("denAbs",&event_[6]);
  tree_->Branch("denTotal",&event_[7]);
  tree_->Branch("gFrac",&event_[8]);
  tree_->Branch("eFrac",&event_[9]);
  tree_->Branch("muFrac",&event_[10]);
  tree_->Branch("hadFrac",&event_[11]);
  tree_->Branch("avgTime",&event_[12]);
  tree_->Branch("nSiHits",&event_[13]);
  tree_->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&hitvec_);

  hitGeom_.PrintGeometry(G4cout);

}

//
EventAction::~EventAction()
{
  outF_->cd();
  //ntuple_->Write();
  tree_->Write();
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
void EventAction::Detect(G4double edep, G4double stepl,G4double globalTime, G4int pdgId, G4VPhysicalVolume *volume, G4ThreeVector position)
{
  for(size_t i=0; i<detector_->size(); i++) (*detector_)[i].add(edep,stepl,globalTime,pdgId,volume,position,i);
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
      event_[13]=(*detector_)[i].getSiHitVec().size();
      hitvec_.clear();
      std::map<unsigned,HGCSSSimHit> lHitMap;
      std::pair<std::map<unsigned,HGCSSSimHit>::iterator,bool> isInserted;

      for (unsigned iSiHit(0); iSiHit<event_[13];++iSiHit){
	G4SiHit lSiHit = (*detector_)[i].getSiHitVec()[iSiHit];
	hitGeom_.SetHit(lSiHit);
	HGCSSSimHit lHit(lSiHit,hitGeom_.cellid());
	isInserted = lHitMap.insert(std::pair<unsigned,HGCSSSimHit>(hitGeom_.cellid(),lHit));
	if (!isInserted.second) isInserted.first->second.Add(lSiHit);
      }
      std::map<unsigned,HGCSSSimHit>::iterator lIter = lHitMap.begin();
      hitvec_.reserve(lHitMap.size());
      for (; lIter != lHitMap.end(); ++lIter){
	hitvec_.push_back(lIter->second);
      }
      //ntuple_->Fill(event_);
      tree_->Fill();
      if(debug) (*detector_)[i].report( (i==0) );
      (*detector_)[i].resetCounters();
    }
}
