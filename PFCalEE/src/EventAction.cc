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
  printModulo = 10;
  outF_=TFile::Open("PFcal.root","RECREATE");
  //ntuple_=new TNtuple("CaloStack","CaloStack","event:volNb:volX0:volX0trans:den:denWeight:denAbs:denTotal:gFrac:eFrac:muFrac:hadFrac:avgTime:nhits");
  // ntuple_->Branch("dendydz",dendydz_,"dendydz[81]/F");
  //cellSize_=4;
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
  tree_->Branch("volLambda",&event_[14]);
  tree_->Branch("neutronFrac",&event_[15]);
  tree_->Branch("nGenParticles",&event_[16]);
  tree_->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&hitvec_);
  tree_->Branch("HGCSSGenParticleVec","std::vector<HGCSSGenParticle>",&genvec_);
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
void EventAction::Detect(G4double edep, G4double stepl,G4double globalTime, 
			 G4int pdgId, G4VPhysicalVolume *volume, const G4ThreeVector & position, 
			 G4int trackID, G4int parentID,
			 const HGCSSGenParticle & genPart)
{
  for(size_t i=0; i<detector_->size(); i++) (*detector_)[i].add(edep,stepl,globalTime,pdgId,volume,position,trackID,parentID,i);
  if (genPart.isIncoming()) genvec_.push_back(genPart);
}

  //void EventAction::Detect(G4double edep, G4double stepl,G4double globalTime, G4int pdgId, G4VPhysicalVolume *volume, int iyiz)
  //{
  //for(size_t i=0; i<detector_->size(); i++) (*detector_)[i].add(edep,stepl,globalTime,pdgId,volume,iyiz);
  //}

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
      event_[14]=(*detector_)[i].getAbsorberLambda();
      event_[15]=(*detector_)[i].getNeutronFraction();
      event_[16]=genvec_.size();
      hitvec_.clear();
      std::map<unsigned,HGCSSSimHit> lHitMap;
      std::pair<std::map<unsigned,HGCSSSimHit>::iterator,bool> isInserted;

      if (i>0) (*detector_)[i].trackParticleHistory((*detector_)[i-1].getSiHitVec());

      for (unsigned iSiHit(0); iSiHit<event_[13];++iSiHit){
	G4SiHit lSiHit = (*detector_)[i].getSiHitVec()[iSiHit];
	HGCSSSimHit lHit(lSiHit);

	//print geometry just once for the record
	if (evtNb_ == 0 && i==0 && iSiHit==0){
	  lHit.PrintGeometry(G4cout);
	}

	isInserted = lHitMap.insert(std::pair<unsigned,HGCSSSimHit>(lHit.cellid(),lHit));
	if (!isInserted.second) isInserted.first->second.Add(lSiHit);
      }
      std::map<unsigned,HGCSSSimHit>::iterator lIter = lHitMap.begin();
      hitvec_.reserve(lHitMap.size());
      for (; lIter != lHitMap.end(); ++lIter){
	(lIter->second).calculateTime();
	hitvec_.push_back(lIter->second);
      }
      //ntuple_->Fill(event_);
      tree_->Fill();

      //for(size_t iyiz=0; iyiz<81; iyiz++) dendydz_[iyiz]=(*detector_)[i].getMeasuredEnergyInPos(iyiz);
      //ntuple_->Fill(event_);

      if(debug) {
	if (i==0) G4cout << " -- Number of truth particles = " << genvec_.size() << G4endl;
	(*detector_)[i].report( (i==0) );
      }

      //if (i==0) G4cout << " ** evt " << evt->GetEventID() << G4endl;

      //reset genvec: want to record only for first layer
      genvec_.clear();

      (*detector_)[i].resetCounters();
    }
}
