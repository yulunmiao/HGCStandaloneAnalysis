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
  tree_->Branch("HGCSSEvent","HGCSSEvent",&event_);
  tree_->Branch("HGCSSSamplingSectionVec","std::vector<HGCSSSamplingSection>",&ssvec_);
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
  hitvec_.clear();

  event_.eventNumber(evtNb_);
  event_.cellSize(CELL_SIZE_X);

  ssvec_.clear();
  ssvec_.reserve(detector_->size());

  for(size_t i=0; i<detector_->size(); i++) 
    {
      HGCSSSamplingSection lSec;
      lSec.volNb(i);
      lSec.volX0((*detector_)[i].getAbsorberX0());
      lSec.volX0trans((*detector_)[i].getX0transversed());
      lSec.volLambdatrans((*detector_)[i].getAbsorberLambda());
      lSec.absorberE((*detector_)[i].getAbsorbedEnergy());

      lSec.measuredE((*detector_)[i].getMeasuredEnergy(false));
      lSec.totalE((*detector_)[i].getTotalEnergy());
      lSec.gFrac((*detector_)[i].getPhotonFraction());
      lSec.eFrac((*detector_)[i].getElectronFraction());
      lSec.muFrac((*detector_)[i].getMuonFraction());
      lSec.neutronFrac((*detector_)[i].getNeutronFraction());
      lSec.hadFrac((*detector_)[i].getHadronicFraction());
      lSec.avgTime((*detector_)[i].getAverageTime());
      lSec.nSiHits((*detector_)[i].getSiHitVec(0).size()+
		   (*detector_)[i].getSiHitVec(1).size()+
		   (*detector_)[i].getSiHitVec(2).size());
      ssvec_.push_back(lSec);

      for (unsigned idx(0); idx<3; ++idx){
	std::map<unsigned,HGCSSSimHit> lHitMap;
	std::pair<std::map<unsigned,HGCSSSimHit>::iterator,bool> isInserted;

	//if (i>0) (*detector_)[i].trackParticleHistory(idx,(*detector_)[i-1].getSiHitVec(idx));

	for (unsigned iSiHit(0); iSiHit<(*detector_)[i].getSiHitVec(idx).size();++iSiHit){
	  G4SiHit lSiHit = (*detector_)[i].getSiHitVec(idx)[iSiHit];
	  HGCSSSimHit lHit(lSiHit,idx);

	  isInserted = lHitMap.insert(std::pair<unsigned,HGCSSSimHit>(lHit.cellid(),lHit));
	  if (!isInserted.second) isInserted.first->second.Add(lSiHit);
	}
	std::map<unsigned,HGCSSSimHit>::iterator lIter = lHitMap.begin();
	hitvec_.reserve(hitvec_.size()+lHitMap.size());
	for (; lIter != lHitMap.end(); ++lIter){
	  (lIter->second).calculateTime();
	  hitvec_.push_back(lIter->second);
	}

      }//loop on 3 sensitive layers

      if(debug) {
	(*detector_)[i].report( (i==0) );
      }
      //if (i==0) G4cout << " ** evt " << evt->GetEventID() << G4endl;
      (*detector_)[i].resetCounters();
    }
  if(debug){
    G4cout << " -- Number of truth particles = " << genvec_.size() << G4endl
	   << " -- Number of simhits = " << hitvec_.size() << G4endl
	   << " -- Number of sampling sections = " << ssvec_.size() << G4endl;
    
  }

  tree_->Fill();
  
  //reset vectors
  genvec_.clear();
  hitvec_.clear();
  ssvec_.clear();
}
