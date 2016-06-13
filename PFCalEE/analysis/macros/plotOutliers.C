#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TH2Poly.h"

int plotOutliers(){//main

  const unsigned nEvts = 10;
  int evtNum[nEvts] = {0,1,2,3,4,5,6,7,8,9};

  TFile *sim = TFile::Open("root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalHexa/gittestHexa/gamma/HGcal__version100_model3_BOFF_et50_eta2.000.root");
  TFile *rec = TFile::Open("root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalHexa/gittestHexa/gamma/DigiIC3_Si2__version100_model3_BOFF_et50_eta2.000.root");

  TCanvas *myc= new TCanvas("myc","myc",1500,1000);
  myc->Divide(3,2);
  TCanvas *mycs= new TCanvas("mycs","mycs",1500,1000);
  mycs->Divide(3,2);

  rec->cd();
  TTree *recTree = (TTree*)gDirectory->Get("RecoTree");
  sim->cd();
  TTree *simTree = (TTree*)gDirectory->Get("HGCSSTree");
  unsigned size=0;
  UInt_t cellid[1000];
  unsigned evt[1000];
  double energy[1000];

  //simTree->SetBranchAddress("@HGCSSSimHitVec.size()",&size);
  simTree->SetBranchAddress("HGCSSSimHitVec.cellid_",&cellid[0]);


  unsigned start = 6;
  /*
  TH2Poly *poly = new TH2Poly();
  poly->Honeycomb(-1000,-1000,6.4,180,208);
  //find id->xy conversion
  TIter next(poly->GetBins());
  TObject *obj=0; 
  TH2PolyBin *polyBin = 0;
  std::map<int,std::pair<double,double> > geom;
  
  while ((obj=next())){
    polyBin=(TH2PolyBin*)obj;
    int id = polyBin->GetBinNumber();
    std::pair<double,double> xy = std::pair<double,double>((polyBin->GetXMax()+polyBin->GetXMin())/2.,(polyBin->GetYMax()+polyBin->GetYMin())/2.);
    geom.insert(std::pair<unsigned,std::pair<double,double> >(id,xy));
    }*/
  
  for (unsigned iE(0); iE<nEvts; ++iE){//loop on events
    TH2F *xvsy[6];
    gStyle->SetOptStat(0);
    for (unsigned iL(0); iL<6; ++iL){
      myc->cd(iL+1);
      gPad->SetLogz(1);

      std::ostringstream lname;
      lname << "xvsy_" << iL+start;
      xvsy[iL] = new TH2F(lname.str().c_str(),";x (mm); y(mm); E(mips)",80,-200,200,80,-200,200);
      
      std::ostringstream var,cut;
      var << "HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_>>" << lname.str();
      cut << "(event_==" << evtNum[iE] << " && HGCSSRecoHitVec.layer_==" << iL+start << ")*(HGCSSRecoHitVec.energy_)";
      
      recTree->Draw(var.str().c_str(),cut.str().c_str(),"colz");

    }
    
    myc->Update();
    std::ostringstream lsave;
    lsave << "Layer6to12_rechits_evt" << evtNum[iE] << ".pdf";
    myc->Print(lsave.str().c_str());


    simTree->GetEntry(evtNum[iE]);
    simTree->Show(iE);
    std::cout << "Number of simhits " << size << " cellid0 = " << cellid[0] << std::endl;

    return 1;

    mycs->cd();    
    mycs->Update();
    lsave.str("");
    lsave << "Layer6to12_simhits_evt" << evtNum[iE] << ".pdf";
    mycs->Print(lsave.str().c_str());

  }//loop on events

  return 0;
}//main
