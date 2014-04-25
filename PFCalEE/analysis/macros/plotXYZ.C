#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TColor.h"

void set_plot_style()
{
    const Int_t NRGBs = 7;
    const Int_t NCont = 64;

    //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    //Double_t stops[NRGBs] = { 0.00, 0.20, 0.30, 0.45, 0.61, 0.84, 1.00 };
    Double_t stops[NRGBs] = { 0.00, 0.10, 0.20, 0.35, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 1.00, 0.00, 0.00, 0.50, 1.00, 0.50 };
    Double_t green[NRGBs] = { 1.00, 1.00, 1.00, 0.81, 0.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.50, 0.00, 0.10, 1.00, 0.50, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

int plotXYZ(){//main  

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 2;
  std::string scenario[nS] = {
    //"quark_u/eta30/",
    //"quark_u/PU/eta30/"
    //"gamma_eta25/GeVCal/",
    //"gamma_eta25/PU/GeVCal/"
    //"VBFH/concept/",
    //"VBFH/PU/concept/"
    //"ZH125/concept/",
    //"ZH125/PU/concept/"
    "GluGlu/concept/",
    "GluGlu/PU/concept/"
  };

  bool doVBF = false;
  bool doGlu = true;

  std::ostringstream ltitleBase;
  //ltitle << "#gamma 200 GeV"
  if (doVBF) ltitleBase << "VBF jet";
  else if (doGlu) ltitleBase << "Gluon jet";
  else ltitleBase << "ZH,H#rightarrow#tau#tau";
  //<< " Event #" << event[ievt]
  ltitleBase << ", E_{1#times1 cm^{2}} > ";

  const unsigned nV = 1;
  TString version[nV] = {"20"};
  //double Emip = 0.0548;//in MeV
  // double Emip[nLayers];
  // for (unsigned iL(0);iL<nLayers;++iL){
  //   if (nLayers < nEcalLayers) Emip[iL] = 0.0548;//in MeV
  //   else if (nLayers < nEcalLayers+24)  Emip[iL] = 0.0849;
  //   else Emip[iL] = 1.49;
  // }

  const unsigned nmips = 3;
  unsigned mipThreshBase[nmips] = {1,5,10};
  bool doThresh = true;

  bool do3by3 = true;

  double minRadius = 0;//316;//mm

  //VBFH
  //const unsigned nEvts = 7;//1000;
  //unsigned event[nEvts]={4,5,6,7,9,11,12};//,6,12};//103,659,875};
  //Htautau
  //const unsigned nEvts = 7;//1000;
  //unsigned event[nEvts]={1,2,5,8,10,11,12};//,6,12};//103,659,875};
  //Htautau 125
  //const unsigned nEvts = 2;//1000;
  //unsigned event[nEvts]={10,14};//,6,12};//103,659,875};
  //gluons
  const unsigned nEvts = 3;//1000;
  unsigned event[nEvts]={7,16,22};//,6,12};//103,659,875};


  //for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
  //  event[ievt] = ievt;
  // }

  const unsigned nLayers = 64;
  const unsigned nEcalLayers = 31;
  unsigned minLayer3by3 = 21;
  if (!do3by3) minLayer3by3 = nLayers;

  double minZ=3170,maxZ=5000;
  unsigned nZ=(maxZ-minZ)/2.;
  //double minX=-150,maxX=150;

  //VBFH
  //double minX[nEvts] = {-700,-600,-400,-400,0,0,0};
  //double maxX[nEvts] = {-100,0,200,100,700,700,500};
  //double minY[nEvts] = {100,-600,-1100,-800,-100,-200,-700};
  //double maxY[nEvts] = {480,0,0,0,400,300,0};

  //double minXecal[nEvts] = {-600,-490,-400,-300,450,320,50};
  //double maxXecal[nEvts] = {-350,-290,200,-100,650,520,350};
  //double minYecal[nEvts] = {200,-490,-1100,-650,100,-100,-580};
  //double maxYecal[nEvts] = {450,-290,0,-450,300,100,-280};

  //Htautau
  /*double minX[nEvts] = {0,-500,0,-400,0,-500,-700};
  double maxX[nEvts] = {600,0,700,500,700,0,600};
  double minY[nEvts] = {0,-700,-100,-700,-400,0,-400};
  double maxY[nEvts] = {500,0,400,900,200,700,700};

  double minXecal[nEvts] = {350,-400,450,-300,400,-400,-700};
  double maxXecal[nEvts] = {550,-200,650,-100,600,-200,600};
  double minYecal[nEvts] = {180,-550,20,600,-200,400,-400};
  double maxYecal[nEvts] = {380,-300,220,800,0,600,700};*/

  //Htautau 125
  /*double minX[nEvts] = {-300,-850,-200,-450,-550};
  double maxX[nEvts] = {100,-450,300,0,-150};
  double minY[nEvts] = {-900,-280,200,-950,-800};
  double maxY[nEvts] = {-500,150,650,-550,-400};

  double minXecal[nEvts] = {-190,-770,-30,-330,-460};
  double maxXecal[nEvts] = {10,-570,170,-130,-260};
  double minYecal[nEvts] = {-800,-160,320,-840,-720};
  double maxYecal[nEvts] = {-600,40,520,-640,-520};

  double minX[nEvts] = {-450,-550};
  double maxX[nEvts] = {50,-150};
  double minY[nEvts] = {-1100,-950};
  double maxY[nEvts] = {-700,-550};

  double minXecal[nEvts] = {-350,-480};
  double maxXecal[nEvts] = {-150,-280};
  double minYecal[nEvts] = {-900,-760};
  double maxYecal[nEvts] = {-700,-560};*/


  //pp to gg
  /*double minX[nEvts] = {-600,-250};
  double maxX[nEvts] = {0,350};
  double minY[nEvts] = {100,250};
  double maxY[nEvts] = {700,850};

  double minXecal[nEvts] = {-600,-250};
  double maxXecal[nEvts] = {0,350};
  double minYecal[nEvts] = {100,250};
  double maxYecal[nEvts] = {700,850};*/

  double minX[nEvts] = {100,-100,100};
  double maxX[nEvts] = {700,500,700};
  double minY[nEvts] = {100,200,100};
  double maxY[nEvts] = {700,800,700};
 
  double minXecal[nEvts] = {100,-100,100};
  double maxXecal[nEvts] = {700,500,700};
  double minYecal[nEvts] = {100,200,100};
  double maxYecal[nEvts] = {700,800,700};

  //double minX=-700,maxX=700;
  //double minY=420,maxY=720;
  //double minY=1150,maxY=1450;
  //double minY=-1100,maxY=500;

  //unsigned nX[nEvts],nY[nEvts];

  bool doAll = false;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
    if (doAll){
      minX[ievt]=-1700;maxX[ievt]=1700;
      minY[ievt]=-1700;maxY[ievt]=1700;
    }
    //nX[ievt]=(maxX[ievt]-minX[ievt])/10;
    //nY[ievt]=(maxY[ievt]-minY[ievt])/10;
  }


  TCanvas *mycECAL = new TCanvas("mycECAL","mycECAL",1200,800);
  TCanvas *mycECALzoom = new TCanvas("mycECALzoom","mycECALzoom",1200,800);
  TCanvas *mycHCAL = new TCanvas("mycHCAL","mycHCAL",1200,800);
  //TCanvas *mycAll = new TCanvas("mycAll","mycAll",1200,800);
  //const unsigned nPads = nEvts>10 ? 10 : nEvts;
  //mycAll->Divide(static_cast<unsigned>(nPads/2.+0.5),nEvts/2<1 ? 1 : 2);

  const unsigned nCanvas = nS;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1300,800);
  }

  std::ostringstream saveName;
  
  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
      
      TString inputDir = "../PLOTS/version_"+version[iV]+"/"+scenario[iS]+"/";
      
      
      bool isRECO = false;
      //if (scenario[iS].find("scenario_") != scenario[iS].npos) isRECO=true;
      
      if (isRECO) inputDir += "Reco/";

      unsigned evtcounter = 0;

      for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
	std::ostringstream lName;
	lName << inputDir << "CalibHistos_E200_evt" << event[ievt] << ".root";
	TFile *inputFile = TFile::Open(lName.str().c_str());
	if (!inputFile) {
	  std::cout << " -- Error, input file " << lName.str() << " cannot be opened. Going to next..." << std::endl;
	  continue;
	  //return 1;
	}
	else {
	  std::cout << " -- Processing event " << event[ievt] << std::endl;
	}
	TH2F *p_xy[nLayers];
	TH2F *p_xytot = 0;

	TString plotDir = inputDir;
	if (doAll) plotDir += "Overview/";

	// lName.str("");
	// lName << "../PLOTS/version_8/gamma_eta17/PU_pipm/GeVCal/CalibHistos_E200_evt" << event[ievt] << ".root";
	// TFile *inputproton = TFile::Open(lName.str().c_str());
	// TH3F *p_xyz_p = (TH3F*)gDirectory->Get("p_xyz")->Clone();
	// lName.str("");
	// lName << "../PLOTS/version_8/gamma_eta17/PU_pi0/GeVCal/CalibHistos_E200_evt" << event[ievt] << ".root";
	// TFile *inputneutron = TFile::Open(lName.str().c_str());
	// TH3F *p_xyz_n = (TH3F*)gDirectory->Get("p_xyz")->Clone();

	inputFile->cd();
	TH3F *p_xyz = 0;
	if (!isRECO) p_xyz = (TH3F*)gDirectory->Get("p_xyz")->Clone();
	else p_xyz = (TH3F*)gDirectory->Get("p_recoxyz")->Clone();
      
	if (!p_xyz) {
	  std::cout << " -- ERROR, pointer for XYZ histogram is null. Exiting..." << std::endl;
	  return 1;
	}
	p_xyz->Sumw2();


	//p_xyz_p->Sumw2();
	//p_xyz_n->Sumw2();
	double EmaxEcal = 0;
	double EmaxHcal = 0;
	//p_xyz->Scale(1./Emip);
	//p_xyz_p->Scale(1./Emip);
	//p_xyz_n->Scale(1./Emip);
	//p_xyz->SetMinimum(100);

	for (unsigned iL(0); iL<nLayers; ++iL){
	  //if (nLayers > nEcalLayers) Emip = 300./200*Emip;
	  lName.str("");
	  if (!isRECO) lName << "p_xy_" << iL;
	  else lName << "p_recoxy_" << iL;
	  p_xy[iL] = (TH2F*)gDirectory->Get(lName.str().c_str());
	  if (!p_xy[iL]) {
	    std::cout << " -- ERROR, pointer for histogram is null for layer: " << iL << ". Exiting..." << std::endl;
	    return 1;
	  }
	  
	  //p_xy[iL]->Scale(1./Emip);
	  if (!isRECO) {
	    p_xy[iL]->RebinX(4);
	    p_xy[iL]->RebinY(4);
	  }
	  
	  double Etot = p_xy[iL]->GetMaximum();
	  if (Etot > EmaxEcal && iL<nEcalLayers) EmaxEcal = Etot;
	  if (Etot > EmaxHcal && iL>=nEcalLayers) EmaxHcal = Etot;
	  
	  if (iL==0){
	    p_xytot = (TH2F*)p_xy[iL]->Clone();
	  }
	  else {
	    if (!doThresh){
	      p_xytot->Add(p_xy[iL]);
	    }
	  }
	}//loop on layers

	std::cout << " -- max E found: ECAL: " <<  EmaxEcal << ", HCAL: " << EmaxHcal << std::endl;

	TH2F *p_xztot = 0;
	TH2F *p_zytot = 0;

	unsigned mipThresh = 1;
	for (unsigned iM(0); iM<nmips;++iM){//loop on mip thresh values
	  
	  set_plot_style();
	  if (iS==1) {
	    mipThresh = mipThreshBase[iM];     
	  }
	  else if (iS==0){
	    if (iM>0) continue;
	  }

	  for (int xb(1); xb<p_xyz->GetNbinsX()+1;++xb){
	    for (int yb(1); yb<p_xyz->GetNbinsY()+1;++yb){
	      double xtmp=p_xyz->GetYaxis()->GetBinCenter(yb);
	      for (int zb(1); zb<p_xyz->GetNbinsZ()+1;++zb){
		double ytmp=p_xyz->GetZaxis()->GetBinCenter(zb);
		double radius = sqrt(xtmp*xtmp+ytmp*ytmp);
		if (radius<minRadius){
		  p_xyz->SetBinContent(xb,yb,zb,0);
		  continue;
		}
		double ltmp = p_xyz->GetBinContent(xb,yb,zb);
		if (ltmp<1) continue;
		//std::cout << xb << " " << yb << " " << zb << " " << p_xyz->GetBinContent(xb,yb,zb) << std::endl;
		if (ltmp < mipThresh) {
		  p_xyz->SetBinContent(xb,yb,zb,0);
		}
		else {
		  p_xyz->SetBinContent(xb,yb,zb,ltmp);
		}
	      }
	    }
	  }
	  if (p_xztot) p_xztot->Delete();
	  p_xztot = new TH2F("p_xztot",";x(mm);z(mm)",
			     p_xyz->GetNbinsY(),p_xyz->GetYaxis()->GetBinLowEdge(1),p_xyz->GetYaxis()->GetBinLowEdge(p_xyz->GetNbinsY()+1),
			     p_xyz->GetNbinsX(),p_xyz->GetXaxis()->GetBinLowEdge(1),p_xyz->GetXaxis()->GetBinLowEdge(p_xyz->GetNbinsX()+1));
	  if (p_zytot) p_zytot->Delete();
	  p_zytot = new TH2F("p_zytot",";z(mm);y(mm)",
			     p_xyz->GetNbinsX(),p_xyz->GetXaxis()->GetBinLowEdge(1),p_xyz->GetXaxis()->GetBinLowEdge(p_xyz->GetNbinsX()+1),
			     p_xyz->GetNbinsZ(),p_xyz->GetZaxis()->GetBinLowEdge(1),p_xyz->GetZaxis()->GetBinLowEdge(p_xyz->GetNbinsZ()+1));


	  if (doThresh){
	    int xbmin = p_xyz->GetYaxis()->FindBin(minX[ievt]);
	    int xbmax = p_xyz->GetYaxis()->FindBin(maxX[ievt]);
	    int ybmin = p_xyz->GetZaxis()->FindBin(minY[ievt]);
	    int ybmax = p_xyz->GetZaxis()->FindBin(maxY[ievt]);
	    int zbmin = p_xyz->GetXaxis()->FindBin(minZ);
	    int zbmax = p_xyz->GetXaxis()->FindBin(maxZ);
	    for (unsigned iL(0); iL<nLayers; ++iL){
	      for (int xb(1); xb<p_xyz->GetNbinsY()+1;++xb){
		double xtmp=p_xyz->GetYaxis()->GetBinCenter(xb);
		for (int yb(1); yb<p_xyz->GetNbinsZ()+1;++yb){
		  double ytmp=p_xyz->GetZaxis()->GetBinCenter(yb);
		  double radius = sqrt(xtmp*xtmp+ytmp*ytmp);
		  for (int zb(1); zb<p_xyz->GetNbinsX()+1;++zb){
		    //std::cout << xb << " " << yb << " " << zb << " " << p_xyz->GetBinContent(xb,yb,zb) << std::endl;
		    if (zb==1) {
		      if (iL==0) {
			p_xytot->SetBinContent(xb,yb,0);
		      }
		      if (radius<minRadius){
			p_xy[iL]->SetBinContent(xb,yb,0);
		      }
		      if (p_xy[iL]->GetBinContent(xb,yb) < mipThresh) {
			p_xytot->SetBinContent(xb,yb,p_xytot->GetBinContent(xb,yb));
		      }
		      else {
			p_xytot->SetBinContent(xb,yb,p_xytot->GetBinContent(xb,yb)+p_xy[iL]->GetBinContent(xb,yb));
		      } 
		    }
		    if (iL==0){
		      if (yb>ybmin && yb < ybmax) {
			p_xztot->SetBinContent(xb,zb,p_xztot->GetBinContent(xb,zb)+p_xyz->GetBinContent(zb,xb,yb));
			if (zb>1 && p_xztot->GetBinContent(xb,zb)>0 && p_xztot->GetBinContent(xb,zb-1)>0){
			  p_xztot->SetBinContent(xb,zb-1,p_xztot->GetBinContent(xb,zb-1)+p_xztot->GetBinContent(xb,zb));
			  p_xztot->SetBinContent(xb,zb,0);
			}
		      }
		      if (xb>xbmin && xb < xbmax) {
			p_zytot->SetBinContent(zb,yb,p_zytot->GetBinContent(zb,yb)+p_xyz->GetBinContent(zb,xb,yb));
			if (zb>1 && p_zytot->GetBinContent(zb,yb)>0 && p_zytot->GetBinContent(zb-1,yb)>0){
			  p_zytot->SetBinContent(zb-1,yb,p_zytot->GetBinContent(zb-1,yb)+p_zytot->GetBinContent(zb,yb));
			  p_zytot->SetBinContent(zb,yb,0);
			}
		      }
		    }
		  }//loop on z
		}//loop on y
	      }//loop on x
	    }//loop on layers
	  }//doThresh

	  gStyle->SetOptStat(0);
	  std::ostringstream ltitle;
	  ltitle << ltitleBase.str() << mipThresh << " Mips" ;
	
	  myc[iS]->Clear();
	  myc[iS]->Divide(2,2);
	  myc[iS]->cd(1);
	  if (!isRECO) {
	    //p_xyz->RebinY(4);
	    //p_xyz->RebinZ(4);
	    //p_xyz_p->RebinY(4);
	    //p_xyz_p->RebinZ(4);
	    //p_xyz_n->RebinY(4);
	    //p_xyz_n->RebinZ(4);
	  }
	  p_xyz->GetXaxis()->SetRangeUser(minZ,maxZ);
	  p_xyz->GetYaxis()->SetRangeUser(minX[ievt],maxX[ievt]);
	  p_xyz->GetZaxis()->SetRangeUser(minY[ievt],maxY[ievt]);

	  p_xyz->GetXaxis()->SetLabelSize(0.05);
	  p_xyz->GetYaxis()->SetLabelSize(0.05);
	  p_xyz->GetZaxis()->SetLabelSize(0.05);
	  p_xyz->GetXaxis()->SetTitleSize(0.05);
	  p_xyz->GetYaxis()->SetTitleSize(0.05);
	  p_xyz->GetZaxis()->SetTitleSize(0.05);
	  p_xyz->SetTitle(ltitle.str().c_str());
	  p_xyz->Draw("");
	  //if (iS==1){
	    //p_xyz_p->SetMarkerColor(2);
	    //p_xyz_p->SetMaximum(p_xyz->GetMaximum());
	    //p_xyz_p->Draw("same");
	    //p_xyz_n->SetMarkerColor(3);
	    //p_xyz_n->SetMaximum(p_xyz->GetMaximum());
	    //p_xyz_n->Draw("same");
	  //}

	  myc[iS]->cd(2);
	  gPad->SetLogz(1);
	  //p_xztot->RebinX(4);
	  //p_xztot->RebinY(2);
	  p_xztot->GetXaxis()->SetLabelSize(0.05);
	  p_xztot->GetYaxis()->SetLabelSize(0.05);
	  p_xztot->GetXaxis()->SetTitleSize(0.05);
	  p_xztot->GetYaxis()->SetTitleSize(0.05);
	  p_xztot->GetZaxis()->SetTitleOffset(-0.5);
	  p_xztot->GetZaxis()->SetTitle("E(MIP)");
	  p_xztot->GetXaxis()->SetRangeUser(minX[ievt],maxX[ievt]);
	  p_xztot->GetYaxis()->SetRangeUser(minZ,maxZ);
	  p_xztot->Draw("colz");
	  //p_xztot->Draw("lego2z");
	  //p_xztot->Draw("");
	
	  myc[iS]->cd(3);
	  gPad->SetLogz(1);
	  //p_zytot->RebinX(2);
	  //p_zytot->RebinY(4);
	  p_zytot->GetXaxis()->SetLabelSize(0.05);
	  p_zytot->GetYaxis()->SetLabelSize(0.05);
	  p_zytot->GetXaxis()->SetTitleSize(0.05);
	  p_zytot->GetYaxis()->SetTitleSize(0.05);
	  p_zytot->GetZaxis()->SetTitleOffset(-0.5);
	  p_zytot->GetZaxis()->SetTitle("E(MIP)");
	  p_zytot->GetXaxis()->SetRangeUser(minZ,maxZ);
	  p_zytot->GetYaxis()->SetRangeUser(minY[ievt],maxY[ievt]);
	  p_zytot->Draw("colz");
	  //p_zytot->Draw("lego2z");
	  //p_zytot->Draw("");
	
	  myc[iS]->cd(4);
	  gPad->SetLogz(1);
	  p_xytot->GetXaxis()->SetLabelSize(0.05);
	  p_xytot->GetYaxis()->SetLabelSize(0.05);
	  p_xytot->GetXaxis()->SetTitleSize(0.05);
	  p_xytot->GetYaxis()->SetTitleSize(0.05);
	  p_xytot->GetZaxis()->SetTitleOffset(-0.5);
	  p_xytot->GetZaxis()->SetTitle("E(MIP)");
	  p_xytot->GetXaxis()->SetRangeUser(minX[ievt],maxX[ievt]);
	  p_xytot->GetYaxis()->SetRangeUser(minY[ievt],maxY[ievt]);
	  //p_xytot->Draw("colz");
	  p_xytot->Draw("lego2z");

	  myc[iS]->Update();
	  saveName.str("");
	  saveName << plotDir << "/xyzSimHits_evt" << event[ievt] << "_mipThresh" << mipThresh;
	  if (!do3by3) saveName << "_1by1";
	  myc[iS]->Print((saveName.str()+".png").c_str());
	  myc[iS]->Print((saveName.str()+".pdf").c_str());
	
	  //return 1;//
	  //continue;

	  // mycAll->cd(evtcounter%nPads+1);
	  // gPad->SetLogz(1);
	  // //TH2D *proj = (TH2D*)p_xyz->Project3D(lLabel);
	  // p_xztot->GetZaxis()->SetLabelSize(0.05);
	  // p_xztot->GetZaxis()->SetTitleSize(0.05);
	  // p_xztot->GetZaxis()->SetTitle("E(MIP)");
	  // p_xztot->GetZaxis()->SetTitleOffset(-0.5);
	  // p_xztot->SetTitle(ltitle.str().c_str());
	  // //p_xztot->SetMinimum(0.1);
	  // //p_xztot->Draw("LEGO2z");
	  // p_xztot->Draw("colz");

	  // if (evtcounter%nPads == nPads-1){
	  //   mycAll->Update();
	  //   saveName.str("");
	  //   saveName << plotDir << "/xzSimHits_" << evtcounter/nPads << "_mipThresh" << mipThresh;
	  //   mycAll->Print((saveName.str()+".png").c_str());
	  //   mycAll->Print((saveName.str()+".pdf").c_str());
	  // }
	  //continue;
	}//loop on mi pthreshold

	mycECAL->Clear();
	unsigned counter = 1;
	mycECAL->Divide(3,3);
	mycECALzoom->Clear();
	mycECALzoom->Divide(3,3);

	mycHCAL->Clear();
	mycHCAL->Divide(5,3);

	unsigned layId = 0;
	for (unsigned iL(1); iL<nLayers; ++iL){//loop on layers
	  //std::cout << " -- Processing layer " << iL << std::endl;
	  if (iL<minLayer3by3){
	    p_xy[iL]->RebinX(2);
	    p_xy[iL]->RebinY(2);
	  }
	  else {
	    p_xy[iL]->RebinX(3);
	    p_xy[iL]->RebinY(3);
	  }
	  if ((iL-1)%3==2 && counter < 10 && iL<nEcalLayers) {
	    std::cout << " -- Ecal layer " << iL << " counter = " << counter << std::endl;
	    mycECAL->cd(counter);
	    counter++;
	    layId = iL;
	  }
	  else if (counter<25 && iL>=nEcalLayers && ((iL<nEcalLayers+24 && iL%2==1)||iL>=nEcalLayers+24))  {
	    std::cout << " -- Hcal layer " << iL << " counter = " << counter << std::endl;
	    mycHCAL->cd(counter-9);
	    layId = nEcalLayers+counter-10;
	    if (iL>=nEcalLayers+24) layId = nEcalLayers+15+counter-10;
	    counter++;
	  }
	  else if (iL==nEcalLayers+24+3){
	    saveName.str("");
	    saveName << plotDir << "/xySimHits_HCAL_evt" << event[ievt] << "_subset1";
	    if (!do3by3) saveName << "_1by1";
	    mycHCAL->Update();
	    mycHCAL->Print((saveName.str()+".png").c_str());
	    mycHCAL->Print((saveName.str()+".pdf").c_str());
	    counter = 10;
	    std::cout << " -- Hcal layer " << iL << " counter = " << counter << std::endl;
	    mycHCAL->Clear();
	    mycHCAL->Divide(2,3);
	    mycHCAL->cd(counter-9);
	    layId = nEcalLayers+15+counter-10;
	    counter++;
	  }
	  else continue;
	  gPad->SetLogz(1);
	  p_xy[iL]->SetMaximum(iL<nEcalLayers ? EmaxEcal : EmaxHcal);
	  p_xy[iL]->GetXaxis()->SetLabelSize(0.05);
	  p_xy[iL]->GetYaxis()->SetLabelSize(0.05);
	  p_xy[iL]->GetXaxis()->SetTitleSize(0.05);
	  p_xy[iL]->GetYaxis()->SetTitleSize(0.05);
	  p_xy[iL]->GetZaxis()->SetTitleOffset(-0.5);
	  p_xy[iL]->GetZaxis()->SetTitle("E(MIP)");

	  //p_xy[iL]->Draw("colz");
	  char buf[500];
	  sprintf(buf,"Layer %d",layId-1);
	  p_xy[iL]->SetTitle(buf);

	  //p_xy[iL]->Draw("LEGO2z");
	  p_xy[iL]->GetXaxis()->SetRangeUser(minX[ievt],maxX[ievt]);
	  p_xy[iL]->GetYaxis()->SetRangeUser(minY[ievt],maxY[ievt]);
	  //p_xy[iL]->Draw("colz");
	  // if (layId-1==8 || layId-1==14 || 
	  //     layId-1==20 || layId-1==26 ||
	  //     (layId>30 && layId%2==0)
	  //     || layId>42
	  //     ) 

	  gPad->SetTheta(60); // default is 30
	  gPad->SetPhi(40); // default is 30
	  gPad->Update();
	  p_xy[iL]->Draw("lego2z");

	}//loop on layers
	saveName.str("");
	saveName << plotDir << "/xySimHits_ECAL_evt" << event[ievt] << "_subset";
	if (!do3by3) saveName << "_1by1";
	mycECAL->Update();
	mycECAL->Print((saveName.str()+".png").c_str());
	mycECAL->Print((saveName.str()+".pdf").c_str());
	saveName.str("");
	saveName << plotDir << "/xySimHits_HCAL_evt" << event[ievt] << "_subset2";
	if (!do3by3) saveName << "_1by1";
	mycHCAL->Update();
	mycHCAL->Print((saveName.str()+".png").c_str());
	mycHCAL->Print((saveName.str()+".pdf").c_str());


	//plot ECAL zoom
	counter=1;
	for (unsigned iL(1); iL<nEcalLayers; ++iL){//loop on ecal layers
          //std::cout << " -- Processing layer " << iL << std::endl;
          if ((iL-1)%3==2 && counter < 10) {
	    std::cout << " -- ZOOM Ecal layer " << iL << " counter = " << counter << std::endl;
            mycECALzoom->cd(counter);
            counter++;
            layId = iL;
          }
          else continue;
          gPad->SetLogz(1);
          p_xy[iL]->SetMaximum(iL<nEcalLayers ? EmaxEcal : EmaxHcal);
          p_xy[iL]->GetXaxis()->SetLabelSize(0.05);
          p_xy[iL]->GetYaxis()->SetLabelSize(0.05);
          p_xy[iL]->GetXaxis()->SetTitleSize(0.05);
          p_xy[iL]->GetYaxis()->SetTitleSize(0.05);
          p_xy[iL]->GetZaxis()->SetTitleOffset(-0.5);
          p_xy[iL]->GetZaxis()->SetTitle("E(MIP)");

          //p_xy[iL]->Draw("colz");
          char buf[500];
          sprintf(buf,"Layer %d",layId-1);
          p_xy[iL]->SetTitle(buf);

          //p_xy[iL]->Draw("LEGO2z");
          p_xy[iL]->GetXaxis()->SetRangeUser(minXecal[ievt],maxXecal[ievt]);
          p_xy[iL]->GetYaxis()->SetRangeUser(minYecal[ievt],maxYecal[ievt]);
          p_xy[iL]->Draw("colz");

        }//loop on layers
        saveName.str("");
        saveName << plotDir << "/xySimHits_ECAL_evt" << event[ievt] << "_subset";
        if (!do3by3) saveName << "_1by1";
	saveName << "_zoom";
        mycECALzoom->Update();
        mycECALzoom->Print((saveName.str()+".png").c_str());
        mycECALzoom->Print((saveName.str()+".pdf").c_str());


	evtcounter++;

	//return 1;

      }//loop on events
 
    }//loop on scenarios

  }//loop on versions
  
  return 0;


}//main
