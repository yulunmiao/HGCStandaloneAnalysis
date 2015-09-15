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

#include "effSigmaMacro.C"

int plotSigmaEffvsVtxPos(){//main

  const unsigned nF = 2;
  TFile *fcalib[nF];
  fcalib[0] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalCracks/gitV06-03-04/version100/model3/gamma/eta20_et60_pu0_IC3_Si2.root");
  fcalib[1] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalCracks/gitV06-03-04/version100/model3/gamma/phi_0.440pi/eta20_et60_pu0_IC3_Si2.root");

  const double eta = 2.0;
  const double pT = 60;
  const double Egen = pT*cosh(eta);

  std::string label[nF] = {"phi90","phi79"};

  TFile *fin[nF];
  fin[0] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalCracks/gitV06-03-04/version100/model4/gamma/eta20_et60_pu0_IC3_Si2.root");
  fin[1] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalCracks/gitV06-03-04/version100/model4/gamma/phi_0.440pi/eta20_et60_pu0_IC3_Si2.root");

  const unsigned nP = 62;
  const double stepsize = 310./nP;
  const double xmin = -100+0.33+2.5;

  TCanvas *mycCalib[nF];
  TCanvas *myc[2*nF];
  for (unsigned iF(0); iF<2*nF;++iF){//loop on files
    std::ostringstream lname;
    lname << "mycCalib" << iF;
    if (iF<2) mycCalib[iF] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
    lname.str("");
    lname << "myc" << iF;
    myc[iF] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
  }

  TGraphErrors *gr[nF];

  TH1F *hE[nF][nP];
  TH1F *hECalib[nF];

  for (unsigned iF(0); iF<nF;++iF){//loop on files
    mycCalib[iF]->cd();
    gStyle->SetOptStat("eMRu");
    gStyle->SetOptFit(1111);
    gStyle->SetStatH(0.2);
    //gStyle->SetStatX(0.15);
    fcalib[iF]->cd("Energies");
    
    TTree *treeCalib = (TTree*)gDirectory->Get("Ereso");
    std::ostringstream lcalib;
    lcalib << "Ecalib_" << label[iF];
    hECalib[iF] = new TH1F(lcalib.str().c_str(),";E (Mips); showers",100,20000,28000);
    treeCalib->Draw(("wgtEtotal>>"+lcalib.str()).c_str());
    hECalib[iF]->Fit("gaus");
    TF1 *fit = (TF1*)hECalib[iF]->GetFunction("gaus");
    if (!fit) continue;
    double normalisation = Egen/fit->GetParameter(1);

    double sigmaeff_ref = effSigmaMacro(hECalib[iF]);

    TLatex lat;
    char buf[100];
    sprintf(buf,"E_{gen} = %3.3f GeV",Egen);
    lat.DrawLatexNDC(0.2,0.8,buf);
    sprintf(buf,"#sigma_{eff} = %3.3f Mips",sigmaeff_ref);
    lat.DrawLatexNDC(0.2,0.7,buf);
    sprintf(buf,"#sigma_{eff} = %3.3f GeV",sigmaeff_ref*normalisation);
    lat.DrawLatexNDC(0.2,0.6,buf);
    sprintf(buf,"#sigma_{eff}/E_{gen} = %3.3f",sigmaeff_ref*normalisation/Egen);
    lat.DrawLatexNDC(0.2,0.5,buf);
    mycCalib[iF]->Update();
    mycCalib[iF]->Print(("PlotsCracks/EtotalCalib_"+label[iF]+".pdf").c_str());



    myc[iF]->Divide(8,8);

    fin[iF]->cd("Energies");
    TTree *tree = (TTree*)gDirectory->Get("Ereso");
    double sigeff[nP];
    double errsig[nP];
    double vtx[nP];
    double errvtx[nP];
    for (unsigned iP(0);iP<nP;++iP){//loop on points
      myc[iF]->cd(iP+1);

      std::ostringstream lname,lvar,lcut;
      lname << "Etotal_" << label[iF] << "_step" << iP;

      hE[iF][iP] = new TH1F(lname.str().c_str(),";E (GeV); showers",150,0,300);

      lvar << "wgtEtotal*" << normalisation << ">>" << lname.str();
      lcut << "vtxX>" << xmin+iP*stepsize 
	   << " && vtxX<" << xmin+(iP+1)*stepsize;
      errvtx[iP] = stepsize/2.;
      vtx[iP] = xmin+(iP+0.5)*stepsize;

      tree->Draw(lvar.str().c_str(),lcut.str().c_str());
      sigeff[iP] = effSigmaMacro(hE[iF][iP])/Egen;
      errsig[iP] = hE[iF][iP]->GetRMSError()/Egen;

    }//loop on points

    myc[iF]->Update();
    myc[iF]->Print(("PlotsCracks/EtotalPerVertexBin_"+label[iF]+".pdf").c_str());

    gr[iF] = new TGraphErrors(nP,vtx,sigeff,errvtx,errsig);
    gr[iF]->SetMarkerStyle(22);
    gr[iF]->SetMinimum(0);
    gr[iF]->SetMaximum(0.06);
    gr[iF]->SetTitle(";vtx x (mm); #sigma_{eff}/E_{gen}");
    myc[2+iF]->cd();
    gStyle->SetOptStat(0);
    gr[iF]->Draw("APE");

    TBox *cr1 = new TBox(-61.7,gr[iF]->GetMinimum(),-51.7,gr[iF]->GetMaximum());
    cr1->SetFillColor(7);
    //cr1->SetFillStyle(3004);
    cr1->Draw();
    TBox *cr2 = new TBox(41.7,gr[iF]->GetMinimum(),51.7,gr[iF]->GetMaximum());
    cr2->SetFillColor(7);
    //cr2->SetFillStyle(3004);
    cr2->Draw();
    TBox *cr3 = new TBox(145,gr[iF]->GetMinimum(),155,gr[iF]->GetMaximum());
    cr3->SetFillColor(7);
    //cr3->SetFillStyle(3004);
    cr3->Draw();

    lat.SetTextColor(9);
    lat.SetTextSize(0.04);
    lat.DrawLatex(-70,0.005,"L10-19 crack");
    lat.DrawLatex(30,0.005,"L20-27 crack");
    lat.DrawLatex(130,0.005,"L0-9 crack");
    lat.SetTextColor(1);
    lat.SetTextSize(0.05);


    TLine *ref = new TLine(-100,sigmaeff_ref*normalisation/Egen,210,sigmaeff_ref*normalisation/Egen);
    ref->SetLineColor(2);
    ref->Draw();
    gr[iF]->Draw("PE");

    sprintf(buf,"Single #gamma, E_{gen} = %3.1f GeV, #eta=%3.1f",Egen,eta);
    lat.DrawLatexNDC(0.3,0.8,buf);

    if (iF==1){
      lat.DrawLatexNDC(0.3,0.7,"3^{o} tilt");
    }

    myc[2+iF]->Update();
    myc[2+iF]->Print(("PlotsCracks/SigmaEffvsVtxPos_"+label[iF]+".pdf").c_str());


  }//loop on files



 return 0;
}//main

