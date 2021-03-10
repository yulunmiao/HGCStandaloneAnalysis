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
#include "TProfile.h"
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


unsigned fitEnergy(TH1F *hist,
		   TPad *pad,
		   std::string unitStr,
		   FitResult & lres,
		   unsigned isr,
		   const bool useSigmaEff){
  
  pad->cd();
  //double eMin = hist->GetMean()-5*hist->GetRMS();
  //double eMax = hist->GetMean()+5*hist->GetRMS();
  //hist->GetXaxis()->SetRangeUser(eMin,eMax);
  hist->Draw("PE");

  double sigmaeff = effSigmaMacro(hist);

  double nRMSm = isr<1? 1 : 2.5;
  double nRMSp = 2.5;
  
  TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  double tmpmean = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
  /*  if (tmpmean<3) {
    //find next maximum
    std::cout << " !!ERROR!! mismatch of events leading to peak at 0 !! " << std::endl;
    double max = 0;
    for (unsigned ix(4); ix<hist->GetNbinsX()+1; ++ix){
      if (hist->GetBinContent(ix)>max) {
	max = hist->GetBinContent(ix);
	tmpmean = hist->GetXaxis()->GetBinCenter(ix);
      }
    }
    }*/
  fitResult->SetParameters(hist->Integral()/10., tmpmean, hist->GetRMS());
  std::cout << "HIST Integral: "  << hist->Integral() << std::endl;
  std::cout << " Initial params: "  << fitResult->GetParameter(0) << " "<< fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
	    << std::endl;

  int status = hist->Fit("fitResult","BL0QEMI","",
			 fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
			 fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
  
  
  std::cout << " First fit: " << status << " " << fitResult->GetParameter(0) << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
	    << std::endl;
  
  if ((status != 0 && status != 4000) || fitResult->GetChisquare()/fitResult->GetNDF()>20){
    std::cout << " -- Bad fit ! Try again..." << std::endl;
    status = hist->Fit("fitResult","BL0QEMI","",
   		       fitResult->GetParameter(1)-1*fitResult->GetParameter(2),
		       fitResult->GetParameter(1)+2*fitResult->GetParameter(2));
    std::cout << " Second fit: " << status << " " << fitResult->GetParameter(0) << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
	    << std::endl;

  }
  fitResult->SetLineColor(2);
  fitResult->Draw("same");
  
  if (status != 0 && status != 4000) {
    std::cout << " Warning! Fit failed with status " << status << "! Please have a look at the verbose output below...." << std::endl;
    hist->Fit("fitResult","BL0EMI","",
	      fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
	      fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
    //totalE for pu140 is expected to be pathological :/
    //if (isr!=7) return 1;
  }

  char buf[500];
  TLatex lat;
  lat.SetTextSize(0.04);
  double latx = hist->GetXaxis()->GetXmin()+(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin())/40.;
  double laty = hist->GetMaximum();
  sprintf(buf,"<Efit> = %3.3f +/- %3.3f %s",fitResult->GetParameter(1),fitResult->GetParError(1),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.9,buf);
  sprintf(buf,"RMSfit = %3.3f +/- %3.3f %s",fitResult->GetParameter(2),fitResult->GetParError(2),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.8,buf);
  sprintf(buf,"#sigma_{eff} = %3.3f +/- %3.3f %s",sigmaeff,hist->GetRMSError(),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.7,buf);
  sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
  lat.DrawLatex(latx,laty*0.6,buf);
  
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
  lat.DrawLatex(latx,laty*0.5,buf);
  
  lres.chi2 = fitResult->GetChisquare();
  lres.ndf = fitResult->GetNDF();
  lres.mean = fitResult->GetParameter(1);
  lres.meanerr = fitResult->GetParError(1);
  lres.sigma = useSigmaEff ? sigmaeff : fitResult->GetParameter(2);
  std::cout << "!!!!!!!!! " << std::endl;
  std::cout << "SIGMA: " << sigmaeff << ", " << fitResult->GetParameter(2) << std::endl;
  std::cout << "!!!!!!!!! " << std::endl;
  lres.sigmaerr = useSigmaEff ? hist->GetRMSError() : fitResult->GetParError(2);

  return 0;
};
