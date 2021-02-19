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

//ORDER MATTERS...
#include "TDRStyle.h"
#include "plotUtilities.C"
#include "fitEnergy.C"
#include "makeBackLeakCor.C"
#include "makePuSigma.C"
#include "makeCalibration.C"
#include "makeResolution.C"

//#include "extractResultsvsEta.C"


int makeEfit(const bool useSigmaEff,
	     const bool dovsE,
	     const bool doRaw,
	     const bool doBackLeakCor,
	     const unsigned nBack,
	     const unsigned pu,
	     const unsigned ICval,
	     const std::string scenario,
	     const TString version,
	     const unsigned nLayers,
	     const unsigned eta,
	     const unsigned pT,
	     const unsigned iSR,
	     const double radius,
	     const double & offset,
	     const double & calib,
	     const double & backLeakCor,
	     TCanvas * mycE,
	     TTree *ltree,
	     TFile *foutEfit,
	     TString plotDir,
	     TGraphErrors *calibRecoFit,
	     TGraphErrors *resoRecoFit
	     ){//main


  double Eval = E(pT,eta);
  double etaval = eta/10.;
  double etaerr = 0;

  const unsigned nEvtMin = 150;

  TString pSuffix = doBackLeakCor?"backLeakCor":"";


  double fitQual = iSR==0?50:30;


  TH1F *p_chi2ndf;
  std::ostringstream label;
  label << "p_chi2ndf_E_" << pT << "_SR" << iSR;
  p_chi2ndf = new TH1F(label.str().c_str(),";#chi^{2}/N;entries",500,0,50);
  p_chi2ndf->StatOverflows();
 
  
  std::string unit = "MIPS";
  if (!doRaw) unit = "GeV";
  
  //identify valid energy values
  
  TH1F *p_Ereco;
  TH2F *p_ErecovsEback;
	      
  mycE->cd();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
		
  std::ostringstream lName;
  std::string lNameTot,lNameBack;

  getTotalEnergyString(nLayers,nBack,lNameTot,lNameBack,iSR);

  //std::cout << " Check totE = " << lNameTot << std::endl;
  //std::cout << " Check backE = " << lNameBack << std::endl;

  lName.str("");
  lName << std::setprecision(6);
  
  if (!doRaw) lName << "(";
  lName << lNameTot;
  if (!doRaw) lName << " - " << offset << ")/" << calib;

  if (doBackLeakCor){// && E>300) {
    lName << " - " << backLeakCor << "*(" << lNameBack << ")/(" << lNameTot << ")";
  }
  
  //cut low outliers....
  std::ostringstream lcut;
  if (iSR>0 && !doRaw){
    lcut << lName.str() << ">0.7*" << Eval;
  }
  
  ltree->Draw(lName.str().c_str(),lcut.str().c_str(),"");
  TH1F  *hTmp =  (TH1F*)gPad->GetPrimitive("htemp");
  if (!hTmp) {
    std::cout << " -- Problem with tree Draw for "
	      <<  lName.str() << " cut " << lcut.str()
	      << " -> not drawn. Exiting..." << std::endl;
    return 1;
  }

  //get mean and rms, and choose binning accordingly then redraw.
  TH1F *hE = new TH1F("hE","",40,hTmp->GetMean()-6*hTmp->GetRMS(),hTmp->GetMean()+6*hTmp->GetRMS());

  ltree->Draw((lName.str()+">>hE").c_str(),lcut.str().c_str(),"");

  lName.str("");
  lName << "energy" << pT << "_SR" << iSR ;
  p_Ereco = (TH1F*)hE->Clone(lName.str().c_str()); // 1D
  
  if (!p_Ereco){
    std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null." << std::endl;
    return 1;
  }
  std::cout << " --- Reco E = entries " << p_Ereco->GetEntries() 
	    << " mean " << p_Ereco->GetMean() 
	    << " rms " << p_Ereco->GetRMS() 
	    << " overflows " << p_Ereco->GetBinContent(p_Ereco->GetNbinsX()+1)
	    << std::endl;
  
  
  p_Ereco->SetTitle((";E ("+unit+");events").c_str());
  
  //take min 20 bins
  //if(p_Ereco->GetNbinsX()>40) p_Ereco->Rebin(2);
  
  //skip data with too little stat
  if (p_Ereco->GetEntries()<nEvtMin) {
    gPad->Clear();
    return 1;
  }
  
  
  TPad *lpad = (TPad*)(mycE->cd());
  FitResult lres;
  if (fitEnergy(p_Ereco,lpad,unit,lres,iSR,useSigmaEff)!=0) return 1;
  lpad->cd();
  char buf[500];
  sprintf(buf,"#gamma p_{T}=%d GeV + PU %d",pT,pu);
  TLatex lat;
  lat.SetTextSize(0.04);
  lat.DrawLatexNDC(0.2,0.96,buf);
  sprintf(buf,"#eta=%3.1f, r = %3.0f mm",etaval,radius);
  lat.SetTextSize(0.04);
  lat.DrawLatexNDC(0.2,0.85,buf);
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
  p_chi2ndf->Fill(lres.chi2/lres.ndf);
  
  //filter out bad points
  if (lres.chi2/lres.ndf > fitQual) {
    std::cout << " --- INFO! Point Egen=" 
	      << pT 
	      << " eta=" << etaval
	      << " pu=" << pu
	      << " skipped, chi2/ndf = "
	      << lres.chi2/lres.ndf
	      << std::endl;
    return 1;
  }

  Int_t np=calibRecoFit->GetN();
  //if (!dovsE) calibRecoFit[iSR]->SetPoint(np,genEn[iE],lres.mean);
  //else 
  calibRecoFit->SetPoint(np,Eval,lres.mean);
  calibRecoFit->SetPointError(np,0.0,lres.meanerr);
  
  //use truth: residual calib won't affect resolution....
  double reso = fabs(lres.sigma/lres.mean);
  //if (pu>0) reso = fabs(lres.sigma/(dovsE?E(genEn[iE],eta[ieta]):genEn[iE]));
  if (pu>0) reso = fabs(lres.sigma/Eval);

  np=resoRecoFit->GetN();
  if (!dovsE) resoRecoFit->SetPoint(np,pT,reso);
  else resoRecoFit->SetPoint(np,Eval,reso);
  double errFit = reso*sqrt(pow(lres.sigmaerr/lres.sigma,2)+pow(lres.meanerr/lres.mean,2));
  resoRecoFit->SetPointError(np,0,errFit);
  
  std::cout << " -- calib/reso value for point " << np << " " << Eval << " " << lres.mean << " +/- " << lres.meanerr 
	    << " sigma " << lres.sigma << " +/- " << lres.sigmaerr
	    << " reso " << reso << " +/- " << errFit
	    << std::endl;


  if (system(TString("mkdir -p ")+plotDir+TString("/Ereco"))) return 1;
  std::ostringstream saveName;
  saveName.str("");
  saveName << plotDir << "/Ereco/Ereco_eta" << eta << "_pu" << pu;
  if (doRaw) saveName << "raw";
  saveName << "_E" << pT << "_SR" << iSR << pSuffix;
  mycE->Update();
  mycE->Print((saveName.str()+".pdf").c_str());
  mycE->Print((saveName.str()+".C").c_str());

  //    drawChi2(p_chi2ndf);
	    
  foutEfit->cd();
  p_chi2ndf->Write();
  p_Ereco->Write();
  
  return 0;
  
};//makeEfit

int plotEGReso(){//main

  SetTdrStyle();
  
  //TString baseDir = "/afs/cern.ch/work/b/bfontana/PFCalEEAna/HGCalTDR/gitv60/";
  //TString saveDir = "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalTDR/gitV08-07-01/";
  TString baseDir = "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalTDR/gitV08-08-00/";
  TString saveDir = "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalTDR/gitV08-08-00/";
  //TString baseDir = "../v60/";
  //TString saveDir = "PlotsEreso/";

  const unsigned nIC = 1;
  const unsigned ICval[nIC] = {3};//0,1,2,3,4,5,10,15,20,50};

  const bool useSigmaEff = true;
  const bool dovsE = true;
  const bool doBackLeakCor = false;
  const bool redoCalib = true;
  const bool redoLeakCor = false;
  const unsigned nBack = 2;

  const unsigned nPu = 1;//4;
  unsigned pu[nPu] = {0};//,140,200};

  const unsigned nS = 1;
  std::string scenario[nS] = {
    //"model2/gamma/",
    "model2/gamma/200u/"
  };
  const unsigned nV = 2;
  TString version[nV] = {"60","68"};

  const unsigned neta = 1;
  unsigned eta[neta]={20};
  //const unsigned neta = 7;
  //unsigned eta[neta]={17,19,21,23,25,27,29};

  //unsigned genEnAll[]={3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};
  unsigned genEnAll[]={5,10,15,20};//,30,40,60,80,100,150,200};
  //unsigned genEnAll[]={5,20,60};
  //unsigned genEnAll[]={7,10,20,30,40};
  const unsigned nGenEnAll=sizeof(genEnAll)/sizeof(unsigned);

  const unsigned nEvtMin = 150;

  const unsigned nSR = 6;
  const unsigned srIdx[nSR] = {0,1,2,3,4,5};
  const double radius[nSR] = {13,15,20,23,26,53};
  //const double radius[nSR] = {25};//{9,14,21,25,37,53};
  double noise100[nSR];
  double noise200[nSR];
  double noise300[nSR];
  unsigned ncellsSmall[nSR] = {7,13,19,31,37,151};
  unsigned ncellsLarge[nSR] = {7,7,13,19,19,85};


  for (unsigned iSR(0); iSR<nSR;++iSR){
    noise100[iSR] = sqrt(pow(sqrt(10*ncellsSmall[iSR])*0.00192,2)+pow(sqrt(10*ncellsSmall[iSR])*0.00241,2)+pow(sqrt(8*ncellsSmall[iSR])*0.00325,2));
    noise200[iSR] = sqrt(pow(sqrt(10*ncellsLarge[iSR])*0.00097,2)+pow(sqrt(10*ncellsLarge[iSR])*0.00121,2)+pow(sqrt(8*ncellsLarge[iSR])*0.00164,2));
    noise300[iSR] = sqrt(pow(sqrt(10*ncellsLarge[iSR])*0.00049,2)+pow(sqrt(10*ncellsLarge[iSR])*0.00062,2)+pow(sqrt(8*ncellsLarge[iSR])*0.00083,2));
    std::cout << " Noise vals for iSR " << iSR << " = " << noise100[iSR] << " " << noise200[iSR] << " " << noise300[iSR] << " GeV." << std::endl;
  }
  //canvas so they are created only once
  TCanvas *mycE[nGenEnAll];
  for (unsigned iE(0); iE<nGenEnAll;++iE){
    std::ostringstream lName;
    lName << "mycE" << genEnAll[iE];
    mycE[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
  }

  TCanvas *mycE2D[nGenEnAll];
  if (redoLeakCor){
    for (unsigned iE(0); iE<nGenEnAll;++iE){
      std::ostringstream lName;
      lName << "mycE2D" << genEnAll[iE];
      mycE2D[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
    }
  }

  for (unsigned ic(0);ic<nIC;++ic){//loop on intercalib
    
    for (unsigned iV(0); iV<nV;++iV){//loop on versions
      //if (iV==nV-1) nLayers = 24;
      const unsigned nLayers = version[iV]=="70" ? 26 : 28;

      for (unsigned iS(0); iS<nS;++iS){//loop on scenarios

	TString inputDir = baseDir+"version"+version[iV]+"/"+TString(scenario[iS])+"/";
	TString plotDir = saveDir+"version"+version[iV]+"/"+TString(scenario[iS])+"/";
	if (system(TString("mkdir -p ")+plotDir)) return 1;


	for (unsigned ieta(0);ieta<neta;++ieta){//loop on eta
	  TTree *ltree[nPu][nGenEnAll];

	  double sigmaStoch[nPu];
	  double sigmaConst[nPu];
	  double sigmaNoise[nPu];
	  double sigmaStochErr[nPu];
	  double sigmaConstErr[nPu];
	  double sigmaNoiseErr[nPu];

	  for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
	    sigmaStoch[ipu] = 0;
	    sigmaConst[ipu] = 0;
	    sigmaNoise[ipu] = 0;
	    sigmaStochErr[ipu] = 0;
	    sigmaConstErr[ipu] = 0;
	    sigmaNoiseErr[ipu] = 0;


	    //loop over energies
	    //check energy point is valid
	    bool skip[nGenEnAll];
	    unsigned nValid = 0;
	    TFile *inputFile[nGenEnAll];
	    for (unsigned iE(0); iE<nGenEnAll; ++iE){
	      skip[iE] = false;
	      inputFile[iE] = 0;
	      ltree[ipu][iE] = 0;
	    }
	    for (unsigned iE(0); iE<nGenEnAll; ++iE){
	      std::ostringstream linputStr;
	      linputStr << inputDir ;
	      //if (eta[ieta]==17) linputStr << "300u/";
	      //else if (eta[ieta]==20) linputStr << "200u/";
	      linputStr << "eta" << eta[ieta] << "_et" << genEnAll[iE];// << "_pu" << pu[ipu];
	      if (pu[ipu]>0) linputStr << "_Pu" << pu[ipu];
	      linputStr << "_IC" << ICval[ic];// << "_Si2";
	      linputStr << ".root";
	      inputFile[iE] = TFile::Open(linputStr.str().c_str());
	      if (!inputFile[iE]) {
		std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Skipping..." << std::endl;
		//	    return 1;
		skip[iE] = true;
	      }
	      else {
		inputFile[iE]->cd("Energies");
		ltree[ipu][iE] = (TTree*)gDirectory->Get("Ereso");
		
		if (!ltree[ipu][iE]){
		  std::cout << " -- File " << inputFile[iE]->GetName() << " sucessfully opened but tree Ereso not found! Skipping." << std::endl;
		  skip[iE] = true;
		} else { 
		  std::cout << " -- File " << inputFile[iE]->GetName() << " sucessfully opened." << std::endl
			    << " ---- Tree found with " << ltree[ipu][iE]->GetEntries() << " entries." << std::endl;
		  if (ltree[ipu][iE]->GetEntries()<nEvtMin) {
		    std::cout << " -- Tree has only " << ltree[ipu][iE]->GetEntries() << " entries, skipping..." << std::endl;
		    skip[iE] = true;
		  } 
		  else nValid++;
		}
	      }
	    }//loop on energies
	    
	    
	    if (nValid < 1) return 1;

	    unsigned newidx = 0;
	    const unsigned nGenEn = nValid;
	    unsigned genEn[nGenEn];
	    unsigned oldIdx[nGenEn];
	    for (unsigned iE(0); iE<nGenEnAll; ++iE){	  
	      if (!skip[iE]) {
		genEn[newidx]=genEnAll[iE];
		oldIdx[newidx]=iE;
		newidx++;
	      }
	    }


	    std::vector<double> noisePu;
	    //loop on SR...
	    //for (unsigned iSR(0); iSR<nSR;++iSR)
	    {//loop on signal region
	      unsigned iSR=srIdx[4];
	      std::cout << " - Processing signal region: " << iSR << " with size " << radius[iSR] << std::endl;


	      double calib = 0;
	      double offset = 0;
	      double calibErr = 0;
	      double offsetErr = 0;
	      
	      for (unsigned iT(redoCalib && pu[ipu]==0 ? 0 : 1); iT<2; ++iT){//loop on calib type
		bool doRaw = iT==0?true:false;
		std::ostringstream loutputStr,linputStr;
		loutputStr << plotDir << "IC" << ICval[ic] << "_pu" << pu[ipu] << "_SR" << iSR << "_Eta" << eta[ieta]; 
		if (dovsE) loutputStr  << "_vsE";
		if (doBackLeakCor) loutputStr  << "_backLeakCor";
		linputStr << loutputStr.str();
		if (doRaw) loutputStr << "_raw";
		loutputStr  << ".root";
		TFile *foutEfit = TFile::Open(loutputStr.str().c_str(),"UPDATE");
		if (!foutEfit){
		  std::cout << " -- Could not create output file " << loutputStr.str() << ". Exit..." << std::endl;
		  return 1;
		}

		TFile *finEfit = 0;
		if (!redoCalib){
		  linputStr << "_raw";
		  linputStr  << ".root";
		  finEfit = TFile::Open(linputStr.str().c_str());
		  if (!finEfit){
		    std::cout << " -- Could not open input file " << linputStr.str() << ". Exit..." << std::endl;
		    return 1;
		  }
		}


		std::cout << " -- Output histograms will be saved in: " << foutEfit->GetName() << std::endl;

		////
		///////// add flag to get graphs from file instead of new
		////
		foutEfit->cd();
		TGraphErrors *corrBackLeakFit = 0;
		if (!doRaw && redoLeakCor) {
		  corrBackLeakFit = new TGraphErrors();
		  corrBackLeakFit->SetName("corrBackLeakFit");
		}
		else if (!doRaw && doBackLeakCor){
		  //get from file
		  foutEfit->cd();
		  corrBackLeakFit = (TGraphErrors*)gDirectory->Get("corrBackLeakFit");
		  if (!corrBackLeakFit) {
		    std::cout << " -- back leak corr TGraph not found! " << std::endl;
		    return 1;
		  }
		}
		/////
		TGraphErrors *calibRecoFit = 0;
		TGraphErrors *calibRecoDelta = 0;
		TGraphErrors *resoRecoFit = 0;
		if (!redoCalib){
		  //get calib and offset from file....
		  std::cout << " -- Reading input calibration histograms from: " << finEfit->GetName() << std::endl;
		  finEfit->cd();
		  calibRecoFit = (TGraphErrors*)gDirectory->Get("calibRecoFitRaw");
		  calibRecoDelta = (TGraphErrors*)gDirectory->Get("calibRecoDeltaRaw");
		  if (!calibRecoFit || !calibRecoDelta){
		    std::cout << " -- Calib TGraphs not found !" << std::endl;
		    return 1;
		  }
		  TF1 *fitFunc = calibRecoFit->GetFunction("calib");
		  calib = fitFunc->GetParameter(1);
		  offset = fitFunc->GetParameter(0);
		  calibErr = fitFunc->GetParError(1);
		  offsetErr = fitFunc->GetParError(0);
		  
		  std::cout << " -- Getting calibration factors from file: " << std::endl
			    << " slope " << calib << " +/- " << calibErr
			    << " offset " << offset << " +/- " << offsetErr
			    << std::endl;
		  calibRecoFit->Delete();
		  calibRecoFit = new TGraphErrors();
		  calibRecoFit->SetName("calibRecoFit");
		  calibRecoDelta->Delete();
		  calibRecoDelta =  new TGraphErrors();
		  calibRecoDelta->SetName("calibRecoDelta");

		}
		else {
		  foutEfit->cd();
		  if (calibRecoFit) calibRecoFit->Delete();
		  calibRecoFit = new TGraphErrors();
		  calibRecoFit->SetName("calibRecoFitRaw");
		  if (calibRecoDelta) calibRecoDelta->Delete();
		  calibRecoDelta =  new TGraphErrors();
		  calibRecoDelta->SetName("calibRecoDeltaRaw");
		}
		
		if (resoRecoFit) resoRecoFit->Delete();
		resoRecoFit = new TGraphErrors();
		if (doRaw) resoRecoFit->SetName("resoRecoFitRaw");
		else resoRecoFit->SetName("resoRecoFit");



		if (doRaw) std::cout << " -- Doing RAW calibration fit " << std::endl;
		else std::cout << " -- Check linearity fit " << std::endl;
		//fit of rawE
		
		
		///////////////////////////////////////////////////////
		////////////// TO BE UPDATED IF PU NEEDED /////////////
		///////////////////////////////////////////////////////
		//fill sigma PU for all SR for 40 GeV E
		if (iT>0 && pu[ipu]>0){// && pT>5 && pT<40) {
		  unsigned iE = 6;
		  double etaval = eta[ieta]/10.;
		  std::vector<double> calib = GetCalib(foutEfit->GetName());
		  bool success = retrievePuSigma(pu[ipu],version[iV],nLayers,
						 ltree[0][oldIdx[iE]], ltree[ipu][oldIdx[iE]],
						 calib,
						 etaval,
						 noisePu);
		  if (!success) return 1;
		}
		///////////////////////////////////////////////////////
		///////////////////////////////////////////////////////
		///////////////////////////////////////////////////////


		for (unsigned iE(0); iE<nGenEn; ++iE){
		  
		  std::cout << " --- Processing energy : " << genEn[iE] 
			    << std::endl;
		  
		  inputFile[iE]->cd("Energies");
		  std::cout << " -- Tree entries for eta=" << eta[ieta] << " pu=" << pu[ipu] << " : " << ltree[ipu][oldIdx[iE]]->GetEntries() << std::endl;

		  double backLeakCor = 0;

		  if (redoLeakCor && !doRaw) {
		    int fail = makeBackLeakCor(nLayers,nBack,
					       iSR,genEn[iE],eta[ieta],pu[ipu],
					       offset,calib,
					       backLeakCor,
					       mycE2D[oldIdx[iE]],
					       ltree[ipu][oldIdx[iE]],
					       foutEfit,
					       corrBackLeakFit,
					       plotDir) ;
		    if (doBackLeakCor && fail){
		      std::cout << " -- Fit failed for back leakage correction, for Egen " << genEn[iE] << std::endl;
		      return 1;
		    }
		  }
		  else {
		    //get it from some saved input file !


		  }
		  
		  int success = makeEfit(useSigmaEff,dovsE,doRaw,
					 doBackLeakCor,nBack,
					 pu[ipu],ICval[ic],
					 scenario[iS],version[iV],nLayers,
					 eta[ieta],genEn[iE],iSR,radius[iSR],
					 offset,calib,backLeakCor,
					 mycE[oldIdx[iE]],
					 ltree[ipu][oldIdx[iE]],
					 foutEfit,
					 plotDir,
					 calibRecoFit,
					 resoRecoFit
					 );
		}//loop on energies

		//get calib
		makeCalibration(doRaw,doBackLeakCor,
				eta[ieta],pu[ipu],iSR,radius[iSR],
				calibRecoFit,calibRecoDelta,
				calib,calibErr,
				offset,offsetErr,
				plotDir,foutEfit
				);

		foutEfit->cd();
		if (redoLeakCor && !doRaw) {
		  corrBackLeakFit->Write();
		  plotBackLeakFit(plotDir,corrBackLeakFit,eta[ieta],pu[ipu]);
		}
		double noise = eta[ieta]==24?noise100[iSR]:eta[ieta]==20?noise200[iSR]:eta[ieta]==17?noise300[iSR]:0 ;
		double noiseRef = pu[ipu]==0? noise : noisePu[iSR];
		
		makeResolution(dovsE,doRaw,doBackLeakCor,
			       eta[ieta],pu[ipu],iSR,radius[iSR],
			       resoRecoFit,
			       sigmaStoch[0],
			       sigmaConst[0],
			       noiseRef,
			       sigmaStoch[ipu],sigmaStochErr[ipu],
			       sigmaConst[ipu],sigmaConstErr[ipu],
			       sigmaNoise[ipu],sigmaNoiseErr[ipu],
			       plotDir,foutEfit
			       );


		foutEfit->Write();
	      }//loop on calib type

	    }//loop on SR
	  }//loop on pu
	}//loop on eta

      }//loop on scenario
    }//loop on version
  }//loop on intercalib


  return 0;

};//main
