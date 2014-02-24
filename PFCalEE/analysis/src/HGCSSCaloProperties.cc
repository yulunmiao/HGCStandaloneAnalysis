#include "HGCSSCaloProperties.hh"
#include "HGCSSSimHit.hh"
#include "RooGaussian.h"
#include "RooPlot.h"
GraphToTF1::GraphToTF1(TString name, TGraph *g){ sp_ = new TSpline3(name,g) ; }
double GraphToTF1::operator()(double *x,double *p) { return sp_->Eval( x[0] ) - p[0]; }

ShowerProfile::ShowerProfile() 
{
  h_rawEn=0;          h_en=0;  h_enFit=0;
  h_enfracVsOverburden=0; h_enVsOverburden=0; h_enVsDistToShowerMax=0;
  gr_raw=0;           gr_centered=0; gr_frac=0;
  gr_unc=0;           gr_relUnc=0;
}

void ShowerProfile::writeTo(TDirectory *dir)
{
  dir->cd();
  if(h_rawEn)               h_rawEn->Clone()->Write();
  if(h_en)                  h_en->Clone()->Write();
  if(h_enFit)               h_enFit->Clone()->Write();
  if(h_enVsOverburden)      h_enVsOverburden->Clone()->Write();
  if(h_enfracVsOverburden)      h_enfracVsOverburden->Clone()->Write();
  if(h_enVsDistToShowerMax) h_enVsDistToShowerMax->Clone()->Write();
  if(gr_raw)                gr_raw->Clone()->Write();
  if(gr_frac)               gr_frac->Clone()->Write();
  if(gr_centered)           gr_centered->Clone()->Write();
  if(gr_unc)                gr_unc->Clone()->Write();
  if(gr_relUnc)             gr_relUnc->Clone()->Write();
}

//
bool ShowerProfile::buildShowerProfile(Float_t eElec, TString version)
{
  TF1 *showerFunc=0;
  
  std::ostringstream indirpath;
  //  indirpath << PATH_TO_G4_FILES << "/"
  //	    << version << "/"
  //	    << PARTICLE_TYPE << "/"
  //	    << "e_" << (Int_t) eElec;
  indirpath << "root://eoscms//eos/cms/store/cmst3/group/hgcal/Geant4/HGcal_" << version << "_e" << (Int_t) eElec << ".root";

  TString inDir(indirpath.str().c_str()); 
  
  //open file
  //TFile *fIn=TFile::Open(inDir+"/PFcal.root");
  TFile *fIn=TFile::Open(inDir);
  if(fIn==0)          return false;
  if(fIn->IsZombie()) return false;
  
  //prepare to run over the ntuple with the information
  TTree *CaloStack=(TTree *)fIn->Get("HGCSSTree");
  if (!CaloStack) return false;

  Float_t event,volNb,volX0,den,denWeight;
  std::vector<HGCSSSimHit> *hitvec = 0;
  CaloStack->SetBranchAddress("event",&event);
  CaloStack->SetBranchAddress("volNb",&volNb);
  CaloStack->SetBranchAddress("volX0trans",&volX0);
  CaloStack->SetBranchAddress("den", &den);
  CaloStack->SetBranchAddress("denWeight", &denWeight);
  CaloStack->SetBranchAddress("HGCSSSimHitVec",&hitvec);
   
  CaloStack->GetEntry(0);
  Float_t curEvent(event);
  
  //energy deposits counter <vol number, < absorber X0, En> >
  std::map<Int_t, std::pair<Float_t,Float_t> > enCounter, enFracCounter;
  Float_t totalRawEn(0), totalEn(0), totalEnFit(0);
  rooRawEn = new RooRealVar("rawEn","rawEn",0,1000000);
  rooEn    = new RooRealVar("En","En",0,1000000);
  data     = new RooDataSet("data","data",RooArgSet(*rooRawEn,*rooEn));

  Float_t mipEn(54.8);
  bool showFit(true);
  float curOverburden(0);
  Float_t refX0(1.0);
  for(Int_t i=0; i<CaloStack->GetEntriesFast(); i++)
    {
      CaloStack->GetEntry(i);
      curOverburden+=volX0;

      if(i==0) refX0=volX0;

      //in case a new event has been found analyze the previous
      if(curEvent!=event){
	curOverburden=0;

	TGraphErrors *showerProf=new TGraphErrors;
	showerProf->SetName("showerprof");
	showerProf->SetMarkerStyle(20);
	Float_t overburden(0);
	for(std::map<Int_t, std::pair<Float_t,Float_t> >::iterator it=enCounter.begin();
	    it!=enCounter.end();
	    it++)
	  {
	    //float ien=it->second.second;
	    float ien=enCounter[it->first].second;

	    overburden += it->second.first;
	    Int_t np=showerProf->GetN();
	    showerProf->SetPoint(np, overburden, ien);
	    showerProf->SetPointError(np,0 ,sqrt(ien) );
	  }
	
	//instantiate fit function if not yet available
	if( showerFunc==0 ) 
	  {
	    showerFunc=new TF1("showerfunc","[0]*pow(x,[1])*exp(-[2]*x)",0,30);//overburden);
	    showerFunc->SetParLimits(0,0.,showerProf->GetMaximum()*1.1);
	    showerFunc->SetParLimits(1,0.,100.);
	    showerFunc->SetParLimits(2,0.,100.);
	    showerFunc->SetLineColor(kBlue);
	 }
	
	//fit for the maximum
	Int_t fitStatus=showerProf->Fit(showerFunc,"RQ+");       
	if(fitStatus==0)
	  {
	    Float_t chi2=showerFunc->GetChisquare();
	    Int_t ndof=showerFunc->GetNDF();
	    Float_t a(showerFunc->GetParameter(1)), b(showerFunc->GetParameter(2));
	    Float_t showerMax(a/b);
	    totalEnFit = showerFunc->Integral(0,28);//overburden);

	    //instantiate shower profile histogram if not yet available
	    if(h_rawEn==0)
	      {
		TString name(""); name+= (Int_t)eElec;
		TString title(""); title+=eElec;  
		h_rawEn=new TH1F("eraw_"+name,title+";Raw energy;Events",400,0,80000);
		h_rawEn->Sumw2();
		h_rawEn->SetDirectory(0);

		h_en=new TH1F("e_"+name,title+";Weighted energy;Events",400,0,80000);
		h_en->SetDirectory(0);

		h_enFit=new TH1F("efit_"+name,title+";Fit energy;Events",400,0,80000);
		h_enFit->SetDirectory(0);

		h_showerMax=new TH1F("smax_"+name,title+";Shower maximum [1/X_{0}];Events",100,0,20);
		h_showerMax->SetDirectory(0);

		h_enVsDistToShowerMax = new TH2F("envsdisttoshowermax_"+name, title+";Distance to shower maximum [1/X_{0}]; Energy; Events", 50,-0.5*overburden,1.5*overburden,100,0,500);
		h_enVsDistToShowerMax->Sumw2();
		h_enVsDistToShowerMax->SetDirectory(0);
		
		h_enVsOverburden = new TH2F("envsoverburden_"+name,title+";Distance transversed [1/X_{0}]; Energy; Events",50,0,overburden,100,0,500);
		h_enVsOverburden->Sumw2();
		h_enVsOverburden->SetDirectory(0);

		h_enfracVsOverburden = new TH2F("enfracvsoverburden_"+name,title+";Distance transversed [1/X_{0}]; Energy fraction; Events",50,0,overburden,100,0,1);
		h_enfracVsOverburden->Sumw2();
		h_enfracVsOverburden->SetDirectory(0);
	      }
	    
	    //energy distributions
	    h_rawEn->Fill( totalRawEn );
	    h_en   ->Fill( totalEn );
	    h_enFit->Fill( totalEnFit );
	    h_showerMax->Fill( showerMax );

	    rooRawEn->setVal(totalRawEn);
	    rooEn->setVal(totalEn);
	    data->add(RooArgSet(*rooRawEn,*rooEn));
	    
	    //shower profiles
	    overburden=0;
	    for(std::map<Int_t, std::pair<Float_t,Float_t> >::iterator it=enCounter.begin();
		it!=enCounter.end();
		it++)
	      {
		overburden += it->second.first;
		Float_t ien=enCounter[it->first].second;
		h_enVsOverburden     ->Fill( overburden,          ien);
		h_enfracVsOverburden     ->Fill( overburden,          enFracCounter[it->first].second);
		h_enVsDistToShowerMax->Fill( overburden-showerMax,ien);
	      }
	    
	    //show fits for 10 events for debug purposes
	    if(showFit){
	      if(event>100)   showFit=false;

	      TCanvas *c=new TCanvas("c","c",500,500);
	      c->SetTopMargin(0.05);
	      c->SetLeftMargin(0.15);
	      c->SetRightMargin(0.05);
	      
	      showerProf->Draw("ap");
	      showerProf->SetMarkerStyle(20);
	      showerProf->GetXaxis()->SetTitle("Transversed thickness [1/X_{0}]");
	      showerProf->GetYaxis()->SetTitle("Energy");
	      showerProf->GetXaxis()->SetLabelSize(0.04);
	      showerProf->GetYaxis()->SetLabelSize(0.04);
	      showerProf->GetXaxis()->SetTitleSize(0.05);
	      showerProf->GetYaxis()->SetTitleSize(0.05);
	      showerProf->GetYaxis()->SetTitleOffset(1.4);
		
	      drawHeader();
	      
	      TPaveText *pt=new TPaveText(0.6,0.6,0.9,0.9,"brNDC");
	      pt->SetBorderSize(0);
	      pt->SetFillStyle(0);
	      pt->SetTextFont(42);
	      pt->SetTextAlign(12);
	      char buf[200];
	      sprintf(buf,"E^{gen}=%3.0f",eElec);
	      pt->AddText(buf);
	      sprintf(buf,"#Sigma E_{i}=%3.0f",totalRawEn);
	      pt->AddText(buf);
	      sprintf(buf,"#Sigma w_{i}E_{i}=%3.0f",totalEn);
	      pt->AddText(buf);
	      sprintf(buf,"Fit E=%3.0f",totalEnFit);
	      pt->AddText(buf);
	      sprintf(buf,"Shower max=%3.1f",showerMax);
	      pt->AddText(buf);
	      sprintf(buf,"#chi^{2}/ndof=%3.0f/%d",chi2,ndof);
	      pt->AddText(buf);
	      pt->Draw();
	   
	      //save 
	      TString name("e_"); name += (Int_t) eElec; name += event;
	      c->SaveAs("PLOTS/"+version+name+"_showerfits.png");

   
	    }
	  }
	   
	//clear energy counters for new event
	enCounter.clear();
	enFracCounter.clear();
	totalRawEn=0;
	totalEn=0;
	totalEnFit=0;
	
	//update event identifier
	curEvent=event;
      }
      
      //account for energy deposit
      //enCounter[ Int_t(volNb) ] = std::pair<Float_t,Float_t>(volX0,den);
      //totalRawEn += den;
      //totalEn    += denWeight;
      
      //account for energy deposit with spatial constraints
      Float_t maxRhoToAcquire(10000);
      //Float_t maxRhoToAcquire(25);
      //Float_t maxRhoToAcquire(10);
      //Float_t maxRhoToAcquire(1.1376*exp(0.1174*volNb));
      //Float_t maxRhoToAcquire(8.705*exp(0.064*volNb));

      Float_t totEnInHits(0);
      for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits 	
	HGCSSSimHit lHit = (*hitvec)[iH];
	lHit.layer(unsigned(volNb));
	double hitEn = lHit.energy()*1e3/mipEn;
	//if(hitEn<1.0) continue;
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	double rho=sqrt(posx*posx+posy*posy);
	if(rho>maxRhoToAcquire) continue;
	totEnInHits += hitEn;
      }

      Float_t weight(volX0/refX0);
      denWeight *= 1.0e3/mipEn;
      den *= 1.0e3/mipEn;
      totalRawEn += totEnInHits;
      if(den>0) totalEn += totEnInHits*weight;
      enCounter[ Int_t(volNb) ]     = std::pair<Float_t,Float_t>(volX0,totEnInHits);
      enFracCounter[ Int_t(volNb) ] = std::pair<Float_t,Float_t>(volX0,den>0 ? totEnInHits/den : 0);
    }
  fIn->Close();
  

  TString name("e_"); name += (Int_t) eElec;

  //draw the centered shower profile
  TCanvas *c=new TCanvas("showerprof","showerprof",500,500);
  c->cd();
  c->SetTopMargin(0.05);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.05);
  drawHeader();
  //h_enVsOverburden->Draw("colz");
  h_enVsDistToShowerMax->Draw("colz");
  c->Modified();
  c->Update();
  c->SaveAs("PLOTS/"+version+name+"_showerprof.png");
  
  //final graphs
  TString title(""); title+=eElec;  
  gr_raw = new TGraphErrors; 
  gr_raw->SetName("rawshowerprof_"+name); 
  gr_raw->SetTitle(title);
  gr_frac     = (TGraphErrors *) gr_raw->Clone("fracshowerprof_"+name);
  gr_centered = (TGraphErrors *) gr_raw->Clone("rawshowerprof_"+name);
  gr_unc      = (TGraphErrors *) gr_raw->Clone("showerprofunc_"+name);
  gr_relUnc   = (TGraphErrors *) gr_raw->Clone("showerprofrelunc_"+name);
  for(Int_t xbin=1; xbin<=h_enVsOverburden->GetXaxis()->GetNbins(); xbin++)
    {
      //centered profile and uncertainties
      TH1D *htemp  = h_enVsDistToShowerMax->ProjectionY("py",xbin,xbin);
      Float_t x    = h_enVsDistToShowerMax->GetXaxis()->GetBinCenter(xbin);
      Float_t xerr = h_enVsDistToShowerMax->GetXaxis()->GetBinWidth(xbin)/2;

      // profile and uncertainties
      //TH1D *htemp  = h_enVsOverburden->ProjectionY("py",xbin,xbin);
      //Float_t x    = h_enVsOverburden->GetXaxis()->GetBinCenter(xbin);
      //Float_t xerr = h_enVsOverburden->GetXaxis()->GetBinWidth(xbin)/2;

      if(htemp->GetMean()>0){
	Int_t np     = gr_centered->GetN();
	gr_centered->SetPoint     (np,x,   htemp->GetMean());
	gr_centered->SetPointError(np,xerr,htemp->GetMeanError());
	gr_unc     ->SetPoint     (np,x,   htemp->GetRMS());
	gr_unc     ->SetPointError(np,xerr,htemp->GetRMSError());
	gr_relUnc  ->SetPoint     (np,x,   htemp->GetRMS()/htemp->GetMean());
	gr_relUnc  ->SetPointError(np,xerr,htemp->GetRMSError()/htemp->GetMean());
      }

     //raw profile
     TH1D *htempraw = h_enVsOverburden->ProjectionY("rawpy",xbin,xbin);
     x              = h_enVsOverburden->GetXaxis()->GetBinCenter(xbin);
     xerr           = h_enVsOverburden->GetXaxis()->GetBinWidth(xbin)/2;
     if(htempraw->GetMean()>0){
       Int_t np     = gr_raw->GetN();
       gr_raw      -> SetPoint     (np,x,   htempraw->GetMean());
       gr_raw      -> SetPointError(np,xerr,htempraw->GetMeanError());
     }

     //energy fraction profile
     TH1D *hfractempraw = h_enfracVsOverburden->ProjectionY("fracpy",xbin,xbin);
     x              = h_enfracVsOverburden->GetXaxis()->GetBinCenter(xbin);
     xerr           = h_enfracVsOverburden->GetXaxis()->GetBinWidth(xbin)/2;
     if(hfractempraw->GetMean()>0){
       Int_t np     = gr_frac->GetN();
       gr_frac      -> SetPoint     (np,x,   hfractempraw->GetMean());
       gr_frac      -> SetPointError(np,xerr,hfractempraw->GetMeanError());
     }

     delete htemp;
     delete htempraw;
     delete hfractempraw;
   }

  //all done here
 return true;
}

CaloProperties::CaloProperties(TString tag) 
{ 
  tag_=tag; gr_showerMax=0; gr_centeredShowerMax=0; 
  genEn_.push_back(5);
  genEn_.push_back(10);
  genEn_.push_back(25);
  genEn_.push_back(50);
  genEn_.push_back(75);
  genEn_.push_back(100);
  genEn_.push_back(200);
  //  genEn_.push_back(300);
  //  genEn_.push_back(500);
}

  
void CaloProperties::writeTo(TDirectoryFile *dir)
{
  if(dir==0) return;
  TDirectory *outDir=dir->mkdir(tag_);
  if(gr_showerMax) gr_showerMax->Clone()->Write();
  if(gr_centeredShowerMax) gr_centeredShowerMax->Clone()->Write();
  for(size_t i=0; i<calibCurve_.size(); i++) calibCurve_[i]->Clone()->Write();
  for(size_t i=0; i<resCurve_.size(); i++)   resCurve_[i]->Clone()->Write();
  for(std::map<Float_t,ShowerProfile>::iterator it=showerProfiles_.begin(); it!=showerProfiles_.end(); it++)
    {
      TString name("e"); name+=Int_t(it->first);
      TDirectory *subDir=outDir->mkdir(name);
      it->second.writeTo(subDir);
    }
}

//
void CaloProperties::characterizeCalo()
{
  gr_showerMax = new TGraph;
  gr_showerMax->SetName("showermax_"+tag_);
  gr_showerMax->SetFillStyle(0);
  gr_showerMax->SetLineColor(kBlue);
  gr_showerMax->SetLineStyle(7);
  gr_showerMax->SetLineWidth(2);
  gr_showerMax->SetMarkerStyle(1);
  gr_showerMax->SetMarkerColor(kBlue);
  gr_centeredShowerMax = (TGraph *) gr_showerMax->Clone("centeredshowermax_"+tag_);

  //build shower profiles
  const Int_t nGenEn=genEn_.size();
  for(Int_t i=0; i<nGenEn; i++)
    {
      Float_t en=genEn_[i];
      ShowerProfile sh;
      sh.buildShowerProfile(en,tag_);
      if(sh.h_rawEn==0) continue;
      showerProfiles_[en]=sh;
      
      GraphToTF1 rawTF1=GraphToTF1( "rawgr2f", sh.gr_raw );
      TF1 *func=new TF1("rawf",rawTF1,sh.gr_raw->GetX()[0],sh.gr_raw->GetX()[ sh.gr_raw->GetN()-1 ],1,"GraphToTF1");
      Double_t xmax=func->GetMaximumX();
      Double_t ymax=func->GetMaximum();
      gr_showerMax->SetPoint(gr_showerMax->GetN(),xmax,ymax);

      GraphToTF1 centeredTF1=GraphToTF1( "centeredgr2f", sh.gr_centered );
      func=new TF1("centeredf",centeredTF1,sh.gr_centered->GetX()[0],sh.gr_centered->GetX()[ sh.gr_centered->GetN()-1 ],1,"GraphToTF1");
      xmax=func->GetMaximumX();
      ymax=func->GetMaximum();
      gr_centeredShowerMax->SetPoint(gr_centeredShowerMax->GetN(),xmax,ymax);
    }
  
  //build the calibration and resolution curves
  for(size_t i=0; i<3; i++)
    {
      TString mode("rawen"),title("raw E");
      if(i==1) { mode="en";    title="weighted E";  }
      if(i==2) { mode="fiten"; title="fitted E";   }
      if(i==3) { mode="showmax"; title="shower max";   }

      if(i!=3){
	calibCurve_.push_back( new TGraphErrors );
	calibCurve_[i]->SetName(mode+"_calib_"+tag_);
	calibCurve_[i]->SetMarkerStyle(20+i);
	calibCurve_[i]->SetTitle(title);
	resCurve_.push_back( (TGraphErrors *) calibCurve_[i]->Clone(mode+"_res_"+tag_) );
      }

      //build the calibration and resolution curves from the energy distributions
      TCanvas *cen=new TCanvas("c"+mode,"c"+mode,500,500);
      TLegend *leg=new TLegend(0.15,0.85,0.9,0.9);
      leg->SetFillStyle(0);
      leg->SetFillColor(0);
      leg->SetBorderSize(0);
      leg->SetTextFont(42);
      Int_t iplot(0);
      for(std::map<Float_t,ShowerProfile>::iterator it=showerProfiles_.begin(); it!=showerProfiles_.end(); it++,iplot++)
	{
	  TH1F *h=it->second.h_rawEn;
	  if(i==1) h=it->second.h_en;
	  if(i==2) h=it->second.h_enFit;
	  if(i==3) h=it->second.h_showerMax;

	  RooRealVar *varToFit=it->second.rooRawEn;
	  if(i==1) varToFit=it->second.rooEn;
	  
	  float unbAvg=it->second.data->mean(*varToFit);
	  
	  varToFit->setRange("fitrange",unbAvg*0.5,unbAvg*2);
	  RooRealVar *unbinnedGauss_mean=new RooRealVar(mode+"mean","mean",unbAvg*0.5,unbAvg*2);
	  RooRealVar *unbinnedGauss_sigma=new RooRealVar(mode+"sigma","sigma",unbAvg*0.001,unbAvg*2);
	  RooGaussian *unbinnedGauss=new RooGaussian(mode+"gauss",title,*varToFit,*unbinnedGauss_mean,*unbinnedGauss_sigma);
	  unbinnedGauss->fitTo( *(it->second.data),RooFit::Range("fitRange") );
	  //it->second.data->printMultiline(std::cout,1000,kTRUE,"\t");

	  // TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",500,500);
	  // RooPlot* frame = varToFit->frame(RooFit::Range("fitRange"));
	  // it->second.data->plotOn(frame);
	  // unbinnedGauss->plotOn(frame);
	  // frame ->Draw();
	  // c->SaveAs(h->GetName()+TString(".png"));

	  Float_t en=it->first;
	  
	  if(i!=3){
	    Int_t np=calibCurve_[i]->GetN();
	    
	    //calibCurve_[i]->SetPoint(np,en,h->GetMean());
	    //calibCurve_[i]->SetPointError(np,0.0,h->GetMeanError());
	    //resCurve_[i]  ->SetPoint(np,1/sqrt(en),h->GetRMS()/h->GetMean());
	    //resCurve_[i]  ->SetPointError(np,0,h->GetRMSError()/h->GetMean());

	    h->Fit("gaus","QR+","same",h->GetMean()-1.5*h->GetRMS(),h->GetMean()+1.5*h->GetRMS());
	    TF1 *gaus=h->GetFunction("gaus");
	    if(gaus){
	      float mean=gaus->GetParameter(1), meanError=gaus->GetParError(1);
	      float sigma=gaus->GetParameter(2), sigmaError=gaus->GetParError(2);
	      calibCurve_[i]->SetPoint(np,en,mean);
	      calibCurve_[i]->SetPointError(np,0.0,meanError);
	      resCurve_[i]  ->SetPoint(np,1/sqrt(en),sigma/mean);
	      resCurve_[i]  ->SetPointError(np,0,sqrt(pow(sigma*meanError,2)+pow(sigmaError*mean,2))/pow(mean,2));
	    }
	    // std::cout << "En=" << en << "GeV mean="<< mean << " sigma=" << sigma << " resol=" << sigma/mean << 
	    //   " -> mean=" << unbinnedGauss_mean->getVal() << " sigma=" << unbinnedGauss_sigma->getVal() << " resol=" << unbinnedGauss_sigma->getVal()/unbinnedGauss_mean->getVal()<< std::endl;
	  }

	  //show it
	  //  if(i!=3) h->Rebin();
	  h->Draw(iplot==0? "hist" : "histsame");
	  leg->AddEntry(h,h->GetTitle(),"f");
	  h->SetLineWidth(1);
	  h->SetLineColor(kCyan-3+iplot);
	  h->SetFillColor(kCyan-3+iplot);
	  h->SetFillStyle(3001);
	  h->SetMarkerStyle(1);
	  h->SetMarkerColor(kCyan-3+iplot);
	  if(iplot>0) continue;
	  h->GetYaxis()->SetTitleOffset(1.6);
	  h->GetYaxis()->SetRangeUser(0,100);
	}
  
      leg->SetNColumns(iplot);
      leg->Draw();
      drawHeader();
      cen->SaveAs("PLOTS/"+tag_+"_c"+mode+".png");
    }
  
  //show the calibration and resolution curves
  for(size_t i=0; i<2; i++)
    {
      TString type(i==0 ? "calib" : "resol");
      std::vector<TGraphErrors *> &grVec=( i==0 ? calibCurve_ : resCurve_);

      TCanvas *c=new TCanvas("c"+type,"c"+type,500,500);
      c->SetTopMargin(0.05);
      c->SetLeftMargin(0.15);
      c->SetRightMargin(0.05);

      TLegend *leg= (i==0 ? 
		     new TLegend(0.15,0.6,0.6,0.9) :
		     new TLegend(0.5,0.2,0.9,0.5) );
      leg->SetFillStyle(0);
      leg->SetFillColor(0);
      leg->SetBorderSize(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.03);
      for(Int_t igr=grVec.size()-1; igr>=0; igr--)
	{
	  if(igr==Int_t(grVec.size()-1)){
	    grVec[igr]->Draw("ap");
	    grVec[igr]->GetYaxis()->SetRangeUser(0,grVec[igr]->GetYaxis()->GetXmax());
	    grVec[igr]->GetYaxis()->SetTitleOffset(1.4);
	    grVec[igr]->GetXaxis()->SetLabelSize(0.04);
	    grVec[igr]->GetYaxis()->SetLabelSize(0.04);
	    grVec[igr]->GetXaxis()->SetTitleSize(0.05);
	    grVec[igr]->GetYaxis()->SetTitleSize(0.05);

	    if(i==0) { grVec[igr]->GetXaxis()->SetTitle("Beam energy [GeV]");                   grVec[igr]->GetYaxis()->SetTitle("Energy estimator"); }
	    else     { grVec[igr]->GetXaxis()->SetTitle("1/#sqrt{Beam energy} [1/#sqrt{GeV}]"); grVec[igr]->GetYaxis()->SetTitle("Relative energy resolution");     }
	  }
	  else{
	    grVec[igr]->Draw("p");
	  }
	  
	  drawHeader();

	  char buf[500];
	  if(i==0) {
	    TString name("calibfunc"); name += grVec[igr]->GetTitle(); name.ReplaceAll(" ","");
	    TF1 *fitFunc=new TF1(name,"[0]+[1]*x",grVec[igr]->GetXaxis()->GetXmin(),grVec[igr]->GetXaxis()->GetXmax());
	    grVec[igr]->Fit(fitFunc,"RME");
	    sprintf(buf,"#splitline{%s}{<E> #propto %3.3f+%3.3f#timesE}",grVec[igr]->GetTitle(),fitFunc->GetParameter(0),fitFunc->GetParameter(1));
	  }else{
	    TString name("resfunc"); name += grVec[igr]->GetTitle(); name.ReplaceAll(" ","");
	    //TF1 *fitFunc=new TF1(name,"sqrt([0]*x*x+[1]+[2]*x*x*x*x)",grVec[igr]->GetXaxis()->GetXmin(),grVec[igr]->GetXaxis()->GetXmax());
	    TF1 *fitFunc=new TF1(name,"sqrt([0]*x*x+[1])",grVec[igr]->GetXaxis()->GetXmin(),grVec[igr]->GetXaxis()->GetXmax());
	    fitFunc->SetParameter(0,0.2);
	    fitFunc->SetParLimits(0,0,1);
	    fitFunc->SetParameter(1,0);
	    fitFunc->SetParLimits(1,0,0.3);
	    //fitFunc->SetParameter(2,0);
	    
	    grVec[igr]->Fit(fitFunc,"RME");
	    Float_t sigmaStoch(sqrt(fitFunc->GetParameter(0)));
	    Float_t sigmaStochErr(fitFunc->GetParError(0)/(2*sigmaStoch));
	    Float_t sigmaConst(sqrt(fitFunc->GetParameter(1)));
	    Float_t sigmaConstErr(fitFunc->GetParError(1)/(2*sigmaConst));
	    //Float_t sigmaNoise(sqrt(fitFunc->GetParameter(2)));
	    //Float_t sigmaNoiseErr(sqrt(fitFunc->GetParameter(2))/(2*sigmaNoise));
	    //sprintf(buf,"#splitline{%s}{#frac{#sigma}{E} #propto #frac{%3.3f}{#sqrt{E}} #oplus %3.3f #oplus #frac{%3.3f}{E}}",grVec[igr]->GetTitle(),sigmaStoch,sigmaConst,sigmaNoise);
	    sprintf(buf,"#splitline{%s}{#frac{#sigma}{E} #propto #frac{%3.3f}{#sqrt{E}} #oplus %3.3f}",grVec[igr]->GetTitle(),sigmaStoch,sigmaConst);
	    stochTerms_.push_back( Measurement_t(sigmaStoch,sigmaStochErr) );
	    constTerms_.push_back( Measurement_t(sigmaConst,sigmaConstErr) );
	  }
	  leg->AddEntry(grVec[igr],buf,"lp");
       }
      leg->Draw();

      c->Modified();
      c->Update();
      c->SaveAs("PLOTS/"+tag_+"_"+type+".png");
      c->SaveAs("PLOTS/"+tag_+"_"+type+".root");
    }
  
  //shower profile
  TCanvas *craw=new TCanvas("rawprof","rawprof",500,500);
  TCanvas *cfrac=new TCanvas("fracprof","rawprof",500,500);
  TCanvas *ccen=new TCanvas("cenprof","cenprof",500,500);
  TCanvas *cunc=new TCanvas("uncprof","uncprof",500,500);
  TCanvas *crelunc=new TCanvas("reluncprof","reluncprof",500,500);
  TLegend *leg=new TLegend(0.2,0.85,0.9,0.94);

  Int_t igr(0);
  for(std::map<Float_t,ShowerProfile>::reverse_iterator it = showerProfiles_.rbegin(); 
      it!=showerProfiles_.rend();
      it++,igr++){
    
    craw->cd();
    it->second.gr_raw->Draw(igr==0? "ae2p" : "e2p");
    it->second.gr_raw->SetLineWidth(2);
    it->second.gr_raw->SetLineColor(igr%2+1);
    it->second.gr_raw->SetFillColor(0);
    it->second.gr_raw->SetFillStyle(0);
    it->second.gr_raw->SetMarkerStyle(20+igr);
    it->second.gr_raw->SetMarkerColor(igr%2+1);
    it->second.gr_raw->GetXaxis()->SetTitle("Transversed thickness [1/X_{0}]");
    it->second.gr_raw->GetYaxis()->SetTitle("Energy");
    // it->second.gr_raw->GetYaxis()->SetRangeUser(0,it->second.gr_raw->GetMaximum()*1.1);
    leg->AddEntry(it->second.gr_raw,it->second.gr_raw->GetTitle(),"p");

    cfrac->cd();
    it->second.gr_frac->Draw(igr==0? "ae2p" : "e2p");
    it->second.gr_frac->SetLineWidth(2);
    it->second.gr_frac->SetLineColor(igr%2+1);
    it->second.gr_frac->SetFillColor(0);
    it->second.gr_frac->SetFillStyle(0);
    it->second.gr_frac->SetMarkerStyle(20+igr);
    it->second.gr_frac->SetMarkerColor(igr%2+1);
    it->second.gr_frac->GetXaxis()->SetTitle("Transversed thickness [1/X_{0}]");
    it->second.gr_frac->GetYaxis()->SetTitle("<Energy fraction>");
    // it->second.gr_raw->GetYaxis()->SetRangeUser(0,it->second.gr_raw->GetMaximum()*1.1);

    ccen->cd();
    it->second.gr_centered->Draw(igr==0? "ap" : "p");
    it->second.gr_centered->SetLineWidth(2);
    it->second.gr_centered->SetLineColor(igr%2+1);
    it->second.gr_centered->SetFillColor(0);
    it->second.gr_centered->SetFillStyle(0);
    it->second.gr_centered->SetMarkerStyle(20+igr);
    it->second.gr_centered->SetMarkerColor(igr%2+1);
    //  it->second.gr_centered->GetYaxis()->SetRangeUser(0,it->second.gr_centered->GetMaximum()*1.1);
    it->second.gr_centered->GetXaxis()->SetTitle("<Transversed thickness - shower max> [1/X_{0}]");
    it->second.gr_centered->GetYaxis()->SetTitle("Energy");
    
    cunc->cd();
    it->second.gr_unc->Draw(igr==0? "ap" : "p");
    it->second.gr_unc->SetLineWidth(2);
    it->second.gr_unc->SetLineColor(igr%2+1);
    it->second.gr_unc->SetFillColor(0);
    it->second.gr_unc->SetFillStyle(0);
    it->second.gr_unc->SetMarkerStyle(20+igr);
    it->second.gr_unc->SetMarkerColor(igr%2+1);
    //  it->second.gr_unc->GetYaxis()->SetRangeUser(0,it->second.gr_unc->GetMaximum()*1.1);
    it->second.gr_unc->GetXaxis()->SetTitle("<Transversed thickness - shower max> [1/X_{0}]");
    it->second.gr_unc->GetYaxis()->SetTitle("RMS(Energy)");

    crelunc->cd();
    it->second.gr_relUnc->Draw(igr==0? "ap" : "p");
    it->second.gr_relUnc->SetLineWidth(2);
    it->second.gr_relUnc->SetLineColor(igr%2+1);
    it->second.gr_relUnc->SetFillColor(0);
    it->second.gr_relUnc->SetFillStyle(0);
    it->second.gr_relUnc->SetMarkerStyle(20+igr);
    it->second.gr_relUnc->SetMarkerColor(igr%2+1);
    // it->second.gr_relUnc->GetYaxis()->SetRangeUser(0,it->second.gr_relUnc->GetMaximum()*1.1);
    it->second.gr_relUnc->GetXaxis()->SetTitle("<Transversed thickness - shower max> [1/X_{0}]");
    it->second.gr_relUnc->GetYaxis()->SetTitle("RMS/<Energy>");
  }
  
  craw->cd();
  gr_showerMax->Draw("l");
  drawHeader();
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetNColumns(5);
  leg->Draw();
  craw->Modified();
  craw->Update();
  craw->SaveAs("PLOTS/"+tag_+"_rawprof.png");

  cfrac->cd();
  drawHeader();
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetNColumns(5);
  leg->Draw();
  cfrac->Modified();
  cfrac->Update();
  cfrac->SaveAs("PLOTS/"+tag_+"_fracprof.png");

  ccen->cd();
  gr_centeredShowerMax->Draw("l");
  drawHeader();
  leg->Draw();
  ccen->Modified();
  ccen->Update();
  ccen->SaveAs("PLOTS/"+tag_+"_cenprof.png");

  cunc->cd();
  drawHeader();
  leg->Draw();
  cunc->Modified();
  cunc->Update();
  cunc->SaveAs("PLOTS/"+tag_+"_cenunc.png");

  crelunc->cd();
  drawHeader();
  leg->Draw();
  crelunc->Modified();
  crelunc->Update();
  crelunc->SaveAs("PLOTS/"+tag_+"_crelunc.png");

  /*
  TFile *fOut=TFile::Open("CaloPerformance_"+tag_+".root","RECREATE");
  fOut->cd();
  for(Int_t i=0; i<erawToSave.GetEntriesFast(); i++) {
    erawToSave.At(i)->Clone()->Write();
    profToSave.At(i)->Clone()->Write();
    profUncToSave.At(i)->Clone()->Write();
  }
  for(Int_t i=0; i<finalToSave.GetEntriesFast(); i++){
    finalToSave.At(i)->Clone()->Write();
  }
  fOut->Close();
*/
  //return props;
}
//
void drawHeader()
{
  TPaveText *pt=new TPaveText(0.15,0.95,0.6,0.99,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->AddText("Geant4 simulation");
  pt->Draw();
}

//
void setStyle()
{
  //configure results style
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.1);
}

