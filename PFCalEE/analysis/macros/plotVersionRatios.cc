#include<iostream>
#include<iomanip>
#include<sstream>
#include<string>
#include<vector>

#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#include "../include/InputParserBase.h"

class InputParserPlotEGResoRatios: public InputParserBase {
public:
  explicit InputParserPlotEGResoRatios(int argc, char** argv): InputParserBase(argc, argv) { run(); }

  std::string tag() const { return tag_; }
  std::vector<float> etas() const { return etas_; }
  std::vector<std::string> versions() const { return versions_; }
  std::vector<unsigned> signalRegions() const { return signalRegions_; }

private:
  std::string tag_;
  std::vector<float> etas_;
  std::vector<std::string> versions_;
  std::vector<unsigned> signalRegions_;

  void set_args_final_()
  {
    tag_ = chosen_args_["--tag"];
    for(auto&& x: chosen_args_v_["--etas"])
      etas_.push_back( std::stof(x) );
    for(auto&& x: chosen_args_v_["--signalRegions"])
      signalRegions_.push_back( static_cast<unsigned>(std::stoi(x)) );
    versions_ = chosen_args_v_["--versions"];
  }
  
  void set_args_options_()
  {
    required_args_ = { "--versions", "--tag", "--etas", "--signalRegions"  };

    valid_args_v_["--versions"] = {"60", "70"};
    valid_args_v_["--signalRegions"] = {"0", "1", "2", "3", "4", "5"};
    free_args_ = {"--tag"};
    free_args_v_ = {"--etas"};
    optional_args_ = {""};
  }
};

int plotVersionRatios(const InputParserPlotEGResoRatios& ip) {
  std::vector<float> etas = ip.etas();
  std::vector<std::string> versions = ip.versions();
  const unsigned etas_s = etas.size();
  const unsigned versions_s = versions.size();

  std::unordered_map<unsigned,float> radius_map;
  radius_map[4] = 26.0;

  TGraphErrors *resoV[etas_s*versions_s];

  std::string dirInBase = "/eos/user/b/bfontana/www/RemoveLayers/" + ip.tag() + "/";

  for(unsigned ieta(0); ieta<etas_s; ++ieta) {
    
    for(unsigned iversion(0); iversion<versions_s; ++iversion) {
      const unsigned idx = iversion + ieta*versions_s;

      std::string dirIn = ( dirInBase + "version" + versions[iversion] + "/model2/gamma/SR4/");
      std::string fileIn = ( dirIn + "IC3_pu0_SR4_Eta" + std::to_string(static_cast<int>(etas[ieta]*10.f))
			     + "_vsE_backLeakCor_raw.root" );

      TFile *fIn = TFile::Open(fileIn.c_str(),"READ");
      if(fIn)
	fIn->cd();
      else {
	std::cout << "Problem reading the file."  << std::endl;
	return 1;
      }
      resoV[idx] = (TGraphErrors*)gDirectory->Get("resoRecoFitRaw");    
    }
  }

  std::unordered_map<std::string, std::string> vmap;
  vmap["60"] = "TDR (no neutron moderator)";
  vmap["70"] = "Scenario 13";

  for (unsigned ieta(0); ieta<etas_s; ++ieta)
    {
      const unsigned idx1 = ieta*versions_s;
      const unsigned idx2 = idx1 + 1;

      int npoints = resoV[idx1]->GetN();      
      int npoints2 = resoV[idx2]->GetN();
      int maxpoints = std::max(npoints,npoints2);
      if(npoints!=npoints2)
	std::cout << "WARNING --- The graphs have a different number of data points: " << npoints << " vs. " << npoints2 << "!" << std::endl;
      std::vector<float> x1v, x2v, y1v, y2v, e1v, e2v;
				 
      for (int j(0), k(0); j<maxpoints and k<maxpoints; j++, k++) {
	if(j>=npoints or k>= npoints2)
	  continue;
	
	double x1, x2, y1, y2;

	bool has_different_abciss = true;
	while(has_different_abciss) {
	    resoV[idx1]->GetPoint(j,x1,y1);
	    resoV[idx2]->GetPoint(k,x2,y2);
	    if(static_cast<int>(x1)!=static_cast<int>(x2)) {
	      if(x1 > x2)
		k++;
	      else
		j++;
	    }
	    else
	      has_different_abciss = false;
	  }

	x1v.push_back(x1);
	x2v.push_back(x2);
	y1v.push_back(y1);
	y2v.push_back(y2);
	e1v.push_back(resoV[idx1]->GetErrorY(j));
	e2v.push_back(resoV[idx2]->GetErrorY(j));
      }
      assert(x1v == x2v);
      assert(x1v.size() == y1v.size());
      
      double xbins[x1v.size()+1];
      for(unsigned q(0); q<x1v.size(); ++q) {
	if(q==0)
	  xbins[q] = x1v[q]-5.;
	else
	  xbins[q] = (x1v[q-1]+x1v[q])/2.;
      }
      xbins[x1v.size()] = x1v.back()+5.;

      std::string name = "cRatio" + std::to_string(ieta);
      std::string title = "resoRatio" + std::to_string(static_cast<int>(etas[ieta]*10.f));
      TCanvas *c = new TCanvas(name.c_str(), title.c_str(), 1000, 600);
      c->SetRightMargin(0.13);
      c->SetLeftMargin(0.05);
      c->SetBottomMargin(0.50);
      c->Draw();
      TPad *p1 = new TPad("p1","p1", 0.1, 0.33, 0.9, 1.);
      p1->SetBottomMargin(0.);
      p1->Draw();
      p1->cd();

      std::string hname = "h" + std::to_string(static_cast<int>(etas[ieta]*10.));
      std::string htitle = "hReso" + std::to_string(static_cast<int>(etas[ieta]*10.));
      TH1F *h1 = new TH1F(hname.c_str(), htitle.c_str(), x1v.size(), xbins);
      TH1F *h2 = new TH1F(hname.c_str(), htitle.c_str(), x1v.size(), xbins);
      for(unsigned q(0); q<x1v.size(); ++q) {
	h1->SetBinContent(q+1,y1v[q]);
	h2->SetBinContent(q+1,y2v[q]);
	h1->SetBinError(q+1,e1v[q]);
	h2->SetBinError(q+1,e2v[q]);
      }
      h1->SetTitle("");
      h2->SetTitle("");
      h1->SetMarkerStyle(1);
      h1->SetMarkerColor(1);
      h1->SetLineColor(1);
      //h1->GetXaxis()->SetRangeUser(-5.,295.);
      //h1->GetYaxis()->SetRangeUser(0.001,0.16);
      h1->GetXaxis()->SetLabelSize(0.05);
      h1->GetYaxis()->SetLabelSize(0.05);
      h1->GetYaxis()->SetTitleSize(0.05);
      h1->GetYaxis()->SetTitle("#sigma/E");
      h1->SetLineColor(kRed);
      h1->Draw();
      h2->SetLineColor(kBlue);
      h2->Draw("SAME");

      std::stringstream etavalstr;
      etavalstr << std::fixed << std::setprecision(1) << etas[ieta];

      char buf1[500];
      sprintf(buf1,"#gamma, PU 0");
      TLatex lat1;
      lat1.SetTextSize(0.06);
      lat1.DrawLatexNDC(0.20,0.8,buf1);
      sprintf(buf1,"r = %3.0f mm", radius_map[ip.signalRegions()[0]]);
      lat1.DrawLatexNDC(0.20,0.74,buf1);
      sprintf(buf1,("|#eta| = " + etavalstr.str()).c_str());
      lat1.DrawLatexNDC(0.20,0.69,buf1);

      lat1.DrawLatexNDC(0.11,0.92,"HGCAL G4 standalone");
      
      auto legend = new TLegend(0.5,0.7,0.9,0.9);
      legend->SetTextSize(0.05);
      legend->AddEntry(h1,vmap[versions[0]].c_str(),"f");
      legend->AddEntry(h2,vmap[versions[1]].c_str(),"f");
      legend->Draw();

      c->cd(0);
      TPad *p2 = new TPad("p2","p2",0.1,0.,0.9,0.33);
      p2->SetTopMargin(0.0);
      p2->Draw();
      p2->cd();
      
      std::string hdivname = "hdiv" + std::to_string(static_cast<int>(etas[ieta]*10.));
      std::string hdivtitle = "hResodiv" + std::to_string(static_cast<int>(etas[ieta]*10.));
      TH1F *div = new TH1F(hdivname.c_str(), hdivtitle.c_str(), x1v.size(), xbins);
      div->Sumw2();
      div->Divide(h1,h2,1.,1.,"b");
      div->SetAxisRange(0.8, 1.1,"Y");
      div->SetTitle("");
      div->GetXaxis()->SetTitle("E (GeV)");
      div->GetXaxis()->SetTitleSize(0.13);
      div->GetXaxis()->SetLabelSize(0.13);
      div->GetYaxis()->SetLabelSize(0.1);
      div->Draw();
      TLine *line = new TLine(x1v[0],1,x1v.back(),1);
      line->SetLineStyle(2);
      line->Draw("same");

      double factor = TMath::Sqrt(26./28.);
      TLine *line2 = new TLine(x1v[0],factor,x1v.back(),factor);
      line2->SetLineColor(kMagenta);
      line2->SetLineStyle(5);
      line2->Draw("same");

      c->SaveAs((dirInBase + title + ".png").c_str());
    }

  return 0;
}

int main(int argc, char** argv)
{
  gStyle->SetOptStat(kFALSE);
  InputParserPlotEGResoRatios ip(argc, argv);

  if(ip.signalRegions().size() != 1)
    std::cout << "This code is not ready to take more than one signal region. " <<  std::endl;

  if(ip.versions().size() == 2)
    plotVersionRatios(ip);
  else
    std::cout << "Please specify the two versions." << std::endl;
  return 0;
}
