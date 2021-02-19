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
#include "TH1F.h"
#include "TLatex.h"

#include "../include/InputParserBase.h"

class InputParserPlotEGResoRatios: public InputParserBase {
public:
  explicit InputParserPlotEGResoRatios(int argc, char** argv): InputParserBase(argc, argv) { run(); }

  std::vector<std::string> tags() const { return tags_; }
  std::vector<float> etas() const { return etas_; }
  std::vector<std::string> versions() const { return versions_; }
  std::vector<unsigned> signalRegions() const { return signalRegions_; }

private:
  std::vector<std::string> tags_;
  std::vector<float> etas_;
  std::vector<std::string> versions_;
  std::vector<unsigned> signalRegions_;

  void set_args_final_()
  {
    for(auto&& x: chosen_args_v_["--tags"])
      tags_.push_back(x);
    for(auto&& x: chosen_args_v_["--etas"])
      etas_.push_back( std::stof(x) );
    for(auto&& x: chosen_args_v_["--signalRegions"])
      signalRegions_.push_back( static_cast<unsigned>(std::stoi(x)) );
    versions_ = chosen_args_v_["--versions"];
  }
  
  void set_args_options_()
  {
    required_args_ = { "--versions", "--tags", "--etas", "--signalRegions"  };
    optional_args_ = {""};
    
    valid_args_v_["--versions"] = {"60", "70"};
    valid_args_v_["--signalRegions"] = {"0", "1", "2", "3", "4", "5"};
    free_args_v_ = {"--etas", "--tags"};
  }
};
int plotRatiosTagSummary(const InputParserPlotEGResoRatios& ip) {
  std::unordered_map<unsigned,float> radius_map;
  radius_map[4] = 26.0;
  
  const std::string dirInBase = "/eos/user/b/bfontana/www/RemoveLayers/";
  const std::string fName = dirInBase + "fOutReso.root";
  std::vector<float> etas = ip.etas();
  std::vector<std::string> tags = ip.tags();
  const unsigned tag_s = tags.size();
  const unsigned etas_s = etas.size();

  TH1F* resoV[etas_s*tag_s];

  TCanvas *c[etas_s];
  for(unsigned ieta(0); ieta<etas_s; ++ieta) {
    std::string cname = "RatiosTagSummary" + std::to_string(static_cast<int>(etas[ieta]*10.));
    c[ieta] = new TCanvas(cname.c_str(), cname.c_str(), 800, 600);
    c[ieta]->SetRightMargin(0.20);
    c[ieta]->SetLeftMargin(0.10);
    c[ieta]->SetBottomMargin(0.10);
    c[ieta]->Draw();
  }

  //access the graphs
  TFile *fIn = TFile::Open(fName.c_str(), "READ");
  if(fIn)
    fIn->cd();
  else {
    std::cout << "Problem reading the file."  << std::endl;
    return 1;
  }

  for(unsigned ieta(0); ieta<etas_s; ++ieta) {
    for(unsigned itag(0); itag<tag_s; ++itag) {    
      {
	const unsigned idx = itag + ieta*tag_s;
	const std::string hName = "hdiv" + std::to_string(static_cast<int>(etas[ieta]*10.f)) + "_" + tags[itag];
	resoV[idx] = (TH1F*)fIn->Get(hName.c_str());
	if(!resoV[idx]) {
	  std::cout << "File " << hName << " not found." << std::endl;
	  return 1;
	}
      }
    }
  }

  //plot the graphs
  std::unordered_map<std::string, std::string> tmap;
  tmap["V08-08-00"] = "Eff Sigma";
  tmap["V08-08-00-noSigmaEff"] = "Fit sigma";

  for(unsigned ieta(0); ieta<etas_s; ++ieta) {
    c[ieta]->cd();
    auto legend = new TLegend(0.68,0.74,0.9,0.9);
    legend->SetTextSize(0.04);

    std::stringstream etavalstr;
    etavalstr << std::fixed << std::setprecision(1) << etas[ieta];

    for(unsigned itag(0); itag<tag_s; ++itag) {
      const unsigned idx = itag + ieta*tag_s;
      resoV[idx]->GetYaxis()->SetRangeUser(.8,1.15);
      //resoV[idx]->GetXaxis()->SetRangeUser(0.,300);
      resoV[idx]->GetXaxis()->SetLabelSize(0.04);
      resoV[idx]->GetYaxis()->SetLabelSize(0.04);
      resoV[idx]->GetXaxis()->SetTitleSize(0.04);
      resoV[idx]->GetYaxis()->SetTitleSize(0.04);
      resoV[idx]->GetYaxis()->SetTitle("#sigma/E ( TDR / Scenario 13 )");
      resoV[idx]->SetLineColor(itag+2);
      std::string opt = itag==0 ? "" : "same";
      resoV[idx]->Draw(opt.c_str());
      legend->AddEntry(resoV[idx],tmap[tags[itag]].c_str(),"f");
    }
    
    legend->Draw();

    char buf1[500];
    sprintf(buf1,"#gamma, PU 0");
    TLatex lat1;
    lat1.SetTextSize(0.06);
    lat1.DrawLatexNDC(0.12,0.85,buf1);
    sprintf(buf1,"r = %3.0f mm", radius_map[ip.signalRegions()[0]]);
    lat1.DrawLatexNDC(0.12,0.79,buf1);
    sprintf(buf1,("|#eta| = " + etavalstr.str()).c_str());
    lat1.DrawLatexNDC(0.12,0.74,buf1);
    lat1.DrawLatexNDC(0.11,0.92,"HGCAL G4 standalone");

    std::string cname = dirInBase + "RatiosTagSummary" + std::to_string(static_cast<int>(etas[ieta]*10.)) + ".png";
    c[ieta]->SaveAs(cname.c_str());
  }
  return 0;
}

int main(int argc, char** argv)
{
  gStyle->SetOptStat(kFALSE);
  InputParserPlotEGResoRatios ip(argc, argv);

  if(ip.tags().size() == 2)
    plotRatiosTagSummary(ip);
  else
    std::cout << "Please specify two tags." << std::endl;

  return 0;
}
