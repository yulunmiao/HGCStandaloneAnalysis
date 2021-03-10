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
#include "TGraphErrors.h"
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
    
    valid_args_v_["--versions"] = {"60", "70", "80"};
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
  std::vector<std::string> versions = ip.versions();
  std::vector<std::string> tags = ip.tags();
  const unsigned tag_s = tags.size();
  const unsigned etas_s = etas.size();

  TGraphErrors *div[etas_s*tag_s];

  TCanvas *c[etas_s];
  for(unsigned ieta(0); ieta<etas_s; ++ieta) {
    std::string cname = "RatiosTagSummary" + std::to_string(static_cast<int>(etas[ieta]*10.));
    c[ieta] = new TCanvas(cname.c_str(), cname.c_str(), 800, 600);
    c[ieta]->SetRightMargin(0.02);
    c[ieta]->SetLeftMargin(0.1);
    c[ieta]->SetBottomMargin(0.1);
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
	const std::string hName = "resoRatio" + std::to_string(static_cast<int>(etas[ieta]*10.f)) + "_" + ip.tags()[itag] + "_" + versions[0] + "_" + versions[1];
	std::cout << fName << std::endl;
	div[idx] = (TGraphErrors*)fIn->Get(hName.c_str());
	if(!div[idx]) {
	  std::cout << "File " << hName << " not found." << std::endl;
	  return 1;
	}
      }
    }
  }

  //plot the graphs
  std::unordered_map<std::string, std::string> tmap;
  tmap["V08-08-00"] = "Eff Sigma";
  tmap["V08-08-00-noSigmaEff"] = "Fit Sigma";
  tmap["V08-08-00_20K"] = "Eff Sigma";
  tmap["V08-08-00_20K-noSigmaEff"] = "Fit Sigma";
  std::unordered_map<std::string, std::string> vmap;
  vmap["60"] = "TDR";
  vmap["70"] = "Scenario 13";
  vmap["80"] = "Uniform";

  for(unsigned ieta(0); ieta<etas_s; ++ieta) {
    c[ieta]->cd();

    std::stringstream etavalstr;
    etavalstr << std::fixed << std::setprecision(1) << etas[ieta];

    float ymax=1.02f, ymin=.9f;
    for(unsigned itag(0); itag<tag_s; ++itag) {
      const unsigned idx = itag + ieta*tag_s;
      ymax = std::max( 0.05f + static_cast<float>(div[idx]->GetYaxis()->GetBinUpEdge(div[idx]->GetYaxis()->GetLast())),
		       ymax );
      ymin = std::min( -0.02f + static_cast<float>(div[idx]->GetYaxis()->GetBinLowEdge(div[idx]->GetYaxis()->GetFirst())),
		       ymin );

      div[idx]->GetYaxis()->SetRangeUser(ymin,ymax);
      std::cout << ymin << ", " << ymax << std::endl;      
      div[idx]->SetMarkerStyle(20);
      div[idx]->SetMarkerSize(1.5);
      div[idx]->SetLineColor(itag+2);
      div[idx]->SetMarkerColor(itag+2);
      std::string axistitle = "Ratio: " + vmap[versions[0]] + " / " + vmap[versions[1]];
      div[idx]->GetYaxis()->SetTitle(axistitle.c_str());
      div[idx]->GetXaxis()->SetTitle("E [GeV]");
      div[idx]->GetXaxis()->SetTitleSize(0.04);
      div[idx]->GetXaxis()->SetLabelSize(0.04);
      div[idx]->GetYaxis()->SetTitleOffset(1.1);
      div[idx]->GetYaxis()->SetTitleSize(0.04);
      div[idx]->GetYaxis()->SetLabelSize(0.04);
      div[idx]->GetXaxis()->SetTitleOffset(1.0);

      std::string opt;
      if(itag==0) opt = "ape";
      else opt = "pe same";
      div[idx]->Draw(opt.c_str());
    }

    auto legend = new TLegend(0.82,0.74,0.98,0.9);
    legend->SetTextSize(0.04);
    legend->AddEntry(div[0],tmap[tags[0]].c_str(),"ep");
    legend->AddEntry(div[1],tmap[tags[1]].c_str(),"ep");
    legend->Draw();

    const unsigned idx1 = ieta*tag_s;
    const unsigned idx2 = ieta*tag_s + 1;
    float bol = std::min( div[idx1]->GetXaxis()->GetBinLowEdge(div[idx1]->GetXaxis()->GetFirst()),
			  div[idx2]->GetXaxis()->GetBinLowEdge(div[idx2]->GetXaxis()->GetFirst()) );
    float eol = std::max( div[idx1]->GetXaxis()->GetBinUpEdge(div[idx1]->GetXaxis()->GetLast()),
			  div[idx2]->GetXaxis()->GetBinUpEdge(div[idx2]->GetXaxis()->GetLast()) );
    TLine *line = new TLine(bol,1,eol,1);
    line->SetLineStyle(2);
    line->Draw();

    char buf1[500];
    sprintf(buf1,"#gamma, PU 0");
    TLatex lat1;
    lat1.SetTextSize(0.045);
    lat1.DrawLatexNDC(0.14,0.85,buf1);
    sprintf(buf1,"r = %3.0f mm", radius_map[ip.signalRegions()[0]]);
    lat1.DrawLatexNDC(0.14,0.79,buf1);
    sprintf(buf1,("|#eta| = " + etavalstr.str()).c_str());
    lat1.DrawLatexNDC(0.14,0.74,buf1);
    lat1.DrawLatexNDC(0.12,0.92,"HGCAL G4 standalone");

    std::string cname = dirInBase + "Summary" + std::to_string(static_cast<int>(etas[ieta]*10.)) + "_" + versions[0] + "_" + versions[1] + ".png";
    c[ieta]->SaveAs(cname.c_str());
  }
  return 0;
}

//plots a comparison between two version ratios, each encoded by a different tag
int main(int argc, char** argv)
{
  gStyle->SetOptStat(kFALSE);
  InputParserPlotEGResoRatios ip(argc, argv);

  if(ip.tags().size() == 2 and ip.versions().size() == 2)
    plotRatiosTagSummary(ip);
  else
    std::cout << "Please specify two tags and/or two versions." << std::endl;

  return 0;
}
