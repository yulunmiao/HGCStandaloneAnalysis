#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
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

    valid_args_v_["--versions"] = {"60", "70", "80"};
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
  std::string dirInBase1 = "/eos/user/b/bfontana/www/RemoveLayers/";
  std::string dirInBase2 = dirInBase1 + ip.tag() + "/";

  std::string fname = dirInBase1 + "fOutReso.root";
  TFile *foutReso = TFile::Open(fname.c_str(),"UPDATE");

  for(unsigned ieta(0); ieta<etas_s; ++ieta) {
    
    for(unsigned iversion(0); iversion<versions_s; ++iversion) {
      const unsigned idx = iversion + ieta*versions_s;

      std::string dirIn = ( dirInBase2 + "version" + versions[iversion] + "/model2/gamma/SR4/");
      std::string fileIn = ( dirIn + "IC3_pu0_SR4_Eta" + std::to_string(static_cast<int>(etas[ieta]*10.f))
			     + "_vsE_backLeakCor.root" );

      TFile *fIn = TFile::Open(fileIn.c_str(),"READ");
      if(fIn)
	fIn->cd();
      else {
	std::cout << "Problem reading the file."  << std::endl;
	return 1;
      }
      resoV[idx] = (TGraphErrors*)gDirectory->Get("resoRecoFit");    
    }
  }

  std::unordered_map<std::string, std::string> vmap;
  vmap["60"] = "TDR";
  vmap["70"] = "Scenario 13";
  vmap["80"] = "Uniform";

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
      
      std::string name = "cRatio" + std::to_string(ieta) + "_" + ip.tag();
      std::string title = "resoRatio" + std::to_string(static_cast<int>(etas[ieta]*10.f)) + "_" + ip.tag() + "_" + versions[0] + "_" + versions[1];
      TCanvas *c = new TCanvas(name.c_str(), title.c_str(), 800, 600);
      c->SetRightMargin(0.1);
      c->SetLeftMargin(0.1);
      c->SetBottomMargin(-0.3);
      c->Draw();
      TPad *p1 = new TPad("p1","p1", 0.01, 0.45, 0.99, 1.);
      p1->SetBottomMargin(0.);
      p1->Draw();
      p1->cd();

      resoV[idx1]->SetTitle("");
      resoV[idx2]->SetTitle("");
      resoV[idx1]->GetYaxis()->SetRangeUser(0.,0.09);
      resoV[idx1]->GetXaxis()->SetLabelSize(0.00);
      resoV[idx1]->GetYaxis()->SetLabelSize(0.06);
      resoV[idx1]->GetYaxis()->SetTitle("#sigma/E");
      resoV[idx1]->GetYaxis()->SetTitleSize(0.07);
      resoV[idx1]->GetYaxis()->SetTitleOffset(0.65);
      resoV[idx1]->SetMarkerStyle(20);
      resoV[idx1]->SetMarkerSize(1.5);
      resoV[idx1]->SetMarkerColor(kRed);
      resoV[idx1]->GetFunction("reso")->SetLineColor(kRed);
      resoV[idx1]->Draw("APE");
      resoV[idx2]->SetMarkerStyle(22);
      resoV[idx2]->SetMarkerSize(1.5);
      resoV[idx2]->SetMarkerColor(kBlue);
      resoV[idx2]->GetFunction("reso")->SetLineColor(kBlue);
      resoV[idx2]->Draw("PE SAME");

      std::stringstream etavalstr;
      etavalstr << std::fixed << std::setprecision(1) << etas[ieta];

      char buf1[500];
      sprintf(buf1,"#gamma, PU 0");
      TLatex lat1;
      lat1.SetTextSize(0.06);
      lat1.DrawLatexNDC(0.13,0.84,buf1);
      sprintf(buf1,"r = %3.0f mm", radius_map[ip.signalRegions()[0]]);
      lat1.DrawLatexNDC(0.13,0.77,buf1);
      sprintf(buf1,("|#eta| = " + etavalstr.str()).c_str());
      lat1.DrawLatexNDC(0.13,0.71,buf1);
      lat1.SetTextSize(0.07);
      lat1.DrawLatexNDC(0.11,0.92,"HGCAL G4 standalone");
      
      c->cd(0);
      TPad *p2 = new TPad("p2","p2",0.01,0.1,0.99,0.45);
      p2->SetTopMargin(0.0);
      p2->Draw();
      p2->cd();

      float xvals[x1v.size()], yvals[y1v.size()];
      float exvals[x1v.size()], eyvals[y1v.size()];
      assert(y1v.size() == y2v.size());
      assert(y1v.size() == x1v.size());
      for(unsigned i(0); i<x1v.size(); ++i) {
	xvals[i] = x1v[i];
	exvals[i] = 0;
	yvals[i] = y1v[i] / y2v[i];
	//error propagation of division
	eyvals[i] = (y1v[i]*y1v[i])/(y2v[i]*y2v[i]) * std::sqrt( (e1v[i]*e1v[i])/(y1v[i]*y1v[i]) + (e2v[i]*e2v[i])/(y2v[i]*y2v[i]) );
      }

      TGraphErrors *div = new TGraphErrors(x1v.size(), xvals, yvals, exvals, eyvals);
      div->SetTitle();
      div->SetName(title.c_str());
      div->GetXaxis()->SetTitle("E [GeV] ");
      div->GetXaxis()->SetTitleSize(0.11);
      div->GetXaxis()->SetLabelSize(0.1);
      div->GetXaxis()->SetTitleOffset(-0.4);
      div->GetYaxis()->SetLabelSize(0.08);
      div->GetYaxis()->SetRangeUser(0.93,1.015);
      div->SetMarkerColor(kGreen);
      div->SetMarkerStyle(20);
      div->SetMarkerSize(1.3);
      div->Draw("ap");
      
      float bol = div->GetXaxis()->GetBinLowEdge(div->GetXaxis()->GetFirst());
      float eol = div->GetXaxis()->GetBinUpEdge(div->GetXaxis()->GetLast());

      TLine *line = new TLine(bol,1,eol,1);
      line->SetLineStyle(2);
      line->Draw();

      TLine *line2 = nullptr;
      std::pair<float,float> legendy;
      if (! (std::find(versions.begin(), versions.end(), "70") != versions.end()
	     and std::find(versions.begin(), versions.end(), "80") != versions.end()) )
	{
	  double factor = TMath::Sqrt(26./28.);
	  line2 = new TLine(bol,factor,eol,factor);
	  line2->SetLineColor(kMagenta);
	  line2->SetLineStyle(5);
	  line2->Draw("same");
	  div->Draw("p"); //draw data points above the lines
	  legendy = std::make_pair(0.9-0.12*(versions_s+1),0.9);
	}
      else
	legendy = std::make_pair(0.9-0.12*versions_s,0.9);
      
      p1->cd();
      auto legend = new TLegend(0.55,legendy.first,0.9,legendy.second); //0.15 heigth per entry
      legend->SetTextSize(0.07);
      legend->AddEntry(resoV[idx1],vmap[versions[0]].c_str(),"ep");
      legend->AddEntry(resoV[idx2],vmap[versions[1]].c_str(),"ep");
      std::string ratiostr = vmap[versions[0]] + " / " + vmap[versions[1]];
      legend->AddEntry(div,ratiostr.c_str(),"ep");
      if ( !(std::find(versions.begin(), versions.end(), "70") != versions.end()
	     and std::find(versions.begin(), versions.end(), "80") != versions.end()) )
	legend->AddEntry(line2,"sqrt(26/28)","l");
      legend->Draw();

      foutReso->cd();
      div->Write();
      foutReso->Write();
      
      c->SaveAs((dirInBase2 + title + ".png").c_str());
      c->SaveAs((dirInBase2 + title + ".pdf").c_str());
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
