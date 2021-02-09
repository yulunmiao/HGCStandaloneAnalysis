#include<iostream>
#include <iomanip>
#include <sstream>
#include<string>
#include<vector>
#include<unordered_map>

#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLine.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#include "../include/InputParserBase.h"

class InputParserPlotEGResoEtas: public InputParserBase {
public:
  explicit InputParserPlotEGResoEtas(int argc, char** argv): InputParserBase(argc, argv) { run(); }

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
    required_args_ = { "--versions", "--tag", "--etas", "--signalRegions" };

    valid_args_v_["--versions"] = {"60", "70"};
    valid_args_v_["--signalRegions"] = {"0", "1", "2", "3", "4", "5"};
    free_args_ = {"--tag"};
    free_args_v_ = {"--etas"};
    optional_args_ = {""};
  }
};

int plotEtas(const InputParserPlotEGResoEtas& ip, std::string thisvers) {
  std::vector<float> etas = ip.etas();
  const unsigned etas_s = etas.size();

  std::unordered_map<unsigned,float> radius_map;
  radius_map[4] = 26.0;

  std::unordered_map<std::string, std::string> vmap;
  vmap["60"] = "TDR (no neutron moderator)";
  vmap["70"] = "Scenario 13";

  TGraphErrors *resoV[etas_s];

  std::string dirIn = ( "/eos/user/b/bfontana/www/RemoveLayers/" + ip.tag() + "/version" +
			thisvers + "/model2/gamma/SR4/" );


    auto legend = new TLegend(0.77,0.74,0.91,0.9);
    legend->SetTextSize(0.05);
  
    std::string name = "c" + thisvers;
    std::string title = "ResoOverlayedEtas_" + thisvers;
    TCanvas *c = new TCanvas(name.c_str(), title.c_str(), 1000, 600);
    c->SetRightMargin(0.09);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->Draw();
    for(unsigned ieta(0); ieta<etas_s; ++ieta) {

      std::string fileIn_ = dirIn + "IC3_pu0_SR4_Eta" ;
      std::string tmp_ =  std::to_string(static_cast<int>(etas[ieta]*10.f));
      std::string fileIn = fileIn_ + tmp_ + "_vsE_backLeakCor_raw.root";

      TFile *fIn = TFile::Open(fileIn.c_str(),"READ");
      if(fIn)
	fIn->cd();
      else {
	std::cout << "Problem reading the file."  << std::endl;
	return 1;
      }
      resoV[ieta] = (TGraphErrors*)gDirectory->Get("resoRecoFitRaw")->Clone(("resoRecoFitRaw_"+std::to_string(ieta)).c_str());
      //resoV[ieta]->GetFunction("resoRaw")->SetBit(TF1::kNotDraw);
      resoV[ieta]->GetFunction("resoRaw")->SetLineColor(ieta+2);
      resoV[ieta]->SetMarkerColor(ieta+2);
      resoV[ieta]->SetMaximum(0.1);
      std::string opt = ieta==0 ? "ap" : "p same";
      resoV[ieta]->Draw(opt.c_str());

      std::stringstream etavalstr;
      etavalstr << std::fixed << std::setprecision(1) << etas[ieta];
      legend->AddEntry(resoV[ieta], ("eta="+etavalstr.str()).c_str(), "p");
      legend->Draw("same");  
    }
  
    char buf1[500];
    sprintf(buf1,"#gamma, PU 0");
    TLatex lat1;
    lat1.SetTextSize(0.04);
    lat1.DrawLatexNDC(0.20,0.85,buf1);
    sprintf(buf1,"r = %3.0f mm", radius_map[ip.signalRegions()[0]]);
    lat1.DrawLatexNDC(0.20,0.80,buf1);
    sprintf(buf1,vmap[thisvers].c_str());
    lat1.DrawLatexNDC(0.20,0.75,buf1);
    lat1.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

    c->SaveAs((dirIn + title + ".png").c_str());

    return 0;
}

int main(int argc, char** argv)
{
  gStyle->SetOptStat(kFALSE);
  InputParserPlotEGResoEtas ip(argc, argv);
  
  if(ip.signalRegions().size() != 1)
    std::cout << "This code is not ready to take more than one signal region. " <<  std::endl;
  
  for(unsigned v(0); v<ip.versions().size(); ++v)
    plotEtas(ip, ip.versions()[v]);

  return 0;
}
