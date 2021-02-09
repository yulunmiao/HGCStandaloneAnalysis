#include<iostream>
#include<iomanip>
#include<sstream>
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

std::string etastr(float eta) {
  return std::to_string(static_cast<int>(eta*10.));
}
  
void check_func(TF1* funcptr) {
  if(!funcptr) {
    std::cout << "No function."  << std::endl;
    std::exit(1);
  }
}

void custom_cd(TFile* fileptr) {
  if(fileptr)
    fileptr->cd();
  else {
    std::cout << "Problem reading the file."  << std::endl;
    std::exit(1);
  }
}

void setCanvas(TCanvas *c) {
  c->SetRightMargin(0.09);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);
  c->Draw();
}

int plotBackLeakComparison(const InputParserPlotEGResoEtas& ip, std::string thisvers) {
  std::vector<float> etas = ip.etas();
  const unsigned etas_s = etas.size();

  std::unordered_map<unsigned,float> radius_map;
  radius_map[4] = 26.0;

  TGraphErrors *reso1[etas_s], *reso2[etas_s];

  std::string dirIn = ( "/eos/user/b/bfontana/www/RemoveLayers/" + ip.tag() + "/version" +
			thisvers + "/model2/gamma/SR4/" );

  std::unordered_map<std::string, std::string> vmap;
  vmap["60"] = "TDR (no neutron moderator)";
  vmap["70"] = "Scenario 13";

  std::string name = "c" + thisvers;    
  TCanvas *c[etas_s];
  TLegend *legend[etas_s];

  for(unsigned ieta(0); ieta<etas_s; ++ieta) {
    std::string title = "ResoOverlayedBackCor_" + thisvers + "_" + etastr(etas[ieta]);
    c[ieta] = new TCanvas((name+"_"+etastr(etas[ieta])).c_str(),
			  title.c_str(), 800, 600);
    setCanvas(c[ieta]);
    legend[ieta] = new TLegend(0.42,0.76,0.91,0.9);
    legend[ieta]->SetTextSize(0.05);

    std::string fileIn_ = dirIn + "IC3_pu0_SR4_Eta" ;
    std::string tmp_ =  std::to_string(static_cast<int>(etas[ieta]*10.f));
    std::string fileIn1 = fileIn_ + tmp_ + "_vsE_backLeakCor_raw.root";
    std::string fileIn2 = fileIn_ + tmp_ + "_vsE_backLeakCor.root";

    TFile *fIn1 = TFile::Open(fileIn1.c_str(),"READ");
    custom_cd(fIn1);
    reso1[ieta] = (TGraphErrors*)gDirectory->Get("resoRecoFitRaw")->Clone(("resoRecoFitRaw_"+std::to_string(ieta)).c_str());
    TF1 *func1 = reso1[ieta]->GetFunction("resoRaw");
    check_func(func1);
    func1->SetLineColor(2);
    reso1[ieta]->SetMarkerColor(2);
    reso1[ieta]->SetMaximum(0.1);
    reso1[ieta]->Draw("ap");

    TFile *fIn2 = TFile::Open(fileIn2.c_str(),"READ");
    custom_cd(fIn2);
    reso2[ieta] = (TGraphErrors*)gDirectory->Get("resoRecoFit")->Clone(("resoRecoFit_"+std::to_string(ieta)).c_str());
      
    TF1 *func2 = reso2[ieta]->GetFunction("reso");
    check_func(func2);
    func2->SetLineColor(3);
    reso2[ieta]->SetMarkerColor(3);
    reso2[ieta]->SetMaximum(0.1);
    reso2[ieta]->Draw("p same");

    std::stringstream etavalstr;
    etavalstr << std::fixed << std::setprecision(1) << etas[ieta];
    legend[ieta]->AddEntry(reso1[ieta], "before back leakage correction", "p");
    legend[ieta]->AddEntry(reso2[ieta], "after back leakage correction", "p");
    legend[ieta]->Draw("same");

    char buf1[500];
    sprintf(buf1,"#gamma, PU 0");
    TLatex lat1;
    lat1.SetTextSize(0.04);
    lat1.DrawLatexNDC(0.20,0.85,buf1);
    sprintf(buf1,"r = %3.0f mm", radius_map[ip.signalRegions()[0]]);
    lat1.DrawLatexNDC(0.20,0.80,buf1);
    sprintf(buf1,("|#eta| = " + etavalstr.str()).c_str());
    lat1.DrawLatexNDC(0.20,0.75,buf1);
    sprintf(buf1,vmap[thisvers].c_str());
    lat1.DrawLatexNDC(0.20,0.7,buf1);
    lat1.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
	    
    c[ieta]->SaveAs((dirIn + title + ".png").c_str());
  }

  for(unsigned ieta(0); ieta<etas_s; ++ieta) {    
    delete legend[ieta];
    delete c[ieta];
  }
    
  return 0;
}

int main(int argc, char** argv)
{
  gStyle->SetOptStat(kFALSE);
  InputParserPlotEGResoEtas ip(argc, argv);
  
  if(ip.signalRegions().size() != 1)
    std::cout << "This code is not ready to take more than one signal region. " <<  std::endl;
  
  for(unsigned v(0); v<ip.versions().size(); ++v)
    plotBackLeakComparison(ip, ip.versions()[v]);

  return 0;
}
