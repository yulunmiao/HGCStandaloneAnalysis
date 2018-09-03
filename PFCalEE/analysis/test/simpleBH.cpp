#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<map>
#include<utility>      // std::pair
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TPaletteAxis.h"
#include "TColor.h"
#include "TROOT.h"
#include "TRint.h"
#include "TMultiGraph.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSPUenergy.hh"
//#include "HGCSSSimpleHit.hh"

#include "PositionFit.hh"
#include "SignalRegion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

double DeltaR(double eta1,double phi1,double eta2,double phi2){
  double dr=99999.;
  double deta=fabs(eta1-eta2);
  double dphi=fabs(phi1-phi2);
  if(dphi>TMath::Pi()) dphi=2.*TMath::Pi()-dphi;
  dr=sqrt(deta*deta+dphi*dphi);
  return dr;
}


int main(int argc, char** argv){//main  


  //Input output and config options
  std::string cfg;
  unsigned pNevts;
  //std::string inFilePath;
  std::string outFilePath;
  std::string filePath;
  std::string digifilePath;
  unsigned nRuns;
  std::string simFileName;
  std::string recoFileName;
  unsigned debug;
  //double etamean;
  //double deta;
  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("outFilePath,o",  po::value<std::string>(&outFilePath)->required())
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("digifilePath", po::value<std::string>(&digifilePath)->default_value(""))
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    //("etamean,e",      po::value<double>(&etamean)->default_value(2.8))
    //("deta",      po::value<double>(&deta)->default_value(0.05))
    ;
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);


  std::string inFilePath = filePath+simFileName;

  std::cout << " -- Input parameters: " << std::endl
            << " -- Input file path: " << filePath << std::endl
            << " -- Digi Input file path: " << digifilePath << std::endl
    	    << " -- Output file path: " << outFilePath << std::endl
    //	    << " -- mean eta: " << etamean 
	    << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events per run." << std::endl;

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;


  //
  // hardcoded
  //


  //global threshold to reduce size of noise hits
  const double threshMin = 0.5;

  std::cout << " ---- Selection settings: ---- " << std::endl
	    << " -------threshMin " << threshMin << std::endl
	    << " ------------------------------------------" << std::endl;



  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////


  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;
  std::ostringstream inputrec;
  if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
  else
    inputrec << digifilePath << "/" << recoFileName;

  std::cout << inputsim.str() << " " << inputrec.str() << std::endl;

  HGCSSInfo * info;

  TChain *lSimTree = new TChain("HGCSSTree");
  TChain *lRecTree = 0;

  TFile * simFile = 0;
  TFile * recFile = 0;

  if (recoFileName.find("Digi") != recoFileName.npos)
    lRecTree = new TChain("RecoTree");
  else lRecTree = new TChain("PUTree");

  if (nRuns == 0){
    if (!testInputFile(inputsim.str(),simFile)) return 1;
    lSimTree->AddFile(inputsim.str().c_str());
    if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
    else {
      std::cout << " -- Error in getting information from simfile!" << std::endl;
      return 1;
    }
    if (!testInputFile(inputrec.str(),recFile)) return 1;
    lRecTree->AddFile(inputrec.str().c_str());
  }
  else {
    for (unsigned i(0);i<nRuns+1;++i){
      std::ostringstream lstrsim;
      std::ostringstream lstrrec;
      lstrsim << inputsim.str() << "_run" << i << ".root";
      if (testInputFile(lstrsim.str(),simFile)){
        if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
        else {
	  std::cout << " -- Error in getting information from simfile!" << std::endl;
          return 1;
        }
      }
      else continue;
      lstrrec << inputrec.str() << "_run" << i << ".root";
      if (!testInputFile(lstrrec.str(),recFile)) continue;
      lSimTree->AddFile(lstrsim.str().c_str());
      lRecTree->AddFile(lstrrec.str().c_str());
    }
  }

  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  //assert(info);

  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  const unsigned shape = info->shape();
  const double cellSize = info->cellSize();
  const double calorSizeXY = info->calorSizeXY();

  bool isTBsetup = (model != 2);
  bool bypassR = false;
  if (isTBsetup) bypassR = true;

  HGCSSDetector & myDetector = theDetector();
  myDetector.buildDetector(versionNumber,true,false,bypassR);

  //corrected for Si-Scint overlap
  const unsigned nLayers = 52;//


  std::cout << " -- Calor size XY = " << calorSizeXY
	    << ", version number = " << versionNumber 
	    << ", model = " << model << std::endl
	    << " -- cellSize = " << cellSize
	    << ", shape = " << shape
	    << ", nLayers = " << nLayers
	    << std::endl;
  HGCSSGeometryConversion geomConv(model,cellSize,bypassR,3);

  geomConv.setXYwidth(calorSizeXY);
  geomConv.setVersion(versionNumber);
  
  if (shape==2) geomConv.initialiseDiamondMap(calorSizeXY,10.);
  else if (shape==3) geomConv.initialiseTriangleMap(calorSizeXY,10.*sqrt(2.));
  else if (shape==1) geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  else if (shape==4) geomConv.initialiseSquareMap(calorSizeXY,10.);

  //square map for BHCAL
  geomConv.initialiseSquareMap1(1.3,3.0,-1.*TMath::Pi(),TMath::Pi(),0.01745);//eta phi segmentation
  geomConv.initialiseSquareMap2(1.3,3.0,-1.*TMath::Pi(),TMath::Pi(),0.02182);//eta phi segmentation
  std::vector<unsigned> granularity;
  granularity.resize(myDetector.nLayers(),1);
  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();




  std::ostringstream label;
  //root doesn't like . in branch names.....
  //label << "hits";
  //miptree->Branch(label.str().c_str(),"std::vector<HGCSSSimpleHit>",&miphitvec);

  double energy_max=-1.;

  TH1F* h_energy = new TH1F("h_energy","hit energy",1000,0.,5.);
  //TH1F* h_z = new TH1F("h_z","z of hit",5000,3100.,5200);
  //TH1F* h_z1 = new TH1F("h_z1","z of hit",5000,3150.,3550);
  //TH1F* h_z2 = new TH1F("h_z2","z of hit",5000,3550.,5200);
  //TH2F* h_xy = new TH2F("h_xy","xy of hit",1000,-2000,2000,1000,-2000,2000);
  //TH2F* h_etaphi = new TH2F("h_etaphi","etaphi of hit",1000,1,3.5,1000,-7,7);
  TH2F* h_getaphi = new TH2F("h_getaphi","gen part etaphi of hit",1000,1,3.5,1000,-7,7);
  //TH1F* h_l = new TH1F("h_l","layer of hit",80,0.,80.);
  //TH1F* h_l2 = new TH1F("h_l2","layer of hit",30,50,80.);
  //TH2F* h_zl = new TH2F("h_zl","z vs l of hit",5000,4300.,5200,25,30.,55.);
  /////////////////////////
  //Bryans analysis stuff//
  /////////////////////////
 

  TH2Poly* map_1 = new TH2Poly();
  //map_1->Honeycomb(-2803.17,-2790.5,6.49635,575,497);
  //map_1 = geomConv.hexagonMap();
  map_1 = geomConv.hexagonMap();
  map_1->SetTitle("map_1");
  map_1->SetName("map_1");

  //TH2Poly* map_1_0 = (TH2Poly*) map_1->Clone("map_1_0"); map_1_0->SetTitle("map_1_0");  map_1_0->SetName("map_1_0"); 
  //TH2Poly* map_1_1 = (TH2Poly*) map_1->Clone("map_1_1"); map_1_1->SetTitle("map_1_1");  map_1_1->SetName("map_1_1");

  char title[100];
  TH2Poly* map_1_layer[52];
  bool develop=false;
  if (develop){
    for (int ilayer=36;ilayer<=51;ilayer++){
      sprintf(title,"map_1_%d",ilayer);
      //map_1_layer[ilayer] = new TH2Poly();
      //map_1_layer[ilayer] = geomConv.hexagonMap();
      map_1_layer[ilayer] = (TH2Poly*)map_1->Clone(title);
      map_1_layer[ilayer]->SetTitle(title);
      map_1_layer[ilayer]->SetName(title);
    }
  }

  //TH2Poly* map_2 = new TH2Poly();
  //map_2 = geomConv.squareMap1();
  //map_2 = geomConv.squareMap1();

  //TH2Poly* map_3 = new TH2Poly();
  //map_3 = geomConv.squareMap2();
  //map_3 = geomConv.squareMap2();

  //---- xmin = 1.4, ymin = -3.14159 side = 0.01745, nx = 91, ny=360
  //---- xmin = 1.4, ymin = -3.14159 side = 0.02182, nx = 73, ny=287
  //double rbins2[92];
  //double rbins3[74];
  //TH2F* map_TH2F_2[4];
  //TH2F* map_TH2F_3[12];
  Float_t z_layer[]={
    3198.0,    3207.1,    3222.4,    3231.5,    3246.8,
    3255.9,    3271.2,    3280.3,    3295.6,    3304.7,
    3320.0,    3329.1,    3344.4,    3353.5,    3368.8,
    3377.9,    3393.2,    3402.3,    3417.6,    3426.7,
    3442.0,    3451.1,    3466.4,    3475.5,    3490.8,
    3499.9,    3515.2,    3524.3,    3577.4,    3626.4,
    3675.4,    3724.4,    3773.4,    3822.4,    3871.4,
    3920.4,    3969.4,    4020.3,    4071.2,    4122.1,
    4206.0,    4289.9,    4373.8,    4457.7,    4541.6,
    4625.5,    4709.4,    4793.3,    4877.2,    4961.1,
    5045.0,    5128.9,       0.0,    3971.2,    4022.1,
    4073.0,    4123.9,    4207.8,    4291.7,    4375.6,
    4459.5,    4543.4,    4627.3,    4711.2,    4795.1,
    4879.0,    4962.9,    5046.8,    5130.7};

  int eta_index[69]={
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,19,24, 25,30,29,32,37,
    41,41,42,43,44, 45,45,46,47};


  double minEta[52] = { 1.461, 1.464, 1.463, 1.466, 1.465, 1.468, 1.467, 1.469, 1.469, 1.471,
			1.471, 1.473, 1.472, 1.475, 1.474, 1.477, 1.476, 1.478, 1.478, 1.480,
			1.479, 1.482, 1.481, 1.483, 1.483, 1.485, 1.484, 1.487,
			1.487, 1.490, 1.494, 1.497, 1.500, 1.503, 1.506, 1.493,
			1.474, 1.455, 1.437, 1.420, 1.403, 1.376, 1.351, 1.327, 1.306, 1.316,
			1.332, 1.348, 1.364, 1.379, 1.395, 1.410 };
  
  int bin_exclude=2;
  double r_cut[69]; // r_cut values to exclude a few boundary hits
  double r_noisecut[69] = { 1567.5, 1567.5, 1575.6, 1575.6, 1583.7, 1583.7, 1591.8, 1591.8, 1599.9, 1599.9,
			    1608.0, 1608.0, 1616.1, 1616.1, 1624.2, 1624.2, 1632.3, 1632.3, 1640.4, 1640.4,
			    1648.5, 1648.5, 1656.6, 1656.6, 1664.7, 1664.7, 1672.8, 1672.8,
			    1696.9, 1713.5, 1730.1, 1746.7, 1763.3, 1779.9, 1796.4, 1844.2,
			    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			    0, 0, 0, 0, 0, 0,
			    0,
			    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			    0, 0, 0, 0, 0, 0 }; //cut values to exlude hits outside of the detector
  // BH fine part
  for (int ilayer=36;ilayer<=39;ilayer++){
    int ilayer_org=ilayer+17;
    //for (int ibin=0;ibin<=91;ibin++){
    //  double eta=1.4+double(ibin)*0.01745;
    //  double z = z_layer[ilayer_org];
    //  rbins2[91-ibin]=z*tan(2.*atan(exp(-eta))); 
    //}    
    sprintf(title,"map_TH2F_2_%d",ilayer);
    //map_TH2F_2[ilayer-36] = new TH2F(title,title,360,-1.*TMath::Pi(),3.1404,91,rbins2);
    double eta_tmp=1.4+double(eta_index[ilayer_org]-bin_exclude)*0.01745;
    double eta_tmp2=1.4+double(eta_index[ilayer_org])*0.01745;// to print out the boundary 
    double z_tmp = z_layer[ilayer_org];
    //r_cut[ilayer_org] = z_tmp*tan(2.*atan(exp(-eta_tmp))); //=====for inner excluded ring boundary
    r_noisecut[ilayer_org] = z_tmp*tan(2.*atan(exp(-(minEta[ilayer])))); // ===== for outer excluded ring boundary
    r_cut[ilayer_org] = z_tmp*tan(2.*atan(exp(-(minEta[ilayer]+bin_exclude*0.01745)))); // ====== for outer excluded ring boundary
    //std::cout<<"r min of scint at layer "<<ilayer_org<<" is: "<<z_tmp*tan(2.*atan(exp(-eta_tmp2)))<<std::endl;
    std::cout<<"r max of scint at layer "<<ilayer_org<<" is: "<< z_tmp*tan(2.*atan(exp(-minEta[ilayer])))<<std::endl;
    std::cout<<"r of excluded outer scint ring at layer "<<ilayer_org<<" is: "<< z_tmp*tan(2.*atan(exp(-(minEta[ilayer]+bin_exclude*0.01745))))<<std::endl;
  }  
  // BH coarse part
  for (int ilayer=40;ilayer<=51;ilayer++){
    int ilayer_org=ilayer+17;
    //for (int ibin=0;ibin<=73;ibin++){
    //  double eta=1.4+double(ibin)*0.02182;
    //  double z = z_layer[ilayer_org];
    //  rbins3[73-ibin]=z*tan(2.*atan(exp(-eta)));    
    //}    
    sprintf(title,"map_TH2F_3_%d",ilayer);
    //map_TH2F_3[ilayer-40] = new TH2F(title,title,287,-1.*TMath::Pi(),1.*TMath::Pi(),73,rbins3);
    double eta_tmp=1.4+double(eta_index[ilayer_org]-bin_exclude)*0.02182;
    double eta_tmp2=1.4+double(eta_index[ilayer_org])*0.02182;// to print out the boundary 
    double z_tmp = z_layer[ilayer_org];
    //r_cut[ilayer_org] = z_tmp*tan(2.*atan(exp(-eta_tmp)));//=====for inner excluded ring boundary
    r_noisecut[ilayer_org] = z_tmp*tan(2.*atan(exp(-(minEta[ilayer])))); // ===== for outer excluded ring boundary
    r_cut[ilayer_org] = z_tmp*tan(2.*atan(exp(-(minEta[ilayer]+bin_exclude*0.02182)))); // ====== for outer excluded ring boundary  
    //std::cout<<"r min of scint at layer "<<ilayer_org<<" is: "<<z_tmp*tan(2.*atan(exp(-eta_tmp2)))<<std::endl;
    std::cout<<"r max of scint at layer "<<ilayer_org<<" is: "<< z_tmp*tan(2.*atan(exp(-minEta[ilayer])))<<std::endl;
    std::cout<<"r of excluded outer scint ring at layer "<<ilayer_org<<" is: "<< z_tmp*tan(2.*atan(exp(-(minEta[ilayer]+bin_exclude*0.02182))))<<std::endl;
  } 
  /*
  for (int i_test = 0 ; i_test<69 ; i_test++){
    std::cout<<r_cut[i_test]<<" for layer "<<i_test<<std::endl;
  };
  */

  TProfile* h_el = new TProfile("h_el","energy per layer",52,0,52);
  TProfile* h_el_off = new TProfile("h_el_off","energy per layer",52,0,52);
  TProfile* h_elw = new TProfile("h_elw","mip per layer",52,0,52);
  TProfile* h_elw_off = new TProfile("h_elw_off","mip per layer",52,0,52);
  double lambins[69] = {
    .246, .062, .032, .062, .032, .062, .032, .062, .032, .062,
    .032, .062, .032, .062, .032, .062, .032, .062, .032, .062,
    .032, .062, .032, .062, .032, .062, .032, .062, .292, .262,
    .262, .262, .262, .262, .262, .262, .262, .262, .262, .262,
    .461, .461, .461, .461, .461, .461, .461, .461, .461, .461,
    .461, .461, .572, .256, .256, .256, .256, .455, .455, .455,
    .455, .455, .455, .455, .455, .455, .455, .455, .455}; // radiation length in lambda per layer 

  double lambda[69] = {0};

  for (int i = 0 ; i <= 68 ; i++){
    lambda[i] = lambins [i];
    if (i > 0) lambins[i] += lambins[i-1];
  }
  
  double z_layer2[53]={
    3190.0,    3200.0,    3209.1,    3224.4,    3233.5,    3248.8,
    3257.9,    3273.2,    3283.3,    3297.6,    3306.7,
    3322.0,    3331.1,    3346.4,    3355.5,    3370.8,
    3379.9,    3395.2,    3404.3,    3419.6,    3428.7,
    3444.0,    3453.1,    3468.4,    3477.5,    3492.8,
    3501.9,    3517.2,    3526.3,    3579.4,    3628.4,
    3677.4,    3725.4,    3775.4,    3824.4,    3874.4,    3922.4,   
    3971.4,    4022.3,    4073.2,    4124.1,    4208.0,    4291.9,    4375.8,    4459.7,   
    4543.6,    4627.5,    4711.4,    4795.3,    4879.2,    4963.1,    5047.0,    5130.9};

  TProfile* h_elam = new TProfile("h_elam","mip vs lam",68,lambins);
  TProfile* h_elam_off = new TProfile("h_elam_off","mip vs lam",68,lambins);
  TProfile* h_elam2 = new TProfile("h_elam2","dE/dlambda vs lambda",68,lambins);

  Int_t binnum = sizeof(z_layer2)/sizeof(double)-1;
  TProfile* h_ezW = new TProfile("h_ezW","mip vs z",binnum,z_layer2);
  TProfile* h_ezW_off = new TProfile("h_ezW_off","mip vs z",binnum,z_layer2);

  //std::cout<<"accumilated radiation length in lambda per layer is: "<<std::endl;
  //for (int i = 1 ; i <= 68 ; i++){
  //  std::cout<<lambins[i]<<" for layer: "<<i<<std::endl;
  //}

  TH1F* h_egenreco_cut = new TH1F("h_egenreco_cut","E reco sum over gen",100,0.,2.);
  TH1F* h_eta = new TH1F("h_eta","eta of gen particle",120,1.6,2.8);
  TProfile* h_lostE = new TProfile("h_lostE","lost energy of eta",24,1.6,2.8,0,100);  
  TProfile* h_lostHit = new TProfile("h_lostHit","lost hits of eta",24,1.6,2.8,0,100000);
  TH1F* h_egenreco_scint = new TH1F("h_egenreco_scint","E reco sum (Scint)",100,0.,0.5);
  TH1F* h_egenreco_scint_cut = new TH1F("h_egenreco_scint_cut","E reco sum (Scint Cut)",100,0.,0.05);
  TH1F* h_egenreco_Si = new TH1F("h_egenreco_Si","E reco sum (Si)",100,0.,2.);
  TH2F* h_e_ScintvsSi = new TH2F("h_e_ScintvsSi","E Scint vs E Si",600,0,600,600,0,600);
  TH1F* h_egenreco_rare = new TH1F("h_egenreco_rare"," E reco sum over gen",100,0,2);
  TH1F* h_egenreco_rare_cut = new TH1F("h_egenreco_rare_cut"," E reco sum over gen",100,0,2);
  TH1F* h_ScintoverSi = new TH1F("h_ScintoverSi","E Scint over E Si",200,0,2);
  TH2F* h_tailEvents = new TH2F("h_tailEvents","ECone/ETot vs ECone/Egen",100,0,.2,100,0,.2);
  TProfile* h_elostoveregen = new TProfile("h_elostoveregen","ELost/Egen vs Eta",24,1.6,2.8,0,100);
  TH2F* h_correction = new TH2F("h_correction","ECone vs E Dead Ring",100,0,200,100,0,200);
  TH2F* h_EtotvsEgen = new TH2F("h_EtotvsEgen","ETot vs EGen",250,0,250,250,0,250);
  TH1F* h_EtotoverEgen = new TH1F("h_EtotoverEgen","Total reco Energy over Egen",100,0,2);
  TH1F* h_dRcounts = new TH1F("h_dRcounts","number rechits vs dR",340,0,3.4); 
  TProfile* h_EvsdeltaR = new TProfile("h_EvsdeltaR","E reco vs delta R",340,0,3.4);
  TProfile* h_EsumvsdeltaR = new TProfile("h_EsumvsdeltaR","E recosum vs delta R",40,0,1);
  TProfile* h_EvsdeltaRtail = new TProfile("h_EvsdeltaRtail","E reco vs delta R in the tail",340,0,3.4);
  TH2F* h_gentracks = new TH2F("h_gentracks","R vs Z",2400,2800,5200,2000,0,2000);
  TH1F* h_RareFraction = new TH1F("h_RareFraction","Rare Fraction of Events",100,0.0,1.0);

  // eta vs reso 
  //bool doetaranges = 0;
  //if (doetaranges)
  //{
  //  TH1F* h_egenreco1617     = new TH1F("h_egenreco1617"," ErecoSum over Egen, eta 1.6-1.7",100,0.,2.);
  //  TH1F* h_egenreco1718     = new TH1F("h_egenreco1718"," ErecoSum over Egen, eta 1.7-1.8",100,0.,2.);
  //  TH1F* h_egenreco1819     = new TH1F("h_egenreco1819"," ErecoSum over Egen, eta 1.8-1.9",100,0.,2.);
  //  TH1F* h_egenreco1920     = new TH1F("h_egenreco1920"," ErecoSum over Egen, eta 1.9-2.0",100,0.,2.);
  //  TH1F* h_egenreco2021     = new TH1F("h_egenreco2021"," ErecoSum over Egen, eta 2.0-2.1",100,0.,2.);
  //  TH1F* h_egenreco2122     = new TH1F("h_egenreco2122"," ErecoSum over Egen, eta 2.1-2.2",100,0.,2.);
  //  TH1F* h_egenreco2223     = new TH1F("h_egenreco2223"," ErecoSum over Egen, eta 2.2-2.3",100,0.,2.);
  //  TH1F* h_egenreco2324     = new TH1F("h_egenreco2324"," ErecoSum over Egen, eta 2.3-2.4",100,0.,2.);
  //  TH1F* h_egenreco2425     = new TH1F("h_egenreco2425"," ErecoSum over Egen, eta 2.4-2.5",100,0.,2.);
  //  TH1F* h_egenreco2526     = new TH1F("h_egenreco2526"," ErecoSum over Egen, eta 2.5-2.6",100,0.,2.);
  //  TH1F* h_egenreco2627     = new TH1F("h_egenreco2627"," ErecoSum over Egen, eta 2.6-2.7",100,0.,2.);
  //  TH1F* h_egenreco2728     = new TH1F("h_egenreco2728"," ErecoSum over Egen, eta 2.7-2.8",100,0.,2.);
  //  
  //  TH1F* h_egenreco1617_cut = new TH1F("h_egenreco1617_cut"," ErecoSum over Egen, eta 1.6-1.7",100,0.,2.);
  //  TH1F* h_egenreco1718_cut = new TH1F("h_egenreco1718_cut"," ErecoSum over Egen, eta 1.7-1.8",100,0.,2.);
  //  TH1F* h_egenreco1819_cut = new TH1F("h_egenreco1819_cut"," ErecoSum over Egen, eta 1.8-1.9",100,0.,2.);
  //  TH1F* h_egenreco1920_cut = new TH1F("h_egenreco1920_cut"," ErecoSum over Egen, eta 1.9-2.0",100,0.,2.);
  //  TH1F* h_egenreco2021_cut = new TH1F("h_egenreco2021_cut"," ErecoSum over Egen, eta 2.0-2.1",100,0.,2.);
  //  TH1F* h_egenreco2122_cut = new TH1F("h_egenreco2122_cut"," ErecoSum over Egen, eta 2.1-2.2",100,0.,2.);
  //  TH1F* h_egenreco2223_cut = new TH1F("h_egenreco2223_cut"," ErecoSum over Egen, eta 2.2-2.3",100,0.,2.);
  //  TH1F* h_egenreco2324_cut = new TH1F("h_egenreco2324_cut"," ErecoSum over Egen, eta 2.3-2.4",100,0.,2.);
  //  TH1F* h_egenreco2425_cut = new TH1F("h_egenreco2425_cut"," ErecoSum over Egen, eta 2.4-2.5",100,0.,2.);
  //  TH1F* h_egenreco2526_cut = new TH1F("h_egenreco2526_cut"," ErecoSum over Egen, eta 2.5-2.6",100,0.,2.);
  //  TH1F* h_egenreco2627_cut = new TH1F("h_egenreco2627_cut"," ErecoSum over Egen, eta 2.6-2.7",100,0.,2.);
  //  TH1F* h_egenreco2728_cut = new TH1F("h_egenreco2728_cut"," ErecoSum over Egen, eta 2.7-2.8",100,0.,2.);
  //}

 
  ///end of mapping per layer///
  /*
  TH2F* h_zx = new TH2F("h_zx","zx of hit",5000,3100,5200,1000,-2000,2000);
  TH2F* h_zx10000 = new TH2F("h_zx10000","zx of hit",10000,3100,5200,1000,-2000,2000);
  TH2F* h_zx1000 = new TH2F("h_zx1000","zx of hit",1000,3100,5200,1000,-2000,2000);
  TH3F* h_xyz = new TH3F("h_xyz","xyz of hit",100,-2000.,2000.,100,-2000.,2000.,500,3100.,5200.); // 3d histo
  TH2F* h_yz = new TH2F("h_yz","yz of hit",5000,3100,5200,1000,-2000,2000);
  TH2F* h_zx1 = new TH2F("h_zx1","zx of hit",5000,3150,3550,4000,-2000,2000);
  TH2F* h_yz1 = new TH2F("h_yz1","yz of hit",5000,3150,3550,4000,-2000,2000);
  TH2F* h_zx2 = new TH2F("h_zx2","zx of hit",5000,3550,5200,4000,-2000,2000);
  TH2F* h_yz2 = new TH2F("h_yz2","yz of hit",5000,3550,5200,4000,-2000,2000);
  */

  
  TH2F* h_zr = new TH2F("h_zr","zr of hit",2100,3100,5200,3000,0,3000);
  TH2F* h_simzr = new TH2F("h_simzr","zr of hit",2100,3100,5200,3000,0,3000);
  TH2F* h_lr = new TH2F("h_lr","lr of hit",55,0,55,3000,0,3000);
  //////layer histos///////
  /*
  TH2F* h_nszxl = new TH2F("h_nszxl","zx of hit not scint (layers)",1400,3800,5200,3000,-2000,2000);
  TH2F* h_szxl = new TH2F("h_szxl","zx of hit scint (layers)",1400,3800,5200,3000,-2000,2000);
  TH2F* h_nszx = new TH2F("h_nszx","zx of hit not scint",1900,3100,5200,4000,-2000,2000);
  TH2F* h_szx = new TH2F("h_szx","zx of hit scint",1900,3100,5200,4000,-2000,2000);
  TH2F* h_nsxy = new TH2F("h_nsxy","xy of hit not scint",2000,-2000,2000,2000,-2000,2000);
  TH2F* h_sxy = new TH2F("h_sxy","xy of hit scint",2000,-2000,2000,2000,-2000,2000);
  TH2F* h_nsxyl = new TH2F("h_nsxyl","xy of hit not scint (layers)",3000,-2000,2000,3000,-2000,2000);
  TH2F* h_sxyl = new TH2F("h_sxyl","xy of hit scint (layers)",3000,-2000,2000,3000,-2000,2000);
  

  TH2F* h_nszx36 = new TH2F("h_nszx36","zx of hit not scint",5000,3100,5200,1000,-1200,1201);
  TH2F* h_nszx37 = new TH2F("h_nszx37","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx38 = new TH2F("h_nszx38","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx39 = new TH2F("h_nszx39","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx40 = new TH2F("h_nszx40","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx41 = new TH2F("h_nszx41","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx42 = new TH2F("h_nszx42","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx43 = new TH2F("h_nszx43","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx44 = new TH2F("h_nszx44","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx45 = new TH2F("h_nszx45","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx46 = new TH2F("h_nszx46","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx47 = new TH2F("h_nszx47","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx48 = new TH2F("h_nszx48","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx49 = new TH2F("h_nszx49","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx50 = new TH2F("h_nszx50","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx51 = new TH2F("h_nszx51","zx of hit not scint",5000,3100,5200,1000,-1200,1200);

  TH2F* h_szx36 = new TH2F("h_szx36","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx37 = new TH2F("h_szx37","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx38 = new TH2F("h_szx38","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx39 = new TH2F("h_szx39","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx40 = new TH2F("h_szx40","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx41 = new TH2F("h_szx41","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx42 = new TH2F("h_szx42","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx43 = new TH2F("h_szx43","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx44 = new TH2F("h_szx44","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx45 = new TH2F("h_szx45","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx46 = new TH2F("h_szx46","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx47 = new TH2F("h_szx47","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx48 = new TH2F("h_szx48","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx49 = new TH2F("h_szx49","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx50 = new TH2F("h_szx50","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx51 = new TH2F("h_szx51","zx of hit scint",5000,3100,5200,1000,-1200,1200);

  
  ///////////END///////////

  TH2F* h_nsxy36 = new TH2F("h_nsxy36","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy37 = new TH2F("h_nsxy37","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy38 = new TH2F("h_nsxy38","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy39 = new TH2F("h_nsxy39","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy40 = new TH2F("h_nsxy40","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy41 = new TH2F("h_nsxy41","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy42 = new TH2F("h_nsxy42","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy43 = new TH2F("h_nsxy43","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy44 = new TH2F("h_nsxy44","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy45 = new TH2F("h_nsxy45","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy46 = new TH2F("h_nsxy46","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy47 = new TH2F("h_nsxy47","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy48 = new TH2F("h_nsxy48","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy49 = new TH2F("h_nsxy49","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy50 = new TH2F("h_nsxy50","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy51 = new TH2F("h_nsxy51","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);

  TH2F* h_sxy36 = new TH2F("h_sxy36","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy37 = new TH2F("h_sxy37","xy of hit scint",1000,-3000,3000,1000,-3000,3000);
  TH2F* h_sxy38 = new TH2F("h_sxy38","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy39 = new TH2F("h_sxy39","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy40 = new TH2F("h_sxy40","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy41 = new TH2F("h_sxy41","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy42 = new TH2F("h_sxy42","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy43 = new TH2F("h_sxy43","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy44 = new TH2F("h_sxy44","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy45 = new TH2F("h_sxy45","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy46 = new TH2F("h_sxy46","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy47 = new TH2F("h_sxy47","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy48 = new TH2F("h_sxy48","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy49 = new TH2F("h_sxy49","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy50 = new TH2F("h_sxy50","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy51 = new TH2F("h_sxy51","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  */
  TH2F* h_Egenreco = new TH2F("h_Egenreco","E reco sum versus gen",1000,0.,1000.,100,0.,20.);
  TH1F* h_egenreco = new TH1F("h_egenreco","E reco sum over gen",100,0.,2.);//changed from 20 to 2
  TH1F* h_egensim = new TH1F("h_egensim","E sim sum over gen",100,0.,2.);//changed from 20 to 2

  TH2F* h_EpCone = new TH2F("h_EpCone","Ereco/gen versus cone size",10,0.,1.,100,0.,2.); // changed from 20 to 2 do to e weighting 
  TH2F* h_EpPhi = new TH2F("h_EpPhi","Ereco/gen versus phi",100,-4.,4.,100,0.,2.); // changed from 20 to 2 
  //TH2F* h_etagenmax= new TH2F("h_etagenmax","eta gen vs max",100,1.,5.,100,1.,5.);
  //TH2F* h_phigenmax= new TH2F("h_phigenmax","phi gen vs max",100,-4,4.,100,-4.,4.);
  //TH1F* h_maxE = new TH1F("h_maxE","energy of highest energy hit",1000,0.,1000.);// changed from 5000 to 1000
  TH1F* h_ECone03 = new TH1F("h_ECone03","Sum energy cone 03",1000,0.,500.);// changed from 50000 to 500
  // histos from sarahs code //
  TH2F* h_banana = new TH2F("h_banana","banana plot",1000,0.,500.,1000,0.,500.);
  TH1F* h_fracBH = new TH1F("h_fracBH","fraction in BH",100,-01.,1.1);
  
  ///////////////////////////////////////////////////////
  //////////////////  start event loop
  //////////////////////////////////////////////////////

  //  const unsigned nEvts = ((pNevts > lRecTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lRecTree->GetEntries()) : pNevts) ;
  
  //std::cout << " -- Processing " << nEvts << " events out of " << 
  //   lRecTree->GetEntries()<< std::endl;


  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;

  std::cout << " -- Processing " << nEvts << " events out of " << lSimTree->GetEntries() << " " << lRecTree->GetEntries() <<"\n"<< std::endl;


  //loop on events
  HGCSSEvent * event = 0;
  HGCSSEvent * eventRec = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;



  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  lRecTree->SetBranchAddress("HGCSSEvent",&eventRec);
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lRecTree->GetBranch("nPuVtx")) lRecTree->SetBranchAddress("nPuVtx",&nPuVtx);

  

  unsigned ievtRec = 0;
  unsigned nSkipped = 0;
  std::vector<double> absW; // from Sarah's code fill with weight
  bool firstEvent = true;

  // ---------- Event loop starts ----------
  
  for (unsigned ievt(0); ievt<nEvts; ++ievt , ++ievtRec){//loop on entries
  //for (unsigned ievt(0);ievt<100;++ievt){ // just lookings at the first 100 events for now
    if (ievtRec>=lRecTree->GetEntries()) continue;


    if (debug && ievt%500 == 0) std::cout << std::endl<<"... Processing entry: " << ievt <<"\n"<< std::endl;
    //else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;


    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievtRec);
    if (nPuVtx>0 && eventRec->eventNumber()==0 && event->eventNumber()!=0) {
      std::cout << " skip !" << ievt << " " << ievtRec << std::endl;
      nSkipped++;
      continue;
    }
    //ievtRec++;
    //if ((*genvec).size() > 1 && false) // skip events that have multiple particles generated before the detector
    //  {
    //	printf("Event %d has %d generated particles before the detector. SKIPPING!!!\n",ievt,(*genvec).size());
    //	continue;
    //  }
    // start getting filling our weighting (from Sarah)

    if(firstEvent) {
      firstEvent=false;
      std::cout<<" size of ssvec of weights is "<<(*ssvec).size()<<std::endl;
      double absweight=0;      
      for (unsigned iL(0); iL<(*ssvec).size();++iL){
	if(iL<((*ssvec).size()-1)) {
	  unsigned next=iL+1;
	  absweight=(((*ssvec)[iL].voldEdx())+((*ssvec)[next].voldEdx()))/2. ;
	} else{
	  absweight+=(*ssvec)[iL].voldEdx()  ;
	}
	absW.push_back(absweight);
	absweight=0;
      }
      std::cout << " -- AbsWeight size: " << absW.size() << std::endl;
      std::cout<<" values are ";
      for (unsigned iL(0); iL<(*ssvec).size();++iL){
	std::cout<<" "<<absW[iL];
      }
      std::cout<<std::endl;
    }
    // end of weightings

    double ptgen=-1.;
    double Egen=-1.;
    double ptgenpx=-1.;
    double ptgenpy=-1.;
    double ptgenpz=-1.;
    double etagen=99999.;
    double phigen=99999.;
    int pidgen=-1;
    double massgen= -1;
    int trackid = -1;
    if((*genvec).size()>0) {
      pidgen=(*genvec)[0].pdgid();
      massgen=(*genvec)[0].mass();
      ptgenpx=(*genvec)[0].px()/1000.;
      ptgenpy=(*genvec)[0].py()/1000.;
      ptgenpz=(*genvec)[0].pz()/1000.;
      //trackid=(*genvec)[0].trackID();
      ptgen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy);
      //double theta=atan(ptgen/ptgenpz);
      //etagen=-log(tan(theta/2));
      etagen=(*genvec)[0].eta();
      //Egen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy+ptgenpz*ptgenpz);
      //phigen=atan2(ptgenpy,ptgenpx);
      phigen=(*genvec)[0].phi();
      //if(phigen<0) phigen=2.*TMath::Pi()+phigen;
    }
    h_getaphi->Fill(etagen,phigen);
    
    h_eta->Fill(etagen); // check statistics by eta distribution

    if(debug) {
      //std::cout<<" gen vec size is "<<(*genvec).size()<<std::endl;
      //std::cout<<" first gen   pt  "<<ptgen<<" egen  "<<Egen<<" pidgen  "<<pidgen<<" etagen  "<<etagen<<" phi gen "<<phigen<<  "mass gen "<< massgen<<std::endl;
      Egen=0.; // added to later sum all egen 
      double egentemp = 0;
      double tempPt   = 0;
      double maxPt    = 0;
      TLorentzVector tlzv;
      TLorentzVector tlzv_tmp;
      TLorentzVector tlzv0;
      bool muonPrint = false;

      tlzv.SetPtEtaPhiE(0.0,0.0,0.0,0.0);
      tlzv_tmp.SetPtEtaPhiE(0.0,0.0,0.0,0.0);
      tlzv0.SetPtEtaPhiE(0.0,0.0,0.0,0.0);
           
      for (unsigned iP(0); iP<(*genvec).size(); ++iP){
        //std::cout<<" gen particle "<<iP<<" is (pdgid) "<<(*genvec)[iP].pdgid()<<std::endl;
	double rtemp = sqrt((*genvec)[iP].x()*(*genvec)[iP].x()+(*genvec)[iP].y()*(*genvec)[iP].y());
	h_gentracks->Fill((*genvec)[iP].z(),rtemp);
	egentemp = sqrt((*genvec)[iP].px()/1000.*(*genvec)[iP].px()/1000.+(*genvec)[iP].py()/1000.*(*genvec)[iP].py()/1000.+(*genvec)[iP].pz()/1000.*(*genvec)[iP].pz()/1000.);
	tempPt = (*genvec)[iP].pt();
	//if ((*genvec)[iP].pdgid() == 13 || (*genvec)[iP].pdgid() == -13) muonPrint = true ;
	if (tempPt > maxPt)
	  {
	    maxPt = tempPt;
	    tlzv0.SetPtEtaPhiE((*genvec)[iP].pt(), (*genvec)[iP].eta(), (*genvec)[iP].phi(), (*genvec)[iP].E());
	  }
	//egentemp = (*genvec)[iP].E()/1000;
	Egen +=  egentemp;
	tlzv_tmp.SetPtEtaPhiE((*genvec)[iP].pt(), (*genvec)[iP].eta(), (*genvec)[iP].phi(), (*genvec)[iP].E());
	tlzv += tlzv_tmp;
	//if ((*genvec).size() > 1 ) std::cout<<"Gen particle "<<iP<<" is (pdgid) "<<(*genvec)[iP].pdgid()<<" with trackID: "<<(*genvec)[iP].trackID()<<" at eta: "<<(*genvec)[iP].eta()<<" with Pt: "<<(*genvec)[iP].pt()/1000<<std::endl;
	//std::cout<<"really real energy "<<Egen<<std::endl;
	//if( (*genvec)[iP].trackID() == 1 )
	//  {
	//    etagen = (*genvec)[iP].eta();
	//    phigen = (*genvec)[iP].phi();
	//  }
      }
      if (muonPrint)
	{
	  std::cout<<"=====MUON DETECTED!!!====="<<std::endl;
	  std::cout<<"Particles in this event: "<<(*genvec).size()<<std::endl;
	  std::cout<<"Max Pt is   "<<tlzv0.Pt()/1000<<std::endl;
	  std::cout<<"Total Pt is "<<tlzv.Pt()/1000<<std::endl;
	}
      //if( etagen == (*genvec)[0].eta() && (*genvec)[0].trackID() != 1 )
      //	{
      //	  std::cout << "NO PARTICLE WITH TRACKID OF 1\nEtagen = "<<etagen<<"\nPhigen = "<<phigen<<std::endl;
      //	}
      //std::cout<<"final egen sum is = "<<Egen<<std::endl;
      if ( true && ((tlzv0.Pt() / tlzv.Pt() < 0.95) || (tlzv0.Pt() / tlzv.Pt() > 1.05))) //switch true to false if doing quarks ==========================
	{
	  printf("======================================================================================================================================================\n");
	  printf("Leading particle in event %d has a total fractional Pt of %f, which contains less than 95 percent of the total transverse momentum. SKIPPING!!!\n",ievt,tlzv0.Pt() / tlzv.Pt());
	  printf("======================================================================================================================================================\n");
	  continue ;
	}
      etagen = tlzv.Eta();// get eta from summed tlorentz vectors 
      phigen = tlzv.Phi();// get phi from summed tlorentz vectors 
      //ptgen = Egen / TMath::CosH(etagen);
      ptgen  = tlzv.Pt()/1000;
      //std::cout<<"Pt for this events is = "<<ptgen<<std::endl<<std::endl;
      //if(etagen < 1.5 || etagen > 2.8)  continue;
      //===================================================================================================================================================================
            
    }


    double ErecoRatio = 0; // for low e reco su,
    double recoNotCone = 0; //for events outside of the cone 
    
   
    bool isScint = false;
    //if (debug) std::cout << " - Event contains " << (*rechitvec).size() << " rechits." << std::endl;

    // make some simple plots about all the rechits
    unsigned iMax=-1;
    double MaxE=-1.;
    
    //find rmin and rmax of rechits 

    //double rmin = 99999;
    //double rmax = -1;

    // ---------- Rechit loop starts ----------

    // for longitudinal energy profile 
   
    //int nHits[80];

    std::map<std::pair<int,int>,float> mymap_rechit;
    mymap_rechit.clear();
    
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      //if (lHit.noiseFraction()>0.5) continue;
      double leta = lHit.eta();
      //double lphi = lHit.phi();
      //std::cout<<"getting layer"<<std::endl;
      unsigned layer = lHit.layer();
      if (lHit.energy()>MaxE) {MaxE=lHit.energy(); iMax=iH;}
      if (debug>20) std::cout << " -- hit " << iH << " eta " << leta << std::endl; 

      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
      isScint = subdet.isScint;
      TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();
       
      //map_1->Fill(lHit.get_x(),lHit.get_y());

      unsigned cellid = map->FindBin(lHit.get_x(),lHit.get_y());

      //if (lHit.energy()>1000.) std::cout << "reco energy"<< lHit.energy() << std::endl;
      //std::cout << "x "<< lHit.get_x() << "\t y "<<lHit.get_y() << "\t z" << lHit.get_z()<< std::endl; // added by Bryan, prints out xyz of each reco hit
      //std::cout<<"reco energy " << lHit.energy()<< " reco eta "<<lHit.eta()<< " reco phi " << lHit.phi() << " reco layer " << lHit.layer()<<" reco noise ratio "<< lHit.noiseFraction()<<std::endl;
      double lenergy = 0;
      double lenergyNoW = 0;
      if (!isScint){
	lenergy = lHit.energy()*absW[layer]/1000.; // Si was weighted too high
	lenergyNoW = lHit.energy();
      }
      if (isScint){
	lenergy = lHit.energy()*absW[layer]/1000.; //Scint was weighted too low
	lenergyNoW = lHit.energy();
      }
      geomConv.fill(lHit.layer(),lenergyNoW,0,cellid,lHit.get_z());
      //double lenergy=lHit.energy()*absW[layer]/1000.; // weight added (from Sarah's code)
      double r_hit = sqrt(lHit.get_x()*lHit.get_x()+lHit.get_y()*lHit.get_y());
      //if (( r_hit > r_noisecut[layer]) && (layer < 36 || layer > 52)) continue;
      
      // printf added by bryan
      /*
      printf("|| reco hit# e, corre, eta, phi, lyr, noiseF: %5d %6.1f %6.1f %6.1f %6.1f %3d %8.3f\n",
	     iH,lHit.energy(),lenergy,lHit.eta(),lHit.phi(),lHit.layer(),lHit.noiseFraction());
      printf("|| reco hit# %d  \t",iH);
      printf("| reco energy = %f \t",lHit.energy());
      printf("| reco weighted E = %f \t", lenergy);
      printf("| reco eta = %f \t",lHit.eta());
      printf("| reco phi = %f \t",lHit.phi());
      ///
      //printf("| reco layer = %d \t",lHit.layer());
      //printf("| reco zhit = %d \t\n",lHit.get_x());

      //printf("| reco noise ratio = %f\t ||\n ",lHit.noiseFraction());
      //printf("| reco cellid = %d \t",cellid);

      */
      // nHits[layer] += 1 ;
      
      
      h_energy->Fill(lenergyNoW);
      //h_z->Fill(lHit.get_z());
      //h_z1->Fill(lHit.get_z());
      //h_z2->Fill(lHit.get_z());
      //h_l->Fill(lHit.layer()+0.5);
      int ixx=lHit.layer();
      if(ixx>52) ixx=ixx-17;
      //h_zl->Fill(lHit.get_z(),ixx);
      //h_l2->Fill(lHit.layer()+0.5);
      //h_xy->Fill(lHit.get_x(),lHit.get_y());
      h_zr->Fill(lHit.get_z(),r_hit); // added by Bryan
      h_lr->Fill(ixx,r_hit);
      //h_zx->Fill(lHit.get_z(),lHit.get_x()); //added by Bryan
      //h_yz->Fill(lHit.get_z(),lHit.get_y()); //added by Bryan
      //h_zx1->Fill(lHit.get_z(),lHit.get_x()); //added by Bryan
      //h_yz1->Fill(lHit.get_z(),lHit.get_y()); //added by Bryan
      //h_zx2->Fill(lHit.get_z(),lHit.get_x()); //added by Bryan
      //h_yz2->Fill(lHit.get_z(),lHit.get_y()); //added by Bryan
      //h_zx10000->Fill(lHit.get_z(),lHit.get_x());//added by Bryan
      //h_zx1000->Fill(lHit.get_z(),lHit.get_x());//added by Bryan
      //h_xyz->Fill(lHit.get_x(),lHit.get_y(),lHit.get_z());// added by Bryan
      
      //h_etaphi->Fill(lHit.eta(),lHit.phi());
      

      //std::cout<<"Layer of hit "<<lHit.layer()<< " at z = "<<lHit.get_z()<<std::endl;

      // mapping stuff

      /*
      if(!isScint && ixx == 37){
	if(r_hit > rmax){
	  rmax = r_hit;
	};
      }

      if(!isScint) {
	map_1->Fill(lHit.get_x(),lHit.get_y());
	if (develop && ixx>=36) map_1_layer[ixx]->Fill(lHit.get_x(),lHit.get_y());
      }
      //if (ixx<=20) map_1_0->Fill(lHit.get_x(),lHit.get_y());
      //else         map_1_1->Fill(lHit.get_x(),lHit.get_y());

      if(isScint){ 
	//if (dR < .3 && r_hit < r_cut[lHit.layer()]){
	//r_hit = 0;
	  //std::cout<<"EXLUDING RecHit is less than "<<r_cut[lHit.layer()]<<" for layer "<<lHit.layer()<<std::endl;
	//};
	if (ixx>=36&&ixx<=39)      map_TH2F_2[ixx-36]->Fill(lHit.phi(),r_hit);
	else if (ixx>=40&&ixx<=51) map_TH2F_3[ixx-40]->Fill(lHit.phi(),r_hit);
      }

      */
      //end of mapping stuff
      /*
      if(isScint)
	{
	  h_sxy->Fill(lHit.get_x(),lHit.get_y());
	  h_szx->Fill(lHit.get_z(),lHit.get_x());
	}
      else
	{
	  h_nsxy->Fill(lHit.get_x(),lHit.get_y());
	  h_nszx->Fill(lHit.get_z(),lHit.get_x());
	}
      */

      //int ilayer = ixx;

      // added in zx comps. Bryan//////
      /*
      if(ilayer==36) {
	if(isScint)
	  {
	    h_sxy36->Fill(lHit.get_x(),lHit.get_y());
	    h_szx36->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else 
	  { 
	    h_nsxy36->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx36->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==37) {
	if(isScint)
	  {
	    h_sxy37->Fill(lHit.get_x(),lHit.get_y());
	    h_szx37->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy37->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx37->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==38) {
	if(isScint)
	  {
	    h_sxy38->Fill(lHit.get_x(),lHit.get_y());
	    h_szx38->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy38->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx38->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==39) {
	if(isScint)
	  {
	    h_sxy39->Fill(lHit.get_x(),lHit.get_y());
	    h_szx39->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy39->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx39->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==40) {
	if(isScint)
	  {
	    h_sxy40->Fill(lHit.get_x(),lHit.get_y());
	    h_szx40->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy40->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx40->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==41) {
	if(isScint)
	  {
	    h_sxy41->Fill(lHit.get_x(),lHit.get_y());
	    h_szx41->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy41->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx41->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==42) {
	if(isScint)
	  {
	    h_sxy42->Fill(lHit.get_x(),lHit.get_y());
	    h_szx42->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy42->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx42->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==43) {
	if(isScint)
	  {
	    h_sxy43->Fill(lHit.get_x(),lHit.get_y());
	    h_szx43->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy43->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx43->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==44) {
	if(isScint)
	  {
	    h_sxy44->Fill(lHit.get_x(),lHit.get_y());
	    h_szx44->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy44->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx44->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==45) {
	if(isScint)
	  {
	    h_sxy45->Fill(lHit.get_x(),lHit.get_y());
	    h_szx45->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy45->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx45->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==46) {
	if(isScint)
	  {
	    h_sxy46->Fill(lHit.get_x(),lHit.get_y());
	    h_szx46->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy46->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx46->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==47) {
	if(isScint)
	  {
	    h_sxy47->Fill(lHit.get_x(),lHit.get_y());
	    h_szx47->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy47->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx47->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==48) {
	if(isScint)
	  {
	    h_sxy48->Fill(lHit.get_x(),lHit.get_y());
	    h_szx48->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy48->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx48->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==49) {
	if(isScint)
	  {
	    h_sxy49->Fill(lHit.get_x(),lHit.get_y());
	    h_szx49->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy49->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx49->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==50) {
	if(isScint)
	  {
	    h_sxy50->Fill(lHit.get_x(),lHit.get_y());
	    h_szx50->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy50->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx50->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==51) {
	if(isScint)
	  {
	    h_sxy51->Fill(lHit.get_x(),lHit.get_y());
	    h_szx51->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy51->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx51->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
v
      */
    }//loop on hits
    
    //fill long energy profile 

    //std::cout<<"nhits in layer 2:"<<nHits[2]<<std::endl;
    //std::cout<<"summed energy in layer 2:"<<penergy[2]<<std::endl;
    //std::cout<<"rmin in layer 37 is:"<<rmin<<std::endl;
    //std::cout<<"rmax in layer 37 is:"<<rmax<<std::endl;


      //if(nHits[ilayer] >= 1000) {

      //std::cout<<"Layer, "<<ilayer<<" has a summed rechit energy of "<<penergy[ilayer]<<std::endl;
      //};
    

    HGCSSRecoHit lHit = (*rechitvec)[iMax];
    double maxeta = lHit.eta();
    double maxphi=lHit.phi();
    //double maxE=lHit.energy();
    //h_maxE->Fill(maxE);
    //h_etagenmax->Fill(maxeta,etagen);
    //h_phigenmax->Fill(maxphi,phigen);
    if(debug>2) {
      std::cout<<" Max hit energy eta phi "<<lHit.energy()<<" "<<lHit.eta()<<" "<<lHit.phi()<<std::endl;
    }

    // make e/p plots for various cones around gen particle
    double rechitsumE03_cut=0; // sum of rechit energy with the exclusion of boundry scint layers 
    double rechitsum[12]={0}; // for eta vs reso plots
    double rechitsum_cut[12]={0}; // for eta vs reso plots with scint cut
    double rechitsum_lost=0;
    double rechitsum_scint=0;
    double rechitsum_Si=0;
    double rechitsumNoW[80] = {0};
    double penergy[80]={0};
    int    lostHit=0;
    double rechitsum_dr[40] = {0};

    double rechitsumE01=0.;
    double rechitsumE02=0.;
    double rechitsumE03=0.;
    double rechitsumE04=0.;
    double rechitsumE05=0.;
    double rechitBHsumE01=0.;// if a rechit is part of the scint 
    double rechitBHsumE02=0.;
    double rechitBHsumE03=0.;
    double rechitBHsumE04=0.;
    double rechitBHsumE05=0.;
    double rechitsumEWoNoise01=0.;
    double rechitsumEWoNoise02=0.;
    double rechitsumEWoNoise03=0.;
    double rechitsumEWoNoise04=0.;
    double rechitsumEWoNoise05=0.;
    double etaW=0.;
    double phiW=0.;
    double norm=0.;
    //double etaaxis=etagen;
    //double phiaxis=phigen;
    double etaaxis=etagen;
    double phiaxis=phigen;
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      //if (lHit.noiseFraction()>0.5) continue;
      double r_hit = sqrt(lHit.get_x()*lHit.get_x()+lHit.get_y()*lHit.get_y());
      unsigned layer = lHit.layer();
      double leta = lHit.eta();
      double lphi = lHit.phi();
      //if(lphi<0) lphi=2.*TMath::Pi()+lphi;
      double lenergy = 0;
      double lenergyNoW = 0;
      if (!isScint){
	lenergy = lHit.energy()*absW[layer]/1000.; // Si was weighted too high
	lenergyNoW = lHit.energy();
      }
      if (isScint){
	lenergy = lHit.energy()*absW[layer]/1000.; //Scint was weighted too low
	lenergyNoW = lHit.energy();
      }
      if (debug>20) std::cout << " -- hit " << iH << " et eta phi " << lenergy<<" "<<leta << " "<< lphi<<std::endl; 
	//clean up rechit collection

      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
      isScint = subdet.isScint;
      

      norm+=lenergy;
      etaW+=leta*lenergy;
      phiW+=lphi*lenergy;

      double dR=DeltaR(etaaxis,phiaxis,leta,lphi);
      //double dR=fabs(etagen-leta);

      double dr_iter = .025 ; double dr_start = 0; double dr_end = 1;
      
      for (int i = 0 ; i*dr_iter+dr_start < dr_end ; i++){
	  if ( dR >= i*dr_iter+dr_start && dR < i*dr_iter+dr_start + dr_iter){
	      rechitsum_dr[i]+=lenergy;
	  }
      }

      TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();

      //unsigned cellid = map->FindBin(lHit.get_x(),lHit.get_y());

      //
      //KH - make a map for selected rechits around gen pions and low noise fraction
      //KH - Warning - cellid not working well for scintilator layers (so not counted)
      //
      if (!isScint && lHit.noiseFraction()<0.5 && dR<0.3){//========================changed 0.4 for quark
	  unsigned cellid = map->FindBin(lHit.get_x(),lHit.get_y());
	  std::map<std::pair<int,int>, float>::iterator it = mymap_rechit.find(std::make_pair(layer,cellid)); 
	  if (it != mymap_rechit.end()) (*it).second += lenergyNoW;
	  else mymap_rechit.insert(std::make_pair(std::make_pair(layer,cellid),lenergyNoW));
      }

      if(debug>20) std::cout<<" dR "<<dR<<" "<<etagen<<" "<<phigen<<" "<<leta<<" "<<lphi<<std::endl;
      if(dR<0.1)
	{
	  rechitsumE01+=lenergy;
	  if (isScint) rechitBHsumE01+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise01+=lenergy;
	}
      if(dR<0.2)
	{
	  rechitsumE02+=lenergy;
	  if (isScint) rechitBHsumE02+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise02+=lenergy;
	}

      ///////////used to make energy reso plots//////////////
      if(dR<0.3)//==========================================change to .4 for quark gun
	{

	  //if (( r_hit > r_noisecut[layer]) && (layer < 36 || layer > 52)) continue;
	  
	  rechitsumE03+=lenergy;
	  rechitsumNoW[layer] += lenergyNoW;
	  penergy[layer]+=lenergy;
	  
	  if (isScint && r_hit <= r_cut[layer]){
	    rechitsum_scint+=lenergy;
	  }
	  if (!isScint){
	    rechitsum_Si+=lenergy;
	  }

	  double eTemp = lenergy;
	  if (isScint && r_hit > r_cut[layer]){
	    //std::cout<<"EXLUDING RecHit is less than "<<r_cut[layer]<<" for layer "<<layer<<std::endl;
	    eTemp = 0;// exclude boundry  scint layers
	    rechitsum_lost+=lenergy;
	    lostHit+=1;
	  }
	  rechitsumE03_cut+=eTemp;
	  

	  double etaStart = 1.6;
	  double etaEnd   = 2.8;
	  double etaEter  = .1;
	  for ( int ieta = 0 ; etaEter*ieta+etaStart < etaEnd ; ++ieta){
	    if ( etagen >= etaEter*ieta+etaStart && etagen <= etaEter*ieta+etaStart+.1){
	      //std::cout<<etaEter*ieta+etaStart<<" < "<<etagen<< " < "<<etaEter*ieta+etaStart+.1<<std::endl;
	      rechitsum[ieta]+=lenergy;
	      if ( r_hit > r_cut[layer] ) lenergy = 0;
	      rechitsum_cut[ieta]+=lenergy;
	    }
	  }
	  

	  if (isScint) rechitBHsumE03+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise03+=lenergy;
	}


      if(dR<0.4)
	{
	  rechitsumE04+=lenergy;
	  if (isScint) rechitBHsumE04+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise04+=lenergy;
	}
      if(dR<0.5)
	{
	  rechitsumE05+=lenergy;
	  if (isScint) rechitBHsumE05+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise05+=lenergy;
	}
      
      h_dRcounts->Fill(dR);
      h_EvsdeltaR->Fill(dR,lenergy);

    }//loop on hits
    //if(debug>1) {
    //std::cout<<" reco gen are "<<rechitsumE01<<" "<<rechitsumE02<<" "<<rechitsumE03<<" "<<rechitsumE04<<" "<<rechitsumE05<<" "<<Egen<<std::endl;
    //}
    //if(debug>5) std::cout<<"weighted eta phi are "<<etaW/norm<<" "<<phiW/norm<<std::endl;
    
    double dr_iter = .025 ; double dr_start = 0; double dr_end = 1;
    
    for (int i = 0 ; i*dr_iter+dr_start < dr_end ; i++){
      h_EsumvsdeltaR->Fill(i*dr_iter+dr_start,rechitsum_dr[i]);
    }

    if (MaxE>energy_max) energy_max=MaxE;
    ErecoRatio  = rechitsumE03/Egen;
    if (ErecoRatio < .1) std::cout<<"E reco sum: ="<<rechitsumE03<<"...E gen: "<<Egen<<std::endl;
    recoNotCone = rechitsumE03/norm;

    h_Egenreco->Fill(Egen,rechitsumE03/Egen);// changed from 5 to 3
    h_egenreco_cut->Fill(rechitsumE03_cut/Egen); // with scint cuts
    h_egenreco->Fill(rechitsumE03/Egen); // changed from 5 to 3
    h_lostE->Fill(etagen,rechitsum_lost); // fill lost energy per eta
    h_elostoveregen->Fill(etagen,rechitsum_lost/Egen);// average ratio of lost energy/gen energy vs eta
    h_lostHit->Fill(etagen,lostHit); // fill lost hits per eta
    h_egenreco_scint->Fill(rechitsum_scint/Egen); //dead scint eta ring 
    h_egenreco_scint_cut->Fill(rechitsum_lost/Egen);
    h_egenreco_Si->Fill(rechitsum_Si/Egen);
    h_e_ScintvsSi->Fill(rechitsum_Si,rechitsum_lost+rechitsum_scint);
    h_ScintoverSi->Fill(rechitBHsumE03/rechitsum_Si);
    h_correction->Fill(rechitsumE03_cut,rechitsumE03);
    
    
    
    if (rechitBHsumE03/(rechitBHsumE03+rechitsum_Si) >= .1){ //tuned---fill only when a substantial amount of energy deposited inside the scint

      h_egenreco_rare->Fill(rechitsumE03/Egen);
      h_egenreco_rare_cut->Fill(rechitsumE03_cut/Egen);
    }

    h_RareFraction->Fill(rechitBHsumE03/(rechitBHsumE03+rechitsum_Si));
    //printf("rechitsum_lost = %f || rechitsum_scint = %f || rechitsum_Si = %f || rechitsumE03 = %f \n",rechitsum_lost,rechitsum_scint,rechitsum_Si,rechitsumE03);
    double tempPenergy[80]= {0};
    double tempWenergy[80] = {0};

    for(int ilayer = 0 ; ilayer < 69 ; ++ilayer){
      h_el->Fill(ilayer+0.5,penergy[ilayer]);
      h_elw->Fill(ilayer+0.5,rechitsumNoW[ilayer]);
      if (ilayer < 53){
	tempPenergy[ilayer] = penergy[ilayer];
	tempWenergy[ilayer] = rechitsumNoW[ilayer];
      }
      if (ilayer >= 53){
	tempPenergy[ilayer-17] += penergy[ilayer];
	tempWenergy[ilayer-17] += rechitsumNoW[ilayer];
      }

    }
    for(int ilayer = 0 ; ilayer < 52 ; ++ilayer){
      h_el_off->Fill(ilayer+0.5,tempPenergy[ilayer]); // tprofile with no offset  
      h_elw_off->Fill(ilayer+0.5,tempWenergy[ilayer]);
      h_elam->Fill(lambins[ilayer],penergy[ilayer]);
      h_elam_off->Fill(lambins[ilayer],tempPenergy[ilayer]);
      h_elam2->Fill(lambins[ilayer],tempPenergy[ilayer]*2./(lambda[ilayer]+lambda[ilayer+1]));
      h_ezW->Fill(z_layer[ilayer],rechitsumNoW[ilayer]);
      h_ezW_off->Fill(z_layer[ilayer],tempWenergy[ilayer]);
      //printf("tempPenergy[%d] = %f \t lambins[%d] = %f \n", ilayer, tempPenergy[ilayer],ilayer,lambins[ilayer]);
    }
    //std::cout<<"at eta gen. "<<etagen<<", the lost energy and hits are "<<rechitsum_lost<<" and "<<lostHit<<" respectively."<<std::endl;
   // if (doetaranges)
   //   {
   //	double etaStart = 1.6;
   //
   //	double etaEter  = .1;
   //	
   //	//for (int index=0;index<12;++index){
   //	//printf("rechitsum[%d]/Egen = %f \n",index,rechitsum[index]/Egen); 
   //	//}
   //	if ( etagen >= etaEter*0+etaStart && etagen <= etaEter*0+etaStart+.1 ){
   //	  h_egenreco1617->Fill(rechitsum[0]/Egen);
   //	  h_egenreco1617_cut->Fill(rechitsum_cut[0]/Egen);
   //	}
   //	if ( etagen >= etaEter*1+etaStart && etagen <= etaEter*1+etaStart+.1 ){
   //	  h_egenreco1718->Fill(rechitsum[1]/Egen);
   //	  h_egenreco1718_cut->Fill(rechitsum_cut[1]/Egen);
   //	}
   //	if ( etagen >= etaEter*2+etaStart && etagen <= etaEter*2+etaStart+.1 ){
   //	  h_egenreco1819->Fill(rechitsum[2]/Egen);
   //	  h_egenreco1819_cut->Fill(rechitsum_cut[2]/Egen);
   //	}
   //	if ( etagen >= etaEter*3+etaStart && etagen <= etaEter*3+etaStart+.1 ){
   //	  h_egenreco1920->Fill(rechitsum[3]/Egen);
   //	  h_egenreco1920_cut->Fill(rechitsum_cut[3]/Egen);
   //	}
   //	if ( etagen >= etaEter*4+etaStart && etagen <= etaEter*4+etaStart+.1 ){
   //	  h_egenreco2021->Fill(rechitsum[4]/Egen);
   //	  h_egenreco2021_cut->Fill(rechitsum_cut[4]/Egen);
   //	}
   //	if ( etagen >= etaEter*5+etaStart && etagen <= etaEter*5+etaStart+.1 ){
   //	  h_egenreco2122->Fill(rechitsum[5]/Egen);
   //	  h_egenreco2122_cut->Fill(rechitsum_cut[5]/Egen);
   //	}
   //	if ( etagen >= etaEter*6+etaStart && etagen <= etaEter*6+etaStart+.1 ){
   //	  h_egenreco2223->Fill(rechitsum[6]/Egen);
   //	  h_egenreco2223_cut->Fill(rechitsum_cut[6]/Egen);
   //	}
   //	if ( etagen >= etaEter*7+etaStart && etagen <= etaEter*7+etaStart+.1 ){
   //	  h_egenreco2324->Fill(rechitsum[7]/Egen);
   //	  h_egenreco2324_cut->Fill(rechitsum_cut[7]/Egen);
   //	}
   //	if ( etagen >= etaEter*8+etaStart && etagen <= etaEter*8+etaStart+.1 ){
   //	  h_egenreco2425->Fill(rechitsum[8]/Egen);
   //	  h_egenreco2425_cut->Fill(rechitsum_cut[8]/Egen);
   //	}
   //	if ( etagen >= etaEter*9+etaStart && etagen <= etaEter*9+etaStart+.1 ){
   //	  h_egenreco2526->Fill(rechitsum[9]/Egen);
   //	  h_egenreco2526_cut->Fill(rechitsum_cut[9]/Egen);
   //	}
   //	if ( etagen >= etaEter*10+etaStart && etagen <= etaEter*10+etaStart+.1 ){
   //	  h_egenreco2627->Fill(rechitsum[10]/Egen);
   //	  h_egenreco2627_cut->Fill(rechitsum_cut[10]/Egen);
   //	}
   //	if ( etagen >= etaEter*11+etaStart && etagen <= etaEter*11+etaStart+.1 ){
   //	  h_egenreco2728->Fill(rechitsum[11]/Egen);
   //	  h_egenreco2728_cut->Fill(rechitsum_cut[11]/Egen);
   //	}
   //   }
   
    
    h_EpPhi->Fill(phigen,rechitsumE03/Egen); // changed from 2 to 3
    h_EpCone->Fill(0.1,rechitsumE01/Egen);
    h_EpCone->Fill(0.2,rechitsumE02/Egen);
    h_EpCone->Fill(0.3,rechitsumE03/Egen);
    h_EpCone->Fill(0.4,rechitsumE04/Egen);
    h_EpCone->Fill(0.5,rechitsumE05/Egen);
    
    h_banana->Fill(rechitsumE03-rechitBHsumE03,rechitBHsumE03); // added from sarah's code
    double frac =-0.05; // from sarah
    //double notBH=rechitsumE03-rechitBHsumE03; // from sarah unused 
    if(rechitsumE03>0) frac=rechitBHsumE03/rechitsumE03; // from sarah
    h_fracBH->Fill(frac); // from sarah

    h_ECone03->Fill(rechitsumE03);
      //miptree->Fill();

    // ---------- Simhit loop starts ----------

    double simhitsumE01=0.;
    double simhitsumE02=0.;
    double simhitsumE03=0.;
    double simhitsumE04=0.;
    double simhitsumE05=0.;
    double simhitBHsumE01=0.;
    double simhitBHsumE02=0.;
    double simhitBHsumE03=0.;
    double simhitBHsumE04=0.;
    double simhitBHsumE05=0.;

    //initialise calibration class
    //KH - myDigitiser
    HGCSSDigitisation myDigitiser;
    //const unsigned interCalib = 3; // check against generation setting    
    //myDigitiser.setIntercalibrationFactor(interCalib);
    const unsigned nSiLayers = 2;  // this is what I see in generation, but why?
    //std::cout << "KHKH: inFilePath,bypassR,nSiLayers: " << inFilePath<<" "<<bypassR<<" "<<nSiLayers << std::endl;
    HGCSSCalibration mycalib(inFilePath,bypassR,nSiLayers);
    mycalib.setVertex(eventRec->vtx_x(),eventRec->vtx_y(),eventRec->vtx_z());

    int nhit=0, nhitscinti=0;
    std::map<std::pair<int,int>,float> mymap_simhit; mymap_simhit.clear();
    std::map<std::pair<int,int>,float> mymap_simhit_x; mymap_simhit_x.clear();
    std::map<std::pair<int,int>,float> mymap_simhit_y; mymap_simhit_y.clear();
    std::map<std::pair<int,int>,float> mymap_simhit_z; mymap_simhit_z.clear();    
    //std::map<std::pair<int,int>,float> mymap_simhit_eta; mymap_simhit_eta.clear();
    //std::map<std::pair<int,int>,float> mymap_simhit_phi; mymap_simhit_phi.clear();
    std::map<std::pair<int,int>,float> mymap_simhit_MeVToMip; mymap_simhit_MeVToMip.clear();

    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];

      unsigned layer = lHit.layer();
      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
      DetectorEnum type = subdet.type;
      isScint = subdet.isScint;
    
      std::pair<double,double> xy = lHit.get_xy(subdet,geomConv,shape);
      double posx = xy.first;//lHit.get_x(cellSize);
      double posy = xy.second;//lHit.get_y(cellSize);
      double posz = lHit.get_z();
      double radius = sqrt(pow(posx,2)+pow(posy,2));
      double energy = lHit.energy()*mycalib.MeVToMip(layer,radius); // if (energy > 0) std::cout << "sim energy = "<<lHit.energy()<<", reco energy = "<<energy<<std::endl;
      double absweight = absW[layer];
      double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);
      
      h_simzr->Fill(posz,radius);

      //KH - myDigitiser
      bool passTime = myDigitiser.passTimeCut(type,realtime);
      if (!passTime) continue;
      if (realtime<-10.) continue;
      if (lHit.time()==0.) continue; // remove bad timing hits
      // This is not a rigorous timing cut, but clean some simhits which won't make visible differences in results

      ROOT::Math::XYZPoint lpos = ROOT::Math::XYZPoint(posx,posy,posz);
      double eta = lpos.eta();
      double phi = lpos.phi();
      double dR=DeltaR(etaaxis,phiaxis,eta,phi);

      if (!isScint){ // si
	//DetectorEnum adet = subdet.type;
	unsigned adc = 0;
	//double digiE = 0;
	adc = myDigitiser.adcConverter(energy,type);
	//digiE = myDigitiser.adcToMIP(adc,type);	
	if (adc<5) continue;                // ADC cut on silicon hits
      } else {
	if (energy<0.5) continue;           // MIP cut (>=0.5MIP) on scintilator hits
      }

      // 
      if (layer<=51&&lHit.silayer()>=geomConv.getNumberOfSiLayers(type,radius)) continue;

      nhit++;
      if (layer>51) nhitscinti++;
      //-----
      /*
      if (layer>51)
      std::cout << "MeVToMip: " << layer << " "
		<< radius << " "
		<< lHit.cellid() << " "
		<< eta << ", " << phi << ", " << dR << " " << energy << " "
		<< mycalib.MeVToMip(layer,radius) << std::endl;
      */
      
      if (!isScint && dR<0.3){
	std::map<std::pair<int,int>, float>::iterator it = mymap_simhit.find(std::make_pair(layer,lHit.cellid())); 
	if (it != mymap_simhit.end()){
	  (*it).second += energy;
	}
	else {
	  mymap_simhit.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),energy));
	  mymap_simhit_MeVToMip.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),mycalib.MeVToMip(layer,radius)));
	  mymap_simhit_x.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),posx));
	  mymap_simhit_y.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),posy));
	  mymap_simhit_z.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),posz));
	}
      }
      
      if(dR<0.1){ 
	simhitsumE01+=energy*absweight/1000.;
	if (isScint) simhitBHsumE01+=energy*absweight/1000.;
      }
      if(dR<0.2){
	simhitsumE02+=energy*absweight/1000.;
	if (isScint) simhitBHsumE02+=energy*absweight/1000.;
      }
      if(dR<0.3){
	simhitsumE03+=energy*absweight/1000.;
	if (isScint) simhitBHsumE03+=energy*absweight/1000.;
      }
      if(dR<0.4){
	simhitsumE04+=energy*absweight/1000.;
	if (isScint) simhitBHsumE04+=energy*absweight/1000.;
      }
      if(dR<0.5){
	simhitsumE05+=energy*absweight/1000.;      
	if (isScint) simhitBHsumE05+=energy*absweight/1000.;
      }
      //std::cout << absWeight << " " << absW[layer] << std::endl;

    }

    /////-----
    /*
    bool print_header=true;
    for (std::map<std::pair<int,int>,float>::iterator it=mymap_simhit.begin(); it!=mymap_simhit.end(); ++it){
      std::pair<int,int> key=it->first;

      std::map<std::pair<int,int>, float>::iterator iit = mymap_rechit.find(std::make_pair(key.first,key.second)); 
      std::map<std::pair<int,int>, float>::iterator iiit = mymap_simhit_MeVToMip.find(std::make_pair(key.first,key.second)); 
      std::map<std::pair<int,int>, float>::iterator ix = mymap_simhit_x.find(std::make_pair(key.first,key.second)); 
      std::map<std::pair<int,int>, float>::iterator iy = mymap_simhit_y.find(std::make_pair(key.first,key.second)); 
      std::map<std::pair<int,int>, float>::iterator iz = mymap_simhit_z.find(std::make_pair(key.first,key.second)); 
      double rr = pow(pow(ix->second,2)+pow(iy->second,2),0.5);
      if (it->second>0.){
      if (print_header){
	printf("mymap_simhit: cellid, layer, sim E (mip) [rec E, rec/sim E, MeVToMip, rr]\n");
	print_header=false;
      }
      if (iit == mymap_rechit.end()) printf("mymap_simhit: %6d %10d %10.4f\n",key.first,key.second,it->second);
      else{ 
	printf("mymap_simhit: %6d %10d %10.4f %10.4f %10.4f %10.4f %8.3f\n",
	       key.first,key.second,
	       it->second,iit->second,iit->second/it->second,iiit->second,rr);      
	if (iit->second/it->second<0.8 && it->second>50.){
	  
	  for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	    HGCSSSimHit lHit = (*simhitvec)[iH];
	    if (lHit.layer()==key.first&&lHit.cellid()==key.second){

	      unsigned layer = lHit.layer();
	      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	      DetectorEnum type = subdet.type;
	      isScint = subdet.isScint;
	      std::pair<double,double> xy = lHit.get_xy(subdet,geomConv,shape);
	      double posx = xy.first;//lHit.get_x(cellSize);
	      double posy = xy.second;//lHit.get_y(cellSize);
	      double posz = lHit.get_z();
	      double radius = sqrt(pow(posx,2)+pow(posy,2));
	      double energy = lHit.energy()*mycalib.MeVToMip(layer,radius);
	      double absweight = absW[layer];
	      double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);
	      ROOT::Math::XYZPoint lpos = ROOT::Math::XYZPoint(posx,posy,posz);
	      double eta = lpos.eta();
	      double phi = lpos.phi();
	      double dR=DeltaR(etaaxis,phiaxis,eta,phi);
	      
	      std::cout << "Details: " << layer << " "
			<< radius << " "
			<< lHit.cellid() << " "
			<< eta << ", " << phi << ", " << dR << " " << energy << " "
			<< realtime << " "
			<< posx << " " << posy << " " << posz << " "
			<< mycalib.MeVToMip(layer,radius) << std::endl;

	    } // strange layer#, cell	    	   
	  }   // loop over simhits
	}     // reco/sim<0.8 and sim energy>50mips
      }       // found reco/sim match
      }

      //printf("mymap_simhit.insert: %6d %10d %10.4f\n",key.first,key.second,it->second);
    }
    std::cout << "mymap_simhit.size: " << mymap_simhit.size() << std::endl;
    /////-----
    print_header=true;
    for (std::map<std::pair<int,int>,float>::iterator it=mymap_rechit.begin(); it!=mymap_rechit.end(); ++it){
      std::pair<int,int> key=it->first;
      std::map<std::pair<int,int>, float>::iterator iit = mymap_simhit.find(std::make_pair(key.first,key.second)); 
      if (it->second>0.){
      if (print_header){
	printf("mymap_rechit: cellid, layer, rec E (mip) [sim E, sim/rec E, MeVToMip, rr]\n");
	print_header=false;
      }
      if (iit == mymap_simhit.end()) printf("mymap_rechit: %6d %10d %10.4f\n",key.first,key.second,it->second);
      else printf("mymap_rechit: %6d %10d %10.4f %10.4f %10.4f\n",key.first,key.second,it->second,iit->second,iit->second/it->second);
      }
    }
    std::cout << "mymap_rechit.size: " << mymap_rechit.size() << std::endl;
    //std::cout << "mymap_rechit: " << it->first << " => " << it->second << '\n';
    */

    h_egensim->Fill(simhitsumE03/Egen);
    /*
    printf("simhitsumE01,2,3:   %8.3f, %8.3f, %8.3f\n",simhitsumE01,simhitsumE02,simhitsumE03);
    printf("rechitsumE01,2,3:   %8.3f, %8.3f, %8.3f\n",rechitsumE01,rechitsumE02,rechitsumE03);
    printf(" (w/o noise):       %8.3f, %8.3f, %8.3f\n",rechitsumE01,rechitsumE02,rechitsumE03);
    printf("rec/sim sumE01,2,3: %8.3f, %8.3f, %8.3f\n",rechitsumE01/simhitsumE01,rechitsumE02/simhitsumE02,
	   rechitsumE03/simhitsumE03);
    printf(" (w/o noise):       %8.3f, %8.3f, %8.3f\n",rechitsumEWoNoise01/simhitsumE01,
	   rechitsumEWoNoise02/simhitsumE02,
	   rechitsumEWoNoise03/simhitsumE03);
    printf("simhitBHsumE01,2,3:   %8.3f, %8.3f, %8.3f\n",simhitBHsumE01,simhitBHsumE02,simhitBHsumE03);
    printf("rechitBHsumE01,2,3:   %8.3f, %8.3f, %8.3f\n",rechitBHsumE01,rechitBHsumE02,rechitBHsumE03);
    printf("rec/sim sumE01,2,3: %8.3f, %8.3f, %8.3f\n",rechitBHsumE01/simhitBHsumE01,
	   rechitBHsumE02/simhitBHsumE02,
	   rechitBHsumE03/simhitBHsumE03);
    std::cout << "nhit: "  << nhit << " / " << (*simhitvec).size() << " " << nhitscinti << std::endl;
    */
    //=========

    geomConv.initialiseHistos();


    // ============== for low Erecosum events (ereco/egen < .1)=================


   
    if(debug && ErecoRatio < .1) {
      std::cout<<"... Processing entry: " << ievt << std::endl;
      std::cout<<"The amount of rechits in this entry is: " <<(*rechitvec).size()<<std::endl;
      std::cout<<"Ereco/Egen is:                                   "<<ErecoRatio<<std::endl;
      std::cout<<"Ereco within the cone / total Ereco of the entry "<<recoNotCone<<std::endl;
      std::cout<<"Ereco is: " <<ErecoRatio*Egen<<std::endl;
      h_tailEvents->Fill(ErecoRatio,recoNotCone); // should be roughly 1 to 1 
      h_EtotvsEgen->Fill(Egen,norm); // should be roughly 1 to 1
      h_EtotoverEgen->Fill(norm/Egen); // should be a gauss distribution centered around 1
      for (unsigned iP(0); iP<(*rechitvec).size(); ++iP)
	{
	  HGCSSRecoHit lHit = (*rechitvec)[iP];
	  unsigned layer = lHit.layer();
	  double lenergy = 0;
	  double lenergyNoW = 0;
	  if (!isScint){
	    lenergy = lHit.energy()*absW[layer]/1000.; // Si was weighted too high
	    lenergyNoW = lHit.energy();
	  }
	  if (isScint){
	    lenergy = lHit.energy()*absW[layer]/1000.; //Scint was weighted too low
	    lenergyNoW = lHit.energy();
	  }
	  double leta = lHit.eta();
	  double lphi = lHit.phi();
	  double dR=DeltaR(etagen,phigen,leta,lphi);
	  h_EvsdeltaRtail->Fill(dR,lenergy);
	}
      //std::cout<<" gen vec size is "<<(*genvec).size()<<std::endl;
      //std::cout<<" first gen   pt  "<<ptgen<<" egen  "<<Egen<<" pidgen  "<<pidgen<<" etagen  "<<etagen<<" phi gen "<<phigen<<  "mass gen "<< massgen<<std::endl;
      Egen=0.; // added to later sum all egen 
      double egentemp = 0; // to temp store the energy of the gen particle
      printf("=========================================================================================\n");
      std::cout<<"eta of cone: "<<etagen<< "....phi of cone: "<<phigen<<std::endl;
      for (unsigned iP(0); iP<(*genvec).size(); ++iP){
        egentemp = sqrt((*genvec)[iP].px()/1000.*(*genvec)[iP].px()/1000.+(*genvec)[iP].py()/1000.*(*genvec)[iP].py()/1000.+(*genvec)[iP].pz()/1000.*(*genvec)[iP].pz()/1000.);
	std::cout<<"Gen particle "<<iP<<" is (pdgid) "<<(*genvec)[iP].pdgid()<<" with trackID: "<<(*genvec)[iP].trackID()<<" at eta: "<<(*genvec)[iP].eta()<<" with energy: "<<egentemp<<std::endl;
	Egen=Egen+egentemp; // sum egen
      }
      std::cout<<"Final egen sum is = "<<Egen<<"\n"<<std::endl;
    }
    // =========== end of low events ===========================
   
    //++ievtRec;
   }//loop on entries
  // ---------- Event loop ends ----------

  //std::cout<< "max energy "<< energy_max << std::endl;

  //std::cout<<" size of gen hits "<< (*genvec).size()<< " size of rechits " << (*rechitvec).size()<<std::endl;

  if(debug) std::cout<<"writing files"<<std::endl;

  // store maps

  if (develop){

    map_1->GetZaxis()->SetRangeUser(1,200);
    map_1->Print();
    map_1->SetMinimum(0.5);
    map_1->Write();
    
    for (int ilayer=36;ilayer<=51;ilayer++){
      //sprintf(title,"map_1_%d",ilayer);
      map_1_layer[ilayer]->Print();
      map_1_layer[ilayer]->Write();
    }
  }

  outputFile->cd();
  outputFile->Write();
  outputFile->Close();

  return 0;
  
}//main
