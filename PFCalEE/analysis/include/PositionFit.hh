#ifndef PositionFit_hh
#define PositionFit_hh

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSPUenergy.hh"
#include "HGCSSGeometryConversion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

class PositionFit{

public:

  PositionFit(const unsigned nSR,const double & residualMax, const unsigned nLayers, const unsigned nSiLayers, const unsigned debug=0);
  ~PositionFit(){};


  void initialise(TFile *outputFile, 
		  const std::string outFolder, 
		  const HGCSSGeometryConversion & geomConv, 
		  const HGCSSPUenergy & puDensity);

  void initialisePositionHistograms();
  void initialiseFitHistograms();

  bool getZpositions();
  void getZpositions(TTree *aSimTree,
		     const unsigned nEvts);

  void getInitialPositions(TTree *simTree, 
			   TTree *recoTree,
			   const unsigned nEvts);
  
  void getGlobalMaximum(std::vector<HGCSSRecoHit> *rechitvec,double & phimax,double & etamax);

  void getTruthPosition(std::vector<HGCSSGenParticle> *genvec,std::vector<ROOT::Math::XYPoint> & truthPos);

  void getMaximumCell(std::vector<HGCSSRecoHit> *rechitvec,const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax);

  void getEnergyWeightedPosition(std::vector<HGCSSRecoHit> *rechitvec,const unsigned nPU, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<ROOT::Math::XYPoint> & recoPos,std::vector<unsigned> & nHits,const bool puSubtracted=true);

  void fillErrorMatrix(const std::vector<ROOT::Math::XYPoint> & recoPos,const std::vector<ROOT::Math::XYPoint> & truthPos, const std::vector<unsigned> & nHits);

  void finaliseErrorMatrix();
  bool fillMatrixFromFile();
  void fillCorrelationMatrix();

  bool getPositionFromFile(const unsigned ievt,
			   std::vector<unsigned> & layerId,
			   std::vector<double> & posx,
			   std::vector<double> & posy,
			   std::vector<double> & posz,
			   std::vector<double> & posxtruth,
			   std::vector<double> & posytruth,
			   bool cutOutliers=false);

  bool performLeastSquareFit(TTree *recoTree, const unsigned nEvts);

  bool fitEvent(const unsigned ievt,
		unsigned & nInvalidFits,
		std::ofstream & fout,
		const bool cutOutliers=false);

  inline void setOutputFile(TFile *outputFile){
    outputFile_ = outputFile;
  };

  inline unsigned nSR() const{
    return nSR_;
  };

  inline double residualMax() const{
    return residualMax_;
  }


private:
  PositionFit(){};

  unsigned nSR_;
  double residualMax_;
  double chi2ndfmax_;
  unsigned nLayers_;
  unsigned nSiLayers_;
  double cellSize_;
  unsigned debug_;

  HGCSSGeometryConversion geomConv_;
  HGCSSPUenergy puDensity_;

  std::vector<double> avgZ_;

  //initialisation for error matrix
  std::vector<double> mean_[2];//sum residuals for x and y
  std::vector<std::vector<double> > sigma_[2];//sum square
  std::vector<unsigned> nL_mean_;//number of valid layers
  std::vector<std::vector<unsigned> > nL_sigma_;
  TMatrixD matrix_;
  TMatrixD corrMatrix_;

  //path for saving data files
  std::string outFolder_;
  TFile *outputFile_;

  std::vector<TH2F *> p_genxy;
  std::vector<TH2F *> p_recoxy;

  TH2F *p_hitPuContrib;

  TH1F *p_residuals_x;
  TH1F *p_residuals_y;
  TH2F *p_etavsphi_max;
  TH1F *p_nLayersFit;
  TH2F *p_recoXvsLayer;
  TH2F *p_recoYvsLayer;
  TH2F *p_recoZvsLayer;
  TH2F *p_truthXvsLayer;
  TH2F *p_truthYvsLayer;
  TH2F *p_fitXvsLayer;
  TH2F *p_fitYvsLayer;
  TH2D *p_errorMatrix;
  TH2D *p_corrMatrix;
  TH1F *p_chi2[2];
  TH1F *p_chi2overNDF[2];
  TH1F *p_impactX[2];
  TH1F *p_impactY[2];
  TH1F *p_tanAngleX[2];
  TH1F *p_tanAngleY[2];
  TH1F *p_positionReso;
  TH1F *p_angularReso;
  TH1F *p_impactX_residual;
  TH1F *p_impactY_residual;
  TH1F *p_tanAngleX_residual;
  TH1F *p_tanAngleY_residual;

};//class

#endif
