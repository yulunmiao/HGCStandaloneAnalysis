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

  PositionFit(const unsigned nSR,const double & residualMax, const unsigned nLayers, const unsigned nSiLayers,const bool applyPuMixFix,const unsigned debug=0);
  ~PositionFit(){};

  std::pair<unsigned, std::pair<double,double> > findMajorityValue(std::vector<std::pair<double,double> > & values) const;

  void initialise(TFile* outputFile,
		  const std::string outputDir, 
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
  
  bool getGlobalMaximum(const unsigned ievt, 
			const unsigned nVtx, 
			std::vector<HGCSSRecoHit> *rechitvec, 
			const ROOT::Math::XYZVector & truthPos0, 
			double & phimax,double & etamax);

  bool getTruthPosition(std::vector<HGCSSGenParticle> *genvec,std::vector<ROOT::Math::XYPoint> & truthPos);

  void getMaximumCellFromGeom(const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax);

  void getMaximumCell(std::vector<HGCSSRecoHit> *rechitvec,const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax);

  void getEnergyWeightedPosition(std::vector<HGCSSRecoHit> *rechitvec,const unsigned nPU, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<ROOT::Math::XYPoint> & recoPos,std::vector<unsigned> & nHits,std::vector<double> & puE, const bool puSubtracted=true);

  void getPuContribution(std::vector<HGCSSRecoHit> *rechitvec, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<double> & puE);

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

  bool initialiseLeastSquareFit();
  //return 1 if no input file or <3 layers
  //return 2 if chi2/ndf>chi2ndfmax_
  //return 0 if success
  unsigned performLeastSquareFit(const unsigned ievt,
			     std::vector<ROOT::Math::XYZVector> & eventPos);
  void finaliseFit();

  //return 1 if no input file or <3 layers
  //return 2 if chi2/ndf>chi2ndfmax_
  //return 0 if success
  unsigned fitEvent(const unsigned ievt,
		    std::vector<ROOT::Math::XYZVector> & eventPos,
		    const bool cutOutliers=false);

  inline void setOutputFile(TFile * outputFile){
    outputFile_ = outputFile;
    outputFile_->mkdir(outputDir_.c_str());
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
  bool useMeanPU_;
  bool fixForPuMixBug_;

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
  
  unsigned nInvalidFits_;
  unsigned nFailedFitsAfterCut_;
  std::ofstream fout_;

  //path for saving data files
  std::string outFolder_;
  std::string outputDir_;
  TFile *outputFile_;

  TH1F *p_nGenParticles;

  std::vector<TH2F *> p_genxy;
  std::vector<TH2F *> p_recoxy;

  TH1F *p_numberOfMaxTried;
  TH1F *p_dRMaxTruth;
  TH2F *p_hitEventPuContrib;
  TH2F *p_hitMeanPuContrib;
  TH1F *p_diffPuContrib;

  TH1F *p_residuals_x;
  TH1F *p_residuals_y;
  TH2F *p_etavsphi;
  TH2F *p_etavsphi_max;
  TH2F *p_etavsphi_truth;

  TH2F *p_yvsx_max;
  TH2F *p_yvsx_truth;

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
