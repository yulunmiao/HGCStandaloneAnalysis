#include "HGCSSGeometryConversion.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSGeometryConversion::HGCSSGeometryConversion(std::string filePath,
						 const HGCSSDetector & aDet){
  width_ = 200;//mm
  if (filePath.find("model1")!=filePath.npos) width_ = 500;
  else if (filePath.find("model2")!=filePath.npos) width_ = 1700*2;
  else if (filePath.find("model3")!=filePath.npos) width_ = 1000;

  cellSize_ = 2.5;

  for (unsigned iS(0); iS<aDet.nSections();++iS){
    bool last = false;
    if (iS==aDet.nSections()-1) last = true;
    detector_.addSubdetector(aDet.getSubDetector(iS),last);
  }

}


HGCSSGeometryConversion::~HGCSSGeometryConversion(){
  std::map<DetectorEnum,std::vector<TH2D *> >::iterator liter =
    2DHistMapE_.begin();
  for (; liter !=2DHistMapE_.end();++liter){
    deleteHistos(liter->second);
  }
  2DHistMapE_.clear();
  std::map<DetectorEnum,std::vector<TH2D *> >::iterator liter =
    2DHistMapTime_.begin();
  for (; liter !=2DHistMapTime_.end();++liter){
    deleteHistos(liter->second);
  }
  2DHistMapTime_.clear();
  std::map<DetectorEnum,std::vector<TH2D *> >::iterator liter =
    2DHistMapZ_.begin();
  for (; liter !=2DHistMapZ_.end();++liter){
    deleteHistos(liter->second);
  }
  2DHistMapZ_.clear();
}


void HGCSSGeometryConversion::deleteHistos(std::vector<TH2D *> & aVec){
  if (aVec.size()!=0){
    for (unsigned iL(0); iL<aVec.size();++iL){
      aVec[iL]->Delete();
    }
  }
}

void HGCSSGeometryConversion::setGranularity(const std::vector<unsigned> & granul){
  granularity_.reserve(granul.size());
  for (unsigned iL(0); iL<granul.size();++iL){
    granularity_.push_back(granul[iL]);
  }
}


 void HGCSSGeometryConversion::initialiseHistos(const bool recreate){

  for (unsigned iS(0); iS<detector_.nSections();++iS){
    std::vector<TH2D *> histVecE;
    resetVector(histVecE,detector_.detName(iS),detector_.subDetector(iS),detector_.nLayers(iS),recreate);
    2DHistMapE_[detector_.detType(iS)]=histVecE;

    std::vector<double> avgvecE;
    avgvecE.resize(detector_.nLayers(iS),0);
    avgMapE_[detector_.detType(iS)]=avgvecE;

    std::vector<TH2D *> histVecTime;
    resetVector(histVecTime,detector_.detName(iS),detector_.subDetector(iS),detector_.nLayers(iS),recreate);
    2DHistMapTime_[detector_.detType(iS)]=histVecTime;
    std::vector<TH2D *> histVecZ;
    resetVector(histVecZ,detector_.detName(iS),detector_.subDetector(iS),detector_.nLayers(iS),recreate);
    2DHistMapZ_[detector_.detType(iS)]=histVecZ;

    std::vector<double> avgvecZ;
    avgvecZ.resize(detector_.nLayers(iS),0);
    avgMapZ_[detector_.detType(iS)]=avgvecZ;
  }
 }

void HGCSSGeometryConversion::fill(const DetectorEnum type,
				   const unsigned newlayer,
				   const double & weightedE,
				   const double & aTime,
				   const double & posx,
				   const double & posy,
				   const double & posz)
{
  2DHistMapE_[type][newlayer]->Fill(posx,posy,weightedE);
  2DHistMapTime_[type][newlayer]->Fill(posx,posy,weightedE*aTime);
  2DHistMapZ_[type][newlayer]->Fill(posx,posy,weightedE*posz);
  avgMapZ_[type][newlayer] += weightedE*posz;
  avgMapE_[type][newlayer] += weightedE;
}

double HGCSSGeometryConversion::getAverageZ(const unsigned layer){
  const HGCSSSubDetector & subdet = detector_.subDetector(layer);
  unsigned newlayer = layer-subdet.layerIdMin;
  double avg = 0;
  if (avgMapE_[subdet.type][newlayer]>0)
    avg =avgMapZ_[subdet.type][newlayer]/avgMapE_[subdet.type][newlayer];
  return avg;
}

TH2F * & HGCSSGeometryConversion::get2DHist(const unsigned layer,std::string name){
  const HGCSSSubDetector & subdet = detector_.subDetector(layer);
  unsigned newlayer = layer-subdet.layerIdMin;
  if (name == "E") return 2DHistMapE_[subdet.type][newlayer];
  else if (name == "Time") return 2DHistMapTime_[subdet.type][newlayer];
  else if (name == "Z") return 2DHistMapZ_[subdet.type][newlayer];
  else {
    std::cerr << " ERROR !! Unknown histogram name. Exiting..." << std::endl;
    exit(1);
  }
}

unsigned HGCSSGeometryConversion::getGranularity(const unsigned aLayer, const HGCSSSubDetector & adet){
  unsigned idx = adet.layerIdMin+aLayer;
  return granularity_[idx];
}

void HGCSSGeometryConversion::resetVector(std::vector<TH2D *> & aVec,
					  std::string aString,
					  const HGCSSSubDetector & aDet,
					  const unsigned nLayers,
					  bool recreate=false)
{
  if (aVec.size()!=0){
    for (unsigned iL(0); iL<aVec.size();++iL){
      aVec[iL]->Reset();
      if (recreate) aVec[iL]->Delete();
    }
  }
  if (aVec.size()==0 || recreate){
    if (nLayers > 0){
      aVec.resize(nLayers,0);
      std::cout << " -- Creating " << nLayers << " 2D histograms for " << aString 
	//<< " with " << nBins << " bins between " << min << " and " << max 
		<< std::endl;
      for (unsigned iL(0); iL<nLayers;++iL){
	std::ostringstream lname;
	lname << "EmipHits_" << aString << "_" << iL ;
	double newcellsize = cellSize_*getGranularity(iL,adet);
	//take smallest pair integer to be sure to fit in the geometry
	//even if small dead area introduced at the edge
	unsigned nBins = static_cast<unsigned>(width_*1./(newcellsize*2.))*2;
	double min = -nBins*newcellsize/2.;
	double max = nBins*newcellsize/2.;
	aVec[iL] = new TH2D(lname.str().c_str(),";x(mm);y(mm)",nBins,min,max,nBins,min,max);
      }
    }
  }
}





