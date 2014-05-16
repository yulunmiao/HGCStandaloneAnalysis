#include "HGCSSGeometryConversion.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSGeometryConversion::HGCSSGeometryConversion(std::string filePath){
  width_ = 200;//mm
  if (filePath.find("model1")!=filePath.npos) width_ = 500;
  else if (filePath.find("model2")!=filePath.npos) width_ = 1700*2;
  else if (filePath.find("model3")!=filePath.npos) width_ = 1000;

  cellSize_ = 2.5;
  
}


HGCSSGeometryConversion::~HGCSSGeometryConversion(){
  std::map<DetectorEnum,std::vector<TH2D *> >::iterator liter =
    HistMapE_.begin();
  for (; liter !=HistMapE_.end();++liter){
    deleteHistos(liter->second);
  }
  HistMapE_.clear();
  liter = HistMapTime_.begin();
  for (; liter !=HistMapTime_.end();++liter){
    deleteHistos(liter->second);
  }
  HistMapTime_.clear();
  liter = HistMapZ_.begin();
  for (; liter !=HistMapZ_.end();++liter){
    deleteHistos(liter->second);
  }
  HistMapZ_.clear();
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

   for (unsigned iS(0); iS<theDetector().nSections();++iS){
    std::vector<TH2D *> histVecE;
    resetVector(histVecE,theDetector().detName(iS),theDetector().subDetector(iS),theDetector().nLayers(iS),recreate);
    HistMapE_[theDetector().detType(iS)]=histVecE;

    std::vector<double> avgvecE;
    avgvecE.resize(theDetector().nLayers(iS),0);
    avgMapE_[theDetector().detType(iS)]=avgvecE;

    std::vector<TH2D *> histVecTime;
    resetVector(histVecTime,theDetector().detName(iS),theDetector().subDetector(iS),theDetector().nLayers(iS),recreate);
    HistMapTime_[theDetector().detType(iS)]=histVecTime;
    std::vector<TH2D *> histVecZ;
    resetVector(histVecZ,theDetector().detName(iS),theDetector().subDetector(iS),theDetector().nLayers(iS),recreate);
    HistMapZ_[theDetector().detType(iS)]=histVecZ;

    std::vector<double> avgvecZ;
    avgvecZ.resize(theDetector().nLayers(iS),0);
    avgMapZ_[theDetector().detType(iS)]=avgvecZ;
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
  HistMapE_[type][newlayer]->Fill(posx,posy,weightedE);
  HistMapTime_[type][newlayer]->Fill(posx,posy,weightedE*aTime);
  HistMapZ_[type][newlayer]->Fill(posx,posy,weightedE*posz);
  avgMapZ_[type][newlayer] += weightedE*posz;
  avgMapE_[type][newlayer] += weightedE;
}

double HGCSSGeometryConversion::getAverageZ(const unsigned layer){
  const HGCSSSubDetector & subdet = theDetector().subDetector(layer);
  unsigned newlayer = layer-subdet.layerIdMin;
  double avg = 0;
  if (avgMapE_[subdet.type][newlayer]>0)
    avg =avgMapZ_[subdet.type][newlayer]/avgMapE_[subdet.type][newlayer];
  return avg;
}

TH2D * HGCSSGeometryConversion::get2DHist(const unsigned layer,std::string name){
  const HGCSSSubDetector & subdet = theDetector().subDetector(layer);
  unsigned newlayer = layer-subdet.layerIdMin;
  if (name == "E") return HistMapE_[subdet.type][newlayer];
  else if (name == "Time") return HistMapTime_[subdet.type][newlayer];
  else if (name == "Z") return HistMapZ_[subdet.type][newlayer];
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
					  bool recreate)
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
	double newcellsize = cellSize_*getGranularity(iL,aDet);
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





