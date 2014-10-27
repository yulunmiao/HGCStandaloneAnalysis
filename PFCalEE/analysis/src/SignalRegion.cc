#include "SignalRegion.hh"

#include "TCanvas.h"
#include "TProfile2D.h"
#include "TStyle.h"

SignalRegion::SignalRegion(const std::string inputFolder, double cellSize, unsigned nLayers, unsigned nevt){

    nevt_ = nevt;
    cellSize_ = cellSize;
    nLayers_ = nLayers;

    double xpos(0),ypos(0),zpos(0),xangle(0),yangle(0);    
    int layerIndex(0),eventIndex(0);

    std::ifstream fzpos;
    std::ostringstream finname;
    finname << inputFolder << "_zPositions.dat";
    fzpos.open(finname.str());
    if (!fzpos.is_open()){
      std::cout << " Cannot open input file " << finname.str() << "! Refilling now..." << std::endl;
    }
    for(unsigned iL(0);iL<nLayers_;iL++){
        fzpos >> layerIndex >> zpos;
        zPos_.push_back(zpos);
    }

    std::ifstream fxypos;
    finname.str("");
    finname << inputFolder << "_accuratePos.dat";
    fxypos.open(finname.str());
    if (!fxypos.is_open()){
      std::cout << " Cannot open input file " << finname.str() << "! Refilling now..." << std::endl;
    }
    for(unsigned ievt(0);ievt<nevt_;ievt++){ 
        std::string fitMatrix[4];
        fxypos >> eventIndex >> xpos >> fitMatrix[0] >> xangle >> fitMatrix[1] >> ypos >> fitMatrix[2] >> yangle >> fitMatrix[3];
        std::vector<ROOT::Math::XYZVector> tmpXYZ;
        for(unsigned iL(0);iL<nLayers_;iL++){
           ROOT::Math::XYZVector position( xpos + xangle*zPos_[iL], ypos + yangle*zPos_[iL], zPos_[iL]);
           tmpXYZ.push_back(position);
        }
        accuratePos_.push_back(tmpXYZ);
    }

}


SignalRegion::~SignalRegion(){
};



void SignalRegion::initialise(std::vector<HGCSSRecoHit> *rechitvec, unsigned ievt){

    if(ievt >= nevt_) return;
    std::vector<ROOT::Math::XYZVector> eventPos = accuratePos_[ievt];
    if(eventPos.size()!=nLayers_) return; 

    const double halfCell = 0.5*cellSize_;

    std::vector<double> accurateX;
    std::vector<double> accurateY;
    std::vector<double> signalSR0, signalSR1, signalSR2, signalSR3, signalSR4;
    signalSR0.resize(eventPos.size(),0);
    signalSR1.resize(eventPos.size(),0);
    signalSR2.resize(eventPos.size(),0);
    signalSR3.resize(eventPos.size(),0);
    signalSR4.resize(eventPos.size(),0);

    // set the accurate coordinate for each layer
    for(unsigned iL(0); iL< eventPos.size(); iL++){
        accurateX.push_back(eventPos[iL].x());
        accurateY.push_back(eventPos[iL].y());
    }

    // Define different signal region and sum over energy
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
        const HGCSSRecoHit & lHit = (*rechitvec)[iH];

        unsigned layer = lHit.layer();
        if (layer >= nLayers_) {
            continue;
        }
        double posx = lHit.get_x();
        double posy = lHit.get_y();
        double posz = lHit.get_z();
        double radius = sqrt(posx*posx+posy*posy);
        double energy = lHit.energy();

        //SR0
        if(fabs(posx + halfCell-accurateX[layer])< halfCell && fabs(posy+ halfCell-accurateY[layer])< halfCell) signalSR0[layer] += energy;
        //SR1
        if(fabs(posx + halfCell-accurateX[layer]) < 2*halfCell && fabs(posy+ halfCell-accurateY[layer]) < 2*halfCell) signalSR1[layer] += energy; 
        //SR2
        if(fabs(posx + halfCell-accurateX[layer]) < 3*halfCell && fabs(posy+ halfCell-accurateY[layer]) < 3*halfCell) signalSR2[layer] += energy;
        //SR3
        if(fabs(posx + halfCell-accurateX[layer]) < 4*halfCell && fabs(posy+ halfCell-accurateY[layer]) < 4*halfCell) signalSR3[layer] += energy;
        //SR4
        if(fabs(posx + halfCell-accurateX[layer]) < 5*halfCell && fabs(posy+ halfCell-accurateY[layer]) < 5*halfCell) signalSR4[layer] += energy;

     }

     energySR0_.push_back(signalSR0); 
     energySR1_.push_back(signalSR1); 
     energySR2_.push_back(signalSR2); 
     energySR3_.push_back(signalSR3); 
     energySR4_.push_back(signalSR4); 

}








    
