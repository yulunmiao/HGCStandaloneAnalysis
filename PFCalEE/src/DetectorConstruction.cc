#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "HGCSSSimHit.hh"

#include <boost/algorithm/string.hpp>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalConstants.hh"


using namespace std;

//
DetectorConstruction::DetectorConstruction(G4int ver, G4int mod,
					   G4int shape,
					   std::string absThickW,
					   std::string absThickPb,
					   std::string dropLayer,
                                           int coarseGranularity) :
  version_(ver), model_(mod), shape_(shape), addPrePCB_(false), m_coarseGranularity(coarseGranularity)
{
  doHF_ = false;

  lastEElayer_ = 9999;
  firstHFlayer_ = 9999;
  firstMixedlayer_ = 9999;
  firstScintlayer_ = 9999;
  firstCoarseScintlayer_ = 9999;


  m_minEta0 = 1.3;
  m_maxEta0 = 3.0;
  // min Eta per layer, 0-27 for EE, 28-35 FH, 36-51 BH
  //double minEta[52] =
  double minEta_tmp[52] = { 1.461, 1.464, 1.463, 1.466, 1.465, 1.468, 1.467, 1.469, 1.469, 1.471,
			1.471, 1.473, 1.472, 1.475, 1.474, 1.477, 1.476, 1.478, 1.478, 1.480,
			1.479, 1.482, 1.481, 1.483, 1.483, 1.485, 1.484, 1.487,
			1.487, 1.490, 1.494, 1.497, 1.500, 1.503, 1.506, 1.493,
			1.474, 1.455, 1.437, 1.420, 1.403, 1.376, 1.351, 1.327, 1.306, 1.316,
			1.332, 1.348, 1.364, 1.379, 1.395, 1.410 };
  //std::cout << "DetectorConstruction constructor:" << std::endl;
  for (int i=0; i<52; i++){
    minEta.push_back(minEta_tmp[i]);
    //std::cout << minEta[i] << std::endl;
  }

  SetWThick(absThickW);
  SetPbThick(absThickPb);
  SetDropLayers(dropLayer);

  //radiation lengths: cf. http://pdg.lbl.gov/2012/AtomicNuclearProperties/
  //W 3.504 mm
  //Pb 5.612 mm
  //Cu 14.36 mm
  switch(version_)
    {
      //cf. http://arxiv.org/abs/0805.4833
    case v_CALICE:
      {
	G4cout << "[DetectorConstruction] starting v_CALICE (10x0.4+10x0.8+10x1.2)X_0 with Tungsten" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(0.4*3.504*mm);lEle.push_back("W");
	lThick.push_back(0.525*mm);lEle.push_back("Si");
	lThick.push_back(1.0*mm);lEle.push_back("PCB");
	lThick.push_back(2.5*mm);lEle.push_back("Air");

	for(unsigned i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 0.8*3.504*mm;
	for(unsigned i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 1.2*3.504*mm;
	for(unsigned i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	break;
      }

    case v_HGCALEE_TB: case v_HGCALEE_TB_gap0: case v_HGCALEE_TB_allW: case v_HGCALEE_TB_samedEdx:
      {
	G4cout << "[DetectorConstruction] starting v_HGCAL for testbeam"<< G4endl;
	G4double airThick = 2*mm;
	if(version_ == v_HGCALEE_TB_gap0) airThick = 0*mm;
	G4double pcbThick = 2*mm;

	std::vector<G4double> lThickL;
	std::vector<std::string> lEleL;
	lThickL.push_back(0.5*mm);lEleL.push_back("Cu");
	lThickL.push_back(0.5*mm);lEleL.push_back("CFMix");
	lThickL.push_back(2.*mm);lEleL.push_back("W");
	lThickL.push_back(0.5*mm);lEleL.push_back("CFMix");
	lThickL.push_back(0.5*mm);lEleL.push_back("Cu");
	lThickL.push_back(airThick);lEleL.push_back("Air");
	lThickL.push_back(pcbThick);lEleL.push_back("PCB");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");

	if (version_ == v_HGCALEE_TB_samedEdx) lThickL[2] = 2.8*mm;

	std::vector<G4double> lThickR;
	std::vector<std::string> lEleR;
	if(version_ != v_HGCALEE_TB_allW){
	  lThickR.push_back(0.6*mm);lEleR.push_back("WCu");
	  lThickR.push_back(6*mm);lEleR.push_back("Cu");
	  lThickR.push_back(0.6*mm);lEleR.push_back("WCu");
	} else {
	  lThickR.push_back(0.5*mm);lEleR.push_back("Cu");
	  lThickR.push_back(2.*mm);lEleR.push_back("W");
	  lThickR.push_back(0.5*mm);lEleR.push_back("Cu");
	}
	lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	lThickR.push_back(pcbThick);lEleR.push_back("PCB");
	lThickR.push_back(airThick);lEleR.push_back("Air");
	if (version_ == v_HGCALEE_TB_samedEdx) {
	  lThickR[0] = 0*mm;
	  lThickR[2] = 0*mm;
	}

	for(unsigned i=0; i<5; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}
	lThickL[2] = 2.8*mm;
	if (version_ == v_HGCALEE_TB_samedEdx) lThickL[2] = 3.7*mm;
	if(version_ != v_HGCALEE_TB_allW){
	  lThickR[0] = 1.2*mm;
	  lThickR[2] = 1.2*mm;
	} else {
	  lThickR[1] = 2.8*mm;
	}
 	if (version_ == v_HGCALEE_TB_samedEdx) {
	  lThickR[0] = 0.5*mm;
	  lThickR[2] = 0.5*mm;
	}
	for(unsigned i=0; i<5; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}
	lThickL[2] = 4.2*mm;
	if (version_ == v_HGCALEE_TB_samedEdx) lThickL[2] = 5.3*mm;
	if(version_ != v_HGCALEE_TB_allW){
	  lThickR[0] = 2.2*mm;
	  lThickR[2] = 2.2*mm;
	} else {
	  lThickR[1] = 4.2*mm;
	}
 	if (version_ == v_HGCALEE_TB_samedEdx) {
	  lThickR[0] = 1.4*mm;
	  lThickR[2] = 1.4*mm;
	}
 	for(unsigned i=0; i<4; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}
	break;
      }//TB setup

    case v_HGCAL_2016TB:
      {
	G4cout << "[DetectorConstruction] starting v_HGCAL for testbeam 2016"<< G4endl;

	G4double absPlateAirGap = 6*mm;

	std::vector<G4double> lThick;
	std::vector<std::string> lEle;

	// One module behind the copper cooling plate with
	// tungsten absorbing plates in front.
	//
	for (int i=0; i<3; i++) {
	  lThick.push_back(4.2*mm);        lEle.push_back("W");
	  lThick.push_back(absPlateAirGap);lEle.push_back("Air");
	}
	for (int i=0; i<2; i++) {
	  lThick.push_back(2.1*mm);        lEle.push_back("W");
	  lThick.push_back(absPlateAirGap);lEle.push_back("Air");
	}
	lThick.push_back(2.1*mm);        lEle.push_back("W");
	lThick.push_back(3.0*mm);        lEle.push_back("Air");

	lThick.push_back(0.6*mm); lEle.push_back("WCu");
	lThick.push_back(6.0*mm); lEle.push_back("Cu");
	lThick.push_back(0.6*mm); lEle.push_back("WCu");
	lThick.push_back(0.01*mm);lEle.push_back("Air"); // kapton
	lThick.push_back(0.1*mm); lEle.push_back("Si");
	lThick.push_back(0.1*mm); lEle.push_back("Si");
	lThick.push_back(0.12*mm); lEle.push_back("Si");

	m_caloStruct.push_back( SamplingSection(lThick,lEle) );

	break;
      }//TB setup

    case v_HGCALEE_v6: case v_HGCAL_v6: case v_HGCALEE_v7: case v_HGCAL_v7: case v_HGCAL_v7_HF: case v_HGCALEE_v624: case v_HGCALEE_v618: case v_HGCAL_v624: case v_HGCAL_v618:
      {
	G4cout << "[DetectorConstruction] starting v_HGCALEE_v6"<< G4endl;
	G4double airThick = 2*mm;
	G4double pcbThick = 2*mm;
        unsigned Nmodule=4;
	G4double wThick = 2.*mm;
	G4double wcuThick = 0.6*mm;

        if(version_ == v_HGCALEE_v624 || version_ == v_HGCAL_v624){
             Nmodule = 3;
	     wThick = 2.6*mm;
	     wcuThick = 1.*mm;
	}
        else if(version_ == v_HGCALEE_v618 || version_ == v_HGCAL_v618){
             Nmodule =2;
	     wThick = 3.5*mm;
	     wcuThick = 1.7*mm;
	}

	std::vector<G4double> lThickL;
	std::vector<std::string> lEleL;
        lThickL.push_back(100.*mm);lEleL.push_back("NeutMod");
	lThickL.push_back(2.*mm);lEleL.push_back("Al");
	lThickL.push_back(16.*mm);lEleL.push_back("Foam");
	lThickL.push_back(2.*mm);lEleL.push_back("Al");
	lThickL.push_back(10.*mm);lEleL.push_back("Air");
	lThickL.push_back(0.5*mm);lEleL.push_back("CFMix");
	lThickL.push_back(wThick);lEleL.push_back("W");
	lThickL.push_back(0.5*mm);lEleL.push_back("CFMix");
	lThickL.push_back(0.5*mm);lEleL.push_back("Cu");
	lThickL.push_back(airThick);lEleL.push_back("Air");
	lThickL.push_back(pcbThick);lEleL.push_back("PCB");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");

	std::vector<G4double> lThickR;
	std::vector<std::string> lEleR;
	lThickR.push_back(wcuThick);lEleR.push_back("WCu");
	lThickR.push_back(6*mm);lEleR.push_back("Cu");
	lThickR.push_back(wcuThick);lEleR.push_back("WCu");
	lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	lThickR.push_back(pcbThick);lEleR.push_back("PCB");
	lThickR.push_back(airThick);lEleR.push_back("Air");

	m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );

	lThickL.clear();
	lEleL.clear();
	lThickL.push_back(0.5*mm);lEleL.push_back("Cu");
	lThickL.push_back(0.5*mm);lEleL.push_back("CFMix");
	lThickL.push_back(wThick);lEleL.push_back("W");
	lThickL.push_back(0.5*mm);lEleL.push_back("CFMix");
	lThickL.push_back(0.5*mm);lEleL.push_back("Cu");
	lThickL.push_back(airThick);lEleL.push_back("Air");
	lThickL.push_back(pcbThick);lEleL.push_back("PCB");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	for(unsigned i=0; i<Nmodule; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}

        Nmodule=5;
	lThickL[2] = 2.8*mm;
	lThickR[0] = 1.2*mm;
	lThickR[2] = 1.2*mm;
        if(version_ == v_HGCALEE_v624 || version_ == v_HGCAL_v624){
            Nmodule=4;
            lThickL[2] = 3.6*mm;
            lThickR[0] = 1.75*mm;
            lThickR[2] = 1.75*mm;
        }
        else if(version_ == v_HGCALEE_v618 || version_ == v_HGCAL_v618){
            Nmodule=3;
            lThickL[2] = 4.9*mm;
            lThickR[0] = 2.7*mm;
            lThickR[2] = 2.7*mm;
        }
	for(unsigned i=0; i<Nmodule; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}

        Nmodule=4;
	lThickL[2] = 4.2*mm;
	lThickR[0] = 2.2*mm;
	lThickR[2] = 2.2*mm;
        if(version_ == v_HGCALEE_v618 || version_ == v_HGCAL_v618){
            Nmodule=3;
            lThickL[2] = 5.8*mm;
            lThickR[0] = 3.3*mm;
            lThickR[2] = 3.3*mm;
        }
	for(unsigned i=0; i<Nmodule; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  //next line moved to HCAL to get proper sampling...
	  //if (i==3) {lThickR.push_back(0.5*mm);lEleR.push_back("Cu");}
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}

	if(version_==v_HGCAL_v6){
	  //add HCAL
	  buildHGCALFHE(6);
	  buildHGCALBHE(6);
	}
	else if(version_==v_HGCAL_v7 || version_==v_HGCAL_v7_HF){
	  //add HCAL
	  buildHGCALFHE(7);
	  buildHGCALBHE(7);
	  if (version_==v_HGCAL_v7_HF) buildHF();
	}
	else if(version_==v_HGCAL_v624){
	  //add HCAL
	  buildHGCALFHE(624);
	  buildHGCALBHE(6);
	}
	else if(version_==v_HGCAL_v618){
	  //add HCAL
	  buildHGCALFHE(618);
	  buildHGCALBHE(6);
	}

	break;
      }

    case v_testCu:
      {
	G4double cuExtra = 1.5*mm;
	std::vector<G4double> lThickL;
	std::vector<std::string> lEleL;

	lThickL.push_back(1.6*mm);lEleL.push_back("PCB");
	lThickL.push_back(0.5*mm);lEleL.push_back("Air");
	lThickL.push_back(cuExtra);lEleL.push_back("CuExtra");
	lThickL.push_back(2*mm);lEleL.push_back("Air");
	lThickL.push_back(1.6*mm);lEleL.push_back("PCB");

	m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	m_minEta.push_back(m_minEta0);m_maxEta.push_back(m_maxEta0);

	break;
      }

    case v_HGCALEE_v8: case v_HGCAL_v8: case v_HGCALEE_v8_air3: case v_HGCALEE_v8_air4: case v_HGCALEE_v8_Cu: case v_HGCALEE_v8_Cu_12: case v_HGCALEE_v8_neutmod: case v_HGCAL_v8_envelope:
      {

	G4cout << "[DetectorConstruction] starting v_HGCAL(EE)_v8"<< G4endl;
	G4double airThick1 = (version_==v_HGCALEE_v8_air4)? 4*mm :(version_==v_HGCALEE_v8_air3)? 3*mm : 1.5*mm;
	G4double airThick2 = 0*mm;
	G4double cuExtra = 0*mm;
	if (version_==v_HGCALEE_v8_Cu || version_==v_HGCALEE_v8_Cu_12){
	  airThick1 = 0.5*mm;
	  cuExtra = 1.5*mm;
	  airThick2 = 2.0*mm;
	}
	G4double pcbThick = 1.6*mm;
	G4double pbThick = 2.1*mm;
	G4double wcuThick = 1.4*mm;
	G4int nK7 = 13;

	if (version_==v_HGCALEE_v8_Cu) {
	  pbThick = 2.*mm;
	  wcuThick = 1.*mm;
	}

	if (version_==v_HGCALEE_v8_Cu_12) {
	  pbThick = 2.4*mm;
	  nK7 = 11;
	}

	std::vector<G4double> lThickL;
	std::vector<std::string> lEleL;


	//in front of first layer
        if(version_!=v_HGCALEE_v8_neutmod) {
          lThickL.push_back(210.5*mm);lEleL.push_back("Air");
        }else {
          lThickL.push_back(2.*mm);lEleL.push_back("Al");
          lThickL.push_back(36.*mm);lEleL.push_back("Foam");
          lThickL.push_back(2.*mm);lEleL.push_back("Al");
          lThickL.push_back(10*mm);lEleL.push_back("Air");
          lThickL.push_back(157.*mm);lEleL.push_back("NeutMod");
          lThickL.push_back(3.5*mm);lEleL.push_back("Air");
        }

	//cassette structure

	lThickL.push_back(0.3*mm);lEleL.push_back("SSteel");
	lThickL.push_back(pbThick);lEleL.push_back("Pb");
	lThickL.push_back(0.3*mm);lEleL.push_back("SSteel");
	lThickL.push_back(0.1*mm);lEleL.push_back("Cu");
	lThickL.push_back(pcbThick);lEleL.push_back("PCB");

	lThickL.push_back(airThick1);lEleL.push_back("Air");
	lThickL.push_back(cuExtra);lEleL.push_back("CuExtra");
	lThickL.push_back(airThick2);lEleL.push_back("Air");
	lThickL.push_back(pcbThick);lEleL.push_back("PCB");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");

	std::vector<G4double> lThickR;
	std::vector<std::string> lEleR;
	lThickR.push_back(wcuThick);lEleR.push_back("WCu");
	lThickR.push_back(6*mm);lEleR.push_back("Cu");
	lThickR.push_back(wcuThick);lEleR.push_back("WCu");
	lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	lThickR.push_back(0.1*mm);lEleR.push_back("Si");

	//first cassette
	m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	m_minEta.push_back(minEta[0]);m_maxEta.push_back(m_maxEta0);
	m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	m_minEta.push_back(minEta[1]);m_maxEta.push_back(m_maxEta0);

	lThickL.clear();
	lEleL.clear();
	lThickL.push_back(pcbThick);lEleL.push_back("PCB");
	lThickL.push_back(airThick2);lEleL.push_back("Air");
	lThickL.push_back(cuExtra);lEleL.push_back("CuExtra");
	lThickL.push_back(airThick1);lEleL.push_back("Air");
	lThickL.push_back(pcbThick);lEleL.push_back("PCB");
	lThickL.push_back(0.1*mm);lEleL.push_back("Cu");
	lThickL.push_back(0.3*mm);lEleL.push_back("SSteel");
	lThickL.push_back(pbThick);lEleL.push_back("Pb");
	lThickL.push_back(0.3*mm);lEleL.push_back("SSteel");
	lThickL.push_back(0.3*mm);lEleL.push_back("SSteel");
	lThickL.push_back(pbThick);lEleL.push_back("Pb");
	lThickL.push_back(0.3*mm);lEleL.push_back("SSteel");
	lThickL.push_back(0.1*mm);lEleL.push_back("Cu");
	lThickL.push_back(pcbThick);lEleL.push_back("PCB");

	lThickL.push_back(airThick1);lEleL.push_back("Air");
	lThickL.push_back(cuExtra);lEleL.push_back("CuExtra");
	lThickL.push_back(airThick2);lEleL.push_back("Air");
	lThickL.push_back(pcbThick);lEleL.push_back("PCB");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	lThickL.push_back(0.1*mm);lEleL.push_back("Si");



	for(G4int i=0; i<nK7; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_minEta.push_back(minEta[2+2*i]);m_maxEta.push_back(m_maxEta0); // minEta index: 2,4,...,26
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	  m_minEta.push_back(minEta[2+2*i+1]);m_maxEta.push_back(m_maxEta0); // minEta index: 3,5,...,27
	}

	lastEElayer_ = m_caloStruct.size()-1;

	if(version_==v_HGCAL_v8 || version_==v_HGCAL_v8_envelope){
	  //add HCAL
	  //FH = FH+BH silicon version = 24 layers
	  buildHGCALFHE(8);
	  //BH = FH+BH scintillator version = 16 layers
	  buildHGCALBHE(8);
	} else {
          //add the back cover which would otherwise be added in buildHGCALFHE
          lThickL.clear(); lEleL.clear();
          lThickL.push_back(1*mm);  lEleL.push_back("Cu");
          lThickL.push_back(45*mm); lEleL.push_back("SSteel");
          m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_minEta.push_back( m_minEta[m_minEta.size()-1] ); m_maxEta.push_back(m_maxEta0);
        }

	break;
      }
      
      //Jan2021 update
    case v_HGCALEE_v9: case v_HGCAL_v9: case v_HGCALEE_v10: case v_HGCAL_v10:
      {

        G4cout << "[DetectorConstruction] starting v_HGCAL(EE)_v9/10"<< G4endl;
	
        G4double z0(2980);

        //hard-point extension
        G4double hpAirThick = 0.1*mm;

        //lead-ss sandwich
        G4double cuThick1 = 0.1*mm;
        G4double inoxThick1 = 0.3*mm;
        G4double leadThick[3] = {3.1*mm,5.2*mm,8.22*mm}; // cassetes 1; 2-9; 10-13;
        if(version_==v_HGCAL_v10 || version_==v_HGCALEE_v10) {
          leadThick[1]=5.70045579;
          leadThick[2]=5.70045579;
        }
        G4double inoxThick2 = 0.3*mm;
        G4double cuThick2 = 0.1*mm;
        G4double airThick = 0.2*mm;

        //motherboard
        G4double mbPCBThick(1.6*mm);
        G4double mbAirThick(0.5*(5.25-mbPCBThick));

        //module
        G4double modPCBThick(1.6*mm);
        G4double modAirThick1(0.4*mm);
        G4double modSiThick(0.1*mm); //x3 below
        G4double modAirThick2(0.3*mm);
        G4double modWCuThick(1.4*mm);
        if(version_==v_HGCAL_v10 || version_==v_HGCALEE_v10) {
          modWCuThick=1.521929934;
        }
        //cooling plate
        G4double coolingCuThick(6.05*mm);

        std::vector<G4double> a_lThick, b_lThick;
        std::vector<std::string> a_lEle, b_lEle;
        
        //cassette structure
        G4double totalThick(0.);
        for(int i=0;i<13; i++) {
          
          a_lThick.clear(); b_lThick.clear();
          a_lEle.clear();   b_lEle.clear();

          //in front of first layer only
          if(i==0) {

            //air just to get first layer in the correct position
            a_lThick.push_back(222.5*mm); a_lEle.push_back("Air");

            /*Chris said to ignore this for the moment
            a_lThick.push_back(2.*mm);   a_lEle.push_back("Al");
            a_lThick.push_back(36.*mm);  a_lEle.push_back("Foam");
            a_lThick.push_back(2.*mm);   a_lEle.push_back("Al");
            a_lThick.push_back(10*mm);   a_lEle.push_back("Air");
            a_lThick.push_back(157.*mm); a_lEle.push_back("NeutMod");
            a_lThick.push_back(3.5*mm);  a_lEle.push_back("Air");
            */
          }

          //A-side of the cassette
          //hard-point extension
          a_lThick.push_back(hpAirThick); a_lEle.push_back("Air");
          
          //lead-ss sandwich
          int Pbidx(0);
          if(i>0 && i<9) Pbidx=1;
          if(i>=9)       Pbidx=2;

          a_lThick.push_back(cuThick1);         a_lEle.push_back("Cu");
          a_lThick.push_back(inoxThick1);       a_lEle.push_back("SSteel");
          a_lThick.push_back(leadThick[Pbidx]); a_lEle.push_back("Pb");
          a_lThick.push_back(inoxThick2);       a_lEle.push_back("SSteel");
          a_lThick.push_back(cuThick2);         a_lEle.push_back("Cu");
          a_lThick.push_back(airThick);         a_lEle.push_back("Air");
          
          //motherboard
          a_lThick.push_back(mbAirThick);    a_lEle.push_back("Air");
          a_lThick.push_back(mbPCBThick);    a_lEle.push_back("PCB");
          a_lThick.push_back(mbAirThick);    a_lEle.push_back("Air");

          //module
          a_lThick.push_back(modPCBThick);   a_lEle.push_back("PCB");
          a_lThick.push_back(modAirThick1);  a_lEle.push_back("Air");
          for(int j=0; j<3; j++){
            a_lThick.push_back(modSiThick);  a_lEle.push_back("Si");
          }
          a_lThick.push_back(modAirThick2);  a_lEle.push_back("Air");
          a_lThick.push_back(modWCuThick);   a_lEle.push_back("WCu");
   
          //B-side of the cassette
          //cooling plate
          b_lThick.push_back(coolingCuThick); b_lEle.push_back("Cu");

          //"inverted" module
          b_lThick.push_back(modWCuThick);   b_lEle.push_back("WCu");
          b_lThick.push_back(modAirThick2);  b_lEle.push_back("Air");
          for(int j=0; j<3; j++){
            b_lThick.push_back(modSiThick);  b_lEle.push_back("Si");
          }
          b_lThick.push_back(modAirThick1);  b_lEle.push_back("Air");
          b_lThick.push_back(modPCBThick);   b_lEle.push_back("PCB");

          //motherboard
          b_lThick.push_back(mbAirThick);    b_lEle.push_back("Air");
          b_lThick.push_back(mbPCBThick);    b_lEle.push_back("PCB");
          b_lThick.push_back(mbAirThick);    b_lEle.push_back("Air");

          //hard-point extension
          b_lThick.push_back(hpAirThick); b_lEle.push_back("Air");

          if(i==12 && (version_==v_HGCALEE_v9 || version_==v_HGCALEE_v10)) {
            //add the back cover which would otherwise be added in buildHGCALFHE
            b_lThick.push_back(1*mm);  b_lEle.push_back("Cu");
            b_lThick.push_back(45*mm); b_lEle.push_back("SSteel");
          }

          m_caloStruct.push_back( SamplingSection(a_lThick,a_lEle) );
          totalThick += m_caloStruct[m_caloStruct.size()-1].Total_thick;
          
          float z=z0+totalThick;
          float rinner=getRealisticRho(z,"inner");
          float router=getRealisticRho(z,"outer");
          m_minEta.push_back(getEtaFromRZ(router,z));
          m_maxEta.push_back(getEtaFromRZ(rinner,z));

          m_caloStruct.push_back( SamplingSection(b_lThick,b_lEle) );
          totalThick += m_caloStruct[m_caloStruct.size()-1].Total_thick;
          z=z0+totalThick;
          rinner=getRealisticRho(z,"inner");
          router=getRealisticRho(z,"outer");
          m_minEta.push_back(getEtaFromRZ(router,z));
          m_maxEta.push_back(getEtaFromRZ(rinner,z));
        }

	lastEElayer_ = m_caloStruct.size()-1;

        //add HCAL if full detector
        if(version_==v_HGCAL_v9 || version_==v_HGCAL_v10){

          //FH = FH+BH silicon version
          buildHGCALFHE(9);
          //BH = FH+BH scintillator version
          buildHGCALBHE(9);
        }

        break;
      }

    case v_HGCALEE_v5: case v_HGCALEE_v5_gap4: case v_HGCAL_v5: case v_HGCAL_v5_gap4:
      {
	G4cout << "[DetectorConstruction] starting v_HCALEE_v5"<< G4endl;
	G4double airThick = 2*mm;
	if(version_==v_HGCALEE_v5_gap4 || version_==v_HGCAL_v5_gap4) airThick = 4*mm;

	G4double pcbThick = 1.2*mm;

	//first and last layers
	std::vector<G4double> lThick1;
	std::vector<std::string> lEle1;
	lThick1.push_back(180.*mm);lEle1.push_back("NeutMod");
	lThick1.push_back(2.*mm);lEle1.push_back("Al");
	lThick1.push_back(26.*mm);lEle1.push_back("Foam");
	lThick1.push_back(2.*mm);lEle1.push_back("Al");
	lThick1.push_back(0*mm);lEle1.push_back("Cu");
	lThick1.push_back(0*mm);lEle1.push_back("W");
	lThick1.push_back(0.5*mm);lEle1.push_back("Cu");
	lThick1.push_back(airThick);lEle1.push_back("Air");
	lThick1.push_back(pcbThick);lEle1.push_back("PCB");
	lThick1.push_back(0.1*mm);lEle1.push_back("Si");
	lThick1.push_back(0.1*mm);lEle1.push_back("Si");
	lThick1.push_back(0.1*mm);lEle1.push_back("Si");
	m_caloStruct.push_back( SamplingSection(lThick1,lEle1) );

	std::vector<G4double> lThickL;
	std::vector<std::string> lEleL;
	std::vector<G4double> lThickR;
	std::vector<std::string> lEleR;
	for (unsigned i=0;i<14;i++){
	  if (dropLayer_[2*i]) {
	    lThickR.clear();lEleR.clear();
	    lThickR.push_back(3*mm);lEleR.push_back("Cu");
	    lThickR.push_back(1*mm);lEleR.push_back("Pb");
	    //reset to 0.5 Cu and no lead for following layers
	    if (i>0) {
	      lThickR[0] = 0.5;
	      lThickR[1] = 0;
	    }
	    lThickR.push_back(absThickW_[i]);lEleR.push_back("W");
	    lThickR.push_back(0.5*mm);lEleR.push_back("Cu");
	    lThickR.push_back(airThick);lEleR.push_back("Air");
	    lThickR.push_back(absThickPb_[i]);lEleR.push_back("Pb");
	    lThickR.push_back(3*mm);lEleR.push_back("Cu");
	    lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	    lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	    lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	    lThickR.push_back(pcbThick);lEleR.push_back("PCB");
	    lThickR.push_back(airThick);lEleR.push_back("Air");
	    m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	  }
	  else {
	    lThickL.push_back(3*mm);lEleL.push_back("Cu");
	    lThickL.push_back(1*mm);lEleL.push_back("Pb");
	    //reset to 0.5 Cu and no lead for following layers
	    if (i>0) {
	      if (lThickL.size()==2){
		lThickL[0] = 0.5;
		lThickL[1] = 0;
	      }
	      else {
		lThickL[3] = 0.5;
		lThickL[4] = 0;
	      }
	    }
	    lThickL.push_back(absThickW_[i]);lEleL.push_back("W");
	    lThickL.push_back(0.5*mm);lEleL.push_back("Cu");
	    lThickL.push_back(airThick);lEleL.push_back("Air");
	    lThickL.push_back(pcbThick);lEleL.push_back("PCB");
	    lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	    lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	    lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	    m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	    lThickL.clear();lEleL.clear();
	    if (dropLayer_[2*i+1]) {
	      lThickL.clear();lEleL.clear();
	      lThickL.push_back(3*mm);lEleL.push_back("Cu");
	      lThickL.push_back(absThickPb_[i]);lEleL.push_back("Pb");
	      lThickL.push_back(airThick);lEleL.push_back("Air");
	    }
	    else {
	      lThickR.clear();lEleR.clear();
	      lThickR.push_back(3*mm);lEleR.push_back("Cu");
	      lThickR.push_back(absThickPb_[i]);lEleR.push_back("Pb");
	      lThickR.push_back(3*mm);lEleR.push_back("Cu");
	      lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	      lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	      lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	      lThickR.push_back(pcbThick);lEleR.push_back("PCB");
	      lThickR.push_back(airThick);lEleR.push_back("Air");
	      m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	    }
	  }
	}

	/*	//lThickR.push_back(0.5*mm);lEleR.push_back("Cu");
	for(unsigned i=0; i<4; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}
	lThickL[2] = 2.8*mm;
	lThickR[1] = 2.1*mm;
	for(unsigned i=0; i<5; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}
	lThickL[2] = 4.2*mm;
	lThickR[1] = 4.4*mm;
	for(unsigned i=0; i<4; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}
	*/

	//last layer: add Cu+W+Cu...
	lThick1[0] = 0.*mm;
	lThick1[1] = 0.*mm;
	lThick1[2] = 0.*mm;
	lThick1[3] = 0.*mm;
	lThick1[4] = 0.5*mm;
	lThick1[5] = absThickW_[14];
	m_caloStruct.push_back( SamplingSection(lThick1,lEle1) );

	if(version_==v_HGCAL_v5 || version_==v_HGCAL_v5_gap4){
	  //add HCAL
	  buildHGCALFHE(5);
	  buildHGCALBHE(5);
	}

	break;
      }


    case v_HGCALEE_Si80: case v_HGCALEE_Si120: case v_HGCALEE_Si200: case v_HGCALEE_Si500: case v_HGCALEE_gap1: case  v_HGCALEE_CALICE: case v_HGCALEE_inverted: case v_HGCALEE_concept: case v_HGCALEE_W: case v_HGCALEE_gap4: case v_HGCALEE_prePCB: case v_HGCAL:
      {
	float siWidth(0.300), gap(2),pad(2);
	float pb1(1.60), pb2(3.30), pb3(5.60), cu(3.0);
	unsigned n1(10),n2(10),n3(10);

	if(version_==v_HGCALEE_Si80)     {siWidth=0.08;                               }
	if(version_==v_HGCALEE_Si120)    {siWidth=0.120;                              }
	if(version_==v_HGCALEE_Si500)    {siWidth=0.500;                              }
	if(version_==v_HGCALEE_gap1)     {gap=1;                                      }
	if(version_==v_HGCALEE_gap4)     {gap=4;                                      }
	if(version_==v_HGCALEE_CALICE)   {pb1=1.07;                                   }
	if(version_==v_HGCALEE_inverted) {pb2=1.63; pb1=3.32;                         }
	if(version_==v_HGCALEE_W)        {pb1=0.70;  pb2=2.10; pb3=3.50;}
	if(version_==v_HGCALEE_prePCB)   {addPrePCB_=true;}
	G4cout << "[DetectorConstruction] starting v_HGCAL* with Si width=" << siWidth << " and gap " << gap << G4endl;

	//add 1um silicon to track incoming particles
	//if(version_==v_HGCALEE_SiDummy) m_caloStruct.push_back( SamplingSection(0,0,0.001*mm,0,0) );
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(pb1);lEle.push_back(version_==v_HGCALEE_W?"W":"Pb");
	lThick.push_back(cu); lEle.push_back("Cu");
	//add PCB to shield from delta-rays?
	if (addPrePCB_) {
	  lThick.push_back(pad);lEle.push_back("PCB");
	}
	lThick.push_back(siWidth/3.);lEle.push_back("Si");
	lThick.push_back(siWidth/3.);lEle.push_back("Si");
	lThick.push_back(siWidth/3.);lEle.push_back("Si");
	lThick.push_back(pad);lEle.push_back("PCB");
	lThick.push_back(gap);lEle.push_back("Air");

	if(version_==v_HGCALEE_concept || version_==v_HGCAL) {
	  lThick[0] = 0;
	  //lThick[1] = 0;
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}

	lThick[0] = pb1;
	lThick[1] = cu;
	for(unsigned i=0; i<n1; i++)         m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = pb2;
	for(unsigned i=0; i<n2; i++)         m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = pb3;
	for(unsigned i=0; i<n3; i++)         m_caloStruct.push_back( SamplingSection(lThick,lEle) );

	if (version_==v_HGCAL){
	  buildHGCALFHE(4);
	  buildHGCALBHE(4);
	}
	break;
      }
    case v_HGCALHE: case v_HGCALHE_CMSSWv4: case v_HGCALHE_v5: case v_HGCALHE_v6: case v_HGCALHE_v7: case v_HGCALHE_v624: case v_HGCALHE_v618:
      {
	//add HCAL
	if (version_== v_HGCALHE_CMSSWv4) {
	  buildHGCALFHE(41);
	  buildHGCALBHE(4);
	}
	else if (version_== v_HGCALHE_v5) {
	  buildHGCALFHE(5);
	  buildHGCALBHE(5);
	}
	else if (version_== v_HGCALHE_v6) {
	  buildHGCALFHE(6);
	  buildHGCALBHE(6);
	}
	else if (version_== v_HGCALHE_v7) {
	  buildHGCALFHE(7);
	  buildHGCALBHE(7);
	}
	else if (version_== v_HGCALHE_v624) {
	  buildHGCALFHE(624);
	  buildHGCALBHE(6);
	}
	else if (version_== v_HGCALHE_v618) {
	  buildHGCALFHE(618);
	  buildHGCALBHE(6);
	}
	else {
	  buildHGCALFHE(4);
	  buildHGCALBHE(4);
	}
	break;
      }
    case v_HGCALHE_v8:
      {
	buildHGCALFHE(8);
	buildHGCALBHE(8);
	break;
      }
    case v_HGCALHEScint:
      {
	buildHGCALBHE(4);
	break;
      }
    case v_HGCALBE_v5:
      {
	buildHGCALBHE(5);
	break;
      }
    case v_HGCALBE_v6:
      {
	buildHGCALBHE(6);
	break;
      }
    case v_HGCALBE_v8:
      {
	buildHGCALBHE(8);
	break;
      }
    case v_HGCALHE_CALICE:
      {
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(16.7*mm);lEle.push_back("Steel");
	lThick.push_back(1.25*mm);lEle.push_back("Air");
	lThick.push_back(2*mm);lEle.push_back("Steel");
	lThick.push_back(1.5*mm);lEle.push_back("CFMix");
	lThick.push_back(1*mm);lEle.push_back("PCB");
	lThick.push_back(0.115*mm);lEle.push_back("Polystyrole");
	lThick.push_back(5*mm);lEle.push_back("Scintillator");
	lThick.push_back(0.115*mm);lEle.push_back("Polystyrole");
	lThick.push_back(2*mm);lEle.push_back("Steel");
	lThick.push_back(1.25*mm);lEle.push_back("Air");
	//for(unsigned i=0; i<12; i++) m_caloStruct.push_back( SamplingSection(52,3.0,0.3,2.0,2.0) );
	//total 5.3 lambda = 22mm brass * 38 layers
	for(unsigned i=0; i<3; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 17.4*mm;
	for(unsigned i=0; i<21; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 17.6*mm;
	for(unsigned i=0; i<14; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	//last absorber chunk ??
	//m_caloStruct.push_back( SamplingSection(20.5,0.,0,0,0) );

	//for(unsigned i=0; i<38; i++) m_caloStruct.push_back( SamplingSection(21,0.,5,0.,0.) );
	lThick[0] = 21*mm;
	for(unsigned i=0; i<9; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 104*mm;
	for(unsigned i=0; i<7; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	break;
      }

    case v_HGCALHF: case v_HGCALCu_HF:
      {
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	//lThick.push_back(1644*mm);lEle.push_back("Cu");
	//lThick.push_back(1067*mm);lEle.push_back("W");
	lThick.push_back(18*mm);lEle.push_back("W");
	for(unsigned i=0; i<60; i++){
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	buildHF();
	break;
      }

    }

  DefineMaterials();
  SetMagField(0);
  m_detectorMessenger = new DetectorMessenger(this);
  UpdateCalorSize();
}

void DetectorConstruction::buildHGCALFHE(const unsigned aVersion){
  std::vector<G4double> lThick;
  std::vector<std::string> lEle;

  // min Eta per layer, 0-27 for EE, 28-35 FH, 36-51 BH
  /*
  double minEta[52] = { 1.461, 1.464, 1.463, 1.466, 1.465, 1.468, 1.467, 1.469, 1.469, 1.471,
			1.471, 1.473, 1.472, 1.475, 1.474, 1.477, 1.476, 1.478, 1.478, 1.480,
			1.479, 1.482, 1.481, 1.483, 1.483, 1.485, 1.484, 1.487,
			1.487, 1.490, 1.494, 1.497, 1.500, 1.503, 1.506, 1.493,
			1.474, 1.455, 1.437, 1.420, 1.403, 1.376, 1.351, 1.327, 1.306, 1.316,
			1.332, 1.348, 1.364, 1.379, 1.395, 1.410 };
  std::cout << "DetectorConstruction::buildHGCALFHE " << std::endl;
  for (int i=0; i<minEta.size(); i++) std::cout << minEta[i] << std::endl;
  */
 
  if (aVersion==8){
 
    bool useRealisticBoundaries( version_==v_HGCAL_v8_envelope);
    if(useRealisticBoundaries) {
      std::cout << "buildHGCALFHE version 8 with realistic boundaries" << std::endl;
    }

    G4double airThick = 1.5*mm;
    G4double pcbThick = 1.6*mm;
    if (version_==v_HGCAL_v8 || version_==v_HGCAL_v8_envelope){
      //back of ecal
      lThick.push_back(pcbThick);lEle.push_back("PCB");
      lThick.push_back(airThick);lEle.push_back("Air");
      lThick.push_back(pcbThick);lEle.push_back("PCB");
      lThick.push_back(0.1*mm);lEle.push_back("Cu");
      lThick.push_back(1*mm);lEle.push_back("SSteel");
    }
    //CE-H layers 1
    lThick.push_back(40*mm);lEle.push_back("SSteel");
    lThick.push_back(6*mm);lEle.push_back("Cu");
    lThick.push_back(1*mm);lEle.push_back("SSteel");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    m_minEta.push_back(minEta[28]);m_maxEta.push_back(m_maxEta0);
    
    //CE-H layers 2-9
    lThick.clear();
    lEle.clear();
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(airThick);lEle.push_back("Air");
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(1*mm);lEle.push_back("Cu");
    lThick.push_back(1*mm);lEle.push_back("Air");
    lThick.push_back(35*mm);lEle.push_back("SSteel");
    lThick.push_back(6*mm);lEle.push_back("Cu");
    lThick.push_back(1*mm);lEle.push_back("SSteel");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");

    //layers 2-8          
    for(unsigned i=0; i<7; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      m_minEta.push_back(minEta[29+i]);m_maxEta.push_back(m_maxEta0); // minEta index: 29,30,...,35      
    }

    //layer 9
    firstMixedlayer_ = m_caloStruct.size();
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    if(useRealisticBoundaries){
      float z=3920.7; //here we latch to whatever was already there
      m_minEta.push_back(getEtaFromRZ(getRealisticRho(z,"mixed"),z));
      m_maxEta.push_back(getEtaFromRZ(getRealisticRho(z,"inner"),z));
      std::cout << z << " " << getRealisticRho(z,"mixed") << " " << getRealisticRho(z,"inner") << std::endl;
    }
    else{
      m_minEta.push_back(getEtaFromRZ(1450,3920.7));m_maxEta.push_back(m_maxEta0);
    }
    //CE-H layers 10-12
    lThick.clear();
    lEle.clear();
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(airThick);lEle.push_back("Air");
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    //added in layer 9 but after Si so belongs to layer10 here...
    lThick.push_back(1.9*mm);lEle.push_back("Air");
    lThick.push_back(1*mm);lEle.push_back("Cu");
    lThick.push_back(1*mm);lEle.push_back("Air");
    lThick.push_back(35*mm);lEle.push_back("SSteel");
    lThick.push_back(6*mm);lEle.push_back("Cu");
    lThick.push_back(1*mm);lEle.push_back("SSteel");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    for(unsigned i=0; i<3; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    }
    if(useRealisticBoundaries){
      float z[]={3969.7,4020.6,4071.5};
      for(size_t i=0; i<sizeof(z)/sizeof(float); i++){
        m_minEta.push_back(getEtaFromRZ(getRealisticRho(z[i],"mixed"),z[i]));
        m_maxEta.push_back(getEtaFromRZ(getRealisticRho(z[i],"inner"),z[i]));
        std::cout << z[i] << " " << getRealisticRho(z[i],"mixed") << " " << getRealisticRho(z[i],"inner") << std::endl;
      }
    }
    else {
      m_minEta.push_back(getEtaFromRZ(1325,3969.7));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(1325,4020.6));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(1225,4071.5));m_maxEta.push_back(m_maxEta0);
    }

    //CE-H layers 13-24
    lThick[6] = 68*mm;
    for(unsigned i=0; i<12; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    }
    if(useRealisticBoundaries){
      //here the last layer is forced to be smaller as in v73
      float z[]={4122.4,4206.3,4290.2,4374.1,4458,4541.9,4625.8,4709.7,4793.6,4877.5,4961.4,5262}; //5045.3};
      for(size_t i=0; i<sizeof(z)/sizeof(float);i++){
        std::cout << z[i] << " " << getRealisticRho(z[i],"mixed") << " " << getRealisticRho(z[i],"inner") << std::endl;
        m_minEta.push_back(getEtaFromRZ(getRealisticRho(z[i],"mixed"),z[i]));
        m_maxEta.push_back(getEtaFromRZ(getRealisticRho(z[i],"inner"),z[i]));
      }
    }
    else{
      m_minEta.push_back(getEtaFromRZ(1110,4122.4));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(1050,4206.3));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(950,4290.2));m_maxEta.push_back(m_maxEta0);
      
      m_minEta.push_back(getEtaFromRZ(900,4374.1));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(900,4458));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(900,4541.9));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(900,4625.8));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(900,4709.7));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(900,4793.6));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(900,4877.5));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(900,4961.4));m_maxEta.push_back(m_maxEta0);
      m_minEta.push_back(getEtaFromRZ(900,5045.3));m_maxEta.push_back(m_maxEta0);      
    }

    //end of last layer
    lThick.clear();
    lEle.clear();
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(airThick);lEle.push_back("Air");
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(1.9*mm);lEle.push_back("Air");
    lThick.push_back(1*mm);lEle.push_back("Cu");
    lThick.push_back(1*mm);lEle.push_back("Air");

    //back disk
    lThick.push_back(93.9*mm);lEle.push_back("SSteel");
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    if(useRealisticBoundaries) {
      float z(5262);
      m_minEta.push_back(getEtaFromRZ(getRealisticRho(z,"outer"),z));
      m_maxEta.push_back(getEtaFromRZ(getRealisticRho(z,"inner"),z));
    }
    else{
      m_minEta.push_back(m_minEta[51]);m_maxEta.push_back(m_maxEta0);
    }
  }

  //Jan2021 version
  else if (aVersion==9){

    G4double z0(2980);
    G4double totalThick(0.);
    for(size_t i=0; i<m_caloStruct.size(); i++)
      totalThick += m_caloStruct[i].Total_thick;  
    std::cout << "HGCALFHE starting at " << z0 << " + " << totalThick << " = " << z0+totalThick << std::endl;
    

    //
    // CE-H fine
    //
    
    //CE-H layer 1
    lThick.clear(); lEle.clear();
    
    //back cover for CEE   
    lThick.push_back(1*mm);  lEle.push_back("Cu");
    lThick.push_back(45*mm); lEle.push_back("SSteel");
    
    //CEH gap
    std::vector<G4double> gap_lThick; 
    std::vector<std::string> gap_lEle;    
    gap_lThick.push_back(0.7*mm);   gap_lEle.push_back("Air");
    gap_lThick.push_back(2.5*mm);   gap_lEle.push_back("Cu");
    gap_lThick.push_back(2.5*mm);   gap_lEle.push_back("Air");
    gap_lThick.push_back(1.6*mm);   gap_lEle.push_back("PCB");
    gap_lThick.push_back(2.5*mm);   gap_lEle.push_back("Air");
    gap_lThick.push_back(1.6*mm);   gap_lEle.push_back("PCB");
    gap_lThick.push_back(2.5*mm);   gap_lEle.push_back("Air");
    for(int i=0; i<3; i++){
      gap_lThick.push_back(0.1*mm); gap_lEle.push_back("Si");
    }
    gap_lThick.push_back(1.0*mm);   gap_lEle.push_back("PCB"); //baseplate
    gap_lThick.push_back(6.35*mm);  gap_lEle.push_back("Cu");  //cooling plate

    lThick.insert(lThick.end(), gap_lThick.begin(), gap_lThick.end());
    lEle.insert(lEle.end(),     gap_lEle.begin(),   gap_lEle.end());
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    totalThick += m_caloStruct[m_caloStruct.size()-1].Total_thick;
    float z=z0+totalThick;
    float rinner=getRealisticRho(z,"inner");
    float router=getRealisticRho(z,"outer");
    float rmixed=getRealisticRho(z,"mixed");
    m_minEta.push_back(getEtaFromRZ(router,z));
    m_maxEta.push_back(getEtaFromRZ(rinner,z));

    //CE-H layers 2-10
    lThick.clear(); lEle.clear();
    lThick.push_back(41.5);  lEle.push_back("SSteel");
    lThick.insert(lThick.end(), gap_lThick.begin(), gap_lThick.end());
    lEle.insert(lEle.end(), gap_lEle.begin(), gap_lEle.end());
    for(unsigned i=0; i<6; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      totalThick += m_caloStruct[m_caloStruct.size()-1].Total_thick;
      z=z0+totalThick;
      rinner=getRealisticRho(z,"inner");
      router=getRealisticRho(z,"outer");      
      m_minEta.push_back(getEtaFromRZ(router,z));
      m_maxEta.push_back(getEtaFromRZ(rinner,z));
    }
    firstMixedlayer_ = m_caloStruct.size();
    for(unsigned i=0; i<4; i++){
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      totalThick += m_caloStruct[m_caloStruct.size()-1].Total_thick;
      z=z0+totalThick;
      rinner=getRealisticRho(z,"inner");
      rmixed=getRealisticRho(z,"mixed");
      m_minEta.push_back(getEtaFromRZ(rmixed,z));
      m_maxEta.push_back(getEtaFromRZ(rinner,z));
    }

    //
    // CE-H coarse
    //
    //CE-H layers 12-22
    lThick.clear(); lEle.clear();
    lThick.push_back(60.7*mm);  lEle.push_back("SSteel");
    lThick.insert(lThick.end(), gap_lThick.begin(), gap_lThick.end());
    lEle.insert(lEle.end(), gap_lEle.begin(), gap_lEle.end());
    for(unsigned i=0; i<10; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      totalThick += m_caloStruct[m_caloStruct.size()-1].Total_thick;
      z=z0+totalThick;
      rinner=getRealisticRho(z,"inner");
      rmixed=getRealisticRho(z,"mixed");
      m_minEta.push_back(getEtaFromRZ(rmixed,z));
      m_maxEta.push_back(getEtaFromRZ(rinner,z));
    }
    
    //end of last layer
    lThick.clear();            lEle.clear();
    lThick.push_back(93.9*mm); lEle.push_back("SSteel");
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    m_minEta.push_back( m_minEta[m_minEta.size()-1] );
    m_maxEta.push_back( m_maxEta[m_maxEta.size()-1] );
  }
  else {
    
    G4double airThick = 2*mm;
    if(version_==v_HGCAL_v5_gap4) airThick = 4*mm;
   
    if(aVersion==6 || aVersion==624 || aVersion==618 || aVersion==7) {
      bool isBrass = aVersion!=7;
      airThick = 2*mm;
      G4double pcbthick = 2*mm;
      G4double brassthick = aVersion==618? 62*mm : 35*mm;
      //putting all absorber in front of each Si layer to have correct reweighting
      lThick.push_back(0.5*mm); lEle.push_back("Cu");
      lThick.push_back(15.*mm);lEle.push_back("SSteel");
      if (isBrass) {lThick.push_back(brassthick);lEle.push_back("Brass");}
      else {lThick.push_back(brassthick);lEle.push_back("SSteel");}
      lThick.push_back(0.5*mm); lEle.push_back("Cu");
      lThick.push_back(airThick);lEle.push_back("Air");
      lThick.push_back(pcbthick);lEle.push_back("PCB");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      
      lThick.clear();
      lEle.clear();
      lThick.push_back(1.);lEle.push_back("CFMix");
      lThick.push_back(6.); lEle.push_back("Cu");
      if (isBrass) {lThick.push_back(brassthick);lEle.push_back("Brass");}
      else {lThick.push_back(brassthick);lEle.push_back("SSteel");}
      lThick.push_back(0.5*mm); lEle.push_back("Cu");
      lThick.push_back(airThick);lEle.push_back("Air");
      lThick.push_back(pcbthick);lEle.push_back("PCB");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      for(unsigned i=0; i<5; i++) {
        m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      }
      unsigned nLay = aVersion==618? 3 : aVersion==624 ? 5 : 6;
      lThick[2] = aVersion==624 ? 55*mm : brassthick;
      for(unsigned i=0; i<nLay; i++) {
        m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      }
      
      
    }
    else {
      G4double pcbthick = (aVersion==4)? 2*mm : 1.2*mm;
      //add last ECAL layer structure
      lThick.push_back(3*mm);lEle.push_back("Cu");
      lThick.push_back(1*mm);lEle.push_back("Pb");
      lThick.push_back(15.*mm);lEle.push_back("SSteel");
      if (aVersion==41) {lThick.push_back(52.*mm);lEle.push_back("Pb");}
      else {lThick.push_back(40.*mm);lEle.push_back("Brass");}
      lThick.push_back(0.5*mm); lEle.push_back("Cu");
      lThick.push_back(airThick);lEle.push_back("Air");
      lThick.push_back(pcbthick);lEle.push_back("PCB");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      
      //next 11 layers
      lThick.clear(); lEle.clear();
      lThick.push_back(3*mm); lEle.push_back("Cu");
      lThick.push_back(1*mm); lEle.push_back("Pb");
      if (aVersion==41) {lThick.push_back(52.*mm);lEle.push_back("Pb");}
      else {lThick.push_back(40.*mm);lEle.push_back("Brass");}
      lThick.push_back(0.5*mm); lEle.push_back("Cu");
      lThick.push_back(airThick);lEle.push_back("Air");
      lThick.push_back(pcbthick);lEle.push_back("PCB");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      lThick.push_back(0.1*mm);lEle.push_back("Si");
      
      for(unsigned i=0; i<11; i++) {
        m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      }
    }
  }
}
//
void DetectorConstruction::buildHGCALBHE(const unsigned aVersion){
  std::vector<G4double> lThick;
  std::vector<std::string> lEle;

  // min Eta per layer, 0-27 for EE, 28-35 FH, 36-51 BH
  /*
  double minEta[52] = { 1.461, 1.464, 1.463, 1.466, 1.465, 1.468, 1.467, 1.469, 1.469, 1.471,
			1.471, 1.473, 1.472, 1.475, 1.474, 1.477, 1.476, 1.478, 1.478, 1.480,
			1.479, 1.482, 1.481, 1.483, 1.483, 1.485, 1.484, 1.487,
			1.487, 1.490, 1.494, 1.497, 1.500, 1.503, 1.506, 1.493,
			1.474, 1.455, 1.437, 1.420, 1.403, 1.376, 1.351, 1.327, 1.306, 1.316,
			1.332, 1.348, 1.364, 1.379, 1.395, 1.410 };
  std::cout << "DetectorConstruction::buildHGCALBHE " << std::endl;
  for (int i=0; i<minEta.size(); i++) std::cout << minEta[i] << std::endl;
  */

  if (aVersion==8){
    firstScintlayer_ = m_caloStruct.size();

    bool useRealisticBoundaries(version_==v_HGCAL_v8_envelope);

    G4double airThick = 1.5*mm;
    G4double pcbThick = 1.6*mm;
    //back of layer 8+layer9 scint.
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(airThick);lEle.push_back("Air");
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(1*mm);lEle.push_back("Cu");
    lThick.push_back(1*mm);lEle.push_back("Air");
    lThick.push_back(35*mm);lEle.push_back("SSteel");
    lThick.push_back(6*mm);lEle.push_back("Cu");
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(1.2*mm);lEle.push_back("Air");
    lThick.push_back(3.0*mm);lEle.push_back("Scintillator");
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    if(useRealisticBoundaries){
      float z(3920.7);
      m_minEta.push_back(getEtaFromRZ(getRealisticRho(z,"outer"),z));
      m_maxEta.push_back(getEtaFromRZ(getRealisticRho(z,"mixed"),z));
      std::cout << z << " " << getRealisticRho(z,"outer") << " " << getRealisticRho(z,"mixed") << std::endl;
    }else{
      m_minEta.push_back(minEta[36]);m_maxEta.push_back(getEtaFromRZ(1450,3920.7));
    }

    //CE-H-scint layers 10-12
    lThick.clear();
    lEle.clear();
    lThick.push_back(0.5*mm);lEle.push_back("Air");
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(1*mm);lEle.push_back("Cu");
    lThick.push_back(1*mm);lEle.push_back("Air");
    lThick.push_back(35*mm);lEle.push_back("SSteel");
    lThick.push_back(6*mm);lEle.push_back("Cu");
    lThick.push_back(pcbThick);lEle.push_back("PCB");
    lThick.push_back(1.2*mm);lEle.push_back("Air");
    lThick.push_back(3.0*mm);lEle.push_back("Scintillator");
    for(unsigned i=0; i<3; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    }
    if(useRealisticBoundaries){
      float z[]={3969.7,4020.6,4071.5};
      for(size_t i=0; i<sizeof(z)/sizeof(float); i++) {
        m_minEta.push_back(getEtaFromRZ(getRealisticRho(z[i],"outer"),z[i]));
        m_maxEta.push_back(getEtaFromRZ(getRealisticRho(z[i],"mixed"),z[i]));
        std::cout << z[i] << " " << getRealisticRho(z[i],"outer") << " " << getRealisticRho(z[i],"mixed") << std::endl;
      }
    }
    else{
      m_minEta.push_back(minEta[37]);m_maxEta.push_back(getEtaFromRZ(1325,3969.7));
      m_minEta.push_back(minEta[38]);m_maxEta.push_back(getEtaFromRZ(1325,4020.6));
      m_minEta.push_back(minEta[39]);m_maxEta.push_back(getEtaFromRZ(1225,4071.5));
    }
    firstCoarseScintlayer_ = m_caloStruct.size();

    //CE-H-scint layers 13-24
    lThick[4] = 68*mm;
    for(unsigned i=0; i<12; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    }
    if(useRealisticBoundaries){
      float z[]={4122.4,4206.3,4290.2,4374.1,4458,4541.9,4625.8,4709.7,4793.6,4877.5,4961.4,5045.3};
      for(size_t i=0; i<sizeof(z)/sizeof(float); i++) {
        m_minEta.push_back(getEtaFromRZ(getRealisticRho(z[i],"outer"),z[i]));
        m_maxEta.push_back(getEtaFromRZ(getRealisticRho(z[i],"mixed"),z[i]));
        std::cout << z[i] << " " << getRealisticRho(z[i],"outer") << " " << getRealisticRho(z[i],"mixed") << std::endl;
      }
    }
    else{
      m_minEta.push_back(minEta[40]);m_maxEta.push_back(getEtaFromRZ(1110,4122.4));
      m_minEta.push_back(minEta[41]);m_maxEta.push_back(getEtaFromRZ(1050,4206.3));
      m_minEta.push_back(minEta[42]);m_maxEta.push_back(getEtaFromRZ(950,4290.2));
   
      m_minEta.push_back(minEta[43]);m_maxEta.push_back(getEtaFromRZ(900,4374.1));
      m_minEta.push_back(minEta[44]);m_maxEta.push_back(getEtaFromRZ(900,4458));
      m_minEta.push_back(minEta[45]);m_maxEta.push_back(getEtaFromRZ(900,4541.9));
      m_minEta.push_back(minEta[46]);m_maxEta.push_back(getEtaFromRZ(900,4625.8));
      m_minEta.push_back(minEta[47]);m_maxEta.push_back(getEtaFromRZ(900,4709.7));
      m_minEta.push_back(minEta[48]);m_maxEta.push_back(getEtaFromRZ(900,4793.6));
      m_minEta.push_back(minEta[49]);m_maxEta.push_back(getEtaFromRZ(900,4877.5));
      m_minEta.push_back(minEta[50]);m_maxEta.push_back(getEtaFromRZ(900,4961.4));
      m_minEta.push_back(minEta[51]);m_maxEta.push_back(getEtaFromRZ(900,5045.3));
    }


    //end of last layer - only if building BH only
    if (version_ == v_HGCALBE_v8){
      lThick.clear();
      lEle.clear();
      lThick.push_back(0.5*mm);lEle.push_back("Air");
      lThick.push_back(pcbThick);lEle.push_back("PCB");
      lThick.push_back(1*mm);lEle.push_back("Cu");
      lThick.push_back(1*mm);lEle.push_back("Air");
      //back disk
      lThick.push_back(93.9*mm);lEle.push_back("SSteel");
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      m_minEta.push_back(minEta[51]);m_maxEta.push_back(m_maxEta0);
    }
  }
  //Jan2021 version
  else if (aVersion==9){

    firstScintlayer_ = m_caloStruct.size();
    G4double z0(2980);
    G4double totalThick(0.);
    for(size_t i=0; i<firstMixedlayer_; i++)
      totalThick += m_caloStruct[i].Total_thick;
    std::cout << "HGCALBHE starting at " << z0 << " + " << totalThick << " = " << z0+totalThick << std::endl;

    //
    // CE-H fine
    //
    
    //CE-H layer 1
    lThick.clear(); lEle.clear();
    
    //CEH gap for scintillator
    std::vector<G4double> gap_lThick; 
    std::vector<std::string> gap_lEle;    
    gap_lThick.push_back(0.7*mm);   gap_lEle.push_back("Air");
    gap_lThick.push_back(2.5*mm);   gap_lEle.push_back("Cu");
    gap_lThick.push_back(5.8*mm);   gap_lEle.push_back("Air");
    gap_lThick.push_back(1.6*mm);   gap_lEle.push_back("PCB");
    gap_lThick.push_back(3*mm);     gap_lEle.push_back("Scintillator");
    gap_lThick.push_back(1.6*mm);   gap_lEle.push_back("PCB");
    gap_lThick.push_back(6.35*mm);  gap_lEle.push_back("Cu");
    
    //fine section
    lThick.push_back(41.5);  lEle.push_back("SSteel");
    lThick.insert(lThick.end(), gap_lThick.begin(), gap_lThick.end());
    lEle.insert(lEle.end(), gap_lEle.begin(), gap_lEle.end());
    for(unsigned i=0; i<4; i++){
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      totalThick += m_caloStruct[m_caloStruct.size()-1].Total_thick;
      float z=z0+totalThick;
      float router=getRealisticRho(z,"outer");
      float rmixed=getRealisticRho(z,"mixed");
      m_minEta.push_back(getEtaFromRZ(router,z));
      m_maxEta.push_back(getEtaFromRZ(rmixed,z));
    }

    firstCoarseScintlayer_ = m_caloStruct.size();
    //
    //coarse section
    //
    //CE-H layers 12-22
    lThick.clear(); lEle.clear();
    lThick.push_back(60.7*mm);  lEle.push_back("SSteel");
    lThick.insert(lThick.end(), gap_lThick.begin(), gap_lThick.end());
    lEle.insert(lEle.end(), gap_lEle.begin(), gap_lEle.end());
    for(unsigned i=0; i<10; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
      totalThick += m_caloStruct[m_caloStruct.size()-1].Total_thick;
      float z=z0+totalThick;
      float router=getRealisticRho(z,"outer");
      float rmixed=getRealisticRho(z,"mixed");
      m_minEta.push_back(getEtaFromRZ(router,z));
      m_maxEta.push_back(getEtaFromRZ(rmixed,z));
    }
    
    //end of last layer
    lThick.clear();            lEle.clear();
    lThick.push_back(93.9*mm); lEle.push_back("SSteel");
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    m_minEta.push_back( m_minEta[m_minEta.size()-1] );
    m_maxEta.push_back( m_maxEta[m_maxEta.size()-1] );
  }
  else {
    //first layer
    bool isBrass = aVersion!=7;

    if (aVersion==6 || aVersion==7){
      lThick.push_back(1.*mm);lEle.push_back("CFMix");
      lThick.push_back(6.*mm); lEle.push_back("Cu");
    } else {
      lThick.push_back(3*mm); lEle.push_back("Cu");
      lThick.push_back(1*mm); lEle.push_back("Pb");
    }
    lThick.push_back(2.*mm);lEle.push_back("Al");
    lThick.push_back(16.*mm);lEle.push_back("Foam");
    lThick.push_back(2.*mm);lEle.push_back("Al");
    lThick.push_back(65.*mm);lEle.push_back("Air");
    if (isBrass) {lThick.push_back(78.*mm);lEle.push_back("Brass");}
    else {lThick.push_back(78.*mm);lEle.push_back("SSteel");}
    if (aVersion==6 || aVersion==7) {
      lThick.push_back(2.6*mm);lEle.push_back("Air");
      lThick.push_back(3.8*mm);lEle.push_back("Scintillator");
      lThick.push_back(2.6*mm);lEle.push_back("Air");
    }
    else {
      lThick.push_back(9.*mm);lEle.push_back("Scintillator");
    }
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );

    //other layers
    lThick.clear();lEle.clear();
    if (isBrass) {lThick.push_back(78.*mm);lEle.push_back("Brass");}
    else {lThick.push_back(78.*mm);lEle.push_back("SSteel");}
    if (aVersion==6 || aVersion==7) {
      lThick.push_back(2.6*mm);lEle.push_back("Air");
      lThick.push_back(3.8*mm);lEle.push_back("Scintillator");
      lThick.push_back(2.6*mm);lEle.push_back("Air");
    }
    else {
      lThick.push_back(9.*mm);lEle.push_back("Scintillator");
    }

    unsigned maxi = (aVersion==4)?9:11;
    for(unsigned i=0; i<maxi; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    }
  }
}

void DetectorConstruction::buildHF(){
  firstHFlayer_ = m_caloStruct.size();
  G4double airThick = 4*mm;
  std::vector<G4double> lThick;
  std::vector<std::string> lEle;
  G4double pcbthick = 2*mm;
  lThick.clear();
  lEle.clear();
  //"ECAL" part
  //first layer without bit after Si
  lThick.push_back(0.); lEle.push_back("Cu");
  lThick.push_back(1.);lEle.push_back("CFMix");
  lThick.push_back(25*mm);lEle.push_back("Pb");
  lThick.push_back(0.5*mm); lEle.push_back("Cu");
  lThick.push_back(airThick);lEle.push_back("Air");
  lThick.push_back(pcbthick);lEle.push_back("PCB");
  lThick.push_back(0.2*mm);lEle.push_back("Si");
  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
  //put back bit after Si in next layer
  lThick[0] = 5*mm;
  for(unsigned i=0; i<5; i++) {
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );
  }
  //"HCAL" part
  lThick[2] = 65*mm;
  lEle[2] = "SSteel";
  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
  lThick.push_back(5.*mm); lEle.push_back("Cu");
  m_caloStruct.push_back( SamplingSection(lThick,lEle) );

  //add HF as of now
  //+11.1m from IP
  lThick.clear();
  lEle.clear();
  lThick.push_back(1*mm);lEle.push_back("Air");
  lThick.push_back(1650*mm); lEle.push_back("Cu");
  m_caloStruct.push_back( SamplingSection(lThick,lEle) );

}//buildHF


//
DetectorConstruction::~DetectorConstruction() { delete m_detectorMessenger;}


double DetectorConstruction::getEtaFromRZ(const double & r, const double & z){
  double theta = atan(r/z);
  return -log(tan(theta/2.));
}

//
void DetectorConstruction::DefineMaterials()
{
  std::cout << "[DetectorConstruction::DefineMaterials]" << std::endl;

  G4NistManager* nistManager = G4NistManager::Instance();
  m_materials["Abs"] = (version_== v_CALICE || version_==v_HGCALEE_W) ?
    nistManager->FindOrBuildMaterial("G4_W",false) :
    nistManager->FindOrBuildMaterial("G4_Pb",false);
  m_materials["Al"] = nistManager->FindOrBuildMaterial("G4_Al",false);
  m_dEdx["Al"] = 0.4358;
  m_materials["W"] = nistManager->FindOrBuildMaterial("G4_W",false);
  m_dEdx["W"] = 2.210;
  //G4cout << " Density of W: " <<  m_materials["W"]->GetDensity()/(g/mm3) << G4endl;
  //G4cout << m_materials["W"] << G4endl;
  m_materials["Pb"] = nistManager->FindOrBuildMaterial("G4_Pb",false);
  m_dEdx["Pb"] = 1.274;
  m_materials["Cu"] = nistManager->FindOrBuildMaterial("G4_Cu",false);
  m_dEdx["Cu"] = 1.257;
  m_materials["CuExtra"] = nistManager->FindOrBuildMaterial("G4_Cu",false);
  m_dEdx["CuExtra"] = 1.257;
  m_materials["Si"] = nistManager->FindOrBuildMaterial("G4_Si",false);
  m_dEdx["Si"] = 0.3876;
  m_materials["Zn"] = nistManager->FindOrBuildMaterial("G4_Zn",false);
  m_dEdx["Zn"] = 1.007;
  m_materials["Air"]=nistManager->FindOrBuildMaterial("G4_AIR",false);
  m_dEdx["Air"] = 0;
  m_materials["Fe"] = nistManager->FindOrBuildMaterial("G4_Fe",false);
  m_dEdx["Fe"] = 1.143;
  m_materials["Mn"] = nistManager->FindOrBuildMaterial("G4_Mn",false);
  m_dEdx["Mn"] = 1.062 ;
  m_materials["C"] = nistManager->FindOrBuildMaterial("G4_C",false);
  m_dEdx["C"] = 0.3952;
  m_materials["H"] = nistManager->FindOrBuildMaterial("G4_H",false);
  m_dEdx["H"] =  0;
  m_materials["Cl"] = nistManager->FindOrBuildMaterial("G4_Cl",false);
  m_dEdx["Cl"] = 0;
  m_materials["Cr"] = nistManager->FindOrBuildMaterial("G4_Cr",false);
  m_dEdx["Cr"] = 1.046;
  m_materials["Ni"] = nistManager->FindOrBuildMaterial("G4_Ni",false);
  m_dEdx["Ni"] = 1.307;
  m_materials["O"] = nistManager->FindOrBuildMaterial("G4_O",false);
  m_materials["Br"] = nistManager->FindOrBuildMaterial("G4_Br",false);

  /*m_materials["PCB"] = new G4Material("G10",1.700*g/cm3,4);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(14), 1);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(8) , 2);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(6) , 3);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(1) , 3);
  m_dEdx["PCB"] = 0;*/
  m_materials["PCB"] = new G4Material("FR4",1.700*g/cm3,5); //from google: 1.850*g/cm3,5);
  m_materials["PCB"]->AddMaterial(m_materials["Si"] , 0.18077359);
  m_materials["PCB"]->AddMaterial(m_materials["O"]  , 0.4056325);
  m_materials["PCB"]->AddMaterial(m_materials["C"]  , 0.27804208);
  m_materials["PCB"]->AddMaterial(m_materials["H"]  , 0.068442752);
  m_materials["PCB"]->AddMaterial(m_materials["Br"] , 0.067109079);
  m_dEdx["PCB"] =  0.399; //matching CMSSW

  m_materials["Brass"]= new G4Material("Brass",8.53*g/cm3,2);
  m_materials["Brass"]->AddMaterial(m_materials["Cu"]  , 70*perCent);
  m_materials["Brass"]->AddMaterial(m_materials["Zn"]  , 30*perCent);

  //G4cout << m_materials["Brass"] << G4endl
  //	 << "-- check density: formula " << 1./(0.7/(m_materials["Cu"]->GetDensity()/(g/mm3))+0.3/(m_materials["Zn"]->GetDensity()/(g/mm3)))
  //	 << " input density " << m_materials["Brass"]->GetDensity()/(g/mm3)
  //	 << G4endl;

  m_dEdx["Brass"] = m_materials["Brass"]->GetDensity()/(g/mm3)*(0.7/(m_materials["Cu"]->GetDensity()/(g/mm3))*m_dEdx["Cu"]+0.3/(m_materials["Zn"]->GetDensity()/(g/mm3))*m_dEdx["Zn"]);

  //G4cout << " Check dEdx: "
  //	 << " -- simple weighted sum (wrong): "<< 0.7*m_dEdx["Cu"]+0.3*m_dEdx["Zn"] << G4endl
  //	 << " -- with Brass density: " << m_dEdx["Brass"] << G4endl
  //	 << " -- Chris's formula: " << (0.7/(m_materials["Cu"]->GetDensity()/(g/mm3))*m_dEdx["Cu"]+0.3/(m_materials["Zn"]->GetDensity()/(g/mm3))*m_dEdx["Zn"])/(0.7/(m_materials["Cu"]->GetDensity()/(g/mm3))+0.3/(m_materials["Zn"]->GetDensity()/(g/mm3))) << G4endl;

  m_materials["Steel"]= new G4Material("Steel",7.87*g/cm3,3);
  m_materials["Steel"]->AddMaterial(m_materials["Fe"]  , 0.9843);
  m_materials["Steel"]->AddMaterial(m_materials["Mn"], 0.014);
  m_materials["Steel"]->AddMaterial(m_materials["C"], 0.0017);

  m_dEdx["Steel"] = m_materials["Steel"]->GetDensity()/(g/mm3)*(0.9843/(m_materials["Fe"]->GetDensity()/(g/mm3))*m_dEdx["Fe"]+0.014/(m_materials["Mn"]->GetDensity()/(g/mm3))*m_dEdx["Mn"]+0.0017/(m_materials["C"]->GetDensity()/(g/mm3))*m_dEdx["C"]);

  m_materials["SSteel"]= new G4Material("SSteel",8.02*g/cm3,4);
  m_materials["SSteel"]->AddMaterial(m_materials["Fe"]  , 0.70);
  m_materials["SSteel"]->AddMaterial(m_materials["Mn"], 0.01);
  m_materials["SSteel"]->AddMaterial(m_materials["Cr"], 0.19);
  m_materials["SSteel"]->AddMaterial(m_materials["Ni"], 0.10);

  //m_dEdx["SSteel"] = 0.7*m_dEdx["Fe"]+0.01*m_dEdx["Mn"]+0.19*m_dEdx["Cr"]+0.1*m_dEdx["Ni"];
  m_dEdx["SSteel"] = m_materials["SSteel"]->GetDensity()/(g/mm3)*(0.7/(m_materials["Fe"]->GetDensity()/(g/mm3))*m_dEdx["Fe"]+0.01/(m_materials["Mn"]->GetDensity()/(g/mm3))*m_dEdx["Mn"]+0.19/(m_materials["Cr"]->GetDensity()/(g/mm3))*m_dEdx["Cr"]+0.1/(m_materials["Ni"]->GetDensity()/(g/mm3))*m_dEdx["Ni"]);

  //G4cout <<  m_materials["SSteel"] << G4endl
  //<< " Check dEdx: " << G4endl
  //<< " -- simple weighted sum (wrong): "<< 0.7*m_dEdx["Fe"]+0.01*m_dEdx["Mn"]+0.19*m_dEdx["Cr"]+0.1*m_dEdx["Ni"] << G4endl
  //<< " -- with SSteel density: " << m_dEdx["SSteel"] << G4endl
  //<< " -- formula: " << (0.7/(m_materials["Fe"]->GetDensity()/(g/mm3))*m_dEdx["Fe"]+0.01/(m_materials["Mn"]->GetDensity()/(g/mm3))*m_dEdx["Mn"]+0.19/(m_materials["Cr"]->GetDensity()/(g/mm3))*m_dEdx["Cr"]+0.1/(m_materials["Ni"]->GetDensity()/(g/mm3))*m_dEdx["Ni"])/(0.7/(m_materials["Fe"]->GetDensity()/(g/mm3))+0.01/(m_materials["Mn"]->GetDensity()/(g/mm3))+0.19/(m_materials["Cr"]->GetDensity()/(g/mm3))+0.1/(m_materials["Ni"]->GetDensity()/(g/mm3))) << G4endl;


  m_materials["AbsHCAL"] = (version_== v_HGCALHE_CALICE) ?
    m_materials["Steel"]:
    m_materials["Brass"];
  m_dEdx["AbsHCAL"] = (version_== v_HGCALHE_CALICE) ?
    m_dEdx["Steel"]:
    m_dEdx["Brass"];
  //m_materials["Scintillator"]= nistManager->FindOrBuildMaterial("G4_POLYSTYRENE",false);
  m_materials["Scintillator"]= new G4Material("Scintillator",1.032*g/cm3,2);
  m_materials["Scintillator"]->AddMaterial(m_materials["C"]  , 91.512109*perCent);
  m_materials["Scintillator"]->AddMaterial(m_materials["H"]  , 8.4878906*perCent);
  m_dEdx["Scintillator"] = m_dEdx["C"];

  m_materials["Polystyrole"]= new G4Material("Polystyrole",1.065*g/cm3,2);
  m_materials["Polystyrole"]->AddMaterial(m_materials["H"]  , 50*perCent);
  m_materials["Polystyrole"]->AddMaterial(m_materials["C"]  , 50*perCent);
  m_dEdx["Polystyrole"] = 0.5*m_dEdx["C"];

  m_materials["PVC"]= new G4Material("PVC",1.350*g/cm3,3);
  m_materials["PVC"]->AddMaterial(m_materials["H"]  , 50*perCent);
  m_materials["PVC"]->AddMaterial(m_materials["C"]  , 33.33*perCent);
  m_materials["PVC"]->AddMaterial(m_materials["Cl"]  , 16.67*perCent);
  m_dEdx["PVC"] = 0.33*m_dEdx["C"];

  /*m_materials["CFMix"]= new G4Material("CFMix",0.120*g/cm3,3);
  m_materials["CFMix"]->AddMaterial(m_materials["Air"]  , 0.009);
  m_materials["CFMix"]->AddMaterial(m_materials["PVC"]  , 0.872);
  m_materials["CFMix"]->AddMaterial(m_materials["Polystyrole"]  , 0.119);
  m_dEdx["CFMix"] = 0;*/

  m_materials["CFMix"]= new G4Material("CFMix",0.845*g/cm3,3);
  m_materials["CFMix"]->AddMaterial(m_materials["C"]  , 0.84491305);
  m_materials["CFMix"]->AddMaterial(m_materials["H"]  , 0.042542086);
  m_materials["CFMix"]->AddMaterial(m_materials["O"]  , 0.11254487);
  m_dEdx["CFMix"] = 0;

  m_materials["Foam"]= new G4Material("Foam",0.0999*g/cm3,2);
  m_materials["Foam"]->AddMaterial(m_materials["C"]  , 0.856);
  m_materials["Foam"]->AddMaterial(m_materials["H"]  , 0.144);
  m_dEdx["Foam"] = 1.749*0.856*0.0999/10.;

  m_materials["WCu"]= new G4Material("WCu",14.979*g/cm3,2);
  m_materials["WCu"]->AddMaterial(m_materials["W"]  , 75*perCent);
  m_materials["WCu"]->AddMaterial(m_materials["Cu"]  , 25*perCent);

  m_dEdx["WCu"] = m_materials["WCu"]->GetDensity()/(g/mm3)*(0.75/(m_materials["W"]->GetDensity()/(g/mm3))*m_dEdx["W"]+0.25/(m_materials["Cu"]->GetDensity()/(g/mm3))*m_dEdx["Cu"]);

  //G4cout << m_materials["WCu"] << G4endl
  //<< "Check dEdx: " << G4endl
  //	 << " -- simple weighted sum (wrong): " << 0.75*m_dEdx["W"]+0.25*m_dEdx["Cu"] << G4endl
  //	 << " -- with WCu density: " << m_dEdx["WCu"] << G4endl
  //	 << " -- Chris's formula: " << (0.75/(m_materials["W"]->GetDensity()/(g/mm3))*m_dEdx["W"]+0.25/(m_materials["Cu"]->GetDensity()/(g/mm3))*m_dEdx["Cu"])/(0.75/(m_materials["W"]->GetDensity()/(g/mm3))+0.25/(m_materials["Cu"]->GetDensity()/(g/mm3)))
  //	 << G4endl;


  m_materials["NeutMod"]= new G4Material("NeutMod",0.950*g/cm3,2);
  m_materials["NeutMod"]->AddMaterial(m_materials["C"]  , 0.85628);
  m_materials["NeutMod"]->AddMaterial(m_materials["H"]  , 0.14372);
  m_dEdx["NeutMod"] = 1.749*0.86*0.950/10.;


  G4cout << m_materials["PCB"] << G4endl;
  G4cout << m_materials["Scintillator"] << G4endl;
  G4cout << m_materials["Si"] << G4endl;
  G4cout << m_materials["Pb"] << G4endl;
  G4cout << m_materials["W"] << G4endl;
  G4cout << m_materials["WCu"] << G4endl;
  G4cout << m_materials["Cu"] << G4endl;
  G4cout << m_materials["SSteel"] << G4endl;
}

//
void DetectorConstruction::UpdateCalorSize(){

  m_CalorSizeZ=0;
  double HFsize = 0;
  for(size_t i=0; i<m_caloStruct.size(); i++){
    if (i>=firstHFlayer_) HFsize += m_caloStruct[i].Total_thick;
    m_CalorSizeZ=m_CalorSizeZ+m_caloStruct[i].Total_thick;
  }
  std::cout << m_CalorSizeZ << std::endl;

  m_z0pos = -m_CalorSizeZ/2.;

  m_nSectors = 1;
  m_interSectorWidth = 0;
  if (model_ == DetectorConstruction::m_SIMPLE_20){
    m_CalorSizeXY=200;
    m_sectorWidth = m_CalorSizeXY;
  }
  else if (model_ == DetectorConstruction::m_SIMPLE_50){
    m_CalorSizeXY=500;
    m_sectorWidth = m_CalorSizeXY;
  }
  else if (model_ == DetectorConstruction::m_SIMPLE_100){
    m_CalorSizeXY=1000;
    m_sectorWidth = m_CalorSizeXY;
  }
  else if (model_ == DetectorConstruction::m_BOXWITHCRACK_100 ){
    m_nSectors = 3;
    m_sectorWidth = 460;
    m_interSectorWidth = 10;
    m_CalorSizeXY=m_nSectors*m_sectorWidth;
  }
  else if (model_ == DetectorConstruction::m_2016TB ){
    m_nSectors    = 1;
    m_sectorWidth = (m_coarseGranularity>0 ? CELL_SIZE_X : m_coarseGranularity<0 ? ULTRAFINE_CELL_SIZE_X : FINE_CELL_SIZE_X) * 2 *11;
    m_interSectorWidth = 0;
    m_CalorSizeXY = m_sectorWidth*m_nSectors;
    m_minRadius   = m_CalorSizeXY/(2*sqrt(3)); // center-to-side radius of hexagon
    m_maxRadius   = m_CalorSizeXY/2.;          // center-to-corner radius of hexagon
  }
  else if (model_ == DetectorConstruction::m_FULLSECTION){
    m_CalorSizeXY=2800*2;//use full length for making hexagon map
    m_minRadius = 150;
    m_maxRadius = m_CalorSizeXY;
    if (version_ != v_HGCAL_v8 && version_ != v_HGCALEE_v8 && version_!=v_HGCALEE_v8_neutmod && version_!=v_HGCALEE_v8_air3 && version_!=v_HGCALEE_v8_air4 && version_!=v_HGCALEE_v8_Cu && version_!=v_HGCALEE_v8_Cu_12 && version_!= v_HGCAL_v8_envelope && version_!=v_HGCAL_v9 && version_!=v_HGCALEE_v9 && version_!=v_HGCAL_v10 && version_!=v_HGCALEE_v10) {
      m_minEta.resize(m_caloStruct.size(),m_minEta0);
      m_maxEta.resize(m_caloStruct.size(),m_maxEta0);
    }
    m_minEtaHF = 3.0;
    m_maxEtaHF = 4.3;
    m_sectorWidth = 2.*pi;//2.*pi/m_nSectors;
    m_interSectorWidth = 0;//0.2*pi/180.;
    m_z0pos = 2990;//3170;
    if (version_ == v_HGCALEE_v5 || version_ == v_HGCAL_v5 || version_ == v_HGCALEE_v5_gap4 || version_ == v_HGCAL_v5_gap4) m_z0pos = 2990;//3170;
    else if (version_ == v_HGCALEE_v6 || version_ == v_HGCAL_v6 || version_ == v_HGCALEE_v7 || version_ == v_HGCAL_v7 || version_ == v_HGCAL_v7_HF ||version_ == v_HGCALEE_v624 || version_ == v_HGCALEE_v618) m_z0pos = 3070;
    else if (version_ == v_HGCALEE_v8 || version_==v_HGCALEE_v8_neutmod || version_ == v_HGCAL_v8 || version_==v_HGCALEE_v8_air3 || version_==v_HGCALEE_v8_air4 || version_==v_HGCALEE_v8_Cu || version_==v_HGCALEE_v8_Cu_12 || version_ == v_HGCAL_v8_envelope ||
             version_ == v_HGCALEE_v9 || version_ == v_HGCAL_v9 || version_==v_HGCAL_v10 || version_==v_HGCALEE_v10) {
      m_z0pos = 2980;
    }
    else if(version_ == v_HGCALBE_v8) {
      m_z0pos=3920.7;
    }

    if (doHF_){
      m_z0HF=11100;
      m_CalorSizeZ=m_z0HF-m_z0pos+HFsize;
    }
    else m_CalorSizeZ+=m_z0pos;
  }
  else {
    m_CalorSizeXY=200;
    m_sectorWidth = m_CalorSizeXY;
  }

  for(size_t i=0; i<m_caloStruct.size(); i++) m_caloStruct[i].setNumberOfSectors(m_nSectors);

  m_WorldSizeZ=m_CalorSizeZ*1.1;
  if (m_nSectors>1) m_WorldSizeXY=(m_CalorSizeXY+2*m_sectorWidth)*1.1;
  else m_WorldSizeXY=m_CalorSizeXY*1.1;

  if (model_ == DetectorConstruction::m_FULLSECTION || model_ == DetectorConstruction::m_2016TB) {
    G4cout << "[DetectorConstruction][UpdateCalorSize] Z x minR * maxR = "
	   << m_CalorSizeZ << " x "
	   << m_minRadius << " x "
	   << m_maxRadius
	   << " mm, eta range "
	   << m_minEta0 << " - "
	   << m_maxEta0 ;
    if (doHF_) G4cout << " eta range HF "
		      << m_minEtaHF << " - "
		      << m_maxEtaHF;
    G4cout << " nsectors = " << m_nSectors
	   << G4endl;
    G4cout << "[DetectorConstruction][UpdateCalorSize] CHECK EE/CE-H boundaries (layer index counting from 0)" << G4endl;
    G4cout << " lastEElayer = " << lastEElayer_   << G4endl;
    G4cout << " firstMixedlayer = " << firstMixedlayer_   << G4endl;
    G4cout << " firstScintlayer = " << firstScintlayer_   << G4endl;
    G4cout << " firstCoarseScintlayer = " << firstCoarseScintlayer_   << G4endl;
    G4cout << "[DetectorConstruction][UpdateCalorSize] Eta ranges: " << G4endl;
    double tmpz=m_z0pos;
    G4cout << "check vector sizes: calostruct " << m_caloStruct.size() << " etamin " << m_minEta.size() << " etamax " << m_maxEta.size() << std::endl;
    for (unsigned il(0); il<m_caloStruct.size(); ++il){
      G4cout << "Layer #" << il << " z0= " << tmpz << " eta " << m_minEta[il] << " " << m_maxEta[il]<< G4endl;
      tmpz += m_caloStruct[il].Total_thick;
    }
  }
  else G4cout << "[DetectorConstruction][UpdateCalorSize] Z x XY = "
	      << m_CalorSizeZ << " x "
	      << m_CalorSizeXY << " mm "
	      << ", nsectors = " << m_nSectors
	      <<  G4endl;

}

//
G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //clean old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //world
  G4double expHall_z = model_ == DetectorConstruction::m_2016TB? 0.5*m : 14*m;
  G4double expHall_x = model_ == DetectorConstruction::m_2016TB? 20*cm : 3*m;
  G4double expHall_y = model_ == DetectorConstruction::m_2016TB? 20*cm : 3*m;

  G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume* experimentalHall_log = new G4LogicalVolume(experimentalHall_box, m_materials["Air"],"expHall_log");
  G4VPhysicalVolume* experimentalHall_phys
    = new G4PVPlacement(0,                       // no rotation
			G4ThreeVector(0.,0.,0.), // translation position
			experimentalHall_log,    // its logical volume
			"expHall",               // its name
			0,                       // its mother volume
			false,                   // no boolean operations
			0);                      // its copy number


  experimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  //detector's World
  G4double pos_x = 0.;
  G4double pos_y = 0.;
  G4double pos_z = 0.;
  if (model_ == DetectorConstruction::m_FULLSECTION){
    pos_x = 0.;
    pos_y = 0.;
    pos_z = m_z0pos+m_CalorSizeZ/2;
  }

  if (model_ == DetectorConstruction::m_FULLSECTION){
    m_solidWorld = new G4Tubs("Wbox",m_minRadius*0.9,m_maxRadius*1.1,m_WorldSizeZ/2,0,2*pi);
  }
  else {
    m_solidWorld = new G4Box("Wbox",m_WorldSizeXY/2,m_WorldSizeXY/2,m_WorldSizeZ/2);
  }
  cout << "m_WorldSizeXY = " << m_WorldSizeXY << ", m_WorldSizeZ = " << m_WorldSizeZ << endl;

  m_logicWorld = new G4LogicalVolume(m_solidWorld, m_materials["Air"], "Wlog");
  m_physWorld = new G4PVPlacement(0, G4ThreeVector(pos_x,pos_y,pos_z), m_logicWorld, "Wphys", experimentalHall_log, false, 0);

  for (unsigned iS(0); iS<m_nSectors; ++iS){
    G4double minL = m_sectorWidth*iS;
    buildSectorStack(iS,minL,m_sectorWidth-m_interSectorWidth);
    if (m_nSectors>1) fillInterSectorSpace(iS,minL+m_sectorWidth-m_interSectorWidth,m_interSectorWidth);
  }
  // Visualization attributes
  //
  m_logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

  //return m_physWorld;
  return experimentalHall_phys;
}

void DetectorConstruction::buildSectorStack(const unsigned sectorNum,
					    const G4double & minL,
					    const G4double & width)
{


  //build the stack
  G4double zOffset(-m_CalorSizeZ/2), zOverburden(0.), zOverburdenRef(0.);
  char nameBuf[100];
  G4VSolid *solid;
  G4VSolid *supportcone;

  G4double totalLengthX0 = 0;
  G4double totalLengthL0 = 0;

  for(size_t i=0; i<m_caloStruct.size(); i++)
    {
      G4double crackOffset = getCrackOffset(i);
      G4double angOffset = getAngOffset(i);

      //std::cout << " sector " << sectorNum << " layer " << i << " offset " << crackOffset  << std::endl;
      if (model_ == DetectorConstruction::m_FULLSECTION) {
	//for HF
	if (i==firstHFlayer_) {
	  zOffset+=m_z0HF-m_z0pos;
	  zOverburden=0;
	}
	crackOffset=0;
	if (i==firstMixedlayer_) zOverburdenRef = zOverburden;
	if (i==firstScintlayer_) zOverburden = zOverburdenRef;
      }
      const unsigned nEle = m_caloStruct[i].n_elements;
      //index for counting Si sensitive layers
      unsigned idx = 0;
      double totalThicknessLayer = 0;
      for (unsigned ie(0); ie<nEle;++ie){
	std::string eleName = m_caloStruct[i].ele_name[ie];
	if (m_nSectors==1) sprintf(nameBuf,"%s%d",eleName.c_str(),int(i+1));
	else sprintf(nameBuf,"%s%d_%d",eleName.c_str(),int(sectorNum),int(i+1));
	if (eleName=="Si") {
	  if (m_nSectors==1) sprintf(nameBuf,"Si%d_%d",int(i+1),idx);
	  else sprintf(nameBuf,"Si%d_%d_%d",int(sectorNum),int(i+1),idx);
	  idx++;
	}
	std::string baseName(nameBuf);
	G4double thick = m_caloStruct[i].ele_thick[ie];
	totalThicknessLayer += thick;
	//
	G4double extraWidth = 0;
	if (m_nSectors>1 && eleName=="W" && model_ != DetectorConstruction::m_FULLSECTION){
	 extraWidth = 10*mm;
	}
	if(thick>0){

#if 0
	  cout << "solid = constructSolid("<<baseName
	       <<",thick="<<thick
	       <<",zOffset+zOverburden="<<zOffset+zOverburden
	       <<",angOffset+minL="<<angOffset+minL
	       <<",width+extraWidth="<<width+extraWidth<<");"<<endl;
#endif
	  //special solid with holes
	  if (eleName=="CuExtra"){
	    solid = constructSolidWithHoles(baseName,thick,width+extraWidth);
	  }
	  else {
	    solid = constructSolid(i,baseName,thick,zOffset+zOverburden,angOffset+minL,width+extraWidth,i>=firstHFlayer_);
	  }

	  G4LogicalVolume *logi = new G4LogicalVolume(solid, m_materials[eleName], baseName+"log");
	  m_caloStruct[i].ele_X0[ie]   = m_materials[eleName]->GetRadlen();
	  m_caloStruct[i].ele_dEdx[ie] = m_dEdx[eleName];
	  m_caloStruct[i].ele_L0[ie]   = m_materials[eleName]->GetNuclearInterLength();
	  if (sectorNum==0 || sectorNum==m_nSectors-1) {
	    G4cout << "************ " << eleName;
	    if (m_nSectors>1) G4cout << " sector " << sectorNum;
	    //G4cout << " layer " << i << " dEdx=" << m_caloStruct[i].ele_dEdx[ie] << " X0=" << m_caloStruct[i].ele_X0[ie]
            //	   << " L0=" << m_caloStruct[i].ele_L0[ie] << " zpos=" ;
            
            G4cout << " layer " << i << " zpos=" ;
            if (i>=firstHFlayer_) G4cout << m_z0HF+zOverburden ;
            else G4cout << m_z0pos+zOverburden ;
            G4cout << "mm w=" << m_caloStruct[i].ele_thick[ie] << "mm";
            
	    //G4cout << " d=" << m_materials[eleName]->GetDensity();
	  //G4cout << G4endl;
	  //G4cout << *(m_materials[eleName]->GetMaterialTable()) << G4endl;
	    totalLengthX0 += m_caloStruct[i].ele_thick[ie]/m_caloStruct[i].ele_X0[ie]; G4cout << " TotX0=" << totalLengthX0;// << G4endl;
	    totalLengthL0 += m_caloStruct[i].ele_thick[ie]/m_caloStruct[i].ele_L0[ie]; G4cout << " TotLambda=" << totalLengthL0 << G4endl;
	  }

	  if (m_caloStruct[i].isSensitiveElement(ie)){
	    if (idx<=1) m_caloStruct[i].sensitiveZ = ((i>=firstHFlayer_)?m_z0HF+zOverburden:m_z0pos+zOverburden);

	    m_logicSi.push_back(logi);
	    //if (i==m_caloStruct.size()-1 && version_ == v_HGCALHF) m_logicSi.push_back(logi);
	  }

	  G4double xpvpos = -m_CalorSizeXY/2.+minL+width/2+crackOffset;
	  if (model_ == DetectorConstruction::m_FULLSECTION) xpvpos=0;
#if 0
	  cout << "m_caloStruct[i].ele_vol[nEle*sectorNum+ie]=new G4PVPlacement(0, G4ThreeVector(xpvpos="<<xpvpos
	       << ",0.,zOffset+zOverburden+thick/2="<<zOffset+zOverburden+thick/2
	       << "), logi,"
	       << baseName+"phys, m_logicWorld, false, 0);" << endl;
#endif
	  m_caloStruct[i].ele_vol[nEle*sectorNum+ie]=
	    new G4PVPlacement(0, G4ThreeVector(xpvpos,0.,zOffset+zOverburden+thick/2), logi, baseName+"phys", m_logicWorld, (eleName=="CuExtra")?true:false, 0);
	  //std::cout << " **** positionning layer " <<  m_caloStruct[i].ele_vol[nEle*sectorNum+ie]->GetName() << " at " << xpvpos << " 0 " << zOffset+zOverburden+thick/2 << std::endl;

	  G4VisAttributes *simpleBoxVisAtt= new G4VisAttributes(m_caloStruct[i].g4Colour(ie));
	  simpleBoxVisAtt->SetVisibility(true);
	  if (eleName=="CuExtra") {
	    simpleBoxVisAtt->SetForceSolid(true);
	    simpleBoxVisAtt->SetForceAuxEdgeVisible(true);
	  }
	  logi->SetVisAttributes(simpleBoxVisAtt);
	  zOverburden = zOverburden + thick;
	  //for sensitive volumes
	  //add region to be able to set specific cuts for it
	  //just for Si
	  if (eleName=="Si"){
	    unsigned nlogicsi = m_logicSi.size();
	    G4Region* aRegion = new G4Region(baseName+"Reg");
	    m_logicSi[nlogicsi-1]->SetRegion(aRegion);
	    aRegion->AddRootLogicalVolume(m_logicSi[nlogicsi-1]);
	  }
	}

      }//loop on elements

      //add support cone, for EE only, and only in full det version
      if (i<=lastEElayer_ && model_==DetectorConstruction::m_FULLSECTION && (version_ == v_HGCALEE_v6 || version_ ==  v_HGCAL_v6 || version_ == v_HGCALEE_v7 || version_ ==  v_HGCAL_v7 || version_ == v_HGCALEE_v8 || version_==v_HGCALEE_v8_neutmod || version_ ==  v_HGCAL_v8 || version_==v_HGCALEE_v8_air3 || version_==v_HGCALEE_v8_air4 || version_==v_HGCALEE_v8_Cu || version_==v_HGCALEE_v8_Cu_12 || version_ == v_HGCAL_v8_envelope || version_==v_HGCAL_v9 || version_==v_HGCALEE_v9 || version_==v_HGCAL_v10 || version_==v_HGCALEE_v10)) {
	//remove support cone for moderator
	//if (i==0) {
	//totalThicknessLayer -= 100;
	//}
	std::string eleName = "SupportCone";
	if (m_nSectors==1) sprintf(nameBuf,"%s%d",eleName.c_str(),int(i+1));
	else sprintf(nameBuf,"%s%d_%d",eleName.c_str(),int(sectorNum),int(i+1));

	std::string baseName(nameBuf);
	G4double extraWidth = 0;
#if 1
	cout << "solid = constructSupportCone("<<baseName
	     <<",thick="<<totalThicknessLayer
	     <<",zOffset+zOverburden="<<zOffset+zOverburden-totalThicknessLayer
	     <<",angOffset+minL="<<angOffset+minL
	     <<",width+extraWidth="<<width+extraWidth<<");"<<endl;
#endif

	supportcone = constructSupportCone(baseName,totalThicknessLayer,zOffset+zOverburden-totalThicknessLayer,angOffset+minL,width+extraWidth);
	G4LogicalVolume *logi = new G4LogicalVolume(supportcone, m_materials["Al"], baseName+"log");
	m_logicAl.push_back(logi);
	G4double xpvpos = -m_CalorSizeXY/2.+minL+width/2+crackOffset;
	if (model_ == DetectorConstruction::m_FULLSECTION) xpvpos=0;
	m_caloStruct[i].supportcone_vol=
	new G4PVPlacement(0, G4ThreeVector(xpvpos,0.,zOffset+zOverburden-totalThicknessLayer/2), logi, baseName+"phys", m_logicWorld, false, 0);
	G4VisAttributes *simpleBoxVisAtt= new G4VisAttributes(G4Colour::Red());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	unsigned nlogical = m_logicAl.size();
	G4Region* aRegion = new G4Region(baseName+"Reg");
	m_logicAl[nlogical-1]->SetRegion(aRegion);
	aRegion->AddRootLogicalVolume(m_logicAl[nlogical-1]);
      }
    }//loop on layers
  std::cout << " Z positions of sensitive layers: " << std::endl;
  for (size_t i=0; i<m_caloStruct.size(); i++) {
    std::cout << "sensitiveZ_[" << i << "] = " << m_caloStruct[i].sensitiveZ << ";"  << std::endl;
  }
  std::cout << " Eta boundary of sensitive layers for HGCSSDetector class: " << std::endl;
  for (size_t i=0; i<m_caloStruct.size(); i++) {
    //<< " maxRadius=" <<
    if (i<firstScintlayer_) std::cout << " etaBoundary_[" << i << "] =" << m_minEta[i] << ";"  << std::endl;
    else std::cout << " etaBoundary_[" << i << "] =" << m_maxEta[i] << ";"  << std::endl;
  }

  //dummy layer to get genparticles
  std::string eleName = "DummyLayer";
  double wDummy = 0.5;
  //-1 to ensure no overlap...
  G4VSolid *dummylayer;
  if (model_ == DetectorConstruction::m_FULLSECTION) 
    dummylayer = constructSolid(eleName,wDummy*2,-m_CalorSizeZ/2-wDummy*2-1,0,width,1.3,5);
  else dummylayer = constructSolid(0,eleName,wDummy*2,-m_CalorSizeZ/2-wDummy*2,0,m_CalorSizeXY);
  G4LogicalVolume *logi = new G4LogicalVolume(dummylayer, m_materials["Air"], eleName+"log");
  G4double xpvpos = 0;//-m_CalorSizeXY/2.;
  if (model_ == DetectorConstruction::m_FULLSECTION) xpvpos=0;
  m_caloStruct[0].dummylayer_vol=
    new G4PVPlacement(0, G4ThreeVector(xpvpos,0.,-m_CalorSizeZ/2-wDummy*2+wDummy-1), logi, eleName+"phys", m_logicWorld, false, 0);
  std::cout << "Adding dummy layer at (local coord): " << xpvpos << " 0 " << -m_CalorSizeZ/2-wDummy*2+wDummy-1 << std::endl;

  G4VisAttributes *simpleBoxVisAtt= new G4VisAttributes(G4Colour::Blue);
  simpleBoxVisAtt->SetVisibility(true);
  simpleBoxVisAtt->SetForceSolid(true);
  //simpleBoxVisAtt->SetForceAuxEdgeVisible(true);
  logi->SetVisAttributes(simpleBoxVisAtt);




}//buildstack

void DetectorConstruction::fillInterSectorSpace(const unsigned sectorNum,
						const G4double & minL,
						const G4double & width)
{

  //build the stack
  G4double zOffset(-m_CalorSizeZ/2), zOverburden(0.), zOverburdenRef(0.);
  char nameBuf[10];
  G4VSolid *solid;

  for(size_t i=0; i<m_caloStruct.size(); i++)
    {
      G4double crackOffset = getCrackOffset(i);
      G4double angOffset = getAngOffset(i);

      if (model_ == DetectorConstruction::m_FULLSECTION) {
	crackOffset=0;
	if (i==firstMixedlayer_) zOverburdenRef = zOverburden;
	if (i==firstScintlayer_) zOverburden = zOverburdenRef;
      }
      const unsigned nEle = m_caloStruct[i].n_elements;
      for (unsigned ie(0); ie<nEle;++ie){

	std::string eleName = m_caloStruct[i].ele_name[ie];
	G4double thick = m_caloStruct[i].ele_thick[ie];
	G4double extraWidth = 0;
	if (eleName=="W" && model_ != DetectorConstruction::m_FULLSECTION){
	 extraWidth = -10.*mm;
	 //std::cout << " -- total width: " << width+extraWidth << " offsets: " << crackOffset << " " << angOffset << std::endl;
	}
	eleName = "CFMix";
	sprintf(nameBuf,"%s%d_%d",eleName.c_str(),int(sectorNum),int(i+1));
	std::string baseName(nameBuf);
	if(thick>0){
	  solid = constructSolid(i,baseName,thick,zOffset+zOverburden,angOffset+minL,width+extraWidth);
	  G4LogicalVolume *logi = new G4LogicalVolume(solid, m_materials[eleName], baseName+"log");
	  G4double xpvpos = -m_CalorSizeXY/2.+minL+width/2+crackOffset;
	  if (model_ == DetectorConstruction::m_FULLSECTION) xpvpos=0;
	  G4PVPlacement *tmp = new G4PVPlacement(0, G4ThreeVector(xpvpos,0.,zOffset+zOverburden+thick/2), logi, baseName+"phys", m_logicWorld, false, 0);
	  //std::cout << "** positionning layer " << baseName << " at " << xpvpos << " 0 " << zOffset+zOverburden+thick/2 << std::endl;

	  G4VisAttributes *simpleBoxVisAtt= new G4VisAttributes(G4Colour::Magenta);
	  simpleBoxVisAtt->SetVisibility(true);
	  simpleBoxVisAtt->SetForceSolid(true);
	  logi->SetVisAttributes(simpleBoxVisAtt);
	  zOverburden = zOverburden + thick;
	}
      }//loop on elements
    }//loop on layers

}//fill intersector space

G4double DetectorConstruction::getCrackOffset(size_t layer){
  //model with 3 cracks identical by block of 10 layers
  //if (m_nSectors>1) return static_cast<unsigned>(layer/10.)*static_cast<unsigned>(m_sectorWidth/30.)*10;
  //with cracks shifted systematically layer-to-layer
  //if (m_nSectors>1) return 10*((7*layer)%31);
  //cracks shifted every two layers by 2cm
  //if (m_nSectors>1) return static_cast<unsigned>(layer/2.)*30;
  //cracks shifted every other layer by 2cm
  if (m_nSectors>1) return static_cast<unsigned>((layer%4)/2.)*30;

  return 0;
}

G4double DetectorConstruction::getAngOffset(size_t layer){
  if (model_ == DetectorConstruction::m_FULLSECTION){
    if (m_nSectors>1) return static_cast<unsigned>(layer/10.)*m_sectorWidth/3.;
    return 0;
  }
  return 0;
}
//
void DetectorConstruction::SetMagField(G4double fieldValue)
{

  if(fieldValue<=0) return;

  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if(m_magField) delete m_magField;                //delete the existing magn field
  m_magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
  fieldMgr->SetDetectorField(m_magField);
  fieldMgr->CreateChordFinder(m_magField);
  fieldMgr->SetDetectorField(m_magField);
}

void DetectorConstruction::SetDetModel(G4int model)
{
  if (model <= 0) return;
  std::cout << " -- Setting detector model to " << model << std::endl;
  model_ = model;
}

void DetectorConstruction::SetWThick(std::string thick)
{
  if (thick.size() <= 0) return;
  std::cout << " -- Setting W thick to " << thick << std::endl;
  std::vector<std::string> vec;
  boost::split(vec, thick, boost::is_any_of(","));
  absThickW_.resize(vec.size(),0);
  for (unsigned iE(0); iE<vec.size(); ++iE){//loop on elements
    std::istringstream(vec[iE])>>absThickW_[iE];
    std::cout << absThickW_[iE] << " ";
  }
  std::cout << std::endl;
}

void DetectorConstruction::SetPbThick(std::string thick)
{
  if (thick.size() <= 0) return;
  std::cout << " -- Setting Pb thick to " << thick << std::endl;
  std::vector<std::string> vec;
  boost::split(vec, thick, boost::is_any_of(","));
  absThickPb_.resize(vec.size(),0);
  for (unsigned iE(0); iE<vec.size(); ++iE){//loop on elements
    std::istringstream(vec[iE])>>absThickPb_[iE];
    std::cout << absThickPb_[iE] << " ";
  }
  std::cout << std::endl;
}

void DetectorConstruction::SetDropLayers(std::string layers)
{
  dropLayer_.resize(lastEElayer_+1,false);
  if (layers.size() <= 0) return;
  std::cout << " -- Dropping layers " << layers << std::endl;
  std::vector<std::string> vec;
  boost::split(vec, layers, boost::is_any_of(","));
  for (unsigned iE(0); iE<vec.size(); ++iE){//loop on elements
    unsigned layerId = 0;
    std::istringstream(vec[iE])>>layerId;
    if (layerId>0 && layerId<lastEElayer_+2) dropLayer_[layerId-1] = true;
    else std::cout << " -- invalid layer to drop, ignoring..." << std::endl;
  }
  for (unsigned iE(0); iE<dropLayer_.size(); ++iE){//loop on elements
    std::cout << dropLayer_[iE] << " ";
  }
  std::cout << std::endl;
}

G4SubtractionSolid *DetectorConstruction::constructSolidWithHoles (std::string baseName, G4double thick, const G4double & width){
  G4VSolid *solid = new G4Box(baseName+"box", width/2, m_CalorSizeXY/2, thick/2 );
  G4SubtractionSolid* result = 0;

  bool isSmall = false;
  if (!m_coarseGranularity) isSmall = true;

  //center of rectangles to remove
  double xc[34] = {
    -24.1,  -2.1, -24.1,  -2.1,     117.1,    -1.1,  57.3,  61.9, -22.1,  58.3,  75.9,  20.9,  16.9,  75.9,  -2.1,  29.3,  81.9,
    -107.9, -59.6,  24.1, -59.6,    100.1,    -26.5,-108.2, -98.0, -25.5,-109.2, -81.6, -66.9, 100.6, -79.6, -54.5, 113.0, -69.0
  };

  double yc[34] = {
    26.8,  26.8, -75.2, -75.2,   -24.2,    -17.2, -64.9,   9.5,  -6.2, -87.2,  22.2,  16.8, -62.9, -21.9,  11.8, -82.3,  -6.1,
    69.8,  26.8,  69.8, -75.2,   -24.2,    80.1, -64.9,   9.5,  56.4, -87.2,  22.2,  82.1,  82.1, -21.9,  62.7,  62.7,  -6.1
  };
  //rotation angles, in pi/3. units
  double rot[34] = {
    0.,0.,0.,0.,  0.,  0.,2.,4.,    0.,2.,4.,   0.,2.,4.,   0.,2.,4.,
    0.,0.,0.,0.,  0.,  2.,2.,4.,    2.,2.,4.,   2.,2.,4.,   2.,2.,4.
  };

  // 1 or 0 depending whether present in 200 and 300µm Si
  unsigned in23[34] = {
    1,0,1,0,   1,  1,1,1,  1,1,1,0,0,0,  1,1,1,
    1,0,1,0,   1,  1,1,1,  1,1,1,0,0,0,  1,1,1
  };

  // 1 or 0 depending on whether a complete hole or just a recess (hole for only 0.85mm)
  unsigned holeD[34] = {
    0,0,0,0,  1,  1,1,1,  0,0,0,0,0,0, 0,0,0,
    0,0,0,0,  1,  1,1,1,  0,0,0,0,0,0, 0,0,0
  };

  //depth of hole depending on type above
  double hDepth[2] = {0.85/1.5*thick,thick};

  unsigned type[34] = {
    0,0,0,0,  1,  2,2,2,  3,3,3,3,3,3, 4,4,4,
    0,0,0,0,  1,  2,2,2,  3,3,3,3,3,3, 4,4,4
  };
  // The hole i has dx, dy dimmensions: dx[type[i]], dy[type[i]]
  double dx[5] = {18.,17., 8., 9.,7.};
  double dy[5] = {18.,85.,32.,18.,7.};

  //replicas to cover the 1*1m^2 area: +/- 1.5 times in x, +/- 2 times in y
  double stepx = 251.2;
  double stepy = 193.3;
  int nrepx = 2;
  int nrepy = 2;
  bool isFirst = true;
  for (int irx(-1); irx<nrepx; irx++){
    for (int iry(-1); iry<nrepy; iry++){
      double newxc[34];
      double newyc[34];
      for (unsigned i(0); i<34; ++i){
	newxc[i] = xc[i]+irx*stepx;
	newyc[i] = yc[i]+iry*stepy;
	//exclude holes that would go beyond the initial volume, with safety margin for rotations...
	double maxxy = std::max(dx[i],dy[i]);
	if (newxc[i]+maxxy > m_CalorSizeXY/2. ||
	    newxc[i]-1.*maxxy < -1.*m_CalorSizeXY/2. ||
	    newyc[i]+maxxy > m_CalorSizeXY/2. ||
	    newyc[i]-1.*maxxy < -1.*m_CalorSizeXY/2.) continue;

	if (!isSmall && !in23[i]) continue;
	std::ostringstream lname;
	lname << baseName << "_" << i << "_" << irx+2 << "_" << iry+2;
	G4VSolid *hole = new G4Box(lname.str(), dx[type[i]]/2., dy[type[i]]/2., hDepth[holeD[i]]/2. );

	G4RotationMatrix zRot;   // Rotates X and Y axes only
	zRot.rotateZ(rot[i]*pi/3.*rad);
	G4ThreeVector  translation(newxc[i],newyc[i],-1.*thick/2.+0.5*hDepth[holeD[i]]);
	G4Transform3D transform(zRot,translation);

	G4SubtractionSolid *tmpSolid = new G4SubtractionSolid(lname.str()+"boxminushole",isFirst?&(*solid):&(*result),&(*hole),transform);
	result = tmpSolid;
	//solid = &tmpSolid;
	isFirst = false;
      }//loop on holes
    }//loop on y replicas
  }//loop on x replicas

  return result;
}

G4VSolid *DetectorConstruction::constructSolid (const unsigned layer, std::string baseName, G4double thick, G4double zpos,const G4double & minL, const G4double & width, const bool isHF){
  G4VSolid *solid;
  if (model_ == DetectorConstruction::m_FULLSECTION){
    if (isHF) solid=constructSolid(baseName,thick,zpos,minL,width,m_minEtaHF,m_maxEtaHF);
    else solid=constructSolid(baseName,thick,zpos,minL,width,m_minEta[layer],m_maxEta[layer]);
  }
  else if (model_ == DetectorConstruction::m_2016TB){
    std::cout << " zpos = " << zpos << "-" << zpos+thick << " hexagon with side " << m_maxRadius << std::endl;
    G4double zPlane[2] = { -thick/2,thick/2};
    G4double rInner[2] = { 0, 0 };
    G4double rOuter[2] = { m_minRadius, m_minRadius }; // feed center-to-side distance to G4Polyhedra
    solid = new G4Polyhedra(baseName+"hexa",0,2*pi,6,2,zPlane,rInner,rOuter); // startphi angle points to a corner
  }
  else{
    solid = new G4Box(baseName+"box", width/2, m_CalorSizeXY/2, thick/2 );
  }
  return solid;
}

G4VSolid *DetectorConstruction::constructSolid (std::string baseName, G4double thick, G4double zpos,const G4double & minL, const G4double & width, const double & etamin, const double & etamax){
  G4VSolid *solid=0;
  if (model_ == DetectorConstruction::m_FULLSECTION){
    double minR = tan(2*atan(exp(-1.*etamax)))*(zpos+m_z0pos+m_CalorSizeZ/2);
    double maxR = tan(2*atan(exp(-1.*etamin)))*(zpos+m_z0pos+m_CalorSizeZ/2);
    //std::cout << baseName << " zpos = " << zpos+m_z0pos+m_CalorSizeZ/2 << " radius range " << minR << " " << maxR << " thick " << thick << " width " << width << std::endl;
    solid = new G4Tubs(baseName+"box",minR,maxR,thick/2,minL,width);
  }
  return solid;
}

G4VSolid *DetectorConstruction::constructSupportCone (std::string baseName, G4double thick, G4double zpos,const G4double & minL, const G4double & width){

 return constructSupportCone(baseName,thick,zpos,minL,width,3.00);

}

G4VSolid *DetectorConstruction::constructSupportCone (std::string baseName, G4double thick, G4double zpos,const G4double & minL, const G4double & width, const double & etamin){

  G4VSolid *solid;
  double maxR = tan(2*atan(exp(-1.*etamin)))*(zpos+m_z0pos+m_CalorSizeZ/2);
  //double minR = tan(2*atan(exp(-1.*etamax)))*(zpos+m_z0pos+m_CalorSizeZ/2);
  //use fix radius but one edge always at eta=3
  double minR = maxR-10;
  std::cout << " SupportCone zpos = " << zpos+m_z0pos+m_CalorSizeZ/2 << " radius range " << minR << " " << maxR << std::endl;
  solid = new G4Tubs(baseName+"box",minR,maxR,thick/2,minL,width);

  return solid;
}

//
float DetectorConstruction::getRealisticRho(float z,std::string boundary){
  if(boundary=="outer") {
    if(z<3870) return 0.3198*z+495.5339;
    if(z>=3870 && z<4582) return 1.2598*z-3142.5478;
    if(z>=4582 && z<5056) return 2630.0000;
    if(z>=5056) return 0.0159*z+2409.7460;
  }
  if(boundary=="inner") {
    if(z<3662) return 325.0000;
    if(z>=3662 && z<4117) return 364.0000;
    if(z>=4117 && z<4369) return 412.0000;
    if(z>=4369) return 511.0000;
  }
  if(boundary=="mixed") {
    if(z<4407) return 1540.0000;
    if(z>=4407 && z<4550) return 1388.0000;
    if(z>=4550 && z<4885) return 1236.0000;
    if(z>=4885) return 1064.0000;
  }

  return 0;
}
