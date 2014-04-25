#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
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
DetectorConstruction::DetectorConstruction(G4int ver, G4int mod) : version_(ver), model_(mod), addPrePCB_(false)
{
  //radiation lengths: cf. http://pdg.lbl.gov/2012/AtomicNuclearProperties/
  //W 3.504 cm
  //Pb 5.612 cm
  //Cu 14.36 cm
  switch(version_)
    {
      //cf. http://arxiv.org/abs/0805.4833
    case v_CALICE:
      {
	G4cout << "[DetectorConstruction] starting v_CALICE (10x0.4+10x0.8+10x1.2)X_0 with Tungsten" << G4endl;
	for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(0.4*3.504*mm,0.0*mm,0.525*mm,1.0*mm,2.5*mm) );
	for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(0.8*3.504*mm,0.0*mm,0.525*mm,1.0*mm,2.5*mm) );
	for(int i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(1.2*3.504*mm,0.0*mm,0.525*mm,1.0*mm,2.5*mm) );
	break;
      }
    case v_HGCALEE_Si80: case v_HGCALEE_Si120: case v_HGCALEE_Si200: case v_HGCALEE_Si500: case v_HGCALEE_gap1: case  v_HGCALEE_CALICE: case v_HGCALEE_inverted: case v_HGCALEE_concept: case v_HGCALEE_W: case v_HGCALEE_gap4: case v_HGCAL:
      {
	float siWidth(0.200), gap(2),pad(1);
	float pb1(1.63), pb2(3.32), pb3(5.56), cu(3.0);
	int n1(10),n2(10),n3(10);
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
	if(version_==v_HGCALEE_concept || version_==v_HGCAL) m_caloStruct.push_back( SamplingSection(0.0*mm,0*mm,siWidth*mm,pad*mm,gap*mm) );
	for(int i=0; i<n1; i++)         m_caloStruct.push_back( SamplingSection(pb1*mm,cu*mm,siWidth*mm,pad*mm,gap*mm) );
	for(int i=0; i<n2; i++)         m_caloStruct.push_back( SamplingSection(pb2*mm,cu*mm,siWidth*mm,pad*mm,gap*mm) );
	for(int i=0; i<n3; i++)         m_caloStruct.push_back( SamplingSection(pb3*mm,cu*mm,siWidth*mm,pad*mm,gap*mm) );

	if (version_==v_HGCAL){
	  for(int i=0; i<12; i++) {
	    //add an intermediate Si layer to study resolution improvement
	    m_caloStruct.push_back( SamplingSection(26*mm,0.*mm,0.3*mm,0.*mm,0.*mm) );
	    m_caloStruct.push_back( SamplingSection(26*mm,3.0*mm,0.3*mm,2.0*mm,2.0*mm) );
	  }
	  //for(int i=0; i<12; i++) m_caloStruct.push_back( SamplingSection(52*mm,3.0*mm,0.3*mm,2.0*mm,2.0*mm) );
	  for(int i=0; i<9; i++) m_caloStruct.push_back( SamplingSection(78*mm,0.*mm,9*mm,0.*mm,0.*mm) );
	}
	break;
      }
    case v_HGCALHE:
      {
	//add HCAL
	for(int i=0; i<12; i++) {
	  //add an intermediate Si layer to study resolution improvement
	  m_caloStruct.push_back( SamplingSection(26*mm,0.*mm,0.3*mm,0.*mm,0.*mm) );
	  m_caloStruct.push_back( SamplingSection(26*mm,3.0*mm,0.3*mm,2.0*mm,2.0*mm) );
	}
	//for(int i=0; i<12; i++) m_caloStruct.push_back( SamplingSection(52*mm,3.0*mm,0.3*mm,2.0*mm,2.0*mm) );
	for(int i=0; i<9; i++) m_caloStruct.push_back( SamplingSection(78*mm,0.*mm,9*mm,0.*mm,0.*mm) );

	break;
      }
    case v_HGCALHEScint:
      {
	//for(int i=0; i<12; i++) m_caloStruct.push_back( SamplingSection(52*mm,3.0*mm,0.3*mm,2.0*mm,2.0*mm) );
	for(int i=0; i<9; i++) m_caloStruct.push_back( SamplingSection(78*mm,0.*mm,9*mm,0.*mm,0.*mm) );

	break;
      }
    case v_HGCALHE_CALICE:
      {
	//for(int i=0; i<12; i++) m_caloStruct.push_back( SamplingSection(52*mm,3.0*mm,0.3*mm,2.0*mm,2.0*mm) );
	//total 5.3 lambda = 22mm brass * 38 layers 
	for(int i=0; i<38; i++) m_caloStruct.push_back( SamplingSection(21*mm,0.*mm,5*mm,0.*mm,0.*mm) );
	for(int i=0; i<9; i++) m_caloStruct.push_back( SamplingSection(21*mm,0.*mm,5*mm,0.*mm,0.*mm) );
	for(int i=0; i<7; i++) m_caloStruct.push_back( SamplingSection(104*mm,0.*mm,5*mm,0.*mm,0.*mm) );
	break;
      }

    }

  DefineMaterials();
  SetMagField(0);
  m_detectorMessenger = new DetectorMessenger(this);
  UpdateCalorSize();
}

//
DetectorConstruction::~DetectorConstruction() { delete m_detectorMessenger;}

//
void DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();
  m_materials["Abs"] = (version_== v_CALICE || version_==v_HGCALEE_W) ? 
    nistManager->FindOrBuildMaterial("G4_W",false) :
    nistManager->FindOrBuildMaterial("G4_Pb",false);
  m_materials["Cu"] = nistManager->FindOrBuildMaterial("G4_Cu",false); 
  m_materials["Si"] = nistManager->FindOrBuildMaterial("G4_Si",false);
  m_materials["Zn"] = nistManager->FindOrBuildMaterial("G4_Zn",false);
  m_materials["PCB"] = new G4Material("G10",1.700*g/cm3,4);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(14), 1);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(8) , 2);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(6) , 3);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(1) , 3);
  m_materials["Air"]=nistManager->FindOrBuildMaterial("G4_AIR",false);
  m_materials["Brass"]= new G4Material("Brass",8.5*g/cm3,2);
  m_materials["Brass"]->AddMaterial(m_materials["Cu"]  , 70*perCent);
  m_materials["Brass"]->AddMaterial(m_materials["Zn"]  , 30*perCent);
  m_materials["AbsHCAL"] = (version_== v_HGCALHE_CALICE) ?
    nistManager->FindOrBuildMaterial("G4_Fe",false) :
    m_materials["Brass"];
  m_materials["C"] = nistManager->FindOrBuildMaterial("G4_C",false); 
  m_materials["H"] = nistManager->FindOrBuildMaterial("G4_H",false); 
  m_materials["Scintillator"]= nistManager->FindOrBuildMaterial("G4_POLYSTYRENE",false); 
  //m_materials["Scintillator"]= new G4Material("Scintillator",1.032*g/cm3,2);
  //m_materials["Scintillator"]->AddMaterial(m_materials["C"]  , 91.512109*perCent);
  //m_materials["Scintillator"]->AddMaterial(m_materials["H"]  , 8.4878906*perCent);
  G4cout << m_materials["Scintillator"] << G4endl;

}

//
void DetectorConstruction::UpdateCalorSize(){  

  m_CalorSizeZ=0;
  for(size_t i=0; i<m_caloStruct.size(); i++)
    m_CalorSizeZ=m_CalorSizeZ+m_caloStruct[i].Total_thick;

  if (model_ == DetectorConstruction::m_SIMPLE_20)
    m_CalorSizeXY=200;
  else if (model_ == DetectorConstruction::m_SIMPLE_50)
    m_CalorSizeXY=500;
  else if (model_ == DetectorConstruction::m_SIMPLE_100)
    m_CalorSizeXY=1000;
  else if (model_ == DetectorConstruction::m_FULLSECTION){
    m_CalorSizeXY=1700;
    m_minRadius = 150;
    m_maxRadius = m_CalorSizeXY;
  }
  else m_CalorSizeXY=200;

  m_WorldSizeZ=m_CalorSizeZ*1.1;  
  m_WorldSizeXY=m_CalorSizeXY*1.1;

  if (model_ == DetectorConstruction::m_FULLSECTION) 
    G4cout << "[DetectorConstruction][UpdateCalorSize] Z x minR * maxR = " 
	   << m_CalorSizeZ << " x " 
	   << m_minRadius << " x " 
	   << m_maxRadius 
	   << " mm " <<  G4endl;
  else G4cout << "[DetectorConstruction][UpdateCalorSize] Z x XY = " 
	      << m_CalorSizeZ << " x " 
	      << m_CalorSizeXY << " mm " <<  G4endl;

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
  G4double expHall_z = 6*m;
  G4double expHall_x = 2*m;
  G4double expHall_y = 2*m;

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

  //detector's World
  G4double pos_x = 0.;
  G4double pos_y = 0.;
  G4double pos_z = 0.;
  if (model_ == DetectorConstruction::m_FULLSECTION){
    pos_x = 0.;
    pos_y = 0.;
    pos_z = 3170*mm+m_CalorSizeZ/2;
  }

  if (model_ == DetectorConstruction::m_FULLSECTION){
    m_solidWorld = new G4Tubs("Wbox",m_minRadius*1.1,m_maxRadius*1.1,m_WorldSizeZ/2,0,2*pi);
  }
  else {
    m_solidWorld = new G4Box("Wbox",m_WorldSizeXY/2,m_WorldSizeXY/2,m_WorldSizeZ/2);
  }
  m_logicWorld = new G4LogicalVolume(m_solidWorld, m_materials["Air"], "Wlog");
  m_physWorld = new G4PVPlacement(0, G4ThreeVector(pos_x,pos_y,pos_z), m_logicWorld, "Wphys", experimentalHall_log, false, 0);



  //build the stack
  G4double zOffset(-m_CalorSizeZ/2), zOverburden(0.);
  char nameBuf[10];
  G4CSGSolid *solid;

  for(size_t i=0; i<m_caloStruct.size(); i++)
    {
      sprintf(nameBuf,"Abs%d",int(i+1)); 
      std::string baseName(nameBuf);
      G4double thick = m_caloStruct[i].Pb_thick;
      if(thick>0){
	solid = constructSolid(baseName,thick);
	G4LogicalVolume *logi;
	if (thick < 20) {
	  //ECAL
	  logi = new G4LogicalVolume(solid, m_materials["Abs"], baseName+"log");
	  m_caloStruct[i].Pb_X0 = m_materials["Abs"]->GetRadlen();
	  m_caloStruct[i].Pb_L0 = m_materials["Abs"]->GetNuclearInterLength();
	  G4cout << "************ Abs " << i << " " << m_caloStruct[i].Pb_X0 << " " << m_caloStruct[i].Pb_thick << " " << m_materials["Abs"]->GetDensity() << G4endl;
	}
	else {
	  //HCAL
	  logi = new G4LogicalVolume(solid, m_materials["AbsHCAL"], baseName+"log");
	  m_caloStruct[i].Pb_X0 = m_materials["AbsHCAL"]->GetRadlen();
	  m_caloStruct[i].Pb_L0 = m_materials["AbsHCAL"]->GetNuclearInterLength();
	  G4cout << "************ Abs " << i << " " << m_caloStruct[i].Pb_X0 << " " << m_caloStruct[i].Pb_thick << " " << m_materials["AbsHCAL"]->GetDensity() << G4endl;
	}
	m_caloStruct[i].Pb_vol= new G4PVPlacement(0, G4ThreeVector(0.,0.,zOffset+zOverburden+thick/2), logi, baseName+"phys", m_logicWorld, false, 0);
	G4VisAttributes *simpleBoxVisAtt= new G4VisAttributes(G4Colour::Gray());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	zOverburden = zOverburden + thick;
      }
      
      sprintf(nameBuf,"Cu%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].Cu_thick;
      if(thick>0){
	solid = constructSolid(baseName,thick);
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["Cu"], baseName+"log");
	//m_logicAbs.push_back(logi);
	m_caloStruct[i].Cu_vol = new G4PVPlacement(0, G4ThreeVector(0.,0.,zOffset+zOverburden+thick/2), logi, baseName+"phys", m_logicWorld, false, 0);
	m_caloStruct[i].Cu_X0 = m_materials["Cu"]->GetRadlen();
	m_caloStruct[i].Cu_L0 = m_materials["Cu"]->GetNuclearInterLength();
	G4cout << "************ Cu " << i << " " << m_caloStruct[i].Cu_X0 << " " << m_caloStruct[i].Cu_thick << " " << m_materials["Cu"]->GetDensity() << G4endl;
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Black());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	zOverburden = zOverburden + thick;
	//add region to be able to set specific cuts for it
	//G4Region* aRegion = new G4Region(baseName+"Reg");
	//m_logicAbs[i]->SetRegion(aRegion);
	//aRegion->AddRootLogicalVolume(m_logicAbs[i]);
      }

      //add PCB to shield from delta-rays?
      if(addPrePCB_)
	{
	  sprintf(nameBuf,"PCB%d",int(i+1));
	  baseName=nameBuf;
	  thick = m_caloStruct[i].PCB_thick;
	  if(thick>0){
	    solid = constructSolid(baseName,thick);
	    G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["PCB"], baseName+"log");
	    m_caloStruct[i].PCB_vol = new G4PVPlacement(0, G4ThreeVector(0.,0.,zOffset+zOverburden+thick/2), logi, baseName+"phys", m_logicWorld, false, 0);
	    m_caloStruct[i].PCB_X0 = m_materials["PCB"]->GetRadlen();
	    G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Blue());
	    simpleBoxVisAtt->SetVisibility(true);
	    logi->SetVisAttributes(simpleBoxVisAtt);
	    zOverburden = zOverburden + thick;
	  }
	}

      thick = m_caloStruct[i].Si_thick;
      if(thick>0){
	G4LogicalVolume *logi;
	if (thick<1){//Si
	  sprintf(nameBuf,"Si%d",int(i+1)); 
	  baseName=nameBuf;
	  solid = constructSolid(baseName,thick);
	  logi  = new G4LogicalVolume(solid, m_materials["Si"], baseName+"log");
	  m_caloStruct[i].Si_X0 = m_materials["Si"]->GetRadlen();
	  G4cout << "************ Si  " << i << " " << m_caloStruct[i].Si_X0 << " " << m_caloStruct[i].Si_thick << G4endl;
	}
	else {//scint
	  sprintf(nameBuf,"Scint%d",int(i+1));
	  baseName=nameBuf;
	  solid = constructSolid(baseName,thick);
	  logi  = new G4LogicalVolume(solid, m_materials["Scintillator"], baseName+"log");
	  m_caloStruct[i].Si_X0 = m_materials["Scintillator"]->GetRadlen();
	  G4cout << "************ Scintillator  " << i << " " << m_caloStruct[i].Si_X0 << " " << m_caloStruct[i].Si_thick << " " << m_materials["Scintillator"]->GetDensity() << G4endl;
	}
	m_logicSi.push_back(logi);
	m_caloStruct[i].Si_vol = new G4PVPlacement(0, G4ThreeVector(0.,0.,zOffset+zOverburden+thick/2), logi, baseName+"phys", m_logicWorld, false, 0);
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::White());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	zOverburden = zOverburden + thick;
	//add region to be able to set specific cuts for it
	G4Region* aRegion = new G4Region(baseName+"Reg");
	m_logicSi[i]->SetRegion(aRegion);
	aRegion->AddRootLogicalVolume(m_logicSi[i]);
      }


      sprintf(nameBuf,"PCB%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].PCB_thick;
      if(thick>0){
	solid = constructSolid(baseName,thick);
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["PCB"], baseName+"log");
	m_caloStruct[i].PCB_vol = new G4PVPlacement(0, G4ThreeVector(0.,0.,zOffset+zOverburden+thick/2), logi, baseName+"phys", m_logicWorld, false, 0);
	m_caloStruct[i].PCB_X0 = m_materials["PCB"]->GetRadlen();
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Blue());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	zOverburden = zOverburden + thick;
      }

      sprintf(nameBuf,"Air%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].Air_thick;
      if(thick>0){
	solid = constructSolid(baseName,thick);
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["Air"], baseName+"log");
	m_caloStruct[i].Air_vol = new G4PVPlacement(0, G4ThreeVector(0.,0.,zOffset+zOverburden+thick/2), logi, baseName+"phys", m_logicWorld, false, 0);
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Cyan());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	zOverburden = zOverburden + thick;
      }
    }
 
  //                                        
  // Visualization attributes
  //
  m_logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  //return m_physWorld;
  return experimentalHall_phys;
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

G4CSGSolid *DetectorConstruction::constructSolid (std::string baseName, G4double thick){
  
  G4CSGSolid *solid;
  if (model_ == DetectorConstruction::m_FULLSECTION){
    solid = new G4Tubs(baseName+"box",m_minRadius,m_maxRadius,thick/2,0,2*pi);
  }
  else{
    solid = new G4Box(baseName+"box", m_CalorSizeXY/2, m_CalorSizeXY/2, thick/2 );
  }
  return solid;
}
