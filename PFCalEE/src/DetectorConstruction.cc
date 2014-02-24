#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
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

using namespace std;

//
DetectorConstruction::DetectorConstruction(G4int ver) : version_(ver), addPrePCB_(false)
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
	G4cout << "[DetectorConstruction] starting v_HGCALEE with Si width=" << siWidth << " and gap " << gap << G4endl;
	
	//add 1um silicon to track incoming particles
	//if(version_==v_HGCALEE_SiDummy) m_caloStruct.push_back( SamplingSection(0,0,0.001*mm,0,0) );
	if(version_==v_HGCALEE_concept) m_caloStruct.push_back( SamplingSection(0.0*mm,cu*mm,siWidth*mm,pad*mm,gap*mm) );
	for(int i=0; i<n1; i++)         m_caloStruct.push_back( SamplingSection(pb1*mm,cu*mm,siWidth*mm,pad*mm,gap*mm) );
	for(int i=0; i<n2; i++)         m_caloStruct.push_back( SamplingSection(pb2*mm,cu*mm,siWidth*mm,pad*mm,gap*mm) );
	for(int i=0; i<n3; i++)         m_caloStruct.push_back( SamplingSection(pb3*mm,cu*mm,siWidth*mm,pad*mm,gap*mm) );

	if (version_==v_HGCAL){
	  //add HCAL
	  for(int i=0; i<12; i++) m_caloStruct.push_back( SamplingSection(52*mm,3.0*mm,0.3*mm,2.0*mm,2.0*mm) );
	  for(int i=0; i<18; i++) m_caloStruct.push_back( SamplingSection(34.5*mm,0*mm,0.3*mm,2.0*mm,2.0*mm) );
	}
	break;
      }
    }

  DefineMaterials();
  UpdateCalorSize();
  SetMagField(0);
  m_detectorMessenger = new DetectorMessenger(this);
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
}

//
void DetectorConstruction::UpdateCalorSize(){  

  m_CalorSizeX=0;
  for(size_t i=0; i<m_caloStruct.size(); i++)
    m_CalorSizeX=m_CalorSizeX+m_caloStruct[i].Total_thick;
  m_CalorSizeYZ=200;
  if (version_==v_HGCAL) m_CalorSizeYZ=500;

  m_WorldSizeX=m_CalorSizeX*1.1;  
  m_WorldSizeYZ=m_CalorSizeYZ*1.1;

  G4cout << "[DetectorConstruction][UpdateCalorSize] X x YZ = " << m_CalorSizeX << " x " << m_CalorSizeYZ << " mm " <<  G4endl;
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
  G4double expHall_x = 5*m;
  G4double expHall_y = 2*m;
  G4double expHall_z = 2*m;
  
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
  m_solidWorld = new G4Box("Wbox",m_WorldSizeX/2,m_WorldSizeYZ/2,m_WorldSizeYZ/2);
  m_logicWorld = new G4LogicalVolume(m_solidWorld, m_materials["Air"], "Wlog");
  m_physWorld = new G4PVPlacement(0, G4ThreeVector(), m_logicWorld, "Wphys", experimentalHall_log, false, 0);



  //build the stack
  char nameBuf[10];
  G4double xOffset(-m_CalorSizeX/2), xOverburden(0.);

  for(size_t i=0; i<m_caloStruct.size(); i++)
    {
      sprintf(nameBuf,"Abs%d",int(i+1)); 
      std::string baseName(nameBuf);
      G4double thick = m_caloStruct[i].Pb_thick;
      if(thick>0){
	G4Box *solid          = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi;
	if (thick < 20) {
	  //ECAL
	  logi = new G4LogicalVolume(solid, m_materials["Abs"], baseName+"log");
	  m_caloStruct[i].Pb_X0 = m_materials["Abs"]->GetRadlen();
	  G4cout << "************ Abs " << i << " " << m_caloStruct[i].Pb_X0 << " " << m_caloStruct[i].Pb_thick << " " << m_materials["Abs"]->GetDensity() << G4endl;
	}
	else {
	  //HCAL
	  logi = new G4LogicalVolume(solid, m_materials["Brass"], baseName+"log");
	  m_caloStruct[i].Pb_X0 = m_materials["Brass"]->GetRadlen();
	  G4cout << "************ Abs " << i << " " << m_caloStruct[i].Pb_X0 << " " << m_caloStruct[i].Pb_thick << " " << m_materials["Brass"]->GetDensity() << G4endl;
	}
	m_caloStruct[i].Pb_vol= new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	G4VisAttributes *simpleBoxVisAtt= new G4VisAttributes(G4Colour::Gray());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
      }
      
      sprintf(nameBuf,"Cu%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].Cu_thick;
      if(thick>0){
	G4Box *solid = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["Cu"], baseName+"log");
	m_caloStruct[i].Cu_vol = new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	m_caloStruct[i].Cu_X0 = m_materials["Cu"]->GetRadlen();
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Black());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
      }

      //add PCB to shield from delta-rays?
      if(addPrePCB_)
	{
	  sprintf(nameBuf,"PCB%d",int(i+1));
	  baseName=nameBuf;
	  thick = m_caloStruct[i].PCB_thick;
	  if(thick>0){
	    G4Box *solid = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	    G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["PCB"], baseName+"log");
	    m_caloStruct[i].PCB_vol = new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	    m_caloStruct[i].PCB_X0 = m_materials["PCB"]->GetRadlen();
	    G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Blue());
	    simpleBoxVisAtt->SetVisibility(true);
	    logi->SetVisAttributes(simpleBoxVisAtt);
	    xOverburden = xOverburden + thick;
	  }
	}


      sprintf(nameBuf,"Si%d",int(i+1)); 

      baseName=nameBuf;
      thick = m_caloStruct[i].Si_thick;
      if(thick>0){
	G4Box *solid = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["Si"], baseName+"log");
	m_caloStruct[i].Si_vol = new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	m_caloStruct[i].Si_X0 = m_materials["Si"]->GetRadlen();
	G4cout << "************ Si  " << i << " " << m_caloStruct[i].Si_X0 << " " << m_caloStruct[i].Si_thick << G4endl;
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::White());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
      }

      sprintf(nameBuf,"PCB%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].PCB_thick;
      if(thick>0){
	G4Box *solid = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["PCB"], baseName+"log");
	m_caloStruct[i].PCB_vol = new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	m_caloStruct[i].PCB_X0 = m_materials["PCB"]->GetRadlen();
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Blue());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
      }

      sprintf(nameBuf,"Air%d",int(i+1)); 
      baseName=nameBuf;
      thick = m_caloStruct[i].Air_thick;
      if(thick>0){
	G4Box *solid = new G4Box(baseName+"box", thick/2, m_CalorSizeYZ/2, m_CalorSizeYZ/2 );
	G4LogicalVolume *logi  = new G4LogicalVolume(solid, m_materials["Air"], baseName+"log");
	m_caloStruct[i].Air_vol = new G4PVPlacement(0, G4ThreeVector(xOffset+xOverburden+thick/2,0,0), logi, baseName+"phys", m_logicWorld, false, 0);
	G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour::Cyan());
	simpleBoxVisAtt->SetVisibility(true);
	logi->SetVisAttributes(simpleBoxVisAtt);
	xOverburden = xOverburden + thick;
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

  //apply a global uniform magnetic field along X axis
  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if(m_magField) delete m_magField;                //delete the existing magn field
  m_magField = new G4UniformMagField(G4ThreeVector(fieldValue,0.,0.));
  fieldMgr->SetDetectorField(m_magField);
  fieldMgr->CreateChordFinder(m_magField);
  fieldMgr->SetDetectorField(m_magField);  
}

