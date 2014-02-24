#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
//for hadronic processes
#include "QGSP_BERT.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{
#ifdef G4VIS_USE
  std::cout << " -- G4VIS_USE is set " << std::endl;
#else
  std::cout << " -- G4VIS_USE is not set " << std::endl;
#endif

  // Choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // User Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //int version=DetectorConstruction::v_UNIFORM;
  //int version=DetectorConstruction::v_PFCALEE;
  int version=DetectorConstruction::v_UNIFORM_03;
  //int version=DetectorConstruction::v_CALICE;
  //int version=DetectorConstruction::v_CALICE_Pb;
  //int version=DetectorConstruction::v_HGCAL;
  if(argc>2) version=atoi(argv[2]);

  std::cout << "-- Running version " << version << std::endl;

  runManager->SetUserInitialization(new DetectorConstruction(version));
  //runManager->SetUserInitialization(new PhysicsList);
  G4VUserPhysicsList* physics = new QGSP_BERT();
  runManager->SetUserInitialization(physics);


  // Set user action classes
  runManager->SetUserAction(new PrimaryGeneratorAction);
  runManager->SetUserAction(new RunAction);
  runManager->SetUserAction(new EventAction);
  runManager->SetUserAction(new SteppingAction);
  
  // Initialize G4 kernel
  runManager->Initialize();
  
  // Initialize visualization
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  /*
  G4TrajectoryDrawByParticleID* model = new G4TrajectoryDrawByParticleID;
  model->SetDefault("yellow");
  model->Set("e+", "green");
  model->Set("e-", "green");
  visManager->RegisterModel(model);
  visManager->SelectTrajectoryModel(model->Name());
  */
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4String fileName;
  if (argc>1) fileName = argv[1];
  if (argc!=1 && fileName != "vis")   // batch mode
    {    
      std::cout << " ====================================== " << std::endl
		<< " ========  Running batch mode ========= " << std::endl
		<< " ====================================== " << std::endl;
      G4String command = "/control/execute ";
      UImanager->ApplyCommand(command+fileName);
    }
  else
    {
      std::cout << " ====================================== " << std::endl
		<< " ====  Running interactive display ==== " << std::endl
		<< " ====================================== " << std::endl;
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac"); 
#endif
      if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");
      ui->SessionStart();
      delete ui;
#endif
    }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

