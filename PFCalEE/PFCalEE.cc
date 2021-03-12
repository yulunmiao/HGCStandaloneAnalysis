#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

void printHelp() {
  std::cout << "===========================================================================" << std::endl
            << "PFCalEE - a standalone simulation of HGCal-type of detectors" << std::endl
            << "---------------------------------------------------------------------------" << std::endl
            << "Execute like" << std::endl
            << "PFCalEE steer.mac [options]" << std::endl
            << "---------------------------------------------------------------------------" << std::endl
            << "The following options can be used" << std::endl
            << "\t--version - the detector version to use" << std::endl
            << "\t--model - the model (stack, endcap etc.)" << std::endl
            << "\t--eta - pseudo-rapidity" << std::endl
            << "\t--shape - cell shape" << std::endl
            << "\t--absThick{W,Pb} - csv list of the thicknesses for W and Pb absorbers" << std::endl
            << "\t--dropLayers - csv list of layers to drop" << std::endl
            << "\t--fineGranularity - use fine granularity cells" << std::endl
            << "\t--ultraFineGranularity - use ultra fine granularity cells" << std::endl
            << "\t--ui - do not run in batch mode" << std::endl
            << "===========================================================================" << std::endl << std::endl;
}


int main(int argc,char** argv)
{

  // Choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // User Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //int version=DetectorConstruction::v_HGCAL_2016TB;
  int version=73;
  //int version=DetectorConstruction::v_HGCALEE_TB;
  int model=DetectorConstruction::m_FULLSECTION;
  int coarseGranularity(1);
  int shape = 4;
  double eta=0;

  std::string absThickW="";//1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2";
  std::string absThickPb="";//1,1,1,1,1,2.1,2.1,2.1,2.1,2.1,4.4,4.4,4.4,4.4";
  std::string dropLayers="";
  bool batchMode(true);

  if (argc<2){
    printHelp();
    return -1;
  }

  //parse command line
  G4String fileName = argv[1];
  for(int i=2;i<argc;i++){
    std::string arg(argv[i]);
    if(arg.find("-h") !=std::string::npos)                       { printHelp(); return -1;} 
    else if(arg.find("--model")!=std::string::npos)              { sscanf(argv[i+1],"%d",&model); i++;}
    else if(arg.find("--version")!=std::string::npos )           { sscanf(argv[i+1],"%d",&version); i++;}
    else if(arg.find("--eta")!=std::string::npos )               { sscanf(argv[i+1],"%lf",&eta); i++;}
    else if(arg.find("--shape")!=std::string::npos )             { sscanf(argv[i+1],"%d",&shape); i++;}
    else if(arg.find("--absThickW")!=std::string::npos)          { absThickW=argv[i+1]; i++;}
    else if(arg.find("--absThickPb")!=std::string::npos)         { absThickPb=argv[i+1]; i++;}
    else if(arg.find("--dropLayers")!=std::string::npos)         { dropLayers=argv[i+1]; i++;}
    else if(arg.find("--fineGranularity")!=std::string::npos)    { coarseGranularity=0;} 
    else if(arg.find("--ultraFineGranularity")!=std::string::npos)    { coarseGranularity=-1;} 
    else if(arg.find("--ui")!=std::string::npos)                 { batchMode=false;} 
  }

  std::cout << "-- Running version=" << version << " model=" << model << " shape=" << shape << std::endl
            << "\teta=" << eta << " coarse granularity=" << coarseGranularity << std::endl
            << "\tabsThickW=" << absThickW << " absThickPb=" << absThickPb << " dropLayers=" << dropLayers << std::endl
            << "\tbatchMode=" << batchMode << std::endl;

  runManager->SetUserInitialization(new DetectorConstruction(version,model,shape,absThickW,absThickPb,dropLayers,coarseGranularity));
  runManager->SetUserInitialization(new PhysicsList);

  // Set user action classes
  runManager->SetUserAction(new PrimaryGeneratorAction(model,eta));
  runManager->SetUserAction(new RunAction);
  runManager->SetUserAction(new EventAction);
  runManager->SetUserAction(new SteppingAction);

  // Initialize G4 kernel
  runManager->Initialize();

  // Initialize visualization if required
  G4VisManager* visManager=0;
  if(!batchMode) {
    visManager = new G4VisExecutive;
    visManager->Initialize();
  }

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if(batchMode)
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
      
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute "+fileName);
      //if (ui->IsGUI())
      // UImanager->ApplyCommand("/control/execute gui.mac");
      ui->SessionStart();
      delete ui;
      
      delete visManager;
    }

  delete runManager;

  return 0;
}
