#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4ios.hh"

#include "npsDetectorConstruction.hh"
#include "npsAddOptics.hh"
#include "QGSP_BERT.hh"
#include "npsPrimaryGeneratorAction.hh"
#include "npsRunAction.hh"
#include "npsEventAction.hh"
#include "npsStackingAction.hh"
#include "npsSteppingAction.hh"
#include "npsSteppingVerbose.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <ostream>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <error.h>
#include <errno.h>

// Random seed,copied from Simon. ---------------------------------------------

int getSeed();

int getSeed()
{
    // open /dev/random
    int fd;
    if ((fd = open("/dev/random", O_RDONLY)) == -1)
    {
        error(1, errno, "Fatal: unable to open /dev/random for reading.");
    }
    // read it
    int32_t value;
    ssize_t len;

    switch (len = read(fd, &value, sizeof(int32_t)))
    {
        case sizeof(int):
            printf("The seed is %i\n", value);
            break;
        case -1:
            error(1, errno, "Fatal: unable to read from /dev/random.");
            break;
        default:
            error(1, 0, "Fatal: %i bytes was read from "
                  "/dev/random instead of %i.", len, sizeof(int32_t));
    }
    // close file
    if (close(fd) == -1)
    {
        error(0, errno, "Warning: unable to close /dev/random properly.");
    }
    return value;
}

// ----------------------------------------------------------------------------

int main(int argc,char** argv)
{
  // Seed the random number generator
  //
  G4long myseed = getSeed();
  CLHEP::HepRandom::setTheSeed(myseed);
  
  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new npsSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // Run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // UserInitialization classes - mandatory
  //

  npsDetectorConstruction* detector = new npsDetectorConstruction;
  runManager-> SetUserInitialization(detector);
  //

  QGSP_BERT * physics = new QGSP_BERT ;


  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  physics->SetCutValue(3.5*cm, "gamma");
  physics->SetCutValue(0.085*mm, "e-");
  physics->SetCutValue(0.085*mm, "e+");

  physics->RegisterPhysics(new npsAddOptics("Cerenkov & Scintillation") );

  runManager-> SetUserInitialization(physics);

  physics->DumpCutValuesTable();

  //  G4VUserPhysicsList* physics = new npsPhysicsList;
  //  runManager-> SetUserInitialization(physics);

#ifdef G4VIS_USE
  // visualization manager
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // UserAction classes
  //

  npsRunAction* run_action = new npsRunAction;
  runManager->SetUserAction(run_action);
  //
  G4VUserPrimaryGeneratorAction* gen_action = new npsPrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);
  //
  npsEventAction* event_action = new npsEventAction(run_action);
  runManager->SetUserAction(event_action);
  //
  G4UserStackingAction* stacking_action = new npsStackingAction;
  runManager->SetUserAction(stacking_action);
  //
  G4UserSteppingAction* stepping_action =
                    new npsSteppingAction(detector, event_action);
  runManager->SetUserAction(stepping_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
    
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UI = G4UImanager::GetUIpointer(); 
   
  if (argc==1)   // Define UI session for interactive mode
    {
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
     UI->ApplyCommand("/control/execute vis.mac");
#endif
     ui->SessionStart();
     delete ui;
#endif

    }
   
  else         // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }
   
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  delete verbosity;

  return 0;
}
