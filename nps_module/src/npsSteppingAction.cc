#include "npsSteppingAction.hh"

#include "npsDetectorConstruction.hh"
#include "npsEventAction.hh"

#include "G4Step.hh"
#include "G4UnitsTable.hh"

#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"

npsSteppingAction::npsSteppingAction(npsDetectorConstruction* det,
                                         npsEventAction* evt)
:detector(det), eventaction(evt)					 
{
  //  G4cout << " ==> npsSteppingAction::npsSteppingAction:" << G4endl;
}


npsSteppingAction::~npsSteppingAction()
{ }


void npsSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //  G4cout << " ==> npsSteppingAction::UserSteppingAction:" << G4endl;

  // ... retrieve the 'pre-step' point
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();

  // ... retrieve a touchable handle and access to the information
  G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
 
  G4VPhysicalVolume* volume = theTouchable->GetVolume();

  //Accrue energy deposition and step length of charged particle
  //in the scintillator block.

  //  if (volume == detector->Get_block()) {
  if (volume->GetName() == "Block_phys") {

    // collect energy and track length step by step
    G4double edep = aStep->GetTotalEnergyDeposit();
  
    G4double stepl = 0.;
    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
      stepl = aStep->GetStepLength();

    eventaction->AddEv(edep,stepl);

    //    if (edep !=0. || stepl != 0.) {
    //      G4int counterNo = theTouchable->GetCopyNumber(2);
    //      G4cout << "counter# " << counterNo << "  ";
    //      G4cout << "   Edep =" << std::setw(6) << G4BestUnit(edep,"Energy") << 
    //	"  stepL =" << G4BestUnit(stepl,"Length") << G4endl;
    //    }

  }

  //PMT hits. Look for boundary process with optical photon. The boundary
  //status must be "detection". The only volume on which it can happen is
  // photocathode.

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4OpBoundaryProcess* boundary=NULL;
  
  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm 
      = aStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }


  G4ParticleDefinition* particleType = aStep->GetTrack()->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){
    //Optical photon only

    //    G4VPhysicalVolume* prevol = 
    //      aStep->GetPreStepPoint()->GetPhysicalVolume();
    //    if (prevol->GetName() == "World") 
    //      aStep->GetTrack()->SetTrackStatus(fStopAndKill);

    boundaryStatus=boundary->GetStatus();
 
    //Check to see if the partcile was actually at a boundary
    //Otherwise the boundary status may not be valid
    //Prior to Geant4.6.0-p1 this would not have been enough to check
    if(aStep->GetPostStepPoint()->GetStepStatus()==fGeomBoundary){

      //      G4VPhysicalVolume* postvol = 
      //      	aStep->GetPostStepPoint()->GetPhysicalVolume();
      //      if (postvol->GetName() == "World") 
      //      	aStep->GetTrack()->SetTrackStatus(fStopAndKill);

      if (boundaryStatus == Detection) {

	//	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
	//	G4TouchableHandle theTouchable = postStepPoint->GetTouchableHandle();

	//	G4VPhysicalVolume* postvol = 
	//	aStep->GetPostStepPoint()->GetPhysicalVolume();
	//	G4VPhysicalVolume* prevol = 
	//		aStep->GetPreStepPoint()->GetPhysicalVolume();
	//G4cout << "SteppingAction: boundaryStatus = Detection" << G4endl;
	//	G4cout << "postvol = " << postvol->GetName()
	//	       << "  prevol = " << prevol->GetName() << G4endl;
	//	getchar();
	
	eventaction->AddNpe();

	////	aStep->GetTrack()->SetTrackStatus(fStopAndKill);
      }

      // if (postvol->GetName() != "Opt_Ins" && postvol->GetName() != "Air_phys"
      //	  && postvol->GetName() != "Block_phys")
      //	{
      //
      //    G4cout << "==> At the boundary of " << volume->GetName() << " - "
      //	  	    	 << postvol->GetName() << G4endl;
      //
      //	  switch(boundaryStatus){
      //	  case Absorption:
      //    G4cout << "    Optical photon absorped at boundary !" << G4endl;
      //	    break;
      //  case Detection: //Note, this assumes that the volume causing detection
      //              //is the photocathode because it is the only one with
      //                      //non-zero efficiency
      //	    {
      //	      G4cout << "    Optical photon Detected !" << G4endl;
      //	      break;
      //	    }
      //	  case FresnelReflection:
      //    G4cout << "    Optical photon Fresnel Reflection !" << G4endl;
      //	    break;
      //	  case TotalInternalReflection:
      //	    G4cout << "    Optical photon total internal reflection !" 
      //	    	   << G4endl;
      //	    break;
      //	  case SpikeReflection:
      //	    G4cout << "    Optical photon spike reflection !" << G4endl;
      //	    break;
      //	  default:
      //	    G4cout << "    Optical photon: boundary status not known !" 
      //	    	   << G4endl;
      //	    break;
      //	  }
      //	}
    }
  }
  //example of saving random number seed of this event, under condition
  //// if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent(); 
}
