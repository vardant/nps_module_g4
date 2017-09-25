// Add additional physics to the Geant4 standard physics lists. See
// http://hypernews.slac.stanford.edu/HyperNews/geant4/get/phys-list/54.html.

#ifndef npsAddOptics_hh
#define npsAddOptics_hh

#include "G4VPhysicsConstructor.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"

class npsAddOptics : public G4VPhysicsConstructor {

public:

  npsAddOptics(G4String aS) : G4VPhysicsConstructor(aS) {};

  virtual void ConstructParticle() 
  {
    // here we do nothing - or may be optical photons?
    G4OpticalPhoton::OpticalPhotonDefinition();
  };

  virtual void ConstructProcess() 
  {
    // here we can construct and register Cherenkov and other stuff.
    // Make sure nothing is there twice, though !!!

    theCerenkovProcess=new G4Cerenkov();
    theScintillationProcess=new G4Scintillation();
    theAbsorptionProcess=new G4OpAbsorption();
    theRayleighScattering=new G4OpRayleigh();
    theBoundaryProcess=new G4OpBoundaryProcess();

    G4ProcessManager * pManager = 0;
    pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
      G4ParticleDefinition* particle = theParticleIterator->value();
      pManager = particle->GetProcessManager();
      if(theCerenkovProcess->IsApplicable(*particle)){
	pManager->AddProcess(theCerenkovProcess);
	pManager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
      }
      if(theScintillationProcess->IsApplicable(*particle)){
      	pManager->AddProcess(theScintillationProcess);
	pManager->SetProcessOrdering(theScintillationProcess,idxPostStep);
      }
      if (particle->GetParticleName() == "opticalphoton") {
	G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
	pManager->AddDiscreteProcess(theAbsorptionProcess);
	pManager->AddDiscreteProcess(theRayleighScattering);
	//	pmanager->AddDiscreteProcess(fMieHGScatteringProcess);
	pManager->AddDiscreteProcess(theBoundaryProcess);
      }

    }

    theCerenkovProcess->SetMaxNumPhotonsPerStep(300);
    theCerenkovProcess->SetTrackSecondariesFirst(true);

    theScintillationProcess->SetTrackSecondariesFirst(true);

  }

private:

  G4Cerenkov* theCerenkovProcess;
  G4Scintillation* theScintillationProcess;
  G4OpAbsorption* theAbsorptionProcess;
  G4OpRayleigh* theRayleighScattering;
  G4OpBoundaryProcess* theBoundaryProcess;

}; 

#endif
