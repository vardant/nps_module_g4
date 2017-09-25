#include "npsStackingAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"



npsStackingAction::npsStackingAction()
: gammaCounter(0)
{}



npsStackingAction::~npsStackingAction()
{}



G4ClassificationOfNewTrack
npsStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    if(aTrack->GetParentID()>0)
    { // particle is secondary
      gammaCounter++;
    }
  }
  return fUrgent;
}



void npsStackingAction::NewStage()
{
  //  G4cout << "Number of optical photons produced in this event : "
  //  << gammaCounter << G4endl;
  G4cout << gammaCounter << " ";
}



void npsStackingAction::PrepareNewEvent()
{ gammaCounter = 0; }
