// Make this appear first!
#include "G4Timer.hh"

#include "npsRunAction.hh"

#include "G4Run.hh"



npsRunAction::npsRunAction()
{
  timer = new G4Timer;
}



npsRunAction::~npsRunAction()
{
  delete timer;
}



void npsRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl; 
  timer->Start();
}



void npsRunAction::EndOfRunAction(const G4Run* aRun)
{   
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() 
         << " " << *timer << G4endl;
}
