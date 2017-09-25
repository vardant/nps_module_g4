#ifndef npsSteppingAction_h
#define npsSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class npsDetectorConstruction;
class npsEventAction;

class npsSteppingAction : public G4UserSteppingAction
{
public:
  npsSteppingAction(npsDetectorConstruction*, npsEventAction*);
  virtual ~npsSteppingAction();

  void UserSteppingAction(const G4Step*);
    
private:
  npsDetectorConstruction* detector;
  npsEventAction*          eventaction;  
};

#endif
