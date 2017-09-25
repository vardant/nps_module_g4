#ifndef npsPrimaryGeneratorMessenger_h
#define npsPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class npsPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;



class npsPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    npsPrimaryGeneratorMessenger(npsPrimaryGeneratorAction*);
   ~npsPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    npsPrimaryGeneratorAction* npsAction;
    G4UIdirectory*               gunDir; 
    G4UIcmdWithADoubleAndUnit*   polarCmd;
};



#endif

