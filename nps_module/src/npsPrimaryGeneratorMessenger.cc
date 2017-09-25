#include "npsPrimaryGeneratorMessenger.hh"

#include "npsPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"



npsPrimaryGeneratorMessenger::npsPrimaryGeneratorMessenger(
                                          npsPrimaryGeneratorAction* npsGun)
:npsAction(npsGun)
{
  gunDir = new G4UIdirectory("/nps/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");
   
  //  polarCmd = new G4UIcmdWithADoubleAndUnit("/nps/gun/optPhotonPolar",this);
  //  polarCmd->SetGuidance("Set linear polarization");
  //  polarCmd->SetGuidance("  angle w.r.t. (k,n) plane");
  //  polarCmd->SetParameterName("angle",true);
  //  polarCmd->SetUnitCategory("Angle");  
  //  polarCmd->SetDefaultValue(-360.0);
  //  polarCmd->SetDefaultUnit("deg");
  //  polarCmd->AvailableForStates(G4State_Idle);
}



npsPrimaryGeneratorMessenger::~npsPrimaryGeneratorMessenger()
{
  //  delete polarCmd;
  delete gunDir;
}



void npsPrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{ 
  //  if( command == polarCmd ) {
  //      G4double angle = polarCmd->GetNewDoubleValue(newValue);
  //      if ( angle == -360.0*deg ) {
  //         npsAction->SetOptPhotonPolar();
  //      } else {
  //         npsAction->SetOptPhotonPolar(angle);
  //      }
  //  }
}
