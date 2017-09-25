#include "npsEventAction.hh"

#include "npsRunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//#include "npsPrimaryGeneratorAction.hh"

G4int const printModulo = 1;

npsEventAction::npsEventAction(npsRunAction* run)
:runAct(run)
//:runAct(run),eventMessenger(0)
{
  //  eventMessenger = new npsEventActionMessenger(this);
}


npsEventAction::~npsEventAction()
{
  //  delete eventMessenger;
}


void npsEventAction::BeginOfEventAction(const G4Event* evt)
{  

 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) { 
   //   G4cout << "\n---> Begin of event: " << evtNb << G4endl;
   //   CLHEP::HepRandom::showEngineStatus();
 }
 
 // initialisation per event
 EdepEv = 0.;
 TrackEv = 0.;
 Npe = 0;
}


void npsEventAction::EndOfEventAction(const G4Event* evt)
{
  //accumulates statistic
  //
  //  runAct->fillPerEvent(EdepEv, TrackEv);
  
  //print per event (modulo n)
  //
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    //    G4cout << "---> End of event: " << evtNb << G4endl;	

    //    G4cout
    //       << "   Total Energy Depostion: " << std::setw(7)
    //                                        << G4BestUnit(EdepEv,"Energy")
    //       << "       Total Track Length: " << std::setw(7)
    //                                        << G4BestUnit(TrackEv,"Length")
    //       << "       Npe total: " << Npe
    //       << G4endl;
  }

  std::ios::fmtflags curr_fmt = G4cout.flags();

  G4cout << std::showpoint << std::setprecision(5)
	 << EdepEv/GeV  << " " << TrackEv/cm << std::setw(5) << Npe << G4endl;
  G4cout.flags(curr_fmt);
}
