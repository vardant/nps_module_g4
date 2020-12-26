#ifndef npsEventAction_h
#define npsEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "npsDetectorConstruction.hh"

class G4Event;

class npsRunAction;
//class npsEventActionMessenger;

class npsEventAction : public G4UserEventAction
{
  public:
    npsEventAction(npsRunAction*);
   ~npsEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  //Accrue energy deposition and track length for the event.
    void AddEv(G4double de, G4double dl) {
      EdepEv += de; 
      TrackEv += dl;
    }

  //Accrue number of detected photoelectrons for the event.
  void AddNpe() {
      Npe += 1; 
    }

  private:
    npsRunAction*  runAct;
   
    G4double  EdepEv, TrackEv;   //Energy deposition and track lenth of event
    G4int Npe;                   //Number of detected photoelectrons of event
                     
  //   npsEventActionMessenger*  eventMessenger;

};


#endif

    
