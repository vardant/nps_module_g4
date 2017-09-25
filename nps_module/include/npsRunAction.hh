#ifndef npsRunAction_h
#define npsRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"



class G4Timer;
class G4Run;

class npsRunAction : public G4UserRunAction
{
  public:
    npsRunAction();
   ~npsRunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

  private:
    G4Timer* timer;
};



#endif /*npsRunAction_h*/
