#ifndef npsStackingAction_H
#define npsStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"



class npsStackingAction : public G4UserStackingAction
{
  public:
    npsStackingAction();
   ~npsStackingAction();

  public:
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    void NewStage();
    void PrepareNewEvent();

  private:
    G4int gammaCounter;
};



#endif

