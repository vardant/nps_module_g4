class npsSteppingVerbose;

#ifndef npsSteppingVerbose_h
#define npsSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"



class npsSteppingVerbose : public G4SteppingVerbose
{
 public:   

   npsSteppingVerbose();
  ~npsSteppingVerbose();

   void StepInfo();
   void TrackingStarted();

};



#endif
