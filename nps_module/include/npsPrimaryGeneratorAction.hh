#ifndef npsPrimaryGeneratorAction_h
#define npsPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class npsPrimaryGeneratorMessenger;



class npsPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    npsPrimaryGeneratorAction();
   ~npsPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);

    void SetOptPhotonPolar();
    void SetOptPhotonPolar(G4double);

  //     G4ParticleGun* GetParticleGun() {return particleGun;};

  private:
    G4ParticleGun* particleGun;
    npsPrimaryGeneratorMessenger* gunMessenger;

  //Primary particle.
  double fXmin, fXmax;   //cm
  double fYmin, fYmax;   //cm
  double fZmin, fZmax;   //cm
  double fPx, fPy, fPz;   //GeV
  std::string fParticle;
  
};



#endif /*npsPrimaryGeneratorAction_h*/
