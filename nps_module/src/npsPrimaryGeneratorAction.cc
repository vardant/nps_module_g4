#include "npsPrimaryGeneratorAction.hh"
#include "npsPrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

npsPrimaryGeneratorAction::npsPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new npsPrimaryGeneratorMessenger(this);
  
  //default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
  //  G4ParticleDefinition* particle = particleTable->FindParticle("chargedgeantino");

  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleTime(0.0*ns);
  particleGun->SetParticleMomentum(50*MeV);

  particleGun->SetParticlePosition(G4ThreeVector(-10.*cm,0.,0.*cm));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

}


npsPrimaryGeneratorAction::~npsPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}



void npsPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //  G4double x = -5.*cm;
  //  G4double x = (G4UniformRand()-0.5)*9.*cm;
  //  G4double y = (G4UniformRand()-0.5)*9.*cm;
  //  G4double z = -25.*cm;
  //  particleGun->SetParticlePosition(G4ThreeVector(x,y,z));

  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,-10.1*cm));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  G4double Pin = particleGun->GetParticleMomentum();
  //  G4cout << "Pin = " << Pin/GeV << " GeV/c" << G4endl;

  G4ThreeVector XYZin =  particleGun->GetParticlePosition();
  //  G4cout << "XYZin (cm): " << XYZin/cm << G4endl;
  //  G4double Xin = XYZin[0];
  G4double Xin = XYZin[0];
  G4double Yin = XYZin[1];
  G4double Zin = XYZin[2];
  //  G4cout << "X, Y, Z (cm): " << Xin/cm << " " << Yin/cm << " " << Zin/cm 
  //	 << G4endl;

  G4int index = particleGun->GetParticleDefinition()->GetPDGEncoding();
  //  G4cout << "PDG encoding: " << index << G4endl;

  std::ios::fmtflags curr_fmt = G4cout.flags();

  G4cout << std::setw(4) << index << " "
	 << std::showpoint << std::setprecision(3)
	 << Pin/MeV << " " << Xin/cm << " " << Yin/cm << " " << Zin/cm << "  ";

  G4cout.flags(curr_fmt);

  particleGun->GeneratePrimaryVertex(anEvent);
}



void npsPrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

void npsPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 if (particleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton")
   {
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
               "the particleGun is not an opticalphoton" << G4endl;
     return;
   }
     	       
 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = particleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton); 
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product; 
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 particleGun->SetParticlePolarization(polar);
}
