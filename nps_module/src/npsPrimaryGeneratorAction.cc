#include "npsPrimaryGeneratorAction.hh"
#include "npsPrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

npsPrimaryGeneratorAction::npsPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new npsPrimaryGeneratorMessenger(this);

  //Input primary particle's coordinates, momenta, name.
  
  ifstream fin;
  fin.open("beam.inp");

  istringstream iss;
  G4String line;
  
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> fXmin >> fXmax;
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> fYmin >> fYmax;
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> fZmin >> fZmax;
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> fPx >> fPy >> fPz;
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> fParticle;

  fXmin *= cm; fXmax *= cm;
  fYmin *= cm; fYmax *= cm;
  fZmin *= cm; fZmax *= cm;
  fPx *= GeV; fPy *= GeV; fPz *= GeV;

  cout << "PrimaryGeneratorAction: --------------------" << endl;
  cout << "Particle: " << fParticle << endl;
  cout << "Px, Py, Pz: " << fPx/GeV << ", " << fPy/GeV << ", " << fPz/GeV
       << " GeV/c" <<endl;
  cout << "Xmin, Xmax: " << fXmin/cm << ", " << fXmax/cm << " cm" << endl;
  cout << "Ymin, Ymax: " << fYmin/cm << ", " << fYmax/cm << " cm" << endl;
  cout << "Zmin, Zmax: " << fZmin/cm << ", " << fZmax/cm << " cm" << endl;
  cout << "-----------------------------------------------" << endl;
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle=particleTable->FindParticle(fParticle.c_str());
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleTime(0.0*ns);
  double P = sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
  particleGun->SetParticleMomentum(P);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPx/P,fPy/P,fPz/P));

}


npsPrimaryGeneratorAction::~npsPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}


void npsPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //Sample input coordinates.
  G4double x = fXmin + (fXmax-fXmin)*G4UniformRand();
  G4double y = fYmin + (fYmax-fYmin)*G4UniformRand();
  G4double z = fZmin + (fZmax-fZmin)*G4UniformRand();
  particleGun->SetParticlePosition(G4ThreeVector(x,y,z));

  //Get coordinates and momnentum of the primary particle and output them.

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

  //Generate primary vertex.

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
