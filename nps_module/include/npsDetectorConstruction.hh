#ifndef npsDetectorConstruction_h
#define npsDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

class npsDetectorConstruction : public G4VUserDetectorConstruction {

public:
  npsDetectorConstruction();
  ~npsDetectorConstruction();

public:
  G4VPhysicalVolume* Construct();

  //  const G4VPhysicalVolume* Get_block() {return block_phys;};

private:

  //Outer sizes of optical insulation.
  G4double tedlar_x;
  G4double tedlar_y;
  G4double tedlar_z;

  //Outer sizes of optical reflector.
  G4double mylar_x;
  G4double mylar_y;
  G4double mylar_z;

  //Thicknesses of components.
  G4double tedlar_thick;
  G4double mylar_thick;
  G4double glue_thick;    //Glue between PMT and scintillator block
  G4double air_gap;       //Air gap between reflector and scintillator block

  //Parameters of reflector to be read from input file (reflector.inp).
  G4int refFlag;
  G4String refName;
  G4int refNumData;
  G4double* refWL;
  G4double* refReIndex;
  G4double* refImIndex;
  G4double subRefrIndex;
  G4double* refRefl;

  //Scintillator block sizes.  
  G4double block_x;
  G4double block_y;
  G4double block_z;

  //PMT sizes.
  G4double PMT_diameter;
  G4double PMTWin_thick;    //Thickness of PMT glass window.
  G4double Cathode_diam;    //Photocathode diameter, thickness.
  G4double Cathode_thick;

  //Outer sizes of the whole assembly (counter).
  G4double counter_x;
  G4double counter_y;
  G4double counter_z;

  G4double expHall_x;
  G4double expHall_y;
  G4double expHall_z;

  G4LogicalVolume* block_log;
  G4LogicalVolume* mylar_log;
  G4LogicalVolume* PMTWin_right_log;
  G4LogicalVolume* PMTHouse_log;
  G4LogicalVolume* Cathode_log;
  G4LogicalVolume* glue_log;
  G4LogicalVolume* tedlar_log;
  G4LogicalVolume* counter_log;
  G4LogicalVolume* expHall_log;
};

#endif /*npsDetectorConstruction_h*/
