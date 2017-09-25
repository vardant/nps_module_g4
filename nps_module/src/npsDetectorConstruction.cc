#include "npsDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

#include "G4NistManager.hh"

#include <fstream>

using namespace std;

// Single module from the Heremes calorimeter.

npsDetectorConstruction::npsDetectorConstruction() {

  ifstream fin;
  fin.open("reflector.inp");

  G4String line;
  getline(fin,line); istringstream iss1(line);
  iss1 >> air_gap;
  getline(fin,line); istringstream iss2(line);
  iss2 >> ReflectorFlag;

  fin.close();

  air_gap *= mm;
    
  G4cout << "npsDetectorConstruction::npsDetectorConstruction: input data:"
	  << G4endl;
  G4cout << "   Air gap = " << air_gap/mm << " mm" << G4endl;
  G4cout << "   ReflectorFlag = " << ReflectorFlag << ", ";
  switch (ReflectorFlag) {
  case 0:
    G4cout << "diffuse reflector";
    break;
  case 1:
    G4cout << "specular reflector with metal facing crystal";
    break;
  case 2:
    G4cout << "specular reflector with substrate facing crystal";
    break;
  default:
    G4cout << "not recognized, assume diffuse reflector";
    ReflectorFlag = 0;
  }
  G4cout << "." << G4endl;
  
  tedlar_thick = 0.040*mm;   //40um Tedlar
  mylar_thick = 0.025*mm;    // + 25um Mylar
  ////  air_gap = 0.035*mm;        //guess
  glue_thick = 0.035*mm;     //guess

  //Photon tracking test
  //  tedlar_thick = 1*mm;   //40um Tedlar
  //  mylar_thick = 1*mm;    // + 25um Mylar
  //  air_gap = 1*mm;        //guess
  //  glue_thick = 1*mm;     //test

  PMT_diameter = 1.86*cm;
  PMTWin_thick = 1*mm;     //??

  Cathode_diam = 1.5*cm;
  Cathode_thick = 0.1*mm;
  
  block_x = 2.05*cm;
  block_y = 2.05*cm;
  block_z = 20*cm;

  mylar_x = block_x + 2*air_gap + 2*mylar_thick;
  mylar_y = block_y + 2*air_gap + 2*mylar_thick;
  mylar_z = block_z + 2*air_gap + 2*mylar_thick;

  tedlar_x = mylar_x + 2*tedlar_thick;
  tedlar_y = mylar_y + 2*tedlar_thick;
  tedlar_z = mylar_z + 2*tedlar_thick;

  counter_x = tedlar_x;
  counter_y = tedlar_y;
  counter_z = tedlar_z + 2*glue_thick +  2*PMTWin_thick;

  expHall_x = counter_x * 1.5;
  expHall_y = counter_y * 1.5;
  expHall_z = counter_z * 1.5;
}

npsDetectorConstruction::~npsDetectorConstruction(){;}



G4VPhysicalVolume* npsDetectorConstruction::Construct()
{

//	------------- Materials -------------

  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);

  //  G4Material* Vac    = man->FindOrBuildMaterial("G4_Galactic");

  G4double density;
  G4int nelements;
  G4int ncomponents;
  //  G4double fractionmass;

  G4Element* H  = man->FindOrBuildElement("H");
  G4Element* Si = man->FindOrBuildElement("Si");
  G4Element* C  = man->FindOrBuildElement("C");
  G4Element* O  = man->FindOrBuildElement("O");
  G4Element* K  = man->FindOrBuildElement("K");
  G4Element* N  = man->FindOrBuildElement("N");
  G4Element* Cs = man->FindOrBuildElement("Cs");

  // Lead tungstate (PbWO4).

  G4Material* PbWO4 = man->FindOrBuildMaterial("G4_PbWO4");

//
// ------------ Generate & Add Material Properties Table ------------
//

  //Wevalengths in nm.
  G4double wlPbWO4[52] = {675.,
			  670.,660.,650.,640.,630.,620.,610.,600.,590.,580.,
			  570.,560.,550.,540.,530.,520.,510.,500.,490.,480.,
			  470.,460.,450.,440.,430.,420.,410.,400.,390.,380.,
			  370.,360.,350.,340.,330.,320.,310.,300.,290.,280.,
			  270.,260.,250.,240.,230.,220.,210.,200.,190.,180.,
			  175.};

  for (G4int i=0; i<52; i++) wlPbWO4[i] *= nanometer;

  const G4double hc = 1.239841857E-6*m*eV;   //(PDG)

  G4double kphotPbWO4[52];   //Momenta of optical photons in eV units.
  for (G4int i=0; i<52; i++) kphotPbWO4[i] = hc/wlPbWO4[i];

  G4double abslength[52] = {
    1400.,
    1400.,1400.,1400.,1400.,1400.,1400.,1400.,933.3,933.3,933.3,
    933.3,933.3,933.3,933.3,933.3,933.3,700.0,700.0,622.2,560.0,
    560.0,466.6,350.0,280.0,233.3,175.0,151.3,112.0,71.79,45.52,
    29.62,17.07,10.17,6.026,3.557,2.092,1.227,0.717,0.418,0.243,
    0.140,0.081,0.047,0.027,0.016,0.009,0.005,0.003,0.002,0.001,
    0.000711281};

  for (G4int i=0; i<52; i++) {
    abslength[i] *= cm;
  };

  G4double rindPbWO4[52];
  for (G4int i=0; i<52; i++) {
    rindPbWO4[i] = 2.2;             //PbWO conventional refractive index
  };

  G4double wlPbWO4_sc_fast[82] = {
    630.,
    626.,622.,618.,614.,610.,606.,602.,598.,594.,590.,
    586.,582.,578.,574.,570.,566.,562.,558.,554.,550.,
    546.,542.,538.,534.,530.,526.,522.,518.,514.,510.,
    506.,502.,498.,494.,490.,486.,482.,478.,474.,470.,
    466.,462.,458.,454.,450.,446.,442.,438.,434.,430.,
    426.,422.,418.,414.,410.,406.,402.,398.,394.,390.,
    386.,382.,378.,374.,370.,366.,362.,358.,354.,350.,
    346.,342.,338.,334.,330.,326.,322.,318.,314.,310.,
    306.};

  G4double wlPbWO4_sc_slow[82];
  for (G4int i=0; i<82; i++) wlPbWO4_sc_slow[i] =  wlPbWO4_sc_fast[i] + 5.;
  
  for (G4int i=0; i<82; i++) wlPbWO4_sc_fast[i] *= nanometer;
  for (G4int i=0; i<82; i++) wlPbWO4_sc_slow[i] *= nanometer;

  G4double kphotPbWO4_sc_fast[82];
  G4double kphotPbWO4_sc_slow[82];
  for (G4int i=0; i<82; i++) kphotPbWO4_sc_fast[i] = hc/wlPbWO4_sc_fast[i];
  for (G4int i=0; i<82; i++) kphotPbWO4_sc_slow[i] = hc/wlPbWO4_sc_slow[i];

  G4double PbWO4_sc_fast[82] = {
    0.,
    0.019,0.045,0.064,0.058,0.058,0.064,0.070,0.064,0.064,0.064,
    0.070,0.070,0.090,0.077,0.096,0.122,0.109,0.141,0.134,0.154,
    0.186,0.166,0.192,0.205,0.218,0.243,0.256,0.269,0.288,0.320,
    0.358,0.390,0.416,0.429,0.467,0.512,0.544,0.589,0.627,0.640,
    0.704,0.730,0.774,0.794,0.838,0.870,0.909,0.928,0.934,0.986,
    0.979,0.998,0.992,0.986,0.973,0.941,0.902,0.870,0.819,0.787,
    0.730,0.691,0.653,0.589,0.538,0.461,0.410,0.326,0.282,0.224,
    0.173,0.102,0.070,0.051,0.013,0.000,0.000,0.000,0.000,0.000,
    0.000};

  G4double PbWO4_sc_slow[82];
  for (G4int i=0; i<82; i++) PbWO4_sc_slow[i] = PbWO4_sc_fast[i];
    
  G4MaterialPropertiesTable *PbWO4MPT = new G4MaterialPropertiesTable();
  
  PbWO4MPT -> AddProperty("RINDEX",kphotPbWO4,rindPbWO4,52);
  PbWO4MPT -> AddProperty("ABSLENGTH",kphotPbWO4,abslength,52);

  PbWO4MPT->AddProperty("FASTCOMPONENT",kphotPbWO4_sc_fast,PbWO4_sc_fast,82);
  PbWO4MPT->AddProperty("SLOWCOMPONENT",kphotPbWO4_sc_slow,PbWO4_sc_slow,82);
  PbWO4MPT->AddConstProperty("SCINTILLATIONYIELD", 40000*0.377/100/MeV);
  PbWO4MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  PbWO4MPT->AddConstProperty("FASTTIMECONSTANT", 10.*ns);
  PbWO4MPT->AddConstProperty("SLOWTIMECONSTANT", 30.*ns);
  PbWO4MPT->AddConstProperty("YIELDRATIO", 0.077/(0.077+0.3));

  PbWO4 -> SetMaterialPropertiesTable(PbWO4MPT);

  //  G4cout << "PbWO4 optical properties:" << G4endl;
  //  for (G4int i=0; i<52; i++) {
  //    G4cout << G4BestUnit(wlPbWO4[i],"Length")
  //	   << G4BestUnit(kphotPbWO4[i],"Energy")
  //	   << G4BestUnit(abslength[i],"Length")
  //	   << rind[i] << G4endl;
  //  }

  // Air
  // 
  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  G4double rindAir[52];
  for (G4int i=0; i<52; i++) {
    rindAir[i] = 1.000293;   //Air @ STP
  };
  G4MaterialPropertiesTable *AirMPT = new G4MaterialPropertiesTable();
  AirMPT -> AddProperty("RINDEX",kphotPbWO4,rindAir,52);
  Air -> SetMaterialPropertiesTable(AirMPT);

  // Glass
  //

  density = 2.23*g/cm3;   //Borosilicate glass (wikipedia)
  G4Material* Glass = new G4Material("Glass", density, ncomponents=2);
  Glass->AddElement(Si, 1);
  Glass->AddElement(O,  2);

  G4double rindGlass[52];
  for (G4int i=0; i<52; i++) {
    rindGlass[i] = 1.525;              //average of 1.51-1.54
  };

  G4MaterialPropertiesTable *GlassMPT = new G4MaterialPropertiesTable();
  GlassMPT -> AddProperty("RINDEX",kphotPbWO4,rindGlass,52);
  Glass -> SetMaterialPropertiesTable(GlassMPT);

  // Optical grease BC630 from Bicron
  //
  density = 1.06*g/cm3;
  G4Material* OpticalGlue = new G4Material("Silgard", density, ncomponents=1);
  OpticalGlue->AddElement(Si, 1); //not known

  G4double rindGlue[52];
  for (G4int i=0; i<52; i++) {
    rindGlue[i] = 1.465;
  };

  G4MaterialPropertiesTable *GlueMPT = new G4MaterialPropertiesTable();
  GlueMPT -> AddProperty("RINDEX",kphotPbWO4,rindGlue,52);
  OpticalGlue -> SetMaterialPropertiesTable(GlueMPT);

  // Optical insulation
  //
  density = 1.5;   //approximately
  G4Material* Polymer = new G4Material("Polymer", density, ncomponents=2);
  Polymer->AddElement(C, 1);
  Polymer->AddElement(H, 1);

  G4Material* Mylar = man->FindOrBuildMaterial("G4_MYLAR");

  if (ReflectorFlag == 2) {
    
    //Mylar refractive index.
    G4double rindMylar[52];
    for (G4int i=0; i<52; i++) {
      rindMylar[i] = 1.65;
    };

    G4MaterialPropertiesTable *MylarMPT = new G4MaterialPropertiesTable();
    MylarMPT -> AddProperty("RINDEX",kphotPbWO4,rindMylar,52);
    Mylar -> SetMaterialPropertiesTable(MylarMPT);

  }

  // Bialcali, the photochathode material
  //

  density = 1*g/cm3;   //Does not matter
  G4Material* Bialcali = new G4Material("Bialcali", density, ncomponents=2);
  Bialcali->AddElement(Cs, 1);
  Bialcali->AddElement(K,  1);

//
//	------------- Volumes --------------

  // PbWO4 block
  //
  G4Box* block_box = new G4Box("Block_box",block_x/2,block_y/2,block_z/2);

  block_log = new G4LogicalVolume(block_box,PbWO4,"Block_log",0,0,0);

  // Optical insulation
  //
  G4Box* tedlar_outer =
    new G4Box("Tedlar_solid",tedlar_x/2,tedlar_y/2,tedlar_z/2);

  G4Box* tedlar_inner = new G4Box("Tedlar_cavity",
				 tedlar_x/2-tedlar_thick,
				 tedlar_y/2-tedlar_thick,
				 tedlar_z/2-tedlar_thick);

  G4SubtractionSolid* tedlar_box = new G4SubtractionSolid("Tedlar",
				  tedlar_outer, tedlar_inner);

  // Make a hole of PMT size
  G4Tubs*  tedlar_hole = new G4Tubs("tedlar_hole",
				    0., PMT_diameter/2, tedlar_thick/2,
  				    0.*deg, 360.*deg);

  // Optical insulation with hole on right side.
  G4RotationMatrix rot;
  G4ThreeVector z_trans_tedlar_hole(0, 0, tedlar_z/2 - tedlar_thick/2);
  G4Transform3D trans_tedlar_hole(rot, z_trans_tedlar_hole);
  G4SubtractionSolid* tedlar_holed = new G4SubtractionSolid("Tedlar_holed",
			      tedlar_box, tedlar_hole, trans_tedlar_hole);

  //Remove front wall of Tedlar
  G4Box* tedlar_front = new G4Box("Tedlar_fr",
				  tedlar_x/2,tedlar_y/2,tedlar_thick/2);
  G4ThreeVector z_trans_tedlar_front(0, 0, -tedlar_z/2 + tedlar_thick/2);
  G4Transform3D trans_tedlar_front(rot, z_trans_tedlar_front);
  G4SubtractionSolid* tedlar_frame = new G4SubtractionSolid("Tedlar",
			      tedlar_holed, tedlar_front, trans_tedlar_front);

  tedlar_log = new G4LogicalVolume(tedlar_frame,Polymer,"Tedlar",0,0,0);

  //Mylar, reflector.
  
  G4Box* mylar_outer = new G4Box("Mylar_solid",mylar_x/2,mylar_y/2,mylar_z/2);

  G4Box* mylar_inner = new G4Box("Mylar_cavity",
				 mylar_x/2-mylar_thick,
				 mylar_y/2-mylar_thick,
				 mylar_z/2-mylar_thick);

  G4SubtractionSolid* mylar_box = new G4SubtractionSolid("Mylar",
				  mylar_outer, mylar_inner);

  G4Tubs*  mylar_hole = new G4Tubs("mylar_hole", 0., PMT_diameter/2,
				   mylar_thick/2, 0.*deg, 360.*deg);

  G4ThreeVector z_trans_mylar_hole(0, 0, mylar_z/2 - mylar_thick/2);
  G4Transform3D trans_mylar_hole(rot, z_trans_mylar_hole);
  G4SubtractionSolid* mylar_holed = new G4SubtractionSolid("Mylar",
				      mylar_box, mylar_hole, trans_mylar_hole);

  //Remove front wall of Mylar
  G4Box* mylar_front = new G4Box("Mylar_fr",mylar_x/2,mylar_y/2,mylar_thick/2);
  G4ThreeVector z_trans_mylar_front(0, 0, -mylar_z/2 + mylar_thick/2);
  G4Transform3D trans_mylar_front(rot, z_trans_mylar_front);
  G4SubtractionSolid* mylar_frame = new G4SubtractionSolid("Mylar_holed",
			      mylar_holed, mylar_front, trans_mylar_front);

  mylar_log=new G4LogicalVolume(mylar_frame,Mylar,"Mylar",0,0,0);

  // PMT Window
  //
  G4Tubs*  PMTWin_tube =
  new G4Tubs("PMTWindow", 0., PMT_diameter/2, PMTWin_thick/2,0.*deg, 360.*deg);

  PMTWin_right_log = new G4LogicalVolume(PMTWin_tube,Glass, "PMTWindow");

  // PMT Housing
  //
  //  G4Tubs*  PMTHouse_tube = new G4Tubs("PMTHouse", PMT_diameter/2,
  //  G4Tubs*  PMTHouse_tube = new G4Tubs("PMTHouse", 0.,
  //		 PMT_diameter/2+1.*mm, PMTWin_thick/2, 0.*deg, 360.*deg);

  //  PMTHouse_log = new G4LogicalVolume(PMTHouse_tube, Polymer, "PMTHousing");

  // Photocathode
  //
  G4Tubs*  Cathode_tube =
  new G4Tubs("Cathode", 0., Cathode_diam/2, Cathode_thick/2,0.*deg, 360.*deg);

  Cathode_log = new G4LogicalVolume(Cathode_tube, Bialcali, "Cathode");

  // Optical glue
  //
  G4Tubs*  glue_tube =
    new G4Tubs("glue", 0., PMT_diameter/2, glue_thick/2, 0.*deg, 360.*deg);

  glue_log = new G4LogicalVolume(glue_tube,OpticalGlue, "Glue");

  // Counter
  //
  G4Box* counter_box = new G4Box("Counter",counter_x/2,counter_y/2,counter_z/2);

  counter_log = new G4LogicalVolume(counter_box,Air,"Counter",0,0,0);

  // The experimental Hall
  //
  G4Box* expHall_box = new G4Box("World",expHall_x/2,expHall_y/2,expHall_z/2);

  expHall_log = new G4LogicalVolume(expHall_box,Air,"World",0,0,0);


  // Place constituents, construct physical volumes.
  //

  G4VPhysicalVolume* expHall_phys =
    new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

  //  G4VPhysicalVolume* counter_phys =
  new G4PVPlacement(0, //no rotation
		    G4ThreeVector(),
		    counter_log, //its logical volume
		    "Counter",   //its name
		    expHall_log,     //its mother  volume
		    false,         //no boolean operation
		    0);  //copy number

  // Insulation for counter
  new G4PVPlacement(0,  //no rotation
  		    G4ThreeVector(),
  		    tedlar_log,      //its logical volume
  		    "Tedlar",          //its name
  		    counter_log,      //its mother  volume
  		    false,              //no boolean operation
  		    0);                 //copy number
  
  //  G4VPhysicalVolume* mylar_phys =
  new G4PVPlacement(0,  //no rotation
		    G4ThreeVector(),
		    mylar_log,        //its logical volume
		    "Mylar_phys",       //its name
		    counter_log,    //its mother  volume
		    false,            //no boolean operation
		    0);               //copy number
  
  // Glass block for counter

  new G4PVPlacement(0,  //no rotation
		    G4ThreeVector(),
		    block_log,        //its logical volume
		    "Block_phys",     //its name
		    counter_log,        //its mother  volume
		    false,            //no boolean operation
		    0);               //copy number

  // Glue, window and cathode for counter.

  G4double x = 0.;
  G4double y = 0.;
  G4double z = block_z/2 + glue_thick/2;
  new G4PVPlacement(0, //no rotation
		    G4ThreeVector(x,y,z),
		    glue_log,    //its logical volume
		    "Glue",      //its name
		    counter_log, //its mother  volume
		    false,         //no boolean operation
		    0);            //copy number

  z = block_z/2 + glue_thick + PMTWin_thick/2;

  ////  new G4PVPlacement(0,    //no rotation
  ////		    G4ThreeVector(x,y,z),
  ////		    PMTHouse_log, //its logical volume
  ////		    "PMTHousing",      //its name
  ////		    counter_log,    //its mother  volume
  ////		    false,            //no boolean oper.
  ////		    0);               //copy number

  //  G4VPhysicalVolume* PMTWin_phys =
  new G4PVPlacement(0,    //no rotation
		    G4ThreeVector(x,y,z),
		    PMTWin_right_log, //its logical volume
		    "PMTWindow",      //its name
		    counter_log,    //its mother  volume
		    false,            //no boolean oper.
		    0);               //copy number

  ////  z = PMTWin_thick/2 - Cathode_thick/2;
  ////  new G4PVPlacement(0, //no rotation
  ////		    G4ThreeVector(x,y,z),
  ////		    Cathode_log,  //its logical volume
  ////		    "Cathode", //its name
  ////		    PMTWin_right_log, //its mother  volume
  ////		    false,       //no boolean operation
  ////		    0);          //copy number

  z = block_z/2 + glue_thick + PMTWin_thick + Cathode_thick/2;
  new G4PVPlacement(0, //no rotation
		    G4ThreeVector(x,y,z),
		    Cathode_log,  //its logical volume
		    "Cathode", //its name
		    counter_log, //its mother  volume
		    false,       //no boolean operation
		    0);          //copy number


//	------------- Surfaces --------------
//

  G4MaterialPropertiesTable* ReflectorMPT = new G4MaterialPropertiesTable();
  G4OpticalSurface* Reflector = new G4OpticalSurface("Reflector");

  if (ReflectorFlag > 0) {

    // Specular reflector, aluminum.
    
    G4double wlAl[35] = {695.32,666.93,639.71,613.59,588.55,564.52,541.47,
			 519.37,498.17,477.83,458.33,439.62,421.67,404.46,
			 387.95,372.11,356.92,342.35,328.37,314.97,302.11,
			 289.78,277.95,266.60,255.72,245.28,235.27,225.66,
			 216.45,207.62,199.14,191.01,183.21,175.73,168.56};

    for (G4int i=0; i<35; i++) wlAl[i] *= nanometer;

    G4double kphotAl[35];   //Momenta of optical photons in eV units.
    for (G4int i=0; i<35; i++) kphotAl[i] = hc/wlAl[i];

    G4double real_index[35] = {
      1.62030,1.42860,1.29520,1.19110,1.100500,1.016200,0.93525,
      0.85697,0.78186,0.71074,0.64439,0.583360,0.527880,0.4779,
      0.43317,0.39330,0.35785,0.32634,0.298320,0.273360,0.25107,
      0.23109,0.21312,0.19687,0.18210,0.168600,0.156170,0.14468,
      0.13398,0.12398,0.11461,0.10582,0.097571,0.089866,0.082703
    };

    G4double imag_index[35] = {
      7.9931,7.6672,7.3576,7.0711,6.8042,6.5519,6.3099,
      6.0750,5.8455,5.6203,5.3994,5.1830,4.9717,4.7659,
      4.5661,4.3726,4.1856,4.0052,3.8313,3.6640,3.503,
      3.3482,3.1992,3.0560,2.9182,2.7856,2.6578,2.5347,
      2.4159,2.3011,2.1901,2.0826,1.9783,1.8771,1.7786
    };
    /*
      G4double wlAg[102] = {
      775.00,
      751.52,729.41,708.57,688.89,670.27,652.63,635.90,620.00,604.88,590.48,
      576.74,563.64,551.11,539.13,527.66,516.67,506.12,496.00,486.27,476.92,
      467.92,459.26,450.91,442.86,435.09,427.59,420.34,413.33,406.56,400.00,
      393.65,387.50,n381.54,375.76,370.15,364.71,359.42,354.29,349.30,344.44,
      339.73,335.14,330.67,326.32,322.08,317.95,313.92,310.00,306.17,302.44,
      298.80,295.24,291.76,288.37,285.06,281.82,278.65,275.56,272.53,269.57,
      266.67,263.83,261.05,258.33,255.67,253.06,250.51,248.00,245.54,243.14,
      240.78,238.46,236.19,233.96,231.78,229.63,227.52,225.45,223.42,221.43,
      219.47,217.54,215.65,213.79,211.97,210.17,208.40,206.67,204.96,203.28,
      201.63,200.00,198.40,196.83,195.28,193.75,192.25,190.77,189.31,187.88,
      175.00
      };

      for (G4int i=0; i<102; i++) wlAg[i] *= nanometer;

      G4double kphotAg[102];   //Momenta of optical photons in eV units.
      for (G4int i=0; i<102; i++) kphotAg[i] = hc/wlAg[i];

      G4double real_index[102] = {
      0.143  ,
      0.1459 ,0.148 ,0.1443 ,0.14 ,0.1401 ,0.14 ,0.1361 ,0.131 ,0.1255 ,0.121 ,
      0.1193 ,0.12 ,0.1244 ,0.129 ,0.1301 ,0.13 ,0.1299 ,0.13 ,0.1303 ,0.132 ,
      0.1373 ,0.144 ,0.1511 ,0.157 ,0.1585 ,0.16 ,0.1667 ,0.173 ,0.1726 ,0.173 ,
      0.182 ,0.192 ,0.1981 ,0.2 ,0.1921 ,0.186 ,0.1948 ,0.209 ,0.2214 ,0.238 ,
      0.2533 ,0.294 ,0.3881 ,0.526 ,0.7191 ,0.932 ,1.1421 ,1.323 ,1.4325 ,1.496 ,
      1.5194 ,1.519 ,1.5136 ,1.502 ,1.4901 ,1.476 ,1.4592 ,1.441 ,1.4223 ,1.404 ,
      1.3875 ,1.372 ,1.3569 ,1.343 ,1.3311 ,1.32 ,1.3086 ,1.298 ,1.2897 ,1.282 ,
      1.2734 ,1.265 ,1.2572 ,1.25 ,1.244 ,1.238 ,1.2307 ,1.223 ,1.2157 ,1.208 ,
      1.1991 ,1.19 ,1.1819 ,1.173 ,1.1614 ,1.149 ,1.1372 ,1.125 ,1.1116 ,1.098 ,
      1.0848 ,1.072 ,1.0596 ,1.048 ,1.0375 ,1.028 ,1.0195 ,1.012 ,1.0043 ,0.995 ,
      0.995
      };

      G4double imag_index[102] = {
      5.09 ,
      4.9081 ,4.74 ,4.5863 ,4.44 ,4.2931 ,4.15 ,4.0106 ,3.88 ,3.7663 ,3.66 ,
      3.5538 ,3.45 ,3.3481 ,3.25 ,3.1594 ,3.07 ,2.9738 ,2.88 ,2.7981 ,2.72 ,
      2.64 ,2.56 ,2.4788 ,2.4 ,2.3294 ,2.26 ,2.1863 ,2.11 ,2.0294 ,1.95 ,
      1.8788 ,1.81 ,1.735 ,1.67 ,1.6419 ,1.61 ,1.5338 ,1.44 ,1.3453 ,1.24 ,
      1.1207 ,0.986 ,0.8186 ,0.663 ,0.5544 ,0.504 ,0.5509 ,0.647 ,0.7611 ,0.882 ,
      0.9888 ,1.08 ,1.143 ,1.19 ,1.2288 ,1.26 ,1.2881 ,1.31 ,1.3219 ,1.33 ,
      1.3413 ,1.35 ,1.3513 ,1.35 ,1.35 ,1.35 ,1.3506 ,1.35 ,1.3456 ,1.34 ,
      1.335 ,1.33 ,1.325 ,1.32 ,1.3147 ,1.31 ,1.3072 ,1.305 ,1.3025 ,1.3 ,
      1.2975 ,1.295 ,1.2928 ,1.29 ,1.2853 ,1.28 ,1.275 ,1.27 ,1.2656 ,1.26 ,
      1.2513 ,1.24 ,1.2256 ,1.21 ,1.1944 ,1.18 ,1.1681 ,1.16 ,1.1494 ,1.13 ,
      1.13
      };
    */
    ReflectorMPT -> AddProperty("REALRINDEX",kphotAl,real_index,35);
    ReflectorMPT -> AddProperty("IMAGINARYRINDEX",kphotAl,imag_index,35);
    ////  ReflectorMPT -> AddProperty("REALRINDEX",kphotAg,real_index,102);
    ////  ReflectorMPT -> AddProperty("IMAGINARYRINDEX",kphotAg,imag_index,102);

    Reflector -> SetType(dielectric_metal);
    Reflector -> SetFinish(polished);
    Reflector -> SetModel(glisur);
  }
  else {
    // Diffuse reflector, PTFE (Teflon).

    G4double wlPTFE[14] = {800.,750.,700.,650.,600.,550.,500.,
			   450.,400.,350.,300.,250.,200.,150.};

    for (G4int i=0; i<14; i++) wlPTFE[i] *= nanometer;

    G4double kphotPTFE[14];
    for (G4int i=0; i<14; i++) kphotPTFE[i] = hc/wlPTFE[i];

    G4double refPTFE[14]= {0.755079,0.793454,0.829571,0.861174,0.889391,
			   0.918736,0.934537,0.941309,0.952596,0.950339,
			   0.936795,0.905192,0.905192,0.905192};

    ReflectorMPT -> AddProperty("REFLECTIVITY",kphotPTFE,refPTFE,14);

    Reflector -> SetModel(unified);
    Reflector -> SetFinish(groundfrontpainted);   //Purely Lambertian reflection
  }

  Reflector -> SetMaterialPropertiesTable(ReflectorMPT);

  G4cout << "===== ReflectorMPT: ============================" << G4endl;
  ReflectorMPT->DumpTable();
  Reflector->DumpInfo();

  //  G4cout << "Dumping optical properties of Reflector ..." << G4endl;
  //  ((G4OpticalSurface*)
  //   (Reflector->GetSurface(tedlar_log)->GetSurfaceProperty()))->DumpInfo();

  if (ReflectorFlag < 2)
    // Reflective front surface of Mylar.
    // Remove refraction index of Mylar (property table) in this case
    // (to make it opaque).
    new G4LogicalSkinSurface("Reflector",mylar_log,Reflector);
  else
    // Reflective back surface of Mylar.
    // Tedlar borders Mylar from back. Making it reflective, makes effectively
    // Mylar back surface reflective.
    // Activate property table of Mylar with its refractive index in this case.
    new G4LogicalSkinSurface("Reflector",tedlar_log,Reflector);

  // Cathode efficiency for Phylips XP3461 PMT.
  //

  G4double wlCat[101] = {675.,670.,665.,660.,655.,650.,645.,640.,635.,630.,
			 625.,620.,615.,610.,605.,600.,595.,590.,585.,580.,
			 575.,570.,565.,560.,555.,550.,545.,540.,535.,530.,
			 525.,520.,515.,510.,505.,500.,495.,490.,485.,480.,
			 475.,470.,465.,460.,455.,450.,445.,440.,435.,430.,
			 425.,420.,415.,410.,405.,400.,395.,390.,385.,380.,
			 375.,370.,365.,360.,355.,350.,345.,340.,335.,330.,
			 325.,320.,315.,310.,305.,300.,295.,290.,285.,280.,
			 275.,270.,265.,260.,255.,250.,245.,240.,235.,230.,
			 225.,220.,215.,210.,205.,200.,195.,190.,185.,180.,
			 175.};

  for (G4int i=0; i<101; i++) {
    wlCat[i] *= nanometer;
  };

  G4double kphotCat[101];   //Momenta of optical photons in eV units.
  for (G4int i=0; i<101; i++) kphotCat[i] = hc/wlCat[i];

  // Hamamatsu R4125 quantum efficiency (bialcali photocathode, borosilicate
  // window). Taken from the Hamamatsu booklet, p.65.
  G4double effCat[101] = {
    0.0030,0.0035,0.0040,0.0046,0.0052,0.0060,0.0068,0.0077,0.0087,0.0099,
    0.0112,0.0126,0.0141,0.0159,0.0177,0.0198,0.0221,0.0245,0.0272,0.0301,
    0.0332,0.0365,0.0401,0.0440,0.0481,0.0525,0.0572,0.0621,0.0673,0.0728,
    0.0785,0.0846,0.0908,0.0973,0.1041,0.1110,0.1181,0.1255,0.1329,0.1405,
    0.1482,0.1560,0.1638,0.1716,0.1793,0.1870,0.1946,0.2020,0.2092,0.2162,
    0.2229,0.2293,0.2354,0.2411,0.2463,0.2511,0.2554,0.2592,0.2625,0.2651,
    0.2673,0.2688,0.2697,0.2700,0.2688,0.2653,0.2595,0.2517,0.2419,0.2305,
    0.2177,0.2038,0.1891,0.1740,0.1586,0.1434,0.1285,0.1141,0.1004,0.0877,
    0.0758,0.0650,0.0553,0.0466,0.0389,0.0322,0.0264,0.0215,0.0173,0.0138,
    0.0110,0.0086,0.0067,0.0052,0.0040,0.0030,0.0023,0.0017,0.0012,0.0009,
    0.0007};

  G4double reflCat[101];
  for (G4int i = 0; i < 101; i++) {
    reflCat[i] = 0.;
  }

  G4OpticalSurface* surfCat = new G4OpticalSurface("Cathode");

  surfCat -> SetType(dielectric_metal);
  surfCat -> SetFinish(polished);
  surfCat -> SetModel(glisur);

  G4MaterialPropertiesTable* surfCatMPT = new G4MaterialPropertiesTable();
  surfCatMPT -> AddProperty("REFLECTIVITY",kphotCat,reflCat,101);
  surfCatMPT -> AddProperty("EFFICIENCY",kphotCat,effCat,101);

  surfCat -> SetMaterialPropertiesTable(surfCatMPT);

  new G4LogicalSkinSurface("Cathode",Cathode_log,surfCat);

  //test. PMT surface, black.
  //
  //  G4double reflPMT[101];
  //  for (G4int i = 0; i < 101; i++) {
  //    reflPMT[i] = 0.;
  //  }
  //
  //  G4OpticalSurface* surfPMT = new G4OpticalSurface("PMTSurface");
  //
  //  surfPMT -> SetType(dielectric_dielectric);
  //  surfPMT -> SetFinish(polished);
  //  surfPMT -> SetModel(glisur);
  //
  //  G4MaterialPropertiesTable* mptPMT = new G4MaterialPropertiesTable();
  //  mptPMT -> AddProperty("REFLECTIVITY",kphotCat,reflPMT,101);
  //
  //  surfPMT -> SetMaterialPropertiesTable(mptPMT);
  //
  //  new G4LogicalBorderSurface("PMTSurface",PMTWin_phys,counter_phys,surfPMT);

// Visualisation attributes
//
  expHall_log->SetVisAttributes (G4VisAttributes::Invisible);
  counter_log->SetVisAttributes (G4VisAttributes::Invisible);
  //  counter_end_log->SetVisAttributes (G4VisAttributes::Invisible);
  //  PMTWin_log->SetVisAttributes (G4VisAttributes::Invisible);
  //  glue_log->SetVisAttributes (G4VisAttributes::Invisible);

  // print the table of materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//
//always return the physical World
//
  return expHall_phys;
}
