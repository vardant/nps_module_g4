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

// Single module of the NPS calorimeter.

npsDetectorConstruction::npsDetectorConstruction() {

  //Read in reflector parameters.

  ifstream fin;
  fin.open("reflector.inp");

  G4String line;
  getline(fin,line); istringstream iss1(line);
  iss1 >> air_gap;
  getline(fin,line); istringstream iss2(line);
  iss2 >> refFlag;
  getline(fin,line); istringstream iss3(line);
  iss3 >> refName;
  getline(fin,line); istringstream iss4(line);
  iss4 >> refNumData;
  refWL = new G4double[refNumData];

  //For the specular reflector, read real and imaginary parts of refractive
  //index. For the diffuse reflector, read reflectivity.
  /*
  if (refFlag!=0) {
    refReIndex = new G4double[refNumData];
    refImIndex = new G4double[refNumData];
    for (G4int i=refNumData-1; i>-1; i--)
      fin >> refWL[i] >> refReIndex[i] >> refImIndex[i];
  }
  else {
    refRefl = new G4double[refNumData];
    for (G4int i=refNumData-1; i>-1; i--)
      fin >> refWL[i] >> refRefl[i];
  }
  */
  if (refFlag==1) {
    refReIndex = new G4double[refNumData];
    refImIndex = new G4double[refNumData];
    for (G4int i=refNumData-1; i>-1; i--)
      fin >> refWL[i] >> refReIndex[i] >> refImIndex[i];
  }
  else {
    refRefl = new G4double[refNumData];
    for (G4int i=refNumData-1; i>-1; i--)
      fin >> refWL[i] >> refRefl[i];
  }
  
  //Read refractive index of substrate of reflector.
  fin >> subRefrIndex;
  fin.close();

  air_gap *= mm;

  for (G4int i=0; i<refNumData; i++) refWL[i] *= nanometer;

  //Print out parameters of reflector.
    
  G4cout << "npsDetectorConstruction::npsDetectorConstruction: input data:"
	  << G4endl;
  G4cout << "   Air gap = " << air_gap/mm << " mm" << G4endl;
  G4cout << "   Reflector: " << refName << ", refFlag = " << refFlag << ", ";
  if (refFlag==0)
    G4cout << "diffuse reflector";
  else
    G4cout << "specular reflector";
  G4cout << "." << G4endl;

  G4cout << "   Reflector data:" << G4endl;
  for (G4int i=refNumData-1; i>-1; i--) {
    G4cout << "   " << refWL[i]/nanometer << " ";
    if (refFlag==1)
      G4cout  << refReIndex[i] << " " << refImIndex[i];
    else
      G4cout << refRefl[i];
    G4cout << " " << i << G4endl;
  };

  G4cout << "   Substrate refr. index = " << subRefrIndex;
  if (subRefrIndex == 0.)
    G4cout << ", no substrate between crystal and reflector.";
  else
    G4cout << ", substrate layer between crystal and reflector.";
  G4cout << G4endl;

  //
    
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

//Construct the module.

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
  G4Element* F  = man->FindOrBuildElement("F");

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

  //PbWO4 measured absorption length.
  /*
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
  */

  //Absorption length from transmittance measurements of 19 SICCAS crystals
  //at CUA. Transmittances from Vladimir.
  /*
  const int nAbsl = 60;
  double wlAbsl[nAbsl] = {800.,790.,780.,770.,760.,750.,740.,730.,720.,710.,
			  700.,690.,680.,670.,660.,650.,640.,630.,620.,610.,
			  600.,590.,580.,570.,560.,550.,540.,530.,520.,510.,
			  500.,490.,480.,470.,460.,450.,440.,430.,420.,410.,
			  400.,390.,380.,370.,360.,350.,340.,330.,320.,310.,
			  300.,290.,280.,270.,260.,250.,240.,230.,220.,210.};
  double abslength[nAbsl] = {440.983, 441.162, 442.037, 447.54, 451.036, 
			     449.096, 440.888, 437.576, 431.158, 434.123, 
			     425.697, 426.832, 428.482, 419.536, 414.457, 
			     416.837, 408.455, 396.617, 398.019, 392.691, 
			     393.791, 384.79 , 373.348, 372.505, 367.286, 
			     348.639, 345.557, 336.343, 324.083, 308.243, 
			     301.369, 289.736, 279.774, 270.555, 256.633, 
			     238.905, 224.355, 206.619, 185.351, 160.239, 
			     134.945, 108.963, 85.5356, 57.4822, 26.29, 
			     //8.28956,  2.3158,      0.,      0., 2.28963, 
			     //2.21313, 2.16771, 2.13827, 2.10315, 1.91101
			     0.     , 0.     , 0.     , 0.     , 0.     , 
			     0.     , 0.     , 0.     , 0.     , 0.     };
  */

  //Absorption length from transmittance measurement of a Crytur crystal
  //at CUA. Transmittances from Vladimir.
  const int nAbsl = 61;
  double wlAbsl[nAbsl] = {800.,790.,780.,770.,760.,750.,740.,730.,720.,710.,
			  700.,690.,680.,670.,660.,650.,640.,630.,620.,610.,
			  600.,590.,580.,570.,560.,550.,540.,530.,520.,510.,
			  500.,490.,480.,470.,460.,450.,440.,430.,420.,410.,
			  400.,390.,380.,370.,360.,350.,340.,330.,320.,310.,
			  300.,290.,280.,270.,260.,250.,240.,230.,220.,210.,
			  200.};
  double abslength[nAbsl] = {796.786, 734.961, 789.161, 747.455, 803.081, 
			     798.21 , 725.257, 815.755, 768.501, 743.501, 
			     756.055, 796.412, 809.239, 720.358, 742.132, 
			     725.029, 751.585, 697.754, 713.542, 709.338, 
			     716.528, 692.502, 675.357, 681.138, 649.767, 
			     645.957, 638.725, 614.09 , 608.476, 606.628, 
			     640.36 , 630.599, 565.443, 564.31 , 596.479, 
			     585.38 , 587.096, 580.228, 521.699, 511.217, 
			     490.883, 430.125, 307.714, 150.343, 46.2254, 
			     //11.4686,   0.   ,   0.   ,   0.   ,  2.75252, 
			     //2.53861, 2.7026 , 2.90355, 2.90012, 2.76918, 
			     //2.41157, 0.     , 0.     , 0.     , 0.     , 
			     11.4686,   0.   ,   0.   ,   0.   ,  0., 
			     0.     ,   0.   ,   0.   ,   0.   ,  0., 
			     0.     ,   0.   ,   0.   ,   0.   ,  0., 
			     0.};

  for (G4int i=0; i<nAbsl; i++) wlAbsl[i] *= nanometer;

  for (G4int i=0; i<nAbsl; i++) {
    abslength[i] *= cm;
  };

  G4double kphotAbsl[nAbsl];   //Momenta of optical photons in eV units.
  for (G4int i=0; i<nAbsl; i++) kphotAbsl[i] = hc/wlAbsl[i];

  //PbWO4 refractive index.

  //  G4double rindPbWO4[52];
  //  for (G4int i=0; i<52; i++) {
  ///    rindPbWO4[i] = 2.2;             //PbWO conventional refractive index
  //    rindPbWO4[i] = 2.3;             //test
  //  };

  //Ordinary refractive index for PbWO4 from R.Chipaux, CMS TN/95-184.
  //Note: it is kept constant below 320 nm.
  G4double rindPbWO4[52] = {
    2.230, 2.231, 2.233, 2.235, 2.237, 2.239, 2.242, 2.244, 2.247, 2.250, 
    2.253, 2.257, 2.261, 2.265, 2.269, 2.274, 2.279, 2.284, 2.290, 2.297, 
    2.303, 2.312, 2.321, 2.330, 2.341, 2.354, 2.368, 2.384, 2.402, 2.424, 
    2.449, 2.479, 2.515, 2.561, 2.618, 2.693, 2.795, 2.795, 2.795, 2.795, 
    2.795, 2.795, 2.795, 2.795, 2.795, 2.795, 2.795, 2.795, 2.795, 2.795,
    2.795, 2.795};


  //PBWO4 scinttilation yields, fast and slow components.

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

  //PbWO4 material properties table.
    
  G4MaterialPropertiesTable *PbWO4MPT = new G4MaterialPropertiesTable();
  
  PbWO4MPT -> AddProperty("RINDEX",kphotPbWO4,rindPbWO4,52);
  ///  PbWO4MPT -> AddProperty("ABSLENGTH",kphotPbWO4,abslength,52);
  PbWO4MPT -> AddProperty("ABSLENGTH",kphotAbsl,abslength,nAbsl);

  PbWO4MPT->AddProperty("FASTCOMPONENT",kphotPbWO4_sc_fast,PbWO4_sc_fast,82);
  PbWO4MPT->AddProperty("SLOWCOMPONENT",kphotPbWO4_sc_slow,PbWO4_sc_slow,82);
  PbWO4MPT->AddConstProperty("SCINTILLATIONYIELD", 40000*0.377/100/MeV);
  PbWO4MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  PbWO4MPT->AddConstProperty("FASTTIMECONSTANT", 10.*ns);
  PbWO4MPT->AddConstProperty("SLOWTIMECONSTANT", 30.*ns);
  PbWO4MPT->AddConstProperty("YIELDRATIO", 0.077/(0.077+0.3));

  PbWO4 -> SetMaterialPropertiesTable(PbWO4MPT);

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
  OpticalGlue->AddElement(Si, 1); //exact composition not known

  G4double rindGlue[52];
  for (G4int i=0; i<52; i++) {
    rindGlue[i] = 1.465;
  };

  G4MaterialPropertiesTable *GlueMPT = new G4MaterialPropertiesTable();
  GlueMPT -> AddProperty("RINDEX",kphotPbWO4,rindGlue,52);
  OpticalGlue -> SetMaterialPropertiesTable(GlueMPT);

  // Optical insulation
  //
  ///  density = 1.5;   //approximately
  ///  G4Material* Polymer = new G4Material("Polymer", density, ncomponents=2);
  ///  Polymer->AddElement(C, 1);
  ///  Polymer->AddElement(H, 1);
  density = 1.38;   //from goodfellow
  G4Material* PVF = new G4Material("Polyvinyl fluoride",density,ncomponents=3);
  PVF->AddElement(C, 2);
  PVF->AddElement(H, 3);
  PVF->AddElement(F, 1);

  //Mylar, reflector substrate material.
  //

  G4Material* Mylar = man->FindOrBuildMaterial("G4_MYLAR");

  if (subRefrIndex != 0.) {
    
    //Mylar refractive index.
    G4double rindMylar[52];
    for (G4int i=0; i<52; i++) {
      rindMylar[i] = subRefrIndex;
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

  tedlar_log = new G4LogicalVolume(tedlar_frame,PVF,"Tedlar",0,0,0);

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

  G4double* refKphot;   //Momenta of optical photons in eV units.
  refKphot = new G4double[refNumData];
  for (G4int i=0; i<refNumData; i++) refKphot[i] = hc/refWL[i];

  /*
  if (refFlag != 0) {
    //Specular reflector.

    ReflectorMPT->AddProperty("REALRINDEX",refKphot,refReIndex,refNumData);
    ReflectorMPT->AddProperty("IMAGINARYRINDEX",refKphot,refImIndex,refNumData);

    Reflector -> SetType(dielectric_metal);
    Reflector -> SetFinish(polished);
    Reflector -> SetModel(glisur);
  }
  else {
    // Diffuse reflector, PTFE (Teflon).

    ReflectorMPT -> AddProperty("REFLECTIVITY",refKphot,refRefl,refNumData);

    Reflector -> SetType(dielectric_dielectric);
    Reflector -> SetModel(unified);
    Reflector -> SetFinish(groundfrontpainted);   //Purely Lambertian reflection
  }
  */

  if (refFlag == 1) {

    ReflectorMPT->AddProperty("REALRINDEX",refKphot,refReIndex,refNumData);
    ReflectorMPT->AddProperty("IMAGINARYRINDEX",refKphot,refImIndex,
			      refNumData);

    Reflector -> SetType(dielectric_metal);
    Reflector -> SetFinish(polished);
    Reflector -> SetModel(glisur);
  }
  else {
    // Non metallic reflector, PTFE (Teflon), VM2000.

    ReflectorMPT -> AddProperty("REFLECTIVITY",refKphot,refRefl,refNumData);
    Reflector -> SetType(dielectric_dielectric);
    Reflector -> SetModel(unified);
    if (refFlag == 0)
      Reflector -> SetFinish(groundfrontpainted); //Purely Lambertian reflection
    else
      Reflector -> SetFinish(polishedfrontpainted); //Purely specular reflection
  }

  Reflector -> SetMaterialPropertiesTable(ReflectorMPT);

  G4cout << "===== ReflectorMPT: ============================" << G4endl;
  ReflectorMPT->DumpTable();
  Reflector->DumpInfo();

  if (subRefrIndex == 0.) {
    // Reflective front surface of Mylar.
    new G4LogicalSkinSurface("Reflector",mylar_log,Reflector);
  }
  else {
    // Reflective back surface of Mylar.
    // Tedlar borders Mylar from back. Making it reflective, makes effectively
    // Mylar back surface reflective.
    new G4LogicalSkinSurface("Reflector",tedlar_log,Reflector);
    G4cout << "   subRefrIndex = " << subRefrIndex
	   << ", substarate between crystal and reflector" << G4endl;
  }
  
  //Quantum efficiency of PMT photocathode.
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
