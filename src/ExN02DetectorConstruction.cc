//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: ExN02DetectorConstruction.cc,v 1.18 2006/06/29 17:48:00 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "ExN02DetectorConstruction.hh"
#include "ExN02DetectorMessenger.hh"
//#include "ExN02ChamberParameterisation.hh"
#include "ExN02MagneticField.hh"
//#include "ExN02TrackerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
//#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN02DetectorConstruction::ExN02DetectorConstruction()
  :solidWorld(0),  logicWorld(0), physiWorld(0), 
   TargetMater(0), fpMagField(0), fWorldLength(0.)
{
  fpMagField = new ExN02MagneticField();
  detectorMessenger = new ExN02DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN02DetectorConstruction::~ExN02DetectorConstruction()
{
  delete fpMagField;
  delete detectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* ExN02DetectorConstruction::Construct()
{
  //--------- Material definition ---------

  G4double a,density, temperature, pressure, abundance,fractionmass;
  G4int nel,z,n,ncomponents;
  G4String name,symbol;
  
  // ---------------
  //    Isotopes
  // ---------------
  // ----
  //  Xe - Stables
  // ----
  G4Isotope* Xe124 = new G4Isotope(name="Xe124",z=54,n=124,a=123.9061*g/mole);
  G4Isotope* Xe126 = new G4Isotope(name="Xe126",z=54,n=126,a=125.9042*g/mole);
  G4Isotope* Xe128 = new G4Isotope(name="Xe128",z=54,n=128,a=127.9035*g/mole);
  G4Isotope* Xe129 = new G4Isotope(name="Xe129",z=54,n=129,a=128.9048*g/mole);
  G4Isotope* Xe130 = new G4Isotope(name="Xe130",z=54,n=130,a=129.9035*g/mole);
  G4Isotope* Xe131 = new G4Isotope(name="Xe131",z=54,n=131,a=130.9051*g/mole);
  G4Isotope* Xe132 = new G4Isotope(name="Xe132",z=54,n=132,a=131.9042*g/mole);
  G4Isotope* Xe134 = new G4Isotope(name="Xe134",z=54,n=134,a=133.9054*g/mole);
  G4Isotope* Xe136 = new G4Isotope(name="Xe136",z=54,n=136,a=135.9072*g/mole);
  
  // ---------------
  //    Elements
  // ---------------
  // ------------
  //  Natural Xe 
  // ------------
  G4Element* elXe = new G4Element(name="Xenon(el)",symbol="Xe",ncomponents=9);
  elXe->AddIsotope(Xe124,abundance=0.096*perCent);
  elXe->AddIsotope(Xe126,abundance=0.090*perCent);
  elXe->AddIsotope(Xe128,abundance=1.92*perCent);
  elXe->AddIsotope(Xe129,abundance=26.44*perCent);
  elXe->AddIsotope(Xe130,abundance=4.08*perCent);
  elXe->AddIsotope(Xe131,abundance=21.18*perCent);
  elXe->AddIsotope(Xe132,abundance=28.69*perCent);
  elXe->AddIsotope(Xe134,abundance=10.44*perCent);
  elXe->AddIsotope(Xe136,abundance=8.87*perCent);

  // ---------------
  //    Materials
  // ---------------
  // -------
  //  Xenon
  // -------
  density = 2.9*g/cm3;
  // Liquid Xe density used in detectors: 2.3, 2.9, 2.95, 3.0, 3.05, 3.1 g/cm3
  G4Material* Xe = new G4Material(name="Xenon",density,nel=1);
  Xe->AddElement(elXe,fractionmass=1.);
  // mass = 131.30 g/mole (?133.67?)

  // Lead
  G4Material* Pb = 
  new G4Material("Lead", z=82, a= 207.19*g/mole, density= 11.35*g/cm3);
    
  // Copper
  G4Material* Cu = 
    new G4Material("Copper", z=29, a=63.55*g/mole, density= 8.96*g/cm3);

  // Aluminium
  G4Material* Al = 
  new G4Material("Aluminium", z=13, a= 26.98*g/mole, density= 2.70*g/cm3);

  // Silicon
  G4Material* Si = 
  new G4Material("Silicon", z=14, a= 28.09*g/mole, density= 2.33*g/cm3);

  // Argon
  G4Material* Ar = 
  new G4Material("Argon", z=18, a= 39.95*g/mole, density= 1.39*g/cm3);

  // Germanium
  G4Material* Ge = 
  new G4Material("Germanium", z=32, a= 72.64*g/mole, density= 5.32*g/cm3);

  // Xenon gas
  G4Material* XeGas = 
  new G4Material("XenonGas", z=54, a=131.29*g/mole, density= 5.458*mg/cm3,
                 kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);
 

  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;


  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World is Target
  //------------------------------ 

  // Size of the principal geometrical components (solids)
  fWorldLength  = 500 * m;
  G4double HalfWorldLength = 0.5*fWorldLength;
  solidWorld = new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  
  // (default) Material of the Target
  TargetMater = Xe;
  logicWorld = new G4LogicalVolume(solidWorld, TargetMater , "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
				 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // copy number

  // Verbosity
  G4cout << "World (Target) is a " << G4BestUnit(fWorldLength,"Length") << "square box of " 
         << TargetMater->GetName() << " (default)" << G4endl;			 
  FILE *fp=fopen("quenching.dat","a");
  fprintf(fp,"# (default) TargetMaterial is %s\n",TargetMater->GetName().data());
  fclose(fp);
  
  /*
  //--------- Visualization attributes -------------------------------

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicWorld  ->SetVisAttributes(BoxVisAtt);  
  logicTarget ->SetVisAttributes(BoxVisAtt);
  logicTracker->SetVisAttributes(BoxVisAtt);
  
  G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  logicChamber->SetVisAttributes(ChamberVisAtt);
  
  //--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in N02PhysicsList how to setup the processes
  // G4StepLimiter or G4UserSpecialCuts).
    
  // Sets a max Step length in the tracker region, with G4StepLimiter
  //
  G4double maxStep = 0.5*ChamberWidth; 
  logicTracker->SetUserLimits(new G4UserLimits(maxStep));
  
  // Set additional contraints on the track, with G4UserSpecialCuts
  //
  // G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  */
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExN02DetectorConstruction::setTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);
  if (pttoMaterial && materialName != TargetMater->GetName()) {
    TargetMater = pttoMaterial;
    logicWorld->SetMaterial(pttoMaterial); // World is Target
    G4cout << "\n------> Geometry (target material) has changed:\n" 
	   << "\tWorld (Target) is now a " << G4BestUnit(fWorldLength,"Length")
	   << "square box of " << materialName << G4endl << G4endl;
    FILE *fp=fopen("quenching.dat","a");
    fprintf(fp,"# TargetMaterial was changed to %s\n",materialName.data());
    fclose(fp);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExN02DetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetFieldValue(fieldValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
