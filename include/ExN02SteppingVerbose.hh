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
// $Id: ExN02SteppingVerbose.hh,v 1.8 2006/06/29 17:47:50 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose.
//   It shows how to extract informations during the tracking of a particle.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN02SteppingVerbose;

#ifndef ExN02SteppingVerbose_h
#define ExN02SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

#include "G4VLowEnergyModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN02SteppingVerbose : public G4SteppingVerbose 
{
public:
   
  ExN02SteppingVerbose();
  ~ExN02SteppingVerbose();
  
  void StepInfo();
  void TrackingStarted();

  void SetNuclearStoppingPowerModel(const G4String& dedxTable)
                               {theNuclearTable = dedxTable;};
  // This method defines the nuclear stopping parametrisation method 
  // via the name of the table. 
  // Default is "ICRU_49".

private:

  // Low Energy Electromagnetic models
  G4VLowEnergyModel* theNuclearStoppingModel;
  G4VLowEnergyModel* theIonEffChargeModel;

  G4String theNuclearTable;  // name of parametrisation table of nuclear stopping power

  G4double factorPDG2AMU;    // Factor to convert PDG mass unit
                             // into AMU

  G4double theZieglerFactor; // Factor to convert the Stopping Power 
                             // unit [ev/(10^15 atoms/cm^2]
                             // into the Geant4 dE/dx unit
  
  G4double MinKineticEnergy; /* MUITO IMPORTANTE!! */

  G4double MinkineE;         /* a energia cinetica minima de Tracking 
				da particula, apos a qual a matamos, 
				tipo User Limits */
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
