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
// $Id: G4Cerenkov.hh,v 1.8 2006/06/29 19:55:31 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// 
////////////////////////////////////////////////////////////////////////
// Cerenkov Radiation Class Definition 
////////////////////////////////////////////////////////////////////////
//
// File:        G4Cerenkov.hh  
// Description:	Continuous Process -- Generation of Cerenkov Photons
// Version:     2.0
// Created:     1996-02-21
// Author:      Juliet Armstrong
// Updated:     2005-07-28 add G4ProcessType to constructor
//              1999-10-29 add method and class descriptors
//              1997-04-09 by Peter Gumplinger
//              > G4MaterialPropertiesTable; new physics/tracking scheme
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef RPNuclearRecoil01_h
#define RPNuclearRecoil01_h 1

/////////////
// Includes
/////////////

#include "G4hLowEnergyIonisation.hh"

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VContinuousProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh" 
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"

// Class Description:
// Continuous Process -- Generation of nuclear recoils
// Class inherits publicly from G4VContinuousProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class RPNuclearRecoil01 : public G4VContinuousProcess  
{

public: // Without description

	////////////////////////////////
	// Constructors and Destructor
	////////////////////////////////

        RPNuclearRecoil01(const G4String& processName = "recoil"); 
  
        ~RPNuclearRecoil01();	

        ////////////
        // Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition&);
        // True for all charged hadrons/ions
    
        G4double GetContinuousStepLimit(const G4Track&,
					G4double,
                                        G4double,
                                        G4double&); 
        // Returns the continuous step limit
  
        G4VParticleChange* AlongStepDoIt(const G4Track& trackData ,
					 const G4Step& stepData ) ;

 
        G4double MinKineticEnergy;
        G4VLowEnergyModel* theNuclearStoppingModel;

};
  
////////////////////
// Inline methods
////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double RPNuclearRecoil01::GetContinuousStepLimit(
                                        const G4Track&,
                                              G4double,
                                              G4double,
                                              G4double&)
{
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool RPNuclearRecoil01::IsApplicable(
                                const G4ParticleDefinition& particle)
{
   return(particle.GetPDGCharge() != 0.0
       && particle.GetPDGMass() > proton_mass_c2*0.1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
