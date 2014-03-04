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
// $Id: ExN02PrimaryGeneratorAction.cc,v 1.7 2006/06/29 17:48:13 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include "G4GeneralParticleSource.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02PrimaryGeneratorAction::ExN02PrimaryGeneratorAction()
{
  // replacing G4ParticleGun with the more powerful G4GeneralParticleSource
  // don't forget to also replace in ExN02PrimaryGeneratorAction.hh
  //  particleGun = new G4GeneralParticleSource();

  // the simple G4ParticleGun
  particleGun = new G4ParticleGun();
  

  ///////////////////////////////////////////////////////////////////////////////////
  // The properties of the incident source can also be defined in the macro file, 
  // overriding the values predefined in the next lines. Consider them default values
  ///////////////////////////////////////////////////////////////////////////////////
 
  
  // predefine the particle type
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // alpha or proton
  //  G4ParticleDefinition* particle = particleTable->FindParticle("alpha");

  // generic ion: Xe: Z=54, A=131; Ge: Z=32, A=73; Ar: Z=18, A=40;
  // Si: Z=14, A=28; Al: Z=13, A=27; N: Z=7, A=14
  G4ParticleDefinition* particle = particleTable->GetIon(54,131,0.0);
  
  particleGun->SetParticleDefinition(particle);

 
  // predefine the particle direction

  // generating with random (isotropic) direction
  particleGun->SetParticleMomentumDirection(G4ThreeVector(G4UniformRand(),G4UniformRand(),G4UniformRand()));
  
  // generating with direction along the z axis
  //  particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));

 
  // predefine the particle initial kinetic energy
  
  particleGun->SetParticleEnergy(100*keV);
   

  // predefine the number of particles to generate
  
  //  G4int n_particle = 1;
  //  particleGun->SetNumberOfParticles(n_particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02PrimaryGeneratorAction::~ExN02PrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  // predefine the particle initial position
  
  // generating at the center of the box,
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

