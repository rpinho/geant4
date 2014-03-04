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
// $Id: ExN02SteppingAction.cc,v 1.9 2006/06/29 17:48:18 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02SteppingAction.hh"
#include "G4SteppingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02SteppingAction::ExN02SteppingAction(ExN02EventAction* event_action)
  :myEA(event_action)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02SteppingAction::UserSteppingAction(const G4Step* fStep)
{
  const G4Track* fTrack = fStep->GetTrack();
  G4int step = fTrack->GetCurrentStepNumber();
  G4int track = fTrack->GetTrackID();
  G4int parent = fTrack->GetParentID();

  // the primary
  if (track == 1) {
    
    // the first step
    if(step == 1) { // don't know why but step == 0 is not working...
      
      // get the initial recoil energy
      G4double kineticEnergy = fTrack->GetKineticEnergy() - fStep->GetDeltaEnergy();
      myEA->InitialKineticEnergy(kineticEnergy);
      
      // get the length of the first step
      myEA->FirstStepLength(fStep->GetStepLength());
    }
    
    // the last step
    if(fTrack->GetKineticEnergy() == 0)
      myEA->PrimaryRange(fTrack->GetTrackLength());
  }

  // the secondaries
  else {

    // the first step
    if(step == 1) {

      // update total number of recoils
      if(track > myEA->TotalNumberOfTracks()) myEA->TotalNumberOfTracks(track);
  
      // update total number of "parents"
      myEA->TotalNumberOfParents(parent);

      // update total number of primary "sons"
      if(parent == 1 && track > myEA->TotalNumberOfPrimarySons()) myEA->TotalNumberOfPrimarySons(track);
    }

    // the last step
    if(fTrack->GetKineticEnergy() == 0)
      myEA->SecondariesRange(fTrack->GetTrackLength());
  }
  
  // update quenching
  myEA->Add2TotalEnergyDeposit( fStep->GetTotalEnergyDeposit() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

