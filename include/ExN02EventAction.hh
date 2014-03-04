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
// $Id: ExN02EventAction.hh,v 1.8 2006/06/29 17:47:35 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef ExN02EventAction_h
#define ExN02EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>

class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN02EventAction : public G4UserEventAction
{
public:
  ExN02EventAction();
  ~ExN02EventAction();
  
public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

  void Add2TotalEnergyDeposit(G4double eloss)
  {TotalEnergyDeposit += eloss;};

  G4double InitialKineticEnergy()
  {return RecoilEnergy;};
  void InitialKineticEnergy(G4double energy)
  {RecoilEnergy = energy;};

  G4int TotalNumberOfTracks()
  {return NumberOfRecoils;};
  void TotalNumberOfTracks(G4int trackID)
  {NumberOfRecoils = trackID;};

  G4int TotalNumberOfParents()
  {return NumberOfParents;};
  void TotalNumberOfParents(G4int trackID)
  {parents.push_back(trackID);};

  G4int TotalNumberOfPrimarySons()
  {return NumberOfPrimarySons;};
  void TotalNumberOfPrimarySons(G4int trackID)
  {NumberOfPrimarySons = trackID;};

  G4double FirstStepLength()
  {return step;};
  void FirstStepLength(G4double length)
  {step = length;};

  G4double PrimaryRange()
  {return range;};
  void PrimaryRange(G4double length)
  {range = length;};

  G4double SecondariesRange()
  {return recoilRange;};
  void SecondariesRange(G4double length)
  {recoilRange += length;};

private:  
  G4double TotalEnergyDeposit;
  G4double RecoilEnergy;
  G4int NumberOfRecoils;
  G4int NumberOfParents;
  G4int NumberOfPrimarySons;
  G4double step;
  G4double range;
  G4double recoilRange;
  FILE *fp;
  std::vector<G4int> parents;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
