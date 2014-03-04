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
// $Id: ExN02EventAction.cc,v 1.11 2006/06/29 17:48:05 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "ExN02EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN02EventAction::ExN02EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN02EventAction::~ExN02EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExN02EventAction::BeginOfEventAction(const G4Event*)
{
  TotalEnergyDeposit=0.0;
  RecoilEnergy=0.0;
  NumberOfRecoils=0;
  NumberOfParents=0;
  NumberOfPrimarySons=0;
  step=0.0;
  range=0.0;
  recoilRange=0.0;
  parents.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExN02EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  
  G4double quenching = TotalEnergyDeposit / RecoilEnergy;
  
  if(NumberOfRecoils) NumberOfRecoils-=1;
  if(NumberOfPrimarySons) NumberOfPrimarySons-=1;
  if(NumberOfRecoils) recoilRange /= NumberOfRecoils;

  std::sort(parents.begin(),parents.end());
  parents.erase(std::unique(parents.begin(),parents.end()),parents.end());

  // periodic printing
  G4int prec = G4cout.precision(3);
  if (event_id < 100 || event_id%100 == 0) {
    G4cout << "\n>>> Event " << evt->GetEventID() << G4endl;
    G4cout << "Initial Kinetic Energy = " << RecoilEnergy/keV << " keV" 
	   << "\tTotal Energy Deposit = " << TotalEnergyDeposit/keV << " keV" 
	   << "\tNuclear Quenching = " << quenching
	   << "\nTotal Number of Recoils = " << NumberOfRecoils
	   << "\tTotal Number of Parents (and sons of the primary recoil) = " << parents.size()
	   << " (" << NumberOfPrimarySons << ")"
	   << "\nPrimary's First Step Length = " << G4BestUnit(step,"Length")
	   << "\tPrimary Range = " << G4BestUnit(range,"Length")
	   << "\tSecondaries Average Range = " << G4BestUnit(recoilRange,"Length")
	   << G4endl << G4endl;
    
    fp=fopen("quenching.dat","a");
    fprintf(fp,"%g\t%lf\t%d\t%d (%u)\t\t%.0lf\t%.0lf\t%.0lf\n",RecoilEnergy/keV,quenching,NumberOfRecoils,parents.size(),NumberOfPrimarySons,step/angstrom,range/angstrom,recoilRange/angstrom);
    fclose(fp);
  }
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
