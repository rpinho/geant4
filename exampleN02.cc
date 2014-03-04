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
// $Id: exampleN02.cc,v 1.12 2006/06/29 17:47:25 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02DetectorConstruction.hh"
#include "ExN02PhysicsList.hh"
#include "ExN02PrimaryGeneratorAction.hh"
#include "ExN02RunAction.hh"
#include "ExN02EventAction.hh"
#include "ExN02SteppingAction.hh"
#include "ExN02SteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // the output data files
  char f[]="dEdx.dat",f1[]="LE_proton.dat",f2[]="LE_alfa.dat",f3[]="LE_ion.dat";
  FILE *fp,*fp1,*fp2,*fp3,*fp4,*fp5;

  // open and close them to delete the already existing ones
  fp=fopen(f,"w");fclose(fp);
  fp1=fopen("nloss.dat","w");fclose(fp1);
  fp2=fopen("eloss.dat","w");fclose(fp2);
  fp3=fopen("quenching.dat","w");fclose(fp3);
  fp4=fopen("dEdx_beta.dat","w");fclose(fp4);
  fp5=fopen("step.dat","w");fclose(fp5);
  
  // choose the Random engine and seed it
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  long seed;
  time(&seed); /* num. de segundos desde 1/1/1970 */
  CLHEP::HepRandom::getTheEngine()->setSeed(seed,1);

  // User Verbose output class
  //
  ExN02SteppingVerbose* verbosity = new ExN02SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // Run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // User Initialization classes (mandatory)
  //
  ExN02DetectorConstruction* detector = new ExN02DetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new ExN02PhysicsList(verbosity);
  runManager->SetUserInitialization(physics);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
   
  // User Action classes
  //
  //G4VUserPrimaryGeneratorAction* gen_action = new ExN02PrimaryGeneratorAction(detector);
  //runManager->SetUserAction(gen_action);
  //
  G4UserRunAction* run_action = new ExN02RunAction;
  runManager->SetUserAction(run_action);
  //
  ExN02EventAction* event_action = new ExN02EventAction;
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action = new ExN02SteppingAction(event_action);
  runManager->SetUserAction(stepping_action);

  // Initialize G4 kernel
  //
  runManager->Initialize();

  //
  G4VUserPrimaryGeneratorAction* gen_action = new ExN02PrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);

  // Get the pointer to the User Interface manager
  //
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  if(argc==1)  // Define (G)UI terminal for interactive mode
  { 
    // G4UIterminal is a (dumb) terminal
    //
    G4UIsession * session = 0;
    session = new G4UIterminal(new G4UItcsh);      
    /*UI->ApplyCommand("/run/verbose 0");
    UI->ApplyCommand("/event/verbose 0");
    UI->ApplyCommand("/tracking/verbose 1");
    for(int i=200;i>0;i--) {
      UI->ApplyCommand("/gun/energy i keV");
      UI->ApplyCommand("/run/beamOn 1");
      }*/
    session->SessionStart();
    delete session;
  }
  else   // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;
  delete verbosity;

  //  printf("\nImpresso o ficheiro %s\n\n",f);

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

