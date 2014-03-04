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
// $Id: ExN02PhysicsList.cc,v 1.22 2006/06/29 17:48:11 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "ExN02PhysicsList.hh"

#include "ExN02SteppingVerbose.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02PhysicsList::ExN02PhysicsList(ExN02SteppingVerbose* verbosity)
  :  G4VUserPhysicsList(),
     mySV(verbosity),
     // default values for printing verbose only
     MinKineticEnergy(10.0*eV),
     theNuclearTable("ICRU_R49"),
     theProtonTable("ICRU_R49p"),
     dRoverRange(0.20)
{
  defaultCutValue = .01*nm;
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02PhysicsList::~ExN02PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructBaryons();
  ConstructIons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PhysicsList::ConstructIons()
{
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
  // mandatory to be able to shoot ions in PrimaryGenerator particleTable->GetIon
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PhysicsList::ConstructLeptons()
{
  G4Electron::ElectronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PhysicsList::ConstructBaryons()
{
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RPNuclearRecoil01.hh"
#include "RPNuclearRecoil02.hh"

#include "G4hIonisation.hh"
#include "G4hLowEnergyIonisation.hh"

#include "G4UserLimits.hh"
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
     
    // User Verbose output class
    //    ExN02SteppingVerbose* verbosity = new ExN02SteppingVerbose;
    
    //////////////////////////////////////////////////////////////////////////////
    // RPNuclearRecoil initialization and behaviour control via its public methods
    //////////////////////////////////////////////////////////////////////////////
    
    RPNuclearRecoil01* nuclearrecoils01 = new RPNuclearRecoil01;
    //    pmanager->AddContinuousProcess(nuclearrecoils01);

    RPNuclearRecoil02* nuclearrecoils02 = new RPNuclearRecoil02;
     
    // Set the threshold ion energy, THE CUT, for both primary and secondaries
    // Default is 10.0*eV.
    // DON'T FORGET TO ALSO CHANGE IN G4hLowEnergyIonisation IF < 10eV !!!
    // Anyway, you should ALWAYS change to keep both processes sinchronized for more precise results
    //    nuclearrecoils02->SetMinKineticEnergy(MinKineticEnergy=15.0*eV);

    // define the nuclear stopping parametrisation method 
    // List of available nuclear stopping power parametrisation tables:
    // "ICRU_R49" (default), "Ziegler1977", "Ziegler1985"
    //    nuclearrecoils02->SetNuclearStoppingPowerModel(theNuclearTable="Ziegler1985");
    // DON'T FORGET TO ALSO CHANGE IN SteppingVerbose just to avoid Energy Conservation Warnings
    mySV->SetNuclearStoppingPowerModel(theNuclearTable);

    // propose a step length corresponding to NumRecoilsProposed recoils per step
    // Default is Off.
    nuclearrecoils02->SetStepLimitOff();
    
    // force NumRecoilsProposed per step independent of step length
    // Default is On. 
    nuclearrecoils02->ForceNumRecoilsPerStepOn();

    // Propose a number of nuclear recoils to be generated during a step
    // Default is 1
    //    nuclearrecoils02->ProposeNumRecoilsPerStep(2);

    // Gaussian sample recoil energy and number (only if NumRecoilsProposed > 1)
    // Default is Off.
    //    nuclearrecoils02->SetNuclearStoppingFluctuationsOn();

    pmanager->AddContinuousProcess(nuclearrecoils02);


    /////////////////////////////////////////////////////////////////////////////////////
    // G4hLowEnergyIonisation initialization and behaviour control via its public methods
    /////////////////////////////////////////////////////////////////////////////////////
    
    G4hLowEnergyIonisation* lowenergy = new G4hLowEnergyIonisation;
    
    /* ELECTRONIC */

    // List of available electronic stopping power parametrisation tables:
    // "ICRU_R49p" (default), "ICRU_R49He", "Ziegler1977p", "Ziegler1985p", "Ziegler1977He", "SRIM2000p"
    //    lowenergy->SetElectronicStoppingPowerModel(particle,theProtonTable="SRIM2000p");

    // Definition of the boundary proton energy.

    // For higher energies Bethe-Bloch formula is used, for lower energies 
    // a parametrisation of the energy losses is performed. 
    // Default is min(100 MeV,theProtonModel->HighEnergyLimit
    //    lowenergy->SetHighEnergyForProtonParametrisation(2.*MeV);

    // For lower energies the Free Electron Gas (Linhard) model is used for the energy losses.
    // Default is 1 keV.
    //    lowenergy->SetLowEnergyForProtonParametrisation (5.0*keV);//4.0);

    /* NUCLEAR */
    
    // List of available nuclear stopping power parametrisation tables:
    // "ICRU_R49" (default), "Ziegler1977", "Ziegler1985"
    //    lowenergy->SetNuclearStoppingPowerModel(theNuclearTable="Ziegler1977");

    // Set the nuclear stopping power off. MANDATORY IF YOU USE THE NUCLEAR RECOILS PROCESS
    // Default is On.
    lowenergy->SetNuclearStoppingOff();
    
    /* GENERAL */

    // print verbosity of G4hLowEnergyIonisation: 3=everything 2=theDEDXpTable
    // Default is 0 (Silent).
    //    lowenergy->SetVerboseLevel(2);

    // approximate maximum step allowed as a fraction of the expected range value
    // Default is 0.20
    lowenergy->SetdRoverRange(dRoverRange=0.01);

    // no fluctuations: continuous slowing down approximation
    // Default is On.
    lowenergy->SetEnlossFluc(false);

    // This method switch off calculation of the Barkas and Bloch effects
    // Default is On.
    //    lowenergy->SetBarkasOff();
    
    pmanager->AddProcess(lowenergy,-1,2,2); // 2 e' a ordem
    

    ///////////////////////
    // standard EM package
    ///////////////////////

    //    pmanager->AddProcess(new G4hIonisation, -1,2,2);
    
    /*
      G4VProcess* theprotoncut = new G4UserLimits();
      theprotoncut->SetUserMinEkine(1.0*MeV);
      pmanager->AddProcess(theprotoncut);
    */
    
    // step limit
    pmanager->AddProcess(new G4StepLimiter,       -1,-1,3);         
    //    pmanager->AddProcess(new G4UserSpecialCuts,   -1,-1,4);  
    
  }
  G4cout << "\n**********************************************************"
	 << "\nhLowEIoni maximum step length = " << dRoverRange*100
	 << "% of the expected range"
	 << "\n**********************************************************"
	 << G4endl;

  FILE *fp=fopen("quenching.dat","a");
  fprintf(fp,"# hLowEIoni:\tdRoverRange = %g\t\ttheProtonTable  = %s\n",dRoverRange,theProtonTable.data());
  fprintf(fp,"# contNucRcls:\tMinKineticEnergy = %g eV \ttheNuclearTable = %s\n",MinKineticEnergy/eV,theNuclearTable.data());
  fclose(fp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PhysicsList::SetCuts()
{
  // G4VUserPhysicsList::SetCutsWithDefault method sets 
  // the default cut value for all particle types 
  //
  //  SetCutsWithDefault();

  G4double lowlimit=100*GeV; // big cut to prevent the production of delta electrons
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100.2*GeV);
     
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
