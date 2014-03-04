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
// $Id: ExN02SteppingVerbose.cc,v 1.12 2006/06/29 17:48:21 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02SteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

#include "G4hLowEnergyIonisation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "G4hIonEffChargeSquare.hh"
#include "G4hNuclearStoppingModel.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include <stdlib.h>
#include <math.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02SteppingVerbose::ExN02SteppingVerbose()
  : theNuclearTable("ICRU_R49")
{
  // Constants - see description in the header file
  factorPDG2AMU    = 1.007276/proton_mass_c2 ;
  theZieglerFactor = eV*cm2*1.0e-15 ; 
  MinKineticEnergy = 10.0*eV ;
  MinkineE         = 0*MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02SteppingVerbose::~ExN02SteppingVerbose()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02SteppingVerbose::StepInfo()
{

  /************** Track and Step Data **************/

  CopyState();

  /*
    G4TrackStatus status = fTrack->GetTrackStatus();
    G4cout << status << G4endl;
    if(kineticEnergy<MinkineE) {
    fTrack->SetTrackStatus(fStopAndKill);             // kill the primary
    fTrack->SetTrackStatus(fKillTrackAndSecondaries); // kill both the primary and the secundaries
    G4cout << "        ********************** "
	   << "A particula vai ser morta agora abaixo de KineE = " << G4BestUnit(MinkineE,"Energy")
	   << " **************************" << G4endl;
	   }
  */
  
  /* a energia no FINAL do step */
  G4double finalT = fTrack->GetKineticEnergy();
  //  G4cout << "finalT = fTrack->GetKineticEnergy() = " << G4BestUnit(finalT,"Energy") << G4endl;

  /* a perca TOTAL de energia durante o step (electronic + nuclear) NOTA: o delta e' negativo */
  G4double ionloss = - fStep->GetDeltaEnergy();
  //  G4cout << "ionloss = - fStep->GetDeltaEnergy() = " << G4BestUnit(ionloss,"Energy") << G4endl;

  /* a energia no INICIO do step */
  G4double kineticEnergy = finalT + ionloss;
  //  G4cout << "kineticEnergy = finalT + ionloss = " << G4BestUnit(kineticEnergy,"Energy") << G4endl;

  /* a DEPOSICAO local de energia durante o step (electronic, se lowenergy->SetNuclearStoppingOff() no PhysicsList) */
  G4double eloss = fStep->GetTotalEnergyDeposit();
  //  G4cout << "eloss = fStep->GetTotalEnergyDeposit() = " << G4BestUnit(eloss,"Energy") << G4endl;

  /* o tamanho do step dx (= PhysicalStep) */
  const G4double step = fStep->GetStepLength();
  //  G4cout << "step = fStep->GetStepLength() = " << G4BestUnit(step,"Length") << G4endl;

  /* ELECTRONIC Stopping Power, NOTA: depende do ProposeLocalEnergyDeposit(eloss) definido no G4hLowEnergyIonisation */
  G4double dEdx = eloss/step;
  //  G4cout << "dEdx = eloss/step = " << G4BestUnit(dEdx,"Energy/Length") << G4endl;

  // get the change in momentum of the primary particle. Note: delta is negative
  G4ThreeVector DeltaMomentum = - fStep->GetDeltaMomentum();

  // get the change in direction of the primary particle
  G4ThreeVector DeltaDirection = DeltaMomentum.unit();
    
  
  /************* get material parameters needed for the energy loss calculation **************/

  const G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  //  G4cout <<  "numOfCouples = " << numOfCouples << G4endl;
  G4double deltaCut = (*(theCoupleTable->GetEnergyCutsVector(1)))[0];
  //  G4cout <<  "deltaCut = " << G4BestUnit(deltaCut,"Energy") << G4endl;
  const G4MaterialCutsCouple* couple = fTrack->GetMaterialCutsCouple();
  //  G4double deltaCutNow = cutForDelta[couple->GetIndex()] ;
  const G4Material* material= couple->GetMaterial();
  G4String materialName = material->GetName();
  //  G4cout << "Material Name = " << materialName << G4endl;
  G4double materialDensity = material->GetDensity();
  //  G4cout << "Material Density = " << G4BestUnit(materialDensity,"Volumic Mass") << G4endl;
  G4double materialAtomicMass = material->GetA();
  //  G4cout << "Material Atomic Mass = " << materialAtomicMass << G4endl;
  G4double TotNbOfAtomsPerVolume = material->GetTotNbOfAtomsPerVolume(); // 13.3007 nm^-3
  //  G4cout << "TotNbOfAtomsPerVolume = " << TotNbOfAtomsPerVolume/(1/(nm*nm*nm)) << G4endl;
  G4double materialAtomicNumber = material->GetZ();
  //  G4cout << "Material Atomic Number = " << materialAtomicNumber << G4endl;
  

  /************* get particle *************/
  
  const G4DynamicParticle* aParticle = fTrack->GetDynamicParticle() ;
  G4ParticleDefinition* particle = fTrack->GetDefinition();
  // == aParticle->GetDefinition()
  G4String particleName = particle->GetParticleName();
  //  G4cout << "Particle Name = " << particleName << G4endl;
  // The mass of the particle, in units of equivalent energy
  G4double particleMass = particle->GetPDGMass();
  //  G4cout << "Particle Mass = " << particleMass << G4endl;
  G4double particleCharge = particle->GetPDGCharge();
  // == aParticle->GetCharge(), i.e., nao e' effective charge 
  // == std::abs((particle->GetPDGCharge())/eplus)
  // == z1
  //  G4cout << "Z = " << particleCharge << G4endl;
    
  
  /************** nuclear stopping power ****************/
  
  //  theNuclearTable = "Ziegler1985"; // ICRU_R49 (default), Ziegler1985 ou Ziegler1977
  theNuclearStoppingModel = new G4hNuclearStoppingModel(theNuclearTable) ;
  
  G4double nloss = theNuclearStoppingModel->TheValue(particle, material, kineticEnergy);
  // != theNuclearStoppingModel->TheValue(G4DynamicParticle* aParticle, material) <=> finalT != kineticEnergy
  // This is the only reason for using one and not the other, as NuclearStopping doesn't use the effective 
  // charge approach, for that would make it velocity dependent, and we use the static approximation - Linhard
  //  G4cout << "nloss = " << G4BestUnit(nloss*step,"Energy") << G4endl;
  
  // Average number of recoils <n> assuming MinKineticEnergy as W, i.e., average energy per recoil
  //  G4cout << "<n> = " << nloss*step/MinKineticEnergy << G4endl;
  
  // Projectile nucleus
  G4double z1 = particleCharge;
  G4double m1 = particleMass*factorPDG2AMU;
  // Target nucleus
  G4double z2 = materialAtomicNumber;
  G4double m2 = materialAtomicMass*mole/g;
  //  G4cout << z1 << " " << m1 << " " << z2 << " " << m2 << G4endl;
    
  // Linhard units
  // energy in keV
  G4double energy = kineticEnergy/keV ;
  // reduced mass
  G4double rm = (m1 + m2) * ( std::pow(z1, .23) + std::pow(z2, .23) ) ;
  //  G4cout << "reduced mass = " << rm << G4endl;
  // reduced energy
  G4double er = 32.536 * m2 * energy / ( z1 * z2 * rm ) ;
  //  G4cout << "reduced energy = " << er << G4endl;
    
    
  /************** effective charge model ****************/
  theIonEffChargeModel = new G4hIonEffChargeSquare("Zigler1988");
  G4double chargeSquare = theIonEffChargeModel->TheValue(particle,material,kineticEnergy);
  
    
  /************** normalizacao para comparar com os graficos do Ziegler ************/
  /*
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* normpart = particleTable->GetIon(2,4,0.0); // Al: Z=13, A=27
    G4double chargeSquareNorm = theIonEffChargeModel->TheValue(normpart,material);
    G4double Z2norm = chargeSquareNorm/chargeSquare;
  */
  
  
  /************** secondary particles *****************/
  // Note: (*fSecondary) is a G4Track
  G4int tN2ndariesTot = fN2ndariesAtRestDoIt + fN2ndariesAlongStepDoIt + fN2ndariesPostStepDoIt;
  G4double t2ndariesEnergy = 0.0;
  if(tN2ndariesTot>0)
    for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; lp1<(*fSecondary).size(); lp1++)
      t2ndariesEnergy += (*fSecondary)[lp1]->GetKineticEnergy();
  
  
  /************** a copy of G4hLowEnergyIonisation::DeltaRaysEnergy **************/
  G4double dloss = 0.0 ;
  G4double electronDensity = material->GetElectronDensity();
  G4double excEnergy = material->GetIonisation()->GetMeanExcitationEnergy();
  //  G4cout << "excEnergy = " << G4BestUnit(excEnergy,"Energy") << G4endl;
  G4double tau = kineticEnergy/particleMass ;
  G4double rateMass = electron_mass_c2/particleMass ;
  // some local variables
  G4double gamma,bg2,beta2,beta,tmax,X ;
  gamma = tau + 1.0 ;
  bg2 = tau*(tau+2.0) ;
  beta2 = bg2/(gamma*gamma) ;
  beta=sqrt(beta2);
  tmax = 2.*electron_mass_c2*bg2/(1.0+2.0*gamma*rateMass+rateMass*rateMass) ;
  // Validity range for delta electron cross section
  deltaCut = std::max(deltaCut, excEnergy);
  if ( deltaCut < tmax) {
    X = deltaCut / tmax ;
    dloss = ( beta2 * (X - 1.0) - std::log(X) ) * twopi_mc2_rcl2 * electronDensity / beta2 ;
  }
  
  
  /************* print data files **************/      
  
  char f[]="dEdx.dat",f1[]="LE_proton.dat",f2[]="LE_alfa.dat",f3[]="LE_ion.dat";
  FILE *fp,*fp1,*fp2,*fp3,*fp4;
  fp=fopen(f,"a");
  fp1=fopen("nloss.dat","a");
  fp2=fopen("eloss.dat","a");
  fp3=fopen("dEdx_beta.dat","a");
  fp4=fopen("step.dat","a");

  /***************** output units ****************************/

  G4double x=kineticEnergy; /* MeV */

  //  x/=keV; // keV
  //  x/=m1;  // energy/amu (energy per nucleon)
  //  x=std::sqrt(er); /* Linhard's reduced energy */
  
  G4double y=dEdx;  /* MeV/mm */
  G4double z=nloss; /* MeV/mm */

  /*  
  y/=(MeV/um);
  z/=(MeV/um);
  */
  /*  
  y/=materialDensity;
  z/=materialDensity;
  */
  /*
  y/=(MeV*cm2/g);
  z/=(MeV*cm2/g);
  */
  /*    
  y/= TotNbOfAtomsPerVolume;      // return to per atom
  z/= TotNbOfAtomsPerVolume;
  y/= theZieglerFactor;           // return to Ziegler's units
  z/= theZieglerFactor;
  */
  /*
  y/= 8.462 * z1 * z2 * m1 / rm ; // return to Linhard's units
  z/= 8.462 * z1 * z2 * m1 / rm ;
  */
  /*****************************************************************/
  
  fprintf(fp,"%g\t%g\n",x,y); /* dE/dx em funcao da energia do iao incidente */
  fprintf(fp1,"%g\t%g\n",x,z); /* nuclear stopping power only */
  fprintf(fp2,"%g\n",eloss/keV); /* electronic energy loss only */
  
  /* para procurar a relacao linear entre beta e Se - Alessio */
  y/=(MeV/cm);
  fprintf(fp3,"%lf\t%lf\n",beta,y/pow(particleCharge,2));
  //  fprintf(fp3,"%lf\t%lf\n",std::sqrt(x),y);
  
  // to study the step size as a function of the ion energy
  fprintf(fp4,"%lf\t%lf\n",x/keV,step/nm);

  /* para tirar o efeito da carga efectiva q e' multiplicada por defeito */
  //  fprintf(fp5,"%lf\t%lf\n",kineticEnergy,dEdx/chargeSquare);
      
  /*
    if(tN2ndariesTot>0)
    for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; lp1<(*fSecondary).size(); lp1++)
    fprintf(fp6,"%lf\t%lf\n",kineticEnergy,(*fSecondary)[lp1]->GetKineticEnergy()/step);
  */
  
  fclose(fp);fclose(fp1);fclose(fp2);fclose(fp3);fclose(fp4);
  
  
  /************** print verbosity output table ****************/
  
  if( verboseLevel >= 1 ){
  
    // change the number of decimal places of G4cout 
    // it returns the previous value to change it back
    G4int prec = G4cout.precision(3);

    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << G4endl;    
      G4cout << std::setw( 5) << "#Step#"     << " "
	     << std::setw( 6) << "X"          << "    "
	     << std::setw( 6) << "Y"          << "    "  
	     << std::setw( 6) << "Z"          << "    "
	     << std::setw( 9) << "KineE"      << " "
	     << std::setw( 9) << "dEStep"     << " "  
	     << std::setw(10) << "StepLeng"     
	     << std::setw(10) << "TrakLeng" 
	     << std::setw(10) << "Volume"    << "  "
	     << std::setw(10) << "Process"   << G4endl;	          
    }

    G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
      /*   
	   << std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	   << std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	   << std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	   << std::setw(7) << fTrack->GetMomentumDirection().x()
	   << std::setw(7) << fTrack->GetMomentumDirection().y()
	   << std::setw(7) << fTrack->GetMomentumDirection().z()
	   << std::setw(7) << DeltaDirection.x()
	   << std::setw(7) << DeltaDirection.y()
	   << std::setw(7) << DeltaDirection.z()
	   << "\n    : "
      */
	   << std::setw(6) << G4BestUnit(kineticEnergy,"Energy")
	   << std::setw(6) << G4BestUnit(finalT,"Energy")
	   << std::setw(6) << G4BestUnit(ionloss,"Energy")
	   << std::setw(6) << G4BestUnit(nloss*step,"Energy")
	   << std::setw(7) << G4BestUnit(t2ndariesEnergy,"Energy")
      //	   << std::setw(4) << G4BestUnit(dloss,"Energy/Length")
	   << std::setw( 9) << beta
	   << std::setw( 7) << chargeSquare //<< " " << chargeSquareNorm
      //	   << std::setw(8) << Z2norm
	   << std::setw( 9) << G4BestUnit(step,"Length")
	   << std::setw( 6) << G4BestUnit(dEdx,"Energy/Length")
	   << std::setw( 6) << G4BestUnit(nloss,"Energy/Length")
	   << std::setw( 9) << G4BestUnit(fTrack->GetTrackLength(),"Length");

     
    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << std::setw(8) << fTrack->GetVolume()->GetName() << " of " << material->GetName();
    } else {
      G4cout << std::setw(10) << "OutOfWorld";
    }
    
    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << "  " 
	     << std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	->GetProcessName();
    } else {
      G4cout << "   UserLimit";
    }
    
    G4cout << G4endl;
      
    
    // secondaries

    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << std::setw(3) << tN2ndariesTot 
	       << "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
	       << " ---------------"
	       << G4endl
	       << "    :  "
	       << std::setw( 6) << "X"          << "    "
	       << std::setw( 6) << "Y"          << "    "  
	       << std::setw( 6) << "Z"          << "    "
	       << std::setw( 6) << "R"          << "    "
	       << std::setw( 7) << "dirx"       << " "
	       << std::setw( 7) << "diry"       << " "
	       << std::setw( 7) << "dirz"       << " "
	       << std::setw( 9) << "kineE"      << " "
	       << std::setw( 9) << "ParticleName"//      << " "
	       << G4endl;
	
	G4double SecondariesEnergy = 0.0;
	G4ThreeVector SecondariesMomentum;
	
	G4String beginWarning = "\t\t\t\t\t *********************** WARNING ********************************";
	G4String symmetricWarning = "\n\t\t\t\t\t * the projectile/target atom combination is NOT symmetrical !! *";
	G4String endWarning = "\n\t\t\t\t\t ****************************************************************\n";

	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; 
	           lp1<(*fSecondary).size(); lp1++){

	  // Verify symmetric projectile/target atom combination
	  G4String secondName = (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	  //	  G4cout << "Secondary Name = " << secondName << G4endl;
	  if(secondName != particleName || secondName(0) != materialName(0) || secondName(1) != materialName(1))
	    G4cout << beginWarning << symmetricWarning << endWarning;
	  
	  G4cout << "    : "
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().mag(),"Length")
		 << std::setw(8)
		 << (*fSecondary)[lp1]->GetMomentumDirection().x()
		 << std::setw(8)
		 << (*fSecondary)[lp1]->GetMomentumDirection().y()
		 << std::setw(8)
		 << (*fSecondary)[lp1]->GetMomentumDirection().z()
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << std::setw(12)
		 << secondName
		 << G4endl;
	  G4cout.precision(prec);
	  
	  SecondariesEnergy += (*fSecondary)[lp1]->GetKineticEnergy();
	  SecondariesMomentum += (*fSecondary)[lp1]->GetMomentum();
	}
	//	G4cout << "SecondariesEnergy = " << G4BestUnit(SecondariesEnergy,"Energy") << G4endl;

	G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << G4endl;


	// Verify energy and momentum conservation

	nloss*=step;
	G4double secEnergyError = std::fabs(nloss - SecondariesEnergy) / nloss;
	//	G4cout << "relative secEnergyError = " << secEnergyError << G4endl;
	G4double deltaEnergyError = std::fabs(eloss+nloss - ionloss) / ionloss;
	//	G4cout << "relative deltaEnergyError = " << deltaEnergyError << G4endl;
	G4double energyTolerance = 2.2e-14;
	G4double momentumTolerance = 0.11;

	G4bool energyConservation = secEnergyError < energyTolerance && deltaEnergyError < energyTolerance;
	G4bool momentumConservation = DeltaMomentum.isNear(SecondariesMomentum,momentumTolerance);

	G4String energyWarning = "\n\t\t\t\t\t **********\tNO ENERGY CONSERVATION !!!\t\t*********";
	energyWarning += "\n\t\t\t\t\t **********\tTRY DECREASING STEP SIZE\t\t*********";

	G4String momentumWarning = "******* WARNING ****** Momentum conservation is being violated ";
	momentumWarning+= "above defined tolerance, try decreasing STEP SIZE\n";

	G4cout.precision(3);
	if(fTrack->GetTrackID() == 1 && fTrack->GetCurrentStepNumber() == 1)
	  G4cout << "----> Nuclear stopping power model defined in SteppingVerbose is "
		 << theNuclearTable
		 << G4endl;
  	G4cout << "Energy conservation? "
	       << (energyConservation?"YES":"NO")
	       << ", by "  << deltaEnergyError
	       << G4endl;
	if(!energyConservation)
	  G4cout << beginWarning << energyWarning << endWarning;
	G4cout << "Momentum conservation? "
	       << (momentumConservation?"YES":"NO")
	       << ", within " 
	       << DeltaMomentum.howNear(SecondariesMomentum) << ", ("
	       << DeltaMomentum.x() << ", "
	       << DeltaMomentum.y() << ", "
	       << DeltaMomentum.z() << ") vs. ("
	       << SecondariesMomentum.x() << ", "
	       << SecondariesMomentum.y() << ", "
	       << SecondariesMomentum.z() << ")"
	       << G4endl;
	if(!momentumConservation) G4cout << momentumWarning;
      }
    }
    G4cout.precision(prec);
  }
  //  exit(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02SteppingVerbose::TrackingStarted()
{
  
  CopyState();

  G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){
    
    G4cout << "Particle initial position = ("
	   << G4BestUnit(fTrack->GetPosition().x(),"Length") << ", "
	   << G4BestUnit(fTrack->GetPosition().y(),"Length") << ", "
	   << G4BestUnit(fTrack->GetPosition().z(),"Length") << ")"
	   << "\tParticle initial direction = ("
	   << fTrack->GetMomentumDirection().x() << ", "
	   << fTrack->GetMomentumDirection().y() << ", "
	   << fTrack->GetMomentumDirection().z() << ")"
	   << G4endl;
    
    G4cout << std::setw( 5) << "Step#"      << " "
      /*  
	   << std::setw( 6) << "X"          << "    "
      	   << std::setw( 6) << "Y"          << "    "  
      	   << std::setw( 6) << "Z"          << "    "
      */
	   << std::setw( 9) << "kineE"      << " "
      	   << std::setw( 9) << "finalT"      << " "
	   << std::setw(10) << "dE Total"     << " "  
	   << std::setw(11) << "dE Nuclear"  << " "  
	   << std::setw( 8) << "T(Sec)"     
      //      	   << std::setw( 9) << "dloss"      << "  "
	   << std::setw( 7) << "beta"       << "  "
      //	   << std::setw(10) << "Zeff n i"
	   << std::setw( 7) << "Zeff^2"
	   << std::setw(11) << "StepLeng"  
	   << std::setw(13) << "dE/dx Elect"      << "  "
	   << std::setw(12) << "dE/dx Nuclear"    << "  "
	   << std::setw(10) << "TrakLeng"
	   << std::setw(12) << "Volume"     << "  "
	   << std::setw(14) << "Process"    << G4endl;	     
  }
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
