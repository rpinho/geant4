// $Id: G4Cerenkov.cc,v 1.21 2006/06/29 19:56:03 gunter Exp $
////////////////////////////////////////////////////////////////////////
// Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:
// Description: Continuous Process -- Generation of 
// Version:
// Created:
// Author:
// Updated:
// mail:
//
////////////////////////////////////////////////////////////////////////

#include "ExN02SteppingVerbose.hh"
#include "G4hLowEnergyIonisation.hh"
#include "RPNuclearRecoil01.hh"

using namespace std;

        /////////////////
        // Constructors
        /////////////////

RPNuclearRecoil01::RPNuclearRecoil01(const G4String& processName)
           : G4VContinuousProcess(processName)
{
	if (verboseLevel>0)
	  G4cout << GetProcessName() << " is created " << G4endl;
	
	MinKineticEnergy = 10.0*eV ;
}

        ////////////////
        // Destructors
        ////////////////

RPNuclearRecoil01::~RPNuclearRecoil01()
{
  if(theNuclearStoppingModel)delete theNuclearStoppingModel;
}

        ////////////
        // Methods
        ////////////

// AlongStepDoIt
// -------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* RPNuclearRecoil01::AlongStepDoIt(
                                           const G4Track& trackData,
                                           const G4Step& stepData)
{
  /* the nuclear recoils, based on the delta electons emission of G4hLowEnergyIonisation::PostStepDoIt */

  size_t totalNumber = 1;
  G4double KineticEnergy,TotalEnergy,TotalMomentum,ParticleMass,
           RecoilKineticEnergy,RecoilTotalMomentum,costheta,sintheta,phi,
           dirx,diry,dirz,finalKineticEnergy,finalPx,finalPy,finalPz,
           Psquare,finalMomentum ;

  aParticleChange.Initialize(trackData) ;

  const G4MaterialCutsCouple* couple = trackData.GetMaterialCutsCouple();
  const G4Material* material = couple->GetMaterial();

  // get the actual (true) Step length from stepData
  const G4double step = stepData.GetStepLength(); //printf("AlongStep GetStepLength = %g\n",step);

  const G4DynamicParticle* aParticle = trackData.GetDynamicParticle() ;
  G4ParticleDefinition* particle = aParticle->GetDefinition();
  KineticEnergy = aParticle->GetKineticEnergy(); //printf("AlongStep kineticEnergy = %g\n",kineticEnergy);

  /************** nuclear stopping power ****************/
  G4String theNuclearTable = "ICRU_R49";
  G4VLowEnergyModel* theNuclearStoppingModel = new G4hNuclearStoppingModel(theNuclearTable) ;
  G4double nloss = theNuclearStoppingModel->TheValue(aParticle, material)*step ;
  //  G4cout <<"AlongStep nloss = " << nloss << G4endl;

  ParticleMass=particle->GetPDGMass();
  TotalEnergy=KineticEnergy + ParticleMass ;
  Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ;
  TotalMomentum = std::sqrt(Psquare) ;
  G4ThreeVector ParticleDirection = aParticle->GetMomentumDirection() ;
  RecoilKineticEnergy = nloss;
  RecoilTotalMomentum = std::sqrt(RecoilKineticEnergy * (RecoilKineticEnergy +
						       2. * ParticleMass));
  costheta = RecoilKineticEnergy * (TotalEnergy + ParticleMass)
           /(RecoilTotalMomentum * TotalMomentum) ;
  G4double  cosTheta = 1-2*RecoilKineticEnergy/KineticEnergy;
  G4double  cosnu = std::sqrt(1-RecoilKineticEnergy/KineticEnergy);
  G4cout << costheta << " " << cosTheta << " " << cosnu << " " << std::sqrt((1+cosTheta)/2) << G4endl;

  //  protection against costheta > 1 or < -1   ---------------
  if ( costheta < -1. )
    costheta = -1. ;
  if ( costheta > +1. )
    costheta = +1. ;
  
  //  direction of the recoiled atom  ........
  phi = twopi * G4UniformRand() ;
  sintheta = std::sqrt(1. - costheta*costheta);
  dirx = sintheta * std::cos(phi) ;
  diry = sintheta * std::sin(phi) ;
  dirz = costheta ;
  
  G4ThreeVector RecoilDirection(dirx,diry,dirz) ;
  RecoilDirection.rotateUz(ParticleDirection) ;
  
  // create G4DynamicParticle object for the recoiled atom
  G4DynamicParticle *theRecoilAtom = new G4DynamicParticle;
  theRecoilAtom->SetKineticEnergy( RecoilKineticEnergy );
  theRecoilAtom->SetMomentumDirection(RecoilDirection.x(),
				    RecoilDirection.y(),
				    RecoilDirection.z());
  theRecoilAtom->SetDefinition(particle);
  
  // Save the recoiled atoms and update the ParticleDirection and Energy, i.e., the aParticleChange method */

  finalKineticEnergy = KineticEnergy - nloss;

  if (finalKineticEnergy > MinKineticEnergy) {
    if (nloss > MinKineticEnergy)
      {
	finalPx = TotalMomentum*ParticleDirection.x()
	  - RecoilTotalMomentum*RecoilDirection.x();
	finalPy = TotalMomentum*ParticleDirection.y()
	  - RecoilTotalMomentum*RecoilDirection.y();
	finalPz = TotalMomentum*ParticleDirection.z()
	  - RecoilTotalMomentum*RecoilDirection.z();
	finalMomentum =
	  std::sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz) ;
	finalPx /= finalMomentum ;
	finalPy /= finalMomentum ;
	finalPz /= finalMomentum ;
	
	const G4ThreeVector NewDirection(finalPx,finalPy,finalPz);
	G4double polarAngle = ParticleDirection.cosTheta(NewDirection);
	G4double escalar = ParticleDirection*NewDirection/(ParticleDirection.mag()*NewDirection.mag());
	G4cout << polarAngle << " " << escalar << G4endl;
	aParticleChange.ProposeMomentumDirection( finalPx,finalPy,finalPz );
	/* = aParticleChange.ProposeMomentumDirection( NewDirection ); */
      }
  }
  else
    {
      finalKineticEnergy = 0.;
      /*aParticleChange.ProposeMomentumDirection(ParticleDirection.x(),
                      ParticleDirection.y(),ParticleDirection.z());
      if(!aParticle->GetDefinition()->GetProcessManager()->
      GetAtRestProcessVector()->size())*/
        aParticleChange.ProposeTrackStatus(fStopAndKill);
	/*else
        aParticleChange.ProposeTrackStatus(fStopButAlive);*/
    }

  aParticleChange.ProposeEnergy( finalKineticEnergy );
  
  if(nloss > MinKineticEnergy) {
    aParticleChange.SetNumberOfSecondaries(totalNumber);
    aParticleChange.AddSecondary(theRecoilAtom);
  }

  return &aParticleChange ;
}
