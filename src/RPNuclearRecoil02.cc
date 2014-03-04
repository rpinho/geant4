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

//#include "ExN02SteppingVerbose.hh"
//#include "G4hLowEnergyIonisation.hh"
#include "RPNuclearRecoil02.hh"
#include "G4UnitsTable.hh"

using namespace std;

/////////////////
// Constructors
/////////////////

RPNuclearRecoil02::RPNuclearRecoil02(const G4String& processName)
  : G4VContinuousProcess(processName),
    // default values, but the user has at his disposal public methods to change them
    theNuclearTable("ICRU_R49"),  // public method SetNuclearStoppingPowerModel
    MinKineticEnergy(10.0*eV),    // public method SetMinKineticEnergy
    NumRecoilsProposed(1),        // public method ProposeNumRecoilsPerStep
    lossFlucFlag(false),          // public methods SetNuclearStoppingFluctuationsOn() and Off()
    stepLimitFlag(true),          // public methods SetStepLimitOn() and Off()
    fNumRecoilsFlag(false),       // public methods ForceNumRecoilsPerStepOn() and Off()
    fTrackSecondariesFirst(false) // public methods SetTrackSecondariesFirstOn() and Off()
{
  // Factor to convert PDG mass unit into AMU
  factorPDG2AMU = 1.007276/proton_mass_c2 ;
  // controle flag for output message
  verboseLevel = 0;
  
  if (verboseLevel>0)
    G4cout << GetProcessName() << " is created " << G4endl;
}

////////////////
// Destructors
////////////////

RPNuclearRecoil02::~RPNuclearRecoil02()
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

G4VParticleChange* RPNuclearRecoil02::AlongStepDoIt(
                                           const G4Track& trackData,
                                           const G4Step& stepData)
{
  // The "continuous" production of nuclear recoils, based on: 
  // -- the delta electrons emission and fluorescence generation of 
  // G4hLowEnergyIonisation::PostStepDoIt and the deexcitation of
  // G4hLowEnergyIonisation::DeexciteAtom in AlongStepDoIt,
  // -- the continous photon generation in G4Cerenkov.cc

  aParticleChange.Initialize(trackData) ;

  // get the material
  const G4MaterialCutsCouple* couple = trackData.GetMaterialCutsCouple();
  const G4Material* material = couple->GetMaterial();

  // get the actual (true) Step length from stepData
  const G4double step = stepData.GetStepLength();
  //  G4cout << "StepLength = " << G4BestUnit(step,"Length") << G4endl;

  // get the incident particle
  const G4DynamicParticle* aParticle = trackData.GetDynamicParticle() ;
  G4ParticleDefinition* particle = aParticle->GetDefinition();
  G4double KineticEnergy = aParticle->GetKineticEnergy();
  //  G4cout << "kineticEnergy = " << KineticEnergy << G4endl;
  G4ThreeVector ParticleDirection = aParticle->GetMomentumDirection() ;

  // The nuclear stopping power parametrization
  theNuclearStoppingModel = new G4hNuclearStoppingModel(theNuclearTable) ;
  G4double nloss = theNuclearStoppingModel->TheValue(aParticle, material)*step ;
  //  G4cout << "nloss = " << nloss << G4endl;
  
  // protection against negative values
  if(nloss < 0.0) nloss = 0.0;

  // The final kinetic energy of the incident particle after the step.
  // Don't forget to set G4hLowEnergyIonisation->SetNuclearStoppingOff()
  // in your Physics List, otherwise nloss will be taken twice
  G4double finalKineticEnergy = KineticEnergy - nloss;

  ////////////////////////////////////////////////////////////////
  // Generation of the recoil atoms above the CUT energy
  if(nloss > MinKineticEnergy) {
    
    // get number of recoils, still per unit length
    G4double MeanNumberOfRecoils = GetAverageNumberOfRecoilsPerUnitLength(material);
    //    G4cout << "MeanNumberOfRecoils = " << MeanNumberOfRecoils << G4endl;

    // protection against negative number of recoils
    if (MeanNumberOfRecoils < 0.0) MeanNumberOfRecoils = 0.0;
    //    G4cout << "MeanNumberOfRecoils = " << MeanNumberOfRecoils << G4endl;

    // the number of recoils in this step
    MeanNumberOfRecoils *= step;
    //    G4cout << "MeanNumberOfRecoils = " << MeanNumberOfRecoils << G4endl;
    
    // to force the number of recoils to be independent of step length
    if(fNumRecoilsFlag) MeanNumberOfRecoils = NumRecoilsProposed;
    //    G4cout << "MeanNumberOfRecoils = " << MeanNumberOfRecoils << G4endl;

    // Straggling
    if(lossFlucFlag && NumRecoilsProposed > 1) {
      G4double sig = GetStragglingFunctionWidth(aParticle,material);
      MeanNumberOfRecoils *= G4RandGauss::shoot(1.0,sig) ;
    }
    //    G4cout << "MeanNumberOfRecoils = " << MeanNumberOfRecoils << G4endl;

    // the number of recoils is, of course, an integer
    G4int NumRecoils = (G4int) rint(MeanNumberOfRecoils);
    //    G4cout << "NumRecoils = " << NumRecoils << G4endl;
    
    // protection against negative number of recoils
    if (NumRecoils < 0) NumRecoils = 0;
    //    G4cout << "NumRecoils = " << NumRecoils << G4endl;

    // a copy of the incident particle for (virtual) momentum updating
    G4DynamicParticle* aParticleCopy = new G4DynamicParticle(*aParticle);

    // get the vector of secondaries
    std::vector<G4DynamicParticle*>* newpart = 0;
    newpart = GenerateRecoils(aParticleCopy, material, nloss, NumRecoils);
    
    // if any secondary was sampled above the CUT energy
    if(newpart != 0) {
    
      // a secondary
      G4DynamicParticle* part = 0;
      size_t nSecondaries = newpart->size();
      aParticleChange.SetNumberOfSecondaries(nSecondaries);
      G4Track* newtrack = 0;
      const G4StepPoint* preStep = stepData.GetPreStepPoint();
      const G4StepPoint* postStep = stepData.GetPostStepPoint();
      G4ThreeVector r0 = preStep->GetPosition();
      G4ThreeVector deltaR = postStep->GetPosition();
      deltaR -= r0; // == stepData.GetDeltaPosition();
      G4double t0 = preStep->GetGlobalTime();
      G4double deltaT = postStep->GetGlobalTime();
      deltaT -= t0; // == stepData.GetDeltaTime();
      G4double time, rand, energy;
      G4ThreeVector position;
      
      for(size_t i=0; i<nSecondaries; i++) {
	
	part = (*newpart)[i];
	if(part) {
	  
	  // protection for energy conservation
	  energy = part->GetKineticEnergy();
	  if(energy <= nloss) {
	    
	    // Generate new G4Track object for the secondary,
	    // distributed evenly along the particle's step segment
	    rand = G4UniformRand();
	    time = t0 + deltaT*rand;
	    position  = r0 + deltaR*rand;
	    newtrack = new G4Track(part, time, position);
	    aParticleChange.AddSecondary(newtrack);
	  }
	  else {
	    G4cout << G4endl << "************************************************"
		   << G4endl << "ERROR!!!!! recoil energy > nloss"
		   << "=> VIOLATION OF ENERGY CONSERVATION!!!!"
		   << G4endl << "************************************************"
		   << G4endl;
	    G4cout << G4endl;
	    delete part;
	  }
	}
      }
      delete newpart;

      // Now really update the particle direction via the aParticleChange method
      if (finalKineticEnergy > MinKineticEnergy) {
	G4ThreeVector finalParticleMomentum = aParticleCopy->GetMomentumDirection();
	aParticleChange.ProposeMomentumDirection(finalParticleMomentum);
      }
    }
    delete aParticleCopy;
  }
  
  //////////////////////////////////////////////////////////////////
  // The finale, independent of the generation or not of secondaries 
  // Keep in mind that the primary particle ALWAYS looses nloss!

  // stop particle if the kinetic energy <= MinKineticEnergy (CUT)
  if (finalKineticEnergy <= MinKineticEnergy) {
    finalKineticEnergy = 0.0;
    aParticleChange.ProposeMomentumDirection(ParticleDirection);
    if(!particle->GetProcessManager()->
       GetAtRestProcessVector()->size())
      aParticleChange.ProposeTrackStatus(fStopAndKill);
    // if there are still AtRest (decay) processes to run keep it alive
    else
      aParticleChange.ProposeTrackStatus(fStopButAlive);
  }

  // update the particle energy via the aParticleChange method */
  aParticleChange.ProposeEnergy( finalKineticEnergy );
  
  return &aParticleChange ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4DynamicParticle*>*
RPNuclearRecoil02::GenerateRecoils(G4DynamicParticle* aParticle,
				   const G4Material* material,
				   G4double nloss,
				   G4int NumRecoils)
{

  if (verboseLevel > 0) {
	G4cout << "GenerateRecoils: cutForRecoils(eV)= " << MinKineticEnergy/eV
               << "  nLoss(keV)= " << nloss/keV
               << G4endl;
  }

  if(nloss <= MinKineticEnergy) return 0;

  // create vector of tracks of secondary particles

  std::vector<G4DynamicParticle*>* partVector =
    new std::vector<G4DynamicParticle*>;
  G4DynamicParticle* aSecondary = 0;

  // sample secondaries

  G4double RecoilKineticEnergy;
  G4double etot = 0.0;

  for (G4int i=0; i<NumRecoils; i++) {

    if(i!=NumRecoils-1) {
      
      // the partitioning of the stopping power per recoil
      RecoilKineticEnergy = nloss / NumRecoils;
      
      // Straggling
      if(lossFlucFlag) {
	G4double sig = GetStragglingFunctionWidth(aParticle,material);
	RecoilKineticEnergy *= G4RandGauss::shoot(1.0,sig) ;
      }
    }
    
    // the energy left for the last (or only) one is always 
    // imposed by energy conservation
    else
      RecoilKineticEnergy = nloss - etot;

    etot += RecoilKineticEnergy;

    // generate secondary only above CUT
    if(RecoilKineticEnergy > MinKineticEnergy) {
      aSecondary = GenerateOneRecoilAtom(aParticle,RecoilKineticEnergy);
      partVector->push_back(aSecondary);
    }
  }

  if(partVector->empty()) {
    if(verboseLevel > 0)
      G4cout << G4endl << "No secondaries were generated in this step" << G4endl;
    delete partVector;
    return 0;
  }
  
  // protection for energy conservation
  if(nloss != etot) {
    G4cout << G4endl << "************************************************"
	   << G4endl << "ERROR!!!!! nloss != etot"
	   << "=> VIOLATION OF ENERGY CONSERVATION!!!!"
	   << G4endl << "************************************************"
	   << G4endl;
    delete partVector;
    return 0;
  }
  
  return partVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticle* 
RPNuclearRecoil02::GenerateOneRecoilAtom(G4DynamicParticle* aParticle,
					 double RecoilKineticEnergy)
{
  G4double KineticEnergy,ParticleMass,TotalEnergy,Psquare,
           TotalMomentum,RecoilTotalMomentum,costheta,phi,
           sintheta,dirx,diry,dirz;

  // get the (relativistic) kinematics of the incoming particle
  G4ParticleDefinition* particle = aParticle->GetDefinition();
  KineticEnergy = aParticle->GetKineticEnergy();
  //  G4cout << "kineticEnergy = " << KineticEnergy << G4endl;
  ParticleMass=particle->GetPDGMass();
  TotalEnergy=KineticEnergy + ParticleMass ; // == aParticle->GetTotalEnergy();
  Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ;
  TotalMomentum = std::sqrt(Psquare) ; // == aParticle->GetTotalMomentum();
  G4ThreeVector ParticleDirection = aParticle->GetMomentumDirection() ;

  RecoilTotalMomentum = std::sqrt(RecoilKineticEnergy * (RecoilKineticEnergy +
							 2. * ParticleMass));
  costheta = RecoilKineticEnergy * (TotalEnergy + ParticleMass)
           /(RecoilTotalMomentum * TotalMomentum) ;
  
  //  protection against costheta > 1 or < -1   ---------------
  if ( costheta < -1. )
    costheta = -1. ;
  if ( costheta > +1. )
    costheta = +1. ;
  
  //  direction of the recoiled atom, uniformly azimuth  ........
  phi = twopi * G4UniformRand() ;
  sintheta = std::sqrt(1. - costheta*costheta);
  dirx = sintheta * std::cos(phi) ;
  diry = sintheta * std::sin(phi) ;
  dirz = costheta ;
  
  // Create recoil momentum direction vector. The momentum 
  // direction is still with respect to the coordinate system 
  // where the primary particle direction is aligned with the z axis  
  G4ThreeVector RecoilDirection(dirx,diry,dirz) ;
  
  // Rotate momentum direction back to global reference system 
  RecoilDirection.rotateUz(ParticleDirection) ;
  
  // create G4DynamicParticle object for the recoiled atom
  G4DynamicParticle *theRecoilAtom = new G4DynamicParticle;
  theRecoilAtom->SetKineticEnergy(RecoilKineticEnergy);
  theRecoilAtom->SetMomentumDirection(RecoilDirection);
  // Don't forget this process is, at the time being, only intended 
  // to symmetric projectile/target atom combinations, i.e., the 
  // incoming particle is identical with the atoms of the substance
  theRecoilAtom->SetDefinition(particle);

  // New primary particle direction by momentum conservation
  G4ThreeVector OldParticleMomentum = TotalMomentum*ParticleDirection; // == aParticle->GetMomentum();
  G4ThreeVector RecoilMomentum = RecoilTotalMomentum*RecoilDirection; // == theRecoilAtom->GetMomentum();
  G4ThreeVector NewParticleMomentum = OldParticleMomentum - RecoilMomentum;

  // Verify energy and momentum conservation by scattering angle
  G4double deltaDirAngle = OldParticleMomentum.cosTheta(NewParticleMomentum);
  G4double labScattAngle = std::sqrt(1-RecoilKineticEnergy/KineticEnergy);
  //  printf("deltaDirAngle=%.10lf\nlabScattAngle=%.10lf\n",deltaDirAngle,labScattAngle);
  // protection from conservation violation
  if( std::fabs(deltaDirAngle - labScattAngle) > 0.000001 ) {
    G4cout << G4endl << "************************************************"
	   << G4endl << "ERROR!!!!! deltaDirAngle != labScattAngle "
	   << "=> VIOLATION OF CONSERVATION LAWS!!!!"
	   << G4endl << "************************************************"
	   << G4endl;
    delete theRecoilAtom;
    return 0;
  }

  // (Virtually) update the incident particle energy and direction ( == momentum)
  aParticle->SetMomentum(NewParticleMomentum);
  
  return theRecoilAtom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double RPNuclearRecoil02::GetStragglingFunctionWidth(const G4DynamicParticle* aParticle,
						       const G4Material* material)
{
  // Projectile nucleus
  G4double KineticEnergy = aParticle->GetKineticEnergy();
  G4ParticleDefinition* particle = aParticle->GetDefinition();
  G4double z1 = particle->GetPDGCharge();
  G4double m1 = particle->GetPDGMass()*factorPDG2AMU;
  // Target nucleus
  G4double z2 = material->GetZ();
  G4double m2 = material->GetA()*mole/g;
  //  G4cout << G4endl << z1 << " " << m1 << " " << z2 << " " << m2 << G4endl;
  /* Linhard units */
  G4double energy = KineticEnergy/keV ;  // energy in keV
  G4double rm = (m1 + m2) * ( std::pow(z1, .23) + std::pow(z2, .23) ) ; // reduced mass
  //  G4cout << G4endl << "reduced mass = " << rm << G4endl;
  G4double er = 32.536 * m2 * energy / ( z1 * z2 * rm ) ;  // reduced energy
  //  G4cout << G4endl << "reduced energy = " << er << G4endl;
  G4double sig = 4.0 * m1 * m2 / ((m1 + m2)*(m1 + m2)*
	        (4.0 + 0.197*std::pow(er,-1.6991)+6.584*std::pow(er,-1.0494))) ;
  //  G4cout << G4endl << "sigma = " << sig << G4endl;
  //  G4cout << "relative width of the straggling function = " << sig/er*100 << "%" << G4endl;

  return sig;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double RPNuclearRecoil02::GetContinuousStepLimit(const G4Track& aTrack,
                                                   G4double  ,
                                                   G4double  ,
                                                   G4double& )
{
  if (!stepLimitFlag) return DBL_MAX;

  // If user has defined an average maximum number of recoils to
  // be generated in a Step, then return the Step length for that 
  // number of recoils. 
  
  if (NumRecoilsProposed <= 0) return DBL_MAX;
  
  const G4Material* material = aTrack.GetMaterial();
  
  G4double MeanNumberOfRecoilsPerUnitLength = GetAverageNumberOfRecoilsPerUnitLength(material);

  if(MeanNumberOfRecoilsPerUnitLength <= 0.0) return DBL_MAX;

  G4double StepLimit = NumRecoilsProposed / MeanNumberOfRecoilsPerUnitLength;

  return StepLimit;
}

// GetAverageNumberOfRecoilsPerUnitLength
// -------------------------
// This routine computes the number of nuclear recoils produced per
// GEANT-unit (millimeter) in the current medium. 
//             ^^^^^^^^^^

G4double RPNuclearRecoil02::GetAverageNumberOfRecoilsPerUnitLength(const G4Material *material) const
{
  /* na falta de melhor, de algo mais fisico para me dar o numero medio de colisoes por step nesta 
     continuous slowing down approximation, uso a distancia interatomica como a distancia media 
     entre diferentes colisoes, definindo portanto, e em media, um recuo nuclear por cada 
     distancia interatomica percorrida.
     desta forma, o numero medio de recoils e' dado por unidade de comprimento, para depois se
     multiplicar pelo devido step length (o menor dos propostos pelos processos fisicos concorrentes).
     nos Linhard nao consigo encontrar nada fisico q m limite o num. de colisoes. vou continuar 'a procura.
     no caso de um processo fisico com as seccoes eficazes diferenciais, ai' sim ja posso definir um step
     length com significado, com uma colisao (um integral) por cada step. ira' ser o proximo a implementar.
  */

  /* n number of atoms per unit volume of the material. Xe = 13.30 nm^-3 */
  G4double TotNbOfAtomsPerVolume = material->GetTotNbOfAtomsPerVolume();
  //  G4cout << "TotNbOfAtomsPerVolume = " << TotNbOfAtomsPerVolume/(1/(nm*nm*nm)) << " nm^-3" << G4endl;

  /* interatomic distance assuming spherical volume around each atom, from n*V = 1. Xe = 5.45 Aº */
  G4double InteratomicDistance = pow((3/(4*pi*TotNbOfAtomsPerVolume)),(1./3));
  //  G4cout << "InteratomicDistance = " << G4BestUnit(InteratomicDistance,"Length") << G4endl;

  G4double NumberOfRecoilsPerUnitLength = 1./(10*InteratomicDistance);
  //  G4cout << "NumberOfRecoilsPerUnitLength = " << NumberOfRecoilsPerUnitLength << G4endl;

  return NumberOfRecoilsPerUnitLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RPNuclearRecoil02::PrintInfoDefinition() const
{
  G4String comments = "The continuous production of nuclear recoils.";
  comments += "\n\t Recoil energy parametrised by Nuclear Stopping Power.";
  comments += "\n\t Only intended for symmetric projectile/target atom combinations:";
  comments += "\n\t the incoming particle is identical with the atoms of the target.";
 
  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n\t MinKineticEnergy is " << MinKineticEnergy / eV << " eV "
	 << "\n\t Nuclear stopping power model is " << theNuclearTable
	 << "\n\t Straggling of recoil energy and number is " << (lossFlucFlag?"On":"Off")
	 << "\n\t Step limiting is " << (stepLimitFlag?"On":"Off")
	 << "\n\t The number of recoils proposed per step is " << NumRecoilsProposed
	 << "\n\t Forcing NumRecoilsProposed recoils per step is " << (fNumRecoilsFlag?"On":"Off")
	 << "\n\t Tracking secondaries first is " << (fTrackSecondariesFirst?"On":"Off")
	 << G4endl << G4endl << G4endl ;
}
