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
// ------------------------------------------------------------
//
// Class Description:
//
// Continuous Process -- Generation of nuclear recoils
// Class inherits publicly from G4VContinuousProcess.
//
// The "continuous" production of nuclear recoils, based on: 
// -- the delta electrons emission and fluorescence generation of 
// G4hLowEnergyIonisation::PostStepDoIt and the deexcitation of
// G4hLowEnergyIonisation::DeexciteAtom in AlongStepDoIt,
// -- the continous photon generation in G4Cerenkov.cc
//
// The user may select parametrisation tables for nuclear stopping powers
// The list of available tables:
// Nuclear stopping powers:    "ICRU_49" (default), "Ziegler1977",
//                             "Ziegler1985"
//
// Don't forget this process is, at the time being, only intended 
// to symmetric projectile/target atom combinations, i.e., the 
// incoming particle is identical with the atoms of the substance
//
// Class Description - End:

////////////////////////////////////////////////////////////////////////

#ifndef RPNuclearRecoil02_h
#define RPNuclearRecoil02_h 1

/////////////
// Includes
/////////////

#include "G4hLowEnergyIonisation.hh"

#include "globals.hh"
#include "Randomize.hh"

#include "G4VContinuousProcess.hh"

#include "G4ProcessManager.hh"

/////////////////////
// Class Definition
/////////////////////

class RPNuclearRecoil02 : public G4VContinuousProcess  
{

public:  // With description

  ////////////////////////////////
  // Constructor and Destructor
  ////////////////////////////////
  
  RPNuclearRecoil02(const G4String& processName = "contNucRcls"); 
  // Constructor: The nuclear recoils process for hadrons/ions 
  // to be include in the UserPhysicsList

  ~RPNuclearRecoil02();	
  // Destructor

  ///////////////////////
  // Mandatory methods
  ///////////////////////
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  // True for all charged hadrons/ions
  
  G4double GetContinuousStepLimit(const G4Track&,
				  G4double,
				  G4double,
				  G4double&); 
  // Returns the continuous step limit
  
  G4VParticleChange* AlongStepDoIt(const G4Track& trackData ,
				   const G4Step& stepData ) ;
  // Simulation of "continuous" nuclear recoils production.

  ///////////////////////////////////////////////////////////////
  // Public methods for the user to set parameters and control
  // this physical process behaviour via his Physics List
  ///////////////////////////////////////////////////////////////

  void SetMinKineticEnergy(G4double energy) 
              {MinKineticEnergy = energy;};
  // Set the threshold ion energy, THE CUT.
  // It is the same for both the primary and the secondary particles 
  // Default is 10.0*eV.
  
  void SetNuclearStoppingPowerModel(const G4String& dedxTable)
                               {theNuclearTable = dedxTable;};
  // This method defines the nuclear stopping parametrisation method 
  // via the name of the table. 
  // Default is "ICRU_49".

  void SetNuclearStoppingFluctuationsOn() {lossFlucFlag = true;}; 
  void SetNuclearStoppingFluctuationsOff(){lossFlucFlag = false;}; 
  // If set, each recoil energy and also the number of recoils per 
  // step (if NumRecoilsProposed > 1) are Gaussian sampled.
  // Default is Off.
  
  void SetStepLimitOn() {stepLimitFlag = true;};
  void SetStepLimitOff(){stepLimitFlag = false;};
  // If set, this process proposes a step length (other than DBL_MAX) 
  // corresponding to NumRecoilsProposed recoils per step.
  // Default is On.

  void ForceNumRecoilsPerStepOn() {fNumRecoilsFlag = true;};
  void ForceNumRecoilsPerStepOff(){fNumRecoilsFlag = false;};
  // If set, it will force NumRecoilsProposed recoils per step 
  // (or its corresponding gaussian sampled number, if straggling is set),
  // independent of step size. Otherwise step limits the NumRecoils.
  // Default is Off. 

  void SetTrackSecondariesFirstOn() {fTrackSecondariesFirst = true;};
  void SetTrackSecondariesFirstOff(){fTrackSecondariesFirst = false;};
  // If set, the primary particle tracking is interrupted and any 
  // produced nuclear recoils are tracked next. When all have 
  // been tracked, the tracking of the primary resumes. 
  // Default is Off.

  void ProposeNumRecoilsPerStep(const G4int NumRecoils){NumRecoilsProposed = NumRecoils;};
  // Propose a number of nuclear recoils to be generated during 
  // a step. Its actual use depends on two other flags: if StepLimit 
  // is set, the (proposed) step is limited by this number; if 
  // ForceNumRecoilsPerStep is set, it will then trully  impose this 
  // number of recoils per step (or its corresponding gaussian 
  // sampled number, if straggling is set and NumRecoils > 1).
  // Default is 1

  //////////////////////////////////
  // Dump out process information
  //////////////////////////////////

  void PrintInfoDefinition() const;
  // Print out the class parameters
  
  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
  {if(&aParticleType == G4Proton::Proton())PrintInfoDefinition();}
  // Messaged by the Particle definition (via the Process manager)
  // whenever cross section tables have to be rebuilt (i.e. if new
  // materials have been defined). 
  // It's only being used to automatically run PrintInfoDefinition()
  
protected:
  
private: // private data members
 
  std::vector<G4DynamicParticle*>* GenerateRecoils(G4DynamicParticle* aParticle,
						   const G4Material* material,
						   G4double nloss,
						   G4int NumRecoils);
  // samples recoil energy and creates vector of tracks of secondary particles

  G4DynamicParticle* GenerateOneRecoilAtom(G4DynamicParticle* aParticle,
					   double RecoilKineticEnergy);
  // samples recoil direction and creates object for the recoiled atom
  
  G4double GetAverageNumberOfRecoilsPerUnitLength(const G4Material *material) const;
  // tries to return something similar to the inverse of a MeanFreePath

  G4double GetStragglingFunctionWidth(const G4DynamicParticle* aParticle,
				      const G4Material* material);
  // taken from G4hICRU49Nuclear.cc or G4hZiegler1985Nuclear.cc

  // the parametrised model
  G4VLowEnergyModel* theNuclearStoppingModel; // the nuclear stopping model
  G4String theNuclearTable; // name of parametrisation table of nuclear stopping power
  
  // the parameters
  G4double MinKineticEnergy;
  G4int  NumRecoilsProposed;

  // Factor to convert PDG mass unit into AMU
  G4double factorPDG2AMU;
  
  // boolean flags
  G4bool lossFlucFlag;
  G4bool stepLimitFlag;
  G4bool fNumRecoilsFlag;
  G4bool fTrackSecondariesFirst;
};

////////////////////
// Inline methods
////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* from G4hLowEnergyIonisation.hh */

inline G4bool RPNuclearRecoil02::IsApplicable(const G4ParticleDefinition& particle)
{
   return(particle.GetPDGCharge() != 0.0
       && particle.GetPDGMass() > proton_mass_c2*0.1);

   /* a teoria de Linhard, i.e., o nuclear stopping power (nuclear quenching) so e' 
      relevante para particulas pesadas, mas nao esquecer que o protao ja e' 
      considerado uma particula pesada..
      as outras condicoes de aplicabilidade envolvem a velocidade (ou energia), uma 
      vez que so estamos interessadas em particulas lentas (ou pouco energeticas), 
      de modo a que o potencial de Coulomb seja screened. no entanto, aqui neste 
      metodo so temos acesso ao ParticleDefinition, que so nos da as propriedades 
      estaticas da particula, nao a velocidade ou energia.
      a condicao mais importante na teoria de Linhard e' v < v1 = Z1^(2/3)*v0.
      Xe: v1 = 14*v0 = 0.10*c => E < E1 = 616 MeV -> ISTO E' UM EXAGERO, TEM Q TAR MAL!!!
   */
}

#endif
