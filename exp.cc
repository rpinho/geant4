  // some kinematics

  ParticleMass=aParticle->GetDefinition()->GetPDGMass(); /* em unidades de energia */
  KineticEnergy=aParticle->GetKineticEnergy();
  TotalEnergy=KineticEnergy + ParticleMass ;
  Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ; /* unidades de c=1 e tendo em atencao que T=E-mc^2 => p^2c^2=E^2-m^2c^4 */
  Esquare=TotalEnergy*TotalEnergy;
  betasquare=Psquare/Esquare; /* unidades de c=1 e tendo em atencao que E=m*gamma => E^2=m^2+p^2 => P^2/E^2=1-1/gamma^2 */
  G4ThreeVector ParticleDirection = aParticle->GetMomentumDirection() ;

  DeltaKineticEnergy = nloss;//x * tmax;

  DeltaTotalMomentum = std::sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
						       2. * ParticleMass));//electron_mass_c2 )) ;
  TotalMomentum = std::sqrt(Psquare) ;
  costheta = DeltaKineticEnergy * (TotalEnergy + ParticleMass/*electron_mass_c2*/)
    /(DeltaTotalMomentum * TotalMomentum) ;
