********* O OUTPUT ******************
========= Table of registered couples ==============================

Index : 0     used in the geometry : Yes     recalculation needed : No 
 Material : Air
 Range cuts        :  gamma 10000 fm     e- 10000 fm     e+ 10000 fm 
 Energy thresholds :  gamma 1 GeV    e- 1 GeV    e+ 1 GeV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 1     used in the geometry : Yes     recalculation needed : No 
 Material : Lead
 Range cuts        :  gamma 10000 fm     e- 10000 fm     e+ 10000 fm 
 Energy thresholds :  gamma 1 GeV    e- 1 GeV    e+ 1 GeV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

====================================================================

************ OS MEUS VALORES ***********
ExN02PhysicsList.cc
ExN02PhysicsList::ExN02PhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = .01*nm;
  SetVerboseLevel(1);
}

void ExN02PhysicsList::SetCuts()
{
  //G4VUserPhysicsList::SetCutsWithDefault method sets 
  //the default cut value for all particle types 

  SetCutsWithDefault();

  G4double lowlimit=1.*GeV; /* cut elevado para nao produzir electroes delta */
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100*GeV);
     
  if (verboseLevel>0) DumpCutValuesTable();
}

************ OS VALORES DEFAULT **********
G4VUserPhysicsList.cc
00068 G4VUserPhysicsList::G4VUserPhysicsList()
00069                    :verboseLevel(1),
(...)
00077 {
00078   // default cut value  (1.0mm)
00079   defaultCutValue = 1.0*mm;
00080 
00081   // pointer to the particle table
00082   theParticleTable = G4ParticleTable::GetParticleTable();
00083   theParticleIterator = theParticleTable->GetIterator();
00084 
00085   // pointer to the cuts table
00086   fCutsTable =  G4ProductionCutsTable::GetProductionCutsTable();
00087 
00088   // set energy range for SetCut calculation
00089   fCutsTable->SetEnergyRange(0.99*keV, 100*TeV);
(...)
00093 }

**************** O CÓDIGO QUE IMPRIME A INFORMAÇÃO NO OUTPUT ************
G4ProductionCutsTable.cc
00283 void G4ProductionCutsTable::DumpCouples() const
00284 {
00285   G4cout << G4endl;
00286   G4cout << "========= Table of registered couples =============================="
00287          << G4endl;
00288   for(CoupleTableIterator cItr=coupleTable.begin();
00289       cItr!=coupleTable.end();cItr++) {
00290     G4MaterialCutsCouple* aCouple = (*cItr);
00291     G4ProductionCuts* aCut = aCouple->GetProductionCuts();
00292     G4cout << G4endl;
00293     G4cout << "Index : " << aCouple->GetIndex() 
00294            << "     used in the geometry : ";
00295     if(aCouple->IsUsed()) G4cout << "Yes";
00296     else                  G4cout << "No ";
00297     G4cout << "     recalculation needed : ";
00298     if(aCouple->IsRecalcNeeded()) G4cout << "Yes";
00299     else                          G4cout << "No ";
00300     G4cout << G4endl;
00301     G4cout << " Material : " << aCouple->GetMaterial()->GetName() << G4endl;
00302     G4cout << " Range cuts        : " 
00303            << " gamma " << G4BestUnit(aCut->GetProductionCut("gamma"),"Length")
00304            << "    e- " << G4BestUnit(aCut->GetProductionCut("e-"),"Length")
00305            << "    e+ " << G4BestUnit(aCut->GetProductionCut("e+"),"Length")
00306            << G4endl;
00307     G4cout << " Energy thresholds : " ;
00308     if(aCouple->IsRecalcNeeded()) {
00309       G4cout << " is not ready to print";
00310     } else {
00311       G4cout << " gamma " << G4BestUnit((*(energyCutTable[0]))[aCouple->GetIndex()],"Energy")
00312              << "    e- " << G4BestUnit((*(energyCutTable[1]))[aCouple->GetIndex()],"Energy")
00313              << "    e+ " << G4BestUnit((*(energyCutTable[2]))[aCouple->GetIndex()],"Energy");
00314     }
00315     G4cout << G4endl;
00316 
00317     if(aCouple->IsUsed()) {
00318       G4cout << " Region(s) which use this couple : " << G4endl;
00319       typedef std::vector<G4Region*>::iterator regionIterator;
00320       for(regionIterator rItr=fG4RegionStore->begin();
00321           rItr!=fG4RegionStore->end();rItr++) {
00322         if (IsCoupleUsedInTheRegion(aCouple, *rItr) ){
00323           G4cout << "    " << (*rItr)->GetName() << G4endl;
00324         }
00325       }
00326     }
00327   }
00328   G4cout << G4endl;
00329   G4cout << "====================================================================" << G4endl;
00330   G4cout << G4endl;
00331 }
