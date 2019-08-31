// Updated on 8/6/2014
// Updated on 9/16/2014, Hexc, Olesya: added CLHEP namespace for g4.10 version
// 
#include "globals.hh"
#include "MAGCOSPhysicsList.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4IonConstructor.hh"
#include "G4ios.hh"
#include "G4UserSpecialCuts.hh"
//#include "MAGCOSPhysicsMessenger.hh"
#include "G4Transportation.hh"
#include "G4PropagatorInField.hh"
//#include "G4TransportationManager.hh"

#include "MYTransportation.hh"
#include "G4ProcessTable.hh"

#ifndef BEFORE_V7
#include "G4StepLimiter.hh"
#endif

#include <iomanip>                

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////////////
//
MAGCOSPhysicsList::MAGCOSPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 3.*mm;
  SetVerboseLevel(1);
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSPhysicsList::~MAGCOSPhysicsList()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 
  
  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
  G4Alpha::AlphaDefinition();
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4He3::He3Definition();
  
  //  generic ion
  G4GenericIon::GenericIonDefinition();
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPhysicsList::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPhysicsList::ConstructBarions()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPhysicsList::ConstructProcess()
{ 
#ifndef BEFORE_V7
//StepLimiter process since version 7.0
  theParticleIterator->reset();
  G4StepLimiter* theStepLimiterProcess = new G4StepLimiter();
  while( (*theParticleIterator)() ){
    	G4ParticleDefinition* particle = theParticleIterator->value();
    	G4ProcessManager* pmanager = particle->GetProcessManager();
    	G4String particleName = particle->GetParticleName();
  	pmanager->AddDiscreteProcess(theStepLimiterProcess);
  }	

  
#endif
  AddTransportation(); 
  AddMyTransportation(); 
}

////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPhysicsList::AddMyTransportation()
{ 
  MYTransportation* theTransportationProcess= new MYTransportation();
  //theTransportationProcess->GetPropagatorInField()->SetMaxLoopCount(1000);
  
  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (!particle->IsShortLived()) {
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
	// Error !! no process manager
	G4cout << " !!!!!!!! G4VUserPhysicsList::AddTransportation : no process manager!" << G4endl;
	//G4Exception("G4VUserPhysicsList::AddTransportation : no process manager!", "Error", 1, "Error");
      } 
      else {
	// add transportation with ordering = ( -1, "first", "first" )
	pmanager->AddProcess(theTransportationProcess);
	pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
	pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
      }
    } 
    else {
      // shortlived particle case
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPhysicsList::SetCuts()
{
  if (verboseLevel >1){
    G4cout << "MAGCOSPhysicsList::SetCuts:";
  }  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutValue(1.*m,"gamma");
  SetCutValue(1.*m,"e-");
  SetCutValue(1.*m,"e+");
  
  SetCutValue(1.*m,"proton");
  SetCutValue(1.*m,"anti_proton");
  
  // Commented out the line for the time being:  8/6/2014
  //  SetCutValueForOthers(defaultCutValue );
}
