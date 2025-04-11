/// \file PhysicsList.hh
/// \brief Implementation of the PhysicsList class (Mandatory)

#include "PhysicsList.hh"

// Include standard Electromagnetic
#include <G4EmStandardPhysics.hh>
// Include header for G4EmLivermorePhysics
#include <G4EmLivermorePhysics.hh>
//Include header for G4EmExtraPhysics
#include <G4EmExtraPhysics.hh>
// Set decay physics list
#include <G4DecayPhysics.hh>
// Needed to implement cuts
#include <G4ProductionCutsTable.hh>

#include  <G4EmLivermorePolarizedPhysics.hh>

#include <G4SystemOfUnits.hh>

#include "G4EmStandardPhysics_option2.hh"

#include "G4EmLowEPPhysics.hh"

#include <G4ParticleTable.hh>
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include <G4HadronPhysicsFTFP_BERT.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4PenelopePhotoElectricModel.hh>

//----------------------------------------

//#include "HadronElasticPhysicsHP.hh"
//#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonPhysics.hh"
#include "G4IonElasticPhysics.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "GammaNuclearPhysics.hh"
#include "ElectromagneticPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "NeutronHPphysics.hh"
#include "G4EmExtraPhysics.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  SetVerboseLevel(1);
  RegisterPhysics(new G4HadronElasticPhysicsHP() );
  RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP());
  RegisterPhysics(new G4IonPhysics(1));
  RegisterPhysics(new GammaNuclearPhysics("gamma"));
  RegisterPhysics(new G4EmStandardPhysics());
}
PhysicsList::~PhysicsList()
{;}


void PhysicsList::ConstructParticle()
{
  G4VModularPhysicsList::ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
  G4VModularPhysicsList::ConstructProcess();
}

void PhysicsList::SetCuts()
{
  //SetCutValue(1 * mm, "alpha");
}
