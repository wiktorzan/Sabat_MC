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

//++++++++++++++++++++++

PhysicsList::PhysicsList() :
    G4VModularPhysicsList()
{

    SetVerboseLevel(1);
 // Hadron Elastic scattering
    RegisterPhysics( new G4HadronElasticPhysicsHP() );
 // RegisterPhysics( new NeutronHPphysics());
 
 // Hadron Inelastic Physics
   // RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(1));   
    RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP());
    // RegisterPhysics( new G4HadronInelasticQBBC(1));
    ////RegisterPhysics( new G4HadronPhysicsINCLXX(verb));
    // RegisterPhysics( new G4HadronPhysicsQGSP_BIC());
  
  
// Ion Physics
     RegisterPhysics( new G4IonPhysics(1));
 // RegisterPhysics( new G4IonINCLXXPhysics(verb));
   
 // Gamma-Nuclear Physics
    RegisterPhysics( new GammaNuclearPhysics("gamma"));
    
// Em interactions    
    RegisterPhysics(new G4EmStandardPhysics());
 
 // Decay
 // RegisterPhysics(new G4DecayPhysics());
 
 // Radioactive decay
 // RegisterPhysics(new G4RadioactiveDecayPhysics());
}
PhysicsList::~PhysicsList()
{    ;  }


void PhysicsList::ConstructParticle()
{
    // Call parent method. can be Replaced, if required
    G4VModularPhysicsList::ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
    // Call parent method. Replace it, if required
    G4VModularPhysicsList::ConstructProcess();
}

void PhysicsList::SetCuts()
{
    // sks - The method SetCuts() is mandatory in the interface. Here, one just use
    // the default SetCuts() provided by the base class.
    // Call parent method. Replace it, if required
   // G4VUserPhysicsList::SetCuts(1.*mm,"e-"); 

/* commented to check for missing tracks - adam stuff
     SetCutValue(50 * mm, "e+");
     SetCutValue(50 * mm, "La139");
     SetCutValue(50 * mm, "Br81");
     SetCutValue(50 * mm, "Br79"); 
     SetCutValue(50 * mm, "La138");
     SetCutValue(50 * mm, "Br80");
     SetCutValue(50 * mm, "Br78");
     SetCutValue(50 * mm, "Br82");
     SetCutValue(50 * mm, "La140");
     SetCutValue(50 * mm, "Se78");
     SetCutValue(50 * mm, "proton");
     SetCutValue(50 * mm, "Ce141");
     SetCutValue(50 * mm, "Ce142");
     SetCutValue(50 * mm, "Ce139");
     SetCutValue(50 * mm, "Ce140");
     SetCutValue(50 * mm, "Ba139");
     SetCutValue(50 * mm, "C13");
     SetCutValue(50 * mm, "deuteron");
  */
  //  SetCutValue(50 * mm, "e-");
  //  SetCutValue(10 * mm, "gamma");


}
