
#include "GammaNuclearPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Processes

#include "G4HadronInelasticProcess.hh"
#include "G4CascadeInterface.hh"

#include "G4SystemOfUnits.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4BetheHeitler5DModel.hh"

#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"

GammaNuclearPhysics::GammaNuclearPhysics(const G4String& name) : G4VPhysicsConstructor(name)
{}

GammaNuclearPhysics::~GammaNuclearPhysics()
{}

void GammaNuclearPhysics::ConstructProcess()
{
   G4ProcessManager* pManager = G4Gamma::Gamma()->GetProcessManager();
   //
  // G4HadronInelasticProcess* process = new G4HadronInelasticProcess("PhotonInelastic", G4Gamma::Gamma());
   //
  // G4CascadeInterface* bertini = new G4CascadeInterface();
  // bertini->SetMaxEnergy(10*GeV);
  // process->RegisterMe(bertini);
   //
  // pManager->AddDiscreteProcess(process);

// Temptorary to fix the compilation problem with missing cross-sections
   G4RayleighScattering* theRayleigh = new G4RayleighScattering();
   pManager->AddDiscreteProcess(theRayleigh);

   G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
   thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
   pManager->AddDiscreteProcess(thePhotoElectricEffect);

   G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
   theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
   pManager->AddDiscreteProcess(theComptonScattering);

   G4GammaConversion* theGammaConversion = new G4GammaConversion();
   theGammaConversion->SetEmModel(new G4BetheHeitler5DModel());
   pManager->AddDiscreteProcess(theGammaConversion);
}
