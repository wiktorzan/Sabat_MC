
#include "GammaNuclearPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Processes
#include "G4PhotoNuclearCrossSection.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4LowEGammaNuclearModel.hh"
#include "G4CascadeInterface.hh"
#include "G4SystemOfUnits.hh"

/*#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4BetheHeitler5DModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"*/

GammaNuclearPhysics::GammaNuclearPhysics(const G4String& name) : G4VPhysicsConstructor(name)
{}

GammaNuclearPhysics::~GammaNuclearPhysics()
{}

void GammaNuclearPhysics::ConstructProcess()
{
   G4HadronInelasticProcess* process = new G4HadronInelasticProcess("photonNuclear", G4Gamma::Definition());
   process->AddDataSet(new G4PhotoNuclearCrossSection);

   // to not register a model, set Emax=0; eg. Emax1 = 0.
   const G4double Emax1 = 200 * MeV, Emax2 = 10 * GeV;

   if (Emax1 > 0.) {  // model 1
      G4LowEGammaNuclearModel* model1 = new G4LowEGammaNuclearModel();
      model1->SetMaxEnergy(Emax1);
      process->RegisterMe(model1);
   }

   if (Emax2 > 0.) {  // model 2
      G4CascadeInterface* model2 = new G4CascadeInterface();
      G4double Emin2 = std::max(Emax1 - 1 * MeV, 0.);
      model2->SetMinEnergy(Emin2);
      model2->SetMaxEnergy(Emax2);
      process->RegisterMe(model2);
   }

 /*  G4ProcessManager* pManager = G4Gamma::Gamma()->GetProcessManager();
   pManager->AddDiscreteProcess(process);

// Temporary to fix the compilation problem with missing cross-sections
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
   pManager->AddDiscreteProcess(theGammaConversion);*/
}
