

#ifndef ElectromagneticPhysics_h
#define ElectromagneticPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class ElectromagneticPhysics : public G4VPhysicsConstructor
{
public:
  ElectromagneticPhysics(const G4String& name = "standard");
  ~ElectromagneticPhysics();

// This method is dummy for physics
  virtual void ConstructParticle() {};
  virtual void ConstructProcess();
};

#endif

