/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"


/// The primary generator action class

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  ///constructor
  PrimaryGeneratorAction();
  ///destructor
  ~PrimaryGeneratorAction();

  G4ParticleGun* GetParticleGun() {return fGun;}

  ///main interface
  void GeneratePrimaries(G4Event* anEvent) override;
private:
  G4ParticleGun* fGun;
};

#endif
