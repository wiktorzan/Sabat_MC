/// \file PhysicsList.hh
/// \brief Definition of the PhysicsList class (Mandatory)

#ifndef PHYSICS_LIST_H
#define PHYSICS_LIST_H 1

#include <G4VModularPhysicsList.hh>
#include "globals.hh"


/// Modular physics list


class PhysicsList : public G4VModularPhysicsList
{
public:
  ///constructor
  PhysicsList();

  /// destructor
  virtual ~PhysicsList();

  /// Builds particles
  void ConstructParticle() override;

  /// Build processes
  void ConstructProcess() override;

  /// Set user cuts

  void SetCuts() override;

};

#endif
