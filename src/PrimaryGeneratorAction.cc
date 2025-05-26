/// \file PrimaryGeneratorAction.hh
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleDefinition.hh"
#include "G4RandomDirection.hh"
#include "G4ParticleTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4Neutron.hh"
// to fill neutron and alpha direction
#include "Analysis.hh"
using namespace CLHEP;


PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction()
{
  fPrimGenMess = new PrimaryGeneratorMessenger(this);
  fGun = new G4ParticleGun();

  fSourcePosition =  G4ThreeVector(0., -55 * cm, 20. * cm);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ThreeVector dirAlp;
  G4double cosTheta = 1. - 2*G4UniformRand(), Phi = CLHEP::twopi*G4UniformRand();
  G4double sinTheta = sqrt(1. - cosTheta*cosTheta);
  G4ThreeVector dirNeu(sinTheta * cos(Phi), sinTheta * sin(Phi), cosTheta);
  G4ThreeVector dirAlpha(-1*sinTheta * cos(Phi), -1*sinTheta * sin(Phi), -1*cosTheta);

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();

  if (fShootNeutron) {
    fGun->SetParticleDefinition(G4Neutron::Definition());
    fGun->SetParticleEnergy(fNeutronEnergy);
    fGun->SetParticlePosition(fSourcePosition);
    fGun->SetParticleMomentumDirection(dirNeu);

    analysis->FillNtupleDColumn(17, dirNeu.theta()*(180/CLHEP::pi));
    analysis->FillNtupleDColumn(18, dirNeu.phi()*(180/CLHEP::pi));
    fGun->GeneratePrimaryVertex(anEvent);
  }

  if (fShootAlpha) {
    G4ParticleDefinition* myParticle = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    fGun->SetParticleDefinition(myParticle);
    fGun->SetParticleEnergy(fAlphaEnergy);
    fGun->SetParticlePosition(fSourcePosition);
    fGun->SetParticleMomentumDirection(dirAlpha);

    analysis->FillNtupleDColumn(22, dirAlpha.theta()*180/CLHEP::pi);
    analysis->FillNtupleDColumn(23, dirAlpha.phi()*180/CLHEP::pi);
    fGun->GeneratePrimaryVertex(anEvent);
  }
}
