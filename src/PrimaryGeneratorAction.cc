/// \file PrimaryGeneratorAction.hh
/// \brief Implementation of the PrimaryGeneratorAction class

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
  fGun = new G4ParticleGun();
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

  for(int i=0; i<=1; i++) {
    if (i == 0) {
      fGun->SetParticleDefinition(G4Neutron::Definition());
      fGun->SetParticleEnergy(14.1*MeV);
      fGun->SetParticlePosition(G4ThreeVector(0,-15*cm,0));

      analysis->FillNtupleDColumn(17, dirNeu.theta()*(180/CLHEP::pi));
      analysis->FillNtupleDColumn(18, dirNeu.phi()*(180/CLHEP::pi));

      fGun->SetParticleMomentumDirection(dirNeu);
    } else if (i==1) {
      G4ParticleDefinition* myParticle = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
      fGun->SetParticleDefinition(myParticle);
      fGun->SetParticleEnergy(3.490*MeV);
      fGun->SetParticlePosition(G4ThreeVector(0,-15*cm,0));
      fGun->SetParticleMomentumDirection(dirAlpha);

      analysis->FillNtupleDColumn(19,dirAlpha.theta()*180/3.1414);
      analysis->FillNtupleDColumn(20,dirAlpha.phi()*180/3.1414);
    }
    fGun->GeneratePrimaryVertex(anEvent);
  }
}
