#include "SensitiveSD.hh"

#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include "SensitiveHit.hh"

#include "G4String.hh"

#include "G4VProcess.hh"   


SensitiveSD::SensitiveSD(G4String name) : G4VSensitiveDetector(name)
{
  collectionName.insert("energy_time");
}

G4bool SensitiveSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /*ROhist*/)
{
  SensitiveHit* hit = new SensitiveHit();
  G4StepPoint * postStepPoint = aStep->GetPostStepPoint();
  G4Track * theTrack = aStep->GetTrack();
  G4TouchableHandle touch1 = postStepPoint->GetTouchableHandle();
  G4int copynumber = touch1->GetCopyNumber();
  G4int ParID = aStep->GetTrack()->GetParentID();
  G4int StepID = theTrack->GetCurrentStepNumber();
  G4double eDep = aStep->GetTotalEnergyDeposit();
  G4double eKin = theTrack->GetKineticEnergy();
  G4double time = aStep->GetPostStepPoint()->GetGlobalTime();
  G4double timeL = aStep->GetPostStepPoint()->GetLocalTime();
  G4ThreeVector position = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector position2 = aStep->GetPostStepPoint()->GetPosition(); // New
  G4String name = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
  G4String Process = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  G4TouchableHandle touchable = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4int CopyNo = touchable->GetVolume()->GetCopyNo();
  G4String IntVol = touchable->GetVolume()->GetName();
  G4int TraID = aStep->GetTrack()->GetTrackID();

  hit->SetTime(time);
  hit->SetTimeL(timeL);
  hit->SetDeltaEnergy(eDep);
  hit->SetKinEnergy(eKin);
  hit->SetPosition(position);
  hit->SetPosition2(position2);
  hit->SetNbCopy(copynumber);
  hit->SetParID(ParID);
  hit->SetTrackID(TraID);
  hit->SetStepID(StepID);
  hit->SetParName(name);
  hit->SetPrcName(Process);
  hit->SetVolName(CopyNo);
  hit->SetVolName2(IntVol);

  fHitsCollection->insert(hit);

  return true;
}

void SensitiveSD::Initialize(G4HCofThisEvent* hcof)
{
  fHitsCollection = new SensitiveHitsCollection(SensitiveDetectorName, collectionName[0]);

  if (fHitsCollectionId < 0) {
    fHitsCollectionId = G4SDManager::GetSDMpointer()->GetCollectionID(GetName() + "/" + collectionName[0]);
  }

  hcof->AddHitsCollection(fHitsCollectionId, fHitsCollection);
}
