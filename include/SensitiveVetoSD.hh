#ifndef SENSITIVEVETOSD_HH
#define SENSITIVEVETOSD_HH

#include <G4VSensitiveDetector.hh>

#include "SensitiveHit.hh"

class SensitiveVetoSD : public G4VSensitiveDetector
{
public:
  SensitiveVetoSD(G4String name);

  void Initialize(G4HCofThisEvent*) override;

protected:
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) override;

private:
  SensitiveHitsCollection* fHitsCollection {nullptr};
  G4int fHitsCollectionId {-1};
};

#endif
