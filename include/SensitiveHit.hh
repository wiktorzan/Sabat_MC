#ifndef SENSITIVEHIT_HH
#define SENSITIVEHIT_HH

#include <G4VHit.hh>
#include <G4THitsMap.hh>
#include <G4ThreeVector.hh>


class SensitiveHit : public G4VHit
{
public:
  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void SetDeltaEnergy(G4double deltaE) { fDeltaEnergy = deltaE; }
  void SetKinEnergy(G4double kinE) { fKEnergy = kinE; }
  void SetTime(G4double time) { fTime = time; }
  void SetTimeL(G4double timeL) { fTimeL = timeL;}
  void SetPosition(G4ThreeVector pos) { fPosition = pos; }
  void SetPosition2(G4ThreeVector pos2) { fPosition2 = pos2; }
  void SetNbCopy(G4int scint) {fscintID = scint;}
  void SetParID(G4int Parentid) {fParID = Parentid;}
  void SetStepID(G4int stepid) {fStepID = stepid;}
  void SetParName(G4String particle) {fParName = particle;}
  void SetPrcName(G4String process) {fPrcName = process;}
  void SetVolName(G4int volume) {fVolName = volume;}
  void SetVolName2(G4String volume2){fVolName2=volume2;}
  void SetTrackID(G4int Trackid) {fTrackID = Trackid;}

  G4double GetDeltaEnergy() const { return fDeltaEnergy; }
  G4double GetKEnergy() const {return fKEnergy;}
  G4double GetTime() const { return fTime; }
  G4double GetTimeL() const {return fTimeL; }
  G4ThreeVector GetPosition() const { return fPosition; }
  G4ThreeVector GetPosition2() const { return fPosition2; }
  G4int GetNbCopy() const { return fscintID;}
  G4int GetParID() const {return fParID;}
  G4int GetStepID() const {return fStepID;}
  G4String GetParName() const {return fParName;}
  G4String GetPrcName() const {return fPrcName;}
  G4int GetVolName() const {return fVolName;}
  G4String GetVolName2() const {return fVolName2;}
  G4int GetTrackID() const {return fTrackID;}
 
private:
  G4double fDeltaEnergy,fKEnergy;
  G4double fTime;
  G4double fTimeL;
  G4ThreeVector fPosition;
  G4ThreeVector fPosition2;
  G4int fscintID;
  G4int fParID;
  G4int fStepID;
  G4String fParName;
  G4String fPrcName,fVolName2;
  G4int fVolName;
  G4int fTrackID;
 // G4ThreeVector fMomentum;
};

using SensitiveHitsCollection = G4THitsCollection<SensitiveHit>;

extern G4ThreadLocal G4Allocator<SensitiveHit> *hitAllocator;

inline void* SensitiveHit::operator new(size_t)
{
  if (!hitAllocator)
  {
      hitAllocator = new G4Allocator<SensitiveHit>;
  }
  return hitAllocator->MallocSingle();
}

inline void SensitiveHit::operator delete(void *aHit)
{
    if (!hitAllocator)
    {
        hitAllocator = new G4Allocator<SensitiveHit>;
    }
    hitAllocator->FreeSingle((SensitiveHit*) aHit);
}

#endif
