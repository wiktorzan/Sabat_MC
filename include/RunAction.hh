#ifndef RUNACTION_HH
#define RUNACTION_HH

#include "G4ParticleDefinition.hh"
#include "G4UserRunAction.hh"
#include "RunMessenger.hh"
#include "G4Run.hh"

class RunAction : public G4UserRunAction
{
public:
  RunAction();
  ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void AddSecondary(const G4ParticleDefinition*, G4double energy);
  void AddTrackLength(G4double length);
  void AddTimeAndSeed() {fOutputAddTimeAndSeed = "t";};
  void SetTimeAndSeed(std::string timeAndSeed) {fTimeAndSeed = timeAndSeed;};
  void RemoveAlphaGen() {fIncludeAlphaDetectorFields = "n";};
  std::string GetFlagForAlphaDetectorFields() {return fIncludeAlphaDetectorFields;};
private:
  RunMessenger*  fRunMess;

  std::string fOutputAddTimeAndSeed = "n";
  std::string fTimeAndSeed = "";
  std::string fIncludeAlphaDetectorFields = "t";
};

#endif
