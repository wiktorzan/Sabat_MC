/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include <G4UserEventAction.hh>
#include <globals.hh>

#include <tuple>
#include <map>

enum InteractionType {
  fNeuCapture, fNeuInelastic, fNeuOther
};

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction();

  void BeginOfEventAction(const G4Event* anEvent) override;
  void EndOfEventAction(const G4Event* anEvent) override;
  void SetAlphaDetectorsFieldsFlag(std::string flag) {fFillAlphaDetectorFields = flag;};
  void AddNewNeutronInteraction(G4double time, std::string inter) {
    if (inter == "nCapture")
      fNeutronTimeVsInteraction.insert(std::pair{time, InteractionType::fNeuCapture});
    else if (inter == "neutronInelastic")
      fNeutronTimeVsInteraction.insert(std::pair{time, InteractionType::fNeuInelastic});
    else
      fNeutronTimeVsInteraction.insert(std::pair{time, InteractionType::fNeuOther});
  };
  void AddParticleEnergyTime(std::string particleName, G4int trackID, G4double time) {
    fParticleTrackTime.push_back(std::make_tuple(particleName, trackID, time));
  };
  void PrepareReactionLabels();
  G4String FindLabel(int trackID);
  G4String GetAddedLabel() {
    G4String temp = "Unknown";
    if (fTimeReactionLabels.size() == 1)
      temp = fTimeReactionLabels.begin()->second;
    return temp;
  };
private:
  G4int fScintillatorId{-1};
  G4int fScintillatorIdVeto{-1};
  G4String fFillAlphaDetectorFields = "t";

  std::map<G4double, InteractionType> fNeutronTimeVsInteraction;
  std::vector<std::tuple<G4String, G4int, G4double>> fParticleTrackTime;
  std::map<G4double, G4String> fTimeReactionLabels;
};

#endif
