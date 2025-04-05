/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include <G4UserEventAction.hh>
#include <globals.hh>

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction();

  void BeginOfEventAction(const G4Event* anEvent) override;
  void EndOfEventAction(const G4Event* anEvent) override;
  void SetAlphaDetectorsFieldsFlag(std::string flag) {fFillAlphaDetectorFields = flag;};
private:
  G4int fScintillatorId{-1};
  G4int fScintillatorIdVeto{-1};
  G4String fFillAlphaDetectorFields = "t";
};

#endif
