#ifndef STACKINGACTION_HH
#define STACKINGACTION_HH

#include <G4UserStackingAction.hh>
#include <G4Track.hh>

class RunAction;

class StackingAction : public G4UserStackingAction
{
public:
  StackingAction(RunAction* const );
  ~StackingAction(){;};

  G4ClassificationOfNewTrack   ClassifyNewTrack (const G4Track*);

private:
  RunAction* fRunAction;
};

#endif
