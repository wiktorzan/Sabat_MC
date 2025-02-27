#include "StackingAction.hh"
#include "RunAction.hh"

#include <G4SystemOfUnits.hh>

StackingAction::StackingAction(RunAction* aRunAction) :
  G4UserStackingAction(),fRunAction(aRunAction)
{;}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack (const G4Track*
 aTrack)
{
   
  return G4UserStackingAction::ClassifyNewTrack(aTrack);
}
