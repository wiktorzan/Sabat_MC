#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include <G4UserSteppingAction.hh>

#include <string>
#include <vector>
#include <tuple>
#include <list>
#include <map>

class RunAction;
class EventAction;

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(RunAction*, EventAction*);

    void UserSteppingAction(const G4Step*) override;

private:
    RunAction* fRunAction;
    EventAction* fEventAction;
    double fCurrNeuTime = -1;
    std::list<int> fUniqueTracks;
};

#endif
