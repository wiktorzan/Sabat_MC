/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"

#include "RunAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization() : G4VUserActionInitialization()
{;}

ActionInitialization::~ActionInitialization()
{;}

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction());

    RunAction* theRunAction = new RunAction();
    theRunAction->SetTimeAndSeed(timeAndSeedAdd);
    SetUserAction(theRunAction);

    EventAction* theEventAction = new EventAction();
    theEventAction->SetAlphaDetectorsFieldsFlag(theRunAction->GetFlagForAlphaDetectorFields());
    SetUserAction(theEventAction);

    SetUserAction(new StackingAction(theRunAction));

    SetUserAction(new SteppingAction(theRunAction));
}

void ActionInitialization::BuildForMaster() const
{
   RunAction* runAction = new RunAction();
   SetUserAction(runAction); 
}
