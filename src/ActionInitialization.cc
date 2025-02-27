/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"

#include "RunAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"

//! Class constructor
ActionInitialization::ActionInitialization() :
  G4VUserActionInitialization()
{;}

//! Class destructor
ActionInitialization::~ActionInitialization()
{;}

void ActionInitialization::Build() const
{
    
    SetUserAction(new PrimaryGeneratorAction());


    RunAction* theRunAction = new RunAction();
    SetUserAction(theRunAction);

    SetUserAction(new EventAction());

    SetUserAction(new StackingAction(theRunAction));

    SetUserAction(new SteppingAction(theRunAction));

}

void ActionInitialization::BuildForMaster() const
{
   RunAction* runAction = new RunAction();
   SetUserAction(runAction); 


}
