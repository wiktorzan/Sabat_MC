/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "PrimaryGeneratorAction.hh"
#include "ActionInitialization.hh"
#include "EventAction.hh"

#include "RunAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* detCons) : G4VUserActionInitialization()
{
    detConsPoint = detCons;
}

ActionInitialization::~ActionInitialization()
{;}

void ActionInitialization::Build() const
{
    PrimaryGeneratorAction* thePrimGenAction = new PrimaryGeneratorAction();
    detConsPoint->SetPrimGen(thePrimGenAction);
    SetUserAction(thePrimGenAction);

    RunAction* theRunAction = new RunAction();
    theRunAction->SetTimeAndSeed(timeAndSeedAdd);
    SetUserAction(theRunAction);

    EventAction* theEventAction = new EventAction();
    theEventAction->SetAlphaDetectorsFieldsFlag(theRunAction->GetFlagForAlphaDetectorFields());
    SetUserAction(theEventAction);

    SetUserAction(new StackingAction(theRunAction));

    SetUserAction(new SteppingAction(theRunAction, theEventAction));
}

void ActionInitialization::BuildForMaster() const
{
   RunAction* runAction = new RunAction();
   SetUserAction(runAction); 
}
