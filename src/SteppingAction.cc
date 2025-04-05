#include "SteppingAction.hh"
#include "RunAction.hh"
#include <G4Step.hh>
#include <G4Electron.hh>
#include <G4Track.hh>
#include <G4SystemOfUnits.hh>
#include "G4HadronicProcess.hh"
#include "EventAction.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include <iostream>
#include <iomanip>
using namespace std;

SteppingAction::SteppingAction(RunAction* runAction) : fRunAction(runAction)
{}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4int Event;
  const G4Event* evnt = G4RunManager::GetRunManager()->GetCurrentEvent();
  Event = evnt->GetEventID();
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchable()->GetVolume();
  G4VPhysicalVolume* volume1 = aStep->GetPostStepPoint()->GetTouchable()->GetVolume();
  G4StepPoint * preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint * postStepPoint = aStep->GetPostStepPoint();
  G4Track * theTrack = aStep->GetTrack();
  G4String ParName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
  G4TouchableHandle touch = postStepPoint->GetTouchableHandle();
}
