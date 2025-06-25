#include "SteppingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"

#include "G4HadronicProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4UnitsTable.hh"
#include <G4Electron.hh>
#include <G4Track.hh>
#include <G4Step.hh>

#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;

SteppingAction::SteppingAction(RunAction* runAction, EventAction* evAction) : fRunAction(runAction), fEventAction(evAction)
{}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4int Event;
  const G4Event* evnt = G4RunManager::GetRunManager()->GetCurrentEvent();
  Event = evnt->GetEventID();
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchable()->GetVolume();
  G4VPhysicalVolume* volume1 = aStep->GetPostStepPoint()->GetTouchable()->GetVolume();
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  G4Track* theTrack = aStep->GetTrack();
  G4int trackID = theTrack->GetTrackID();
  G4String parName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
  G4TouchableHandle touch = postStepPoint->GetTouchableHandle();

  if (preStepPoint->GetGlobalTime() == fCurrNeuTime) {
    if (std::find(fUniqueTracks.begin(), fUniqueTracks.end(), trackID) == fUniqueTracks.end()) {
      fEventAction->AddParticleEnergyTime(parName, trackID, preStepPoint->GetGlobalTime());
      fUniqueTracks.push_back(trackID);
    }
  }

  if (parName == "neutron") {
    if (theTrack->GetKineticEnergy() == 0) {
      fEventAction->AddNewNeutronInteraction(postStepPoint->GetGlobalTime(), postStepPoint->GetProcessDefinedStep()->GetProcessName());
      fCurrNeuTime = postStepPoint->GetGlobalTime();
      fUniqueTracks.clear();
    }
  }
}
