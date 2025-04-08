/// \file EventAction.cc
/// \brief Implementation of the EventAction class
#include "G4SystemOfUnits.hh"
#include "SensitiveHit.hh"
#include "G4RunManager.hh"
#include "EventAction.hh"
#include "G4SDManager.hh"
#include "InitConfig.hh"
#include "G4THitsMap.hh"
#include "G4Event.hh"
#include "Analysis.hh"

using namespace std;

EventAction::EventAction() : G4UserEventAction()
{;}

EventAction::~EventAction()
{;}

void EventAction::BeginOfEventAction(const G4Event* /*anEvent*/)
{;}

void EventAction::EndOfEventAction(const G4Event* anEvent)
{
  G4SDManager* sdm = G4SDManager::GetSDMpointer();

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->SetDefaultFileType("root");

  G4HCofThisEvent* hcofEvent = anEvent->GetHCofThisEvent();

  if (!hcofEvent)
    return;

  if (fScintillatorId<0) {
    fScintillatorId = sdm->GetCollectionID("Detector/LaBr_det");
    G4cout << "Event action : Scintillator id is = " << fScintillatorId << G4endl;

    fScintillatorIdVeto = sdm->GetCollectionID("Veto/Veto_det");
    G4cout << "Event action : Veto scintillator id is = " << fScintillatorIdVeto << G4endl;
  }

  const G4Event* evnt = G4RunManager::GetRunManager()->GetCurrentEvent();
  SensitiveHitsCollection* hitsColl=0;
  SensitiveHitsCollection* hitsVetoColl=0;

  if(hcofEvent) {
    hitsColl = dynamic_cast<SensitiveHitsCollection*>(hcofEvent->GetHC(fScintillatorId));
    hitsVetoColl = dynamic_cast<SensitiveHitsCollection*>(hcofEvent->GetHC(fScintillatorIdVeto));
  }

  if (hitsColl) {
    int numberHits = hitsColl->entries(), tempEvent=0;

    double totEdep = 0.; // in MeV, but trying to set as double to fix the issue with zero TotEDep

    for(int i1=0; i1<numberHits; i1++) {
      auto hit = (*hitsColl)[i1];
      G4String processCheck = hit->GetPrcName();

      double energy_step=0.;
      energy_step = hit->GetDeltaEnergy()/MeV;
      totEdep += energy_step;

      analysis->FillNtupleDColumn(0, energy_step);
      analysis->FillNtupleDColumn(1, hit->GetTime() / us);
      G4ThreeVector position = hit->GetPosition();
      analysis->FillNtupleDColumn(2, position.getX() / cm);
      analysis->FillNtupleDColumn(3, position.getY() / cm);
      analysis->FillNtupleDColumn(4, position.getZ() / cm);
      analysis->FillNtupleIColumn(5, hit->GetNbCopy());
      analysis->FillNtupleIColumn(6, hit->GetParID());
      analysis->FillNtupleIColumn(7, hit->GetStepID());
      analysis->FillNtupleSColumn(8, hit->GetParName());
      analysis->FillNtupleSColumn(9, hit->GetPrcName());
      analysis->FillNtupleDColumn(10, hit->GetTimeL() / us);
      analysis->FillNtupleDColumn(11, hit->GetKEnergy()/MeV);
      analysis->FillNtupleIColumn(12, hit->GetVolName());
      analysis->FillNtupleSColumn(13, hit->GetVolName2());
      analysis->FillNtupleIColumn(14, hit->GetTrackID());
      analysis->FillNtupleIColumn(15, evnt->GetEventID());
      G4ThreeVector positionPostStep = hit->GetPosition2();
      analysis->FillNtupleDColumn(19, positionPostStep.getX() / cm);
      analysis->FillNtupleDColumn(20, positionPostStep.getY() / cm);
      analysis->FillNtupleDColumn(21, positionPostStep.getZ() / cm);

      analysis->AddNtupleRow();

      if (energy_step > 14.1) {
        G4cout << " Too high energy step in event = " << evnt->GetEventID() << G4endl;
        tempEvent = evnt->GetEventID();
      }
    }
        
    if (totEdep > 0)
      analysis->FillNtupleDColumn(16, totEdep);

    if (tempEvent == evnt->GetEventID())
      G4cout << "////////////////////////////////////////////" << G4endl;
  }

  if (hitsVetoColl && fFillAlphaDetectorFields == "t") {
    int numberHits = hitsVetoColl->entries();

    for(int i1=0; i1<numberHits; i1++) {
      auto hit = (*hitsVetoColl)[i1];
      G4String processCheck = hit->GetPrcName();

      double energy_step=0.;
      energy_step = hit->GetDeltaEnergy()/MeV;

      analysis->FillNtupleIColumn(15, evnt->GetEventID());
      analysis->FillNtupleDColumn(24, energy_step);
      analysis->FillNtupleDColumn(25, hit->GetTime() / us);
      G4ThreeVector position = hit->GetPosition();
      analysis->FillNtupleDColumn(26, position.getX() / cm);
      analysis->FillNtupleDColumn(27, position.getY() / cm);
      analysis->FillNtupleDColumn(28, position.getZ() / cm);
      analysis->FillNtupleSColumn(29, hit->GetParName());

      analysis->AddNtupleRow();
    }
  }
}
