#include "G4SystemOfUnits.hh"
#include "SensitiveHit.hh"
#include "G4RunManager.hh"
#include "EventAction.hh"
#include "G4SDManager.hh"
#include "InitConfig.hh"
#include "G4THitsMap.hh"
#include "G4Event.hh"
#include "Analysis.hh"

//Progress bar -WZ
#include "G4Run.hh"

using namespace std;

EventAction::EventAction() : G4UserEventAction()
{;}

EventAction::~EventAction()
{;}

void EventAction::BeginOfEventAction(const G4Event* anEvent)
{
    // will use later for more information - SKS

    //Progress bar -WZ
    G4int eventID=anEvent->GetEventID();
    G4Run* run = static_cast<G4Run*>( G4RunManager::GetRunManager()->GetNonConstCurrentRun() );
    G4int nOfEvents = run->GetNumberOfEventToBeProcessed();
    G4double perCent = 1.; // status increment in percent
  
    if(fmod(eventID,double(nOfEvents*perCent*0.01))==0)
    {
      time_t my_time = time(NULL);
      tm *ltm = localtime(&my_time);
      G4double status=(100*(eventID/double(nOfEvents))); 
      std::cout << "=> Event " << eventID << " starts ("<< status << "%, "<< ltm->tm_hour << ":" <<  ltm->tm_min << ":" << ltm->tm_sec << ")" << std::endl;
    }
}

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

    if (numberHits)
      PrepareReactionLabels();

    for(int i1=0; i1<numberHits; i1++) {
      auto hit = (*hitsColl)[i1];
      G4String processCheck = hit->GetPrcName();
      if (processCheck == "Transportation")
        continue;

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
      G4String hitLabel = FindLabel(hit->GetTrackID());
      analysis->FillNtupleSColumn(22, hitLabel);

      analysis->AddNtupleRow();

      if (energy_step > 14.1) {
        G4cout << " Too high energy step in event = " << evnt->GetEventID() << G4endl;
        tempEvent = evnt->GetEventID();
      }
    }

    if (totEdep > 0)
      analysis->FillNtupleDColumn(16, totEdep);

    if (hitsVetoColl && fFillAlphaDetectorFields == "t") {
      int numberHits = hitsVetoColl->entries();

      for(int i1=0; i1<numberHits; i1++) {
        auto hit = (*hitsVetoColl)[i1];
        G4String processCheck = hit->GetPrcName();

        double energy_step=0.;
        energy_step = hit->GetDeltaEnergy()/MeV;

        analysis->FillNtupleIColumn(15, evnt->GetEventID());
        analysis->FillNtupleDColumn(25, energy_step);
        analysis->FillNtupleDColumn(26, hit->GetTime() / us);
        G4ThreeVector position = hit->GetPosition();
        analysis->FillNtupleDColumn(27, position.getX() / cm);
        analysis->FillNtupleDColumn(28, position.getY() / cm);
        analysis->FillNtupleDColumn(29, position.getZ() / cm);
        analysis->FillNtupleSColumn(30, hit->GetParName());

        analysis->AddNtupleRow();
      }
    }

    if (tempEvent == evnt->GetEventID())
      G4cout << "////////////////////////////////////////////" << G4endl;
  }

  fNeutronTimeVsInteraction.clear();
  fParticleTrackTime.clear();
}

void EventAction::PrepareReactionLabels()
{
  fTimeReactionLabels.clear();
  G4String label;
  for (auto it=fNeutronTimeVsInteraction.begin(); it!=fNeutronTimeVsInteraction.end(); it++) {
    label = "";
    double time = it->first;

    if (it->second == InteractionType::fNeuCapture)
      label = "nC";
    else if (it->second == InteractionType::fNeuInelastic)
      label = "nI";
    else
      label = "nO";

    int noOfGamma = 0;
    for (unsigned i=0; i<fParticleTrackTime.size(); i++) {
      if (std::get<2>(fParticleTrackTime.at(i)) == time) {
        if (std::get<0>(fParticleTrackTime.at(i)) == "gamma")
          noOfGamma++;
        else
          label += "_" + std::get<0>(fParticleTrackTime.at(i));
      }
    }
    fTimeReactionLabels.insert(std::pair{time, label});
  }
}

G4String EventAction::FindLabel(int trackID) {
  for (unsigned i=0; i<fParticleTrackTime.size(); i++) {
    if (std::get<1>(fParticleTrackTime.at(i)) == trackID) {
      auto it = fTimeReactionLabels.find(std::get<2>(fParticleTrackTime.at(i)));
      return it->second;
    }
  }
  G4String labelForOthers;
  if (trackID == 1)
    labelForOthers = "primaryN";
  else
    labelForOthers = "secondary";
  return labelForOthers;
}
