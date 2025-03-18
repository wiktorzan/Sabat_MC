#include "RunAction.hh"
#include <G4Gamma.hh>
#include <G4Electron.hh>
#include <G4AccumulableManager.hh>
#include <G4SystemOfUnits.hh>
#include "Analysis.hh"
#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hh"
#include <iomanip>
//#include "HistoManager.hh"
#include <G4Alpha.hh>
#include <G4Proton.hh>


RunAction::RunAction()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFirstNtupleId(0);

  analysisManager->CreateNtuple("EventTree", "Tracked Events");
  analysisManager->CreateNtupleDColumn("Energy_Deposition");
  analysisManager->CreateNtupleDColumn("Time");
  analysisManager->CreateNtupleDColumn("Hit_X");
  analysisManager->CreateNtupleDColumn("Hit_Y");
  analysisManager->CreateNtupleDColumn("Hit_Z");
  analysisManager->CreateNtupleIColumn("Count_rate");
  analysisManager->CreateNtupleIColumn("Parent_ID");
  analysisManager->CreateNtupleIColumn("Step_ID");
  analysisManager->CreateNtupleSColumn("Particles");
  analysisManager->CreateNtupleSColumn("Process");
  analysisManager->CreateNtupleDColumn("Time_L");
  analysisManager->CreateNtupleDColumn("Kinetic_Energy");
  analysisManager->CreateNtupleIColumn("Volume");
  analysisManager->CreateNtupleSColumn("Volume2");
  analysisManager->CreateNtupleIColumn("Track_ID");
  analysisManager->CreateNtupleIColumn("Event_ID");
  analysisManager->CreateNtupleDColumn("Total_Gamma_Energy_Deposition");
  analysisManager->CreateNtupleDColumn("Neutron_Theta");
  analysisManager->CreateNtupleDColumn("Neutron_Phi");
  analysisManager->CreateNtupleDColumn("Hit_X2");
  analysisManager->CreateNtupleDColumn("Hit_Y2");
  analysisManager->CreateNtupleDColumn("Hit_Z2");
  analysisManager->CreateNtupleDColumn("Alpha_Theta");
  analysisManager->CreateNtupleDColumn("Alpha_Phi");
  analysisManager->CreateNtupleDColumn("Veto_Energy_Deposition");
  analysisManager->CreateNtupleDColumn("Veto_Time");
  analysisManager->CreateNtupleDColumn("Vet_Hit_X");
  analysisManager->CreateNtupleDColumn("Veto_Hit_Y");
  analysisManager->CreateNtupleDColumn("Veto_Hit_Z");
  analysisManager->FinishNtuple();
}

RunAction::~RunAction()
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
}

void RunAction::AddSecondary(const G4ParticleDefinition* particle, G4double energy)
{
  return;
}
