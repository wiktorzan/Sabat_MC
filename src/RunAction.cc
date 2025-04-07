#include "PrimaryGeneratorAction.hh"
#include "G4AccumulableManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Electron.hh"
#include "InitConfig.hh"
#include "RunAction.hh"
#include "Analysis.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4Alpha.hh"
#include "iostream"
#include "iomanip"
#include "string"


RunAction::RunAction()
{
  InitConfig* init = InitConfig::getInstance();
  fOutputAddTimeAndSeed = init->GetVariable("filenameAddTimeAndSeed");
  if (fOutputAddTimeAndSeed == "t" || fOutputAddTimeAndSeed == "true") {
    fTimeAndSeed = init->GetTimeAndSeed();
  }
  fIncludeAlphaDetectorFields = init->GetVariable("includeAlphaDetection");
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
  if (fIncludeAlphaDetectorFields == "t") {
    analysisManager->CreateNtupleDColumn("Veto_Energy_Deposition");
    analysisManager->CreateNtupleDColumn("Veto_Time");
    analysisManager->CreateNtupleDColumn("Veto_Hit_X");
    analysisManager->CreateNtupleDColumn("Veto_Hit_Y");
    analysisManager->CreateNtupleDColumn("Veto_Hit_Z");
    analysisManager->CreateNtupleSColumn("Veto_Particles");
  }
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

  if (fTimeAndSeed != "") {
    std::string oldName = analysisManager->GetFileName();
    int pos = oldName.find_first_of('.');
    oldName = oldName.substr(0, pos);
    //Change dir to ../Data
    std::string path = "../DATA/" + oldName + "_" + fTimeAndSeed + ".root";
    // std::string path = oldName + "_" + fTimeAndSeed + ".root";
    std::cout << path << std::endl;
    analysisManager->SetFileName(path);
  }else
  {
    std::cout << "Brak seed+time" << std::endl;
  }

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
