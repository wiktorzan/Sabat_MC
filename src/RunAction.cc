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

//automatic name creation
#include <chrono>


RunAction::RunAction()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFirstNtupleId(0);
  analysisManager->CreateNtuple("ekin_time", "Energy and time");
  analysisManager->CreateNtupleDColumn("EnergyDeposit");
  analysisManager->CreateNtupleDColumn("Time");
  analysisManager->CreateNtupleDColumn("X");
  analysisManager->CreateNtupleDColumn("Y");
  analysisManager->CreateNtupleDColumn("Z");
  analysisManager->CreateNtupleIColumn("count_rate");
  analysisManager->CreateNtupleIColumn("ParentID");
  analysisManager->CreateNtupleIColumn("StepID");
  analysisManager->CreateNtupleSColumn("Particles");
  analysisManager->CreateNtupleSColumn("Process");
  analysisManager->CreateNtupleDColumn("TimeL");
  analysisManager->CreateNtupleDColumn("KEnergy");
  analysisManager->CreateNtupleIColumn("volume");
  analysisManager->CreateNtupleSColumn("volume2");
  analysisManager->CreateNtupleIColumn("TrackID");
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleDColumn("TotalEdepGamma");
  analysisManager->CreateNtupleDColumn("NeutronTheta");
  analysisManager->CreateNtupleDColumn("NeutronPhi");
  analysisManager->CreateNtupleDColumn("AlphaTheta");
  analysisManager->CreateNtupleDColumn("AlphaPhi");
  analysisManager->CreateNtupleDColumn("X2");
  analysisManager->CreateNtupleDColumn("Y2");
  analysisManager->CreateNtupleDColumn("Z2");
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

  //automatic name creation

  std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
  std::time_t in_time_t = std::chrono::system_clock::to_time_t(now);
  std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

  std::tm* time_info = std::localtime(&in_time_t);
  std::ostringstream oss;
  oss << std::put_time(time_info, "%Y-%m-%d_%H_%M_%S");

  std::string filename = "output_";
  filename += oss.str();
  filename += "_";
  filename += std::to_string(ms.count());
  filename += ".root";
  
  std::cout << "Filename: " << filename << std::endl;

  analysisManager->OpenFile(filename);

  // analysisManager->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
}

void RunAction::AddSecondary(const G4ParticleDefinition* particle, G4double energy)
{
  return;
}
