#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include <Randomize.hh>

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PhysicsList.hh"

#include <unistd.h>
#include <chrono>

int main(int argc, char** argv)
{
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  int mask = 01001010;
  using namespace std::chrono;
  system_clock::time_point now = system_clock::now();
  long seed = (unsigned int)(system_clock::to_time_t(now)) * 677 * ::getpid();
  G4Random::setTheSeed(seed);

  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(new ActionInitialization());
  runManager->Initialize();

    // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

    // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

    // Execute visualization macro (if in interactive mode)
  if (ui) {
    UImanager->ApplyCommand("/control/execute macros/vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  delete visManager;
  delete runManager;

  return 0;
}