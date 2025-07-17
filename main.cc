#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include <Randomize.hh>

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PhysicsList.hh"
#include "InitConfig.hh"

#include <sstream>
#include <unistd.h>
#include <chrono>

int main(int argc, char** argv)
{
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  int mask = 01001010;
  using namespace std::chrono;
  system_clock::time_point now = system_clock::now();
  long seed = (unsigned int)(system_clock::to_time_t(now)) * 677 * ::getpid();
  G4Random::setTheSeed(seed);

  std::time_t tt = std::chrono::system_clock::to_time_t(now);
  std::tm tm = *std::gmtime(&tt); //GMT (UTC)
  std::stringstream ss;
  const std::string& format = "%m%d%H%M%S";
  ss << std::put_time( &tm, format.c_str() );

  std::string seedAndTime = ss.str() + "_" + std::to_string(seed);

  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization(new PhysicsList());
  DetectorConstruction* detCons = new DetectorConstruction();
  runManager->SetUserInitialization(detCons);
  ActionInitialization* actInit = new ActionInitialization(detCons);
  actInit->SetTimeAndSeedAdd(seedAndTime);
  runManager->SetUserInitialization(actInit);
//  runManager->Initialize(); // initialization is in the macro

// Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
// Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

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