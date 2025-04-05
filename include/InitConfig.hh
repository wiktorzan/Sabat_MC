#ifndef InitConfig_h
#define InitConfig_h 1

#include <G4SystemOfUnits.hh>
#include <globals.hh>
#include <chrono>
#include <vector>
#include <map>

#include <iostream>
#include <fstream>

class InitConfig
{

public:
  InitConfig();
  ~InitConfig();

  static InitConfig* getInstance();
  std::string GetVariable(std::string key) {
    std::string empty = "";
    if (auto search = variables.find(key); search != variables.end())
      return variables.at(key);
    else
      return empty;
  };

  void SetTimeAndSeed(std::string seedPlusTime) {fSeedPlusTime = seedPlusTime;}
  void SetFileName(std::string fileName) {initName = fileName;};
  void Initialization();
  void Read();

  std::string GetTimeAndSeed() {return fSeedPlusTime;};

private:
  static InitConfig* fInit;
  std::string initName = "initConfig.dat";
  std::map<G4String, G4String> variables;

  std::string fSeedPlusTime;
};

#endif
