#ifndef InitConfig_h
#define InitConfig_h 1

#include <G4SystemOfUnits.hh>
#include <globals.hh>

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

  void SetFileName(std::string fileName) {initName = fileName;};
  void Initialization();
  void Read();

private:
  static InitConfig* fInit;
  std::string initName = "initConfig.dat";
  std::map<G4String, G4String> variables;
};

#endif
