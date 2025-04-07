#include <sys/stat.h>
#include <functional>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <sstream>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>

#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

#include "Histo_Collection.h"

//Smearing parameters
double a = 2.0*pow(10, -4); // in MeV
double b = 2.22*pow(10, -2);
double c = 0.5;

double LaBrDecayTime = 0.025; // in us
double SiDecayTime = 0.09; // in us for silicon drift detectors

bool FileCheck(const std::string& NameOfFile);
double SmearEnergy(double energy);
void AnalyzeFile(std::string NameOfFile, HistCollection histo);

int main(int argc, char* argv[])
{
//--------------------------------------------------------
//Reading parameters
  if (argc < 2) {
    std::cout << "Not enough arguments. Up to two needed: name of the file | name of the second file ..." << std::endl;
    return 0;
  }

  std::string fileOrPattern = argv[1];

//--------------------------------------------------------
//Setting environment for results
  Double_t minEnergy = 0, maxEnergy = 24, binSize = 0.01; // in MeV
  Double_t minTime = 0, timeSeparator = 1, maxTime = 1500, binSizeSmall = 0.001, binSizeLarge = 1; // in us
  HistCollection hist;
  hist.CreateEnergyHistos(minEnergy, maxEnergy, binSize);
  hist.CreateTimingHistos(minTime, timeSeparator, maxTime, binSizeSmall, binSizeLarge);
  TString outputName = "";
//--------------------------------------------------------
//Analysis
  std::vector<std::string> filesToAnalyze;

  if (argc > 2) {
    outputName = "Out_LastFile_" + fileOrPattern;
    filesToAnalyze.push_back(fileOrPattern);
    for (unsigned i=2; i<argc; i++) {
      fileOrPattern = argv[i];
      filesToAnalyze.push_back(fileOrPattern);
    }
    for (unsigned fileNo = 0; fileNo < filesToAnalyze.size(); fileNo++) {
      AnalyzeFile(filesToAnalyze.at(fileNo), hist);
    }
  } else {
    TString root_file = fileOrPattern;
    TString star = "*";
    TString slash = "/";
    std::size_t starPlace = fileOrPattern.rfind(star);
    std::size_t slashPlace = fileOrPattern.rfind(slash);

    if (root_file[strlen(root_file) - 1] == star) {
      std::string pattern;
      std::string dirName;
      if (slashPlace > 0) {
        dirName = fileOrPattern.substr(0, slashPlace+1);
        pattern = fileOrPattern.substr(slashPlace+1, starPlace);
      }
      else {
        dirName = "";
        pattern = fileOrPattern.substr(0, starPlace);

      }
      std::cout << "Getting files from directory: " << dirName << " and pattern " << pattern << std::endl;
      TString dirNameRoot = dirName;
      if (dirNameRoot == "")
        dirNameRoot = ".";
      TSystemDirectory dir(dirNameRoot, dirNameRoot);
      TList *files = dir.GetListOfFiles();

      if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
          fname = file->GetName();

          if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith(pattern)) {
            std::cout << "Adding file to the analysis queue " << fname.Data() << std::endl;
            filesToAnalyze.push_back(dirName + fname.Data());
          }
        }
      }
      outputName = "Out_Patt_" + pattern + ".root";
    } else {
      filesToAnalyze.push_back(fileOrPattern);
      outputName = "Out_" + fileOrPattern;
    }

    for (unsigned fileNo = 0; fileNo < filesToAnalyze.size(); fileNo++) {
      AnalyzeFile(filesToAnalyze.at(fileNo), hist);
    }
  }
  hist.SaveHistos(outputName);

  return 0;
}

void AnalyzeFile(std::string NameOfFile, HistCollection histo)
{
  TString fileName = NameOfFile;
  std::cout << " Reading file " << fileName << std::endl;
  TFile* hfile = new TFile(fileName, "READ");
  TTree *ntuple = (TTree *) hfile->Get("EventTree");
  Double_t Energy_Deposit, Neutron_Theta, Neutron_Phi, Time, Hit_X, Hit_Y, Hit_Z;
  Double_t Veto_Energy_Deposit, Veto_Time, Veto_Hit_X, Veto_Hit_Y ,Veto_Hit_Z;
  Int_t Event, EventVeto, Parent_ID;
  Char_t Particles[6], Veto_Particles[6], Process[22];
  Int_t Volume;
  Char_t Volume2[13];
  Int_t Event_ID;
  ntuple->SetBranchAddress("Energy_Deposition", &Energy_Deposit);
  ntuple->SetBranchAddress("Time", &Time);
  ntuple->SetBranchAddress("Hit_X", &Hit_X);
  ntuple->SetBranchAddress("Hit_Y", &Hit_Y);
  ntuple->SetBranchAddress("Hit_Z", &Hit_Z);
  ntuple->SetBranchAddress("Parent_ID", &Parent_ID);
  ntuple->SetBranchAddress("Particles", Particles);
  ntuple->SetBranchAddress("Process", Process);
  ntuple->SetBranchAddress("Volume", &Volume);
  ntuple->SetBranchAddress("Volume2", Volume2);
  ntuple->SetBranchAddress("Event_ID", &Event_ID);
  ntuple->SetBranchAddress("Neutron_Theta", &Neutron_Theta);
  ntuple->SetBranchAddress("Neutron_Phi", &Neutron_Phi);
  ntuple->SetBranchAddress("Veto_Energy_Deposition", &Veto_Energy_Deposit);
  ntuple->SetBranchAddress("Veto_Time", &Veto_Time);
  ntuple->SetBranchAddress("Veto_Hit_X", &Veto_Hit_X);
  ntuple->SetBranchAddress("Veto_Hit_Y", &Veto_Hit_Y);
  ntuple->SetBranchAddress("Veto_Hit_Z", &Veto_Hit_Z);
  ntuple->SetBranchAddress("Veto_Particles", Veto_Particles);

  Int_t nentries = (Int_t)ntuple->GetEntries();

  Double_t EnergyDepositFinal = 0, EnergyDepositVetoFinal = 0;
  Double_t FirstTime, FirstTimeVeto;
  for (Int_t i=0; i<nentries; i++) {
    ntuple->GetEntry(i);

 //   std::cout << i << " " << Event_ID << std::endl;

    if (Time > 0) {
      histo.FillTimeLaBr(Time);
      if (Event != Event_ID) {
//std::cout << "Time LaBr " << Time << " and Energy " << EnergyDepositFinal << " " << Event_ID << std::endl;
//std::cin >> Parent_ID;
        if (EnergyDepositFinal > 0) {
          histo.FillEnergyDeposition(EnergyDepositFinal);
          histo.FillEnergyDepositionSmeared(EnergyDepositFinal + SmearEnergy(EnergyDepositFinal));

          if (Veto_Time > 0) {
            histo.FillTimeDifference(FirstTimeVeto - FirstTime);
          }
        }

        FirstTime = Time;
        EnergyDepositFinal = 0;
        Event = Event_ID;
      } else {
        if (strcmp(Volume2, "DetectorLaBr") == 0 && (strcmp(Particles, "gamma") == 0 || strcmp(Particles, "e-") == 0))
          EnergyDepositFinal = EnergyDepositFinal + Energy_Deposit;
      }
    }

    if (Veto_Time > 0) {

      histo.FillTimeVeto(Veto_Time);
      if (EventVeto != Event_ID) {
//std::cout << "Time Veto " << Veto_Time << " and Energy " << EnergyDepositVetoFinal << " " << Event_ID << std::endl;
//std::cin >> Parent_ID;
        if (EnergyDepositVetoFinal > 0) {
          histo.FillEnergyDepositionVeto(EnergyDepositVetoFinal/100);
        }

        FirstTimeVeto = Veto_Time;
        EnergyDepositVetoFinal = 0;
        EventVeto = Event_ID;
      } else {
        if ((strcmp(Veto_Particles, "alpha") == 0 || strcmp(Veto_Particles, "e-") == 0))
          EnergyDepositVetoFinal = EnergyDepositVetoFinal + Veto_Energy_Deposit;
      }
    }
  }
}

double SmearEnergy(double energy)
{
  return gRandom->Gaus(0,1)*(a + b*sqrt(energy + c*pow(energy, 2)))/(2.35482004503);
}

bool FileCheck(const std::string& NameOfFile)
{
    struct stat buffer;
    return (stat(NameOfFile.c_str(), &buffer) == 0);
}
