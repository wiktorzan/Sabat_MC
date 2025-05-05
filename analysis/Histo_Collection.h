#ifndef Histo_Collection_h
#define Histo_Collection_h

#include "TLatex.h"
#include "TAxis.h"
#include "TFile.h"
#include "TLine.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

class HistCollection
{
public:
  HistCollection() {;};
  ~HistCollection() {;};

  void CreateEnergyHistos(double minEnergy, double maxEnergy, double binSize)
  {
    EnergyDeposition = new TH1D("EnergyDeposition","Energy Deposition; Energy [MeV]; Counts", (maxEnergy - minEnergy)/binSize,
                                minEnergy - 0.5*binSize, maxEnergy - 0.5*binSize);
    EnergyDepositionSmeared = new TH1D("EnergyDepositionSmeared","Energy Deposition Smeared; Energy [MeV]; Counts", (maxEnergy - minEnergy)/binSize,
                                       minEnergy - 0.5*binSize, maxEnergy - 0.5*binSize);
    EnergyDepositionVeto = new TH1D("EnergyDepositionVeto","Energy Deposition Veto; Energy [MeV]; Counts", (5*maxEnergy - minEnergy)/binSize/5,
                                minEnergy - 0.5*binSize*5, 5*maxEnergy - 0.5*binSize*5);

    EnergyDepositionWithVeto = new TH1D("EnergyDepositionWithVeto","Energy Deposition; Energy [MeV]; Counts", (maxEnergy - minEnergy)/binSize,
                                minEnergy - 0.5*binSize, maxEnergy - 0.5*binSize);
    EnergyDepositionWithVetoSmeared = new TH1D("EnergyDepositionWithVetoSmeared","Energy Deposition Smeared; Energy [MeV]; Counts", (maxEnergy - minEnergy)/binSize,
                                       minEnergy - 0.5*binSize, maxEnergy - 0.5*binSize);
  }

  void CreateTimingHistos(double minTime, double timeSeparator, double maxTime, double binSizeSmall, double binSizeLarge)
  {
    TimeLaBr = new TH1D("TimeLaBr","Time; Time [us]; Counts", (timeSeparator - minTime)/binSizeSmall,
                                   minTime - 0.5*binSizeSmall, timeSeparator - 0.5*binSizeSmall);
    TimeVeto = new TH1D("TimeVeto","Time; Time [us]; Counts", (timeSeparator - minTime)/binSizeSmall,
                        minTime - 0.5*binSizeSmall, timeSeparator - 0.5*binSizeSmall);
    TimeDifferenceSmall = new TH1D("TimeDifferenceSmall","Time Difference; Time [us]; Counts", (timeSeparator - minTime)/binSizeSmall,
                                minTime - 0.5*binSizeSmall, timeSeparator - 0.5*binSizeSmall);
    TimeDifferenceLarge = new TH1D("TimeDifferenceLarge","Time Difference; Time [us]; Counts", (maxTime - timeSeparator)/binSizeLarge,
                                       timeSeparator - 0.5*binSizeLarge, maxTime - 0.5*binSizeLarge);
    timeSep = timeSeparator;
  }

  void CreateEnergyVsTimingHistos(double minEnergy, double maxEnergy, double binSize, double minTime, double maxTime, double binSizeLarge)
  {
    double eneBinSize = 2*binSize;
    double timeBinSize = 5*binSizeLarge;
    EnergyDepositionVsTimeDiff = new TH2D("EnergyDepositionVsTimeDiff","Energy Deposition; Energy [MeV]; Time [us]", (maxEnergy - minEnergy)/eneBinSize,
                                          minEnergy - 0.5*eneBinSize, maxEnergy - 0.5*eneBinSize, (maxTime - minTime)/timeBinSize,
                                          minTime - 0.5*timeBinSize, maxTime - 0.5*timeBinSize);
    EnergyDepositionVsTimeDiffSmeared = new TH2D("EnergyDepositionWithVetoSmeared","Energy Deposition Smeared; Energy [MeV]; Time [us]",
                                                 (maxEnergy - minEnergy)/eneBinSize, minEnergy - 0.5*eneBinSize, maxEnergy - 0.5*eneBinSize,
                                                 (maxTime - minTime)/timeBinSize, minTime - 0.5*timeBinSize, maxTime - 0.5*timeBinSize);
  }

// If more histos create table/vector/container and add references to the common functions
  void FillEnergyDeposition(double energy) {EnergyDeposition->Fill(energy);}
  void FillEnergyDepositionSmeared(double energy) {EnergyDepositionSmeared->Fill(energy);}
  void FillEnergyDepositionVeto(double energy) {EnergyDepositionVeto->Fill(energy);}
  void FillEnergyDepositionWithVeto(double energy) {EnergyDepositionWithVeto->Fill(energy);}
  void FillEnergyDepositionWithVetoSmeared(double energy) {EnergyDepositionWithVetoSmeared->Fill(energy);}
  void FillTimeLaBr(double time) {TimeLaBr->Fill(time);}
  void FillTimeVeto(double time) {TimeVeto->Fill(time);}
  void FillTimeDifference(double time)
  {
    if (time < timeSep)
      TimeDifferenceSmall->Fill(time);
    else
      TimeDifferenceLarge->Fill(time);
  }
  void FillEnergyDepositionVsTimeDiff(double energy, double time) {EnergyDepositionVsTimeDiff->Fill(energy, time);}
  void FillEnergyDepositionVsTimeDiffSmeared(double energy, double time) {EnergyDepositionVsTimeDiffSmeared->Fill(energy, time);}

  void SaveHistos(TString output)
  {
    TFile* outfile = new TFile(output, "RECREATE");
    outfile->cd();
    EnergyDeposition->Write("EnergyDeposition");
    EnergyDepositionSmeared->Write("EnergyDepositionSmeared");
    EnergyDepositionVeto->Write("EnergyDepositionVeto");
    EnergyDepositionWithVeto->Write("EnergyDepositionWithVeto");
    EnergyDepositionWithVetoSmeared->Write("EnergyDepositionWithVetoSmeared");
    TimeLaBr->Write("TimeLaBr");
    TimeVeto->Write("TimeVeto");
    TimeDifferenceSmall->Write("TimeDifferenceSmall");
    TimeDifferenceLarge->Write("TimeDifferenceLarge");
    EnergyDepositionVsTimeDiff->Write("EnergyDepositionVsTimeDiff");
    EnergyDepositionVsTimeDiffSmeared->Write("EnergyDepositionVsTimeDiffSmeared");
    outfile->Close();
  }

private:
  double timeSep;

  TH1D *EnergyDeposition;
  TH1D *EnergyDepositionSmeared;
  TH1D *EnergyDepositionVeto;
  TH1D *EnergyDepositionWithVeto;
  TH1D *EnergyDepositionWithVetoSmeared;
  TH1D *TimeLaBr;
  TH1D *TimeVeto;
  TH1D *TimeDifferenceSmall;
  TH1D *TimeDifferenceLarge;

  TH2D *EnergyDepositionVsTimeDiff;
  TH2D *EnergyDepositionVsTimeDiffSmeared;
};

#endif
