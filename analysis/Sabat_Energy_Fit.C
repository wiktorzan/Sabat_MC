#include <sys/stat.h>
#include <functional>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <string>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <map>
#include <time.h>
#include <unistd.h>

#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMinuit.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TKey.h>
#include <TF12.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TF1.h>
#include <TF2.h>

enum TypeOfIntensityEstimation {
  Gauss, Integral
};

int typeOfBackground = 3;

int sizeOfHalfWindow = 5;
int orderOfPolynomialToInterpolate = 2;

std::vector<double> desiredLines = {0.512, 0.695, 0.93, 1.95, 2.738, 3.84, 4.44, 5.11, 5.74, 6.125, 6.67, 7.11, 7.79, 8.2, 8.58/*, 10.8*/}; // in MeV
double nitrogenOffset = 9.2; // in MeV. From which energy one should calculate background from nitrogen
std::vector<std::vector<double>> complexLines = {{0.785, 0.846},{1.17, 1.24}, {1.61, 1.72}, {2.12, 2.23, 2.31, 2.435}}; // in MeV

bool FileCheck(const std::string& NameOfFile);
std::string NumberToChar(double number, int precision);

void PlotDerivatives(TH1D* histo, std::string nameOfInput);
void EstimateBGfromMinima(TH1D* histoToFit, std::string nameOfInput);
void EstimateBGfromDerivative(TH1D* histoToFit, std::string nameOfInput);
void FitSpectrum(TH1D* histoToFit, std::string nameOfInput);
std::pair<double, double> FitSingleLine(TH1D* histoToFit, double energy, TypeOfIntensityEstimation type);
std::vector<std::pair<double, double>> FitMultipleLines(TH1D* histoToFit, std::vector<double> energies, TypeOfIntensityEstimation type);
void FitLineLocally(TH1D* histoToFit, std::string nameOfInput, TypeOfIntensityEstimation type);

Double_t eneLineFunction(Double_t *A, Double_t *P);
Double_t eneMultLinesFunction(Double_t *A, Double_t *P);
Double_t bgFunctionWithLogGauss(Double_t *A, Double_t *P);
double GaussDistr(double x, double mean, double sigma);
double LogGaussDistr(double x, double mean, double sigma);
std::pair<int, int> FindLocalEx(TH1D* derivHisto, int centerBin);
double GramPolynomial(int i, int m, int k, int s);
double GenFactor(int a, int b);
double CalcWeight(int i, int t, int m, int n, int s);
std::vector<double> ComputeWeights(int m, int t, int n, int s);

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cout << "Not enough arguments. Provide the name of the input root file" << std::endl;
    return 0;
  }
  std::string nameOfInput = argv[1];

  int checkOfExistence = FileCheck(nameOfInput);

  if (!checkOfExistence) {
    std::cout << "No such file: " << nameOfInput << std::endl;
    return 0;
  }

  TString fileToAnalyze = nameOfInput;

  std::cout << " Reading file " << fileToAnalyze << std::endl;
  TFile* inputFile = new TFile(fileToAnalyze, "READ");

  TIter next(inputFile->GetListOfKeys());
  TKey *key;
  TH1D* histoToFit;
  TH1D* temp1D;
  TH2D* temp2D;
  std::vector<TH1D*> HistogramsToPotentialFit1D;
  std::vector<TH2D*> HistogramsToPotentialFit2D;

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());

    if (cl->InheritsFrom("TH1D")) {
      temp1D = (TH1D*) key->ReadObj()->Clone();
      temp1D->SetDirectory(0);
      std::cout << "1D Histogram with name: " << key->GetName() << " have ID equal to " << HistogramsToPotentialFit1D.size() << std::endl;
      HistogramsToPotentialFit1D.push_back(temp1D);
    } else if (cl->InheritsFrom("TH2D")) {
      temp2D = (TH2D*) key->ReadObj()->Clone();
      temp2D->SetDirectory(0);
      std::cout << "2D Histogram with name: " << key->GetName() << " have ID equal to " << HistogramsToPotentialFit2D.size() << std::endl;
      HistogramsToPotentialFit2D.push_back(temp2D);
    }
  }

  TString histoName;

  if (HistogramsToPotentialFit1D.size() + HistogramsToPotentialFit2D.size() == 0) {
    std::cout << "No histogram to fit" << std::endl;
    return 0;
  } else {
    int IDofHistogram;
    TString histoDimension;
    std::cout << "Provide the ID of the histogram to fit:" << std::endl;
    std::cin >> IDofHistogram;
    std::cout << "Do you want to fit 1D or 2D histogram? Write 1D or 2D" << std::endl;
    std::cin >> histoDimension;

    if (histoDimension == "1D" || histoDimension == "1") {
      histoToFit = HistogramsToPotentialFit1D.at(IDofHistogram);
      histoName = histoToFit->GetName();
    } else if (histoDimension == "2D" || histoDimension == "2") {
      int yBinMin, yBinMax;
      std::cout << "From which bin for projection of 2D histogram you want to fit. Starting from 1" << std::endl;
      std::cin >> yBinMin;
      std::cout << "Until which bin for projection of 2D histogram you want to fit. Less than " << HistogramsToPotentialFit2D.at(IDofHistogram)->GetNbinsY() << std::endl;
      std::cin >> yBinMax;

      if (yBinMin < 1 || yBinMin > HistogramsToPotentialFit2D.at(IDofHistogram)->GetNbinsY()) {
        std::cout << "Provded bin " << yBinMin << " is less than 1 or greater than number of bins equal to ";
        std::cout << HistogramsToPotentialFit2D.at(IDofHistogram)->GetNbinsY() << std::endl;
        return 0;
      } else if (yBinMax < 1 || yBinMax > HistogramsToPotentialFit2D.at(IDofHistogram)->GetNbinsY()) {
        std::cout << "Provded bin " << yBinMax << " is less than 1 or greater than number of bins equal to ";
        std::cout << HistogramsToPotentialFit2D.at(IDofHistogram)->GetNbinsY() << std::endl;
        return 0;
      }

      std::cout << "Getting projection for Y min: " << HistogramsToPotentialFit2D.at(IDofHistogram)->GetYaxis()->GetBinCenter(yBinMin) << std::endl;
      std::cout << "Getting projection for Y max: " << HistogramsToPotentialFit2D.at(IDofHistogram)->GetYaxis()->GetBinCenter(yBinMax) << std::endl;
      histoToFit = (TH1D*)HistogramsToPotentialFit2D.at(IDofHistogram)->ProjectionX("projX",yBinMin,yBinMax);
      histoToFit->SetDirectory(0);
      histoName = HistogramsToPotentialFit2D.at(IDofHistogram)->GetName();
      TString tempAdd = "for_bin_";
      TString binAdd = NumberToChar(yBinMin, 0) + " " + NumberToChar(yBinMax, 0);
      histoName = histoName + tempAdd + binAdd;
      if (histoToFit->GetEntries() == 0) {
        std::cout << "Histogram " << HistogramsToPotentialFit2D.at(IDofHistogram)->GetName() << " for bins " << yBinMin << " " << yBinMax << " have no entries" << std::endl;
        return 0;
      }
    } else {
      std::cout << "Wrong option given. Type 1D or 2D " << std::endl;
      return 0;
    }
  }
  inputFile->Close();

//  PlotDerivatives(histoToFit, nameOfInput);
//  EstimateBGfromMinima(histoToFit, nameOfInput);
//  EstimateBGfromDerivative(histoToFit, nameOfInput);
//  FitSpectrum(histoToFit, nameOfInput);
  FitLineLocally(histoToFit, nameOfInput, TypeOfIntensityEstimation::Gauss);

  return 0;
}

void PlotDerivatives(TH1D* histoToFit, std::string nameOfInput)
{
  std::vector<double> weightsForSmooth = ComputeWeights(sizeOfHalfWindow, 0, orderOfPolynomialToInterpolate, 0);
  // 2nd arg - central point, 4th - order of interpolation funtion, where 0 is smooth and 1 is derivative
  std::vector<double> weightsForDerivative = ComputeWeights(sizeOfHalfWindow, 0, orderOfPolynomialToInterpolate, 1);
  std::vector<double> weightsFor2ndDerivative = ComputeWeights(sizeOfHalfWindow, 0, orderOfPolynomialToInterpolate, 2);

  TH1D* smoothHisto = (TH1D*)histoToFit->Clone("smoothHisto");
  TH1D* derivHisto = (TH1D*)histoToFit->Clone("derivHisto");
  TH1D* secDerivHisto = (TH1D*)histoToFit->Clone("secDerivHisto");
  int endBin = histoToFit->GetNbinsX();
  double sizeOfBin = 1.;
  //  if (endBin > 3)
  //    sizeOfBin = derivHisto->GetXaxis()->GetBinCenter(2) - derivHisto->GetXaxis()->GetBinCenter(1);

  for (unsigned i=sizeOfHalfWindow+2; i<endBin-sizeOfHalfWindow-2; i++) {
    double weight = 0;
    for (int j=-sizeOfHalfWindow; j<=sizeOfHalfWindow; j++) {
      weight += weightsForDerivative.at(j+sizeOfHalfWindow)*histoToFit->GetBinContent(i+j);
    }
    derivHisto->SetBinContent(i, weight/sizeOfBin);

    weight = 0;
    for (int j=-sizeOfHalfWindow; j<=sizeOfHalfWindow; j++) {
      weight += weightsFor2ndDerivative.at(j+sizeOfHalfWindow)*histoToFit->GetBinContent(i+j);
    }
    secDerivHisto->SetBinContent(i, weight/sizeOfBin);

    weight = 0;
    for (int j=-sizeOfHalfWindow; j<=sizeOfHalfWindow; j++) {
      weight += weightsForSmooth.at(j+sizeOfHalfWindow)*histoToFit->GetBinContent(i+j);
    }
    smoothHisto->SetBinContent(i, weight/sizeOfBin);
  }
  for (unsigned i=0; i<sizeOfHalfWindow+1+1; i++) {
    derivHisto->SetBinContent(1+i, 0);
    derivHisto->SetBinContent(endBin-1-i, 0);
    secDerivHisto->SetBinContent(1+i, 0);
    secDerivHisto->SetBinContent(endBin-1-i, 0);
  }

  TString outputName = "Derivatives_" + nameOfInput;
  TFile* outputFile = new TFile(outputName, "RECREATE");

  histoToFit->Write("RawEnergy");
  smoothHisto->Write("SmoothEnergy");
  derivHisto->Write("Derivative");
  secDerivHisto->Write("2ndDerivative");

  outputFile->Close();
}

void EstimateBGfromMinima(TH1D* histoToFit, std::string nameOfInput)
{
  double minRange = 0.17;
  double maxRange = 9;

  double restrictionStartRange = 2.2;

  int maxSize = 100;

  std::vector<double> weightsForSmooth = ComputeWeights(sizeOfHalfWindow, 0, orderOfPolynomialToInterpolate, 0);
  TH1D* smoothHisto = (TH1D*)histoToFit->Clone("smoothHisto");
  int endBin = histoToFit->GetNbinsX();
  double sizeOfBin = 1.;
  for (unsigned i=sizeOfHalfWindow+2; i<endBin-sizeOfHalfWindow-2; i++) {
    double weight = 0;
    for (int j=-sizeOfHalfWindow; j<=sizeOfHalfWindow; j++) {
      weight += weightsForSmooth.at(j+sizeOfHalfWindow)*histoToFit->GetBinContent(i+j);
    }
    smoothHisto->SetBinContent(i, weight/sizeOfBin);
  }

  Double_t args[maxSize], values[maxSize];
  int iterator = 0;
  for (unsigned i=1; i<smoothHisto->GetNbinsX() && iterator<maxSize; i++) {
    double binCenter = smoothHisto->GetBinCenter(i);
    if (binCenter >= minRange && binCenter <= maxRange) {
      double binValue = smoothHisto->GetBinContent(i);
      double binValueM2 = smoothHisto->GetBinContent(i-2);
      double binValueM1 = smoothHisto->GetBinContent(i-1);
      double binValueP1 = smoothHisto->GetBinContent(i+1);
      double binValueP2 = smoothHisto->GetBinContent(i+2);

      if (binCenter <= restrictionStartRange) {
        if (binValue < binValueM2 && binValue < binValueM1 && binValue < binValueP1 && binValue < binValueP2) {
          args[iterator] = binCenter;
          values[iterator] = binValue;
          iterator++;
        }
      } else {
        double binValueM3 = smoothHisto->GetBinContent(i-3);
        double binValueP3 = smoothHisto->GetBinContent(i+3);
        if (binValue < binValueM3 && binValue < binValueM2 && binValue < binValueM1 && binValue < binValueP1 && binValue < binValueP2 && binValue < binValueP3) {
          args[iterator] = binCenter;
          values[iterator] = binValue;
          iterator++;
        }
      }
    }
  }

  if (iterator < maxSize) {
    for (unsigned i=iterator; i<maxSize; i++) {
      args[i] = 12 + i*0.1;
      values[i] = 0;
    }
  }

  TCanvas *c1 = new TCanvas("c1", "", 710, 500);
  c1->SetHighLightColor(2);
  c1->SetFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameLineWidth(2);
  c1->SetLeftMargin(0.1162571);
  c1->SetRightMargin(0.0831758);

  TGraph* graph = new TGraph(maxSize, args, values);
  graph->SetMarkerColor(kMagenta);
  graph->SetMarkerSize(1);
  graph->SetMarkerStyle(21);

  std::vector<double> backgroundGaussesOff = {0.5, 1.1, 1.5, 3.2, 5.6, 7.2};
 // std::vector<double> backgroundGaussesOff = {2, 4, 5, 6, 7, 8};
  int numOfGauss = backgroundGaussesOff.size();
  TF1* graphFit = new TF1("gFit", eneMultLinesFunction, minRange, maxRange, 4 + 3*numOfGauss);
  graphFit->SetParameter(0, 10);
  graphFit->SetParName(0, "Background constant");
  graphFit->SetParLimits(0, 0, 0.01*smoothHisto->GetMaximum());
  graphFit->SetParameter(1, smoothHisto->GetBinContent((smoothHisto->FindBin(minRange))));
  graphFit->SetParName(1, "Background amplitude");
  graphFit->SetParLimits(1, 0, 10*smoothHisto->GetBinContent((smoothHisto->FindBin(minRange))));
  graphFit->SetParameter(2, smoothHisto->GetMean());
  graphFit->SetParLimits(2, 0, 5*smoothHisto->GetMean());
  graphFit->SetParName(2, "Background mean");
/*  TF1* graphFit = new TF1("gFit", bgFunctionWithLogGauss, minRange, maxRange, 4 + 3*numOfGauss);
  graphFit->SetParameter(0, 10);
  graphFit->SetParName(0, "Background constant");
  graphFit->SetParLimits(0, 0, 0.01*smoothHisto->GetMaximum());
  graphFit->SetParameter(1, 0.5);
  graphFit->SetParName(1, "Background sigma");
  graphFit->SetParLimits(1, 0, 10*smoothHisto->GetBinContent((smoothHisto->FindBin(minRange))));
  graphFit->FixParameter(2, 0);
  graphFit->SetParName(2, "Background mean");*/
  graphFit->FixParameter(3, numOfGauss);
  double freedomFraction = 0.75;
  for (unsigned i=0; i<numOfGauss; i++) {
    graphFit->SetParameter(4+3*i, 10);
    graphFit->SetParName(4+3*i, "Intensity");
    graphFit->SetParLimits(4+3*i, 0, std::pow(50,3));

    graphFit->SetParameter(5+3*i, 4);
    graphFit->SetParName(5+3*i, "Sigma");
    graphFit->SetParLimits(5+3*i, 0.5, 20);

    double offset = backgroundGaussesOff.at(i);
    graphFit->SetParameter(6+3*i, offset);
    graphFit->SetParName(6+3*i, "Offset");
    graphFit->SetParLimits(6+3*i, (1-freedomFraction)*offset, (1+freedomFraction)*offset);
  }
  graph->Fit(graphFit, "RM");

  smoothHisto->Draw();
  graph->Draw("Psame");

  TString outputName = "BGEst_" + nameOfInput;
  TFile* outputFile = new TFile(outputName, "RECREATE");

  histoToFit->Write("Raw");
  smoothHisto->Write("Smooth");
  graph->Write("BG");
  c1->Write("All");

  outputFile->Close();
}

void EstimateBGfromDerivative(TH1D* histoToFit, std::string nameOfInput)
{
  double minRange = 0.17;
  double maxRange = 9;

  double restrictionStartRange = 2.2;

  int maxSize = 50;

  std::vector<double> weightsForDerivative = ComputeWeights(sizeOfHalfWindow, 0, orderOfPolynomialToInterpolate, 1);
  TH1D* derivHisto = (TH1D*)histoToFit->Clone("derivHisto");
  int endBin = histoToFit->GetNbinsX();
  double sizeOfBin = 1.;
  for (unsigned i=sizeOfHalfWindow+2; i<endBin-sizeOfHalfWindow-2; i++) {
    double weight = 0;
    for (int j=-sizeOfHalfWindow; j<=sizeOfHalfWindow; j++) {
      weight += weightsForDerivative.at(j+sizeOfHalfWindow)*histoToFit->GetBinContent(i+j);
    }
    derivHisto->SetBinContent(i, weight/sizeOfBin);
  }

  Double_t args[maxSize], values[maxSize];
  int iterator = 0;
  for (unsigned i=1; i<derivHisto->GetNbinsX() && iterator<maxSize; i++) {
    double binCenter = derivHisto->GetBinCenter(i);
    if (binCenter >= minRange && binCenter <= maxRange) {
      double binValueM3 = derivHisto->GetBinContent(i-3);
      double binValueM2 = derivHisto->GetBinContent(i-2);
      double binValueM1 = derivHisto->GetBinContent(i-1);
      double binValue = derivHisto->GetBinContent(i);
      double binValueP1 = derivHisto->GetBinContent(i+1);
      double binValueP2 = derivHisto->GetBinContent(i+2);
      double binValueP3 = derivHisto->GetBinContent(i+3);

      if (binValueM3 > binValueM2 && binValueM2 > binValueM1 && binValueM1 > binValue && binValue > binValueP1 && binValueP1 > binValueP2 && binValueP2 > binValueP3) {
        std::pair<int, int> extremaBins = FindLocalEx(derivHisto, i);
        std::cout << derivHisto->GetBinCenter(extremaBins.first) << " " << derivHisto->GetBinCenter(extremaBins.second) << std::endl;
        i += 5;
      }
    }
  }


}

void FitSpectrum(TH1D* histoToFit, std::string nameOfInput)
{
  std::cout << "Fitting background separately " << (typeOfBackground == 1 ? "exponential type " : "power type ") << std::endl;
  std::vector<double> backgroundGaussesOff = {0.5, 1.1, 1.5, 3.2, 5.6, 7.2};
  std::vector<double> predefGaussOff = {0.512, 0.695, 0.785, 0.846, 0.93, 1.17, 1.24, 1.61, 1.95, 2.12, 2.23, 2.31, 2.435, 2.738, 4.44, 5.11, 6.13, 6.63, 7.12, 7.65, 7.79, 8.2, 8.58, 10.8};
  int noOfPredefinedGauss = predefGaussOff.size(), noOfRandomGauss = backgroundGaussesOff.size();
  unsigned noOfParameters = 3 + 1 + 3*(noOfRandomGauss + noOfPredefinedGauss);
  Double_t startFit = 0.45, endFit = 11;
  int numberOfPointsForDrawing = (endFit - startFit)*500;
/*
Parameters
0 -2 -> BG as exp/poly2 + const
3 -> noOfRandomGausses
the rest are intensity + sigma + offset for each gaussian
*/
  TF1* fitFunc = new TF1("fitFunc", eneMultLinesFunction, startFit, endFit, noOfParameters);
  fitFunc -> SetNpx(numberOfPointsForDrawing);

  double min = histoToFit->GetMinimum();
  double max = histoToFit->GetMaximum();
  if (typeOfBackground == 1) {
    fitFunc->SetParameter(0, 10);
    fitFunc->SetParName(0, "Background constant");
    fitFunc->SetParLimits(0, 0, 0.01*max);
    fitFunc->SetParameter(1, histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
    fitFunc->SetParName(1, "Background amplitude");
    fitFunc->SetParLimits(1, 0, 10*histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
    fitFunc->SetParameter(2, histoToFit->GetMean());
    fitFunc->SetParLimits(2, 0, 5*histoToFit->GetMean());
    fitFunc->SetParName(2, "Background mean");
  } else if (typeOfBackground == 2) {
    double amplitudeLimit = 10*max;
    fitFunc->SetParameter(0, histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
    fitFunc->SetParName(0, "Background amplitude");
    fitFunc->SetParLimits(0, 0, amplitudeLimit);
    fitFunc->SetParameter(1, -1.2);
    fitFunc->SetParName(1, "Background power");
    fitFunc->SetParLimits(1, -10, 0);
    fitFunc->SetParameter(2, 1);
    fitFunc->SetParName(2, "Background constant");
    fitFunc->SetParLimits(2, 0, amplitudeLimit);
  } else {
    double amplitudeLimit = (max - min)/4E8;
    fitFunc->SetParameter(0, 0.5*amplitudeLimit);
    fitFunc->SetParName(0, "Background amplitude");
    fitFunc->SetParLimits(0, 0, amplitudeLimit);
    fitFunc->SetParameter(1, 0.5*(endFit - startFit));
    fitFunc->SetParName(1, "Background minimum horizontal shift");
    fitFunc->SetParLimits(1, 0, endFit);
    double minLimit = histoToFit->GetBinContent((histoToFit->FindBin(endFit)));
    fitFunc->SetParameter(2, minLimit);
    fitFunc->SetParName(2, "Background minimum vertical shift");
    fitFunc->SetParLimits(2, 0, 1.2*minLimit);
  }
  fitFunc->FixParameter(3, noOfRandomGauss + noOfPredefinedGauss);
  fitFunc->SetParName(3, "Number of Gaussians");

  double veryLargeIntensity = max*10;
  double freedomFraction = 0.01;
  for (unsigned i=0; i<noOfPredefinedGauss; i++) {
    Double_t offset = predefGaussOff.at(i);
    Double_t intens = 0.001*veryLargeIntensity/std::pow(offset, 2);
    Double_t sigma = 0.02*offset;

    fitFunc->SetParameter(4+3*i, intens);
    fitFunc->SetParName(4+3*i, "Intensity");
    fitFunc->SetParLimits(4+3*i, 0, veryLargeIntensity);

    fitFunc->SetParameter(5+3*i, sigma);
    fitFunc->SetParName(5+3*i, "Sigma");
    fitFunc->SetParLimits(5+3*i, 0.5*sigma, 1.5*sigma);

    fitFunc->SetParameter(6+3*i, offset);
    fitFunc->SetParName(6+3*i, "Offset");
    fitFunc->SetParLimits(6+3*i, (1-freedomFraction)*offset, (1+freedomFraction)*offset);
  }
  for (unsigned i=noOfPredefinedGauss; i<noOfRandomGauss + noOfPredefinedGauss; i++) {
    Double_t offset = backgroundGaussesOff.at(i-noOfPredefinedGauss);//startFit + (i+1 - noOfPredefinedGauss)*(endFit-startFit)/noOfRandomGauss;

    fitFunc->SetParameter(4+3*i, 50);
    fitFunc->SetParName(4+3*i, "Intensity");
    fitFunc->SetParLimits(4+3*i, 0, veryLargeIntensity);

    fitFunc->SetParameter(5+3*i, 0.5);
    fitFunc->SetParName(5+3*i, "Sigma");
    fitFunc->SetParLimits(5+3*i, 0.2, 5);

    fitFunc->SetParameter(6+3*i, offset);
    fitFunc->SetParName(6+3*i, "Offset");
    fitFunc->SetParLimits(6+3*i, (1-20*freedomFraction)*offset, (1+20*freedomFraction)*offset);
  }

  TGraphAsymmErrors *graphToFit = new TGraphAsymmErrors(histoToFit);
  Int_t np = graphToFit->GetN();
  for (Int_t i=0; i<np; i++) {
    graphToFit->SetPointEXhigh(i,0.);
    graphToFit->SetPointEXlow(i,0.);
  }
  graphToFit->Fit(fitFunc,"RM");

  const Int_t n = noOfPredefinedGauss;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  std::cout << "Energy  Intensity        Error" << std::endl;
  int it = 0;
  for (unsigned i=0; i<noOfPredefinedGauss; i++) {
    std::cout << NumberToChar(fitFunc->GetParameter(6+3*i), 3) << "\t" << NumberToChar(fitFunc->GetParameter(4+3*i), 3) << "\t" << NumberToChar(fitFunc->GetParError(4+3*i), 3) << std::endl;
    x[i] = fitFunc->GetParameter(6+3*i);
    ex[i] = 0.;
    y[i] = fitFunc->GetParameter(4+3*i);
    ey[i] = fitFunc->GetParError(4+3*i);
  }

  auto gr = new TGraphErrors(n,x,y,ex,ey);
  gr->SetTitle("TGraphErrors");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);

  TF1* fitFuncClone = new TF1("fitFuncClone", eneMultLinesFunction, startFit, endFit, noOfParameters);
  fitFuncClone->SetParameters(fitFunc->GetParameters());
  for (unsigned i=0; i<noOfPredefinedGauss; i++) {
    fitFuncClone->SetParameter(4+3*i, 0);
  }

  TF1* backgroundClone = new TF1("backgroundClone", eneMultLinesFunction, startFit, endFit, noOfParameters);
  backgroundClone->SetParameters(fitFunc->GetParameters());
  TF1* backgroundClone2 = new TF1("backgroundClone2", eneMultLinesFunction, startFit, endFit, noOfParameters);
  backgroundClone2->SetParameters(fitFunc->GetParameters());
  backgroundClone2->SetParameter(0, 0);
  backgroundClone2->SetParameter(1, 0);
  for (unsigned i=0; i<noOfPredefinedGauss; i++) {
    backgroundClone->SetParameter(4+3*i, 0);
    backgroundClone2->SetParameter(4+3*i, 0);
  }
  for (unsigned i=noOfPredefinedGauss; i<noOfRandomGauss + noOfPredefinedGauss; i++) {
    backgroundClone->SetParameter(4+3*i, 0);
  }

  TCanvas *c1 = new TCanvas("c1", "", 710, 500);
  c1->SetHighLightColor(2);
  c1->SetFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameLineWidth(2);
  c1->SetLeftMargin(0.1162571);
  c1->SetRightMargin(0.0831758);
  graphToFit->Draw();
  //fitFunc->Draw("same");
  fitFuncClone->SetLineColor(kGreen);
  fitFuncClone->Draw("same");
  backgroundClone->SetLineColor(kBlack);
  backgroundClone->Draw("same");
  backgroundClone2->SetLineColor(kCyan);
  backgroundClone2->Draw("same");

  std::cout << "Fitted background parameters " << (typeOfBackground == 1 ? "exponential type " : "power type ");
  std::cout << "and parameters y0, Amplitude and Mean [MeV] ";
  std::cout << fitFunc->GetParameter(0) << " " << fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2) << std::endl;

  TString add = "Exp_";
  if (typeOfBackground == 2)
    add = "Pow_";
  else if (typeOfBackground == 3)
    add = "Lin_";

  TString outputName = "Fitted_" + add + nameOfInput;
  TFile* outputFile = new TFile(outputName, "RECREATE");

  histoToFit->Write("Raw");
  graphToFit->Write("Fitted");
  c1->Write("All");
  gr->Write("ResTogether");

  outputFile->Close();
}

std::pair<double, double> FitSingleLine(TH1D* histoToFit, double energy, TypeOfIntensityEstimation type)
{
  double a = 2.0*pow(10, -4); // in MeV
  double b = 2.22*pow(10, -2);
  double c = 0.5;
  double sigma = (a + b*sqrt(energy + c*pow(energy, 2)))/(2.35482004503);
  double sigmaFracExtend = 5;
  Double_t startFit = energy - sigmaFracExtend*sigma, endFit = energy + sigmaFracExtend*sigma;
  if (type == TypeOfIntensityEstimation::Gauss) {
    int numberOfPointsForDrawing = (endFit - startFit)*500;
    TF1* fitFunc = new TF1("fitFunc", eneLineFunction, startFit, endFit, 6);
    fitFunc -> SetNpx(numberOfPointsForDrawing);
    double min = histoToFit->GetMinimum();
    double max = histoToFit->GetMaximum();
    if (typeOfBackground == 1) {
      fitFunc->SetParameter(0, 0.5*histoToFit->GetBinContent((histoToFit->FindBin(endFit))));
      fitFunc->SetParName(0, "Background constant");
      fitFunc->SetParLimits(0, 0, max);
      double mean = (histoToFit->GetBinCenter((histoToFit->FindBin(endFit))) - histoToFit->GetBinCenter((histoToFit->FindBin(startFit))));
      mean = mean*log(histoToFit->GetBinContent((histoToFit->FindBin(startFit))) / histoToFit->GetBinContent((histoToFit->FindBin(endFit))));
      double amp = histoToFit->GetBinContent((histoToFit->FindBin(startFit)));
      fitFunc->SetParameter(1, amp);
      fitFunc->SetParName(1, "Background amplitude");
      fitFunc->SetParLimits(1, 0, 10*histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
      fitFunc->SetParameter(2, mean);
      fitFunc->SetParLimits(2, 0, 5*histoToFit->GetMean());
      fitFunc->SetParName(2, "Background mean");
    } else if (typeOfBackground == 2) {
      double amplitudeLimit = 10*max;
      fitFunc->SetParameter(0, histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
      fitFunc->SetParName(0, "Background amplitude");
      fitFunc->SetParLimits(0, 0, amplitudeLimit);
      fitFunc->SetParameter(1, -1.2);
      fitFunc->SetParName(1, "Background power");
      fitFunc->SetParLimits(1, -10, 0);
      fitFunc->SetParameter(2, 1);
      fitFunc->SetParName(2, "Background constant");
      fitFunc->SetParLimits(2, 0, amplitudeLimit);
    } else {
      double amplitudeLimit = 10*max;
      double slope = (histoToFit->GetBinContent((histoToFit->FindBin(startFit))) / histoToFit->GetBinContent((histoToFit->FindBin(endFit))));
      slope = slope/(histoToFit->GetBinCenter((histoToFit->FindBin(startFit))) / histoToFit->GetBinCenter((histoToFit->FindBin(endFit))));
      double constant = histoToFit->GetBinContent((histoToFit->FindBin(startFit))) + slope*histoToFit->GetBinCenter((histoToFit->FindBin(startFit)));
      fitFunc->SetParameter(0, constant);
      fitFunc->SetParName(0, "Background constant");
      fitFunc->SetParLimits(0, histoToFit->GetBinContent((histoToFit->FindBin(startFit))), amplitudeLimit);
      fitFunc->SetParameter(1, slope);
      fitFunc->SetParName(1, "Background slope");
      fitFunc->SetParLimits(1, 0, 2000);
      fitFunc->SetParameter(2, 1);
      fitFunc->SetParName(2, "Empty parameter");
    }
    double veryLargeIntensity = max*10;
    double freedomFraction = 0.025;

    Double_t offset = energy;
    Double_t intens = 0.001*veryLargeIntensity/std::pow(offset, 2);

    fitFunc->SetParameter(3, intens);
    fitFunc->SetParName(3, "Intensity");
    fitFunc->SetParLimits(3, 0, veryLargeIntensity);

    fitFunc->SetParameter(4, sigma);
    fitFunc->SetParName(4, "Sigma");
    fitFunc->SetParLimits(4, 0.5*sigma, 2*sigma);

    fitFunc->SetParameter(5, offset);
    fitFunc->SetParName(5, "Offset");
    fitFunc->SetParLimits(5, (1-freedomFraction)*offset, (1+freedomFraction)*offset);

    TGraphAsymmErrors *graphToFit = new TGraphAsymmErrors(histoToFit);
    Int_t np = graphToFit->GetN();
    for (Int_t i=0; i<np; i++) {
      graphToFit->SetPointEXhigh(i,0.);
      graphToFit->SetPointEXlow(i,0.);
    }
    graphToFit->Fit(fitFunc,"RM");

    fitFunc->Draw("same");

    TF1* fitFuncExp = new TF1("fitFuncExp", eneLineFunction, startFit, endFit, 6);
    fitFuncExp->SetParameters(fitFunc->GetParameters());
    fitFuncExp->SetParameter(3,0);
    fitFuncExp->SetLineColor(kGreen);
    fitFuncExp->Draw("same");

    std::cout << "Fitted background parameters " << (typeOfBackground == 1 ? "exponential type " : "other type ");
    std::cout << "and parameters y0, Amplitude and Mean [MeV] " << std::endl;
    std::cout << fitFunc->GetParameter(0) << " " << fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2) << std::endl;
    std::cout << "Gauss parameters [Intensity] [Sigma] [Offset]" << std::endl;
    std::cout << fitFunc->GetParameter(3) << " " << fitFunc->GetParameter(4) << " " << fitFunc->GetParameter(5) << std::endl;

    std::pair<double, double> res = {fitFunc->GetParameter(3), fitFunc->GetParError(3)};
    return res;
  } else {
    Int_t minBin = histoToFit->FindBin(startFit);
    Int_t maxBin = histoToFit->FindBin(endFit);

    double sum = 0, error = 0;
    for (unsigned i=minBin; i<=maxBin; i++) {
      sum += histoToFit->GetBinContent(i);
      error += histoToFit->GetBinError(i);
    }
    std::cout << "Sum: " << sum << " Err " << error << std::endl;

    Double_t binSize = histoToFit->GetBinCenter(minBin+1) - histoToFit->GetBinCenter(minBin);
    std::cout << "Integral: " << binSize*histoToFit->Integral(minBin, maxBin) << std::endl;
    Double_t startCenter = histoToFit->GetBinCenter(minBin), endCenter = histoToFit->GetBinCenter(maxBin);
    Double_t startCount = histoToFit->GetBinContent(minBin), endCount = histoToFit->GetBinContent(maxBin);
    double aParam = (endCount - startCount)/(endCenter - startCenter);
    double bParam = startCount - startCenter*aParam;

    TF1* fitFunc = new TF1("fitFunc", "[0] + [1]*x", startFit, endFit);
    fitFunc->SetParameters(bParam,aParam);
    fitFunc->Draw("same");

    double BGsum = 0.;
    for (unsigned i=minBin; i<=maxBin; i++) {
      BGsum += binSize*(bParam+aParam*histoToFit->GetBinCenter(i));
    }
    std::cout << std::endl;
    std::cout << "Lin Background: " << BGsum << std::endl;
    double res1 = binSize*histoToFit->Integral(minBin, maxBin) - BGsum;
    std::cout << "Diff: " << res1 << std::endl;
    double res2 = 0; // to do: calculate error

    std::pair<double, double> res = {res1, res2};
    return res;
  }
}

std::vector<std::pair<double, double>> FitMultipleLines(TH1D* histoToFit, std::vector<double> energies, TypeOfIntensityEstimation type)
{
  double a = 2.0*pow(10, -4); // in MeV
  double b = 2.22*pow(10, -2);
  double c = 0.5;
  std::vector<double> sigmas;
  int numOfLines = energies.size();
  for (unsigned i=0; i<numOfLines; i++) {
    sigmas.push_back((a + b*sqrt(energies.at(i) + c*pow(energies.at(i), 2)))/(2.35482004503));
  }
  double sigmaFracExtend = 5;
  Double_t startFit = energies.at(0) - sigmaFracExtend*sigmas.at(0);
  Double_t endFit = energies.at(numOfLines-1) + sigmaFracExtend*sigmas.at(numOfLines-1);
  if (type == TypeOfIntensityEstimation::Gauss) {
    int numberOfPointsForDrawing = (endFit - startFit)*500;
    TF1* fitFunc = new TF1("fitFunc", eneMultLinesFunction, startFit, endFit, 3+1+3*numOfLines);
    fitFunc -> SetNpx(numberOfPointsForDrawing);
    double min = histoToFit->GetMinimum();
    double max = histoToFit->GetMaximum();
    if (typeOfBackground == 1) {
      fitFunc->SetParameter(0, 0.5*histoToFit->GetBinContent((histoToFit->FindBin(endFit))));
      fitFunc->SetParName(0, "Background constant");
      fitFunc->SetParLimits(0, 0, max);
      double mean = (histoToFit->GetBinCenter((histoToFit->FindBin(endFit))) - histoToFit->GetBinCenter((histoToFit->FindBin(startFit))));
      mean = mean*log(histoToFit->GetBinContent((histoToFit->FindBin(startFit))) / histoToFit->GetBinContent((histoToFit->FindBin(endFit))));
      double amp = histoToFit->GetBinContent((histoToFit->FindBin(startFit)));
      fitFunc->SetParameter(1, amp);
      fitFunc->SetParName(1, "Background amplitude");
      fitFunc->SetParLimits(1, 0, 10*histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
      fitFunc->SetParameter(2, mean);
      fitFunc->SetParLimits(2, 0, 5*histoToFit->GetMean());
      fitFunc->SetParName(2, "Background mean");
    } else if (typeOfBackground == 2) {
      double amplitudeLimit = 10*max;
      fitFunc->SetParameter(0, histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
      fitFunc->SetParName(0, "Background amplitude");
      fitFunc->SetParLimits(0, 0, amplitudeLimit);
      fitFunc->SetParameter(1, -1.2);
      fitFunc->SetParName(1, "Background power");
      fitFunc->SetParLimits(1, -10, 0);
      fitFunc->SetParameter(2, 1);
      fitFunc->SetParName(2, "Background constant");
      fitFunc->SetParLimits(2, 0, amplitudeLimit);
    } else {
      double amplitudeLimit = 10*histoToFit->GetMaximum();
      double slope = (histoToFit->GetBinContent((histoToFit->FindBin(startFit))) / histoToFit->GetBinContent((histoToFit->FindBin(endFit))));
      slope = slope/(histoToFit->GetBinCenter((histoToFit->FindBin(startFit))) / histoToFit->GetBinCenter((histoToFit->FindBin(endFit))));
      double constant = histoToFit->GetBinContent((histoToFit->FindBin(startFit))) + slope*histoToFit->GetBinCenter((histoToFit->FindBin(startFit)));
      fitFunc->SetParameter(0, constant);
      fitFunc->SetParName(0, "Background constant");
      fitFunc->SetParLimits(0, histoToFit->GetBinContent((histoToFit->FindBin(startFit))), amplitudeLimit);
      fitFunc->SetParameter(1, slope);
      fitFunc->SetParName(1, "Background slope");
      fitFunc->SetParLimits(1, 0, 2000);
      fitFunc->SetParameter(2, 1);
      fitFunc->SetParName(2, "Empty parameter");
    }
    fitFunc->FixParameter(3, numOfLines);

    double veryLargeIntensity = histoToFit->GetMaximum()*10;
    double freedomFraction = 0.025;

    for (unsigned i=0; i<numOfLines; i++) {
      Double_t offset = energies.at(i);
      Double_t intens = 0.001*veryLargeIntensity/std::pow(offset, 2);

      fitFunc->SetParameter(4+3*i, intens);
      fitFunc->SetParName(4+3*i, "Intensity");
      fitFunc->SetParLimits(4+3*i, 0, veryLargeIntensity);

      fitFunc->SetParameter(5+3*i, sigmas.at(i));
      fitFunc->SetParName(5+3*i, "Sigma");
      fitFunc->SetParLimits(5+3*i, 0.5*sigmas.at(i), 2*sigmas.at(i));

      fitFunc->SetParameter(6+3*i, energies.at(i));
      fitFunc->SetParName(6+3*i, "Offset");
      fitFunc->SetParLimits(6+3*i, (1-freedomFraction)*energies.at(i), (1+freedomFraction)*energies.at(i));
    }

    TGraphAsymmErrors *graphToFit = new TGraphAsymmErrors(histoToFit);
    Int_t np = graphToFit->GetN();
    for (Int_t i=0; i<np; i++) {
      graphToFit->SetPointEXhigh(i,0.);
      graphToFit->SetPointEXlow(i,0.);
    }
    graphToFit->Fit(fitFunc,"RM");

    fitFunc->Draw("same");

    TF1* fitFuncExp = new TF1("fitFuncExp", eneLineFunction, startFit, endFit, 6);
    fitFuncExp->SetParameter(0, fitFunc->GetParameter(0));
    fitFuncExp->SetParameter(1, fitFunc->GetParameter(1));
    fitFuncExp->SetParameter(2, fitFunc->GetParameter(2));
    fitFuncExp->SetParameter(3, 0);
    fitFuncExp->SetParameter(4, 1);
    fitFuncExp->SetParameter(5, 1);
    fitFuncExp->SetLineColor(kGreen);
    fitFuncExp->Draw("same");

    std::cout << "Fitted background parameters " << (typeOfBackground == 1 ? "exponential type " : "other type ");
    std::cout << "and parameters y0, Amplitude and Mean [MeV] " << std::endl;
    std::cout << fitFunc->GetParameter(0) << " " << fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2) << std::endl;
    std::vector<std::pair<double, double>> res;
    for (unsigned i=0; i<numOfLines; i++) {
      std::cout << "Gauss parameters [Intensity] [Sigma] [Offset]" << std::endl;
      std::cout << fitFunc->GetParameter(4+3*i) << " " << fitFunc->GetParameter(5+3*i) << " " << fitFunc->GetParameter(6+3*i) << std::endl;
      std::pair<double, double> temp = {fitFunc->GetParameter(4+3*i), fitFunc->GetParError(4+3*i)};
      res.emplace_back(temp);
    }

    return res;
  } else {
    Int_t minBin = histoToFit->FindBin(startFit);
    Int_t maxBin = histoToFit->FindBin(endFit);

    double sum = 0, error = 0;
    for (unsigned i=minBin; i<=maxBin; i++) {
      sum += histoToFit->GetBinContent(i);
      error += histoToFit->GetBinError(i);
    }
    std::cout << "Sum: " << sum << " Err " << error << std::endl;

    Double_t binSize = histoToFit->GetBinCenter(minBin+1) - histoToFit->GetBinCenter(minBin);
    std::cout << "Integral: " << binSize*histoToFit->Integral(minBin, maxBin) << std::endl;
    Double_t startCenter = histoToFit->GetBinCenter(minBin), endCenter = histoToFit->GetBinCenter(maxBin);
    Double_t startCount = histoToFit->GetBinContent(minBin), endCount = histoToFit->GetBinContent(maxBin);
    double aParam = (endCount - startCount)/(endCenter - startCenter);
    double bParam = startCount - startCenter*aParam;

    TF1* fitFunc = new TF1("fitFunc", "[0] + [1]*x", startFit, endFit);
    fitFunc->SetParameters(bParam,aParam);
    fitFunc->Draw("same");

    double BGsum = 0.;
    for (unsigned i=minBin; i<=maxBin; i++) {
      BGsum += binSize*(bParam+aParam*histoToFit->GetBinCenter(i));
    }
    std::cout << std::endl;
    std::cout << "Lin Background: " << BGsum << std::endl;
    double res1 = binSize*histoToFit->Integral(minBin, maxBin) - BGsum;
    std::cout << "Diff: " << res1 << std::endl;
    double res2 = 0; // to do: calculate error

    std::vector<std::pair<double, double>> res;
    std::pair<double, double> temp = {res1, res2};
    res.emplace_back(temp);
    return res;
  }
}

void FitLineLocally(TH1D* histoToFit, std::string nameOfInput, TypeOfIntensityEstimation type)
{
  TCanvas *c1 = new TCanvas("c1", "", 710, 500);
  c1->SetHighLightColor(2);
  c1->SetFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameLineWidth(2);
  c1->SetLeftMargin(0.1162571);
  c1->SetRightMargin(0.0831758);

  histoToFit->Draw();

  std::vector<std::pair<double, double>> results;
  for (unsigned i=0; i<desiredLines.size(); i++) {
    results.emplace_back(FitSingleLine(histoToFit, desiredLines.at(i), type));
  }

  std::vector<std::pair<double, double>> temp;
  for (unsigned i=0; i<complexLines.size(); i++) {
    temp = FitMultipleLines(histoToFit, complexLines.at(i), type);
    for (unsigned j=0; j<temp.size(); j++) {
      desiredLines.push_back(complexLines.at(i).at(j));
      results.emplace_back(temp.at(j));
    }
  }

  const Int_t n = desiredLines.size();
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  std::cout << "Energy  Intensity        Error" << std::endl;
  int it = 0;
  for (auto iter : results) {
    std::cout << NumberToChar(desiredLines.at(it), 3) << "\t" << NumberToChar(iter.first, 6) << "\t" << NumberToChar(iter.second, 6) << std::endl;
    x[it] = desiredLines.at(it);
    ex[it] = 0.;
    y[it] = iter.first;
    ey[it] = iter.second;
    it++;
  }

  auto gr = new TGraphErrors(n,x,y,ex,ey);
  gr->SetTitle("TGraphErrors");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);

  TString outputName;
  if (type == TypeOfIntensityEstimation::Gauss) {
    if (typeOfBackground == 1)
      outputName = "FitLocGaussExp_" + nameOfInput;
    else if (typeOfBackground == 2)
      outputName = "FitLocGaussPow_" + nameOfInput;
    else
      outputName = "FitLocGaussLin_" + nameOfInput;
  } else
    outputName = "FitLocInteg_" + nameOfInput;

  TFile* outputFile = new TFile(outputName, "RECREATE");
  c1->Write("Fitted");
  gr->Write("ResTogether");
  outputFile->Close();
}

Double_t eneLineFunction(Double_t *A, Double_t *P)
{
  Double_t result = 0;
  switch (typeOfBackground) {
    case 1:
      result += P[0] + P[1]*exp(-A[0]/P[2]);
      break;
    case 2:
      result += P[0]*std::pow(A[0], P[1]) + P[2];
      break;
    case 3:
      result += P[0] - P[1]*A[0];
      break;
  }

  result += P[3]*GaussDistr(A[0], P[5], P[4]);

  return result;
}

Double_t eneMultLinesFunction(Double_t *A, Double_t *P)
{
  Double_t result = 0;
  switch (typeOfBackground) {
    case 1:
      result += P[0] + P[1]*exp(-A[0]/P[2]);
      break;
    case 2:
      result += P[0]*std::pow(A[0], P[1]) + P[2];
      break;
    case 3:
      result += P[0] - P[1]*A[0];
      break;
  }

  Int_t noOfRandomGauss = P[3];
  for (unsigned i=0; i<noOfRandomGauss; i++) {
    result += P[4+3*i]*GaussDistr(A[0], P[6+3*i], P[5+3*i]);
  }

  return result;
}

Double_t bgFunctionWithLogGauss(Double_t *A, Double_t *P)
{
  Double_t result = 0;

  result += P[0]*LogGaussDistr(A[0], P[2], P[1]);

  Int_t noOfRandomGauss = P[3];
  for (unsigned i=0; i<noOfRandomGauss; i++) {
    result += P[4+3*i]*LogGaussDistr(P[6+3*i]-
    A[0], P[6+3*i], P[5+3*i]);
  }

  return result;
}

double GaussDistr(double x, double mean, double sigma)
{
  return 1/sqrt(2*M_PI)/sigma*exp(-0.5*std::pow((x-mean)/sigma, 2));
}

double LogGaussDistr(double x, double mean, double sigma)
{
  if (x > 0)
    return 1/x/sqrt(2*M_PI)/sigma*exp(-0.5*std::pow((log(x)-mean)/sigma, 2));
  else
    return 0;
}

std::pair<int, int> FindLocalEx(TH1D* derivHisto, int centerBin)
{
  int binMax = 0, binMin = 0;
  int rangeForEx = 3;
  for (unsigned i=0; i<centerBin-rangeForEx && i<derivHisto->GetNbinsX()-centerBin-rangeForEx; i++) {
    double binValuePlus = derivHisto->GetBinContent(centerBin + i);
    double binValueMinus = derivHisto->GetBinContent(centerBin - i);

    if (binMax == 0) {
      int test = 0;
      for (unsigned j=0; j<rangeForEx; j++) {
        if (derivHisto->GetBinContent(centerBin - i - j - 1) < binValueMinus)
          test++;
      }
      if (test == rangeForEx)
        binMax = centerBin - i;
    }
    if (binMin == 0) {
      int test = 0;
      for (unsigned j=0; j<rangeForEx; j++) {
        if (derivHisto->GetBinContent(centerBin + i + j + 1) > binValuePlus)
          test++;
      }
      if (test == rangeForEx)
        binMin = centerBin + i;
    }
    if (binMax*binMin != 0)
      break;
  }

  std::pair<int, int> res = {binMax, binMin};
  return res;
}

double GramPolynomial(int i, int m, int k, int s)
{
  if (k>0) {
    return (4.*k - 2.)/(k*(2.*m - k + 1.))*(i*GramPolynomial(i, m, k-1, s) + s*GramPolynomial(i, m, k-1, s-1)) - ((k-1.)*(2.*m + k))/(k*(2.*m - k + 1.))*GramPolynomial(i, m, k-2, s);
  } else {
    if (k==0 && s==0)
      return 1.;
    else
      return 0.;
  }
}

double GenFactor(int a, int b)
{
  double gf = 1.;
  for (int j=(a-b)+1; j<=a; j++) {
    gf*=j;
  }
  return gf;
}

double CalcWeight(int i, int t, int m, int n, int s)
{
  double w = 0;
  for (int k=0; k<=n; ++k) {
    w = w + (2*k + 1)*(GenFactor(2*m, k)/GenFactor(2*m + k + 1, k + 1))*GramPolynomial(i, m, k, 0)*GramPolynomial(t, m, k, s);
  }
  return w;
}

// m -> size of half window - how many points from the central are taken
// t -> number of data point (give 0)
// n -> order of interpolation polynomial -> 1 - line, 2 - quadratic ...
// s -> order of function, 0 - smooth, 1 - derivative, 2 - second derivative, ...
std::vector<double> ComputeWeights(int m, int t, int n, int s)
{
  std::vector<double> weights(2*m + 1);
  for (int i=0; i<2*m+1; ++i) {
    weights.at(i) = CalcWeight(i-m, t, m, n, s);
  }

  std::cout << "Calculated weights ";
  for (unsigned i=0; i<weights.size(); i++)
    std::cout << weights.at(i) << " ";
  std::cout << std::endl;

  return weights;
}

bool FileCheck(const std::string& NameOfFile)
{
    struct stat buffer;
    return (stat(NameOfFile.c_str(), &buffer) == 0);
}

std::string NumberToChar(double number, int precision)
{
    std::ostringstream conv;
    conv << std::fixed << std::setprecision(precision);
    conv << number;
    return conv.str();
}
