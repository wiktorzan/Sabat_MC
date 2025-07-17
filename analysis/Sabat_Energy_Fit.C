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

int typeOfBackground = 1;

int sizeOfHalfWindow = 5;
int orderOfPolynomialToInterpolate = 2;

bool FileCheck(const std::string& NameOfFile);
std::string NumberToChar(double number, int precision);

void PlotDerivatives(TH1D* histo, std::string nameOfInput);
void FitSpectrum(TH1D* histoToFit, std::string nameOfInput);
std::string FitLine(TH1D* histoToFit, double energy);

Double_t eneLineFunction(Double_t *A, Double_t *P);
Double_t bgFunction(Double_t *A, Double_t *P);
double GaussDistr(double x, double mean, double sigma);
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
//  FitSpectrum(histoToFit, nameOfInput);
  double desiredEnergy = 1.;
  std::string test = "y";
  while (test == "y") {
    std::cout << "What energy line to estimate in [MeV]: " << std::endl;
    std::cin >> desiredEnergy;
    test = FitLine(histoToFit, desiredEnergy);
  }

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

void FitSpectrum(TH1D* histoToFit, std::string nameOfInput)
{
  std::cout << "Fitting background separately " << (typeOfBackground == 1 ? "exponential type " : "quadratic type ") << std::endl;
  std::vector<double> backgroundGaussesOff = {0.5, 1.5, 5.5, 8, 9};
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
  TF1* fitFunc = new TF1("fitFunc", bgFunction, startFit, endFit, noOfParameters);
  fitFunc -> SetNpx(numberOfPointsForDrawing);

  double min = histoToFit->GetMinimum();
  if (typeOfBackground == 1) {
    fitFunc->SetParameter(0, 10);
    fitFunc->SetParName(0, "Background constant");
    fitFunc->SetParLimits(0, 0, min);
    fitFunc->SetParameter(1, histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
    fitFunc->SetParName(1, "Background amplitude");
    fitFunc->SetParLimits(1, 0, 10*histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
    fitFunc->SetParameter(2, histoToFit->GetMean());
    fitFunc->SetParLimits(2, 0, 5*histoToFit->GetMean());
    fitFunc->SetParName(2, "Background mean");
  } else {
    double amplitudeLimit = (histoToFit->GetMaximum() - min)/4E8;
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

  double veryLargeIntensity = histoToFit->GetMaximum()*10;
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
    fitFunc->SetParLimits(5+3*i, 0.5*sigma, 5*sigma);

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
    fitFunc->SetParLimits(6+3*i, (1-10*freedomFraction)*offset, (1+10*freedomFraction)*offset);
  }

  TGraphAsymmErrors *graphToFit = new TGraphAsymmErrors(histoToFit);
  Int_t np = graphToFit->GetN();
  for (Int_t i=0; i<np; i++) {
    graphToFit->SetPointEXhigh(i,0.);
    graphToFit->SetPointEXlow(i,0.);
  }
  graphToFit->Fit(fitFunc,"RM");

  TCanvas *c1 = new TCanvas("c1", "", 710, 500);
  c1->SetHighLightColor(2);
  c1->SetFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameLineWidth(2);
  c1->SetLeftMargin(0.1162571);
  c1->SetRightMargin(0.0831758);
  graphToFit->Draw();
  fitFunc->Draw("same");

  std::cout << "Fitted background parameters " << (typeOfBackground == 1 ? "exponential type " : "quadratic type ");
  std::cout << "and parameters y0, Amplitude and Mean [MeV] ";
  std::cout << fitFunc->GetParameter(0) << " " << fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2) << std::endl;

  TString outputName = "Fitted_" + nameOfInput;
  TFile* outputFile = new TFile(outputName, "RECREATE");

  histoToFit->Write("Raw");
  graphToFit->Write("Fitted");

  outputFile->Close();
}

std::string FitLine(TH1D* histoToFit, double energy)
{
  double a = 2.0*pow(10, -4); // in MeV
  double b = 2.22*pow(10, -2);
  double c = 0.5;
  double sigma = (a + b*sqrt(energy + c*pow(energy, 2)))/(2.35482004503);
  double sigmaFracExtend = 5;
  Double_t startFit = energy - sigmaFracExtend*sigma, endFit = energy + sigmaFracExtend*sigma;
  int numberOfPointsForDrawing = (endFit - startFit)*500;
  TF1* fitFunc = new TF1("fitFunc", eneLineFunction, startFit, endFit, 6);
  fitFunc -> SetNpx(numberOfPointsForDrawing);
  double min = histoToFit->GetMinimum();
  if (typeOfBackground == 1) {
    fitFunc->SetParameter(0, 10);
    fitFunc->SetParName(0, "Background constant");
    fitFunc->SetParLimits(0, 0, min);
    fitFunc->SetParameter(1, histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
    fitFunc->SetParName(1, "Background amplitude");
    fitFunc->SetParLimits(1, 0, 10*histoToFit->GetBinContent((histoToFit->FindBin(startFit))));
    fitFunc->SetParameter(2, histoToFit->GetMean());
    fitFunc->SetParLimits(2, 0, 5*histoToFit->GetMean());
    fitFunc->SetParName(2, "Background mean");
  } else {
    double amplitudeLimit = (histoToFit->GetMaximum() - min)/4E8;
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
  double veryLargeIntensity = histoToFit->GetMaximum()*10;
  double freedomFraction = 0.05;

  Double_t offset = energy;
  Double_t intens = 0.001*veryLargeIntensity/std::pow(offset, 2);

  fitFunc->SetParameter(3, intens);
  fitFunc->SetParName(3, "Intensity");
  fitFunc->SetParLimits(3, 0, veryLargeIntensity);

  fitFunc->SetParameter(4, sigma);
  fitFunc->SetParName(4, "Sigma");
  fitFunc->SetParLimits(4, 0.5*sigma, 5*sigma);

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

  TCanvas *c1 = new TCanvas("c1", "", 710, 500);
  c1->SetHighLightColor(2);
  c1->SetFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameLineWidth(2);
  c1->SetLeftMargin(0.1162571);
  c1->SetRightMargin(0.0831758);
  graphToFit->Draw();
  fitFunc->Draw("same");

  std::cout << "Fitted background parameters " << (typeOfBackground == 1 ? "exponential type " : "quadratic type ");
  std::cout << "and parameters y0, Amplitude and Mean [MeV] " << std::endl;
  std::cout << fitFunc->GetParameter(0) << " " << fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2) << std::endl;
  std::cout << "Gauss parameters [Intensity] [Sigma] [Offset]" << std::endl;
  std::cout << fitFunc->GetParameter(4) << " " << fitFunc->GetParameter(5) << " " << fitFunc->GetParameter(6) << std::endl;

  Int_t minBin = histoToFit->FindBin(startFit);
  Int_t maxBin = histoToFit->FindBin(endFit);

  Double_t binSize = histoToFit->GetBinCenter(minBin+1) - histoToFit->GetBinCenter(minBin);
  std::cout << "Integral: " << binSize*histoToFit->Integral(minBin, maxBin) << std::endl;
  Double_t startCenter = histoToFit->GetBinCenter(minBin), endCenter = histoToFit->GetBinCenter(maxBin);
  Double_t startCount = histoToFit->GetBinContent(minBin), endCount = histoToFit->GetBinContent(maxBin);
  double aParam = (endCount - startCount)/(endCenter - startCenter);
  double bParam = startCount - startCenter*aParam;
  double BGsum = 0.;
  for (unsigned i=minBin; i<=maxBin; i++) {
    BGsum += binSize*(bParam+aParam*histoToFit->GetBinCenter(i));
  }
  std::cout << std::endl;
  std::cout << "Lin Background: " << BGsum << std::endl;
  std::cout << "Diff: " << binSize*histoToFit->Integral(minBin, maxBin) - BGsum << std::endl;

  std::string test;
  std::cout << "Next line? [y/n]" << std::endl;
  std::cin >> test;

  return test;
}

Double_t eneLineFunction(Double_t *A, Double_t *P)
{
  Double_t result = 0;
  if (typeOfBackground == 1)
    result += P[0] + P[1]*exp(-A[0]/P[2]);
  else
    result += P[0]*std::pow(A[0] - P[1], 2) + P[2];

  result += P[3]*GaussDistr(A[0], P[5], P[4]);

  return result;
}

Double_t bgFunction(Double_t *A, Double_t *P)
{
  Double_t result = 0;
  if (typeOfBackground == 1)
    result += P[0] + P[1]*exp(-A[0]/P[2]);
  else
    result += P[0]*std::pow(A[0] - P[1], 2) + P[2];

  Int_t noOfRandomGauss = P[3];
  for (unsigned i=0; i<noOfRandomGauss; i++) {
    result += P[4+3*i]*GaussDistr(A[0], P[6+3*i], P[5+3*i]);
  }

  return result;
}

double GaussDistr(double x, double mean, double sigma)
{
  return 1/sqrt(2*M_PI)/sigma*exp(-0.5*std::pow((x-mean)/sigma, 2));
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
