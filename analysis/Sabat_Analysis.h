#ifndef Sabat_Analysis_h
#define Sabat_Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Sabat_Analysis {
public :
  TTree          *fChain;
  Int_t           fCurrent;

  Double_t        EnergyDeposit, a, b, c, EnergyDepositFinal=0, EnergyDepositFinal_smr=0, EngSmr=0, TotalEdepGamma, DistanceBtAlphaGamma;
  Double_t        NeutronTheta, NeutronPhi, AlphaTheta, AlphaPhi;
  Double_t        Time;
  Double_t        X;
  Double_t        Y;
  Double_t        Z;
  Int_t           count_rate, Event=0;
  Int_t           ParentID;
  Int_t           StepID;
  Char_t          Particles[6];
  Char_t          Process[22];
  Double_t        TimeL;
  Double_t        KEnergy;
  Int_t           volume;
  Char_t          volume2[13];
  Int_t           TrackID;
  Int_t           EventID;

   // Branches - 
  TBranch        *b_EnergyDeposit;
  TBranch        *b_Time;
  TBranch        *b_X;
  TBranch        *b_Y;
  TBranch        *b_Z;
  TBranch        *b_count_rate;
  TBranch        *b_ParentID;
  TBranch        *b_StepID;
  TBranch        *b_Particles;
  TBranch        *b_Process;
  TBranch        *b_TimeL;
  TBranch        *b_KEnergy;
  TBranch        *b_volume;
  TBranch        *b_volume2;
  TBranch        *b_TrackID;
  TBranch        *b_EventID;
  TBranch        *b_TotalEdepGamma;
  TBranch        *b_NeutronTheta;
  TBranch        *b_NeutronPhi;
  TBranch        *b_AlphaTheta;
  TBranch        *b_AlphaPhi;

  Sabat_Analysis(TTree *tree=0);
  Sabat_Analysis(TString input);
  virtual ~Sabat_Analysis();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

private:
  TString inputName;
};

#endif

#ifdef Sabat_Analysis_cxx
Sabat_Analysis::Sabat_Analysis(TTree *tree) : fChain(0) 
{
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("neutron_merged.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("test.root");
    }
  f->GetObject("ekin_time", tree);
  }
  Init(tree);
}

Sabat_Analysis::Sabat_Analysis(TString input) : fChain(0)
{
  inputName = input;
  TTree *tree = 0;
  TFile *f = new TFile(input);
  if (f != 0) {
    f->GetObject("ekin_time", tree);
    if (tree) {
      std::cout << input << std::endl;
      Init(tree);
      Loop();
    }
  }
}

Sabat_Analysis::~Sabat_Analysis()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t Sabat_Analysis::GetEntry(Long64_t entry)
{
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t Sabat_Analysis::LoadTree(Long64_t entry)
{
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void Sabat_Analysis::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("EnergyDeposit", &EnergyDeposit, &b_EnergyDeposit);
  fChain->SetBranchAddress("Time", &Time, &b_Time);
  fChain->SetBranchAddress("X", &X, &b_X);
  fChain->SetBranchAddress("Y", &Y, &b_Y);
  fChain->SetBranchAddress("Z", &Z, &b_Z);
  fChain->SetBranchAddress("count_rate", &count_rate, &b_count_rate);
  fChain->SetBranchAddress("ParentID", &ParentID, &b_ParentID);
  fChain->SetBranchAddress("StepID", &StepID, &b_StepID);
  fChain->SetBranchAddress("Particles", Particles, &b_Particles);
  fChain->SetBranchAddress("Process", Process, &b_Process);
  fChain->SetBranchAddress("TimeL", &TimeL, &b_TimeL);
  fChain->SetBranchAddress("KEnergy", &KEnergy, &b_KEnergy);
  fChain->SetBranchAddress("volume", &volume, &b_volume);
  fChain->SetBranchAddress("volume2", volume2, &b_volume2);
  fChain->SetBranchAddress("TrackID", &TrackID, &b_TrackID);
  fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
  fChain->SetBranchAddress("TotalEdepGamma", &TotalEdepGamma, &b_TotalEdepGamma);
  fChain->SetBranchAddress("TotalEdepGamma", &TotalEdepGamma, &b_TotalEdepGamma);
  fChain->SetBranchAddress("NeutronTheta", &NeutronTheta, &b_NeutronTheta);
  fChain->SetBranchAddress("NeutronPhi",   &NeutronPhi,   &b_NeutronPhi);
  fChain->SetBranchAddress("AlphaTheta",   &AlphaTheta,   &b_AlphaTheta);
  fChain->SetBranchAddress("AlphaPhi",     &AlphaPhi,     &b_AlphaPhi);
  Notify();
}

Bool_t Sabat_Analysis::Notify()
{
  return kTRUE;
}

void Sabat_Analysis::Show(Long64_t entry)
{
// Print contents of entry.
  if (!fChain)
    return;
  fChain->Show(entry);
}

Int_t Sabat_Analysis::Cut(Long64_t entry)
{
  return 1;
}

#endif // #ifdef Sabat_Analysis_cxx
