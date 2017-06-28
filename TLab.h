#ifndef TLab_h
#define TLab_h

#include "./includes.h"
#include "TTheory.h"
#include "TSim.h"

//------------------------------------------------------------------------------------------------

class TLab : public TObject{ 
  
 public :
  TLab( ); 
  TLab( TString); 
  TLab( TString, TString);
  TLab( TString, TString, TString); 
  ~TLab();
  
  //======================
  //======================
  //== Member Functions ==
  //======================
  //======================
  
  void SetFilenames(TString);
  void SetFilenames(TString,TString,TString);

  Bool_t RawROOTFileExists();

  Bool_t CalibratedROOTFileExists();
  
  void MakeRawDataTreeFile();
  
  void MakeCalibratedDataTreeFile();
  
  void SetPedestals();
  Float_t GetPedestal(Int_t);
  
  void SetPhotopeaks();
  Float_t GetPhotopeak(Int_t);
    
  Float_t GetEnergy(Float_t, Int_t);

  Bool_t GoodTiming(Float_t);
  Bool_t GoodTheta(Float_t);

  Float_t ElectronEnergyToTheta(Float_t);
  Float_t PhotonEnergyToTheta(Float_t);

  void CalculateAsymmetry(Int_t,Float_t,Float_t);
  
  void GraphAsymmetry(Char_t);
  
  void SetStyle();
  
  //======================
  //======================
  //==== Data Members ====
  //======================
  //======================
  
  static const Int_t nChannels = 10;
  
  Int_t runNumberInt;

  ifstream *inData;

  TFile  *rootFileRawData;
  TTree  *rawDataTree;
  
  TFile  *rootFileCalData;
  TTree  *calDataTree;
  
  TString textFileName;
  
  TString rootFileRawName;
  TString rootFileCalName;
  
  TCanvas *canvas1;
  TCanvas *canvas2;
  
  // Raw data
  TH1F   *hQ[nChannels];
  TH1F   *hT[nChannels];

  Long64_t eventNumber;

  Float_t Q[nChannels];
  
  Float_t T[nChannels];

  Float_t pedQ[nChannels];
  Float_t phoQ[nChannels];

  // Cal data
  TH1F   *hEA[nChannels];
  TH1F   *hEB[nChannels];
  
  Float_t QA[nChannels];
  Float_t QB[nChannels];

  Float_t EA[nChannels];
  Float_t EB[nChannels];
  
  Float_t TA[nChannels];
  Float_t TB[nChannels];
  
  Float_t tHA[nChannels];
  Float_t tHB[nChannels];
  
  Float_t Asym;
  Float_t AsymErr;
  
  Double_t R000;
  Double_t R090;
  
  Float_t AR;
  Float_t BR;

  Char_t type;

  TSim * simData;

  Int_t npeaks;
  
  ClassDef(TLab,1);
};

#endif

// ------------------------------------------------------------------------------------------------
