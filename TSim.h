#ifndef TSim_h
#define TSim_h

#include "./includes.h"
#include "TComplex.h"
#include "TRandom.h"

//--------------------------------------------------------------------------------------

class TSim : public TObject{ 
  
 public :
  TSim( ); 
  TSim( TString); 
  ~TSim();

  //======================
  //======================
  //== Member Functions ==
  //======================
  //======================
  
  //==================================
  // convert raw root file to root file
  // containing sorted events
  Int_t SortEvents(TString);
  Bool_t SortedROOTFileExists();

  //==================================
  // 
  Int_t GetAsymmetry(TString);
  
  //==================================
    
  Float_t PhotonEnergyToTheta(Double_t);  
  Float_t ElectronEnergyToTheta(Double_t);
  Float_t ThetaToPhotonEnergy(Float_t);
  Float_t ThetaToElectronEnergy(Float_t);

  //=====================
  //====== Setters ======
  
  void    SetStyle();

  //=====================
  //====== Getters ======

  Int_t   GetThetaBin(Float_t);
  void    GetThetaBinValues();
  
  Float_t CrystalToPhi(Int_t); 
  //=====================
  //==== conditions =====
  
  Bool_t GoodTheta(Float_t);

  //====================================
  // overload operator to use with TF1
  double operator() (double *v, double *p) {
    Float_t theta  = DegToRad()*v[0];
    Float_t sintsq = Sin(theta)*Sin(theta);
    Float_t gamma  = 2 - Cos(theta) + 1/(2-Cos(theta));
    return p[0] + 2*sintsq*sintsq/(gamma*gamma - 2*gamma*sintsq);
  }


  void SetAsymmetry(TString);
  
  Float_t GetAsymm(Int_t bin){
    return 0 /*asymmCr90[bin]*/;
  }
  
  Float_t GetAsymErr(Int_t bin){
    return 0 /* asyerCr90[bin]*/;
  }
  
  //======================
  //======================
  //==== Data Members ====
  //======================
  //======================

  //bin sizes
  static const Int_t nCrystals = 9;

  Double_t EA[nCrystals];
  Double_t EB[nCrystals];
  
  Float_t etHA[nCrystals];
  Float_t etHB[nCrystals];
  Float_t ltHA[nCrystals];
  Float_t ltHB[nCrystals];

  Float_t mintHAErr[nCrystals];
  Float_t mintHBErr[nCrystals];
  Float_t maxtHAErr[nCrystals];
  Float_t maxtHBErr[nCrystals];

  Float_t simtHA[nCrystals];
  Float_t simtHB[nCrystals];
  Int_t centrFirst[nCrystals];

  static const Int_t nThbins = 8;
  static const Int_t nPhibins = 4;
  
  Float_t ThMin[nThbins];
  Float_t ThMax[nThbins];
  Float_t plotTheta[nThbins]; 
  Float_t AsymMatrix[nThbins][nPhibins];

  Int_t n000;
  Int_t n090;
  Int_t n180;
  Int_t n270;

  Double_t sigmaA[nCrystals];
  Double_t sigmaB[nCrystals];
  
  // filename without extensions 
  // or folder
  TString rootFileRawName;
  TString rootFileSortName;
  
  // corresponding file
  TFile  *theFile;
  
  // corresponding tree
  TTree *sortDataTree;
  TTree  *simDataTree;
  
  // 
  void Initialise();
  
  // 
  void Loop();
  
 // Declaration of leaf types
  Double_t        crystal0;
  Double_t        crystal1;
  Double_t        crystal2;
  Double_t        crystal3;
  Double_t        crystal4;
  Double_t        crystal5;
  Double_t        crystal6;
  Double_t        crystal7;
  Double_t        crystal8;
  Double_t        crystal9;
  Double_t        crystal10;
  Double_t        crystal11;
  Double_t        crystal12;
  Double_t        crystal13;
  Double_t        crystal14;
  Double_t        crystal15;
  Double_t        crystal16;
  Double_t        crystal17;
  Int_t           TrueEvent;
  Double_t        cosA;
  Double_t        cosB;

  // List of branches
  TBranch        *b_crystal0;   
  TBranch        *b_crystal1;   
  TBranch        *b_crystal2;   
  TBranch        *b_crystal3;   
  TBranch        *b_crystal4;
  TBranch        *b_crystal5;   
  TBranch        *b_crystal6;   
  TBranch        *b_crystal7;   
  TBranch        *b_crystal8;   
  TBranch        *b_crystal9;   
  TBranch        *b_crystal10;   
  TBranch        *b_crystal11;   
  TBranch        *b_crystal12;   
  TBranch        *b_crystal13;   
  TBranch        *b_crystal14;   
  TBranch        *b_crystal15;   
  TBranch        *b_crystal16;   
  TBranch        *b_crystal17;   
  TBranch        *b_TrueEvent;   
  TBranch        *b_cosA;   
  TBranch        *b_cosB;   


   ClassDef(TSim,1);
};



#endif

//--------------------------------------------------------------------------------------
