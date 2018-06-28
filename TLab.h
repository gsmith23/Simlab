#ifndef TLab_h
#define TLab_h

#include "./includes.h"
#include "TRunInfo.h"

//------------------------------------------------------------------------------------------------

class TLab : public TObject{ 
  
 public :
  TLab( ); 
  TLab( TString); 
  TLab( TString, TString);
  TLab( TString, TString, TString); 
  TLab( TString, TString, TString, TString); 
  ~TLab();
  
  //======================
  //======================
  //== Member Functions ==
  //======================
  //======================
  
  void SetInfo(TString);

  Bool_t RawTextFileExists();
  Bool_t RawROOTFileExists();

  Bool_t CalibratedROOTFileExists();
  
  void MakeRawDataTreeFile();

  void MakeCalibratedDataTreeFile();
  
  void    SetPedestals();
  Float_t GetPedestal(Int_t);
  Int_t   DefaultPedestalPart();
  
  void    FillQSumHistos();
  
  void    SetPhotopeaks();
  void    InitPhotopeaks();
  
  Int_t   GetCentralIndex(Int_t);

  Bool_t  photopeaksInitByChan = kFALSE;

  Int_t   GetMinQ();
  Int_t   GetMaxQ();
  
  Bool_t  DoFitPhotopeaks();
  void    FitPhotopeaks();
  
  Float_t GetPhotopeak(Int_t);
  Int_t   DefaultPhotopeakPart(Int_t);
  
  void    HWHM_Q_To_E();
  
  Bool_t QIsInComptonRange(Float_t, Int_t);

  Bool_t GoodTiming(Float_t);
  Bool_t GoodTheta(Float_t);

  Float_t ElectronEnergyToTheta(Float_t);
  Float_t PhotonEnergyToTheta(Float_t);

  Float_t ThetaToThetaError(Float_t, Int_t);
  Float_t ThetaToPhotonEnergy(Float_t);
  Float_t ThetaToElectronEnergy(Float_t);

  Int_t Chan2ArrayA(Int_t channel);
  Int_t Chan2ArrayB(Int_t channel);
  
  void CalculateAsymmetry();

  void GetThetaBinValues();
  
  void GraphAsymmetry(Char_t);
  
  void SetStyle();

  Float_t RandomLabPhi();
  Bool_t  RandomGoodLabPhi(Float_t, Int_t);
  
  void  SetEventNumbers(Int_t);

  //======================
  //======================
  //==== Data Members ====
  //======================
  //======================

  
  // OR, AND, OR
  const static Int_t nParts = 3;

  Bool_t onePart = kFALSE;
    
  // crystals per array
  static const Int_t nCrystals = 9;
  
  // only five per array are recorded
  //static const Int_t nChannels = 18;
  static const Int_t nChannels = 10;
  
  // For Graphing
  static const Int_t nPhiBins = 4;

  static const Int_t nThBins  = 7;
  //static const Int_t nThBins  = 1;

/*   Float_t thetaLowEdge  = 10.0; */
/*   Float_t thetaHighEdge = 170.; */
  
  // 39, 54, 69, 84, 99, 114, 129
  Float_t thetaLowEdge  = 31.5;
  Float_t thetaHighEdge = 136.5;

  Long64_t nOR1;
  Long64_t nAND;
  Long64_t nOR2;
  Long64_t eventSum;
  
  Float_t ThMin[nThBins];
  Float_t ThMax[nThBins];
  Float_t plotTheta[nThBins]; 
  Float_t AsymMatrix[nThBins][nPhiBins];
  
  Int_t   runNumberInt;
  
  TString simRun;
  TString simRunU;
  TString simRunP;
  TString simRunE;
  
  ifstream *inData;
  
  TRunInfo * runInfo = nullptr;

  TFile  *rootFileRawData;
  TTree  *rawDataTree;
  
  TFile  *rootFileCalData;
  TTree  *calDataTree;
  
  TString textFileName;
  TString rootFileRawName;
  TString rootFileCalName;

  TCanvas *canvas1;
  TCanvas *canvas2;
  
  // pre-run/part OR data
  TH1F   *hQ_0[nChannels];
  // main run/part
  TH1F   *hQ_1[nChannels];
  // post-run/part OR data
  TH1F   *hQ_2[nChannels];

  TH1F   *hT[nChannels];
  
  // main run/part outer summed with inner
  // (Compton region only)
  TH1F   *hQQ_1[nChannels];

  Long64_t eventNumber;

  Float_t Q[nChannels];
  
  Float_t T[nChannels];

  // Fit Results
  Float_t pedQ[nChannels][nParts];
  Float_t phoQ[nChannels][nParts];
  Float_t HWHM[nChannels][nParts];

  // Cal data
  TH1F   *hEA[nCrystals];
  TH1F   *hEB[nCrystals];
  
  Float_t QA[nCrystals];
  Float_t QB[nCrystals];

  Float_t EA[nCrystals];
  Float_t EB[nCrystals];
  
  Float_t TA[nCrystals];
  Float_t TB[nCrystals];
  
  Float_t tHA[nCrystals];
  Float_t tHB[nCrystals];

  Float_t tHAErr[nCrystals];
  Float_t tHBErr[nCrystals];

  Float_t Asym;
  Float_t AsymErr;
  
  Float_t AsymPhi[nPhiBins];
  Float_t AsymPhiErr[nPhiBins];
  
  Double_t R000;
  Double_t R090;
  
  Float_t AR;
  Float_t BR;

  Char_t type;

  TSim *simData;

  TH1F *hr;
  
  ClassDef(TLab,1);
};

#endif

// ------------------------------------------------------------------------------------------------
