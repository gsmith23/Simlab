#ifndef TSim_h
#define TSim_h

#include "./includes.h"
//#include "TTheory.h"


//--------------------------------------------------------------------------------------

class TSim : public TObject{ 
  
 public :
  TSim( ); 
  TSim( TString);
  TSim( TString, TString);
  ~TSim();

  //======================
  //======================
  //== Member Functions ==
  //======================
  //======================
  
  //==================================
  // convert raw root file to root file
  // containing sorted events
  Int_t  SortEvents(TString);
  Bool_t SortedROOTFileExists();

  //==================================
  // 

  
  Int_t CalculateAsymmetrySim(TString);
  Int_t CalculateAsymmetrySimScattered(TString,
				       Float_t,
				       Float_t);
  
  void CalculateABC();
  void CalculateABC_Lab();
  void CalculateABC_True();
  
  void GraphAsymmetrySim(TString, 
			 TString, 
			 Int_t   nThSBins = 0,
			 Float_t thSMin = 0.,
			 Float_t thSMax = 0.);
  
  Int_t CalculateAsymmetryLab(TString);
  void  GraphAsymmetryLab(TString,
			  TString file2);
  
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

  Float_t GetAsymLab(Int_t,Int_t);
  Float_t GetAsymLabErr(Int_t,Int_t);
  Float_t GetAsymLabTrue(Int_t,Int_t);
  Float_t GetAsymLabTrueErr(Int_t,Int_t);
  
  //=====================
  //==== conditions =====
  
  Bool_t GoodTheta(Float_t);
  Bool_t CentralYZ(Double_t);
  Bool_t CentralZ(Double_t);
  Bool_t CentralXA(Double_t);
  Bool_t CentralXB(Double_t);
  //====================================
  // overload operator to use with TF1
  double operator() (double *v, double *p) {
    Float_t theta  = DegToRad()*v[0];
    Float_t sintsq = Sin(theta)*Sin(theta);
    Float_t gamma  = 2 - Cos(theta) + 1/(2-Cos(theta));
    return p[0] + 2*sintsq*sintsq/(gamma*gamma - 2*gamma*sintsq);
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

  Double_t simtHA[nCrystals];
  Double_t simtHB[nCrystals];
  Double_t simPhiA[nCrystals];
  Double_t simPhiB[nCrystals];

  Int_t nb_ComptA[nCrystals];
  Int_t nb_ComptB[nCrystals];

  static const Int_t nThbins = 8;
  static const Int_t nPhibins = 4;
  // can be changed (to a multiple of 4)
  static const Int_t nPhibinsSim = 8; 

  Int_t bin000 = 0;
  Int_t bin090 = nPhibinsSim*1/4;
  Int_t bin180 = nPhibinsSim*2/4;
  Int_t bin270 = nPhibinsSim*3/4;

  Double_t XposA[nCrystals];
  Double_t YposA[nCrystals];
  Double_t ZposA[nCrystals];

  Double_t XposB[nCrystals];
  Double_t YposB[nCrystals];
  Double_t ZposB[nCrystals];
  
  Float_t ThMin[nThbins];
  Float_t ThMax[nThbins];
  Float_t plotTheta[nThbins]; 
  
  Float_t AsymMatrix_sim[nThbins][nPhibinsSim];
  Float_t AsymMatrix[nThbins][nPhibins];
  Float_t AsymTrue[nThbins][nPhibins];

  Float_t AsPhiDiff[nThbins];
  Float_t AePhiDiff[nThbins];
  Float_t AsTrue[nThbins];
  Float_t AeTrue[nThbins];

  Float_t fABC[nThbins][3];
  Float_t pABC[nThbins][3];
  Float_t pA[nThbins];
  Float_t pB[nThbins];
  Float_t pC[nThbins];
  
  Float_t fABC_Lab[nThbins][3];
  Float_t pABC_Lab[nThbins][3];
  Float_t pA_Lab[nThbins];
  Float_t pB_Lab[nThbins];
  Float_t pC_Lab[nThbins];

  Float_t fABC_True[nThbins][3];
  Float_t pABC_True[nThbins][3];
  Float_t pA_True[nThbins];
  Float_t pB_True[nThbins];
  Float_t pC_True[nThbins];

  Int_t n000;
  Int_t n090;
  Int_t n180;
  Int_t n270;

  Double_t P_lab[nThbins];
  Double_t Pq_lab[nThbins];

  Double_t P_true_lab[nThbins];
  Double_t Pq_true_lab[nThbins];
  
  Double_t sigmaA[nCrystals];
  Double_t sigmaB[nCrystals];
  
  // filename without extensions 
  // or folder
  TString rootFileRawName;
  TString rootFileSortName;

  TString rootFileRawName1;
  TString rootFileSortName1;
  TString rootFileRawName2;
  TString rootFileSortName2;
  
    
  TH1F * hDPhi[nPhibinsSim]; 
  
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
  Double_t        edep0;
  Double_t        edep1;
   Double_t        edep2;
   Double_t        edep3;
   Double_t        edep4;
   Double_t        edep5;
   Double_t        edep6;
   Double_t        edep7;
   Double_t        edep8;
   Double_t        edep9;
   Double_t        edep10;
   Double_t        edep11;
   Double_t        edep12;
   Double_t        edep13;
   Double_t        edep14;
   Double_t        edep15;
   Double_t        edep16;
   Double_t        edep17;
   Double_t        edepColl1;
   Double_t        edepColl2;
   Int_t           nb_Compt0;
   Int_t           nb_Compt1;
   Int_t           nb_Compt2;
   Int_t           nb_Compt3;
   Int_t           nb_Compt4;
   Int_t           nb_Compt5;
   Int_t           nb_Compt6;
   Int_t           nb_Compt7;
   Int_t           nb_Compt8;
   Int_t           nb_Compt9;
   Int_t           nb_Compt10;
   Int_t           nb_Compt11;
   Int_t           nb_Compt12;
   Int_t           nb_Compt13;
   Int_t           nb_Compt14;
   Int_t           nb_Compt15;
   Int_t           nb_Compt16;
   Int_t           nb_Compt17;
   Double_t        XposA_1st;
   Double_t        YposA_1st;
   Double_t        ZposA_1st;
   Double_t        XposA_2nd;
   Double_t        YposA_2nd;
   Double_t        ZposA_2nd;
   Double_t        XposB_1st;
   Double_t        YposB_1st;
   Double_t        ZposB_1st;
   Double_t        XposB_2nd;
   Double_t        YposB_2nd;
   Double_t        ZposB_2nd;
   Double_t        ThetaA_1st;
   Double_t        PhiA_1st;
   Double_t        ThetaB_1st;
   Double_t        PhiB_1st;
   Double_t        dPhi_1st;
   Double_t        ThetaA_2nd;
   Double_t        PhiA_2nd;
   Double_t        ThetaB_2nd;
   Double_t        PhiB_2nd;
   Double_t        dPhi_2nd;

   // List of branches
   TBranch        *b_edep0;   //!
   TBranch        *b_edep1;   //!
   TBranch        *b_edep2;   //!
   TBranch        *b_edep3;   //!
   TBranch        *b_edep4;   //!
   TBranch        *b_edep5;   //!
   TBranch        *b_edep6;   //!
   TBranch        *b_edep7;   //!
   TBranch        *b_edep8;   //!
   TBranch        *b_edep9;   //!
   TBranch        *b_edep10;   //!
   TBranch        *b_edep11;   //!
   TBranch        *b_edep12;   //!
   TBranch        *b_edep13;   //!
   TBranch        *b_edep14;   //!
   TBranch        *b_edep15;   //!
   TBranch        *b_edep16;   //!
   TBranch        *b_edep17;   //!
   TBranch        *b_edepColl1;   //!
   TBranch        *b_edepColl2;   //!
   TBranch        *b_nb_Compt0;   //!
   TBranch        *b_nb_Compt1;   //!
   TBranch        *b_nb_Compt2;   //!
   TBranch        *b_nb_Compt3;   //!
   TBranch        *b_nb_Compt4;   //!
   TBranch        *b_nb_Compt5;   //!
   TBranch        *b_nb_Compt6;   //!
   TBranch        *b_nb_Compt7;   //!
   TBranch        *b_nb_Compt8;   //!
   TBranch        *b_nb_Compt9;   //!
   TBranch        *b_nb_Compt10;   //!
   TBranch        *b_nb_Compt11;   //!
   TBranch        *b_nb_Compt12;   //!
   TBranch        *b_nb_Compt13;   //!
   TBranch        *b_nb_Compt14;   //!
   TBranch        *b_nb_Compt15;   //!
   TBranch        *b_nb_Compt16;   //!
   TBranch        *b_nb_Compt17;   //!
   TBranch        *b_XposA_1st;   //!
   TBranch        *b_YposA_1st;   //!
   TBranch        *b_ZposA_1st;   //!
   TBranch        *b_XposA_2nd;   //!
   TBranch        *b_YposA_2nd;   //!
   TBranch        *b_ZposA_2nd;   //!
   TBranch        *b_XposB_1st;   //!
   TBranch        *b_YposB_1st;   //!
   TBranch        *b_ZposB_1st;   //!
   TBranch        *b_XposB_2nd;   //!
   TBranch        *b_YposB_2nd;   //!
   TBranch        *b_ZposB_2nd;   //!
   TBranch        *b_ThetaA_1st;   //!
   TBranch        *b_PhiA_1st;   //!
   TBranch        *b_ThetaB_1st;   //!
   TBranch        *b_PhiB_1st;   //!
   TBranch        *b_dPhi_1st;   //!
   TBranch        *b_ThetaA_2nd;   //!
   TBranch        *b_PhiA_2nd;   //!
   TBranch        *b_ThetaB_2nd;   //!
   TBranch        *b_PhiB_2nd;   //!
   TBranch        *b_dPhi_2nd;   //!


   ClassDef(TSim,1);
};



#endif

//--------------------------------------------------------------------------------------
