#ifndef TSim_h
#define TSim_h

#include "./includes.h"

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
  // convert text file to root file
  // containing tree
  Int_t Hits2Tree(TString,Int_t);
  
  //==================================
  // convert raw root file to root file
  // containing sorted events
  Int_t SortEvents(TString,Int_t);

  //==================================
  // 
  Int_t GetAsymmetry(TString);
  
  //==================================
    
  Float_t PhotonEnergyToTheta(Float_t);  
  Float_t ElectronEnergyToTheta(Float_t);
  Float_t ThetaToPhotonEnergy(Float_t);
  Float_t ThetaToElectronEnergy(Float_t);

  //=====================
  //====== Setters ======
  
  void    SetStyle();

  //=====================
  //====== Getters ======

  Int_t   GetDPhiBin(Float_t, Float_t);
  Int_t   GetThetaBin(Float_t,Int_t);
  
  Float_t GetAverageEnergy();
  
  Float_t GetHitAngle(Float_t,Float_t);  
  //=====================
  //==== conditions =====
  
  Bool_t GoodNumHits(Int_t);
  Bool_t GoodLOR(Float_t,Char_t);
  Bool_t GoodHitAngle(Float_t,Float_t);
  Bool_t GoodHitSeparation(Float_t);
  Bool_t GoodDTheta(Float_t,Float_t,Float_t);
  Bool_t GoodE(UInt_t,Float_t,Float_t);
  Bool_t GoodScatterDistance(Float_t);
  Bool_t GoodScatterDistances(Float_t,Float_t);
  Bool_t BadDPhi(Float_t);
  Bool_t GoodTheta(Float_t);

  //====================================
  // overload operator to use with TF1
  double operator() (double *v, double *p) {
    Float_t theta  = DegToRad()*v[0];
    Float_t sintsq = Sin(theta)*Sin(theta);
    Float_t gamma  = 2 - Cos(theta) + 1/(2-Cos(theta));
    return p[0] + 2*sintsq*sintsq/(gamma*gamma - 2*gamma*sintsq);
  }
  
  void PlotTF1();

  void SetAsymmetry(TString);
  
  Float_t GetAsymm(Int_t bin){
    return asymmCr90[bin];
  }
  
  Float_t GetAsymErr(Int_t bin){
    return asyerCr90[bin];
  }
  
  //======================
  //======================
  //==== Data Members ====
  //======================
  //======================
  
  // filename without extensions 
  // or folder
  TString rootFile;
  
  // corresponding file
  TFile  *theFile;
  
  // corresponding tree
  TTree  *theTree;
  
  // 
  void Initialise();
  
  // 
  void Loop();
  
  // Branch Leaves
  Int_t           nHits;
  Int_t           nParticles;
  Long64_t        eventNumber[32];
  Char_t          detector[32];
  Int_t           detectorUnitID[32];
  Float_t         energy[32]; 
  Float_t         x[32];
  Float_t         y[32];
  Float_t         z[32];
  Float_t         minTime[32];
  Float_t         maxTime[32];
  Int_t           nOriginalTracks[32];  
  Int_t           originalParticleNumber[32]; 
  Int_t           nTracks[32];
  Int_t           particleNumber[32][16];
  Int_t           nHitsC;
  Int_t           nHitsW;
  Int_t           eventType;
  Float_t         iC[32];
  Float_t         iW[32];
  Int_t           i1P1;
  Int_t           i1P2;
  Float_t         sinogramR;
  Float_t         sinogramPhi;
  Float_t         meanHitAngle;
  Float_t         energy1P1;
  Float_t         energy1P2;
  Int_t           i2P1;
  Int_t           i2P2;
  Float_t         energy2P1;
  Float_t         energy2P2;
  Float_t         thetaP1E;
  Float_t         thetaP2E;
  Float_t         thetaP1V;
  Float_t         thetaP2V;
  Float_t         thetaP1VTrue;
  Float_t         thetaP2VTrue;
  Float_t         phiP1;
  Float_t         phiP2;
  Float_t         dPhi;

  TVector3        p3P1Lab;
  TVector3        p3P2Lab;
  TVector3        LOR;
  TVector3        LORXY;
  TVector3        v3P1Lab;
  TVector3        v3P2Lab;
  TVector2        midpointXY;
  TVector3        midpoint;
  TVector3        v3ScP1Lab;
  TVector3        v3ScP2Lab;
  TVector3        v3P1Event;
  TVector3        v3P2Event;
  
  // List of branches
  TBranch        *b_nHits;
  TBranch        *b_nParticles;
  TBranch        *b_eventNumber;
  TBranch        *b_detector;
  TBranch        *b_detectorUnitID;
  TBranch        *b_energy;
  TBranch        *b_x;
  TBranch        *b_y;
  TBranch        *b_z;
  TBranch        *b_nOriginalTracks;
  TBranch        *b_originalParticleNumber;
  TBranch        *b_nTracks;
  TBranch        *b_nHitsC;
  TBranch        *b_nHitsW;
  TBranch        *b_eventType; 
  TBranch        *b_iC; 
  TBranch        *b_iW;
  TBranch        *b_i1P1;
  TBranch        *b_i1P2;
  TBranch        *b_sinogramR;
  TBranch        *b_sinogramPhi;
  TBranch        *b_meanHitAngle;
  TBranch        *b_energy1P1;
  TBranch        *b_energy1P2;
  TBranch        *b_i2P1;   
  TBranch        *b_i2P2;   
  TBranch        *b_energy2P1;
  TBranch        *b_energy2P2;
  TBranch        *b_thetaP1E; 
  TBranch        *b_thetaP2E; 
  TBranch        *b_thetaP1V; 
  TBranch        *b_thetaP2V;
  TBranch        *b_thetaP1VTrue;
  TBranch        *b_thetaP2VTrue;
  TBranch        *b_phiP1;   
  TBranch        *b_phiP2;  
  TBranch        *b_dPhi;   

  //  
  Float_t asymm90[32];
  Float_t asyer90[32];

  Float_t asymmCr90[32];
  Float_t asyerCr90[32];

   ClassDef(TSim,1);
};



#endif

//--------------------------------------------------------------------------------------
