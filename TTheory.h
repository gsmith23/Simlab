#ifndef TTheory_h
#define TTheory_h

#include "./includes.h"
//--------------------------------------------------------------------------------------

class TTheory : public TObject{ 
  
 public :
  TTheory( ); 
  ~TTheory();

  //======================
  //======================
  //== Member Functions ==
  //======================
  //======================

  //===================================
  // Theoretical curves for
  // finite geometry 
  void GraphFiniteAsymmetry(Int_t,
			    Float_t,
			    Char_t);
  
  Float_t X(Float_t);
  Float_t J(Float_t,Float_t);
  Float_t JLim(Float_t);
  Float_t Jdash(Float_t,Float_t);
  Float_t JdashLim(Float_t);
  Float_t rho1(Float_t, Float_t);
  
  Float_t u(Float_t);
  Float_t w(Float_t);
  Float_t Z(Float_t);
  Float_t rho2(Float_t,Float_t,Float_t);

  Float_t ThetaToPhotonEnergy(Float_t);
  Float_t ThetaToElectronEnergy(Float_t);

  //======================
  //======================
  //==== Data Members ====
  //======================
  //======================
  
  ClassDef(TTheory,1);
};

#endif
