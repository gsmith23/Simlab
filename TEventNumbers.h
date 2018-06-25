#ifndef TEventNumbers_h
#define TEventNumbers_h

#include "TObject.h"

class TEventNumbers : public TObject{

 public :
  TEventNumbers();
  TEventNumbers(int);
  ~TEventNumbers();
  
  void SetEventNumbers(int run);
  Bool_t QIsInComptonRange(Float_t Q,
			   Int_t ch);
  
  Bool_t   GetOnePart();
  
  Long64_t GetnOR1();
  Long64_t GetnAND();
  Long64_t GetnOR2();
  Long64_t GetEventSum();
  
  Int_t GetMinQ();
  Int_t GetMaxQ();
  
  Float_t GetPhotoStartVal(Int_t ch, Int_t part, Bool_t onePart);
  Float_t GetPedStartVal();

  Bool_t onePart = kFALSE;
  
  Long64_t nOR1;
  Long64_t nAND;
  Long64_t nOR2;
  Long64_t eventSum;

  Int_t   runNumber;

  static const Int_t nParts    = 3;
  static const Int_t nChannels = 18;

  Float_t phoQStVal[nChannels][nParts];
  Float_t pedQStVal[nChannels][nParts];
  
  ClassDef(TEventNumbers,1);
};

#endif
