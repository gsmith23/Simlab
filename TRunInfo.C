#include "TRunInfo.h"
#include "iostream"

using namespace std;

#if !defined(__CINT__)
ClassImp(TRunInfo)
#endif

TRunInfo::TRunInfo(){

}

TRunInfo::TRunInfo(int run){
  
  runNumber = run;
  SetEventNumbers(run); // and onePart

}

TRunInfo::~TRunInfo(){}

void TRunInfo::SetEventNumbers(int run){

  nOR1 = 0;
  nAND = 0;
  nOR2 = 0;
  eventSum = 0.;
  
  if     (run == 0){
    cout << endl;
    cout << " Is your mind entangled ?? " << endl; 
  }
  else if(run == 423){
    // Runs: 423 - testing code
    nOR1 = 0; 
    nAND = 210910;
    nOR2 = 0;  
  }
  // ------------------------------------------
  // NB from here all thresholds were at 300 mV
  // ------------------------------------------
  else if(run == 400){ // collimator run
    // Runs: 451 (OR), 452 (AND), 450 (OR)
    nOR1 = 1019283; // 451 - OR collimators
    nAND = 6675453; // 452 - AND collimators
    nOR2 = 890963;  // 450 - closest previous no collimators run
  }
  else if(run == 449){
    // Runs: 449 (OR + AND), 450 (OR)
    nOR1 = 500000;  // start of run 449 was OR
    nAND = 6580429; // rest of 449 was AND
    nOR2 = 890963;  // 450 OR 
  }
  else if(run == 454){
    // Runs: 451 (OR), 452 (AND), 453 (OR) 
    nOR1 = 1019283; // 451 - OR collimators (unclear photopeaks)
    nAND = 6675453; // 452 - AND collimators
    nOR2 = 1360796; // 453 - OR collimators (unclear photopeaks)
  }
  else if(run == 1454){ // longest collimator run
    // Runs: 453 (OR), 454 (AND), 456 (OR)
    nOR1 = 1360796;  // 453 - OR, collimators
    nAND = 38528184; // 454 - AND, collimators, 0.1 mV HV drift
    nOR2 = 1459920;  // 456 - OR, closest no collimator run after
  }
  // ------------------------------------------
  // NB Now crystals - source = 3.0 cm 
  // Collimators were also removed
  // check logbook for previous values
  // ------------------------------------------
  else if(run == 1457 ){ 
    // Runs: 456 (OR), 457 (AND), 458 (OR)
    nOR1 = 1524550;  // 456 - OR, 
    nAND = 65734400; // 457 - AND
    nOR2 = 1572198;  // 458 - OR
  }
  // ------------------------------------------
  // ------------------------------------------
  // NB Best runs from here on. 
  // Central channels thesholds 
  // set to 100 mV 
  // and crystals at 3.0 cm
  // ------------------------------------------
  // ------------------------------------------
  else if(run == 1460){ 
    // Runs: 459 (or), 460 (AND), 461 (OR)
    nOR1 = 1395877;
    nAND = 47142893; 
    nOR2 = 1228535;
  }
  else if(run == 14600){ 
    // Runs: 459 (or), 460 (AND), 461 (OR)
    nOR1 = 1395877;
    nAND = 3111673; 
    nOR2 = 0;
  }
  else if(run == 1470){
    // Runs: 467 (or), 469 (AND), 470 (OR)
    nOR1 = 1121370;   // 467 - OR prior to power up/down 
    nAND = 170907253; // 469 - AND (468 was interupted)
    nOR2 = 1269750;   // 470 - OR
  }
  else if(run == 14700){
    // Runs: 467 (or), 469 (AND), 470 (OR)
    nOR1 = 100000; // 467 - OR prior to power up/down 
    nAND = 1000000; // 469 - AND (468 was interupted)
    nOR2 = 100000; // 470 - OR
  }
  else if(run == 123){
    nOR1 = 11111111;
    nAND = 22222222;
    nOR2 = 33333333;
  }
  //-------
  // NB New set up with active corner crystals
  //-------
  //-------
  // Test with two files
  //-------
  else if(run == 49801){
    nOR1 = 0;
    nAND = 29824032; //AND on
    nOR2 = 0;
    onePart = kTRUE;
  }
  else if(run == 4985){
    nOR1 = 0;
    nAND = 4999; //AND on
    nOR2 = 0;
    onePart = kTRUE;
  }
  else if(run == 4980){
    nOR1 = 0;
    nAND = 500000; 
    nOR2 = 0;
    onePart = kTRUE;
  }
  else if(run == 5010){
    nOR1 = 0;
    nAND = 1000000; 
    nOR2 = 0;
    onePart = kTRUE;
  }
  else if(run == 501500){
    nOR1 = 0;
    nAND = 500; 
    nOR2 = 0;
    onePart = kTRUE;
  }
  else{

    cout << endl;
    cout << "----------------------------------------------------" << endl; 
    cout << " Warning: no event number information for this run. " << endl; 
    cout << "----------------------------------------------------" << endl; 
    cout << endl;
  }

  eventSum = nOR1 + nAND + nOR2;

  cout << endl;
  cout << " \t\t\t Events " << endl; 
  cout << " First  OR  run \t" 
       << nOR1 << endl;
  cout << " Main   AND run \t" 
       << nAND << endl;
  cout << " Second OR  run \t" 
       << nOR2  << endl;
  cout << " Total          \t" 
       << eventSum  << endl;
  
}



Bool_t TRunInfo::QIsInComptonRange(Float_t Q,
					Int_t ch){
  
  Float_t lowQ  = Q;
  Float_t highQ = 0;
  
  // mimic 5010 to data with no outer ch photopeaks
  if(runNumber == 5010)
    highQ = Q;

  if(onePart){
    // run 501 (should be okay for 498)
    switch(ch){
    case 0  : if( lowQ > 1100 && highQ < 2200 ) return kTRUE; break;
    case 1  : if( lowQ > 1050 && highQ < 2200 ) return kTRUE; break;
    case 2  : if( lowQ >  900 && highQ < 2100 ) return kTRUE; break;
    case 3  : if( lowQ > 1000 && highQ < 2400 ) return kTRUE; break;
    case 4  : if( lowQ > 1100 && highQ < 2200 ) return kTRUE; break;
    case 5  : if( lowQ > 1100 && highQ < 2300 ) return kTRUE; break;
    case 6  : if( lowQ > 1000 && highQ < 2400 ) return kTRUE; break;
    case 7  : if( lowQ > 900  && highQ < 2200 ) return kTRUE; break;
    case 8  : if( lowQ > 975  && highQ < 2200 ) return kTRUE; break;
    case 9  : if( lowQ > 1000 && highQ < 2300 ) return kTRUE; break;
    }
  }
  else if(runNumber == 1470  ||
	  runNumber == 14700 ){
    switch(ch){
    case 0  : if( lowQ > 1100 && highQ < 2200 ) return kTRUE; break;
    case 1  : if( lowQ > 1050 && highQ < 2200 ) return kTRUE; break;
    case 2  : if( lowQ >  900 && highQ < 2100 ) return kTRUE; break;
    case 3  : if( lowQ > 1000 && highQ < 2400 ) return kTRUE; break;
    case 4  : if( lowQ > 1100 && highQ < 2200 ) return kTRUE; break;
    case 5  : if( lowQ > 1100 && highQ < 2300 ) return kTRUE; break;
    case 6  : if( lowQ > 1000 && highQ < 2400 ) return kTRUE; break;
    case 7  : if( lowQ > 900  && highQ < 2200 ) return kTRUE; break;
    case 8  : if( lowQ > 975  && highQ < 2200 ) return kTRUE; break;
    case 9  : if( lowQ > 1000 && highQ < 2300 ) return kTRUE; break;
    }
  }
  return kFALSE;
    
}

Int_t TRunInfo::GetMinQ(){

  Int_t minQ = 2200;
  if    ( runNumber == 1460 ||
	  runNumber == 14600){
    minQ = 2700;
  }
  else if( runNumber == 1470  ||
	   runNumber == 14700){
    minQ = 2500;
  }
  return minQ;
}

Int_t TRunInfo::GetMaxQ(){

  Int_t maxQ = 3600;

  if( runNumber == 1460  ||
      runNumber == 14600 ||
      runNumber == 1470  ||
      runNumber == 14700 ){
    maxQ = 3900;
  }
  
  return maxQ;
}

Float_t TRunInfo::GetPhotoStartVal(int ch,
					int part,
					Bool_t onePart){

    for (Int_t ch = 0 ; ch < nChannels ; ch++)
      for (Int_t part = 0 ; part < nParts ; part++)
	phoQStVal[ch][part] =  2600.;
  
  if(runNumber==1460||
     runNumber==14600 ){
    phoQStVal[0][0] = 3340., phoQStVal[1][0] = 3420.;
    phoQStVal[2][1] = 3300., phoQStVal[3][0] = 3140.;
    phoQStVal[4][0] = 3340., phoQStVal[5][0] = 3500.;
    phoQStVal[6][0] = 3440., phoQStVal[7][1] = 3660.; 
    phoQStVal[8][0] = 3410., phoQStVal[9][0] = 3050.; 
    
  }
  if(runNumber==1470 || 
     runNumber==14700 ){
  
    // channel 9 drifted from previous run
    phoQStVal[0][2] = 3350., phoQStVal[1][2] = 3420.;
    phoQStVal[2][1] = 3300., phoQStVal[3][2] = 3140.;
    phoQStVal[4][2] = 3340., phoQStVal[5][2] = 3500.;
    phoQStVal[6][2] = 3440., phoQStVal[7][1] = 3660.; 
    phoQStVal[8][2] = 3410., phoQStVal[9][2] = 2800.; 
    
  }
  else if( onePart ){
    cout << endl;
    cout << " Using run 501 values " << endl;
    
    phoQStVal[0][1] = 2596., phoQStVal[1][1] = 2625.;
    phoQStVal[2][1] = 2763., phoQStVal[3][1] = 2902.;
    phoQStVal[4][1] = 2728., phoQStVal[5][1] = 2800.;
    phoQStVal[6][1] = 2656., phoQStVal[7][1] = 2629.; 
    phoQStVal[8][1] = 2588., phoQStVal[9][1] = 2742.; 
    
  }
  
  return phoQStVal[ch][part];  
}

Float_t TRunInfo::GetPedStartVal(){  
  return 600.0;
}

Long64_t TRunInfo::GetnOR1(){  
  return nOR1;
}  

Long64_t TRunInfo::GetnAND(){  
  return nAND;
}  

Long64_t TRunInfo::GetnOR2(){  
  return nOR2;
}  

Long64_t TRunInfo::GetEventSum(){  
  return eventSum;
}  

Bool_t TRunInfo::GetOnePart(){  
  return onePart;
}  
