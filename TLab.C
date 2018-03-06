#include "TLab.h"

#if !defined(__CINT__)
ClassImp(TLab)
#endif

// Default Constructor
TLab::TLab( ) {
  cout <<  endl;
  cout << " Default constructor. " << endl;
}

// option for use with one raw file 
TLab::TLab(TString runNumber) {
  SetFilenames(runNumber);
}

//option for use with one raw file and one sim file
TLab::TLab(TString runNumber,
	   TString simNumber) {
  simRun = simNumber;
  SetFilenames(runNumber);
}

// option for use with two simulated files
// one being unpolarised
TLab::TLab(TString runNumber,
	   TString simNumber,
	   TString simNumberU) {
  
  simRun  = simNumber;
  simRunU = simNumberU;
  SetFilenames(runNumber);
  
}

//-------------------------------------------------

TLab::~TLab()
{
  cout << endl;
  cout << " TLab object has been destructed!" << endl;
}

// ------------------------------------------------

void TLab::SetFilenames(TString runNumber){

  cout << endl;
  cout << " TLab object has been created " << endl;
  
  runNumberInt = runNumber.Atoi();
  
  cout << endl;
  cout << " Run Number = " << runNumberInt << endl;

  textFileName    = runNumber + ".txt";
  rootFileRawName = runNumber + ".root";  
  rootFileCalName = runNumber + ".root";
  
  textFileName    = "../Data/run" + textFileName;
  rootFileRawName = "../Data/run" + rootFileRawName;
  rootFileCalName = "../Data/cal" + rootFileCalName;

}


/** Public member functions *********/

Bool_t TLab::RawTextFileExists(){

  inData = new ifstream(textFileName);
  
  return inData;
}

Bool_t TLab::RawROOTFileExists(){
  
  TFile *file = TFile::Open(rootFileRawName);
    
  return file;
}


void TLab::SetEventNumbers(Int_t run){

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
  else if(run == 1470){
    // Runs: 467 (or), 469 (AND), 470 (OR)
    nOR1 = 1121370;   // 467 - OR prior to power up/down 
    nAND = 170907253; // 469 - AND (468 was interupted)
    nOR2 = 1269750;   // 470 - OR
  }
  else if(run == 123){
    nOR1 = 11111111;
    nAND = 22222222;
    nOR2 = 33333333;
  }
  

  eventSum = nOR1 + nAND + nOR2;

  cout << endl;
  cout << " The first calibration run had " 
       << nOR1 << " events " << endl;
  cout << " The main run had " 
       << nAND  << " events " << endl;
  cout << " The second calibration run had " 
       << nOR2 << " events " << endl;
  cout << " In total there are  " 
       << eventSum << " events " << endl;
  
}

void TLab::MakeRawDataTreeFile(){

  SetEventNumbers(runNumberInt);

  TString textFile = "??";
  TString rootFile = "??";
  
  eventNumber = 0;
  
  for(Int_t i = 0 ; i < nChannels ; i++ ){
    Q[i] = 0.0;
    T[i] = 0.0;
  }
  
  textFile = textFileName;
  rootFile = rootFileRawName;
    
  inData = new ifstream(textFile);
  
  rootFileRawData = new TFile(rootFile,"RECREATE","Raw LYSO data");
  rawDataTree     = new TTree("rawDataTree",
			      "LYSO data QDC and TDC values");

  rawDataTree->Branch("eventNumber",
		      &eventNumber,
		      "eventNumber/L");
  
  TString tempString = ""; 
  
  tempString.Form("Q[%d]/F",nChannels);
  rawDataTree->Branch("Q",Q,tempString);
  
  tempString.Form("T[%d]/F",nChannels);
  rawDataTree->Branch("T",T,tempString);
  
  TString nameHist;
  TString titleHist;
  
  for(Int_t i = 0 ; i < nChannels ; i++ ){

    // usually three runs: OR, AND, OR
    for(Int_t run = 0 ; run < nRuns ; run++ ){    

      nameHist.Form("hQ%d_%d",i,run);
      titleHist.Form("hQ%d_%d;QDC bin;Counts",i,run);
      hQ[i][run] = new TH1F(nameHist,titleHist,4096,0,4096);
      
    }

    nameHist.Form("hT%d",i);
    titleHist.Form("hT%d;TDC bin;Counts",i);
    hT[i] = new TH1F(nameHist,titleHist,5200,0,5200);
  }
  cout << endl;
  cout << " Making tree " << endl;
  
  eventNumber = 0;
  
  Int_t index = 0;

  while(*inData 
	>> Q[0+index] 
	>> Q[1+index] 
	>> Q[2+index] 
	>> Q[3+index]  
	>> Q[4+index]  
	>> T[0+index] 
	>> T[1+index] 
	>> T[2+index] 
	>> T[3+index]  
	>> T[4+index] 
	){

    if(eventNumber==0){
      
      if(index==0){
	cout << endl;
	cout << " Line format: " << endl;
      }
      cout << " "
	   << Q[0+index] << " " 
	   << Q[1+index] << " " 
	   << Q[2+index] << " " 
	   << Q[3+index] << " " 
	   << Q[4+index] << " " 
	   << T[0+index] << " " 
	   << T[1+index] << " " 
	   << T[2+index] << " " 
	   << T[3+index] << " " 
	   << T[4+index] << " " 
	   << endl;
    }

    // raw text file has event data over two lines  
    // with five channels per line
    for(Int_t i = index ; i < (index+5) ; i++ ){
    
      // one histogram per channel and per run
      if      ( eventNumber < nOR1 )
	hQ[i][0]->Fill(Q[i]);
      else if ( eventNumber > (nOR1+nAND)){
	hQ[i][2]->Fill(Q[i]);
      }
      else
      hQ[i][1]->Fill(Q[i]);
      
      hT[i]->Fill(T[i]);
      
    } // end of: for(Int_t i = 0 ; i < 16 ; i++ ...
    
    if(index==5){
      rawDataTree->Fill();
      eventNumber++;
    }
    
    if(index==0)
      index = 5;
    else
      index = 0;
    
  }// end of: while( .....
  
  cout << endl;
  cout << " Making " << rootFile << endl;
  cout << " " << eventNumber << " events " << endl;
  

  rawDataTree->Write();
  
  // writes the histograms
  rootFileRawData->Write();
  
  rootFileRawData->Close();
  
}

Bool_t TLab::CalibratedROOTFileExists(){
  
  TFile *file = TFile::Open(rootFileCalName);
  return file;
}

Int_t TLab::Chan2ArrayA(Int_t channel){
  
  Int_t crystal = -1;


  if     (channel == 0)
    crystal = 1;
  else if(channel == 1)
    crystal = 3;
  else if(channel == 2)
    crystal = 4;
  else if(channel == 3)
    crystal = 5;
  else if(channel == 4)
    crystal = 7;
  
  return crystal;
}

Int_t TLab::Chan2ArrayB(Int_t channel){
  
  Int_t crystal = -1;
  
  if     (channel == 5)
    crystal = 1;
  else if(channel == 6)
    crystal = 5;
  else if(channel == 7)
    crystal = 4;
  else if(channel == 8)
    crystal = 3;
  else if(channel == 9)
    crystal = 7;
  
  return crystal;
}

void TLab::MakeCalibratedDataTreeFile(){
  
  cout << endl;
  cout << " Making calibrated data tree " << endl;

  SetPedestals();
  FitPhotopeaks();

  Float_t temp_phoQ = 0.;
  
  // HWHM QDC to energy
  for (Int_t run = 0 ; run < nRuns ; run++){
      cout << endl;
      
      for (Int_t i=0; i<nChannels; i++){
      
      // Set to use OR data after AND data
      temp_phoQ = GetPhotopeak(i) - GetPedestal(i);
      HWHM[i][run] = HWHM[i][run]*511./temp_phoQ;

      cout << "HWHM["<< i 
	   << "]["   << run 
	   << "] = " << HWHM[i][run] 
	   << endl;
    }
  }
  
  cout << endl;
  cout << "Post photopeak setting rootFileRawName is " 
       << rootFileRawName << endl;
  
  // read from this
  rootFileRawData = new TFile(rootFileRawName);
  rawDataTree     = (TTree*)rootFileRawData->Get("rawDataTree");
  
  // write to this
  rootFileCalData = new TFile(rootFileCalName,
			      "RECREATE","Calibrated LYSO data");
  
  calDataTree     = new TTree("calDataTree",
			      "LYSO data QDC and TDC values");
  
  TString tempString = ""; 
  
  rawDataTree->SetBranchAddress("Q",Q);
  rawDataTree->SetBranchAddress("T",T);
  rawDataTree->SetBranchAddress("eventNumber",&eventNumber);
  
  calDataTree->Branch("eventNumber",
		      &eventNumber,
		      "eventNumber/L");
  
  tempString.Form("EA[%d]/F",nCrystals);
  calDataTree->Branch("EA",EA,tempString);

  tempString.Form("EB[%d]/F",nCrystals);
  calDataTree->Branch("EB",EB,tempString);

  tempString.Form("tHA[%d]/F",nCrystals);
  calDataTree->Branch("tHA",tHA,tempString);

  tempString.Form("tHAErr[%d]/F",nCrystals);
  calDataTree->Branch("tHAErr",tHAErr,tempString);
  
  tempString.Form("tHB[%d]/F",nCrystals);
  calDataTree->Branch("tHB",tHB,tempString);

  tempString.Form("tHBErr[%d]/F",nCrystals);
  calDataTree->Branch("tHBErr",tHBErr,tempString);
  
  tempString.Form("TA[%d]/F",nCrystals);
  calDataTree->Branch("TA",TA,tempString);

  tempString.Form("TB[%d]/F",nCrystals);
  calDataTree->Branch("TB",TB,tempString);
  
  TString plotNameA  = "";
  TString plotTitleA = "";
  TString plotNameB  = "";
  TString plotTitleB = "";

  Float_t min = 0;
  Float_t max = 1300;

  for( Int_t i = 0 ; i < nCrystals ; i++ ){
    plotNameA.Form("hEA%d",i);
    plotTitleA.Form("hEA%d;Energy (keV);Counts",i);
    hEA[i] = new TH1F(plotNameA,plotTitleA,512,min,max);
    
    plotNameB.Form("hEB%d",i);
    plotTitleB.Form("hEB%d;QDC;Counts",i);
    hEB[i] = new TH1F(plotNameB,plotTitleB,512,min,max);
  }

  cout << " rawDataTree->GetEntries() = " << 
    rawDataTree->GetEntries() << endl;
  
  Long64_t maxEntries = rawDataTree->GetEntries();

  for (Int_t i = 0 ; i < nCrystals ; i++){
    EA[i]  = 0.;
    EB[i]  = 0.;
    tHA[i] = 0.;
    tHAErr[i] = 0.;
    tHB[i] = 0.;
    tHBErr[i] = 0.;
    TA[i]  = 0.;
    TB[i]  = 0.;
    
    // not used
    QA[i]  = 0.;
    QB[i]  = 0.;
  }
  
  Int_t chaA, chaB, cryA, cryB;
  
  Float_t QA_temp    = 0.;
  Float_t QB_temp    = 0.;
      
  Float_t phoQA_temp = 0.;
  Float_t phoQB_temp = 0.;
  
  // Calculate E,T,theta
  for( Long64_t i = 0 ; i <  maxEntries ; i++ ){
    
    rawDataTree->GetEntry(i);
    
    for (Int_t k = 0 ; k < 5 ; k++){ 

      QA_temp    = 0.;
      QB_temp    = 0.;
      
      phoQA_temp = 0.;
      phoQB_temp = 0.;
      
      // channels for A go from 0 - 4
      chaA = k;
      // channels for B go from 5 - 9
      chaB = (k+5);
      
      // crytals for A 
      cryA = Chan2ArrayA(chaA);
      // crytals for B 
      cryB = Chan2ArrayB(chaB);
      
      // to do: time calibration
      TA[cryA]  = T[chaA];
      TB[cryB]  = T[chaB];
      
      // Energy Arrays 
      // Set to use OR data after AND data
      QA_temp   = Q[chaA]-GetPedestal(chaA);
      QB_temp   = Q[chaB]-GetPedestal(chaB);

      // pedestal subtracted
      QA[cryA]  = QA_temp;
      QB[cryB]  = QB_temp;
      
      // Set to use OR data after AND data
      phoQA_temp = GetPhotopeak(chaA) - GetPedestal(chaA);
      phoQB_temp = GetPhotopeak(chaB) - GetPedestal(chaB);
      
      EA[cryA]  = QA_temp * 511./phoQA_temp;
      EB[cryB]  = QB_temp * 511./phoQB_temp;
      
      // We presume the photons interacted 
      // in the central crystal first
      // for all apart from centre crystal

      tHA[cryA] = PhotonEnergyToTheta(EA[cryA]);
      tHB[cryB] = PhotonEnergyToTheta(EB[cryB]);
      
      tHAErr[cryA] = ThetaToThetaError(tHA[cryA],chaA);
      tHBErr[cryB] = ThetaToThetaError(tHB[cryB],chaB);
    }

    // central crystals
    tHA[4] = ElectronEnergyToTheta(EA[4]);
    tHB[4] = ElectronEnergyToTheta(EB[4]);

    tHAErr[4] = ThetaToThetaError(tHA[4],2);
    tHBErr[4] = ThetaToThetaError(tHB[4],7);
    
    // Create Energy Histos
    for( Int_t j = 0 ; j < nCrystals ; j++ ) {
      
      hEA[j]->Fill(EA[j]);
      hEB[j]->Fill(EB[j]);
      
    }
    
    calDataTree->Fill();
  }
  
  calDataTree->Write();
  
  rootFileCalData->Write();
  
}

void TLab::SetPedestals(){

  // This method takes the value 
  // in the maximum bin which works
  // for data taken with OR trigger
  
  // Calibrations use the
  // second OR run pedestals
  
  cout << endl;
  cout << " Setting Pedestals " << endl;
  
  TString rawFileName;
  
  cout << endl;
  cout << " Reading " << rootFileRawName << endl;
  rootFileRawData = new TFile(rootFileRawName);
  
  TString histName = "";
  
  cout << endl;
  for(Int_t run = 0 ; run < nRuns ; run++ ){    
    for( Int_t i = 0 ; i < nChannels ; i++ ){
      histName.Form("hQ%d_%d",i,run);
      hQ[i][run] = (TH1F*)rootFileRawData->Get(histName);
      pedQ[i][run] = hQ[i][run]->GetXaxis()->
	GetBinCenter(hQ[i][run]->GetMaximumBin());
    }
  }
  
  cout << endl;
  for(Int_t run = 0 ; run < nRuns ; run++ ){
    cout << endl;
    for( Int_t i = 0 ; i < nChannels ; i++ )
      cout << " pedQ["<< i << "][" << run 
	   << "] =  " << pedQ[i][run] << endl;
  }
  cout << endl;

  rootFileRawData->Close();
  
}

Float_t TLab::GetPedestal(Int_t channel){
  return pedQ[channel][DefaultPedestalRun()]; 
}

void TLab::FitPhotopeaks(){

  cout << endl;
  cout << " Setting Photopeaks " << endl;
  
  TString rawFileName;
  
  cout << endl;
  cout << " Reading " << rootFileRawName << endl;
  rootFileRawData = new TFile(rootFileRawName);
  
  SetStyle();
	       
  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);
  TString histName = "";
  Char_t  plotName[128];
  
  Double_t maxv = 0.;
  
  TF1 *phoQfit = nullptr;
  
  for(Int_t run = 0 ; run < nRuns ; run++){
    for( Int_t i = 0 ; i < nChannels ; i++ ){
      
      histName.Form("hQ%d_%d",i,run);
      hQ[i][run] = (TH1F*)rootFileRawData->Get(histName);
      hQ[i][run]->GetXaxis()->SetRangeUser(2400,4000);
      maxv = hQ[i][run]->GetXaxis()->
	GetBinCenter(hQ[i][run]->GetMaximumBin());
      
      phoQfit = new TF1("phoQfit",
			"[0]*exp(-0.5*(((x-[1])/[2])^2))",
			maxv-250,maxv+250);
      
      phoQfit->SetLineColor(2);
      phoQfit->SetParameters(10.,3000.,100.);

      if      (i == 9)
	phoQfit->SetParameters(10.,2800.,100.);
      else if (i == 2)
	phoQfit->SetParameters(10.,3200.,100.);
      
      //phoQfit->SetParLimits(1.,2700.,3700.);
      //phoQfit->SetParLimits(2.,100.,300.);
      
      hQ[i][run]->Fit("phoQfit","RQ");
      
      sprintf(plotName,"../Plots/%d_hQ%d_%d.pdf",
	      runNumberInt,
	      i,run);
      
      canvas->SaveAs(plotName);
      
      phoQ[i][run] = phoQfit->GetParameter(1.);
      HWHM[i][run] = (phoQfit->GetParameter(2.))*Sqrt(Log(2.));

    }
  }
  
  
  for( Int_t run = 0 ; run < nRuns ; run++ ){
    cout << endl;  
    for( Int_t i = 0 ; i < nChannels ; i++ ){
      cout << " phoQ["<< i << "][" << run 
	   << "] =  " << phoQ[i][run] << endl;
    }
  }

  cout << endl;

  rootFileRawData->Close();
  
}

Float_t TLab::GetPhotopeak(Int_t channel){
  return phoQ[channel][DefaultPhotopeakRun(channel)]; 
}

Int_t TLab::DefaultPedestalRun(){  
  return 2;
}

Int_t TLab::DefaultPhotopeakRun(Int_t channel){  
  
  if(channel == 2 || channel == 7 )
    return 1;
  else  
    return 2;
}

Float_t TLab::ThetaToPhotonEnergy(Float_t theta){
  return (511./(2 - Cos(TMath::DegToRad()*theta)));
}

Float_t TLab::ThetaToElectronEnergy(Float_t theta){
  return (511. - (511./(2. - Cos(TMath::DegToRad()*theta))));
}

Float_t TLab::ThetaToThetaError(Float_t theta,
				Int_t channel){
  
  Float_t EnergyRes = 0.;
  EnergyRes = HWHM[channel][DefaultPhotopeakRun(channel)];
  EnergyRes = EnergyRes/Sqrt(511.);
  
  Float_t energy     = 0.;
  Float_t energy_max = 0.;
  Float_t energy_min = 0.;
  Float_t theta_min  = 0.;
  Float_t theta_max  = 0.;
  
  // Channels 2 and 7 correspond to central crystals
  
  // This function calculates the theta difference
  // at energy +- HWHM and returns the largest 
  // value as theta error
  if ( (channel!=2) &&
       (channel!=7) ){
    energy     = ThetaToPhotonEnergy(theta);
    energy_max = energy + EnergyRes*Sqrt(energy);
    theta_min  = theta  - PhotonEnergyToTheta(energy_max);
    energy_min = energy - EnergyRes*Sqrt(energy);
    theta_max  = PhotonEnergyToTheta(energy_min) - theta;
  }
  else{
    energy     = ThetaToElectronEnergy(theta);
    energy_max = energy + EnergyRes*Sqrt(energy);
    theta_max  = ElectronEnergyToTheta(energy_max) - theta;
    energy_min = energy - EnergyRes*Sqrt(energy);
    theta_min  = theta - ElectronEnergyToTheta(energy_min);
  }
  return (theta_max>theta_min? theta_max : theta_min);
}

Bool_t TLab::GoodTiming(Float_t time){
  
  Bool_t goodTime = kFALSE;
  
  // Unit is TDC bins
  Float_t timingRange[2];
  timingRange[0] =  800.;
  timingRange[1] = 1200.;
  
  if( time > timingRange[0] && 
      time < timingRange[1]   )
    goodTime = kTRUE;
  
  return goodTime;
}


Bool_t TLab::GoodTheta(Float_t theta){
  
  Bool_t goodTheta = kFALSE;
  
  // Asymmety is minimal below 
  // 30 degrees but allow
  // for energy resolution
  Float_t thetaRange[2];
  thetaRange[0] =  10.;
  thetaRange[1] =  170.;
  
  if( theta > thetaRange[0]  && 
      theta < thetaRange[1] )
    goodTheta = kTRUE;
  
  return goodTheta;
}

void TLab::CalculateAsymmetry(){
  
  cout << endl;
  cout << " Calculating Asymmetry  " << endl;
  
  GetThetaBinValues();
  
  // this function will calculate values
  // for the following data member:
  for(Int_t j = 0 ; j <nThBins; j++){
    for(Int_t k = 0 ; k < nPhiBins; k++){
      AsymMatrix[j][k] = 0;
    }
 }
  
  rootFileCalData = new TFile(rootFileCalName);
  
  calDataTree = (TTree*)rootFileCalData->Get("calDataTree");
  
  calDataTree->SetBranchAddress("EA",EA);
  calDataTree->SetBranchAddress("EB",EB);
  
  calDataTree->SetBranchAddress("tHA",tHA);
  calDataTree->SetBranchAddress("tHB",tHB);
  
  calDataTree->SetBranchAddress("tHAErr",tHAErr);
  calDataTree->SetBranchAddress("tHBErr",tHBErr);
  
  calDataTree->SetBranchAddress("TA",&TA);
  calDataTree->SetBranchAddress("TB",&TB);
  
  Bool_t    A[nCrystals],  B[nCrystals];
  Long64_t nA[nCrystals], nB[nCrystals];
  
  for ( Int_t i = 0 ; i < nCrystals ; i++ ){
    A[i]  = kFALSE;
    B[i]  = kFALSE;
    nA[i] = 0;
    nB[i] = 0;
  }
  
  Bool_t AB000 = kFALSE, AB090 = kFALSE,
    AB180 = kFALSE, AB270 = kFALSE;
  
  Long64_t maxEntry = calDataTree->GetEntries();
  
  Int_t nDuplicates = 0;
  
  Float_t phiB;
  Bool_t  randomisePhiB = kFALSE;
  
  Float_t minE = 450.;
  Float_t maxE = 550.;

  Float_t totEA = 0.;
  Float_t totEB = 0.;
  Int_t   thBin = -1.;
    
  for(Long64_t nEntry = 0 ; nEntry < maxEntry; nEntry++ ){
    
    calDataTree->GetEvent(nEntry);
    
    totEA = 0.;
    totEB = 0.;
    thBin = -1.;
    
    for (Int_t i = 0 ; i < nCrystals; i++){
      totEA += EA[i];
      totEB += EB[i];
    }
    
    //---------------------------
    // Event Selection
    
    if ((totEA < minE) ||
	(totEA > maxE) ||
	(totEB < minE) ||
	(totEB > maxE))
      continue;
    
    A[4] = kFALSE;
    B[4] = kFALSE;
    
    // Central crystal selections first
    // as they are essential 
    
    // assign theta bin to array A  
    // central crystal
    for (Int_t k = 0 ; k < nThBins; k++){
      if( ( GoodTheta(tHA[4]) ) &&
	  ( tHA[4] > ThMin[k] ) &&
	  ( tHA[4] < ThMax[k] )){
	A[4] = kTRUE;
	nA[4]++;
	thBin = k;
      }
    }
    
    // assign theta bin to array B
    // central crystal 
    if((A[4]) &&
       ( GoodTheta(tHB[4]) ) &&
       ( tHB[4] > ThMin[thBin] ) &&
       ( tHB[4] < ThMax[thBin] )){
      B[4] = kTRUE;
      nB[4]++;
    }
    
    // Central Crystals are always required
    if( !A[4] || !B[4] )
      continue;
    
    //  Outer crystal selections
    for (Int_t j = 0 ; j < nCrystals ; j++){
      
      if(j!=4){
	A[j] = kFALSE;  
	B[j] = kFALSE;
	
	// Conditions on summed central 
	// and summed crystal energies
	// to assign phi
	// theta bin must match central
	// crystals
	if( 
	   ( GoodTheta(tHA[j]) )     &&
	   // ( (EA[j]+EA[4]) > minE )  &&
	   // ( (EA[j]+EA[4]) < maxE )  &&
	   ( tHA[j] > ThMin[thBin] ) &&
	   ( tHA[j] < ThMax[thBin] )
	    ){
	  A[j] = kTRUE;
	  nA[j]++;
	}
	
	if( 
	   ( GoodTheta(tHB[j]) )     &&
	   // ( (EB[j]+EB[4]) > minE )  &&
	   // ( (EB[j]+EB[4]) < maxE )  &&
	   ( tHB[j] > ThMin[thBin] ) &&
	   ( tHB[j] < ThMax[thBin] )
	    ){
	  B[j] = kTRUE;
	  nB[j]++;
	}
	
      }
    } // end of: for (Int_t j = 0 ; j < nCrystals  
    
    // optional check for bugs
    if( randomisePhiB ){
      phiB = RandomLabPhi();
      
      for(Int_t r = 0 ; r < nCrystals ; r++)
	B[r] = RandomGoodLabPhi(phiB,r);
    }
    
    AB000 = kFALSE;
    AB090 = kFALSE;
    AB180 = kFALSE; 
    AB270 = kFALSE;
    
    // determine delta phi
    // for this event
    if((A[1]&&B[1])||
       (A[3]&&B[3])||
       (A[5]&&B[5])||
       (A[7]&&B[7]))
      AB000 = kTRUE;
    
    if((A[1]&&B[3])||
       (A[3]&&B[7])||
       (A[7]&&B[5])||
       (A[5]&&B[1]))
      AB090 = kTRUE;
    
    if((A[1]&&B[7])||
       (A[3]&&B[5])||
       (A[7]&&B[1])||
       (A[5]&&B[3]))
      AB180 = kTRUE;
    
    if((A[1]&&B[5])||
       (A[5]&&B[7])||
       (A[7]&&B[3])||
       (A[3]&&B[1]))
      AB270 = kTRUE;
    
    // check that only one combination
    // is satisfied
    if( (AB000 && AB090) || 
	(AB000 && AB180) ||
	(AB000 && AB270) ||
	(AB090 && AB180) ||
	(AB090 && AB270) ||
	(AB180 && AB270) ){
      nDuplicates++;
      continue;
    }
    
    if     (AB000)
      AsymMatrix[thBin][0]+=1.;
    else if(AB090)
      AsymMatrix[thBin][1]+=1.;
    else if(AB180)
      AsymMatrix[thBin][2]+=1.;
    else if(AB270)
      AsymMatrix[thBin][3]+=1.;
    
  } // end of : for(Int_t i = 0 ; i < calDa...

  if(nDuplicates!=0)
    cout << " nDuplicates = " << nDuplicates << endl;
  
  cout << endl;
  cout << " Asymmetry"<<endl;
  cout << " theta \t" << "dPhi=0 \t" 
       << "90 \t" << "180 \t" << "270" << endl;
  for (Int_t i = 0 ; i < nThBins ; i++)
    cout << " " << plotTheta[i]     << "\t" 
	 << " " << AsymMatrix[i][0] << "\t" 
	 << " " << AsymMatrix[i][1] << "\t"
	 << " " << AsymMatrix[i][2] << "\t"
	 << " " << AsymMatrix[i][3] << endl;
    
}

Float_t TLab::RandomLabPhi(){
  
  TRandom1 * rand1 = new TRandom1(); 
  Float_t phi = rand1->Uniform()*360;
  
  if     ( phi < 90 ){
    phi = 0.;
  }
  else if( phi >= 90  &&
	   phi < 180){
    phi = 90.;
  }
  else if( phi >= 180 &&
	   phi < 270){
    phi = 180.;
  }
  else if( phi >= 270 &&
	   phi < 360.){
    phi = 270.;
  }
  
  return phi;
}

Bool_t  TLab::RandomGoodLabPhi(Float_t phi, 
			       Int_t crystal){
  
  Bool_t  good = kFALSE;

  if      (phi     ==  0.  && 
	   crystal == 1){
    good = kTRUE;
  }
  else if (phi     == 90.  &&
	   crystal == 3){
    good = kTRUE;
  }
  else if (phi     == 180. &&
	   crystal == 5){
    good = kTRUE;
  }
  else if (phi     == 270. &&
	   crystal == 7){
    good = kTRUE;	
  }
  
  return good;
}


Float_t TLab::ElectronEnergyToTheta(Float_t energy){

  Float_t theta = 0.0;

  energy = 511. - energy;
  
  const Float_t m = 511.0;

  // cos(theta)
  theta = 2 - m/energy;

  theta = ACos(theta);
  theta = RadToDeg()*theta;

  return theta;
}

Float_t TLab::PhotonEnergyToTheta(Float_t energy){

  Float_t theta = 0.0;

  const Float_t m = 511.0;

  // cos(theta)
  theta = 2 - m/energy;
  theta = ACos(theta);
  theta = RadToDeg()*theta;
  
  return theta;
}

void TLab::GetThetaBinValues(){
  
  Float_t thetaBinWidth = (thetaHighEdge - thetaLowEdge)/(Float_t)nThBins;

  for (Int_t i = 0 ; i < nThBins ; i++){
    ThMin[i] = thetaLowEdge + i*thetaBinWidth;
    ThMax[i] = ThMin[i] + thetaBinWidth;
    plotTheta[i] = ThMin[i] + (ThMax[i] - ThMin[i])/2.;
  }
}

void TLab::GraphAsymmetry(Char_t option){

  SetStyle();
  GetThetaBinValues();
  
  Float_t  AsPhiDiff[nThBins];
  Float_t  AePhiDiff[nThBins];
  
  Bool_t divByUnPol = kFALSE;
  Bool_t divByF     = kFALSE;
  Bool_t correctA   = kFALSE;
  
  if     (option=='d'){
    correctA   = kTRUE;
    divByUnPol = kTRUE;
  }
  else if(option=='f'){
    correctA = kTRUE;
    divByF   = kTRUE;
  }
  
  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = N(90)/N(0) 
  Int_t   dPhiDiff = 90;

  Float_t  phi[nPhiBins];

  for (Int_t i = 0 ; i < nPhiBins ; i++)
    phi[i] = i*90.;
  
  // for delta phi graph
  Float_t phiLowEdge  = -45.0;
  Float_t phiHighEdge = 315.0;

  Float_t thetaBinWidth = (thetaHighEdge - thetaLowEdge)/(Float_t)nThBins;
  
  // Asymmetry plot range
  Float_t maxY = 3.5;
  Float_t minY = 0.5;

  if( dPhiDiff==180 )
    maxY = 6.0;
  
  // Theory curve
  Float_t aTheory[nThBins];
  Float_t aTheory1[nThBins];

  // 35.0 is result from Chloe Schoolings fits
  // perhaps an underestimate due to distribution
  // wings
  Float_t alpha1   = DegToRad()*35.0*2.355/2.;
    
  // theta half width for theory
  Float_t semiSpan = thetaBinWidth/2.*DegToRad();
  
  cout << endl;
  cout << " semiSpan = " << semiSpan*RadToDeg() << endl;
  cout << " alpha1   = " << alpha1*RadToDeg()   << endl;

  TTheory * theory =  new TTheory();
  
  // Simulation results
  Float_t aSim[nThBins]={0};
  Float_t aSimE[nThBins]={0};
  Float_t aSimTrue[nThBins]={0};
  Float_t aSimTrueE[nThBins]={0};
  
  Float_t f_aSim[nThBins]={0};
  Float_t f_aSimE[nThBins]={0};
  // Float_t f_aSimTrue[nThBins]={0};
  // Float_t f_aSimTrueE[nThBins]={0};
  
  Float_t aSimU[nThBins]={0};
  Float_t aSimUE[nThBins]={0};
  Float_t aSimUTrue[nThBins]={0};
  Float_t aSimUTrueE[nThBins]={0};
  
  Float_t AsPhiDiffD[nThBins]={0};
  Float_t AePhiDiffD[nThBins]={0};

  Float_t AsPhiDiffF[nThBins]={0};
  Float_t AePhiDiffF[nThBins]={0};

  Float_t mu[nThBins]={0};
  // To Do:
  Float_t muE[nThBins]={0};
  
  // use 90 and 270 for 90 degrees?
  Bool_t use270 = kTRUE;
  
  //  lab calculation (not theory only)
  if( option!='t' ){
    
    cout << endl;
    cout << " Calculating asymmetry values and " << endl;
    cout << " associated errors for lab data ... " << endl;
    
    CalculateAsymmetry();
    
    for (Int_t i = 0 ; i < nThBins ; i++){
      
      if (AsymMatrix[i][0] != 0){
	
	mu[i] = (AsymMatrix[i][1] - AsymMatrix[i][0]);
	mu[i] = mu[i]/(AsymMatrix[i][1] + AsymMatrix[i][0]);

	if (dPhiDiff  == 90){
	  
	  if(!use270){
	    AsPhiDiff[i] =
	      AsymMatrix[i][1]/AsymMatrix[i][0];
	    AePhiDiff[i] =
	      AsPhiDiff[i]*Sqrt((1/AsymMatrix[i][1])+(1/AsymMatrix[i][0]));
	  }
	  else{
	    AsPhiDiff[i] =
	      (AsymMatrix[i][1]+AsymMatrix[i][3])/(2*AsymMatrix[i][0]);
	    AePhiDiff[i] =
	      AsPhiDiff[i]*Sqrt((1/(AsymMatrix[i][1]+AsymMatrix[i][3]))+(1/AsymMatrix[i][0]));
	  }
	}
	if (dPhiDiff  == 180){
	  AsPhiDiff[i] =
	    AsymMatrix[i][2]/AsymMatrix[i][0];
	  AePhiDiff[i] =
	    AsPhiDiff[i]*Sqrt((1/AsymMatrix[i][2])+(1/AsymMatrix[i][0]));
	}
	if (dPhiDiff  == 270){
	  AsPhiDiff[i] =
	    AsymMatrix[i][3]/AsymMatrix[i][0];
	  AePhiDiff[i] =
	    AsPhiDiff[i]*Sqrt((1/AsymMatrix[i][3])+(1/AsymMatrix[i][0]));
	}
      }
      
    // Only plot points in axis range
      if(AsPhiDiff[i] < 0.5 || AsPhiDiff[i] > maxY)
      	AsPhiDiff[i] = 0.0; 
      
    }
 
  }// end of: if( option!='t' &...


  // theory and lab
  if( option=='t' ||
      option=='b'  ){
    
    cout << endl;
    cout << " Calculating theory curve ... " << endl;
    cout << endl;
    cout << " semiSpan = " << semiSpan*RadToDeg() << endl;
    cout << " alpha1   = " << alpha1*RadToDeg()   << endl;
    
    for(Int_t i = 0 ; i < nThBins ; i++){
      
      if( dPhiDiff == 180 ){
	aTheory[i] = 1.0;
	continue;
      }
      plotTheta[i] = plotTheta[i]*DegToRad();
      aTheory[i] = theory->rho2(plotTheta[i],semiSpan,alpha1);
      aTheory1[i] = theory->rho1(plotTheta[i],semiSpan);
      plotTheta[i] = plotTheta[i]*RadToDeg();
 
    }
    
  }// end of:  if( option=='t'


  // sim
  if( option=='a' || 
      option=='d' || 
      option=='f' || 
      option=='s' || // s,c: relics?
      option=='c'  ){
    
    cout << endl;
    cout << " Calculating theory curve " << endl;
    cout << " and simulation results ... " << endl;
    TSim *simData ;
      
    if(divByUnPol){
      
      // two argument option requires some work
      // but is functional 
      // (no sort file creations yet)
      simData = new TSim(simRun,simRunU);
      
      // unpolarised (second file) first
      simData->CalculateAsymmetryLab(simRunU);
    }
    else{
      simData = new TSim(simRun);
      simData->CalculateAsymmetryLab(simRun);
    }
    
    for(Int_t i = 0 ; i < nThBins ; i++){
      plotTheta[i] = plotTheta[i]*DegToRad();
      
      aTheory[i]   = theory->rho2(plotTheta[i],semiSpan,alpha1);            
      aTheory1[i]  = theory->rho1(plotTheta[i],semiSpan);
      plotTheta[i] = plotTheta[i]*RadToDeg();

      aSim[i]        = simData->GetAsymLab(dPhiDiff,i);
      aSimE[i]       = simData->GetAsymLabErr(dPhiDiff,i);
      aSimTrue[i]    = simData->GetAsymLabTrue(dPhiDiff,i);
      aSimTrueE[i]   = simData->GetAsymLabTrueErr(dPhiDiff,i);
      
      // divide lab data by unpolarised sim
      AsPhiDiffD[i] = (AsPhiDiff[i]/aSim[i]);
      AePhiDiffD[i] = AsPhiDiffD[i] * 
	Sqrt(AePhiDiff[i]*AePhiDiff[i]/
	     (AsPhiDiff[i]*AsPhiDiff[i]) + 
	     aSimE[i]*aSimE[i]/(aSim[i]*aSim[i]) );
      
      if(divByF){
	
	// error for acceptance from theory
	f_aSimE[i]     = aSimE[i]/aTheory1[i];
	//	f_aSimTrueE[i] = aSimTrueE[i]/aTheory1[i];
	
	// acceptance from theory & simulation
	// divide lab data by this
	// as alternative to the unpolarised
	// simulated data
	f_aSim[i]     = aSim[i]/aTheory1[i];
	//	f_aSimTrue[i] = aSimTrue[i]/aTheory1[i];
	
	AsPhiDiffF[i] = AsPhiDiff[i]/f_aSim[i];
	
	AePhiDiffF[i] = AsPhiDiffF[i] * 
	  Sqrt(AePhiDiff[i]*AePhiDiff[i]/
	       (AsPhiDiff[i]*AsPhiDiff[i]) + 
	       f_aSimE[i]*f_aSimE[i]/
	       (f_aSim[i]*f_aSim[i]) );

	cout << endl;
	cout << " AsPhiDiffF[" << i << "] = " 
	     << AsPhiDiffF[i] << endl;
	cout << " AePhiDiffF[" << i << "] = " 
	     << AePhiDiffF[i] << endl;
	  


      }
      else if(divByUnPol){
	// same values as aSim[i] etc
	// - they will be overwritten below
	aSimU[i]      = aSim[i];
	aSimUE[i]     = aSimE[i];
	aSimUTrue[i]  = aSimTrue[i];
	aSimUTrueE[i] = aSimTrueE[i];
      }
    }
    
    // divide simulated entangled/polarised 
    // by simulated unpolarised
    if(divByUnPol){
      
      // now the entangled/polarised data (first file)
      simData->CalculateAsymmetryLab(simRun);
      
      for(Int_t i = 0 ; i < nThBins ; i++){
	aSim[i]        = simData->GetAsymLab(dPhiDiff,i);
	aSimE[i]       = simData->GetAsymLabErr(dPhiDiff,i);
	aSimTrue[i]    = simData->GetAsymLabTrue(dPhiDiff,i);
	aSimTrueE[i]   = simData->GetAsymLabTrueErr(dPhiDiff,i);
	
	// ----------------------
	// Acceptance correction

	// calculate the errors first
	
	// errors for dividing by unpolarised
	aSimE[i] = (aSim[i]/aSimU[i])*
	  Sqrt( aSimUE[i]*aSimUE[i]/(aSimU[i]*aSimU[i]) + 
		aSimE[i]*aSimE[i]/(aSim[i]*aSim[i]));
	
	aSimTrueE[i] = (aSimTrue[i]/aSimUTrue[i]) *
	  Sqrt( aSimUTrueE[i]*aSimUTrueE[i]/
		( aSimUTrue[i]*aSimUTrue[i]) + 
		aSimTrueE[i]*aSimTrueE[i]/
		(aSimTrue[i]*aSimTrue[i]));
	
	
	// acceptance from unpolarised simulation
	// divide by unpolarised 
	aSim[i]      = aSim[i]/aSimU[i];
	aSimTrue[i]  = aSimTrue[i]/aSimUTrue[i];
	
	
      }
    }

  }// end of:  if( option=='t' || option=='T'
  
  cout << endl;
  cout << " ... done.  " << endl;
  cout << " Plotting the results.  " << endl;
    
  Float_t AsInt[nPhiBins], AeInt[nPhiBins];
  Float_t maxCountsPhi = 0;
  
  for ( Int_t j = 0 ; j < nPhiBins ; j++){
    
    AsInt[j] = 0.0;
    AeInt[j] = 0.0;
      
    for ( Int_t i = 0 ; i < nThBins ; i++ ){
      
      AsInt[j] += AsymMatrix[i][j];
      AeInt[j] += AsymMatrix[i][j];  
      
      }
    
    AeInt[j] = Sqrt(AeInt[j]);
    
    cout << endl;
    cout << " AsInt[" << j << "] = " << AsInt[j] << endl;
    cout << " AeInt[" << j << "] = " << AeInt[j] << endl;
    
    if( (AsInt[j]*1.1) > maxCountsPhi)
      maxCountsPhi = AsInt[j]*1.1;
  }
  TCanvas *canvas1 = new TCanvas("canvas","canvas",
				10,10,1200,800);
    
  // Axis
  TH1F *hr;

  // Graph counts vs dPhi integrated over theta
  TGraphErrors *grDPhi = new TGraphErrors(nPhiBins,
					  phi,
					  AsInt,
					  0,
					  AeInt);
  
  Char_t plotName[128];
  // dPhi Plot
  hr = canvas1->DrawFrame(phiLowEdge,0,
			  phiHighEdge,maxCountsPhi);
  hr->GetXaxis()->SetTitle("#phi (deg)");
  hr->GetYaxis()->SetTitle("N(#Delta#phi)");
  
  grDPhi->Draw("P E");
  
  sprintf(plotName,"../Plots/%d_DeltaPhi_%c.pdf",
	  runNumberInt,option);
  canvas1->SaveAs(plotName);
 
  // Asymmetry

  TGraphErrors * grAsym[5];
  
  TGraphErrors * grMu = new TGraphErrors(nThBins,plotTheta,mu,
					 0,muE);


  if( correctA ){
    if     (divByF){
      grAsym[0] = new TGraphErrors(nThBins,plotTheta,AsPhiDiff,
				   0,AePhiDiff);
      grAsym[1] = new TGraphErrors(nThBins,plotTheta,aTheory1,
				   0,0);
    }
    else if(divByUnPol){
      grAsym[0] = new TGraphErrors(nThBins,plotTheta,AsPhiDiffD,
				   0,AePhiDiffD);
      grAsym[1] = new TGraphErrors(nThBins,plotTheta,aTheory,
				   0,0);
    }
  }
  else{
    grAsym[0] = new TGraphErrors(nThBins,plotTheta,AsPhiDiff,
				 0,AePhiDiff);
    grAsym[1] = new TGraphErrors(nThBins,plotTheta,aTheory,
				 0,0);
  }
  for (Int_t k = 0; k<nThBins; k++)
    plotTheta[k] -= 1.;
  
  grAsym[2] = new TGraphErrors(nThBins,plotTheta,aSim,0,aSimE);
  
  for (Int_t k = 0; k<nThBins; k++)
    plotTheta[k] += 2.;
  
  if(option!='f')
    grAsym[3] = new TGraphErrors(nThBins,plotTheta,aSimTrue,
				 0,aSimTrueE);
  else{
    grAsym[3] = new TGraphErrors(nThBins,plotTheta,AsPhiDiffF,
				 0,AePhiDiffF);
    grAsym[4] = new TGraphErrors(nThBins,plotTheta,f_aSim,
				 0,f_aSimE);
  }

    
  grAsym[0]->SetLineColor(kBlue);
  grAsym[0]->SetMarkerColor(kBlue);
  grAsym[0]->SetFillColor(kBlue);
  grAsym[0]->SetFillStyle(3000);
  
  grAsym[1]->SetLineColor(kRed);
  grAsym[1]->SetMarkerColor(kRed);
  grAsym[1]->SetFillColor(kRed);
  
  
  grAsym[2]->SetLineColor(kGreen+2.7);
  grAsym[2]->SetMarkerColor(kGreen+2.7);
  grAsym[2]->SetFillColor(kGreen+2.7);
  
  if(option!='f'){
    grAsym[3]->SetLineColor(kGreen);
    grAsym[3]->SetMarkerColor(kGreen);
    grAsym[3]->SetFillColor(kGreen);
    grAsym[3]->SetFillStyle(3001);
  }
  else{
    grAsym[3]->SetLineColor(kBlue+2.7);
    grAsym[3]->SetMarkerColor(kBlue+2.7);
    grAsym[4]->SetLineColor(kOrange);
    grAsym[4]->SetMarkerColor(kOrange);
    grAsym[4]->SetFillColor(kOrange);
    grAsym[4]->SetFillStyle(3002);
  }
  
  cout << " here " << endl;      

  TLegend *leg =  new TLegend(0.6,0.75,0.9,0.85);
  
  TString theoryLegendTitle = " ";
  
  alpha1 = alpha1*RadToDeg();
  
  theoryLegendTitle.Form("theory curve #alpha_{#phi} = %.1f^{o}", alpha1);
  
  if(divByF)
    theoryLegendTitle.Form("theory curve ");
  
  Char_t yAxis[128];
  
  // Ratio Plot
  hr = canvas1->DrawFrame(thetaLowEdge,minY,thetaHighEdge,maxY);
  
  hr->GetXaxis()->SetTitle("#theta (deg)");

  sprintf(yAxis,"P(%d^{o})/P(0^{o})",dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);
    
  if     (option=='b'){
    leg->AddEntry(grAsym[0],"laboratory","E P");
    leg->AddEntry(grAsym[1],
		  theoryLegendTitle,"L P");
    grAsym[0]->Draw("P E");
    grAsym[1]->Draw("same P L");
  }
  else if(option=='l'){
    leg->AddEntry(grAsym[0],"laboratory","E P");
    grAsym[0]->Draw("P E");
  }
  else if(option == 't'){
    leg->AddEntry(grAsym[1],
		  theoryLegendTitle,"L P");
    grAsym[1]->Draw("L P");
  }
  else if(option=='a' ||
	  option=='d' ||
	  option=='f'){
    leg->AddEntry(grAsym[0],
		  "lab","E P");
    leg->AddEntry(grAsym[1],
		  theoryLegendTitle,"L P");
    leg->AddEntry(grAsym[2],
		  "sim ","E P");
    if(option!='f')
      leg->AddEntry(grAsym[3],
		    "sim (extra)","E P");
    else{
      leg->AddEntry(grAsym[3],
		    "lab data / f ","E P");
      leg->AddEntry(grAsym[4],
		    " f = sim data / theory","3");
    }
    // grAsym[0]->Draw("L P E");
//     grAsym[1]->Draw("same L P");
//     grAsym[2]->Draw("same L P E");
//     grAsym[3]->Draw("same L P E");

    // f
    if(option=='f'){
      grAsym[4]->Draw("3");
      // lab
      grAsym[0]->Draw("same P E");
    }
    else 
      grAsym[0]->Draw("P E");
    
    // theory
    grAsym[1]->Draw("same LP");
    // sim
    grAsym[2]->Draw("same EP");
    // lab/f
    grAsym[3]->Draw("same EP ");
    
  }
  else if(option=='c' || option=='C'){
    leg->AddEntry(grAsym[0],
		  "laboratory","E P");
    leg->AddEntry(grAsym[2],
		  "simulation (no entanglement)","P E");
    grAsym[0]->Draw("P E");
    grAsym[2]->Draw("same P");
  }
  else if(option=='s' || option=='S'){
    leg->AddEntry(grAsym[1],
		  theoryLegendTitle,"L P");
    leg->AddEntry(grAsym[2],
		  "simulation (no entanglement)","P E");
    grAsym[1]->Draw("P L");
    grAsym[2]->Draw("same P L");
  }
  
  leg->Draw();
  
  sprintf(plotName,"../Plots/%d_A_%d_%c.pdf",
	  runNumberInt,dPhiDiff,option);
  
  canvas1->SaveAs(plotName);

  // mu plot
  hr = canvas1->DrawFrame(thetaLowEdge,-0.2,
			  thetaHighEdge,0.5);
  
  hr->GetXaxis()->SetTitle("#theta (deg)");

  sprintf(yAxis,"(N(%d^{o})-N(0^{o})/(N(%d^{o})+N(0^{o})",
	  dPhiDiff,dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);

  sprintf(plotName,"../Plots/%d_Mu_%d_%c.pdf",
	  runNumberInt,dPhiDiff,option);
  
  grMu->Draw("P E");
  
  canvas1->SaveAs(plotName);

  // beep
  cout << '\a' << endl;
  
  }

void TLab::SetStyle(){
  
  TStyle     *garyStyle  = new TStyle("garyStyle","My Root Styles");
  
  const Int_t NCont = 255;
  const Int_t NRGBs = 5;
  
  // Color scheme for 2D plotting with a better defined scale 
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };          
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  
  garyStyle->SetNumberContours(NCont);
  
  // General
  // OPTIONS - FILL LINE TEXT MARKER

  garyStyle->SetFillColor(0);
  garyStyle->SetTextSize(0.05);
  
  //-----------  Canvas
  
  garyStyle->SetCanvasBorderMode(0);
  garyStyle->SetCanvasColor(kWhite);
  
  //------------- Pad
  
  garyStyle->SetPadBorderMode(0); 
  garyStyle->SetPadColor(kWhite);
  
  //Make more room for X and Y titles
  // one pad
  garyStyle->SetPadRightMargin(0.05);  //percentage
  garyStyle->SetPadLeftMargin(0.1);    //percentage
  garyStyle->SetPadBottomMargin(0.12); //percentage

  // six sub-pads
  // garyStyle->SetPadRightMargin(0.16);  //percentage
  // garyStyle->SetPadLeftMargin(0.2);    //percentage
  // garyStyle->SetPadBottomMargin(0.14); //percentage

  //----------- Histogram
  
  //Histos
  garyStyle->SetHistLineWidth(2);
  garyStyle->SetMarkerStyle(20);
  
  //  FILL CONTOURS LINE BAR 
  //  Frames
  garyStyle->SetFrameBorderMode(0);
  
  //  FILL BORDER LINE
  //  Graphs
  //  LINE ERRORS
  
  //---------  Axis 
  
  garyStyle->SetLabelFont(132,"XYZ"); 
  garyStyle->SetLabelSize(0.04,"XYZ");
  garyStyle->SetLabelOffset(0.01 ,"Y");
  
  //---------  Title
  garyStyle->SetOptTitle(0);
  
  garyStyle->SetTitleFont(132,"XYZ"); 
  garyStyle->SetTitleSize(0.05,"XYZ");

  garyStyle->SetTitleOffset(1.0,"XYZ");
  
  // 6 sub-pads
  //garyStyle->SetTitleOffset(1.6,"Y");
  
  //----------  Stats
  garyStyle->SetOptStat(0);
  garyStyle->SetOptFit(1);

  //----------  Legend
  garyStyle->SetLegendBorderSize(0);
  //garyStyle->SetLegendFont(132);
  
  gROOT->SetStyle("garyStyle");
  gROOT->ForceStyle();

}

// ------------------------------------------------------------------------------------------------
