#include "TLab.h"
#include "./includes.h"

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

TLab::TLab(TString runNumber,
	   TString simNumber) {
  simData = new TSim(simNumber);
  SetFilenames(runNumber);
}

// option for use with multiple raw files
// per run
TLab::TLab(TString runNumber,
	   TString fileNumStart,
	   TString fileNumFinish) {
  
  cout << endl;
  cout << " fileNumStart  = " << fileNumStart  << endl;
  cout << " fileNumFinish = " << fileNumFinish << endl;
  
  fileNumStart  = "_" + fileNumStart;
  fileNumFinish = "_" + fileNumFinish;
  
  SetFilenames(runNumber,
	       fileNumStart,
	       fileNumFinish);

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

void TLab::SetFilenames(TString runNumber,
			TString fileNumStart,
			TString fileNumFinish) {

  cout << endl;
  cout << " TLab object has been created " << endl;

  runNumberInt = runNumber.Atoi();
  
  cout << endl;
  cout << " Run Number    = " << runNumberInt  << endl;

  textFileName    = runNumber + fileNumStart;
  rootFileRawName = runNumber + fileNumStart;  
  rootFileCalName = runNumber + fileNumStart;

  textFileName    = textFileName    + fileNumFinish;
  rootFileRawName = rootFileRawName + fileNumFinish;  
  rootFileCalName = rootFileCalName + fileNumFinish;

  textFileName    = textFileName    + ".txt";
  rootFileRawName = rootFileRawName + ".root";  
  rootFileCalName = rootFileCalName + ".root";
  
  textFileName    = "../Data/run" + textFileName;
  rootFileRawName = "../Data/run" + rootFileRawName;
  rootFileCalName = "../Data/cal" + rootFileCalName;

  cout << " rootFileRawName = " << rootFileRawName << endl;
}

/** Public member functions *********/

Bool_t TLab::RawROOTFileExists(){
  
  TFile *file = TFile::Open(rootFileRawName);
    
  return file;
}

void TLab::MakeRawDataTreeFile(){

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
    
    for(Int_t i = index ; i < (index+5) ; i++ ){
      
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
  
  // To Do 
  // - automate the fitting,
  // at present all the photopeaks 
  // are hardcoded from manual fits
  
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
  //maxEntries = 10000;

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
  // in the maximum bin which will work
  // for data taken with OR trigger
  
  // To Do: take pedestals for AND trigger
  // data from OR data file/s
  
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
  return pedQ[channel][2]; 
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
  Char_t plotName[128];

  Double_t maxv = 0.;

  TF1 *phoQfit = nullptr;
  
  for(Int_t run = 0 ; run < nRuns ; run++){
    for( Int_t i = 0 ; i < nChannels ; i++ ){
      
      histName.Form("hQ%d_%d",i,run);
      hQ[i][run] = (TH1F*)rootFileRawData->Get(histName);
      
      hQ[i][run]->GetXaxis()->SetRangeUser(2500,4000);
      
      maxv = hQ[i][run]->GetXaxis()->
	GetBinCenter(hQ[i][run]->GetMaximumBin());
      
      phoQfit = new TF1("phoQfit",
			"[0]*exp(-0.5*(((x-[1])/[2])^2))",
			maxv-200,maxv+200);
      
      phoQfit->SetLineColor(2);
      phoQfit->SetParameters(10.,3000.,100.,0.,0.);
      phoQfit->SetParLimits(1.,2700.,3700.);
      //phoQfit->SetParLimits(2.,100.,300.);
      
      hQ[i][run]->Fit("phoQfit","RQ");
    
      sprintf(plotName,"../Plots/Run_%d_hQ%d_%d.pdf",
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
  
  // !! Use the post AND data result for now
  // ie use HWHM[channel][2]
  //
  // Note that since we have 
  // using namespace TMath in
  // include.h, TMath:: can be omitted
  Float_t EnergyRes = HWHM[channel][2]/Sqrt(511.);
  Float_t energy = 0.;
  Float_t energy_max = 0.;
  Float_t energy_min = 0.;
  Float_t theta_min = 0.;
  Float_t theta_max = 0.;
  
  //Channels 2 and 7 correspond to central crystals
  if ((channel!=2)&&(channel!=7)){
    energy = ThetaToPhotonEnergy(theta);
    energy_max = energy + EnergyRes*TMath::Sqrt(energy);
    theta_min = theta -  PhotonEnergyToTheta(energy_max);
    energy_min = energy - EnergyRes*TMath::Sqrt(energy);
    theta_max = PhotonEnergyToTheta(energy_min)- theta;}
  else{
    energy = ThetaToElectronEnergy(theta);
    energy_max = energy + EnergyRes*TMath::Sqrt(energy);
    theta_max = ElectronEnergyToTheta(energy_max) - theta;
    energy_min = energy - EnergyRes*TMath::Sqrt(energy);
    theta_min = theta - ElectronEnergyToTheta(energy_min);
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

void TLab::CalculateAsymmetry(Int_t   dPhi,
			      Float_t minTh,
			      Float_t maxTh){
  
  // this function will calculate values
  // for the following data members:
  Asym    = 1.0;
  AsymErr = 0.0;
  
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
  
  // 
  Bool_t A[nCrystals];
  Bool_t B[nCrystals];
  
  Long64_t nA[nCrystals], nB[nCrystals];
  
  for ( Int_t i = 0 ; i < nCrystals ; i++ ){
    A[i]  = kFALSE;
    B[i]  = kFALSE;
    nA[i] = 0.;
    nB[i] = 0.;
  }
  Bool_t AB000 = kFALSE, AB090 = kFALSE,
    AB180 = kFALSE, AB270 = kFALSE;
  
  Double_t n000 = 0., n090 = 0.,
    n180 = 0, n270 = 0;
  
  Long64_t maxEntry = calDataTree->GetEntries();

  cout << endl;
  cout << " Calculating Asymmetry  " << endl;
  cout << " A(" << dPhi << ") in the range " 
       << minTh << " < #theta < " << maxTh << endl;
  
  // !!!???!!!
  Float_t thRes = 0.0;

  Int_t nDuplicates = 0;
  
  for(Long64_t i = 0 ; i < maxEntry; i++ ){
    
    calDataTree->GetEvent(i);
    
    AB000 = kFALSE, AB090 = kFALSE, 
      AB180 = kFALSE, AB270 = kFALSE;
    
    for (Int_t j = 0 ; j < nCrystals ; j++){
      
      A[j] = kFALSE;  
      B[j] = kFALSE;
      
      if( ( GoodTheta(tHA[j]) ) &&
	  ( tHA[j] > minTh - thRes  ) &&
	  ( tHA[j] < maxTh + thRes )){
	A[j] = kTRUE;
	
	nA[j]++;
      }
      
      if( ( GoodTheta(tHB[j]) ) &&
	  ( tHB[j] > minTh - thRes ) &&
	  ( tHB[j] < maxTh + thRes )){
	B[j] = kTRUE;
	nB[j]++;
	
      }
    } // end of: for (Int_t j = 0 ; j < nCrystals  
    
    // Central Crystals are always required
    if( !A[4] || !B[4] )
      continue;
    
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

    // To Do
    // check that 511 keV was 
    // deposited in two crytals
    // per array
    
    if     (AB000)
      n000++;
    else if(AB090)
      n090++;
    else if(AB180)
      n180++;
    else if(AB270)
      n270++;
    
  } // end of : for(Int_t i = 0 ; i < calDa...
  
  AsymPhi[0]    = n000;
  AsymPhiErr[0] = Sqrt(n000);
  
  AsymPhi[1]    = n090;
  AsymPhiErr[1] = Sqrt(n090);
  
  AsymPhi[2]    = n180;
  AsymPhiErr[2] = Sqrt(n180);
  
  AsymPhi[3]    = n270;
  AsymPhiErr[3] = Sqrt(n270);
  
  if     (dPhi==90){
    Asym = (Float_t)n090/n000;
    AsymErr = (Float_t)Asym*Sqrt((1/n090)+(1/n000));
  }
  else if(dPhi==180){
    Asym = (Float_t)n180/n000;
    AsymErr = (Float_t)Asym*Sqrt((1/n180)+(1/n000));
  }
  else if(dPhi==270){
    Asym = (Float_t)n270/n000;
    AsymErr = (Float_t)Asym*Sqrt((1/n270)+(1/n000));
  }
  
  cout << endl;
  cout << " n000 = " << n000 <<  endl;
  cout << " n090 = " << n090 <<  endl;
  cout << " n180 = " << n180 <<  endl;
  cout << " n270 = " << n270 <<  endl;
  cout << endl;
  cout << " nDuplicates = " << nDuplicates << endl;
  
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



void TLab::GraphAsymmetry(Char_t option){

  SetStyle();

  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);

  // old less-general method
  Float_t  theta[nThBins];
  Float_t  thetaRange[nThBins][2];
  
  Float_t  AsPhiDiff[nThBins];
  Float_t  AePhiDiff[nThBins];
  
  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = P(90)/P(0) 
  Int_t   dPhiDiff = 90;
  
  // new more-general method
  Float_t  As[nThBins][nPhiBins];
  Float_t  Ae[nThBins][nPhiBins];

  Float_t  phi[nPhiBins];

  for (Int_t i = 0 ; i < nPhiBins ; i++)
    phi[i] = i*90.;
  
  // Theta range 
  Float_t thetaLowEdge  = 10.;
  Float_t thetaHighEdge = 170.;
  
  Float_t phiLowEdge  = -45.0;
  Float_t phiHighEdge = 315.0;

  Float_t thetaBinWidth = (thetaHighEdge - thetaLowEdge)/(Float_t)nThBins;
  
  // Axis
  TH1F * hr;
  
  Float_t maxY = 2.5;
  Float_t minY = 0.5;

  if(dPhiDiff==180)
    maxY = 6.0;
  
  
  // Theory curve
  Float_t aTheory[nThBins];
  // half resolution in dPhi 
  // !!to do - access alpha1 from user input
  Float_t alpha1   = DegToRad()*22.5;
  alpha1   = DegToRad()*37.0;
  // half resolution in theta
  
  //!!!!!!!
  //Float_t thetaBinWidthTemp = thetaBinWidth + 3.0;
  Float_t thetaBinWidthTemp = thetaBinWidth;
  Float_t semiSpan = thetaBinWidthTemp/2.*DegToRad();
  
  cout << endl;
  cout << " semiSpan = " << semiSpan*RadToDeg() << endl;
  cout << " alpha1   = " << alpha1*RadToDeg()   << endl;

  TTheory * theory =  new TTheory();
  
  // Simulation results
  
  Float_t aSim[nThBins];
  Float_t aSimE[nThBins];
  
  for(Int_t i = 0 ; i < nThBins ; i++){
    thetaRange[i][0]  = thetaLowEdge + thetaBinWidth*i;
    thetaRange[i][1]  = thetaRange[i][0] + thetaBinWidth;
    theta[i] = (thetaRange[i][0] + thetaRange[i][1])/2.;
  }

  // lab only
  if( option!='t' && option!='T'){
    
    cout << endl;
    cout << " Calculating asymmetry values and " << endl;
    cout << " associated errors for lab data ... " << endl;
    
    for(Int_t i = 0 ; i < nThBins ; i++){
      
      CalculateAsymmetry(dPhiDiff,thetaRange[i][0],
			 thetaRange[i][1]);
      
      AsPhiDiff[i] = Asym;
      AePhiDiff[i] = AsymErr;
      
      // Only plot points in axis range
      if(AsPhiDiff[i] < 0.5 || AsPhiDiff[i] > maxY)
	AsPhiDiff[i] = 0.0;
      
      // Counts vs dPhi 
      for (Int_t j = 0 ; j < nPhiBins ;  j++){
	
	As[i][j] = AsymPhi[j];
	Ae[i][j] = AsymPhiErr[j];
	
      }   
      
    }
  }// end of: if( option!='t' &...

  // theory and lab
  if( option=='t' || option=='T' ||
      option=='b' || option=='B' ){
    
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
      
      theta[i] = theta[i]*DegToRad();
      
      aTheory[i] = theory->rho2(theta[i],semiSpan,alpha1);

      theta[i] = theta[i]*RadToDeg();
      
      cout << endl;
      cout << " theta[" << i << "] = " << theta[i]  << endl;
      
    }
    
  }// end of:  if( option=='t' || option=='T'

  // sim
  if( option=='a' || option=='A' ||
      option=='s' || option=='S' ||
      option=='c' || option=='C' ){
    
    cout << endl;
    cout << " Calculating theory curve " << endl;
    cout << " and simulation results ... " << endl;
    
    for(Int_t i = 0 ; i < nThBins ; i++){
      theta[i] = theta[i]*DegToRad();
      
      aTheory[i] = theory->rho2(theta[i],semiSpan,alpha1);
      
      aSim[i]    = simData->GetAsymm(i);
      aSimE[i]   = simData->GetAsymErr(i);
            
      theta[i] = theta[i]*RadToDeg();
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
      
      AsInt[j] += As[i][j];
      AeInt[j] += Ae[i][j]*Ae[i][j];  
      
      }
    
    AeInt[j] = Sqrt(AeInt[j]);
    
    cout << endl;
    cout << " AsInt[" << j << "] = " << AsInt[j] << endl;
    cout << " AeInt[" << j << "] = " << AeInt[j] << endl;
    
    if( (AsInt[j]*1.1) > maxCountsPhi)
      maxCountsPhi = AsInt[j]*1.1;;
    
  }

  // Graph counts vs dPhi integrated over theta
  TGraphErrors * grDPhi = new TGraphErrors(nPhiBins,
					   phi,
					   AsInt,
					   0,
					   AeInt);
  
  // dPhi Plot
  hr = canvas->DrawFrame(phiLowEdge,0,phiHighEdge,maxCountsPhi);
  hr->GetXaxis()->SetTitle("#phi (deg)");
  hr->GetYaxis()->SetTitle("N(#Delta#phi)");
  
  grDPhi->Draw("P E");
  
  Char_t plotName[128];
  
  sprintf(plotName,"../Plots/DeltaPhi_%d.pdf",runNumberInt);
  canvas->SaveAs(plotName);
  
  
  TGraphErrors *grAsym[3];
  grAsym[0] =  new TGraphErrors(nThBins,theta,AsPhiDiff,0,AePhiDiff);
  grAsym[1] =  new TGraphErrors(nThBins,theta,aTheory,0,0);
  
  grAsym[2] =  new TGraphErrors(nThBins,theta,aSim,0,aSimE);
  
  grAsym[0]->SetLineColor(kBlue);
  grAsym[0]->SetMarkerColor(kBlue);
  grAsym[1]->SetLineColor(kRed);
  grAsym[1]->SetMarkerColor(kRed);
  grAsym[2]->SetLineColor(kGreen+2);
  grAsym[2]->SetMarkerColor(kGreen+2);

  TLegend * leg =  new TLegend(0.6,0.8,0.9,0.9);
 
  TString theoryLegendTitle = " ";
  
  alpha1 = alpha1*RadToDeg();
  theoryLegendTitle.Form("theory curve #alpha_{#Delta#phi} = %f^{o}", alpha1);
  
  Char_t yAxis[128];
  
  // Ratio Plot
  hr = canvas->DrawFrame(thetaLowEdge,minY,thetaHighEdge,maxY);
  hr->GetXaxis()->SetTitle("#theta (deg)");

  sprintf(yAxis,"P(%d^{o})/P(0^{o})",dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);

  
  if     (option=='b' || option=='B'){
    leg->AddEntry(grAsym[0],"laboratory","E P");
    leg->AddEntry(grAsym[1],
		  theoryLegendTitle,"L P");
    grAsym[0]->Draw("P E");
    grAsym[1]->Draw("same P L");
  }
  else if(option=='l' || option=='L'){
    leg->AddEntry(grAsym[0],"laboratory","E P");
    grAsym[0]->Draw("P E");
  }
  else if(option == 't' || option == 'T'){
    leg->AddEntry(grAsym[1],
		  theoryLegendTitle,"L P");
    grAsym[1]->Draw("L P");
  }
  else if(option=='a' || option=='A'){
    leg->AddEntry(grAsym[0],
		  "laboratory","E P");
    leg->AddEntry(grAsym[1],
		  theoryLegendTitle,"L P");
    leg->AddEntry(grAsym[2],
		  "simulation (no entanglement)","P E");
    grAsym[0]->Draw("P E");
    grAsym[1]->Draw("same L P");
    grAsym[2]->Draw("same P");
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
  
  sprintf(plotName,"../Plots/A_%d_%d.pdf",runNumberInt,dPhiDiff);
  
  canvas->SaveAs(plotName);
  
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
