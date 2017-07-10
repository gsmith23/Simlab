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
      else if ( eventNumber > (nOR1+nAND))
	hQ[i][2]->Fill(Q[i]);
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
  
  SetPhotopeaks();
  
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
  
  tempString.Form("tHB[%d]/F",nCrystals);
  calDataTree->Branch("tHB",tHB,tempString);

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
    tHB[i] = 0.;
    TA[i]  = 0.;
    TB[i]  = 0.;
    
    // not used
    QA[i]  = 0.;
    QB[i]  = 0.;
  }
  
  Int_t chaA, chaB, cryA, cryB;
  
  // Calculate E,T,theta
  for( Long64_t i = 0 ; i <  maxEntries ; i++ ){
    
    rawDataTree->GetEntry(i);
    
    // To Do: map channels to 
    // crystal number scheme
    // for ease of use in 
    // asymmetry calculation
    
    for (Int_t k = 0 ; k < 5 ; k ++){ 

      // channels for A go from 0 - 4
      chaA = k;
      // channels for B go from 5 - 9
      chaB = (k+5);
      
      // crytals for A 
      cryA = Chan2ArrayA(chaA);
      // crytals for B 
      cryB = Chan2ArrayB(chaB);
      
      // pedestal subtracted

      QA[cryA]  = Q[chaA] - pedQ[chaA][0];
      QB[cryB]  = Q[chaB] - pedQ[chaB][0];

      
      // to do: time calibration
      TA[cryA]  = T[chaA];
      TB[cryB]  = T[chaB];
      
      // Energy Arrays 

      EA[cryA]  = (Q[chaA]-pedQ[chaA][0])*511./(phoQ[chaA]-pedQ[chaA][0]) ;
      EB[cryB]  = (Q[chaB]-pedQ[chaB][0])*511./(phoQ[chaB]-pedQ[chaB][0]) ; ;
      
      if ( i == 100 ){
	cout << endl;
	cout << "Q[" << chaA <<"]    = " << Q[chaA] << endl;
	cout << "pedQ[" << chaA <<"] = " << pedQ[chaA][0] << endl;
	cout << "phoQ[" << chaA <<"] = " << phoQ[chaA] << endl;
	cout << "EA[" << chaA <<"]   = " << EA[chaA] << endl;

      }

      // We presume the photons interacted 
      // in the central crystal first
      
      // for all apart from centre crystal
      tHA[cryA] = PhotonEnergyToTheta(EA[cryA]);
      tHB[cryB] = PhotonEnergyToTheta(EB[cryB]);
<<<<<<< HEAD
      
      // central crystals
      tHA[4] = ElectronEnergyToTheta(EA[4]);
      tHB[4] = ElectronEnergyToTheta(EB[4]);
      
=======
>>>>>>> upstream/master
    }
    
    // central crystals
    tHA[4] = ElectronEnergyToTheta(EA[4]);
    tHB[4] = ElectronEnergyToTheta(EB[4]);
    
    // Create Energy Histos
    for( Int_t j = 0 ; j < nCrystals ; j++ ) {
      
      hEA[j]->Fill(EA[j]);
      hEB[j]->Fill(EB[j]);
      
    }
    
    if(EA[4] > .200 && EB[4] > .200)
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
  for( Int_t i = 0 ; i < nChannels ; i++ ){

    for(Int_t run = 0 ; run < nRuns ; run++ ){    
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
      cout << " pedQ["<< i << "][" << run << "] =  " << pedQ[i][run] << endl;
  }
  cout << endl;

  rootFileRawData->Close();
  
}

Float_t TLab::GetPedestal(Int_t channel){
  return pedQ[channel][0]; 
}

void TLab::SetPhotopeaks(){

   cout << endl;
  cout << " Setting Photopeaks " << endl;
  
  TString rawFileName;
  
  cout << endl;
  cout << " Reading " << rootFileRawName << endl;
  rootFileRawData = new TFile(rootFileRawName);
  
  TString histName = "";
  Double_t HWHM[10]={0.};
  Double_t maxv = 0.;
  
  cout << endl;
  for( Int_t i = 0 ; i < nChannels ; i++ ){
    histName.Form("hQ%d_0",i);
    hQ[i][0] = (TH1F*)rootFileRawData->Get(histName);
    hQ[i][0]->GetXaxis()->SetRangeUser(3000,4000);
    maxv = hQ[i][0]->GetXaxis()->GetBinCenter(hQ[i][0]->GetMaximumBin());
    TF1 *phoQfit = new TF1("phoQfit","[0]*exp(-0.5*(((x-[1])/[2])^2))",maxv-300,maxv+300);
    phoQfit->SetParameters(10.,3000.,100.,0.,0.);
    phoQfit->SetParLimits(1.,3000.,3600.);
    phoQfit->SetParLimits(2.,100.,300.);
    hQ[i][0]->Fit("phoQfit","R");
    TF1 *fit = hQ[i][0]->GetFunction("phoQfit");
    phoQ[i] = fit->GetParameter(1.);
    HWHM[i]=(fit->GetParameter(2.))*TMath::Sqrt(TMath::Log(2.));
    }
  
  cout << endl;
  for( Int_t i = 0 ; i < nChannels ; i++ )
    cout << " phoQ["<< i << "] =  " << phoQ[i] << endl;

  cout << endl;

  rootFileRawData->Close();
  
   Float_t phoQ_temp[10] = {3347.,
			    3515.,
			    3225., //central crystal A
			    3161.,
			    3395.,
			    3148.,
			    3309.,
			    3458., //central crystal B
			    3412.,
			    3263.};
   

  
  cout << endl;
  for( Int_t i = 0 ; i < nChannels ; i++ ){
    phoQ[i] = phoQ_temp[i];
     cout << " phoQ["<< i << "] =  " << phoQ[i] << endl;
   }
   cout << endl;
}

Float_t TLab::GetPhotopeak(Int_t channel){
  return phoQ[channel]; 
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
			      Float_t maxTh
			      ){
  
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
  
  Bool_t A000 = kFALSE, A090 = kFALSE, A180 = kFALSE, A270 = kFALSE,
    B000 = kFALSE, B090 = kFALSE, B180 = kFALSE, B270  = kFALSE;
  
  Bool_t AB000 = kFALSE, AB090 = kFALSE, AB180 = kFALSE;
  
  Double_t n000 = 0., n090 = 0., n180 = 0.;

  
  Long64_t maxEntry = calDataTree->GetEntries();
  
  //  maxEntry = 1000000;
  
  // !!calculate
  Float_t thRes = 10.;
  
  minTh = minTh - thRes;
  maxTh = maxTh + thRes;
  
  cout << endl;
  cout << " Calculating Asymmetry  " << endl;
  cout << " A(" << dPhi << ") in the range " 
       << minTh << " < #theta < " << maxTh << endl;
  cout << endl;

  for(Long64_t i = 0 ; i < maxEntry; i++ ){
    
    calDataTree->GetEvent(i);
  
    A000 = kFALSE, A090 = kFALSE, A180 = kFALSE, A270 = kFALSE;
    B000 = kFALSE, B090 = kFALSE, B180 = kFALSE, B270  = kFALSE;
    AB000 = kFALSE, AB090 = kFALSE, AB180 = kFALSE;
    
    for (Int_t j = 0 ; j < nCrystals ; j++){
      
      A[j] = kFALSE;  
      B[j] = kFALSE;  
      
      if( ( tHA[j] > minTh  ) &&
	  ( tHA[j] < maxTh  )){
	A[j] = kTRUE;
	nA[j]++;
      }
      
      if( ( tHB[j] > minTh  ) &&
	  ( tHB[j] < maxTh  )){
	B[j] = kTRUE;
	nB[j]++;
      }
    } // end of: for (Int_t j = 0 ; j < nCrystals  
    
    // Central Crystals are always required
    
    if( !A[4] || !B[4])
      continue;
    
    if((A[1]&&B[1])||
       (A[3]&&B[3])||
       (A[5]&&B[5])||
       (A[7]&&B[7])){
      AB000 = kTRUE;
      n000++;
    }
    
    if((A[1]&&B[3])||
       (A[3]&&B[7])||
       (A[7]&&B[5])||
       (A[5]&&B[1])){
      AB090 = kTRUE;
      n090++;
    }
    
    if((A[1]&&B[7])||
       (A[3]&&B[5])||
       (A[7]&&B[1])||
       (A[5]&&B[3])){
      AB180 = kTRUE;
      n180++;
    }
    
    if(AB000 && AB090){
      cout << endl;
      cout << " AB000 && AB090 " << endl;
      cout << endl;
    }
    
    
  } // end of : for(Int_t i = 0 ; i < calDa...
  
  
  Asym = (Float_t)n090/n000;
  
  cout << endl;
  cout << " n000 = " << n000 <<  endl;
  cout << " n090 = " << n090 <<  endl;
  cout << " n180 = " << n180 <<  endl;
  
  
}

Float_t TLab::ElectronEnergyToTheta(Float_t energy){

  Float_t theta = 0.0;

  energy = 511. - energy;
  
  const Float_t m = 511.0;

  // cos(theta)
  theta = 2 - m/energy;

  theta = TMath::ACos(theta);
  theta = TMath::RadToDeg()*theta;

  return theta;
}

Float_t TLab::PhotonEnergyToTheta(Float_t energy){

  Float_t theta = 0.0;

  const Float_t m = 511.0;

  // cos(theta)
  theta = 2 - m/energy;

  theta = TMath::ACos(theta);
  theta = TMath::RadToDeg()*theta;

  return theta;
}



void TLab::GraphAsymmetry(Char_t option){

  SetStyle();

  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);

  //const Int_t nBins = 10;
  const Int_t nBins = 5;

  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = P(90)/P(0) 
  Int_t   dPhiDiff = 90;
  
  Float_t  theta[nBins];
  Float_t  thetaRange[nBins][2];
  Float_t  As090[nBins];
  Float_t  Ae090[nBins];
  
  // Theta range 
  Float_t thetaLowEdge  = 50.;
  //thetaLowEdge  = 0.;
  Float_t thetaHighEdge = 150.;
  //thetaHighEdge = 180.;
  
  Float_t thetaBinWidth = (thetaHighEdge - thetaLowEdge)/(Float_t)nBins;
  
  // Axis
  TH1F * hr;
  hr = canvas->DrawFrame(thetaLowEdge,0.5,thetaHighEdge,3.5);
  hr->GetXaxis()->SetTitle("#theta (deg)");
  
  hr->GetYaxis()->SetTitle("P(90^{o})/P(0^{o})");
  
  // Theory curve
  Float_t aTheory[nBins];
  // half resolution in dPhi 
  // !!to do - access alpha1 from user input
  Float_t alpha1   = DegToRad()*22.5;
  // half resolution in theta
  Float_t semiSpan = thetaBinWidth/2.*DegToRad();
  
  cout << endl;
  cout << " semiSpan = " << semiSpan*RadToDeg() << endl;
  cout << " alpha1   = " << alpha1*RadToDeg()   << endl;

  TTheory * theory =  new TTheory();
  
  // Simulation results
  
  Float_t aSim[nBins];
  Float_t aSimE[nBins];
  
  for(Int_t i = 0 ; i < nBins ; i++){
    thetaRange[i][0]  = thetaLowEdge + thetaBinWidth*i;
    thetaRange[i][1]  = thetaRange[i][0] + thetaBinWidth;
    theta[i] = (thetaRange[i][0] + thetaRange[i][1])/2.;
  }

  // lab only
  if( option!='t' && option!='T'){
    
    cout << endl;
    cout << " Calculating asymmetry values and " << endl;
    cout << " associated errors for lab data ... " << endl;
    
    for(Int_t i = 0 ; i < nBins ; i++){
      
      CalculateAsymmetry(dPhiDiff,thetaRange[i][0],
			 thetaRange[i][1]);
      
      As090[i] = Asym;
	
      Ae090[i] = AsymErr;

      cout << endl;
      cout << " A(90,"<< theta[i] <<")/A(0,"<< theta[i] <<") = "  
	   << As090[i] << " ± " << Ae090[i] << endl;
    }   
    
  }// end of: if( option!='t' &...

  // theory and lab
  if( option=='t' || option=='T' ||
      option=='b' || option=='B' ){
    
    cout << endl;
    cout << " Calculating theory curve ... " << endl;
    
    for(Int_t i = 0 ; i < nBins ; i++){
      theta[i] = theta[i]*DegToRad();
      
      
      aTheory[i] = theory->rho2(theta[i],semiSpan,alpha1);

      theta[i] = theta[i]*RadToDeg();
      
      cout << " theta[" << i << "] = " << theta[i]  << endl;
      cout << " semiSpan = " << semiSpan*RadToDeg() << endl;
      cout << " alpha1   = " << alpha1*RadToDeg()   << endl;
    
    }
    
  }// end of:  if( option=='t' || option=='T'

  // sim
  if( option=='a' || option=='A' ||
      option=='s' || option=='S' ||
      option=='c' || option=='C' ){
    
    cout << endl;
    cout << " Calculating theory curve " << endl;
    cout << " and simulation results ... " << endl;
    
    for(Int_t i = 0 ; i < nBins ; i++){
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
    
  TGraphErrors *grAsym[3];
  grAsym[0] =  new TGraphErrors(nBins,theta,As090,0,Ae090);
  grAsym[1] =  new TGraphErrors(nBins,theta,aTheory,0,0);
  
  grAsym[2] =  new TGraphErrors(nBins,theta,aSim,0,aSimE);
  
  grAsym[0]->SetLineColor(kBlue);
  grAsym[0]->SetMarkerColor(kBlue);
  grAsym[1]->SetLineColor(kRed);
  grAsym[1]->SetMarkerColor(kRed);
  grAsym[2]->SetLineColor(kGreen+2);
  grAsym[2]->SetMarkerColor(kGreen+2);

  TLegend * leg =  new TLegend(0.6,0.8,0.9,0.85);
 
  TString theoryLegendTitle = "theory curve #alpha_{#Delta#phi} = 22.5^{o}";
  
  if     (option=='b' || option=='B'){
    leg->AddEntry(grAsym[0],"laboratory","E P");
    leg->AddEntry(grAsym[1],
		  "theoryLegendTitle","L P");
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
  
  Char_t plotName[128];
  
  sprintf(plotName,"../Plots/A_%d.pdf",runNumberInt);
  
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
  garyStyle->SetOptFit(0);

  //----------  Legend
  garyStyle->SetLegendBorderSize(0);
  //garyStyle->SetLegendFont(132);
  
  gROOT->SetStyle("garyStyle");
  gROOT->ForceStyle();

}

// ------------------------------------------------------------------------------------------------
