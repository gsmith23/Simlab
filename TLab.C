#include "TLab.h"
#include "./includes.h"

#if !defined(__CINT__)
ClassImp(TLab)
#endif

// Default Constructor
TLab::TLab( ) {
  cout <<  endl;
  cout << " Default constructor. " << endl;
	cout << " lalalalala " << endl;
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
    
    nameHist.Form("hQ%d",i);
    titleHist.Form("hQ%d;QDC bin;Counts",i);
    hQ[i] = new TH1F(nameHist,titleHist,4096,0,4096);
    
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
      hQ[i]->Fill(Q[i]);
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

void TLab::MakeCalibratedDataTreeFile(){
  
  cout << endl;
  cout << " Making calibrated data tree " << endl;

  SetPedestals();
  
  // To Do 
  // - fit photopeaks and hardcode mean
  // in the function, or even better...
  // - automate the fitting
  // at present all the photopeaks are set 
  // to 3000.
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
  
  tempString.Form("EA[%d]/F",nChannels);
  calDataTree->Branch("EA",EA,tempString);

  tempString.Form("EB[%d]/F",nChannels);
  calDataTree->Branch("EB",EB,tempString);

  tempString.Form("tHA[%d]/F",nChannels);
  calDataTree->Branch("tHA",tHA,tempString);
  
  tempString.Form("tHB[%d]/F",nChannels);
  calDataTree->Branch("tHB",tHB,tempString);

  tempString.Form("TA[%d]/F",nChannels);
  calDataTree->Branch("TA",TA,tempString);

  tempString.Form("TB[%d]/F",nChannels);
  calDataTree->Branch("TB",TB,tempString);
  
  TString plotNameA  = "";
  TString plotTitleA = "";
  TString plotNameB  = "";
  TString plotTitleB = "";

  Float_t min = 0;
  Float_t max = 1300;

  for( Int_t i = 0 ; i < nChannels ; i++ ){
    plotNameA.Form("hEA%d",i);
    plotTitleA.Form("hEA%d;Energy (keV);Counts",i);
    hEA[i] = new TH1F(plotNameA,plotTitleA,512,min,max);
    
    plotNameB.Form("hEB%d",i);
    plotTitleB.Form("hEB%d;QDC;Counts",i);
    hEB[i] = new TH1F(plotNameB,plotTitleB,512,min,max);
  }

  cout << " rawDataTree->GetEntries() = " << 
    rawDataTree->GetEntries() << endl;
  

  
  // Calculate E,T,theta
  for( Int_t i = 0 ; i < rawDataTree->GetEntries() ; i++ ){
    
    rawDataTree->GetEntry(i);
    
    // To Do: map channels to 
    // crystal number scheme
    // for ease of use in 
    // asymmetry calculation
    
    Int_t kA, kB;
    
    for (Int_t k = 0 ; k < 5 ; k ++){ 
      
      kA = k;
      kB = (k+5);
      
      QA[k]  = Q[kA];
      QB[k]  = Q[kB];
      
      TA[k]  = T[kA];
      TB[k]  = T[kB];
      
      // Array A
      EA[k]  = (Q[kA]-pedQ[kA])*511./(phoQ[kA]-pedQ[kA]) ;
      EB[k]  = (Q[kB]-pedQ[kB])*511./(phoQ[kB]-pedQ[kB]) ; ;
      
      // for all apart from centre crystal
      tHA[k] = PhotonEnergyToTheta(EA[kA]);
      tHB[k] = PhotonEnergyToTheta(EB[kB]);
      
      // !!!!!!  check these    !!!!!!!!
      // !!!!!! are the correct !!!!!!!!
      // !!!!!!   channels     !!!!!!!!
      // central crystals
      tHA[2] = ElectronEnergyToTheta(EA[2]);
      tHB[2] = ElectronEnergyToTheta(EB[2]);
      
    }
    
    // Create Energy Histos
    for( Int_t j = 0 ; j < nChannels ; j++ ) {
      
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
  for( Int_t i = 0 ; i < nChannels ; i++ ){
    histName.Form("hQ%d",i);
    hQ[i] = (TH1F*)rootFileRawData->Get(histName);
    pedQ[i] = hQ[i]->GetXaxis()->GetBinCenter(hQ[i]->GetMaximumBin());
  }
  
  for( Int_t i = 0 ; i < nChannels ; i++ )
    cout << " pedQ["<< i << "] =  " << pedQ[i] << endl;
  
  rootFileRawData->Close();
  
}

Float_t TLab::GetPedestal(Int_t channel){
  return pedQ[channel]; 
}

void TLab::SetPhotopeaks(){

 // for (Int_t i = 0 ; i < nChannels ; i++)
 //   phoQ[i] = 3000.;
Double_t phoQ[10]={2871.,
		   3076.,
		   3034., //central crystal A
      		   2489.,
      		   2570.,
      		   2741.,
      		   2917.,
                   3011., //central crystal B
                   2475.,
                   3161.};

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
  
  // example variables to use for checking
  // hits 
  Bool_t A[nChannels];
  Bool_t B[nChannels];
  
  for ( Int_t i = 0 ; i < nChannels ; i++ ){
    A[i] = kFALSE;
    B[i] = kFALSE;
  }
  
  
  
  // To Do:
  // Implement routine to calculate 
  // Asymmetry A(dPhi) = N(dPhi)/N(0)
  // in a given theta range
  // by counting four hit combinations
  
  
  cout << endl;
  cout << " entries = " << calDataTree->GetEntries() << endl;
  
  Int_t nA4 = 0; 

  for(Int_t i = 0 ; i < calDataTree->GetEntries(); i++ ){
    
    calDataTree->GetEvent(i);
    
    A[4] = kFALSE;  
    
    //cout << " EA[4] = " << EA[4] << endl;
    
    if( ( EA[4] > 200.  ) &&
	( EA[4] < 400.  )){
      
      A[4] = kTRUE;
      
      nA4++;
      
      //cout << " A[4] = " << A[4] << endl;    
    }
    
  }
  
  cout << " nA4 = " << nA4 << endl;     

  cout << endl;
  cout << " Here is where the Asymmetry be calculated " << endl;
  cout << " A(" << dPhi << ") in the range " 
       << minTh << " < #theta < " << maxTh << endl;
  cout << endl;
  
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
  const Int_t nBins = 4;

  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = P(90)/P(0) 
  Int_t   dPhiDiff = 90;
  
  Float_t  theta[nBins];
  Float_t  thetaRange[nBins][2];
  Float_t  As090[nBins];
  Float_t  Ae090[nBins];
  
  // Theta range 
  Float_t thetaLowEdge  = 30.;
  thetaLowEdge  = 0.;
  Float_t thetaHighEdge = 130.;
  thetaHighEdge = 360.;
  
  Float_t thetaBinWidth = (thetaHighEdge - thetaLowEdge)/(Float_t)nBins;
  
  // Axis
  TH1F * hr;
  hr = canvas->DrawFrame(thetaLowEdge,0.5,thetaHighEdge,2.5);
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
