#include "TSim.h"
#include "./includes.h"

#if !defined(__CINT__)
ClassImp(TSim)
#endif

TSim::TSim(){
}

TSim::TSim(TString fileNumber){
  
  cout << endl;
  cout << " Constructing TSim object " << endl;
    
  rootFileRawName = "../Data/sim" + fileNumber;
  rootFileSortName = "../Data/sort" + fileNumber;
  
  rootFileRawName = rootFileRawName + ".root";
  rootFileSortName = rootFileSortName + ".root";

  Initialise();
  
  //!! temporary 
  // SetAsymmetry(fileNumber);
  

}

TSim::~TSim(){
}

//----------------------------------------------

/** Public member functions *********/


Bool_t TSim::SortedROOTFileExists(){
  TFile *file = TFile::Open(rootFileSortName);

  return file;
}

void TSim::Loop()
{
  
   Long64_t nentries = simDataTree->GetEntriesFast();

   for (Long64_t i = 0 ; i < nentries ; i++) {
     simDataTree->GetEntry(i);
           
   }
}


Float_t TSim::ElectronEnergyToTheta(Double_t nrg){
  return RadToDeg()*ACos(2. - 511./(511. - nrg));
}

Float_t TSim::PhotonEnergyToTheta(Double_t nrg){
  return RadToDeg()*ACos(2. - 511./nrg);;
}

Float_t TSim::ThetaToPhotonEnergy(Float_t theta){
  return (511./(2 - Cos(TMath::DegToRad()*theta)));
}

Float_t TSim::ThetaToElectronEnergy(Float_t theta){
  return (511. - (511./(2. - Cos(TMath::DegToRad()*theta))));
}



Int_t TSim::GetThetaBin(Float_t theta){
  
  Int_t bin = -1;

  for (Int_t i = 0 ; i < nThbins ; i++){
    if ((theta > ThMin[i])&&(theta < ThMax[i])){
      bin = i;
    }
  }
  
  return bin;
}




void TSim::Initialise(){
  
  cout << endl;
  cout << " Connecting Branches " << endl;
  
  theFile = new TFile(rootFileRawName);
  //name of the Tree?
  simDataTree = (TTree*)theFile->Get("Tangle2");
  
  /*simDataTree->SetBranchAddress("edep0", &edep0, &b_edep0);
  simDataTree->SetBranchAddress("edep1", &edep1, &b_edep1);
  simDataTree->SetBranchAddress("edep2", &edep2, &b_edep2);
  simDataTree->SetBranchAddress("edep3", &edep3, &b_edep3);
  simDataTree->SetBranchAddress("edep4", &edep4, &b_edep4);
  simDataTree->SetBranchAddress("edep5", &edep5, &b_edep5);
  simDataTree->SetBranchAddress("edep6", &edep6, &b_edep6);
  simDataTree->SetBranchAddress("edep7", &edep7, &b_edep7);
  simDataTree->SetBranchAddress("edep8", &edep8, &b_edep8);
  simDataTree->SetBranchAddress("edep9", &edep9, &b_edep9);
  simDataTree->SetBranchAddress("edep10", &edep10, &b_edep10);
  simDataTree->SetBranchAddress("edep11", &edep11, &b_edep11);
  simDataTree->SetBranchAddress("edep12", &edep12, &b_edep12);
  simDataTree->SetBranchAddress("edep13", &edep13, &b_edep13);
  simDataTree->SetBranchAddress("edep14", &edep14, &b_edep14);
  simDataTree->SetBranchAddress("edep15", &edep15, &b_edep15);
  simDataTree->SetBranchAddress("edep16", &edep16, &b_edep16);
  simDataTree->SetBranchAddress("edep17", &edep17, &b_edep17);
  simDataTree->SetBranchAddress("nb_Compt0", &nb_Compt0, &b_nb_Compt0);
  simDataTree->SetBranchAddress("nb_Compt1", &nb_Compt1, &b_nb_Compt1);
  simDataTree->SetBranchAddress("nb_Compt2", &nb_Compt2, &b_nb_Compt2);
  simDataTree->SetBranchAddress("nb_Compt3", &nb_Compt3, &b_nb_Compt3);
  simDataTree->SetBranchAddress("nb_Compt4", &nb_Compt4, &b_nb_Compt4);
  simDataTree->SetBranchAddress("nb_Compt5", &nb_Compt5, &b_nb_Compt5);
  simDataTree->SetBranchAddress("nb_Compt6", &nb_Compt6, &b_nb_Compt6);
  simDataTree->SetBranchAddress("nb_Compt7", &nb_Compt7, &b_nb_Compt7);
  simDataTree->SetBranchAddress("nb_Compt8", &nb_Compt8, &b_nb_Compt8);
  simDataTree->SetBranchAddress("nb_Compt9", &nb_Compt9, &b_nb_Compt9);
  simDataTree->SetBranchAddress("nb_Compt10", &nb_Compt10, &b_nb_Compt10);
  simDataTree->SetBranchAddress("nb_Compt11", &nb_Compt11, &b_nb_Compt11);
  simDataTree->SetBranchAddress("nb_Compt12", &nb_Compt12, &b_nb_Compt12);
  simDataTree->SetBranchAddress("nb_Compt13", &nb_Compt13, &b_nb_Compt13);
  simDataTree->SetBranchAddress("nb_Compt14", &nb_Compt14, &b_nb_Compt14);
  simDataTree->SetBranchAddress("nb_Compt15", &nb_Compt15, &b_nb_Compt15);
  simDataTree->SetBranchAddress("nb_Compt16", &nb_Compt16, &b_nb_Compt16);
  simDataTree->SetBranchAddress("nb_Compt17", &nb_Compt17, &b_nb_Compt17);
  simDataTree->SetBranchAddress("XposA_1st", &XposA_1st, &b_XposA_1st);
  simDataTree->SetBranchAddress("YposA_1st", &YposA_1st, &b_YposA_1st);
  simDataTree->SetBranchAddress("ZposA_1st", &ZposA_1st, &b_ZposA_1st);
  simDataTree->SetBranchAddress("XposA_2nd", &XposA_2nd, &b_XposA_2nd);
  simDataTree->SetBranchAddress("YposA_2nd", &YposA_2nd, &b_YposA_2nd);
  simDataTree->SetBranchAddress("ZposA_2nd", &ZposA_2nd, &b_ZposA_2nd);
  simDataTree->SetBranchAddress("XposB_1st", &XposB_1st, &b_XposB_1st);
  simDataTree->SetBranchAddress("YposB_1st", &YposB_1st, &b_YposB_1st);
  simDataTree->SetBranchAddress("ZposB_1st", &ZposB_1st, &b_ZposB_1st);
  simDataTree->SetBranchAddress("XposB_2nd", &XposB_2nd, &b_XposB_2nd);
  simDataTree->SetBranchAddress("YposB_2nd", &YposB_2nd, &b_YposB_2nd);
  simDataTree->SetBranchAddress("ZposB_2nd", &ZposB_2nd, &b_ZposB_2nd);
  simDataTree->SetBranchAddress("ThetaA", &ThetaA, &b_ThetaA);
  simDataTree->SetBranchAddress("PhiA", &PhiA, &b_PhiA);
  simDataTree->SetBranchAddress("ThetaB", &ThetaB, &b_ThetaB);
  simDataTree->SetBranchAddress("PhiB", &PhiB, &b_PhiB);
  simDataTree->SetBranchAddress("dPhi", &dPhi, &b_dPhi);
  */
  simDataTree->SetBranchAddress("crystal0", &crystal0, &b_crystal0);
  simDataTree->SetBranchAddress("crystal1", &crystal1, &b_crystal1);
  simDataTree->SetBranchAddress("crystal2", &crystal2, &b_crystal2);
  simDataTree->SetBranchAddress("crystal3", &crystal3, &b_crystal3);
  simDataTree->SetBranchAddress("crystal4", &crystal4, &b_crystal4);
  simDataTree->SetBranchAddress("crystal5", &crystal5, &b_crystal5);
  simDataTree->SetBranchAddress("crystal6", &crystal6, &b_crystal6);
  simDataTree->SetBranchAddress("crystal7", &crystal7, &b_crystal7);
  simDataTree->SetBranchAddress("crystal8", &crystal8, &b_crystal8);
  simDataTree->SetBranchAddress("crystal9", &crystal9, &b_crystal9);
  simDataTree->SetBranchAddress("crystal10", &crystal10, &b_crystal10);
  simDataTree->SetBranchAddress("crystal11", &crystal11, &b_crystal11);
  simDataTree->SetBranchAddress("crystal12", &crystal12, &b_crystal12);
  simDataTree->SetBranchAddress("crystal13", &crystal13, &b_crystal13);
  simDataTree->SetBranchAddress("crystal14", &crystal14, &b_crystal14);
  simDataTree->SetBranchAddress("crystal15", &crystal15, &b_crystal15);
  simDataTree->SetBranchAddress("crystal16", &crystal16, &b_crystal16);
  simDataTree->SetBranchAddress("crystal17", &crystal17, &b_crystal17);
  simDataTree->SetBranchAddress("TrueEvent", &TrueEvent, &b_TrueEvent);
  simDataTree->SetBranchAddress("cosA", &cosA, &b_cosA);
  simDataTree->SetBranchAddress("cosB", &cosB, &b_cosB);
}

Int_t  TSim::SortEvents(TString fileNumber){

  //=========================
  // Initialise variables &
  // connect tree and leaves
  
  TString rootFileRawName;
  TString rootFileSortName;
  
  rootFileRawName = "../Data/sim" + fileNumber;
  rootFileSortName = "../Data/sort" + fileNumber;
  
  rootFileRawName = rootFileRawName + ".root";
  rootFileSortName = rootFileSortName + ".root";
  
  TFile* inputFile = new TFile(rootFileRawName);
  TFile* outputFile = new TFile(rootFileSortName,"recreate");

  cout << endl;
  cout << " SortEvents: " << endl;
  cout << endl;
  cout << " as you requested, I will now sort the  " << endl;
  cout << " events in the root file for you and create" << endl;
  cout << " a new root file containing additional variables " << endl;
  cout << " to use in the PET analysis. " << endl;
  cout << endl;
  
  cout << endl;
  cout << " Input File : " << rootFileRawName  << endl;
  cout << " Output File: " << rootFileSortName << endl;
  
  //==========================
  // read in output from 
  // Geant4 simulation

  TTree* simDataTree =(TTree*)inputFile->Get("Tangle2");

  //write to

  sortDataTree = new TTree("sortDataTree",
			   "Sorted simulation events");

  //============================================
  //--------------------------------------------
  // Variable declarations

  //Energy drom the simulation will be smeared with a Gaussian
  //Set the HWHM for each crystal
  
  Double_t HWHM = 25.; //average from lab data

  //energy resolution coefficient k, such thatt sigma(Energy) = k*Sqrt(Energy)
  for (Int_t i = 0; i < nCrystals; i++){
    sigmaA[i] = HWHM/(Sqrt(2*Log(2.)*511));
    sigmaB[i] = HWHM/(Sqrt(2*Log(2.)*511));
  }
  
  //-------------------------------
  // Initial array initialisations
  for (Int_t i = 0; i < nCrystals; i++){

    //smeared energy
    EA[i] = 0; 
    EB[i] = 0;

    //angle calculated from exact energy deposit = 'energy theta'
    etHA[i] = 0;
    //angle calculated from the smeared energy deposit = 'lab theta'
    ltHA[i] = 0;
    //angle calculated from exact energy - a*sigma
    mintHAErr[i] = 0.;
    //angle calculated from exact energy + a*sigma
    maxtHAErr[i] = 0.;
    
    etHB[i] = 0.;
    ltHB[i] = 0.;
    mintHBErr[i] = 0.;
    maxtHBErr[i] = 0.;

    //angle after first Compton scattering retrieved from the simulation = 'sim theta'
    simtHA[i] = 180;
    simtHB[i] = 180;

    //1 if central crystal was hit first, 0 otherwise
    centrFirst[i] = 0;
  }


  simDataTree->SetBranchAddress("crystal0", &crystal0, &b_crystal0);
  simDataTree->SetBranchAddress("crystal1", &crystal1, &b_crystal1);
  simDataTree->SetBranchAddress("crystal2", &crystal2, &b_crystal2);
  simDataTree->SetBranchAddress("crystal3", &crystal3, &b_crystal3);
  simDataTree->SetBranchAddress("crystal4", &crystal4, &b_crystal4);
  simDataTree->SetBranchAddress("crystal5", &crystal5, &b_crystal5);
  simDataTree->SetBranchAddress("crystal6", &crystal6, &b_crystal6);
  simDataTree->SetBranchAddress("crystal7", &crystal7, &b_crystal7);
  simDataTree->SetBranchAddress("crystal8", &crystal8, &b_crystal8);
  simDataTree->SetBranchAddress("crystal9", &crystal9, &b_crystal9);
  simDataTree->SetBranchAddress("crystal10", &crystal10, &b_crystal10);
  simDataTree->SetBranchAddress("crystal11", &crystal11, &b_crystal11);
  simDataTree->SetBranchAddress("crystal12", &crystal12, &b_crystal12);
  simDataTree->SetBranchAddress("crystal13", &crystal13, &b_crystal13);
  simDataTree->SetBranchAddress("crystal14", &crystal14, &b_crystal14);
  simDataTree->SetBranchAddress("crystal15", &crystal15, &b_crystal15);
  simDataTree->SetBranchAddress("crystal16", &crystal16, &b_crystal16);
  simDataTree->SetBranchAddress("crystal17", &crystal17, &b_crystal17);
  simDataTree->SetBranchAddress("TrueEvent", &TrueEvent, &b_TrueEvent);
  simDataTree->SetBranchAddress("cosA", &cosA, &b_cosA);
  simDataTree->SetBranchAddress("cosB", &cosB, &b_cosB);


  //============================================
  // New variables
  TString tempString = "";

  tempString.Form("EA[%d]/D", nCrystals);
  sortDataTree->Branch("EA", EA, tempString);

  tempString.Form("EB[%d]/D", nCrystals);
  sortDataTree->Branch("EB", EB, tempString);

  tempString.Form("etHA[%d]/F",nCrystals);
  sortDataTree->Branch("etHA",etHA,tempString);

  tempString.Form("ltHA[%d]/F",nCrystals);
  sortDataTree->Branch("ltHA",ltHA,tempString);

  tempString.Form("maxtHAErr[%d]/F",nCrystals);
  sortDataTree->Branch("maxtHAErr",maxtHAErr,tempString);
  tempString.Form("mintHAErr[%d]/F",nCrystals);
  sortDataTree->Branch("mintHAErr",mintHAErr,tempString);
  
  tempString.Form("etHB[%d]/F",nCrystals);
  sortDataTree->Branch("etHB",etHB,tempString);

  tempString.Form("ltHB[%d]/F",nCrystals);
  sortDataTree->Branch("ltHB",ltHB,tempString);

  tempString.Form("maxtHBErr[%d]/F",nCrystals);
  sortDataTree->Branch("maxtHBErr",maxtHBErr,tempString);
  tempString.Form("mintHBErr[%d]/F",nCrystals);
  sortDataTree->Branch("mintHBErr",mintHBErr,tempString);

  tempString.Form("simtHA[%d]/F",nCrystals);
  sortDataTree->Branch("simtHA", simtHA, tempString);

  tempString.Form("simtHB[%d]/F",nCrystals);
  sortDataTree->Branch("simtHB", simtHB, tempString);

  tempString.Form("centrFirst[%d]/F",nCrystals);
  sortDataTree->Branch("centrFirst", centrFirst, tempString);

TRandom *rand = new TRandom();

  //============
  // EVENT LOOP
 Float_t sigmaPar = 2.5;

  Long64_t nEvents = simDataTree->GetEntries();
  
  for(Int_t i = 0 ; i < nEvents ; i++){ 
    simDataTree->GetEvent(i);
    
    Double_t CrystEnergyDep[18] = {crystal0,crystal1,crystal2,
				   crystal3,crystal4,crystal5,
				   crystal6,crystal7,crystal8,
				   crystal9,crystal10,crystal11,
				   crystal12,crystal13,crystal14,
				   crystal15,crystal16,crystal17};

    for (Int_t j = 0 ; j < nCrystals ; j++){
      
      EA[j] = rand->Gaus(CrystEnergyDep[j],sigmaA[j]*Sqrt(CrystEnergyDep[j]));
      EB[j] = rand->Gaus(CrystEnergyDep[j+nCrystals], sigmaB[j]*Sqrt(CrystEnergyDep[j+nCrystals]));
      
      etHA[j] = PhotonEnergyToTheta(CrystEnergyDep[j]);
      etHB[j] = PhotonEnergyToTheta(CrystEnergyDep[j+nCrystals]);
      ltHA[j] = PhotonEnergyToTheta(EA[j]);
      ltHB[j] = PhotonEnergyToTheta(EB[j]);

      // max/min defined as sigmaPar* sigma away from the mean
      maxtHAErr[j] = PhotonEnergyToTheta(CrystEnergyDep[j] - (sigmaPar*sigmaA[j]*Sqrt(CrystEnergyDep[j])));
      mintHAErr[j] = PhotonEnergyToTheta(CrystEnergyDep[j] + (sigmaPar*sigmaA[j]*Sqrt(CrystEnergyDep[j])));
      maxtHBErr[j] = PhotonEnergyToTheta(CrystEnergyDep[j+nCrystals] - (sigmaPar*sigmaB[j]*Sqrt(CrystEnergyDep[j+nCrystals])));
      mintHBErr[j] = PhotonEnergyToTheta(CrystEnergyDep[j+nCrystals] + (sigmaPar*sigmaB[j]*Sqrt(CrystEnergyDep[j+nCrystals])));
 

      simtHA[j] = ACos(cosA)*RadToDeg();
      simtHB[j] = ACos(cosB)*RadToDeg();

      centrFirst[j] = TrueEvent;
    }

    etHA[4] = ElectronEnergyToTheta(CrystEnergyDep[4]);
    etHB[4] = ElectronEnergyToTheta(CrystEnergyDep[13]);
    ltHA[4] = ElectronEnergyToTheta(EA[4]);
    ltHB[4] = ElectronEnergyToTheta(EB[4]);
    maxtHAErr[4] = ElectronEnergyToTheta(CrystEnergyDep[4] + (sigmaPar*sigmaA[4]*Sqrt(CrystEnergyDep[4])));
    mintHAErr[4] = ElectronEnergyToTheta(CrystEnergyDep[4] - (sigmaPar*sigmaA[4]*Sqrt(CrystEnergyDep[4])));
    maxtHBErr[4] = ElectronEnergyToTheta(CrystEnergyDep[13] + (sigmaPar*sigmaB[4]*Sqrt(CrystEnergyDep[13])));
    mintHBErr[4] = ElectronEnergyToTheta(CrystEnergyDep[13] - (sigmaPar*sigmaB[4]*Sqrt(CrystEnergyDep[13])));				 

    // introduce additional conditions
    /*Int_t nHits = 0;

    for (Int_t j = 0 ; j < nCrystals; j++){
      if (EA[j] > 0) nHits++;
      if (EB[j] > 0) nHits++;}

    if (nHits > 3){
      sortDataTree->Fill();
      }*/

    sortDataTree->Fill();
  } // end of: for(Int_t i = 0 ; i < nEvents     
    
  //=========================
  sortDataTree->Write();
      
  
  outputFile->Write();
  
  cout << endl;
  cout << " the new root file has been written " << endl;
  cout << endl;
  
  outputFile->Close();

  return 0;
  
} //Closes the program



Bool_t TSim::GoodTheta(Float_t theta){
  
  Bool_t goodTheta = kFALSE;
  
  if( theta >= 10. &&
      theta <  170.)
    goodTheta = kTRUE;
  
  return goodTheta;
}

Float_t TSim::CrystalToPhi(Int_t crystal){
  
  Float_t crystalToPhi[9] = {-1.,0.,-1.,270.,-1.,90.,-1.,180.,-1.};

  Float_t phi = crystalToPhi[crystal];

  return phi;
}

void TSim::SetAsymmetry(TString inputFileNumber){
  
  cout << endl;
  cout << " Setting asymmetry " << endl;
  
  this->GetAsymmetry(inputFileNumber);
  
}

void TSim::GetThetaBinValues(){
  Float_t thetaLowEdge  = 10.;
  Float_t thetaHighEdge = 170.;

  Float_t thetaBinWidth = (thetaHighEdge - thetaLowEdge)/(Float_t)nThbins;

  for (Int_t i = 0 ; i < nThbins ; i++){
    ThMin[i] = thetaLowEdge + i* thetaBinWidth;
    ThMax[i] = thetaLowEdge + (i+1)* thetaBinWidth;
    plotTheta[i] = thetaLowEdge + (i+0.5)* thetaBinWidth;
  }
}

Int_t TSim::GetAsymmetry(TString inputFileNumber){

  GetThetaBinValues();

  //difference 'sim theta' - 'energy theta'
  TH2F* histThD = new TH2F("histThD","theta difference",100,50,150,250,-120,130);

  //difference 'lab theta' - 'energy theta'
  TH2F* hLEdiff = new TH2F("hLEdiff","lab - energy theta",100,50,150,120,-50,70);

  //difference 'max/min lab theta' - 'energy theta'
  TH2F* hMEdiff = new TH2F("hMEdiff","max/min - energy theta",100,50,150,120,-50,70);

  cout << endl;
  cout << " Getting asymmetry " << endl;
  
  this->SetStyle();

  Int_t inputFileInt = inputFileNumber.Atoi();
  
  TString plotName;
  plotName = "../Plots/Asym_" + inputFileNumber;

  plotName = plotName + ".pdf";

  TString inputFileName = "../Data/sort" + inputFileNumber;
  inputFileName  = inputFileName + ".root";
    
  cout << endl;
  cout << " Input File : " << inputFileName  << endl;
  cout << endl;

  TFile* inputFile = new TFile(inputFileName);

  TTree* sortDataTree=(TTree*)inputFile->Get("sortDataTree");

  sortDataTree->SetBranchAddress("EA",EA);
  sortDataTree->SetBranchAddress("EB",EB);
  
  sortDataTree->SetBranchAddress("etHA",etHA);
  sortDataTree->SetBranchAddress("etHB",etHB);
  sortDataTree->SetBranchAddress("ltHA",ltHA);
  sortDataTree->SetBranchAddress("ltHB",ltHB);
  
  sortDataTree->SetBranchAddress("mintHAErr",mintHAErr);
  sortDataTree->SetBranchAddress("mintHBErr",mintHBErr);
  sortDataTree->SetBranchAddress("maxtHAErr",maxtHAErr);
  sortDataTree->SetBranchAddress("maxtHBErr",maxtHBErr);
  
  sortDataTree->SetBranchAddress("simtHA",&simtHA);
  sortDataTree->SetBranchAddress("simtHB",&simtHB);
  sortDataTree->SetBranchAddress("centrFirst", &centrFirst);

  n000 = 0;
  n090 = 0;
  n180 = 0;
  n270 = 0;

  for(Int_t j = 0 ; j <nThbins; j++){
    for(Int_t k = 0 ; k < 4; k++){
      AsymMatrix[j][k] =0;}
  }

  Long64_t nEntries = sortDataTree->GetEntries();

  //Event loop
  for (Int_t ientry = 0 ; ientry < nEntries ; ientry++){
    sortDataTree->GetEvent(ientry);

    Int_t thBin = -1;
    Float_t phiA = -1;
    Float_t phiB = -1;
    Int_t indexA = -1;
    Int_t indexB = -1;

    
    if (GetThetaBin(etHA[4]) == GetThetaBin(etHB[4])){
      thBin = GetThetaBin(etHA[4]);

      for (Int_t i = 0 ; i < nCrystals ; i++){
	if ((i!=4)&&(thBin>-1)&&(etHA[i]<ThMax[thBin])&&(etHA[i]>ThMin[thBin])){
	  phiA = CrystalToPhi(i);
	  indexA = i ;}
	if ((i!=4)&&(thBin>-1)&&(etHB[i]<ThMax[thBin])&&(etHB[i]>ThMin[thBin])){
	  phiB = CrystalToPhi(i);
	  indexB = i;}	
      }
    }

    if ((phiA > -1)&&(phiB > -1)){
      
      histThD->Fill(etHA[4],simtHA[4]-etHA[4]);
      histThD->Fill(etHB[4],simtHB[4]-etHB[4]);
      
      hLEdiff->Fill(etHA[4],ltHA[4] - etHA[4]);
      hLEdiff->Fill(etHA[indexA],ltHA[indexA] - etHA[indexA]);
      hLEdiff->Fill(etHB[4],ltHB[4] - etHB[4]);
      hLEdiff->Fill(etHB[indexB],ltHB[indexB] - etHB[indexB]);
      
      hMEdiff->Fill(etHA[4],mintHAErr[4] - etHA[4]);
      hMEdiff->Fill(etHA[indexA],mintHAErr[indexA] - etHA[indexA]);
      hMEdiff->Fill(etHB[4],mintHBErr[4] - etHB[4]);
      hMEdiff->Fill(etHB[indexB],mintHBErr[indexB] - etHB[indexB]);
      hMEdiff->Fill(etHA[4],maxtHAErr[4] - etHA[4]);
      hMEdiff->Fill(etHA[indexA],maxtHAErr[indexA] - etHA[indexA]);
      hMEdiff->Fill(etHB[4],maxtHBErr[4] - etHB[4]);
      hMEdiff->Fill(etHB[indexB],maxtHBErr[indexB] - etHB[indexB]);
      
      Float_t phiDiff = phiB - phiA;
      if (phiDiff == 0.) {
	AsymMatrix[thBin][0] += 1;
	n000++;}
      if ((phiDiff == 90)||(phiDiff == -270)){
	AsymMatrix[thBin][1] += 1;
	n090++;}
      if ((phiDiff == 180)||(phiDiff == -180)){
	AsymMatrix[thBin][2] += 1;
	n180++;}
      if ((phiDiff == -90)||(phiDiff == 270)){
	AsymMatrix[thBin][3] += 1;
	n270++;}
    }
    
  }//end of: for(Int_t ientry...
  
  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = P(90)/P(0) 
  Int_t   dPhiDiff = 90;

  //Calculating ratios for desired dPhiDiff
  Float_t AsPhiDiff[nThbins] = {0.};
  Float_t AePhiDiff[nThbins] = {0.};

  
  for (Int_t i = 0 ; i < nThbins ; i++){
    if (AsymMatrix[i][0] != 0){
      if (dPhiDiff  == 90){
	AsPhiDiff[i] = (AsymMatrix[i][1]+AsymMatrix[i][3])/(2*AsymMatrix[i][0]);
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/(AsymMatrix[i][1]+AsymMatrix[i][3]))+(1/AsymMatrix[i][0]));
      }
      if (dPhiDiff  == 180){
	AsPhiDiff[i] = AsymMatrix[i][2]/AsymMatrix[i][0];
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix[i][2])+(1/AsymMatrix[i][0]));
      }
      if (dPhiDiff  == 270){
	AsPhiDiff[i] = AsymMatrix[i][3]/AsymMatrix[i][0];
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix[i][3])+(1/AsymMatrix[i][0]));
      }
    }
  }
  
 
  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);

  TGraphErrors *grAsym =
    new TGraphErrors(nThbins,plotTheta,AsPhiDiff,0,AePhiDiff);

  grAsym->SetLineColor(kBlue);
  grAsym->SetMarkerColor(kBlue);


  TH1F *hr;

  hr = canvas->DrawFrame(10,0.5,170,3);
  hr->GetXaxis()->SetTitle("#theta (deg)");

  Char_t plotN[128];
  Char_t yAxis[128];
  
  sprintf(yAxis,"P(%d^{o})/P(0^{o})",dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);
  grAsym->Draw("P E");

  sprintf(plotN,"../Plots/A_%d_%d.pdf",inputFileInt,dPhiDiff);
  
  canvas->SaveAs(plotN);

  hLEdiff->Draw("colz");
  hLEdiff->GetXaxis()->SetTitle("energy #theta (deg)");
  hLEdiff->GetYaxis()->SetTitle("lab #theta - energy #theta (deg)");

  sprintf(plotN,"../Plots/histLEdiff_%d.pdf",inputFileInt);
  canvas->SaveAs(plotN);
   
  histThD->Draw("colz");
  histThD->GetXaxis()->SetTitle("energy #theta (deg)");
  histThD->GetYaxis()->SetTitle("sim #theta - energy #theta (deg)");
  
  sprintf(plotN,"../Plots/histThDiff_%d.pdf",inputFileInt);
  canvas->SaveAs(plotN);

  hMEdiff->Draw("colz");
  hMEdiff->GetXaxis()->SetTitle("energy #theta (deg)");
  hMEdiff->GetYaxis()->SetTitle("max/min lab #theta - energy #theta (deg)");

  sprintf(plotN,"../Plots/histME_%d.pdf",inputFileInt);
  canvas->SaveAs(plotN);
  return 0;

}

void TSim::SetStyle(){

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
  
  //-----------  Canvas
  
  garyStyle->SetCanvasBorderMode(0);
  garyStyle->SetCanvasColor(kWhite);
  
  //------------- Pad
  
  garyStyle->SetPadBorderMode(0); 
  garyStyle->SetPadColor(kWhite);
  
  //Make more room for X and Y titles
  garyStyle->SetPadRightMargin(0.16);  //percentage
  garyStyle->SetPadLeftMargin(0.2);    //percentage
  garyStyle->SetPadBottomMargin(0.14); //percentage
  
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
  garyStyle->SetTitleOffset(1.6,"Y");
  
  //----------  Stats
  garyStyle->SetOptStat(0);
  garyStyle->SetOptFit(0);
  
  //----------  Legend

  gROOT->SetStyle("garyStyle");
  gROOT->ForceStyle();

}

// ------------------------------------------------------------------------------------------------
