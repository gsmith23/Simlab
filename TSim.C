#include "TSim.h"

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
  
}

//option for use with two input files
TSim::TSim(TString fileNumber1, TString fileNumber2){
  
  cout << endl;
  cout << " Constructing TSim object " << endl;
    
  rootFileRawName1 = "../Data/sim" + fileNumber1;
  rootFileSortName1 = "../Data/sort" + fileNumber1;
  
  rootFileRawName1 = rootFileRawName1 + ".root";
  rootFileSortName1 = rootFileSortName1 + ".root";

  rootFileRawName2 = "../Data/sim" + fileNumber2;
  rootFileSortName2 = "../Data/sort" + fileNumber2;
  
  rootFileRawName2 = rootFileRawName2 + ".root";
  rootFileSortName2 = rootFileSortName2 + ".root";

  //  Initialise();
  
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
    if ( 
	(theta > ThMin[i])  &&
	(theta < ThMax[i])
	 ){
      bin = i;
    }
  }
  
  return bin;
}


void TSim::Initialise(){
  
  cout << endl;
  cout << " Connecting Branches " << endl;
  
  theFile = new TFile(rootFileRawName);
  //  name of the Tree?
  simDataTree = (TTree*)theFile->Get("Tangle2");

  simDataTree->SetBranchAddress("edep0", &edep0, &b_edep0);
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
  simDataTree->SetBranchAddress("ThetaA_1st", &ThetaA_1st, &b_ThetaA_1st);
  simDataTree->SetBranchAddress("PhiA_1st", &PhiA_1st, &b_PhiA_1st);
  simDataTree->SetBranchAddress("ThetaB_1st", &ThetaB_1st, &b_ThetaB_1st);
  simDataTree->SetBranchAddress("PhiB_1st", &PhiB_1st, &b_PhiB_1st);
  simDataTree->SetBranchAddress("dPhi_1st", &dPhi_1st, &b_dPhi_1st);
  simDataTree->SetBranchAddress("ThetaA_2nd", &ThetaA_2nd, &b_ThetaA_2nd);
  simDataTree->SetBranchAddress("PhiA_2nd", &PhiA_2nd, &b_PhiA_2nd);
  simDataTree->SetBranchAddress("ThetaB_2nd", &ThetaB_2nd, &b_ThetaB_2nd);
  simDataTree->SetBranchAddress("PhiB_2nd", &PhiB_2nd, &b_PhiB_2nd);
  
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

  //TTree* simDataTree =(TTree*)inputFile->Get("Tangle2");
  simDataTree =(TTree*)inputFile->Get("Tangle2");
  simDataTree->SetBranchAddress("edep0", &edep0, &b_edep0);
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
  simDataTree->SetBranchAddress("ThetaA_1st", &ThetaA_1st, &b_ThetaA_1st);
  simDataTree->SetBranchAddress("PhiA_1st", &PhiA_1st, &b_PhiA_1st);
  simDataTree->SetBranchAddress("ThetaB_1st", &ThetaB_1st, &b_ThetaB_1st);
  simDataTree->SetBranchAddress("PhiB_1st", &PhiB_1st, &b_PhiB_1st);
  simDataTree->SetBranchAddress("dPhi_1st", &dPhi_1st, &b_dPhi_1st);
  simDataTree->SetBranchAddress("ThetaA_2nd", &ThetaA_2nd, &b_ThetaA_2nd);
  simDataTree->SetBranchAddress("PhiA_2nd", &PhiA_2nd, &b_PhiA_2nd);
  simDataTree->SetBranchAddress("ThetaB_2nd", &ThetaB_2nd, &b_ThetaB_2nd);
  simDataTree->SetBranchAddress("PhiB_2nd", &PhiB_2nd, &b_PhiB_2nd);
  
  //write to

  sortDataTree = new TTree("sortDataTree",
			   "Sorted simulation events");

  //============================================
  //--------------------------------------------
  // Variable declarations

  //Energy drom the simulation will be smeared with a Gaussian
  //Set the HWHM for each crystal
  
  Double_t HWHM = 25.; //average from lab data

  // energy resolution coefficient k 
  // such that sigma(Energy) = k*Sqrt(Energy)
  // sigma = DeltaE_511 / sqrt(511) * 1 / (2sqrt(2ln2)
  // IE: k = DeltaE_511 / sqrt(511) = 25/sqrt(511)
  
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
    
    nb_ComptA[i] = 0;
    nb_ComptB[i] = 0;
    
  }
  
  for(Int_t i = 0 ; i < 2 ; i++){
    // angle after first Compton scattering 
    // retrieved from the simulation = 'sim theta'
    simtHA[i] = 180;
    simtHB[i] = 180;
    
    simPhiA[i] = 999;
    simPhiB[i] = 999;
  }
  
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

  
  tempString.Form("nb_ComptA[%d]/I",nCrystals);
  sortDataTree->Branch("nb_ComptA", nb_ComptA, tempString);

  tempString.Form("nb_ComptB[%d]/I",nCrystals);
  sortDataTree->Branch("nb_ComptB", nb_ComptB, tempString);

  tempString.Form("simtHA[%d]/D",2);
  sortDataTree->Branch("simtHA", simtHA, tempString);

  tempString.Form("simtHB[%d]/D",2);
  sortDataTree->Branch("simtHB", simtHB, tempString);

  tempString.Form("simPhiA[%d]/D",2);
  sortDataTree->Branch("simPhiA", simPhiA, tempString);

  tempString.Form("simPhiB[%d]/D",2);
  sortDataTree->Branch("simPhiB", simPhiB, tempString);


  tempString.Form("XposA[%d]/D",2);
  sortDataTree->Branch("XposA", XposA, tempString);

  tempString.Form("YposA[%d]/D",2);
  sortDataTree->Branch("YposA", YposA, tempString);

  tempString.Form("ZposA[%d]/D",2);
  sortDataTree->Branch("ZposA", ZposA, tempString);
  
  tempString.Form("XposB[%d]/D",2);
  sortDataTree->Branch("XposB", XposB, tempString);

  tempString.Form("YposB[%d]/D",2);
  sortDataTree->Branch("YposB", YposB, tempString);

  tempString.Form("ZposB[%d]/D",2);
  sortDataTree->Branch("ZposB", ZposB, tempString);

  // // per event variables
  // sortDataTree->Branch("dPhi_1st",&dPhi_1st, "dPhi_1st/D");
  
  TRandom *rand = new TRandom();

  //============
  // EVENT LOOP
  Float_t sigmaPar = 2.5;

  Long64_t nEvents = simDataTree->GetEntries();
  
  for(Int_t i = 0 ; i < nEvents ; i++){ 
    simDataTree->GetEvent(i);
    
    Double_t CrystEnergyDep[18] = {edep0,edep1,edep2,
				   edep3,edep4,edep5,
				   edep6,edep7,edep8,
				   edep9,edep10,edep11,
				   edep12,edep13,edep14,
				   edep15,edep16,edep17};
    for (Int_t k=0; k<18; k++)
      CrystEnergyDep[k] = 1000*CrystEnergyDep[k];

    Int_t nb_Compt[18] = {nb_Compt0,nb_Compt1,nb_Compt2,
			  nb_Compt3,nb_Compt4,nb_Compt5,
			  nb_Compt6,nb_Compt7,nb_Compt8,
			  nb_Compt9,nb_Compt10,nb_Compt11,
			  nb_Compt12,nb_Compt13,nb_Compt14,
			  nb_Compt15,nb_Compt16,nb_Compt17};

    for (Int_t j = 0 ; j < nCrystals ; j++){
      
      // sigmaE = sigma511 * sqrt( E / 511) 
      
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
      
      nb_ComptA[j] = nb_Compt[j];
      nb_ComptB[j] = nb_Compt[j + nCrystals];
      
    }
    
    simtHA[0] = ThetaA_1st;
    simtHB[0] = ThetaB_1st;
    simtHA[1] = ThetaA_2nd;
    simtHB[1] = ThetaB_2nd;
    
    simPhiA[0] = PhiA_1st;
    simPhiB[0] = PhiB_1st;
    simPhiA[1] = PhiA_2nd;
    simPhiB[1] = PhiB_2nd;
    
    XposA[0] = XposA_1st;
    YposA[0] = YposA_1st;
    ZposA[0] = ZposA_1st;
    XposB[0] = XposB_1st;
    YposB[0] = YposB_1st;
    ZposB[0] = ZposB_1st;

    XposA[1] = XposA_2nd;
    YposA[1] = YposA_2nd;
    ZposA[1] = ZposA_2nd;
    XposB[1] = XposB_2nd;
    YposB[1] = YposB_2nd;
    ZposB[1] = ZposB_2nd;
    
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


Bool_t TSim::CentralXA(Double_t posXA){
  
  Bool_t centralXA = kFALSE;
  
  if( posXA > 0. && 
      posXA < 54.25 )
    centralXA = kTRUE;

//   if( posXA > 35. && 
//       posXA < 45. )
//     centralXA = kTRUE;

  return centralXA;
}

Bool_t TSim::CentralXB(Double_t posXB){
  
  Bool_t centralXB = kFALSE;
  
  if( posXB < 0.  && 
      posXB > -54.25 )
    centralXB = kTRUE;
  
//   if( posXB < -35.  && 
//       posXB > -45. )
//     centralXB = kTRUE;
  
  return centralXB;
}


Bool_t TSim::CentralYZ(Double_t posYZ){
  
  Bool_t centralYZ = kFALSE;
  
  //!!
  Float_t crystalHalfSizeYZ = 2.0;
  
  posYZ = Abs(posYZ);

  if(posYZ < crystalHalfSizeYZ)
    centralYZ = kTRUE;
  
  return centralYZ;
}

Bool_t TSim::OuterYZ(Double_t posYZ){
  
  Bool_t outerYZ = kFALSE;
  
  Float_t crystalHalfSizeYZ = 2.0;
  
  posYZ = Abs(posYZ);
  
  if(posYZ > crystalHalfSizeYZ )
    outerYZ = kTRUE;

  return outerYZ;
}

// original simulation had 3 x 4 x 22 mm crystals
Bool_t TSim::CentralZ(Double_t posZ){
  
  Bool_t centralZ = kFALSE;
  
  //! works only for this particular detector geometry
  posZ = Abs(posZ);

  if( posZ < 1.5)
    centralZ = kTRUE;

  return centralZ;
}

Float_t TSim::CrystalToPhi(Int_t crystal){

//   Float_t crystalToPhi[9] = { -1.  ,   0. , -1.,
// 			      270. ,  -1. , 90.,
// 			      -1.  , 180. , -1.};
  
  // !! Use corner crystals
  Float_t crystalToPhi[9] = {  0.  ,  -1  , 90. ,
			       -1. ,  -1. , -1. ,
			       270.,  -1. , 180. };
  
  Float_t phi = crystalToPhi[crystal];

  return phi;
}

Float_t TSim::GetAsymLab(Int_t dPhiDiff, Int_t i){

  Float_t AsPhiDiff = 0;
  if (AsymMatrix[i][0] != 0){
      if (dPhiDiff  == 90){
	AsPhiDiff = AsymMatrix[i][1]/AsymMatrix[i][0];
	//using average 90 and 270
	AsPhiDiff = (AsymMatrix[i][1]+AsymMatrix[i][3])/(2*AsymMatrix[i][0]);
      }
      if (dPhiDiff  == 180)
	AsPhiDiff = AsymMatrix[i][2]/AsymMatrix[i][0];
      if (dPhiDiff  == 270)
	AsPhiDiff = AsymMatrix[i][3]/AsymMatrix[i][0];
    }
  return AsPhiDiff;
}

Float_t TSim::GetAsymLabErr(Int_t dPhiDiff, Int_t i){

  Float_t AsPhiDiff = 0;
  Float_t AePhiDiff = 0;
  if (AsymMatrix[i][0] != 0){
      if (dPhiDiff  == 90){
	AsPhiDiff = AsymMatrix[i][1]/AsymMatrix[i][0];
	AePhiDiff = AsPhiDiff*Sqrt((1/AsymMatrix[i][1])+(1/AsymMatrix[i][0]));
	//using average 90 and 270
	AsPhiDiff = (AsymMatrix[i][1]+AsymMatrix[i][3])/(2*AsymMatrix[i][0]);
	AePhiDiff = AsPhiDiff*Sqrt((1/(AsymMatrix[i][1]+AsymMatrix[i][3]))+(1/AsymMatrix[i][0]));
      }
      if (dPhiDiff  == 180){
	AsPhiDiff = AsymMatrix[i][2]/AsymMatrix[i][0];
	AePhiDiff = AsPhiDiff*Sqrt((1/AsymMatrix[i][2])+(1/AsymMatrix[i][0]));
      }
      if (dPhiDiff  == 270){
	AsPhiDiff = AsymMatrix[i][3]/AsymMatrix[i][0];
	AePhiDiff = AsPhiDiff*Sqrt((1/AsymMatrix[i][3])+(1/AsymMatrix[i][0]));
      }
    }
  
  return AePhiDiff;
}

Float_t TSim::GetAsymLabTrue(Int_t dPhiDiff, Int_t i){

  Float_t AsTrue = 0;
  if (AsymTrue[i][0] != 0){
      if (dPhiDiff  == 90){
	AsTrue = AsymTrue[i][1]/AsymTrue[i][0];
	//using average 90 and 270
	AsTrue = (AsymTrue[i][1]+AsymTrue[i][3])/(2*AsymTrue[i][0]);
      }
      if (dPhiDiff  == 180)
	AsTrue = AsymTrue[i][2]/AsymTrue[i][0];
      if (dPhiDiff  == 270)
	AsTrue = AsymTrue[i][3]/AsymTrue[i][0];
    }
  return AsTrue;
}

Float_t TSim::GetAsymLabTrueErr(Int_t dPhiDiff, Int_t i){

  Float_t AsTrue = 0;
  Float_t AeTrue = 0;
  if (AsymTrue[i][0] != 0){
      if (dPhiDiff  == 90){
	AsTrue = AsymTrue[i][1]/AsymTrue[i][0];
	AeTrue = AsTrue*Sqrt((1/AsymTrue[i][1])+(1/AsymTrue[i][0]));
	//using average 90 and 270
	AsTrue = (AsymTrue[i][1]+AsymTrue[i][3])/(2*AsymTrue[i][0]);
	AeTrue = AsTrue*Sqrt((1/(AsymTrue[i][1]+AsymTrue[i][3]))+(1/AsymTrue[i][0]));
      }
      if (dPhiDiff  == 180){
	AsTrue = AsymTrue[i][2]/AsymTrue[i][0];
	AeTrue = AsTrue*Sqrt((1/AsymTrue[i][2])+(1/AsymTrue[i][0]));
      }
      if (dPhiDiff  == 270){
	AsTrue = AsymTrue[i][3]/AsymTrue[i][0];
	AeTrue = AsTrue*Sqrt((1/AsymTrue[i][3])+(1/AsymTrue[i][0]));
      }
    }
  
  return AeTrue;
}

void TSim::GetThetaBinValues(){
  Float_t thetaLowEdge  = 10.;
  Float_t thetaHighEdge = 170.;

  Float_t thetaBinWidth = (thetaHighEdge - thetaLowEdge)/(Float_t)nThbins;

  cout << endl;
  for (Int_t i = 0 ; i < nThbins ; i++){
    ThMin[i] = thetaLowEdge + i * thetaBinWidth;
    ThMax[i] = thetaLowEdge + (i+1)* thetaBinWidth;
    plotTheta[i] = thetaLowEdge + (i+0.5)* thetaBinWidth;
    
    //cout << " plotTheta[" << i << "] = " << plotTheta[i] << endl;
  }
}

Int_t TSim::CalculateAsymmetryLab(TString inputFileNumber){

  GetThetaBinValues();

  //difference 'sim theta' - 'energy theta'
  TH2F* histThD = new TH2F("histThD",
			   "sim - energy #theta difference",
			   160,10,170,250,-120,130);

  //difference 'lab theta' - 'energy theta'
  TH2F* hLEdiff = new TH2F("hLEdiff","lab - energy theta",
			   160,10,170,120,-50,70);

  //difference 'max/min lab theta' - 'energy theta'
  TH2F* hMEdiff = new TH2F("hMEdiff","max/min - energy theta",
			   160,10,170,120,-50,70);

  
  TH1F * hThRes00[nThbins];
  TH1F * hThRes90[nThbins];

  TH1F * hThRes00_TL[nThbins];
  TH1F * hThRes90_TL[nThbins];

  TH1F * hDPhiRes00[nThbins];
  TH1F * hDPhiRes90[nThbins];
  
  TH1F * hDPhiRes00_TL[nThbins];
  TH1F * hDPhiRes90_TL[nThbins];
  
  TH1F * hBeta[nThbins];
  
  TString hTitle = "hThRes00";
  
  for( Int_t th = 0 ; th < nThbins ; th++){
    
    hTitle.Form("hThRes00_%d",th);
    hThRes00[th] = new TH1F(hTitle,hTitle,
			    128, -10., 180.);
    
    hTitle.Form("hThRes90_%d",th);
    hThRes90[th] = new TH1F(hTitle,hTitle,
			    128, -10., 180.);
    
    hTitle.Form("hThRes00_TL_%d",th);
    hThRes00_TL[th] = new TH1F(hTitle,hTitle,
			       128, -10., 180.);
    
    hTitle.Form("hThRes90_TL_%d",th);
    hThRes90_TL[th] = new TH1F(hTitle,hTitle,
			       128, -10., 180.);
    
    hTitle.Form("hDPhiRes00_%d",th);
    hDPhiRes00[th] = new TH1F(hTitle,hTitle,
			      32, -180.,180.);
    
    hTitle.Form("hDPhiRes90_%d",th);
    hDPhiRes90[th] = new TH1F(hTitle,hTitle,
			      32, -180.,180.);
    
    hTitle.Form("hDPhiRes00_TL_%d",th);
    hDPhiRes00_TL[th] = new TH1F(hTitle,hTitle,
				 32, -180.,180.);
    
    hTitle.Form("hDPhiRes90_TL_%d",th);
    hDPhiRes90_TL[th] = new TH1F(hTitle,hTitle,
				 32, -180.,180.);
  
    hBeta[th] = new TH1F(hTitle,hTitle,
			 32,-0.0, 10.0);
  }
  
  cout << endl;
  cout << " Getting asymmetry " << endl;
  
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
  
  sortDataTree=(TTree*)inputFile->Get("sortDataTree");

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
  sortDataTree->SetBranchAddress("simPhiA",&simPhiA);
  sortDataTree->SetBranchAddress("simPhiB",&simPhiB);

  sortDataTree->SetBranchAddress("nb_ComptA",&nb_ComptA);
  sortDataTree->SetBranchAddress("nb_ComptB",&nb_ComptB);
  
  sortDataTree->SetBranchAddress("XposA",&XposA);
  sortDataTree->SetBranchAddress("YposA",&YposA);
  sortDataTree->SetBranchAddress("ZposA",&ZposA);

  sortDataTree->SetBranchAddress("XposB",&XposB);
  sortDataTree->SetBranchAddress("YposB",&YposB);
  sortDataTree->SetBranchAddress("ZposB",&ZposB);
  
  //sortDataTree->SetBranchAddress("dPhi_1st",&dPhi_1st);
  
  n000 = 0;
  n090 = 0;
  n180 = 0;
  n270 = 0;

  for(Int_t j = 0 ; j <nThbins; j++){
    for(Int_t k = 0 ; k < 4; k++){
      AsymMatrix[j][k] = 0;}
  }

  for(Int_t j = 0 ; j < nThbins; j++){
    for(Int_t k = 0 ; k < 4; k++){
      AsymTrue[j][k] = 0;}
  }
  
  Long64_t nEntries = sortDataTree->GetEntries();

  for(Int_t j = 0 ; j <nThbins; j++){
    P_lab[j]  = 0.;
    Pq_lab[j] = 0.;
    P_true_lab[j]  = 0.;
    Pq_true_lab[j] = 0.;
  }

  Int_t   thBin  = -1;
  Float_t phiA   = -1;
  Float_t phiB   = -1;
  Int_t   indexA = -1;
  Int_t   indexB = -1;
  Float_t totEA  = 0;
  Float_t totEB  = 0;
  
  // Float_t phiB_prev = -1;
  // Float_t phiB_this = -1;

  Float_t phiDiff = -10.;
  Float_t dPhiXact = -999.;
  
  //Event loop
  for (Int_t ientry = 0 ; ientry < nEntries ; ientry++){
    sortDataTree->GetEvent(ientry);
      
   thBin  = -1;
   phiA   = -1;
   phiB   = -1;
   indexA = -1;
   indexB = -1;
   totEA  = 0;
   totEB  = 0;
   
   for (Int_t k = 0; k < nCrystals; k++){
      totEA += EA[k];
      totEB += EB[k];
   }
   
   // Select same theta bins and total energy
   if (GetThetaBin(ltHA[4]) == GetThetaBin(ltHB[4]) &&
       (totEA > 450) && (totEA < 550)               &&
       (totEB > 450) && (totEB < 550)){
     
     thBin = GetThetaBin(ltHA[4]);
     
     for (Int_t i = 0 ; i < nCrystals ; i++){
       if ( (i!=4)                  && 
	    (thBin>-1)              && 
	    (ltHA[i]<ThMax[thBin])  &&
	    (ltHA[i]>ThMin[thBin]) ){
	 phiA = CrystalToPhi(i);
	 indexA = i;
       }
       if ( (i!=4)                  &&
	    (thBin>-1)              &&
	    (ltHB[i]<ThMax[thBin])  &&
	    (ltHB[i]>ThMin[thBin]) ){
	 phiB = CrystalToPhi(i);
	 indexB = i;
       }	
     }
   } // end of : if (GetThetaBin(ltHA
   
   // randomise phiB
   // TRandom1 * rand1 = new TRandom1(); 
    // phiB = rand1->Uniform()*360;
    // if     ( phiB < 90 ){
    //   phiB = 0.;
    // }
    // else if( phiB >= 90  &&
    // 	     phiB < 180){
    //   phiB = 90.;
    // }
    // else if( phiB >= 180 &&
    // 	     phiB < 270){
    //   phiB = 180.;
    // }
    // else if( phiB >= 270 &&
    // 	     phiB < 360.){
    //   phiB = 270.;
    // }
    
    // make phiB value from previous event
    // phiB_this = phiB;
    // phiB = phiB_prev;
    // phiB_prev = phiB_this;
        
    // cout << endl;
    // cout << " phiB = " << phiB << endl;

    if ( (phiA > -1) && 
	 (phiB > -1) ){
      
      histThD->Fill(etHA[4],simtHA[0]-etHA[4]);
      histThD->Fill(etHB[4],simtHB[0]-etHB[4]);
      
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
      
      phiDiff   = phiB - phiA;
      
      dPhiXact  = simPhiA[0] + simPhiB[0];
      if(dPhiXact < 0)
	dPhiXact = dPhiXact + 360.;
      
      if (phiDiff == 0.) {
       
	if (dPhiXact > 180.)
	  dPhiXact = dPhiXact - 360. ;
	
	AsymMatrix[thBin][0] += 1;
	n000++;
	
	hThRes00[thBin]->Fill(simtHA[0]);
	hDPhiRes00[thBin]->Fill(dPhiXact);	
      }
      if ((phiDiff == 90)||(phiDiff == -270)){
	
	AsymMatrix[thBin][1] += 1;
	n090++;
	
	hThRes90[thBin]->Fill(simtHA[0]);
	
	dPhiXact = dPhiXact-90.;
	
	if(dPhiXact < 0)
	  dPhiXact = dPhiXact + 360.;
	
	if (dPhiXact > 180.)
	  dPhiXact =  dPhiXact - 360. ;
	
	hDPhiRes90[thBin]->Fill(dPhiXact);
	
      }
      if ((phiDiff == 180)||(phiDiff == -180)){
	AsymMatrix[thBin][2] += 1;
	n180++;}
      if ((phiDiff == -90)||(phiDiff == 270)){
	AsymMatrix[thBin][3] += 1;
	n270++;}
      
//       cout << endl;
//       cout << endl;
//       cout << " YposA[1] = " << YposA[1] << endl;
//       cout << endl;
//       cout << endl;
      
      // 'true lab' analysis
      if (
	  // simtHA[0]<ThMax[thBin] && 
// 	  simtHA[0]>ThMin[thBin] &&
// 	  simtHB[0]<ThMax[thBin] && 
// 	  simtHB[0]>ThMin[thBin] &&
	  CentralYZ(YposA[0])    && 
	  CentralYZ(ZposA[0])    &&
	  CentralYZ(YposB[0])    && 
	  CentralYZ(ZposB[0])   
	  // nb_ComptA[4] == 1.     && 
// 	  nb_ComptB[4] == 1.    
	  
	  ){
	
	if (phiDiff == 0.){
	  AsymTrue[thBin][0] += 1;
	  
	  hThRes00_TL[thBin]->Fill(simtHA[0]);	
	  hDPhiRes00_TL[thBin]->Fill(dPhiXact);	
	  
	}
	if ((phiDiff == 90) || (phiDiff == -270)){
	  AsymTrue[thBin][1] += 1;
	  
	  hThRes90_TL[thBin]->Fill(simtHA[0]);
		  
	  if      (phiDiff == 90)
	    hDPhiRes90_TL[thBin]->Fill(dPhiXact);
 	  else if (phiDiff == -270)
 	    hDPhiRes90_TL[thBin]->Fill(dPhiXact);
	  
	}
	if ((phiDiff == 180) || (phiDiff == -180))
	  AsymTrue[thBin][2] += 1;
	if ((phiDiff == -90) || (phiDiff == 270))
	  AsymTrue[thBin][3] += 1;
		
      }
    }
  }
  
  cout<<"Lab Asymmetry"<<endl;
  cout<< "theta -" << " dPhi=0 -" << " dPhi=90 -" << " dPhi=180 -" << " dPhi=270" << endl;
  cout<<plotTheta[0]<<" - "<<AsymMatrix[0][0]<<" - "<<AsymMatrix[0][1]<<" - "<<AsymMatrix[0][2]<<" - "<<AsymMatrix[0][3]<<endl;
  cout<<plotTheta[1]<<" - "<<AsymMatrix[1][0]<<" - "<<AsymMatrix[1][1]<<" - "<<AsymMatrix[1][2]<<" - "<<AsymMatrix[1][3]<<endl;
  cout<<plotTheta[2]<<" - "<<AsymMatrix[2][0]<<" - "<<AsymMatrix[2][1]<<" - "<<AsymMatrix[2][2]<<" - "<<AsymMatrix[2][3]<<endl;
  cout<<plotTheta[3]<<" - "<<AsymMatrix[3][0]<<" - "<<AsymMatrix[3][1]<<" - "<<AsymMatrix[3][2]<<" - "<<AsymMatrix[3][3]<<endl;
  cout<<plotTheta[4]<<" - "<<AsymMatrix[4][0]<<" - "<<AsymMatrix[4][1]<<" - "<<AsymMatrix[4][2]<<" - "<<AsymMatrix[4][3]<<endl;
  cout<<plotTheta[5]<<" - "<<AsymMatrix[5][0]<<" - "<<AsymMatrix[5][1]<<" - "<<AsymMatrix[5][2]<<" - "<<AsymMatrix[5][3]<<endl;
  cout<<plotTheta[6]<<" - "<<AsymMatrix[6][0]<<" - "<<AsymMatrix[6][1]<<" - "<<AsymMatrix[6][2]<<" - "<<AsymMatrix[6][3]<<endl;
  cout<<plotTheta[7]<<" - "<<AsymMatrix[7][0]<<" - "<<AsymMatrix[7][1]<<" - "<<AsymMatrix[7][2]<<" - "<<AsymMatrix[7][3]<<endl;

  cout <<"True Asymmetry"<<endl;
  cout<< "theta -" << " dPhi=0 -" << " dPhi=90 -" << " dPhi=180 -" << " dPhi=270" << endl;
  cout<<plotTheta[0]<<" - "<<AsymTrue[0][0]<<" - "<<AsymTrue[0][1]<<" - "<<AsymTrue[0][2]<<" - "<<AsymTrue[0][3]<<endl;
  cout<<plotTheta[1]<<" - "<<AsymTrue[1][0]<<" - "<<AsymTrue[1][1]<<" - "<<AsymTrue[1][2]<<" - "<<AsymTrue[1][3]<<endl;
  cout<<plotTheta[2]<<" - "<<AsymTrue[2][0]<<" - "<<AsymTrue[2][1]<<" - "<<AsymTrue[2][2]<<" - "<<AsymTrue[2][3]<<endl;
  cout<<plotTheta[3]<<" - "<<AsymTrue[3][0]<<" - "<<AsymTrue[3][1]<<" - "<<AsymTrue[3][2]<<" - "<<AsymTrue[3][3]<<endl;
  cout<<plotTheta[4]<<" - "<<AsymTrue[4][0]<<" - "<<AsymTrue[4][1]<<" - "<<AsymTrue[4][2]<<" - "<<AsymTrue[4][3]<<endl;
  cout<<plotTheta[5]<<" - "<<AsymTrue[5][0]<<" - "<<AsymTrue[5][1]<<" - "<<AsymTrue[5][2]<<" - "<<AsymTrue[5][3]<<endl;
  cout<<plotTheta[6]<<" - "<<AsymTrue[6][0]<<" - "<<AsymTrue[6][1]<<" - "<<AsymTrue[6][2]<<" - "<<AsymTrue[6][3]<<endl;
  cout<<plotTheta[7]<<" - "<<AsymTrue[7][0]<<" - "<<AsymTrue[7][1]<<" - "<<AsymTrue[7][2]<<" - "<<AsymTrue[7][3]<<endl;
    
  TCanvas *canvas = new TCanvas("canvas","canvas",
				1500,1000);

  Char_t plotN[128];

  hLEdiff->Draw("colz");
  hLEdiff->GetXaxis()->SetTitle("energy #theta (deg)");
  hLEdiff->GetYaxis()->SetTitle("lab #theta - energy #theta (deg)");

  sprintf(plotN,"../Plots/histLEdiff_%d.pdf",inputFileInt);
  
  plotName = "../Plots/histLEdiff_" + inputFileNumber;
  plotName += ".pdf";

  canvas->SaveAs(plotName);
   
  histThD->Draw("colz");
  histThD->GetXaxis()->SetTitle("#theta (deg)");
  histThD->GetYaxis()->SetTitle("sim #theta - energy #theta (deg)");
  
  sprintf(plotN,"../Plots/histThDiff_%d.pdf",inputFileInt);

  plotName = "../Plots/histThDiff_" + inputFileNumber;
  plotName += ".pdf";

  canvas->SaveAs(plotName);

  hMEdiff->Draw("colz");
  hMEdiff->GetXaxis()->SetTitle("energy #theta (deg)");
  hMEdiff->GetYaxis()->SetTitle("max/min lab #theta - energy #theta (deg)");
  
  sprintf(plotN,"../Plots/histME_%d.pdf",inputFileInt);
  
  plotName = "../Plots/histME_" + inputFileNumber;
  plotName += ".pdf";

  canvas->SaveAs(plotName);
  
  canvas->Clear();
  
  canvas->Divide(4,2);
  
  for( Int_t th = 0 ; th < nThbins ; th++){

    hThRes00[th]->SetLineColor(kBlue);
    hDPhiRes00[th]->SetLineColor(kBlue);

    hThRes90[th]->SetLineColor(kRed);
    hDPhiRes90[th]->SetLineColor(kRed);
    
    hThRes00_TL[th]->SetLineColor(kBlue);
    hDPhiRes00_TL[th]->SetLineColor(kBlue);
    
    hThRes90_TL[th]->SetLineColor(kRed);
    hDPhiRes90_TL[th]->SetLineColor(kRed);
    
    hThRes00[th]->GetXaxis()->SetTitle("#theta_{A} (deg)");
    hDPhiRes00[th]->GetXaxis()->SetTitle("#Delta#phi (deg)");

    hThRes90[th]->GetXaxis()->SetTitle("#theta_{A} (deg)");
    hDPhiRes90[th]->GetXaxis()->SetTitle("#Delta#phi (deg)");

    hThRes00_TL[th]->GetXaxis()->SetTitle("#theta_{A} (deg)");
    hDPhiRes00_TL[th]->GetXaxis()->SetTitle("#Delta#phi (deg)");
    
    hThRes90_TL[th]->GetXaxis()->SetTitle("#theta_{A} (deg)");
    hDPhiRes90_TL[th]->GetXaxis()->SetTitle("#Delta#phi (deg)");
    
    hThRes00[th]->GetYaxis()->SetTitle("Counts");
    hDPhiRes00[th]->GetYaxis()->SetTitle("Counts");
    
    hThRes90[th]->GetYaxis()->SetTitle("Counts");
    hDPhiRes90[th]->GetYaxis()->SetTitle("Counts");
      
    hThRes00_TL[th]->GetYaxis()->SetTitle("Counts");
    hDPhiRes00_TL[th]->GetYaxis()->SetTitle("Counts");
    
    hThRes90_TL[th]->GetYaxis()->SetTitle("Counts");
    hDPhiRes90_TL[th]->GetYaxis()->SetTitle("Counts");
    
  }    
  
  for( Int_t th = 0 ; th < nThbins ; th++){
    canvas->cd(th+1);
    hThRes90[th]->Draw("");
    hThRes00[th]->Draw("same");  
  }
  plotName = "../Plots/hThRes_" + inputFileNumber;
  plotName += ".pdf";
  canvas->SaveAs(plotName);
  
  for( Int_t th = 0 ; th < nThbins ; th++){
    canvas->cd(th+1);
    hThRes90_TL[th]->Draw("");
    hThRes00_TL[th]->Draw("same");  
  }
  plotName = "../Plots/hThRes_TL_" + inputFileNumber;
  plotName += ".pdf";
  canvas->SaveAs(plotName);

  for( Int_t th = 0 ; th < nThbins ; th++){
    canvas->cd(th+1);
    hDPhiRes90[th]->Draw("");
    hDPhiRes00[th]->Draw("same");  
  }
  plotName = "../Plots/hDPhiRes_" + inputFileNumber;
  plotName += ".pdf";
  canvas->SaveAs(plotName);

  for( Int_t th = 0 ; th < nThbins ; th++){
    canvas->cd(th+1);
    hDPhiRes90_TL[th]->Draw("");
    hDPhiRes00_TL[th]->Draw("same");  
  }
  plotName = "../Plots/hDPhiRes_TL_" + inputFileNumber;
  plotName += ".pdf";
  canvas->SaveAs(plotName);
  
  return 0;
}

void TSim::GraphAsymmetryLab(TString inputFileNumber1,
			     TString inputFileNumber2) {
  
  cout << endl;
  cout << "--------------------- " << endl;
  cout << "  GraphAsymmetryLab   " << endl;
  
  this->SetStyle();

  // GetThetaBinValues() is implemented in
  // CalculateAsymmetryLab
    
  CalculateAsymmetryLab(inputFileNumber1);
  
  Int_t inputFileInt1 = inputFileNumber1.Atoi();
  
  Int_t inputFileInt2 = 0;
  
  Int_t nFiles = 2;

  if(inputFileNumber2=="??")
    nFiles = 1;
  else 
    inputFileInt2 = inputFileNumber2.Atoi();
  
  cout << "  " << nFiles << " Files " << endl;

  TString plotName;
  plotName = "../Plots/Asym_" + inputFileNumber1;
  plotName = plotName + ".pdf";

  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = P(90)/P(0) 
  Int_t   dPhiDiff = 90;
  
  Float_t thetaLowEdge  = 10.;
  Float_t thetaHighEdge = 170.;
  
  CalculateABC_Lab();
  CalculateABC_True();
  
  // To Do:
  //CalculateABC_Lab_Err()
  

  for (Int_t i = 0 ; i < nThbins ; i++){
    pA_Lab[i] = pABC_Lab[i][0];
    pB_Lab[i] = pABC_Lab[i][1];
    pC_Lab[i] = pABC_Lab[i][2];
    
    pA_True[i] = pABC_True[i][0];
    pB_True[i] = pABC_True[i][1];
    pC_True[i] = pABC_True[i][2];
  }
 
  TGraphErrors *grB_Lab;
  TGraphErrors *grC_Lab;
  TGraphErrors *grB_True;
  TGraphErrors *grC_True;
  
  grB_Lab     = new TGraphErrors(nThbins,plotTheta,pB_Lab,0,0);
  grC_Lab     = new TGraphErrors(nThbins,plotTheta,pC_Lab,0,0);
  grB_True    = new TGraphErrors(nThbins,plotTheta,pB_True,0,0);
  grC_True    = new TGraphErrors(nThbins,plotTheta,pC_True,0,0);
  
  //Calculating ratios for desired dPhiDiff
  Float_t AsPhiDiff[nThbins] = {0.};
  Float_t AePhiDiff[nThbins] = {0.};

  Float_t AsTrue[nThbins] = {0.};
  Float_t AeTrue[nThbins] = {0.};

  Float_t AsPhiDiff1[nThbins] = {0.};
  Float_t AePhiDiff1[nThbins] = {0.};
  Float_t AsPhiDiffR[nThbins] = {0.};
  Float_t AePhiDiffR[nThbins] = {0.};

  Float_t AsTrue1[nThbins] = {0.};
  Float_t AeTrue1[nThbins] = {0.};
  Float_t AsTrueR[nThbins] = {0.};
  Float_t AeTrueR[nThbins] = {0.};
  
  Float_t mu[nThbins]={0};
  Float_t muE[nThbins]={0};

  for (Int_t file = 0 ; file < nFiles ; file++){
    
    if (file==1)
      CalculateAsymmetryLab(inputFileNumber2);

    for (Int_t i = 0 ; i < nThbins ; i++){
      if (AsymMatrix[i][0] != 0){
	
	mu[i] = (AsymMatrix[i][1] - AsymMatrix[i][0]);
	mu[i] = mu[i]/(AsymMatrix[i][1] + AsymMatrix[i][0]);
	
	if (dPhiDiff  == 90){
	  AsPhiDiff[i] = AsymMatrix[i][1]/AsymMatrix[i][0];
	  AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix[i][1])+(1/AsymMatrix[i][0]));
	  //using average 90 and 270
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
    
      if     (file==0){
	AsPhiDiff1[i] = AsPhiDiff[i];
	AePhiDiff1[i] = AePhiDiff[i];
	
      }
      else if(file==1){
	AsPhiDiffR[i] = AsPhiDiff1[i]/AsPhiDiff[i] ;
	
	AePhiDiffR[i] = AsPhiDiffR[i] * Sqrt( AePhiDiff[i]*AePhiDiff[i]/(AsPhiDiff[i]*AsPhiDiff[i]) + AePhiDiff1[i]*AePhiDiff1[i]/(AsPhiDiff1[i]*AsPhiDiff1[i]));

      }
         
    
    }
    
    for (Int_t i = 0 ; i < nThbins ; i++){
      if (AsymTrue[i][0] != 0){
	if (dPhiDiff  == 90){
	  AsTrue[i] = AsymTrue[i][1]/AsymTrue[i][0];
	  AeTrue[i] = AsTrue[i]*Sqrt((1/AsymTrue[i][1])+(1/AsymTrue[i][0]));
	  //using average 90 and 270
	  AsTrue[i] = (AsymTrue[i][1]+AsymTrue[i][3])/(2*AsymTrue[i][0]);
	  AeTrue[i] = AsTrue[i]*Sqrt((1/(AsymTrue[i][1]+AsymTrue[i][3]))+(1/AsymTrue[i][0]));
	}
	if (dPhiDiff  == 180){
	  AsTrue[i] = AsymTrue[i][2]/AsymTrue[i][0];
	  AeTrue[i] = AsTrue[i]*Sqrt((1/AsymTrue[i][2])+(1/AsymTrue[i][0]));
	}
	if (dPhiDiff  == 270){
	  AsTrue[i] = AsymTrue[i][3]/AsymTrue[i][0];
	  AeTrue[i] = AsTrue[i]*Sqrt((1/AsymTrue[i][3])+(1/AsymTrue[i][0]));
	}
      }
    
      if     (file==0){
	AsTrue1[i] = AsTrue[i];
	AeTrue1[i] = AeTrue[i];
	
      }
      else if(file==1){
	AsTrueR[i] = AsTrue1[i]/AsTrue[i] ;
	
	AeTrueR[i] = AsTrueR[i] * Sqrt( AeTrue[i]*AeTrue[i]/(AsTrue[i]*AsTrue[i]) + AeTrue1[i]*AeTrue1[i]/(AsTrue1[i]*AsTrue1[i]));
      }
      
    } // end of: for (Int_t i = 0 ; i < nThbins ...
  
  } //end of: for (Int_t file = 0 ; 
  
  //theory curve
  Float_t aTheory[nThbins];

  //half resolution in dPhi

  Float_t alpha1 = DegToRad()*26.0*2.34/2.;


  //half resolution in theta
  Float_t semiSpan = DegToRad()*(ThMax[0] - ThMin[0])/2.;

  TTheory *theory = new TTheory();

  //calculating theory curve
  cout << endl;
  cout << " Calculating theory curve ... " << endl;
  cout << endl;
  cout << " semiSpan = " << semiSpan*RadToDeg() << endl;
  cout << " alpha1   = " << alpha1*RadToDeg()   << endl;
  
  for (Int_t i = 0; i<nThbins; i++){
    if( dPhiDiff == 180 ){
	aTheory[i] = 1.0;
	continue;
      }
    plotTheta[i] = plotTheta[i]*DegToRad();
    aTheory[i] = theory->rho2(plotTheta[i],semiSpan,alpha1);
    plotTheta[i] = plotTheta[i]*RadToDeg();
  }
    
 
  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);

  TH1F *hr;
  
  Float_t maxY = 4.0;

  if(dPhiDiff==180){
    maxY = 6.0;
    if(nFiles==2){
	maxY = 2.0;
      }
  }
  
   
  hr = canvas->DrawFrame(10,0.5,170,maxY);
  hr->GetXaxis()->SetTitle("#theta (deg)");

  TGraphErrors *grAsym[3];
  
  TGraphErrors * grMu = new TGraphErrors(nThbins,plotTheta,mu,0,muE);

  if     (nFiles==1)
    grAsym[0] = new TGraphErrors(nThbins,plotTheta,AsPhiDiff,0,AePhiDiff);
  else if(nFiles == 2)
    grAsym[0] = new TGraphErrors(nThbins,plotTheta,AsPhiDiffR,0,AePhiDiffR);
  
  grAsym[1] = new TGraphErrors(nThbins,plotTheta,aTheory,0,0);

  for(Int_t i = 0; i < nThbins; i++)
    plotTheta[i] += 2.;

  if     (nFiles==1)
    grAsym[2] = new TGraphErrors(nThbins,plotTheta,AsTrue,0,AeTrue);
  else if(nFiles == 2)
    grAsym[2] = new TGraphErrors(nThbins,plotTheta,AsTrueR,0,AeTrueR);

  grAsym[0]->SetLineColor(kGreen+2.7);
  grAsym[0]->SetMarkerColor(kGreen+2.7);
  grAsym[1]->SetLineColor(kRed);
  grAsym[1]->SetMarkerColor(kRed);
  grAsym[2]->SetLineColor(kGreen);
  grAsym[2]->SetMarkerColor(kGreen);

  TLegend *leg = new TLegend(0.6,0.75,0.9,0.85);
  
  
  
  TString theoryLegendTitle = " ";
  alpha1 = alpha1*RadToDeg();
  theoryLegendTitle.Form("theory curve #alpha_{#phi} = %.1f^{o}", alpha1);
  
  Char_t  plotN[128];
  //TString plotName = "";
  Char_t  yAxis[128];
  
  if(nFiles==1)
    sprintf(yAxis,"N(%d^{o})/N(0^{o})",dPhiDiff);
  else if(nFiles == 2)
    sprintf(yAxis,"N^{E}(%d^{o})/N^{E}(0^{o}) * N^{U}(0^{o})/N^{U}(%d^{o}) ",dPhiDiff,dPhiDiff);
  
  hr->GetYaxis()->SetTitle(yAxis);
  leg->AddEntry(grAsym[1],
		  theoryLegendTitle,"L P");
  leg->AddEntry(grAsym[0],"simulated lab", "E P");
  leg->AddEntry(grAsym[2],"true simulated lab","E P");
  grAsym[0]->Draw("P E");
  grAsym[1]->Draw("P L SAME");
  grAsym[2]->Draw("P E");
  leg->Draw();

  if     (nFiles==1){
    sprintf(plotN,"../Plots/A_%d_%d.pdf",inputFileInt1,dPhiDiff);
    plotName = "../Plots/A_" + inputFileNumber1;
  }
  else if(nFiles==2){
    sprintf(plotN,"../Plots/A_%d_%d_%d.pdf",inputFileInt1,inputFileInt2,dPhiDiff);
    plotName  =  "../Plots/A_" + inputFileNumber1;
    plotName += "_";
    plotName += inputFileNumber2;
  }  
  
  plotName.Form(plotName + "_%d",dPhiDiff);
  plotName += ".pdf";
  
  canvas->SaveAs(plotName);
  
  // ------------------------------------------------
  //  Fourier Coefficients Graph

  leg->Clear();
  
  maxY =  1.0;
  Float_t minY = -1.2;
    
  hr = canvas->DrawFrame(10,minY,170,maxY);
  hr->GetXaxis()->SetTitle("#theta (deg)");
  hr->GetYaxis()->SetTitle("cos(#Delta#phi) coefficient");
  hr->GetYaxis()->SetTitleOffset(0.7);

  grB_Lab->Draw("L P");
  
  grB_True->SetLineColor(kGreen);
  grB_True->Draw("L P same");
  
  plotName = "../Plots/B_coeff_";
  plotName = plotName + inputFileNumber1;
  plotName = plotName + ".pdf";
    
  canvas->SaveAs(plotName);
  
  maxY =  0.6;
  minY = -0.4;
    
  hr = canvas->DrawFrame(thetaLowEdge,minY,thetaHighEdge,maxY);
  hr->GetXaxis()->SetTitle("#theta (deg)");
  hr->GetYaxis()->SetTitle("cos(2#Delta#phi) coefficient");
  hr->GetYaxis()->SetTitleOffset(0.7);
  
  grC_Lab->Draw("L P");
  
  grC_True->SetLineColor(kGreen);
  grC_True->Draw("L P same");
  
  plotName = "../Plots/C_coeff_";
  plotName = plotName + inputFileNumber1;
  plotName = plotName + ".pdf";
    
  canvas->SaveAs(plotName);
  
    // mu plot
  hr = canvas->DrawFrame(thetaLowEdge,-0.2,
			  thetaHighEdge,0.5);
  
  hr->GetXaxis()->SetTitle("#theta (deg)");

  sprintf(yAxis,"(N(%d^{o})-N(0^{o})/(N(%d^{o})+N(0^{o})",
	  dPhiDiff,dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);

  plotName.Form("../Plots/Mu_%d_%d.pdf",inputFileInt1,dPhiDiff);

  grMu->Draw("P E");

  canvas->SaveAs(plotName);
  
  cout << "  GraphAsymmetryLab   " << endl;
  cout << "--------------------- " << endl;
  cout << endl;
  
}
//-----------------------------------------------------------------------

void TSim::CalculateABC(){
  
  for(Int_t i = 0 ; i < nThbins ; i++){
    
    for(Int_t j = 0 ; j < 3 ; j++){
      fABC[i][j] = 0.;
      pABC[i][j] = 0.;
    }
    
    if (AsymMatrix_sim[i][bin000]== 0) continue;
    
    fABC[i][0]  =  AsymMatrix_sim[i][bin000];
    fABC[i][1]  =  AsymMatrix_sim[i][bin000];
    fABC[i][2]  =  AsymMatrix_sim[i][bin000];
    
    fABC[i][0] +=  AsymMatrix_sim[i][bin090];
    fABC[i][1] +=  AsymMatrix_sim[i][bin000];
    fABC[i][2] -=  AsymMatrix_sim[i][bin090];
    
    fABC[i][0] +=  AsymMatrix_sim[i][bin180];
    fABC[i][1] -=  AsymMatrix_sim[i][bin180];
    fABC[i][2] +=  AsymMatrix_sim[i][bin180];
    
    fABC[i][0] +=  AsymMatrix_sim[i][bin270];
    fABC[i][1] -=  AsymMatrix_sim[i][bin180];
    fABC[i][2] -=  AsymMatrix_sim[i][bin270];
    
    for(Int_t j = 0 ; j < 3 ; j++){
      fABC[i][j] = fABC[i][j]/4.;
      pABC[i][j] = fABC[i][j]/fABC[i][bin000];
      
//       cout << " pABC["<< i << "][" << j << "] = " 
// 	   <<  pABC[i][j] << endl;
    }
  }
  
}

void TSim::CalculateABC_True(){
  
  for(Int_t i = 0 ; i < nThbins ; i++){
    
    for(Int_t j = 0 ; j < 3 ; j++){
      fABC_True[i][j] = 0.;
      pABC_True[i][j] = 0.;
    }
    
    if (AsymTrue[i][bin000]== 0) continue;
    
    fABC_True[i][0]  =  AsymTrue[i][0];
    fABC_True[i][1]  =  AsymTrue[i][0];
    fABC_True[i][2]  =  AsymTrue[i][0];
    
    fABC_True[i][0] +=  AsymTrue[i][1];
    fABC_True[i][1] +=  AsymTrue[i][0];
    fABC_True[i][2] -=  AsymTrue[i][1];
    
    fABC_True[i][0] +=  AsymTrue[i][2];
    fABC_True[i][1] -=  AsymTrue[i][2];
    fABC_True[i][2] +=  AsymTrue[i][2];
    
    fABC_True[i][0] +=  AsymTrue[i][3];
    fABC_True[i][1] -=  AsymTrue[i][2];
    fABC_True[i][2] -=  AsymTrue[i][3];
    
    cout << endl;
    for(Int_t j = 0 ; j < 3 ; j++){
      fABC_True[i][j] = fABC_True[i][j]/4.;
      pABC_True[i][j] = fABC_True[i][j]/fABC_True[i][bin000];
            
//       cout << " pABC_True["<< i << "][" << j << "] = " 
// 	   <<  pABC_True[i][j] << endl;
    }
    
  } // end of: for(Int_t i = 0 ; i < nThbins ;
  
}

void TSim::CalculateABC_Lab(){
  
  for(Int_t i = 0 ; i < nThbins ; i++){
    
    for(Int_t j = 0 ; j < 3 ; j++){
      fABC_Lab[i][j] = 0.;
      pABC_Lab[i][j] = 0.;
    }
    
    if (AsymMatrix[i][0]== 0) continue;
    
    fABC_Lab[i][0]  =  AsymMatrix[i][0];
    fABC_Lab[i][1]  =  AsymMatrix[i][0];
    fABC_Lab[i][2]  =  AsymMatrix[i][0];
        
    fABC_Lab[i][0] +=  AsymMatrix[i][1];
    fABC_Lab[i][1] +=  AsymMatrix[i][0];
    fABC_Lab[i][2] -=  AsymMatrix[i][1];
    
    fABC_Lab[i][0] +=  AsymMatrix[i][2];
    fABC_Lab[i][1] -=  AsymMatrix[i][2];
    fABC_Lab[i][2] +=  AsymMatrix[i][2];
    
    fABC_Lab[i][0] +=  AsymMatrix[i][3];
    fABC_Lab[i][1] -=  AsymMatrix[i][2];
    fABC_Lab[i][2] -=  AsymMatrix[i][3];
    
//     cout << endl;
//     for(Int_t j = 0 ; j < 3 ; j++){
//       //fABC_Lab[i][j] = 1.0*fABC_Lab[i][j]/4.;
//       pABC_Lab[i][j] = 1.0*fABC_Lab[i][j]/fABC_Lab[i][0];
      
//       cout << " pABC_Lab["<< i << "][" << j << "] = " 
// 	   <<  pABC_Lab[i][j] << endl;
//     }
    
//     for(Int_t j = 0 ; j < 4 ; j++){
//       cout << " AsymMatrix["<< i << "][" << j << "] = " 
// 	   <<  AsymMatrix[i][j] << endl;
//     }
    
  }
}

Int_t TSim::CalculateAsymmetrySim(TString inputFileNumber){

  GetThetaBinValues();

  cout << endl;
  cout << " Getting asymmetry sim " << endl;
  
  TString plotName;
  plotName = "../Plots/Asym_sim_" + inputFileNumber;
  plotName = plotName + ".pdf";

  TString inputFileName = "../Data/sim" + inputFileNumber;
  inputFileName  = inputFileName + ".root";
  
  cout << endl;
  cout << " Input File : " << inputFileName  << endl;
  cout << endl;

  TFile* inputFile = new TFile(inputFileName);
  
  //TTree* simDataTree=(TTree*)inputFile->Get("Tangle2");
  simDataTree=(TTree*)inputFile->Get("Tangle2");
  simDataTree->SetBranchAddress("ThetaA_1st",&ThetaA_1st);
  simDataTree->SetBranchAddress("ThetaA_2nd",&ThetaA_2nd);
  simDataTree->SetBranchAddress("ThetaB_1st",&ThetaB_1st);
  simDataTree->SetBranchAddress("ThetaB_2nd",&ThetaB_2nd);
  simDataTree->SetBranchAddress("PhiA_1st",&PhiA_1st);
  simDataTree->SetBranchAddress("PhiA_2nd",&PhiA_2nd);
  simDataTree->SetBranchAddress("PhiB_1st",&PhiB_1st);
  simDataTree->SetBranchAddress("PhiB_2nd",&PhiB_2nd);
  
  for(Int_t j = 0 ; j <nThbins; j++)
    for(Int_t k = 0 ; k < nPhibinsSim; k++)
      AsymMatrix_sim[j][k] = 0;
  
  Float_t  halfBinSize = 180./nPhibinsSim; 

  Long64_t nEntries = simDataTree->GetEntries();
  for (Int_t ientry = 0 ; ientry < nEntries ; ientry++){
    simDataTree->GetEvent(ientry);
    
    Float_t dPhi_1st = PhiA_1st + PhiB_1st;
    if(dPhi_1st<0)
      dPhi_1st = dPhi_1st + 360; 
        
    Int_t thBin = -1;
    if (GetThetaBin(ThetaA_1st) == 
	GetThetaBin(ThetaB_1st)){
      thBin = GetThetaBin(ThetaA_1st);
    }
    
    if(thBin < 0) continue;      
    // dphi=0 bin
    if( (dPhi_1st <   halfBinSize     ) || 
	(dPhi_1st > (360-halfBinSize) ) )
      AsymMatrix_sim[thBin][bin000] += 1;
    // fill the rest 
    for (Int_t i = 1 ; i < nPhibinsSim ; i++){
      if( (dPhi_1st > (halfBinSize*(2*i - 1)) ) &&
	  (dPhi_1st < (halfBinSize*(2*i + 1)) ) )
	AsymMatrix_sim[thBin][i] += 1;
    }
  }//end of: for(Int_t ientry...
    
  // //printing out the asym matrix 
  
  Bool_t comments = kFALSE;
  if(comments)
    for(Int_t j = 0 ; j <nThbins; j++)
      for(Int_t k = 0 ; k < nPhibinsSim; k++)
	cout << " assym matrix for theta bin " 
	     << j 
	     << " and phi bin "               
	     << k
	     << " is "
	     << AsymMatrix_sim[j][k] 
	     << endl;
  
  return 0;

} //end of CalculateAsymmetrySim


//-----------------------------------------------------------------------
Int_t TSim::CalculateAsymmetrySimScattered(TString inputFileNumber,
					   Float_t thetaS,
					   Float_t thetaSHalf){
  
  GetThetaBinValues();
  
  cout << endl;
  cout << " Getting asymmetry sim scattered " << endl;
  
  cout << " thetaS     = " << thetaS     << endl;
  cout << " thetaSHalf = " << thetaSHalf << endl;

  TString plotName;
  plotName = "../Plots/Asym_sim_" + inputFileNumber;
  
  plotName = plotName + ".pdf";
  
  TString inputFileName = "../Data/sim" + inputFileNumber;
  inputFileName  = inputFileName + ".root";
    
  cout << endl;
  cout << " Input File : " << inputFileName  << endl;
  cout << endl;

  TFile* inputFile = new TFile(inputFileName);

  //Tree* simDataTree=(TTree*)inputFile->Get("Tangle2");
  simDataTree=(TTree*)inputFile->Get("Tangle2");
  
  simDataTree->SetBranchAddress("ThetaA_1st",&ThetaA_1st);
  simDataTree->SetBranchAddress("ThetaA_2nd",&ThetaA_2nd);
  simDataTree->SetBranchAddress("ThetaB_1st",&ThetaB_1st);
  simDataTree->SetBranchAddress("ThetaB_2nd",&ThetaB_2nd);

  simDataTree->SetBranchAddress("PhiA_1st",&PhiA_1st);
  simDataTree->SetBranchAddress("PhiA_2nd",&PhiA_2nd);
  simDataTree->SetBranchAddress("PhiB_1st",&PhiB_1st);
  simDataTree->SetBranchAddress("PhiB_2nd",&PhiB_2nd);
  
  simDataTree->SetBranchAddress("XposA_1st", &XposA_1st);
  simDataTree->SetBranchAddress("YposA_1st", &YposA_1st);
  simDataTree->SetBranchAddress("ZposA_1st", &ZposA_1st);

  simDataTree->SetBranchAddress("XposB_1st", &XposB_1st);
  simDataTree->SetBranchAddress("YposB_1st", &YposB_1st);
  simDataTree->SetBranchAddress("ZposB_1st", &ZposB_1st);
  
  
  for(Int_t j = 0 ; j <nThbins; j++){
    for(Int_t k = 0 ; k < nPhibinsSim; k++)
      AsymMatrix_sim[j][k] = 0;
  }
  
  Long64_t nEntries = simDataTree->GetEntries();
  
  
  Float_t halfBinSize = 180/nPhibinsSim; 
  
  // Event loop - scattering study
  
  Float_t dPhi;
  
  // thetas used in asymmetry calc
  Float_t thetaA  = 0.;
  Float_t thetaB  = 0.;
  // scattering angle in data
  Float_t thetaABS = 0.;

  Char_t scatterArray = 'N';
    
  for (Int_t ientry = 0 ; ientry < nEntries ; ientry++){
    
    scatterArray = 'N';
    
    thetaA  = 0.;
    thetaB  = 0.;
    thetaABS = 0.;
    
    dPhi = -999.;
    
    simDataTree->GetEvent(ientry);
    
    // scattering in array A
    if ( PhiA_2nd   < 499. && 
	 ThetaA_2nd < 499. ){
      
      scatterArray = 'A';
      
      thetaA   = ThetaA_2nd;
      thetaB   = ThetaB_1st;
      thetaABS = ThetaA_1st;
      //!!!!!
      // warning - check it is never negative
      dPhi =  PhiB_1st + PhiA_2nd;
    }
    
    // scattering in array B
    if( PhiB_2nd   < 499. &&
	ThetaB_2nd < 499.){
      
      if(scatterArray=='A'){
	// scattering in both
	scatterArray='C'; 
      }
      else{ 
	scatterArray = 'B';
	
	thetaA   = ThetaA_1st;
	thetaB   = ThetaB_2nd;
	thetaABS = ThetaB_1st;
	
	//!!!!!
	// warning - check it is never negative
	dPhi = PhiB_2nd + PhiA_1st;	
      }
    }
    
    if(scatterArray=='C') {
      //scattering in both
      thetaA   = ThetaA_2nd;
      thetaB   = ThetaB_2nd;
      continue; // for now
    }
    
    //cout << " scatterArray = " << scatterArray << endl;
    
    if( thetaA < 1.  ||
	thetaB < 1.
	)
      continue;
    
    if(dPhi < 0)
      dPhi = dPhi + 360; 
    
    Int_t thBin = -1;
    
    // theta1 == theta2 ? 
    if (GetThetaBin(thetaA) == GetThetaBin(thetaB)){
      thBin = GetThetaBin(thetaA);
    }
    
    if(thBin < 0) continue;      

    // scattering in both
    // if (GetThetaBin(ThetaB_1st) != GetThetaBin(ThetaA_1st) ||
    // 	GetThetaBin(ThetaB_1st) == -1)
    //       continue;
    
    // scattering angle in range and
    // hit in central crystal
    if( (thetaABS > (thetaS - thetaSHalf) ) &&
	(thetaABS < (thetaS + thetaSHalf) )){
      
	// fill dphi=0 bin
	if( (dPhi < halfBinSize       ) || 
	    (dPhi > (360-halfBinSize) )
	    )
	  AsymMatrix_sim[thBin][bin000] += 1;
	
	// fill the rest 
	// !! applying condition that scattering is
	// in central conditions !!
	for (Int_t i = 1 ; i < nPhibinsSim ; i++){
	  if( dPhi > (halfBinSize*(2*i - 1))  &&
	      dPhi < (halfBinSize*(2*i + 1))  
	      ){
	    AsymMatrix_sim[thBin][i] += 1;
	  }
	}
      }
    
  }//end of: for(Int_t ientry...
  
  
  Bool_t comments = kFALSE;
  
  if(comments){
    for(Int_t j = 0 ; j <nThbins; j++)
      for(Int_t k = 0 ; k < nPhibinsSim; k++)
	cout << " assym matrix for theta bin " 
	     << j 
	     << " and phi bin "               
	     << k
	     << " is "
	     << AsymMatrix_sim[j][k] 
	     << endl;
  }
  
  return 0;
} //end of CalculateAsymmetrySimScattered


void TSim::GraphAsymmetrySim(TString inputFileNumber1,
			     TString inputFileNumber2,
			     Int_t   nThetaSBins,
			     Float_t thetaSMin,
			     Float_t thetaSMax			      
			     ){
  
  cout << "-------------------" << endl;
  cout << " GraphAsymmetrySim " << endl;
  cout << "                   " << endl;

  this->SetStyle();
  
  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = P(90)/P(0) 
  Int_t  dPhiDiff = 90;
  
  Bool_t drawTheory = kFALSE;
  
  // Set Data types for files/analysis
  // One and Two
  Bool_t entangled[2];
  Bool_t polarised[2];
  
  // To Do: set file type via user input 
  //        implement in simlab.cc
  
  entangled[0] = kFALSE;
  entangled[1] = kFALSE;
  
  entangled[0] = kTRUE;
  //  entangled[1] = kTRUE;
  
  polarised[0] = kFALSE;
  polarised[1] = kFALSE;

  if(!entangled[0])
    polarised[0] = kTRUE;
  
  if(!entangled[1])
    polarised[1] = kTRUE;
  
  //Calculating ratios for desired dPhiDiff
  Float_t AsPhiDiff[nThbins] = {0.};
  Float_t AePhiDiff[nThbins] = {0.};
  
  Float_t AsPhiDiff1[nThbins] = {0.};
  Float_t AePhiDiff1[nThbins] = {0.};

  Float_t AsPhiDiffR[nThbins] = {0.};
  Float_t AePhiDiffR[nThbins] = {0.};
  
  Int_t inputFileInt1 = inputFileNumber1.Atoi();
  Int_t inputFileInt2 = inputFileNumber2.Atoi();
  
  TCanvas *canvas = new TCanvas("canvas","canvas",
				200,200,1400,1000);
  
  
  Double_t halfBinSize = 180.0/nPhibinsSim;
  
  Float_t thetaSBinWidth = 0.;
  thetaSBinWidth = (thetaSMax  - thetaSMin);
  thetaSBinWidth = thetaSBinWidth/(Float_t)nThetaSBins;
  
  Float_t thetaSHalf = thetaSBinWidth/2.;
  Float_t thetaS_arr[nThetaSBins];
  Float_t energyS_arr[nThetaSBins];
    
  for (Int_t i = 0 ; i < nThetaSBins ; i++){
    thetaS_arr[i]  = thetaSMin + (i+0.5)*thetaSBinWidth;
    energyS_arr[i] = ThetaToPhotonEnergy(thetaS_arr[i]); 
  }
  
  Int_t   nGraphs = 2;
  if( nThetaSBins!=0 )
    nGraphs = nThetaSBins + 1;
  
  TGraphErrors *grAsym[nGraphs];
  TGraphErrors *grAsymR;
  
  TGraphErrors *grA[nGraphs];
  TGraphErrors *grB[nGraphs];
  TGraphErrors *grC[nGraphs];
  
  
  for (Int_t g = 0 ; g < nGraphs ; g++){
    
    if     ( g == 0 ){
      // Set: AsymMatrix_sim[thBin][]
      CalculateAsymmetrySim(inputFileNumber1);
      //cout << " skipping " << endl;
    }
    else if( g >= 1 ){
      // Standard analysis, second file (option 4) 
      if(nThetaSBins == 0){
	CalculateAsymmetrySim(inputFileNumber2);
      }
      // Scattered analysis same file (option 5)
      else
	CalculateAsymmetrySimScattered(inputFileNumber1,
				       thetaS_arr[g-1],
				       thetaSHalf);
    }
    
    // now that AsymMatrix_sim[thBin][phiBin]
    // is set ...
    CalculateABC();
    
    // To Do:
    //CalculateABC_Err()
    
    for (Int_t i = 0 ; i < nThbins ; i++){
      pA[i] = pABC[i][0];
      pB[i] = pABC[i][1];
      pC[i] = pABC[i][2];
      
    }

    // Set: N(dPhi)/N(0)
    for (Int_t i = 0 ; i < nThbins ; i++){
      if (AsymMatrix_sim[i][0]== 0) continue;
      
      if (dPhiDiff  == 90){
	AsPhiDiff[i] = AsymMatrix_sim[i][bin090];
	AsPhiDiff[i] = AsPhiDiff[i]+AsymMatrix_sim[i][bin270];
	AsPhiDiff[i] = AsPhiDiff[i]/(2.*AsymMatrix_sim[i][0]);
	
	AePhiDiff[i] = AsymMatrix_sim[i][bin090];
	AePhiDiff[i] = AsymMatrix_sim[i][bin270];
	AePhiDiff[i] = 1./AePhiDiff[i];
	AePhiDiff[i] = AePhiDiff[i] + (1./AsymMatrix_sim[i][0]);
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt(AePhiDiff[i]);
      }
      else if (dPhiDiff  == 180){
	AsPhiDiff[i] = AsymMatrix_sim[i][bin180];
	AsPhiDiff[i] = AsPhiDiff[i]/AsymMatrix_sim[i][0];
	
	AePhiDiff[i] = (1./AsymMatrix_sim[i][bin180]);
	AePhiDiff[i] = AePhiDiff[i] +(1./AsymMatrix_sim[i][0]);
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt(AePhiDiff[i]);
      }
      else if (dPhiDiff  == 270){
	AsPhiDiff[i] = AsymMatrix_sim[i][bin270];
	AsPhiDiff[i] = AsPhiDiff[i]/AsymMatrix_sim[i][0];
	
	AePhiDiff[i] = 1./AsymMatrix_sim[i][bin270];
	AePhiDiff[i] = AePhiDiff[i] + 1./AsymMatrix_sim[i][0];
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt(AePhiDiff[i]);
      }	
      
      if(nThetaSBins==0){
	if     (g == 0){
	  AsPhiDiff1[i] = AsPhiDiff[i];
	  AePhiDiff1[i] = AePhiDiff[i];
	}
	else if(g == 1){
	  AsPhiDiffR[i] = AsPhiDiff1[i]/AsPhiDiff[i];
	  AePhiDiffR[i] = Sqrt(AePhiDiff[i]*AePhiDiff[i] + AePhiDiff1[i]*AePhiDiff1[i]);
	}
      }
      
    } // end of: for (Int_t i = 0 ; i < nThbins ; i+ 
    
    
    grAsym[g] = new TGraphErrors(nThbins,plotTheta,AsPhiDiff,0,AePhiDiff);
    
    if(g==1)
      grAsymR   = new TGraphErrors(nThbins,plotTheta,AsPhiDiffR,0,AePhiDiffR);
    
    grA[g]    = new TGraphErrors(nThbins,plotTheta,pA,0,0);
    grB[g]    = new TGraphErrors(nThbins,plotTheta,pB,0,0);
    grC[g]    = new TGraphErrors(nThbins,plotTheta,pC,0,0);
    
  } //end of : for (Int_t g = 0 ; g < nGr
  
  for( Int_t g = 0 ; g < nGraphs ; g++ ){
    grAsym[g]->SetLineColor( g+1 );
    grAsym[g]->SetMarkerColor( g+1 );

    grAsymR->SetLineColor( g+1 );
    grAsymR->SetMarkerColor( g+1 );
    
    grA[g]->SetLineColor( g+1 );
    grA[g]->SetLineStyle( g );
    grA[g]->SetMarkerColor( g+1 );
    
    grB[g]->SetLineColor( g+1 );
    grB[g]->SetLineStyle( g+1 );
    grB[g]->SetMarkerColor( g+1 );
    
    grC[g]->SetLineColor( g+1 );
    grC[g]->SetLineStyle( g+2 );
    grC[g]->SetMarkerColor( g+1 );
  }
  
  // Theory Curve
  Float_t aTheory[nThbins];
  
  // half resolution in phi: from delta phi binning
  Float_t alpha1 = DegToRad()*halfBinSize/Sqrt(2);
  
  // half resolution in theta
  Float_t semiSpan = DegToRad()*(ThMax[0] - ThMin[0])/2.;
  
  Char_t theoryLegendTitle[128];
  
  Char_t  plotN[128];
  TString plotName;
  
  Char_t yAxis[128];
  
  alpha1 = RadToDeg()*alpha1;
  
  sprintf(theoryLegendTitle,
	  "theory curve #alpha_{#phi} = %.1f^{o}",
	  alpha1);  
  
  Float_t maxY = 1.8;
  Float_t minY = 0.8;

  if( entangled[0] || 
      entangled[1]){
    minY =   0.0;
    maxY =   3.5;
  }
  
  TH1F *hr;
  hr = canvas->DrawFrame(10,minY,170,maxY);
  hr->GetXaxis()->SetTitle("#theta (deg)");
  sprintf(yAxis,"P(%d^{o})/P(0^{o})",dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);
  hr->GetYaxis()->SetTitleOffset(0.7);
  
  TLegend *leg = new TLegend(0.50,0.75,0.8,0.89);
  
  TTheory *theory = new TTheory();
  cout << endl;
  cout << " Calculating theory curve ... " << endl;
  cout << endl;
  cout << " semiSpan = " << semiSpan*RadToDeg() << endl;
  cout << " alpha1   = " << alpha1*RadToDeg()   << endl;
  
  for (Int_t i = 0; i < nThbins; i++){
    if( dPhiDiff == 180 ){
      aTheory[i] = 1.0;
      continue;
    }
    plotTheta[i] = plotTheta[i]*DegToRad();
    aTheory[i] = theory->rho2(plotTheta[i],semiSpan,alpha1);
    plotTheta[i] = plotTheta[i]*RadToDeg();
  }
  
  TGraphErrors* grThe = new TGraphErrors(nThbins,plotTheta,aTheory,0,0);
  grThe->SetLineColor(kRed);
  grThe->SetMarkerColor(kRed);
  
  if(drawTheory)
    leg->AddEntry(grThe,theoryLegendTitle,"L P");
  
  TString gLegTitle;
  gLegTitle = "";

  // first file
  if     (entangled[0])
    gLegTitle = "Entangled," + gLegTitle;
  else if(polarised[0])
    gLegTitle = "Polarised,"   + gLegTitle;
  else 
    gLegTitle = "UnPolarised," + gLegTitle;
  
  leg->AddEntry(grAsym[0],gLegTitle,"E P");
  
  grAsym[0]->Draw("P L E");
  
  // second+ file/s
  for (Int_t g = 1 ; g < nGraphs ; g++ ){
    gLegTitle = "";
    
      if     (entangled[1])
	gLegTitle = "Entangled,"   + gLegTitle;
      else if(polarised[1])
	gLegTitle = "Polarised,"   + gLegTitle;
      else 
	gLegTitle = "UnPolarised," + gLegTitle;
      
      if(nThetaSBins > 0.){
	gLegTitle.Form(gLegTitle + " #theta_{S} = %.0f^{o} (%.0f keV)",
		       thetaS_arr[g-1],
		       energyS_arr[g-1]);
	
      }
      
      leg->AddEntry(grAsym[g],gLegTitle,"E P");
      grAsym[g]->Draw("same P L E");
      
  }
  
  gStyle->SetTitleH(0.1);
  
  if(drawTheory)
    grThe->Draw("same P L");
  
  leg->Draw();
  
  if(nThetaSBins==0)
    sprintf(plotN,"../Plots/A_%d_Files%d_%d_%dPhiBins.pdf", 
	    dPhiDiff, inputFileInt1, inputFileInt2, nPhibinsSim);
  else
    sprintf(plotN,"../Plots/A_%d_Scattered_%d_%dPhiBins.pdf", 
	    dPhiDiff, inputFileInt1, nPhibinsSim);  
  
  canvas->SaveAs(plotN);
  
  //-----------------------------------------------------------------
  // Graphing Ratio of Asymms for Entangled,Polarised to Unpolarised
  
  leg->Clear();
  //leg->AddEntry(grAsymR,gLegTitle,"E P");
  
  hr = canvas->DrawFrame(10,minY,170,maxY);
  hr->GetXaxis()->SetTitle("#theta (deg)");
  sprintf(yAxis,"A^{E,P}/A^{E,U} #Delta #phi = %d",dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);
  hr->GetYaxis()->SetTitleOffset(0.7);
  
  grAsymR->Draw("P L E");

  plotName =  "../Plots/A_" + inputFileNumber1;
  plotName += "_";
  plotName += inputFileNumber2;
  plotName += ".pdf";
  
  canvas->SaveAs(plotName);
  
  // ------------------------------------------------
  //  Fourier Coefficients Graph
  
  leg->Clear();
  
  maxY =  0.3;
  minY = -0.3;
  if(entangled[0] ||
     entangled[1]){
    maxY =  0.5;
    minY = -0.5;
  }
  
  hr = canvas->DrawFrame(10,minY,170,maxY);
  hr->GetXaxis()->SetTitle("#theta (deg)");
  hr->GetYaxis()->SetTitle("cos(2#Delta#phi) coefficient");
  hr->GetYaxis()->SetTitleOffset(0.7);
  
  gLegTitle = "";

  // first file
  if     (entangled[0])
    gLegTitle = "Entangled,"   + gLegTitle;
  else if(polarised[0])
    gLegTitle = "Polarised,"   + gLegTitle;
  else 
    gLegTitle = "UnPolarised," + gLegTitle;
  
  leg->AddEntry(grC[0],gLegTitle,"L");
  
  //  grA[0]->Draw("L");
  //   grB[0]->Draw("L");
  grC[0]->Draw("L P");
  
  // second+ file/s
  for (Int_t g = 1 ; g < nGraphs ; g++ ){
    gLegTitle = "";
    if     (entangled[1])
      gLegTitle = "Entangled,"   + gLegTitle;
    else if(polarised[1])
      gLegTitle = "Polarised,"   + gLegTitle;
    else 
      gLegTitle = "UnPolarised," + gLegTitle;
    
    if(nThetaSBins > 0.){
      
      gLegTitle.Form(gLegTitle + " #theta_{S} = %.0f^{o} (%.0f keV)",
		     thetaS_arr[g-1],energyS_arr[g-1]);
      
    }
    
    leg->AddEntry(grC[g],gLegTitle,"L");
    //grA[g]->Draw("same P E");
    //grB[g]->Draw("same P E");
    grC[g]->Draw("same L P");
    
  }
  
  if(nThetaSBins==0)
    sprintf(plotN,"../Plots/FourierCoeffs_%d_%d_%dPhiBins.pdf", 
	    inputFileInt1, inputFileInt2, nPhibinsSim);
  else
    sprintf(plotN,"../Plots/FourierCoeffs_Scattered_%d_%dPhiBins.pdf", 
	    inputFileInt1, nPhibinsSim);  
  

  leg->Draw();
  
  canvas->SaveAs(plotN);
  
  cout << "                   " << endl;
  cout << " GraphAsymmetrySim " << endl;
  cout << "-------------------" << endl;
 

}//end of GraphAsymmetrySim



    
//-----------------------------------------------------------------------

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
  garyStyle->SetFrameFillColor(0);

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
  
  garyStyle->SetLegendBorderSize(0);
  garyStyle->SetLegendFillColor(0);
  garyStyle->SetLegendFont(42);
  garyStyle->SetLegendTextSize(0.);  
  
  gROOT->SetStyle("garyStyle");
  gROOT->ForceStyle();

}

// ------------------------------------------------------------------------------------------------
