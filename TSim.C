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

  
  //!! temporary 
  // SetAsymmetry(fileNumber);
  
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

  TTree* simDataTree =(TTree*)inputFile->Get("Tangle2");
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

    simPhiA[i] = 999;
    simPhiB[i] = 999;

    nb_ComptA[i] = 0;
    nb_ComptB[i] = 0;


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

  tempString.Form("simtHA[%d]/D",nCrystals);
  sortDataTree->Branch("simtHA", simtHA, tempString);

  tempString.Form("simtHB[%d]/D",nCrystals);
  sortDataTree->Branch("simtHB", simtHB, tempString);

  tempString.Form("simPhiA[%d]/D",nCrystals);
  sortDataTree->Branch("simPhiA", simPhiA, tempString);

  tempString.Form("simPhiB[%d]/D",nCrystals);
  sortDataTree->Branch("simPhiB", simPhiB, tempString);

  tempString.Form("nb_ComptA[%d]/I",nCrystals);
  sortDataTree->Branch("nb_ComptA", nb_ComptA, tempString);

  tempString.Form("nb_ComptB[%d]/I",nCrystals);
  sortDataTree->Branch("nb_ComptB", nb_ComptB, tempString);

  tempString.Form("XposA[%d]/D",nCrystals);
  sortDataTree->Branch("XposA", XposA, tempString);

  tempString.Form("YposA[%d]/D",nCrystals);
  sortDataTree->Branch("YposA", YposA, tempString);

  tempString.Form("ZposA[%d]/D",nCrystals);
  sortDataTree->Branch("ZposA", ZposA, tempString);

  tempString.Form("XposB[%d]/D",nCrystals);
  sortDataTree->Branch("XposB", XposB, tempString);

  tempString.Form("YposB[%d]/D",nCrystals);
  sortDataTree->Branch("YposB", YposB, tempString);

  tempString.Form("ZposB[%d]/D",nCrystals);
  sortDataTree->Branch("ZposB", ZposB, tempString);

 

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
 
      simtHA[j] = ThetaA_1st;
      simtHB[j] = ThetaB_1st;

      simPhiA[j] = PhiA_1st;
      simPhiB[j] = PhiB_1st;

      nb_ComptA[j] = nb_Compt[j];
      nb_ComptB[j] = nb_Compt[j + nCrystals];

      XposA[j] = XposA_1st;
      YposA[j] = YposA_1st;
      ZposA[j] = ZposA_1st;

      XposB[j] = XposB_1st;
      YposB[j] = YposB_1st;
      ZposB[j] = ZposB_1st;


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

Bool_t TSim::CentralY(Double_t posY){
  
  Bool_t centralY = kFALSE;
  
  //! works only for this particular detector geometry
  posY = Abs(posY);

  if( posY < 2.)
    centralY = kTRUE;

  return centralY;
}

Bool_t TSim::CentralZ(Double_t posZ){
  
  Bool_t centralZ = kFALSE;
  
  //! works only for this particular detector geometry
  posZ = Abs(posZ);

  if( posZ < 1.5)
    centralZ = kTRUE;

  return centralZ;
}

Float_t TSim::CrystalToPhi(Int_t crystal){

  //!! using only 5 crystals
  //works only for a particular detector geometry
  Float_t crystalToPhi[9] = {-1.,0.,-1.,270.,-1.,90.,-1.,180.,-1.};

  Float_t phi = crystalToPhi[crystal];

  return phi;
}

void TSim::SetAsymmetry(TString inputFileNumber){
  
  cout << endl;
  cout << " Setting asymmetry " << endl;
  
  this->CalculateAsymmetryLab(inputFileNumber);
  
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

  for (Int_t i = 0 ; i < nThbins ; i++){
    ThMin[i] = thetaLowEdge + i* thetaBinWidth;
    ThMax[i] = thetaLowEdge + (i+1)* thetaBinWidth;
    plotTheta[i] = thetaLowEdge + (i+0.5)* thetaBinWidth;
  }
}

Int_t TSim::CalculateAsymmetryLab(TString inputFileNumber){

  GetThetaBinValues();

  //difference 'sim theta' - 'energy theta'
  TH2F* histThD = new TH2F("histThD","sim - energy #theta difference",160,10,170,250,-120,130);

  //difference 'lab theta' - 'energy theta'
  TH2F* hLEdiff = new TH2F("hLEdiff","lab - energy theta",160,10,170,120,-50,70);

  //difference 'max/min lab theta' - 'energy theta'
  TH2F* hMEdiff = new TH2F("hMEdiff","max/min - energy theta",160,10,170,120,-50,70);

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

  n000 = 0;
  n090 = 0;
  n180 = 0;
  n270 = 0;

  for(Int_t j = 0 ; j <nThbins; j++){
    for(Int_t k = 0 ; k < 4; k++){
      AsymMatrix[j][k] = 0;}
  }

  for(Int_t j = 0 ; j <nThbins; j++){
    for(Int_t k = 0 ; k < 4; k++){
      AsymTrue[j][k] = 0;}
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
    Float_t totEA = 0;
    Float_t totEB = 0;

    for (Int_t k=0; k<nCrystals; k++){
      totEA += EA[k];
      totEB += EB[k];
    }
    
    if (GetThetaBin(ltHA[4]) == GetThetaBin(ltHB[4])&&
	(totEA > 450)&&(totEA < 550)&&(totEB > 450)&&(totEB < 550)){
      thBin = GetThetaBin(ltHA[4]);

      for (Int_t i = 0 ; i < nCrystals ; i++){
	if ((i!=4)&&(thBin>-1)&&(ltHA[i]<ThMax[thBin])&&(ltHA[i]>ThMin[thBin])){
	  phiA = CrystalToPhi(i);
	  indexA = i ;}
	if ((i!=4)&&(thBin>-1)&&(ltHB[i]<ThMax[thBin])&&(ltHB[i]>ThMin[thBin])){
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

      if ((simtHA[4]<ThMax[thBin])&&(simtHA[4]>ThMin[thBin])&&
	  (simtHB[4]<ThMax[thBin])&&(simtHB[4]>ThMin[thBin])&&
	  (CentralY(YposA[4]))&&(CentralZ(ZposA[4]))&&
	  (CentralY(YposB[4]))&&(CentralZ(ZposB[4]))&&
	  (nb_ComptA[4] == 1.)&&(nb_ComptB[4] == 1.)){
	if (phiDiff == 0.)
	  AsymTrue[thBin][0] += 1;
	if ((phiDiff == 90)||(phiDiff == -270))
	  AsymTrue[thBin][1] += 1;
	if ((phiDiff == 180)||(phiDiff == -180))
	  AsymTrue[thBin][2] += 1;
	if ((phiDiff == -90)||(phiDiff == 270))
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
				  10,10,1200,800);

  Char_t plotN[128];

  hLEdiff->Draw("colz");
  hLEdiff->GetXaxis()->SetTitle("energy #theta (deg)");
  hLEdiff->GetYaxis()->SetTitle("lab #theta - energy #theta (deg)");

  sprintf(plotN,"../Plots/histLEdiff_%d.pdf",inputFileInt);
  canvas->SaveAs(plotN);
   
  histThD->Draw("colz");
  histThD->GetXaxis()->SetTitle("#theta (deg)");
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

void TSim::GraphAsymmetryLab(TString inputFileNumber){

  GetThetaBinValues();

  this->SetStyle();

  Int_t inputFileInt = inputFileNumber.Atoi();
  
  TString plotName;
  plotName = "../Plots/Asym_" + inputFileNumber;

  plotName = plotName + ".pdf";

  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = P(90)/P(0) 
  Int_t   dPhiDiff = 90;

  //Calculating ratios for desired dPhiDiff
  Float_t AsPhiDiff[nThbins] = {0.};
  Float_t AePhiDiff[nThbins] = {0.};

  Float_t AsTrue[nThbins] = {0.};
  Float_t AeTrue[nThbins] = {0.};

  
  for (Int_t i = 0 ; i < nThbins ; i++){
    if (AsymMatrix[i][0] != 0){
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
  }

  //theory curve
  Float_t aTheory[nThbins];

  //half resolution in dPhi
  Float_t alpha1 = DegToRad()*35.0;

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

  TGraphErrors *grAsym[3];
  grAsym[0] = new TGraphErrors(nThbins,plotTheta,AsPhiDiff,0,AePhiDiff);
  grAsym[1] = new TGraphErrors(nThbins,plotTheta,aTheory,0,0);

  for(Int_t i = 0; i < nThbins; i++)
    plotTheta[i] +=2.;

  grAsym[2] = new TGraphErrors(nThbins,plotTheta,AsTrue,0,AeTrue);
  
  grAsym[0]->SetLineColor(kGreen+2.7);
  grAsym[0]->SetMarkerColor(kGreen+2.7);
  grAsym[1]->SetLineColor(kRed);
  grAsym[1]->SetMarkerColor(kRed);
  grAsym[2]->SetLineColor(kGreen);
  grAsym[2]->SetMarkerColor(kGreen);

  TLegend *leg = new TLegend(0.6,0.75,0.9,0.85);
  
  
  TH1F *hr;

  hr = canvas->DrawFrame(10,0.5,170,3);
  hr->GetXaxis()->SetTitle("#theta (deg)");
  
  TString theoryLegendTitle = " ";
  alpha1 = alpha1*RadToDeg();
  theoryLegendTitle.Form("theory curve #alpha_{#Delta#phi} = %.1f^{o}", alpha1);
  
  Char_t plotN[128];
  Char_t yAxis[128];
  
  sprintf(yAxis,"P(%d^{o})/P(0^{o})",dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);
  leg->AddEntry(grAsym[1],
		  theoryLegendTitle,"L P");
  leg->AddEntry(grAsym[0],"simulated lab", "E P");
  leg->AddEntry(grAsym[2],"true simulated lab","E P");
  grAsym[0]->Draw("P E");
  grAsym[1]->Draw("P L SAME");
  grAsym[2]->Draw("P E");
  leg->Draw();

  sprintf(plotN,"../Plots/A_%d_%d.pdf",inputFileInt,dPhiDiff);
  
  canvas->SaveAs(plotN);

}
//-----------------------------------------------------------------------
Int_t TSim::CalculateAsymmetrySim(TString inputFileNumber){

  GetThetaBinValues();

  cout << endl;
  cout << " Getting asymmetry sim " << endl;
  
  this->SetStyle();

  Int_t inputFileInt = inputFileNumber.Atoi();
  
  TString plotName;
  plotName = "../Plots/Asym_sim_" + inputFileNumber;

  plotName = plotName + ".pdf";

  TString inputFileName = "../Data/sim" + inputFileNumber;
  inputFileName  = inputFileName + ".root";
    
  cout << endl;
  cout << " Input File : " << inputFileName  << endl;
  cout << endl;

  TFile* inputFile = new TFile(inputFileName);

  TTree* simDataTree=(TTree*)inputFile->Get("Tangle2");

  simDataTree->SetBranchAddress("ThetaA_1st",&ThetaA_1st);
  simDataTree->SetBranchAddress("ThetaB_1st",&ThetaB_1st);
  simDataTree->SetBranchAddress("PhiA_1st",&PhiA_1st);
  simDataTree->SetBranchAddress("PhiB_1st",&PhiB_1st);

 

   for(Int_t j = 0 ; j <nThbins; j++){
    for(Int_t k = 0 ; k < nPhibinsSim; k++){
      AsymMatrix_sim[j][k] = 0;}
  }

   Long64_t nEntries = simDataTree->GetEntries();

Float_t binsize = 45;//180/nPhibinsSim; //45
  //Event loop
  for (Int_t ientry = 0 ; ientry < nEntries ; ientry++){
    simDataTree->GetEvent(ientry);
     Float_t dPhi_1st = PhiA_1st + PhiB_1st;
  if(dPhi_1st<0){
    dPhi_1st = dPhi_1st + 360; }
    

    Int_t thBin = -1;
    
    if (GetThetaBin(ThetaA_1st) == GetThetaBin(ThetaB_1st)){
	thBin = GetThetaBin(ThetaA_1st);
      }


      if(thBin>-1){
	
	//fill matrix for dphi=0 bin
	if((dPhi_1st<binsize)||(dPhi_1st>360-binsize)){
	 
	    AsymMatrix_sim[thBin][0] += 1;}
	//fill the rest of the matrix 
	for (Int_t i = 1 ; i < nPhibinsSim ; i++){
	  if((dPhi_1st>2*i*binsize-binsize)&&(dPhi_1st<2*i*binsize+binsize)){
	    AsymMatrix_sim[thBin][i] += 1;}
	}
      }

      }//end of: for(Int_t ientry...
  
  //printing out the asym matrix 
    for(Int_t j = 0 ; j <nThbins; j++){
      for(Int_t k = 0 ; k < nPhibinsSim; k++){
	cout<<"assym matrix for theta bin "<<j<<" and phi bin "<<k<<" is "<<AsymMatrix_sim[j][k]<<endl;}
    }
  

} //end of CalculateAsymmetrySim




Int_t TSim::GraphAsymmetrySim(TString inputFileNumber1, TString inputFileNumber2){
  
  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = P(90)/P(0) 
  Int_t   dPhiDiff = 90;
    
  //Calculating ratios for desired dPhiDiff
  Float_t AsPhiDiff[nThbins] = {0.};
  Float_t AePhiDiff[nThbins] = {0.};

  Int_t inputFileInt1 = inputFileNumber1.Atoi();
  Int_t inputFileInt2 = inputFileNumber2.Atoi();
  
  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);

  int bin90 = nPhibinsSim/4;
  int bin180 = nPhibinsSim/4;
  int bin270 = 3*nPhibinsSim/4;
  
  
  
  CalculateAsymmetrySim(inputFileNumber1);

  for (Int_t i = 0 ; i < nThbins ; i++){
    if (AsymMatrix_sim[i][0] != 0){
      if (dPhiDiff  == 90){
	//using average 90 and 270
	AsPhiDiff[i] = (AsymMatrix_sim[i][bin90]+AsymMatrix_sim[i][bin270])/(2*AsymMatrix_sim[i][0]);
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/(AsymMatrix_sim[i][bin90]+AsymMatrix_sim[i][bin270]))+(1/AsymMatrix_sim[i][0]));
	  	}
      
      if (dPhiDiff  == 180){
	AsPhiDiff[i] = AsymMatrix_sim[i][bin180]/AsymMatrix_sim[i][0];
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix_sim[i][bin180])+(1/AsymMatrix_sim[i][0]));
      }
      if (dPhiDiff  == 270){
	AsPhiDiff[i] = AsymMatrix_sim[i][bin270]/AsymMatrix_sim[i][0];
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix_sim[i][bin270])+(1/AsymMatrix_sim[i][0]));
	
      }	
    }
  }
  
  
  TGraphErrors *grAsym1 = new TGraphErrors(nThbins,plotTheta,AsPhiDiff,0,AePhiDiff);
 

  
  CalculateAsymmetrySim(inputFileNumber2);
  
  for (Int_t i = 0 ; i < nThbins ; i++){
    if (AsymMatrix_sim[i][0] != 0){
      if (dPhiDiff  == 90){
	//using average 90 and 270
	AsPhiDiff[i] = (AsymMatrix_sim[i][bin90]+AsymMatrix_sim[i][bin270])/(2*AsymMatrix_sim[i][0]);
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/(AsymMatrix_sim[i][bin90]+AsymMatrix_sim[i][bin270]))+(1/AsymMatrix_sim[i][0]));
	  	}
      
      if (dPhiDiff  == 180){
	AsPhiDiff[i] = AsymMatrix_sim[i][bin180]/AsymMatrix_sim[i][0];
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix_sim[i][bin180])+(1/AsymMatrix_sim[i][0]));
      }
      if (dPhiDiff  == 270){
	AsPhiDiff[i] = AsymMatrix_sim[i][bin270]/AsymMatrix_sim[i][0];
	AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix_sim[i][bin270])+(1/AsymMatrix_sim[i][0]));
	
      }	
    }
    cout<<AePhiDiff[i]<<endl;
  }
   
  TGraphErrors *grAsym2 =
    new TGraphErrors(nThbins,plotTheta,AsPhiDiff,0,AePhiDiff);
   

    //theory curve
  Float_t aTheory[nThbins];

  //half resolution in dPhi
  Float_t alpha1 = DegToRad()*40.0;

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

    TGraphErrors* grThe = new TGraphErrors(nThbins,plotTheta,aTheory,0,0);


  grAsym1->SetLineColor(kBlue);
  grAsym1->SetMarkerColor(kBlue);

  grAsym2->SetLineColor(kGreen);
  grAsym2->SetMarkerColor(kGreen);

  grThe->SetLineColor(kRed);
  grThe->SetMarkerColor(kRed);


  TH1F *hr;

  Char_t plotN[128];
  Char_t yAxis[128];

  hr = canvas->DrawFrame(10,0.5,170,3);
  hr->GetXaxis()->SetTitle("#theta (deg)");
  
  sprintf(yAxis,"P(%d^{o})/P(0^{o})",dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);


  grAsym1->Draw("P E");
  grAsym2->Draw("same P E");
  grThe->Draw("same P L");

  sprintf(plotN,"../Plots/A_%d_%d_%d.pdf",inputFileInt1, inputFileInt2, dPhiDiff);
  
  canvas->SaveAs(plotN);
 
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
