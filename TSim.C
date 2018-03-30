#include "TSim.h"

#if !defined(__CINT__)
ClassImp(TSim)
#endif

TSim::TSim(){
}

TSim::TSim(TString fileNumber){
  
  //!!!!
  InitHs();
  
  cout << endl;
  cout << " Constructing TSim object " << endl;
  
  rootFileRawName = "../Data/sim" + fileNumber;
  rootFileSortName = "../Data/sort" + fileNumber;
  
  rootFileRawName = rootFileRawName + ".root";
  rootFileSortName = rootFileSortName + ".root";
  
}

TSim::TSim(TString fileNumber1, TString fileNumber2){
  
  //!!
  InitHs();
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
  
}


TSim::~TSim(){
}

//----------------------------------------------

/** Public member functions *********/


Bool_t TSim::SortedROOTFileExists(){
  TFile *file = TFile::Open(rootFileSortName);
  
  return file;
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
    if ( (theta > ThMin[i])  &&
	 (theta < ThMax[i]) )
      bin = i;
  }
  
  return bin;
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


   //============================================
   // write to
   sortDataTree = new TTree("sortDataTree",
			    "Sorted simulation events");

   // New variables
   TString tempString = "";

   tempString.Form("EA[%d]/D", nCrystals);
   sortDataTree->Branch("EA", EA, tempString);

   tempString.Form("EB[%d]/D", nCrystals);
   sortDataTree->Branch("EB", EB, tempString);

   tempString.Form("EAX[%d]/D", nCrystals);
   sortDataTree->Branch("EAX", EAX, tempString);

   tempString.Form("EBX[%d]/D", nCrystals);
   sortDataTree->Branch("EBX", EBX, tempString);


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

     // smeared energy
     EA[i] = 0; 
     EB[i] = 0;

     // exact energy    
     EAX[i] = 0; 
     EBX[i] = 0;

     // angle calculated from exact energy deposit = 'energy theta'
     etHA[i] = 0;
     etHB[i] = 0.;

     // angle calculated from the smeared energy deposit = 'lab theta'
     ltHA[i] = 0;
     ltHB[i] = 0.;

     // angle calculated from exact energy - a*sigma
     mintHAErr[i] = 0.;
     mintHBErr[i] = 0.;

     // angle calculated from exact energy + a*sigma
     maxtHAErr[i] = 0.;
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

   //============
   // EVENT LOOP
   Float_t sigmaPar = 2.5;

   Long64_t nEvents = simDataTree->GetEntries();

   TRandom3 * rand3 = new TRandom3(); 

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

     Int_t k = 0;
     
     for (Int_t j = 0 ; j < nCrystals ; j++){

       // sigmaE = sigma511 * sqrt( E / 511) 
       
       // index for array B
       k = j+nCrystals;

       EAX[j] = CrystEnergyDep[j];
       EBX[j] = CrystEnergyDep[k];

       EA[j] = rand3->Gaus(CrystEnergyDep[j],
			   sigmaA[j]*Sqrt(CrystEnergyDep[j]));
       EB[j] = rand3->Gaus(CrystEnergyDep[k],
			   sigmaB[j]*Sqrt(CrystEnergyDep[k]));

       etHA[j] = PhotonEnergyToTheta(CrystEnergyDep[j]);
       etHB[j] = PhotonEnergyToTheta(CrystEnergyDep[k]);
       ltHA[j] = PhotonEnergyToTheta(EA[j]);
       ltHB[j] = PhotonEnergyToTheta(EB[j]);

       // max/min defined as sigmaPar* sigma away from the mean
       maxtHAErr[j] = PhotonEnergyToTheta(CrystEnergyDep[j] - 
					  (sigmaPar*sigmaA[j]*
					   Sqrt(CrystEnergyDep[j])));
       mintHAErr[j] = PhotonEnergyToTheta(CrystEnergyDep[j] + 
					  (sigmaPar*sigmaA[j]*
					   Sqrt(CrystEnergyDep[j])));
       maxtHBErr[j] = PhotonEnergyToTheta(CrystEnergyDep[k] - 
					  (sigmaPar*sigmaB[j]*
					   Sqrt(CrystEnergyDep[k])));
       mintHBErr[j] = PhotonEnergyToTheta(CrystEnergyDep[k] +
					  (sigmaPar*sigmaB[j]*
					   Sqrt(CrystEnergyDep[k])));

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

     maxtHAErr[4] = ElectronEnergyToTheta(CrystEnergyDep[4] +
					  (sigmaPar*sigmaA[4]*
					   Sqrt(CrystEnergyDep[4])));
     mintHAErr[4] = ElectronEnergyToTheta(CrystEnergyDep[4] -
					  (sigmaPar*sigmaA[4]*
					   Sqrt(CrystEnergyDep[4])));
     maxtHBErr[4] = ElectronEnergyToTheta(CrystEnergyDep[13] +
					  (sigmaPar*sigmaB[4]*
					   Sqrt(CrystEnergyDep[13])));
     mintHBErr[4] = ElectronEnergyToTheta(CrystEnergyDep[13] - 
					  (sigmaPar*sigmaB[4]*
					   Sqrt(CrystEnergyDep[13])));				 

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

 } 

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

   Float_t crystalHalfSizeYZ =  1.0;

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

   Float_t crystalToPhi[9] = { -1.  ,   0. , -1.,
			       90. ,  -1. , 270.,
			       -1.  , 180. , -1.};

 //    !!Temp - Use corner crystals
 //   Float_t crystalToPhi[9] = {  0.  ,  -1  , 90. ,
 //   			       -1. ,  -1. , -1. ,
 //   			       270.,  -1. , 180. };

   return crystalToPhi[crystal];
 }

 // Transformation below matches the
 // simulation with beam in X and
 // arbitrary reference in Z (lab 'up' direction)
 // If using full PET array then use
 // arbitrary reference as scanner axis
 // so it is never in beam direction

Float_t TSim::CrystalToPhiA(Int_t crystal){
  
  Float_t crystalToPhiA[9] = { -1.  ,  0.  ,  -1.,
			       -90. ,  -1. ,  90.,
			       -1.  , 180. , -1.};
  
 //   // corner crystals
 //   Float_t crystalToPhiA[9] = { 45.  , -1.  , -45.,
 // 			       -1.  , -1. ,  -1.,
 // 			       135. , -1. , -135.};

return crystalToPhiA[crystal];
 }

Float_t TSim::CrystalToPhiB(Int_t crystal){

   Float_t crystalToPhiB[9] = { -1. ,  0.  , -1.,
				90. ,  -1. , -90.,
				-1. , 180. , -1.};
   

 //   Float_t crystalToPhiB[9] = { -45. , -1. , 45.,
 // 			       -1.  , -1. , -1.,
 // 			       -135., -1. , 135.};

   return crystalToPhiB[crystal];
 }

 Float_t TSim::GetAsymLab(Int_t dPhiDiff, Int_t i){

   Float_t AsPhiDiff = 0;
   if (AsymMatrix[i][0] != 0){
       if (dPhiDiff  == 90){
	 AsPhiDiff = AsymMatrix[i][1]/AsymMatrix[i][0];
	 //using average 90 and 270
	 //AsPhiDiff = (AsymMatrix[i][1]+AsymMatrix[i][3])/(2*AsymMatrix[i][0]);
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

	 // using average 90 and 270
	 //AsPhiDiff = (AsymMatrix[i][1]+AsymMatrix[i][3])/(2*AsymMatrix[i][0]);
	 //AePhiDiff = AsPhiDiff*Sqrt((1/(AsymMatrix[i][1]+AsymMatrix[i][3]))+(1/AsymMatrix[i][0]));
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
	 //AsTrue = (AsymTrue[i][1]+AsymTrue[i][3])/(2*AsymTrue[i][0]);
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
	 // AsTrue = (AsymTrue[i][1]+AsymTrue[i][3])/(2*AsymTrue[i][0]);
 // 	AeTrue = AsTrue*Sqrt((1/(AsymTrue[i][1]+AsymTrue[i][3]))+(1/AsymTrue[i][0]));
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

//!!! temporary
void TSim::InitHs(){


   for( Int_t th = 0 ; th < nThbins ; th++){

     TString hTitle = "";

     hTitle.Form("hDPhi_%d_0",th);
     hDPhi[th][0] = new TH1F(hTitle,hTitle,
			     32, 0.,360.);
     hTitle.Form("hDPhi_%d_1",th);
     hDPhi[th][1] = new TH1F(hTitle,hTitle,
			     32, 0.,360.);

     hTitle.Form("hDPhi_2_%d_0",th);
     hDPhi_2[th][0] = new TH1F(hTitle,hTitle,
				32, 0.,360.);

     hTitle.Form("hDPhi_2_%d_1",th);
     hDPhi_2[th][1] = new TH1F(hTitle,hTitle,
				32, 0.,360.);
   }


 }

Float_t TSim::CalculateR(Float_t y,
			 Float_t z){
  return Sqrt(y*y+z*z);
}

Float_t TSim::CalculateBeta(Float_t x,
			    Float_t y,
			    Float_t z){
  return (ATan(CalculateR(y,z)/x)*RadToDeg());
}

Float_t TSim::GetDeltaPhi(Float_t phiA,
			 Float_t phiB){
  
  Float_t dPhi = phiB + phiA;
  
  // shift (-360,0) to (0,360.)
  if(dPhi < 0.)
    dPhi += 360.;
  
  // 360 appears in one lablike combination
  // phiA = 180, phiB = 180
  if(dPhi==360)
    dPhi = 0;
  
  // for QETlab simulation:
  //... see Tangle2SteppingAction.cc

  return dPhi;
}

Float_t TSim::GetDeltaPhiRes(Float_t dPhiX,
			     Float_t dPhiL){

  // centre on 0 deg
  Float_t dPhiRes = dPhiX - dPhiL;
  
  if     (dPhiRes <= -180.0 )
    dPhiRes += 360.;
  else if(dPhiRes > 180.0 )
    dPhiRes -= 360.;
  
  if( dPhiRes <  -180.0 ||
      dPhiRes >=  180.) {
    cout << endl;
    cout << " dPhiX   = " <<  dPhiX   <<  endl;
    cout << " dPhiL   = " <<  dPhiL   <<  endl;
    cout << " dPhiRes = " <<  dPhiRes <<  endl;
  }
  
  return dPhiRes;
}


Float_t TSim::GetPhiRes(Float_t phiX,
			Float_t phiL){

  Float_t phiRes = phiX - phiL;

  if     (phiRes < -180.0 )
    phiRes += 360.;
  else if(phiRes >= 180.0 )
    phiRes -= 360.;
  
  if( phiRes < -180 ||
      phiRes > 180 ){
    cout << endl;
    cout << " phiX   =  " << phiX   << endl;
    cout << " phiL   =  " << phiL   << endl;
    cout << " phiRes =  " << phiRes << endl;
  }
    
  return phiRes;
}


Float_t TSim::ShiftPhi(Float_t phiX,
		       Float_t lowEdge){

  if( phiX < lowEdge)
    phiX += 360.0;
  
  return phiX;
}


Int_t TSim::InvestigateAcceptance(TString inputFileNumber){

   GetThetaBinValues();

   //------------------------------
   // Acceptance and Bad Event 
   // Reconstruction Investigation
   
   // compare exact theta to 
   // theta from unsmeared energy
   // inner crystal / theta
   
   TH2F* hThExaSubThEI = new TH2F("hThExaSubThEI",
				  "hThExaSubThEI",
				  160,10,170,
				  250,-120,130);

   //   as above with extra conditions
   TH2F* hThExaSubThEI_2 = new TH2F("hThExaSubThEI_2",
				    "hThExaSubThEI_2",
				    160,10,170,
				    250,-120,130);

   //   outer crystal / theta
   TH2F* hThExaSubThEO = new TH2F("hThExaSubThEO",
				  "hThExaSubThEO",
				  160,10,170,
				  250,-120,130);
   
   TH2F* hThExaSubThEO_2 = new TH2F("hThExaSubThEO_2",
				    "hThExaSubThEO_2",
				    160,10,170,
				    250,-120,130);
   
   // compare exact theta to 
   // 'lablike' theta from smeared energy
   // inner crystal / theta
   TH2F* hThExaSubThLI = new TH2F("hThExaSubThLI",
				  "hThExaSubThLI",
				  160,10,170,
				  250,-120,130);

   TH2F* hThExaSubThLI_2 = new TH2F("hThExaSubThLI_2",
				    "hThExaSubThLI_2",
				    160,10,170,
				    250,-120,130);
   
   // compare energy calculated theta to 
   // 'lablike' theta from smeared energy
   // inner crystal / theta
   TH2F* hThESubThLI = new TH2F("hThESubThLI",
				"hThESubThLI",
				160,10,170,
				120,-50,70);

   
   // 'max/min lab theta' - 'energy theta'
   TH2F* hMaxThdiff = new TH2F("hMaxThdiff",
			      "max/min - energy theta",
			      160,10,170,
			      120,-50,70);
   
   // beta is beam angle wrt x-axis
   TH2F* hBetaVsDPhi = new TH2F("hBetaVsDPhi",
				"hBetaVsDPhi",
				32,-10.0,370,
				32,0.,10.);
   //  R = sqrt( Y*Y + Z*Z )
   TH2F* hRVsDPhi = new TH2F("hRVsDPhi",
			     "hRVsDPhi",
			     32,-10.0,370,
			     32,-0.5,6.5);
   
   TH2F * hYZ[nThbins], * hYZ_2[nThbins];
   TH2F * hXY[nThbins], * hXY_2[nThbins];

   TH2F * hYZ_000[nThbins];
   TH2F * hYZ_090[nThbins];
   TH2F * hYZ_180[nThbins];
   TH2F * hYZ_270[nThbins];

   TH2F * hYZ_dPhi[4], * hYZ_dPhi_2[4];
   
   TH1F * hBeta[nThbins];
   TH1F * hBeta_2[nThbins];
   
   //  exact phi for final sample
   TH1F * hPhiA;
   TH1F * hPhiB;
   
   TH1F * hDPhiRes00[nThbins];
   TH1F * hDPhiRes90[nThbins];
   
   TH1F * hDPhiRes00_2[nThbins];
   TH1F * hDPhiRes90_2[nThbins];

   TH1F * hPhi000Res[nThbins];
   TH1F * hPhi090Res[nThbins];
   TH1F * hPhi180Res[nThbins];
   TH1F * hPhi270Res[nThbins];

   TH1F * hPhi000Res_2[nThbins];
   TH1F * hPhi090Res_2[nThbins];
   TH1F * hPhi180Res_2[nThbins];
   TH1F * hPhi270Res_2[nThbins];

   TH1F * hPhiRes000[nThbins];
   TH1F * hPhiRes090[nThbins];
   TH1F * hPhiRes180[nThbins];
   TH1F * hPhiRes270[nThbins];

   TH1F * hPhiRes000_2[nThbins];
   TH1F * hPhiRes090_2[nThbins];
   TH1F * hPhiRes180_2[nThbins];
   TH1F * hPhiRes270_2[nThbins];

   TH1F * hThRes00[nThbins];
   TH1F * hThRes90[nThbins];

   TH1F * hThRes00_2[nThbins];
   TH1F * hThRes90_2[nThbins];


   TH1F * hBeta000[nThbins], * hBeta000_2[nThbins];
   TH1F * hBeta090[nThbins], * hBeta090_2[nThbins];
   TH1F * hBeta180[nThbins], * hBeta180_2[nThbins];
   
   //------------------------
   // Initialise Histos

   hPhiA = new TH1F("hPhiA","hPhiA",
		    72, -135.,225.);
   
   hPhiB = new TH1F("hPhiB","hPhiB",
		    72, -135.,225.);

   TString hTitle = "";

   Int_t   nYZBins = 25;
   Float_t maxYZ   = 6.0*(1. + 2./(nYZBins-1));
   Float_t minYZ   = -maxYZ; 
   
   for( Int_t th = 0 ; th < nThbins ; th++){

     hTitle.Form("hDPhiRes00_%d",th);
     hDPhiRes00[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);

     hTitle.Form("hDPhiRes90_%d",th);
     hDPhiRes90[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);

     hTitle.Form("hDPhiRes00_2_%d",th);
     hDPhiRes00_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);

     hTitle.Form("hDPhiRes90_2_%d",th);
     hDPhiRes90_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);
     
     hTitle.Form("hPhiRes000_%d",th);
     hPhiRes000[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);

     hTitle.Form("hPhiRes090_%d",th);
     hPhiRes090[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);

     hTitle.Form("hPhiRes180_%d",th);
     hPhiRes180[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);

     hTitle.Form("hPhiRes270_%d",th);
     hPhiRes270[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);

     hTitle.Form("hPhiRes000_2_%d",th);
     hPhiRes000_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);

     hTitle.Form("hPhiRes090_2_%d",th);
     hPhiRes090_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);

     hTitle.Form("hPhiRes180_2_%d",th);
     hPhiRes180_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);

     hTitle.Form("hPhiRes270_2_%d",th);
     hPhiRes270_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);

     hTitle.Form("hPhi000Res_%d",th);
     hPhi000Res[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);

     hTitle.Form("hPhi090Res_%d",th);
     hPhi090Res[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);
     
     hTitle.Form("hPhi180Res_%d",th);
     hPhi180Res[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);

     hTitle.Form("hPhi270Res_%d",th);
     hPhi270Res[th] = new TH1F(hTitle,hTitle,
			       32, -180.,180.);

     hTitle.Form("hPhi000Res_2_%d",th);
     hPhi000Res_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);

     hTitle.Form("hPhi090Res_2_%d",th);
     hPhi090Res_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);

     hTitle.Form("hPhi180Res_2_%d",th);
     hPhi180Res_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);

     hTitle.Form("hPhi270Res_2_%d",th);
     hPhi270Res_2[th] = new TH1F(hTitle,hTitle,
				  32, -180.,180.);
     
     hTitle.Form("hThRes00_%d",th);
     hThRes00[th] = new TH1F(hTitle,hTitle,
			     128, -10., 180.);

     hTitle.Form("hThRes90_%d",th);
     hThRes90[th] = new TH1F(hTitle,hTitle,
			     128, -10., 180.);

     hTitle.Form("hThRes00_2_%d",th);
     hThRes00_2[th] = new TH1F(hTitle,hTitle,
				128, -10., 180.);

     hTitle.Form("hThRes90_2_%d",th);
     hThRes90_2[th] = new TH1F(hTitle,hTitle,
				128, -10., 180.);
     
     hTitle.Form("hBeta_%d",th);
     hBeta[th] = new TH1F(hTitle,hTitle,
			  32,-0.0, 10.0);
     
     hTitle.Form("hBeta_2_%d",th);
     hBeta_2[th] = new TH1F(hTitle,hTitle,
			     32,-0.0, 10.0);

     hTitle.Form("hBeta000_%d",th);
     hBeta000[th] = new TH1F(hTitle,hTitle,
			     32,-0.0, 10.0);

     hTitle.Form("hBeta090_%d",th);
     hBeta090[th] = new TH1F(hTitle,hTitle,
			     32,-0.0, 10.0);

     hTitle.Form("hBeta180_%d",th);
     hBeta180[th] = new TH1F(hTitle,hTitle,
			     32,-0.0, 10.0);

     hTitle.Form("hBeta000_2_%d",th);
     hBeta000_2[th] = new TH1F(hTitle,hTitle,
				32,-0.0, 10.0);

     hTitle.Form("hBeta090_2_%d",th);
     hBeta090_2[th] = new TH1F(hTitle,hTitle,
				32,-0.0, 10.0);

     hTitle.Form("hBeta180_2_%d",th);
     hBeta180_2[th] = new TH1F(hTitle,hTitle,
				32,-0.0, 10.0);

     hTitle.Form("hXY_%d",th);
     hXY[th] = new TH2F(hTitle,hTitle,
			20, 25.0, 55.0,
			nYZBins, minYZ, maxYZ);

     hTitle.Form("hYZ_%d",th);
     hYZ[th] = new TH2F(hTitle,hTitle,
			nYZBins, minYZ, maxYZ,
			nYZBins, minYZ, maxYZ);

     hTitle.Form("hYZ_000_%d",th);
     hYZ_000[th] = new TH2F(hTitle,hTitle,
			    nYZBins, minYZ, maxYZ,
			    nYZBins, minYZ, maxYZ);

     hTitle.Form("hYZ_090_%d",th);
     hYZ_090[th] = new TH2F(hTitle,hTitle,
			    nYZBins, minYZ, maxYZ,
			    nYZBins, minYZ, maxYZ);

     hTitle.Form("hYZ_180_%d",th);
     hYZ_180[th] = new TH2F(hTitle,hTitle,
			    nYZBins, minYZ, maxYZ,
			    nYZBins, minYZ, maxYZ);

     hTitle.Form("hYZ_270_%d",th);
     hYZ_270[th] = new TH2F(hTitle,hTitle,
			    nYZBins, minYZ, maxYZ,
			    nYZBins, minYZ, maxYZ);
     
     hTitle.Form("hXY_2_%d",th);
     hXY_2[th] = new TH2F(hTitle,hTitle,
			  10, 25.0, 55.0,
			  nYZBins, minYZ, maxYZ);
     
     hTitle.Form("hYZ_2_%d",th);
     hYZ_2[th] = new TH2F(hTitle,hTitle,
			  nYZBins, minYZ, maxYZ,
			  nYZBins, minYZ, maxYZ);


   }

   for (Int_t p = 0 ; p < 4 ; p ++){

     Int_t dp = p*90;

     hTitle.Form("hYZ_dPhi_%d",dp);
     hYZ_dPhi[p] = new TH2F(hTitle,hTitle,
			    41,-4.1, 4.1,
			    41,-4.1, 4.1);
     
     hTitle.Form("hYZ_dPhi_2_%d",dp);
     hYZ_dPhi_2[p] = new TH2F(hTitle,hTitle,
			      16,-1.55, 1.55,
			      16,-1.55, 1.55);
   }

   cout << endl;
   cout << " Getting asymmetry " << endl;
   
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

   sortDataTree->SetBranchAddress("EAX",EAX);
   sortDataTree->SetBranchAddress("EBX",EBX);

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

   for(Int_t j = 0 ; j <nThbins; j++){
     for(Int_t k = 0 ; k < 4; k++){
       AsymMatrix[j][k] = 0;}
   }

   for(Int_t j = 0 ; j < nThbins; j++){
     for(Int_t k = 0 ; k < 4; k++){
       AsymTrue[j][k] = 0;}
   }

   Long64_t nEntries = sortDataTree->GetEntries();

   Int_t   thBin  = -1;
   Int_t   thBinA = -1;
   Int_t   thBinB = -1;
   Float_t phiA   = -1;
   Float_t phiB   = -1;
   Int_t   indexA = -1;
   Int_t   indexB = -1;
   Float_t totEA  = 0;
   Float_t totEB  = 0;

   Float_t totEAX  = 0;
   Float_t totEBX  = 0;

   Float_t dPhiLab = -999.;
   
   Float_t dPhiXact    = -999.;
   Float_t dPhiXactRes = -999.;
   
   Float_t betaA = -99;
   Float_t betaB = -99;
   Float_t rA    = -99;

   Float_t simPhiARes = -999.;

   Int_t nDuplicatesAB = 0;
   Int_t nDuplicatesB  = 0;
   Int_t nDuplicatesA  = 0;

   for (Int_t i = 0 ; i < nThbins ; i++){
     hPhi000Res[i]->SetMinimum(0);
     hPhi090Res[i]->SetMinimum(0);
   }
   
   hPhiA->SetMinimum(0);
   hPhiB->SetMinimum(0);

   // Some histograms
   // have multiples
   // per file e.g.
   // polarised/unpolarised
   fileNum++;
   cout << endl;
   cout << " fileNum = " << fileNum << endl;
   
   // Event loop
   for (Int_t ientry = 0 ; ientry < nEntries ; ientry++){
     
     cout << " start " << endl;
     
     sortDataTree->GetEvent(ientry);
     
     totEA  = 0;
     totEB  = 0;

     totEAX  = 0;
     totEBX  = 0;

     for (Int_t k = 0; k < nCrystals; k++){
       totEA  += EA[k];
       totEB  += EB[k];
       totEAX += EAX[k];
       totEBX += EBX[k];
     }

     thBinA = -1;
     thBinB = -1;
     thBin  = -1;
     phiA   = -1;
     phiB   = -1;
     indexA = -1;
     indexB = -1;
     betaA = -99;
     betaB = -99;
     rA    = -99;
     
     //----------------------------------------
     //----------------------------------------
     //----------------------------------------
     // True Data Analysis
     //----------------------------------------
     //----------------------------------------
     //----------------------------------------
     
     //thBinA = GetThetaBin(simtHA[0]);
     //thBinB = GetThetaBin(simtHB[0]);
     
     thBinA = GetThetaBin(etHA[4]);
     thBinB = GetThetaBin(etHB[4]);
     
     // To do: tidy up this to test
     // exact variables distributions
     if( (totEA > 450 )  &&
	 (totEA < 550 )  &&
	 (thBinA != -1 ) &&
	 (totEB > 450 )  &&
	 (totEB < 550 )  &&
	 (thBinB != -1 ) &&
	 (thBinA == thBinB) 
	 ){

       for (Int_t i = 0 ; i < nCrystals ; i++){

	 if( i == 4 || i == 0 || i == 2 ||
	     i == 6 || i == 8 ) continue;

	 if ( ( etHA[i] < ThMax[thBinA])  &&
	      ( etHA[i] > ThMin[thBinA])  
	      ){
	   hPhi000Res[thBinA]->Fill(simPhiA[0]);
	 }

	 if ( ( etHB[i] < ThMax[thBinB])  &&
	      ( etHB[i] > ThMin[thBinB]) 
	      ){
	   hPhi090Res[thBinB]->Fill(simPhiB[0]);
	 }

       } //end of :for (Int_t i = 0 ; i < nCr...
     }
     
     //----------------------------------------
     //----------------------------------------
     //----------------------------------------
     // Lablike Analysis
     //----------------------------------------
     //----------------------------------------
     //----------------------------------------
     
      
     
     //----------------------------------------
     // only interested in events with
     // full energy deposited
     if( (totEA < 450) || 
	 (totEA > 550) ||
	 (totEB < 450) ||
	 (totEB > 550) )
       continue;

     //----------------------------------------
     // Select same theta bins
     if ( GetThetaBin(ltHA[4]) == GetThetaBin(ltHB[4]) &&
	  GetThetaBin(ltHA[4]) != -1 ){ 

       thBin = GetThetaBin(ltHA[4]);
       
       nDuplicatesA  = 0;
       nDuplicatesB  = 0;
       nDuplicatesAB = 0;

       // cout << " before assign phiA and phiB " << endl;

       // assign phiA and phiB
       for (Int_t i = 0 ; i < nCrystals ; i++){
	 
	 // phi determined by outer crystals
	 if( i == 4 ) continue;
	 
	 // array B phi 
	 if ( ( ltHB[i]         < ThMax[thBin])  &&
	      ( ltHB[i]         > ThMin[thBin])  &&
	      ( (EB[i] + EB[4]) > 450.        )  &&
	      ( (EB[i] + EB[4]) < 550.        )  
	      ) {
	   phiB = CrystalToPhiB(i);
	   indexB = i;
	   nDuplicatesB++;
	 }

	 // array A phi
	 if ( ( ltHA[i]         < ThMax[thBin])  &&
	      ( ltHA[i]         > ThMin[thBin])  &&
	      ( (EA[i] + EA[4]) > 450.        )  &&
	      ( (EA[i] + EA[4]) < 550.        )
	      ) {
	   phiA = CrystalToPhiA(i);
	   indexA = i;
	   nDuplicatesA++;
	 }

       } // end of: for (Int_t i = 0 ; i < nCryst
       
       nDuplicatesAB = nDuplicatesA + nDuplicatesB;

     } // end of : if (GetThetaBin(ltHA
     
     if( nDuplicatesAB > 2 ){
	cout << " nDuplicatesA  = " << nDuplicatesA  << endl;
	cout << " nDuplicatesB  = " << nDuplicatesB  << endl;
	cout << " nDuplicatesAB = " << nDuplicatesAB << endl;
	continue;
      }

     if ( (phiA  == -1) || 
	  (phiB  == -1) ||
	  (thBin == -1) )
       continue;
     
     if( nDuplicatesAB != 2 ){
       cout << " nDuplicatesA  = " << nDuplicatesA  << endl;
       cout << " nDuplicatesB  = " << nDuplicatesB  << endl;
       cout << " nDuplicatesAB = " << nDuplicatesAB << endl;
       continue;
     }
     
     //--------------------------------------
     //--------------------------------------
     // Lablike events have been selected
     
     //------------------------
     // Delta phi and phi res
     
     //Lablike delta phi in [0,360) range
     dPhiLab   = GetDeltaPhi(phiA,phiB);
     
     //Exact delta phi in [0,360) range
     dPhiXact  = GetDeltaPhi(simPhiA[0],simPhiB[0]);
     
     
     // dPhis all shifted to 0 for resolution plots
     dPhiXactRes = GetDeltaPhiRes(dPhiXact,dPhiLab);
     
     // first hit distance from x axis
     rA = CalculateR(YposA[0],
		     ZposA[0]);
     
     // beam angle wrt x axis
     betaA = CalculateBeta(XposA[0],
			   YposA[0],
			   ZposA[0]);
     
     //------------------------------------
     //------------------------------------
     //------------------------------------
     // Iterate the counting variables
     // used to form lablike asymmetries
     //------------------------------------
     
     // Plot Lab Delta Phi depending variables
     if (dPhiLab == 0) {
       hDPhiRes00[thBin]->Fill(dPhiXactRes);	
       hThRes00[thBin]->Fill(simtHA[0]);
       hBeta000[thBin]->Fill(betaA);

     }
     if (dPhiLab == 90){
       hDPhiRes90[thBin]->Fill(dPhiXactRes);
       hThRes90[thBin]->Fill(simtHA[0]);
       hBeta090[thBin]->Fill(betaA);

     }
     if (dPhiLab == 180 ){
       hBeta180[thBin]->Fill(betaA);

     }
     if (dPhiLab == 270){
       hBeta090[thBin]->Fill(betaA);
     }
     
     //--------------------------------------
     // Theta angle study
     
     hThExaSubThEI->Fill(etHA[4],simtHA[0]-etHA[4]);
     hThExaSubThEI->Fill(etHB[4],simtHB[0]-etHB[4]);

     hThExaSubThEO->Fill(etHA[4],simtHA[0]-etHA[indexA]);
     hThExaSubThEO->Fill(etHB[4],simtHB[0]-etHB[indexB]);

     if(dPhiLab!=180){
       hThExaSubThLI->Fill(ltHA[4],simtHA[0]-ltHA[4]);
       hThExaSubThLI->Fill(ltHB[4],simtHB[0]-ltHB[4]);
     }
     
     hThESubThLI->Fill(etHA[4],ltHA[4] - etHA[4]);

     hMaxThdiff->Fill(etHA[4],mintHAErr[4] - etHA[4]);
     hMaxThdiff->Fill(etHA[indexA],
		      mintHAErr[indexA] - etHA[indexA]);
     hMaxThdiff->Fill(etHB[4],mintHBErr[4] - etHB[4]);
     hMaxThdiff->Fill(etHB[indexB],
		      mintHBErr[indexB] - etHB[indexB]);
     hMaxThdiff->Fill(etHA[4],maxtHAErr[4] - etHA[4]);
     hMaxThdiff->Fill(etHA[indexA],
		      maxtHAErr[indexA] - etHA[indexA]);
     hMaxThdiff->Fill(etHB[4],maxtHBErr[4] - etHB[4]);
     hMaxThdiff->Fill(etHB[indexB],
		      maxtHBErr[indexB] - etHB[indexB]);
     
     //--------------------------------------
     // First hit position / 
     // incident angle study
     
     hXY[thBin]->Fill(Abs(XposA[0]),YposA[0]);
     hYZ[thBin]->Fill(YposA[0],ZposA[0]);
     
     hRVsDPhi->Fill(dPhiXact,rA);
     
     hBetaVsDPhi->Fill(dPhiXact,betaA);
     hBeta[thBin]->Fill(betaA);
     
     //exact angles delta phi
     //this will be acceptance corrected 
     //for two file (pol/unpol) case 
     hDPhi[thBin][fileNum]->Fill(dPhiXact);

     //first hit positions distributions 
     //depend on delta phi
     hYZ_dPhi[(Int_t)(dPhiLab/90.)]->Fill(YposA[0],ZposA[0]);
     

     //------------
     // Notes on phi
     //------------
     // When the beam is not in a fixed (say x) 
     // direction the lablike phi distribution 
     // can be in a different plane to the
     // real phi and therefore should
     // be projected before investigating 
     // biases
      
     // raw data delta phi is biased
     // for isotropic beam due to thresholds
     // on central crystals - clearly noticeable
     // when first hits are in outer crystals
     
     simPhiARes = GetPhiRes(simPhiA[0],phiA);
     
     //-------------------------
     // delta phi selections
     
     // phi minus phiLab (ie plotted around 0)
     // grouped in delta phi
     if     ( dPhiLab == 0 ){
       hPhiRes000[thBin]->Fill(simPhiARes);
       hYZ_000[thBin]->Fill(YposA[0],ZposA[0]);
     }
     else if( dPhiLab == 90 ){
       hPhiRes090[thBin]->Fill(simPhiARes);
       hYZ_090[thBin]->Fill(YposA[0],ZposA[0]);
     }
     else if     ( dPhiLab == 180 ){
       hPhiRes180[thBin]->Fill(simPhiARes);
       hYZ_180[thBin]->Fill(YposA[0],ZposA[0]);
     }
     else if     ( dPhiLab == 270 ){
       hPhiRes270[thBin]->Fill(simPhiARes);
       hYZ_270[thBin]->Fill(YposA[0],ZposA[0]);
     }
     // phi angle plots - shift again so
     // that no peaks are at boundary
     hPhiA->Fill(ShiftPhi(simPhiA[0],-135.0));
     hPhiB->Fill(ShiftPhi(simPhiB[0],-135.0));
     
     // phi minus phiLab (ie plotted around 0)
     // grouped in phiA (by crystals)
     if      (phiA == 0)
       hPhi000Res[thBin]->Fill(simPhiARes);
     else if(phiA == 90)
       hPhi090Res[thBin]->Fill(simPhiARes);
     else if(phiA == 180)
       hPhi180Res[thBin]->Fill(simPhiARes);
     else if(phiA == -90)
       hPhi270Res[thBin]->Fill(simPhiARes);

     if      (phiB==0)
       hPhi000Res[thBin]->Fill(simPhiB[0]);
     else if(phiB==90)
       hPhi090Res[thBin]->Fill(simPhiB[0]);
     else if(phiB==180)
       hPhi180Res[thBin]->Fill(simPhiB[0]);
     else if(phiB==-90)
       hPhi270Res[thBin]->Fill(simPhiB[0]);
     
     // Lablike events Selected
     //--------------------------------------
     
     //--------------------------------------
     // Extra selection criteria
     // hit in central crystal
     if ( !CentralYZ(YposA[0]) ||
	  !CentralYZ(ZposA[0]) ||
	  !CentralYZ(YposB[0]) || 
	  !CentralYZ(ZposB[0]) )
       continue;
     
     //------------------------------------
     // Iterate the counting variables
     // used to form lablike asymmetries
     // with extra conditions
     
     if (dPhiLab == 0.){
       hDPhiRes00_2[thBin]->Fill(dPhiXactRes);	
       hThRes00_2[thBin]->Fill(simtHA[0]);	
       hBeta000_2[thBin]->Fill(betaA);
     }
     if ( dPhiLab == 90 ){
       hDPhiRes90_2[thBin]->Fill(dPhiXactRes);
       hThRes90_2[thBin]->Fill(simtHA[0]);
       hBeta090_2[thBin]->Fill(betaA);	
     }
     if ( dPhiLab == 180 ){
       hBeta180_2[thBin]->Fill(betaA);	
     }
     if (dPhiLab == 270){
     }
     
     
     // evident for an isotropic source
     // is that for delta phi = 180 there is a
     // theta difference band that deviates far from zero
     // due to the wrong ordering of hits
     hThExaSubThEI_2->Fill(etHA[4],simtHA[0]-etHA[4]);
     hThExaSubThEI_2->Fill(etHB[4],simtHB[0]-etHB[4]);
     
     hThExaSubThEO_2->Fill(etHA[4],simtHA[0]-etHA[indexA]);
     hThExaSubThEO_2->Fill(etHB[4],simtHB[0]-etHB[indexB]);
     
     hThExaSubThLI_2->Fill(ltHA[4],simtHA[0]-ltHA[4]);
     hThExaSubThLI_2->Fill(ltHB[4],simtHB[0]-ltHB[4]);
     
     hDPhi_2[thBin][fileNum]->Fill(dPhiXact);     
     hYZ_dPhi_2[(Int_t)(dPhiLab/90.)]->Fill(YposA[0],ZposA[0]);
     hBeta_2[thBin]->Fill(betaA);
     hXY_2[thBin]->Fill(Abs(XposA[0]),YposA[0]);
     hYZ_2[thBin]->Fill(YposA[0],ZposA[0]);
     
     
     
     // phi minus phiLab (ie plotted around 0)
     // grouped in delta phi
     if     ( dPhiLab == 0 )
       hPhiRes000_2[thBin]->Fill(simPhiARes);
     else if( dPhiLab == 90 )
       hPhiRes090_2[thBin]->Fill(simPhiARes);
     if     ( dPhiLab == 180 )
       hPhiRes180_2[thBin]->Fill(simPhiARes);
     if     ( dPhiLab == 270 )
       hPhiRes270_2[thBin]->Fill(simPhiARes);
     
     if      (phiA==0)
       hPhi000Res_2[thBin]->Fill(simPhiA[0]);
     else if(phiA==90)
       hPhi090Res_2[thBin]->Fill(simPhiA[0]);
     else if(phiA==180)
       hPhi180Res_2[thBin]->Fill(simPhiA[0]);
     else if(phiA==270)
       hPhi270Res_2[thBin]->Fill(simPhiA[0]);
     
     
   } // end of : for (Int_t ientry = 0 ; ien.....
   
   //cout << " end " << endl;
   
   cout << endl;
   cout << " Lab Asymmetry"<<endl;
   cout << " theta \t" << "dPhi=0 \t" 
	<< "90 \t" << "180 \t" << "270" << endl;
   for (Int_t i = 0 ; i < 8 ; i++){
     cout << " " << plotTheta[i]     << "\t"
	  << " " << AsymMatrix[i][0] << "\t" 
	  << " " << AsymMatrix[i][1] << "\t"
	  << " " << AsymMatrix[i][2] << "\t" 
	  << " " << AsymMatrix[i][3] << endl;
   }
   
   cout << endl;
   cout << " True Asymmetry"<<endl;
   cout << " theta \t" << "dPhi=0 \t" 
	<< "90 \t" << "180 \t" << "270" << endl;
   for (Int_t i = 0 ; i < 8 ; i++){
     cout << " " << plotTheta[i]   << "\t "
	  << " " << AsymTrue[i][0] << "\t "
	  << " " << AsymTrue[i][1] << "\t " 
	  << " " << AsymTrue[i][2] << "\t " 
	  << " " << AsymTrue[i][3] << endl;
   }
   cout << endl;
   
   TCanvas *canvas = new TCanvas("canvas","canvas",
				 1500,1000);
   
   Bool_t plotAll = kFALSE;
   
   hThExaSubThEI->Draw("colz");
   hThExaSubThEI->GetXaxis()->SetTitle("#theta_{exact} (deg)");
   hThExaSubThEI->GetYaxis()->SetTitle("#theta_{exact} - #theta_{E} (deg)");

   plotName = "../Plots/hThExaSubThEI_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);

   hThExaSubThEI_2->Draw("colz");
   hThExaSubThEI_2->GetXaxis()->SetTitle("#theta_{exact} (deg)");
   hThExaSubThEI_2->GetYaxis()->SetTitle("#theta_{exact} - #theta_{lab} (deg)");

   plotName = "../Plots/hThExaSubThEI_2_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   hThExaSubThEO->Draw("colz");
   hThExaSubThEO->GetXaxis()->SetTitle("#theta_{exact} (deg)");
   hThExaSubThEO->GetYaxis()->SetTitle("#theta_{exact} - #theta_{E} (deg)");

   plotName = "../Plots/hThExaSubThEO_" + inputFileNumber;
   plotName += ".pdf";

   if(plotAll)
     canvas->SaveAs(plotName);

   hThExaSubThEO_2->Draw("colz");
   hThExaSubThEO_2->GetXaxis()->SetTitle("#theta_{exact} (deg)");
   hThExaSubThEO_2->GetYaxis()->SetTitle("#theta_{exact} - #theta_{E} (deg)");
   
   plotName = "../Plots/hThExaSubThEO_2_" + inputFileNumber;
   plotName += ".pdf";

   if(plotAll)
     canvas->SaveAs(plotName);
   
   hThExaSubThLI->Draw("colz");
   hThExaSubThLI->GetXaxis()->SetTitle("#theta_{exact} (deg)");
   hThExaSubThLI->GetYaxis()->SetTitle("#theta_{exact} - #theta_{lab} (deg)");

   plotName = "../Plots/hThExaSubThLI_" + inputFileNumber;
   plotName += ".pdf";

   canvas->SaveAs(plotName);

   
   hThExaSubThLI_2->Draw("colz");
   hThExaSubThLI_2->GetXaxis()->SetTitle("#theta_{exact} (deg)");
   hThExaSubThLI_2->GetYaxis()->SetTitle("#theta_{exact} - #theta_{E} (deg)");

   plotName = "../Plots/hThExaSubThLI_2_" + inputFileNumber;
   plotName += ".pdf";

   canvas->SaveAs(plotName);

   hThESubThLI->Draw("colz");
   hThESubThLI->GetXaxis()->SetTitle("energy #theta (deg)");
   hThESubThLI->GetYaxis()->SetTitle("lab #theta - energy #theta (deg)");
   
   plotName = "../Plots/hThESubThLI_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   hMaxThdiff->Draw("colz");
   hMaxThdiff->GetXaxis()->SetTitle("#theta_{E} (deg)");
   hMaxThdiff->GetYaxis()->SetTitle("max/min #theta_{Lab} - #theta_{E} (deg)");
   
   plotName = "../Plots/hMaxThdiff_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   hBetaVsDPhi->Draw("colz");
   hBetaVsDPhi->GetXaxis()->SetTitle("#Delta #phi (deg)");
   hBetaVsDPhi->GetYaxis()->SetTitle("#beta (deg)");

   plotName = "../Plots/hBetaVsDPhi_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   hRVsDPhi->Draw("colz");
   hRVsDPhi->GetXaxis()->SetTitle("#Delta #phi (deg)");
   hRVsDPhi->GetYaxis()->SetTitle("(y^{2} + z^{2})^{1/2}");
   
   plotName = "../Plots/hRVsDPhi_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   canvas->Clear();

   Int_t nRows = 1; 
   Int_t nColumns = 2;
   
   canvas->Divide(nColumns,nRows);

   canvas->cd(1);
   
   hPhiA->GetYaxis()->SetTitle("Counts");
   hPhiA->GetXaxis()->SetTitle("#phi_{A}");
   hPhiA->SetLineColor(kBlue);    
   
   hPhiA->Draw();
   
   canvas->cd(2);
   
   hPhiB->GetYaxis()->SetTitle("Counts");
   hPhiB->GetXaxis()->SetTitle("#phi_{B}");
   hPhiB->SetLineColor(kRed);
   
   hPhiB->Draw();

   plotName = "../Plots/hPhiAB_" + inputFileNumber;
   plotName += ".pdf";
   
   canvas->SaveAs(plotName);
   
   canvas->Clear();

   nRows = 2; 
   nColumns = nThbins/nRows;
   
   canvas->Divide(nColumns,nRows);
   
   for( Int_t th = 0 ; th < nThbins ; th++){

     hDPhi[th][fileNum]->SetLineColor(kBlue);
     hDPhi[th][fileNum]->SetMinimum(0.0);

     hDPhi_2[th][fileNum]->SetLineColor(kBlue);
     hDPhi_2[th][fileNum]->SetMinimum(0.0);

     if(fileNum==1){

       hDPhi[th][1]->Divide(hDPhi[th][0]);
       hDPhi_2[th][1]->Divide(hDPhi_2[th][0]);
       
     }
     
     hThRes00[th]->SetLineColor(kBlue);
     hDPhiRes00[th]->SetLineColor(kBlue);
     
     hThRes90[th]->SetLineColor(kRed);
     hDPhiRes90[th]->SetLineColor(kRed);
     
     hThRes00_2[th]->SetLineColor(kBlue);
     hDPhiRes00_2[th]->SetLineColor(kBlue);
     
     hThRes90_2[th]->SetLineColor(kRed);
     hDPhiRes90_2[th]->SetLineColor(kRed);

     hThRes00[th]->GetXaxis()->SetTitle("#theta_{A} (deg)");
     hDPhiRes00[th]->GetXaxis()->SetTitle("#Delta#phi (deg)");

     hThRes90[th]->GetXaxis()->SetTitle("#theta_{A} (deg)");
     hDPhiRes90[th]->GetXaxis()->SetTitle("#Delta#phi (deg)");

     hThRes00_2[th]->GetXaxis()->SetTitle("#theta_{A} (deg)");
     hDPhiRes00_2[th]->GetXaxis()->SetTitle("#Delta#phi (deg)");

     hThRes90_2[th]->GetXaxis()->SetTitle("#theta_{A} (deg)");
     hDPhiRes90_2[th]->GetXaxis()->SetTitle("#Delta#phi (deg)");

     hThRes00[th]->GetYaxis()->SetTitle("Counts");
     hDPhiRes00[th]->GetYaxis()->SetTitle("Counts");

     hThRes90[th]->GetYaxis()->SetTitle("Counts");
     hDPhiRes90[th]->GetYaxis()->SetTitle("Counts");

     hThRes00_2[th]->GetYaxis()->SetTitle("Counts");
     hDPhiRes00_2[th]->GetYaxis()->SetTitle("Counts");

     hThRes90_2[th]->GetYaxis()->SetTitle("Counts");
     hDPhiRes90_2[th]->GetYaxis()->SetTitle("Counts");

     hPhiRes000[th]->GetYaxis()->SetTitle("Counts");
     hPhiRes000[th]->GetXaxis()->SetTitle("#phi");
     hPhiRes000[th]->SetLineColor(kBlue);    

     hPhiRes000_2[th]->GetYaxis()->SetTitle("Counts");
     hPhiRes000_2[th]->GetXaxis()->SetTitle("#phi");
     hPhiRes000_2[th]->SetLineColor(kBlue);    

     hPhiRes090[th]->GetYaxis()->SetTitle("Counts");
     hPhiRes090[th]->GetXaxis()->SetTitle("#phi");
     hPhiRes090[th]->SetLineColor(kRed);    

     hPhiRes090_2[th]->GetYaxis()->SetTitle("Counts");
     hPhiRes090_2[th]->GetXaxis()->SetTitle("#phi");
     hPhiRes090_2[th]->SetLineColor(kRed);    

     hPhiRes180[th]->GetYaxis()->SetTitle("Counts");
     hPhiRes180[th]->GetXaxis()->SetTitle("#phi");
     hPhiRes180[th]->SetLineColor(kGreen);    

     hPhiRes180_2[th]->GetYaxis()->SetTitle("Counts");
     hPhiRes180_2[th]->GetXaxis()->SetTitle("#phi");
     hPhiRes180_2[th]->SetLineColor(kGreen);    

     hPhiRes270[th]->GetYaxis()->SetTitle("Counts");
     hPhiRes270[th]->GetXaxis()->SetTitle("#phi");
     hPhiRes270[th]->SetLineColor(kMagenta);    

     hPhiRes270_2[th]->GetYaxis()->SetTitle("Counts");
     hPhiRes270_2[th]->GetXaxis()->SetTitle("#phi");
     hPhiRes270_2[th]->SetLineColor(kMagenta);    

     hPhi000Res[th]->GetYaxis()->SetTitle("Counts");
     hPhi000Res[th]->GetXaxis()->SetTitle("#phi");
     hPhi000Res[th]->SetLineColor(kBlue);    
     
     hPhi090Res[th]->GetYaxis()->SetTitle("Counts");
     hPhi090Res[th]->GetXaxis()->SetTitle("#phi");
     hPhi090Res[th]->SetLineColor(kRed);    

     hPhi180Res[th]->GetYaxis()->SetTitle("Counts");
     hPhi180Res[th]->GetXaxis()->SetTitle("#phi");
     hPhi180Res[th]->SetLineColor(kGreen);    
     
     hPhi270Res[th]->GetYaxis()->SetTitle("Counts");
     hPhi270Res[th]->GetXaxis()->SetTitle("#phi");
     hPhi270Res[th]->SetLineColor(kMagenta);    
     
     hPhi000Res_2[th]->GetYaxis()->SetTitle("Counts");
     hPhi000Res_2[th]->GetXaxis()->SetTitle("#phi");
     hPhi000Res_2[th]->SetLineColor(kBlue);    

     hPhi090Res_2[th]->GetYaxis()->SetTitle("Counts");
     hPhi090Res_2[th]->GetXaxis()->SetTitle("#phi");
     hPhi090Res_2[th]->SetLineColor(kRed);    
     
     hPhi180Res_2[th]->GetYaxis()->SetTitle("Counts");
     hPhi180Res_2[th]->GetXaxis()->SetTitle("#phi");
     hPhi180Res_2[th]->SetLineColor(kGreen);    

     hPhi270Res_2[th]->GetYaxis()->SetTitle("Counts");
     hPhi270Res_2[th]->GetXaxis()->SetTitle("#phi");
     hPhi270Res_2[th]->SetLineColor(kMagenta);    
     
     hBeta000[th]->SetLineColor(kBlue);
     hBeta090[th]->SetLineColor(kRed);
     hBeta180[th]->SetLineColor(kGreen);
     
     hBeta000_2[th]->SetLineColor(kBlue);
     hBeta090_2[th]->SetLineColor(kRed);
     hBeta180_2[th]->SetLineColor(kGreen);
     
     hBeta090[th]->Scale(1./2);
     hBeta090_2[thBin]->Scale(1./2);
     
     hBeta_2[th]->SetLineColor(kGreen);
     
   }

   canvas->Clear();
   canvas->Divide(nColumns,nRows);
   
   for( Int_t th = 0 ; th < nThbins ; th++){ 
     canvas->cd(th+1);
     
     cout << " th+1 " << (th+1) << endl;
     
     hPhi000Res[th]->Draw("");
     hPhi090Res[th]->Draw("same");
          
     hPhi180Res[th]->Draw("same");
     hPhi270Res[th]->Draw("same");
   }
   plotName = "../Plots/hPhiXXXRes_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
      
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
     hThRes90_2[th]->Draw("");
     hThRes00_2[th]->Draw("same");  
   }
   plotName = "../Plots/hThRes_2_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hDPhiRes90[th]->Draw("HIST");
     hDPhiRes00[th]->Draw("same");  
   }
   plotName = "../Plots/hDPhiRes_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);

   for( Int_t th = 0 ; th < nThbins ; th++){ 
     canvas->cd(th+1);
     hPhiRes090[th]->Draw("");
     hPhiRes270[th]->Draw("same");
     hPhiRes180[th]->Draw("same");
     hPhiRes000[th]->Draw("same");
   }
   plotName = "../Plots/hPhiResXXX_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){ 
     canvas->cd(th+1);
     
     hPhiRes090_2[th]->Draw("");
     hPhiRes270_2[th]->Draw("same");
     hPhiRes000_2[th]->Draw("same");
     hPhiRes180_2[th]->Draw("same");
     
   }
   plotName = "../Plots/hPhiResXXX_2_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){ 
     canvas->cd(th+1);
     hPhi000Res_2[th]->Draw("");
     hPhi090Res_2[th]->Draw("same");     
     hPhi180Res_2[th]->Draw("same");
     hPhi270Res_2[th]->Draw("same");
     
   }
   plotName = "../Plots/hPhiXXXRes_2_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hDPhiRes90_2[th]->Draw("HIST");
     hDPhiRes00_2[th]->Draw("same");  
   }
   plotName = "../Plots/hDPhiRes_2_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hXY[th]->Draw("colz");
   }
   plotName = "../Plots/hXY_" + inputFileNumber;
   plotName += ".pdf";

   if(plotAll)
     canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hXY_2[th]->Draw("colz");
   }
   plotName = "../Plots/hXY_2_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hYZ[th]->Draw("colz");
   }
   plotName = "../Plots/hYZ_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hYZ_2[th]->Draw("colz");
   }
   plotName = "../Plots/hYZ_2_" + inputFileNumber;
   plotName += ".pdf";

   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hYZ_000[th]->Draw("colz");
   }
   plotName = "../Plots/hYZ_000_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);

   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hYZ_090[th]->Draw("colz");
   }
   plotName = "../Plots/hYZ_090_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hYZ_180[th]->Draw("colz");
   }
   plotName = "../Plots/hYZ_180_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
 
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hYZ_270[th]->Draw("colz");
   }

   plotName = "../Plots/hYZ_270_" + inputFileNumber;
   plotName += ".pdf";
   
   canvas->SaveAs(plotName);
   
   for( Int_t p = 0 ; p < 4 ; p++){
     canvas->cd(p+1);
     hYZ_dPhi[p]->Draw("colz");
   }
   plotName = "../Plots/hYZ_dPhi_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   for( Int_t p = 0 ; p < 4 ; p++){
     canvas->cd(p+5);
     hYZ_dPhi_2[p]->Draw("colz");
   }
   plotName = "../Plots/hYZ_dPhi_2_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hBeta[th]->Draw("");
   }
   plotName = "../Plots/hBeta_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hBeta_2[th]->Draw("same");
   }
   plotName = "../Plots/hBeta_2_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     
     hBeta180[th]->Draw();
     hBeta090[th]->Draw("same HIST");
     hBeta000[th]->Draw("same");
   }
   plotName = "../Plots/hBetaX_" + inputFileNumber;
   plotName += ".pdf";

   if(plotAll)
     canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     
     hBeta090_2[th]->Draw("HIST");
     hBeta180_2[th]->Draw("same");
     hBeta000_2[th]->Draw("same");
   }
   plotName = "../Plots/hBetaX_2_" + inputFileNumber;
   plotName += ".pdf";
   
   if(plotAll)
     canvas->SaveAs(plotName);
   
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hDPhi[th][fileNum]->Draw("");
   }
   plotName = "../Plots/hDPhi_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
   
   for( Int_t th = 0 ; th < nThbins ; th++){
     canvas->cd(th+1);
     hDPhi_2[th][fileNum]->Draw("");
   }
   plotName = "../Plots/hDPhi_2_" + inputFileNumber;
   plotName += ".pdf";
   canvas->SaveAs(plotName);
   
   return 0;
   
} // End of: InvestigateAcceptance()


Int_t TSim::CalculateAsymmetryLab(TString inputFileNumber){

   GetThetaBinValues();
   
   cout << endl;
   cout << " Getting asymmetry " << endl;
   
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

   sortDataTree->SetBranchAddress("EAX",EAX);
   sortDataTree->SetBranchAddress("EBX",EBX);

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
   
   for(Int_t j = 0 ; j <nThbins; j++){
     for(Int_t k = 0 ; k < 4; k++){
       AsymMatrix[j][k] = 0;}
   }

   for(Int_t j = 0 ; j < nThbins; j++){
     for(Int_t k = 0 ; k < 4; k++){
       AsymTrue[j][k] = 0;}
   }

   Long64_t nEntries = sortDataTree->GetEntries();

   Int_t   thBin  = -1;
   Float_t phiA   = -1;
   Float_t phiB   = -1;
   Int_t   indexA = -1;
   Int_t   indexB = -1;
   Float_t totEA  = 0;
   Float_t totEB  = 0;

   Float_t totEAX  = 0;
   Float_t totEBX  = 0;

   Float_t dPhiLab = -999.;
   
   Float_t dPhiXact    = -999.;
  
   Int_t nDuplicatesAB = 0;
   Int_t nDuplicatesB  = 0;
   Int_t nDuplicatesA  = 0;
   
   // Event loop
   for (Int_t ientry = 0 ; ientry < nEntries ; ientry++){
     
     cout << " start " << endl;
     
     sortDataTree->GetEvent(ientry);
     
     totEA  = 0;
     totEB  = 0;
     totEAX  = 0;
     totEBX  = 0;

     for (Int_t k = 0; k < nCrystals; k++){
       totEA  += EA[k];
       totEB  += EB[k];
       totEAX += EAX[k];
       totEBX += EBX[k];
     }

     
     //----------------------------------------
     // only interested in events with
     // full energy deposited
     if( (totEA < 450) || 
	 (totEA > 550) ||
	 (totEB < 450) ||
	 (totEB > 550) )
       continue;

     thBin  = -1;
     
     //----------------------------------------
     // Select theta bin 
     // it must be the same for A and B
     if ( GetThetaBin(ltHA[4]) == GetThetaBin(ltHB[4]) &&
	  GetThetaBin(ltHA[4]) != -1 ){ 

       thBin = GetThetaBin(ltHA[4]);

       phiA   = -1;
       phiB   = -1;
       indexA = -1;
       indexB = -1;
       nDuplicatesA  = 0;
       nDuplicatesB  = 0;
       nDuplicatesAB = 0;

       // assign phiA and phiB
       for (Int_t i = 0 ; i < nCrystals ; i++){
	 
	 // phi determined by outer crystals
	 if( i == 4 ) continue;
	 
	 // array B phi 
	 if ( ( ltHB[i]         < ThMax[thBin])  &&
	      ( ltHB[i]         > ThMin[thBin])  &&
	      ( (EB[i] + EB[4]) > 450.        )  &&
	      ( (EB[i] + EB[4]) < 550.        )  
	      ) {
	   phiB = CrystalToPhiB(i);
	   indexB = i;
	   nDuplicatesB++;
	 }
	 
	 // array A phi
	 if ( ( ltHA[i]         < ThMax[thBin])  &&
	      ( ltHA[i]         > ThMin[thBin])  &&
	      ( (EA[i] + EA[4]) > 450.        )  &&
	      ( (EA[i] + EA[4]) < 550.        )
	      ) {
	   phiA = CrystalToPhiA(i);
	   indexA = i;
	   nDuplicatesA++;
	 }

       } // end of: for (Int_t i = 0 ; i < nCryst
       
       nDuplicatesAB = nDuplicatesA + nDuplicatesB;

     } // end of : if (GetThetaBin(ltHA
     
     if( nDuplicatesAB > 2 ){
       cout << " nDuplicatesA  = " << nDuplicatesA  << endl;
       cout << " nDuplicatesB  = " << nDuplicatesB  << endl;
       cout << " nDuplicatesAB = " << nDuplicatesAB << endl;
       continue;
     }

     if ( (phiA  == -1) || 
	  (phiB  == -1) ||
	  (thBin == -1) )
       continue;
     
     if( nDuplicatesAB != 2 ){
       cout << " nDuplicatesA  = " << nDuplicatesA  << endl;
       cout << " nDuplicatesB  = " << nDuplicatesB  << endl;
       cout << " nDuplicatesAB = " << nDuplicatesAB << endl;
       continue;
     }
     
     //--------------------------------------
     //--------------------------------------
     // Lablike events have been selected
     
     //------------------------
     // Delta phi and phi res
     
     //Lablike delta phi in [0,360) range
     dPhiLab   = GetDeltaPhi(phiA,phiB);
     
     //Exact delta phi in [0,360) range
     dPhiXact  = GetDeltaPhi(simPhiA[0],simPhiB[0]);
     
     
     //------------------------------------
     //------------------------------------
     //------------------------------------
     // Iterate the counting variables
     // used to form lablike asymmetries
     //------------------------------------

     if      (dPhiLab == 0) 
       AsymMatrix[thBin][0]+= 1;
     else if (dPhiLab == 90)
       AsymMatrix[thBin][1] += 1;
     else if (dPhiLab == 180 )
       AsymMatrix[thBin][2] += 1;
     else if (dPhiLab == 270)
       AsymMatrix[thBin][3] += 1;
     
     // Lablike events Selected
     //--------------------------------------
     
     //--------------------------------------
     // Extra selection criteria
     // hit in central crystal
     if ( !CentralYZ(YposA[0]) ||
	  !CentralYZ(ZposA[0]) ||
	  !CentralYZ(YposB[0]) || 
	  !CentralYZ(ZposB[0]) )
       continue;
     
     //------------------------------------
     // Iterate the counting variables
     // used to form lablike asymmetries
     // with extra conditions
     
     if      (dPhiLab == 0.)
       AsymTrue[thBin][0] += 1;
     else if ( dPhiLab == 90 )
       AsymTrue[thBin][1] += 1;
     else if ( dPhiLab == 180 )
       AsymTrue[thBin][2] += 1;
     else if (dPhiLab == 270)
       AsymTrue[thBin][3] += 1;
     
     
   } // end of : for (Int_t ientry = 0 ; ien.....

   
   cout << endl;
   cout << " Lab Asymmetry"<<endl;
   cout << " theta \t" << "dPhi=0 \t" 
	<< "90 \t" << "180 \t" << "270" << endl;
   for (Int_t i = 0 ; i < 8 ; i++){
     cout << " " << plotTheta[i]     << "\t"
	  << " " << AsymMatrix[i][0] << "\t" 
	  << " " << AsymMatrix[i][1] << "\t"
	  << " " << AsymMatrix[i][2] << "\t" 
	  << " " << AsymMatrix[i][3] << endl;
   }
   
   cout << endl;
   cout << " True Asymmetry"<<endl;
   cout << " theta \t" << "dPhi=0 \t" 
	<< "90 \t" << "180 \t" << "270" << endl;
   for (Int_t i = 0 ; i < 8 ; i++){
     cout << " " << plotTheta[i]   << "\t "
	  << " " << AsymTrue[i][0] << "\t "
	  << " " << AsymTrue[i][1] << "\t " 
	  << " " << AsymTrue[i][2] << "\t " 
	  << " " << AsymTrue[i][3] << endl;
   }
   cout << endl;
   
   return 0;

} //end of: TSim::CalculateAsymmetryLa

void TSim::GraphAsymmetryLab(TString inputFileNumber1,
			     TString inputFileNumber2) {
  
  cout << endl;
  cout << "--------------------- " << endl;
  cout << "  GraphAsymmetryLab   " << endl;
  
  this->SetStyle();
  
  //  GetThetaBinValues() is implemented in
  //  CalculateAsymmetryLab
  
  // entangled/polarised first
  // To do: change to have this second
  // with the unpolarised calculated 
  // first as is done in TLab.C
  
  // Iterate the counting variables
  CalculateAsymmetryLab(inputFileNumber1);
  
  Int_t nFiles = 2;

  if(inputFileNumber2=="??")
    nFiles = 1;
  
  cout << endl;
  cout << " There are " << nFiles << " Files " << endl;
  
  if(nFiles==2)
    cout << " Dividing by unpolarised data " << endl;

  TString plotName;
  plotName = "../Plots/Asym_" + inputFileNumber1;
  plotName = plotName + ".pdf";

  // The ratio to be calculated for the
  // lab data:  90 e.g corresponds to 
  // A(90) = N(90)/N(0) 
  Int_t   dPhiDiff = 90;
  
  Float_t thetaLowEdge  = 10.;
  Float_t thetaHighEdge = 170.;
  
  // Consider removing
  CalculateABC_Lab();
  CalculateABC_True();
   
  TGraphErrors *grB_Lab;
  TGraphErrors *grC_Lab;
  TGraphErrors *grB_True;
  TGraphErrors *grC_True;
  
  grB_Lab     = new TGraphErrors(nThbins,plotTheta,pB_Lab,0,0);
  grC_Lab     = new TGraphErrors(nThbins,plotTheta,pC_Lab,0,0);
  grB_True    = new TGraphErrors(nThbins,plotTheta,pB_True,0,0);
  grC_True    = new TGraphErrors(nThbins,plotTheta,pC_True,0,0);
  
  for (Int_t th = 0 ; th < nThbins ; th++){
    AsPhiDiff[th] = 0.;
    AePhiDiff[th] = 0.;
    AsTrue[th] = 0.;
    AeTrue[th] = 0.;
  }
  
  // temporary variables
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

  // use [ N(90) + N(270) ]/2 for N(90)
  Bool_t use270 = kFALSE;

  for (Int_t file = 0 ; file < nFiles ; file++){
    
    if (file==1)
      CalculateAsymmetryLab(inputFileNumber2);
    
    for (Int_t i = 0 ; i < nThbins ; i++){
      if (AsymMatrix[i][0] != 0){
	
	mu[i] = (AsymMatrix[i][1] - AsymMatrix[i][0]);
	mu[i] = mu[i]/(AsymMatrix[i][1] + AsymMatrix[i][0]);
	
	if (dPhiDiff  == 90){
	  
	  if(!use270){
	    AsPhiDiff[i] = AsymMatrix[i][1]/AsymMatrix[i][0];
	    AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix[i][1])+
					     (1/AsymMatrix[i][0]));
	  }
	  else{
	    AsPhiDiff[i] = (AsymMatrix[i][1]+
			    AsymMatrix[i][3])/(2*AsymMatrix[i][0]);
	    AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/(AsymMatrix[i][1]+
						 AsymMatrix[i][3]))+
					     (1/AsymMatrix[i][0]));
	  }
	}
	if (dPhiDiff  == 180){
	  AsPhiDiff[i] = AsymMatrix[i][2]/AsymMatrix[i][0];
	  AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix[i][2])+
					   (1/AsymMatrix[i][0]));
	}
	if (dPhiDiff  == 270){
	  AsPhiDiff[i] = AsymMatrix[i][3]/AsymMatrix[i][0];
	  AePhiDiff[i] = AsPhiDiff[i]*Sqrt((1/AsymMatrix[i][3])+
					   (1/AsymMatrix[i][0]));
	}
      }
    
      // duplicate for later use in two file case
      if     (file==0){
	AsPhiDiff1[i] = AsPhiDiff[i];
	AePhiDiff1[i] = AePhiDiff[i];
	
      }
      else if(file==1){
	AsPhiDiffR[i] = AsPhiDiff1[i]/AsPhiDiff[i] ;
	AePhiDiffR[i] = AsPhiDiffR[i] * 
	  Sqrt( AePhiDiff[i]*AePhiDiff[i]/(AsPhiDiff[i]*AsPhiDiff[i]) +
		AePhiDiff1[i]*AePhiDiff1[i]/(AsPhiDiff1[i]*AsPhiDiff1[i]));
      }
             
    }
    
    for (Int_t i = 0 ; i < nThbins ; i++){
      if (AsymTrue[i][0] != 0){
	if (dPhiDiff  == 90){
	  
	  if(!use270){
	    AsTrue[i] = AsymTrue[i][1]/AsymTrue[i][0];
	    AeTrue[i] = AsTrue[i]*Sqrt((1/AsymTrue[i][1])+
				       (1/AsymTrue[i][0]));
	  }
	  else{
	    AsTrue[i] = (AsymTrue[i][1]+AsymTrue[i][3])/(2*AsymTrue[i][0]);
	    AeTrue[i] = AsTrue[i]*Sqrt((1/(AsymTrue[i][1]+
					   AsymTrue[i][3]))+
				       (1/AsymTrue[i][0]));
	  }
	}
	if (dPhiDiff  == 180){
	  AsTrue[i] = AsymTrue[i][2]/AsymTrue[i][0];
	  AeTrue[i] = AsTrue[i]*Sqrt((1/AsymTrue[i][2])+
				     (1/AsymTrue[i][0]));
	}
	if (dPhiDiff  == 270){
	  AsTrue[i] = AsymTrue[i][3]/AsymTrue[i][0];
	  AeTrue[i] = AsTrue[i]*Sqrt((1/AsymTrue[i][3])+
				     (1/AsymTrue[i][0]));
	}
      }
      
      if     (file==0){
	AsTrue1[i] = AsTrue[i];
	AeTrue1[i] = AeTrue[i];
	
      }
      else if(file==1){
	AsTrueR[i] = AsTrue1[i]/AsTrue[i] ;
	AeTrueR[i] = AsTrueR[i] * 
	  Sqrt( AeTrue[i]*AeTrue[i]/(AsTrue[i]*AsTrue[i]) + 
		AeTrue1[i]*AeTrue1[i]/(AsTrue1[i]*AsTrue1[i]));
      }
      
    } // end of: for (Int_t i = 0 ; i < nThbins ...
  
  } // end of: for (Int_t file = 0 ; 
   
  //use rho1 (detectors only finite in theta)
  Float_t aTheoryX[nThbins];
  
  //theory curve with half angles

  Float_t aTheory1[nThbins];
  Float_t aTheory[nThbins];
  Float_t aTheory2[nThbins];
  Float_t aTheoryE[nThbins];

  // half resolution in dPhi
  // 35.0 is result from Chloe Schoolings fits
  // sigma
  Float_t alpha1  = 1.0*DegToRad()*35.0;
  Float_t alpha   = 1.5*DegToRad()*35.0;
  Float_t alpha2  = 2.0*DegToRad()*35.0;


  // half resolution in theta
  Float_t semiSpan = DegToRad()*(ThMax[0] - ThMin[0])/2.;
    
  TTheory *theory = new TTheory();

  // calculating theory curves
  cout << endl;
  cout << " Calculating theory curve/s ... " << endl;
  cout << endl;
  cout << " semiSpan = " << semiSpan*RadToDeg() << endl;
  cout << " alpha1   = " << alpha1*RadToDeg()   << endl;
  cout << " alpha    = " << alpha*RadToDeg()    << endl;
  cout << " alpha2   = " << alpha2*RadToDeg()   << endl;
  cout << endl;
  
  for (Int_t i = 0; i<nThbins; i++){
    plotTheta[i] = plotTheta[i]*DegToRad();
    
    aTheoryX[i] = theory->rho1(plotTheta[i],semiSpan);

    aTheory[i]  = theory->rho2(plotTheta[i],semiSpan,alpha);
    aTheory1[i] = theory->rho2(plotTheta[i],semiSpan,alpha1);
    aTheory2[i] = theory->rho2(plotTheta[i],semiSpan,alpha2);
    
    aTheory[i] = 0.7*aTheory1[i] + 0.3*aTheory2[i];

    aTheoryE[i] = 0.;Abs(aTheory[i]-aTheory1[i]); 
    
    if( dPhiDiff == 180 ){
      aTheoryX[i] = 1.0;
      aTheory1[i] = 1.0;
      aTheory[i]  = 1.0;
      aTheory2[i] = 1.0;
    }
    
    // acceptance
    f_AsPhiDiff[i] = AsPhiDiff[i]/aTheoryX[i];
    // error on acceptance
    f_AePhiDiff[i] = AePhiDiff[i]/aTheoryX[i];
    
    plotTheta[i] = plotTheta[i]*RadToDeg();
  }
  
  
  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);
  
  TH1F *hr;
  
  Float_t maxY = 2.8;
  Float_t minY = 0.8;

  
  maxY = 1.4;
  minY = 0.9;
  
  if(dPhiDiff==180){
    maxY = 6.0;
    if(nFiles==2){
	maxY = 2.0;
      }
  }
  
  hr = canvas->DrawFrame(10,minY,170,maxY);
  hr->GetXaxis()->SetTitle("#theta (deg)");

  TGraphErrors *grAsym[3];
  
  TGraphErrors * grMu = new TGraphErrors(nThbins,plotTheta,mu,0,muE);

  if     (nFiles==1)
    grAsym[0] = new TGraphErrors(nThbins,plotTheta,AsPhiDiff,0,AePhiDiff);
  else if(nFiles == 2)
    grAsym[0] = new TGraphErrors(nThbins,plotTheta,AsPhiDiffR,0,AePhiDiffR);
  
  grAsym[1] = new TGraphErrors(nThbins,plotTheta,aTheory,0,aTheoryE);

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
  grAsym[1]->SetFillColor(kRed);
  grAsym[1]->SetFillStyle(3003);
  grAsym[2]->SetLineColor(kGreen);
  grAsym[2]->SetMarkerColor(kGreen);

  TLegend *leg = new TLegend(0.6,0.75,0.9,0.85);
  
  TString theoryLegendTitle = " ";
  alpha1 = alpha1*RadToDeg();
  theoryLegendTitle.Form("theory #alpha_{#phi} = %.1f^{o}", alpha1);
  
  Char_t  yAxis[128];
  
  if(nFiles==1)
    sprintf(yAxis,"N(%d^{o})/N(0^{o})",dPhiDiff);
  else if(nFiles == 2)
    sprintf(yAxis,"N^{E}(%d^{o})/N^{E}(0^{o}) * N^{U}(0^{o})/N^{U}(%d^{o}) ",dPhiDiff,dPhiDiff);
  
  hr->GetYaxis()->SetTitle(yAxis);
  leg->AddEntry(grAsym[1],theoryLegendTitle,"P L");
  leg->AddEntry(grAsym[0],"simulated lab", "E P");
  leg->AddEntry(grAsym[2],"true simulated lab","E P");
  grAsym[0]->Draw("P E");

  // plot theory curve
  grAsym[1]->Draw("3 P L E SAME");
  grAsym[2]->Draw("P E SAME");
  leg->Draw();

  if     (nFiles==1){
    plotName = "../Plots/A_" + inputFileNumber1;
  }
  else if(nFiles==2){
    plotName  =  "../Plots/A_" + inputFileNumber1;
    plotName += "_";
    plotName += inputFileNumber2;
  }  
  
  plotName.Form(plotName + "_%d",dPhiDiff);
  plotName += ".pdf";
  
  canvas->SaveAs(plotName);
  
  //------------------------------------------------
  // Fourier Coefficients Graph

  Bool_t plotFCs = kFALSE;
     
  leg->Clear();
  
  maxY =  1.0;
  minY = -1.2;
    
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
    
  if(plotFCs) 
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
    
  if(plotFCs)
    canvas->SaveAs(plotName);
  
  //mu plot
  hr = canvas->DrawFrame(thetaLowEdge,-0.2,
			  thetaHighEdge,0.5);
  
  hr->GetXaxis()->SetTitle("#theta (deg)");

  sprintf(yAxis,"(N(%d^{o})-N(0^{o})/(N(%d^{o})+N(0^{o})",
	  dPhiDiff,dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);
  
  plotName = "../Plots/Mu_";
  
  plotName += inputFileNumber1;

  plotName += "_%d.pdf";
  
  plotName.Form(plotName,dPhiDiff);

  grMu->Draw("P E");
  
  if(plotFCs)
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
      
      cout << " pABC["<< i << "][" << j << "] = " 
	   <<  pABC[i][j] << endl;
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
    

    for(Int_t j = 0 ; j < 3 ; j++){
      fABC_True[i][j] = fABC_True[i][j]/4.;
      pABC_True[i][j] = fABC_True[i][j]/fABC_True[i][bin000];
      
      cout << endl;
      cout << " pABC_True["<< i << "][" << j << "] = " 
	   <<  pABC_True[i][j] << endl;
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

    Bool_t comments = kFALSE;
    if(comments){    
      cout << endl;
      for(Int_t j = 0 ; j < 3 ; j++){
	//fABC_Lab[i][j] = 1.0*fABC_Lab[i][j]/4.;
	pABC_Lab[i][j] = 1.0*fABC_Lab[i][j]/fABC_Lab[i][0];
	
	cout << " pABC_Lab["<< i << "][" << j << "] = " 
	     <<  pABC_Lab[i][j] << endl;
      }
      
      for(Int_t j = 0 ; j < 4 ; j++){
	cout << " AsymMatrix["<< i << "][" << j << "] = " 
	     <<  AsymMatrix[i][j] << endl;
      }
    }  
  
  } // end of: for(Int_t i = 0 ; i < nTh....
  
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
  
  for(Int_t j = 0 ; j <nThbins; j++)
    for(Int_t k = 0 ; k < nPhibinsSim; k++)
      AsymMatrix_sim[j][k] = 0;
  
  Float_t  halfBinSize = 180./nPhibinsSim; 

  Long64_t nEntries = simDataTree->GetEntries();
  for (Int_t ientry = 0 ; ientry < nEntries ; ientry++){
    simDataTree->GetEvent(ientry);
    
    // Ensure first hits are in central crystals
    // We know there was an energy deposit there
    if( !CentralYZ(YposA_1st) ||
	!CentralYZ(ZposA_1st) ||
	!CentralYZ(YposB_1st) ||
	!CentralYZ(ZposB_1st) )
      continue;
    
    Float_t dPhi_1st = PhiA_1st + PhiB_1st;
    if(dPhi_1st < 0)
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
    
   //printing out the asym matrix 
  
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
    
    // Ensure first hits are in central crystals
    // We know there was an energy deposit there
    if( !CentralYZ(YposA_1st) ||
	!CentralYZ(ZposA_1st) ||
	!CentralYZ(YposB_1st) ||
	!CentralYZ(ZposB_1st) )
      continue;
    

    // scattering in array A
    if ( PhiA_2nd   < 499. && 
	 ThetaA_2nd < 499. ){
      
      scatterArray = 'A';
      
      thetaA   = ThetaA_2nd;
      thetaB   = ThetaB_1st;
      thetaABS = ThetaA_1st;

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
	thetaB < 1. )
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
  //entangled[1] = kTRUE;
  
  polarised[0] = kFALSE;
  polarised[1] = kFALSE;

  if(!entangled[0])
    polarised[0] = kTRUE;
  
  // if(!entangled[1])
  //   polarised[1] = kTRUE;
   

  // Calculating ratios for desired dPhiDiff
  Float_t AsPhiDiff[nThbins] = {0.};
  Float_t AePhiDiff[nThbins] = {0.};
  
  Float_t AsPhiDiff1[nThbins] = {0.};
  Float_t AePhiDiff1[nThbins] = {0.};

  Float_t AsPhiDiffR[nThbins] = {0.};
  Float_t AePhiDiffR[nThbins] = {0.};
  
  // Entangled, Polarised
  Float_t AsMx_simEP[nThbins][nPhibinsSim] = {{0.},{0.}};
  Float_t AsMx_simEPInt[nThbins] = {0.}; 
  
  // Entangled, Unpolarised
  Float_t AsMx_simEU[nThbins][nPhibinsSim] = {{0.},{0.}};
  Float_t AsMx_simEUInt[nThbins] = {0.}; 
  
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
  
  // nGraphs = 1 for single file analysis
  // nGraphs = 2 for two file analysis
  // nGraphs = 2 for single file scattered analysis
  Int_t   nGraphs = 2;
  if( nThetaSBins!=0 ){
    nGraphs = nThetaSBins + 1;
  }
  else if(inputFileNumber1 == inputFileNumber2 &&
	  nThetaSBins==0 ){
    nGraphs = 1;
  }
  
  cout << endl;
  cout <<  " There are " << nGraphs << " graphs " << endl;
    
  Bool_t subUnPol = kFALSE;

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
      
      if(subUnPol){
	//-------------------
	//!!!!!!!!
	// Unpolarised subtraction method 
	// NB this method is incorrect but may 
	// be a useful framework for implementing 
	// method of dividing by unpolarised data
	
	// sum the contribrutions
	for (Int_t p = 0 ; p < nPhibinsSim ; p++ ){
	  if     (g == 0){
	    AsMx_simEP[i][p]  = AsymMatrix_sim[i][p];
	    AsMx_simEPInt[i] += AsMx_simEP[i][p];
	  }
	  else if(g == 1){
	    AsMx_simEU[i][p]  = AsymMatrix_sim[i][p];
	    AsMx_simEUInt[i] += AsMx_simEU[i][p];
	  }
	}
        
	if     (g == 0){
	  cout << endl;
	  cout << " AsMx_simEPInt[" << i 
	       << "] = "
	       << AsMx_simEPInt[i] 
	       << endl;
	}
	else if(g == 1){
	  cout << endl;
	  cout << " AsMx_simEUInt[" << i 
	       << "] = "
	       << AsMx_simEUInt[i] 
	       << endl;
	}
	
	cout << endl;
	// normalise
	for (Int_t p = 0 ; p < nPhibinsSim ; p++ ){
	  if     (g == 0 && 
		  AsMx_simEPInt[i] != 0 ){
	    AsMx_simEP[i][p]  = AsMx_simEP[i][p]/AsMx_simEPInt[i];
	    
	    cout << " AsMx_simEP[" << i 
		 << "][" << p << "] = "
		 << AsMx_simEP[i][p] 
		 << endl;
	    
	    
	  }
	  else if(g == 1 && 
		  AsMx_simEUInt[i] != 0 ){
	    AsMx_simEU[i][p]  = AsMx_simEU[i][p]/AsMx_simEUInt[i];
	    
	    AsymMatrix_sim[i][p] = AsMx_simEP[i][p] - AsMx_simEU[i][p] + 1./nPhibinsSim;
	    
	    cout << " AsymMatrix_sim[" << i 
		 << "][" << p << "] = "
		 << AsymMatrix_sim[i][p] 
		 << endl;
	  }
	}
      
      } //end of: if(subUnPol){
   
      // cout << endl;
//       cout << " Calculating Ratios " << endl;
      
      if (dPhiDiff  == 90){
	
	AsPhiDiff[i] = AsymMatrix_sim[i][bin090];
	AsPhiDiff[i] = AsPhiDiff[i]/(AsymMatrix_sim[i][0]);
	
// 	cout << endl;
// 	cout << " AsPhiDiff[" << i 
// 	     << "] = " << AsPhiDiff[i] << endl;
	
	//AsPhiDiff[i] = AsPhiDiff[i]+AsymMatrix_sim[i][bin270];
	//AsPhiDiff[i] = AsPhiDiff[i]/(2.*AsymMatrix_sim[i][0]);
	
	AePhiDiff[i] = AsymMatrix_sim[i][bin090];
	//AePhiDiff[i] = AsymMatrix_sim[i][bin270];
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
    
    grAsymR   = new TGraphErrors(nThbins,plotTheta,AsPhiDiffR,0,AePhiDiffR);
    
    grA[g]    = new TGraphErrors(nThbins,plotTheta,pA,0,0);
    grB[g]    = new TGraphErrors(nThbins,plotTheta,pB,0,0);
    grC[g]    = new TGraphErrors(nThbins,plotTheta,pC,0,0);
    
  } //end of : for (Int_t g = 0 ; g < nGr
  
  cout << endl;
  cout << " graphs have been made" << endl;
  
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
  
  TString plotName = "";
  
  Char_t yAxis[128];
  
  alpha1 = RadToDeg()*alpha1;

  sprintf(theoryLegendTitle,
	  "theory curve #alpha_{#phi} = %.1f^{o}",
	  alpha1);  
  
  alpha1 = DegToRad()*alpha1;

  Float_t maxY = 1.8;
  Float_t minY = 0.8;

  if( entangled[0] || 
      entangled[1] ){
    minY =   0.0;
    maxY =   3.5;
  }
  
  TH1F *hr;
  hr = canvas->DrawFrame(10,minY,170,maxY);
  hr->GetXaxis()->SetTitle("#theta (deg)");
  sprintf(yAxis,"P(%d^{o})/P(0^{o})",dPhiDiff);
  hr->GetYaxis()->SetTitle(yAxis);
  hr->GetYaxis()->SetTitleOffset(0.7);
  
  
  TLegend *leg = new TLegend(0.6,0.75,0.9,0.85);
  
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
  
  plotName = "../Plots/" ;
  plotName += inputFileNumber1;
  plotName += "_";
  
  if(nThetaSBins==0){
    plotName += inputFileNumber2;
    plotName += "_A_%d_";  
    plotName += "_%d";
    plotName += "_%d";
    plotName.Form(plotName,nPhibinsSim);
  }
  else{
    plotName += "_A_%d_";  
    plotName += "Scattered_%d";
    plotName.Form(plotName,dPhiDiff,nPhibinsSim);

  }
  plotName += "dPhiBins.pdf";
  canvas->SaveAs(plotName);
  
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
  
  if(nThetaSBins==0){
    plotName =  "../Plots/FourierCoeffs_";
    plotName += inputFileNumber1;
    plotName += inputFileNumber2;
    plotName += "_%dPhiBins.pdf";
    plotName.Form(plotName,nPhibinsSim);
  }
  else{
    plotName = "../Plots/FourierCoeffs_Scattered_";
    plotName += inputFileNumber1;
    plotName += "_%dPhiBins.pdf";
    plotName.Form(plotName,nPhibinsSim);
  }

  leg->Draw();
  
  canvas->SaveAs(plotName);
  
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
  garyStyle->SetPadRightMargin(0.05);  //percentage
  garyStyle->SetPadLeftMargin(0.12);    //percentage
  garyStyle->SetPadBottomMargin(0.12); //percentage
  
  //----------- Histogram
  
  //Histos
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
  //garyStyle->SetTitleOffset(1.6,"Y");
  
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
