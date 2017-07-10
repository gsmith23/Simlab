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
    
  rootFile = "../Data/hitsLab" + fileNumber ;;
  
  rootFile = rootFile + ".pet.root";;

  Initialise();
  
  //!! temporary 
  // SetAsymmetry(fileNumber);

}

TSim::~TSim(){
}

//----------------------------------------------

/** Public member functions *********/

void TSim::Loop()
{
  
   Long64_t nentries = theTree->GetEntriesFast();

   for (Long64_t i = 0 ; i < nentries ; i++) {
     theTree->GetEntry(i);
           
   }
}

Float_t TSim::GetAverageEnergy(){
  
  Float_t averageEnergy = 0.;
  
  for( Int_t i = 0 ; i < theTree->GetEntries() ; i++ ){
    theTree->GetEvent(i);
  
    averageEnergy += energy1P1;
  
  }
  return averageEnergy/theTree->GetEntries();
}

void TSim::PlotTF1(){
  
  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);
  
  TH1F * hr;
  hr = canvas->DrawFrame(0.0,0.5,180.0,3.0);
    
  TF1 * fA90dA0 = new TF1("fA90dA0",this,0.,180.,1,"TSim"); 
  
  // Dummy
  fA90dA0->SetParameter(0,1);
  
  fA90dA0->GetYaxis()->SetTitle("#rho");
  fA90dA0->GetXaxis()->SetTitle("#theta (deg.)");
  fA90dA0->Draw();

  canvas->SaveAs("../Plots/fA90dA0.pdf");
  
}

Float_t TSim::ElectronEnergyToTheta(Float_t nrg){
  return RadToDeg()*ACos(2. - 511./(511. - nrg));
}

Float_t TSim::PhotonEnergyToTheta(Float_t nrg){
  return RadToDeg()*ACos(2. - 511./nrg);;
}

Float_t TSim::ThetaToPhotonEnergy(Float_t theta){
  return (511./(2 - Cos(TMath::DegToRad()*theta)));
}

Float_t TSim::ThetaToElectronEnergy(Float_t theta){
  return (511. - (511./(2. - Cos(TMath::DegToRad()*theta))));
}

    
Int_t TSim::GetDPhiBin(Float_t dP, Float_t rangeDPhi){
  
  Int_t bin = -1;
  
  if     ( dP > (0.0  - rangeDPhi) && dP < ( 0.0 + rangeDPhi) ){ 
    bin = 0;
  }
  else if( dP > (90.0 - rangeDPhi) && dP < (90.0 + rangeDPhi) ){
    bin = 90;
  }
  else if( dP > (180. - rangeDPhi) || dP < (-180. + rangeDPhi) ){
    bin = 180;
  }
  else if( dP > (-90. - rangeDPhi) && dP < (-90.0 + rangeDPhi) ){
    bin = 270;
  }
  
  return bin;
}


Int_t TSim::GetThetaBin(Float_t theta,
			      Int_t nbins){
  
  Int_t bin = -1;
  
  if(nbins==9){
    if     (theta >=   0.0 && theta <  20.){
      bin = 0;
    }
    else if(theta >=  20.0 && theta <  40.){
      bin = 1;
    }
    else if(theta >=  40.0 && theta <  60.){
      bin = 2;
    }
    else if(theta >=  60.0 && theta <  80.){
      bin = 3;
    }
    else if(theta >=  80.0 && theta < 100.){
      bin = 4;
    }
    else if(theta >= 100.0 && theta < 120.){
      bin = 5;
    }
    else if(theta >= 120.0 && theta < 140.){
      bin = 6;
    }
    else if(theta >= 140.0 && theta < 160.){
      bin = 7;
    }
    else if(theta >= 160.0 && theta < 180.){
      bin = 8;
    }
    else{
      bin = 9;
    }
  }
  else if(nbins==4){
    
    if     ( theta >= 0.0 && theta < 45 ){
      bin = 0;
    }
    else if( theta >= 45.0 && theta < 90. ){
      bin = 1;
    }
    else if( theta >= 90.0 && theta < 135. ){
      bin = 2;
    }
    else if( theta >= 135.0 && theta < 180. ){
      bin = 3;
    }
    
  }
  
  return bin; 
}




void TSim::Initialise(){
  
  cout << endl;
  cout << " Connecting Branches " << endl;
  
  theFile = new TFile(rootFile);
  theTree = (TTree*)theFile->Get("phiTree");
  
  theTree->SetBranchAddress("nHits", &nHits, &b_nHits);
  theTree->SetBranchAddress("nParticles", &nParticles, &b_nParticles);
  theTree->SetBranchAddress("eventNumber", eventNumber, &b_eventNumber);
  theTree->SetBranchAddress("detector", detector, &b_detector);
  theTree->SetBranchAddress("detectorUnitID", detectorUnitID, &b_detectorUnitID);
  theTree->SetBranchAddress("energy", energy, &b_energy);
  theTree->SetBranchAddress("x", x, &b_x);
  theTree->SetBranchAddress("y", y, &b_y);
  theTree->SetBranchAddress("z", z, &b_z);
  theTree->SetBranchAddress("nOriginalTracks", nOriginalTracks, &b_nOriginalTracks);
  theTree->SetBranchAddress("originalParticleNumber", originalParticleNumber, &b_originalParticleNumber);
  theTree->SetBranchAddress("nTracks", nTracks, &b_nTracks);
  theTree->SetBranchAddress("nHitsC", &nHitsC, &b_nHitsC);
  theTree->SetBranchAddress("nHitsW", &nHitsW, &b_nHitsW);
  theTree->SetBranchAddress("eventType", &eventType, &b_eventType);
  theTree->SetBranchAddress("iC", iC, &b_iC);
  theTree->SetBranchAddress("iW", &iW, &b_iW);
  theTree->SetBranchAddress("i1P1", &i1P1, &b_i1P1);
  theTree->SetBranchAddress("i1P2", &i1P2, &b_i1P2);
  theTree->SetBranchAddress("sinogramR", &sinogramR, &b_sinogramR);
  theTree->SetBranchAddress("sinogramPhi", &sinogramPhi, &b_sinogramPhi);
  theTree->SetBranchAddress("meanHitAngle", &meanHitAngle, &b_meanHitAngle);
  theTree->SetBranchAddress("energy1P1", &energy1P1, &b_energy1P1);
  theTree->SetBranchAddress("energy1P2", &energy1P2, &b_energy1P2);
  theTree->SetBranchAddress("i2P1", &i2P1, &b_i2P1);
  theTree->SetBranchAddress("i2P2", &i2P2, &b_i2P2);
  theTree->SetBranchAddress("energy2P1", &energy2P1, &b_energy2P1);
  theTree->SetBranchAddress("energy2P2", &energy2P2, &b_energy2P2);
  theTree->SetBranchAddress("thetaP1E", &thetaP1E, &b_thetaP1E);
  theTree->SetBranchAddress("thetaP2E", &thetaP2E, &b_thetaP2E);
  theTree->SetBranchAddress("thetaP1V", &thetaP1V, &b_thetaP1V);
  theTree->SetBranchAddress("thetaP2V", &thetaP2V, &b_thetaP2V);
  theTree->SetBranchAddress("thetaP1VTrue", &thetaP1VTrue, &b_thetaP1VTrue);
  theTree->SetBranchAddress("thetaP2VTrue", &thetaP2VTrue, &b_thetaP2VTrue);
  theTree->SetBranchAddress("phiP1", &phiP1, &b_phiP1);
  theTree->SetBranchAddress("phiP2", &phiP2, &b_phiP2);
  theTree->SetBranchAddress("dPhi", &dPhi, &b_dPhi);
  
}

Int_t TSim::SortEvents(TString inputFileName,
		       Int_t   comments){

  //=========================
  // Initialise variables &
  // connect tree and leaves
  
  inputFileName = "../Data/hitsLab" + inputFileName;
  
  TString outputFileName;
    
  outputFileName = inputFileName + ".pet.root";
  inputFileName  = inputFileName + ".root";
  
  TFile* inputFile = new TFile(inputFileName);
  TFile* outputFile = new TFile(outputFileName,"recreate");

  cout << endl;
  cout << " SortEvents: " << endl;
  cout << endl;
  cout << " as you requested, I will now sort the  " << endl;
  cout << " events in the root file for you and create" << endl;
  cout << " a new root file containing additional variables " << endl;
  cout << " to use in the PET analysis. " << endl;
  cout << endl;
  cout << " I will disgard the unwanted events " << endl;
  cout << " ie those not with 4 hits or " << endl;
  cout << " an unsuitable LOR  " << endl;
  
  cout << endl;
  cout << " Input File : " << inputFileName  << endl;
  cout << " Output File: " << outputFileName << endl;
  
  //==========================
  // read in output from 
  // Hits2Tree.C

  TTree* tree=(TTree*)inputFile->Get("tree");

  //============================================
  //--------------------------------------------
  // Variable declarations
  
  //-------------------------------
  // Initial array initialisations
  for (Int_t i = 0; i < 32; i++){
    
    eventNumber[i] = 0;
    detectorUnitID[i] = 0;
    
    detector[i] = '?';
    
    energy[i] = -999.0;

    x[i] = -999.9;
    y[i] = -999.9;
    z[i] = -999.9;

    nOriginalTracks[i] = 0;
    originalParticleNumber[i] = 0;
  
    nTracks[i] = 0;
    
  }

  tree->SetBranchAddress("nHits"         ,&nHits);
  tree->SetBranchAddress("nParticles"    ,&nParticles);
  
  tree->SetBranchAddress("eventNumber"   ,eventNumber);

  tree->SetBranchAddress("detector"      ,detector);
  tree->SetBranchAddress("detectorUnitID",detectorUnitID);
  tree->SetBranchAddress("energy"        ,energy);
  
  tree->SetBranchAddress("x"             ,x);
  tree->SetBranchAddress("y"             ,y);
  tree->SetBranchAddress("z"             ,z);

  tree->SetBranchAddress("nOriginalTracks",nOriginalTracks);
  tree->SetBranchAddress("originalParticleNumber",originalParticleNumber);
  tree->SetBranchAddress("nTracks",nTracks);


  //============================================
  // New variables
  Long64_t nEvents = 0;

  Float_t         newX[32];
  Float_t         newY[32];
  Float_t         newZ[32];
  
  for (Int_t i = 0; i < 32; i++){
    iC[i] = -1;
    iW[i] = -1;

    newX[i] = 0.;
    newY[i] = 0.;
    newZ[i] = 0.;
  }
  
  TVector3 crossVector(0.,0.,0.);
  TVector3 perpVector(0.,0.,0.);

  TVector3 xVector(1.,0.,0.);
  TVector3 yVector(0.,1.,0.);
  TVector3 zVector(0.,0.,1.);
  TVector3 mixedVector(1.,1.,1.);

  Float_t  rotationArray[3][3];
  
  Float_t hitAngleP1;
  Float_t hitAngleP2;

  // Type 0 - neither of below (1,3,5+ hits)
  // Type 1 - candidate PET event        (2 hits)
  // Type 2 - candidate scattered event  (4 hits)
  
  TTree *outputTree = new TTree("phiTree","Output from SortEvents.C");

  outputTree->Branch("nHits",&nHits,"nHits/I");
  outputTree->Branch("nParticles",&nParticles,"nParticles/I");
  
  outputTree->Branch("eventNumber",&eventNumber,"eventNumber[nHits]/L");
  outputTree->Branch("detector",&detector,"detector[nHits]/C");
  outputTree->Branch("detectorUnitID",&detectorUnitID,"detectorUnitID[nHits]/I");
  outputTree->Branch("energy",&energy,"energy[nHits]/F");
  
  outputTree->Branch("x",&x,"x[nHits]/F");
  outputTree->Branch("y",&y,"y[nHits]/F");
  outputTree->Branch("z",&z,"z[nHits]/F");
  
  outputTree->Branch("nOriginalTracks",&nOriginalTracks,
	       "nOriginalTracks[nHits]/I");
  outputTree->Branch("originalParticleNumber",&originalParticleNumber,
	       "originalParticleNumber[nHits]/I");
  outputTree->Branch("nTracks",&nTracks,
	       "nTracks[nHits]/I");

  // New Variables
  
  // used for all PET events
  outputTree->Branch("nHitsC",&nHitsC,"nHitsC/I");
  outputTree->Branch("nHitsW",&nHitsW,"nHitsW/I");

  outputTree->Branch("eventType",&eventType,"eventType/I");

  outputTree->Branch("iC",&iC,"iC[nHitsC]/F");
  outputTree->Branch("iW",&iW,"iW[nHitsW]/F");

  outputTree->Branch("i1P1",&i1P1,"i1P1/I");
  outputTree->Branch("i1P2",&i1P2,"i1P2/I");

  outputTree->Branch("p3P1Lab","TVector3",&p3P1Lab);
  outputTree->Branch("p3P2Lab","TVector3",&p3P2Lab);

  outputTree->Branch("LOR","TVector3",&LOR);
  outputTree->Branch("LORXY","TVector3",&LORXY);

  outputTree->Branch("v3P1Lab","TVector3",&v3P1Lab);
  outputTree->Branch("v3P2Lab","TVector3",&v3P2Lab);

  outputTree->Branch("midpoint","TVector3",&midpoint);
  outputTree->Branch("midpointXY","TVector2",&midpointXY);

  outputTree->Branch("sinogramR",&sinogramR,"sinogramR/F");
  outputTree->Branch("sinogramPhi",&sinogramPhi,"sinogramPhi/F");

  outputTree->Branch("meanHitAngle",&meanHitAngle,"meanHitAngle/F");
  
  outputTree->Branch("energy1P1",&energy1P1,"energy1P1/F");
  outputTree->Branch("energy1P2",&energy1P2,"energy1P2/F");

  //used for scattered PET events only

  outputTree->Branch("i2P1",&i2P1,"i2P1/I");
  outputTree->Branch("i2P2",&i2P2,"i2P2/I");

  outputTree->Branch("energy2P1",&energy2P1,"energy2P1/F");
  outputTree->Branch("energy2P2",&energy2P2,"energy2P2/F");
  
  outputTree->Branch("thetaP1E",&thetaP1E,"thetaP1E/F");
  outputTree->Branch("thetaP2E",&thetaP2E,"thetaP2E/F");

  outputTree->Branch("v3ScP1Lab","TVector3",&v3ScP1Lab);
  outputTree->Branch("v3ScP2Lab","TVector3",&v3ScP2Lab);

  outputTree->Branch("v3P1Event","TVector3",&v3P1Event);
  outputTree->Branch("v3P2Event","TVector3",&v3P2Event);

  outputTree->Branch("thetaP1V",&thetaP1V,"thetaP1V/F");
  outputTree->Branch("thetaP2V",&thetaP2V,"thetaP2V/F");

  outputTree->Branch("thetaP1VTrue",&thetaP1VTrue,"thetaP1VTrue/F");
  outputTree->Branch("thetaP2VTrue",&thetaP2VTrue,"thetaP2VTrue/F");
  
  outputTree->Branch("phiP1",&phiP1,"phiP1/F");
  outputTree->Branch("phiP2",&phiP2,"phiP2/F");
  outputTree->Branch("dPhi",&dPhi,"dPhi/F");
  
  //---------------------
  //    Histograms

  TH1I* hEventTypes = new TH1I("hEventTypes","hEventTypes",3,-0.5,2.5);
  
  nEvents = tree->GetEntries();
  
  cout << endl;
  cout << " " << nEvents << " events in this tree " <<  endl;
    

  //============
  // EVENT LOOP
  
  for(Int_t i = 0 ; i < nEvents ; i++){ 
      
    //-------------------------------------------
    // Re-initialise variables
    
    nHits      = 0;
    nHitsC     = 0;
    nHitsW     = 0;
    
    nParticles = 0;
        
    i1P1 = -1;
    i2P1 = -1;
    i1P2 = -1;
    i2P2 = -1;
    
    for (Int_t j = 0; j < 32; j++){

      eventNumber[j] = 0;
      detectorUnitID[j] = 0;
    
      detector[j] = '?';
    
      energy[j] = -999.;
      
      x[j] = -999.;
      y[j] = -999.;
      z[j] = -999.;
    
      nOriginalTracks[j] = 0;
      originalParticleNumber[j] = 0;
  
      nTracks[j] = 0;
    
    } // end of: for (Int_t j = 0; j < 32; j++

    for (Int_t j = 0; j < 32; j++){
      iC[j] = -1;
      iW[j] = -1;
      
      newX[j] = -999.9;
      newY[j] = -999.9;
      newZ[j] = -999.9;
    }

    phiP1 = -999.9;
    phiP2 = -999.9;
    dPhi = -999.9;
    thetaP1E  = -999.9;
    thetaP2E  = -999.9;
    thetaP1V = -999.9;
    thetaP2V = -999.9;

    thetaP1VTrue = -999.9;
    thetaP2VTrue = -999.9;

    energy1P1 = -999.0;
    energy1P2 = -999.0;
    energy2P1 = -999.0;
    energy2P2 = -999.0;

    hitAngleP1      = -999.9;
    hitAngleP2      = -999.9;
    meanHitAngle    = -999.9;
    
	
    LOR.SetXYZ(-999.9,-999.9,-999.9);
    LORXY.SetXYZ(-999.9,-999.9,-999.9);
            
    crossVector.SetXYZ(-999.9,-999.9,-999.9);
    perpVector.SetXYZ(-999.9,-999.9,-999.9);
 
    v3P1Lab.SetXYZ(-999.9,-999.9,-999.9);
    v3P2Lab.SetXYZ(-999.9,-999.9,-999.9);
    
    p3P1Lab.SetXYZ(-999.9,-999.9,-999.9);
    p3P2Lab.SetXYZ(-999.9,-999.9,-999.9);

    v3ScP1Lab.SetXYZ(-999.9,-999.9,-999.9);
    v3ScP2Lab.SetXYZ(-999.9,-999.9,-999.9);
  
    v3P1Event.SetXYZ(-999.9,-999.9,-999.9);
    v3P2Event.SetXYZ(-999.9,-999.9,-999.9);

    for( Int_t j = 0 ; j < 3 ; j++)
      for (Int_t k = 0 ; k < 3 ; k++)
	rotationArray[j][k] = -999.9;
	     
    midpointXY.Set(-999.9,-999.9);
    midpoint.SetXYZ(-999.9,-999.9,-999.9);
    
    sinogramR   = -999.9;
    sinogramPhi = -999.9;
    
    eventType = 0;
    //------------------------------------------------------
    
    tree->GetEvent(i);
    
    if(comments > 1){
      cout << endl;
      cout << " event   " << i   << endl;
      cout << " nHits = " << nHits << endl;
    }

    for(Int_t j = 0; j < nHits; j++){ 
      
      // change to keV;
      energy[j] = energy[j]*1000.;
      
      if(detectorUnitID[j] == 1 ){
	iW[nHitsW] = j;
	nHitsW++;
      }
      else if(detectorUnitID[j] > 1){
	iC[nHitsC] = j;
	nHitsC++;
      }
      else{
	cout << endl;
	cout << " !!! Shouldny !!! " << endl;
	cout << " detectorUnitID["<< j <<"] = " << detectorUnitID[j] << endl;
	cout << " detector["<< j <<"] = " << detector[j] << endl;
      }
      
     
      // // only want one original track
      // if(nOriginalTracks[j]!=1){
	
      // 	if(comments > 1)
      // 	  cout << " more than one original track " << endl;
	
      // 	// this will force a continue in the next step
      // 	break;
      // }
            
    } //end of: for(Int_t j = 0; j < nHits; j++){


    if(comments > 1 ){
      
      cout << endl;
      cout << " i               = " << i << endl;
      cout << " eventNumber[0]  = " << eventNumber[0] << endl;
      cout << " eventNumber[1]  = " << eventNumber[1] << endl;
      cout << " nHits  = " << nHits << endl;
      cout << " nHitsC = " << nHitsC << endl;
      cout << " nHitsW = " << nHitsW << endl;
    }

    if( (nHitsC+nHitsW)!=nHits ){
      cout << endl;
      cout << " (nHitsC+nHitsW)!=nHits " << endl;
      cout << " i               = " << i << endl;
      cout << " eventNumber[0]  = " << eventNumber[0] << endl;
      cout << " eventNumber[1]  = " << eventNumber[1] << endl;
      cout << " nHits  = " << nHits << endl;
      cout << " nHitsC = " << nHitsC << endl;
      cout << " nHitsW = " << nHitsW << endl;
      
      return 3;
    }
    
    // total number of events in tree 
    // (with x-hits per event)

    // ------------------------------
    // CONDITION 
    // require 2 [4] hits per events
    // for PET [scattered] events 
    // ------------------------------
    
    if ( nHitsC!=4 ){
      if( comments > 1 ){
	cout << " there are neither three nor four calorimeter hits in this event " 
	     << endl;
      }
      
      hEventTypes->Fill(eventType);
      
      continue;
    }
    
    if      (nHitsC==4) eventType = 2;
    else                eventType = 0; // shouldn't ever happen
    
    if(comments == 1){
      cout << endl;
      cout << " eventType = " << eventType << endl;
    }
    
    // total number of events in tree with 
    // the desired number of hits
    
    // Indices for hits in event
    // particle one hit one (i1P1)
    // particle one hit two (i2P1)
    
    // particle two hit one (i1P2)
    // particle two hit two (i2P2)

    if     (eventType==1){
      i1P1 = iC[0];
      i1P2 = iC[1];
    }
    else if(eventType==2){
      i1P1 = iC[0];
      i2P1 = iC[1];
      i1P2 = iC[2];
      i2P2 = iC[3];
    }
    else if(eventType==3){
      i1P1 = iC[0];
      i2P1 = iC[1];
      i1P2 = iC[2];
    }

    
    // position vector for photon one
    p3P1Lab.SetX(x[i1P1]);   
    p3P1Lab.SetY(y[i1P1]);
    p3P1Lab.SetZ(z[i1P1]);	

    // position vector for photon 2
    p3P2Lab.SetX(x[i1P2]);   
    p3P2Lab.SetY(y[i1P2]);
    p3P2Lab.SetZ(z[i1P2]);	

    LOR = p3P2Lab - p3P1Lab;
    LORXY.SetXYZ( x[i1P2] - x[i1P1], y[i1P2] - y[i1P1], 0.);

        
    // if( nHits==4 ){
      
    //   cout << endl;
    //   for (Int_t j = 0 ; j < nHits ; j++)
    // 	cout << " iC[" << j << "] = " <<  iC[j] << endl;
      
    // }
    
    // cout << " nHitsC    = " << nHitsC    << endl;
    // cout << " LOR.Mag() = " << LOR.Mag() << endl;
    // cout << " x[i1P2]   = " << x[i1P2]   << endl;
    // cout << " i1P1      = " << i1P1      << endl;
    // cout << " i1P12     = " << i1P2      << endl; 
    
    // Events where LOR is 
    // very small likely have
    // first hit indices 
    // corresponding to scattering
    // of single photon
    if( !GoodLOR(LOR.Mag(),'l')){
      hEventTypes->Fill(0); 
      continue;
    }
    

    
    hEventTypes->Fill(eventType);
    
    midpointXY.Set( (x[i1P1]+x[i1P2])/2.,(y[i1P1]+y[i1P2])/2. );

    midpoint.SetXYZ( (x[i1P1]+x[i1P2])/2.,(y[i1P1]+y[i1P2])/2., (z[i1P1]+z[i1P2])/2.  );

    // hit vector for photon 1
    v3P1Lab.SetX(x[i1P1]-midpoint.X());   
    v3P1Lab.SetY(y[i1P1]-midpoint.Y());
    v3P1Lab.SetZ(z[i1P1]-midpoint.Z());	
      
    // hit vector for photon 2
    v3P2Lab.SetX(x[i1P2]-midpoint.X());   
    v3P2Lab.SetY(y[i1P2]-midpoint.Y());
    v3P2Lab.SetZ(z[i1P2]-midpoint.Z());
    
    sinogramR     = TMath::Sqrt( TMath::Power(midpointXY.X(),2) + TMath::Power(midpointXY.Y(),2)  );
    
    sinogramPhi   = LORXY.Angle(yVector)*TMath::RadToDeg();

    hitAngleP1    = 180. - LOR.Angle(p3P1Lab)*TMath::RadToDeg();
    hitAngleP2    = LOR.Angle(p3P2Lab)*TMath::RadToDeg();
    meanHitAngle  = (hitAngleP1 + hitAngleP2)/2;
    
    // cout << endl;
    // cout << " hitAngleP1   = " << hitAngleP1 << endl;
    // cout << " hitAngleP2   = " << hitAngleP2 << endl;
    // cout << " meanHitAngle = " << meanHitAngle << endl;

    if(comments > 1 ){
      cout << energy[i1P1] << " is energy[i1P1]" << endl;
      cout << energy[i1P2] << " is energy[i1P2]" << endl;
      
      cout << " hitAngleP1 = " << hitAngleP1 << endl;
      cout << " hitAngleP2 = " << hitAngleP2 << endl;
      
      if(eventType == 2){
	cout << energy[i2P1] << " is energy[i2P1]" << endl;
	cout << energy[i2P2] << " is energy[i2P2]" << endl;
      }
      
    }
	
    if(comments > 1){
      cout << endl;
      cout << " theta P1 = " << thetaP1E << " degrees " << endl;
      cout << " theta P2 = " << thetaP2E << " degrees " << endl;
    }
    
    energy1P1 = energy[i1P1];
    energy1P2 = energy[i1P2];
    
    //==============================
    // Variables for scattered events 
    if(eventType == 2 ){
      
      energy2P1 = energy[i2P1];
      energy2P2 = energy[i2P2];
      thetaP1E  = PhotonEnergyToTheta(energy[i2P1]);
      thetaP2E  = PhotonEnergyToTheta(energy[i2P2]);

      v3ScP1Lab.SetX(x[i2P1] - x[i1P1]);   
      v3ScP1Lab.SetY(y[i2P1] - y[i1P1]);
      v3ScP1Lab.SetZ(z[i2P1] - z[i1P1]);	

      v3ScP2Lab.SetX(x[i2P2] - x[i1P2]);   
      v3ScP2Lab.SetY(y[i2P2] - y[i1P2]);
      v3ScP2Lab.SetZ(z[i2P2] - z[i1P2]);

      thetaP1V = v3P1Lab.Angle(v3ScP1Lab)*TMath::RadToDeg();
      thetaP2V = v3P2Lab.Angle(v3ScP2Lab)*TMath::RadToDeg();
      
      thetaP1VTrue = p3P1Lab.Angle(v3ScP1Lab)*TMath::RadToDeg();
      thetaP2VTrue = p3P2Lab.Angle(v3ScP2Lab)*TMath::RadToDeg();

      //crossVector = (LOR.Unit()).Cross(mixedVector.Unit());
      
      crossVector = (LOR.Unit()).Cross(zVector.Unit());
      perpVector  = (LOR.Unit()).Cross(crossVector.Unit());

      //========================================================
      // Matrix Rotation of the lab frame into the event frame

      rotationArray[0][0] =  xVector.Unit().Dot(perpVector.Unit());
      rotationArray[1][0] =  yVector.Unit().Dot(perpVector.Unit());
      rotationArray[2][0] =  zVector.Unit().Dot(perpVector.Unit());
      rotationArray[0][1] =  xVector.Unit().Dot(crossVector.Unit());
      rotationArray[1][1] =  yVector.Unit().Dot(crossVector.Unit());
      rotationArray[2][1] =  zVector.Unit().Dot(crossVector.Unit());
      rotationArray[0][2] =  xVector.Unit().Dot(LOR.Unit());
      rotationArray[1][2] =  yVector.Unit().Dot(LOR.Unit());
      rotationArray[2][2] =  zVector.Unit().Dot(LOR.Unit());

      // Hit position in event co-ordinate system
      for(Int_t j = 0 ; j < nHits ; j++){
	newX[j] = (x[j]*rotationArray[0][0]) + (y[j]*rotationArray[1][0]) + (z[j]*rotationArray[2][0]);
	newY[j] = (x[j]*rotationArray[0][1]) + (y[j]*rotationArray[1][1]) + (z[j]*rotationArray[2][1]);
	newZ[j] = (x[j]*rotationArray[0][2]) + (y[j]*rotationArray[1][2]) + (z[j]*rotationArray[2][2]);
      }
    
      v3P1Event.SetX(newX[i2P1] - newX[i1P1]);   
      v3P1Event.SetY(newY[i2P1] - newY[i1P1]);
      v3P1Event.SetZ(newZ[i2P1] - newZ[i1P1]);	
      
      phiP1 = v3P1Event.Phi() * TMath::RadToDeg();
  
      v3P2Event.SetX(newX[i2P2] - newX[i1P2]);   
      v3P2Event.SetY(newY[i2P2] - newY[i1P2]);
      v3P2Event.SetZ(newZ[i2P2] - newZ[i1P2]);

      phiP2 =  v3P2Event.Phi()* TMath::RadToDeg();
  
      dPhi = v3P1Event.DeltaPhi(v3P2Event)*TMath::RadToDeg();
            
      if(comments > 1) 

	cout << " delta phi " << dPhi << endl; 
      
	
    } //end of: if(eventType == 2)       
  
  
  //=========================
  //=========================
  
    outputTree->Fill();
    
  } // end of: for(Int_t i = 0 ; i < nEvents     
    
  //=========================
  //=========================
  
  cout << endl;  
  cout << " " << hEventTypes->GetBinContent(3) << 
    " entries with 4 hits " << endl;
      
  if(comments > 2)
    outputTree->Print();
  
  outputFile->Write();
  
  cout << endl;
  cout << " the new root file has been written " << endl;
  cout << endl;
  
  outputFile->Close();

  return 0;
  
} //Closes the program

Int_t TSim::Hits2Tree(TString inputFileName,
			    Int_t   comments){
  
  TString outputFileName;
  outputFileName = inputFileName + ".root";
    
  inputFileName  = inputFileName + ".text.out";
    
  cout << endl;
  cout << " Hits2Tree: " << endl;
  cout << endl;
  cout << " as you requested, I will now create a  " << endl;
  cout << " root file from list mode data for you " << endl;
  cout << endl;
  cout << " text file : " << inputFileName << endl;
  cout << " root rile : " << outputFileName   << endl;
  
  FILE  *inputFile = fopen(inputFileName,"r");
  
  // for reading in from text file
  Char_t line[128] = "??";
  
  TFile  outputFile(outputFileName,
		    "RECREATE",
		    "Root File with TTree of simulated hits");

  Int_t  currentEventNumber  = -999;
  Int_t  previousEventNumber = -999;
  
  TTree *tree = new TTree("tree",
			  "A beautiful birch");

  tree->Branch("nHits",&nHits,"nHits/I");
  tree->Branch("nParticles",&nParticles,"nParticles/I");
  
  tree->Branch("eventNumber",&eventNumber,"eventNumber[nHits]/L");
  tree->Branch("detector",&detector,"detector[nHits]/C");
  tree->Branch("detectorUnitID",&detectorUnitID,"detectorUnitID[nHits]/I");
  tree->Branch("energy",&energy,"energy[nHits]/F");
  tree->Branch("minTime",&minTime,"minTime[nHits]/F");
  tree->Branch("maxTime",&maxTime,"maxTime[nHits]/F");
  tree->Branch("x",&x,"x[nHits]/F");
  tree->Branch("y",&y,"y[nHits]/F");
  tree->Branch("z",&z,"z[nHits]/F");
  tree->Branch("nOriginalTracks",&nOriginalTracks,
	       "nOriginalTracks[nHits]/I");
  tree->Branch("originalParticleNumber",&originalParticleNumber,
	       "originalParticleNumber[nHits]/I");
  tree->Branch("nTracks",&nTracks,
	       "nTracks[nHits]/I");
  tree->Branch("particleNumber",&particleNumber,
	       "particleNumber[nHits][16]/I");
  
  // Initial initialisations
  for (Int_t i = 0; i < 16; i++){
    eventNumber[i]     = 0; 
    detectorUnitID[i]  = 0;
    energy[i]  = 0.;
    minTime[i] = 0.;
    maxTime[i] = 0.;
    x[i]       = 0.;
    y[i]       = 0.;
    z[i]       = 0.;
    
    nOriginalTracks[i]        = 0;
    originalParticleNumber[i] = 0;
    nTracks[i]                = 0;
    
    detector[i] = '?';

    for (Int_t j = 0; j < 16 ; j++)
      particleNumber[i][j] = 0;
    
  } // end of :  for (Int_t i = 0; i < 16 ...

  //_______________________________________
  //=======================================
  // Read from text file and write to TTree
  
  // Each original hit for each event
  // is written into a new array element
  
  // Multiple tracks per original hit 
  // are also written using particle
  // numbers
  
  Long64_t nTotalHits = 0;

  while (fgets(line,200,inputFile)) {
    
    nTotalHits++;

    if(comments){
      cout << endl;    
      cout << " line " << " = " << line << endl;
    }

    sscanf(&line[0],"%i",&currentEventNumber);
    
    if(comments){
      cout << endl;
      cout << " next event will be   " << currentEventNumber << endl;
    }
    
    // If the next line to be read in will be the very first, 
    // or will correspond to the same event as the previous line,
    // then fill the next array positions and continue
    // (or if next line will be a new event then fill tree
    // with previous event info, then read in next the line)
    
    if(currentEventNumber==previousEventNumber || 
       previousEventNumber==-999){
      
      nHits++;
      
      sscanf(&line[0],
	     "%lli %c %i %f %f %f %f %f %f %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i", 
	     &eventNumber[nHits-1], &detector[nHits-1],
	     &detectorUnitID[nHits-1], &energy[nHits-1],
	     &minTime[nHits-1], &maxTime[nHits-1],
	     &x[nHits-1], &y[nHits-1], &z[nHits-1],
	     &nOriginalTracks[nHits-1],
	     &originalParticleNumber[nHits-1],
	     &nTracks[nHits-1],
	     &particleNumber[nHits-1][0],
	     &particleNumber[nHits-1][1],
	     &particleNumber[nHits-1][2],
	     &particleNumber[nHits-1][3],
	     &particleNumber[nHits-1][4],
	     &particleNumber[nHits-1][5],
	     &particleNumber[nHits-1][6],
	     &particleNumber[nHits-1][7],
	     &particleNumber[nHits-1][8],
	     &particleNumber[nHits-1][9],
	     &particleNumber[nHits-1][10],
	     &particleNumber[nHits-1][11],
	     &particleNumber[nHits-1][12],
	     &particleNumber[nHits-1][13],
	     &particleNumber[nHits-1][14],
	     &particleNumber[nHits-1][15]);
	     
      if(comments){
      	cout << endl;
	cout << " previousEventNumber       = " << previousEventNumber     << endl;
	cout << endl;
	cout << " once per event  " << endl;
	cout << " (last value per event is stored) " << endl;
	cout << " nHits                     = " << nHits                   << endl;
	cout << endl;
	cout << " once per line " << endl;
	cout << " (all values stored in array)" << endl;
	cout << " eventNumber["<<(nHits-1)<<"]            = " << eventNumber[nHits-1]        << endl;
	cout << " energy["<<(nHits-1)<<"]                 = " << energy[nHits-1]        << endl;
	cout << " nOriginalTracks["<<(nHits-1)<<"] = " << nOriginalTracks[nHits-1] << endl;
	cout << " originalParticleNumber["<<(nHits-1)<<"] = " << originalParticleNumber[nHits-1] << endl;
	cout << " nTracks["<<(nHits-1)<<"]         = " << nTracks[nHits-1] << endl;
      }

      // set nParticles to the largest particle number in the event
      for (Int_t i = 0 ; i < 16 ; i++){
	if( particleNumber[nHits-1][i]!=0){
	  
	  if(comments) 
	    cout << " particleNumber[" <<(nHits-1)<< "," << i << "]    = " 
		 << particleNumber[nHits-1][i] << endl;
	  
	  if( particleNumber[nHits-1][i] > nParticles )
	    nParticles = particleNumber[nHits-1][i];
	}
	else{
	  break;
	}
      } // end of: for (Int_t i = 0 ; i < 16 ; 

      if( comments ){
	cout << endl;
	cout << " once per event  " << endl;
	cout << " (last value per event is stored) " << endl;
	cout << " nParticles                = " << nParticles              << endl;
      }
      
    } // end of :  if(currentEventNumber==previousEventNumber ||... ...previousEventNumber==-999){
    // if next line will be a new event then fill tree
    // with previous event info, then read in next the line
    else{       
      
      if(comments){
	cout << endl;
	cout << " Filling tree with event " << previousEventNumber << endl;
      }
      
      tree->Fill();

      // Re-initialise 
      nHits      = 1;
      nParticles = 1;
      for (Int_t i = 0; i < 16 ; i++){
	eventNumber[i]  = 0; 
	detectorUnitID[i]  = 0;
	energy[i]  = 0.;
	minTime[i] = 0.;
	maxTime[i] = 0.;
	x[i]    = -99.9;
	y[i]    = -99.9;
	z[i]    = -99.9;
	detector[i] = '?';

	nOriginalTracks[i] = 0;
	originalParticleNumber[i] = 0;
	nTracks[i] = 0;
    
	for (Int_t j = 0; j < 16; j++)
	  particleNumber[i][j] = 0;
      
      } 

      // read in next line, which will correspond 
      // to the first hit in a new event
      sscanf(&line[0],
	     "%lli %c %i %f %f %f %f %f %f %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i", 
	     &eventNumber[nHits-1], &detector[nHits-1],
	     &detectorUnitID[nHits-1], &energy[nHits-1],
	     &minTime[nHits-1], &maxTime[nHits-1],
	     &x[nHits-1], &y[nHits-1], &z[nHits-1],
	     &nOriginalTracks[nHits-1],
	     &originalParticleNumber[nHits-1],
	     &nTracks[nHits-1],
	     &particleNumber[nHits-1][0],
	     &particleNumber[nHits-1][1],
	     &particleNumber[nHits-1][2],
	     &particleNumber[nHits-1][3],
	     &particleNumber[nHits-1][4],
	     &particleNumber[nHits-1][5],
	     &particleNumber[nHits-1][6],
	     &particleNumber[nHits-1][7],
	     &particleNumber[nHits-1][8],
	     &particleNumber[nHits-1][9],
	     &particleNumber[nHits-1][10],
	     &particleNumber[nHits-1][11],
	     &particleNumber[nHits-1][12],
	     &particleNumber[nHits-1][13],
	     &particleNumber[nHits-1][14],
	     &particleNumber[nHits-1][15]);
      
      if(comments){
      	cout << endl;
	cout << " previousEventNumber       = " << previousEventNumber     << endl;
	cout << endl;
	cout << " once per event  " << endl;
	cout << " (last value per event is stored) " << endl;
	cout << " nHits                     = " << nHits                   << endl;
	cout << endl;
	cout << " once per line " << endl;
	cout << " (all values stored in array)" << endl;
	cout << " eventNumber["<<(nHits-1)<<"]            = " << eventNumber[nHits-1]        << endl;
	cout << " energy["<<(nHits-1)<<"]                 = " << energy[nHits-1]        << endl;
	cout << " nOriginalTracks["<<(nHits-1)<<"] = " << nOriginalTracks[nHits-1] << endl;
	cout << " originalParticleNumber["<<(nHits-1)<<"] = " << originalParticleNumber[nHits-1] << endl;
	cout << " nTracks["<<(nHits-1)<<"]         = " << nTracks[nHits-1] << endl;
      }
      
      for (Int_t i = 0 ; i < 16 ; i++){
	if( particleNumber[nHits-1][i]!=0){
	  
	  if(comments) 
	    cout << " particleNumber[" <<(nHits-1)<< "," << i << "]    = " 
		 << particleNumber[nHits-1][i] << endl;
	  
	  if( particleNumber[nHits-1][i] > nParticles )
	    nParticles = particleNumber[nHits-1][i];
	}
	else break;
     
      }
      if( comments){
	cout << endl;
	cout << " once per event  " << endl;
	cout << " (last value per event is stored) " << endl;
	cout << " nParticles                = " << nParticles              << endl;
      }
           
    } // end of: else{
    
    previousEventNumber = currentEventNumber;  
    
  } // end of: while (fgets(line,10....

  // The very last event
  // in the text file
  tree->Fill();
  
  if(comments) tree->Print();
  
  outputFile.Write();

  cout << endl;
  cout << " the root file has been written, there  " << endl;
  cout << " were " << nTotalHits << " hits in total " << endl;
  cout << endl;
  
  fclose(inputFile);

  outputFile.Close();
  
  return 0;
}



Float_t TSim::GetHitAngle(Float_t hitSeparation,
				Float_t distancetoOrigin){
  
  hitSeparation = hitSeparation/2.;
  
  Float_t angle = RadToDeg()*ACos(hitSeparation/distancetoOrigin);
  
  return angle;
}

Bool_t TSim::GoodNumHits(Int_t hits){
  
  Bool_t goodEvent = kFALSE;
  
  if (hits==2 || hits==4 )
    goodEvent = kTRUE;
  
  return goodEvent;
}

Bool_t TSim::GoodLOR(Float_t lor,
			   Char_t dataType){
  
  Bool_t goodEvent = kFALSE;
  
  Float_t minLOR = 200;
  
  if(dataType == 'l')    
    minLOR = 40;
  
  if (lor > minLOR )
    goodEvent = kTRUE;
  
  return goodEvent;
}

Bool_t TSim::GoodHitAngle(Float_t hitAngle,
				Float_t maxHitAngle){
  
  Bool_t goodEvent = kFALSE;
  
  if ( hitAngle < maxHitAngle )
    goodEvent = kTRUE;
  
  return goodEvent;
}


Bool_t TSim::GoodHitSeparation(Float_t hitSeparation){
  
  Bool_t goodEvent = kFALSE;
  
  if ( hitSeparation > 800. )
    goodEvent = kTRUE;
  
  return goodEvent;
}

Bool_t TSim::GoodDTheta(Float_t thetaV,
			      Float_t thetaE,
			      Float_t maxAngle){
  
  Bool_t goodEvent = kFALSE;
  
  if( Abs(thetaV-thetaE) < maxAngle &&
      thetaV < 179. &&
      thetaE < 179.
      ) 
    goodEvent = kTRUE;
  
  return goodEvent;
}

Bool_t TSim::GoodE(UInt_t  eventtype,
			 Float_t energy1,
			 Float_t energy2){
  
  Bool_t goodEvent = kFALSE;
  
  Float_t totalEnergy = energy1 + energy2;
  
  Float_t photopeakMin = 435.0;
  Float_t photopeakMax = 520.0;

  //------------------------------------
  // energy limits for scattered events
  //
  
  // Compton edge at (2/3 * 511) = 340.7
  Float_t ComptonEdge  = 340.7;
  
  ComptonEdge  = 345.0;
  
  // 170 is energy of scattered photon at 180 deg.
  // for energies below this, the first deposit cant have been 
  // a single Compton scatter as 180 is clearly the physical maximum
  // Float_t ComptonMin   = 170.0;

  // fluorescence peaks at 53, 54, 61.5 and 63. keV ?
  // Float_t safeMin      = 70.0;
  
  //
  //-------------------------------- 
  
  if     (eventtype==1){
    if( energy1 > photopeakMin && 
	energy1 < photopeakMax)      
      goodEvent = kTRUE;
  }
  else if(eventtype==2){
    if( energy1     > 100.0        &&
	energy1     < ComptonEdge  &&
	energy2     > 100.0        &&
	totalEnergy < photopeakMax 
	)
      goodEvent = kTRUE;
  }
  return goodEvent;
  
}

Bool_t TSim::GoodScatterDistance(Float_t scatterDistance){
  
  Bool_t goodEvent = kFALSE;
  
  if(scatterDistance < 4.1 && 
     scatterDistance > 3.9) 
    goodEvent = kTRUE;
  
  return goodEvent;
}

Bool_t TSim::GoodScatterDistances(Float_t scatterDistance1,
					Float_t scatterDistance2 
					){
  
  Bool_t goodEvent = kFALSE;
  
  if(scatterDistance1 < 0.0 && 
     scatterDistance1 > 1000.0  && 
     scatterDistance2 < 0.0 &&  
     scatterDistance2 > 1000.0 ) 
  goodEvent = kTRUE;
  
  return goodEvent;
}

Bool_t TSim::BadDPhi(Float_t dP){
  
  Bool_t badEvent = kFALSE;
  
  if(dP < -0.1 || dP > 0.1)
    badEvent = kTRUE;
  
  return badEvent;
}


Bool_t TSim::GoodTheta(Float_t theta){
  
  Bool_t goodTheta = kFALSE;
  
  if( theta >= 0.0 &&
      theta <  180.)
    goodTheta = kTRUE;
  
  return goodTheta;
}


void TSim::SetAsymmetry(TString inputFileNumber){
  
  cout << endl;
  cout << " Setting asymmetry " << endl;
  
  this->GetAsymmetry(inputFileNumber);
  
}

Int_t TSim::GetAsymmetry(TString inputFileNumber){
  
  cout << endl;
  cout << " Getting asymmetry " << endl;
  
  this->SetStyle();

  TString plotName;
  plotName = "../Plots/Asym_" + inputFileNumber;

  plotName = plotName + ".pdf";

  TString inputFileName = "../Data/hitsLab" + inputFileNumber;
  
  TString outputFileName;
  outputFileName = inputFileName + ".goodPet.root";
    
  inputFileName  = inputFileName + ".pet.root";
    
  cout << endl;
  cout << " Input File : " << inputFileName  << endl;
  cout << " Output File: " << outputFileName << endl;
  cout << endl;

  TFile* inputFile = new TFile(inputFileName);
  
  Float_t theta1P1;
  Float_t theta1P2;
  Float_t theta2P1;
  Float_t theta2P2;

  //-------------------------------
  // Initial array initialisations
  for (Int_t i = 0; i < 32; i++){
    
    eventNumber[i] = 0;
    detectorUnitID[i] = 0;
    
    detector[i] = '?';
    
    energy[i] = -999.0;

    x[i] = -999.9;
    y[i] = -999.9;
    z[i] = -999.9;

    nOriginalTracks[i] = 0;
    originalParticleNumber[i] = 0;
  
    nTracks[i] = 0;
    
  }
  
  energy1P1 = -999.0;
  energy1P2 = -999.0;
  energy2P1 = -999.0;
  energy2P2 = -999.0;
  
  theta1P1  = -999.0;
  theta1P2  = -999.0;
  theta2P1  = -999.0;
  theta2P2  = -999.0;

  thetaP1E   = -999.9;
  thetaP2E   = -999.9;

  thetaP1V   = -999.9;
  thetaP2V   = -999.9;

  thetaP1VTrue = -999.9;
  thetaP2VTrue = -999.9;

  Long64_t nEvents = 0;

  TTree* tree=(TTree*)inputFile->Get("phiTree");

  // Variables from Hits2Tree.C
  tree->SetBranchAddress("nHits",&nHits);
  tree->SetBranchAddress("nHitsC",&nHitsC);
  tree->SetBranchAddress("nHitsW",&nHitsW);
  
  tree->SetBranchAddress("nParticles",&nParticles);
  
  tree->SetBranchAddress("eventNumber",&eventNumber);
  tree->SetBranchAddress("detector",&detector);
  tree->SetBranchAddress("detectorUnitID",&detectorUnitID);
  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("x",&x);
  tree->SetBranchAddress("y",&y);
  tree->SetBranchAddress("z",&z);
  
  tree->SetBranchAddress("nOriginalTracks",&nOriginalTracks);
  tree->SetBranchAddress("originalParticleNumber",&originalParticleNumber);
  tree->SetBranchAddress("nTracks",&nTracks);

  // Variables from SortEvents.C
  tree->SetBranchAddress("eventType",&eventType);

  tree->SetBranchAddress("i1P1",&i1P1);
  tree->SetBranchAddress("i2P1",&i2P1);

  tree->SetBranchAddress("sinogramR",&sinogramR);
  tree->SetBranchAddress("sinogramPhi",&sinogramPhi);

  tree->SetBranchAddress("meanHitAngle",&meanHitAngle);

  tree->SetBranchAddress("energy1P1",&energy1P1);
  tree->SetBranchAddress("energy1P2",&energy1P2);

  //used for scattered PET events only

  tree->SetBranchAddress("i1P2",&i1P2);
  tree->SetBranchAddress("i2P2",&i2P2);

  tree->SetBranchAddress("energy2P1",&energy2P1);
  tree->SetBranchAddress("energy2P2",&energy2P2);

  tree->SetBranchAddress("thetaP1E",&thetaP1E);
  tree->SetBranchAddress("thetaP2E",&thetaP2E);

  tree->SetBranchAddress("thetaP1V",&thetaP1V);
  tree->SetBranchAddress("thetaP2V",&thetaP2V);

  tree->SetBranchAddress("thetaP1VTrue",&thetaP1VTrue);
  tree->SetBranchAddress("thetaP2VTrue",&thetaP2VTrue);

  tree->SetBranchAddress("phiP1",&phiP1);
  tree->SetBranchAddress("phiP2",&phiP2);
  tree->SetBranchAddress("dPhi",&dPhi);
  
  nEvents = (UInt_t)tree->GetEntries();

  TTree *outputTree = new TTree("goodPhiTree","Output from DeltaPhi.C");
    
  // Variables from Hits2Tree.C
  outputTree->Branch("nHits",&nHits,"nHits/I");
  outputTree->Branch("nHitsC",&nHitsC,"nHitsC/I");
  outputTree->Branch("nHitsW",&nHitsW,"nHitsW/I");
  
  outputTree->Branch("nParticles",&nParticles,"nParticles/I");
  
  outputTree->Branch("eventNumber",&eventNumber,"eventNumber[nHits]/L");
  outputTree->Branch("detector",&detector,"detector[nHits][128]/C");
  outputTree->Branch("detectorUnitID",&detectorUnitID,"detectorUnitID[nHits]/I");
  outputTree->Branch("energy",&energy,"energy[nHits]/F");
  outputTree->Branch("x",&x,"x[nHits]/F");
  outputTree->Branch("y",&y,"y[nHits]/F");
  outputTree->Branch("z",&z,"z[nHits]/F");
  
  outputTree->Branch("nOriginalTracks",&nOriginalTracks,
	       "nOriginalTracks[nHits]/I");
  outputTree->Branch("originalParticleNumber",&originalParticleNumber,
	       "originalParticleNumber[nHits]/I");
  outputTree->Branch("nTracks",&nTracks,
	       "nTracks[nHits]/I");

  // New Variables
  
  // used for all PET events
  outputTree->Branch("eventType",&eventType,"eventType/I");
  
  outputTree->Branch("i1P1",&i1P1,"i1P1/I");
  outputTree->Branch("i1P2",&i1P2,"i1P2/I");

  outputTree->Branch("p3P1Lab","TVector3",&p3P1Lab);
  outputTree->Branch("p3P2Lab","TVector3",&p3P2Lab);

  outputTree->Branch("LOR","TVector3",&LOR);
  outputTree->Branch("LORXY","TVector3",&LORXY);

  outputTree->Branch("v3P1Lab","TVector3",&v3P1Lab);
  outputTree->Branch("v3P2Lab","TVector3",&v3P2Lab);

  outputTree->Branch("midpoint","TVector3",&midpoint);
  outputTree->Branch("midpointXY","TVector2",&midpointXY);

  outputTree->Branch("sinogramR",&sinogramR,"sinogramR/F");
  outputTree->Branch("sinogramPhi",&sinogramPhi,"sinogramPhi/F");

  outputTree->Branch("meanHitAngle",&meanHitAngle,"meanHitAngle/F");
  
  outputTree->Branch("energy1P1",&energy1P1,"energy1P1/F");
  outputTree->Branch("energy1P2",&energy1P2,"energy1P2/F");

  //used for scattered PET events only

  outputTree->Branch("i2P1",&i2P1,"i2P1/I");
  outputTree->Branch("i2P2",&i2P2,"i2P2/I");

  outputTree->Branch("energy2P1",&energy2P1,"energy2P1/F");
  outputTree->Branch("energy2P2",&energy2P2,"energy2P2/F");
  
  outputTree->Branch("thetaP1E",&thetaP1E,"thetaP1E/F");
  outputTree->Branch("thetaP2E",&thetaP2E,"thetaP2E/F");

  outputTree->Branch("v3ScP1Lab","TVector3",&v3ScP1Lab);
  outputTree->Branch("v3ScP2Lab","TVector3",&v3ScP2Lab);

  outputTree->Branch("v3P1Event","TVector3",&v3P1Event);
  outputTree->Branch("v3P2Event","TVector3",&v3P2Event);

  outputTree->Branch("thetaP1V",&thetaP1V,"thetaP1V/F");
  outputTree->Branch("thetaP2V",&thetaP2V,"thetaP2V/F");
  
  outputTree->Branch("phiP1",&phiP1,"phiP1/F");
  outputTree->Branch("phiP2",&phiP2,"phiP2/F");
  outputTree->Branch("dPhi",&dPhi,"dPhi/F");
    
  Bool_t blockCentres = kFALSE;
  
  Bool_t blockAZero  = kFALSE;
  Bool_t blockANinty = kFALSE;
  Bool_t blockAPi    = kFALSE;
  Bool_t blockA270   = kFALSE;
  
  Bool_t blockBZero  = kFALSE;
  Bool_t blockBNinty = kFALSE;
  Bool_t blockBPi    = kFALSE;
  Bool_t blockB270   = kFALSE;
  
  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);
 
  canvas->Divide(3,2);
  
  canvas->cd(1);

  gPad->SetTickx();
  gPad->SetTicky(1);
  
  TGraphErrors *grAsym[3];
  TGraphErrors *grAsymCr[3];
  TGraphErrors *grCoefCr[3];
  TGraphErrors *grCoef[3];
  
  TH1F *deltaPhi = new TH1F("deltaPhi",
			    "deltaPhi;#Delta#phi (deg);relative counts",
			    16,-180.,180.);

 //   TH1F *dPFitErr = new TH1F("dPFitErr",
// 			    "dPFitErr;#Delta#phi;Counts",
// 			    3,-180.,180.);
  
  const Int_t nBins = 2;
  
  Int_t nZero[nBins];
  Int_t nNinty[nBins];
  Int_t nPi[nBins];
  Int_t n270[nBins];

  Int_t nZeroCr[nBins];
  Int_t nNintyCr[nBins];
  Int_t nPiCr[nBins];
  Int_t n270Cr[nBins];
  
  for( Int_t i = 0 ; i < nBins ; i++){
    nZero[i]   = 0;
    nNinty[i]  = 0;
    nPi[i]     = 0;
    n270[i]    = 0;
    
    nZeroCr[i]  = 0;
    nNintyCr[i] = 0;
    nPiCr[i]    = 0;
    n270Cr[i]   = 0;
    
  } 
  Float_t theta[nBins];
  
  Float_t asymm180[nBins];
  Float_t asyer180[nBins];
  
  Float_t asymm270[nBins];
  Float_t asyer270[nBins];
  
  Float_t asymmCr180[nBins];
  Float_t asyerCr180[nBins];

  Float_t asymmCr270[nBins];
  Float_t asyerCr270[nBins];
  
  // Function is
  // A + Bcos(phi) + Ccos2(phi)
  
  Float_t A[nBins];
  Float_t B[nBins];
  Float_t C[nBins];

  Float_t ACrerr[nBins];
  Float_t BCrerr[nBins];
  Float_t CCrerr[nBins];
  
  Float_t ACr[nBins];
  Float_t BCr[nBins];
  Float_t CCr[nBins];
  
  Float_t Aerr[nBins];
  Float_t Berr[nBins];
  Float_t Cerr[nBins];
  
  // Amplitudes at:
  // dPhi = 0
  Float_t ApBpC[nBins];
  // dPhi = 90
  Float_t AmC[nBins];
  // dPhi = 180
  Float_t AmBpC[nBins];
  // dPhi = 0
  Float_t ApBpCCr[nBins];
  // dPhi = 90
  Float_t AmCCr[nBins];
  // dPhi = 180
  Float_t AmBpCCr[nBins];

  for(Int_t i = 0 ; i < nBins ; i++){
    //!!
    //theta[i]  = 10. + 20*i;
    theta[i] = 90.;
    
    A[i] = 0.;
    B[i] = 0.;
    C[i] = 0.;
    
    Aerr[i] = 0.0;
    Berr[i] = 0.0;
    Cerr[i] = 0.0;

    ACr[i] = 0.;
    BCr[i] = 0.;
    CCr[i] = 0.;
    
    ACrerr[i] = 0.0;
    BCrerr[i] = 0.0;
    CCrerr[i] = 0.0;
    
    ApBpC[i] = 0.;
    AmC[i]   = 0.;
    AmBpC[i] = 0.;

    ApBpCCr[i] = 0.;
    AmCCr[i]   = 0.;
    AmBpCCr[i] = 0.;
    
    asymm90[i]   = 0.;
    asyer90[i]   = 0.;
    asymm180[i]  = 0.;
    asyer180[i]  = 0.;
    asymm270[i]  = 0.;
    asyer270[i]  = 0.;
    
    asymmCr90[i]  = 0.;
    asyerCr90[i]  = 0.;
    asymmCr180[i] = 0.;
    asyerCr180[i] = 0.;
    asymmCr270[i] = 0.;
    asyerCr270[i] = 0.;
  }
     
  //--------------------------------------//
  //------------- event loop -------------//
  //--------------------------------------//

  Long64_t events[8];
  
  for (Int_t i = 0 ; i < 8 ; i++)
    events[i] = 0;
    
  const Float_t rangeDPhi = 45.;

  Int_t dPhiCr = 0;
  
  for (Int_t i = 0 ; i < nEvents ; i++){
    
    blockCentres = kFALSE;
    blockAZero   = kFALSE;
    blockANinty  = kFALSE;
    blockAPi     = kFALSE;
    blockA270    = kFALSE;
    blockBZero   = kFALSE;
    blockBNinty  = kFALSE;
    blockBPi     = kFALSE;
    blockB270    = kFALSE;
    dPhiCr = 0;

    events[0]++;
    
    tree->GetEvent(i);
    
    //============================================
    // only retain events 
    // containing four hits
    
    if( nHitsC!=4 ) continue;
    
    events[1]++;
    
    //============================================
    
    theta1P1 = ElectronEnergyToTheta(energy1P1); 
    theta2P1 = PhotonEnergyToTheta(energy2P1);

    theta1P2 = ElectronEnergyToTheta(energy1P2); 
    theta2P2 = PhotonEnergyToTheta(energy2P2);

    //============================================
    // only retain events in which 
    // central crystals were hit

    if( (detectorUnitID[0] == 30004  && 
	 detectorUnitID[2] == 40004) ||
	(detectorUnitID[2] == 30004  && 
	 detectorUnitID[0] == 40004)
	){
      
      for(Int_t j = 0; j < nHits ; j++)
	detectorUnitID[j] = detectorUnitID[j]-20000;
      
      blockCentres = kTRUE;
    }
    else if( (detectorUnitID[0] == 50004  && 
	      detectorUnitID[2] == 60004) ||
	     (detectorUnitID[2] == 50004  && 
	      detectorUnitID[0] == 60004)
	     ){
      
      for(Int_t j = 0; j < nHits ; j++)
	detectorUnitID[j] = detectorUnitID[j]-40000;
      
      blockCentres = kTRUE;
    }
    else if( (detectorUnitID[0] == 10004  && 
	      detectorUnitID[2] == 20004) ||
	     (detectorUnitID[2] == 10004  && 
	      detectorUnitID[0] == 20004)
	     ){
      blockCentres = kTRUE;
    }
    
    if( !blockCentres ) continue;
    
    // put first particle in first block
    if( detectorUnitID[0] == 20004 ){
      detectorUnitID[0] = detectorUnitID[0] - 10000;
      detectorUnitID[1] = detectorUnitID[1] - 10000;
      detectorUnitID[2] = detectorUnitID[2] + 10000;
      detectorUnitID[3] = detectorUnitID[3] + 10000;
    }
        
    events[2]++;
    
    //============================================

    //============================================
    // Match theta for two particles

    if( !GoodTheta(theta1P1) || !GoodTheta(theta1P2) )
      continue;
    
    if( GetThetaBin(theta1P1,nBins) != GetThetaBin(theta1P2,nBins) )
      continue;

    // if( Abs(theta1P1 - theta1P2 ) > 10. )
    //   continue;

    events[3]++;

    // deltaPhi plot for one bin
    if( theta1P1 >= 0.0 && theta1P1 < 180.0)
      deltaPhi->Fill(dPhi);
        
    //============================================
    //--------------------------------------------
    // hit position asymmetry analysis
    
    if      (GetDPhiBin(dPhi,rangeDPhi) == 0  ) nZero[GetThetaBin(theta1P1,nBins)]++;
    else if (GetDPhiBin(dPhi,rangeDPhi) == 90 ) nNinty[GetThetaBin(theta1P1,nBins)]++;
    else if (GetDPhiBin(dPhi,rangeDPhi) == 180) nPi[GetThetaBin(theta1P1,nBins)]++;
    else if (GetDPhiBin(dPhi,rangeDPhi) == 270) n270[GetThetaBin(theta1P1,nBins)]++;

    
    //--------------------------------------------
    //============================================

    
    //============================================
    //In which crystals did the second hits occur?

    //============================================
    // where is second hit in the first block ?
    if     (detectorUnitID[1] == 10001)
      blockAZero   = kTRUE;
    else if(detectorUnitID[1] == 10003)
      blockANinty = kTRUE;
    else if(detectorUnitID[1] == 10007)
      blockAPi = kTRUE;
    else if(detectorUnitID[1] == 10005)
      blockA270 = kTRUE;
    else continue;
    
    events[4]++;
    
    // where is second hit in the second block ?
    if     (detectorUnitID[3] == 20001)
      blockBZero   = kTRUE;
    else if(detectorUnitID[3] == 20003)
      blockBNinty = kTRUE;
    else if(detectorUnitID[3] == 20007)
      blockBPi = kTRUE;
    else if(detectorUnitID[3] == 20005)
      blockB270 = kTRUE;
    else continue;

    events[5]++;
            
    if      ( (blockAZero  && blockBZero)  ||
	      (blockANinty && blockBNinty) ||
	      (blockAPi    && blockBPi)    ||
	      (blockA270   && blockB270) 
	      ) dPhiCr = 0;
    else if ( (blockAZero  && blockBNinty) ||
	      (blockANinty && blockBPi)    ||
	      (blockAPi    && blockB270)   ||
	      (blockA270   && blockBZero) 
	      ) dPhiCr = 90;
    else if ( (blockAZero  && blockBPi)    ||
	      (blockANinty && blockB270)   ||
	      (blockAPi    && blockBZero)  ||
	      (blockA270   && blockBNinty) 
	      ) dPhiCr = 180;
    else if ( (blockAZero  && blockB270)   ||
	      (blockANinty && blockBZero)  ||
	      (blockAPi    && blockBNinty) ||
	      (blockA270   && blockBPi) 
	      ) dPhiCr = 270;

    if ( dPhiCr == 0 ){
      nZeroCr[GetThetaBin(theta1P1,nBins)]++; 
    }
    else if( dPhiCr == 90 ){
      nNintyCr[GetThetaBin(theta1P1,nBins)]++; 
    }
    else if( dPhiCr == 180 ){
      nPiCr[GetThetaBin(theta1P1,nBins)]++; 
    }
    else if( dPhiCr == 270 ){
      n270Cr[GetThetaBin(theta1P1,nBins)]++; 
    }
    else{
      cout << endl;
      cout << " NO!! " << endl;
    }
      
  } // end of: for (Int_t i = 0 ; i < nEve.....
  //========================================================
  //========================================================
  
  cout << endl;
  for( Int_t i = 0 ; i < 6 ; i++ )
    cout << " " << events[i] << "\t events after cut " << i << endl;
  
    
  for(Int_t i = 0 ; i < nBins ; i++ ){

    asymm90[i] = (Float_t)(nNinty[i])/(nZero[i]);
    asyer90[i] = asymm90[i]*Sqrt((nNinty[i]+nZero[i])/(Float_t)(nNinty[i]*nZero[i]));
    
    asymm180[i] = (Float_t)nPi[i]/nZero[i];
    asyer180[i] = asymm180[i]*Sqrt((nPi[i]+nZero[i])/(Float_t)(nPi[i]*nZero[i]));
    
    asymm270[i] = (Float_t)n270[i]/nZero[i];
    asyer270[i] = asymm270[i]*Sqrt((n270[i]+nZero[i])/(Float_t)(n270[i]*nZero[i]));

    asymmCr90[i] = (Float_t)(nNintyCr[i])/(nZeroCr[i]);
    asyerCr90[i] = asymmCr90[i]*Sqrt((nNintyCr[i]+nZeroCr[i])/(Float_t)(nNintyCr[i]*nZeroCr[i]));
    
    asymmCr180[i] = (Float_t)nPiCr[i]/nZeroCr[i];
    asyerCr180[i] = asymmCr180[i]*Sqrt((nPiCr[i]+nZeroCr[i])/(Float_t)(nPiCr[i]*nZeroCr[i]));

    asymmCr270[i] = (Float_t)n270Cr[i]/nZeroCr[i];
    asyerCr270[i] = asymmCr270[i]*Sqrt((n270Cr[i]+nZeroCr[i])/(Float_t)(n270Cr[i]*nZeroCr[i]));

    // error is Sqrt(2/nZero[i]);
    ApBpC[i] = 1.;
    
    // error is asyer90[i]
    AmC[i]   = asymm90[i];
    
    // error is  = asyer180[i]
    AmBpC[i] = asymm180[i];
   
    B[i]    = 1/2.*(ApBpC[i] - AmBpC[i]);
    Berr[i] = 1/2.*Sqrt( asyer180[i]*asyer180[i] );
    
    C[i]    = 1/4.*(ApBpC[i] - 2*AmC[i] + AmBpC[i]);
    Cerr[i] = 1/4.*Sqrt( 4*asyer90[i] + asyer180[i]*asyer180[i]);
    
    A[i] = ApBpC[i] - B[i] - C[i];
    Aerr[i] = Sqrt( Berr[i]*Berr[i] + Cerr[i]*Cerr[i] );
          
    ApBpCCr[i] = 1.;
    AmCCr[i]   = asymmCr90[i];
    AmBpCCr[i] = asymmCr180[i];
    
    BCr[i] = 1/2.*(ApBpCCr[i] - AmBpCCr[i]);
    BCrerr[i] = 1/2.*Sqrt( asyerCr180[i]*asyerCr180[i] );
    
    CCr[i] = 1/4.*(ApBpCCr[i] - 2*AmCCr[i] + AmBpCCr[i]);
    CCrerr[i] = 1/4.*Sqrt( 2./nZeroCr[i] + 4*asyerCr90[i] + asyerCr180[i]*asyerCr180[i]);
    
    ACr[i] = ApBpCCr[i] - BCr[i] - CCr[i];
    ACrerr[i] = Sqrt( BCrerr[i]*BCrerr[i] + CCrerr[i]*CCrerr[i] );
    
  }

  grAsym[0] =  new TGraphErrors(nBins,theta,asymm90,0,asyer90);
  grAsym[1] =  new TGraphErrors(nBins,theta,asymm180,0,asyer180);
  grAsym[2] =  new TGraphErrors(nBins,theta,asymm270,0,asyer270);

  grAsym[0]->SetMarkerSize(0.5);
  grAsym[0]->SetMarkerStyle(20);
  grAsym[0]->SetMarkerColor(kRed);
  grAsym[0]->SetLineColor(kRed);
  grAsym[0]->SetFillColor(kRed-2);
  grAsym[0]->SetFillStyle(3003);
  grAsym[0]->SetName("grAsym[0]");

  grAsym[1]->SetMarkerSize(0.5);
  grAsym[1]->SetMarkerStyle(20);
  grAsym[1]->SetMarkerColor(kBlue);
  grAsym[1]->SetLineColor(kBlue);
  grAsym[1]->SetFillColor(kBlue-2);
  grAsym[1]->SetFillStyle(3004);
  grAsym[1]->SetName("grAsym[1]");

  grAsym[2]->SetMarkerSize(0.5);
  grAsym[2]->SetMarkerStyle(20);
  grAsym[2]->SetMarkerColor(kGreen+2);
  grAsym[2]->SetLineColor(kGreen+2);
  grAsym[2]->SetFillColor(kGreen+2);
  grAsym[2]->SetFillStyle(3006);
  grAsym[2]->SetName("grAsym[2]");
  
  grAsym[0]->Draw("AP L E3");
  grAsym[0]->GetXaxis()->SetTitle("#theta (deg)");

  grAsym[0]->GetYaxis()->SetTitle("A(#Delta#phi) - hit positions");
  grAsym[0]->GetYaxis()->SetRangeUser(0.7,1.9);
  grAsym[0]->GetXaxis()->SetRangeUser(0.0,180.);
  grAsym[0]->Draw("AP L E3");
  
  grAsym[1]->Draw("same P L E3");
  grAsym[2]->Draw("same P L E3");
  
  TLegend *leg = new TLegend(0.25,0.7,0.4,0.85);
  leg->SetTextFont(132);
  leg->SetTextSize(0.03);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetBorderSize(1);
  leg->AddEntry("grAsym[0]","90","p");
  leg->AddEntry("grAsym[1]","180","p");
  leg->AddEntry("grAsym[2]","270","p");
  leg->Draw();

  grAsymCr[0] =  new TGraphErrors(nBins,theta,asymmCr90,0,asyerCr90);
  grAsymCr[1] =  new TGraphErrors(nBins,theta,asymmCr180,0,asyerCr180);
  grAsymCr[2] =  new TGraphErrors(nBins,theta,asymmCr270,0,asyerCr270);
  
  grAsymCr[0]->SetMarkerSize(0.5);
  grAsymCr[0]->SetMarkerStyle(20);
  grAsymCr[0]->SetMarkerColor(kRed);
  grAsymCr[0]->SetLineColor(kRed);
  grAsymCr[0]->SetFillStyle(3003);
  grAsymCr[0]->SetFillColor(kRed-2);
  grAsymCr[0]->SetName("grAsymCr[0]");

  grAsymCr[1]->SetMarkerSize(0.5);
  grAsymCr[1]->SetMarkerStyle(20);
  grAsymCr[1]->SetMarkerColor(kBlue);
  grAsymCr[1]->SetLineColor(kBlue);
  grAsymCr[1]->SetFillStyle(3004);
  grAsymCr[1]->SetFillColor(kBlue-2);
  grAsymCr[1]->SetName("grAsymCr[1]");

  grAsymCr[2]->SetMarkerSize(0.5);
  grAsymCr[2]->SetMarkerStyle(20);
  grAsymCr[2]->SetMarkerColor(kGreen+2);
  grAsymCr[2]->SetLineColor(kGreen+2);
  grAsymCr[2]->SetFillStyle(3006);
  grAsymCr[2]->SetFillColor(kGreen+2);
  grAsymCr[2]->SetName("grAsymCr[2]");

  canvas->cd(2);
  
  gPad->SetTickx();
  gPad->SetTicky(1);

  grAsymCr[0]->Draw("AP");
  grAsymCr[0]->GetXaxis()->SetTitle("#theta (deg)");
  grAsymCr[0]->GetYaxis()->SetTitle("A(#Delta#phi) - crystal positions");
  grAsymCr[0]->GetXaxis()->SetRangeUser(0.0,180.);
  grAsymCr[0]->GetYaxis()->SetRangeUser(0.7,1.9);
  grAsymCr[0]->Draw("AP L E3");

  
  grAsymCr[1]->Draw("same P L E3");
  grAsymCr[2]->Draw("same P L E3");

  
  TLegend *leg2 = new TLegend(0.25,0.7,0.4,0.85);
  leg2->SetTextFont(132);
  leg2->SetTextSize(0.03);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->SetBorderSize(1) ;
  leg2->AddEntry("grAsymCr[0]","90","p");
  leg2->AddEntry("grAsymCr[1]","180","p");
  leg2->AddEntry("grAsymCr[2]","270","p");
  leg2->Draw();

  //================ 
 
  canvas->cd(3);

  gPad->SetTickx();
  gPad->SetTicky(1);
  
  deltaPhi->Sumw2();
  deltaPhi->Scale(2./(deltaPhi->GetBinContent(8)+deltaPhi->GetBinContent(9)));
  deltaPhi->SetMinimum(0.8);
  deltaPhi->SetMaximum(1.8);
  
  deltaPhi->SetLineWidth(1);
  deltaPhi->Draw("hist");

  canvas->cd(4);
  
  gPad->SetTickx();
  gPad->SetTicky(1);

  grCoef[0]  =  new TGraphErrors(nBins,theta,A,0,Aerr);
  grCoefCr[0]  =  new TGraphErrors(nBins,theta,ACr,0,ACrerr);
  
  grCoefCr[0]->SetMarkerSize(0.5);
  grCoefCr[0]->SetMarkerStyle(20);
  grCoefCr[0]->SetMarkerColor(kRed);
  grCoefCr[0]->SetLineColor(kRed);
  grCoefCr[0]->SetFillColor(kRed-2);
  grCoefCr[0]->SetFillStyle(3003);
  grCoefCr[0]->Draw("AP");
  grCoefCr[0]->GetXaxis()->SetTitle("#theta (deg)");
  grCoefCr[0]->GetYaxis()->SetTitle("constant coefficient");
  grCoefCr[0]->GetXaxis()->SetRangeUser(0.0,180.);
  grCoefCr[0]->GetYaxis()->SetRangeUser(0.5,1.5);
  grCoefCr[0]->Draw("AP L E3 ");

  grCoef[0]->SetMarkerSize(0.5);
  grCoef[0]->SetMarkerStyle(20);
  grCoef[0]->SetMarkerColor(kBlue);
  grCoef[0]->SetLineColor(kBlue);
  grCoef[0]->SetFillColor(kBlue-2);
  grCoef[0]->SetFillStyle(3004);
  grCoef[0]->Draw("E3 P L same");

  TLegend *leg3 = new TLegend(0.25,0.75,0.4,0.85);
  leg3->SetTextFont(132);
  leg3->SetTextSize(0.03);
  leg3->SetFillColor(kWhite);
  leg3->SetLineColor(kWhite);
  leg3->SetBorderSize(1) ;
  leg3->AddEntry(grCoefCr[0],"crystal positions","p");
  leg3->AddEntry(grCoef[0],"hit positions","p");
  leg3->Draw();
  
  
  canvas->cd(5);
  
  gPad->SetTickx();
  gPad->SetTicky(1);

  grCoef[1]  =  new TGraphErrors(nBins,theta,B,0,Berr);
  grCoefCr[1]  =  new TGraphErrors(nBins,theta,BCr,0,BCrerr);

  grCoefCr[1]->SetMarkerSize(0.5);
  grCoefCr[1]->SetMarkerStyle(20);
  grCoefCr[1]->SetMarkerColor(kRed);
  grCoefCr[1]->SetLineColor(kRed);
  grCoefCr[1]->SetFillColor(kRed-2);
  grCoefCr[1]->SetFillStyle(3003);
  grCoefCr[1]->Draw("AP");
  grCoefCr[1]->GetXaxis()->SetTitle("#theta (deg)");
  grCoefCr[1]->GetYaxis()->SetTitle("cos(#Delta#phi) coefficient");
  grCoefCr[1]->GetXaxis()->SetRangeUser(0.0,180.);
  grCoefCr[1]->GetYaxis()->SetRangeUser(-0.5,0.5);
  grCoefCr[1]->Draw("AP L E3");
  
  grCoef[1]->SetMarkerSize(0.5);
  grCoef[1]->SetMarkerStyle(20);
  grCoef[1]->SetMarkerColor(kBlue);
  grCoef[1]->SetLineColor(kBlue);
  grCoef[1]->SetFillColor(kBlue-2);
  grCoef[1]->SetFillStyle(3004);
  grCoef[1]->Draw("P L E3 same");

  
  TLegend *leg4 = new TLegend(0.25,0.75,0.4,0.85);
  leg4->SetTextFont(132);
  leg4->SetTextSize(0.03);
  leg4->SetFillColor(kWhite);
  leg4->SetLineColor(kWhite);
  leg4->SetBorderSize(1) ;
  leg4->AddEntry(grCoefCr[0],"crystal positions","p");
  leg4->AddEntry(grCoef[0],"hit positions","p");
  leg4->Draw();

  canvas->cd(6);
  
  gPad->SetTickx();
  gPad->SetTicky(1);
  
  grCoef[2]  =  new TGraphErrors(nBins,theta,C,0,Cerr);
  grCoefCr[2]  =  new TGraphErrors(nBins,theta,CCr,0,CCrerr);
  
  grCoefCr[2]->Draw("AP");
  grCoefCr[2]->GetXaxis()->SetTitle("#theta (deg)");
  grCoefCr[2]->GetYaxis()->SetTitle("cos(2#Delta#phi) coefficient");
  grCoefCr[2]->GetXaxis()->SetRangeUser(0.0,180.);
  grCoefCr[2]->GetYaxis()->SetRangeUser(-0.5,0.5);
  
  grCoefCr[2]->SetMarkerSize(0.5);
  grCoefCr[2]->SetMarkerStyle(20);
  grCoefCr[2]->SetMarkerColor(kRed);
  grCoefCr[2]->SetLineColor(kRed);
  grCoefCr[2]->SetFillColor(kRed-3);
  grCoefCr[2]->SetFillStyle(3003);
  grCoefCr[2]->Draw("AP L E3");
  
  grCoef[2]->SetMarkerSize(0.5);
  grCoef[2]->SetMarkerStyle(20);
  grCoef[2]->SetMarkerColor(kBlue);
  grCoef[2]->SetLineColor(kBlue);
  grCoef[2]->SetFillColor(kBlue-2);
  grCoef[2]->SetFillStyle(3004);
  
  grCoef[2]->Draw("P L E3 same");

  TLegend *leg5 = new TLegend(0.25,0.75,0.4,0.85);
  leg5->SetTextFont(132);
  leg5->SetTextSize(0.03);
  leg5->SetFillColor(kWhite);
  leg5->SetLineColor(kWhite);
  leg5->SetBorderSize(1) ;
  leg5->AddEntry(grCoefCr[0],"crystal positions","p");
  leg5->AddEntry(grCoef[0],"hit positions","p");
  leg5->Draw();

  
  canvas->cd(3);
  
  TLegend *leg6 = new TLegend(0.25,0.70,0.4,0.85);
  leg6->SetTextSize(0.03);
  leg6->SetTextFont(132);
  leg6->SetFillColor(kWhite);
  leg6->SetLineColor(kWhite);
  leg6->SetBorderSize(1) ;

  leg6->SetName("Bin 4: (60^{o} - 80^{o}");
  leg6->AddEntry(deltaPhi,"directly from hit positions","l");
  leg6->AddEntry(grCoefCr[0],"from coefs using crystal positions","p");
  leg6->AddEntry(grCoef[0],"from coefs using hit positions","p");


  leg6->Draw();
  
  TF1 *fun1 = new TF1("fun1","[0] + [1]*cos(TMath::DegToRad()*x) + [2]*cos(2*TMath::DegToRad()*x)",-180,180);
  
  fun1->SetLineColor(kBlue);
  fun1->SetMarkerColor(kBlue);
  fun1->SetMarkerStyle(2);
  fun1->SetLineWidth(1);
  
  //!!!
  const UShort_t index = 0;

  fun1->SetParameter(0,A[index]);
  fun1->SetParameter(1,B[index]);
  fun1->SetParameter(2,C[index]);
  
  fun1->SetParError(0,Aerr[index]);
  fun1->SetParError(1,Berr[index]);
  fun1->SetParError(2,Cerr[index]);
  
  fun1->Draw("same");
  
  const Int_t nBins2 = 5;

  Float_t dPhi1[nBins2];

  Float_t coefs[nBins2];
  Float_t coeferrs[nBins2];

  dPhi1[0]    = -180.0;
  coefs[0]    = AmBpC[index];
  coeferrs[0] = asyer180[index];
  
  dPhi1[1]    = -90.0;
  coefs[1]    = AmC[index];
  coeferrs[1] = asyer90[index];
  
  dPhi1[2]    = 0.0;  
  coefs[2]    = ApBpC[index];
  coeferrs[2] = Sqrt(2./nZero[index]);

  dPhi1[3]    = 90.0;
  coefs[3]    = AmC[index];
  coeferrs[3] = asyer90[index];

  dPhi1[4]    = 180.0;  
  coefs[4]    = AmBpC[index];
  coeferrs[4] = asyer180[index];

  TGraphErrors * fun1Err = new TGraphErrors(nBins2,dPhi1,coefs,0,coeferrs);

  fun1Err->SetFillColor(kBlue-2);
  fun1Err->SetMarkerColor(kBlue);
  fun1Err->SetMarkerStyle(2);
  fun1Err->SetMarkerSize(0.5);
  fun1Err->SetFillStyle(3003);
  fun1Err->SetLineWidth(2);
  fun1Err->SetLineColor(kBlue);
  
  fun1Err->Draw("E1 P same");

  Float_t coefsCr[nBins2];
  Float_t coeferrsCr[nBins2];

  coefsCr[0]    = AmBpCCr[index];
  coeferrsCr[0] = asyerCr180[index];
  
  coefsCr[1]    = AmCCr[index];
  coeferrsCr[1] = asyerCr90[index];
  
  coefsCr[2]    = ApBpCCr[index];
  coeferrsCr[2] = Sqrt(2./nZeroCr[index]);

  coefsCr[3]    = AmCCr[index];
  coeferrsCr[3] = asyerCr90[index];

  coefsCr[4]    = AmBpCCr[index];
  coeferrsCr[4] = asyerCr180[index];
  
  
  TGraphErrors * fun1ErrCr = new TGraphErrors(nBins2,dPhi1,coefsCr,0,coeferrsCr);

  fun1ErrCr->SetFillColor(kRed-2);
  fun1ErrCr->SetMarkerColor(kRed);
  fun1ErrCr->SetLineColor(kRed);
  fun1ErrCr->SetLineWidth(2);
  fun1ErrCr->SetMarkerStyle(20);
  fun1ErrCr->SetMarkerSize(0.5);
  fun1ErrCr->SetFillStyle(3001);
  fun1ErrCr->Draw("E1 P same");
  
  TF1 *fun2 = new TF1("fun2","[0] + [1]*cos(TMath::DegToRad()*x) + [2]*cos(2*TMath::DegToRad()*x)",-180,180);
  fun2->SetParameter(0,ACr[index]);
  fun2->SetParameter(1,BCr[index]);
  fun2->SetParameter(2,CCr[index]);
  fun2->SetLineColor(kRed);
  fun2->SetLineWidth(1);
  fun2->Draw("same");
  
    
  canvas->SaveAs(plotName);

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
