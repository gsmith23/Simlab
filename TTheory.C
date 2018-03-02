#include "TTheory.h"
//#include "./includes.h"

#if !defined(__CINT__)
ClassImp(TTheory)
#endif

// Default constructor
TTheory::TTheory(){
}

// Destructor
TTheory::~TTheory(){
}

//----------------------------------------------
// Finite Geometry Asymmetry Calculation
// Snyder 1948 
// Angular Correlation of Scattered Annihilation Radiation
// http://link.aps.org/doi/10.1103/PhysRev.73.440

// Member functions below are used to calculate
// rho = N(90)/N(0)

Float_t TTheory::X(Float_t theta){
  return ( 2. - Cos(theta) );
}

Float_t TTheory::JLim(Float_t theta){
  return ( Log(X(theta)) - 1./( 2.*X(theta)*X(theta) ) );
}

Float_t TTheory::J(Float_t theta, Float_t semiSpan){
  return ( JLim(theta+semiSpan) - JLim(theta-semiSpan) ) ;
}

Float_t TTheory::JdashLim(Float_t theta){
  return ( -X(theta) + 4.*Log(X(theta)) + 3./( X(theta) ) );  
}

Float_t TTheory::Jdash(Float_t theta, Float_t semiSpan){
  return ( JdashLim(theta+semiSpan) - JdashLim(theta-semiSpan) ) ;
}

// Asymmetry for finite theta only
// semi-span is half the size of the detector in theta
Float_t TTheory::rho1(Float_t theta, Float_t semiSpan){
  return ( 1. + 1./( 1./2.*( J(theta,semiSpan)/Jdash(theta,semiSpan) )*( J(theta,semiSpan)/Jdash(theta,semiSpan) ) - ( J(theta,semiSpan)/Jdash(theta,semiSpan) ) ) );
}

Float_t TTheory::u(Float_t alpha){
  return ( 2*alpha*alpha - 1./2.*Sin(2*alpha)*Sin(2*alpha) );
}

Float_t TTheory::w(Float_t alpha){
  return ( 2*alpha*alpha + 1./2.*Sin(2*alpha)*Sin(2*alpha) );
}

Float_t TTheory::Z(Float_t alpha){
  return (u(alpha)/w(alpha));
}

// Asymmetry for finite theta and finite phi
// alpha is half the size of the detector in phi
Float_t TTheory::rho2(Float_t theta, Float_t semiSpan, Float_t alpha){
  return ( (Z(alpha) + rho1(theta,semiSpan) )/( 1 + Z(alpha)*rho1(theta,semiSpan)) ); 
}  

//working stage - no resolution
Float_t TTheory::modFactor(Float_t theta1, Float_t theta2){
  Float_t m1 = Power(Sin(theta1),2.)*(2. - Cos(theta1))/(2. + Power((1.- Cos(theta1)),3.));
  Float_t m2 = Power(Sin(theta2),2.)*(2. - Cos(theta2))/(2. + Power((1.- Cos(theta2)),3.));
  return m1*m2;
}
  
//--------------------------------------------------

// Compton scattering variable conversions
Float_t TTheory::ThetaToPhotonEnergy(Float_t theta){
  return (511./(2 - Cos(TMath::DegToRad()*theta)));
}

Float_t TTheory::ThetaToElectronEnergy(Float_t theta){
  return (511. - (511./(2. - Cos(TMath::DegToRad()*theta))));
}

//--------------------------------------------------
// Plotting routine 

void TTheory::GraphFiniteAsymmetry(Int_t   nBins,
				   Float_t semiSpan,
				   Char_t  xVariable){
  
  TCanvas *canvas = new TCanvas("canvas","canvas",
				10,10,1200,800);
  
  // This is a dummy histogram which is 
  // used for setting the axis on the TGraph
  TH1F * hr = new TH1F();
  
  // Plot as a function of energy or theta
  if     ( xVariable=='t') 
    hr = canvas->DrawFrame(0.0,0.5,180.0,3.0);
  else if( xVariable=='e'){ 
    hr = canvas->DrawFrame(0.0,0.5, 511.,3.0);
  }

  TGraphErrors * gr[3];

  // this second graph array is used for simultaneously
  // plotting as a function of photon energy when plotting
  // as a function of electron energy (energy deposited
  // in 1st crystal)
  TGraphErrors * gr2[3];
  
  Float_t theta[nBins];
  
  Float_t elecNRG[nBins];
  Float_t gammaNRG[nBins];
  
  Float_t asymm1[nBins];
  Float_t asymm2[nBins];
  Float_t asymm3[nBins];
  
  // To Do: determine resolution
  // of detector system using simulation
  Float_t alpha1   = 20.;
  Float_t alpha2   = 30.;
  Float_t alpha3   = 40.;
  
  TLegend * leg =  new TLegend(0.6,0.7,0.8,0.85);
  
  TString legTit;
  legTit.Form("semi span #theta  = %.1f ^{o}",
	      semiSpan);
  
  leg->AddEntry((TObject*)0,legTit, "");

  TString legStr[3];
  
  legStr[0].Form("#alpha_{#phi} = %.1f ^{o}",
		 alpha1);
  legStr[1].Form("#alpha_{#phi} = %.1f ^{o}",
		 alpha2);
  legStr[2].Form("#alpha_{#phi} = %.1f ^{o}",
		 alpha3);
  
  
  alpha1 = alpha1 * DegToRad();
  alpha2 = alpha2 * DegToRad();
  alpha3 = alpha3 * DegToRad();
  
  semiSpan = DegToRad()*semiSpan;
  
  Float_t centreTh = 94.;
  
  for(Int_t i = 0 ; i < nBins ; i++){
    
    // symmetric bins around centreTh
    theta[i]    = centreTh - (nBins-1)*semiSpan*RadToDeg();
    theta[i]    = theta[i] +  i*2*semiSpan*RadToDeg();
    elecNRG[i]  = ThetaToElectronEnergy(theta[i]);
    gammaNRG[i] = ThetaToPhotonEnergy(theta[i]);
    
    theta[i]  = theta[i]*DegToRad();
    
    asymm1[i] = rho2(theta[i],semiSpan,alpha1);
    asymm2[i] = rho2(theta[i],semiSpan,alpha2);
    asymm3[i] = rho2(theta[i],semiSpan,alpha3);
    
    theta[i] = theta[i]*RadToDeg();
  }


  if     ( xVariable=='t'){
    hr->GetXaxis()->SetTitle("#theta (deg)");
    gr[0] = new TGraphErrors(nBins,theta,asymm1,0,0);
    gr[1] = new TGraphErrors(nBins,theta,asymm2,0,0);
    gr[2] = new TGraphErrors(nBins,theta,asymm3,0,0);    
  }
  else if( xVariable=='e'){
    hr->GetXaxis()->SetTitle("energy (keV): electron (dashed), photon (solid)");
    gr[0]  = new TGraphErrors(nBins,elecNRG,asymm1,0,0);
    gr[1]  = new TGraphErrors(nBins,elecNRG,asymm2,0,0);
    gr[2]  = new TGraphErrors(nBins,elecNRG,asymm3,0,0);    
    gr2[0] = new TGraphErrors(nBins,gammaNRG,asymm1,0,0);
    gr2[1] = new TGraphErrors(nBins,gammaNRG,asymm2,0,0);
    gr2[2] = new TGraphErrors(nBins,gammaNRG,asymm3,0,0);    
  }

  gr[0]->SetLineColor(kRed+1);
  gr[0]->SetMarkerColor(kRed+1);
  gr[0]->SetMarkerStyle(20);
  gr[0]->SetLineStyle(2);
  
  if( xVariable=='e'){  
    gr2[0]->SetLineColor(kRed-2);
    gr2[0]->SetMarkerColor(kRed-2);
    gr2[0]->SetMarkerStyle(20);
  }
  
  gr[1]->SetLineColor(kBlue+1);
  gr[1]->SetMarkerColor(kBlue+1);
  gr[1]->SetMarkerStyle(20);
  gr[1]->SetLineStyle(2);
  
  if( xVariable=='e'){
    gr2[1]->SetLineColor(kBlue-2);
    gr2[1]->SetMarkerColor(kBlue-2);
    gr2[1]->SetMarkerStyle(20);
  }

  gr[2]->SetLineColor(kGreen+3);
  gr[2]->SetMarkerColor(kGreen+3);
  gr[2]->SetMarkerStyle(20);
  gr[2]->SetLineStyle(2);
  
  if( xVariable=='e'){
    gr2[2]->SetLineColor(kGreen+1);
    gr2[2]->SetMarkerColor(kGreen+1);
    gr2[2]->SetMarkerStyle(20);
  }
  
  TString plotStyle = "PL";
  
  gr[0]->Draw(plotStyle);
  
  plotStyle = plotStyle + "same";
  
  gr[1]->Draw(plotStyle);
  gr[2]->Draw(plotStyle);

  if( xVariable=='e'){
    gr2[0]->Draw(plotStyle);
    gr2[1]->Draw(plotStyle);
    gr2[2]->Draw(plotStyle);
  }
  
  hr->GetYaxis()->SetTitle("P(#Delta#phi = 90)/P(#Delta#phi = 0)");
  
  Char_t plotName[128]; 
    
  if     ( xVariable=='t')
    sprintf(plotName,"../Plots/A_Theory_%d_bins_theta.pdf", nBins);
  else if( xVariable=='e')
    sprintf(plotName,"../Plots/A_Theory_%d_bins_energy.pdf", nBins);
  

  leg->AddEntry(gr[0],legStr[0],plotStyle);
  leg->AddEntry(gr[1],legStr[1],plotStyle);
  leg->AddEntry(gr[2],legStr[2],plotStyle);
  
  leg->Draw();

  canvas->SaveAs(plotName);
  
}

//--------------------------------------------------
