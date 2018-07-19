#include "Compton.h" 
#include "TMath.h" 

using namespace TMath;

float ElectronEnergyToTheta(float energy){
  return ACos(2. - (511./(511. - energy)))*RadToDeg();
}

float PhotonEnergyToTheta(float energy){
  return ACos(2. - (511./energy))* RadToDeg();;
}

float ThetaToPhotonEnergy(float theta){
  return (511./(2 - Cos(DegToRad()*theta)));
}

float ThetaToElectronEnergy(float theta){
  return (511. - (511./(2. - Cos(DegToRad()*theta))));
}
