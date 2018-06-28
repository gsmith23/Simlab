#include "TLab.h"
#include "TSim.h"
#include "TTheory.h"
#include "Messages.h"
#include "./includes.h"

Int_t main(int argc, char **argv){ 

  Intro_Message();
  
  if(!Check_Arguments(argc,argv))
    return 1;
  

  //////////////////////////////////////////////////////////////////
  //-----------------------------------------------------------------
  //-----------------------------------------------------------------
  //-----------------------------------------------------------------
  //-----                    THEORY CURVE                       -----

  if( strcmp(argv[1],"6")==0 ){
  
    cout << endl;
    cout << "       -------------------------" << endl; 
    cout << "       | theory curve plotting |" << endl; 
    cout << "       -------------------------" << endl; 
    cout << endl;
  
    TTheory * theory = new TTheory();
    
    Char_t  xVariable = 't';
    
    cout << endl;
    cout << " Graphing Asymmetry " << endl;
    cout << endl;
    cout << " Enter plot type: "       << endl;
    cout << " e - by energy "          << endl;
    cout << " t - by theta (default) " << endl;
    cin  >> xVariable;
    
    if     (xVariable=='T')
      xVariable = 't';
    else if(xVariable=='E')
      xVariable = 'e';
    
    if(xVariable !='e' && 
       xVariable !='t' )
      xVariable = 't';
    
    cout << endl;
    Int_t   nBins = 6;
    Float_t semiSpan = 10;
    
    cout << " Plotting in "             << nBins 
	 << " "                         << (2*semiSpan)
	 << " degrees wide theta bins " << endl;
    
    if     (xVariable=='t')
      cout << " as a function of theta " << endl;
    else if(xVariable=='e')
      cout << " as a function of energy " << endl;
    cout << endl;
    
    theory->GraphFiniteAsymmetry(nBins,
				 semiSpan,
				 xVariable);
    
    
  }
  
  //////////////////////////////////////////////////////////////////
  //-----------------------------------------------------------------
  //-----------------------------------------------------------------
  //-----------------------------------------------------------------
  //-----                 LAB DATA ANALYSIS                     -----
  
  
  if( strcmp(argv[1],"1")==0 || 
      strcmp(argv[1],"0")==0 ||
      strcmp(argv[1],"8")==0 ||
      strcmp(argv[1],"9")==0 ){
    
    cout << endl;
    cout << "       -----------------------" << endl; 
    cout << "       | laboratory analysis |" << endl; 
    cout << "       -----------------------" << endl; 
    cout << endl;
    
    TLab* data;
    data = new TLab(argv[2]);
    
    Char_t overwrite = 'n';
    Char_t option    = 'b';
    
    cout << endl;
    cout << " Checking if Text file exists " << endl; 
      
    if(!(data->RawTextFileExists())){
      cout << endl;
      cout << " ...                                   " << endl;
      cout << " Text File of raw data does not exist. " << endl;
      
      if     (strcmp(argv[1],"9")==0 ||
	      strcmp(argv[1],"0")==0)   {
	
	cout << " Option " << argv[1] 
	     << " requires the raw text file." << endl;
	
	return 1;
      }
    }
    else{
      cout << endl;
      cout << " Text file of raw data does exist   " << endl;
    }
    
    cout << endl;
    cout << " Checking if Raw ROOT file exists " << endl; 
    
    if(!(data->RawROOTFileExists())){
      cout << endl;
      cout << " ...                                   " << endl;
      cout << " ROOT file of raw data does not exist. " << endl;
      cout << " Shall I create one ?                  " << endl;
      cout << " [answer (y/n)]  (default y)           " << endl;
      cout << endl;
      
      Char_t makeRAWROOTFile = 'y';
      cin >> makeRAWROOTFile;
      
      if(makeRAWROOTFile!='n' &&
	 makeRAWROOTFile!='N')
	data->MakeRawDataTreeFile();
    }
    else{
      cout << endl;
      cout << " ROOT file of raw data does exist   " << endl;
      
      if(strcmp(argv[1],"0")==0){
	cout << " Would you like me to overwrite it? " << endl;
	cout << " [answer (y/n)]                     " << endl;
	cout << endl;
	cout << " ";
	cin  >> overwrite;
	
	if( overwrite=='y' || overwrite=='Y'){
	  cout << endl;
	  cout << " Okay, I will overwrite the file. " << endl;
	  data->MakeRawDataTreeFile();
	}
	
      } // end of: if(strcmp(argv[1],"0")=...
    } // end of: else{....
    
    
    if( strcmp(argv[1],"1")==0 || 
	strcmp(argv[1],"0")==0 ){
      
      //////////////////////////////////////////
      // Create calibrated ROOT file
      
      cout << endl;
      cout << "    Calibrate data " << endl; 
	
      cout << endl;
      cout << " Checking for calibrated ROOT file  " << endl; 
      if(!(data->CalibratedROOTFileExists()))
	data->MakeCalibratedDataTreeFile();
      else{
	cout << endl;
	cout << " ROOT file of calibrated data does exist." << endl;
	  
	if(strcmp(argv[1],"0")==0){
	  cout << " Would you like me to overwrite it?      " << endl;
	  cout << " [answer (y/n)]                          " << endl;
	  cout << endl;
	  cout << " ";
	  cin  >> overwrite;
	    
	  if( overwrite=='y' || overwrite=='Y'){
	    cout << endl;
	    cout << " Okay, I will overwrite the file  " << endl;
	    data->MakeCalibratedDataTreeFile();
	  }
	}// end of: if(strcmp(argv[1],"0")...
      }// end of: else{...
	
      cout << endl;
      cout << " Graphing Asymmetry " << endl;
      cout << endl;
      cout << " Enter plot type:   " << endl;
      cout << " b - Theory and Lab (default) " << endl;
      cout << " l - Lab only " << endl;
      cout << " t - Theory only " << endl;
      cout << " ";
      cin  >> option; 
      
      if     (option == 'l' ||
	      option == 'L')
	option = 'l';
      else if(option == 't' ||
	      option == 'T')
	option = 't';
      else if(option == 'b' ||
	      option == 'B')
	option = 'b';
      else{
	cout << endl;
	cout << " invalid choice, setting to default (b) " << endl;
	option = 'b';
      }
      
      if(option!='l' &&
	 option!='t'  ){
	cout << endl;
	cout << " plotting lab results and theory curve " << endl;
	cout << " ..... " << endl;
      }
      else if(option=='l'){
	cout << endl;
	cout << " plotting lab results only " << endl;
	cout << " ..... " << endl;
      }
      else if(option=='t'){
	cout << endl;
	cout << " plotting theory results only " << endl;
	cout << " ..... " << endl;
      }
      
      data->GraphAsymmetry(option);
	
    }// end of: if( strcmp(argv[1],"1")==0 || ......
      
    //delete data;
      
  }// end of: if( strcmp(argv[1],"1")==0 || ......
  
  /////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //-----------------------------------------------------------------
  //-----------------------------------------------------------------
  //-----------------------------------------------------------------
  //-----               SIMULATION DATA ANALYSIS                -----
  //
  else if( strcmp(argv[1],"2")==0 ) {
    cout << endl;
    cout << "     ---------------------------" << endl; 
    cout << "     | simulation lab analysis |" << endl; 
    cout << "     ---------------------------" << endl; 
    cout << endl;
    
    cout << endl;
    cout << " Analysing : " << argv[2] << endl;
    
    // To do extend to two files/argument case
    TSim * simData = new TSim(argv[2]);

    // To do - add check for sorted file/s
    Char_t sort = 'n';
    cout << endl;
    cout << " Sorting ROOT file " << endl;
    cout << " Sort/Re-sort ROOT file ? " << endl;
    cout << " n (default)/ y " << endl;
    cout << " ";
    //!!!!
    cin  >> sort;

    if (sort == 'y' || sort == 'Y' ){
      cout << endl;
      cout << " Sorting data " << endl;
      // sort argv[2]
      simData->SortEvents();
    }
    
    cout << endl;
    cout << " Calculating & Plotting Asymmetry " << endl;

    // second file unpolarised data
    if     (argc==3)
       simData->GraphAsymmetryLab(argv[2],"??");
    else if(argc==4)
      simData->GraphAsymmetryLab(argv[2],argv[3]);
    
    simData->InvestigateAcceptance(argv[2]);
    
    delete simData;
  
  }/////////////////////////////////////////////////////////////////
  else if( strcmp(argv[1],"3")==0 ) {
    
    cout << endl;
    cout << "     --------------------" << endl;
    cout << "    | lab and simulation |" << endl;
    cout << "     --------------------" << endl;
    
    TLab * data = nullptr;
    
    cout << endl;
    if     (argc == 4){
      cout << " Analysing: " << argv[2] << " and " << argv[3] << endl;
      data = new TLab(argv[2],argv[3]);
    }
    else if(argc == 5){
      
      cout << " Analysing: " 
	   << argv[2] << ", " << argv[3] 
	   << " and " << argv[4] << endl;
      data = new TLab(argv[2],argv[3],argv[4]);
    }
    
    Char_t overwrite = 'n';
    Char_t option    = 'b';
    
    cout << endl;
    cout << " Checking if ROOT file exists " << endl; 
    
    if(!(data->RawROOTFileExists())){
      cout << endl;
      cout << " ...                                   " << endl;
      cout << " ROOT file of raw data does not exist. " << endl;
      cout << " Shall I create one ?                  " << endl;
      cout << " [answer (y/n)]  (default y)           " << endl;
      cout << endl;
      
      Char_t makeRAWROOTFile = 'y';
      cin >> makeRAWROOTFile;
      
      if(makeRAWROOTFile!='n' &&
	 makeRAWROOTFile!='N')
	data->MakeRawDataTreeFile();

    }
    else{
      cout << endl;
      cout << " ROOT file of raw data does exist   " << endl;
      
      if(strcmp(argv[1],"0")==0){
	cout << " Would you like me to overwrite it? " << endl;
	cout << " [answer (y/n)]                     " << endl;
	cout << endl;
	cout << " ";
	cin  >> overwrite;
	
	if( overwrite=='y' || overwrite=='Y'){
	  cout << endl;
	  cout << " Okay, I will overwrite the file. " << endl;
	  data->MakeRawDataTreeFile();
	}
	
      } // end of: if(strcmp(argv[1],"0")=...
    } // end of: else{....

    cout << endl;
    cout << " Graphing Asymmetry " << endl;
    cout << endl;
    cout << " Enter plot type:   " << endl;
    cout << " a - Lab, Theory and Simulation (default)" << endl;
    cout << " d - Lab, Theory and Simulation: \n " 
	 <<  "    Divide Lab Asym by Unpol Sim Asym " << endl;
    cout << " s - Lab, Theory and Simulation: \n " 
	 <<  "    Subtract Unpol Sim Asym from Asym " << endl;
    cout << " f - Lab, Theory and Simulation: \n " 
	 <<  "    Divide Lab Asym by f = (sim/theory) " << endl;
    cin  >> option;
    
    if     (option == 'a' || 
	    option == 'A')
      option = 'a';
    else if(option == 'd' ||
	    option == 'D'){
      option = 'd';
      
      if(argc != 5){
	cout << endl;
	cout << " you must use two simulated data files for this option " << endl;
	cout << endl;
	return 1;
      }
    }
    else if(option == 's' ||
	    option == 'S'){
      option = 's';
      
      if(argc != 5){
	cout << endl;
	cout << " you must use two simulated data files for this option " << endl;
	cout << endl;
	return 1;
      }
    }
    else if(option == 'f' ||
	    option == 'F')
    
      option = 'f';
    
    else{
      cout << endl;
      cout << " invalid choice, setting to default (a) " << endl;
      option = 'a';
    }
    
    data->GraphAsymmetry(option);
    
    delete data;
  } 
  else if( strcmp(argv[1],"4")==0 ) {
    
    cout << endl;
    cout << "       ----------------------" << endl; 
    cout << "       | simulation analysis |" << endl;
    cout << "       | for two simulations |" << endl;
    cout << "       ----------------------" << endl; 
    cout << endl;
    
    TSim * simData = NULL;

    cout << endl;
    if     (argc==3){
      cout << " Analysing : " << argv[2]  << endl;
      simData = new TSim(argv[2], argv[2]);
      simData->GraphAsymmetrySim(argv[2], argv[2]);
    }
    else if(argc==4){
      cout << " Analysing : " << argv[2] << " and "<< argv[3] << endl;
      simData = new TSim(argv[2], argv[3]);
      simData->GraphAsymmetrySim(argv[2], argv[3]);
    }
    
    cout << endl;
    cout << " Calculating & Plotting Asymmetry " << endl;
    
    
        
    delete simData;
    
  }
  else if( strcmp(argv[1],"5")==0 ) {
    
    cout << endl;
    cout << "       ----------------------" << endl; 
    cout << "       | simulation analysis |" << endl;
    cout << "       |  scattering study   |" << endl;
    cout << "       ----------------------" << endl; 
    cout << endl;

    cout << endl;
    cout << " Analysing : " << argv[2] << endl;
    
    TSim * simData = new TSim(argv[2]);
    
    cout << endl;
    cout << " Calculating & Plotting Asymmetry " << endl;
    
    Int_t   nThetaBins = 3;
    Float_t thetaSMin  = 40.;
    Float_t thetaSMax  = 70.;
    
    cout << endl;
    cout << "  "
	 << thetaSMin  
	 << " < thetaS < " 
	 << thetaSMax 
	 << endl;
    cout << "  "
	 << nThetaBins 
	 << " theta_{S} bins " 
	 << endl;
    
    simData->GraphAsymmetrySim(argv[2], argv[2],
			       nThetaBins,
			       thetaSMin,
			       thetaSMax);
    
    delete simData;
    
   }
  
  cout << endl;
  cout << " -----------------------------------" << endl;
  cout << " --------------- End  --------------" << endl; 
  cout << " -----------------------------------" << endl;
  
  return 0;  
}
