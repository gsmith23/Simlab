#include "Messages.h" 
#include "iostream"

using namespace std;

bool Check_Arguments(int argc, char **argv){
  
  if( ( argc < 2 ) || ( argc > 6 ) ||
      ( strcmp(argv[1],"0")==0 && argc != 3)     ||
      ( strcmp(argv[1],"1")==0 && argc != 3)     ||
      (
       ( strcmp(argv[1],"2")==0 && argc != 3) &&
       ( strcmp(argv[1],"2")==0 && argc != 4) 
       )                                         ||
      (
       ( strcmp(argv[1],"3")==0 && argc != 4) &&
       ( strcmp(argv[1],"3")==0 && argc != 5) 
       )                                         ||
      (
       ( strcmp(argv[1],"4")==0 && argc != 3)    &&
       ( strcmp(argv[1],"4")==0 && argc != 4)    
       )                                         ||
      ( strcmp(argv[1],"5")==0 && argc != 3)     ||
      ( strcmp(argv[1],"6")==0 && argc != 2)     ||
      ( strcmp(argv[1],"8")==0 && argc != 5)     ||
      ( strcmp(argv[1],"9")==0 && argc != 3) 
      ){
    Argument_Options_Message();
    return false;
  }
  
  return true;
}

void Intro_Message(){

  cout << endl;
  cout << endl;
  cout << endl;
  cout << " -----------------------------------" << endl;
  cout << " ---------  QET Analysis -----------" << endl; 
  cout << " -----------------------------------" << endl;
  cout << " ---------- Version 7.0 ------------" << endl; 
  cout << " -----------------------------------" << endl;
  cout << " ---- Lab and/or Simulated data ----" << endl;
  cout << " --------  and/or Theory  ----------" << endl; 
  cout << " -----------------------------------" << endl;
  cout << " ------ Ten Crystal Analysis -------" << endl;
  cout << " -----------------------------------" << endl;

}

void Argument_Options_Message(){
  
  cout << endl;
  cout << " ------------------------------------- " << endl;
  cout << " First argument: Option                " << endl; 
  cout << " Options:                              " << endl; 
  cout << " 1 - lab data analysis                 " << endl; 
  cout << " 2 - sim data analysis                 " << endl; 
  cout << " 3 - lab and sim data analysis         " << endl;
  cout << " 4 - exact angles sim data analysis    " << endl;
  cout << " 5 - exact angles scattered analysis   " << endl;
  cout << " 6 - plot theory curve only            " << endl; 
  cout << " 0 - lab data analysis (overwrite)     " << endl; 
  cout << " 9 - lab data - make raw trees only    " << endl; 
  cout <<                                              endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Option 1 - one further argument       " << endl;
  cout << "        raw lab run number             " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Examples:                             " << endl; 
  cout << " ./simLab 1 1470                       " << endl; 
  cout << " ------------------------------------- " << endl;
  cout <<                                              endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Option 2 - one/two further arguments  " << endl; 
  cout << " simulated data file number, or        " << endl;
  cout << " entangled/polarised then unpolarised  " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Examples:                             " << endl; 
  cout << " ./simLab 2 001132000                  " << endl;
  cout << " ./simLab 2 001132000 000132000        " << endl;
  cout << " ------------------------------------- " << endl;
  cout <<                                              endl;
  cout << " ------------------------------------- " << endl;
  cout << " Option 3 - two/three further args.    " << endl; 
  cout << "    lab then sim run numbers           " << endl;
  cout <<                                              endl;
  cout << "   lab, tangled sim then unpol sim     " << endl;
  cout << "  (lab & sim data normalised to unpol) " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Examples:                             " << endl; 
  cout << " ./simLab 3 1470 001132000             " << endl;
  cout << " ./simLab 3 1470 001132000 000132000   " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Option 4 - two further arguments      " << endl;
  cout << " simulated data file number (entangled)" << endl;
  cout << " simulated data file number (polarised)" << endl;
  cout << " OR, for single file analysis..        " << endl; 
  cout << " simulated data file number            " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Examples:                             " << endl; 
  cout << " ./simLab 4 001132000 000132000;       " << endl;
  cout << " ./simLab 4 001132000                  " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Option 5 - two further arguments      " << endl;
  cout << " simulated data file number            " << endl;
  cout << " same simulated data file number       " << endl;
  cout << " ( modify simLab.cc to set scattering  " << endl;
  cout << "  angles and number of  bins)          " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Examples:                             " << endl; 
  cout << " ./simLab 5 001132000 001132000;       " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Option 0 - one further argument       " << endl;
  cout << "        raw lab run number             " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Examples:                             " << endl; 
  cout << " ./simLab 0 1470;                      " << endl; 
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Option 9 - one further argument       " << endl;
  cout << "        raw lab run number             " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " ------------------------------------- " << endl;
  cout << " Examples:                             " << endl; 
  cout << " ./simLab 9 1470;                      " << endl; 
  cout << " ------------------------------------- " << endl;
  cout <<                                              endl;
  cout <<                                              endl;
  cout << endl;

}
