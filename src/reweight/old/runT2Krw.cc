//////////////////////////////////////////////////
//
// T2K reweighting program
//
//    Author: Patrick de Perio
//
//    Purpose: Tack on event weights to NEUT ntuple
//
//    Usage: ./runT2Krw -i <input ROOT file> -d [sk,nd280numc,piscat] -f nintfsi_pars.dat
//
//    Options:
//
//    	-i : ROOT tree must have the NEWORK/VCWORK and FSIHIST blocks included
//
//    	-d : Select type of simulation
//
//    	-f : Input parameter variation file (must follow format of e.g. 
//    	     $NEUT_ROOT/src/reweight/nintfsi_pars.dat
//
//      -v : Verbosity (0: low, 1: medium, 2: high, 3: fortran output)
//
//////////////////////////////////////////////////

#include <iostream> 
#include <stdlib.h>

#include "t2kneutreweight.h"

using namespace std;

int getArgs(int argc, char* argv[]);
void loadParFile(char *filename);
void checkROOTfile(char *filename);
Int_t getMCID(char *cMCID);

char *infilename ="";
char *neutparfile = "";
char *cMCID = "";
int verbose = 0;

Int_t    nParSets = 0;
Double_t **pars = NULL;

int main(int argc, char *argv[])
{


  // process the arguments
  int args = getArgs(argc, argv);
  if(args != 0){
    cerr << "Usage " << endl;
    return 1;
  }

  Int_t MCID = getMCID(cMCID);

  checkROOTfile(infilename);
  loadParFile(neutparfile);

  TString theFile = infilename;

  t2kneutreweight ana(theFile,MCID,nParSets,pars,verbose);
  
  ana.appendWeightsToTree();

  // Need to delete allocated memory
  //delete pars;
}




int getArgs(int argc, char* argv[]){

  while( (argc > 1) && (argv[1][0] == '-') ){
    switch(argv[1][1]){
      
    case 'i': 
      infilename = argv[2];
      ++argv; --argc;
      break;
    case 'f': 
      neutparfile = argv[2];
      ++argv; --argc;
      break;
    case 'd':
      cMCID = argv[2];
      ++argv; --argc;
      break;
    case 'v':
      verbose = atoi(argv[2]);
      ++argv; --argc;
      break;

    }
    


    ++argv; --argc;
  }

  return 0;

}

void checkROOTfile(char *filename) {
  string str_filename(filename);
  size_t file_extpos = str_filename.find_last_of(".");
  if (file_extpos >= string::npos) {
    cerr << "Malformed input filename: \"" << filename << "\"" << endl;
    exit(-1);
  }
  else {
    string filetype = str_filename.substr(file_extpos+1);

    if (filetype.compare("root")) {
      cerr << "Input file must be a ROOT file: " << filename << endl;
      exit(-1);
    }
  }

  cout << "Loading input ROOT file: " << filename << endl << endl;
}


void loadParFile(char *filename) {
  // Number of variable NEUT parameters 
  // Must change when adding more variable paramters
  const int nPars = 14; 

  // Check FSI reweighting parameters
  const int maxCharacters = 100;  //GUESS FROM LOOKING AT FILE
  char header[maxCharacters];
  FILE *in;
  if (strcmp(filename,"")) {
    in = fopen(filename,"r");
    if (in == NULL) { 
      cerr << endl << filename << " FSI Reweighting file does not exist" << endl << endl;
      exit(-1);
    }
    else {
      // Read file title
      fgets(header, maxCharacters, in);
      if (strcmp(header,"NEUT REWEIGHTING PARAMETER FILE\n")) {
	cerr << "Check FSI Reweighting Parameter File Format: " << header << endl;
	exit(-1);
      }
      
      // Get number of parameter sets
      fscanf(in,"%d", &nParSets);
      //cout << "Number of NEUT parameter sets = " << nParSets << endl;

      // Read column headers (Change this as parameter set definition changes)
      fgets(header, maxCharacters, in);  // next line
      fgets(header, maxCharacters, in);
      if (strcmp(header,"QEMA QEMV QEKAPP SPIMA SPIMV COHMA DISPDF DISBOD FSIQE FSIQEH FSIINEL FSIABS FSICX FSICXH\n")) {
	cerr << "Check FSI Reweighting Parameter File Columns: " << endl;
	cerr << header << endl;
	exit(-1);
      }
     
      
      pars = new Double_t *[nParSets];
      for (int i=0; i<nParSets; i++) {
	pars[i] = new Double_t[nPars];
	for (int ipar=0; ipar<nPars; ipar++) {
	  fscanf(in,"%lf", &pars[i][ipar]);
	  //cout << pars[i][ipar] << " ";
	}
	//cout << endl;
      }
      
      cout << "Loaded NEUT Parameter file: " << filename << endl << endl;  

      cout << "Par. Set " << header << endl;
    }
  } else {

    cerr << "Must specify NEUT parameter file with '-f <filename>'" << endl;
    exit(-1);
  }
  
  
}


int getMCID(char *cMCID) {

  if (!strcmp(cMCID,"sk")) {
    cout << "SK input file specified" << endl << endl;
    return 0;
  }
  else if (!strcmp(cMCID,"nd280numc")) {
    cout << "ND280-numc file specified" << endl << endl;
    return 1;
  }
  else if (!strcmp(cMCID,"nd280anal")) {
    cout << "ND280-anal file specified" << endl << endl;
    return 2;
  }
  else if (!strcmp(cMCID,"piscat")) {
    cout << "Pion scattering file specified" << endl << endl;
    return 3;
  }
  else {
    cerr << "Unknown MCID = '" << cMCID << "'. Must be one of -d {sk, nd280numc, nd280anal, piscat}." << endl;
    exit(-1);
  }
}
