#include <iostream> 
#include <stdlib.h>
#include <string.h>

#include "piscatana.h"

using namespace std;

int getArgs(int argc, char* argv[]);

string infilename = "";
string legNames="";
int ExpNumEvents=-1;
int drawType=1;

int main(int argc, char *argv[])
{

  // process the arguments
  int args = getArgs(argc, argv);
  if(args != 0){
    cerr << "Usage " << endl;
    return 1;
  }
  
  piscatana ana(infilename,legNames,ExpNumEvents);

  gROOT->ProcessLine(".x ~/.rootlogon.C");  

  if (drawType==1)
    ana.plot_sep_pion_nuc_iint();
  else if (drawType==2)
    ana.plot_sep_nuc();
  else if (drawType==3)
    ana.plot_sep_pion_iint();
    
}

int getArgs(int argc, char* argv[]){

  while( (argc > 1) && (argv[1][0] == '-') ){
    switch(argv[1][1]){

    // Get ROOT input root file
    case 'i': 
      infilename = argv[2];
      ++argv; --argc;
      break;

    case 'l': 
      legNames = argv[2];
      ++argv; --argc;
      break;

    case 'e': 
      ExpNumEvents = atoi(argv[2]);
      ++argv; --argc;
      break;

    case 'd': 
      drawType = atoi(argv[2]);
      ++argv; --argc;
      break;

    }
    ++argv; --argc;
  }

  return 0;

}
