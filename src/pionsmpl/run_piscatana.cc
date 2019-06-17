#include <iostream> 
#include <stdlib.h>
#include <string.h>

#define MAIN
#include "pars_piscat.h"
#include "piscatana.h"

using namespace std;

int getArgs(int argc, char* argv[]);

string infilename = "";
string outfilename = "";

int main(int argc, char *argv[])
{

  // process the arguments
  int args = getArgs(argc, argv);
  if(args != 0){
    cerr << "Usage " << endl;
    return 1;
  }
  
  piscatana ana(infilename,outfilename);
  ana.Loop();
}

int getArgs(int argc, char* argv[]){

  while( (argc > 1) && (argv[1][0] == '-') ){
    switch(argv[1][1]){

    // Get ROOT input root file
    case 'i': 
      infilename = argv[2];
      ++argv; --argc;
      break;

    // Get ROOT output root file
    case 'o': 
      outfilename = argv[2];
      ++argv; --argc;
      break;

    }
    
    ++argv; --argc;
  }

  return 0;

}
