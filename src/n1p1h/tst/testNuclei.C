#include <iostream>
#include <vector>
#include <math.h>
#include "f2cdeclarations.h" 
#include "Nucleus.h" 

int main(int argc,char **argv ){

 
  if( argc < 4 ) {
    std::cout << " Command line format:   " << argv[0] << "  cost Enu  tmu" << std::endl; 
    exit (0); 
  }

  double cost = atof(argv[1]); 
  double eingev  = atof(argv[2]); 
  double tmugev  = atof(argv[3]); 
  int NRC = 0; 
  int NPC = 0;
  double R; 
  double dP; 
  double vc; 

  double xmag     = 1.05; 
  int    irpa     = 1; //1->NUCLEAR EFFECT
  int    anti     = 1;//1->neutrino;-1->antineutrinos
  double type     = 12;//12,16,40.;
  int contrlepton = 1;//1->muons,0->electron
  double rm; 

  initialization_(&xmag,&irpa,&rm,&anti,&type,&contrlepton);

  for(int i = 0; i < 10; i++ ) 
    std::cout << i << "   " << vcfto_.vcd[i] << std::endl; 

  std::cout << " Orig C neutrino " << std::endl; 
  std::cout << sigthnacho_(&eingev,&tmugev,&cost) << std::endl; 

  int ipol = 1; 

  constantsinitialization_(&xmag,&ipol);

  Nucleus *C = new Nucleus(12); 
  Nucleus *O = new Nucleus(16); 
  Nucleus *Fe = new Nucleus(56);
  Nucleus *Ar = new Nucleus(40);
  Nucleus *Pb = new Nucleus(208);
  Nucleus *Al = new Nucleus(27); 
  Nucleus *Si = new Nucleus(28);

  std::cout << " Nuclei defined " << std::endl; 

  int lepton = 1;
  SelectLepton(lepton); 



  C->CopyNucleiStructure(1); 

  for(int i = 0; i < 10; i++ ) 
    std::cout << i << "   " << vcfto_.vcd[i] << std::endl; 

  std::cout << " C neutrino " << sigthnacho_(&eingev,&tmugev,&cost) << std::endl; 


  C->CopyNucleiStructure(-1); 

  for(int i = 0; i < 10; i++ ) 
    std::cout << i << "   " << vcfto_.vcd[i] << std::endl; 

  
  std::cout << " C antineutrino " << sigthnacho_(&eingev,&tmugev,&cost) << std::endl; 

  O->CopyNucleiStructure(1); 

  for(int i = 0; i < 10; i++ ) 
    std::cout << i << "   " << vcfto_.vcd[i] << std::endl; 


  std::cout << " O neutrino " << sigthnacho_(&eingev,&tmugev,&cost) << std::endl;  

  O->CopyNucleiStructure(-1); 
  
  std::cout << " O antineutrino " << sigthbruno_(&eingev,&tmugev,&cost,&NRC,&NPC,&R,&dP,&vc)<< std::endl;  

   C->CopyNucleiStructure(-1); 
  
  std::cout << " C antineutrino " << sigthnacho_(&eingev,&tmugev,&cost) << std::endl; 

} 
