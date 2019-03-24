#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h> 
#include "HadronTensor.h"


// Reads the Tensor from a file. 

HadronTensor::HadronTensor(){
};

HadronTensor::HadronTensor(std::string filename ){
  Initialize(filename);
}

void
HadronTensor::Initialize(std::string filename ){
  std::string line;
  std::ifstream myfile;
  myfile.open (filename.c_str());

  std::cout << " Reading Tensor File " << filename << std::endl; 

  nbinq3 = nbinq0 = inucleus = icode = -1;
  q3max = q0max = -1.; 

  isPN = false; 

  if (myfile.is_open()){
    while( getline (myfile,line) ) {
      std::cout << line << std::endl; 
      if( !line.compare(0,7,"Nucleus") ) 
	sscanf(line.c_str(),"Nucleus %d ",&inucleus);
      else if( !line.compare(0,4,"isPN") )
	isPN = true;
      else if ( !line.compare(0,5,"Q3Bin") ) 
	sscanf(line.c_str(),"Q3Bin %d ",&nbinq3);
      else if ( !line.compare(0,5,"Q0Bin") ) 
	sscanf(line.c_str(),"Q0Bin %d ",&nbinq0);
      else if ( !line.compare(0,5,"Q3Max") ) {
	sscanf(line.c_str(),"Q3Max  %lf ",&q3max);
      }
      else if ( !line.compare(0,5,"Q0Max") ) 
	sscanf(line.c_str(),"Q0Max %lf ",&q0max);
      else if ( !line.compare(0,4,"Code") )
	sscanf(line.c_str(),"Code %d ",&icode);
      else if( !line.compare(0,5,"-----") ) {
	std::cout << " Reading " << 5 * nbinq0 * nbinq3 << " entries " << std::endl;
	Tensor = new double[5 * nbinq0 * nbinq3];
	int ii = 0;
	double d; 
	while( myfile >> d  ) {
	  Tensor[ii] = d; 
#if 0 
	  if( ii < 10 || 5 * nbinq0 * nbinq3 - ii < 10 ) 
	    std::cout << ii << "  " << Tensor[ii] << std::endl;
#endif 
	  ii++;
	  
	  if( ii > 5 * nbinq0 * nbinq3 ) {
	    std::cout << " Error, too many input in Tensor " << std::endl;
	    exit(1); 
	  }
	}
      }
    }
  }

  std::cout << " Nucleus " << inucleus << std::endl;
  std::cout << " Code    " << icode << std::endl;
  if( isPN ) std::cout << " This is a pn cross-section " << std::endl; 
  std::cout << " Q3 bins " << nbinq3 << " from  0 to  " << q3max <<  std::endl; 
  std::cout << " Q0 bins " << nbinq0 << " from  0 to  " << q0max <<  std::endl; 

  if( nbinq3 > 0 && nbinq0 > 0 ) {
    q3bin = q3max/(double)nbinq3;
    q0bin = q0max/(double)nbinq0;
  }
  else {
    std::cout << " Wrong number of bins " << nbinq3 << "  " << nbinq0 << std::endl; 
  }
    
  myfile.close(); 
}

double HadronTensor::Interpolate(int i, int j, double q3,double q0){

  int ic = -1;

  if( i == 0 && j == 0 ) 
    ic = 0;
  else if ( i == 0 && j == 3 )   
    ic = 1;
  else if ( i == 1 && j == 1 )   
    ic = 2;
  else if ( i == 1 && j == 2 )   
    ic = 3;
  else if ( i == 3 && j == 3 )   
    ic = 4;
  else {
    std::cout << " Wrong tensor code " << i << "  " << j << std::endl;
    return 0.; 
  }

  return Interpolate(ic,q3,q0); 
}


double HadronTensor::htfull(int i, int iq0, int iq3) {

  if( iq3 >= nbinq3 || iq3 < 0.0 ) std::cout << " Error in bin q3 " << iq3 << std::endl; 

  if( iq0 >= nbinq0 || iq0 < 0.0 ) std::cout << " Error in bin q0 " << iq0 << std::endl; 

  if( i >= 5 || i < 0.0 ) std::cout << " Error in tensor index " << i << std::endl;

  int indx = iq3*5*nbinq0 + iq0*5 + i; 

  if( indx > 5 * nbinq0 * nbinq3 ) {
    std::cout << " Error in tensor absolute index " << indx << std::endl;
    return 0; 
  }
  else 
    return Tensor[indx]; 
}

double HadronTensor::Interpolate(int ic, double q3,double q0){
  // Compute the q3 bin.  
  double xiq=q3/q3bin;
  int iq = (int)xiq;
  double sq = xiq-(double) iq; 

  // CHANGE xiq and xiq0
  double xiq0= q0/q0bin; 
  int    iq0 = (int)xiq0;
  double sq0 = xiq0-(double)iq0;

  if( sq0 > 1. || sq > 1. ) std::cout << " ERROR in Interpolation " << sq0 << "  " << sq << std::endl; 

  if( ic < 0 || ic > 4 ) { std::cout << " ERROR " << ic << std::endl; return 0.0; }

  // std::cout << " q3 " << iq << "( " << nbinq3 << " ) q0  " << iq0 << " ( " << nbinq0 << " ) " << std::endl; 
  
  // care with limit cases
  if( (iq >= 0 && iq < nbinq3-1 ) && (iq0 >= 0 && iq0 < nbinq0-1 ) ) {
  
    // interpolating
    double interpolation =  htfull(ic,iq0,iq)*(1.-sq)*(1.-sq0)+
      htfull(ic,iq0,iq+1)*sq*(1.-sq0) +
      htfull(ic,iq0+1,iq)*sq0*(1.-sq) + 
      htfull(ic,iq0+1,iq+1)*sq*sq0;

    return interpolation; 
  }
  else 
    return 0.0;

}
