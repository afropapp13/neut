#include "Riostream.h"

void check(){
  
  ifstream in[2];  

  Int_t   IP, ITYPE,i;
  Double_t etbl[38];
  Double_t xsec[2][4][2][7][38];
  Double_t ratio[7];

  in[0].open("spi_nue_xsec_ma1.2.dat");
  in[1].open("../../crsdat/spi_nue_xsec_ma1.2.dat");
  
  TFile *f = new TFile("nue_crs.root","RECREATE");

  TNtuple *nt1 = new TNtuple("nt1","nue_new",
							 "E:Itype:crs1:crs2:crs3:crs4:crs5:crs6:crs7");

  
  Int_t ie,idx,ipar,ip,itype;
  Double_t e,crs1,crs2,crs3,crs4,crs5,crs6,crs7;
  
  Char_t tmpstr[1024];
  in[0].getline(tmpstr,sizeof(tmpstr));
  //  cout << tmpstr << endl;
  in[1].getline(tmpstr,sizeof(tmpstr));
  //  cout << tmpstr << endl;
  ie = 0;
  idx = 0;
  while(idx < 2){
	in[idx] >> ipar >> itype >> e 
	   >> crs1 >> crs2 >> crs3 >> crs4 >> crs5 >> crs6 >> crs7;
	//	cout << ie << " " << ipar << " " << e << " " << itype << endl;
	if (!in[idx].good()){
	  //	  cout << "file1 finished" << endl;
	  idx++;
	  continue;
	  if (idx == 2) break;
	}
	if ( ipar > 0 ){
	  ip = 0;
	}else{
	  ip = 1;
	}

	itype = itype - 1;
	if (e < 0.38){
	  ie = 0;
	}
	//	cout << idx << " " << ip << " " << itype << " " << ie << endl;

	etbl[ie] = e;
	xsec[idx][ip][itype][0][ie] = crs1;
	xsec[idx][ip][itype][1][ie] = crs2;
	xsec[idx][ip][itype][2][ie] = crs3;
	xsec[idx][ip][itype][3][ie] = crs4;
	xsec[idx][ip][itype][4][ie] = crs5;
	xsec[idx][ip][itype][5][ie] = crs6;
	xsec[idx][ip][itype][6][ie] = crs7;
	
	ie = ie + 1;
  }	
  //  for ( idx = 0 ; idx < 2 ; idx++ ){
	for ( ip = 0 ; ip < 2 ; ip++ ){
	  for ( itype = 0 ; itype < 2 ; itype++ ){
		for ( ie = 0 ; ie < 38 ; ie++ ){
		  if (etbl[ie] > 0){
			if (xsec[1][ip][itype][0][ie]>0){
			  for ( i = 0 ; i < 7 ; i++ ){
				ratio[i]=xsec[0][ip][itype][i][ie]/xsec[1][ip][itype][i][ie]*100.-100.;
			  }
			  printf("%3d%3d%8.2f%7.1f%7.1f%7.1f%7.1f%7.1f%7.1f%7.1f\n",
					 ip,itype,etbl[ie],ratio[0],ratio[1],ratio[2],ratio[3],
					 ratio[4],ratio[5],ratio[6]);
			}
		  }
		}
	  }
	}
	//  }
  
}
