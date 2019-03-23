
#include <cstdlib>


#ifndef _N_MODEDEFN_H_
#define _N_MODEDEFN_H_

namespace neut {
  namespace rew   {

    class NModeDefn
    { 
    public:
  
      NModeDefn();
      //NModeDefn(const NModeDefn & err);
      ~NModeDefn();

      bool isCC(int imode) {
	if ( 1<=abs(imode)&&abs(imode)<=26 ) return true;
	else return false;
      }

      bool isNC(int imode) {
	if ( 31<=abs(imode)&&abs(imode)<=52 ) return true;
	else return false;
      }

      bool isNCEL(int imode) {
	if (abs(imode)==51 || abs(imode)==52) return true;
	else return false;
      }

      bool isCCQE(int imode) {
	if (abs(imode)==1) return true;
	else return false;
      }

      bool isCCRES(int imode) {
	if ( (11<=abs(imode)&&abs(imode)<=13) || abs(imode)==17 || abs(imode)==22 || abs(imode)==23  ) return true;
	else return false;
      }

      bool isCC1PI(int imode) {
	if ( (11<=abs(imode)&&abs(imode)<=13) ) return true;
	else return false;
      }

      bool isCC1PI0(int imode) {
	if ( abs(imode)==12 ) return true;
	else return false;
      }

      bool isCC1PIP(int imode) {
	if ( (11==abs(imode)||abs(imode)==13) ) return true;
	else return false;
      }

      bool isCCCOH(int imode) {
	if ( abs(imode)==16 ) return true;
	else return false;
      }

      bool isNCRES(int imode) {
	if ( (31<=abs(imode)&&abs(imode)<=34) || abs(imode)==38 || abs(imode)==39 || (42<=abs(imode)&&abs(imode)<=45)  ) return true;
	else return false;
      }

      bool isNC1PI(int imode) {
	if ( (31<=abs(imode)&&abs(imode)<=34) ) return true;
	else return false;
      }

      bool isNC1PI0(int imode) {
	if ( (31<=abs(imode)&&abs(imode)<=32) ) return true;
	else return false;
      }

      bool isNCCOH(int imode) {
	if ( abs(imode)==36 ) return true;
	else return false;
      }

      bool isCOH(int imode) {
	if ( abs(imode)==16 || abs(imode)==36 ) return true;
	else return false;
      }

      bool isDIS(int imode) {
	if ( abs(imode)==26 || abs(imode)==46 ) return true;
	else return false;
      }

      bool isMPI(int imode) {
	if ( abs(imode)==21 || abs(imode)==41 ) return true;
	else return false;
      }

      bool isRES(int imode) {
	if ( isCCRES(imode) || isNCRES(imode) ) return true;
	else return false;
      }
      
      bool is1PI(int imode) {
	if ( isCC1PI(imode) || isNC1PI(imode) ) return true;
	else return false;
      }

    };

  } // rew   namespace
} // neut namespace

#endif

