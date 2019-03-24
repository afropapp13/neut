***********************************************************************
*  Common block for spectral functions
*
* sFArray - 2D array of spectral function values in (p, ETilde)
* p_Array - 1D array of momentum values corresponding to sFArray rows
* E_Array - 1D array of ETIlde values corresponding to sFArray columns
* nBinsP - number of bins used in momentum coordinate
* nBinsE - number of bins used in energy coordinate
* ETilde - currently selected value of ETilde
* fermiMomentum - fermi momentum being used for current nucleus
*               - note, the fermiMomentum ONLY controls pauli blocking
*
***********************************************************************
      double precision p_Array(12, 651), SFArray(12, 651)
     & , ETilde, fermiMomentum
      integer nBinsP(12)
  
      common /effspecfunc/ p_Array, SFArray, nBinsP, 
     & ETilde, fermiMomentum
