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
      double precision sFArray(9, 200, 200), sFArray1D(9, 200),
     & p_Array(9, 200), E_Array(9, 200)
     & , ETilde, fermiMomentum
      integer nBinsP(9), nBinsE(9)
      
      common /specfunc/ sFArray, p_Array, E_Array, sfArray1D,
     &  nBinsP, nBinsE, ETilde, fermiMomentum
