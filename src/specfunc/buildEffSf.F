* loads SF from grid files and builds an array of the values
* This file contains the subroutines to build these arrays for several nuclei
* All subroutines are called buildSF_chemicalsymbol_atomicmass

***********************************************************************
* buildSF reads in all tables *****************************************
* Currently only carbon, iron, oxygen *********************************
***********************************************************************

      subroutine buildEFFSF      
        implicit none
        integer i_builtSf
        character*1024 inputfile
        save i_builtSf
        data i_builtSf/0/

        if (i_builtSf .EQ. 0) then
          write (*,*) "Building Effective SF tables"
          inputfile = 'qelSfData/eff_sf/H2_effective_sf.grid'
          call buildEFFSF_NUC(1, inputfile)
          inputfile = 'qelSfData/eff_sf/He3_effective_sf.grid'
          call buildEFFSF_NUC(2, inputfile)
          inputfile = 'qelSfData/eff_sf/He4_effective_sf.grid'
          call buildEFFSF_NUC(3, inputfile)
          inputfile = 'qelSfData/eff_sf/C12_effective_sf.grid'
          call buildEFFSF_NUC(4, inputfile)
          inputfile = 'qelSfData/eff_sf/O16_effective_sf.grid'
          call buildEFFSF_NUC(5, inputfile)
          inputfile = 'qelSfData/eff_sf/Ne20_effective_sf.grid'
          call buildEFFSF_NUC(6, inputfile)
          inputfile = 'qelSfData/eff_sf/Al27_effective_sf.grid'
          call buildEFFSF_NUC(7, inputfile)
          inputfile = 'qelSfData/eff_sf/Ar40_effective_sf.grid'
          call buildEFFSF_NUC(8, inputfile)
          inputfile = 'qelSfData/eff_sf/Fe56_effective_sf.grid'
          call buildEFFSF_NUC(9, inputfile)
          inputfile = 'qelSfData/eff_sf/Cu63_effective_sf.grid'
          call buildEFFSF_NUC(10, inputfile)
          inputfile = 'qelSfData/eff_sf/Zn64_effective_sf.grid'
          call buildEFFSF_NUC(11, inputfile)
          inputfile = 'qelSfData/eff_sf/Pb_effective_sf.grid'
          call buildEFFSF_NUC(12, inputfile)
          i_builtSf = 1
        endif
      return
      end

      subroutine buildEFFSF_NUC(nuc, inputFile)
        implicit none
        character*1024 inputFile
        integer nuc       
#include "effspecfunc.h"
#include "neutfilepath.h"
        double precision pval, thisp
        integer p_iterator
        integer pRes
        integer pMin
        integer pMax

        integer*4      LENSTR
        external       lenstr
        integer*4      lenpath
      
        lenpath = lenstr(CRSTBLPATH)
c load file and fill array with values
        inputfile = CRSTBLPATH(1:lenpath)//inputfile

#ifdef __GFORTRAN__
        open(99, file=inputFile, status='old', action='read')        
#else
        open(99, file=inputFile, status='old')
#endif        
        read(99, *) pRes
        read(99, *) pMin
        read(99, *) pMax

        do p_iterator = 1, pRes
          read (99, *) thisp, pval
          p_Array(nuc, p_iterator) = thisp
          SFArray(nuc, p_iterator) = pVal
       enddo

c set all values outside defined range to zero
        do p_iterator = 1, 651
          if (p_iterator .GT. pRes) then
             p_Array(nuc, p_iterator) = (pMax/pRes)*p_iterator
             SFArray(nuc, p_iterator) = 0
          endif
        enddo

        nBinsP(nuc) = pRes

        return
        end

