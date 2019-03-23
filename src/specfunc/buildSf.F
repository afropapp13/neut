* loads SF from grid files and builds an array of the values
* This file contains the subroutines to build these arrays for several nuclei
* All subroutines are called buildSF_chemicalsymbol_atomicmass

***********************************************************************
* buildSF reads in all tables *****************************************
* Currently only carbon, iron, oxygen *********************************
***********************************************************************

      subroutine buildSF      
        implicit none
        integer i_builtSf
        save i_builtSf
        data i_builtSf/0/

        if (i_builtSf .EQ. 0) then
          write (*,*) "Building Spectral function tables"
          call buildSF_C12
          call buildSF_O16
          call buildSF_Fe56
          i_builtSf = 1
        endif
      return
      end


***********************************************************************
*  Carbon - 12  *******************************************************
***********************************************************************
      subroutine buildSF_C12
        implicit none
c include sf arrays, variables etc
c sfArray built as (p, E)
#include "specfunc.h"
#include "neutfilepath.h"

        double precision p, E1, E2, E3, E4
        double precision val1, val2, val3, val4
        integer p_iterator, E_iterator
        integer pRes, ERes
        integer pMin, EMin
        integer pMax, EMax
        integer E_iteratorMax

        CHARACTER*1024 FNAME
        integer*4      LENSTR
        external       lenstr
        integer*4      lenpath
      
        lenpath = lenstr(CRSTBLPATH)
c load file and fill array with values
        FNAME = 'qelSfData/sf/pke12_tot.grid'
        FNAME = CRSTBLPATH(1:lenpath)//FNAME
        open(99, file=FNAME, status='old')
        
        read(99, *) ERes, pRes
        read(99, *) EMin, pMin
        read(99, *) EMax, pMax

c if the last line isn't full, do one more loop over E to catch it
c due to integer division failure
        if ( mod(ERes, 4) .EQ. 0) then
          E_iteratorMax = ERes/4 - 1
        else
          E_iteratorMax = ERes/4
        endif
        
        do p_iterator = 1, pRes
c initialise 1D array as it will only be added to        
          sFArray1D(1,p_iterator) = 0
          read (99, *) p
          p_Array(1,p_iterator) = p
          do E_iterator = 0, E_iteratorMax
            if (E_iterator .EQ. E_iteratorMax .AND. 
     +          E_iteratorMax .EQ. ERes/4) then
              read(99, *) E1, val1, E2, val2
              val3 = 0
              val4 = 0
              E3 = EMax + EMax/ERes
              E4 = E3 + EMax / ERes
            else
              read(99, *) E1, val1, E2, val2, E3, val3, E4, val4
            endif
            if ((p_iterator .EQ. 150).AND.(E_iterator .EQ. 33)) then
              write (*,*) val2
            endif
            sFArray(1, p_iterator, 4*E_iterator + 1) = val1
     +      * p**2
            E_Array(1, 4*E_iterator + 1) = E1
            sFArray(1, p_iterator, 4*E_iterator + 2) = val2
     +      * p**2
            E_Array(1, 4*E_iterator + 2) = E2
            sFArray(1, p_iterator, 4*E_iterator + 3) = val3
     +      * p**2
            E_Array(1, 4*E_iterator + 3) = E3
            sFArray(1, p_iterator, 4*E_iterator + 4) = val4
     +      * p**2
            E_Array(1, 4*E_iterator + 4) = E4
c ----- 1D momentum array            
            sFArray1D(1, p_iterator) = sFArray1D(1, p_iterator) + 
     +      (val1 + val2 + val3 + val4)*p**2
          enddo
        enddo

c set all values outside defined range to zero
        do p_iterator = 1, 200
          if (p_iterator .GT. pRes) then
            p_Array(1, p_iterator) = (pMax / pRes) * p_iterator
            sFArray1D(1, p_iterator) = 0
          endif
          do E_iterator = 1, 200
            if (E_iterator .GE. ERes) then
              sFArray(1, p_iterator, E_iterator) = 0
              E_Array(1, E_iterator) = (EMax / ERes) * E_iterator
            endif
          enddo
        enddo

        nBinsP(1) = pRes
        nBinsE(1) = eRes


        return
        end


***********************************************************************
*  Oxygen - 16  *******************************************************
***********************************************************************
      subroutine buildSF_O16
        implicit none
c include sf arrays, variables etc
c sfArray built as (p, E)
#include "specfunc.h"
#include "neutfilepath.h"

        double precision p, E1, E2, E3, E4
        double precision val1, val2, val3, val4
        integer p_iterator, E_iterator
        integer pRes, ERes
        integer pMin, EMin
        integer pMax, EMax
        integer E_iteratorMax

        CHARACTER*1024 FNAME
        integer*4      LENSTR
        external       lenstr
        integer*4      lenpath
      
        lenpath = lenstr(CRSTBLPATH)
c load file and fill array with values
        FNAME = 'qelSfData/sf/pke16.grid'
        FNAME = CRSTBLPATH(1:lenpath)//FNAME
        open(99, file=FNAME, status='old')
        
        read(99, *) ERes, pRes
        read(99, *) EMin, pMin
        read(99, *) EMax, pMax

c        write (*,*) ERes, ERes/4 - 1

        write (*,*) "reading in values"
        
c if the last line isn't full, do one more loop over E to catch it
c due to integer division failure
        if ( mod(ERes, 4) .EQ. 0) then
          E_iteratorMax = ERes/4 - 1
        else
          E_iteratorMax = ERes/4
        endif
        
        do p_iterator = 1, pRes
          sFArray1D(2 ,p_iterator) = 0
          read (99, *) p
          p_Array(2, p_iterator) = p
          do E_iterator = 0, E_iteratorMax
            if (E_iterator .EQ. E_iteratorMax .AND. 
     +          E_iteratorMax .EQ. ERes/4) then
              read(99, *) E1, val1, E2, val2
              val3 = 0
              val4 = 0
              E3 = EMax + EMax/ERes
              E4 = E3 + EMax / ERes
            else
              read(99, *) E1, val1, E2, val2, E3, val3, E4, val4
            endif
            if ((p_iterator .EQ. 150).AND.(E_iterator .EQ. 33)) then
c              write (*,*) val2
            endif
            sFArray(2, p_iterator, 4*E_iterator + 1) = val1
     +      * p**2
            E_Array(2, 4*E_iterator + 1) = E1
            sFArray(2, p_iterator, 4*E_iterator + 2) = val2
     +      * p**2
            E_Array(2, 4*E_iterator + 2) = E2
            sFArray(2, p_iterator, 4*E_iterator + 3) = val3
     +      * p**2
            E_Array(2, 4*E_iterator + 3) = E3
            sFArray(2, p_iterator, 4*E_iterator + 4) = val4
     +      * p**2
            E_Array(2, 4*E_iterator + 4) = E4
c ----- 1D array - integrated over energy            
            sFArray1D(2, p_iterator) = sFArray1D(2, p_iterator)
     +      + (val1 + val2 + val3 + val4)*p**2
c          write (*,*) p, E_iterator
c            write (*,*) E_iterator, E4, p, val4
          enddo
        enddo

c set all values outside defined range to zero
        do p_iterator = 1, 200
          if (p_iterator .GT. pRes) then
            p_Array(2, p_iterator) = (pMax / pRes) * p_iterator
          endif
          do E_iterator = 1, 200
            if (E_iterator .GE. ERes) then
              sFArray(2, p_iterator, E_iterator) = 0
              E_Array(2, E_iterator) = (EMax / ERes) * E_iterator
            endif
          enddo
        enddo

        nBinsP(2) = pRes
        nBinsE(2) = eRes


        return
        end


***********************************************************************
*  Iron - 56  *********************************************************
***********************************************************************
      subroutine buildSF_Fe56
        implicit none
c include sf arrays, variables etc
c sfArray built as (p, E)
#include "specfunc.h"
#include "neutfilepath.h"

        double precision p, E1, E2, E3, E4
        double precision val1, val2, val3, val4
        integer p_iterator, E_iterator
        integer pRes, ERes
        integer pMin, EMin
        integer pMax, EMax
        integer E_iteratorMax

        CHARACTER*1024 FNAME
        integer*4      LENSTR
        external       lenstr
        integer*4      lenpath
      
        lenpath = lenstr(CRSTBLPATH)
c load file and fill array with values
        FNAME = 'qelSfData/sf/pke56_tot.grid'
        FNAME = CRSTBLPATH(1:lenpath)//FNAME
        open(99, file=FNAME, status='old')
        
        read(99, *) ERes, pRes
        read(99, *) EMin, pMin
        read(99, *) EMax, pMax

c if the last line isn't full, do one more loop over E to catch it
c due to integer division failure
        if ( mod(ERes, 4) .EQ. 0) then
          E_iteratorMax = ERes/4 - 1
        else
          E_iteratorMax = ERes/4
        endif
        
        do p_iterator = 1, pRes
c initialise 1D array as it will only be added to        
          sFArray1D(3, p_iterator) = 0
          read (99, *) p
          p_Array(3, p_iterator) = p
          do E_iterator = 0, E_iteratorMax
            if (E_iterator .EQ. E_iteratorMax .AND. 
     +          E_iteratorMax .EQ. ERes/4) then
              read(99, *) E1, val1, E2, val2
              val3 = 0
              val4 = 0
              E3 = EMax + EMax/ERes
              E4 = E3 + EMax / ERes
            else
              read(99, *) E1, val1, E2, val2, E3, val3, E4, val4
            endif
            if ((p_iterator .EQ. 150).AND.(E_iterator .EQ. 33)) then
              write (*,*) val2
            endif
            sFArray(3, p_iterator, 4*E_iterator + 1) = val1
     +      * p**2
            E_Array(3, 4*E_iterator + 1) = E1
            sFArray(3, p_iterator, 4*E_iterator + 2) = val2
     +      * p**2
            E_Array(3, 4*E_iterator + 2) = E2
            sFArray(3, p_iterator, 4*E_iterator + 3) = val3
     +      * p**2
            E_Array(3, 4*E_iterator + 3) = E3
            sFArray(3, p_iterator, 4*E_iterator + 4) = val4
     +      * p**2
            E_Array(3, 4*E_iterator + 4) = E4
c ----- 1D momentum array            
            sFArray1D(3, p_iterator) = sFArray1D(3, p_iterator) + 
     +      (val1 + val2 + val3 + val4)*p**2
          enddo
        enddo

c set all values outside defined range to zero
        do p_iterator = 1, 200
          if (p_iterator .GT. pRes) then
            p_Array(3, p_iterator) = (pMax / pRes) * p_iterator
            sFArray1D(3, p_iterator) = 0
          endif
          do E_iterator = 1, 200
            if (E_iterator .GE. ERes) then
              sFArray(3, p_iterator, E_iterator) = 0
              E_Array(3, E_iterator) = (EMax / ERes) * E_iterator
            endif
          enddo
        enddo

        nBinsP(3) = pRes
        nBinsE(3) = eRes


        return
        end



