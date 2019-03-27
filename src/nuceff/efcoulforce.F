*     -----------------------------------------------------
      SUBROUTINE EFCOULFORCE(P, IP, X, R, dt, iNuc)
*     -----------------------------------------------------
*
*     ( purpose )
C        Apply momentum change due to Coulomb force of
C        nucleus (called from eftrace.F)
*
*     ( input )
C         P   : Momentum vector (MeV/c)
C         IP  : PDG Particle ID
C         X   : Position vector in nucleus (fm)
C         R   : Radial position in nucleus (fm)
C         dt  : Step time (seconds * speed of light)
*
*     ( output )
C         P   : New momentum (MeV/c)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C#define DEBUG

      implicit none

#include <efcoulforce.h>

      real P(3), X(3), R, dt
      integer IP, charge_sign, iNuc, iRad, i

C     Step size and max radius defined in src/neututils/calc_nuc_dens/neutnuclei.h
      real r_step, max_r
      parameter (r_step=0.1)
      parameter (max_r=14.9)

C     For calculating Coulomb constant
      real intDensSlope, coulombConstant


C     Only apply to charged pions for now
      IF (ABS(IP).NE.211) THEN
         RETURN
      END IF

C      write(*,*) "Pos = ", X(1), X(2), X(3)
C      write(*,*) "Mom before = ", P(1), P(2), P(3)

      charge_sign = IP/abs(IP)

C      write(*,*) ip, charge_sign
C      write(*,*) iNuc, nucZ(iNuc)

C     Get iRad index

C     Use full spherical volume if given radius is too large
C     (This should not occur if the table was defined up to large enough r)
      if (R.ge.max_r) then
         iRad = max_r/r_step+1
         coulombConstant = integratedDensity(iRad,iNuc)
         write(*,*) "EFCOULFORCE Warning: R >= max_r"

C     Linear interpolation of table
      else
         iRad = R/r_step+1

         intDensSlope = ( integratedDensity(iRad+1,iNuc) -
     &                    integratedDensity(iRad,iNuc) ) / r_step

         coulombConstant = intDensSlope * (R-radialPosition(iRad)) +
     &                     integratedDensity(iRad,iNuc)

C         write(*,*) R, r_step, iRad, integratedDensity(iRad+1,iNuc),
C     &  integratedDensity(iRad,iNuc), intDensSlope, radialPosition(iRad)
C     &  ,coulombConstant

      end if

      coulombConstant = coulombConstant / R**2

      do i = 1, 3
         P(i) = P(i) + charge_sign*coulombConstant*X(i)/R
      end do

C      write(*,*) "Mom after = ", P(1), P(2), P(3)


      END


*     -----------------------------------------------------
      SUBROUTINE loadIntDensTable()
*     -----------------------------------------------------

      implicit none

#include <efcoulforce.h>

      integer iNuc, iRad

      character*28 filename
      logical first /.true./

      nNuclei = nNucleiPar
      nRadPoints = nRadPointsPar

C     Load integrated density table (normalized with Coulomb factor)
      filename='nucdens_int_coulomb.dat'
      if (first) then
         write (*,*) "Loading integrated density table: ",filename
         open (unit=10,file=filename,status='old')
         read (10,*) nucZ

#ifdef DEBUG
         do iNuc = 1, nNuclei
            write(*,*) nucZ(iNuc)
         end do
#endif

         do iRad = 1, nRadPoints
            read (10,*) radialPosition(iRad), (integratedDensity(iRad,iNuc),iNuc=1,nNuclei)

#ifdef DEBUG
            write (*,*) 'RADIUS #', iRad
            write (*,*) radialPosition(iRad)

            do iNuc = 1, nNuclei
               write (*,*) integratedDensity(iRad,iNuc)
            end do
#endif

         end do

         close (10)
         first = .false.

      end if


      END
