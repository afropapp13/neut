***********************************************************************
*
      SUBROUTINE RSTCRSRD(LUNI,IPIN,IEMAX,ETBL,XNUCRS)
*
*     ( purpose )
*       read cross-section table from the file (LUNI)
*
*     ( input )
*       LUNI        : Logical Unit number for the cross-section table
*       IPIN        : Particle code (12,14,16)
*       IEMAX       : Number of energy bins (Number of lines)
*       ETBL(IEMAX) : Energy table ( reference )
*
*     ( output )
*       XNUCRS(17,2,IEMAX) : Cross-section table
*
*     ( creation date and author )
*       2010.Oct. ; use file for the cross-section table by Y.Hayato
*      
*     ( comment )
*       for single pi prod.
*
      IMPLICIT NONE

      INTEGER*4 LUNI,IEMAX,IPIN
      REAL*4    ETBL(IEMAX),XNUCRS(17,2,IEMAX)
      
      INTEGER*4 IP,I33FLG
      INTEGER*4 IPTMP,I33TMP

      INTEGER*4 I,J,IOFST

      REAL*4    ETMP,EDIFF

C     IPIN        : Particle code (12,14,16)
      DO 10 IP=IPIN,-1*IPIN,-2*IPIN
C     Read the total and delta cross-section
         DO 20 I33FLG=1,2
C     Read the full energy table
            DO 30 I=1,IEMAX
               IOFST = 0
C     For anti-neutrino (ip < 0) all the modes are offset by 10
               if (ip.lt.0) then
                  IOFST=10
               endif
C     Read the 7 different interaction modes (3 CC, 4 NC)
               do 40 J=1,7
                  XNUCRS(J+IOFST,I33FLG,I)=0.
 40            continue
C     Read the pre-calculated value and store it in XNUCRS
               read (LUNI,1000) IPTMP,I33TMP,ETMP,
     $              XNUCRS(1+IOFST,I33FLG,I),XNUCRS(2+IOFST,I33FLG,I),
     $              XNUCRS(3+IOFST,I33FLG,I),XNUCRS(4+IOFST,I33FLG,I),
     $              XNUCRS(5+IOFST,I33FLG,I),XNUCRS(6+IOFST,I33FLG,I),
     $              XNUCRS(7+IOFST,I33FLG,I)
C     The different in energy between the requrested energy and the
C     stored energy in the table (which the table was calculated with)
               EDIFF = abs(ETMP-ETBL(I))
               if ((IPTMP.ne.IP).or.(I33TMP.ne.I33FLG).
     $              or.(EDIFF.gt.0.01)) 
     $              then
                  write(*,*) 'Corrupted? Expected IP,I33FLG,E=',
     $                 IP,',',I33FLG,',',ETBL(I),
     $                 ' but IP,I33FLG,E=',IPTMP,',',I33TMP,',',ETMP
                  stop
               endif
 30         continue
 20      continue
 10   continue

C     This specifies the format of the input file
C     Integer of 4 lenght, integer of 2 length, float of 8
C     lenght with 2 decimal places, 7 floats of 13 length with
C     6 decimal places
 1000 format (I4,I2,F8.2,7E13.6)
            
      return
      end

      
