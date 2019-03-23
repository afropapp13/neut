***********************************************************************
*
      SUBROUTINE RSTCRSETRD(LUNI,IPIN,IEMAX,ETBL,XNUCRS)
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
      REAL*4    ETBL(IEMAX),XNUCRS(13,IEMAX)
      
      INTEGER*4 IP,I33FLG
      INTEGER*4 IPTMP,I33TMP

      INTEGER*4 I,J,IOFST

      REAL*4    ETMP,EDIFF

      DO 10 IP=IPIN,-1*IPIN,-2*IPIN
         DO 20 I33FLG=1,1
            DO 30 I=1,IEMAX
               IOFST = 0
               if (ip.lt.0) then
                  IOFST=10
               endif
               do 40 J=1,3
                  XNUCRS(J+IOFST,I)=0.
 40            continue
               read (LUNI,1000) IPTMP,I33TMP,ETMP,
     $              XNUCRS(1+IOFST,I),XNUCRS(2+IOFST,I),
     $              XNUCRS(3+IOFST,I)
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

 1000 format (I4,I2,F8.2,7E13.6)
            
      return
      end

      
