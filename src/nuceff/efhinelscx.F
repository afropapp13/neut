*     -----------------------------------------------------
      SUBROUTINE EFHINELSCX(IP,INUC,IPF,INF)      
*     -----------------------------------------------------
*     
*     ( purpose )
*       Charge exchange for high energy in(quasi)elastic scattering
*
*     ( input ) 
C          IP   : INITIAL KIND OF PION
C          IN   : TARGET NUCLEON TYPE
*
*     ( output ) 
C          IPF  : FINAL KIND OF PION
C          INF  : FINAL NUCLEON TYPE
*
c     ( creation & updates)
C          2010.06: P. de Perio
C
C

      INTEGER*4 IP,INUC,IPF,INF

      IF (ABS(IP).eq.211) THEN
         IF (IP.eq.211 .and. INUC.eq.2112) THEN
            INF=2212
            IPF=111
         ELSE IF (IP.eq.-211 .and. INUC.eq.2212) THEN
            INF=2112
            IPF=111
         ELSE
            WRITE(6,*) "Warning in efhinels: Illegal High energy 
     $CX, pi=",IP,"IN=",INUC
            
         END IF
         
      ELSE IF (IP.eq.111) THEN
         IF (INUC.eq.2212) THEN
            IPF=211
            INF=2112
         ELSE IF (INUC.eq.2112) THEN
            IPF=-211
            INF=2212
         ELSE
            WRITE(6,*) "Warning efhinels: Undefined nucleon 
     $for H.E. CX"
         END IF
         
      ELSE
         WRITE(6,*) "Warning efhinels: Undefined pion type in HE.CX" 
         RETURN
         
      END IF

      RETURN
      END

