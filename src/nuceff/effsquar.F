*     -----------------------------------------------------
      FUNCTION EFFSQUAR(MODE,SK,SKD,GK)
*     -----------------------------------------------------
*
*     ( purpose )
*        CALCULATE !F!**2  -- SCATTERING RATE
*
*     ( input )
*      MODE :   INTERACTION MODE
*            =1    PI+,P ---> PI+,P
*            =2    PI+,N ---> PI+,N
*            =3    PI+,N ---> PI0,P
*            =4    PI0,P ---> PI0,P
*            =5    PI0,N ---> PI0,N
*            =6    PI0,P ---> PI+,N
*            =7    PI0,N ---> PI-,P
*            =8    PI-,P ---> PI-,P
*            =9    PI-,N ---> PI-,N
*            =10   PI-,P ---> PI0,N
*       SK   :   INITIAL PI-N RELATIVE MOMENTUM IN C-M SYSTEM(MEV)
*       SKD  :   FINAL   PI-N RELATIVE MOMENTUM IN C-M SYSTEM(MEV)
*       GK   :   PI-N TOTAL MOMENTUM IN LAB. SYSTEM(MEV)
*
*
*     ( output )
*       EFFSQUAR :   !F!**2
*
*     ( update )
*       2010.05    ; P.de Perio add CX scaling factor

#include "neutparams.h"

      DIMENSION SK(3),SKD(3),GK(3)
      COMPLEX FS1,FS3,FP1,FP3,GP1,GP3,FD1,GD1,CTEMP1,CTEMP2
      COMMON/PSHIFT/FS1,FS3,FP1,FP3,GP1,GP3,FD1,GD1
C
      EFFSQUAR=0.
      SKABS=SQRT(SK(1)**2+SK(2)**2+SK(3)**2)
      SKDABS=SQRT(SKD(1)**2+SKD(2)**2+SKD(3)**2)
      SKSKDI=SK(1)*SKD(1)+SK(2)*SKD(2)+SK(3)*SKD(3)
      COST=SKSKDI/SKABS/SKDABS
      SSINT=1.-COST**2
      IF(SSINT.LE.0.)GO TO 1
      SINT=SQRT(SSINT)
      GO TO 2
    1 SINT=0.
    2 SKSKDO=SINT*SKABS*SKDABS
      GO TO (10,20,30,40,50,60,70,80,90,100),MODE
      RETURN
C
C -- PI+,P ---> PI+,P
C
   10 CTEMP1=FS3+SKSKDI*FP3
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(GP3)
      EFFSQUAR=TEMP1**2+SKSKDO**2*TEMP2**2
      RETURN
C
C -- PI+,N ---> PI+,N
C
   20 CTEMP1=2./3.*(FS1+SKSKDI*FP1+FD1*(3*COST**2-1.)/2.)
     &+1./3.*(FS3+SKSKDI*FP3)
      CTEMP2=2./3.*(GP1+GD1*3.*COST/SKABS**2)+1./3.*GP3
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(CTEMP2)
      EFFSQUAR=TEMP1**2+SKSKDO**2*TEMP2**2
      RETURN
C
C -- PI+,N ---> PI0,P
C
   30 CTEMP1=1.414214/3.*(-FS1-SKSKDI*FP1-FD1*(3*COST**2-1.)/2.
     &+FS3+SKSKDI*FP3)
      CTEMP2=1.414214/3.*(-GP1-GD1*3.*COST/SKABS**2+GP3)
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(CTEMP2)
      EFFSQUAR=(TEMP1**2+SKSKDO**2*TEMP2**2)*FEFCX
      RETURN
C
C -- PI0,P ---> PI0,P
C
   40 CTEMP1=1./3.*(FS1+SKSKDI*FP1+FD1*(3*COST**2-1.)/2.)
     &+2./3.*(FS3+SKSKDI*FP3)
      CTEMP2=1./3.*(GP1+GD1*3.*COST/SKABS**2)+2./3.*GP3
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(CTEMP2)
      EFFSQUAR=TEMP1**2+SKSKDO**2*TEMP2**2
      RETURN
C
C -- PI0,N ---> PI0,N
C
   50 CTEMP1=1./3.*(FS1+SKSKDI*FP1+FD1*(3*COST**2-1.)/2.)
     &+2./3.*(FS3+SKSKDI*FP3)
      CTEMP2=1./3.*(GP1+GD1*3.*COST/SKABS**2)+2./3.*GP3
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(CTEMP2)
      EFFSQUAR=TEMP1**2+SKSKDO**2*TEMP2**2
      RETURN
C
C -- PI0,P ---> PI+,N
C
   60 CTEMP1=1.414214/3.*(-FS1-SKSKDI*FP1-FD1*(3*COST**2-1.)/2.
     &+FS3+SKSKDI*FP3)
      CTEMP2=1.414214/3.*(-GP1-GD1*3.*COST/SKABS**2+GP3)
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(CTEMP2)
      EFFSQUAR=(TEMP1**2+SKSKDO**2*TEMP2**2)*FEFCX
      RETURN
C
C -- PI0,N ---> PI-,P
C
   70 CTEMP1=1.414214/3.*(-FS1-SKSKDI*FP1-FD1*(3*COST**2-1.)/2.
     &+FS3+SKSKDI*FP3)
      CTEMP2=1.414214/3.*(-GP1-GD1*3.*COST/SKABS**2+GP3)
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(CTEMP2)
      EFFSQUAR=(TEMP1**2+SKSKDO**2*TEMP2**2)*FEFCX
      RETURN
C
C -- PI-,P ---> PI-,P
C
   80 CTEMP1=2./3.*(FS1+SKSKDI*FP1+FD1*(3*COST**2-1.)/2.)
     &+1./3.*(FS3+SKSKDI*FP3)
      CTEMP2=2./3.*(GP1+GD1*3.*COST/SKABS**2)+1./3.*GP3
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(CTEMP2)
      EFFSQUAR=TEMP1**2+SKSKDO**2*TEMP2**2
      RETURN
C
C -- PI-,N ---> PI-,N
C
   90 CTEMP1=FS3+SKSKDI*FP3
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(GP3)
      EFFSQUAR=TEMP1**2+SKSKDO**2*TEMP2**2
      RETURN
C
C -- PI-,P ---> PI0,N
C
  100 CTEMP1=1.414214/3.*(-FS1-SKSKDI*FP1-FD1*(3*COST**2-1.)/2.
     &+FS3+SKSKDI*FP3)
      CTEMP2=1.414214/3.*(-GP1-GD1*3.*COST/SKABS**2+GP3)
      TEMP1=CABS(CTEMP1)
      TEMP2=CABS(CTEMP2)
      EFFSQUAR=(TEMP1**2+SKSKDO**2*TEMP2**2)*FEFCX
      RETURN
      END
