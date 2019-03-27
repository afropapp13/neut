*     -----------------------------------------------------
      SUBROUTINE EFCALAMP(SK,GK,PPI,radius)
*     -----------------------------------------------------
*
*     ( purpose )
*       CALCULATE AMPLITUDE OF EACH PHASE
*
*     ( input )
*       SK       : RELATIVE MOMENTUM VECTOR IN C-M SYSTEM
*       GK       : PI-N TOTAL MOMENTUM IN LAB. SYSTEM
*       PPI      : PION MOMENTUM IN LAB. SYSTEM
*       radius   : RADIUS
*
*     ( output )
*       COMMON/PSHIFT/FS1,FS3,FP1,FP3,GP1,GP3,FD1,GD1
*       amplitude of each phase

      DIMENSION SK(3),GK(3),PPI(3)
      COMPLEX FS1,FS3,FP1,FP3,GP1,GP3,FD1,GD1
      COMPLEX CTEMP1,CTEMP2
      COMPLEX EFS11MC,EFS31MC,EFP11MC,EFP13MC,EFP31MC,EFP33MC,EFD13
      COMPLEX EFD15
      COMPLEX BP11,BP13,BP31,BP33,BD13,BD15
      COMMON/PSHIFT/FS1,FS3,FP1,FP3,GP1,GP3,FD1,GD1
C
      Q=SQRT(SK(1)**2+SK(2)**2+SK(3)**2)
      FS1=EFS11MC(SK,GK,PPI,radius)
      FS3=EFS31MC(SK,GK,PPI,radius)
      BP11=EFP11MC(SK,GK,PPI,radius)
      BP13=EFP13MC(SK,GK,PPI,radius)
      BP31=EFP31MC(SK,GK,PPI,radius)
      BP33=EFP33MC(SK,GK,PPI,radius)
      BD13=EFD13(Q)
      BD15=EFD15(Q)
      FP1=BP11+2.*BP13
      FP3=BP31+2.*BP33
      FD1=2.*BD13+3.*BD15
      GP1=BP11-BP13
      GP3=BP31-BP33
      GD1=BD13-BD15
      RETURN
      END
