*     -----------------------------------------------------
      SUBROUTINE EFPNCRSH_PNSEP(PAIMOM,NQPAI,NQN,
     $                    SIGEL,SIGINEL,SIGCX,IERROR)
*     -----------------------------------------------------
*
*     ( purpose )
C        CROSS SECTION OF PAI-N INTERACTION AT HIGH ENERGIES (>800MEV).
C                    T.KAJITA  25-JUN-88
*
*     ( input )
C           PAIMOM : MOMENTUM OF PION IN MEV/C
C           NQPAI  : CHARGE OF PION
C           NQN    : CHARGE OF NUCLEON
*
*     ( output )
C           SIGEL : (QUASI)ELASTIC CROSS SECTION IN MB.
C           SIGINEL : INELASTIC CROSS SECTION IN MB
C           SIGCX : CHARGE EXCHANGE CROSS SECTION IN MB
C           IERROR : ERROR CODE  =0 -- NO ERROR
C                                =1 -- ERROR ( PAIMOM < 20MEV )
C
C     (ASSUMPTIONS)
C  (1) SIGMA(PAI(0)-N) is average(SIGMA(PAI(+/-)-p))
C  (2) ISOSPIN SYMMETRY
C           SIGMA(PAI(-)-n) IS SAME AS SIGMA(PAI(+)-p)
C           SIGMA(PAI(+)-n) IS SAME AS SIGMA(PAI(-)-p)
C
C   2010/06 Output Elastic, Inelastic, CX cross sections
C           Use SAID PWA of free pion-nucleon cross sections
C             (http://gwdac.phys.gwu.edu/analysis/pin_analysis.html)
C                         - P.de Perio
C
C
C
      IMPLICIT NONE
      REAL PAIMOM,SIGTOTP,SIGTOTN,SIGEL,SIGINEL,SIGCX
      REAL SIGELP,SIGINELP,SIGELN,SIGINELN,SIGCXP,SIGCXN
      INTEGER NQPAI,NQN, IERROR
      INTEGER IE, IEL, IEH
      REAL DEL, DEH
      REAL PMEV,SEPPP,SIPPP,SEPMP,SIPMP,SCPMP
      DIMENSION PMEV(121),SEPPP(121),SIPPP(121)
      DIMENSION SEPMP(121),SIPMP(121),SCPMP(121)
C
C     -- PION MOMENTUM IN MEV/C
      DATA PMEV  /
     $     20,30,40,50,60,70,80,
     $     90,100,110,120,130,140,150,
     $     160,170,180,190,200,210,220,
     $     230,240,250,260,270,280,290,
     $     300,310,320,330,340,350,360,
     $     370,380,390,400,410,420,425,
     $     455,485,515,545,575,605,635,
     $     665,695,725,755,785,815,845,
     $     875,905,935,965,995,1025,1055,
     $     1085,1115,1145,1175,1205,1235,1265,
     $     1295,1325,1355,1385,1415,1445,1475,
     $     1505,1535,1565,1595,1625,1655,1685,
     $     1715,1745,1775,1805,1835,1865,1895,
     $     1925,1955,1985,2015,2045,2075,2105,
     $     2135,2165,2195,2225,2255,2285,2315,
     $     2345,2375,2400,2660,2920,3180,3440,
     $     3700,3960,4220,4480,4740,5000,10000,
     $     20000,50000/
C
C     -- QE Scatter CROSS SECTION pip-P
      DATA SEPPP /
     $     1.02,1.28,1.51,1.79,2.17,2.68,3.37,
     $     4.31,5.54,7.14,9.2,11.8,15.2,19.4,
     $     24.7,31.4,39.7,50,62.8,78.1,96.3,
     $     117,140,162,182,197,204,202,
     $     195,182,167,152,137,123,110,
     $     99.3,89.5,80.9,73.4,66.9,61.1,58.5,
     $     45.7,36.7,30.1,25.1,21.2,18.1,15.6,
     $     13.6,11.9,10.4,9.16,8.13,7.46,7.43,
     $     8.15,9.31,10.5,11.3,11.9,12.2,12.5,
     $     12.7,12.9,13.1,13.5,14,14.6,15.3,
     $     16.1,16.8,17.5,17.9,18,17.8,17.3,
     $     16.7,16,15.2,14.5,13.8,13.1,12.5,
     $     12,11.6,11.2,10.8,10.5,10.2,9.95,
     $     9.72,9.52,9.34,9.17,9.02,8.89,8.77,
     $     8.67,8.57,8.49,8.41,8.35,8.3,8.25,
     $     8.21,8.19,8.79,8.55,8.01,7.5,7.14,
     $     6.92,6.63,6.3,6.02,5.8,5.63,4.75,
     $     3.5,3.29/
C
C     -- Inelastic CROSS SECTION pip-P
      DATA SIPPP /
     $     0,0,0,0,0,0,0,
     $     0,0,0,0,0,0,0,
     $     0,0,0,0,0,0,0,
     $     0,0,0,0,0,0.0126,0.077,
     $     0.138,0.192,0.237,0.275,0.305,0.331,0.353,
     $     0.372,0.39,0.409,0.428,0.449,0.473,0.486,
     $     0.59,0.752,0.982,1.29,1.68,2.17,2.77,
     $     3.48,4.31,5.29,6.43,7.77,9.27,10.7,
     $     11.8,12.4,12.6,12.6,12.6,12.6,12.8,
     $     13.2,13.7,14.3,15.2,16.2,17.3,18.5,
     $     19.7,20.8,21.8,22.6,23.1,23.3,23.2,
     $     22.9,22.5,22,21.5,21.1,20.6,20.3,
     $     20,19.7,19.6,19.4,19.4,19.3,19.3,
     $     19.3,19.4,19.5,19.6,19.7,19.9,20,
     $     20.2,20.4,20.5,20.7,20.9,21.1,21.3,
     $     21.5,21.7,21.8,22.3,21.7,21.3,21.2,
     $     21.1,21.2,21.1,21.4,21,20.9,20.1,
     $     20.1,19.8/
C
C     -- QE Scatter CROSS SECTION pim-P
      DATA SEPMP /
     $     2.57,2.3,2.13,2.01,1.91,1.83,1.77,
     $     1.74,1.73,1.76,1.85,1.99,2.21,2.52,
     $     2.96,3.54,4.32,5.32,6.59,8.16,10.1,
     $     12.3,14.7,17.1,19.2,20.9,21.8,21.9,
     $     21.4,20.4,19.2,17.8,16.4,15.2,14,
     $     13,12.1,11.4,10.8,10.3,10,9.88,
     $     9.67,10.1,10.8,11.3,11.7,12.3,13.4,
     $     15.4,18.5,19.5,18.2,15.7,13.8,13,
     $     13.6,15.6,19,23,25.3,24.7,22.3,
     $     19.5,17.2,15.5,14.3,13.4,12.9,12.6,
     $     12.3,12.2,12,11.9,11.7,11.5,11.3,
     $     11.1,10.9,10.7,10.5,10.4,10.2,10.1,
     $     9.97,9.88,9.8,9.73,9.67,9.62,9.57,
     $     9.52,9.48,9.44,9.4,9.36,9.31,9.26,
     $     9.21,9.15,9.09,9.03,8.97,8.9,8.83,
     $     8.76,8.69,8.62,8.06,7.7,7.35,7.05,
     $     6.8,6.53,6.23,5.97,5.75,5.57,4.7,
     $     4.2,3.45/
C
C     -- Inelastic CROSS SECTION pim-P
      DATA SIPMP /
     $     0.00129,0.00257,0.00467,0.00731,0.0107,0.0147,0.0195,
     $     0.0252,0.0317,0.0392,0.0478,0.0576,0.069,0.0821,
     $     0.0972,0.115,0.135,0.158,0.185,0.215,0.248,
     $     0.284,0.319,0.351,0.377,0.393,0.404,0.443,
     $     0.477,0.509,0.543,0.583,0.634,0.697,0.775,
     $     0.872,0.989,1.13,1.3,1.5,1.74,1.87,
     $     2.86,4.16,5.55,6.82,7.96,9.15,10.6,
     $     12.5,15.7,18.1,18.3,17.5,17.2,17.6,
     $     18.9,20.7,23,25.3,26.7,26.5,25.1,
     $     23.6,22.3,21.6,21.2,21,21.1,21.2,
     $     21.4,21.7,21.9,22.1,22.2,22.3,22.4,
     $     22.4,22.4,22.5,22.5,22.6,22.8,22.9,
     $     23.1,23.4,23.6,23.9,24.1,24.4,24.7,
     $     24.9,25.2,25.4,25.6,25.8,25.9,26,
     $     26.1,26.1,26,26,25.8,25.7,25.5,
     $     25.4,25.2,25,24.4,24.5,24,24.5,
     $     24.2,24,23.8,24.1,23.6,23.5,22.1,
     $     20.8,20.6/
C
C     -- Single CX CROSS SECTION pim-P
      DATA SCPMP /
     $     5.76,4.97,4.63,4.49,4.45,4.5,4.64,
     $     4.86,5.17,5.6,6.15,6.87,7.77,8.91,
     $     10.3,12.1,14.2,16.9,20.1,23.8,28.2,
     $     32.9,37.8,42.3,46,48.3,48.8,47.8,
     $     45.4,42.2,38.7,35.1,31.8,28.8,26.1,
     $     23.8,21.8,20.1,18.6,17.4,16.4,15.9,
     $     13.9,12.7,11.7,10.7,9.93,9.43,9.31,
     $     9.64,9.63,8.67,7.35,6.07,5.22,4.94,
     $     5.25,6.1,7.22,8.05,8.16,7.43,6.21,
     $     5,4.07,3.44,3.04,2.83,2.74,2.74,
     $     2.78,2.86,2.92,2.96,2.94,2.88,2.76,
     $     2.62,2.45,2.29,2.13,1.98,1.85,1.73,
     $     1.63,1.54,1.46,1.39,1.32,1.27,1.21,
     $     1.16,1.11,1.07,1.02,0.978,0.935,0.894,
     $     0.854,0.815,0.777,0.741,0.706,0.674,0.644,
     $     0.616,0.592,0.556,0.459,0.371,0.278,0.213,
     $     0.163,0.133,0.11,0.0913,0.0778,0.0685,0.00685,
     $     0.000685,0/
C
C
C
      SIGINEL=0.
      SIGEL=0.
      SIGCX=0.
      IERROR=1
      IF(PAIMOM.LT.20.)RETURN
      IERROR=0
C
      DO 10 IE=2,121
         IF(PMEV(IE).GT.PAIMOM) THEN
            IEH=IE
            IEL=IE-1
            DEH=ABS(PMEV(IEH)-PAIMOM)
            DEL=ABS(PAIMOM-PMEV(IEL))
            GO TO 11
         END IF
   10 CONTINUE
      IEH=121
      IEL=120
      DEH=0.
      DEL= ABS (PMEV(IEH)-PMEV(IEL))
   11 CONTINUE
C
C
      IF(NQPAI.EQ.0) THEN
         SIGEL = ( (DEL*SEPMP(IEH) + DEH*SEPMP(IEL)) +
     $  (DEL*SEPPP(IEH) + DEH*SEPPP(IEL)) ) / (DEL+DEH)/2
         SIGINEL = ( (DEL*SIPMP(IEH) + DEH*SIPMP(IEL)) +
     $  (DEL*SIPPP(IEH) + DEH*SIPPP(IEL)) ) / (DEL+DEH)/2
         SIGELP = SIGEL
         SIGELN = SIGEL
         SIGINELP = SIGINEL
         SIGINELN = SIGINEL
         SIGCXP = (DEL*SCPMP(IEH) + DEH*SCPMP(IEL))/(DEL+DEH)/2
         SIGCXN = (DEL*SCPMP(IEH) + DEH*SCPMP(IEL))/(DEL+DEH)/2
      ELSE IF(NQPAI.EQ.1) THEN
         SIGELP = ( DEL*SEPPP(IEH) + DEH*SEPPP(IEL) ) / ( DEL+DEH )
         SIGINELP = ( DEL*SIPPP(IEH) + DEH*SIPPP(IEL) ) / ( DEL+DEH )
         SIGELN = ( DEL*SEPMP(IEH) + DEH*SEPMP(IEL) ) / ( DEL+DEH )
         SIGINELN = ( DEL*SIPMP(IEH) + DEH*SIPMP(IEL) ) / ( DEL+DEH )
         SIGCXP = 0
         SIGCXN = (DEL*SCPMP(IEH) + DEH*SCPMP(IEL)) / (DEL+DEH)
      ELSE IF(NQPAI.EQ.-1) THEN
         SIGELP = ( DEL*SEPMP(IEH) + DEH*SEPMP(IEL) ) / ( DEL+DEH )
         SIGINELP = ( DEL*SIPMP(IEH) + DEH*SIPMP(IEL) ) / ( DEL+DEH )
         SIGELN = ( DEL*SEPPP(IEH) + DEH*SEPPP(IEL) ) / ( DEL+DEH )
         SIGINELN = ( DEL*SIPPP(IEH) + DEH*SIPPP(IEL) ) / ( DEL+DEH )
         SIGCXP = (DEL*SCPMP(IEH) + DEH*SCPMP(IEL)) / (DEL+DEH)
         SIGCXN = 0
      ELSE
         PRINT *, 'Error in pion charge in nuclear effect'
         STOP
      ENDIF

C      SIGTOTP = SIGELP+SIGINELP+SIGCXP
C      SIGTOTN = SIGELN+SIGINELN+SIGCXN

      IF(NQN.EQ.1) THEN
         SIGEL = SIGELP
         SIGINEL = SIGINELP
         SIGCX = SIGCXP
      ELSE IF(NQN.EQ.0) THEN
         SIGEL = SIGELN
         SIGINEL = SIGINELN
         SIGCX = SIGCXN
      ELSE
         PRINT *, 'Error in nucleon charge in nuclear effect'
         STOP
      ENDIF

C
C
      RETURN
      END
