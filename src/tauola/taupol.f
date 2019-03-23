*  **************************************************************************
*  *                                                                        *
*  *                           TAU POLARIZATION                             *
*  *                                  IN                                    *
*  *                   TAU-NEUTRINO NUCLEON SCATTERING  v0.95               *
*  *                                                                        *
*  *        K. Hagiwara (KEK, JAPAN),  K. Mawatari (Kobe U., JAPAN),        *
*  *               and H. Yokoya (Hiroshima U. and RIKEN, JAPAN)            *
*  *                                                                        *
*  *            published in Nucl.Phys.B668 (2003) 364.                     *
*  *                                    [arXiv : hep-ph/0305324]            *
*  *                                                                        *
*  *                                                                        * 
*  *      This program provides the cross section and spin polarization     *
*  *      of the tau-lepton, produced by the tau-neutrino nucleon           *
*  *      scattering. Detail definitions are found in our paper.            *
*  *                                                                        *
*  *                                                                        *
*  *   BASIC INPUT                                                          *
*  *     + Particle Data : the mass of particles, Fermi constant            *
*  *       and Cabibbo angle are according to the values of PDG,            *
*  *       K.Hagiwara et al., Phys.Rev.D66(2002)010001.                     *
*  *                                                                        *
*  *     + Form Factors : original sets are according to our paper.         *
*  *       see also our another paper about pseudoscalar form factors,      *
*  *       [arXiv:hep-ph/0403076].                                          *
*  *                                                                        *
*  *     + Parton Distribution Functions : we use MRST2002 NLO set          * 
*  *       of A.D.Martin et al., Eur.Phys.J.C28(2003)455,                   *
*  *       [arXiv:hep-ph/021180]. you can change if you want.               *
*  *       "http://durpdg.dur.ac.uk/hepdata/pdf.html" is useful.            *
*  *                                                                        *
*  *     + Technical Parameter                                              *
*  *         * DE   : bin-size of QE histogram, correspond to               *
*  *                  the resolution of ETAU parameter.                     *
*  *                  QE region is extended to outside of the               *
*  *                  boundary by DE.                                       *
*  *         * WCUT : artificial boundary between RES and DIS.              *
*  *                                                                        *
*  *     + Isoscalar target : nucleon is set to isoscalar target,           *
*  *       (proton + nucleon) / 2, as default.                              * 
*  *                                                                        *  
*  *   INPUT PARAMETER                                                      *
*  *                                                                        *
*  *     + SIGN (Double Precision) : this value decide the incoming         *
*  *       neutrino type. SIGN = 1.D0 is for neutrino,                      *
*  *       SIGN = -1.D0 is for anti-neutrino.                               *
*  *                                                                        *
*  *     + ENU (GeV, Double Precision) : incoming neutrino energy           *
*  *       in laboratory frame.                                             *
*  *                                                                        *
*  *     + ETAU (GeV, Double Precision) : outgoing tau-lepton energy        *
*  *       in laboratory frame.                                             *
*  *                                                                        *
*  *     + THETA (DEG, Double Precision) : scattering angle of lepton       *
*  *       in laboratory frame.                                             *
*  *                                                                        *
*  *   OUTPUT PARAMETER                                                     *
*  *                                                                        *
*  *     + SIGMA (pb/GEV, Double Precision) : the cross section of          *
*  *       tau-lepton production, differentiated by ETAU and COS(THETA)     *
*  *       of laboratory frame.                                             *
*  *                                                                        *
*  *     + SX, SY, SZ (Double Precision) : spin polarization vector of      *
*  *       produced tau-lepton. defined in tau's rest frame in which        *
*  *       the z-axis is taken along its momentum direction in the          *
*  *       laboratory frame.                                                *
*  *                                                                        *
*  *     + FLAG (Integer) : this value tells which process should occur     *
*  *       for given input parameter sets. FLAG = 1 for QE,                 *
*  *       FLAG = 2 for RES and FLAG = 3 for DIS.                           *
*  *       FLAG = 0 for out of collision, return with zero-value output.    *
*  *                                                                        *
*  *                                                                        *
*  *   Common blocks /MASS/ and /CABIBBO/ are used.                         *
*  *                                                                        *
*  *   Subroutine MASSSET determine the mass of tau, nucleon, pion,         *
*  *    charm quark, delta resonance and W-boson.                           *
*  *                                                                        *
*  *   Subroutine EDGE gives minimal COS(CMIN) for given ENU, and           *
*  *    maximal/minimal ETAU(EMAX/EMIN) for given ENU and THETA.            *
*  *                                                                        *
*  *   Subroutine KINEMA choice the process which should occur by given     *
*  *    hadronic invariant mass and  return the FLAG                        *
*  *                                (1 for QE, 2 for RES and 3 for DIS).    *
*  *                                                                        *
*  *   Subroutine QE, RES and DIS give the structure functions W_{1~5}      *
*  *    for each processes.                                                 *
*  *                                                                        *
*  *                                                                        *
*  *   If you have any problem or question, please contact to H.Yokoya.     *
*  *           <yokoya@theo.phys.sci.hiroshima-u.ac.jp or yokoya@bnl.gov>   * 
*  *                                                                        *
*  **************************************************************************
      SUBROUTINE TAUPOL(SIGN,ENU,ETAU,THETA,SIGMA,SX,SY,SZ,FLAG,MODENE)

      IMPLICIT DOUBLE PRECISION (A-Z)
C

C     MODENE added by CWW (05/09/2004) mode from neut. Use this instead of 
C     mode calculated by KINEMA if NEUT says the event is QE.
      INTEGER MODENE

      INTEGER FLAG
      DIMENSION F(5)
C
      COMMON /MASS/ MTAU, MNUCL, MPI, MC, M33, MW
      COMMON /CABIBBO/ COSC, SINC
C
C--------------------- BASIC INPUT ---------------------
C
      PI = DACOS(-1.D0)
C
      CALL MASSSET (X) ! masses are determined in subroutine MASSSET
C     
C-- Fermi constant 
C
      GF = 1.16639D-5
C
C-- Cabibbo angle
C
      COSC = 0.975D0
      SINC = DSQRT(1.D0-COSC**2)
C
C-- Delta Etau at QE region.
C  To make histogram of QE precess, we introduce the resolution
C  of ETAU parameter, which is equal to the bin-size of histogram. 
C  Notice, QE region is extended to the outside of boundary!
C
C      DE = 0.03D0 * ENU

C     Changed by C.W.W to .05
      DE = 0.03D0 * ENU
C
C-- WCUT 
C  An artificial boundary between RES and DIS 
C  to avoid the double counting.
C
      WCUT = 1.4D0
C
C--------------------------------------------
C
      COS = DCOS(THETA*PI/180.D0)
      SIN = DSIN(THETA*PI/180.D0)
C
      PTAU = DSQRT(ETAU**2-MTAU**2)
C
      Q2 = 2.D0 * ENU * (ETAU-PTAU*COS) - MTAU**2
      PQ = MNUCL * (ENU - ETAU)
      W2 = MNUCL**2 + 2.D0*PQ - Q2
C
      CALL EDGE (ENU,COS,CMIN,EMIN,EMAX,FLAG)
      CALL KINEMA(ENU,ETAU,COS,W2,WCUT,DE,FLAG) 
C

C     ADDED BY CWW (05/09/2004) TO USE NEUT'S QE KINEMATICS
C
      IF (ABS(MODENE).EQ.1) THEN
         FLAG  = 1
         write(*,*) "ETAU:", etau, " EMIN ", emin-etau , " emax ", emax-etau

         IF ( ABS(EMAX-ETAU).LT.ABS(EMIN-ETAU)) THEN
            WRITE(*,*) "USING FORWARD"
            FLAG = 1 
         ELSE
            WRITE(*,*) "USING BACKWARD"
            FLAG = -1
         ENDIF

      ENDIF

      IF (FLAG.EQ.1) THEN
         CALL QE(F,DE,ENU,EMAX,PTAU,COS,Q2)
         ETAO = EMAX
      ELSEIF (FLAG.EQ.-1) THEN
         CALL QE(F,DE,ENU,EMIN,PTAU,COS,Q2)
         ETAO = EMIN
      ELSEIF (FLAG.EQ.2) THEN
         CALL RES(F,Q2,PQ,W2)
         ETAO = ETAU
      ELSEIF (FLAG.EQ.3) THEN
         CALL DIS(F,Q2,PQ,SIGN)
         ETAO = ETAU
      ELSEIF (FLAG.EQ.0) THEN
         SIGMA = 0.D0 ! out of kinetically arrowed region
         SX = 0.D0
         SY = 0.D0
         SZ = 0.D0
         GOTO 100
      ELSE
         WRITE(6,*) "ERROR FLAG"
         STOP
      ENDIF
C
      K = MW**2 / (Q2+MW**2)
      FAC = GF**2 * K**2 / 2.D0 / PI * PTAU / MNUCL
C
      FF = 
     1     (2.D0*F(1)+MTAU**2/MNUCL**2*F(4))*(ETAO-PTAU*COS)
     2     + F(2)*(ETAO+PTAU*COS)
     3     + SIGN*F(3)/MNUCL*(ENU*ETAO+PTAU**2-(ENU+ETAO)*PTAU*COS)
     4     - MTAU**2/MNUCL*F(5)
C
      SIGMA = FAC * FF * 0.389D9 ! pb
C
      SX = - SIGN*MTAU*SIN/2.D0*(2.D0*F(1)-F(2)+SIGN*ENU/MNUCL*F(3)-
     x     MTAU**2/MNUCL**2*F(4)+ETAO/MNUCL*F(5)) / FF
      SY = 0.D0
      SZ = - SIGN/2.D0*((2.D0*F(1)-MTAU**2/MNUCL**2*F(4))*
     z     (PTAU-ETAO*COS) + F(2)*(PTAU+ETAO*COS)
     z     + SIGN*F(3)/MNUCL*((ENU+ETAO)*PTAU-(ENU*ETAO+PTAU**2)*COS)
     z     - MTAU**2/MNUCL*F(5)*COS) / FF
C
      FLAG = ABS(FLAG)
C
 100  CONTINUE
C
      RETURN
C
      END
C
C--------------------------------------
C     This subroutine determine mass of elementary particles.
C
      SUBROUTINE MASSSET (X)
C
      IMPLICIT DOUBLE PRECISION (A-Z)
C
      COMMON /MASS/ MTAU, MNUCL, MPI, MC, M33, MW
C
      MTAU = 1.78D0             ! tau mass
      MNUCL = 0.938D0           ! nucleon mass
      MPI = 0.138D0             ! pion mass
      MC = 1.5D0                ! charm quark mass
      M33 = 1.232D0             ! delta resonance mass
      MW = 80.4D0               ! W-boson mass
C
      RETURN
C
      END
C
C----------------
C
C  This subroutine gives the kinematical edge 
C  of COS and ETAU. CMIN is the maximal angle for 
C  given neutrino energy. EMAX and EMIN are the 
C  maximal and minimal tau-lepton energy
C  for given neutrino energy and scattering angle.
C
      SUBROUTINE EDGE(ENU,COS,CMIN,EMIN,EMAX,FLAG)
C
      IMPLICIT DOUBLE PRECISION (A-Z)
C
      INTEGER FLAG
C
      COMMON /MASS/ MTAU, MNUCL, MPI, MC, M33, MW
C
      CALL MASSSET (X)
C
      ETHRE = (2.D0*MNUCL*MTAU+MTAU**2)/2.D0/MNUCL 
             ! threshold neutrino energy to produce tau-lepton.
C
      IF (ENU.LT.ETHRE) THEN
         FLAG = -1 ! too small incoming neutrino energy.
         RETURN
      ENDIF
C
      CMIN = DSQRT(1.D0 + MNUCL/ENU + MNUCL**2/ENU**2 -
     $     MNUCL**2/MTAU**2 - MTAU**2/4.D0/ENU**2)
C
      A = (ENU+MNUCL)**2 - (ENU*COS)**2
      B = (ENU + MNUCL)*(2.D0*MNUCL*ENU + MTAU**2)
      C = (ENU*MTAU*COS)**2 + (MNUCL*ENU + MTAU**2/2.D0)**2
      D = B**2 - 4.D0*A*C
C
      IF(COS.LT.CMIN) THEN
         FLAG = -2 ! too large angle
         RETURN
      ELSEIF(D.LT.0.D0) THEN
         FLAG = 0 ! out of kinematically arrowed tau-energy
         RETURN
      ENDIF
C
      EMIN = (B-DSQRT(D))/2.D0/A
      EMAX = (B+DSQRT(D))/2.D0/A
C
      RETURN
C
      END
C
C----------------------------------------------
C
C  This subroutine distinguishs the kinematics, which process 
C  occur or not, for given ENU, ETAU and COS.   
C    QE:1, RES:2, DIS:3, NOT ARROWED:0
C
      SUBROUTINE KINEMA(ENU,ETAU,COS,W2,WCUT,DE,FLAG)
C
      IMPLICIT DOUBLE PRECISION (A-Z)
C
      INTEGER FLAG
C
      COMMON /MASS/ MTAU, MNUCL, MPI, MC, M33, MW
C
      W = DSQRT(W2)
C
      IF (FLAG.LT.0) THEN
         FLAG = 0
         RETURN
      ENDIF
C
      IF (W.GE.WCUT) THEN
         FLAG = 3 ! DIS process
         RETURN
      ELSEIF (W.GT.(MNUCL+MPI).AND.W.LT.WCUT) THEN
         FLAG = 2 ! RES process
         RETURN
      ELSE 
         CALL EDGE(ENU,COS,CMIN,EMIN,EMAX,FLAG)
         FLAG = 0
         IF (ETAU.GT.EMIN-DE.AND.ETAU.LE.EMIN) FLAG = -1 
C                                  ! QE backward process
         IF (ETAU.LT.EMAX+DE.AND.ETAU.GE.EMAX) FLAG = 1
C                                  ! QE forward process         
         RETURN
      ENDIF
C
      RETURN
C
      END
C
C----------------------------------------------------------
C This subroutine gives the QE structure functions.
C 
      SUBROUTINE QE(F,DE,ENU,ETAU,PTAU,COS,Q2)

      IMPLICIT DOUBLE PRECISION (A-Z)

      DIMENSION F(5)

      COMMON /MASS/ MTAU, MNUCL, MPI, MC, M33, MW
      COMMON /CABIBBO/ COSC, SINC
C
      PTAU = DSQRT(ETAU**2-MTAU**2)
      Q2 = 2.D0 * ENU * (ETAU-PTAU*COS) - MTAU**2
      PQ = MNUCL * (ENU - ETAU)
C
      FDASH = DABS(2.D0*ENU*ETAU/PTAU*COS-2.D0*(ENU+MNUCL))
C                     ! Jacobian from delta(W2-M2) to delta(ETAU-E) 
C
      TAR = 0.5D0     ! isoscalar target correction
C
      FAC = COSC**2 / FDASH * PQ / DE * TAR 
C                        ! delta(ETAU-E) is replaced by 1./DE
C
C-- Mass parameters in QE.
C
      MV = 0.84D0
      MA = 1.2D0
C
C-- Vector Form Factors
C
      ZAI = 3.7D0
      FA0 = -1.267D0
C
      GE = 1.D0 / (1.D0 + Q2/MV**2)**2
      GM = (1.D0 + ZAI) * GE
C
      FV1 = (GE + Q2*GM/(4.D0*MNUCL**2)) / (1.D0+Q2/(4.D0*MNUCL**2))
C
      ZFV2 = (GM - GE) / (1.D0 + Q2/(4.D0*MNUCL**2))
C
C-- Axial Form Factor 
C
      FA = FA0 / (1.D0+Q2/MA**2)**2
C
C---- Pseudoscalar Form Factor
C
      FP = 2.D0 * (MNUCL**2) / (MPI**2+Q2)
     p     * FA0 / (1.D0+Q2/MA**2)**2 ! choose suppresion power
C
C-- Structure functions 
C 
      F(1) = FAC*((FV1+ZFV2)**2 + FA**2*(1.D0+2.D0*MNUCL**2/PQ))
C
      F(2) = FAC*(2.D0*(FV1+ZFV2)**2 + 2.D0*FA**2 + 2.D0*(-ZFV2)**2
     2     * (1.D0+PQ/(2.D0*MNUCL**2)) + 4.D0*(FV1+ZFV2)*(-ZFV2))
     2     * MNUCL**2 / PQ
C
      F(3) = FAC*(- 4.D0*(FV1+ZFV2)*FA) * MNUCL**2 / PQ
C
      F(4) = FAC*((-ZFV2)**2*(1.D0+PQ/2.D0/MNUCL**2)/2.D0
     4     + (FV1+ZFV2)*(-ZFV2) + FP**2*PQ/MNUCL**2
     4     - 2.D0 * FA*FP) * MNUCL**2 / PQ
C
      F(5) = FAC*(2.D0*(FV1+ZFV2)**2 + 2.D0*FA**2 + 2.D0*(-ZFV2)**2
     5     * (1.D0+PQ/2.D0/MNUCL**2) + 4.D0*(FV1+ZFV2)*(-ZFV2))
     5     * MNUCL**2 / PQ
C
      RETURN
C
      END
C
C-----------------------------------------------------------
C  This subroutine gives RES structure functions
C
      SUBROUTINE RES(F,Q2,PQ,W2)
C
      IMPLICIT DOUBLE PRECISION (A-Z)
C
      DIMENSION F(5)
      COMMON /MASS/ MTAU, MNUCL, MPI, MC, M33, MW
      COMMON /CABIBBO/ COSC, SINC
C
      W = DSQRT(W2)
C
      GAM = 0.12D0 * M33 / W
     1     * (DSQRT((W2-MNUCL**2-MPI**2)**2-4.D0*MNUCL**2*MPI**2))
     2     / (DSQRT((M33**2-MNUCL**2-MPI**2)**2-4.D0*MNUCL**2*MPI**2))
C
      BW = 1.D0 / DACOS(-1.D0) * W*GAM / ((W2-M33**2)**2+W2*GAM**2)
C
      TAR = 2.d0 ! isoscalar target 
C
      FAC = COSC**2 / 4.D0 * BW * TAR 
C
C-- Mass parameters in RES
C
      MV = 0.735D0
      MA = 1.D0
C
C-- Vector part form factors
C
      CV3 = 2.05D0 / (1.D0+Q2/MV**2)**2
*     3     / (1.D0+Q2/(4.D0*MV**2)) !-- PSY's factor 
C
      CV4 = - MNUCL / M33 * CV3
C
      CV5 = 0.D0
C
C-- Axial part form factors
C
      CA3 = 0.D0
C
      CA5 = 1.2D0 / (1.D0+Q2/MA**2)**2
*     5     / (1.D0+Q2/(3.D0*MA**2)) !-- PSY's factor
     5     * (1.D0-1.21D0*Q2/(2.D0+Q2)) !-- Adler's factor
C
      CA4 = - CA5 / 4.D0
C
C-- Pseudoscalar form factor
C
      CA6 = MNUCL**2 / (MPI**2+Q2) 
*     6     * 1.2D0 / (1.D0+Q2/MA**2)**2 ! can change suppression power
     6     * CA5    ! propotional to CA5. 
C
C -- Structure functions
C
      R1 =  (8.*(CA4**2*W**3*MNUCL*PQ**2 + 
     a     CV3*CV4*W**3*MNUCL*PQ**2 - 
     b     CV4**2*W**3*MNUCL*PQ**2 + 
     c     CV3*CV5*W**3*MNUCL*PQ**2 - 
     d     2*CV4*CV5*W**3*MNUCL*PQ**2 - 
     e     CV5**2*W**3*MNUCL*PQ**2 + 
     f     CA4**2*W**2*MNUCL**2*PQ**2 + 
     g     CV3**2*W**2*MNUCL**2*PQ**2 - 
     h     2*CV3*CV4*W**2*MNUCL**2*PQ**2 + 
     i     CV4**2*W**2*MNUCL**2*PQ**2 - 
     j     2*CV3*CV5*W**2*MNUCL**2*PQ**2 + 
     k     2*CV4*CV5*W**2*MNUCL**2*PQ**2 + 
     l     CV5**2*W**2*MNUCL**2*PQ**2 + 
     m     CV3*CV4*W*MNUCL**3*PQ**2 + 
     n     CV3*CV5*W*MNUCL**3*PQ**2 + CV3**2*MNUCL**4*PQ**2 + 
     o     CA4**2*W**2*PQ**3 + CV4**2*W**2*PQ**3 + 
     p     2*CV4*CV5*W**2*PQ**3 + CV5**2*W**2*PQ**3 + 
     q     CV3*CV4*W*MNUCL*PQ**3 + CV3*CV5*W*MNUCL*PQ**3 + 
     r     CV3**2*MNUCL**2*PQ**3 + 
     s     CA5**2*W**2*MNUCL**4*(W*MNUCL + MNUCL**2 + PQ) + 
     t     CA5*W*MNUCL**2*
     u     (CA3*MNUCL*(W**2*PQ + 2*W*MNUCL*(PQ - Q2) + 
     v     (MNUCL**2 + PQ)*(PQ - Q2)) + 
     w     2*CA4*W*(W*MNUCL + MNUCL**2 + PQ)*(PQ - Q2)) + 
     x     CA3*CA4*W*MNUCL*
     y     (W**2*PQ + 2*W*MNUCL*(PQ - Q2) + 
     z     (MNUCL**2 + PQ)*(PQ - Q2))*(PQ - Q2) + 
     a     CV3**2*W**3*MNUCL**3*Q2 + 
     b     CV3**2*W**2*MNUCL**4*Q2 - 
     c     2*CA4**2*W**3*MNUCL*PQ*Q2 - 
     d     CV3*CV4*W**3*MNUCL*PQ*Q2 + 
     e     2*CV4**2*W**3*MNUCL*PQ*Q2 + 
     f     2*CV4*CV5*W**3*MNUCL*PQ*Q2 - 
     g     2*CA4**2*W**2*MNUCL**2*PQ*Q2 + 
     h     4*CV3*CV4*W**2*MNUCL**2*PQ*Q2 - 
     i     2*CV4**2*W**2*MNUCL**2*PQ*Q2 + 
     j     2*CV3*CV5*W**2*MNUCL**2*PQ*Q2 - 
     k     2*CV4*CV5*W**2*MNUCL**2*PQ*Q2 - 
     l     2*CV3*CV4*W*MNUCL**3*PQ*Q2 - 
     m     CV3*CV5*W*MNUCL**3*PQ*Q2 - 
     n     2*CV3**2*MNUCL**4*PQ*Q2 - 2*CA4**2*W**2*PQ**2*Q2 - 
     o     2*CV4**2*W**2*PQ**2*Q2 - 
     p     2*CV4*CV5*W**2*PQ**2*Q2 - 
     q     2*CV3*CV4*W*MNUCL*PQ**2*Q2 - 
     r     CV3*CV5*W*MNUCL*PQ**2*Q2 - 
     s     2*CV3**2*MNUCL**2*PQ**2*Q2 + 
     t     CA4**2*W**3*MNUCL*Q2**2 - 
     u     CV4**2*W**3*MNUCL*Q2**2 + 
     v     CA4**2*W**2*MNUCL**2*Q2**2 - 
     w     2*CV3*CV4*W**2*MNUCL**2*Q2**2 + 
     x     CV4**2*W**2*MNUCL**2*Q2**2 + 
     y     CV3*CV4*W*MNUCL**3*Q2**2 + CV3**2*MNUCL**4*Q2**2 + 
     z     CA4**2*W**2*PQ*Q2**2 + CV4**2*W**2*PQ*Q2**2 + 
     a     CV3*CV4*W*MNUCL*PQ*Q2**2 + 
     b     CV3**2*MNUCL**2*PQ*Q2**2 + 
     c     CA3**2*MNUCL**2*
     d     ((MNUCL**2 + PQ)*(PQ - Q2)**2 - W**3*MNUCL*Q2 + 
     e     W**2*(PQ**2 + MNUCL**2*Q2))))/
     f     (3.*W**2*MNUCL**4)
C
      R2 =   (8.*PQ*(CA5**2*MNUCL**4*(W*MNUCL + MNUCL**2 + PQ) + 
     a     CA3*CA5*W*MNUCL**3*Q2 + 
     b     Q2*(CV3*CV4*W**3*MNUCL - CV4**2*W**3*MNUCL + 
     c     CV3*CV5*W**3*MNUCL - 2*CV4*CV5*W**3*MNUCL - 
     d     CV5**2*W**3*MNUCL + CV3**2*W**2*MNUCL**2 - 
     e     2*CV3*CV4*W**2*MNUCL**2 + 
     f     CV4**2*W**2*MNUCL**2 - 
     g     2*CV3*CV5*W**2*MNUCL**2 + 
     h     2*CV4*CV5*W**2*MNUCL**2 + 
     i     CV5**2*W**2*MNUCL**2 + CV3*CV4*W*MNUCL**3 + 
     j     CV3*CV5*W*MNUCL**3 + CV3**2*MNUCL**4 + 
     k     CV4**2*W**2*PQ + 2*CV4*CV5*W**2*PQ + 
     l     CV5**2*W**2*PQ + CV3*CV4*W*MNUCL*PQ + 
     m     CV3*CV5*W*MNUCL*PQ + CV3**2*MNUCL**2*PQ + 
     n     CA3**2*MNUCL**2*(W**2 + MNUCL**2 + PQ) + 
     o     CA4**2*W**2*(W*MNUCL + MNUCL**2 + PQ) + 
     p     CA3*CA4*W*MNUCL*
     q     (W**2 + 2*W*MNUCL + MNUCL**2 + PQ) + 
     r     CV3*CV5*W*MNUCL*Q2 - CV5**2*W*MNUCL*Q2 + 
     s     CV5**2*MNUCL**2*Q2 + CV5**2*PQ*Q2)))/
     t     (3.*W**2*MNUCL**4)
C
      R3 =   (16.*PQ*(CA5*W*MNUCL**2*
     a     (W*(CV5*PQ + CV4*(PQ - Q2)) + 
     b     CV3*MNUCL*(2*W**2 + 2*W*MNUCL - PQ + Q2)) + 
     c     CA4*W*(PQ - Q2)*
     d     (W*(CV5*PQ + CV4*(PQ - Q2)) + 
     e     CV3*MNUCL*(2*W**2 + 2*W*MNUCL - PQ + Q2)) - 
     f     CA3*MNUCL*(CV5*W*PQ*
     g    (-2*W**2 + 2*W*MNUCL + PQ - Q2) + 
     h     CV4*W*(PQ - Q2)*
     i     (-2*W**2 + 2*W*MNUCL + PQ - Q2) + 
     j     CV3*MNUCL*(2*(PQ - Q2)**2 + W**2*(-4*PQ + 3*Q2)))
     k     ))/(3.*W**2*MNUCL**4)
C
      R4 =   (8.*PQ*(-2*CA5*CA6*W**3*MNUCL**3 - 
     a     CV3**2*W**3*MNUCL**3 - 2*CA5*CA6*W**2*MNUCL**4 - 
     b     CV3**2*W**2*MNUCL**4 + CA5**2*W*MNUCL**5 + 
     c     CA5**2*MNUCL**6 + CV3*CV4*W**3*MNUCL*PQ - 
     d     2*CV4**2*W**3*MNUCL*PQ - 
     e     2*CV4*CV5*W**3*MNUCL*PQ - 
     f     2*CA5*CA6*W**2*MNUCL**2*PQ - 
     g     4*CV3*CV4*W**2*MNUCL**2*PQ + 
     h     2*CV4**2*W**2*MNUCL**2*PQ - 
     i     2*CV3*CV5*W**2*MNUCL**2*PQ + 
     j     2*CV4*CV5*W**2*MNUCL**2*PQ + 
     k     2*CA5*CA6*W*MNUCL**3*PQ + 
     l     2*CV3*CV4*W*MNUCL**3*PQ + CV3*CV5*W*MNUCL**3*PQ + 
     m     CA5**2*MNUCL**4*PQ + 2*CA5*CA6*MNUCL**4*PQ + 
     n     2*CV3**2*MNUCL**4*PQ + 2*CV4**2*W**2*PQ**2 + 
     o     2*CV4*CV5*W**2*PQ**2 + CA6**2*W*MNUCL*PQ**2 + 
     p     2*CV3*CV4*W*MNUCL*PQ**2 + 
     q     2*CV3*CV5*W*MNUCL*PQ**2 - CV5**2*W*MNUCL*PQ**2 + 
     r     2*CA5*CA6*MNUCL**2*PQ**2 + CA6**2*MNUCL**2*PQ**2 + 
     s     2*CV3**2*MNUCL**2*PQ**2 + CV5**2*MNUCL**2*PQ**2 + 
     t     CA6**2*PQ**3 + CV5**2*PQ**3 + 
     u     2*CA4*W**2*(W*MNUCL + MNUCL**2 + PQ)*
     v     (CA5*MNUCL**2 - CA6*PQ) - 
     w     CA3**2*MNUCL**2*
     x     (-(W**3*MNUCL) + W**2*MNUCL**2 - 
     y     (MNUCL**2 + PQ)*(2*PQ - Q2)) + 
     z     CA4**2*W**2*(W*MNUCL + MNUCL**2 + PQ)*
     a     (2*PQ - Q2) + CA6**2*W**3*MNUCL*Q2 + 
     b     CV4**2*W**3*MNUCL*Q2 + CA6**2*W**2*MNUCL**2*Q2 + 
     c     2*CV3*CV4*W**2*MNUCL**2*Q2 - 
     d     CV4**2*W**2*MNUCL**2*Q2 - 
     e     2*CA5*CA6*W*MNUCL**3*Q2 - CV3*CV4*W*MNUCL**3*Q2 - 
     f     2*CA5*CA6*MNUCL**4*Q2 - CV3**2*MNUCL**4*Q2 + 
     g     CA6**2*W**2*PQ*Q2 - CV4**2*W**2*PQ*Q2 - 
     h     2*CA6**2*W*MNUCL*PQ*Q2 - CV3*CV4*W*MNUCL*PQ*Q2 - 
     i     2*CA5*CA6*MNUCL**2*PQ*Q2 - 2*CA6**2*MNUCL**2*PQ*Q2 - 
     j     CV3**2*MNUCL**2*PQ*Q2 - 2*CA6**2*PQ**2*Q2 + 
     k     CA6**2*W*MNUCL*Q2**2 + CA6**2*MNUCL**2*Q2**2 + 
     l     CA6**2*PQ*Q2**2 + 
     m     CA3*W*MNUCL*(CA5*MNUCL**2*
     n     (2*W*MNUCL + MNUCL**2 + 2*PQ) + 
     o     CA4*(W**2*PQ + W*MNUCL*(4*PQ - 2*Q2) + 
     p     (MNUCL**2 + PQ)*(2*PQ - Q2)) - 
     q     CA6*PQ*(W**2 + 2*W*MNUCL + MNUCL**2 + Q2))))/
     r     (3.*W**2*MNUCL**4)
C
      R5 =   (8.*PQ*(2*CA3**2*MNUCL**2*PQ*(W**2 + MNUCL**2 + PQ) + 
     a     CA3*W*MNUCL*(2*CA4*PQ*
     b     (W**2 + 2*W*MNUCL + MNUCL**2 + PQ) - 
     c     CA6*Q2*(W**2 + 2*W*MNUCL + MNUCL**2 + Q2) + 
     d     CA5*MNUCL**2*
     e     (W**2 + 2*W*MNUCL + MNUCL**2 + 2*PQ + Q2)) + 
     f     2*(CA5**2*MNUCL**4*(W*MNUCL + MNUCL**2 + PQ) + 
     g     CA4**2*W**2*PQ*(W*MNUCL + MNUCL**2 + PQ) + 
     h     CA5*CA6*MNUCL**2*(W*MNUCL + MNUCL**2 + PQ)*
     i     (PQ - Q2) + 
     j     CA4*W**2*(W*MNUCL + MNUCL**2 + PQ)*
     k     (CA5*MNUCL**2 - CA6*Q2) + 
     l     PQ*(CV3**2*MNUCL**2*(W**2 + MNUCL**2 + PQ) - 
     m     (W*MNUCL - MNUCL**2 - PQ)*
     n     (CV4**2*W**2 + 2*CV4*CV5*W**2 + 
     o     CV5**2*(W**2 + Q2)) + 
     p     CV3*W*MNUCL*
     r     (CV4*(W**2 - 2*W*MNUCL + MNUCL**2 + PQ) + 
     s     CV5*(W**2 - 2*W*MNUCL + MNUCL**2 + PQ + 
     t     Q2))))))/(3.*W**2*MNUCL**4)
C
      F(1) = FAC * R1
      F(2) = FAC * R2 * MNUCL**2 / PQ
      F(3) = FAC * R3 * MNUCL**2 / PQ
      F(4) = FAC * R4 * MNUCL**2 / PQ
      F(5) = FAC * R5 * MNUCL**2 / PQ
C
      RETURN
C
      END
C
C-----------------------------------------------------------
C  This subroutine gives DIS structure functions
C
      SUBROUTINE DIS(F,Q2,PQ,SIGN)
C
      IMPLICIT DOUBLE PRECISION (A-Z)

      INTEGER MODE

      DIMENSION F(5)

      COMMON /MASS/ MTAU, MNUCL, MPI, MC, M33, MW
      COMMON /CABIBBO/ COSC, SINC
C 
      MODE = 1 ! mode of called pdf data.
C
      X = Q2 / (2.D0*PQ)   ! Bjorken x
      XC = X * (1.D0+MC**2/Q2)   ! slow rescaling
C
      Q = DSQRT(Q2)
      IF(Q2.LT.1.25D0) Q = DSQRT(1.25D0)
C
C-- Light quark production
C     parton distribution functions for neutrino scatterings 
C     * isoscalar target : (proton + neutron ) / 2 
C
      CALL MRST2002(X,Q,MODE,UPV,DNV,US,DS,STR,CHM,BOT,GLU)
*      call mrs99(x,q,mode,upv,dnv,us,ds,str,chm,bot,glu)
C
      IF(SIGN.EQ.1.D0) THEN
         PDF = ((UPV+US+DNV+DS)*COSC**2 + 2.D0*STR*SINC**2) 
     -        / 2.D0 / X
         PDFBAR = (US + DS) / 2.D0 / X
      ELSEIF(SIGN.EQ.-1.D0) THEN
         PDF = (UPV+US+DNV+DS) / 2.D0 / X
         PDFBAR = ((US+DS)*COSC**2 + 2.D0*STR*SINC**2)
     -        / 2.D0 / X
      ENDIF
C
      F1DIS = (PDF + PDFBAR)
     1     * (1.D0+X*MNUCL**2/PQ) ! our modification, see Eq.(55).
      F2DIS = 2.D0 * X * (PDF + PDFBAR)
      F3DIS = 2.D0 * (PDF - PDFBAR)
      F4DIS = 0.D0
      F5DIS = 2.D0 * (PDF + PDFBAR)
C
C-- Heavy quark production
C
      IF(XC.GT.1.D0) THEN
         MPDF = 0.D0
         MPDFBAR = 0.D0
      ELSE
C
         CALL MRST2002(XC,Q,MODE,UPV,DNV,US,DS,STR,CHM,BOT,GLU)
*         call mrs99(xc,q,mode,upv,dnv,us,ds,str,chm,bot,glu)
C
         IF(SIGN.EQ.1.D0) THEN
            MPDF = (2.D0*STR*COSC**2 + (UPV+US+DNV+DS)*SINC**2)
     -           / 2.D0 / XC
            MPDFBAR = 0.D0
         ELSEIF(SIGN.EQ.-1.D0) THEN
            MPDF = 0.D0
            MPDFBAR = ((US+DS)*SINC**2 + 2.D0*STR*COSC**2)
     -           / 2.D0 / XC
         ENDIF
      ENDIF
C
      MF1DIS = (MPDF + MPDFBAR)
     1     * (1.D0+XC*MNUCL**2/PQ) ! our modification factor
      MF2DIS = 2.D0 * XC * (MPDF + MPDFBAR)
      MF3DIS = 2.D0 * (MPDF - MPDFBAR) 
      MF4DIS = 0.D0
      MF5DIS = 2.D0 * (MPDF + MPDFBAR) 
C
      F(1) = F1DIS + MF1DIS
      F(2) = (F2DIS + MF2DIS) * MNUCL**2 / PQ
      F(3) = (F3DIS + MF3DIS) * MNUCL**2 / PQ
      F(4) = (F4DIS + MF4DIS) * MNUCL**2 / PQ
      F(5) = (F5DIS + MF5DIS) * MNUCL**2 / PQ
C
      RETURN
C
      END
C
C-----------------------------------------------------------------
C END of program 
