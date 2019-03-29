C
C     COMMON BLOCK WHICH ARE USED IN TAUOLA
C
      REAL GFERMI, GV, GA, CCABIB, SCABIB, GAMEL

      COMMON /DECPAR/ GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
C
C     GFERMI : Fermi coupling.  G=1.16637E-5 (GeV^-2)
C     GV     : Tau vector coupling. GV=1 in the standard model
C     GA     : Tau axial coupling. GA=-1 in the standard model
C     CCABIB : Cosine of Cabibbo angle
C     SCABIB : Sine of Cabibo angle
C     GAMEL  : Branching ratio for tau -> neu neu_bar e
C
      REAL AMTAU, AMNUTA, AMEL, AMNUE, AMMU, AMNUMU, AMPIZ, AMPI, AMRO,
     &   GAMRO, AMA1, GAMA1, AMK, AMKZ, AMKST, GAMKST

      COMMON /PARMAS/ AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU,
     &                AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1,
     &                AMK,AMKZ,AMKST,GAMKST
C
C     AMTAU  : tau mass
C     AMNUTA : neu-tau mass, non-zero mass required for numerical stability
C     AMEL   : e mass
C     AMNUE  : neu-e mass ( dummy parameter )
C     AMMU   : mu mass
C     AMNUMU : neu-mu mass ( dummy parameter )
C     AMPIZ  : pi0 mass
C     AMPI   : pi+- mass
C     AMRO   : ro mas, for crude MC distributions.
C              See Comp. Phys. Comm. 64 (1991) 275-299
C     GAMRO  : ro width, for crude MC distributions.
C     AMA1   : a1 mass
C     GAMA1  : a1 width
C     AMK    : K+- mass
C     AMKZ   : K0 mass
C     AMKST  : K* mass
C     GAMKST : K* width
C
      REAL    GAMPRT
      INTEGER JLIST, NCHAN

      COMMON /TAUBRA/ GAMPRT(30),JLIST(30),NCHAN
C
C     NCHAN  : Number of decay channels. At present NCHAN=8. The channels
C              appropriately ordered are: e neu neu-b, mu neu nue-b, pi neu,
C              rho neu, a1 neu, K neu, K* neu, N pi neu.
C     JLIST(I)  : Number of the decay channel according to the above list.
C                 If decay modes are ordered as above JLIST(I)=I
C     GAMPRT(I) : Branching ratio for the JLIST(I) decay mode. Arbitary units.
C                 These parameters define actual proportion of decays in the
C                 sample to be generated.
C
      REAL    CBRNPI, AMAS
      INTEGER KPI, MULT
      COMMON /TAUNPI/ CBRNPI(4),AMAS(6,4),KPI(6,4),MULT(4)
C
C     CBRNP(I)  : Individual branching ratios of different type of multipion
C                 final states relative to the total multipion branching ratio.
C     AMAS(J,I) : Mass of the J-th pi for the I-th type of multipion final
C                 states
C     KPI(J,I)  : Type of the J-th pi for the I-th type of multipion final
C                 states. KPI(J,I)=-17,17,23,0 denote respectively pi of
C                 same/opposite charge as the mother tau, pi0 and no pi at all
C     MULT(I)   : Multiplicity of the I-th type of multipion decay mode.
C
      INTEGER JAK1, JAK2, JAKP, JAKM, KTOM
      COMMON /JAKI/   JAK1,JAK2,JAKP,JAKM,KTOM
C
C     JAK1   : Type of the tau+ decay mode according to the list given in the
C              previous table.
C              JAK1 = 0 denotes inclusive tau decay,
C              JAK1=-1 no decay at all.
C     JAK2   : The same as JAK1 but for tau -
C     JAKP   : Internal variable, does not require initialization.
C     JAKM   : Internal variable, does not require initialization.
C     KTOM   : Internal variable, does not require initialization.
C
      INTEGER IDFF
      COMMON /IDFC/   IDFF
C
C     IDFF   : Lund indentifier for tau+, should be set to IDFF=-15
C
      INTEGER INUT, IOUT
      COMMON /INOUT/  INUT,IOUT
C
C     INUT   : LUN for input ( dummy )
C     IOUT   : LUN for output
C
      INTEGER IA1
      COMMON /IDPART/ IA1
C
C     Lund type identifier for a1
C

      INTEGER ITDKRC
      REAL*8            XK0DEC
      COMMON / TAURAD / XK0DEC,ITDKRC
C
C     XK0DEC : Soft-photon cut parameter ( 10**-4 - 10**-3 )
C     ITDKRC : QED correction switch ( 1 : on, 0 : off )
C
      REAL*8           ALFINV,ALFPI,XK0
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
C
C     ALFINV : 1/alpha
C     ALFPI  : alpha/pi
C     XK0    :
C
      REAL BRA1, BRK0, BRK0B, BRKS
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
C
C     Decay paramater of a1, K0, K0B and K*
C     BRA1  : 0.5, (0,1), if 0, a1 -> 2PI0,PI+(-)
C                         if 1, a1 -> three charged pions
C     BRK0  : 0.5, (0,1), if 0, K0 -> Kl
C                         if 1, K0 -> Ks
C     BRK0B : 0.5, (0,1), if 0, K0B -> Kl
C                         if 1, K0B -> Ks
C     BRK0S : 0.667, (0,1), if 0, K* -> K+(-),PI0
C                           if 1, K* -> K0,PI+(-)
C
      INTEGER NMODE, NM1, NM2, NM3, NM4, NM5, NM6, IDFFIN, MULPIK
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31
C
C
C
      INTEGER NEVDEC
      REAL GAMPMC,GAMPER
      COMMON /TAUBMC/ GAMPMC(30),GAMPER(30),NEVDEC(30)
C
C     This COMMON is output
C
C     NEVDEC(I) : Number of generated decays for the I-th channel.
C     GAMPMC(I) : Branching ratio for the i-th decay mode.
C     GAMPER(I) : Relative statistical error for the branching ratio
C                 calculated in the Monte Carlo out of the actual
C                 matrix element in the generator.
C
