***********************************************************************
*
*     rscons.h
*  
*     ( purpose )
*       set constants file for Rein-Sehgal simulation
*  
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgals code)
*       1997.Dec    for eta
*       1998.Feb    for K + LAMBDA
*                   change Wmax 1.4 -> 2.0 GeV
*       2001?     ; MA=1.02 to 1.01
*       2002.09.02; MA=1.01 to 1.21
*       2007.12.02; Selection of MA is added, 1.11 or 1.21
*       2009.02.09; Z = 0.74 to 0.76 (equivalent to g_A = 1.23 to 1.267)
*       2010.10.10; XMA and XMV now set in necard.F(nefillmodel.F)
*       2010.10.13; Modify masses to correspond to mcmass.F
*
***********************************************************************
      REAL  XMN,XMN2,XMPI,XMPI2,XMP,XMP2,XMNE,XMNE2

      REAL  Z
      REAL  PI
      REAL  OMEG
      REAL  XW
      REAL  EPSI
      REAL  XME,XMMU,XMTAU
      REAL  WMAX 
      REAL  XMKAON,XMKAON2,XMLAMD,XMLAMD2
      REAL  XMETA,XMETA2
      REAL  XMGAM,XMGAM2

C      COMMON /RSCONS/Z,PI,MN,MN2,MPI,MPI2,MV,MV2,MA,MA2,OMEG,XW,EPSI,
C     $               ME,MMU,MTAU

c      PARAMETER (Z    = .74)
      PARAMETER (Z    = .76)
      PARAMETER (PI   = 3.1415926)
      PARAMETER (XMP = .938272)
      PARAMETER (XMP2 = XMP**2)
      PARAMETER (XMNE = .939566)
      PARAMETER (XMNE2 = XMNE**2)
      PARAMETER (XMN  = 0.938919)
      PARAMETER (XMN2 = XMN**2)
      PARAMETER (XMPI = .138)
      PARAMETER (XMPI2= XMPI**2)
      PARAMETER (XMETA = .5488)
      PARAMETER (XMETA2= XMETA**2)
      PARAMETER (XMKAON = .496)
      PARAMETER (XMKAON2= XMKAON**2)
      PARAMETER (XMLAMD = 1.116)
      PARAMETER (XMLAMD2= XMLAMD**2)
      PARAMETER (XMGAM = 0.0E0)
      PARAMETER (XMGAM2= 0.0E0)
      PARAMETER (OMEG = 1.05)
      PARAMETER (XW   = .22)
      PARAMETER (EPSI = .00001)
      PARAMETER (XME  = 0.000511)
      PARAMETER (XMMU = 0.105658)
      PARAMETER (XMTAU= 1.77699)
      PARAMETER (WMAX = 2.000)

C--RESONANCE PARAMETERS--------------------------------------------
C-------XEE & FBWNO0 is for SINGLE PION DECAY----------------------

      REAL XMRR(31),BRR(31)
      REAL XEE(31),FBWNOO(31)
      REAL XEEE(31),FBWNOOE(31)
      REAL XEEK(31),FBWNOOK(31)
      REAL XEEG(31),FBWNOOG(31)

      DATA XMRR/1.232,1.52,1.52,1.535,1.535,1.62,1.65,
     1     1.65,1.7,1.7,1.675,1.675,1.7,1.44,1.44,
     2     1.6,1.68,1.68,1.71,1.71,1.72,1.72,1.91,
     3     1.905,1.95,1.92,1.99,1.99,.94,.94,.94/

      DATA BRR/.115,.125,.125,.15,.15,.14,.15,.15,.1,.1,
     1     .155,.155,.25,.2,.2,.37,.125,.125,.11,.11,
     2     .2,.2,.22,.3,.24,.25,.325,.325,3.,3.,3./
cc

ccc for gamma          05/04/30
      DATA XEEG   /0.006,0.005,0.005,0.001,0.001,0.000,0.000,0.000,
     $             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     $             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     $             0.000,0.000,0.000,0.000,0.000,0.000,0.000/
      
      DATA FBWNOOG/1.080,1.014,1.014,1.253,1.253,1.210,1.215,1.215,
     $             1.172,1.172,1.025,1.025,0.842,1.085,1.085,0.930,
     $             0.942,0.942,1.349,1.349,1.246,1.246,1.272,0.700,
     $             0.807,1.242,0.702,0.702,0.000,0.000,0.000/

ccc for K + LAMBDA        98/02/25 J.K.

      DATA XEEK/.0,.0,.0,.0,.0,.0,.11,.11,
     &          .03,.03,.0,.0,.0,.0,.0,.0,
     &          .0,.0,.25,.25,.0,.0,.0,.0,
     &          .0,.0,.0,.0,.0,.0,.0/

      DATA FBWNOOK/1.,1.,1.,1.,1.,1.32,1.24,
     1     1.24,0.586,0.586,0.37,0.37,0.337,1.,1.,
     2     1.,0.353,0.353,0.947,0.947,0.686,0.686,1.0,
     3     0.457,0.558,1.,0.497,0.497,1.,1.,1./

ccc  for eta        98/02/25 J.K.
      DATA XEEE/.0,.0,.0,.4,.4,.0,.0,.0,.0,.0,.0,.0,
     1     .0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     2     .0,.0,.0,.0,.0,.0/

      DATA FBWNOOE/1.,0.277,0.277,1.194,1.194,1.1729,1.126,
     1     1.126,0.917,0.917,0.677,0.677,0.543,1.,1.,
     2     0.45,0.61,0.61,1.25,1.25,1.0,1.0,1.17,
     3     0.548,0.655,1.12,0.574,0.574,1.0,1.,1./

ccc for pion

      DATA XEE/1.,.55,.55,.43,.43,.3,.6,.6,.1,.1,.35,.35,
     1     .15,.6,.6,.1,.6,.6,.15,.15,.15,.15,.22,.12,.4,
     2     .17,.01,.01,1.,1.,1./

      DATA FBWNOO/.96,.98,.98,1.06,1.06,1.05,1.05,1.05,
     1     1.15,1.15,.99,.99,.81,1.04,1.04,.89,.91,.91,
     2     1.32,1.32,1.21,1.21,1.2,.65,.74,1.17,.63,.63,
     3     1.,1.,1./


C---------- NOT USED -----------------------------------------
C      DATA ETAA/1.03,.3,.3,.34,.34,1.6,1.7,1.7,.53,.53,
C     1     .59,.59,.48,.49,.49,.3,.38,.38,.029,.029,
C     2     .054,.054,.072,.54,.22,.125,.31,.31,1.,1.,1./



