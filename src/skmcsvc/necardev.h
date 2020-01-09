************************************************************************
*     ----------
*     NECARDVC.H
*     ----------
*
*     (Purpose)
*       COMMON for CARD on EVCT to NEUT
*
*     (Variables)
*       NEVTEVCD   : ( EVCT-NEV  ) Number of Events
*
*       IDPTEVCD   : ( EVCT-IDPT ) Particle code (neutrino)
*
*       MPOSEVCT   : ( EVCT-MPOS ) Flag on Vertex
*                     1 : Fixed
*                     2 : Random
*       POSEVCT(3) : ( EVCT-POS ) Vertex Position,
*                          used when MPOS = 1
*       RADEVCT    : ( EVCT-RAD ) maximum Radius of vertex,
*                          used when MPOS = 2
*
*       MDIREVCT   : ( EVCT-MDIR ) Flag on Direction
*                     1 : Fixed
*                     2 : Random
*       DIREVCT(3) : ( EVCT-DIR ) Direction if MDIREVCT = 1
*
*       MPVEVCT    : ( EVCT-MPV ) Flag on Momentum
*                     1 : Fixed  ( PEVCT(1) )
*                     2 : Random ( Momentum = [PEVCT(1),PEVCT(2)] )
*       PVEVCT(2)  : ( EVCT-PV  ) Momentum
*
*     (Creation Date and Author)
*       2007.01.02 ; Y.Hayato
*
************************************************************************
      INTEGER         NEVTEVCT,IDPTEVCT,MPOSEVCT,MDIREVCT,MPVEVCT,
     $                INMEVEVCT, OUTPUTFORMAT
      REAL            POSEVCT(3),RADEVCT,DIREVCT(3),PVEVCT(2), ANGDEG,ANGWIDTH
      character*80    FILENMEVCT,HISTNMEVCT

      COMMON/NEVCCARD/NEVTEVCT, IDPTEVCT,
     $                MPOSEVCT, POSEVCT, RADEVCT,
     $                MDIREVCT, DIREVCT,
     $                MPVEVCT,  PVEVCT, ANGDEG,ANGWIDTH,
     $                FILENMEVCT, HISTNMEVCT, INMEVEVCT, OUTPUTFORMAT
