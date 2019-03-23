************************************************************************
*     ------------------
*     INCLUDE 'skhead.h'
*     ------------------
*
*       NRUNSK ; run #
*       NSUBSK ; subrun #
*       NEVSK  ; ev # 
*
*       NDAYSK(1) ; year
*             (2) ; month 
*             (3) ; day
*       NTIMSK(1) ; hour
*             (2) ; minute
*             (3) ; second
*             (4) ; 1/100 sec
*
*       NT48SK(1)-(3) ; 48bit clock
*
*       MDRNSK ; run mode
*                0  Monte Carlo
*                1  Normal run
*                2  Laser calib.
*                3  Pedestal data
*                4  Xe lamp calib.
*                5  Nickel calib.
*                6  random trigger run
*                7  Linac calib.
*       IDTGSK ; trigger ID
*                0x01         Low Energy trigger       
*                0x02         High Energy trigger
*                0x04         Super Low Energy trigger (In Normal Run)
*                0x08         Outer Detector trigger
*                             Fission trigger (In Nickel Run)
*                0x10         Periodical trigger firing one of 
*                               a. nothing (null trigger)
*                               b. TQ map laser
*                               c. water attenuation meas. laser
*                               d. Xe ball
*                0x20         After trigger (In Normal Run)
*                             Calibration trigger (In Calibration Run)
*                               Laser trigger
*                               Xe trigger
*                               Ni trigger
*                               Linac trigger
*                0x40         veto start
*                0x80         vato end
*       IFEVSK ; event status flag
*                ATM                  0x00000001
*                TRG                  0x00000002
*                SMP REGISTER         0x00000004
*                SCALER               0x00000008
*                PEDESTAL START       0x00000010
*                PEDESTAL DATA(ATM)   0x00000020
*                PEDESTAL HISTOGRAM   0x00000040
*                PEDESTAL END         0x00000080
*
*                END OF RUN           0x00000100
*                PEDESTAL(ON)         0x00000200
*                GPS DATA             0x00000800
*
*                CAMAC ADC            0x00001000
*                ANTI DATA            0x00002000
*                INNER SLOW DATA      0x00004000
*                RUN INFORMATION      0x00008000
*
*                ERROR (TKO-PS)       0x00010000
*                ERROR (HV-PS)        0x00020000
*                ERROR (TEMPERARTURE) 0x00040000
*
*                UNCOMPLETED ATM DATA 0x00100000
*                INVALID     ATM DATA 0x00200000
*
*                ERROR (DATA)         0x01000000
*                UNREASONABLE DATA    0x02000000
*                LED BURST ON         0x04000000
*
*                INNER DETECTOR OFF   0x10000000
*                ANTI  DETECTOR OFF   0x20000000
*                TRG IS AVAILABLE     0x80000000
*               
*      ---------
*      Contents of /SKHEADA/:
*      --------- 
*      LTCGPS ; Local time clock at last GPS time
*      NSGPS  ; GPS time (sec)
*      NUSGPS ; GPS time (usec)
*      LTCTRG ; Local time clock at TRG
*      LTCBIP ; Local time clock at end of BIP
*      ITDCT0(I); TDC T0 (TRG) time for hut I (I = 1,4)
*      IFFSCC ; FSCC busy flags
*      ICALVA ; Calibration constant version.
*     
*      ---------
*      Contents of /SKHEADG/:
*      --------- 
*      SK_GEOMETRY ; Control geometry definition used
*      (INTEGER*4)         = SK_I  (=1) SUPER-KAMIOKANDE I
*                          = SK_II (=2) SUPER-KAMIOKANDE II
*                          = SK_III(=3) SUPER-KAMIOKANDE III?? 
*
*     (Creation Date and Author)
*       1992.08.26 ; First version by K.S.Hirata
*         95.11.19 ;      modified by Y.Koshio
*         96.01.16 ;      modified by Y.Hayato(add comments:IFEVSK)
*         96.03.30 ;      modified by J. Flanagan (add /SKHEADA/)
*         96.04.12 ;      modified by Y.Hayato(add new comments:IFEVSK)
*         96.05.02 ;      modified by Y.Hayato(add new comments:IFEVSK)
*         03.01.17 ; Added SK geometry version
*
************************************************************************

      INTEGER NRUNSK, NSUBSK, NEVSK, NDAYSK, NTIMSK, NT48SK, MDRNSK,
     &   IDTGSK, IFEVSK
      INTEGER    LTCGPS, NSGPS, NUSGPS, LTCTRG, LTCBIP, IFFSCC,
     &           ITDCT0, ICALVA
      COMMON /SKHEAD/ NRUNSK, NSUBSK, NEVSK,
     &                NDAYSK(3), NTIMSK(4), NT48SK(3),
     &                MDRNSK, IDTGSK, IFEVSK

      COMMON /SKHEADA/ LTCGPS, NSGPS, NUSGPS, 
     &                 LTCTRG, LTCBIP, 
     &                 ITDCT0(4), IFFSCC, ICALVA

      INTEGER*4 SK_GEOMETRY
      COMMON /SKHEADG/ SK_GEOMETRY
