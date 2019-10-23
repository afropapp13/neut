************************************************************************
*     ------------------
*     INCLUDE 'dsbnkf.h'
*     ------------------
*
*     Flags for disk-writing of individual banks
*     and logical unit numbers for inputs and outputs
*
*     (Creation Date and Author)
*       1993.04.07 ; First version by K.S.Hirata
*
************************************************************************

      INTEGER  NBNKDS, LIN, LOUT, IFLGV

      COMMON /DSBANK/ NBNKDS(6)
      LOGICAL IMCPHI, IMCRI, IMCRAI, IMCTQI, ITQI, ITQAI
      COMMON/INBANK/ IMCPHI, IMCRI, IMCRAI, IMCTQI, ITQI, ITQAI

      COMMON /DSLUNS/LIN, LOUT, IFLGV
