***********************************************************************
*
*     ------------------------
      SUBROUTINE KZLIMIT(IDIM)
*     ------------------------
*
*     (Purpose)
*       Create a bank structure ZBS (Zebra Bank Structure),
*              variable bank size
*
*     (Input)
*       IDIM : Size of KZBS common block
*
*     (History)
*       Sep. 25, 1998.    First version created.        M.Nakahata
*
***********************************************************************


#include "zbs.h"
#include "zbssys.h"

*     Default values and constants for VDS parameters
*          LBDIR0 ; Initial size of bank directory table
*          LBDIR1 ; Secondary size of bank directory table
*          LINDHD ; Length of index header
*          LINDSG ; Length of index segment
*          LDTRC0 ; Initial size of data record
*          FDTRC1 ; Factor of secondary size of data record
*          NSEGM0 ; Initial no. of segments
*          FSEGM1 ; Factor of secondary no. of segments

      PARAMETER (NPSTRL=1000,NPREFL=100,LLIMD2=10000)
      PARAMETER (LREFL=NPSTRL+1,LWORK=NPSTRL+NPREFL+1)

      PARAMETER (NLSYSD=2000,NMSYSD=10000)

*     NPBANK ; Maximum number of 'banks' in one events
*     NPEVNT ; Maximum number of events

      PARAMETER (NPBANK=200,NPEVNT=100)

      CHARACTER*5 CNZBS/'/ZBS/'/
      CHARACTER*6 CNSYSD/'SYSTEM'/

      CHARACTER*4 CNDIR/'*DIR'/
      PARAMETER (NLDIR=LLNKHD,NSDIR=LLNKHD,NDDIR=LDIRHD)

*          LUN0   ; Default value for LUN for print out

      PARAMETER (LUN0 = 6)

      PARAMETER (NPNAME=8)
      CHARACTER*4 CNNAM(NPNAME)
     &    /'*NA1','*NA2','*NA3','*NA4',
     &     '*NA5','*NA6','*NA7','*NA8'/
      PARAMETER (NLNAM= 0,NSNAM= 0,NDNAM=NPBANK,IONAM=5)

      CHARACTER*4 CNPNT/'*PNT'/
      PARAMETER (NLPNT= 0,NSPNT= 0,NDPNT=NPBANK,IOPNT=2)

      CHARACTER*4 CNEXF/'*EXF'/
      PARAMETER (NLEXF= 0,NSEXF= 0,NDEXF=NPBANK,IOEXF=2)

      CHARACTER*4 CNIOF/'*IOF'/
      PARAMETER (NLIOF= 0,NSIOF= 0,NDIOF=NPBANK,IOIOF=2)

      CHARACTER*4 CNSRT/'*SRT'/
      PARAMETER (NLSRT= 0,NSSRT= 0,NDSRT=NPBANK,IOSRT=2)

      CHARACTER*4 CNIOX/'*IOX'/
      PARAMETER (NLIOX= 0,NSIOX= 0,NDIOX=NPBANK,IOIOX=2)

      CHARACTER*4 CNEVH/'*EVH'/
      PARAMETER (NLEVH= 0,NSEVH= 0,NDEVH=NPBANK,IOEVH=2)

      CHARACTER*4 CNEVC/'*EVC'/
      PARAMETER (NLEVC= 2,NSEVC= 1,NDEVC=1,IOEVC=2)

      PARAMETER (NSYSB=2)

      INTEGER*4 LIST(4)

      INTEGER IDEFLOG
      COMMON /KZLOGL_INTERNAL/ IDEFLOG
      SAVE /KZLOGL_INTERNAL/

C ***** Calculate NPSTOR

      NPSTEX = IDIM - 1 - NPDIV - 1


      LIST(1) = 3
      LIST(2) = IDEFLOG
      LIST(3) = -1
      LIST(4) = -1

      IF(IFZEBI.EQ.1) CALL MZEBRA(LIST)

***** Create ZBS bank structure

      CALL MZSTOR(IXSTOR,CNZBS,' ',IFENCE,LZBS,LZBS(LREFL),LZBS(LWORK),
     &            LZBS(LLIMD2),LZBS(NPSTEX))
*
***** CALL DZSTOR('KZINIT-01',IXSTOR)
*
      IXDIV(1) = IXSTOR+1
      IXDIV(2) = IXSTOR+2

***** Create system division (division 3) for Bank Control

      CALL MZDIV(IXSTOR,IXDIV(3),CNSYSD,NLSYSD,NMSYSD,'L')

      IXSYSD = IXDIV(3)

***** Create a directory bank.

      CALL MZFORM(CNDIR,'8I 2F',IODIR)

      CALL MZBOOK(IXSYSD,LZBS(IPDIRT),LZBS(IPDIRT),1,
     &            CNDIR,NLDIR,NSDIR,NDDIR,IODIR,0)

      LDIR=LZBS(IPDIRT)

      IZBS(LDIR+INBTAB)=NPBANK
      IZBS(LDIR+INBANK)=0
      IZBS(LDIR+INLUN )=LUN0

*     ***** Space for bank name
      DO 10 INAME=1,NPNAME
      CALL MZBOOK(IXSYSD,LNAME,LZBS(IPDIRT),-INAME,CNNAM(INAME),
     &            NLNAM,NSNAM,NDNAM,IONAM,0)
   10 CONTINUE

*     ***** Space for bank pointer#
      CALL MZBOOK(IXSYSD,LPNT,LZBS(IPDIRT),-ILBPNT,
     &            CNPNT,NLPNT,NSPNT,NDPNT,IOPNT,0)

*     ***** Store pointer#
      DO 20 I = 1, NPBANK
      IZBS(LPNT + I) = I + NSYSB
   20 CONTINUE

*     ***** Space for bank existing flag
      CALL MZBOOK(IXSYSD,LEXF,LZBS(IPDIRT),-ILBEXF,
     &            CNEXF,NLEXF,NSEXF,NDEXF,IOEXF,0)

*     ***** Space for output flag
      CALL MZBOOK(IXSYSD,LIOF,LZBS(IPDIRT),-ILBIOF,
     &            CNIOF,NLIOF,NSIOF,NDIOF,IOIOF,0)

*     ***** Space for binary sorting
      CALL MZBOOK(IXSYSD,LSRT,LZBS(IPDIRT),-ILBSRT,
     &            CNSRT,NLSRT,NSSRT,NDSRT,IOSRT,0)

*     ***** Space for I/O format word
      CALL MZBOOK(IXSYSD,LIOX,LZBS(IPDIRT),-ILBIOX,
     &            CNIOX,NLIOX,NSIOX,NDIOX,IOIOX,0)

*     ***** Space for event header flag word
      CALL MZBOOK(IXSYSD,LEVH,LZBS(IPDIRT),-ILBEVH,
     &            CNEVH,NLEVH,NSEVH,NDEVH,IOEVH,0)

***** Create event control bank in division 2.

      IXEVTD = IXDIV(2)

      CALL MZBOOK(IXEVTD,LZBS(IPEVC),LZBS(IPEVC),1,
     &            CNEVC,NLEVC,NSEVC,NDEVC,IOEVC,0)

      RETURN
      END

      BLOCK DATA

#include "zbs.h"
      COMMON /KZLOGL_INTERNAL/ IDEFLOG

      DATA IFZEBI/ 1/
      DATA IDEFLOG / 0/

      END
