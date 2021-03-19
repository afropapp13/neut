************************************************************************
*     ----------------------------
      SUBROUTINE KZREAD (LUN,IERR)
*     ----------------------------
*
*     (Purpose)
*       Read one event from sequential file and restore it in the
*       Venus Data Structure.
*
*     (Input)
*       LUN    : Logical unit# of the input device.
*
*     (Output)
*       IERR =  1 : Find EOF while reading data.
*       IERR =  2 : Ilegal sequence of record occured while reading.
*
*     (Comment)
*       For the format of data in a seqential file, see the comment
*       in SUBROUTINE ZVWRIT.
*
*     (History)
*       1991.12.25 : First version by J.Kanzaki
*       1996.04.25 :  J.Kanzaki, Add declarations.
*       1998.01.14 : K.Inoue for rfa
*
************************************************************************

#include "zbsdeclr.h"
#include "zbssys.h"
#include "zbs.h"
#include "zbsinf.h"

      INTEGER LUN,IERR

      INTEGER NUH0
      PARAMETER (NUH0=100)
      INTEGER IUHEAD(NUH0)

#ifdef SITE_HAS_RFA
      integer ichange_file,ichange
      character*160 fname,unit,lgname,access,stg,multi,comment,hname
      character*160 mode,file_name_now
      common /FILE_NAME/ichange,file_name_now
      data ichange_file/1/
      data ichange/0/
      save ichange_file
#endif

      INTEGER IXEVTD
      INTEGER LEVC,LEMB,LLAST

      INTEGER NUH

      INTEGER LZLAST

#ifdef SITE_HAS_RFA
      integer nptcom,icncom,iftype,itpcom,iszcom,iret
      common/luncom/nptcom(99),icncom(99),iftype(99),
     &              itpcom(99),iszcom(99)
#endif

***** Initial set up.
      IERR   = 0
#ifdef SITE_HAS_RFA
      ichange=ichange-1
#endif

 1    IXEVTD = IXDIV(2)

***** Get pointer to the event control bank
      LEVC = LZBS(IPEVC)
      LLAST = LZLAST(IXSTOR,LZBS(LEVC-ILEMLS))

      NUH = NUH0

      IF(LLAST.NE.0) THEN

         CALL FZIN(LUN,IXEVTD,LLAST,0,' ',NUH,IUHEAD)
         LEMB = LZBS(LLAST)

      ELSE

         CALL FZIN(LUN,IXEVTD,LEVC,-1,' ',NUH,IUHEAD)
         LEMB = LZBS(LEVC-ILEMLS)

      END IF

* Save event size information (Added by M.Nakahata)

      ISZEVT = IQUEST(14)

#ifdef SITE_HAS_RFA

      if(iquest(1).lt.0)then
         ierr=2
         return
      endif
      if(iquest(1).eq.2)then
        ichange_file=1
        goto 1
      endif
      if(iquest(1).ge.5)then
         if(icncom(lun).eq.1)then
           call fzendi(lun,'TX')
           call skopenf(lun,0,'Z',iret)
           if(iret.eq.0)goto 1
         endif
         ierr=1
         return
      endif
      if(iquest(1).ne.0)goto 1

      if(ichange_file.eq.1)then
        ichange_file=0
        if(mod(iftype(lun),100).ne.0)then
          call rgetv(lun,nptcom(lun),fname,unit,lgname,
     &           access,stg,multi,comment,hname,mode)
        else
          call rgetv(lun,0,fname,unit,lgname,
     &           access,stg,multi,comment,hname,mode)
        endif
        file_name_now=fname
        ichange=20
      endif
#else
      IF(IQUEST(1).GE.2) THEN
         IERR = 1
         RETURN
      END IF
#endif

      LZBS(LEVC-ILEMCR) = LEMB

C     CALL KZCLDT

      CALL KZMKDT

      RETURN
      END
