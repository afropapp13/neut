      subroutine fort_fopen(lun,fname,access,ftype)


      integer lun
      character*(*) fname,access
      character*1   ftype

      common /QUEST/iquest(100)
      common /STORE/ifence(2),lev,ilvect(10000)


C      call mzebra(-1)
C      call mzstor(ixstor,'//',' ',ifence,lev,ilvect(1),ilvect(1),
C     &            ilvect(2000),ilvect(10000))

c
C      CALL CFOPEN(IQUEST(1),0,5670,'w ',0,fname,ISTAT)
     
C      CALL FZFILE(LUN, 5670,'LO')

      if (ACCESS.eq.'WRT') then
         CALL SET_RFLIST(LUN,fname,'LOCAL',' ',ACCESS,' ',' ',
     $                   'recl=5670 status=new','','')
      else
         CALL SET_RFLIST(LUN,fname,'LOCAL',' ',ACCESS,' ',' ',
     $                   'recl=5670 status=old','','')
      endif
     
      CALL SKOPENF(LUN, 1 , ftype, IRET)
      IF (IRET.lt.0) THEN
        write(*,*) "reformat:***ERROR*** CAN NOT OPEN OUTPUT FILE ***"
        call exit(150)
      ENDIF
      call reset_rflist()
       
      return
      end
      
