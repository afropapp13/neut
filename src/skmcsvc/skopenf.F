      subroutine skopenf(lun,npt0,flg,ihndl)
c
c  lun    logical unit number linked to RFLIST
c           for fortran open, lun is fortran-LUN
c
c  npt    open npt-th file listed in RFLIST as specified LUN
c         if zero, file concatenation works with rread
c
c  flg = 'f' or 'F'     open for fortran access
c        'c' or 'C'     open for C read/write
c        'r' or 'R'     open for rread/rwrite
c        ' ','z' or 'Z' open for ZEBRA file
c
c ((return value))
c  ihndl  file hundler if specified C read/write
c            (if zero or negative that means open error)
c         if file not specified in RFLIST return -1001
c
      character*160 fname,unit,lgname,access,stg,multi,
     &     comment,hname,mode
      character*160 comnt,f_acc,f_form,f_recl,f_stat,f_unit,f_type
      character*1  flg
      character*1  amode,umode,iof
      external iolib_get_tapemode
      integer*4  iolib_get_tapemode
      integer*4  itapemode
      character*1  tapemode
      common/quest/iquest(100)
      common/luncom/nptcom(99),icncom(99),iftype(99),
     &              itpcom(99),iszcom(99)
      data icall/0/
c
      if(lun.le.0 .or. lun.ge.100)then
        write(6,*)'skopenf: bad lun',lun,npt
        stop
      endif
c
      if(icall.eq.0)then
        icall=1
        do 10 i=1,99
          nptcom(i)=0
          icncom(i)=0
          itpcom(i)=0
          iszcom(i)=0
   10   continue
      endif
c
      npt=npt0
      izflg=0
      ihndl=0
c
      if(npt.lt.0)then
        write(6,*)'skopenf: ',lun,npt,' ??'
        stop
      endif
      if(nptcom(lun).gt.0)then
        if(npt.ne.0 .or. icncom(lun).eq.0)then
          write(6,*)'skopenf: LUN already open',lun,npt
          stop
        endif
      endif
      if(npt.eq.0)then
        nptcom(lun)=nptcom(lun)+1
        npt=nptcom(lun)
        icncom(lun)=1
        icont=1
      else
        nptcom(lun)=npt
        icont=0
      endif
c
      call rgetv(lun,npt,fname,unit,lgname,
     &           access,stg,multi,comment,hname,mode)
c     write(6,*)fname,unit,lgname,access,stg,multi,comment,hname,mode
      if(fname.eq.' ')then
c file not found in RFLIST
        ihndl=-1001
        return
      endif
      call ucase(unit)
      call ucase(access)
      call ucase(stg)
      call ucase(multi)
      comnt=comment
      call lcase(comnt)
C     disable remote/rfa access   2019.2.14  KO 
C      if(unit.eq.'LOCAL')then
C        if(hname.eq.'rfa')then
C          ihost=2                 ! local w/ rfa
C        else
C          ihost=1                 ! local
C        endif
C      else
C          ihost=0                 ! remote access
C      endif
      ihost=1
       
      call getvar(comnt,'unit',f_unit)
      call getvar(comnt,'recl',f_recl)
      nrecl=0
      if(f_recl.ne.' ')read(f_recl,*)nrecl
      call getvar(comnt,'access',f_acc)
      call getvar(comnt,'status',f_stat)
c
c FORTRAN OPEN
      if(flg.eq.'f' .or. flg.eq.'F')goto 1001
c C OPEN
      if(flg.eq.'c' .or. flg.eq.'C')goto 1002
c ROPEN
      if(flg.eq.'r' .or. flg.eq.'R')goto 1003
c ZEBRA
      if(flg.eq.' ' .or. flg.eq.'z' .or. flg.eq.'Z')goto 1004
      write(6,*)'skopenf: ',flg,' ??',lun,npt
      stop
c
c ======================================================================
c  actual open for fortran file
 1001 continue
      ihost=1
      if(icont.eq.1 .and. izflg.eq.0)then
        write(6,*)'skopenf: npt=0 with F-opt ??',lun,npt0
        stop
      endif
      call getvar(comnt,'form',f_form)
c     On Alpha crashes at run-time when status or access are ' '.
c     'UNKNOWN' and 'SEQUENTIAL' are the defaults from documentation.
c     --Tomasz.
      if (f_stat.eq.' ') f_stat='UNKNOWN'
      if (f_acc .eq.' ') f_acc ='SEQUENTIAL'
      if(f_acc.eq.'direct')then
        open(lun,file=fname,access=f_acc,
     &       form=f_form,recl=nrecl,status=f_stat)
      else
        open(lun,file=fname,access=f_acc,form=f_form,status=f_stat)
      endif
      if(izflg.eq.1)goto 2000
      iftype(lun)=1
      return
c  actual open for C file
 1002 continue
      ihost=1
      if(icont.eq.1)then
        write(6,*)'skopenf: npt=0 with C-opt ??',lun,npt0
        stop
      endif
      call ccopnf(fname,f_stat,f_acc,f_unit,ihndl)
      iftype(lun)=10
      return
c  actual open for rfile
 1003 continue
      if(icont.eq.1)
     &  write(6,*)'skopenf: file concatenation for rread',lun,npt0
      call skropen(lun,npt0,icode)
      if(icode.ne.0)stop
      iftype(lun)=100
      return
c  actual open for ZEBRA file
 1004 continue
      izflg=1
      call getvar(comnt,'type',f_type)
c  branch for various file type
      if(f_type.eq.'fortran')then
        if(icont.eq.1)
     &     write(6,*)'skopenf: file concatenation for F-ZEBRA',lun,npt
        goto 1001
      endif
      if(f_type.eq.'remote')goto 3001
      if(f_type.eq.' ' .or. f_type.eq.'c')goto 3002
      write(6,*)'skopenf: bad file_type',lun,npt
      stop
c return from fortran ZEBRA open
 2000 continue
      call fzfile(lun,nrecl,'XIO')
      iftype(lun)=1001
      if(f_stat.eq.'new' .or. f_stat.eq.'unknown' .or.
     &   f_acc.eq.'append')iftype(lun)=11001
      return
c
c ------------------------------------------------------------------------
c calling modified ZEBRA open
c
c for remote ZEBRA
 3001 continue
      medium=3
      umode='K'
      if(access.eq.'RED')then
        amode='r'
        iof='I'
      else
        amode='w'
        iof='O'
      endif
      if(icont.eq.1)then
        write(6,*)'skopenf: file concatenation for remote-ZEBRA',lun,npt0
        icncom(lun)=0
      endif
      iuni=npt0*100+lun
      call cfopen(iuni,medium,nrecl,amode,0,fname,ist)
      if(ist.ne.0)stop
      iquest(1)=iuni
      call fzfile(lun,nrecl,umode//'TX'//iof)
      iftype(lun)=1100
      if(iof.eq.'O')iftype(lun)=11100
      return
c for C ZEBRA
 3002 continue
      if(ihost.eq.0 .or. ihost.eq.2)then
        write(6,*)'skopenf: changed to remote file access ',lun,npt
        goto 3001
      endif
      if(icont.eq.1)
     &   write(6,*)'skopenf: file concatenation for C-ZEBRA',lun,npt
c  check read/write
      amode=' '
      umode='L'
      if(access.eq.'RED' .or. f_stat.eq.'old')then
        amode='r'
        iof='I'
      elseif(access.eq.'WRT'     .or. access.eq.'RPL' .or.
     &       f_stat.eq.'unknown' .or. f_stat.eq.'new')then
        amode='w'
        iof='O'
      endif
      if(access.eq.'APP' .or. f_acc.eq.'append')then
        amode='a'
        iof='O'
      endif
      if(access.eq.' ' .and. f_stat.eq.' ' .and. f_acc.ne.'append')then
        amode='w'
        iof='O'
      endif
      if(amode.eq.' ')then
        write(6,*)'skopenf: cannot identify read/write',lun,npt
        stop
      endif
      medium=0
      call cfopen(iuni,medium,nrecl,amode,0,fname,ist)
      if(ist.ne.0)stop
      iquest(1)=iuni
      itapemode=iolib_get_tapemode(fname)
      call itoch(itapemode,tapemode,101)
 101  continue
      call fzfile(lun,nrecl,umode//'X'//iof//tapemode)
      iftype(lun)=1010
      if(iof.eq.'O')iftype(lun)=11010
      return
      end
