C     $Id: getvar.f,v 1.6 1998/04/04 20:57:11 bviren Exp $
C
C     $Log: getvar.f,v $
C     Revision 1.6  1998/04/04 20:57:11  bviren
C     Iolib overhaul:
C     Some Linux incompatibilities and generic bugs fixed (for the 2nd or
C     3rd time). Many missing header files added (for 2nd time). Buffer
C     overrun bugs in charterm.c and gethostname.c fixed. Some unused
C     variables removed.
C
C
      subroutine getvar(a,var,val)
      character*1024 b
      character*(*) a,var,val
      character*1024 tmp
      if(len(a).gt.1024)then
        write(6,*)'warning length greater than 1024 in getvar'
      endif
      val=' '
      b=a
      call lcase(b)
      j=1
      i=1
      goto 200
   10 continue
      do 100 i=j,1024
        if(b(i:i).eq.' ')goto 200
  100 continue
      goto 1000
  200 continue
      do 300 j=i,1024
        if(b(j:j).ne.' ')goto 400
  300 continue
      goto 1000
  400 continue
      tmp(1:lenchr(var))=var(1:lenchr(var))
      tmp(lenchr(var)+1:lenchr(var)+1)='='
      if(b(j:j+lenchr(var)).ne.tmp(1:1+lenchr(var)))goto 10
      do 500 i=j+lenchr(var)+1,1024
        if(b(i:i).eq.' ')goto 1000
        val(i-j-lenchr(var):i-j-lenchr(var))=a(i:i)
  500 continue
 1000 continue
      return
      end
