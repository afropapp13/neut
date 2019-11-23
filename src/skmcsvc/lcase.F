      subroutine lcase(a)
      character*(*) a
      do 100 i=1,len(a)
        if(ichar(a(i:i)).ge.ichar('A') .and.
     &     ichar(a(i:i)).le.ichar('Z'))then
           a(i:i)=char(ichar(a(i:i))+ichar('a')-ichar('A'))
        endif
  100 continue
      return
      end
