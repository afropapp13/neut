      subroutine ucase(a)
      character*(*) a
      do 100 i=1,len(a)
        if(ichar(a(i:i)).ge.ichar('a') .and.
     &     ichar(a(i:i)).le.ichar('z'))then
           a(i:i)=char(ichar(a(i:i))+ichar('A')-ichar('a'))
        endif
  100 continue
      return
      end
