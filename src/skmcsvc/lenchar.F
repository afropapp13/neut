      integer function lenchr(a)
      character*(*) a
      do 10 lenchr=len(a),1,-1
        if(a(lenchr:lenchr).ne.' ')goto 20
   10 continue
      lenchr=0
   20 continue
      return
      end
