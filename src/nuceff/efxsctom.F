C==================================================================
      real function efxsctom(absp,type)
C==================================================================
C
C     calcurate cross section for each process of omega meson
C
C--------------------------------------------------------------------
      implicit none

C - argument input
      real absp
      integer type

* ==================================================================
* type indicate interaction type
* ------------------------------------------------------------------
*  1. omN - piN
*  2. omN - rhoN
*  3. omN - rhopiN
*  4. omN - pipiN
*  5. omN - omN
*  6. omN - sigmaN
* ==================================================================
      integer i
      integer dnum(6),dd
      data dnum/21,12,12,13,21,21/
      real px
      real crosi1(21),pwi1(21)
      real crosi2(12),pwi2(12)
      real crosi3(12),pwi3(12)
      real crosi4(13),pwi4(13)
      real crosi5(21),pwi5(21)
      real crosi6(21),pwi6(21)
      real pw(21),cros(21)

      data pwi1/10,20,30,40,50,60,70,80,90,100,200,300,400,
     &     500,600,700,800,900,1000,2000,3000/
      data crosi1/238.7,122.7,81.49,59.95,48.85,40.84,35.03,
     &     30.05,26.44,23.87,11.66,7.356,5.275,3.981,3.082,
     &     2.448,1.896,1.585,1.359,0.4642,0.2268/

      data pwi2/100,200,300,400,500,600,700,800,900,1000,2000,3000/
      data crosi2/17.34,18.46,20.41,21.19,20.66,19.66,18.23,16.70,
     &     15.11,13.85,5.766,3.008/

      data pwi3/800,1000,1200,1400,1600,1800,2000,2200,2400,
     &     2600,2800,3000/
      data crosi3/0.1,2.,8.112,10.72,11.50,11.50,11.10,10.00,9.657,
     &     9.326,9.326,9.326/

      data pwi4/10,20,30,40,50,60,70,80,90,100,200,300,400/
      data crosi4/34.69 ,17.43,11.41,8.645,6.904, 5.736,5.025,
     &     4.402,3.857, 3.469, 1.814,1.286, 1.054/

      data pwi5/10,20,30,40,50,60,70,80,90,100,200,300,400,
     &     500,600,700,800,900,1000,2000,3000/
      data crosi5/15.06,15.06,15.06,15.06,15.06,15.06,15.06,
     &     15.06,15.06,14.68,14.13,13.59,13.25,12.75,12.27,
     &     11.81,11.51,11.08,10.80,8.468,6.310/

      data pwi6/10,20,30,40,50,60,70,80,90,100,200,300,400,
     &     500,600,700,800,900,1000,2000,3000/
      data crosi6/30,15.28,9.412, 7.613, 5.975, 4.982, 4.281,
     &     3.793, 3.360, 3.068, 1.528, 1.000,0.7386, 0.5796,
     &     0.4833, 0.3910, 0.3360,0.2801, 0.2336, 0.06543, 0.03360/

*
*------- program start
      dd=dnum(type)

      px=log10(absp)
      do i=1,dd
         if(type.eq.1)then
            cros(i)=log10(crosi1(i))
            pw(i)=log10(pwi1(i))
         else if(type.eq.2)then
            cros(i)=log10(crosi2(i))
            pw(i)=log10(pwi2(i))
         else if(type.eq.3)then
            cros(i)=log10(crosi3(i))
            pw(i)=log10(pwi3(i))
         else if(type.eq.4)then
            cros(i)=log10(crosi4(i))
            pw(i)=log10(pwi4(i))
         else if(type.eq.5)then
            cros(i)=log10(crosi5(i))
            pw(i)=log10(pwi5(i))
         else if(type.eq.6)then
            cros(i)=log10(crosi6(i))
            pw(i)=log10(pwi6(i))
         endif
      enddo

      if(px .lt. pw(1)) then
         efxsctom=cros(1)
         goto 22
      endif

      do i=1,dd-1
         if((px .ge. pw(i)) .and. (px.lt.pw(i+1)))then
            efxsctom=(cros(i+1)-cros(i))/(pw(i+1)-pw(i))*(px-pw(i))
     &           +cros(i)
            goto 22
         endif
      enddo

      if(px .ge. pw(dd)) then
         efxsctom=cros(dd)
         goto 22
      endif

      print *,'error'
      efxsctom=999

 22   continue
      efxsctom=10**efxsctom

      return
      end


