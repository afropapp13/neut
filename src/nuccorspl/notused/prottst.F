      program prottst
      common /nucstack/ntrack,
     $     pptrack(4,20),chtrack(20),amass(20),stpt(3,20)

      common /nuchad/ihap,nel(3)
      common /hadron/inot,itrans,idec,ine1,ine2,iel,iabs,icx
      dimension stpi(3),pi(4)
      dimension po(4,20),ido(20),stpo(3,20)

      integer npawc
      parameter (npawc=50000)
      real h
      common/pawc/h(npawc)
      call hlimit(npawc)
      call hbook1(1,'iel',  20,0.,10.,0.)
      call hbook1(2,'spp',  20,0.,10.,0.)
      call hbook1(3,'dpp',  20,0.,10.,0.)
      call hbook1(4,'no',   20,0.,10.,0.)
      call hbook1(5,'trans',20,0.,10.,0.)
      call hbook1(6,'ntrac',20,0.,10.,0.)
      call hbook1(7,'ido'  ,20,0.,20.,0.)

      call hbook1(10,'ep'   ,100,0.,1000.,0.)
      call hbook1(11,'en'   ,100,0.,1000.,0.)
      call hbook1(12,'epi0' ,100,0.,1000.,0.)
      call hbook1(13,'epi+' ,100,0.,1000.,0.)
      call hbook1(14,'epi-' ,100,0.,1000.,0.)

      call hbook1(15,'pang',360,0.,360.,0.)

      call hbook1(16,'pangm',360,0.,360.,0.)
      call hbook1(17,'epmax',100,0.,1000.,0.)
      

      inot=0
      izer =0
      icx=0
      nel1=0
      nel2=0
      nel3=0
      anuc = 16.
      imode = 2
      icont = 0
      do i=1,3
         stpi(i)=0
      end do
      idmc = 2212               ! 2112 for n
      print '(a,$)','proton momentum MeV   '
      read(5,*)pp

      pi(1)=0.
      pi(2)=0.
      pi(3)=pp

      pi(4)= sqrt(939.**2+pi(1)**2+pi(2)**2+pi(3)**2)

      do i=1,100000
         if (mod(i,10000).eq.0) print *,' got ',i,' events'
         call nrprton(anuc,stpi,pi,idmc,
     1                 po,ido,stpo,no,
     2                        imode,icont)

         if(iabs.eq.-1) izer = izer+1
         if(nel(1)+nel(2)+nel(3).eq.0) icx = icx+1
         if(nel(1).ne.0)nel1=nel1+1
         if(nel(2).ne.0)nel2=nel2+1
         if(nel(3).ne.0)nel3=nel3+1

         if(nel(1).ne.0) call hfill(1,float(nel(1)),0.,1.)
         if(nel(2).ne.0) call hfill(2,float(nel(2)),0.,1.)
         if(nel(3).ne.0) call hfill(3,float(nel(3)),0.,1.)
         call hfill(4,float(no),0.,1.)
         if(itrans.ne.0) call hfill(5,float(itrans),0.,1.)
         call hfill(6,float(ntrack),0.,1.)
         ppmax=0
         itmax=0
         do j=1,no
            if(ido(j).eq.2212 .and. po(4,j).gt.ppmax) then
               ppmax = po(4,j)
               itmax = j
            end if
         end do
         do j=1,no
            if(ido(j).eq.2212) then
               nid = 15
               call hfill(10,po(4,j)-939.,0.,1.)
               pnorm = sqrt(po(1,j)**2+po(2,j)**2+po(3,j)**2)
               ang = acos(po(3,j)/pnorm)/0.0174533
               if(ang.gt.1.) call hfill(15,ang,0.,1.)
               if(j.eq.itmax) then
                  call hfill(17,po(4,j)-939.,0.,1.)
                  if(ang.gt.1.) call hfill(16,ang,0.,1.)
               end if
            else if(ido(j).eq.2112) then
               nid = 10
               call hfill(11,po(4,j)-939.,0.,1.)
            else if(ido(j).eq.111) then
               nid = 5
               call hfill(12,po(4,j)-139.6,0.,1.)
            else if(ido(j).eq.211) then
               nid = 6
               call hfill(13,po(4,j)-139.6,0.,1.)
            else if(ido(j).eq.-211) then
               nid = 4
               call hfill(14,po(4,j)-139.6,0.,1.)
            else
               print *, 'unknown particle ',ido(j)
            end if
            call hfill(7,float(nid),0.,1.)
         end do
      end do
      call hrput(0, 'prttst.hbk', ' ')
      print *,inot,icx,izer
      print *,nel1,nel2,nel3
      end
