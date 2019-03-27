C************************************************
C  $Id: partnuc.F,v 1.1 2007-01-20 07:13:07 skrep Exp $
C
C  This is the kinematics used by IMB (west cost).
C
C  $Log: not supported by cvs2svn $
C  Revision 1.11  2000/03/01 08:55:37  cmauger
C  Merge ichikilo and main branches.  All changes to the main branch
C  are cosmetic.
C
C  Revision 1.10  1999/08/13 21:30:26  bviren
C  Change `type' to `print' and sind()/cosd() to sin()/cos() to appease g77.
C
C  Revision 1.9.4.1  1999/04/16 05:43:43  mcgrew
C  "PRINT *" is fortran, "TYPE *" is not.
C
C  Revision 1.9  1998/03/19 03:08:57  cmauger
C  Introduced a fudge factor of 3.0 in the calculation of the cross
C  sections in nuckaon.  The cross sections now somewhat match what
C  is in the particle data book.
C
C  Revision 1.8  1998/03/18 14:17:57  cmauger
C  A couple of bug fixes.  I had incorrectly entered the cms energy.  Also,
C  there was a value of pi missing from the cross section calculation.
C
C  Revision 1.7  1998/03/18 11:08:28  cmauger
C  Added subroutines nucknknp0 and nucknknp1 which are for kaons with
C  strangeness +1.  Added the necessary switch in the subroutine nuckaon.
C  The partial-wave amplitudes come from Hyslop et al. 1992.  Phys. Rev. D
C  Volume 46, No. 3.  p. 961.
C  This is only valid up to CMS energies of 1.8 GeV or so.
C  There is still only elastic scattering and charge exchange.
C
C  Revision 1.6  1998/03/04 20:47:47  cmauger
C  Added conditional assignments to proton and kaon masses in nuckaon
C  so nuckaon can be called without using mcinit.  It is a workaround
C  in order to use nuckaon in the old-onsite MC.
C
C  Revision 1.5  1998/03/02 21:44:36  bviren
C  Somehow this mod never made it into the repo. This adds the ability to
C  track which interactions occured via the `partnucint' common. See
C  partnucint.h for C interface.
C
C  Revision 1.4  1997/10/26 09:43:55  bviren
C  The pizero absorbtion is bassed on charged pion absorbtion on
C  deuteron. This crossection is supposed to be adjusted by the ratio of
C  neutral to charged pion absorbtion due to isospin symmetry (0.5) as
C  well as by the ratio of oxygen nuclear density to deuteron density
C  (RHON/RHOD). Somehow this adjustment was lost. It is not known, when.
C
C  Revision 1.3  1997/10/17 13:53:59  habig
C  Changes to compile under absoft f77.  Mostly syntactical, such as moving data
C  statements to the end of the declarations.
C
C  Revision 1.2  1997/06/02 22:08:24  mcgrew
C  Changes to sync skrep with sksuny.
C
C  Revision 1.1  1996/10/06 01:11:57  fgoebel
C  This is basically the IMB (west) kinematics program. A few
C  things have been changed to make it compatible with the SuperK
C  analysis software.
C
C************************************************


C #define sind(x) (sin((x)*0.01745))
C #define cosd(x) (cos((x)*0.01745))


C**************************************************************8
C     @(#)partnuc.f     1.8 modified on 1/28/93
C
#define JPION 0
#define JETA 1
#define JOMEGA 2
#define JKAON 3
#define JNUCLEON 100
#define sind(x) sin((x)*0.0174532)
#define cosd(x) cos((x)*0.0174532)

      subroutine partnuc ( isetn,anuc,
     $                     stpt,
     $                     ei,pi,ui,chi,ni,
     $                     eo,po,uo,cho,
     $                     start,decp,
     $                     imode,icont)
      integer i1, i2, i3, stepnum
      integer isetn
      real anuc
C**** The input energy, momentum, mass and charge
      real ei(20),pi(3,20),ui(20),chi(20)
C**** The number of input particles.
      integer ni
C**** The Number of output particles.
      integer no
C**** The output energy, momentum, mass, and charge.
      real eo(20),po(3,20),uo(20),cho(20)
C**** The starting point for the hadron.  This is an output variable for
C     icont equal 0 and an input variable if icont is 1.
      real start(3)
C**** The final position of the hadron.  This is used to restart a new
C     hadron from the stopping point after a decay.
      real decp(3)
C**** The starting point for each hadron in the stack.  This should be
C     zero unless icont is 1
      real stpt(3,20)
C**** The density mode of the model
      integer imode
C**** Flag if this hadron is being restarted
      integer icont
      real xtr(3)
C**** The nuclear interaction parameters
      real pin(3), upin(3)
      common /nucint/uin,chint,ein,pin,upin
C**** The nuclear center of mass parameters.
      real up4cms(3), beta(3)
      logical prot
      common /nuccms/rhon,pfermi,unucl,ecms2,up4cms,beta,prot
      common /nuchad/ihap,nel,ncx
#include "hadron.h"
#include "nucnewprd.h"
#include "partnucint.h"
C**** The mc step, and the nuclear maximum radius.
      real step, rmax
      parameter (step=.2, rmax=6.)
C**** A flag that the routine is called.
      integer called
      data called/0/

      if(called.eq.0)then
         called=1
         call nucinit(anuc)
      end if

C**** Make sure there are hadrons to track.
      if (ni.lt.1) then
         write (6,*) 'partnuc_: initialized'
         call flush(6)
         return
      end if

C**** Generate a starting position.
      call nucstp(start,icont,imode)

      i1=0
      nx=ni
      do while (i1.lt.nx)
         i1=i1+1

         stepnum = 0

         if(i1.le.ni)then
            ihap=-1
            nel=-1
            ncx=-1
         end if

         iin=0
         uo(i1)=ui(i1)
         uin=ui(i1)*1000.
         ein=ei(i1)*1000.
         chint=chi(i1)
         pinss=0
         do i2=1,3
            if(i1.gt.ni)start(i2)=0.
            xtr(i2)=stpt(i2,i1)+start(i2)
            pin(i2)=pi(i2,i1)*1000.
            pinss=pinss+pin(i2)**2
         end do

         jpar=-1
         if(abs(uin-137.).lt.5.)jpar=JPION !pion
         if(abs(uin-549.).lt.5.)jpar=JETA !eta
         if(abs(uin-783.).lt.5.)jpar=JOMEGA !omega
         if(abs(uin-939.).lt.5.)jpar=JNUCLEON !nucleon
         if(abs(uin-495.).lt.5) jpar=JKAON !kaon
         if(jpar.lt.JPION)go to 99
         if(jpar.eq.JNUCLEON.and.ein.lt.1100.)go to 99

         pinm=sqrt(pinss)
         do i2=1,3
            upin(i2)=pin(i2)/pinm
         end do
         loop=0

C initialize storing which interactions happened.
            do i2=1,100
               pninter(i2,i1) = 0
            end do
            pnnumber(i1)=0



c****         write(6,'('' % tracking '',i3,i3,5f8.1,f5.1)')i1,jpar,
c****     $        ein,pin(1),pin(2),pin(3),uin,chint

 200     continue               !****start tracking***********

         r=0
         do i2=1,3
            xtr(i2)=xtr(i2)+upin(i2)*step
            r=r+xtr(i2)**2
         end do
         r=sqrt(r)

         stepnum = stepnum + 1

C**** Check if the hadron is still inside the nucleus.
         if(imode.eq.1)then
            if(iin.eq.0.and.r.gt.rmax)then
               loop=loop+1
               if(loop.gt.3)go to 99
               az=rmax**2-xtr(1)**2-xtr(2)**2
               if(az.lt.0.)go to 99
               xtr(3)=-sqrt(az)-0.18+0.05*loop
               go to 200
            end if
            if(r.le.rmax)iin=1
         end if
         if(r.gt.rmax)then
            if(i1.le.ni)then
               itrans=itrans+1
               if(ihap.lt.0)inot=inot+1
               if(nel.gt.0)iel=iel+1
               if(ncx.gt.0)icx=icx+1
            end if
            go to 99
         end if

         call nucfermi(r)

         pninttype=0
         ides=i1*10
         if(jpar.eq.JNUCLEON)then
            call nucnuc(ides)
         else if((jpar.eq.JETA).or.(jpar.eq.JOMEGA))then
            icont=1
            call nucvmes(jpar,ides)
         else if (jpar.eq.JKAON) then
            call nuckaon(ides,step)
         else if (jpar.eq.JPION) then
            call nucpion(ides)
         else
            print *, 'PARTNUC: bogus particle in partnuc'
            stop
         end if

         call printinttype

#ifdef DEBUG
         write (6,*) i1,stepnum,pninttype,chint,uin,ides
#endif
         pninter(stepnum,i1) = pninttype
         pnnumber(i1) = stepnum


         if(ides.eq.1)go to 200         ! Step finished w/nothing special

         if(ides.lt.0)then              ! Second particle created
            ia=abs(ides)
            do l=1,ia
               nx=nx+1
               ei(nx)=piso(4,2*l)/1000.
               ui(nx)=uiso(2*l)/1000.
               chi(nx)=chiso(2*l)
               do l1=1,3
                  pi(l1,nx)=piso(l1,2*l)/1000.
                  stpt(l1,nx)=xtr(l1)
               end do
            end do
            if(jpar.eq.JNUCLEON.and.ein.lt.1100.)then
               itrans=itrans+1
               if(ihap.lt.0)inot=inot+1
               if(nel.gt.0)iel=iel+1
               if(ncx.gt.0)icx=icx+1
               go to 99
            end if
            go to 200
         end if

 99      continue               ! Come here to finish particle
         do i2=1,3
            decp(i2)=xtr(i2)
            po(i2,i1)=pin(i2)/1000.
         end do
         eo(i1)=ein/1000.
         cho(i1)=chint
         uo(i1)=uin/1000.
#ifdef DEBUG
         write (6,*)'Loop end: ',i1,stepnum,pninttype,chint,uin
#endif
      enddo                     ! Finish loop over particles

      ni=nx
      call flush(6)
      return
      end

C ****************************************************************

      subroutine nucfermi(r)
      common /nucdis/probmi(31),probdi(31),rhotab(30),pnorm
      common /nucint/uin,chint,ein,pin(3),upin(3)
      common /nuccms/ rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
      logical prot
      dimension rec(3),p4cms(3)

      parameter (toten2=931.**2)
      parameter (rbin=.2, rmin=0., rmax=6., nbin=(rmax-rmin)/rbin)

      common /nucpart/prat
      data prat/0.5/

c     decide if interaction is on proton (prot=.true.) or neutron (.false.)
      if(ranf().lt.prat)then
         prot=.true.
      else
         prot=.false.
      end if
C     NUCLEAR DENSITY AT THIS VALUE OF R
C     INCLUDE 0.1 FACTOR is TO CHANGE FROM MB TO FMSQ (cms is in fm!)
      call nuchis(r-rmin,rbin,30,irlo,irhi,rrem)
      RHON=(RHOTAB(IRHI)*RREM+RHOTAB(IRLO)*(1.0-RREM))*0.1
      RHONF=RHON/PNORM

C     FERMI MOMENTUM AT THIS VALUE OF R
      PFERMI=270.0*((RHONF*10.)**0.3333)
      IF(PFERMI.LE.0.0) PFERMI=1.0

      AREC=PFERMI*ranf()**0.3333
      call ranve(arec,rec,2.,-1.)
      unucl=sqrt(toten2-arec**2)

C     ENERGY AND 3 MOMENTUM PARTS OF CMS 4 MOMENTUM
      E4CMS=EIN+931.
      DO i2=1,3
         P4CMS(i2)=PIN(i2)+REC(i2)
      END DO
      P4CMSM=SQRT(P4CMS(1)**2+P4CMS(2)**2+P4CMS(3)**2)
      DO i2=1,3
         UP4CMS(i2)=P4CMS(i2)/P4CMSM
         BETA(i2)= -P4CMS(i2)/E4CMS
      END DO
C     SQUARE OF CMS ENERGY
      ECMS2=E4CMS**2-p4cmsm**2

      return
      end

C Table look up allowing simple interpolation in a hand done manner
      subroutine nuchis(x,rbin,num,ilo,ihi,erem)
      ilo=(x)/rbin+1
      ihi=ilo+1
      erem=(x-float(ilo-1)*rbin)/rbin
      if(ilo.ge.num)then
         ilo=num
         ihi=num
         erem=1
      end if
      if(ilo.lt.1)then
         ilo=1
         ihi=2
         erem=0.
      end if
      return
      end

C One time initialization
      subroutine nucinit(anuc)
      common /nucdis/ probmi(31),probdi(31),rhotab(30),pnorm
      dimension probm(30),probd(30)


C   DOVER ET AL DISTRIBUTION
      DATA PROBD/
     1 0.0,0.2,4.7,20.2,67.6,186.2,384.2,802.4,1409.0,2285.8
     1 ,3264.0,3913.1,4665.6,4901.0,4170.9,4050.0,2973.6,1976.8,1296.0,985.5
     1 ,704.0,370.4,232.3,177.7,80.6,35.7,20.0,10.0,5.0,2.5/


C   SET UP MATTER DISTRIBUTION IN NUCLEUS
C   UNITS ARE IN FERMIS
      PROBS=0.0
      RBIN=0.2
      RMIN=0.0
      RMAX=6.0
      NBIN=(RMAX-RMIN)/RBIN
      call nucrodis(rmin,rbin,nbin,probm,rhotab)
      CALL CFD(PROBMI,NBIN,PROBM)
      CALL CFD(PROBDI,NBIN,PROBD)

      DO I1=1,NBIN
       PROBS=PROBS+PROBM(I1)
      END DO
      PNORM=ANUC/PROBS
      DO i1=1,NBIN
       RHOTAB(i1)=PNORM*RHOTAB(i1)
      END DO

      return
      end


      subroutine nuckin(vin,pin,am,pm,pout,pn,eout)
      dimension vin(3),pout(3)
      common/nuccms/rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
      logical prot
      dimension p1in(3),p2in(3),pdum(3)
      parameter (bm2=939.**2)
      do i1=1,3
         p1in(i1)=pin*vin(i1)
         p2in(i1)=-p1in(i1)
      end do
      e1=sqrt(pin**2+am**2)
      e2=sqrt(pin**2+bm2)
      call lloren(pout,eout,p1in,e1,beta)
      call lloren(pdum,edum,p2in,e2,beta)
      pn=sqrt(eout**2-am**2)
      pm=sqrt(edum**2-bm2)
      return
      end


      subroutine nucnuc(idec)
      common/nucint/uin,chint,ein,pin(3),upin(3)
      common/nuccms/rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
#include "nucnewprd.h"
      common/nuchad/ihap,nel(2)
      logical prot

      parameter (step=0.2, fac=1e-2, amn=939., amn2=amn**2)

      dimension xspp(158),xspn(158),xdpp(130),xdpn(130)
      data xspp/
     $     13.,  26.,  50.,  90., 126., 175., 235., 300., 380., 480.,
     $     585., 710., 860.,1035.,1242.,1436.,1580.,1695.,1785.,1866.,
     $     1937.,2000.,2050.,2090.,2128.,2160.,2188.,2208.,2228.,2240.,
     $     2248.,2252.,2253.,2253.,2253.,2252.,2250.,2248.,2245.,2240.,
     $     2238.,2234.,2230.,2225.,2220.,2215.,2210.,2205.,2200.,2195.,
     $     2190.,2185.,2178.,2172.,2168.,2160.,2154.,2150.,2140.,2135.,
     $     2128.,2120.,2112.,2102.,2096.,2089.,2079.,2070.,2061.,2054.,
     $     2042.,2036.,2028.,2016.,2005.,1995.,1985.,1975.,1968.,1959.,
     $     1949.,1939.,1928.,1919.,1908.,1898.,1888.,1878.,1868.,1858.,
     $     1848.,1838.,1828.,1817.,1807.,1797.,1786.,1777.,1765.,1758.,
     $     1748.,1737.,1728.,1716.,1705.,1696.,1685.,1675.,1665.,1655.,
     $     1645.,1635.,1627.,1615.,1605.,1595.,1584.,1575.,1563.,1555.,
     $     1542.,1532.,1520.,1512.,1500.,1492.,1480.,1470.,1460.,1450.,
     $     1440.,1430.,1419.,1408.,1398.,1387.,1377.,1367.,1356.,1344.,
     $     1335.,1324.,1315.,1305.,1295.,1285.,1275.,1265.,1257.,1245.,
     $     1237.,1228.,1218.,1208.,1198.,1188.,1178.,1169./

      data xspn/
     $     18.,  28.,  48.,  75., 117., 159., 208., 270., 340., 440.,
     $     560., 720., 910.,1060.,1170.,1260.,1325.,1375.,1390.,1396.,
     $     1397.,1396.,1393.,1390.,1385.,1382.,1380.,1376.,1372.,1366.,
     $     1360.,1357.,1352.,1346.,1341.,1338.,1332.,1325.,1320.,1316.,
     $     1310.,1304.,1300.,1295.,1288.,1281.,1278.,1270.,1264.,1259.,
     $     1254.,1245.,1240.,1235.,1230.,1223.,1218.,1212.,1203.,1199.,
     $     1193.,1184.,1180.,1175.,1170.,1162.,1158.,1152.,1146.,1140.,
     $     1136.,1130.,1122.,1120.,1112.,1108.,1100.,1096.,1090.,1082.,
     $     1078.,1070.,1065.,1060.,1056.,1050.,1045.,1040.,1034.,1028.,
     $     1021.,1018.,1010.,1005.,1000., 995., 988., 980., 977., 970.,
     $     962., 958., 951., 945., 940., 935., 930., 922., 918., 910.,
     $     904., 900., 895., 889., 882., 878., 872., 865., 860., 855.,
     $     850., 842., 838., 832., 826., 820., 817., 810., 804., 800.,
     $     797., 790., 783., 780., 774., 768., 761., 758., 752., 744.,
     $     740., 735., 730., 722., 720., 715., 710., 705., 698., 694.,
     $     690., 683., 680., 675., 670., 662., 660., 653./

      data xdpp/
     $     0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   8.,
     $     9.,   9.,  11.,  16.,  21.,  38.,  60.,  90., 112., 130.,
     $     146., 160., 177., 190., 200., 215., 228., 240., 252., 265.,
     $     280., 292., 305., 320., 335., 349., 361., 379., 392., 405.,
     $     420., 435., 450., 460., 479., 490., 504., 520., 535., 550.,
     $     565., 580., 595., 610., 623., 640., 655., 670., 680., 698.,
     $     710., 726., 740., 756., 770., 788., 800., 820., 834., 850.,
     $     860., 880., 893., 910., 923., 940., 957., 970., 986.,1002.,
     $     1020.,1035.,1050.,1066.,1080.,1098.,1112.,1128.,1141.,1160.,
     $     1175.,1190.,1208.,1221.,1239.,1255.,1270.,1285.,1300.,1320.,
     $     1338.,1355.,1370.,1388.,1400.,1420.,1436.,1450.,1466.,1480.,
     $     1496.,1510.,1525.,1540.,1558.,1573.,1588.,1600.,1620.,1637.,
     $     1652.,1670.,1682.,1700.,1718.,1731.,1748.,1762.,1780.,1798./

      data xdpn/
     $     0.,  40.,  80., 128., 170., 230., 300., 370., 450., 520.,
     $     560., 600., 630., 660., 690., 715., 740., 760., 785., 805.,
     $     828., 850., 870., 888., 908., 928., 943., 962., 980., 998.,
     $     1017.,1030.,1045.,1063.,1080.,1095.,1110.,1125.,1140.,1156.,
     $     1170.,1185.,1200.,1215.,1230.,1240.,1254.,1266.,1280.,1293.,
     $     1305.,1320.,1333.,1345.,1360.,1373.,1382.,1397.,1410.,1420.,
     $     1437.,1448.,1460.,1473.,1488.,1498.,1510.,1521.,1537.,1548.,
     $     1560.,1572.,1582.,1598.,1610.,1620.,1635.,1645.,1658.,1670.,
     $     1680.,1696.,1708.,1720.,1732.,1743.,1757.,1770.,1780.,1797.,
     $     1808.,1820.,1834.,1844.,1860.,1870.,1880.,1898.,1910.,1920.,
     $     1938.,1950.,1960.,1978.,1990.,2002.,2020.,2032.,2048.,2060.,
     $     2077.,2090.,2104.,2120.,2130.,2150.,2162.,2180.,2195.,2210.,
     $     2225.,2240.,2257.,2270.,2286.,2300.,2315.,2330.,2342.,2360./

      idec=1
      if(ecms2.lt.4070710.)go to 200
      if(ecms2.le.(amn+unucl)**2)go to 200
      pcmsm=sqrt((ecms2-(amn-unucl)**2)*(ecms2-(amn+unucl)**2)/(4*ecms2))
      eineq=(ecms2-unucl**2-amn**2)/(2*unucl)
      emin=amn+360.
      if(eineq.lt.emin)go to 200
      call nuchis(eineq-emin,20.,158,ilo,ihi,erem)
      call nuchis(eineq-amn-920.,20.,130,ilo1,ihi1,erem1)
      if(prot)then
         ka=1
      else
         ka=0
      end if
      ichint=chint
      if(ka-ichint.eq.0)then
         psp=xspp(ihi)*erem+xspp(ilo)*(1-erem)
         pdp=xdpp(ihi1)*erem1+xdpp(ilo1)*(1-erem1)
      else
         psp=xspn(ihi)*erem+xspn(ilo)*(1-erem)
         pdp=xdpn(ihi1)*erem1+xdpn(ilo1)*(1-erem1)
      end if
      if(ecms2.lt.4653520.)pdp=0.
      psp=psp*fac
      pdp=pdp*fac
      ptot=psp+pdp
      if(ptot.eq.0)go to 200
      xr=ranf()+1.e-6
      yr=ranf()*ptot
      ptot=ptot*rhon
      dist=-log(xr)/ptot
      if(step.lt.dist)go to 200
      if(yr.lt.psp)then
         iqq=1                  !single production
      else
         iqq=2                  !double production
      end if
      call nucprd(iqq,eineq,ka,piso,chiso,uiso)

      p1=sqrt(piso(1,1)**2+piso(2,1)**2+piso(3,1)**2)
      if(p1.lt.pfermi)go to 200
      p2=sqrt(piso(1,3)**2+piso(2,3)**2+piso(3,3)**2)
      if(p2.lt.pfermi)go to 200
      idec=-iqq
      ihap=1
      nel(iqq)=1
      if(p1.gt.p2)then
         iqq=1
         pp=p1
      else
         iqq=3
         pp=p2
      end if
      ein=piso(4,iqq)
      chint=chiso(iqq)
      do i1=1,3
         pin(i1)=piso(i1,iqq)
         upin(i1)=pin(i1)/pp
      end do

 200  return
      end

      subroutine nucpion(ides)
      common /nucint/uin,chint,ein,pin(3),upin(3)
      common /nuccms/ rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
#include "nucnewprd.h"
      common /nuchad/ ihap,nel,ncx
#include "hadron.h"
#include "partnucint.h"
      dimension xabs(21),xpeln(126),xpelp(126),xpelm(126),xpchx(126),
     $     xspp(117),xspn(117),xspm(117),
     $     dang(5,21,4)
      dimension angld(18),angldi(19),aa(5),b1(51),b2(51),b3(51)
      dimension pout1l(3),upoutc(3)
      logical prot,notuj

      parameter (upi=139.6, amn=939., amn2=amn**2, step=0.2 )
      parameter (rhod=9./(2.*3.141*9.5**3) ,fac=1e-2)
      parameter (nangl =18)
      integer called /0/
      real junk

      dimension phs31(21),php31(21),php33(21)
      dimension phs11(21),php11(21),php13(21)
      dimension pep11(21),phd13(21),phd15(21)
      dimension angdata(51,428)
      dimension gm1(4),gm3(4)
      data gm1/0.,0.66667,-0.4714,0.33333/
      data gm3/1.,0.33333, 0.4714,0.66667/

      data phs11/
     $     0.,   4.93, 6.51, 7.46, 8.08, 8.50, 8.79, 9.00, 9.17, 9.33,
     $     9.49, 9.69, 9.93,10.23,10.61,11.06,11.62,12.28,13.05,13.93,
     $     14.94/
      data php11/
     $     0.00,-0.68,-1.49,-2.06,-2.29,-2.18,-1.78,-1.13,-0.26, 0.82,
     $     2.12, 3.71, 5.67, 8.18,11.41,15.58,20.87,23.92,28.09,32.26,
     $     36.43/
      data php13/
     $     0.,  -0.15,-0.40,-0.68,-0.97,-1.26,-1.53,-1.79,-2.04,-2.28,
     $     -2.52,-2.78,-3.05,-3.34,-3.69,-4.08,-4.55,-5.10,-5.76,-6.54,
     $     -7.46/
      data phs31/
     $     0.00,-2.85,  -4.28, -5.80, -7.30,- 8.80,-10.23,-11.71,-13.16,
     $     -14.57,-15.93,-17.22,-18.44,-19.56,-20.59,-21.51,-22.30,-23.02,
     $     -23.61,-24.09,-24.46/
      data php31/
     $     0.00,-0.22,-0.60,-1.07,-1.60,-2.17,-2.77,-3.39,-4.01,-4.63,
     $     -5.24,-5.84,-6.41,-6.96,-7.48,-7.96,-8.41,-9.16,-9.71,-10.25,
     $     -10.80/
      data php33/
     $     0.00,  1.29,  3.95,  7.95, 13.62, 21.56, 32.53, 47.17,
     $     66.10, 81.98, 97.51,109.47,118.42,125.19,130.44,134.63,
     $     138.07,140.97,143.45,145.62,147.54/
      data pep11/
     $     1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,
     $     1.   ,0.996,0.986,0.970,0.947,0.908,0.861,0.810,0.762,0.710,
     $     0.660/
      data phd13/
     $     0.  , 0.00, 0.01, 0.04, 0.08, 0.15, 0.25, 0.38, 0.56, 0.78,
     $     1.05, 1.38, 1.77, 2.24, 2.78, 3.41, 4.13, 4.94, 5.86, 6.89,
     $     8.04/
      data phd15/
     $     0.  , 0.00, 0.02, 0.06, 0.12, 0.20, 0.30, 0.42, 0.57, 0.73,
     $     0.91, 1.10, 1.31, 1.52, 1.73, 1.95, 2.16, 2.36, 2.54, 2.71,
     $     2.85/

C     PARTIAL AND       TOTAL CROSS SECTIONS (*100 mb)
c     pi0-p elastic
      data xpeln/
     $     0.,  90., 258., 745.,1575.,2760.,4240.,6190.,7900.,8930.,
     $     7430.,6480.,5540.,4633.,3725.,2985.,2490.,2160.,1908.,1743.,
     $     1565.,1465.,1380.,1300.,1260.,1228.,1200.,1200.,1205.,1235.,
     $     1238.,1140.,1050., 980., 925., 930., 930., 950., 980.,1000.,
     $     1040.,1095.,1170.,1280.,1413.,1568.,1580.,1545.,1513.,1495.,
     $     1445.,1410.,1380.,1350.,1320.,1285.,1280.,1270.,1270.,1280.,
     $     1295.,1320.,1340.,1370.,1390.,1405.,1393.,1360.,1340.,1310.,
     $     1280.,1250.,1220.,1190.,1180.,1163.,1150.,1150.,1150.,1140.,
     $     1125.,1125.,1120.,1120.,1110.,1110.,1100.,1090.,1080.,1050.,
     $     1018., 990., 980., 950., 930., 883., 880., 870., 860., 840.,
     $     835., 820., 810., 800., 790., 800., 780., 775., 770., 765.,
     $     765., 750., 740., 740., 740., 738., 735., 730., 725., 725.,
     $     725., 715., 710., 710., 710., 710./

c     pi-p charge exchange
      data xpchx/
     $     390., 480., 620., 875.,1240.,1675.,2200.,3100.,4350.,4460.,
     $     4000.,3400.,2950.,2625.,2300.,2000.,1650.,1485.,1375.,1300.,
     $     1215.,1150.,1080.,1030., 975., 920., 870., 815., 765., 718.,
     $     675., 630., 595., 555., 515., 490., 470., 451., 450., 468.,
     $     500., 550., 625., 680., 675., 575., 500., 430., 375., 345.,
     $     300., 285., 265., 250., 240., 225., 210., 200., 195., 185.,
     $     178., 175., 170., 165., 160., 155., 152., 150., 149., 148.,
     $     145., 140., 135., 130., 125., 120., 118., 115., 110., 107.,
     $     104., 100.,  99.,  98.,  95.,  90.,  85.,  80.,  78.,  77.,
     $     75.,  73.,  72.,  69.,  67.,  65.,  63.,  61.,  58.,  55.,
     $     52.,  51.,  50.,  50.,  50.,  50.,  50.,  49.,  48.,  47.,
     $     46.,  44.,  42.,  39.,  37.,  35.,  34.,  33.,  32.,  31.,
     $     29.,  27.,  25.,  24.,  24.,  23./

c     pi-p elastic
      data xpelm/
     $     125., 150., 175., 245., 380., 600., 970.,1500.,2140.,2310.,
     $     2295.,2070.,1795.,1550.,1360.,1230.,1130.,1070.,1050.,1070.,
     $     1120.,1175.,1235.,1300.,1400.,1500.,1600.,1700.,1835.,1970.,
     $     2050.,1915.,1770.,1650.,1570.,1520.,1510.,1525.,1550.,1600.,
     $     1685.,1800.,2000.,2230.,2475.,2635.,2510.,2300.,2140.,2000.,
     $     1870.,1750.,1670.,1585.,1505.,1440.,1395.,1340.,1299.,1260.,
     $     1215.,1175.,1140.,1099.,1060.,1040.,1010., 999., 990., 975.,
     $     960., 955., 945., 935., 925., 925., 935., 965., 985.,1000.,
     $     1015.,1030.,1060.,1080.,1095.,1100.,1095.,1090.,1070.,1035.,
     $     1000., 960., 905., 855., 820., 800., 785., 780., 775., 770.,
     $     765., 760., 759., 757., 756., 755., 750., 745., 740., 735.,
     $     730., 728., 725., 725., 724., 723., 722., 722., 722., 721.,
     $     720., 715., 710., 710., 710., 710./

c     pi+p elastic
      data xpelp/
     $     180., 400., 990.,2170.,4000.,6580.,9680.,13920.,18000.,20000.,
     $     16550.,14200.,12250.,10320.,8400.,6725.,5510.,4725.,4130.,3690.,
     $     3230.,2885.,2600.,2300.,2090.,1875.,1675.,1500.,1340.,1200.,
     $     1100., 998., 920., 860., 820., 810., 810., 825., 850., 875.,
     $     900., 940., 975.,1000.,1030.,1075.,1130.,1200.,1275.,1330.,
     $     1350.,1335.,1330.,1330.,1345.,1355.,1380.,1400.,1460.,1500.,
     $     1555.,1625.,1700.,1800.,1875.,1920.,1925.,1890.,1830.,1790.,
     $     1725.,1690.,1640.,1600.,1550.,1505.,1475.,1430.,1400.,1365.,
     $     1335.,1300.,1280.,1250.,1225.,1205.,1195.,1175.,1150.,1135.,
     $     1105.,1095.,1080.,1060.,1030.,1020.,1005., 990., 980., 970.,
     $     960., 940., 930., 920., 900., 898., 890., 880., 870., 860.,
     $     840., 830., 820., 810., 800., 795., 785., 780., 775., 765.,
     $     760., 755., 750., 745., 740., 735./

c     pi+p single pi production
      data xspp/
     $     0.,   0.,   1.,   3.,   5.,   9.,  13.,  21.,  30.,  41.,
     $     55.,  75.,  98., 122., 153., 183., 218., 255., 293., 335.,
     $     383., 445., 530., 650., 780., 880., 990.,1082.,1170.,1251.,
     $     1321.,1380.,1425.,1460.,1490.,1512.,1520.,1500.,1445.,1395.,
     $     1330.,1330.,1388.,1450.,1521.,1595.,1690.,1755.,1842.,1930.,
     $     2012.,2090.,2145.,2180.,2180.,2160.,2155.,2170.,2195.,2220.,
     $     2222.,2220.,2205.,2178.,2140.,2090.,2030.,1990.,1955.,1925.,
     $     1900.,1890.,1862.,1850.,1838.,1830.,1820.,1822.,1823.,1828.,
     $     1833.,1840.,1850.,1859.,1870.,1883.,1900.,1910.,1928.,1943.,
     $     1960.,1980.,2000.,2018.,2038.,2058.,2080.,2098.,2119.,2139.,
     $     2158.,2180.,2198.,2215.,2235.,2252.,2270.,2278.,2286.,2295.,
     $     2300.,2305.,2313.,2319.,2321.,2325.,2330./

c     pi0p single pi production
      data xspn/
     $     0.,   0.,   1.,   8.,  23.,  42.,  69., 100., 135., 175.,
     $     221., 270., 320., 372., 430., 492., 560., 640., 755.,1085.,
     $     1210.,1252.,1350.,1405.,1395.,1365.,1330.,1365.,1455.,1610.,
     $     1755.,1928.,2105.,2240.,2260.,2219.,2170.,2115.,2045.,1965.,
     $     1879.,1825.,1775.,1750.,1752.,1790.,1835.,1885.,1938.,1995.,
     $     2055.,2125.,2172.,2210.,2230.,2235.,2248.,2265.,2288.,2302.,
     $     2310.,2317.,2312.,2300.,2280.,2259.,2230.,2199.,2170.,2145.,
     $     2121.,2100.,2081.,2067.,2055.,2045.,2041.,2041.,2048.,2062.,
     $     2085.,2110.,2140.,2170.,2200.,2225.,2240.,2252.,2265.,2278.,
     $     2290.,2305.,2319.,2332.,2348.,2360.,2372.,2388.,2400.,2412.,
     $     2425.,2438.,2448.,2458.,2468.,2475.,2485.,2488.,2490.,2489.,
     $     2488.,2485.,2482.,2479.,2475.,2470.,2465./

c     pi-p single pi production
      data xspm/
     $     0.,   5.,  10.,  20.,  39.,  68., 110., 160., 230., 305.,
     $     382., 470., 548., 630., 712., 805., 905.,1021.,1230.,1820.,
     $     2035.,2060.,2170.,2160.,2040.,1840.,1670.,1635.,1740.,1965.,
     $     2200.,2435.,2740.,3020.,3030.,2920.,2815.,2730.,2639.,2530.,
     $     2420.,2315.,2190.,2060.,1992.,1979.,1985.,2008.,2045.,2085.,
     $     2128.,2160.,2205.,2241.,2275.,2308.,2340.,2358.,2378.,2392.,
     $     2400.,2410.,2420.,2422.,2425.,2425.,2425.,2410.,2380.,2358.,
     $     2335.,2310.,2295.,2282.,2272.,2265.,2260.,2265.,2279.,2302.,
     $     2335.,2375.,2425.,2471.,2512.,2545.,2580.,2590.,2602.,2612.,
     $     2622.,2630.,2643.,2650.,2659.,2665.,2670.,2679.,2685.,2692.,
     $     2698.,2700.,2703.,2705.,2702.,2701.,2700.,2692.,2688.,2679.,
     $     2671.,2660.,2651.,2640.,2628.,2613.,2600./

C     PI DEUTERON       ABSORPTION (in mb (without factor 100))
      DATA XABS/
     $     2., 3., 5.,12.,18.,27.,37.,50.,60.,45.,
     $     30.,20., 10.,5.,3.,2.,2.,1.,1.,1.,0./


C *** one time initialization
      if(called.eq.0)then
         called=1
         do  i1=1,126
            xpeln(i1)=xpeln(i1)*fac
            xpelp(i1)=xpelp(i1)*fac
            xpelm(i1)=xpelm(i1)*fac
            xpchx(i1)=xpchx(i1)*fac
            if(i1.le.117)then
               xspp(i1)=xspp(i1)*fac
               xspn(i1)=xspn(i1)*fac
               xspm(i1)=xspm(i1)*fac
            end if
         end do
         do i1=1,20
            tl=i1*20.
            el=tl+upi
            pl=sqrt(el**2-upi**2)
            pc=pl/sqrt(1.+2*el/amn+(upi/amn)**2)/197.329
            pk=pc**2
            call nuchis(tl,20.,21,ilo,ihi,erem)

            s11=phs11(ihi)*erem+phs11(ilo)*(1.-erem)
            s31=phs31(ihi)*erem+phs31(ilo)*(1.-erem)

            p11=php11(ihi)*erem+php11(ilo)*(1.-erem)
            p13=php13(ihi)*erem+php13(ilo)*(1.-erem)
            p31=php31(ihi)*erem+php31(ilo)*(1.-erem)
            p33=php33(ihi)*erem+php33(ilo)*(1.-erem)

            ep11=pep11(ihi)*erem+pep11(ilo)*(1.-erem)
            d13=phd13(ihi)*erem+phd13(ilo)*(1.-erem)
            d15=phd15(ihi)*erem+phd15(ilo)*(1.-erem)

            w01=      (sind(s31)**2+sind(p31-p33)**2) !*g3**2

            w02= (sind(s11)**2-sind(s31-s11)**2-3/2.*sind(s31)**2
     $           -ep11*(sind(p31-p11)**2-sind(p33-p11)**2)
     $           +sind(p31-p13)**2-sind(p33-p13)**2
     $           +sind(s31-d13)**2+3/2.*sind(s31-d15)**2
     $           -sind(d13)**2-3/2.*sind(d15)**2) !*g1*g3

            w03=(-3/2.*sind(s11)**2+ep11*sind(p11-p13)**2+(ep11-1)**2/4.
     $           +sind(s11-d13)**2+3/2.*sind(s11-d15)**2
     $           -3/2.*sind(d13-d15)**2+3/2.*sind(d13)**2
     $           +9/4.*sind(d15)**2) !*g1**2

            w11=2*(  sind(s31)*sind(p31)*cosd(s31-p31)
     $           +2*sind(s31)*sind(p33)*cosd(s31-p33)) !*g3**2

            w12=(ep11*(sind(p11)**2-sind(s31-p11)**2)
     $           +sind(p31)**2-sind(s11-p31)**2
     $           +2*(sind(p33)**2+sind(p13)**2-sind(s31-p13)**2
     $           -sind(s11-p33)**2)
     $           +3*(sind(s11)**2+sind(s31)**2))
     $           +(-2*sind(p31-d13)**2+9/2.*sind(p31-d15)**2-5/2.*sind(p31)**2
     $           +5*sind(p33-d13)**2-5*sind(p33)**2-3*sind(d13)**2
     $           -9/2.*sind(d15)**2) !*g1*g3

            w13=(ep11*(sind(p11)**2-sind(s11-p11)**2)
     $           +2*(sind(p13)**2-sind(s11-p13)**2)+3*sind(s11)**2
     $           +ep11*(-2*sind(p11-d13)**2+9/2.*sind(p11-d15)**2
     $           -5/2.*sind(p11)**2)
     $           +5*sind(p13-d13)**2-5*sind(p13)**2-3*sind(d13)**2
     $           -9/2.*sind(d15)**2) !*g1**2

            w21=3*(sind(p33)**2+2*sind(p31)*sind(p33)*cosd(p31-p33)) !*g3**2

            w22=     +3*(ep11*(sind(p11)**2-sind(p33-p11)**2)
     $           +    2*(sind(p13)**2+sind(p33)**2)+sind(p31)**2
     $           -      (sind(p33-p13)**2+sind(p31-p13)**2)
     $           -sind(s31-d13)**2-3/2.*sind(s31-d15)**2+5/2.*sind(s31)**2
     $           +sind(d13)**2+3/2.*sind(d15)**2) !*g1*g3

            w23=   +  3*(2*sind(p13)**2+ep11*(sind(p11)**2-sind(p11-p13)**2)
     $           -sind(s11-d13)**2-3/2.*sind(s11-d15)**2+5/2.*sind(s11)**2
     $           +6*sind(d13-d15)**2-4*sind(d13)**2-6*sind(d15)**2) !*g1**2

            w32=(-15/2.*sind(p31-d15)**2+15/2.*sind(p31)**2-9*sind(p33-d13)**2
     $           -6*sind(p33-d15)**2+15*sind(p33)**2+9*sind(d13)**2
     $           +27/2.*sind(d15)**2) !*g1*g3

            w33= +(15/2.*ep11*(-sind(p11-d15)**2+sind(p11)**2)
     $           -9*sind(p13-d13)**2-6*sind(p13-d15)**2+15*sind(p13)**2
     $           +9*sind(d13)**2+27/2.*sind(d15)**2) !*g1**2

            w43= 45/2.*(sind(d13)**2-sind(d13-d15)**2)+135/4.*sind(d15)**2 !*g1**2
            do i2=1,4
               g1=gm1(i2)
               g3=gm3(i2)
               i3=i1+1
               dang(1,i3,i2)=(w01*g3**2+w02*g1*g3+w03*g1**2)/pk
               dang(2,i3,i2)=(w11*g3**2+w12*g1*g3+w13*g1**2)/pk
               dang(3,i3,i2)=(w21*g3**2+w22*g1*g3+w23*g1**2)/pk
               dang(4,i3,i2)=(          w32*g1*g3+w33*g1**2)/pk
               dang(5,i3,i2)=(                    w43*g1**2)/pk
            end do
         end do

#if defined(VMS) &&    !defined(USE_TEXT_FILE)
         OPEN(UNIT=33, FILE='/u/sbd0v3/fgoebel/nucang.dat', STATUS='OLD',
     $        readonly,shared)

         do i1=1,428
            read(33,*)(angdata(i2,i1),i2=1,51)
         end do

c     changed by Florian Goebel on 11/15/95

c     open(33,status='old',form='unformatted',access='direct',
c     $    organization='relative',recl=51,readonly,shared,
c     $    file='sys$usrdevice:[gajewski.dat]nucang.dat')
c     do i1=1,428
c     read(33,rec=i1)(angdata(i2,i1),i2=1,51)
c     end do
#else
C     open(33,status='old',form='unformatted',
C     $    recl=51,readonly,shared,convert='vaxd',
C     $    file='~/dat/nucang.dat')
c     read(33)angdata
#endif
         close(33)

      end if
C***** End one time initialization

      if(ides.eq.10)then
         notuj=.true.
      else
         notuj=.false.
      end if
      ides=1

C     MAGNITUDE OF CMS MOMENTUM
      IF(ECMS2.LE.(UPI+UNUCL)**2) then
         pninttype=0
         GO TO 200
      endif
      PCMSM=SQRT((ECMS2-(UPI-UNUCL)**2)*(ECMS2-(UPI+UNUCL)**2)/(4.0*ECMS2))

C     EQUIVALENT PION   ENERGY FOR STATIONARY TARGET
      EINEQ=(ECMS2-UNUCL**2-UPI**2)/(2.0*UNUCL)
      IF(EINEQ.LT.UPI) then
         pninttype=0
         GO TO 200
      endif
      t=eineq-upi

C     ABSORPTION CROSS SECTION WEIGHTED BY RATIO OF LOCAL
C     NUCLEON DENSITY TO DEUTERON DENSITY.AT THIS STAGE RHON
C     CONTAINS 0.1 FACTOR FOR UNITS.THE RHOD FACTOR ADJUSTMENT
C     THROUGH PUTTING THE "DEUTERON RADIUS" AT 9.5 TAKES THIS INTO
C     ACCOUNT.PUTTING IN THE 0.1 FACTOR INTO RHOD IMPLIES THAT 9.5
C     CORRESPONDS TO A SO CALLED "QUASI DEUTERON" RADIUS OF 4.4 FM.

      call nuchis(t,20.,21,ilo,ihi,erem)
c      PPABS=(XABS(ihi)*erem+xabs(ilo)*(1.-erem))*0.5*RHON/RHOD
      PPABS=(XABS(ihi)*erem+xabs(ilo)*(1.-erem))

      call nuchis(t,20.,126,ilo,ihi,erem)
      call nuchis(t-180.,20.,117,ilo1,ihi1,erem1)

      ichint=chint
      if(ichint)1,2,3
 1    continue                  !pi-
      if(prot)then
         ppel =(xpelm(ihi)*erem+xpelm(ilo)*(1-erem))
         ppchx=(xpchx(ihi)*erem+xpchx(ilo)*(1-erem))
         psp  =(xspm(ihi1)*erem1+xspm(ilo1)*(1-erem1))
         ka=2
      else
         ppel =(xpelp(ihi)*erem+xpelp(ilo)*(1-erem))
         psp  =(xspp(ihi1)*erem1+xspp(ilo1)*(1-erem1))
         ppchx=0.
         ka=1
      end if
      go to 4
 2    continue                  !pi0
      ppel =(xpeln(ihi)*erem+xpeln(ilo)*(1-erem))
      ppchx=(xpchx(ihi)*erem+xpchx(ilo)*(1-erem))
      psp  =(xspn(ihi1)*erem1+xspn(ilo1)*(1-erem1))
      ppabs=0.7*ppabs
      ka=4
      go to 4
 3    continue                  !pi+
      if(prot)then
         ppel =(xpelp(ihi)*erem+xpelp(ilo)*(1-erem))
         ppchx=0.
         psp  =(xspp(ihi1)*erem1+xspp(ilo1)*(1-erem1))
         ka=1
      else
         ppel =(xpelm(ihi)*erem+xpelm(ilo)*(1-erem))
         ppchx=(xpchx(ihi)*erem+xpchx(ilo)*(1-erem))
         psp  =(xspm(ihi1)*erem1+xspm(ilo1)*(1-erem1))
         ka=2
      end if
 4    continue
      if(ecms2.lt.1488400.)psp=0.

      PTOT=PPEL+PPCHX+PPABS+psp
      junk=ptot                         ! just used to print later

C********CHOOSE INTERACTION DISTANCE
      XR=ranf()+1.e-6
      YR=ranf()*PTOT
C********HERE WE NEED TO NORMALIZE
      PTOT=PTOT*RHON
      if(ptot.eq.0.) then
         pninttype=0
         go to 200
      endif
      DIST=-LOG(XR)/ptot

C********IF NOTHING HAPPENS
      IF(STEP.LT.DIST) then
         pnintttype=0
         GO TO 200
      endif


      write(6,'('' RAND '',6F8.3)') yr,junk,ppabs,psp,ppel,ppchx

      call hf1(1001,ppabs/junk,1.0)
      call hf1(1002,psp/junk,1.0)
      call hf1(1003,ppel/junk,1.0)
      call hf1(1004,ppchx/junk,1.0)
      call hfill(1011,pcmsm,ppabs/junk,1.0)
      call hfill(1012,pcmsm,psp/junk,1.0)
      call hfill(1013,pcmsm,ppel/junk,1.0)
      call hfill(1014,pcmsm,ppchx/junk,1.0)
      call hfill(1021,pcmsm,ppabs,1.0)
      call hfill(1022,pcmsm,psp,1.0)
      call hfill(1023,pcmsm,ppel,1.0)
      call hfill(1024,pcmsm,ppchx,1.0)

C********IS IT ABSORBED BY DI NUCLEON?
      PSUM=PPABS
      IF(YR.Le.PSUM) THEN
         pninttype=1
         if(notuj)IABS=IABS+1
         Uin=1.e-2
         CHINT=0.0
         ihap=1
         ides=0
         go to 200
      END IF
c*********does pion produce an additional pion ?
      psum=psum+psp
      if(yr.lt.psum) then
         pninttype=2
         go to 300
      endif

C*********DOES  PION SCATTER OFF  FERMI NUCLEON?
      PSUM=PSUM+PPEL
      IF(YR.LT.PSUM) THEN
         pninttype=3
         IPATH=0
      else
C*********PION CHARGE EXCHANGE ON FERMI NUCLEON
         pninttype=4
         IPATH=1
         ka=3
      END IF


C*********COME HERE FOR PI-P KINEMATICS WITH FERMI MOTION
 104  CONTINUE

      if(t.le.400.)then
         do k=1,5
            aa(k)=dang(k,ihi,ka)*erem+dang(k,ilo,ka)*(1-erem)
         end do

         DO K=1,18
            angb=(k-1)*10.+5.
            cab=cosd(angb)
            ang=0
            do k1=5,2,-1
               ang=(ang+aa(k1))*cab
            end do
            ANGLD(K)=sind(ANGB)*(aa(1)+ang)
         END DO
         CALL CFD(ANGLDI,NANGL,ANGLD)

C     DO SCATTERING IN CMS
         xXR=ranf()
         ANG=PRBIN(ANGLDI,NANGL,0.,10.,xXR)
         CANGL=cosd(ANG)
      else
         call nuchis(t-400.,20.,107,ilo2,ihi2,erem2)
         krec1=(ka-1)*107+ilo2
         krec2=krec1-ilo2+ihi2
         if(krec1.ne.krecsav)then
            krecsav=krec1

            do k=1,51
               b1(k)=angdata(k,krec1)
               b2(k)=angdata(k,krec2)
            end do

c     read(33,rec=krec1)b1
c     read(33,rec=krec2)b2

            do k=1,51
               b3(k)=b1(k)*erem2+b2(k)*(1-erem2)
            end do
         end if
         xxr=ranf()
         call nuchis(xxr,0.02,51,ill,ihh,errm)
         cangl=b3(ill)*errm+b3(ihh)*(1-errm)
      end if
      PHI=-1
      CALL CONE(UP4CMS,UPOUTC,CANGL,PHI)
      call nuckin(upoutc,pcmsm,upi,pout2lm,pout1l,pout1lm,eout1l)

      IF(POUT2LM.LT.PFERMI) then
         pninttype=0
         GO TO 200
      endif

      DO i2=1,3
         PIN(i2)=POUT1L(i2)
         UPIN(i2)=PIN(i2)/POUT1LM
      END DO
      EIN=EOUT1L

      IF(IPATH.ne.0)THEN
         if(chint.eq.0.)then
            IF(ranf().LT.0.5) THEN
               CHINT=1.0
            ELSE
               CHINT=-1.0
            END IF
         else
            CHINT=0.0
         end if
         ncx=1
      else
         nel=1
      END IF
      ihap=1

 200  return

 300  continue                  !pion production
      call nucpiprd(eineq,ka,piso,chiso,uiso)
      pout2lm=sqrt(piso(4,1)**2-amn2)
      IF(POUT2LM.LT.PFERMI) then
         pninttype=0
         GO TO 200
      endif
      pout1lm=sqrt(piso(1,3)**2+piso(2,3)**2+piso(3,3)**2)
      DO i2=1,3
         PIN(i2)=piso(i2,3)
         UPIN(i2)=PIN(i2)/POUT1LM
      END DO
      EIN=piso(4,3)
      chint=chiso(3)
      ides=-1
      ihap=1
      if(notuj)inel=inel+1
      do i1=1,3
         write(6,'('' % products '',i3,f5.1,5f8.1)')i1,chiso(i1),
     $        piso(1,i1),piso(2,i1),piso(3,i1),piso(4,i1),uiso(i1)
      end do
      go to 200
      end

      subroutine nucpiprd(e,ka,piso,chiso,uiso)
      dimension piso(4,3),chiso(3),uiso(3)
c
c       ka=      1 pi+ p, pi- n  ??? I FUCKED THIS DOCUMENTATION UP - CDM
c                2 pi+ n, pi- p
c                3 pi0 p, pi0 n
c
      common /nuccms/rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
      logical prot

      dimension fripn(117),fmxsp(117),pnmi(101)
      dimension pcfsl1(117),pcfsl2(117),pnfsl1(117),pnfsl2(117)

      dimension dldir(3),pdl(4),pdldec(4,3),rc(5)
      parameter (upi=139.6, amnucl=939.,upo=134.9)
      data fripn/
     $    3.,  14.,  39.,  90., 181., 326., 540., 838.,1236.,1726.,
     $ 2259.,2732.,3135.,3484.,3784.,4036.,4242.,4410.,4549.,4669.,
     $ 4775.,4869.,4952.,5025.,5092.,5152.,5207.,5258.,5305.,5348.,
     $ 5388.,5426.,5461.,5496.,5530.,5564.,5599.,5636.,5674.,5714.,
     $ 5755.,5798.,5843.,5889.,5937.,5986.,6036.,6088.,6140.,6193.,
     $ 6246.,6299.,6353.,6406.,6459.,6513.,6566.,6621.,6676.,6733.,
     $ 6791.,6851.,6912.,6975.,7040.,7107.,7175.,7246.,7319.,7393.,
     $ 7467.,7542.,7617.,7691.,7765.,7838.,7909.,7979.,8048.,8114.,
     $ 8179.,8242.,8304.,8365.,8424.,8482.,8538.,8594.,8648.,8702.,
     $ 8754.,8806.,8857.,8907.,8957.,9007.,9057.,9106.,9154.,9202.,
     $ 9250.,9298.,9346.,9393.,9440.,9487.,9534.,9581.,9628.,9674.,
     $ 9721.,9767.,9814.,9860.,9907.,9954.,10000./

      data fmxsp/
     $   53., 106., 150., 192., 233., 275., 318., 360., 403., 447.,
     $  491., 536., 581., 627., 672., 718., 765., 811., 858., 905.,
     $  953.,1000.,1048.,1096.,1144.,1192.,1240.,1288.,1336.,1385.,
     $ 1433.,1482.,1530.,1579.,1628.,1677.,1725.,1774.,1823.,1872.,
     $ 1921.,1970.,2019.,2068.,2116.,2165.,2214.,2263.,2312.,2361.,
     $ 2410.,2459.,2508.,2557.,2606.,2655.,2704.,2753.,2802.,2850.,
     $ 2899.,2948.,2997.,3046.,3095.,3144.,3192.,3241.,3290.,3339.,
     $ 3387.,3436.,3485.,3534.,3582.,3631.,3680.,3728.,3777.,3826.,
     $ 3874.,3923.,3971.,4020.,4069.,4117.,4166.,4214.,4263.,4311.,
     $ 4360.,4408.,4457.,4505.,4553.,4602.,4650.,4699.,4747.,4795.,
     $ 4844.,4892.,4940.,4989.,5037.,5085.,5134.,5182.,5230.,5278.,
     $ 5327.,5375.,5423.,5471.,5519.,5568.,5616./

      data pnmi/
     $ 1078.6,1134.,1149.,1159.,1166.,1173.,1178.,1183.,1188.,1192.,
     $ 1195.,1199.,1203.,1206.,1208.,1211.,1214.,1217.,1220.,1222.,
     $ 1225.,1227.,1230.,1233.,1236.,1238.,1241.,1244.,1247.,1251.,
     $ 1254.,1257.,1261.,1265.,1269.,1273.,1277.,1282.,1286.,1291.,
     $ 1297.,1303.,1309.,1316.,1324.,1333.,1342.,1353.,1366.,1380.,
     $ 1395.,1413.,1434.,1458.,1486.,1519.,1551.,1580.,1606.,1629.,
     $ 1651.,1671.,1691.,1710.,1729.,1748.,1766.,1784.,1801.,1817.,
     $ 1832.,1846.,1860.,1873.,1886.,1898.,1910.,1923.,1935.,1948.,
     $ 1960.,1974.,1987.,2002.,2016.,2031.,2047.,2062.,2079.,2095.,
     $ 2112.,2129.,2146.,2163.,2180.,2197.,2215.,2232.,2249.,2265.,
     $ 2282./

      data pcfsl1/
     $ 9000.,1800.,  56.,1430.,4860.,6480.,7079.,7263.,7425.,7455.,
     $ 7268.,7104.,6772.,6363.,5685.,4962.,4709.,3776.,3904.,5056.,
     $ 5131.,4892.,4933.,4781.,4417.,3887.,3309.,3140.,3507.,4142.,
     $ 4671.,5104.,5554.,5893.,5927.,5828.,5730.,5651.,5544.,5406.,
     $ 5261.,5106.,4964.,4748.,4665.,4662.,4768.,4886.,5009.,5149.,
     $ 5263.,5406.,5518.,5630.,5714.,5779.,5860.,5911.,5981.,6023.,
     $ 6070.,6124.,6155.,6191.,6204.,6238.,6236.,6245.,6232.,6208.,
     $ 6210.,6203.,6199.,6211.,6218.,6241.,6258.,6297.,6329.,6263.,
     $ 6453.,6524.,6606.,6691.,6744.,6797.,6854.,6897.,6942.,6972.,
     $ 7001.,7043.,7068.,7106.,7135.,7158.,7192.,7247.,7300.,7300.,
     $ 7335.,7358.,7389.,7420.,7444.,7481.,7505.,7522.,7547.,7565.,
     $ 7595.,7612.,7627.,7650.,7665.,7687.,7759./

      data pcfsl2/
     $ 11000.,6200.,4926.,5820.,8100.,9190.,9605.,9721.,9825.,9838.,
     $ 9703.,9576.,9333.,9045.,8574.,8079.,7898.,7268.,7367.,8186.,
     $ 8232.,8045.,8044.,7886.,7560.,7113.,6613.,6431.,6665.,7124.,
     $ 7514.,7836.,8183.,8445.,8459.,8367.,8280.,8217.,8148.,8053.,
     $ 7958.,7830.,7678.,7461.,7346.,7298.,7327.,7386.,7439.,7506.,
     $ 7561.,7636.,7706.,7780.,7851.,7917.,7986.,8020.,8064.,8083.,
     $ 8117.,8162.,8190.,8226.,8251.,8297.,8320.,8337.,8334.,8322.,
     $ 8326.,8317.,8321.,8329.,8335.,8352.,8366.,8394.,8420.,8382.,
     $ 8519.,8575.,8643.,8711.,8753.,8793.,8832.,8863.,8890.,8907.,
     $ 8921.,8942.,8953.,8973.,8986.,8995.,9013.,9039.,9055.,9070.,
     $ 9089.,9098.,9112.,9127.,9135.,9153.,9162.,9168.,9180.,9186.,
     $ 9201.,9208.,9211.,9220.,9225.,9234.,9276./

      data pnfsl1/
     $ 7000.,   0.,   0.,   0.,2146.,4013.,4700.,4934.,5134.,5181.,
     $ 4979.,4820.,4484.,4062.,3366.,2619.,2392.,1426.,1499.,2559.,
     $ 2658.,2496.,2642.,2676.,2588.,2389.,2236.,2300.,2621.,3042.,
     $ 3394.,3685.,3978.,4207.,4253.,4214.,4169.,4118.,4022.,3911.,
     $ 3783.,3702.,3702.,3672.,3724.,3804.,3953.,4069.,4200.,4335.,
     $ 4442.,4564.,4647.,4723.,4763.,4782.,4821.,4855.,4907.,4945.,
     $ 4976.,5004.,5018.,5029.,5022.,5022.,4994.,4989.,4972.,4950.,
     $ 4950.,4949.,4939.,4947.,4951.,4966.,4975.,5002.,5021.,4964.,
     $ 5096.,5140.,5188.,5241.,5274.,5310.,5352.,5385.,5421.,5447.,
     $ 5474.,5511.,5533.,5569.,5595.,5617.,5647.,5696.,5749.,5738.,
     $ 5768.,5789.,5816.,5842.,5865.,5896.,5917.,5933.,5953.,5970.,
     $ 5994.,6008.,6023.,6043.,6055.,6074.,6127./

      data pnfsl2/
     $ 13000.,   0.,   0.,   0.,4996.,7858.,8937.,9253.,9533.,9571.,
     $ 9221.,8905.,8307.,7607.,6483.,5310.,4920.,3441.,3606.,5405.,
     $ 5537.,5177.,5284.,5125.,4714.,4141.,3621.,3559.,4019.,4729.,
     $ 5340.,5857.,6414.,6855.,6902.,6773.,6650.,6549.,6411.,6235.,
     $ 6046.,5864.,5730.,5530.,5487.,5521.,5669.,5810.,5957.,6117.,
     $ 6243.,6396.,6512.,6625.,6708.,6772.,6855.,6907.,6980.,7023.,
     $ 7073.,7130.,7162.,7201.,7217.,7256.,7257.,7268.,7254.,7229.,
     $ 7232.,7223.,7220.,7232.,7241.,7267.,7285.,7329.,7366.,7293.,
     $ 7508.,7591.,7689.,7790.,7854.,7917.,7983.,8034.,8084.,8117.,
     $ 8147.,8190.,8210.,8253.,8282.,8303.,8338.,8392.,8436.,8447.,
     $ 8482.,8502.,8532.,8560.,8581.,8615.,8635.,8649.,8671.,8686.,
     $ 8713.,8727.,8737.,8756.,8767.,8785.,8855./

      t=e-upi
      w=sqrt(ecms2)
      call nuchis(t-180.,20.,117,ilo,ihi,erem)
      r1  =(fripn(ihi)*erem+fripn(ilo)*(1-erem))*1e-4
      phsm=(fmxsp(ihi)*erem+fmxsp(ilo)*(1-erem))*100.
 1    r=r1*ranf()
      call nuchis(r*100.,1.,101,ilo1,ihi1,erem1)
      amd=pnmi(ihi1)*erem1+pnmi(ilo1)*(1-erem1)
      ed=(amd**2+ecms2-upi**2)/(2*w)
      pd2=ed**2-amd**2
      if(pd2.lt.0.)go to 1
      pdlt=sqrt(pd2)
      phs=(w-ed)*ed*pdlt/w
      if(ranf().gt.phs/phsm)go to 1

      if(ka.eq.1)then           !pi+p,pi-n
         if(ranf().gt.0.75)then !backward 25%
            cdl=-0.9999995
         else                   !isotropic 75%
            cdl=2*ranf()-1
         end if
      else                      !pi+n,pi-p,pi0p,pi0n
         if(ranf().gt.0.80)then
            cdl=-0.9999995      !e>500 -backward 20%
            if(t.lt.500.)cdl=-cdl !e<500 -forward  20%
         else
            cdl=2*ranf()-1      !isotropic 80%
         end if
      end if
      call cone(up4cms,dldir,cdl,-1.)
      do i1=1,3
         pdl(i1)=pdlt*dldir(i1) !delta momentum
         pdldec(i1,3)=-pdl(i1)  !pion momentum (leading pi)
      end do
      pdl(4)=ed
      pdldec(4,3)=w-ed
      call twobiso(amd,amnucl,upi,pdl,pdldec) !delta decay isotropic 1=nuc
      do i1=1,3                 !go back to lab
         call lloren(piso(1,i1),piso(4,i1),pdldec(1,i1),pdldec(4,i1),beta)
      end do
      uiso(1) = amnucl
      uiso(2) = upi
      uiso(3) = upi

c     ********* distribute charges ******
      r=ranf()
      if(ka.eq.1)then
         rc(1)=0.1333
         rc(2)=0.4
         rc(3)=1.
         do k=1,3
            if(r.le.rc(k))go to 101
         end do
         k=3
 101     continue
         if(prot)then           !pi+p
            go to (102,103,104)k
 102        chiso(1)=0          !pi+pi+n
            chiso(2)=1
            chiso(3)=1
            go to 105
 103        chiso(1)=1          !pi+pi0p
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=1
            go to 105
 104        chiso(1)=1          !pi0pi+p
            chiso(2)=1
            chiso(3)=0
            uiso(3) = upo
 105        continue
         else                   !pi-n
            go to (106,107,108)k
 106        chiso(1)=1          !pi-pi-p
            chiso(2)=-1
            chiso(3)=-1
            go to 109
 107        chiso(1)=0          !pi-pi0n
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=-1
            go to 109
 108        chiso(1)=0          !pi0pi-n
            chiso(2)=-1
            chiso(3)=0
            uiso(3) = upo
 109        continue
         end if
      else if(ka.eq.2)then
         v1=(pcfsl1(ihi)*erem+pcfsl1(ilo)*(1-erem))*1e-4
         v2=(pcfsl2(ihi)*erem+pcfsl2(ilo)*(1-erem))*1e-4
         v3=(1-v1)*v2/3
         v4=(1-v1)/3-v3
         rc(1)=v1
         rc(2)=rc(1)+v4
         rc(3)=rc(2)+v3
         rc(4)=rc(3)+v4*2
         rc(5)=1
         do k=1,5
            if(r.le.rc(k))go to 201
         end do
         k=5
 201     continue
         if(prot)then           !pi-p
            go to (202,203,204,205,206)k
 202        chiso(1)=0          !pi+pi-n
            chiso(2)=-1
            chiso(3)=1
            go to 208
 203        chiso(1)=0          !pi-pi+n
            chiso(2)=1
            chiso(3)=-1
            go to 208
 204        chiso(1)=1          !pi0pi-p
            chiso(2)=-1
            chiso(3)=0
            uiso(3) = upo
            go to 208
 205        chiso(1)=1          !pi-pi0p
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=-1
            go to 208
 206        chiso(1)=0          !pi0pi0n
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=0
            uiso(3) = upo
 208        continue
         else                   !pi+n
            go to (212,213,214,215,216)k
 212        chiso(1)=1          !pi-pi+p
            chiso(2)=1
            chiso(3)=-1
            go to 218
 213        chiso(1)=1          !pi+pi-p
            chiso(2)=-1
            chiso(3)=1
            go to 218
 214        chiso(1)=0          !pi0pi+n
            chiso(2)=1
            chiso(3)=0
            uiso(3) = upo
            go to 218
 215        chiso(1)=0          !pi+pi0n
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=1
            go to 218
 216        chiso(1)=1          !pi0pi0p
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=0
            uiso(3) = upo
 218        continue
         end if
      else
         v1=(pnfsl1(ihi)*erem+pnfsl1(ilo)*(1-erem))*1e-4
         v2=(pnfsl2(ihi)*erem+pnfsl2(ilo)*(1-erem))*1e-4
         v3=(1-v1)*v2/3
         v4=(1-v1)-3*v3
         rc(1)=v1/3
         rc(2)=rc(1)+v4
         rc(3)=rc(2)+2*rc(1)
         rc(4)=rc(3)+v3
         rc(5)=1
         do k=1,5
            if(r.lt.rc(k))go to 301
         end do
         k=5
 301     continue
         if(prot)then           !pi0p
            go to (302,303,304,305,306)k
 302        chiso(1)=1          !pi+pi-p
            chiso(2)=-1
            chiso(3)=1
            go to 308
 303        chiso(1)=1          !pi-pi+p
            chiso(2)=1
            chiso(3)=-1
            go to 308
 304        chiso(1)=0          !pi+pi0n
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=1
            go to 308
 305        chiso(1)=0          !pi0pi+n
            chiso(2)=1
            chiso(3)=0
            uiso(3) = upo
            go to 308
 306        chiso(1)=1          !pi0pi0p
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=0
            uiso(3) = upo
 308        continue
         else                   !pi0n
            go to (312,313,314,315,316)k
 312        chiso(1)=0          !pi-pi+n
            chiso(2)=1
            chiso(3)=-1
            go to 318
 313        chiso(1)=0          !pi+pi-n
            chiso(2)=-1
            chiso(3)=1
            go to 318
 314        chiso(1)=1          !pi-pi0p
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=-1
            go to 318
 315        chiso(1)=1          !pi0pi-p
            chiso(2)=-1
            chiso(3)=0
            uiso(3) = upo
            go to 318
 316        chiso(1)=0          !pi0pi0n
            chiso(2)=0
            uiso(2) = upo
            chiso(3)=0
            uiso(3) = upo
 318        continue
         end if
      end if
      return
      end

      subroutine nucprd(no,e,ka,piso,chiso,uiso)
      dimension piso(4,4),chiso(4),uiso(4)
c     ka=       0       primary neutron
c               1       primary proton
      common /nuccms/rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
      logical prot

      dimension dldir(3),pdl(4),pdldec(4,4)

      parameter (upi=139.6, amnucl=939., am2n2=4*amnucl**2)

      dimension frinn(161),fmxsm(161),dmin(101),fmxdn(130),fsln(176)
      data frinn/
     $     2.,   7.,  16.,  30.,  52.,  86., 134., 201., 290., 406.,
     $     551., 727., 945.,1204.,1502.,1832.,2179.,2510.,2806.,3075.,
     $     3321.,3546.,3751.,3935.,4100.,4246.,4374.,4486.,4585.,4674.,
     $     4755.,4830.,4898.,4962.,5022.,5076.,5127.,5174.,5219.,5260.,
     $     5300.,5337.,5372.,5405.,5437.,5467.,5496.,5523.,5549.,5574.,
     $     5599.,5623.,5646.,5670.,5694.,5718.,5742.,5768.,5794.,5820.,
     $     5848.,5877.,5906.,5936.,5967.,5999.,6032.,6065.,6100.,6134.,
     $     6170.,6206.,6242.,6279.,6317.,6354.,6392.,6430.,6469.,6507.,
     $     6545.,6583.,6622.,6660.,6699.,6737.,6777.,6816.,6857.,6898.,
     $     6940.,6983.,7026.,7071.,7116.,7162.,7209.,7258.,7307.,7358.,
     $     7409.,7462.,7516.,7570.,7624.,7679.,7734.,7790.,7845.,7900.,
     $     7954.,8008.,8062.,8115.,8167.,8218.,8268.,8318.,8366.,8414.,
     $     8461.,8507.,8552.,8597.,8641.,8684.,8727.,8769.,8811.,8852.,
     $     8892.,8932.,8972.,9011.,9049.,9088.,9126.,9164.,9202.,9239.,
     $     9277.,9314.,9351.,9388.,9424.,9461.,9497.,9534.,9570.,9606.,
     $     9642.,9678.,9714.,9750.,9786.,9821.,9857.,9893.,9929.,9964.,
     $     10000./

      data fmxsm/
     $     341., 598., 776., 923.,1052.,1168.,1276.,1377.,1472.,1563.,
     $     1651.,1735.,1816.,1896.,1973.,2048.,2122.,2194.,2265.,2335.,
     $     2404.,2471.,2538.,2604.,2669.,2733.,2797.,2860.,2922.,2984.,
     $     3045.,3106.,3167.,3227.,3286.,3345.,3404.,3462.,3520.,3578.,
     $     3635.,3692.,3749.,3806.,3862.,3918.,3974.,4030.,4085.,4140.,
     $     4195.,4250.,4305.,4359.,4413.,4468.,4522.,4575.,4629.,4683.,
     $     4736.,4789.,4842.,4895.,4948.,5001.,5054.,5106.,5158.,5211.,
     $     5263.,5315.,5367.,5419.,5471.,5523.,5574.,5626.,5678.,5729.,
     $     5780.,5832.,5883.,5934.,5985.,6036.,6087.,6138.,6189.,6240.,
     $     6290.,6341.,6392.,6442.,6493.,6543.,6594.,6644.,6694.,6745.,
     $     6795.,6845.,6895.,6945.,6995.,7045.,7095.,7145.,7195.,7245.,
     $     7295.,7344.,7394.,7444.,7493.,7543.,7593.,7642.,7692.,7741.,
     $     7791.,7840.,7889.,7939.,7988.,8038.,8087.,8136.,8185.,8235.,
     $     8284.,8333.,8382.,8431.,8480.,8529.,8578.,8627.,8676.,8725.,
     $     8774.,8823.,8872.,8921.,8970.,9019.,9068.,9116.,9165.,9214.,
     $     9263.,9312.,9360.,9409.,9458.,9506.,9555.,9604.,9652.,9701.,
     $     9749./

      data dmin/
     $     1078.6,1133.,1148.,1158.,1166.,1172.,1178.,1183.,1187.,1191.,
     $     1195.,1198.,1202.,1205.,1208.,1211.,1213.,1216.,1219.,1221.,
     $     1224.,1226.,1229.,1231.,1234.,1236.,1239.,1242.,1245.,1248.,
     $     1251.,1255.,1258.,1262.,1265.,1269.,1273.,1277.,1282.,1286.,
     $     1291.,1296.,1302.,1308.,1315.,1323.,1331.,1341.,1351.,1363.,
     $     1376.,1391.,1408.,1427.,1449.,1475.,1505.,1538.,1568.,1594.,
     $     1618.,1640.,1660.,1680.,1699.,1718.,1736.,1754.,1772.,1789.,
     $     1805.,1820.,1835.,1849.,1862.,1875.,1887.,1899.,1911.,1923.,
     $     1935.,1947.,1960.,1973.,1986.,2000.,2014.,2029.,2044.,2059.,
     $     2075.,2091.,2107.,2124.,2140.,2157.,2174.,2190.,2207.,2224.,
     $     2240./

      data fmxdn/
     $     2221.,2297.,2373.,2446.,2519.,2590.,2660.,2730.,2798.,2865.,
     $     2932.,2998.,3063.,3128.,3192.,3255.,3318.,3380.,3442.,3504.,
     $     3565.,3625.,3685.,3745.,3804.,3863.,3922.,3981.,4039.,4096.,
     $     4154.,4211.,4268.,4325.,4382.,4438.,4494.,4550.,4606.,4661.,
     $     4717.,4772.,4827.,4882.,4936.,4991.,5045.,5100.,5154.,5208.,
     $     5261.,5315.,5369.,5422.,5475.,5529.,5582.,5635.,5687.,5740.,
     $     5793.,5845.,5898.,5950.,6003.,6055.,6107.,6159.,6211.,6263.,
     $     6315.,6366.,6418.,6470.,6521.,6573.,6624.,6675.,6727.,6778.,
     $     6829.,6880.,6931.,6982.,7033.,7084.,7135.,7185.,7236.,7287.,
     $     7337.,7388.,7438.,7489.,7539.,7590.,7640.,7690.,7740.,7791.,
     $     7841.,7891.,7941.,7991.,8041.,8091.,8141.,8191.,8241.,8291.,
     $     8341.,8390.,8440.,8490.,8539.,8589.,8639.,8688.,8738.,8787.,
     $     8837.,8886.,8936.,8985.,9035.,9084.,9133.,9183.,9232.,9281./

      data fsln/
     $     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,
     $     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,
     $     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,
     $     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,
     $     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5047.,
     $     5047.,5043.,5040.,5038.,5033.,5031.,5030.,5030.,5035.,5048.,
     $     5061.,5106.,5162.,5237.,5285.,5323.,5356.,5379.,5403.,5428.,
     $     5443.,5463.,5483.,5497.,5515.,5531.,5550.,5565.,5582.,5602.,
     $     5620.,5636.,5658.,5672.,5688.,5701.,5718.,5733.,5747.,5761.,
     $     5771.,5794.,5807.,5821.,5836.,5851.,5866.,5879.,5893.,5905.,
     $     5919.,5933.,5946.,5959.,5968.,5978.,5992.,6005.,6015.,6026.,
     $     6035.,6051.,6061.,6075.,6087.,6096.,6109.,6118.,6126.,6136.,
     $     6149.,6160.,6168.,6181.,6191.,6199.,6213.,6221.,6230.,6241.,
     $     6249.,6258.,6265.,6274.,6281.,6290.,6300.,6306.,6318.,6325.,
     $     6335.,6339.,6348.,6357.,6362.,6374.,6381.,6387.,6394.,6389.,
     $     6416.,6417.,6422.,6427.,6431.,6436.,6439.,6445.,6451.,6455.,
     $     6462.,6462.,6467.,6473.,6476.,6480.,6485.,6489.,6495.,6499.,
     $     6500.,6507.,6510.,6515.,6518.,6525./

      iraz=0
      w=sqrt(ecms2)
      tt=e-amnucl
      if(prot)then
         ka1=1
         ka2=0
         ka3=1
         ka4=-1
      else
         ka1=0
         ka2=1
         ka3=-1
         ka4=1
      end if
      if(no.eq.1)then
         t=tt
      else
         t=((w-upi)**2-am2n2)/2/amnucl
      end if
 1    call nuchis(t-300.,20.,161,ilo,ihi,erem)
      rf1=(frinn(ihi)*erem+frinn(ilo)*(1-erem))*1e-4
      phsm=(fmxsm(ihi)*erem+fmxsm(ilo)*(1-erem))*100.

 2    rf=rf1*ranf()
      call nuchis(rf*100,1.,101,ilo1,ihi1,erem1)
      amd1=dmin(ihi1)*erem1+dmin(ilo1)*(1-erem1)

      if(no.eq.1)then
c     ******** single production ***********
         ed=(amd1**2+ecms2-amnucl**2)/2/w
         pd2=ed**2-amd1**2
         if(pd2.lt.0.)go to 2
         pdlt=sqrt(pd2)
         phs=(w-ed)*ed*pdlt/w
         if(ranf().gt.phs/phsm)go to 2

         if(tt.lt.500.)then
            aiso=1.
         else if(tt.lt.1000.)then
            aiso=0.75
         else if(tt.lt.1300.)then
            aiso=0.50
         else if(tt.lt.2500.)then
            aiso=0.25
         else
            aiso=0
         end if
         if(ranf().lt.aiso)then
            cdl=2*ranf()-1.
         else
            cdl=-.9999995
            if(ranf().lt.0.5)cdl=-cdl
         end if
         call cone(up4cms,dldir,cdl,-1.0)
         do i1=1,3
            pdl(i1)=pdlt*dldir(i1)
            pdldec(i1,3)=-pdl(i1)
         end do
         pdl(4)=ed
         pdldec(4,3)=w-ed
         call twobiso(amd1,amnucl,upi,pdl,pdldec)
         do i1=1,3
            call lloren(piso(1,i1),piso(4,i1),
     &           pdldec(1,i1),pdldec(4,i1),beta)
         end do
         uiso(1) = amnucl
         uiso(2) = upi
         uiso(3) = amd1
         rpi=ranf()
         rnu=ranf()
         if((ka-ka1).eq.0)then  !p-p n-n
            if(rpi.lt.0.75)then
               chiso(1)=ka1
               chiso(2)=ka3
               chiso(3)=ka2
            else if(rpi.lt.0.9167)then
               chiso(1)=ka1
               chiso(2)=0
               chiso(3)=ka1
            else
               chiso(1)=ka2
               chiso(2)=ka3
               chiso(3)=ka1
            end if
         else                   !p-n n-p
            if(rpi.lt.0.66667)then
               chiso(2)=0
               if(rnu.lt.0.5)then
                  chiso(1)=ka
                  chiso(3)=ka1
               else
                  chiso(1)=ka1
                  chiso(3)=ka
               end if
            else
               if(rnu.lt.0.5)then
                  chiso(1)=1
                  chiso(2)=-1
                  chiso(3)=1
               else
                  chiso(1)=0
                  chiso(2)=1
                  chiso(3)=0
               end if
            end if
         end if
c     ****************** end of single production ****************
      else
         if(iraz.eq.0)then      !two deltas
            iraz=1
            amd2=amd1
            rf2=rf1
            t=((w+amnucl-amd1)**2-am2n2)/2/amnucl
            go to 1
         else
            ed=(amd1**2-amd2**2+ecms2)/2/w
            pd2=ed**2-amd1**2
            if(pd2.lt.0)then
               iraz=0
               rf1=rf2
               go to 2
            end if
            pdlt=sqrt(pd2)
            phs=(w-ed)*ed*pdlt/w
            call nuchis(tt-920.,20.,130,ilo1,ihi1,erem1)
            phsm=(fmxdn(ihi1)*erem1+fmxdn(ilo1)*(1-erem1))*100.
            if(ranf().gt.phs*rf1/(phsm*rf2))then
               iraz=0
               rf1=rf2
               go to 2
            end if

            cdl=-.9999995
            if(ranf().lt.0.5)cdl=-cdl
            call cone(up4cms,dldir,cdl,-1.0)
            do i1=1,3
               pdl(i1)=pdlt*dldir(i1)
            end do
            pdl(4)=ed
            call twobiso(amd1,amnucl,upi,pdl,pdldec)
            do i1=1,3
               pdl(i1)=-pdl(i1)
            end do
            pdl(4)=w-ed
            call twobiso(amd2,amnucl,upi,pdl,pdldec(1,3))
            do i1=1,4
               call lloren(piso(1,i1),piso(4,i1),
     &              pdldec(1,i1),pdldec(4,i1),beta)
            end do
            uiso(1) = amnucl
            uiso(2) = upi
            uiso(3) = amnucl
            uiso(4) = upi
            rpi=ranf()
            rnu=ranf()
            if((ka-ka1).eq.0)then !p-p n-n
               if(rpi.lt.0.6)then
                  chiso(1)=ka1
                  chiso(2)=ka3
                  if(rnu.lt.0.3333)then
                     chiso(3)=ka1
                     chiso(4)=ka4
                  else
                     chiso(4)=0
                     chiso(3)=ka2
                  end if
               else if(rpi.lt.0.8667)then
                  chiso(1)=ka1
                  chiso(2)=0
                  if(rnu.lt.0.3333)then
                     chiso(3)=ka2
                     chiso(4)=ka3
                  else
                     chiso(3)=ka1
                     chiso(4)=0
                  end if
               else
                  chiso(1)=ka2
                  chiso(2)=ka3
                  if(rnu.lt.0.3333)then
                     chiso(3)=ka2
                     chiso(4)=ka3
                  else
                     chiso(3)=ka1
                     chiso(4)=0
                  end if
               end if
            else                !p-n n-p
               call nuchis(tt,20.,176,ilo,ihi,erem)
               rat=(fsln(ihi)*erem+fsln(ilo)*(1-erem))*1e-4
               if(rpi.lt.rat)then
                  chiso(1)=ka1
                  chiso(2)=ka3
                  chiso(3)=ka2
                  chiso(4)=ka4
               else if(rpi.lt.(1.-rat)*0.6667)then
                  chiso(1)=ka1
                  chiso(2)=0
                  if(rnu.lt.0.3333)then
                     chiso(3)=ka1
                     chiso(4)=ka4
                  else
                     chiso(3)=ka2
                     chiso(4)=0
                  end if
               else
                  chiso(1)=ka2
                  chiso(2)=ka3
                  if(rnu.lt.0.3333)then
                     chiso(3)=ka1
                     chiso(4)=ka4
                  else
                     chiso(3)=ka2
                     chiso(4)=0
                  end if
               end if
            end if
         end if
      end if
      return
      end

      SUBROUTINE nucRODIS(RMIN,DR,NBIN,PROBM,ROT)
      DIMENSION ROT(NBIN),PROBM(NBIN)
      R=RMIN-DR
      DO i1=1,NBIN
         R=R+DR
         ROT(i1)=ROXY(R)
         PROBM(i1)=4*3.1416*((R+DR)**3-R**3)/3
      END DO
      i2=NBIN-1
      DO i1=1,i2
         PROBM(i1)=PROBM(i1)*(ROT(i1)+ROT(i1+1))/2
      END DO
      PROBM(NBIN)=PROBM(NBIN)*ROT(NBIN)
      RETURN
      END

C*******************************************************
C     Find a starting point inside the nucleus for the hadron.
      subroutine nucstp(start,icont,mode)
      dimension start(3)
      common /nucdis/ probmi(31),probdi(31),rhotab(30),pnorm
      parameter ( rbin=.2, rmin=0, rmax=5.99, nbin=30)

      if(icont.eq.1)then
         do i1=1,3
            start(i1)=0
         end do
      else if(icont.eq.0)then
         xx=ranf()
         if(mode.eq.1)then      !surphace mode
            r=rmax*sqrt(xx)
            phi=360*ranf()
            start(1)=r*cosd(phi)
            start(2)=r*sind(phi)
            start(3)=-rmax
         else if(mode.eq.2)then !volume mode
            r=prbin(probmi,nbin,rmin,rbin,xx)
            call ranve(r,start,2.,-1.)
         else                   !Dover mode
            r=prbin(probdi,nbin,rmin,rbin,xx)
            call ranve(r,start,2.,-1.)
         end if
      end if
      return
      end

      subroutine nucvmes(jpar,ides)
      common /nucint/uin,chint,ein,pin(3),upin(3)
      common /nuccms/ rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
      common /nuchad/ihap,nel,ncx
#include "hadron.h"
      logical prot
      dimension pout1l(3),poutel(3),rv(3)
      parameter (upi=139.6,  step=0.2)

      real xeta(20),xome(20)
C     PARTIAL AND TOTAL CROSS SECTIONS
C     PI-P TO ETA+N CROSS SECTION (MB) AT 20*20MEV/C VALUES
C     OF THE CMS MOM. OF THE ETA
      DATA XETA/
     $     0.,.1,.25,.55,1.,1.25,1.4,1.6,1.7,1.75,
     $     1.8,1.8,1.75,1.7,1.6,1.4,1.3,1.15,1.0,.9/
C     PI-P TO OME+N CROSS       SECTION (MB) AT 20*20MEV/C VALUES
C     OF THE CMS MOM. OF THE OME
      DATA XOME/
     $     0.,.1,.15,.28,.45,.65,.90,1.05,1.20,1.35,
     $     1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85/


      ides=1
      IF(ECMS2.LE.(UPI+UNUCL)**2) GO TO 200
      IF(ECMS2.LE.(uin+unucl)**2) GO TO 200

      PCMSM2=SQRT((ECMS2-(UPI-UNUCL)**2)*(ECMS2-(UPI+UNUCL)**2)/(4.0*ECMS2))
      PCMSM1=SQRT((ECMS2-(UIn-UNUCL)**2)*(ECMS2-(UIn+UNUCL)**2)/(4.0*ECMS2))
      pinm=sqrt(pin(1)**2+pin(2)**2+pin(3)**2)

      call nuchis(pcmsm1,20.,20,ilo,ihi,prem)
      if(jpar.eq.JETA)then
         PINT=1.5/.71*
     $        (XETA(IHI)*PREM+XETA(ILO)*(1.-PREM))*RHON*(PCMSM2/PCMSM1)**2
         PDEC=UIn/(PINM*230000.0)
         if(pcmsm1.lt.100.)then
            pela=.43*pcmsm1
         else if(pcmsm1.gt.270.)then
            pela=156.-.42*pcmsm1
            if(pela.lt.0.)pela=0.
         else
            pela=43.
         end if
         pela=pela*rhon
         PTOT=PDEC+PINT+pela
      else if (jpar.eq.JOMEGA) then
         pint=0.5*
     $        (xome(ihi)*prem+xome(ilo)*(1.-prem))*RHON*(PCMSM2/PCMSM1)**2
         pdec=uin/(20.*pinm)
         pela=0
         ptot=pdec+pint
      else
         print *, 'Not a vector meson', jpar
      end if

      FINCH=0.6666*pint
C     RHON FACTOR ALREADY IN

      XR=ranf()+1.e-6
      YR=ranf()*ptot
      DIST=-LOG(XR)/ptot

C**********IFNOTHING HAPPENS
      IF(STEP.LT.DIST)GO TO 200
C**********IFSOMETHING HAPPENS
c
C**********DOESIT DECAY?
      FSUM=pDEC
      IF(YR.LT.FSUM)THEN
         ihap=1
         IDEC=IDEC+1
         GO TO 99
      END IF

C     VECTOR MESON MUST HAVE INTERACTED

      call ranve(1.,rv,2.,-1.)

C****************DOES VECTOR    MESON ELASTICALLY SCATTER?
      FSUM=FSUM+pELA
      IF(YR.LT.FSUM)THEN

         call nuckin(rv,pcmsm1,uin,poutnlm,poutel,poutelm,eoutel)
         IF(POUTNLM.LT.PFERMI) go to 200
         DO i2=1,3
            PIN(i2)=POUTEL(i2)
            UPIN(i2)=PIN(i2)/POUTELM
         END DO
         EIN=EOUTEL
         ihap=1
         nel=1
         GO TO 200
      END IF

C****************VECTOR MESON MAKES A PION*********

      call nuckin(rv,pcmsm2,upi,pout2lm,pout1l,pout1lm,eout1l)
      IF(POUT2LM.LT.PFERMI)GO TO 200
      ihap=1
      INEL=INEL+1
C*****IS SCATTERING OFF P OR N?
      FSUM=FSUM+FINCH
      if(prot)then
C*****PROTON
         IF(YR.LT.FSUM)THEN
            CHINT=1.0
         ELSE
            CHINT=0.0
         END IF
      ELSE
C*****NEUTRON
         IF(YR.LT.FSUM) THEN
            CHINT=-1.0
         ELSE
            CHINT=0.0
         END IF
      END IF
      DO i2=1,3
         PIN(i2)=POUT1L(i2)
         UPIN(i2)=PIN(i2)/POUT1LM
      END DO
      EIN=EOUT1L
      UIn=upi
 99   ides=0
 200  return
      end

      subroutine AZIO(COSINE,SINE)
    5 r1=ranf()
      R1SQ=R1*R1
C     XSQ
      r2=ranf()
      R2SQ=R2*R2
C     YSQ
      SUM=R1SQ+R2SQ
      IF(SUM.GT.1.0)GOTO5
      SUM=SUM*0.5
C     (XSQ+YSQ)/2
      COSINE=(SUM-R1SQ)/SUM
C     (YSQ-XSQ)/(XSQ+YSQ)
      SINE=(R1*R2)/SUM
C     (2*X*Y)/(XSQ+YSQ)
      r1=ranf()
      IF(R1.LT.0.5)GOTO10
      SINE=-SINE
   10 RETURN
      END

C**************************************************************
C     Discribe the volume mode of the nucleus.  The hadron is assumed
C     to originate from a nucleon inside the nucleus.
      REAL FUNCTION ROXY(R)
      parameter (w=-.05, c=2.61, z=.513, d=.102, g=2.76, e=.35)
      RO=(1+W*(R/C)**2)/(1+EXP((R-C)/Z))
      IF(R.EQ.0)THEN
         DRO=1
      ELSE
         GR=G*R
         DRO=SIN(GR)/GR
      END IF
      DRO=DRO*D*EXP(-(E*R)**2)
      ROXY=RO-DRO
      RETURN
      END

C Decay atomic mass AMD with momentum PD into masses AM1 and AM2 with
C momentum P1(,1) and P1(,2). Then transform the P1's to lab frame.
C (To read this comment, please add 1 grain, salt. -bviren)
      SUBROUTINE TWOBiso(AMD,AM1,AM2,PD,P1)
      DIMENSION PD(4),P1(4,2),PTEM(4,2),dir(3)
      AMD2=AMD**2
      A2M1=AM1**2
      A2M2=AM2**2
      E1=(AMD2+A2M1-A2M2)/AMD/2.
      E2=(AMD2+A2M2-A2M1)/AMD/2.
      P=E1**2-A2M1
      P=SQRT(P)
      call azio(cs,ss)                          ! Get random cos/sin
      pdd=sqrt(pd(1)**2+pd(2)**2+pd(3)**2)
      do i1=1,3
         dir(i1)=pd(i1)/pdd
      end do
      CALL cone(dir,ptem,cs,-1.0)
      do i1=1,3
         ptem(i1,1)=ptem(i1,1)*p
         ptem(i1,2)=-ptem(i1,1)
      end do
      PTEM(4,1)=E1
      PTEM(4,2)=E2
      CALL LOREN(PD,PTEM,P1,2,-1)
      RETURN
      END

C****************************************************
C   Handle kaon nucleon interactions in the range 1.5-2Gev Center of Mass
C     Energy.  This is only set up to handle K+N->m+b where m is a meson and
C     b is a baryon.
C**** The charge of the kaon controls if it is a k+ k0 k0b or a k-
C     k+ > 0.9
C     k0 > 0.0
C     k0b > -0.9
C     k- < -0.9
      subroutine nuckaon(ides,step)
      implicit none
C**** Take one step in the nucleus.
      integer i1, i2, i3
      real r1, r2, r3
C**** The index of the incoming particle.  On output this is:
C     1 - the particle leaves the step.
C     0 - the particle was absorbed.
C    -1 - the particle created a second particle.
      integer ides
C**** The step size for this simulation in fermi.  Usually this is .2 fm.
      real step
C**** The data about the interaction particle
      real uin, chint, ein, pin(3), upin(3)
      common /nucint/uin,chint,ein,pin,upin
#include "nucnewprd.h"
#include "masses.h"
C**** The nuclear center of mass parameters.
      real up4cms(3), beta(3)
      logical prot
      real rhon, pfermi, unucl, ecms2
      common /nuccms/rhon,pfermi,unucl,ecms2,up4cms,beta,prot
C**** The center of mass energy
      real Ecm
C**** The center of mass momentum of the Meson.
      real Pmes(4), PKaon
C**** The mass of the meson produced.
      real Mmes
C**** The charge of the meson produced.
      real Qmes
C**** The center of mass momentum of the baryon.
      real Pbar(4)
C**** The mass of the baryon produced.
      real Mbar
C**** The charge of the baryon produced.
      real Qbar
C**** The number of partial waves in the expansion.
      integer nwaves
      parameter (nwaves=8)
C**** The total partial wave coefficients for kaon+nuc -> kaon+nuc
C     note that <1,0|1,0> = <1,1|1,1> = <1,-1|1,-1> = knkn1
      complex knkn0(nwaves), knkn1(nwaves)
C**** The cross sections.
      integer ncross
      parameter (ncross=3)
      real cross(ncross)
C**** The partial wave coefficients for the interaction.
      complex tab(nwaves,ncross)
C**** The total cross section.
      real tcross
C**** The minimum free path found.
      real free
C**** The index of the minimum free path.
      integer ifree
      real ranf
C**** The beta to go from the lab to the center of mass.
      real bcms(3)

      ides = 1

      if (chint.eq.0) then
C**** Make zero charge kaons into k0 not k0b
         chint = 0.1
      endif

C**** Find the center of mass energy from the energy squared found by
C     nucfermi.
      Ecm = sqrt(Ecms2)/1000.0

C**** Find the momentum of the kaon in the center of mass.
      do i1 = 1,3
         Pmes(i1) = pin(i1)/1000.0
      enddo
      Pmes(4) = ein/1000.0

      do i1 = 1,3
         bcms(i1) = -beta(i1)
      enddo
      call lloren(pmes,pmes(4),pmes,pmes(4),bcms)
      Pkaon = sqrt(Pmes(1)**2+Pmes(2)**2+Pmes(3)**2)

C**** Find the partial wave coefficients for each isospin
C**** Need switch for kaons <---> antikaons
      if(chint.gt.0) then
         call nucknknp0(Ecm,knkn0)
         call nucknknp1(Ecm,knkn1)
      else
         call nucknkn0(Ecm,knkn0)
         call nucknkn1(Ecm,knkn1)
      endif
C**** Find the <11|11> scattering matrix
C     This is for:
C          K+    + p -> K+    + p
C          K0    + n -> K0    + n
C          K0bar + p -> K0bar + p
C          K-    + n -> K-    + n
      do i1 = 1,nwaves
         tab(i1,1) = knkn1(i1)
      enddo

C**** Find the .5(<10|10>+<00|00>) cross section.
C     This is for:
C          K+    + n -> K+    + n
C          K0    + p -> K0    + p
C          K0bar + n -> K0bar + n
C          K-    + p -> K-    + p
      do i1 = 1,nwaves
         tab(i1,2) = 0.5*(knkn1(i1)+knkn0(i1))
      enddo

c**** Find the 0.5(<10|10>-<00|00>) cross section.
C     This is for K0+n -> K-+p and K-+p -> K0+n
C          K+    + n -> K0    + p
C          K0    + p -> K+    + n
C          K0bar + n -> K-    + p
C          K-    + p -> K0bar + n
      do i1 = 1,nwaves
         tab(i1,3) = 0.5*(knkn1(i1)-knkn0(i1))
      enddo

C**** Turn the partial wave amplitudes in to total cross sections.
      do i1 = 1,ncross
         call nuckcrs(PKaon,tab(1,i1),cross(i1))
      enddo

C**** Convert the cross sections into free paths and find the one that is
C     the shortest.
      free = 99999999.9
      ifree = 0
      do i1 = 1,ncross
         r1 = ranf()
         if (r1.gt.1.0E-7) then
            r1 = -log(r1)/cross(i1)/rhon
         else
            r1 = free + 1.0
         endif
         if (r1.lt.free) then
            free = r1
            ifree = i1
         endif
      enddo

      if ((ifree.eq.0).or.(free.gt.step)) then
         ides = 1
         return
      endif
      if(amkc.ne.0.493646) then
         amkc = 0.493646
      endif
      if(amko.ne.0.497671) then
         amko = 0.497671
      endif
      if(amp.ne.0.9382798) then
         amp = 0.9382798
      endif
      if(amn.ne.0.9395731) then
         amn = 0.9395731
      endif



C**** Figure out the final products.
      if (ifree.eq.1) then
         if (chint.gt.0.9) then
C**** This is a Kplus
            Qmes = 1.0
            Mmes = amkc
            Qbar = 1.0
            Mbar = amp
         else if (chint.gt.0) then
C**** This is a K0
            Qmes = 0.1
            Mmes = amko
            Qbar = 0.0
            Mbar = amn
         else if (chint.gt.-0.9) then
C**** This is a K0bar
            Qmes = -0.1
            Mmes = amko
            Qbar = 1.0
            Mbar = amp
         else
C**** This is a Kminus
            Qmes = -1.0
            Mmes = amkc
            Qbar = 0.0
            Mbar = amn
         endif
      else if (ifree.eq.2) then
         if (chint.gt.0.9) then
C**** This is a Kplus
            Qmes = 1.0
            Mmes = amkc
            Qbar = 0.0
            Mbar = amn
         else if (chint.gt.0) then
C**** This is a K0
            Qmes = 0.1
            Mmes = amko
            Qbar = 1.0
            Mbar = amp
         else if (chint.gt.-0.9) then
C**** This is a K0bar
            Qmes = -0.1
            Mmes = amko
            Qbar = 0.0
            Mbar = amn
         else
C**** This is a Kminus
            Qmes = -1.0
            Mmes = amkc
            Qbar = 1.0
            Mbar = amp
         endif
      else if (ifree.eq.3) then
         if (chint.gt.0.9) then
C**** This is a Kplus
            Qmes = 0.1
            Mmes = amko
            Qbar = 1.0
            Mbar = amp
         else if (chint.gt.0) then
C**** This is a K0
            Qmes = 1.0
            Mmes = amkc
            Qbar = 0.0
            Mbar = amn
         else if (chint.gt.-0.9) then
C**** This is a K0bar
            Qmes = -1.0
            Mmes = amkc
            Qbar = 1.0
            Mbar = amp
         else
C**** This is a Kminus
            Qmes = -0.1
            Mmes = amko
            Qbar = 0.0
            Mbar = amn
         endif
      else
C**** No interaction this step.
         print *, 'INVALID CROSS SECTION IN NUCKAON.'
         stop
      endif

C**** Find the direction of the meson after the scatter.  Notice that
C     Pmes goes in a momentum and comes out a direction.
      call nuckdir(Pmes,Pkaon,tab,Pmes)

C**** Convert Pmes back into a 4 momentum.
      Pkaon = ecm**2 - Mmes**2 - Mbar**2
      if (Pkaon.lt.0.0) then
         ides = 1
         stop
      endif
      Pkaon = sqrt(Pkaon)/2.0

      do i1 = 1,3
         PMes(i1) = PMes(i1)*PKaon
      enddo
      PMes(4) = sqrt(Mmes**2+PKaon**2)

      call lloren(Pmes,Pmes(4),Pmes,Pmes(4),beta)

      uin = Mmes*1000.0
      ein = Pmes(4)*1000.0
      chint = Qmes
      r1 = 0.0
      do i1 = 1,3
         pin(i1) = Pmes(i1)*1000.0
         r1 = r1 + pin(i1)**2
      enddo
      r1 = sqrt(r1)
      do i1 = 1,3
         upin(i1) = pin(i1)/r1
      enddo

      return
      end

C**** Return the complex legendre coefficients for center of mass energy enr.
      subroutine nucclxint(enr,tab,energy,coeffs,nenr)
      implicit none
      integer i1
C**** The energy to get the coefficients for.
      real enr
C**** The return variable for the coefficients.
      complex tab(8)
C**** The fraction distance in the energy bin.
      real frac
C**** The energy bin below the energy.
      integer low
C**** The number of energy bins.
      integer nenr
C**** The energy bins and the coefficients.
      real energy(nenr)
      complex coeffs(8,nenr)

C**** Do a simple minded search for the correct energy bin.  If this runs
C     too slow then here is where to speed it up.
      if (enr.le.energy(1)) then
         low = 1
         frac = 0.0
      else if (enr.ge.energy(nenr)) then
         low = nenr-1
         frac = 1.0
      else
         low = 1
         do while (enr.gt.energy(low))
            low = low + 1
         enddo
         low = low - 1
         frac = (enr-energy(low))/(energy(low+1)-energy(low))
      endif

C**** Do a simple linear interpolation between the bins.
      do i1 = 1,8
         tab(i1) = (1.0-frac)*(coeffs(i1,low)) + frac*(coeffs(i1,low+1))
      enddo

      return
      end

      subroutine nuckcrs(PKaon,tab,cros)
      implicit none
      integer i1
C**** The center of mass momentum.
      real PKaon
C**** The partial wave ampitudes for the cross section.
      complex tab(*)
C**** The return value for the cross section.
      real cros
C*** The factor of 3.0 next to the pi is a fudge factor
C*** Shouldn't it be 4.0???
C*** This matches the cross sections.
      cros = 0.0
      do i1 = 1,8
         cros = cros + 3.0*3.14159*(2*(i1-1)+1)*tab(i1)*conjg(tab(i1))
      enddo
      cros = cros/PKaon/PKaon

      return
      end

      subroutine nuckdir(PMes,PKaon,Tab,Pout)
      implicit none
      integer i1, i2
      real r1, r2
C**** The 4 momentum of the kaon in the center of mass frame.
      real Pmes(4)
C**** The momentum of the kaon in the center of mass frame.
      real PKaon
C**** The partial wave amplitudes for the cross section.
      complex Tab(*)
C**** The output momentum in the center of mass frame.
      real Pout(4)
C**** The table used to invert the differential cross section.
      integer ncrs
      parameter (ncrs=20)
      real Dcrs(ncrs)
      real ranf
C**** The cosine
      real css

C**** Find the differential cross section at ncrs different grid points.
      do i1 = 1,ncrs
         r1 = i1*2.0/ncrs - 1.0
         call nucdcrs(PKaon,tab,r1,Dcrs(i1))
      enddo

C**** Sum the cross section to get the total.  Notice that this is not
C     normalized to the total cross section, but it doesn't matter since
C     we have already decided that this interaction is taking place and
C     are just generating the angular distribution.
      do i1 = 2,ncrs
         Dcrs(i1) = Dcrs(i1) + Dcrs(i1-1)
      enddo

C**** Find the cosine
      r1 = 0.999*ranf()*Dcrs(Ncrs)
      i1 = 1
      do while ((r1.lt.Dcrs(i1)).and.(i1.lt.Ncrs))
         i1 = i1 + 1
      enddo

C**** Find the fractional diffenence from the cosine point found and the
C     one wanted.
      r1 = r1 - Dcrs(i1)
      if (i1.eq.1) then
         r1 = r1/Dcrs(i1)
      else
         r1 = r1/(Dcrs(i1)-Dcrs(i1-1))
      endif

      css = (i1+r1)*2.0/ncrs - 1.0

      call cone(Pmes,Pout,css,-1.0)

      return
      end

      subroutine nucdcrs(PKaon,tab,dcos,cros)
C**** Return the differential cross section as a function of cosine.
      implicit none
      integer i1
C**** The center of mass momentum.
      real PKaon
C**** The partial wave amplitudes for the cross section.
      complex tab(*)
C**** The cosine to return the differential cross section for.
      real dcos
C**** The cross section.
      real cros
C**** The wave function
      complex ff
C**** The legendre polynomial
      real plgndr

      ff = (0.0,0.0)
      do i1 = 1,8
         call nucplgndr(i1-1,dcos,plgndr)
         ff = ff + (2.0*(i1-1.0)+1.0)*tab(i1)*plgndr
      enddo
      ff = ff/PKaon

      cros = abs(ff)**2

      return
      end

C***************************************
C   From Numerical Recipes ch 6.6
C
      subroutine NUCPLGNDR(L,X,plgndr)
C**** Computes the associated legendre polynomial P(l)(x) where l is an
c     integer and -1<x<1.
      if (l.lt.0.or.abs(x).gt.1.) then
         print *, 'bad arguments to plgndr',l,m,x
         stop
      endif
      pmm = 1.0
      if (l.eq.0) then
         plgndr = pmm
      else
         pmmp1 = x*1.*pmm
         if (l.eq.1) then
            plgndr = pmmp1
         else
            do ll = 2,l
               pll = (x*(2*ll-l)*pmmp1-(ll-1)*pmm)/(ll)
               pmm = pmmp1
               pmmp1 = pll
            enddo
            plgndr = pll
         endif
      endif

      return
      end

C**** Return the complex legendre coefficients for center of mass energy enr.
C     B.R.Martin and M.K.Pidcock, Nuc.Phys. B126 (1977) 266-284
C     B.R.Martin and M.K.Pidcock, Nuc.Phys. B126 (1977) 285-297
C        Kbar-N Interactions in the Resonance Region. (I) Analysis of data
C             and (II) Amplitudes
      subroutine nucknkn0(enr,tab)
      implicit none
C**** The energy to get the coefficients for.
      real enr
C**** The return variable for the coefficients.
      complex tab(8)
C**** The number of energy bins.
      integer nenr
      parameter (nenr=36)
C**** The energy bins and the coefficients.
      real energy(nenr)
      complex coeffs(8,nenr)
      data energy /
     $     1.550, 1.560, 1.570, 1.580, 1.590,
     $     1.600, 1.610, 1.620, 1.630, 1.640,
     $     1.650, 1.660, 1.670, 1.680, 1.690,
     $     1.700, 1.710, 1.720, 1.730, 1.740,
     $     1.750, 1.760, 1.770, 1.780, 1.790,
     $     1.800, 1.820, 1.840, 1.860, 1.880,
     $     1.900, 1.920, 1.940, 1.960, 1.980,
     $     2.000 /
      data coeffs /
     $     ( 0.247, 0.434),( 0.233, 0.236),( 0.117, 0.022),(-0.353, 0.173),
     $     ( 0.025, 0.026),( 0.003, 0.001),( 0.002, 0.001),( 0.007, 0.000),
     $     ( 0.275, 0.452),( 0.230, 0.261),( 0.124, 0.025),(-0.337, 0.155),
     $     ( 0.026, 0.027),( 0.005, 0.001),( 0.002, 0.001),( 0.008, 0.000),
     $     ( 0.296, 0.472),( 0.223, 0.284),( 0.131, 0.027),(-0.321, 0.139),
     $     ( 0.026, 0.028),( 0.007, 0.001),( 0.002, 0.001),( 0.008, 0.000),
     $     ( 0.312, 0.495),( 0.214, 0.305),( 0.138, 0.030),(-0.304, 0.125),
     $     ( 0.026, 0.028),( 0.010, 0.001),( 0.002, 0.001),( 0.009, 0.000),
     $     ( 0.325, 0.519),( 0.202, 0.323),( 0.146, 0.034),(-0.287, 0.113),
     $     ( 0.026, 0.029),( 0.013, 0.001),( 0.003, 0.001),( 0.009, 0.000),
     $     ( 0.335, 0.546),( 0.190, 0.339),( 0.153, 0.037),(-0.269, 0.102),
     $     ( 0.026, 0.029),( 0.016, 0.002),( 0.003, 0.001),( 0.010, 0.000),
     $     ( 0.342, 0.575),( 0.177, 0.351),( 0.161, 0.041),(-0.249, 0.092),
     $     ( 0.026, 0.030),( 0.019, 0.002),( 0.003, 0.001),( 0.010, 0.000),
     $     ( 0.345, 0.606),( 0.164, 0.361),( 0.169, 0.045),(-0.227, 0.084),
     $     ( 0.027, 0.030),( 0.023, 0.003),( 0.003, 0.001),( 0.011, 0.001),
     $     ( 0.345, 0.639),( 0.152, 0.368),( 0.177, 0.050),(-0.202, 0.080),
     $     ( 0.027, 0.030),( 0.028, 0.003),( 0.003, 0.001),( 0.011, 0.001),
     $     ( 0.337, 0.675),( 0.142, 0.374),( 0.185, 0.055),(-0.173, 0.080),
     $     ( 0.027, 0.031),( 0.033, 0.004),( 0.003, 0.001),( 0.012, 0.001),
     $     ( 0.315, 0.712),( 0.132, 0.379),( 0.194, 0.061),(-0.139, 0.090),
     $     ( 0.027, 0.031),( 0.039, 0.005),( 0.003, 0.001),( 0.013, 0.001),
     $     ( 0.244, 0.711),( 0.124, 0.383),( 0.203, 0.067),(-0.103, 0.114),
     $     ( 0.028, 0.031),( 0.046, 0.006),( 0.003, 0.001),( 0.013, 0.001),
     $     ( 0.339, 0.664),( 0.118, 0.387),( 0.212, 0.074),(-0.076, 0.162),
     $     ( 0.028, 0.031),( 0.054, 0.008),( 0.004, 0.001),( 0.014, 0.001),
     $     ( 0.348, 0.747),( 0.113, 0.391),( 0.221, 0.081),(-0.077, 0.224),
     $     ( 0.028, 0.031),( 0.063, 0.010),( 0.004, 0.001),( 0.015, 0.001),
     $     ( 0.329, 0.797),( 0.109, 0.395),( 0.231, 0.089),(-0.130, 0.285),
     $     ( 0.029, 0.032),( 0.073, 0.013),( 0.004, 0.001),( 0.015, 0.001),
     $     ( 0.305, 0.836),( 0.105, 0.400),( 0.240, 0.099),(-0.200, 0.294),
     $     ( 0.029, 0.032),( 0.085, 0.017),( 0.004, 0.001),( 0.016, 0.001),
     $     ( 0.279, 0.870),( 0.102, 0.406),( 0.250, 0.109),(-0.252, 0.267),
     $     ( 0.030, 0.032),( 0.099, 0.022),( 0.004, 0.001),( 0.017, 0.001),
     $     ( 0.249, 0.901),( 0.099, 0.413),( 0.260, 0.121),(-0.279, 0.230),
     $     ( 0.031, 0.032),( 0.115, 0.029),( 0.005, 0.002),( 0.018, 0.001),
     $     ( 0.216, 0.928),( 0.096, 0.420),( 0.271, 0.135),(-0.288, 0.195),
     $     ( 0.032, 0.032),( 0.134, 0.039),( 0.005, 0.002),( 0.019, 0.002),
     $     ( 0.181, 0.951),( 0.093, 0.428),( 0.280, 0.150),(-0.288, 0.165),
     $     ( 0.033, 0.032),( 0.156, 0.053),( 0.005, 0.002),( 0.020, 0.002),
     $     ( 0.142, 0.969),( 0.089, 0.437),( 0.290, 0.167),(-0.282, 0.141),
     $     ( 0.035, 0.032),( 0.181, 0.073),( 0.005, 0.002),( 0.021, 0.002),
     $     ( 0.102, 0.983),( 0.085, 0.446),( 0.298, 0.186),(-0.274, 0.121),
     $     ( 0.037, 0.033),( 0.210, 0.102),( 0.005, 0.002),( 0.022, 0.002),
     $     ( 0.060, 0.991),( 0.080, 0.456),( 0.306, 0.208),(-0.264, 0.104),
     $     ( 0.039, 0.034),( 0.241, 0.144),( 0.006, 0.002),( 0.023, 0.002),
     $     ( 0.017, 0.994),( 0.074, 0.467),( 0.311, 0.233),(-0.254, 0.090),
     $     ( 0.042, 0.035),( 0.269, 0.206),( 0.006, 0.002),( 0.024, 0.003),
     $     (-0.025, 0.992),( 0.067, 0.477),( 0.313, 0.260),(-0.242, 0.077),
     $     ( 0.046, 0.039),( 0.283, 0.294),( 0.006, 0.003),( 0.025, 0.003),
     $     (-0.066, 0.984),( 0.059, 0.488),( 0.312, 0.289),(-0.230, 0.066),
     $     ( 0.049, 0.046),( 0.261, 0.406),( 0.006, 0.003),( 0.027, 0.004),
     $     (-0.142, 0.953),( 0.040, 0.509),( 0.295, 0.352),(-0.204, 0.047),
     $     ( 0.041, 0.065),( 0.038, 0.582),( 0.007, 0.004),( 0.030, 0.004),
     $     (-0.204, 0.904),( 0.017, 0.529),( 0.253, 0.409),(-0.174, 0.033),
     $     ( 0.022, 0.063),(-0.207, 0.506),( 0.007, 0.005),( 0.034, 0.005),
     $     (-0.247, 0.843),(-0.010, 0.545),( 0.191, 0.447),(-0.139, 0.024),
     $     ( 0.017, 0.055),(-0.291, 0.357),( 0.008, 0.006),( 0.038, 0.007),
     $     (-0.269, 0.778),(-0.040, 0.558),( 0.123, 0.457),(-0.098, 0.023),
     $     ( 0.016, 0.049),(-0.293, 0.256),( 0.009, 0.008),( 0.044, 0.009),
     $     (-0.268, 0.714),(-0.072, 0.567),( 0.067, 0.442),(-0.053, 0.035),
     $     ( 0.016, 0.046),(-0.277, 0.197),( 0.010, 0.012),( 0.051, 0.012),
     $     (-0.249, 0.661),(-0.105, 0.570),( 0.027, 0.416),(-0.007, 0.064),
     $     ( 0.017, 0.044),(-0.259, 0.162),( 0.012, 0.016),( 0.059, 0.017),
     $     (-0.215, 0.623),(-0.138, 0.569),( 0.005, 0.387),( 0.027, 0.114),
     $     ( 0.017, 0.043),(-0.244, 0.141),( 0.015, 0.020),( 0.070, 0.025),
     $     (-0.176, 0.604),(-0.169, 0.564),(-0.006, 0.362),( 0.037, 0.177),
     $     ( 0.018, 0.042),(-0.232, 0.129),( 0.015, 0.026),( 0.083, 0.038),
     $     (-0.138, 0.604),(-0.199, 0.554),(-0.010, 0.343),( 0.018, 0.238),
     $     ( 0.018, 0.041),(-0.226, 0.123),( 0.012, 0.034),( 0.099, 0.061),
     $     (-0.108, 0.618),(-0.225, 0.542),(-0.010, 0.329),(-0.021, 0.280),
     $     ( 0.019, 0.040),(-0.226, 0.120),( 0.005, 0.040),( 0.110, 0.103) /

      call nucclxint(enr,tab,energy,coeffs,nenr)

      return
      end

C**** Return the complex legendre coefficients for center of mass energy enr.
C     B.R.Martin and M.K.Pidcock, Nuc.Phys. B126 (1977) 266-284
C     B.R.Martin and M.K.Pidcock, Nuc.Phys. B126 (1977) 285-297
C        Kbar-N Interactions in the Resonance Region. (I) Analysis of data
C             and (II) Amplitudes
      subroutine nucknkn1(enr,tab)
      implicit none
C**** The energy to get the coefficients for.
      real enr
C**** The return variable for the coefficients.
      complex tab(8)
C**** The number of energy bins.
      integer nenr
      parameter (nenr=23)
C**** The energy bins and the coefficients.
      real energy(nenr)
      complex coeffs(8,nenr)
      data energy /
     $     1.550, 1.570, 1.590, 1.610, 1.630,
     $     1.650, 1.670, 1.690, 1.710, 1.730,
     $     1.750, 1.770, 1.790, 1.810, 1.830,
     $     1.850, 1.870, 1.890, 1.910, 1.930,
     $     1.950, 1.970, 1.990/
      data coeffs /
     $     ( 0.150, 0.101),( 0.079, 0.206),( 0.016, 0.022),(-0.074, 0.025),
     $     ( 0.054, 0.007),( 0.008, 0.001),( 0.017, 0.001),( 0.006, 0.000),
     $     ( 0.171, 0.125),( 0.041, 0.230),( 0.016, 0.024),(-0.071, 0.025),
     $     ( 0.063, 0.010),( 0.009, 0.001),( 0.019, 0.001),( 0.007, 0.000),
     $     ( 0.193, 0.152),(-0.001, 0.237),( 0.017, 0.027),(-0.065, 0.024),
     $     ( 0.074, 0.013),( 0.010, 0.001),( 0.021, 0.002),( 0.007, 0.001),
     $     ( 0.213, 0.184),(-0.037, 0.227),( 0.019, 0.029),(-0.055, 0.023),
     $     ( 0.086, 0.018),( 0.011, 0.002),( 0.023, 0.002),( 0.008, 0.001),
     $     ( 0.231, 0.221),(-0.061, 0.206),( 0.020, 0.032),(-0.038, 0.022),
     $     ( 0.101, 0.026),( 0.012, 0.002),( 0.026, 0.003),( 0.009, 0.001),
     $     ( 0.245, 0.264),(-0.071, 0.184),( 0.021, 0.035),(-0.009, 0.035),
     $     ( 0.118, 0.037),( 0.014, 0.003),( 0.029, 0.003),( 0.010, 0.001),
     $     ( 0.254, 0.312),(-0.071, 0.165),( 0.022, 0.038),(-0.008, 0.086),
     $     ( 0.139, 0.053),( 0.016, 0.004),( 0.031, 0.004),( 0.012, 0.001),
     $     ( 0.255, 0.366),(-0.063, 0.151),( 0.023, 0.041),(-0.037, 0.087),
     $     ( 0.163, 0.080),( 0.017, 0.004),( 0.034, 0.005),( 0.013, 0.002),
     $     ( 0.247, 0.424),(-0.050, 0.143),( 0.023, 0.044),(-0.039, 0.074),
     $     ( 0.188, 0.123),( 0.019, 0.006),( 0.038, 0.006),( 0.014, 0.002),
     $     ( 0.226, 0.485),(-0.034, 0.144),( 0.023, 0.047),(-0.030, 0.067),
     $     ( 0.201, 0.194),( 0.022, 0.007),( 0.042, 0.008),( 0.016, 0.003),
     $     ( 0.181, 0.545),(-0.018, 0.153),( 0.023, 0.050),(-0.017, 0.065),
     $     ( 0.171, 0.293),( 0.024, 0.009),( 0.046, 0.009),( 0.018, 0.003),
     $     ( 0.126, 0.582),(-0.004, 0.171),( 0.023, 0.052),(-0.001, 0.068),
     $     ( 0.063, 0.370),( 0.027, 0.012),( 0.050, 0.011),( 0.020, 0.004),
     $     ( 0.067, 0.603),( 0.003, 0.198),( 0.023, 0.055),( 0.015, 0.075),
     $     (-0.072, 0.351),( 0.030, 0.015),( 0.055, 0.014),( 0.023, 0.006),
     $     ( 0.016, 0.606),(-0.002, 0.232),( 0.024, 0.059),( 0.030, 0.087),
     $     (-0.145, 0.269),( 0.033, 0.020),( 0.061, 0.017),( 0.026, 0.007),
     $     (-0.023, 0.607),(-0.022, 0.267),( 0.023, 0.063),( 0.043, 0.106),
     $     (-0.159, 0.194),( 0.036, 0.027),( 0.067, 0.021),( 0.029, 0.010),
     $     (-0.061, 0.610),(-0.057, 0.293),( 0.022, 0.068),( 0.047, 0.130),
     $     (-0.148, 0.142),( 0.038, 0.036),( 0.074, 0.026),( 0.033, 0.013),
     $     (-0.101, 0.608),(-0.101, 0.304),( 0.020, 0.072),( 0.039, 0.152),
     $     (-0.130, 0.110),( 0.037, 0.047),( 0.082, 0.033),( 0.038, 0.019),
     $     (-0.139, 0.600),(-0.143, 0.298),( 0.017, 0.075),( 0.021, 0.164),
     $     (-0.110, 0.089),( 0.031, 0.060),( 0.090, 0.042),( 0.043, 0.036),
     $     (-0.174, 0.587),(-0.178, 0.281),( 0.014, 0.078),( 0.002, 0.161),
     $     (-0.090, 0.078),( 0.019, 0.071),( 0.099, 0.054),( 0.048, 0.039),
     $     (-0.204, 0.570),(-0.203, 0.257),( 0.011, 0.081),(-0.006, 0.147),
     $     (-0.070, 0.073),( 0.003, 0.076),( 0.108, 0.070),( 0.048, 0.056),
     $     (-0.230, 0.551),(-0.219, 0.232),( 0.008, 0.083),(-0.003, 0.132),
     $     (-0.050, 0.075),(-0.014, 0.074),( 0.115, 0.091),( 0.038, 0.079),
     $     (-0.251, 0.531),(-0.228, 0.209),( 0.005, 0.086),( 0.008, 0.120),
     $     (-0.031, 0.086),(-0.027, 0.065),( 0.118, 0.119),( 0.013, 0.096),
     $     (-0.268, 0.511),(-0.232, 0.188),( 0.001, 0.087),( 0.023, 0.113),
     $     (-0.015, 0.109),(-0.034, 0.055),( 0.112, 0.154),(-0.019, 0.094) /

      call nucclxint(enr,tab,energy,coeffs,nenr)

      return
      end

C**** Return the complex legendre coefficients for center of mass energy enr.
C     Hyslop, et al. Phys. Rev. D, Vol. 46, No. 3, 1 August 1992
C     p.961
C     Partial-wave analysis of K+-nucleon scattering
C
      subroutine nucknknp0(enr,tab)
      implicit none
C**** The energy to get the coefficients for.
      real enr
C**** The return variable for the coefficients.
      complex tab(8)
C**** The number of energy bins.
      integer nenr
      parameter (nenr=22)
C**** The energy bins and the coefficients.
      real energy(nenr)
      complex coeffs(8,nenr)
      data energy /
     $     1.466375,1.498051,1.529070,
     $     1.559473,1.589294,1.618566,
     $     1.647318,1.675576,1.703366,
     $     1.730710,1.757628,1.784140,
     $     1.810264,1.836016,1.861412,
     $     1.886467,1.911192,1.935602,
     $     1.959708,1.983521,2.007051,
     $     2.030309 /
      data coeffs /
     $     (-0.054520, 0.002981), (0.039056, 0.001528),
     $     (-0.034182, 0.001170), (0.016752, 0.000281),
     $     (0.006458, 0.000042), (-0.000175, 0.000000),
     $     (-0.000349, 0.000000), (0.000175, 0.000000),
     $     (-0.120622, 0.014768), (0.120791, 0.014810),
     $     (-0.069068, 0.004793), (0.055561, 0.003097),
     $     (0.019368, 0.000375), (-0.000524, 0.000000),
     $     (-0.001920, 0.000004), (0.001396, 0.000002),
     $     (-0.178043, 0.032773), (0.225594, 0.053786),
     $     (-0.095404, 0.009186), (0.091289, 0.008404),
     $     (0.028956, 0.000839), (-0.000698, 0.000000),
     $     (-0.004712, 0.000022), (0.004189, 0.000018),
     $     (-0.223879, 0.052923), (0.330263, 0.124598),
     $     (-0.113835, 0.013131), (0.113326, 0.013012),
     $     (0.032092, 0.001031), (-0.000175, 0.000000),
     $     (-0.008203, 0.000067), (0.007330, 0.000054),
     $     (-0.258715, 0.072137), (0.414519, 0.220403),
     $     (-0.126541, 0.016278), (0.121746, 0.015307),
     $     (0.029653, 0.000880), (0.001396, 0.000002),
     $     (-0.012216, 0.000149), (0.010297, 0.000106),
     $     (-0.283778, 0.088333), (0.468111, 0.325724),
     $     (-0.134460, 0.018419), (0.121057, 0.015908),
     $     (0.022856, 0.000523), (0.003840, 0.000015),
     $     (-0.016054, 0.000258), (0.012391, 0.000154),
     $     (-0.300768, 0.100577), (0.493321, 0.426624),
     $     (-0.139160, 0.019756), (0.116973, 0.017229),
     $     (0.013263, 0.000176), (0.007505, 0.000056),
     $     (-0.019891, 0.000396), (0.013437, 0.000181),
     $     (-0.310847, 0.108370), (0.497259, 0.515279),
     $     (-0.141003, 0.020294), (0.115293, 0.021244),
     $     (0.002094, 0.000004), (0.012216, 0.000149),
     $     (-0.023728, 0.000563), (0.013082, 0.000421),
     $     (-0.315067, 0.111757), (0.487246, 0.589953),
     $     (-0.140836, 0.020245), (0.120756, 0.031046),
     $     (-0.009599, 0.000092), (0.018138, 0.000579),
     $     (-0.027214, 0.000741), (0.012036, 0.000395),
     $     (-0.314389, 0.111207), (0.468636, 0.652087),
     $     (-0.139160, 0.019756), (0.136501, 0.049592),
     $     (-0.020938, 0.000439), (0.025271, 0.001140),
     $     (-0.030698, 0.000943), (0.010112, 0.000603),
     $     (-0.309752, 0.107503), (0.443451, 0.703964),
     $     (-0.136308, 0.018939), (0.162538, 0.080082),
     $     (-0.031744, 0.001009), (0.033783, 0.001895),
     $     (-0.033833, 0.001146), (0.007838, 0.001063),
     $     (-0.301743, 0.101313), (0.412335, 0.747951),
     $     (-0.132610, 0.017906), (0.194147, 0.126092),
     $     (-0.041317, 0.001710), (0.043794, 0.003430),
     $     (-0.036793, 0.001356), (0.005739, 0.001786),
     $     (-0.291203, 0.093551), (0.372290, 0.784430),
     $     (-0.128229, 0.016722), (0.218925, 0.191258),
     $     (-0.049486, 0.002455), (0.055282, 0.005587),
     $     (-0.039577, 0.001569), (0.003992, 0.002774),
     $     (-0.279017, 0.085091), (0.319790, 0.808601),
     $     (-0.123499, 0.015492), (0.210199, 0.263745),
     $     (-0.056081, 0.003155), (0.068786, 0.009318),
     $     (-0.042361, 0.001798), (0.003459, 0.004533),
     $     (-0.266142, 0.076717), (0.261038, 0.807586),
     $     (-0.118589, 0.014267), (0.176578, 0.298580),
     $     (-0.060588, 0.003684), (0.083689, 0.015494),
     $     (-0.044795, 0.002011), (0.004300, 0.007322),
     $     (-0.253769, 0.069185), (0.214924, 0.782945),
     $     (-0.113496, 0.013052), (0.153387, 0.305778),
     $     (-0.063186, 0.004009), (0.098793, 0.024616),
     $     (-0.047228, 0.002235), (0.006811, 0.012195),
     $     (-0.242862, 0.062944), (0.187700, 0.750537),
     $     (-0.108390, 0.011890), (0.139260, 0.302879),
     $     (-0.063705, 0.004075), (0.113475, 0.035916),
     $     (-0.049312, 0.002438), (0.010569, 0.019492),
     $     (-0.234581, 0.058444), (0.174278, 0.720515),
     $     (-0.103444, 0.010818), (0.129089, 0.294851),
     $     (-0.061974, 0.003856), (0.127420, 0.049707),
     $     (-0.051396, 0.002649), (0.015440, 0.029615),
     $     (-0.229945, 0.056012), (0.168024, 0.695622),
     $     (-0.098486, 0.009796), (0.119268, 0.284303),
     $     (-0.057815, 0.003354), (0.140125, 0.065126),
     $     (-0.053479, 0.002868), (0.020776, 0.042487),
     $     (-0.229945, 0.056012), (0.164861, 0.675559),
     $     (-0.094033, 0.008922), (0.108832, 0.271624),
     $     (-0.051223, 0.002631), (0.151560, 0.082232),
     $     (-0.055214, 0.003058), (0.026239, 0.058342),
     $     (-0.235506, 0.058936), (0.162869, 0.659604),
     $     (-0.089744, 0.008120), (0.096756, 0.257847),
     $     (-0.042187, 0.001783), (0.161451, 0.100396),
     $     (-0.057122, 0.003274), (0.031198, 0.077180),
     $     (-0.247122, 0.065339), (0.160723, 0.646349),
     $     (-0.085964, 0.007445), (0.082437, 0.244160),
     $     (-0.030524, 0.000933), (0.169685, 0.119596),
     $     (-0.058682, 0.003456), (0.035554, 0.098458) /

      call nucclxint(enr,tab,energy,coeffs,nenr)

      return
      end

C**** Return the complex legendre coefficients for center of mass energy enr.
C     Hyslop, et al. Phys. Rev. D, Vol. 46, No. 3, 1 August 1992
C     p.961
C     Partial-wave analysis of K+-nucleon scattering
C
      subroutine nucknknp1(enr,tab)
      implicit none
C**** The energy to get the coefficients for.
      real enr
C**** The return variable for the coefficients.
      complex tab(8)
C**** The number of energy bins.
      integer nenr
      parameter (nenr=22)
C**** The energy bins and the coefficients.
      real energy(nenr)
      complex coeffs(8,nenr)
      data energy /
     $     1.466375,1.498051,1.529070,
     $     1.559473,1.589294,1.618566,
     $     1.647318,1.675576,1.703366,
     $     1.730710,1.757628,1.784140,
     $     1.810264,1.836016,1.861412,
     $     1.886467,1.911192,1.935602,
     $     1.959708,1.983521,2.007051,
     $     2.030309 /
      data coeffs /
     $     (-0.224659, 0.053314), (-0.044447, 0.001979),
     $     (0.021112, 0.000446), (-0.000175, 0.000000),
     $     (-0.001222, 0.000001), (-0.000175, 0.000000),
     $     (0.000000, 0.000000), (0.000000, 0.000000),
     $     (-0.299512, 0.099634), (-0.094719, 0.009054),
     $     (0.048097, 0.002319), (-0.000873, 0.000001),
     $     (-0.004363, 0.000019), (-0.001571, 0.000002),
     $     (0.000524, 0.000000), (0.000349, 0.000000),
     $     (-0.348457, 0.141423), (-0.137651, 0.019321),
     $     (0.073559, 0.005441), (-0.002618, 0.000007),
     $     (-0.008203, 0.000067), (-0.004014, 0.000016),
     $     (0.001222, 0.000001), (0.001222, 0.000001),
     $     (-0.384253, 0.180079), (-0.172813, 0.030814),
     $     (0.096261, 0.009354), (-0.005236, 0.000027),
     $     (-0.011867, 0.000141), (-0.006981, 0.000049),
     $     (0.002269, 0.000005), (0.002269, 0.000005),
     $     (-0.411766, 0.216365), (-0.201293, 0.042309),
     $     (0.116553, 0.013774), (-0.008726, 0.000076),
     $     (-0.015182, 0.000231), (-0.010122, 0.000102),
     $     (0.003491, 0.000012), (0.003316, 0.000011),
     $     (-0.433622, 0.251058), (-0.223543, 0.053594),
     $     (0.135232, 0.018895), (-0.012908, 0.000417),
     $     (-0.017799, 0.000317), (-0.013437, 0.000181),
     $     (0.004887, 0.000024), (0.004189, 0.000018),
     $     (-0.450991, 0.284114), (-0.239911, 0.063603),
     $     (0.153358, 0.024625), (-0.017607, 0.000811),
     $     (-0.019717, 0.000389), (-0.016578, 0.000275),
     $     (0.006283, 0.000039), (0.005061, 0.000026),
     $     (-0.464888, 0.315937), (-0.250048, 0.073121),
     $     (0.171477, 0.031923), (-0.022613, 0.002015),
     $     (-0.021112, 0.000446), (-0.019359, 0.000625),
     $     (0.007854, 0.000062), (0.005585, 0.000031),
     $     (-0.475666, 0.346730), (-0.253606, 0.083546),
     $     (0.190347, 0.042269), (-0.027874, 0.004548),
     $     (-0.021810, 0.000476), (-0.021962, 0.000983),
     $     (0.009599, 0.000092), (0.005757, 0.000283),
     $     (-0.483694, 0.376348), (-0.251362, 0.097424),
     $     (0.209651, 0.058077), (-0.034215, 0.008239),
     $     (-0.021635, 0.000468), (-0.024178, 0.002089),
     $     (0.011169, 0.000125), (0.005751, 0.000784),
     $     (-0.489159, 0.405271), (-0.247176, 0.117045),
     $     (0.227089, 0.081753), (-0.041780, 0.012656),
     $     (-0.020927, 0.000688), (-0.026370, 0.003457),
     $     (0.012914, 0.000167), (0.005223, 0.001279),
     $     (-0.492260, 0.433443), (-0.245186, 0.140300),
     $     (0.236841, 0.114116), (-0.050389, 0.017340),
     $     (-0.019697, 0.000889), (-0.028694, 0.005352),
     $     (0.014477, 0.000460), (0.004348, 0.001772),
     $     (-0.493183, 0.461012), (-0.246241, 0.164046),
     $     (0.233876, 0.150640), (-0.059470, 0.022545),
     $     (-0.017598, 0.001061), (-0.031000, 0.007264),
     $     (0.016221, 0.000513), (0.003126, 0.002516),
     $     (-0.492034, 0.487631), (-0.249479, 0.186137),
     $     (0.216980, 0.183109), (-0.068980, 0.028310),
     $     (-0.014970, 0.001476), (-0.033442, 0.009458),
     $     (0.017955, 0.000823), (0.001388, 0.003011),
     $     (-0.488946, 0.513486), (-0.253387, 0.206240),
     $     (0.191811, 0.205537), (-0.078709, 0.034645),
     $     (-0.011646, 0.002140), (-0.035842, 0.011927),
     $     (0.019697, 0.000889), (-0.000520, 0.003764),
     $     (-0.483517, 0.538223), (-0.257449, 0.223920),
     $     (0.165425, 0.217362), (-0.088610, 0.041581),
     $     (-0.007806, 0.003070), (-0.038386, 0.014432),
     $     (0.021429, 0.001211), (-0.002767, 0.004528),
     $     (-0.476127, 0.561668), (-0.260892, 0.239836),
     $     (0.142134, 0.221526), (-0.098632, 0.049143),
     $     (-0.003459, 0.004533), (-0.040711, 0.017203),
     $     (0.023158, 0.001539), (-0.005010, 0.005051),
     $     (-0.466385, 0.583579), (-0.264032, 0.253958),
     $     (0.122761, 0.222188), (-0.108413, 0.057285),
     $     (0.001205, 0.006798), (-0.043175, 0.020015),
     $     (0.024885, 0.001873), (-0.007590, 0.005842),
     $     (-0.454410, 0.603739), (-0.266203, 0.266792),
     $     (0.107211, 0.220705), (-0.118215, 0.066093),
     $     (0.006337, 0.009378), (-0.045583, 0.023112),
     $     (0.026611, 0.002213), (-0.010334, 0.006651),
     $     (-0.440709, 0.621556), (-0.267669, 0.278407),
     $     (0.094484, 0.219246), (-0.127678, 0.075502),
     $     (0.011397, 0.012794), (-0.047929, 0.026499),
     $     (0.028160, 0.002549), (-0.013234, 0.007735),
     $     (-0.424911, 0.636750), (-0.268664, 0.288885),
     $     (0.083948, 0.217683),(-0.136901, 0.085563),
     $     (0.016527, 0.017064),(-0.050237, 0.029919),
     $     (0.029866, 0.003152),(-0.016293, 0.008844),
     $     (-0.407827, 0.648920),(-0.268958, 0.298353),
     $     (0.074901, 0.216910),(-0.145683, 0.096235),
     $     (0.021207, 0.021931),(-0.052639, 0.033653),
     $     (0.031411, 0.003499),(-0.019339, 0.009973) /

      call nucclxint(enr,tab,energy,coeffs,nenr)

      return
      end



      subroutine printinttype
      common /nucint/uin,chint,ein,pin(3),upin(3)
#include "partnucint.h"
      real chint
      if(pninttype.eq.1) then
         write (6,'('' % interaction '',i8,'' absorbtion'')')pninttype
      else if(pninttype.eq.2) then
         write (6,'('' % interaction '',i8,'' production'')')pninttype
      else if(pninttype.eq.3) then
         write (6,'('' % interaction '',i8,'' scatter'')')pninttype
      else if(pninttype.eq.4) then
         write (6,'('' % interaction '',i8,'' exchange'')')pninttype
      end if
C      write(6,'('' %% c,m,E '',f4.0,2f10.1)')chint,uin,ein

      return
      end

