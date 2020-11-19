      subroutine nrprd(no,e,ka,piso,chiso)
      implicit none
      real*4 piso(4,4),chiso(4) !momentum and charge of products
c     C     interaction type, total energy, charge, momentum, charge

c     ka= 0   primary neutron
c     1   primary proton
c     no=     1       elastic
c     2       single pi production
c     3       double pi production
c
      integer no,ka
      real*4 e                  !proton energy in resting target

#include "nrcms.h"
C     common /nuccms/rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
C     real*4 rhon,pfermi,unucl,ecms2
C     real*4 up4cms,beta
C     logical prot

      real*4 frinn(161),fmxsm(161),dmin(101),fmxdn(130),fsln(176)

C     the data as in MECC7 and GCALOR  (BERTINI's)

      data frinn/
     1     2.,   7.,  16.,  30.,  52.,  86., 134., 201., 290., 406.,
     2     551., 727., 945.,1204.,1502.,1832.,2179.,2510.,2806.,3075.,
     3     3321.,3546.,3751.,3935.,4100.,4246.,4374.,4486.,4585.,4674.,
     4     4755.,4830.,4898.,4962.,5022.,5076.,5127.,5174.,5219.,5260.,
     5     5300.,5337.,5372.,5405.,5437.,5467.,5496.,5523.,5549.,5574.,
     6     5599.,5623.,5646.,5670.,5694.,5718.,5742.,5768.,5794.,5820.,
     7     5848.,5877.,5906.,5936.,5967.,5999.,6032.,6065.,6100.,6134.,
     8     6170.,6206.,6242.,6279.,6317.,6354.,6392.,6430.,6469.,6507.,
     9     6545.,6583.,6622.,6660.,6699.,6737.,6777.,6816.,6857.,6898.,
     1     6940.,6983.,7026.,7071.,7116.,7162.,7209.,7258.,7307.,7358.,
     2     7409.,7462.,7516.,7570.,7624.,7679.,7734.,7790.,7845.,7900.,
     3     7954.,8008.,8062.,8115.,8167.,8218.,8268.,8318.,8366.,8414.,
     4     8461.,8507.,8552.,8597.,8641.,8684.,8727.,8769.,8811.,8852.,
     5     8892.,8932.,8972.,9011.,9049.,9088.,9126.,9164.,9202.,9239.,
     6     9277.,9314.,9351.,9388.,9424.,9461.,9497.,9534.,9570.,9606.,
     7     9642.,9678.,9714.,9750.,9786.,9821.,9857.,9893.,9929.,9964.,
     8     10000./

      data fmxsm/
     1     341., 598., 776., 923.,1052.,1168.,1276.,1377.,1472.,1563.,
     2     1651.,1735.,1816.,1896.,1973.,2048.,2122.,2194.,2265.,2335.,
     3     2404.,2471.,2538.,2604.,2669.,2733.,2797.,2860.,2922.,2984.,
     4     3045.,3106.,3167.,3227.,3286.,3345.,3404.,3462.,3520.,3578.,
     5     3635.,3692.,3749.,3806.,3862.,3918.,3974.,4030.,4085.,4140.,
     6     4195.,4250.,4305.,4359.,4413.,4468.,4522.,4575.,4629.,4683.,
     7     4736.,4789.,4842.,4895.,4948.,5001.,5054.,5106.,5158.,5211.,
     8     5263.,5315.,5367.,5419.,5471.,5523.,5574.,5626.,5678.,5729.,
     9     5780.,5832.,5883.,5934.,5985.,6036.,6087.,6138.,6189.,6240.,
     1     6290.,6341.,6392.,6442.,6493.,6543.,6594.,6644.,6694.,6745.,
     2     6795.,6845.,6895.,6945.,6995.,7045.,7095.,7145.,7195.,7245.,
     3     7295.,7344.,7394.,7444.,7493.,7543.,7593.,7642.,7692.,7741.,
     4     7791.,7840.,7889.,7939.,7988.,8038.,8087.,8136.,8185.,8235.,
     5     8284.,8333.,8382.,8431.,8480.,8529.,8578.,8627.,8676.,8725.,
     6     8774.,8823.,8872.,8921.,8970.,9019.,9068.,9116.,9165.,9214.,
     7     9263.,9312.,9360.,9409.,9458.,9506.,9555.,9604.,9652.,9701.,
     8     9749./

      data dmin/
     1     1078.6,1133.,1148.,1158.,1166.,1172.,1178.,1183.,1187.,1191.,
     2     1195.,1198.,1202.,1205.,1208.,1211.,1213.,1216.,1219.,1221.,
     3     1224.,1226.,1229.,1231.,1234.,1236.,1239.,1242.,1245.,1248.,
     4     1251.,1255.,1258.,1262.,1265.,1269.,1273.,1277.,1282.,1286.,
     5     1291.,1296.,1302.,1308.,1315.,1323.,1331.,1341.,1351.,1363.,
     6     1376.,1391.,1408.,1427.,1449.,1475.,1505.,1538.,1568.,1594.,
     7     1618.,1640.,1660.,1680.,1699.,1718.,1736.,1754.,1772.,1789.,
     8     1805.,1820.,1835.,1849.,1862.,1875.,1887.,1899.,1911.,1923.,
     9     1935.,1947.,1960.,1973.,1986.,2000.,2014.,2029.,2044.,2059.,
     1     2075.,2091.,2107.,2124.,2140.,2157.,2174.,2190.,2207.,2224.,
     2     2240./

      data fmxdn/
     1     2221.,2297.,2373.,2446.,2519.,2590.,2660.,2730.,2798.,2865.,
     2     2932.,2998.,3063.,3128.,3192.,3255.,3318.,3380.,3442.,3504.,
     3     3565.,3625.,3685.,3745.,3804.,3863.,3922.,3981.,4039.,4096.,
     4     4154.,4211.,4268.,4325.,4382.,4438.,4494.,4550.,4606.,4661.,
     5     4717.,4772.,4827.,4882.,4936.,4991.,5045.,5100.,5154.,5208.,
     6     5261.,5315.,5369.,5422.,5475.,5529.,5582.,5635.,5687.,5740.,
     7     5793.,5845.,5898.,5950.,6003.,6055.,6107.,6159.,6211.,6263.,
     8     6315.,6366.,6418.,6470.,6521.,6573.,6624.,6675.,6727.,6778.,
     9     6829.,6880.,6931.,6982.,7033.,7084.,7135.,7185.,7236.,7287.,
     1     7337.,7388.,7438.,7489.,7539.,7590.,7640.,7690.,7740.,7791.,
     2     7841.,7891.,7941.,7991.,8041.,8091.,8141.,8191.,8241.,8291.,
     3     8341.,8390.,8440.,8490.,8539.,8589.,8639.,8688.,8738.,8787.,
     4     8837.,8886.,8936.,8985.,9035.,9084.,9133.,9183.,9232.,9281./

      data fsln/
     1     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,
     2     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,
     3     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,
     4     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,
     5     5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5050.,5047.,
     6     5047.,5043.,5040.,5038.,5033.,5031.,5030.,5030.,5035.,5048.,
     7     5061.,5106.,5162.,5237.,5285.,5323.,5356.,5379.,5403.,5428.,
     8     5443.,5463.,5483.,5497.,5515.,5531.,5550.,5565.,5582.,5602.,
     9     5620.,5636.,5658.,5672.,5688.,5701.,5718.,5733.,5747.,5761.,
     1     5771.,5794.,5807.,5821.,5836.,5851.,5866.,5879.,5893.,5905.,
     2     5919.,5933.,5946.,5959.,5968.,5978.,5992.,6005.,6015.,6026.,
     3     6035.,6051.,6061.,6075.,6087.,6096.,6109.,6118.,6126.,6136.,
     4     6149.,6160.,6168.,6181.,6191.,6199.,6213.,6221.,6230.,6241.,
     5     6249.,6258.,6265.,6274.,6281.,6290.,6300.,6306.,6318.,6325.,
     6     6335.,6339.,6348.,6357.,6362.,6374.,6381.,6387.,6394.,6389.,
     7     6416.,6417.,6422.,6427.,6431.,6436.,6439.,6445.,6451.,6455.,
     8     6462.,6462.,6467.,6473.,6476.,6480.,6485.,6489.,6495.,6499.,
     9     6500.,6507.,6510.,6515.,6518.,6525./

      real*4 dldir(3),pdl(4),pdldec(4,4)
      real*4 upi,amnucl,am2n2
      parameter (upi=139.6, amnucl=939., am2n2=4*amnucl**2)
      real*4 pcmsm,t,tt,w
      integer iraz,ka1,ka2,ka3,ka4,i
      real*4 rf,rf1,ranf,erem,amd1,phsm,ed,pd2,pdlt,phs,aiso,cdl
      real*4 rpi,rnu,amd2,rf2,erem1,rat
      integer ilo,ihi,ilo1,ihi1
      real*4 ekin,nrcst
      external nrcst


      if(no.eq.1) then          !elastic
         pcmsm=sqrt((ecms2-(amnucl-unucl)**2)*
     &        (ecms2-(amnucl+unucl)**2)/(4*ecms2)) ! cms momentum
         ekin=e-amnucl
         call nrranve(pcmsm,pdldec(1,1),nrcst(ekin,ka,prot),-1.)
         do i=1,3
            pdldec(i,3)=-pdldec(i,1)
         end do
         pdldec(4,1)=sqrt(pcmsm**2+amnucl**2)
c     nrllor performs a lorentz shift from first argument is prime of third
c     second argument is prime of fourth
         call nrllor(piso(1,1),piso(4,1),pdldec(1,1),pdldec(4,1),beta)
         call nrllor(piso(1,3),piso(4,3),pdldec(1,3),pdldec(4,1),beta)
         chiso(1)=ka
         if(prot) then
            chiso(3)=1
         else
            chiso(3)=0
         end if
         return
      end if

C     all below is Lindenbaum's Isobar model for pion production
C     Phys Rev 105,1874(1957) taken with some fixup from MECC7

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

      if(no.eq.2)then
         t=tt
      else if(no.eq.3) then
         t=((w-upi)**2-am2n2)/2/amnucl
      end if
C     call nrhis(t-300.,20.,161,ilo,ihi,erem)
 1    call nrhis(t-300.,20.,161,ilo,ihi,erem)
      rf1=(frinn(ihi)*erem+frinn(ilo)*(1-erem))*1e-4
      phsm=(fmxsm(ihi)*erem+fmxsm(ilo)*(1-erem))*100.

 2    rf=rf1*ranf()
      call nrhis(rf*100,1.,101,ilo1,ihi1,erem1)
      amd1=dmin(ihi1)*erem1+dmin(ilo1)*(1-erem1)
      if(amd1.lt.dmin(1)) amd1=dmin(1)
      if(no.eq.2)then
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
         call nrcone(up4cms,dldir,cdl,-1.)
         do i=1,3
            pdl(i)=pdlt*dldir(i)
            pdldec(i,3)=-pdl(i)
         end do
         pdl(4)=ed
         pdldec(4,3)=w-ed
         call nr2iso(amd1,amnucl,upi,pdl,pdldec)
         do i=1,3
            call nrllor(piso(1,i),piso(4,i),
     $           pdldec(1,i),pdldec(4,i),beta)
         end do
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
            call nrhis(tt-920.,20.,130,ilo1,ihi1,erem1)
            phsm=(fmxdn(ihi1)*erem1+fmxdn(ilo1)*(1-erem1))*100.
            if(ranf().gt.phs*rf1/(phsm*rf2))then
               iraz=0
               rf1=rf2
               go to 2
            end if

            cdl=-.9999995
            if(ranf().lt.0.5)cdl=-cdl
            call nrcone(up4cms,dldir,cdl,-1.)
            do i=1,3
               pdl(i)=pdlt*dldir(i)
            end do
            pdl(4)=ed
            call nr2iso(amd1,amnucl,upi,pdl,pdldec)
            do i=1,3
               pdl(i)=-pdl(i)
            end do
            pdl(4)=w-ed
            call nr2iso(amd2,amnucl,upi,pdl,pdldec(1,3))
            do i=1,4
               call nrllor(piso(1,i),piso(4,i),
     $              pdldec(1,i),pdldec(4,i),beta)
            end do
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
               call nrhis(tt,20.,176,ilo,ihi,erem)
               rat=(fsln(ihi)*erem+fsln(ilo)*(1-erem))*1e-4
               if(rpi.lt.rat)then
                  chiso(1)=ka1
                  chiso(2)=ka3
                  chiso(3)=ka2
                  chiso(4)=ka4
C               else if(rpi.lt.(1.-rat)*0.6667)then
C     - Bug fix (R. Tacik July 11, 2011)
               else if (rpi.lt. (rat+(1.-rat)*0.6667) ) then
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
