      subroutine nrnuc(idec)
      implicit none
#include "nrcard.h"
C     COMMON/NUCRES/NUCRESCAT,XNUCFACT,XNUCELAFACT,XNUCSPIFACT,XNUCDPIFACT,NUCRESFLG
#include "nrint.h"
C     common/nucint/uin,chint,ein,pin(3),upin(3),ptotnuc
#include "nrcms.h"
C     common/nuccms/rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
C     logical prot
#include "nrnewprd.h"        
C     common/nucnewprd/piso(4,4),chiso(4)
#include "nrhad.h"
C     common/nuchad/ihap,nel(3)
#include "nucleonfsihist.h"

c     A single trip to this subroutine corresponds to a single 
c     step (not vertex) in a nucleon cascade (regen or reweight)

      real*4    xspp(158),xspn(158),xdpp(130),xdpn(130)
      real*4    xepp(176),xepn(176)
C     all the data the same as in MECC7 or GCALOR (BERTINI's)
c     bunch of data tables from which we will throw to calculate interaction probabilities
c     p-p elastic
      data xepp/
     $     67500.,15500.,6750.,4420.,3230.,2800.,2520.,2370.,2300.,2275.,
     $     2260.,2260.,2260.,2260.,2270.,2280.,2295.,2300.,2310.,2330.,
     $     2350.,2380.,2395.,2420.,2460.,2485.,2500.,2530.,2565.,2600.,
     $     2620.,2640.,2660.,2675.,2690.,2700.,2705.,2710.,2715.,2720.,
     $     2725.,2725.,2720.,2715.,2710.,2700.,2695.,2680.,2670.,2660.,
     $     2640.,2625.,2605.,2590.,2570.,2545.,2525.,2500.,2480.,2470.,
     $     2450.,2430.,2410.,2395.,2370.,2360.,2340.,2325.,2305.,2290.,
     $     2275.,2270.,2260.,2250.,2230.,2225.,2210.,2200.,2195.,2190.,
     $     2175.,2165.,2150.,2140.,2125.,2120.,2105.,2100.,2090.,2075.,
     $     2065.,2055.,2045.,2030.,2020.,2005.,2000.,1995.,1980.,1975.,
     $     1965.,1960.,1945.,1935.,1925.,1910.,1900.,1895.,1885.,1875.,
     $     1865.,1860.,1850.,1830.,1825.,1810.,1800.,1795.,1785.,1775.,
     $     1765.,1755.,1745.,1730.,1720.,1715.,1700.,1695.,1680.,1675.,
     $     1665.,1650.,1645.,1630.,1615.,1605.,1600.,1590.,1580.,1565.,
     $     1555.,1550.,1530.,1520.,1510.,1500.,1490.,1475.,1465.,1455.,
     $     1445.,1440.,1425.,1415.,1405.,1398.,1390.,1380.,1370.,1360.,
     $     1350.,1335.,1325.,1315.,1300.,1295.,1280.,1270.,1260.,1250.,
     $     1240.,1225.,1210.,1200.,1195.,1175./


c     p-n elastic
      data xepn/
     $     196500.,47500.,22000.,13000.,9180.,7300.,6030.,5180.,4680.,
     $     4320.,4080.,3910.,3760.,3650.,3550.,3480.,3415.,3370.,3325.,
     $     3290.,3275.,3250.,3255.,3275.,3285.,3275.,3220.,3150.,3075.,
     $     2990.,2875.,2775.,2695.,2630.,2590.,2565.,2560.,2560.,2560.,
     $     2565.,2570.,2575.,2578.,2580.,2585.,2580.,2575.,2560.,2540.,
     $     2505.,2470.,2425.,2375.,2315.,2275.,2230.,2200.,2175.,2155.,
     $     2145.,2130.,2125.,2115.,2105.,2100.,2095.,2090.,2080.,2070.,
     $     2060.,2050.,2045.,2040.,2030.,2025.,2020.,2015.,2010.,2005.,
     $     2002.,2000.,1999.,1990.,1985.,1975.,1970.,1965.,1960.,1950.,
     $     1945.,1940.,1925.,1920.,1915.,1910.,1900.,1898.,1895.,1890.,
     $     1880.,1875.,1870.,1865.,1855.,1850.,1845.,1830.,1825.,1820.,
     $     1810.,1805.,1800.,1795.,1790.,1785.,1775.,1770.,1765.,1750.,
     $     1745.,1740.,1725.,1720.,1715.,1705.,1699.,1690.,1685.,1675.,
     $     1670.,1650.,1640.,1630.,1620.,1610.,1600.,1595.,1580.,1575.,
     $     1560.,1550.,1530.,1520.,1510.,1500.,1490.,1485.,1470.,1460.,
     $     1450.,1440.,1435.,1425.,1415.,1410.,1400.,1390.,1385.,1375.,
     $     1365.,1350.,1340.,1330.,1320.,1310.,1300.,1290.,1280.,1270.,
     $     1260.,1250.,1240.,1230.,1220.,1210.,1200./

c     p-p single pi production
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

c     p-n single pi production
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

c     p-p double pi production
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

c     p-n double pi production
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
      real*4    step,fac,amn,amn2
      parameter (step=0.2, fac=1e-2, amn=939., amn2=amn**2)
      
      real*4    ampi
      parameter (ampi=139.6)
      
      integer*4 idec
      
      real*4    eineq,ekin,psp,pdp,erem1,erem2,ptot,xr,yr,dist,p1,p2
      real*4    erem,pel
      
      integer*4 ka,ichint,ilo,ihi,ilo1,ihi1,ilo2,ihi2,iqq
      
      real*4    ranf
      external  ranf
      
      real*4 A,B,C
      common /martindebug/  ptot,A,B,C

c     REWEIGHTING: ratio of pcascprob in regen and reweight gives the weight
c     pcascprob is the product of the probability of each step
c     REWEIGHTING: NFiflagStep gives the interaction type of a step:
c     5 = non event, 4 = survival, 3 = douple pi, 2 = single pi, 1 = elastic



c     not sure what 1 is (I'm guessing it will turn out to be survival)
      idec=1
c     REWEIGHTING: if reweighting, and the interaction type flag is 5, it's a non event
      if(NReweightFlag.eq.1) then

         if( NFiflagStep(RW_step_i).eq.5) then
            

            goto 200
         endif
      end if
      
      


      if(NReweightFlag.eq.0) then
         if(ecms2.le.(2.*amn)**2) then
         

            iqq = 4
            NFiflagStep(NFNStep) = 5;
C            print *,"NFiFlagstep recorded at ",NFNstep

c     treat kinematic failure as survival

            goto 200            ! check if kinematically possible    
         endif
      
         

C      endif
      
c     if not kinematically possible QUIT and return 1

c     calculate total and kinetic energy
      eineq=(ecms2-unucl**2-amn**2)/(2*unucl) ! unucl at rest
      ekin = eineq - amn
      
c     prot is target and ka is its charge
c     chint is for incoming
c     set ka for proton or neutron
      if(prot)then
         ka=1
      else
         ka=0
      end if
c     set charge from chint from nrprton.F
      ichint=chint

C      if(NReweightFlag.eq.0) then
         
c     if the kinetic energy is less than 360
c     then we will not produce pions (why 360?)

         if( ekin  .lt. 360. .or.
     $        ecms2 .lt. (amn+amn+ampi)**2 ) then
            psp = 0.
            pdp = 0.
            call nrhis(ekin,20.,176,ilo,ihi,erem)
c     set up the histogram of kinetic energy
c     ilo, ihi, erem get set within nrhis
c     if proton or neutral? prob of elastic taken from proton table
            if(ka-ichint.eq.0) then
               pel = xepp(ihi)*erem + xepp(ilo)*(1-erem)
c     otherwise otherwise prob of elastic taken from neutront table
            else
               pel = xepn(ihi)*erem + xepn(ilo)*(1-erem)
            end if
c     next kinetic range less than 920 or ecms2 is less than 2*mass of nucleon + 2*mass of pion squared
c     now single pion production, but no two pion production
         else if( ekin .lt. 920. .or. 
     $           ecms2 .lt. (amn+amn+ampi+ampi)**2 ) then   
            pdp = 0.
c     set up two hists
            call nrhis(ekin,20.,176,ilo,ihi,erem)
            call nrhis(ekin-360.,20.,158,ilo1,ihi1,erem1)
c     p-p or n-p? calculate pel and psp
            if(ka-ichint.eq.0)then
               pel = xepp(ihi) *erem  + xepp(ilo) *(1-erem)
               psp = xspp(ihi1)*erem1 + xspp(ilo1)*(1-erem1)
            else
               pel = xepn(ihi) *erem  + xepn(ilo) *(1-erem)
               psp = xspn(ihi1)*erem1 + xspn(ilo1)*(1-erem1)
            end if
         else
c     Now also with double pion production
c     calculate pel, psp, pdp
            call nrhis(ekin,20.,176,ilo,ihi,erem)
            call nrhis(ekin-360.,20.,158,ilo1,ihi1,erem1) 
            call nrhis(ekin-920.,20.,130,ilo2,ihi2,erem2)
            if(ka-ichint.eq.0)then
               pel = xepp(ihi) *erem  + xepp(ilo) *(1-erem)
               psp = xspp(ihi1)*erem1 + xspp(ilo1)*(1-erem1)
               pdp = xdpp(ihi2)*erem2 + xdpp(ilo2)*(1-erem2)
            else
               pel = xepn(ihi) *erem  + xepn(ilo) *(1-erem)
               psp = xspn(ihi1)*erem1 + xspn(ilo1)*(1-erem1)
               pdp = xdpn(ihi2)*erem2 + xdpn(ilo2)*(1-erem2)
            end if
         end if
      end if

c     apply nucleon fsi dials
c     REWEIGHTING: get elastic, single pi, and double pi probs from generated
      if(NReweightFlag.eq.1) then

         pel= nfipel(RW_step_i)
         psp=nfipsp(RW_step_i)
         pdp=nfipdp(RW_step_i)


      endif


c     REWEIGHTING: set elastic, single pi, and double pi probs during generation 
      if(NReweightFlag.eq.0) then

         nfipel(nfnstep) = pel
         nfipsp(nfnstep)=psp
         nfipdp(nfnstep)=pdp
      endif


      

c     REWEIGHTING: reweight or generate apply the dial values to the probabilities
C  Toby Nonnenmacher change april 2018
C      pel=pel*fac*xnucfact
      pel=pel*fac*(xnucfact * xnucelafact) ! xnucfact changed to xnucfact * xnucelafact with default xnucfact of 1
C      psp=psp*fac*xnucfact
      psp=psp*fac*(xnucfact * xnucspifact) ! xnucfact changed to xnucfact * xnucspifact with default xnucfact of 1 
C      pdp=pdp*fac*xnucfact
      pdp=pdp*fac*(xnucfact * xnucdpifact) ! xnucfact changed to xnucfact * xnucdpifact with default xnucfact of 1 

 

      ptot=pel+psp+pdp
c     calculate ptot for normalisation



c     if ptot is 0, have a break down, QUIT, and return idec 1 
      if(ptot.eq.0) then

         goto 200
      endif
c     set random xr between 0 and 1 and y between 0 and ptot
      xr=ranf()+1.e-6
      yr=ranf()*ptot
c     ptot multiplied by rhon (Which is set in nrfermi.F)
c     REWEIGHTING: Store RHON for reweighting purposes
c     This lets us skip nrfermi.F in the reweight
      if(NReweightFlag.eq.0) then
      nfirhon(nfnstep) = rhon
      ENDIF
      if(NReweightFlag.eq.1) then
      rhon = nfirhon(RW_step_i)
      ENDIF
c     apply rhon to the probability
      ptot=ptot*rhon


c     set mean free path type distance according to how cascades are calculated
      dist=-log(xr)/ptot
c     count number of times we go through here
c      pnuccounter = pnuccounter +1;
c     caluclate probability of survival


      psurvnuc = exp(-ptot*step);

      ptotnuc = ptotnuc*ptot
c     ptotnuc multiplied by ptot at every single step


c     normalise pel-,psp-,pdp- nuc to survival probability
      pelnuc = (1-psurvnuc)*(pel/(pel+psp+pdp))
      pspnuc = (1-psurvnuc)*(psp/(pel+psp+pdp))
      pdpnuc = (1-psurvnuc)*(pdp/(pel+psp+pdp))
c     ptotalprob is now normalised and should equal 1
      ptotalprob = psurvnuc + pelnuc + pspnuc + pdpnuc
C     Still need to do everything above this point for reweighting, which will depend on table
C     dial values, and energies: (maybe kinetic, total, centre of mass squared?)


      if(NReweightFlag.eq.0) then


c     if surival;
      If(step.lt.dist) then
c     increment pcascprob by psurvnuc for this step


         pcascprob = pcascprob * psurvnuc
c     count a survival
c         survivalcounter = survivalcounter + 1
         NFiflagStep(NFNStep) = 4;





      end if
c     then QUIT
      if(step.lt.dist) goto 200 ! no interaction

c     if non-survival interaction

c     determine the event, using the random number
c     REWEIGHTING: determine the event during generation and record the interaction type
      if(yr.le.pel) then
         iqq=1                  !elastic
         NFiflagStep(NFNStep) = 1

      else if(yr.le.pel+psp) then
         iqq=2                  !single production
         NFiflagStep(NFNStep) = 2

      else
         iqq=3                 !double production
         NFiflagStep(NFNStep) = 3

      end if

      
      endif
      
c     iqq has been set to 1,2,3 for our different interaction types
      
C     Depending which interaction occurs, pcascprob gets multiplied by the now normalised factor
c     we also keep a record of which event occurred

c     REWEIGHTING: depending which event type occurred, change pcascprob by the right factor
      if(NReweightFlag.eq.1) then


         if(NFiflagStep(RW_step_i).eq.4) then
            iqq =4              !survival
            pcascprob = pcascprob * psurvnuc
c            survivalcounter = survivalcounter + 1


            goto 200
         elseif(NFiflagStep(RW_step_i).eq.1) then
            iqq =1              !elastic
         elseif(NFiflagStep(RW_step_i).eq.2) then
            iqq =2              !single pi
         elseif(NFiflagStep(RW_step_i).eq.3) then
            iqq =3              !double pi 
         else

         endif
c     if suriva;                                                                                                                                                                            

c     if non-survival interaction                                                                                                                                                            
c     determine the event, using the random number                                                                                                                                           
 


      endif

c     REWEIGHTING: as above multiply pcascprob by the right factor depending on 
c     which interaction type occurred 
      If(iqq==1) then 
         pcascprob = pcascprob * pelnuc;
c         elasticcounter = elasticcounter +1;
      else If(iqq==2) then
         pcascprob = pcascprob * pspnuc;
c         singlepicounter = singlepicounter +1;
      else If(iqq==3) then
         pcascprob = pcascprob * pdpnuc;
c         twopicounter = twopicounter +1;
      else
         unknowntype = unknowntype +1;

      endif



c     make kinematics and final state

C     interaction type, total energy, charge, momentum, charge

      if(NReweightFlag.eq.0) then

      

      call nrprd(iqq,eineq,ichint,piso,chiso) 
      endif


c     only piso and chiso get set

c     chiso(1,2,etc) gets set to charge from ka1,ka2,etc depening on p-p, p-n, ...
c     piso just set from lorentz transforms of pdeldec

      ! check for Pauli blocking
c     not sure what this bit does
c     p1 and p2 are total momenta (for nucleon and other nucleon respectively?)
      p1=sqrt(piso(1,1)**2+piso(2,1)**2+piso(3,1)**2)
      p2=sqrt(piso(1,3)**2+piso(2,3)**2+piso(3,3)**2)
c     if these are smaller than pfermi, set idec to 11,12,13 amd QUIT

      if(p1.lt.pfermi.or.p2.lt.pfermi) then
        idec=10+iqq             ! return code for interaction that would have taken place, needed for reweighting systematics


        goto 200
      end if
c     REWEIGHTING: End without checking nel(iqq) because iqq can be 4 now - would give an error
      if(NReweightFlag.eq.1) then
         goto 200
      endif
c     otherwise set idec to -1,-2,-3 and increment nel for the interaction type we did      
      idec=-iqq                 !tell what took place
c      print *,"iqq = ",iqq      
      ihap=1                    !make some bookkeep
      nel(iqq)=nel(iqq)+1
c     idec only set to 1 if survival, I think


 200  return

       end
      
      
