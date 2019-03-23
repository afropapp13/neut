	  subroutine nrpion(ides)
	  implicit none
	  integer*4 ides

#include "nrint.h"
C     common /nucint/uin,chint,ein,pin(3),upin(3)
#include "nrcms.h"
C     common /nuccms/ rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
#include "nrnewprd.h"
C     common /nucnewprd/piso(4,4),chiso(4)
#include "nrhad.h"
C     common /nuchad/ ihap,nel(3)
#include "nrhadr.h"
C     common /hadron/inot,itrans,idec,ine1,ine2,iel,iabs,icx

C     dimension xabs(21),xpeln(126),xpelp(126),xpelm(126),xpchx(126),
C     1	  xspp(117),xspn(117),xspm(117),
C     2	  dang(5,21,4)
	  real*4 xabs(21),xpeln(126),xpelp(126),xpelm(126),xpchx(126),
	1	   xspp(117),xspn(117),xspm(117),
	2	   dang(5,21,4)


C     dimension angld(18),angldi(19),aa(5),b1(51),b2(51),b3(51)
	  real*4 angld(18),angldi(19),aa(5),b1(51),b2(51),b3(51)
C     dimension pout1l(3),upoutc(3)
	  real*4 pout1l(3),upoutc(3)
C     logical prot,notuj
	  logical notuj

	  real*4    phs31(21),php31(21),php33(21)
C	  dimension phs31(21),php31(21),php33(21)
	  real*4    phs11(21),php11(21),php13(21)
C	  dimension phs11(21),php11(21),php13(21)
	  real*4    pep11(21),phd13(21),phd15(21)
C	  dimension pep11(21),phd13(21),phd15(21)
	  real*4     gm1(4),gm3(4)
C	  dimension gm1(4),gm3(4)
	  real nucprbin

	  data gm1/0.,0.66667,-0.4714,0.33333/
	  data gm3/1.,0.33333, 0.4714,0.66667/

	  real*4 angdata(51,428)

	  data phs11/
	1	   0.,   4.93, 6.51, 7.46, 8.08, 8.50, 8.79, 9.00, 9.17, 9.33,
	2	   9.49, 9.69, 9.93,10.23,10.61,11.06,11.62,12.28,13.05,13.93,
	3	   14.94/
	  data php11/
	1	   0.00,-0.68,-1.49,-2.06,-2.29,-2.18,-1.78,-1.13,-0.26, 0.82,
	2	   2.12, 3.71, 5.67, 8.18,11.41,15.58,20.87,23.92,28.09,32.26,
	3	   36.43/
	  data php13/
	1	   0.,  -0.15,-0.40,-0.68,-0.97,-1.26,-1.53,-1.79,-2.04,-2.28,
	2	   -2.52,-2.78,-3.05,-3.34,-3.69,-4.08,-4.55,-5.10,-5.76,-6.54,
	3	   -7.46/
	  data phs31/
	1	   0.00,-2.85,  -4.28, -5.80, -7.30,- 8.80,-10.23,-11.71,-13.16,
	2	   -14.57,-15.93,-17.22,-18.44,-19.56,-20.59,-21.51,-22.30,-23.02,
	3	   -23.61,-24.09,-24.46/
	  data php31/
	1	   0.00,-0.22,-0.60,-1.07,-1.60,-2.17,-2.77,-3.39,-4.01,-4.63,
	2	   -5.24,-5.84,-6.41,-6.96,-7.48,-7.96,-8.41,-9.16,-9.71,-10.25,
	3	   -10.80/
	  data php33/
	1	   0.00,  1.29,  3.95,  7.95, 13.62, 21.56, 32.53, 47.17,
	2	   66.10, 81.98, 97.51,109.47,118.42,125.19,130.44,134.63,
	3	   138.07,140.97,143.45,145.62,147.54/
	  data pep11/
	1	   1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,1.   ,
	2	   1.   ,0.996,0.986,0.970,0.947,0.908,0.861,0.810,0.762,0.710,
	3	   0.660/
	  data phd13/
	1	   0.  , 0.00, 0.01, 0.04, 0.08, 0.15, 0.25, 0.38, 0.56, 0.78,
	2	   1.05, 1.38, 1.77, 2.24, 2.78, 3.41, 4.13, 4.94, 5.86, 6.89,
	3	   8.04/
	  data phd15/
	1	   0.  , 0.00, 0.02, 0.06, 0.12, 0.20, 0.30, 0.42, 0.57, 0.73,
	2	   0.91, 1.10, 1.31, 1.52, 1.73, 1.95, 2.16, 2.36, 2.54, 2.71,
	3	   2.85/

C     PARTIAL AND	TOTAL CROSS SECTIONS (*100 mb)
c     pi0-p elastic
	  data xpeln/     
	1	   0.,  90., 258., 745.,1575.,2760.,4240.,6190.,7900.,8930.,
	2	   7430.,6480.,5540.,4633.,3725.,2985.,2490.,2160.,1908.,1743.,
	3	   1565.,1465.,1380.,1300.,1260.,1228.,1200.,1200.,1205.,1235.,
	4	   1238.,1140.,1050., 980., 925., 930., 930., 950., 980.,1000.,
	5	   1040.,1095.,1170.,1280.,1413.,1568.,1580.,1545.,1513.,1495.,
	6	   1445.,1410.,1380.,1350.,1320.,1285.,1280.,1270.,1270.,1280.,
	7	   1295.,1320.,1340.,1370.,1390.,1405.,1393.,1360.,1340.,1310.,
	8	   1280.,1250.,1220.,1190.,1180.,1163.,1150.,1150.,1150.,1140.,
	9	   1125.,1125.,1120.,1120.,1110.,1110.,1100.,1090.,1080.,1050.,
	1	   1018., 990., 980., 950., 930., 883., 880., 870., 860., 840.,
	2	   835., 820., 810., 800., 790., 800., 780., 775., 770., 765.,
	3	   765., 750., 740., 740., 740., 738., 735., 730., 725., 725.,
	4	   725., 715., 710., 710., 710., 710./

c     pi-p charge exchange
	  data xpchx/
	1	   390., 480., 620., 875.,1240.,1675.,2200.,3100.,4350.,4460.,
	2	   4000.,3400.,2950.,2625.,2300.,2000.,1650.,1485.,1375.,1300.,
	3	   1215.,1150.,1080.,1030., 975., 920., 870., 815., 765., 718.,
	4	   675., 630., 595., 555., 515., 490., 470., 451., 450., 468.,
	5	   500., 550., 625., 680., 675., 575., 500., 430., 375., 345.,
	6	   300., 285., 265., 250., 240., 225., 210., 200., 195., 185.,
	7	   178., 175., 170., 165., 160., 155., 152., 150., 149., 148.,
	8	   145., 140., 135., 130., 125., 120., 118., 115., 110., 107.,
	9	   104., 100.,  99.,  98.,  95.,  90.,  85.,  80.,  78.,  77.,
	1	   75.,  73.,  72.,  69.,  67.,  65.,  63.,  61.,  58.,  55.,
	2	   52.,  51.,  50.,  50.,  50.,  50.,  50.,  49.,  48.,  47.,
	3	   46.,  44.,  42.,  39.,  37.,  35.,  34.,  33.,  32.,  31.,
	4	   29.,  27.,  25.,  24.,  24.,  23./

c     pi-p elastic
	  data xpelm/
	1	   125., 150., 175., 245., 380., 600., 970.,1500.,2140.,2310.,
	2	   2295.,2070.,1795.,1550.,1360.,1230.,1130.,1070.,1050.,1070.,
	3	   1120.,1175.,1235.,1300.,1400.,1500.,1600.,1700.,1835.,1970.,
	4	   2050.,1915.,1770.,1650.,1570.,1520.,1510.,1525.,1550.,1600.,
	5	   1685.,1800.,2000.,2230.,2475.,2635.,2510.,2300.,2140.,2000.,
	6	   1870.,1750.,1670.,1585.,1505.,1440.,1395.,1340.,1299.,1260.,
	7	   1215.,1175.,1140.,1099.,1060.,1040.,1010., 999., 990., 975.,
	8	   960., 955., 945., 935., 925., 925., 935., 965., 985.,1000.,
	9	   1015.,1030.,1060.,1080.,1095.,1100.,1095.,1090.,1070.,1035.,
	1	   1000., 960., 905., 855., 820., 800., 785., 780., 775., 770.,
	2	   765., 760., 759., 757., 756., 755., 750., 745., 740., 735.,
	3	   730., 728., 725., 725., 724., 723., 722., 722., 722., 721.,
	4	   720., 715., 710., 710., 710., 710./

c     pi+p elastic
	  data xpelp/
	1	   180., 400., 990.,2170.,4000.,6580.,9680.,13920.,18000.,20000.,
	2	   16550.,14200.,12250.,10320.,8400.,6725.,5510.,4725.,4130.,3690.,
	3	   3230.,2885.,2600.,2300.,2090.,1875.,1675.,1500.,1340.,1200.,
	4	   1100., 998., 920., 860., 820., 810., 810., 825., 850., 875.,
	5	   900., 940., 975.,1000.,1030.,1075.,1130.,1200.,1275.,1330.,
	6	   1350.,1335.,1330.,1330.,1345.,1355.,1380.,1400.,1460.,1500.,
	7	   1555.,1625.,1700.,1800.,1875.,1920.,1925.,1890.,1830.,1790.,
	8	   1725.,1690.,1640.,1600.,1550.,1505.,1475.,1430.,1400.,1365.,
	9	   1335.,1300.,1280.,1250.,1225.,1205.,1195.,1175.,1150.,1135.,
	1	   1105.,1095.,1080.,1060.,1030.,1020.,1005., 990., 980., 970.,
	2	   960., 940., 930., 920., 900., 898., 890., 880., 870., 860.,
	3	   840., 830., 820., 810., 800., 795., 785., 780., 775., 765.,
	4	   760., 755., 750., 745., 740., 735./

c     pi+p single pi production
	  data xspp/
	1	   0.,   0.,   1.,   3.,   5.,   9.,  13.,  21.,  30.,  41.,
	2	   55.,  75.,  98., 122., 153., 183., 218., 255., 293., 335.,
	3	   383., 445., 530., 650., 780., 880., 990.,1082.,1170.,1251.,
	4	   1321.,1380.,1425.,1460.,1490.,1512.,1520.,1500.,1445.,1395.,
	5	   1330.,1330.,1388.,1450.,1521.,1595.,1690.,1755.,1842.,1930.,
	6	   2012.,2090.,2145.,2180.,2180.,2160.,2155.,2170.,2195.,2220.,
	7	   2222.,2220.,2205.,2178.,2140.,2090.,2030.,1990.,1955.,1925.,
	8	   1900.,1890.,1862.,1850.,1838.,1830.,1820.,1822.,1823.,1828.,
	9	   1833.,1840.,1850.,1859.,1870.,1883.,1900.,1910.,1928.,1943.,
	1	   1960.,1980.,2000.,2018.,2038.,2058.,2080.,2098.,2119.,2139.,
	2	   2158.,2180.,2198.,2215.,2235.,2252.,2270.,2278.,2286.,2295.,
	3	   2300.,2305.,2313.,2319.,2321.,2325.,2330./

c     pi0p single pi production
	  data xspn/
	1	   0.,   0.,   1.,   8.,  23.,  42.,  69., 100., 135., 175.,
	2	   221., 270., 320., 372., 430., 492., 560., 640., 755.,1085.,
	3	   1210.,1252.,1350.,1405.,1395.,1365.,1330.,1365.,1455.,1610.,
	4	   1755.,1928.,2105.,2240.,2260.,2219.,2170.,2115.,2045.,1965.,
	5	   1879.,1825.,1775.,1750.,1752.,1790.,1835.,1885.,1938.,1995.,
	6	   2055.,2125.,2172.,2210.,2230.,2235.,2248.,2265.,2288.,2302.,
	7	   2310.,2317.,2312.,2300.,2280.,2259.,2230.,2199.,2170.,2145.,
	8	   2121.,2100.,2081.,2067.,2055.,2045.,2041.,2041.,2048.,2062.,
	9	   2085.,2110.,2140.,2170.,2200.,2225.,2240.,2252.,2265.,2278.,
	1	   2290.,2305.,2319.,2332.,2348.,2360.,2372.,2388.,2400.,2412.,
	2	   2425.,2438.,2448.,2458.,2468.,2475.,2485.,2488.,2490.,2489.,
	3	   2488.,2485.,2482.,2479.,2475.,2470.,2465./

c     pi-p siggle pi production
	  data xspm/
	1	   0.,   5.,  10.,  20.,  39.,  68., 110., 160., 230., 305.,
	2	   382., 470., 548., 630., 712., 805., 905.,1021.,1230.,1820.,
	3	   2035.,2060.,2170.,2160.,2040.,1840.,1670.,1635.,1740.,1965.,
	4	   2200.,2435.,2740.,3020.,3030.,2920.,2815.,2730.,2639.,2530.,
	5	   2420.,2315.,2190.,2060.,1992.,1979.,1985.,2008.,2045.,2085.,
	6	   2128.,2160.,2205.,2241.,2275.,2308.,2340.,2358.,2378.,2392.,
	7	   2400.,2410.,2420.,2422.,2425.,2425.,2425.,2410.,2380.,2358.,
	8	   2335.,2310.,2295.,2282.,2272.,2265.,2260.,2265.,2279.,2302.,
	9	   2335.,2375.,2425.,2471.,2512.,2545.,2580.,2590.,2602.,2612.,
	1	   2622.,2630.,2643.,2650.,2659.,2665.,2670.,2679.,2685.,2692.,
	2	   2698.,2700.,2703.,2705.,2702.,2701.,2700.,2692.,2688.,2679.,
	3	   2671.,2660.,2651.,2640.,2628.,2613.,2600./

C     PI DEUTERON	ABSORPTION (in mb (without factor 100))
	  DATA XABS/
	1	   2., 3., 5.,12.,18.,27.,37.,50.,60.,45.,
	2	   30.,20.,10., 5., 3., 2., 2., 1., 1., 1.,0./
c     1  0.5,1.0,2.0,4.0,7.0,15.0,23.0,36.0,60.0,50.0,
c     2 12.0,4.0,1.0,0.2,0.2, 0.2, 0.2, 0.2, 0.2, 0.2,0/

	  integer*4 iraz
	  data iraz/0/
	  
	  real*4 upi,amn,amn2,step
	  real*4 rhod,fac
	  real*4 nangl

	  parameter (upi=139.6, amn=939., amn2=amn**2, step=0.2 )
	  parameter (rhod=9./(2.*3.141*9.5**3) ,fac=1e-2)
	  parameter (nangl =18)

	  integer*4 i,j,k
	  real*4 t,tl,el,pl,pc,pk,ppel,ppchx,psp,erem,erem1,ppabs
	  real*4 s11,s31,p31,p11,p13,p33,ep11,d13,d15
	  real*4 w01,w02,w03,w11,w12,w13,w21,w22,w23,w32,w33,w43
	  real*4 g1,g3
	  real*4 pcmsm,eineq,dist,psum,erem2,ptot,xr,yr
	  real*4 xr,yr,angb,cab,ang,xxr,cangl,errm,phi
	  real*4 pout1lm,pout2lm,eout1l
	  integer*4 ilo,ihi,i1,ilo1,ihi1,ilo2,ihi2,ill,ihh,ncx,inel
	  integer*4 ichint,ka,k1,ipath,krec1,krec2,krecsav
	  

	  real*4   ranf
	  external ranf


	  if(iraz.eq.0)then
		 iraz=1
		 DO  I=1,126
			xpeln(i)=xpeln(i)*fac
			xpelp(i)=xpelp(i)*fac
			xpelm(i)=xpelm(i)*fac
			xpchx(i)=xpchx(i)*fac
			if(i.le.117)then
			   xspp(i)=xspp(i)*fac
			   xspn(i)=xspn(i)*fac
			   xspm(i)=xspm(i)*fac
			end if
		 end do
		 do i=1,20
			tl=i*20.
			el=tl+upi
			pl=sqrt(el**2-upi**2)
			pc=pl/sqrt(1.+2*el/amn+(upi/amn)**2)/197.329
			pk=pc**2
			call nrhis(tl,20.,21,ilo,ihi,erem)

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
	1			 -ep11*(sind(p31-p11)**2-sind(p33-p11)**2)
	2			 +sind(p31-p13)**2-sind(p33-p13)**2
	3			 +sind(s31-d13)**2+3/2.*sind(s31-d15)**2
	4			 -sind(d13)**2-3/2.*sind(d15)**2) !*g1*g3

			w03=(-3/2.*sind(s11)**2+ep11*sind(p11-p13)**2+(ep11-1)**2/4.
	1			 +sind(s11-d13)**2+3/2.*sind(s11-d15)**2
	2			 -3/2.*sind(d13-d15)**2+3/2.*sind(d13)**2
	3			 +9/4.*sind(d15)**2) !*g1**2

			w11=2*(  sind(s31)*sind(p31)*cosd(s31-p31)
	1			 +2*sind(s31)*sind(p33)*cosd(s31-p33)) !*g3**2

			w12=(ep11*(sind(p11)**2-sind(s31-p11)**2)
	1			 +sind(p31)**2-sind(s11-p31)**2
	2			 +2*(sind(p33)**2+sind(p13)**2-sind(s31-p13)**2
	3			 -sind(s11-p33)**2)
	4			 +3*(sind(s11)**2+sind(s31)**2))
	5			 +(-2*sind(p31-d13)**2+9/2.*sind(p31-d15)**2-5/2.*sind(p31)**2
	6			 +5*sind(p33-d13)**2-5*sind(p33)**2-3*sind(d13)**2
	7			 -9/2.*sind(d15)**2) !*g1*g3

			w13=(ep11*(sind(p11)**2-sind(s11-p11)**2)
	1			 +2*(sind(p13)**2-sind(s11-p13)**2)+3*sind(s11)**2
	2			 +ep11*(-2*sind(p11-d13)**2+9/2.*sind(p11-d15)**2
	3			 -5/2.*sind(p11)**2)
	4			 +5*sind(p13-d13)**2-5*sind(p13)**2-3*sind(d13)**2
	5			 -9/2.*sind(d15)**2) !*g1**2

			w21=3*(sind(p33)**2+2*sind(p31)*sind(p33)*cosd(p31-p33)) !*g3**2

			w22=     +3*(ep11*(sind(p11)**2-sind(p33-p11)**2)
	1			 +    2*(sind(p13)**2+sind(p33)**2)+sind(p31)**2
	2			 -      (sind(p33-p13)**2+sind(p31-p13)**2)
	3			 -sind(s31-d13)**2-3/2.*sind(s31-d15)**2+5/2.*sind(s31)**2
	4			 +sind(d13)**2+3/2.*sind(d15)**2) !*g1*g3

			w23=   +  3*(2*sind(p13)**2+ep11*(sind(p11)**2-sind(p11-p13)**2)
	1			 -sind(s11-d13)**2-3/2.*sind(s11-d15)**2+5/2.*sind(s11)**2
	2			 +6*sind(d13-d15)**2-4*sind(d13)**2-6*sind(d15)**2) !*g1**2

			w32=(-15/2.*sind(p31-d15)**2+15/2.*sind(p31)**2-9*sind(p33-d13)**2
	1			 -6*sind(p33-d15)**2+15*sind(p33)**2+9*sind(d13)**2
	2			 +27/2.*sind(d15)**2) !*g1*g3

			w33= +(15/2.*ep11*(-sind(p11-d15)**2+sind(p11)**2)
	1			 -9*sind(p13-d13)**2-6*sind(p13-d15)**2+15*sind(p13)**2
	2			 +9*sind(d13)**2+27/2.*sind(d15)**2) !*g1**2

			w43= 45/2.*(sind(d13)**2-sind(d13-d15)**2)+135/4.*sind(d15)**2 !*g1**2
			do j=1,4
			   g1=gm1(j)
			   g3=gm3(j)
			   i1=i+1
			   dang(1,i1,j)=(w01*g3**2+w02*g1*g3+w03*g1**2)/pk
			   dang(2,i1,j)=(w11*g3**2+w12*g1*g3+w13*g1**2)/pk
			   dang(3,i1,j)=(w21*g3**2+w22*g1*g3+w23*g1**2)/pk
			   dang(4,i1,j)=(          w32*g1*g3+w33*g1**2)/pk
			   dang(5,i1,j)=(	 		  w43*g1**2)/pk
			end do
		 end do

		 open(33,status='old',readonly,shared,
	1		  file='nucang.dat')

		 do i=1,428
			read(33,777)(angdata(j,i),j=1,51)
		 end do
 777	 format(10f10.7)
		 close(33)

	  end if					!end of initialization

	  if(ides.eq.10)then
		 notuj=.true.
	  else
		 notuj=.false.
	  end if
	  ides=1

C     MAGNITUDE OF CMS MOMENTUM
	  IF(ECMS2.LE.(UPI+UNUCL)**2) GO TO 200
	  PCMSM=SQRT((ECMS2-(UPI-UNUCL)**2)*(ECMS2-(UPI+UNUCL)**2)/(4.0*ECMS2))

C     EQUIVALENT PION	ENERGY FOR STATIONARY TARGET
	  EINEQ=(ECMS2-UNUCL**2-UPI**2)/(2.0*UNUCL)
	  IF(EINEQ.LT.UPI)	GO TO 200
	  t=eineq-upi

C     ABSORPTION	CROSS SECTION WEIGHTED BY RATIO	OF LOCAL
C     NUCLEON DENSITY TO	DEUTERON DENSITY.AT THIS STAGE RHON
C     CONTAINS 0.1 FACTOR FOR UNITS.THE RHOD FACTOR ADJUSTMENT
C     THROUGH PUTTING THE "DEUTERON RADIUS" AT 9.5 TAKES	THIS INTO
C     ACCOUNT.PUTTING IN	THE 0.1	FACTOR INTO RHOD IMPLIES THAT 9.5
C     CORRESPONDS TO A SO CALLED	"QUASI DEUTERON" RADIUS	OF 4.4 FM.

	  call nrhis(t,20.,21,ilo,ihi,erem)
	  PPABS=XABS(ihi)*erem+xabs(ilo)*(1.-erem) !*0.5*rhon/rhod

	  call nrhis(t,20.,126,ilo,ihi,erem)
	  call nrhis(t-180.,20.,117,ilo1,ihi1,erem1)

	  ichint=chint
	  if(ichint)1,2,3
 1	  continue					!pi-
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
 2	  continue					!pi0
	  ppel =(xpeln(ihi)*erem+xpeln(ilo)*(1-erem))
	  ppchx=(xpchx(ihi)*erem+xpchx(ilo)*(1-erem))
	  psp  =(xspn(ihi1)*erem1+xspn(ilo1)*(1-erem1))
	  ppabs=0.7*ppabs
	  ka=4
	  go to 4
 3	  continue					!pi+
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
 4	  continue
	  if(ecms2.lt.1488400.)psp=0.

	  PTOT=PPEL+PPCHX+PPABS+psp

C********CHOOSE INTERACTION DISTANCE
	  XR=ranf()+1.e-6
	  YR=ranf()*PTOT
C********HERE WE NEED TO NORMALIZE
	  PTOT=PTOT*RHON
	  if(ptot.eq.0.)go to 200
	  DIST=-LOG(XR)/ptot

C********IF NOTHING HAPPENS
	  IF(STEP.LT.DIST) GO TO 200

C********IS IT ABSORBED BY DI NUCLEON?
	  PSUM=PPABS
	  IF(YR.Le.PSUM) THEN

		 if(notuj)IABS=IABS+1
		 Uin=1.e-2
		 CHINT=0.0
		 ihap=1
		 ides=0
		 go to 200
	  END IF
c*********does pion produce an additional pion ?
	  psum=psum+psp
	  if(yr.lt.psum)go to 300

C*********DOES  PION SCATTER	OFF  FERMI NUCLEON?
	  PSUM=PSUM+PPEL
	  IF(YR.LT.PSUM) THEN
		 IPATH=0
	  else
C*********PION CHARGE EXCHANGE ON FERMI NUCLEON
		 IPATH=1
		 ka=3
	  END IF


C*********COME HERE FOR PI-P KINEMATICS WITH	FERMI MOTION
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
			ANGLD(K)=SINd(ANGB)*(aa(1)+ang)
		 END DO
		 CALL nrCFD(ANGLDI,NANGL,ANGLD)

C     DO SCATTERING IN CMS
		 xXR=ranf()
		 ANG=nucPRBIN(ANGLDI,NANGL,0.,10.,xXR)
		 CANGL=COSd(ANG)
	  else
		 call nrhis(t-400.,20.,107,ilo2,ihi2,erem2)
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
		 call nrhis(xxr,0.02,51,ill,ihh,errm)
		 cangl=b3(ill)*errm+b3(ihh)*(1-errm)
	  end if
	  PHI=-1
	  CALL nrCONE(UP4CMS,UPOUTC,CANGL,PHI)
	  call nrkin(upoutc,pcmsm,upi,pout2lm,pout1l,pout1lm,eout1l)

	  IF(POUT2LM.LT.PFERMI) GO TO 200

	  DO J=1,3
		 PIN(J)=POUT1L(J)
		 UPIN(J)=PIN(J)/POUT1LM
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
C		 nel=1
		 nel(1)=1
	  END IF
	  ihap=1

 200  return

 300  continue					!pion production
	  call nrpiprd(eineq,ka,piso,chiso)
	  pout2lm=sqrt(piso(4,1)**2-amn2)
	  IF(POUT2LM.LT.PFERMI) GO TO 200
	  pout1lm=sqrt(piso(1,3)**2+piso(2,3)**2+piso(3,3)**2)
	  DO J=1,3
		 PIN(J)=piso(j,3)
		 UPIN(J)=PIN(J)/POUT1LM
	  END DO
	  EIN=piso(4,3)
	  chint=chiso(3)
	  ides=-1
	  ihap=1
	  if(notuj)inel=inel+1
	  go to 200

	  end
	  
