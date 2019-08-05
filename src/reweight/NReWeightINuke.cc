//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightINuke

\brief    Reweighting NEUT Nucleon FSI Cascade Model

\author   Martin Hierholzer <martin@hierholzer.info>
          University of Bern

\created  Nov 04, 2014
*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "NReWeightControls.h"
#include "NReWeightINuke.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

using namespace neut;
using namespace neut::rew;

using std::cout;
using std::endl;

extern "C" {
  float effrmgas_(float *dum1, float *dum2, float *r);
}


//
// constants and tables, need to be in sync with nrprton.F and nrnuc.F
//
const double step = 0.2;            // step size from nrprton.F


//     p-p elastic
const double xepp[176] = {
          67500.,15500.,6750.,4420.,3230.,2800.,2520.,2370.,2300.,2275.,
          2260.,2260.,2260.,2260.,2270.,2280.,2295.,2300.,2310.,2330.,
          2350.,2380.,2395.,2420.,2460.,2485.,2500.,2530.,2565.,2600.,
          2620.,2640.,2660.,2675.,2690.,2700.,2705.,2710.,2715.,2720.,
          2725.,2725.,2720.,2715.,2710.,2700.,2695.,2680.,2670.,2660.,
          2640.,2625.,2605.,2590.,2570.,2545.,2525.,2500.,2480.,2470.,
          2450.,2430.,2410.,2395.,2370.,2360.,2340.,2325.,2305.,2290.,
          2275.,2270.,2260.,2250.,2230.,2225.,2210.,2200.,2195.,2190.,
          2175.,2165.,2150.,2140.,2125.,2120.,2105.,2100.,2090.,2075.,
          2065.,2055.,2045.,2030.,2020.,2005.,2000.,1995.,1980.,1975.,
          1965.,1960.,1945.,1935.,1925.,1910.,1900.,1895.,1885.,1875.,
          1865.,1860.,1850.,1830.,1825.,1810.,1800.,1795.,1785.,1775.,
          1765.,1755.,1745.,1730.,1720.,1715.,1700.,1695.,1680.,1675.,
          1665.,1650.,1645.,1630.,1615.,1605.,1600.,1590.,1580.,1565.,
          1555.,1550.,1530.,1520.,1510.,1500.,1490.,1475.,1465.,1455.,
          1445.,1440.,1425.,1415.,1405.,1398.,1390.,1380.,1370.,1360.,
          1350.,1335.,1325.,1315.,1300.,1295.,1280.,1270.,1260.,1250.,
          1240.,1225.,1210.,1200.,1195.,1175. };

//     p-n elastic
const double xepn[176] = {
          196500.,47500.,22000.,13000.,9180.,7300.,6030.,5180.,4680.,
          4320.,4080.,3910.,3760.,3650.,3550.,3480.,3415.,3370.,3325.,
          3290.,3275.,3250.,3255.,3275.,3285.,3275.,3220.,3150.,3075.,
          2990.,2875.,2775.,2695.,2630.,2590.,2565.,2560.,2560.,2560.,
          2565.,2570.,2575.,2578.,2580.,2585.,2580.,2575.,2560.,2540.,
          2505.,2470.,2425.,2375.,2315.,2275.,2230.,2200.,2175.,2155.,
          2145.,2130.,2125.,2115.,2105.,2100.,2095.,2090.,2080.,2070.,
          2060.,2050.,2045.,2040.,2030.,2025.,2020.,2015.,2010.,2005.,
          2002.,2000.,1999.,1990.,1985.,1975.,1970.,1965.,1960.,1950.,
          1945.,1940.,1925.,1920.,1915.,1910.,1900.,1898.,1895.,1890.,
          1880.,1875.,1870.,1865.,1855.,1850.,1845.,1830.,1825.,1820.,
          1810.,1805.,1800.,1795.,1790.,1785.,1775.,1770.,1765.,1750.,
          1745.,1740.,1725.,1720.,1715.,1705.,1699.,1690.,1685.,1675.,
          1670.,1650.,1640.,1630.,1620.,1610.,1600.,1595.,1580.,1575.,
          1560.,1550.,1530.,1520.,1510.,1500.,1490.,1485.,1470.,1460.,
          1450.,1440.,1435.,1425.,1415.,1410.,1400.,1390.,1385.,1375.,
          1365.,1350.,1340.,1330.,1320.,1310.,1300.,1290.,1280.,1270.,
          1260.,1250.,1240.,1230.,1220.,1210.,1200. };

//     p-p single pi production
const double xspp[158] = {
          13.,  26.,  50.,  90., 126., 175., 235., 300., 380., 480.,
          585., 710., 860.,1035.,1242.,1436.,1580.,1695.,1785.,1866.,
          1937.,2000.,2050.,2090.,2128.,2160.,2188.,2208.,2228.,2240.,
          2248.,2252.,2253.,2253.,2253.,2252.,2250.,2248.,2245.,2240.,
          2238.,2234.,2230.,2225.,2220.,2215.,2210.,2205.,2200.,2195.,
          2190.,2185.,2178.,2172.,2168.,2160.,2154.,2150.,2140.,2135.,
          2128.,2120.,2112.,2102.,2096.,2089.,2079.,2070.,2061.,2054.,
          2042.,2036.,2028.,2016.,2005.,1995.,1985.,1975.,1968.,1959.,
          1949.,1939.,1928.,1919.,1908.,1898.,1888.,1878.,1868.,1858.,
          1848.,1838.,1828.,1817.,1807.,1797.,1786.,1777.,1765.,1758.,
          1748.,1737.,1728.,1716.,1705.,1696.,1685.,1675.,1665.,1655.,
          1645.,1635.,1627.,1615.,1605.,1595.,1584.,1575.,1563.,1555.,
          1542.,1532.,1520.,1512.,1500.,1492.,1480.,1470.,1460.,1450.,
          1440.,1430.,1419.,1408.,1398.,1387.,1377.,1367.,1356.,1344.,
          1335.,1324.,1315.,1305.,1295.,1285.,1275.,1265.,1257.,1245.,
          1237.,1228.,1218.,1208.,1198.,1188.,1178.,1169. };

//     p-n single pi production
const double xspn[158] = {
          18.,  28.,  48.,  75., 117., 159., 208., 270., 340., 440.,
          560., 720., 910.,1060.,1170.,1260.,1325.,1375.,1390.,1396.,
          1397.,1396.,1393.,1390.,1385.,1382.,1380.,1376.,1372.,1366.,
          1360.,1357.,1352.,1346.,1341.,1338.,1332.,1325.,1320.,1316.,
          1310.,1304.,1300.,1295.,1288.,1281.,1278.,1270.,1264.,1259.,
          1254.,1245.,1240.,1235.,1230.,1223.,1218.,1212.,1203.,1199.,
          1193.,1184.,1180.,1175.,1170.,1162.,1158.,1152.,1146.,1140.,
          1136.,1130.,1122.,1120.,1112.,1108.,1100.,1096.,1090.,1082.,
          1078.,1070.,1065.,1060.,1056.,1050.,1045.,1040.,1034.,1028.,
          1021.,1018.,1010.,1005.,1000., 995., 988., 980., 977., 970.,
          962., 958., 951., 945., 940., 935., 930., 922., 918., 910.,
          904., 900., 895., 889., 882., 878., 872., 865., 860., 855.,
          850., 842., 838., 832., 826., 820., 817., 810., 804., 800.,
          797., 790., 783., 780., 774., 768., 761., 758., 752., 744.,
          740., 735., 730., 722., 720., 715., 710., 705., 698., 694.,
          690., 683., 680., 675., 670., 662., 660., 653. };

//     p-p double pi production
const double xdpp[130] = {
          0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   8.,
          9.,   9.,  11.,  16.,  21.,  38.,  60.,  90., 112., 130.,
          146., 160., 177., 190., 200., 215., 228., 240., 252., 265.,
          280., 292., 305., 320., 335., 349., 361., 379., 392., 405.,
          420., 435., 450., 460., 479., 490., 504., 520., 535., 550.,
          565., 580., 595., 610., 623., 640., 655., 670., 680., 698.,
          710., 726., 740., 756., 770., 788., 800., 820., 834., 850.,
          860., 880., 893., 910., 923., 940., 957., 970., 986.,1002.,
          1020.,1035.,1050.,1066.,1080.,1098.,1112.,1128.,1141.,1160.,
          1175.,1190.,1208.,1221.,1239.,1255.,1270.,1285.,1300.,1320.,
          1338.,1355.,1370.,1388.,1400.,1420.,1436.,1450.,1466.,1480.,
          1496.,1510.,1525.,1540.,1558.,1573.,1588.,1600.,1620.,1637.,
          1652.,1670.,1682.,1700.,1718.,1731.,1748.,1762.,1780.,1798. };

//     p-n double pi production
const double xdpn[130] = {
          0.,  40.,  80., 128., 170., 230., 300., 370., 450., 520.,
          560., 600., 630., 660., 690., 715., 740., 760., 785., 805.,
          828., 850., 870., 888., 908., 928., 943., 962., 980., 998.,
          1017.,1030.,1045.,1063.,1080.,1095.,1110.,1125.,1140.,1156.,
          1170.,1185.,1200.,1215.,1230.,1240.,1254.,1266.,1280.,1293.,
          1305.,1320.,1333.,1345.,1360.,1373.,1382.,1397.,1410.,1420.,
          1437.,1448.,1460.,1473.,1488.,1498.,1510.,1521.,1537.,1548.,
          1560.,1572.,1582.,1598.,1610.,1620.,1635.,1645.,1658.,1670.,
          1680.,1696.,1708.,1720.,1732.,1743.,1757.,1770.,1780.,1797.,
          1808.,1820.,1834.,1844.,1860.,1870.,1880.,1898.,1910.,1920.,
          1938.,1950.,1960.,1978.,1990.,2002.,2020.,2032.,2048.,2060.,
          2077.,2090.,2104.,2120.,2130.,2150.,2162.,2180.,2195.,2210.,
          2225.,2240.,2257.,2270.,2286.,2300.,2315.,2330.,2342.,2360. };





//_______________________________________________________________________________________
NReWeightINuke::NReWeightINuke()
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightINuke::~NReWeightINuke()
{
}
//_______________________________________________________________________________________
void NReWeightINuke::Init(void)
{
  fMFPDef = 1.0;

  fMFPCurr  = fMFPDef;

  fMFPDial = fMFPDef;

}
//_______________________________________________________________________________________
bool NReWeightINuke::IsHandled(NSyst_t syst)
{
   bool handle;

   switch(syst) {
     case ( kINukeTwkDial_MFP_N ) :
          handle = true;
          break;

     default:
          handle = false;
	  break;
   }

   return handle;
}
//_______________________________________________________________________________________
void NReWeightINuke::SetSystematic(NSyst_t syst, double twk_dial)
{
  switch(syst) {
  case ( kINukeTwkDial_MFP_N ) :
    fMFPDial   =  twk_dial ;  
    break;

  default:
    break;
  }
}
//_______________________________________________________________________________________
void NReWeightINuke::Reset(void)
{
  fMFPCurr  = fMFPDef;

  fMFPDial  = fMFPDef;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightINuke::Reconfigure(void)
{
  fMFPCurr  = fMFPDial;
}
//_______________________________________________________________________________________
double NReWeightINuke::CalcWeight() 
{ 
    bool tweaked = ( fMFPCurr != fMFPDef );  
    if(!tweaked) return 1.0;
    
    // check whether block size limits were exceeded
    if(nucleonfsihist_.nfnvert >= MAXNUCLEONVERT) {
      std::cout << "NReWeightINuke::CalcWeight(): NFnvert hits the limit, not reweighting this event." << std::endl;
      return 1.0;
    }
    if(nucleonfsihist_.nfnstep >= MAXNUCLEONSTEP) {
      std::cout << "NReWeightINuke::CalcWeight(): NFnstep hits the limit, not reweighting this event." << std::endl;
      return 1.0;
    }
    

    // setup nuclear properties
    A = neuttarget_.numatom;
    Z = neuttarget_.numbndp;
    if(Z <= 1) return 1.0;                        // Hydrogen has no FSI, so skip it
    nrsettarg();

    // probabilities
    double p0 = 1.0;        // probability for unchanged parameters
    double p1 = 1.0;        // probability for changed parameters
    
    for(int i=0; i<nucleonfsihist_.nfnvert; i++) {
      int P = nucleonfsihist_.nfiflag[i]%10;
      if(P != 0) {
        cout << "First vertex of track has wrong flag: " << nucleonfsihist_.nfiflag[i] << " but ``P'' should be 0." << endl;
        exit(1);
      }
      while(P != 4) {  // follow track until it leaves the nucleus

        TVector3 start(nucleonfsihist_.nfx[i],nucleonfsihist_.nfy[i],nucleonfsihist_.nfz[i]);   // determine starting point
        TVector3 mom(nucleonfsihist_.nfpx[i],nucleonfsihist_.nfpy[i],nucleonfsihist_.nfpz[i]);  // determine momentum
        double cmom[3];
        mom.GetXYZ(cmom);                       // convert momentum into C array
        double ein = nucleonfsihist_.nfe[i];                    // determine energy
        int ifirststep = nucleonfsihist_.nffirststep[i]-1;      // index of first step in NFecms2 array (convert from Fortran to C convention)
        if(i+1 >= nucleonfsihist_.nfnvert) {                    // safety check array bounds
          cout << "Last vertex in stack with ID=" << i << " has wrong flag: " << nucleonfsihist_.nfiflag[i] << " but should be 4: event weight will be 1." << endl;
          return 1.;
        }
        i++;
        TVector3 stop(nucleonfsihist_.nfx[i],nucleonfsihist_.nfy[i],nucleonfsihist_.nfz[i]);    // determine stopping point
        P = nucleonfsihist_.nfiflag[i]%10;                      // determine interaction type
        int ilaststep = nucleonfsihist_.nffirststep[i]-1;       // index of last step in NFecms2 array (convert from Fortran to C convention)
        int N = (nucleonfsihist_.nfiflag[i]/100)%10;            // determine nucleon type being tracked
        int T = (nucleonfsihist_.nfiflag[i]/10)%10;             // determine target nucleon type being hit
        TVector3 dir = stop-start;              // determine direction
        double dist = dir.Mag();                // distance between start and stop
        dir.SetMag(1.);
        
        // do steps along the track, determine no-interaction probability
        int nsteps = TMath::Nint(dist/step);
        TVector3 pos = start;

        for(int is=0; is<nsteps-1; is++) {
          pos += step*dir;
          double r = pos.Mag();      // radius
          int idx = ifirststep+is+1;
          nrfermi(r, cmom, N);      // determine local nuclear properties (including randomly drawn CMS energy)
          ecms2 = TMath::Abs(nucleonfsihist_.nfecms2[idx]);   // CMS energy
          if(nucleonfsihist_.nfecms2[idx] > 0) {   // proton target
            ichtrgt = 1;
          }
          else {                   // neutron target
            ichtrgt = 0;
          }
          nrnuc(N);                      // determine probability of interaction (filled to member variables)

          p0 *= 1. - ptot;      // no interaction
          p1 *= 1. - ptot_s;    // no interaction (changed by tweak dial)

        }
        
        // final step: probability for particular interaction
        double r = stop.Mag();       // radius
        nrfermi(r, cmom, N);          // determine local nuclear properties (excluding CMS energy, as it is provided by the stack)
        if(ilaststep >= 0 && ilaststep < MAXNUCLEONSTEP) {
          ecms2 = TMath::Abs(nucleonfsihist_.nfecms2[ilaststep]);                // determine CMS energy squared as calculated from NEUT
        }
        else {
          ecms2 = 0;
        }
        ichtrgt = T;
        nrnuc(N);                      // determine probability of interaction (filled to member variables)
      

        if(P >= 1 && P <= 3) {  // any interaction
          p0 *= ptot;
          p1 *= ptot_s;
        }
/*
        if(P == 1 && N == T) {        // elastic p-p
          p0 *= ppel;// * (1.-(ptot-ppel));
          p1 *= ppel_s;// * (1.-(ptot_s-ppel_s));
        }
        else if(P == 1 && N != T) {   // elastic p-n
          p0 *= pnel;// * (1.-(ptot-pnel));
          p1 *= pnel_s;// * (1.-(ptot_s-pnel_s));
        }
        else if(P == 2 && N == T) {   // single pion p-p
          p0 *= ppsp;// * (1.-(ptot-ppsp));
          p1 *= ppsp_s;// * (1.-(ptot_s-ppsp_s));
        }
        else if(P == 2 && N != T) {   // single pion p-n
          p0 *= pnsp;// * (1.-(ptot-pnsp));
          p1 *= pnsp_s;// * (1.-(ptot_s-pnsp_s));
        }
        else if(P == 3 && N == T) {   // double pion p-p
          p0 *= ppdp;// * (1.-(ptot-ppdp));
          p1 *= ppdp_s;// * (1.-(ptot_s-ppdp_s));
        }
        else if(P == 3 && N != T) {   // double pion p-n
          p0 *= pndp;// * (1.-(ptot-pndp));
          p1 *= pndp_s;// * (1.-(ptot_s-pndp_s));
        }
*/        
      }

    }

    if(p0 < 1e-30) return 1.;
    
    return p1/p0;
}
//_______________________________________________________________________________________
double NReWeightINuke::CalcChisq(void)
{
  double chisq = 0.;
/*
  // how to do this for absolute change?

  // For 1-sigma change
  chisq += TMath::Power(fOverallXsecDial, 2.);
*/
  return chisq;
}
//_______________________________________________________________________________________
vector<double> NReWeightINuke::GetCurrParVals(void)
{
  vector<double> parVals;
  // This must be in the same order as ${T2KREWEIGHT}/src/T2KSyst.h
  parVals.push_back(fMFPCurr);   
  return parVals;  
}
//_______________________________________________________________________________________
void NReWeightINuke::nrsettarg() {
      if (Z == 5) { // 11B (estimate)
         NRRMSRAD= 2.42;
         NRC     = 2.00;
         NRCNN   = 2.00;
         NRAF    = 0.51;
         NRWPARM = 0.00;
         NRDFACT = 1.0;
         NRPFSURF = 0.200;
      }
      else if (Z == 6) { // 12C
         NRRMSRAD= 2.455;
         NRC     = 2.355;
         NRCNN   = 2.355;
         NRAF    = 0.5224;
         NRWPARM =-0.149;
         NRDFACT = 1.0;
         NRPFSURF = 0.217;
      }
      else if (Z == 7) { // 14N
         NRRMSRAD= 2.524;
         NRC     = 2.570;
         NRCNN   = 2.570;
         NRAF    = 0.5052;
         NRWPARM =-0.180;
         NRDFACT = 1.0;
         NRPFSURF = 0.221;
      }
      else if (Z == 8) { // 16O
         NRRMSRAD= 2.730;
         NRC     = 2.69;
         NRCNN   = 2.69;
         NRAF    = 0.40961;
         NRWPARM = 0.;
         NRDFACT = 0.9985962;
         NRPFSURF = 0.225;
      }
      else if (Z == 9) { // 19F
         NRRMSRAD= 2.900;
         NRC     = 2.59;
         NRCNN   = 2.59;
         NRAF    = 0.564;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.226;
      }
      else if (Z == 11) { // 23Na (estimate)
         NRRMSRAD= 2.94;
         NRC     = 2.70;
         NRCNN   = 2.70;
         NRAF    = 0.56;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.229;
      }
      else if (Z == 13) { // 27Al
         NRRMSRAD= 3.05;
         NRC     = 2.84;
         NRCNN   = 2.84;
         NRAF    = 0.569;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.233;
      }
      else if (Z == 14) { // 28Si (2pF)
         NRRMSRAD= 3.15;
         NRC     = 3.14;
         NRCNN   = 3.14;
         NRAF    = 0.537;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.234;
      }
      else if (Z == 16) { // 32S
         NRRMSRAD= 3.243;
         NRC     = 3.503;
         NRCNN   = 3.503;
         NRAF    = 0.633;
         NRWPARM = -0.250;
         NRDFACT = 1.0;
         NRPFSURF = 0.235;
      }
      else if (Z == 17) { // 35Cl
         NRRMSRAD= 3.388;
         NRC     = 3.476;
         NRCNN   = 3.476;
         NRAF    = 0.599;
         NRWPARM =-0.10;
         NRDFACT = 1.0;
         NRPFSURF = 0.236;
      }
      else if (Z == 18) { // 40Ar
         NRRMSRAD= 3.48;
         NRC     = 3.73;
         NRCNN   = 3.73;
         NRAF    = 0.62;
         NRWPARM =-0.19;
         NRDFACT = 1.0;
         NRPFSURF = 0.237;
      }
      else if (Z == 20) { // 40Ca
         NRRMSRAD= 3.482;
         NRC     = 3.766;
         NRCNN   = 3.766;
         NRAF    = 0.586;
         NRWPARM =-0.161;
         NRDFACT = 1.0;
         NRPFSURF = 0.241;
      }
      else if (Z == 22) { // 48Ti
         NRRMSRAD= 3.59;
         NRC     = 3.75;
         NRCNN   = 3.75;
         NRAF    = 0.567;
         NRWPARM = 0.0;
         NRDFACT = 1.0;
         NRPFSURF = 0.243;
      }
      else if (Z == 26) { // 56Fe
         NRRMSRAD= 3.787;
         NRC     = 3.971;
         NRCNN   = 3.971;
         NRAF    = 0.5935;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.250;
      }
      else if (Z == 27) { // 59Co
         NRRMSRAD= 3.80;
         NRC     = 4.08;
         NRCNN   = 4.08;
         NRAF    = 0.569;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.250;
      }
      else if (Z == 28) { // 58Ni
         NRRMSRAD= 3.823;
         NRC     = 4.14 ;
         NRCNN   = 4.14 ;
         NRAF    = 0.56 ;
         NRWPARM = 0.   ;
         NRDFACT = 1.0  ;
         NRPFSURF = 0.250;
      }
      else if (Z == 29) { // 63Cu
         NRRMSRAD= 3.925;
         NRC     = 4.214;
         NRCNN   = 4.214;
         NRAF    = 0.586;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.250;
      }
      else if (Z == 30) { // 64Zn
         NRRMSRAD= 3.965;
         NRC     = 4.285;
         NRCNN   = 4.285;
         NRAF    = 0.584;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.250;
      }
      else if (Z == 38) { // 88Sr
         NRRMSRAD= 4.17;
         NRC     = 4.83;
         NRCNN   = 4.83;
         NRAF    = 0.496;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.250;
      }
      else if (Z == 40) { // 90Zr (Estimate)
         NRRMSRAD= 4.274 ;
         NRC     = 4.8   ;
         NRCNN   = 4.8   ;
         NRAF    = 0.57  ;
         NRWPARM = 0.    ;
         NRDFACT = 1.0   ;
         NRPFSURF = 0.250;
      }
      else if (Z == 41) { // 93Nb
         NRRMSRAD= 4.31;  
         NRC     = 4.87 ; 
         NRCNN   = 4.87  ;
         NRAF    = 0.573 ;
         NRWPARM = 0.    ;
         NRDFACT = 1.0 ;
         NRPFSURF = 0.250;
      }
      else if (Z == 45) { // 103Rh (Estimate from 110Cd)
         NRRMSRAD= 4.6   ;
         NRC     = 5.3   ;
         NRCNN   = 5.3   ;
         NRAF    = 0.55  ;
         NRWPARM = 0.    ;
         NRDFACT = 1.0  ;
         NRPFSURF = 0.250;
      }
      else if (Z == 50) { // 120Sn
         NRRMSRAD= 4.64;  
         NRC     = 5.315; 
         NRCNN   = 5.315 ;
         NRAF    = 0.576 ;
         NRWPARM = 0.    ;
         NRDFACT = 1.0   ;
         NRPFSURF = 0.250;
      }
      else if (Z == 67) { // 165Ho
         NRRMSRAD= 5.23; 
         NRC     = 6.18 ;
         NRCNN   = 6.18 ;
         NRAF    = 0.57 ;
         NRWPARM = 0.   ;
         NRDFACT = 1.0 ;
         NRPFSURF = 0.250;
      }
      else if (Z == 79) { // 197Au
         NRRMSRAD= 5.33  ;
         NRC     = 6.38  ;
         NRCNN   = 6.38  ;
         NRAF    = 0.535 ;
         NRWPARM = 0.    ;
         NRDFACT = 1.0  ;
         NRPFSURF = 0.250;
      }
      else if (Z == 82) { // 208Pb
         NRRMSRAD= 5.521;
         NRC     = 6.624;
         NRCNN   = 6.624;
         NRAF    = 0.549;
         NRWPARM = 0.;
         NRDFACT = 1.;
         NRPFSURF = 0.245;
      }
      else if (Z == 83) { // 209Bi
         NRRMSRAD= 5.51;
         NRC     = 6.75;
         NRCNN   = 6.75;
         NRAF    = 0.468;
         NRWPARM = 0.;
         NRDFACT = 1.0;
         NRPFSURF = 0.245;
      }
      else {
         cout << "This target (A = " << A << ", Z = " << Z << ") is not supported yet, see nrsettarg()." << endl;
         exit(1);
      }



      // this code is from nrinit.F
      
      const double RMSRADOXY = 2.730;

      double rbin = 0.2;
      double rmin = 0.0;
      double rmax = 6.*(NRRMSRAD/RMSRADOXY);
      //double pfermi=NRPFSURF*1000.; // NRPFSURF is GeV unit
	  nbin = TMath::Nint((rmax-rmin)/rbin);

      double probm[200];

	  // call nrrodis(RMIN=rmin,DR=rbin,NBIN=nbin,PROBM=probm,ROT=rhotab)  ==>
	  double R = rmin-rbin;
	  for(int i=0; i<nbin; i++) {
		 R += rbin;
		 rhotab[i] = nrROXY(R);
		 probm[i] = 4./3. * 3.1416 * (pow(R+rbin,3)-pow(R,3));
	  }
	  for(int i=0; i<nbin-1; i++) {
		 probm[i] = probm[i]*(rhotab[i]+rhotab[i+1])/2.;
	  }
	  probm[nbin-1] *= rhotab[nbin-1];
      // end nrodis

      double probs=0.0;
	  for(int i=0; i<nbin; i++) probs += probm[i];
	  pnorm = A/probs;
	  for(int i=0; i<nbin; i++) rhotab[i] *= pnorm;
}
//_______________________________________________________________________________________
double NReWeightINuke::nrROXY(double r) {

	  double w,c,z, gr,ro,dro;
	  const double d=.102, g=2.76, e=.35;

	  w = NRWPARM;
	  c = NRC;
	  z = NRAF;
	  ro = (1+w*pow(r/c,2))/(1+exp((r-c)/z));
	  if(r == 0) {
		 dro = 1;
	  }
	  else {
		 gr = g*r;
		 dro = sin(gr)/gr;
	  }
	  dro *= d*exp(-pow(e*r,2));
	  double result = ro-dro;
	  if(result >= 1.) result = 0.9999999;
	  if(result <= 0.) result = 1.e-7;
	  return result;

}
//_______________________________________________________________________________________
void NReWeightINuke::nrfermi(double r, double *pin, int N) {

      const double toten2 = pow(931.,2);        // atomic mass (in u) to MeV conversion factor

      prat = (double)Z/(double)A;
      ichtrgt = -1;

      double rbin = 0.2;
      double rmin = 0.;
      double pfermi = NRPFSURF*1000.; // NRPFSURF is GeV unit
      if(abs(nework_.modene) == 2) {    // For MEC events fermi surface momentum is function of radial distance r (Asmita R., 22-07-2013)
        cout << "HIER EFFRMGAS" << endl;
        float rfloat = r;
      	pfermi=effrmgas_(0,0,&rfloat);
      }


//    NUCLEAR DENSITY AT THIS VALUE OF R
//    INCLUDE 0.1 FACTOR is TO CHANGE FROM MB TO FMSQ (cms is in fm!)
      int irlo,irhi;
      double rrem;
      nrhis(r-rmin,rbin,nbin,irlo,irhi,rrem);
      rhon = (rhotab[irhi]*rrem + rhotab[irlo]*(1.0-rrem))*0.1;
      //rhonf = rhon/pnorm;

      unucl = sqrt(toten2-pow(pfermi,2));

}
//_______________________________________________________________________________________
void NReWeightINuke::nrhis(double x, double rbin, int num, int &ilo, int &ihi, double &erem) {

	  ilo = TMath::Floor(x/rbin+0.01);
	  erem = (x-(double)ilo*rbin)/rbin;
	  ihi = ilo+1;
	  if(ilo >= num-1) {
		 ilo = num-1;
		 ihi = num-1;
		 erem = 1.;
	  }
	  if(ilo < 0) {
		 ilo = 0;
		 ihi = 0;
		 erem = 1.;
	  }

}
//_______________________________________________________________________________________
void NReWeightINuke::nrnuc(int ichint) {

      const double fac=1e-2, amn=939., amn2=amn*amn, ampi=139.6;
      
      ptot = 0;
      ppel = 0;
      ppsp = 0;
      ppdp = 0;
      pnel = 0;
      pnsp = 0;
      pndp = 0;
      ptot_s = 0;
      ppel_s = 0;
      ppsp_s = 0;
      ppdp_s = 0;
      pnel_s = 0;
      pnsp_s = 0;
      pndp_s = 0;
      if(ecms2 <= 4.*amn2) return;     // kinematically forbidden, no interaction
      
      double eineq = (ecms2 - pow(unucl,2) - pow(amn,2))/(2.*unucl);
      double ekin = eineq - amn;
      
      int ihi,ilo,ihi1,ilo1,ihi2,ilo2;
      double erem,erem1,erem2;
      
      double prat;
      if(ichtrgt < 0) {
        if(ichint == 1) {
          prat = prat;
        }
        else {
          prat = 1.-prat;
        }
      }
      else {
        if(ichint == ichtrgt) {
          prat = 1.;
        }
        else {
          prat = 0.;
        }
      }
      
      //cout << "HIER " << ekin << " " << ecms2 << " " << rhon << " " << prat << endl;

      if(ekin < 360. || ecms2 < pow(amn+amn+ampi,2) ) {            // below pion-production threshold
         nrhis(ekin,20.,176,ilo,ihi,erem);
         // probability for p-p and p-n interactions
         ppel = xepp[ihi]*erem + xepp[ilo]*(1-erem);
         pnel = xepn[ihi]*erem + xepn[ilo]*(1-erem);
      }
      else if(ekin < 920. || ecms2 < pow(amn+amn+ampi+ampi,2) ) {  // below double-pion-production threshold
         nrhis(ekin,20.,176,ilo,ihi,erem);
         nrhis(ekin-360.,20.,158,ilo1,ihi1,erem1);
         // probability for p-p and p-n interactions
         ppel = xepp[ihi] *erem  + xepp[ilo] *(1-erem);
         pnel = xepn[ihi] *erem  + xepn[ilo] *(1-erem);
         ppsp = xspp[ihi1]*erem1 + xspp[ilo1]*(1-erem1);
         pnsp = xspn[ihi1]*erem1 + xspn[ilo1]*(1-erem1);
      }
      else {
         nrhis(ekin,20.,176,ilo,ihi,erem);
         nrhis(ekin-360.,20.,158,ilo1,ihi1,erem1);
         nrhis(ekin-920.,20.,130,ilo2,ihi2,erem2);
         // probability for p-p and p-n interactions
         ppel = xepp[ihi] *erem  + xepp[ilo] *(1-erem);
         pnel = xepn[ihi] *erem  + xepn[ilo] *(1-erem);
         ppsp = xspp[ihi1]*erem1 + xspp[ilo1]*(1-erem1);
         pnsp = xspn[ihi1]*erem1 + xspn[ilo1]*(1-erem1);
         ppdp = xdpp[ihi2]*erem2 + xdpp[ilo2]*(1-erem2);
         pndp = xdpn[ihi2]*erem2 + xdpn[ilo2]*(1-erem2);
      }
      
      //cout << ppel << " " << ppsp << " " << ppdp << " " << pnel << " " << pnsp << " " << pndp << endl;

      ppel *= fac*rhon;
      ppsp *= fac*rhon;
      ppdp *= fac*rhon;
      pnel *= fac*rhon;
      pnsp *= fac*rhon;
      pndp *= fac*rhon;

      ptot = (ppel+ppsp+ppdp)*prat + (pnel+pnsp+pndp)*(1.-prat);

      if(ptot < 0 || ppel < 0 || ppsp < 0 || ppdp < 0 || pnel < 0 || pnsp < 0 || pndp < 0) {
        cout << "negative probability!" << endl;
        return;
      }

      double factor = fMFPDial;
      ptot_s = 1. - exp(-step*ptot*factor);
      ppel_s = 1. - exp(-step*ppel*factor);
      ppsp_s = 1. - exp(-step*ppsp*factor);
      ppdp_s = 1. - exp(-step*ppdp*factor);
      pnel_s = 1. - exp(-step*pnel*factor);
      pnsp_s = 1. - exp(-step*pnsp*factor);
      pndp_s = 1. - exp(-step*pndp*factor);

      
      ptot = 1. - exp(-step*ptot);
      ppel = 1. - exp(-step*ppel);
      ppsp = 1. - exp(-step*ppsp);
      ppdp = 1. - exp(-step*ppdp);
      pnel = 1. - exp(-step*pnel);
      pnsp = 1. - exp(-step*pnsp);
      pndp = 1. - exp(-step*pndp);
}
