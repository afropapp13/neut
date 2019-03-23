************************************************************************
*     ---------------------------------
      REAL FUNCTION DCOHMXNEREIN(E)
*     ---------------------------------
*
*     (Purpose)
C     THIS SUBROUTINE IS FOR 'COHVCT' THAT IS COHERENT PI0 PRODUCTION
C     PROGRAM. HERE,WE GIVES THE MAX VALUE OF THE FUNCTION 'DCOH'.
*
*     (Input)
C       E     : INCIDENT NEUTRINO ENERGY ( GEV )      
C
*     (Output)
C       Maximum value of function DCOHNE
*      
*     (Creation Date and Author)
*       1994.??.?? ; K.Okumura
*       1996.??.?? ; K.Kaneyuki   COHCRS -> COHCRSNE
*       1996.??.?? ; K.Okumura    add cross section from 100GeV to 1000GeV
*       2007.01.08 ; G.Mitsuka    Rename from DCOHMXNE to DCOHMXNEREIN
*       2007.11.05 ; G.Mitsuka    support target nucleus besides 16O
*
*     (Comment)
*
*
************************************************************************

#include "necard.h"
C-Mod 97.03.22 K.O  extension to 1000GeV
C      REAL SIGMAX(100),SIG1(50),SIG2(50)
      REAL SIGMAX(109),SIG1(50),SIG2(50),SIG3(9)
      EQUIVALENCE (SIGMAX( 1),SIG1(1))
      EQUIVALENCE (SIGMAX(51),SIG2(1))
C-Add 97.03.22 K.O
      EQUIVALENCE (SIGMAX(101),SIG3(1))
      
      DATA SIG1/.91041E+02,.29533E+03,.50238E+03,.11129E+04,.19779E+04,
     &          .31529E+04,.51310E+04,.77925E+04,.10388E+05,.13143E+05,
     &          .16261E+05,.20892E+05,.26309E+05,.32594E+05,.39814E+05,
     &          .47954E+05,.57206E+05,.67541E+05,.78972E+05,.91552E+05,
     &          .10550E+06,.12094E+06,.13771E+06,.15576E+06,.17560E+06,
     &          .19665E+06,.21977E+06,.24442E+06,.27051E+06,.29916E+06,
     &          .32941E+06,.36123E+06,.39446E+06,.43133E+06,.46916E+06,
     &          .50981E+06,.55188E+06,.59736E+06,.64482E+06,.69468E+06,
     &          .74272E+06,.79984E+06,.85710E+06,.91888E+06,.98136E+06,
     &          .10450E+07,.11140E+07,.11844E+07,.12573E+07,.13307E+07/
      DATA SIG2/.14226E+07,.15033E+07,.15822E+07,.16745E+07,.17683E+07,
     &          .18697E+07,.19740E+07,.20710E+07,.21676E+07,.22708E+07,
     &          .24169E+07,.25300E+07,.26457E+07,.27778E+07,.29023E+07,
     &          .30574E+07,.31822E+07,.33210E+07,.34842E+07,.36061E+07,
     &          .37860E+07,.39241E+07,.40978E+07,.42526E+07,.44150E+07,
     &          .45855E+07,.48086E+07,.49541E+07,.51745E+07,.53100E+07,
     &          .55221E+07,.57462E+07,.59881E+07,.61494E+07,.64161E+07,
     &          .66489E+07,.68628E+07,.70893E+07,.73916E+07,.76727E+07,
     &          .78610E+07,.81928E+07,.83352E+07,.87537E+07,.89878E+07,
     &          .93418E+07,.95590E+07,.97356E+07,.10161E+08,.10457E+08/
      DATA SIG3/.90738E+08,.33152E+09,.88628E+09,.32364E+10,.55872E+10,
     &          .89624E+10,.13347E+11,.18967E+11,.25976E+11/


      IE = INT(AINT(E))+1
C-Add 97.03.16 Y.H
      IF (IE.gt.100) THEN
C-Mod 97.03.22 K.O
C         write(*,*) "E is OUT OF RANGE...(E=",E,")"
C         write(*,*) "(Approx by E=100.)"
         IE = 100 + AINT(E/100.)
         IF (IE.GT.109) THEN
            write(*,*) "E is OUT OF RANGE...(E=",E,")"
            write(*,*) "(Approx by E=1000.)"
            IE = 109
         ENDIF
      ENDIF
*     scaling factor is applied for other target than 16O
      DCOHMXNEREIN = SIGMAX(IE)*(float(NUMATOM)/16.)**(1./3.)

      RETURN
      END
