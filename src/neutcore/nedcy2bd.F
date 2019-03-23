**********************************************************************
*     -----------------------------------
      SUBROUTINE NEDCY2BD (AM0,AM1,AM2,P)
*     -----------------------------------
*     
*     ( purpose )
*     GET MOMENTUM OF DECAYED PARTICLE IN C.M.S. ONLY 2 BODY DECAY IS
*     CONSIDERED.
*
*     ( input ) 
*       AM0 : MASS OF DECAYED PARTICLE
*       AM1 :    "    DECAY PRODUCT 1
*       AM2 :    "          "       2
*
*     ( output ) 
*       P   : MOMENTUM
*
*     ( creation date and author )
*       1985.10.19 ; First version by T.KAJITA
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    DCY2BD -> NEDCY2BD
*
**********************************************************************
      IF (AM0 .GE. (AM1+AM2) ) THEN
         P=SQRT((AM0**2-(AM1+AM2)**2)*(AM0**2-(AM1-AM2)**2)/
     &     (4.*AM0**2)  )
      ELSE
         WRITE(6,600)AM0,AM1,AM2
 600     FORMAT(' ','** ERROR ** ',F7.1,'<',F7.1,'+',F7.1)
      ENDIF

      RETURN
      END
