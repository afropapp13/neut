      program main
      INTEGER*4 IP,MODE
      REAL*4 ENU,DIRNU(3)
      INTEGER*4 IERR
      
#include "neutinc.h"

      IBOUND = 0 
      NUMATOM = 12
      
      PFMAX = 0.225
      
      DIRNU(1) = 0.
      DIRNU(2) = 0.
      DIRNU(3) = 1.
      
      crstblpath = "./HT2p2h/"
      
      DO I = 1, 1000

         E = 1.
         
        PRINT*," CROSS-SECTION ",FN_2P2H_HT(E,14,MODE) 
        
        AA  = NE_2p2h_HT(14,MODE,E,DIRNU,IERR)
           
        PRINT*,MODENE
        PRINT*,NUMNE
        DO K = 1 , NUMNE
           PRINT*,"---"
           PRINT*,IPNE(K),IORGNE(K),IFLGNE(K),ICRNNE(I) 
           PRINT*,"(",PNE(1,K),",",PNE(2,K),",",PNE(3,K),")"
        ENDDO
        
        PRINT*,"---------------------------------------"
        
      ENDDO
      

      return
      end 
