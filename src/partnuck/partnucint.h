C  Holds what type of interaction occurred in the last call to
C  PARTNUC.
C
C  pninttype is:
C       0:      No interaction
C       1:      Pion Absorbtion
C       2:      Pion production
C       3:      Pion scatter
C       4:      Pion charge exchange
C
C So far this only supports interaction in nucpion.
C
        integer pninttype       ! last interaction type
        integer pnnumber(5)     ! number of interactions
        integer pninter(100,5)  ! The list of interactions at each step
        common /partnucint/ pninttype,pnnumber,pninter
