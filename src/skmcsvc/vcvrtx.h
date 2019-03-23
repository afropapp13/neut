C-------------------------------------------------------------
C
C     INCLUDE FILE FOR VECTOR GENERATOR ( vcvrtx.h )
C
C     VERTEX INFORMATION
C
C-------------------------------------------------------------
      INTEGER MAXVX
      PARAMETER(MAXVX=100)
C
C     NVTXVC      : # OF VERTEX 
C     PVTXVC(3,I) : POSITION OF VERTEX
C     IFLVVC(I)   : KIND OF VERTEX
C                   0 : DETERMINED LATER PROCEDURE
C                   1 : PRIMARY POSITION
C                   2 : DECAY
C                   3 : STOP
C     IPARVC(I)   : PARENT PARTICLE
C     TIMVVC(I)   : TIME OF VERTEX
C
      INTEGER NVTXVC, IFLVVC, IPARVC
      REAL    PVTXVC, TIMVVC
      COMMON /VCVRTX/ NVTXVC,PVTXVC(3,MAXVX),IFLVVC(MAXVX),
     $    IPARVC(MAXVX),TIMVVC(MAXVX)


