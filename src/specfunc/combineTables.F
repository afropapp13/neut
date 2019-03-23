******Subroutine combineTables combines individual nucleon tables 
******into one table for reading in by fnelspau_sf
******takes input argument either CC or NC
******for NC, only nue tables are read in, as all flavours have the same cross-section
      
      subroutine combineTables(mode)
      implicit none
      character(len=2) mode
      character(len=80) outputFileName
      character(len=80) inputFileName

c     ----- logical needed to check whether files exist
      logical Iexist

c     ----- integers corresponding to in and out files
      integer outI, inI

      integer pdg, i, j, k, nBins
      character(len=1) sgn
      double precision EE(210)
      double precision EM(210)
      double precision ET(210)
      double precision CRSLE(9,210)
      double precision CRSRE(9,210)
      double precision CRSLM(9,210)
      double precision CRSRM(9,210)
      double precision CRSLT(9,210)
      double precision CRSRT(9,210)
      double precision CRSNCNU(2,9,210) ! nc neutrino
      double precision CRSNCAN(2,9,210) ! nc antineutrino

      integer table

      data outI/80/
      data inI/90/

      do table = 1, 2 ! Loop through total and maxDiff

         if (table .eq. 1) then
            outputFileName = "qelSfData/totXsec/ALL_1000000000_"
         else 
            outputFileName = "qelSfData/maxDiff/ALL_1000000000_"
         endif

         if (mode.EQ. "CC") then
            outputFileName = outputFileName(1:33) // "cc.csv"
         else if (mode.EQ. "NC")then
            outputFileName = outputFileName(1:33) // "nc.csv"
         endif


*     * Read in tables one-by-one
*     * if no table exists, set all xsec values to 0.0

*     loop through neutrino types, and nuclei, for cc and nc
*     first - nu_mu, C12
         do i = 1, 6 ! neutrino type
            if (i .EQ. 1) then
               pdg = 12
            else if (i .EQ. 2) then
               pdg = -12
            else if (i .EQ. 3) then
               pdg = 14
            else if (i .EQ. 4) then
               pdg = -14
            else if (i .EQ. 5) then
               pdg = 16
            else if (i .EQ. 6) then
               pdg = -16
            endif

            if (pdg .LT. 0) then
               sgn = "-"
            else
               sgn = "+"
            endif
            do j = 1, 9 ! nucleus loop
               

               if (table .eq. 1) then
                  write(inputFileName,"(A18,A1,I2)")
     &                 "qelSfData/totXsec/",sgn,abs(pdg)
               else
                  write(inputFileName,"(A18,A1,I2)")
     &                 "qelSfData/maxDiff/",sgn,abs(pdg)
               endif
c     write (*,*) inputFileName
c     inputFileName = inputFileName(1:18) // sgn
c     write (*,*) inputFileName
c     inputFileName = inputFileName(1:19) // char(abs(pdg))
               
               if (j .EQ. 1) then ! C12
                  inputFileName = inputFileName(1:21) 
     $                 // "_1000060120_"
               else if (j .EQ. 2) then ! O16
                  inputFileName = inputFileName(1:21) 
     $                 // "_1000080160_"
               else if (j .EQ. 3) then ! Fe56
                  inputFileName = inputFileName(1:21) 
     $                 // "_1000260560_"
               else if (j .GT. 3) then ! no more nuclei
                  inputFileName = inputFileName(1:21) 
     $                 // "_1000000000_"
                  if (abs(i) .LT. 15) then
                     do k = 1, 210
                        CRSLE(j,k)=0
                        CRSRE(j,k)=0
                        CRSLM(j,k)=0
                        CRSRM(j,k)=0
                     enddo
                  else
                     do k = 1, 125
                        CRSLT(j,k)=0
                        CRSRT(j,k)=0
                     enddo
                  endif
               endif

               if (mode.EQ. "CC") then
                  inputFileName = inputFileName(1:33) // "cc.csv"
               else if (mode.EQ. "NC")then
                  inputFileName = inputFileName(1:33) // "nc.csv"
               endif
               
               inquire(file=inputFileName, exist=Iexist)
               if (Iexist .EQV. .TRUE.) then
                  open(inI, file=inputFileName, status="old")
                  read(inI, *) nBins
                  do k = 1, nBins ! loop through energies
c     ----- fill the correct neutrino/nucleus array 
                     if (mode .EQ. "CC") then
                        if (i .EQ. 1) then
                           read(inI, *) EE(k), CRSLE(j, k)
                        else if (i .EQ. 2) then
                           read(inI, *) EE(k), CRSRE(j, k)
                        else if (i .EQ. 3) then
                           read(inI, *) EM(k), CRSLM(j, k)
                        else if (i .EQ. 4) then
                           read(inI, *) EM(k), CRSRM(j, k)
                        else if (i .EQ. 5) then
                           read(inI, *) ET(k), CRSLT(j, k)
                        else if (i .EQ. 6) then
                           read(inI, *) ET(k), CRSRT(j, k)
                        endif
                     else if (mode .EQ. "NC") then
                        if (i .EQ. 1) then
                           read(inI, *) EE(k), CRSNCNU(1,j,k),
     $                          CRSNCNU(2,j,k)
                        else if (i.EQ.2) then
                           read(inI, *) EE(k), CRSNCAN(1,j,k),
     $                          CRSNCAN(2,j,k)
                        endif
                     endif
                  enddo ! end loop through individual table
                  close (inI)
               else ! file doesn't exist: set to zero
                  if (j .LT. 4 .AND. mode .NE. "NC")then
                     if (table .eq. 1) then
                        write (*,"(A24,I3,A2,A2,A9,I1)") 
     &                       "no totXsec File for pdg ", pdg,", 
     $                       ", mode, ", nucleus ",j
                     else
                        write (*,"(A24,I3,A2,A2,A9,I1)") 
     &                       "no maxDiff File for pdg ", pdg,", 
     $                       ", mode, ", nucleus ",j
                     endif
                     write(*,*) "1-C12, 2-O16, 3-Fe56"
                  endif
                  if (mode .EQ."CC")then
                     if (pdg .EQ. 12) then
                        do k = 1, 210
                           CRSLE(j,k)=0
                        enddo
                     else if (pdg .EQ. -12) then
                        do k = 1, 210
                           CRSRE(j,k)=0
                        enddo
                     else if (pdg .EQ. 14) then
                        do k = 1,210
                           CRSLM(j,k)=0
                        enddo
                     else if (pdg .EQ. -14) then
                        do k = 1, 210
                           CRSRM(j,k)=0
                        enddo
                     else if (pdg .EQ. 16) then
                        do k = 1, 210
                           CRSLT(j,k)=0
                        enddo
                     else
                        do k = 1,210
                           CRSRT(j,k)=0
                        enddo
                     endif
                  endif ! end mode
               endif ! end file exists check

c     ----- ending loop over nuclei
            enddo

c     ----- ending loop over PDGs
         enddo

*******now, write new File with values read in

c     ---------- pdg, energy, xsecs by nucleus (9 in total)
 7       format(I3,F12.6,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,
     &        F12.4,F12.4,F12.4)

 14      format(I3,F12.6,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,
     &        F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,
     &        F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)
         
         open(outI, file=outputFileName, status="new")
         if (mode .EQ. "CC")then
            do i = 1, 1090
               if (i.le.210) then ! nue
                  pdg = 12
                  write(outI, 7) pdg, EE(i),
     &                 CRSLE(1,i),CRSLE(2,i),CRSLE(3,i),CRSLE(4,i),
     &                 CRSLE(5,i),CRSLE(6,i),CRSLE(7,i),CRSLE(8,i),
     $                 CRSLE(9,i)
                  
               else if (i.gt.210 .and. i.le.420) then ! nue-bar
                  pdg = -12
                  write(outI, 7) pdg, EE(i-210),
     &                 CRSRE(1,i-210),CRSRE(2,i-210),CRSRE(3,i-210),
     $                 CRSRE(4,i-210),
     &                 CRSRE(5,i-210),CRSRE(6,i-210),CRSRE(7,i-210),
     $                 CRSRE(8,i-210),
     &                 CRSRE(9,i-210)
                  
               else if (i.gt.420 .and. i .le.630)then ! numu
                  pdg = 14
                  write(outI, 7) pdg, EM(i-420),
     &                 CRSLM(1,i-420),CRSLM(2,i-420),CRSLM(3,i-420),
     $                 CRSLM(4,i-420),
     &                 CRSLM(5,i-420),CRSLM(6,i-420),CRSLM(7,i-420),
     $                 CRSLM(8,i-420),
     &                 CRSLM(9,i-420)

               else if (i.gt.630 .and. i.le.840) then ! numu-bar
                  pdg = -14
                  write(outI, 7) pdg, EM(i-630),
     &                 CRSRM(1,i-630),CRSRM(2,i-630),CRSRM(3,i-630),
     $                 CRSRM(4,i-630),
     &                 CRSRM(5,i-630),CRSRM(6,i-630),CRSRM(7,i-630),
     $                 CRSRM(8,i-630),
     &                 CRSRM(9,i-630)

               else if (i.gt.840 .and. i.le.965) then ! nutau
                  pdg = 16
                  write(outI, 7) pdg, ET(i-840),
     &                 CRSLT(1,i-840),CRSLT(2,i-840),CRSLT(3,i-840),
     $                 CRSLT(4,i-840),
     &                 CRSLT(5,i-840),CRSLT(6,i-840),CRSLT(7,i-840),
     $                 CRSLT(8,i-840),
     &                 CRSLT(9,i-840)
                  
               else if (i.gt.965 .and. i.le.1090) then ! nutau-bar
                  pdg = -16
                  write(outI, 7) pdg, ET(i-965),
     &                 CRSRT(1,i-965),CRSRT(2,i-965),CRSRT(3,i-965),
     $                 CRSRT(4,i-965),
     &                 CRSRT(5,i-965),CRSRT(6,i-965),CRSRT(7,i-965),
     $                 CRSRT(8,i-965),
     &                 CRSRT(9,i-965)
               endif
            enddo
         else if (mode .EQ."NC")then
            do i = 1, 800
               if (i.LE.200) then
                  pdg = 12
                  write (outI,14) pdg, EE(i),
     &                 CRSNCNU(1,1,i), CRSNCNU(1,2,i),CRSNCNU(1,3,i),
     $                 CRSNCNU(1,4,i),
     &                 CRSNCNU(1,5,i), CRSNCNU(1,6,i),CRSNCNU(1,7,i),
     $                 CRSNCNU(1,8,i),
     &                 CRSNCNU(1,9,i)
               else if (i.GT.200 .AND. i.LE.400) then
                  pdg = 12
                  write (outI,14) pdg, EE(i-200),
     &                 CRSNCNU(2,1,i-200),
     &                 CRSNCNU(2,2,i-200),
     &                 CRSNCNU(2,3,i-200),
     &                 CRSNCNU(2,4,i-200),
     &                 CRSNCNU(2,5,i-200),
     &                 CRSNCNU(2,6,i-200),
     &                 CRSNCNU(2,7,i-200),
     &                 CRSNCNU(2,8,i-200),
     &                 CRSNCNU(2,9,i-200)
               else if (i .GT. 400 .AND. i .LE. 600) then
                  pdg = -12
                  write (outI,14) pdg, EE(i-400),
     &                 CRSNCAN(1,1,i-400),
     &                 CRSNCAN(1,2,i-400),
     &                 CRSNCAN(1,3,i-400),
     &                 CRSNCAN(1,4,i-400),
     &                 CRSNCAN(1,5,i-400),
     &                 CRSNCAN(1,6,i-400),
     &                 CRSNCAN(1,7,i-400),
     &                 CRSNCAN(1,8,i-400),
     &                 CRSNCAN(1,9,i-400)
               else if (i .GT. 600 .AND. i .LE. 800) then
                  pdg = -12
                  write (outI,14) pdg, EE(i-600),
     &                 CRSNCAN(2,1,i-600),
     &                 CRSNCAN(2,2,i-600),
     &                 CRSNCAN(2,3,i-600),
     &                 CRSNCAN(2,4,i-600),
     &                 CRSNCAN(2,5,i-600),
     &                 CRSNCAN(2,6,i-600),
     &                 CRSNCAN(2,7,i-600),
     &                 CRSNCAN(2,8,i-600),
     &                 CRSNCAN(2,9,i-600)
               endif ! end neutrino type
            enddo
         endif ! end NC/CC

c     --- close the output file
         close (outI)
      enddo ! end loop through total and maxdiff

      return
      end
      
