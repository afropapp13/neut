*----------------------------------------------------------------------
*             MACRO ZBSSYS
*       << Include file for system constants >>
*----------------------------------------------------------------------
*
* <<  Position of the pointer to the Bank Directory Table >>
      INTEGER IPDIRT
      PARAMETER (IPDIRT=1)
* <<  Position of the pointer to the Event Control Bank >>
      INTEGER IPEVC
      PARAMETER (IPEVC=2)

* <<  Length of the header of Bank Directory Table >>
      INTEGER LDIRHD
      PARAMETER (LDIRHD=10)
* <<  Position of parameters in the header >>
*     Current size of the Table
      INTEGER INBTAB
      PARAMETER (INBTAB= 1)
*     No. of currently registered banks
      INTEGER INBANK
      PARAMETER (INBANK= 2)
*     Logical unit number for print out
      INTEGER INLUN
      PARAMETER (INLUN = 3)
* <<  No. of links of the Bank Directory Table >>
      INTEGER LLNKHD
      PARAMETER (LLNKHD=14)
*     Pointer table
      INTEGER ILBPNT
      PARAMETER (ILBPNT= 9)
*     Existing flag
      INTEGER ILBEXF
      PARAMETER (ILBEXF=10)
*     Output flag
      INTEGER ILBIOF
      PARAMETER (ILBIOF=11)
*     Sort table
      INTEGER ILBSRT
      PARAMETER (ILBSRT=12)
*     I/O Format word
      INTEGER ILBIOX
      PARAMETER (ILBIOX=13)
*     event header flag word
      INTEGER ILBEVH
      PARAMETER (ILBEVH=14)

*     Initial size of the Bank Directory Table
      INTEGER LBDIR1
      PARAMETER (LBDIR1=50)
*     Initial size of the Data Record
      INTEGER LDTRC0
      PARAMETER (LDTRC0=500)
*     Initial factor for calculating secondary size of the Data Record
      REAL FDTRC1
      PARAMETER (FDTRC1=0.4)
*     Initial no. of segments
      INTEGER NSEGM0
      PARAMETER (NSEGM0=50)
*     Initial factor for calculating secondary no of segments
      REAL FSEGM1
      PARAMETER (FSEGM1=0.4)

***** Parameters for Event Control Bank
*  << Link pointer positions in the Event Control Bank >>
*     Link to linear structure for the Event Mother banks
      INTEGER ILEMLS
      PARAMETER (ILEMLS=1)
*     Link to the Event Mother bank of the current event
      INTEGER ILEMCR
      PARAMETER (ILEMCR=2)

***** Parameters for Event Mother Bank
*  << Length of the data of Event Mother Bank >>
      INTEGER LEMBDT
      PARAMETER (LEMBDT=1)
*  << Position of parameter in the data of Event Mother Bank >>
*     Flag for event mark
      INTEGER IEVMRK
      PARAMETER (IEVMRK=1)
*  << Link pointer positions in the data of Event Mother Bank >>
*     Link to linear structure for the Event Header banks
      INTEGER ILEHLS
      PARAMETER (ILEHLS=1)
*     Link to the Event Header bank of the current sub-event
      INTEGER ILEHCR
      PARAMETER (ILEHCR=2)

***** Parameters for Event Header Bank
*  << Link pointer positions in the data of Event Header Bank >>
*     Link to linear structure for the Bank Index Records
      INTEGER ILBILS
      PARAMETER (ILBILS=1)

***** Parameters for Bank Index Record
*  << Length of the header of Bank Index Record >>
      INTEGER LINDHD
      PARAMETER (LINDHD=21)
*  << Position of parameters in the header of Bank Index Record >>
*     Type of the bank
      INTEGER IBKTYP
      PARAMETER (IBKTYP= 9)
*     Current total length of the data record
      INTEGER ILDTRC
      PARAMETER (ILDTRC=10)
*     Currently filled length of the data record
      INTEGER ILFILL
      PARAMETER (ILFILL=11)
*     No. of expansion of the data record
      INTEGER INEXDT
      PARAMETER (INEXDT=12)
*     Current secondary size of the data record
      INTEGER ILDTR1
      PARAMETER (ILDTR1=13)
*     Current no. of segments
      INTEGER INSEGM
      PARAMETER (INSEGM=14)
*     Maximum filled segment number
      INTEGER INSGMX
      PARAMETER (INSGMX=15)
*     Lastly filled segment number
      INTEGER IISLST
      PARAMETER (IISLST=16)
*     No. of expansion of the segment table
      INTEGER INEXSG
      PARAMETER (INEXSG=17)
*     Current secondary size of the segment table
      INTEGER INSGM1
      PARAMETER (INSGM1=18)
*     Following items for fixed type banks
*     No. of words for master segment
      INTEGER INMAST
      PARAMETER (INMAST=19)
*     No. of data segments
      INTEGER INSEGF
      PARAMETER (INSEGF=20)
*     Length of data segments
      INTEGER ILSEGF
      PARAMETER (ILSEGF=21)
*  << Link pointer positions in the Bank Index Record >>
*     Pointer to the data record
      INTEGER ILDATR
      PARAMETER (ILDATR=1)
*   
*     
*     Length of the segments in the Bank Index Records
      INTEGER LINDSG
      PARAMETER (LINDSG=2)
C------------------------ END OF MACRO ZBSSYS ----------------------
