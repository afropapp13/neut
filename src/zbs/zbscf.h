/*------------------------------------------------------------------
fortran filename   : kzaddn.f
------------------------------------------------------------------*/
/*
#define kzaddn_ELEMS_2          ZTRINGV_NUM(1)
#define kzaddn_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZADDN(A1,A2,A3)  CCALLSFSUB3(KZADDN,kzaddn,INT,STRING,PINT,A1,A2,A3)

/*------------------------------------------------------------------
fortran filename   : kzbcr0.f
------------------------------------------------------------------*/
/*
#define kzbcr0_ELEMS_1          ZTRINGV_NUM(1)
#define kzbcr0_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZBCR0(A1,A2)  CCALLSFSUB2(KZBCR0,kzbcr0,STRING,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzbcr1.f
------------------------------------------------------------------*/
/*
#define kzbcr1_ELEMS_1          ZTRINGV_NUM(1)
#define kzbcr1_ELEMLEN_1        ZTRINGV_NUM(255)
#define kzbcr1_ELEMS_3          ZTRINGV_NUM(1)
#define kzbcr1_ELEMLEN_3        ZTRINGV_NUM(255)
#define kzbcr1_ELEMS_6          ZTRINGV_NUM(1)
#define kzbcr1_ELEMLEN_6        ZTRINGV_NUM(255)
*/

#define KZBCR1(A1,A2,A3,A4,A5,A6,A7)  CCALLSFSUB7(KZBCR1,kzbcr1,STRING,INT,STRING,INT,INT,STRING,PINT,A1,A2,A3,A4,A5,A6,A7)

/*------------------------------------------------------------------
fortran filename   : kzbcr2.f
------------------------------------------------------------------*/
/*
#define kzbcr2_ELEMS_1          ZTRINGV_NUM(1)
#define kzbcr2_ELEMLEN_1        ZTRINGV_NUM(255)
#define kzbcr2_ELEMS_3          ZTRINGV_NUM(1)
#define kzbcr2_ELEMLEN_3        ZTRINGV_NUM(255)
#define kzbcr2_ELEMS_6          ZTRINGV_NUM(1)
#define kzbcr2_ELEMLEN_6        ZTRINGV_NUM(255)
*/

#define KZBCR2(A1,A2,A3,A4,A5,A6,A7)  CCALLSFSUB7(KZBCR2,kzbcr2,STRING,INT,STRING,INT,INT,STRING,PINT,A1,A2,A3,A4,A5,A6,A7)

/*------------------------------------------------------------------
fortran filename   : kzbcr3.f
------------------------------------------------------------------*/
/*
#define kzbcr3_ELEMS_1          ZTRINGV_NUM(1)
#define kzbcr3_ELEMLEN_1        ZTRINGV_NUM(255)
#define kzbcr3_ELEMS_4          ZTRINGV_NUM(1)
#define kzbcr3_ELEMLEN_4        ZTRINGV_NUM(255)
#define kzbcr3_ELEMS_7          ZTRINGV_NUM(1)
#define kzbcr3_ELEMLEN_7        ZTRINGV_NUM(255)
*/

#define KZBCR3(A1,A2,A3,A4,A5,A6,A7,A8)  CCALLSFSUB8(KZBCR3,kzbcr3,STRING,INT,INT,STRING,INT,INT,STRING,PINT,A1,A2,A3,A4,A5,A6,A7,A8)

/*------------------------------------------------------------------
fortran filename   : kzbcre.f
------------------------------------------------------------------*/
/*
#define kzbcre_ELEMS_1          ZTRINGV_NUM(1)
#define kzbcre_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZBCRE(A1,A2)  CCALLSFSUB2(KZBCRE,kzbcre,STRING,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzbdel.f
------------------------------------------------------------------*/
/*
#define kzbdel_ELEMS_1          ZTRINGV_NUM(1)
#define kzbdel_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZBDEL(A1)  CCALLSFSUB1(KZBDEL,kzbdel,STRING,A1)

/*------------------------------------------------------------------
fortran filename   : kzbloc.f
------------------------------------------------------------------*/
/*
#define kzbloc_ELEMS_1          ZTRINGV_NUM(1)
#define kzbloc_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZBLOC(A1,A2)  CCALLSFSUB2(KZBLOC,kzbloc,STRING,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzbtyp.f
------------------------------------------------------------------*/
/*
#define kzbtyp_ELEMS_1          ZTRINGV_NUM(1)
#define kzbtyp_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZBTYP(A1,A2)  CCALLSFSUB2(KZBTYP,kzbtyp,STRING,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzchg0.f
------------------------------------------------------------------*/
/*
#define kzchg0_ELEMS_1          ZTRINGV_NUM(1)
#define kzchg0_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZCHG0(A1,A2,A3,A4)  CCALLSFSUB4(KZCHG0,kzchg0,STRING,INT,INT,INT,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : kzcldt.f
------------------------------------------------------------------*/

#define KZCLDT() CCALLSFSUB0(KZCLDT,kzcldt)

/*------------------------------------------------------------------
fortran filename   : kzcmpn.f
------------------------------------------------------------------*/

 PROTOCCALLSFFUN3(INT,KZCMPN,kzcmpn,INTV,INT,INTV)
#define KZCMPN(A2,A3,A4)  CCALLSFFUN3(KZCMPN,kzcmpn,INTV,INT,INTV,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : kzcrb0.f
------------------------------------------------------------------*/
/*
#define kzcrb0_ELEMS_2          ZTRINGV_NUM(1)
#define kzcrb0_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZCRB0(A1,A2,A3)  CCALLSFSUB3(KZCRB0,kzcrb0,INT,STRING,PINT,A1,A2,A3)

/*------------------------------------------------------------------
fortran filename   : kzcrb1.f
------------------------------------------------------------------*/
/*
#define kzcrb1_ELEMS_2          ZTRINGV_NUM(1)
#define kzcrb1_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZCRB1(A1,A2,A3,A4,A5)  CCALLSFSUB5(KZCRB1,kzcrb1,INT,STRING,INT,INT,INT,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : kzcrb2.f
------------------------------------------------------------------*/
/*
#define kzcrb2_ELEMS_2          ZTRINGV_NUM(1)
#define kzcrb2_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZCRB2(A1,A2,A3,A4,A5)  CCALLSFSUB5(KZCRB2,kzcrb2,INT,STRING,INT,INT,INT,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : kzcrb3.f
------------------------------------------------------------------*/
/*
#define kzcrb3_ELEMS_2          ZTRINGV_NUM(1)
#define kzcrb3_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZCRB3(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(KZCRB3,kzcrb3,INT,STRING,INT,INT,INT,INT,A1,A2,A3,A4,A5,A6)

/*------------------------------------------------------------------
fortran filename   : kzcrd0.f
------------------------------------------------------------------*/
/*
#define kzcrd0_ELEMS_2          ZTRINGV_NUM(1)
#define kzcrd0_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZCRD0(A1,A2,A3,A4)  CCALLSFSUB4(KZCRD0,kzcrd0,INT,STRING,INT,PINT,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : kzcrd1.f
------------------------------------------------------------------*/
/*
#define kzcrd1_ELEMS_2          ZTRINGV_NUM(1)
#define kzcrd1_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZCRD1(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(KZCRD1,kzcrd1,INT,STRING,INT,INT,INT,INT,A1,A2,A3,A4,A5,A6)

/*------------------------------------------------------------------
fortran filename   : kzcrd2.f
------------------------------------------------------------------*/
/*
#define kzcrd2_ELEMS_2          ZTRINGV_NUM(1)
#define kzcrd2_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZCRD2(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(KZCRD2,kzcrd2,INT,STRING,INT,INT,INT,INT,A1,A2,A3,A4,A5,A6)

/*------------------------------------------------------------------
fortran filename   : kzcrd3.f
------------------------------------------------------------------*/
/*
#define kzcrd3_ELEMS_2          ZTRINGV_NUM(1)
#define kzcrd3_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZCRD3(A1,A2,A3,A4,A5,A6,A7)  CCALLSFSUB7(KZCRD3,kzcrd3,INT,STRING,INT,INT,INT,INT,INT,A1,A2,A3,A4,A5,A6,A7)

/*------------------------------------------------------------------
fortran filename   : kzcreh.f
------------------------------------------------------------------*/

#define KZCREH(A1,A2)  CCALLSFSUB2(KZCREH,kzcreh,INT,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzcrem.f
------------------------------------------------------------------*/

#define KZCREM(A1,A2)  CCALLSFSUB2(KZCREM,kzcrem,INT,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzcrmh.f
------------------------------------------------------------------*/

#define KZCRMH() CCALLSFSUB0(KZCRMH,kzcrmh)

/*------------------------------------------------------------------
fortran filename   : kzctoi.f
------------------------------------------------------------------*/
/*
#define kzctoi_ELEMS_1          ZTRINGV_NUM(1)
#define kzctoi_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZCTOI(A1,A2,A3)  CCALLSFSUB3(KZCTOI,kzctoi,STRING,PINT,INT,A1,A2,A3)

/*------------------------------------------------------------------
fortran filename   : kzddiv.f
------------------------------------------------------------------*/

#define KZDDIV(A1)  CCALLSFSUB1(KZDDIV,kzddiv,INT,A1)

/*------------------------------------------------------------------
fortran filename   : kzdevt.f
------------------------------------------------------------------*/

#define KZDEVT(A1)  CCALLSFSUB1(KZDEVT,kzdevt,INT,A1)

/*------------------------------------------------------------------
fortran filename   : kzdsev.f
------------------------------------------------------------------*/

#define KZDSEV(A1)  CCALLSFSUB1(KZDSEV,kzdsev,INT,A1)

/*------------------------------------------------------------------
fortran filename   : kzdsto.f
------------------------------------------------------------------*/

#define KZDSTO(A1)  CCALLSFSUB1(KZDSTO,kzdsto,INT,A1)

/*------------------------------------------------------------------
fortran filename   : kzdtab.f
------------------------------------------------------------------*/

#define KZDTAB(A1)  CCALLSFSUB1(KZDTAB,kzdtab,INT,A1)

/*------------------------------------------------------------------
fortran filename   : kzdump.f
------------------------------------------------------------------*/
/*
#define kzdump_ELEMS_1          ZTRINGV_NUM(1)
#define kzdump_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZDUMP(A1,A2)  CCALLSFSUB2(KZDUMP,kzdump,STRING,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzeclr.f
------------------------------------------------------------------*/

#define KZECLR() CCALLSFSUB0(KZECLR,kzeclr)

/*------------------------------------------------------------------
fortran filename   : kzecur.f
------------------------------------------------------------------*/

#define KZECUR(A1,A2)  CCALLSFSUB2(KZECUR,kzecur,INT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzefst.f
------------------------------------------------------------------*/

#define KZEFST() CCALLSFSUB0(KZEFST,kzefst)

/*------------------------------------------------------------------
fortran filename   : kzejmp.f
------------------------------------------------------------------*/

#define KZEJMP(A1,A2)  CCALLSFSUB2(KZEJMP,kzejmp,INT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzelst.f
------------------------------------------------------------------*/

#define KZELST() CCALLSFSUB0(KZELST,kzelst)

/*------------------------------------------------------------------
fortran filename   : kzemrg.f
------------------------------------------------------------------*/

#define KZEMRG() CCALLSFSUB0(KZEMRG,kzemrg)

/*------------------------------------------------------------------
fortran filename   : kzemrk.f
------------------------------------------------------------------*/

#define KZEMRK() CCALLSFSUB0(KZEMRK,kzemrk)

/*------------------------------------------------------------------
fortran filename   : kzend.f
------------------------------------------------------------------*/

#define KZEND(A1)  CCALLSFSUB1(KZEND,kzend,INT,A1)

/*------------------------------------------------------------------
fortran filename   : kzenew.f
------------------------------------------------------------------*/

#define KZENEW() CCALLSFSUB0(KZENEW,kzenew)

/*------------------------------------------------------------------
fortran filename   : kzenum.f
------------------------------------------------------------------*/

#define KZENUM(A1,A2)  CCALLSFSUB2(KZENUM,kzenum,PINT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzenxt.f
------------------------------------------------------------------*/

#define KZENXT(A1)  CCALLSFSUB1(KZENXT,kzenxt,PINT,A1)

/*------------------------------------------------------------------
fortran filename   : kzeprv.f
------------------------------------------------------------------*/

#define KZEPRV(A1)  CCALLSFSUB1(KZEPRV,kzeprv,PINT,A1)

/*------------------------------------------------------------------
fortran filename   : kzexb0.f
------------------------------------------------------------------*/

#define KZEXB0(A1,A2)  CCALLSFSUB2(KZEXB0,kzexb0,INT,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzexb1.f
------------------------------------------------------------------*/

#define KZEXB1(A1,A2)  CCALLSFSUB2(KZEXB1,kzexb1,INT,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzexb2.f
------------------------------------------------------------------*/

#define KZEXB2(A1,A2)  CCALLSFSUB2(KZEXB2,kzexb2,INT,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzexdi.f
------------------------------------------------------------------*/

#define KZEXDI(A1)  CCALLSFSUB1(KZEXDI,kzexdi,PINT,A1)

/*------------------------------------------------------------------
fortran filename   : kzexdr.f
------------------------------------------------------------------*/

#define KZEXDR(A1,A2)  CCALLSFSUB2(KZEXDR,kzexdr,INT,INT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzfile.f
------------------------------------------------------------------*/

#define KZFILE(A1)  CCALLSFSUB1(KZFILE,kzfile,INT,A1)

/*------------------------------------------------------------------
fortran filename   : kzform.f
------------------------------------------------------------------*/
/*
#define kzform_ELEMS_1          ZTRINGV_NUM(1)
#define kzform_ELEMLEN_1        ZTRINGV_NUM(255)
#define kzform_ELEMS_2          ZTRINGV_NUM(1)
#define kzform_ELEMLEN_2        ZTRINGV_NUM(255)
#define kzform_ELEMS_3          ZTRINGV_NUM(1)
#define kzform_ELEMLEN_3        ZTRINGV_NUM(255)
*/

#define KZFORM(A1,A2,A3,A4,A5)  CCALLSFSUB5(KZFORM,kzform,STRING,STRING,STRING,INT,PINT,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : kzget.f
------------------------------------------------------------------*/
/*
#define kzget_ELEMS_1          ZTRINGV_NUM(1)
#define kzget_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZGET(A1,A2,A3,A4)  CCALLSFSUB4(KZGET,kzget,STRING,INT,PINT,INTV,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : kzget0.f
------------------------------------------------------------------*/
/*
#define kzget0_ELEMS_1          ZTRINGV_NUM(1)
#define kzget0_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZGET0(A1,A2,A3,A4)  CCALLSFSUB4(KZGET0,kzget0,STRING,INT,PINT,PINT,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : kzget1.f
------------------------------------------------------------------*/
/*
#define kzget1_ELEMS_1          ZTRINGV_NUM(1)
#define kzget1_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZGET1(A1,A2,A3,A4)  CCALLSFSUB4(KZGET1,kzget1,STRING,INT,PINT,PINT,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : kzget2.f
------------------------------------------------------------------*/
/*
#define kzget2_ELEMS_1          ZTRINGV_NUM(1)
#define kzget2_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZGET2(A1,A2,A3,A4,A5)  CCALLSFSUB5(KZGET2,kzget2,STRING,INT,INT,PINT,PINT,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : kzgeth.f
------------------------------------------------------------------*/

#define KZGETH(A1)  CCALLSFSUB1(KZGETH,kzgeth,PINT,A1)

/*------------------------------------------------------------------
fortran filename   : kzinit.f
------------------------------------------------------------------*/
/* [00756]SYNTAXIC-WARNING (kzinit.f) line 226
    ==>CHARACTER*4 CNNAM(NPNAME)/'*NA1','*NA2','*NA3','*NA4','*NA5','*NA6','*NA7','*NA8'/ */

#define KZINIT() CCALLSFSUB0(KZINIT,kzinit)

/*------------------------------------------------------------------
fortran filename   : kziofm.f
------------------------------------------------------------------*/
/*
#define kziofm_ELEMS_2          ZTRINGV_NUM(1)
#define kziofm_ELEMLEN_2        ZTRINGV_NUM(1)
*/

 PROTOCCALLSFFUN1(INT,KZIOFM,kziofm,STRING)
#define KZIOFM(A2)  CCALLSFFUN1(KZIOFM,kziofm,STRING,A2)

/*------------------------------------------------------------------
fortran filename   : kzitoc.f
------------------------------------------------------------------*/
/*
#define kzitoc_ELEMS_2          ZTRINGV_NUM(1)
#define kzitoc_ELEMLEN_2        ZTRINGV_NUM(255)
*/

#define KZITOC(A1,A2,A3)  CCALLSFSUB3(KZITOC,kzitoc,INTV,PSTRING,INT,A1,A2,A3)

/*------------------------------------------------------------------
fortran filename   : kzldat.f
------------------------------------------------------------------*/
/*
#define kzldat_ELEMS_1          ZTRINGV_NUM(1)
#define kzldat_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZLDAT(A1,A2)  CCALLSFSUB2(KZLDAT,kzldat,STRING,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzleng.f
------------------------------------------------------------------*/
/*
#define kzleng_ELEMS_2          ZTRINGV_NUM(1)
#define kzleng_ELEMLEN_2        ZTRINGV_NUM(255)
*/

 PROTOCCALLSFFUN1(INT,KZLENG,kzleng,STRING)
#define KZLENG(A2)  CCALLSFFUN1(KZLENG,kzleng,STRING,A2)

/*------------------------------------------------------------------
fortran filename   : kzlseg.f
------------------------------------------------------------------*/
/*
#define kzlseg_ELEMS_1          ZTRINGV_NUM(1)
#define kzlseg_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZLSEG(A1,A2,A3)  CCALLSFSUB3(KZLSEG,kzlseg,STRING,INT,PINT,A1,A2,A3)

/*------------------------------------------------------------------
fortran filename   : kzmcur.f
------------------------------------------------------------------*/

#define KZMCUR(A1,A2)  CCALLSFSUB2(KZMCUR,kzmcur,INT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzmfst.f
------------------------------------------------------------------*/

#define KZMFST() CCALLSFSUB0(KZMFST,kzmfst)

/*------------------------------------------------------------------
fortran filename   : kzmjmp.f
------------------------------------------------------------------*/

#define KZMJMP(A1,A2)  CCALLSFSUB2(KZMJMP,kzmjmp,INT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzmkdt.f
------------------------------------------------------------------*/

#define KZMKDT() CCALLSFSUB0(KZMKDT,kzmkdt)

/*------------------------------------------------------------------
fortran filename   : kzmlst.f
------------------------------------------------------------------*/

#define KZMLST() CCALLSFSUB0(KZMLST,kzmlst)

/*------------------------------------------------------------------
fortran filename   : kzmnum.f
------------------------------------------------------------------*/

#define KZMNUM(A1,A2)  CCALLSFSUB2(KZMNUM,kzmnum,PINT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzmnxt.f
------------------------------------------------------------------*/

#define KZMNXT(A1)  CCALLSFSUB1(KZMNXT,kzmnxt,PINT,A1)

/*------------------------------------------------------------------
fortran filename   : kzmprv.f
------------------------------------------------------------------*/

#define KZMPRV(A1)  CCALLSFSUB1(KZMPRV,kzmprv,PINT,A1)

/*------------------------------------------------------------------
fortran filename   : kzname.f
------------------------------------------------------------------*/
/*
#define kzname_ELEMS_2          ZTRINGV_NUM(1)
#define kzname_ELEMLEN_2        ZTRINGV_NUM(255)
*/

 PROTOCCALLSFFUN1(INT,KZNAME,kzname,STRING)
#define KZNAME(A2)  CCALLSFFUN1(KZNAME,kzname,STRING,A2)

/*------------------------------------------------------------------
fortran filename   : kznseg.f
------------------------------------------------------------------*/
/*
#define kznseg_ELEMS_1          ZTRINGV_NUM(1)
#define kznseg_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZNSEG(A1,A2)  CCALLSFSUB2(KZNSEG,kznseg,STRING,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kznsg0.f
------------------------------------------------------------------*/
/*
#define kznsg0_ELEMS_1          ZTRINGV_NUM(1)
#define kznsg0_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZNSG0(A1,A2)  CCALLSFSUB2(KZNSG0,kznsg0,STRING,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kznssg.f
------------------------------------------------------------------*/
/*
#define kznssg_ELEMS_1          ZTRINGV_NUM(1)
#define kznssg_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZNSSG(A1,A2,A3)  CCALLSFSUB3(KZNSSG,kznssg,STRING,INT,PINT,A1,A2,A3)

/*------------------------------------------------------------------
fortran filename   : kzoptn.f
------------------------------------------------------------------*/
/*
#define kzoptn_ELEMS_1          ZTRINGV_NUM(1)
#define kzoptn_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZOPTN(A1)  CCALLSFSUB1(KZOPTN,kzoptn,STRING,A1)

/*------------------------------------------------------------------
fortran filename   : kzpdir.f
------------------------------------------------------------------*/
/*
#define kzpdir_ELEMS_1          ZTRINGV_NUM(1)
#define kzpdir_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZPDIR(A1,A2)  CCALLSFSUB2(KZPDIR,kzpdir,STRING,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzput0.f
------------------------------------------------------------------*/
/*
#define kzput0_ELEMS_1          ZTRINGV_NUM(1)
#define kzput0_ELEMLEN_1        ZTRINGV_NUM(255)
#define kzput0_ELEMS_4          ZTRINGV_NUM(1)
#define kzput0_ELEMLEN_4        ZTRINGV_NUM(1)
*/

#define KZPUT0(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(KZPUT0,kzput0,STRING,INT,INT,STRING,INT,INTV,A1,A2,A3,A4,A5,A6)

/*------------------------------------------------------------------
fortran filename   : kzput1.f
------------------------------------------------------------------*/
/*
#define kzput1_ELEMS_1          ZTRINGV_NUM(1)
#define kzput1_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZPUT1(A1,A2,A3,A4,A5)  CCALLSFSUB5(KZPUT1,kzput1,STRING,INT,INT,INT,INTV,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : kzput2.f
------------------------------------------------------------------*/
/*
#define kzput2_ELEMS_1          ZTRINGV_NUM(1)
#define kzput2_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZPUT2(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(KZPUT2,kzput2,STRING,INT,INT,INT,INT,INTV,A1,A2,A3,A4,A5,A6)

/*------------------------------------------------------------------
fortran filename   : kzputh.f
------------------------------------------------------------------*/

#define KZPUTH(A1)  CCALLSFSUB1(KZPUTH,kzputh,INT,A1)

/*------------------------------------------------------------------
fortran filename   : kzread.f
------------------------------------------------------------------*/

#define KZREAD(A1,A2)  CCALLSFSUB2(KZREAD,kzread,INT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzreda.f
------------------------------------------------------------------*/

#define KZREDA(A1,A2)  CCALLSFSUB2(KZREDA,kzreda,INT,PINT,A1,A2)

/*------------------------------------------------------------------
fortran filename   : kzrep0.f
------------------------------------------------------------------*/
/*
#define kzrep0_ELEMS_1          ZTRINGV_NUM(1)
#define kzrep0_ELEMLEN_1        ZTRINGV_NUM(255)
#define kzrep0_ELEMS_3          ZTRINGV_NUM(1)
#define kzrep0_ELEMLEN_3        ZTRINGV_NUM(1)
*/

#define KZREP0(A1,A2,A3,A4,A5)  CCALLSFSUB5(KZREP0,kzrep0,STRING,INT,STRING,INT,INTV,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : kzrep1.f
------------------------------------------------------------------*/
/*
#define kzrep1_ELEMS_1          ZTRINGV_NUM(1)
#define kzrep1_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZREP1(A1,A2,A3,A4)  CCALLSFSUB4(KZREP1,kzrep1,STRING,INT,INT,INTV,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : kzrep2.f
------------------------------------------------------------------*/
/*
#define kzrep2_ELEMS_1          ZTRINGV_NUM(1)
#define kzrep2_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZREP2(A1,A2,A3,A4,A5)  CCALLSFSUB5(KZREP2,kzrep2,STRING,INT,INT,INT,INTV,A1,A2,A3,A4,A5)

/*------------------------------------------------------------------
fortran filename   : kzwpos.f
------------------------------------------------------------------*/

#define KZWPOS(A1,A2,A3,A4)  CCALLSFSUB4(KZWPOS,kzwpos,INT,INT,INT,PINT,A1,A2,A3,A4)

/*------------------------------------------------------------------
fortran filename   : kzwps2.f
------------------------------------------------------------------*/
/*
#define kzwps2_ELEMS_1          ZTRINGV_NUM(1)
#define kzwps2_ELEMLEN_1        ZTRINGV_NUM(255)
*/

#define KZWPS2(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(KZWPS2,kzwps2,STRING,INT,INT,INT,INT,PINT,A1,A2,A3,A4,A5,A6)

/*------------------------------------------------------------------
fortran filename   : kzwrit.f
------------------------------------------------------------------*/

#define KZWRIT(A1)  CCALLSFSUB1(KZWRIT,kzwrit,INT,A1)

