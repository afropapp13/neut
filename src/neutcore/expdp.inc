*
* $Id: expdp.inc,v 1.1.1.1 2016/11/04 20:31:28 mdunkman Exp $
*
* $Log: expdp.inc,v $
* Revision 1.1.1.1  2016/11/04 20:31:28  mdunkman
* Import of NEUT 5.3.3 v1r27 for use with T2KReWeight v1r27 (Winter OA 2017)
*
* Revision 1.1  2002/08/29 13:55:47  saji
* Forgot to add 'expdp.inc'
*
* Revision 1.1.1.2  1996/10/30 08:27:10  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:28:48  plothow
* Version 7.01
*
*
*
* expdp.inc
*
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION
#endif
#if defined(CERNLIB_SINGLE)
      REAL
#endif
