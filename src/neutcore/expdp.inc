*
* $Id: expdp.inc,v 1.1 2007-01-20 07:10:40 skrep Exp $
*
* $Log: not supported by cvs2svn $
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
