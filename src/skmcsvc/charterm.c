/* 
   $Id: charterm.c,v 1.2 1998/04/04 20:57:00 bviren Exp $

   Terminate a FORTRAN character string at the first space or at the
   1024th character which ever comes first

   $Log: charterm.c,v $
   Revision 1.2  1998/04/04 20:57:00  bviren
   Iolib overhaul:
   Some Linux incompatibilities and generic bugs fixed (for the 2nd or
   3rd time). Many missing header files added (for 2nd time). Buffer
   overrun bugs in charterm.c and gethostname.c fixed. Some unused
   variables removed.

 */
#include <stdio.h>

void charterm(buffer)
      char buffer[1024];
{
      int i;

      for (i=0;buffer[i]!=' ' && buffer[i]!='\0' && i<1023;i++)
	  ;

      buffer[i]='\0';
      return;
}

