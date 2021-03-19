/* $Id: iolib_get_tapemode.c,v 1.1 2000/10/13 14:23:33 sharkey Exp $ */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cfortran/cfortran.h>

int iolib_get_tapemode(char *filename) {

   struct stat buf;
   int rv;

   rv=stat(filename,&buf);

   if ((rv==0)&&(S_ISREG(buf.st_mode)))
     return(' ');
   else
     return('T');
}

FCALLSCFUN1(INT,iolib_get_tapemode,IOLIB_GET_TAPEMODE,iolib_get_tapemode,STRING)
