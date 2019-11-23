/* 
   $Id: ccopen.c,v 1.4 1998/04/04 20:56:53 bviren Exp $

   Fortran interface to open

   $Log: ccopen.c,v $
   Revision 1.4  1998/04/04 20:56:53  bviren
   Iolib overhaul:
   Some Linux incompatibilities and generic bugs fixed (for the 2nd or
   3rd time). Many missing header files added (for 2nd time). Buffer
   overrun bugs in charterm.c and gethostname.c fixed. Some unused
   variables removed.

 */
#include <fcntl.h>
#include <string.h>
#include <stdio.h>

void charterm(char *);

void ccopen_(fname, ironly, iwonly, irdwr, icreat, iappnd, iexcl, fi, lcom)
      char *fname;
      int  *ironly,*iwonly,*irdwr,*icreat,*iappnd,*iexcl;
      int  lcom;
      int  *fi;
{
      int  oflag;
      mode_t MODE;
      char buffer[1024];
      strncpy(buffer, fname, lcom);
      buffer[lcom] = '\0';
      charterm(buffer);

      oflag = 0;

      if(*ironly == 1){
        oflag = (oflag | O_RDONLY);
      }

      if(*iwonly == 1){
        oflag = (oflag | O_WRONLY);
      }

      if(*irdwr == 1){
        oflag = (oflag | O_RDWR  );
      }

      if(*icreat == 1){
        oflag = (oflag | O_CREAT );
      }

      if(*iappnd == 1){
        oflag = (oflag | O_APPEND);
      }

      if(*iexcl == 1){
        oflag = (oflag | O_EXCL  );
      }

      MODE=00777;
      if ((*fi = open(buffer,oflag,MODE)) == -1) {
        if(*ironly == 1 || *iwonly == 1){
	    fprintf(stderr," open err (ccopen) %s \n",buffer);
        }
      }

      return;
}
