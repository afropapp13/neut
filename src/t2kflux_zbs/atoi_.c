#include <stdlib.h>

int
atoi_(char *arg, int len)
{
   int ret;

   ret = (int) strtol(arg, (char **)NULL, 10);
   return ret;

}
