#include <stdlib.h>

int
atoi_(char *str, int len)
{
  int retval;

  retval = strtol(str, (char **)(NULL), 0);

  return retval;

}



