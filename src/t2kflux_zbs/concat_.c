#include <strings.h>
#include <stdio.h>
void
concat_(char *s1, char *s2, char *s3, int l1, int l2, int l3)
{
  char *s_end;

  memset(s3,0,l3);

  s_end = strstr(s1," ");
  if (s_end != NULL)  *s_end='\0';

  snprintf(s3,l3,"%s%s",s1,s2);

  s_end = strstr(s1," ");
  if (s_end != NULL)  *s_end='\0';
}
