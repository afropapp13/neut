#include <NeutRootHandlers.h>

extern "C"
{
  int rootopen_(char * /* filename */, char * /* option */, int , int);
  int rootwrite_(char * /* filename */, int);
  int rootclose_(char * /* filename */, int);
};

int
rootopen_(char *filename, char *opt, int i, int j)
{

  NeutRootHandlers nrh;
  int ret;

  nrh.nulltermstr(filename,i);
  nrh.nulltermstr(opt,j);
  
  ret = nrh.open(filename, opt);
  return ret;
}

int
rootwrite_(char *filename, int i)
{
  NeutRootHandlers nrh;
  int ret;

  nrh.nulltermstr(filename,i);
  ret = nrh.write(filename);
  return ret;
}

int
rootclose_(char *filename, int i)
{
  NeutRootHandlers nrh;
  int ret;

  nrh.nulltermstr(filename,i);
  ret = nrh.close(filename);
  return ret;
}

