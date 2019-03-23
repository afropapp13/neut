#ifndef __T2KNEUT_INCLUDED__

#include <TRandom3.h>

#ifdef MAIN
TRandom3 *global_tr3;
#else
extern TRandom3 *global_tr3;
#endif


#define __T2KNEUT_INCLUDED__
#endif
