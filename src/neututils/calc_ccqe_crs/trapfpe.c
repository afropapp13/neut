#define _GNU_SOURCE 1
#include <signal.h>
#include <fenv.h>
#include <fpu_control.h>
static void __attribute__ ((constructor))
trapfpe ()
{
  /* Enable some exceptions.  At startup all exceptions are masked.  */

#ifdef ARCH_AARM64
   fpu_control_t cw = _FPU_DEFAULT;
#else
   fpu_control_t cw = (_FPU_DEFAULT & ~_FPU_EXTENDED) | _FPU_DOUBLE;
#endif
   _FPU_SETCW(cw);

  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  signal(SIGFPE, SIG_DFL);
  signal(SIGBUS, SIG_DFL);
  signal(SIGSEGV, SIG_DFL);

}

trapfpe_()
{
  trapfpe();
}
