#include <signal.h>

int
reset_signal_()
{

  signal(SIGFPE,SIG_DFL);
  signal(SIGILL,SIG_DFL);
  signal(SIGSEGV,SIG_DFL);

}
