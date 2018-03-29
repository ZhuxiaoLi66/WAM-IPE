// Need to define _GNU_SOURCE sometimes to get value of CPU_SETSIZE
#define _GNU_SOURCE
#include <stdio.h>
#include <sched.h>

int set_affinity_ (int *core)
{
  int retval = 0;
  cpu_set_t setmask;

  {
    // Set affinity so MPI task rank runs on a specific available core
    // Use the following to pin to a specific core
    CPU_ZERO (&setmask);
    CPU_SET (*core, &setmask);
    if (sched_setaffinity (0, sizeof (setmask), &setmask) < 0) {
      printf ("setaff: bad return from sched_setaffinity %d \n",*core);
      retval = -1;
    }
  }
  return retval;
}
