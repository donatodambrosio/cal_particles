#include <sim_stop.h>
#include <ep_utils.h>

CALbyte caminalu(struct CALModel3D* modello)
{
  CALint step = a_simulazioni->step;

  if (step <= STEPS)
    {
      printSummary(modello);
      return CAL_FALSE;
    }

  return CAL_TRUE;
}
