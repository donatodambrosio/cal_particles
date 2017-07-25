#include <sim_stop.h>
#include <ep_utils.h>

CALbyte caminalu(struct CALModel3D* modello)
{
  if (step < STEPS)
      return CAL_FALSE;

  return CAL_TRUE;
}
