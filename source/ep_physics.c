#include <ep_physics.h>
#include <math.h>
#include <stdlib.h>

void applyForce(CALreal* F, CALreal* p0, CALreal* v0, CALreal m, CALreal t, CALreal* pf, CALreal* vf)
{
  CALreal a[3];
//  CALreal F[3];

  // F[0] =  0;
  // F[1] =  0;
  // F[2] = -m*G;

//  F[0] =  m*0 - 6*M_PI*AIR_VISCOSITY*PARTICLE_RADIUS*v0[0];
//  F[1] =  m*0 - 6*M_PI*AIR_VISCOSITY*PARTICLE_RADIUS*v0[1];
//  F[2] = -m*G - 6*M_PI*AIR_VISCOSITY*PARTICLE_RADIUS*v0[2];

  a[0] = F[0]/m;
  a[1] = F[1]/m;
  a[2] = F[2]/m;

  for (int i=0; i<3; i++)
    {
      vf[i] = v0[i]+a[i]*t;
      pf[i] = p0[i] + v0[i]*t + 0.5*a[i]*t*t;

      if (fabs(pf[i]-p0[i]) >= CELL_SIDE)
        {
          pf[i] = p0[i];
          vf[i] = v0[i];

#ifdef VERBOSE
          printf("ERROR: a particle displacemnt is greater than CELL_SIDE.\n");
          exit(EXIT_FAILURE);
#endif
        }
    }
}
