#include <ep_physics.h>
#include <ep_utils.h>
#include <math.h>
#include <stdlib.h>

// #define AIR_VISCOSITY 1.81e-5

void applyForce(CALreal* F, CALreal* p0, CALreal* v0, CALreal m, CALreal t, CALreal* pf, CALreal* vf)
{
//  CALreal F[3];
//  F[0] =  m*0 - 6*M_PI*AIR_VISCOSITY*PARTICLE_RADIUS*v0[0];
//  F[1] =  m*0 - 6*M_PI*AIR_VISCOSITY*PARTICLE_RADIUS*v0[1];
//  F[2] = -m*G - 6*M_PI*AIR_VISCOSITY*PARTICLE_RADIUS*v0[2];

  CALreal a[3];

  a[0] = F[0]/m;
  a[1] = F[1]/m;
  a[2] = F[2]/m;

  for (int i=0; i<3; i++)
    {
      vf[i] = v0[i]+a[i]*t;
      //pf[i] = p0[i] + v0[i]*t + 0.5*a[i]*t*t;
      pf[i] = p0[i] + vf[i]*t;
    }

  CALreal displacement = distance(p0, pf);
  if (displacement >= PARTICLE_RADIUS)
    {
      printf("ERROR: a particle displacemnt is greater than CELL_SIDE.\n");
#ifdef VERBOSE
      printf("F = (%f, %f, %f)\n", F[0], F[1], F[2]);
      printf("a = (%f, %f, %f)\n", a[0], a[1], a[2]);
      printf("p0 = (%f, %f, %f)\n", p0[0], p0[1], p0[2]);
      printf("pf = (%f, %f, %f)\n", pf[0], pf[1], pf[2]);
      printf("particle displacemnt = %f\n", displacement);
      printf("CELL_SIDE = %f \n", CELL_SIDE);
#endif
      exit(EXIT_FAILURE);
    }
}

void resetF(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal F[3];
  CALreal gravity_force = -PARTICLE_MASS*G;

  F[0] = 0.0;
  F[1] = 0.0;
  F[2] = 0.0;
#ifdef GRAVITY
  F[2] = gravity_force;
#endif

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    //if (calGet3Di(ca, Q.ID[slot],cell_x,cell_y,cell_z) > NULL_ID)
    if (calGetBuffer3DElement(Q_ID_current[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z) > NULL_ID)
      {
        //calSet3Dr(ca, Q.Fx[slot],cell_x,cell_y,cell_z,F[0]);
        //calSet3Dr(ca, Q.Fy[slot],cell_x,cell_y,cell_z,F[1]);
        //calSet3Dr(ca, Q.Fz[slot],cell_x,cell_y,cell_z,F[2]);
        calSetBuffer3DElement(Q_Fx_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,F[0]);
        calSetBuffer3DElement(Q_Fy_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,F[1]);
        calSetBuffer3DElement(Q_Fz_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,F[2]);
      }
}
