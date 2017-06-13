#include <ep_collision.h>
#include <ep_utils.h>
#include <math.h>

#include <model.h>

CALbyte pointPlaneCollision (CALreal* p, CALreal* v, CALreal* B, CALreal* N, CALreal particle_radius)
{
  CALint i = -1;
  CALreal d;

  // x = CELL_SIDE plane OR x = (X_CELLS-1)*CELL_SIDE plane
  if (B[0] == CELL_SIDE || B[0] == (X_CELLS-1)*CELL_SIDE)
    i = 0;

  // y = CELL_SIDE plane OR y = (Y_CELLS-1)*CELL_SIDE plane
  if (B[1] == CELL_SIDE || B[1] == (Y_CELLS-1)*CELL_SIDE)
    i = 1;

  // z = CELL_SIDE plane OR z = (Z_CELLS-1)*CELL_SIDE plane
  if (B[2] == CELL_SIDE || B[2] == (Z_CELLS-1)*CELL_SIDE)
    i = 2;

  if (i != -1 && v[i]*N[i] < 0.0)
    {
      d = sqrt((p[i]-B[i])*(p[i]-B[i]));
      if (d < particle_radius)
        {
          CALreal p_vel = v[i];
          v[i] = -p_vel;
          return CAL_TRUE;
        }
    }

  return CAL_FALSE;
}

void collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal F[3], p[3], v[3];
  CALreal B[3], N[3];
  CALreal Fi[3], pi[3], vi[3];
  CALreal Fn[3], pn[3], vn[3];

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
      {
        F[0] = calGet3Dr(ca, Q.Fx[slot],cell_x,cell_y,cell_z);
        F[1] = calGet3Dr(ca, Q.Fy[slot],cell_x,cell_y,cell_z);
        F[2] = calGet3Dr(ca, Q.Fz[slot],cell_x,cell_y,cell_z);

        p[0] = calGet3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z);
        p[1] = calGet3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z);
        p[2] = calGet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z);

        v[0] = calGet3Dr(ca, Q.vx[slot],cell_x,cell_y,cell_z);
        v[1] = calGet3Dr(ca, Q.vy[slot],cell_x,cell_y,cell_z);
        v[2] = calGet3Dr(ca, Q.vz[slot],cell_x,cell_y,cell_z);

        // particle-plane collision
        for (int n=1; n<ca->sizeof_X; n++)
          if (calGetX3Di(ca, Q.imove[0], cell_x,cell_y, cell_z, n) == PARTICLE_BORDER)
            {
              B[0] = calGetX3Dr(ca, Q.px[0],cell_x,cell_y,cell_z,n);
              B[1] = calGetX3Dr(ca, Q.py[0],cell_x,cell_y,cell_z,n);
              B[2] = calGetX3Dr(ca, Q.pz[0],cell_x,cell_y,cell_z,n);

              N[0] = calGetX3Dr(ca, Q.vx[0],cell_x,cell_y,cell_z,n);
              N[1] = calGetX3Dr(ca, Q.vy[0],cell_x,cell_y,cell_z,n);
              N[2] = calGetX3Dr(ca, Q.vz[0],cell_x,cell_y,cell_z,n);

              if(pointPlaneCollision(p, v, B, N, PARTICLE_RADIUS))
                {
                  calSet3Dr(ca, Q.vx[slot],cell_x,cell_y,cell_z,v[0]);
                  calSet3Dr(ca, Q.vy[slot],cell_x,cell_y,cell_z,v[1]);
                  calSet3Dr(ca, Q.vz[slot],cell_x,cell_y,cell_z,v[2]);
                }
            }

        // particle-particle collision
        CALreal k = 0.0001;
        for (int inner_slot=slot+1; inner_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; inner_slot++)
          if (calGet3Di(ca, Q.imove[inner_slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
            {
              pi[0] = calGet3Dr(ca, Q.px[inner_slot],cell_x,cell_y,cell_z);
              pi[1] = calGet3Dr(ca, Q.py[inner_slot],cell_x,cell_y,cell_z);
              pi[2] = calGet3Dr(ca, Q.pz[inner_slot],cell_x,cell_y,cell_z);

              if (distance(p, pi) < 2*PARTICLE_RADIUS)
                {
                  Fi[0] = calGet3Dr(ca, Q.Fx[inner_slot],cell_x,cell_y,cell_z);
                  Fi[1] = calGet3Dr(ca, Q.Fy[inner_slot],cell_x,cell_y,cell_z);
                  Fi[2] = calGet3Dr(ca, Q.Fz[inner_slot],cell_x,cell_y,cell_z);

                  vi[0] = calGet3Dr(ca, Q.vx[inner_slot],cell_x,cell_y,cell_z);
                  vi[1] = calGet3Dr(ca, Q.vy[inner_slot],cell_x,cell_y,cell_z);
                  vi[2] = calGet3Dr(ca, Q.vz[inner_slot],cell_x,cell_y,cell_z);

                  CALreal delta_n[3];
                  for (int i=0; i<3; i++)
                    delta_n[i] = (2*PARTICLE_RADIUS - (p[i] - pi[i])) / PARTICLE_RADIUS;

                  calSet3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z)+k*delta_n[0]);
                  calSet3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z)+k*delta_n[1]);
                  calSet3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z)+k*delta_n[2]);

                  calSet3Dr(ca, Q.Fx[inner_slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[inner_slot], cell_x,cell_y,cell_z)-k*delta_n[0]);
                  calSet3Dr(ca, Q.Fy[inner_slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[inner_slot], cell_x,cell_y,cell_z)-k*delta_n[1]);
                  calSet3Dr(ca, Q.Fz[inner_slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[inner_slot], cell_x,cell_y,cell_z)-k*delta_n[2]);
                }
            }
        for (int n = 1; n<ca->sizeof_X; n++)
          for (int other_slot=0; other_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; other_slot++)
            if (calGetX3Di(ca, Q.imove[other_slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
              {
                pn[0] = calGetX3Dr(ca, Q.px[other_slot],cell_x,cell_y,cell_z,n);
                pn[1] = calGetX3Dr(ca, Q.py[other_slot],cell_x,cell_y,cell_z,n);
                pn[2] = calGetX3Dr(ca, Q.pz[other_slot],cell_x,cell_y,cell_z,n);

                if (distance(p, pn) < 2*PARTICLE_RADIUS)
                  {
                    Fn[0] = calGet3Dr(ca, Q.Fx[other_slot],cell_x,cell_y,cell_z);
                    Fn[1] = calGet3Dr(ca, Q.Fy[other_slot],cell_x,cell_y,cell_z);
                    Fn[2] = calGet3Dr(ca, Q.Fz[other_slot],cell_x,cell_y,cell_z);

                    vi[0] = calGet3Dr(ca, Q.vx[other_slot],cell_x,cell_y,cell_z);
                    vi[1] = calGet3Dr(ca, Q.vy[other_slot],cell_x,cell_y,cell_z);
                    vi[2] = calGet3Dr(ca, Q.vz[other_slot],cell_x,cell_y,cell_z);

                    CALreal delta_n[3];
                    for (int i=0; i<3; i++)
                      delta_n[i] = (2*PARTICLE_RADIUS - (p[i] - pi[i])) / PARTICLE_RADIUS;

                    calSet3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z)+k*delta_n[0]);
                    calSet3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z)+k*delta_n[1]);
                    calSet3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z)+k*delta_n[2]);

                    calSet3Dr(ca, Q.Fx[other_slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[other_slot], cell_x,cell_y,cell_z)-k*delta_n[0]);
                    calSet3Dr(ca, Q.Fy[other_slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[other_slot], cell_x,cell_y,cell_z)-k*delta_n[1]);
                    calSet3Dr(ca, Q.Fz[other_slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[other_slot], cell_x,cell_y,cell_z)-k*delta_n[2]);
                  }
              }
      }

}

