#include <ep_collision.h>
#include <ep_utils.h>
#include <math.h>

#include <model.h>

CALbyte pointPlaneCollision (CALreal* p, CALreal* v, CALreal* pB, CALreal* nB, CALreal particle_radius)
{
  CALint i = -1;
  CALreal d;

  // x = CELL_SIDE plane OR x = (X_CELLS-1)*CELL_SIDE plane
  if (pB[0] == CELL_SIDE || pB[0] == (X_CELLS-1)*CELL_SIDE)
    i = 0;

  // y = CELL_SIDE plane OR y = (Y_CELLS-1)*CELL_SIDE plane
  if (pB[1] == CELL_SIDE || pB[1] == (Y_CELLS-1)*CELL_SIDE)
    i = 1;

  // z = CELL_SIDE plane OR z = (Z_CELLS-1)*CELL_SIDE plane
  if (pB[2] == CELL_SIDE || pB[2] == (Z_CELLS-1)*CELL_SIDE)
    i = 2;

  if (i != -1 /*&& v[i] != 0*/)
    {
      d = sqrt((p[i]-pB[i])*(p[i]-pB[i]));
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
  CALreal p[3], v[3], B[3], N[3];

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
    {
      p[0] = calGet3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z);
      p[1] = calGet3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z);
      p[2] = calGet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z);

      v[0] = calGet3Dr(ca, Q.vx[slot],cell_x,cell_y,cell_z);
      v[1] = calGet3Dr(ca, Q.vy[slot],cell_x,cell_y,cell_z);
      v[2] = calGet3Dr(ca, Q.vz[slot],cell_x,cell_y,cell_z);

      for (int n=1; n<ca->sizeof_X; n++)
        {
          // particle-plane collision
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

        }
    }

}

