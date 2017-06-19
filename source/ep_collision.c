#include <ep_collision.h>
#include <ep_utils.h>
#include <math.h>

#include <model.h>

void pointPlaneCollisionDetect (CALreal* rj, CALint* possible_collisions_along_directions)
{
  possible_collisions_along_directions[0] =
      possible_collisions_along_directions[1] =
      possible_collisions_along_directions[2] = 0;

  // x = CELL_SIDE plane OR x = (X_CELLS-1)*CELL_SIDE plane
  if (rj[0] == CELL_SIDE || rj[0] == (X_CELLS-1)*CELL_SIDE)
    possible_collisions_along_directions[0] = 1;

  // y = CELL_SIDE plane OR y = (Y_CELLS-1)*CELL_SIDE plane
  if (rj[1] == CELL_SIDE || rj[1] == (Y_CELLS-1)*CELL_SIDE)
    possible_collisions_along_directions[1] = 1;

  // z = CELL_SIDE plane OR z = (Z_CELLS-1)*CELL_SIDE plane
  if (rj[2] == CELL_SIDE || rj[2] == (Z_CELLS-1)*CELL_SIDE)
    possible_collisions_along_directions[2] = 1;
}

void collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal kn = KN;
  CALreal etha = ETHA;
  CALreal ri[3], vi[3], delta_Fi[3];
  CALreal rj[3], vj[3], delta_Fj[3];
  CALreal rij[3], dij, enij[3], vnij[3];
  CALreal delta_n;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
      {
        ri[0] = calGet3Dr(ca,Q.rx[slot],cell_x,cell_y,cell_z);
        ri[1] = calGet3Dr(ca,Q.ry[slot],cell_x,cell_y,cell_z);
        ri[2] = calGet3Dr(ca,Q.rz[slot],cell_x,cell_y,cell_z);
#ifdef VISCOUS
        vi[0] = calGet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z);
        vi[1] = calGet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z);
        vi[2] = calGet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z);
#endif
        delta_Fi[0] = 0.0;
        delta_Fi[1] = 0.0;
        delta_Fi[2] = 0.0;

        // inner particle-particle collision
        for (int inner_slot=slot+1; inner_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; inner_slot++)
          if (calGet3Di(ca, Q.imove[inner_slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
            {
              rj[0] = calGet3Dr(ca,Q.rx[inner_slot],cell_x,cell_y,cell_z);
              rj[1] = calGet3Dr(ca,Q.ry[inner_slot],cell_x,cell_y,cell_z);
              rj[2] = calGet3Dr(ca,Q.rz[inner_slot],cell_x,cell_y,cell_z);
#ifdef VISCOUS
              vj[0] = calGet3Dr(ca,Q.vx[inner_slot],cell_x,cell_y,cell_z);
              vj[1] = calGet3Dr(ca,Q.vy[inner_slot],cell_x,cell_y,cell_z);
              vj[2] = calGet3Dr(ca,Q.vz[inner_slot],cell_x,cell_y,cell_z);
#endif
              delta_Fj[0] = 0.0;
              delta_Fj[1] = 0.0;
              delta_Fj[2] = 0.0;

              dij = distance(ri, rj);
              if (dij < 2*PARTICLE_RADIUS)
                {
                  for (int k=0; k<3; k++)
                    {
                      rij[k] = rj[k] - ri[k];
                      enij[k] = rij[k] / dij;
                      delta_n = 2*PARTICLE_RADIUS - dij;
                      delta_Fi[k] += -kn * delta_n * enij[k];
                      delta_Fj[k] +=  kn * delta_n * enij[k];
#ifdef VISCOUS
                      vnij[k] = vj[k] - vi[k];
                      delta_Fi[k] += -etha * vnij[k];
                      delta_Fj[k] +=  etha * vnij[k];
#endif
                    }

                  // update phase
                  calSet3Dr(ca,Q.Fx[inner_slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca,Q.Fx[inner_slot],cell_x,cell_y,cell_z) + delta_Fj[0]);
                  calSet3Dr(ca,Q.Fy[inner_slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca,Q.Fy[inner_slot],cell_x,cell_y,cell_z) + delta_Fj[1]);
                  calSet3Dr(ca,Q.Fz[inner_slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca,Q.Fz[inner_slot],cell_x,cell_y,cell_z) + delta_Fj[2]);
                }
            }

        // outer particle-particle collision
        for (int n = 1; n<ca->sizeof_X; n++)
          for (int outer_slot=0; outer_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; outer_slot++)
            if (calGetX3Di(ca, Q.imove[outer_slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
              {
                rj[0] = calGetX3Dr(ca,Q.rx[outer_slot],cell_x,cell_y,cell_z,n);
                rj[1] = calGetX3Dr(ca,Q.ry[outer_slot],cell_x,cell_y,cell_z,n);
                rj[2] = calGetX3Dr(ca,Q.rz[outer_slot],cell_x,cell_y,cell_z,n);
#ifdef VISCOUS
                vj[0] = calGetX3Dr(ca,Q.vx[outer_slot],cell_x,cell_y,cell_z,n);
                vj[1] = calGetX3Dr(ca,Q.vy[outer_slot],cell_x,cell_y,cell_z,n);
                vj[2] = calGetX3Dr(ca,Q.vz[outer_slot],cell_x,cell_y,cell_z,n);
#endif
                dij = distance(ri, rj);
                if (dij < 2*PARTICLE_RADIUS)
                  {
                    for (int k=0; k<3; k++)
                      {
                        rij[k] = rj[k] - ri[k];
                        enij[k] = rij[k] / dij;
                        delta_n = 2*PARTICLE_RADIUS - dij;
                        delta_Fi[k] += -kn * delta_n * enij[k];
#ifdef VISCOUS
                        vnij[k] = vj[k] - vi[k];
                        delta_Fi[k] += -etha * vnij[k];
#endif
                      }
                  }
              }

        // particle-plane collision
        for (int n=1; n<ca->sizeof_X; n++)
          if (calGetX3Di(ca, Q.imove[0], cell_x,cell_y, cell_z, n) == PARTICLE_BORDER)
            {
              rj[0] = calGetX3Dr(ca,Q.rx[0],cell_x,cell_y,cell_z,n);
              rj[1] = calGetX3Dr(ca,Q.ry[0],cell_x,cell_y,cell_z,n);
              rj[2] = calGetX3Dr(ca,Q.rz[0],cell_x,cell_y,cell_z,n);
#ifdef VISCOUS
              vj[0] = calGetX3Dr(ca,Q.vx[0],cell_x,cell_y,cell_z,n);
              vj[1] = calGetX3Dr(ca,Q.vy[0],cell_x,cell_y,cell_z,n);
              vj[2] = calGetX3Dr(ca,Q.vz[0],cell_x,cell_y,cell_z,n);
#endif
              CALint e[3];
              pointPlaneCollisionDetect(rj, e);
              for (int k=0; k<3; k++)
                {
                  dij = sqrt((rj[k]-ri[k])*(rj[k]-ri[k]));
                  if (dij < PARTICLE_RADIUS)
                    {
                      rij[k] = rj[k] - ri[k];
                      enij[k] = rij[k] / dij;
                      delta_n = PARTICLE_RADIUS - dij;
                      delta_Fi[k] += -kn * delta_n * enij[k] * e[k];
#ifdef VISCOUS
                      vnij[k] = vj[k] - vi[k];
                      delta_Fi[k] += -etha * vnij[k] * e[k];
#endif
                    }
                }
            }// end of particle-plane collision


        // update phase
        calSet3Dr(ca, Q.Fx[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z) + delta_Fi[0]);
        calSet3Dr(ca, Q.Fy[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z) + delta_Fi[1]);
        calSet3Dr(ca, Q.Fz[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z) + delta_Fi[2]);

      } // end of for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
}
