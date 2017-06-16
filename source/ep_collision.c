#include <ep_collision.h>
#include <ep_utils.h>
#include <math.h>

#include <model.h>

void pointPlaneCollisionDetect (CALreal* rj, CALint* collisions_along_directions)
{
  collisions_along_directions[0] = collisions_along_directions[1] = collisions_along_directions[2] = 0;

  // x = CELL_SIDE plane OR x = (X_CELLS-1)*CELL_SIDE plane
  if (rj[0] == CELL_SIDE || rj[0] == (X_CELLS-1)*CELL_SIDE)
    collisions_along_directions[0] = 1;

  // y = CELL_SIDE plane OR y = (Y_CELLS-1)*CELL_SIDE plane
  if (rj[1] == CELL_SIDE || rj[1] == (Y_CELLS-1)*CELL_SIDE)
    collisions_along_directions[1] = 1;

  // z = CELL_SIDE plane OR z = (Z_CELLS-1)*CELL_SIDE plane
  if (rj[2] == CELL_SIDE || rj[2] == (Z_CELLS-1)*CELL_SIDE)
    collisions_along_directions[2] = 1;
}

void collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal kn = KN;
  CALreal ri[3];
  CALreal rj[3];
  CALreal rij[3], dij, enij[3];
  CALreal delta_n;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
      {
        ri[0] = calGet3Dr(ca, Q.rx[slot],cell_x,cell_y,cell_z);
        ri[1] = calGet3Dr(ca, Q.ry[slot],cell_x,cell_y,cell_z);
        ri[2] = calGet3Dr(ca, Q.rz[slot],cell_x,cell_y,cell_z);

        // particle-plane collision
        for (int n=1; n<ca->sizeof_X; n++)
          if (calGetX3Di(ca, Q.imove[0], cell_x,cell_y, cell_z, n) == PARTICLE_BORDER)
            {
              rj[0] = calGetX3Dr(ca, Q.rx[0],cell_x,cell_y,cell_z,n);
              rj[1] = calGetX3Dr(ca, Q.ry[0],cell_x,cell_y,cell_z,n);
              rj[2] = calGetX3Dr(ca, Q.rz[0],cell_x,cell_y,cell_z,n);

              CALint e[3];
              CALreal delta_n[3];
              pointPlaneCollisionDetect(rj, e);
              for (int k=0; k<3; k++)
                {
                  dij = sqrt((rj[k]-ri[k])*(rj[k]-ri[k]));
                  if (dij < PARTICLE_RADIUS)
                    {
                      rij[k] = rj[k] - ri[k];
                      enij[k] = rij[k] / dij;
                      delta_n[k] = PARTICLE_RADIUS - dij;
                    }
                  else
                    {
                      enij[k] = 0.0;
                      delta_n[k] = 0.0;
                    }
                }

              calSet3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z)-kn*delta_n[0]*enij[0] * e[0]);
              calSet3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z)-kn*delta_n[1]*enij[1] * e[1]);
              calSet3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z)-kn*delta_n[2]*enij[2] * e[2]);
            }// end of particle-plane collision

        //continue;

        // inner particle-particle collision
        for (int inner_slot=slot+1; inner_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; inner_slot++)
          if (calGet3Di(ca, Q.imove[inner_slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
            {
              rj[0] = calGet3Dr(ca, Q.rx[inner_slot],cell_x,cell_y,cell_z);
              rj[1] = calGet3Dr(ca, Q.ry[inner_slot],cell_x,cell_y,cell_z);
              rj[2] = calGet3Dr(ca, Q.rz[inner_slot],cell_x,cell_y,cell_z);

              dij = distance(ri, rj);
              if (dij < 2*PARTICLE_RADIUS)
                {
                  for (int k=0; k<3; k++)
                    {
                      rij[k] = rj[k] - ri[k];
                      enij[k] = rij[k] / dij;
                    }
                  delta_n = 2*PARTICLE_RADIUS - dij;

                  calSet3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z)-kn*delta_n*enij[0]);
                  calSet3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z)-kn*delta_n*enij[1]);
                  calSet3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z)-kn*delta_n*enij[2]);

                  calSet3Dr(ca, Q.Fx[inner_slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[inner_slot], cell_x,cell_y,cell_z)+kn*delta_n*enij[0]);
                  calSet3Dr(ca, Q.Fy[inner_slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[inner_slot], cell_x,cell_y,cell_z)+kn*delta_n*enij[1]);
                  calSet3Dr(ca, Q.Fz[inner_slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[inner_slot], cell_x,cell_y,cell_z)+kn*delta_n*enij[2]);
                }
            }

        // outer particle-particle collision
        for (int n = 1; n<ca->sizeof_X; n++)
          for (int outer_slot=0; outer_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; outer_slot++)
            if (calGetX3Di(ca, Q.imove[outer_slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
              {
                rj[0] = calGetX3Dr(ca, Q.rx[outer_slot],cell_x,cell_y,cell_z,n);
                rj[1] = calGetX3Dr(ca, Q.ry[outer_slot],cell_x,cell_y,cell_z,n);
                rj[2] = calGetX3Dr(ca, Q.rz[outer_slot],cell_x,cell_y,cell_z,n);

                dij = distance(ri, rj);
                if (dij < 2*PARTICLE_RADIUS)
                  {
                    for (int k=0; k<3; k++)
                      {
                        rij[k] = rj[k] - ri[k];
                        enij[k] = rij[k] / dij;
                      }
                    delta_n = 2*PARTICLE_RADIUS - dij;

                    calSet3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z)-kn*delta_n*enij[0]);
                    calSet3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z)-kn*delta_n*enij[1]);
                    calSet3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z)-kn*delta_n*enij[2]);
                  }
              }

      } // end of for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
}
