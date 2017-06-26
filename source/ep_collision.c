#include <ep_collision.h>
#include <ep_utils.h>
#include <math.h>

#include <model.h>

void inner_collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal kn = KN;
  CALreal etha = ETHA;
  CALreal ri[3], vi[3], delta_Fi[3];
  CALreal rj[3], vj[3], delta_Fj[3];
  CALreal rij[3], dij, enij[3], vij[3], vnij;
  CALreal delta_n;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.ID[slot],cell_x,cell_y,cell_z) > NULL_ID)
      {
        ri[0] = calGet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z);
        ri[1] = calGet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z);
        ri[2] = calGet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z);
#ifdef VISCOELASTIC
        vi[0] = calGet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z);
        vi[1] = calGet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z);
        vi[2] = calGet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z);
#endif
        delta_Fi[0] = 0.0;
        delta_Fi[1] = 0.0;
        delta_Fi[2] = 0.0;

        // inner particle-particle collision
        for (int inner_slot=slot+1; inner_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; inner_slot++)
          if (calGet3Di(ca, Q.ID[inner_slot],cell_x,cell_y,cell_z) > NULL_ID )
            {
              rj[0] = calGet3Dr(ca,Q.px[inner_slot],cell_x,cell_y,cell_z);
              rj[1] = calGet3Dr(ca,Q.py[inner_slot],cell_x,cell_y,cell_z);
              rj[2] = calGet3Dr(ca,Q.pz[inner_slot],cell_x,cell_y,cell_z);
#ifdef VISCOELASTIC
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
#ifdef VISCOELASTIC
                      vij[k] = vi[k] - vj[k];
#endif
                    }
#ifdef VISCOELASTIC
                  vnij = 0.0;
                  for (int k=0; k<3; k++)
                    vnij += vij[k] * enij[k];

                  for (int k=0; k<3; k++)
                    {
                      delta_Fi[k] += -etha * vnij * enij[k];
                      delta_Fj[k] +=  etha * vnij * enij[k];
                    }
#endif
                  // update phase
                  calSet3Dr(ca,Q.Fx[inner_slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca,Q.Fx[inner_slot],cell_x,cell_y,cell_z) + delta_Fj[0]);
                  calSet3Dr(ca,Q.Fy[inner_slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca,Q.Fy[inner_slot],cell_x,cell_y,cell_z) + delta_Fj[1]);
                  calSet3Dr(ca,Q.Fz[inner_slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca,Q.Fz[inner_slot],cell_x,cell_y,cell_z) + delta_Fj[2]);
                } //if dij < 2*PARTICLE_RADIUS
            } //for inner_slot AND if
        // update phase
        calSet3Dr(ca, Q.Fx[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z) + delta_Fi[0]);
        calSet3Dr(ca, Q.Fy[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z) + delta_Fi[1]);
        calSet3Dr(ca, Q.Fz[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z) + delta_Fi[2]);
      }
}

void outer_collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal kn = KN;
  CALreal etha = ETHA;
  CALreal pi[3], vi[3], delta_Fi[3];
  CALreal pj[3], vj[3];
  CALreal rij[3], dij, enij[3], vij[3], vnij;
  CALreal delta_n;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.ID[slot],cell_x,cell_y,cell_z) > NULL_ID)
      {
        pi[0] = calGet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z);
        pi[1] = calGet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z);
        pi[2] = calGet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z);
#ifdef VISCOELASTIC
        vi[0] = calGet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z);
        vi[1] = calGet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z);
        vi[2] = calGet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z);
#endif
        delta_Fi[0] = 0.0;
        delta_Fi[1] = 0.0;
        delta_Fi[2] = 0.0;

        // outer particle-particle collision
        for (int n = 1; n<ca->sizeof_X; n++)
          {
            for (int outer_slot=0; outer_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; outer_slot++)
              if (calGetX3Di(ca, Q.ID[outer_slot],cell_x,cell_y,cell_z,n) > NULL_ID)
                {
                  pj[0] = calGetX3Dr(ca,Q.px[outer_slot],cell_x,cell_y,cell_z,n);
                  pj[1] = calGetX3Dr(ca,Q.py[outer_slot],cell_x,cell_y,cell_z,n);
                  pj[2] = calGetX3Dr(ca,Q.pz[outer_slot],cell_x,cell_y,cell_z,n);
#ifdef VISCOELASTIC
                  vj[0] = calGetX3Dr(ca,Q.vx[outer_slot],cell_x,cell_y,cell_z,n);
                  vj[1] = calGetX3Dr(ca,Q.vy[outer_slot],cell_x,cell_y,cell_z,n);
                  vj[2] = calGetX3Dr(ca,Q.vz[outer_slot],cell_x,cell_y,cell_z,n);
#endif
                  dij = distance(pi, pj);
                  if (dij < 2*PARTICLE_RADIUS)
                    {
                      for (int k=0; k<3; k++)
                        {
                          rij[k] = pj[k] - pi[k];
                          enij[k] = rij[k] / dij;
                          delta_n = 2*PARTICLE_RADIUS - dij;
                          delta_Fi[k] += -kn * delta_n * enij[k];
#ifdef VISCOELASTIC
                          vij[k] = vi[k] - vj[k];
#endif
                        }
#ifdef VISCOELASTIC
                      vnij = 0.0;
                      for (int k=0; k<3; k++)
                        vnij += vij[k] * enij[k];

                      for (int k=0; k<3; k++)
                        delta_Fi[k] += -etha * vnij * enij[k];
#endif
                    } //if dij < 2*PARTICLE_RADIUS
                } // for outer_slot AND if
          } //for n = 1 ...

        // update phase
        calSet3Dr(ca, Q.Fx[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z) + delta_Fi[0]);
        calSet3Dr(ca, Q.Fy[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z) + delta_Fi[1]);
        calSet3Dr(ca, Q.Fz[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z) + delta_Fi[2]);
      }
}

void boundary_collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal kn = KN;
  CALreal etha = ETHA;
  CALreal pi[3], vi[3], delta_Fi[3];
  CALreal pj[3], vj[3], Nj[3];
  CALreal rij[3], dij, enij[3], vij[3], vnij;
  CALreal delta_n;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.ID[slot],cell_x,cell_y,cell_z) > NULL_ID)
      {
        pi[0] = calGet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z);
        pi[1] = calGet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z);
        pi[2] = calGet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z);
#ifdef VISCOELASTIC
        vi[0] = calGet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z);
        vi[1] = calGet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z);
        vi[2] = calGet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z);
#endif

        delta_Fi[0] = 0.0;
        delta_Fi[1] = 0.0;
        delta_Fi[2] = 0.0;
        for (int n=1; n<ca->sizeof_X; n++)
          {
            if (calGetX3Di(ca, Q.ID[0], cell_x,cell_y, cell_z, n) == BORDER_ID)
              {
                pj[0] = calGetX3Dr(ca,Q.px[0],cell_x,cell_y,cell_z,n);
                pj[1] = calGetX3Dr(ca,Q.py[0],cell_x,cell_y,cell_z,n);
                pj[2] = calGetX3Dr(ca,Q.pz[0],cell_x,cell_y,cell_z,n);

                Nj[0] = calGetX3Dr(ca,Q.vx[0],cell_x,cell_y,cell_z,n);
                Nj[1] = calGetX3Dr(ca,Q.vy[0],cell_x,cell_y,cell_z,n);
                Nj[2] = calGetX3Dr(ca,Q.vz[0],cell_x,cell_y,cell_z,n);

#ifdef VISCOELASTIC
                vj[0] = 0.0;
                vj[1] = 0.0;
                vj[2] = 0.0;
#endif
                dij = pointPlaneDistance(pi, Nj, pj);
                if (dij < PARTICLE_RADIUS)
                  {
                    orthogonalProjectedPointToPlane(pi, pj, Nj, pj);
                    for (int k=0; k<3; k++)
                      {
                        rij[k] = pj[k] - pi[k];
                        enij[k] = rij[k] / dij;
                        delta_n = PARTICLE_RADIUS - dij;
                        delta_Fi[k] += -kn * delta_n * enij[k];
#ifdef VISCOELASTIC
                        vij[k] = vi[k] - vj[k];
#endif
                      }
#ifdef VISCOELASTIC
                    vnij = 0.0;
                    for (int k=0; k<3; k++)
                      vnij += vij[k] * enij[k];

                    for (int k=0; k<3; k++)
                      delta_Fi[k] += -etha * vnij * enij[k];
#endif
                  } //if dij < PARTICLE_RADIUS
              } // if BORDER
          } //for n = 1 ...

        // update phase
        calSet3Dr(ca, Q.Fx[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fx[slot], cell_x,cell_y,cell_z) + delta_Fi[0]);
        calSet3Dr(ca, Q.Fy[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fy[slot], cell_x,cell_y,cell_z) + delta_Fi[1]);
        calSet3Dr(ca, Q.Fz[slot],cell_x,cell_y,cell_z,calGetNext3Dr(ca, Q.Fz[slot], cell_x,cell_y,cell_z) + delta_Fi[2]);
      }
}
