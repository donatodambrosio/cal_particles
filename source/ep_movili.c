#include <ep_movili.h>
#include <ep_physics.h>
#include <ep_utils.h>

void movili(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  if (!ncestiArmenuNaParticella(ca, cell_x, cell_y, cell_z, 0))
    return;

  CALreal F0[3];
  CALreal p0[3];
  CALreal v0[3];
  CALreal pf[3];
  CALreal vf[3];

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.ID[slot],cell_x,cell_y,cell_z) > NULL_ID)
      {
        F0[0] = calGet3Dr(ca,Q.Fx[slot],cell_x,cell_y,cell_z);
        F0[1] = calGet3Dr(ca,Q.Fy[slot],cell_x,cell_y,cell_z);
        F0[2] = calGet3Dr(ca,Q.Fz[slot],cell_x,cell_y,cell_z);

        p0[0] = pf[0] = calGet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z);
        p0[1] = pf[1] = calGet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z);
        p0[2] = pf[2] = calGet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z);

        v0[0] = vf[0] = calGet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z);
        v0[1] = vf[1] = calGet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z);
        v0[2] = vf[2] = calGet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z);

        applyForce(F0, p0, v0, PARTICLE_MASS, DELTA_T, pf, vf);

        calSet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z,pf[0]);
        calSet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z,pf[1]);
        calSet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z,pf[2]);

        calSet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z,vf[0]);
        calSet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z,vf[1]);
        calSet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z,vf[2]);
      }
}
