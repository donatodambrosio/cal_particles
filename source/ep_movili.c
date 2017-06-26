#include <ep_movili.h>
#include <ep_physics.h>
#include <ep_utils.h>

void movili(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  if (!ncestiArmenuNaParticella(ca, cell_x, cell_y, cell_z, 0))
    return;

  CALreal F0[3];
  CALreal r0[3];
  CALreal v0[3];
  CALreal rf[3];
  CALreal vf[3];

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.ID[slot],cell_x,cell_y,cell_z) > NULL_ID)
      {
        F0[0] = calGet3Dr(ca,Q.Fx[slot],cell_x,cell_y,cell_z);
        F0[1] = calGet3Dr(ca,Q.Fy[slot],cell_x,cell_y,cell_z);
        F0[2] = calGet3Dr(ca,Q.Fz[slot],cell_x,cell_y,cell_z);

        r0[0] = rf[0] = calGet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z);
        r0[1] = rf[1] = calGet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z);
        r0[2] = rf[2] = calGet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z);

        v0[0] = vf[0] = calGet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z);
        v0[1] = vf[1] = calGet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z);
        v0[2] = vf[2] = calGet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z);

        applyForce(F0, r0, v0, PARTICLE_MASS, DELTA_T, rf, vf);

        calSet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z,rf[0]);
        calSet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z,rf[1]);
        calSet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z,rf[2]);

        calSet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z,vf[0]);
        calSet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z,vf[1]);
        calSet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z,vf[2]);
      }
}
