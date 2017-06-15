#include <ep_init.h>
#include <ep_utils.h>
#include <stdlib.h>

void mmiscali_nta_cella(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  //if (x_cell < TOP_LAYERS || calGet3Di(ca, Q.imove[0],x_cell,y_cell,z_cell) == PARTICLE_BORDER)
  //if (y_cell < TOP_LAYERS || calGet3Di(ca, Q.imove[0],x_cell,y_cell,z_cell) == PARTICLE_BORDER)
  //if (z_cell < TOP_LAYERS || calGet3Di(ca, Q.imove[0],x_cell,y_cell,z_cell) == PARTICLE_BORDER)
  if (cell_x < TOP_LAYERS || cell_z < TOP_LAYERS || calGet3Di(ca, Q.imove[0],cell_x,cell_y,cell_z) == PARTICLE_BORDER)
    return;

  CALreal c;
  for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      c = (CALreal)rand()/(CALreal)(RAND_MAX); // 0 <= c <= 1
      if (c < CELL_FILL_RATE)
        {
          CALreal cx = (CALreal)rand()/(CALreal)(RAND_MAX); // 0 <= c <= 1
          CALreal cy = (CALreal)rand()/(CALreal)(RAND_MAX); // 0 <= c <= 1
          CALreal cz = (CALreal)rand()/(CALreal)(RAND_MAX); // 0 <= c <= 1

          CALreal px = CELL_SIDE * (cell_x + cx);
          CALreal py = CELL_SIDE * (cell_y + cy);
          CALreal pz = CELL_SIDE * (cell_z + cz);

          calInit3Dr(ca, Q.Fx[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Dr(ca, Q.Fy[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Dr(ca, Q.Fz[slot],cell_x,cell_y,cell_z,-PARTICLE_MASS*G);
          calInit3Dr(ca, Q.rx[slot],cell_x,cell_y,cell_z,px);
          calInit3Dr(ca, Q.ry[slot],cell_x,cell_y,cell_z,py);
          calInit3Dr(ca, Q.rz[slot],cell_x,cell_y,cell_z,pz);
          calInit3Dr(ca, Q.vx[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Dr(ca, Q.vy[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Dr(ca, Q.vz[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,PARTICLE_PRESENT);
        }
    }
}


// ATTENZIONE, QUESTO PROCESSO ELEMENTARE DEVE AVVENIRE IN SEQUENZIALE
void cancella_particelle_in_urto(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal p[3], pi[3], pn[3];

  for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca,Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
      {
        p[0] = calGet3Dr(ca, Q.rx[slot],cell_x,cell_y,cell_z);
        p[1] = calGet3Dr(ca, Q.ry[slot],cell_x,cell_y,cell_z);
        p[2] = calGet3Dr(ca, Q.rz[slot],cell_x,cell_y,cell_z);

        for (int inner_slot=slot+1; inner_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; inner_slot++)
          if (calGet3Di(ca, Q.imove[inner_slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
          {
            pi[0] = calGet3Dr(ca, Q.rx[inner_slot],cell_x,cell_y,cell_z);
            pi[1] = calGet3Dr(ca, Q.ry[inner_slot],cell_x,cell_y,cell_z);
            pi[2] = calGet3Dr(ca, Q.rz[inner_slot],cell_x,cell_y,cell_z);

            if (distance(p, pi) < 2*PARTICLE_RADIUS)
              {
                //printf("particle (%d,%d,%d,%d) is in inner collision state\n", cell_x, cell_y, cell_z, slot);
                calInit3Dr(ca, Q.Fx[slot],   cell_x,cell_y,cell_z,0.0);
                calInit3Dr(ca, Q.Fy[slot],   cell_x,cell_y,cell_z,0.0);
                calInit3Dr(ca, Q.Fz[slot],   cell_x,cell_y,cell_z,0.0);
                calInit3Dr(ca, Q.rx[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                calInit3Dr(ca, Q.ry[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                calInit3Dr(ca, Q.rz[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                calInit3Dr(ca, Q.vx[slot],   cell_x,cell_y,cell_z,0);
                calInit3Dr(ca, Q.vy[slot],   cell_x,cell_y,cell_z,0);
                calInit3Dr(ca, Q.vz[slot],   cell_x,cell_y,cell_z,0);
                calInit3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,PARTICLE_ABSENT);
              }
          }
        for (int n = 1; n<ca->sizeof_X; n++)
          for (int other_slot=0; other_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; other_slot++)
            if (calGetX3Di(ca, Q.imove[other_slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
            {
              pn[0] = calGetX3Dr(ca, Q.rx[other_slot],cell_x,cell_y,cell_z,n);
              pn[1] = calGetX3Dr(ca, Q.ry[other_slot],cell_x,cell_y,cell_z,n);
              pn[2] = calGetX3Dr(ca, Q.rz[other_slot],cell_x,cell_y,cell_z,n);

              if (distance(p, pn) < 2*PARTICLE_RADIUS)
                {
                  //printf("particle (%d,%d,%d,%d) is in outher collision state\n", cell_x, cell_y, cell_z, slot);
                  calInit3Dr(ca, Q.Fx[slot],   cell_x,cell_y,cell_z,0.0);
                  calInit3Dr(ca, Q.Fy[slot],   cell_x,cell_y,cell_z,0.0);
                  calInit3Dr(ca, Q.Fz[slot],   cell_x,cell_y,cell_z,0.0);
                  calInit3Dr(ca, Q.rx[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                  calInit3Dr(ca, Q.ry[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                  calInit3Dr(ca, Q.rz[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                  calInit3Dr(ca, Q.vx[slot],   cell_x,cell_y,cell_z,0);
                  calInit3Dr(ca, Q.vy[slot],   cell_x,cell_y,cell_z,0);
                  calInit3Dr(ca, Q.vz[slot],   cell_x,cell_y,cell_z,0);
                  calInit3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,PARTICLE_ABSENT);
                }
            }
      }
}
