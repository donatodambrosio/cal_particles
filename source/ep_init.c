#include <ep_init.h>
#include <ep_utils.h>

void mmiscali_nta_cella(struct CALModel3D* ca, int x_cell, int y_cell, int z_cell)
{
  //if (x_cell < TOP_LAYERS || calGet3Di(ca, Q.imove[0],x_cell,y_cell,z_cell) == PARTICLE_BORDER)
  //if (y_cell < TOP_LAYERS || calGet3Di(ca, Q.imove[0],x_cell,y_cell,z_cell) == PARTICLE_BORDER)
  //if (z_cell < TOP_LAYERS || calGet3Di(ca, Q.imove[0],x_cell,y_cell,z_cell) == PARTICLE_BORDER)
  if (x_cell < TOP_LAYERS || z_cell < TOP_LAYERS || calGet3Di(ca, Q.imove[0],x_cell,y_cell,z_cell) == PARTICLE_BORDER)
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

          CALreal px = CELL_SIDE * (x_cell + cx);
          CALreal py = CELL_SIDE * (y_cell + cy);
          CALreal pz = CELL_SIDE * (z_cell + cz);

          calInit3Dr(ca, Q.px[slot],x_cell,y_cell,z_cell,px);
          calInit3Dr(ca, Q.py[slot],x_cell,y_cell,z_cell,py);
          calInit3Dr(ca, Q.pz[slot],x_cell,y_cell,z_cell,pz);
          calInit3Dr(ca, Q.vx[slot],x_cell,y_cell,z_cell,0.0);
          calInit3Dr(ca, Q.vy[slot],x_cell,y_cell,z_cell,0.0);
          calInit3Dr(ca, Q.vz[slot],x_cell,y_cell,z_cell,0.0);
          calInit3Di(ca, Q.imove[slot],x_cell,y_cell,z_cell,PARTICLE_PRESENT);
        }
    }
}


// ATTENZIONE, QUESTO PROCESSO ELEMENTARE DEVE AVVENIRE IN SEQUENZIALE
void cancella_particelle_in_urto(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal p0[3], pi[3], pn[3];

  for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca,Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
      {
        p0[0] = calGet3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z);
        p0[1] = calGet3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z);
        p0[2] = calGet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z);

        for (int inner_slot=slot+1; inner_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; inner_slot++)
          {
            pi[0] = calGet3Dr(ca, Q.px[inner_slot],cell_x,cell_y,cell_z);
            pi[1] = calGet3Dr(ca, Q.py[inner_slot],cell_x,cell_y,cell_z);
            pi[2] = calGet3Dr(ca, Q.pz[inner_slot],cell_x,cell_y,cell_z);

            if (distance(p0, pi) < 2*PARTICLE_RADIUS)
              {
                //printf("particle (%d,%d,%d,%d) is in inner collision state\n", cell_x, cell_y, cell_z, slot);
                calInit3Dr(ca, Q.px[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                calInit3Dr(ca, Q.py[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                calInit3Dr(ca, Q.pz[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                calInit3Dr(ca, Q.vx[slot],   cell_x,cell_y,cell_z,0);
                calInit3Dr(ca, Q.vy[slot],   cell_x,cell_y,cell_z,0);
                calInit3Dr(ca, Q.vz[slot],   cell_x,cell_y,cell_z,0);
                calInit3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,PARTICLE_ABSENT);
              }
          }
        for (int n = 1; n<ca->sizeof_X; n++)
          for (int outher_slot=0; outher_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; outher_slot++)
            {
              pn[0] = calGetX3Dr(ca, Q.px[outher_slot],cell_x,cell_y,cell_z,n);
              pn[1] = calGetX3Dr(ca, Q.py[outher_slot],cell_x,cell_y,cell_z,n);
              pn[2] = calGetX3Dr(ca, Q.pz[outher_slot],cell_x,cell_y,cell_z,n);

              if (distance(p0, pn) < 2*PARTICLE_RADIUS)
                {
                  //printf("particle (%d,%d,%d,%d) is in outher collision state\n", cell_x, cell_y, cell_z, slot);
                  calInit3Dr(ca, Q.px[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                  calInit3Dr(ca, Q.py[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                  calInit3Dr(ca, Q.pz[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
                  calInit3Dr(ca, Q.vx[slot],   cell_x,cell_y,cell_z,0);
                  calInit3Dr(ca, Q.vy[slot],   cell_x,cell_y,cell_z,0);
                  calInit3Dr(ca, Q.vz[slot],   cell_x,cell_y,cell_z,0);
                  calInit3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,PARTICLE_ABSENT);
                }
            }
      }
}
