#include <ep_init.h>
#include <ep_utils.h>
#include <stdlib.h>

void mmiscali_nta_cella(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  /*
  CALreal v = 0.1;
  //if (cell_x == 2 && cell_y == 2 && (cell_z == TOP_LAYERS - 1 || cell_z == TOP_LAYERS - 2))
  if (cell_x == 2 && cell_y == 2 && cell_z == TOP_LAYERS - 4)
  {
      CALreal px = CELL_SIDE * (cell_x);
      CALreal py = CELL_SIDE * (cell_y);
      CALreal pz = CELL_SIDE * (cell_z);

      calInit3Dr(ca, Q.Fx[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q.Fy[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q.Fz[0],cell_x,cell_y,cell_z,-PARTICLE_MASS*G);
      calInit3Dr(ca, Q.rx[0],cell_x,cell_y,cell_z,px);
      calInit3Dr(ca, Q.ry[0],cell_x,cell_y,cell_z,py);
      calInit3Dr(ca, Q.rz[0],cell_x,cell_y,cell_z,pz);
      calInit3Dr(ca, Q.vx[0],cell_x,cell_y,cell_z, v);
      calInit3Dr(ca, Q.vy[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q.vz[0],cell_x,cell_y,cell_z,0.0);
      calInit3Di(ca, Q.imove[0],cell_x,cell_y,cell_z,PARTICLE_PRESENT);
  }

  if (cell_x == 3 && cell_y == 2 && cell_z == TOP_LAYERS - 4)
  {
      CALreal px = CELL_SIDE * (cell_x);
      CALreal py = CELL_SIDE * (cell_y);
      CALreal pz = CELL_SIDE * (cell_z);

      calInit3Dr(ca, Q.Fx[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q.Fy[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q.Fz[0],cell_x,cell_y,cell_z,-PARTICLE_MASS*G);
      calInit3Dr(ca, Q.rx[0],cell_x,cell_y,cell_z,px);
      calInit3Dr(ca, Q.ry[0],cell_x,cell_y,cell_z,py);
      calInit3Dr(ca, Q.rz[0],cell_x,cell_y,cell_z,pz);
      calInit3Dr(ca, Q.vx[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q.vy[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q.vz[0],cell_x,cell_y,cell_z,0.0);
      calInit3Di(ca, Q.imove[0],cell_x,cell_y,cell_z,PARTICLE_PRESENT);
  }
  return;
  */

  //if (cell_x < TOP_LAYERS || calGet3Di(ca, Q.imove[0],cell_x,cell_y,cell_z) == PARTICLE_BORDER)
  //if (cell_y < TOP_LAYERS || calGet3Di(ca, Q.imove[0],cell_x,cell_y,cell_z) == PARTICLE_BORDER)
  if (cell_z < TOP_LAYERS || calGet3Di(ca, Q.imove[0],cell_x,cell_y,cell_z) == PARTICLE_BORDER)
  //if (cell_x < TOP_LAYERS || cell_y < TOP_LAYERS || cell_z < TOP_LAYERS || calGet3Di(ca, Q.imove[0],cell_x,cell_y,cell_z) == PARTICLE_BORDER)
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
          calInit3Dr(ca, Q.Fz[slot],cell_x,cell_y,cell_z,0.0);
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
          if (calGetX3Di(ca, Q.imove[0],cell_x,cell_y,cell_z,n) == PARTICLE_BORDER)
            {
              pn[0] = calGetX3Dr(ca, Q.rx[0],cell_x,cell_y,cell_z,n);
              pn[1] = calGetX3Dr(ca, Q.ry[0],cell_x,cell_y,cell_z,n);
              pn[2] = calGetX3Dr(ca, Q.rz[0],cell_x,cell_y,cell_z,n);

              if (distance(p, pn) < 3*PARTICLE_RADIUS)
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
          else
            for (int outer_slot=0; outer_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; outer_slot++)
              if (calGetX3Di(ca, Q.imove[outer_slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
                {
                  pn[0] = calGetX3Dr(ca, Q.rx[outer_slot],cell_x,cell_y,cell_z,n);
                  pn[1] = calGetX3Dr(ca, Q.ry[outer_slot],cell_x,cell_y,cell_z,n);
                  pn[2] = calGetX3Dr(ca, Q.rz[outer_slot],cell_x,cell_y,cell_z,n);

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
