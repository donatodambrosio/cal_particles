#include <ep_init.h>
#include <ep_utils.h>
#include <stdlib.h>

unsigned int seed = 1;

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

      calInit3Dr(ca,Q.Fx[0],cell_x,cell_y,cell_z,PARTICLE_MASS*G);
      calInit3Dr(ca,Q.Fy[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca,Q.Fz[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca,Q.px[0],cell_x,cell_y,cell_z,px);
      calInit3Dr(ca,Q.py[0],cell_x,cell_y,cell_z,py);
      calInit3Dr(ca,Q.pz[0],cell_x,cell_y,cell_z,pz);
      calInit3Dr(ca,Q.vx[0],cell_x,cell_y,cell_z, v);
      calInit3Dr(ca,Q.vy[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca,Q.vz[0],cell_x,cell_y,cell_z,0.0);
      calInit3Di(ca,Q.ID[0],cell_x,cell_y,cell_z,  DEFAULT_PARTICLE_ID);
  }

  if (cell_x == 3 && cell_y == 2 && cell_z == TOP_LAYERS - 4)
  {
      CALreal px = CELL_SIDE * (cell_x);
      CALreal py = CELL_SIDE * (cell_y);
      CALreal pz = CELL_SIDE * (cell_z);

      calInit3Dr(ca,Q.Fx[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca,Q.Fy[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca,Q.Fz[0],cell_x,cell_y,cell_z,-PARTICLE_MASS*G);
      calInit3Dr(ca,Q.px[0],cell_x,cell_y,cell_z,px);
      calInit3Dr(ca,Q.py[0],cell_x,cell_y,cell_z,py);
      calInit3Dr(ca,Q.pz[0],cell_x,cell_y,cell_z,pz);
      calInit3Dr(ca,Q.vx[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca,Q.vy[0],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca,Q.vz[0],cell_x,cell_y,cell_z,0.0);
      calInit3Di(ca,Q.ID[0],cell_x,cell_y,cell_z,  DEFAULT_PARTICLE_ID);
  }
  return;
*/
  //if (cell_x < TOP_LAYERS || calGet3Di(ca, Q.ID[0],cell_x,cell_y,cell_z) == BORDER)
  //if (cell_y < TOP_LAYERS || calGet3Di(ca, Q.ID[0],cell_x,cell_y,cell_z) == BORDER)
  if (cell_z > TOP_LAYERS || calGet3Di(ca, Q.ID[0],cell_x,cell_y,cell_z) == BORDER_ID)
  //if (cell_x < TOP_LAYERS || cell_y < TOP_LAYERS || cell_z < TOP_LAYERS || calGet3Di(ca, Q.ID[0],cell_x,cell_y,cell_z) == BORDER)
    return;

  CALreal c;
  for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      c = (CALreal)rand_r(&seed)/(CALreal)(RAND_MAX); // 0 <= c <= 1
      if (c < CELL_FILL_RATE)
        {
          CALreal cx = (CALreal)rand_r(&seed)/(CALreal)(RAND_MAX); // 0 <= c <= 1
          CALreal cy = (CALreal)rand_r(&seed)/(CALreal)(RAND_MAX); // 0 <= c <= 1
          CALreal cz = (CALreal)rand_r(&seed)/(CALreal)(RAND_MAX); // 0 <= c <= 1

          CALreal px = CELL_SIDE * (cell_x + cx);
          CALreal py = CELL_SIDE * (cell_y + cy);
          CALreal pz = CELL_SIDE * (cell_z + cz);

          calInit3Dr(ca,Q.Fx[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Dr(ca,Q.Fy[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Dr(ca,Q.Fz[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z,px);
          calInit3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z,py);
          calInit3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z,pz);
          calInit3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z,0.0);
          calInit3Di(ca,Q.ID[slot],cell_x,cell_y,cell_z,DEFAULT_PARTICLE_ID);
        }
    }
}

void mmiscali_nta_cella_seriale(struct CALModel3D* ca)
{
  int i, j, k;
  for (i = 0; i < ca->rows; i++)
    for (j = 0; j < ca->columns; j++)
      for (k = 0; k < ca->slices; k++)
        mmiscali_nta_cella(ca, i, j, k);
}

void pezzialaMo(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z, int slot)
{
  calInit3Dr(ca,Q.Fx[slot],cell_x,cell_y,cell_z,0.0);
  calInit3Dr(ca,Q.Fy[slot],cell_x,cell_y,cell_z,0.0);
  calInit3Dr(ca,Q.Fz[slot],cell_x,cell_y,cell_z,0.0);
  calInit3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z,0.0);
  calInit3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z,0.0);
  calInit3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z,0.0);
  calInit3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z,0.0);
  calInit3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z,0.0);
  calInit3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z,0.0);
  calInit3Di(ca,Q.ID[slot],cell_x,cell_y,cell_z,NULL_ID);
}

void cancella_particelle_in_urto(struct CALModel3D* ca)
{
  CALreal p[3], pj[3], pn[3], Nn[3];
  CALbyte particle_OK = CAL_TRUE;
  CALint particle_id = 0;

  for (int cell_x=0; cell_x<ca->rows; cell_x++)
    for (int cell_y=0; cell_y<ca->columns; cell_y++)
      for (int cell_z = 0; cell_z<ca->slices; cell_z++)
        for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
          if (calGet3Di(ca,Q.ID[slot],cell_x,cell_y,cell_z) > NULL_ID)
            {
              particle_OK = CAL_TRUE;
              p[0] = calGet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z);
              p[1] = calGet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z);
              p[2] = calGet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z);

              for (int inner_slot=slot+1; inner_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; inner_slot++)
                if (calGet3Di(ca, Q.ID[inner_slot],cell_x,cell_y,cell_z) > NULL_ID)
                  {
                    pj[0] = calGet3Dr(ca,Q.px[inner_slot],cell_x,cell_y,cell_z);
                    pj[1] = calGet3Dr(ca,Q.py[inner_slot],cell_x,cell_y,cell_z);
                    pj[2] = calGet3Dr(ca,Q.pz[inner_slot],cell_x,cell_y,cell_z);

                    if (distance(p, pj) <= 2.0*PARTICLE_RADIUS)
                      {
                        pezzialaMo(ca,cell_x,cell_y,cell_z,slot);
                        particle_OK = CAL_FALSE;
                      }
                  }

              for (int n = 1; n<ca->sizeof_X; n++)
                if (calGetX3Di(ca, Q.ID[0],cell_x,cell_y,cell_z,n) == BORDER_ID)
                  {
                    pn[0] = calGetX3Dr(ca,Q.px[0],cell_x,cell_y,cell_z,n);
                    pn[1] = calGetX3Dr(ca,Q.py[0],cell_x,cell_y,cell_z,n);
                    pn[2] = calGetX3Dr(ca,Q.pz[0],cell_x,cell_y,cell_z,n);

                    Nn[0] = calGetX3Dr(ca,Q.vx[0],cell_x,cell_y,cell_z,n);
                    Nn[1] = calGetX3Dr(ca,Q.vy[0],cell_x,cell_y,cell_z,n);
                    Nn[2] = calGetX3Dr(ca,Q.vz[0],cell_x,cell_y,cell_z,n);

                    if (pointPlaneDistance(p, pn, Nn) < PARTICLE_RADIUS)
                      {
                        pezzialaMo(ca,cell_x,cell_y,cell_z,slot);
                        particle_OK = CAL_FALSE;
                      }
                  }
                else
                  for (int outer_slot=0; outer_slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; outer_slot++)
                    if (calGetX3Di(ca, Q.ID[outer_slot],cell_x,cell_y,cell_z,n) > NULL_ID)
                      {
                        pn[0] = calGetX3Dr(ca,Q.px[outer_slot],cell_x,cell_y,cell_z,n);
                        pn[1] = calGetX3Dr(ca,Q.py[outer_slot],cell_x,cell_y,cell_z,n);
                        pn[2] = calGetX3Dr(ca,Q.pz[outer_slot],cell_x,cell_y,cell_z,n);

                        if (distance(p, pn) <= 2.0*PARTICLE_RADIUS)
                          {
                            pezzialaMo(ca,cell_x,cell_y,cell_z,slot);
                            particle_OK = CAL_FALSE;
                          }
                      }

              if (particle_OK)
                {
                  particle_id++;
                  calInit3Di(ca,Q.ID[slot],cell_x,cell_y,cell_z,particle_id);
                }
            }

  initial_nummber_of_particles = particle_id;

  return;
}
