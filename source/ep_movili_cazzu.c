#include <ep_movili_cazzu.h>
#include <ep_utils.h>
#include <stdlib.h>

void pezziala(int slot, struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  calSet3Dr(ca,Q.Fx[slot],cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca,Q.Fy[slot],cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca,Q.Fz[slot],cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z,0.0);
  calSet3Di(ca,Q.ID[slot],cell_x,cell_y,cell_z,NULL_ID);
}

void sucala(int destination_slot, int source_slot, struct CALModel3D* ca, int cell_x, int cell_y, int cell_z, int n)
{
  calSet3Dr(ca,Q.Fx[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.Fx[source_slot],cell_x,cell_y,cell_z,n));
  calSet3Dr(ca,Q.Fy[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.Fy[source_slot],cell_x,cell_y,cell_z,n));
  calSet3Dr(ca,Q.Fz[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.Fz[source_slot],cell_x,cell_y,cell_z,n));
  calSet3Dr(ca,Q.px[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.px[source_slot],cell_x,cell_y,cell_z,n));
  calSet3Dr(ca,Q.py[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.py[source_slot],cell_x,cell_y,cell_z,n));
  calSet3Dr(ca,Q.pz[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.pz[source_slot],cell_x,cell_y,cell_z,n));
  calSet3Dr(ca,Q.vx[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.vx[source_slot],cell_x,cell_y,cell_z,n));
  calSet3Dr(ca,Q.vy[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.vy[source_slot],cell_x,cell_y,cell_z,n));
  calSet3Dr(ca,Q.vz[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.vz[source_slot],cell_x,cell_y,cell_z,n));
  calSet3Di(ca,Q.ID[destination_slot],cell_x,cell_y,cell_z, calGetX3Di(ca,Q.ID[source_slot],cell_x,cell_y,cell_z,n));
}

void moviliCazzu(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal px, py, pz;
  CALint  new_cell_x, new_cell_y, new_cell_z;

  //pezziali
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.ID[slot],cell_x,cell_y,cell_z) > NULL_ID)
      {
        px = calGet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z);
        py = calGet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z);
        pz = calGet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z);

        new_cell_x = px/CELL_SIDE;
        new_cell_y = py/CELL_SIDE;
        new_cell_z = pz/CELL_SIDE;

        if ((cell_x != new_cell_x) || (cell_y != new_cell_y) || (cell_z != new_cell_z))
          pezziala(slot, ca,cell_x,cell_y,cell_z);
      }

  //sucali
  for (int n=1; n<ca->sizeof_X; n++)
    for (int source_slot = 0; source_slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; source_slot++)
      if (calGetX3Di(ca,Q.ID[source_slot],cell_x,cell_y,cell_z,n) > NULL_ID)
        {
          px = calGetX3Dr(ca,Q.px[source_slot],cell_x,cell_y,cell_z,n);
          py = calGetX3Dr(ca,Q.py[source_slot],cell_x,cell_y,cell_z,n);
          pz = calGetX3Dr(ca,Q.pz[source_slot],cell_x,cell_y,cell_z,n);

          new_cell_x = px/CELL_SIDE;
          new_cell_y = py/CELL_SIDE;
          new_cell_z = pz/CELL_SIDE;

          if ((cell_x == new_cell_x) && (cell_y == new_cell_y) && (cell_z == new_cell_z))
            {
              CALbyte sucked = CAL_FALSE;
              int destination_slot;
              for (destination_slot = 0; destination_slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; destination_slot++)
                if (calGetNext3Di(ca,Q.ID[destination_slot],cell_x,cell_y,cell_z) == NULL_ID)
                  {
                    sucala(destination_slot,source_slot,ca,cell_x,cell_y,cell_z,n);
                    sucked = CAL_TRUE;
                    break;
                  }
              if (sucked == CAL_FALSE)
                {
                  printf("ERROR: unable to suck a particle.\n");
#ifdef VERBOSE
                  printf("cell_capacity: %d\n", MAX_NUMBER_OF_PARTICLES_PER_CELL);
                  printf("current_step: %d\n", a_simulazioni->step);
                  printf("source_slot: %d\n", source_slot);
                  printf("destination_slot: %d\n", destination_slot);
#endif
                  exit(EXIT_FAILURE);
                }
            }
        }
}
