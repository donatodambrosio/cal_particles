#include <ep_movili_cazzu.h>
#include <ep_utils.h>

void pezziala(int slot, struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  calSet3Dr(ca, Q.px[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
  calSet3Dr(ca, Q.py[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
  calSet3Dr(ca, Q.pz[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
  calSet3Dr(ca, Q.vx[slot],   cell_x,cell_y,cell_z,0);
  calSet3Dr(ca, Q.vy[slot],   cell_x,cell_y,cell_z,0);
  calSet3Dr(ca, Q.vz[slot],   cell_x,cell_y,cell_z,0);
  calSet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,PARTICLE_ABSENT);
}

void sucala(int detination_slot, int source_sloot, struct CALModel3D* ca, int cell_x, int cell_y, int cell_z, int n)
{
  calSet3Dr(ca, Q.px[detination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.px[source_sloot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.py[detination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.py[source_sloot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.pz[detination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.pz[source_sloot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.vx[detination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.vx[source_sloot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.vy[detination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.vy[source_sloot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.vz[detination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.vz[source_sloot],   cell_x,cell_y,cell_z,n));
  calSet3Di(ca, Q.imove[detination_slot],cell_x,cell_y,cell_z,calGetX3Di(ca,Q.imove[source_sloot],cell_x,cell_y,cell_z,n));
}

void moviliCazzu(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALint a = 0;
  CALreal px, py, pz;
  CALint  new_cell_x, new_cell_y, new_cell_z;

  if (ncestiArmenuNaParticella(ca, cell_x, cell_y, cell_z, 0))
    {
      //pezziali
      for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
        if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
          {
            px = calGet3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z);
            py = calGet3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z);
            pz = calGet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z);

            new_cell_x = px/CELL_SIDE;
            new_cell_y = py/CELL_SIDE;
            new_cell_z = pz/CELL_SIDE;

            if ((cell_x != new_cell_x) || (cell_y != new_cell_y) || (cell_z != new_cell_z))
                pezziala(slot, ca,cell_x,cell_y,cell_z);
          }
    }

  //sucali
  for (int n=1; n<ca->sizeof_X; n++)
    if (ncestiArmenuNaParticella(ca, cell_x, cell_y, cell_z, n))
      for (int source_slot = 0; source_slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; source_slot++)
        if (calGetX3Di(ca,Q.imove[source_slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
          {
            px = calGetX3Dr(ca, Q.px[source_slot],cell_x,cell_y,cell_z,n);
            py = calGetX3Dr(ca, Q.py[source_slot],cell_x,cell_y,cell_z,n);
            pz = calGetX3Dr(ca, Q.pz[source_slot],cell_x,cell_y,cell_z,n);

            new_cell_x = px/CELL_SIDE;
            new_cell_y = py/CELL_SIDE;
            new_cell_z = pz/CELL_SIDE;

            if ((cell_x == new_cell_x) && (cell_y == new_cell_y) && (cell_z == new_cell_z))
              for (int destination_slot = 0; destination_slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; destination_slot++)
                if (calGetNext3Di(ca,Q.imove[destination_slot],cell_x,cell_y,cell_z) == PARTICLE_ABSENT)
                  {
                    sucala(destination_slot,source_slot,ca,cell_x,cell_y,cell_z,n);
                    break;
                  }
          }
}
