#include <ep_movili_cazzu.h>
#include <ep_utils.h>

void pezziala(int slot, struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  calSet3Dr(ca, Q.Fx[slot],   cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca, Q.Fy[slot],   cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca, Q.Fz[slot],   cell_x,cell_y,cell_z,0.0);
  calSet3Dr(ca, Q.px[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
  calSet3Dr(ca, Q.py[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
  calSet3Dr(ca, Q.pz[slot],   cell_x,cell_y,cell_z,PARTICLE_NODATA);
  calSet3Dr(ca, Q.vx[slot],   cell_x,cell_y,cell_z,0);
  calSet3Dr(ca, Q.vy[slot],   cell_x,cell_y,cell_z,0);
  calSet3Dr(ca, Q.vz[slot],   cell_x,cell_y,cell_z,0);
  calSet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,PARTICLE_ABSENT);
}

void sucala(int destination_slot, int source_slot, struct CALModel3D* ca, int cell_x, int cell_y, int cell_z, int n)
{
  calSet3Dr(ca, Q.Fx[destination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.Fx[source_slot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.Fy[destination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.Fy[source_slot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.Fz[destination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.Fz[source_slot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.px[destination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.px[source_slot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.py[destination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.py[source_slot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.pz[destination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.pz[source_slot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.vx[destination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.vx[source_slot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.vy[destination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.vy[source_slot],   cell_x,cell_y,cell_z,n));
  calSet3Dr(ca, Q.vz[destination_slot],   cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.vz[source_slot],   cell_x,cell_y,cell_z,n));
  calSet3Di(ca, Q.imove[destination_slot],cell_x,cell_y,cell_z,calGetX3Di(ca,Q.imove[source_slot],cell_x,cell_y,cell_z,n));
}

void moviliCazzu(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal px, py, pz;
  CALint  new_cell_x, new_cell_y, new_cell_z;

  if (calGet3Di(ca, Q.imove[0],cell_x,cell_y,cell_z) == PARTICLE_BORDER)
    return;

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
              {
                CALbyte sucked = CAL_FALSE;
                int destination_slot;
                for (destination_slot = 0; destination_slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; destination_slot++)
                  if (calGetNext3Di(ca,Q.imove[destination_slot],cell_x,cell_y,cell_z) == PARTICLE_ABSENT)
                    {
                      sucala(destination_slot,source_slot,ca,cell_x,cell_y,cell_z,n);
                      sucked = CAL_TRUE;
                      break;
                    }
                if (sucked == CAL_FALSE)
                  {
                    CALint capacity = MAX_NUMBER_OF_PARTICLES_PER_CELL;
                    CALint dest_slot = destination_slot;
                    break;
                  }
              }
          }
}

//void moviliCazzuSucali(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
//{
//  CALreal px, py, pz;
//  CALint  new_cell_x, new_cell_y, new_cell_z;

///*
//  if (ncestiArmenuNaParticella(ca, cell_x, cell_y, cell_z, 0))
//    {
//      //pezziali
//      for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
//        if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
//          {
//            px = calGet3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z);
//            py = calGet3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z);
//            pz = calGet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z);

//            new_cell_x = px/CELL_SIDE;
//            new_cell_y = py/CELL_SIDE;
//            new_cell_z = pz/CELL_SIDE;

//            if ((cell_x != new_cell_x) || (cell_y != new_cell_y) || (cell_z != new_cell_z))
//                pezziala(slot, ca,cell_x,cell_y,cell_z);
//          }
//    }
//*/

//  //sucali
//  for (int n=1; n<ca->sizeof_X; n++)
//    if (ncestiArmenuNaParticella(ca, cell_x, cell_y, cell_z, n))
//      for (int source_slot = 0; source_slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; source_slot++)
//        if (calGetX3Di(ca,Q.imove[source_slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
//          {
//            px = calGetX3Dr(ca, Q.px[source_slot],cell_x,cell_y,cell_z,n);
//            py = calGetX3Dr(ca, Q.py[source_slot],cell_x,cell_y,cell_z,n);
//            pz = calGetX3Dr(ca, Q.pz[source_slot],cell_x,cell_y,cell_z,n);

//            new_cell_x = px/CELL_SIDE;
//            new_cell_y = py/CELL_SIDE;
//            new_cell_z = pz/CELL_SIDE;

//            if ((cell_x == new_cell_x) && (cell_y == new_cell_y) && (cell_z == new_cell_z))
//              {
//                CALbyte sucked = CAL_FALSE;
//                int destination_slot;
//                for (destination_slot = 0; destination_slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; destination_slot++)
//                  if (calGetNext3Di(ca,Q.imove[destination_slot],cell_x,cell_y,cell_z) == PARTICLE_ABSENT)
//                    {
//                      sucala(destination_slot,source_slot,ca,cell_x,cell_y,cell_z,n);
//                      sucked = CAL_TRUE;
//                      break;
//                    }
//                  else
//                    {
//                      CALreal p0[3], pn[3];
//                      pn[0] = px;
//                      pn[1] = py;
//                      pn[2] = pz;
//                      for (int i = 0; i < MAX_NUMBER_OF_PARTICLES_PER_CELL; i++)
//                        if (calGet3Di(ca, Q.imove[i],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
//                          {
//                            p0[0] = calGet3Dr(ca, Q.px[i],cell_x,cell_y,cell_z);
//                            p0[1] = calGet3Dr(ca, Q.py[i],cell_x,cell_y,cell_z);
//                            p0[2] = calGet3Dr(ca, Q.pz[i],cell_x,cell_y,cell_z);

//                            CALreal d = distance(p0,pn);

//                            CALbyte collision = CAL_FALSE;
//                            if (d < 2*PARTICLE_RADIUS)
//                              collision = CAL_TRUE;
//                          }
//                    }

//                if (sucked == CAL_FALSE)
//                  {
//                    CALint current_step = a_simulazioni->step;
//                    CALint capacity = MAX_NUMBER_OF_PARTICLES_PER_CELL;
//                    CALint src_slot = source_slot;
//                    CALint dest_slot = destination_slot;
//                    break;
//                  }
//              }
//          }
//}
