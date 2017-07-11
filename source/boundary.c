#include <boundary.h>

void mb(struct CALModel3D* ca, struct Substates *Q, int cell_x, int cell_y, int cell_z)
{
  CALreal x, y, z, nx, ny, nz;

  if (cell_x == 0)
    {
      x = PARTICLE_RADIUS;
      y = (cell_y * CELL_SIDE) + CELL_SIDE/2;
      z = (cell_z * CELL_SIDE) + CELL_SIDE/2;
      nx = 1;
      ny = 0;
      nz = 0;

      for (int boundary_slot = MAX_NUMBER_OF_PARTICLES_PER_CELL - 1; boundary_slot >= 0; boundary_slot--)
        if (calGet3Di(ca, Q->ID[boundary_slot], cell_x, cell_y, cell_z) != BORDER_ID )
          {
            calInit3Dr(ca, Q->Fx[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fy[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fz[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->px[boundary_slot],cell_x,cell_y,cell_z,x);
            calInit3Dr(ca, Q->py[boundary_slot],cell_x,cell_y,cell_z,y);
            calInit3Dr(ca, Q->pz[boundary_slot],cell_x,cell_y,cell_z,z);
            calInit3Dr(ca, Q->vx[boundary_slot],cell_x,cell_y,cell_z,nx);
            calInit3Dr(ca, Q->vy[boundary_slot],cell_x,cell_y,cell_z,ny);
            calInit3Dr(ca, Q->vz[boundary_slot],cell_x,cell_y,cell_z,nz);
            calInit3Di(ca, Q->ID[boundary_slot],cell_x,cell_y,cell_z,BORDER_ID);
            break;
          }
    }

  if (cell_x == X_CELLS-1)
    {
      x = X_CELLS * CELL_SIDE - PARTICLE_RADIUS;
      y = (cell_y * CELL_SIDE) + CELL_SIDE/2;
      z = (cell_z * CELL_SIDE) + CELL_SIDE/2;
      nx = -1;
      ny = 0;
      nz = 0;

      for (int boundary_slot = MAX_NUMBER_OF_PARTICLES_PER_CELL - 1; boundary_slot >= 0; boundary_slot--)
        if (calGet3Di(ca, Q->ID[boundary_slot], cell_x, cell_y, cell_z) != BORDER_ID )
          {
            calInit3Dr(ca, Q->Fx[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fy[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fz[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->px[boundary_slot],cell_x,cell_y,cell_z,x);
            calInit3Dr(ca, Q->py[boundary_slot],cell_x,cell_y,cell_z,y);
            calInit3Dr(ca, Q->pz[boundary_slot],cell_x,cell_y,cell_z,z);
            calInit3Dr(ca, Q->vx[boundary_slot],cell_x,cell_y,cell_z,nx);
            calInit3Dr(ca, Q->vy[boundary_slot],cell_x,cell_y,cell_z,ny);
            calInit3Dr(ca, Q->vz[boundary_slot],cell_x,cell_y,cell_z,nz);
            calInit3Di(ca, Q->ID[boundary_slot],cell_x,cell_y,cell_z,BORDER_ID);
            break;
          }
    }

  if (cell_y == 0)
    {
      x = (cell_x * CELL_SIDE) + CELL_SIDE/2;
      y = PARTICLE_RADIUS;
      z = (cell_z * CELL_SIDE) + CELL_SIDE/2;
      nx = 0;
      ny = 1;
      nz = 0;

      for (int boundary_slot = MAX_NUMBER_OF_PARTICLES_PER_CELL - 1; boundary_slot >= 0; boundary_slot--)
        if (calGet3Di(ca, Q->ID[boundary_slot], cell_x, cell_y, cell_z) != BORDER_ID )
          {
            calInit3Dr(ca, Q->Fx[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fy[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fz[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->px[boundary_slot],cell_x,cell_y,cell_z,x);
            calInit3Dr(ca, Q->py[boundary_slot],cell_x,cell_y,cell_z,y);
            calInit3Dr(ca, Q->pz[boundary_slot],cell_x,cell_y,cell_z,z);
            calInit3Dr(ca, Q->vx[boundary_slot],cell_x,cell_y,cell_z,nx);
            calInit3Dr(ca, Q->vy[boundary_slot],cell_x,cell_y,cell_z,ny);
            calInit3Dr(ca, Q->vz[boundary_slot],cell_x,cell_y,cell_z,nz);
            calInit3Di(ca, Q->ID[boundary_slot],cell_x,cell_y,cell_z,BORDER_ID);
            break;
          }
    }

  if (cell_y == Y_CELLS-1)
    {
      x = (cell_x * CELL_SIDE) + CELL_SIDE/2;
      y = Y_CELLS * CELL_SIDE - PARTICLE_RADIUS;
      z = (cell_z * CELL_SIDE) + CELL_SIDE/2;
      nx = 0;
      ny = -1;
      nz = 0;

      for (int boundary_slot = MAX_NUMBER_OF_PARTICLES_PER_CELL - 1; boundary_slot >= 0; boundary_slot--)
        if (calGet3Di(ca, Q->ID[boundary_slot], cell_x, cell_y, cell_z) != BORDER_ID )
          {
            calInit3Dr(ca, Q->Fx[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fy[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fz[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->px[boundary_slot],cell_x,cell_y,cell_z,x);
            calInit3Dr(ca, Q->py[boundary_slot],cell_x,cell_y,cell_z,y);
            calInit3Dr(ca, Q->pz[boundary_slot],cell_x,cell_y,cell_z,z);
            calInit3Dr(ca, Q->vx[boundary_slot],cell_x,cell_y,cell_z,nx);
            calInit3Dr(ca, Q->vy[boundary_slot],cell_x,cell_y,cell_z,ny);
            calInit3Dr(ca, Q->vz[boundary_slot],cell_x,cell_y,cell_z,nz);
            calInit3Di(ca, Q->ID[boundary_slot],cell_x,cell_y,cell_z,BORDER_ID);
            break;
          }
    }

  if (cell_z == 0)
    {
      x = (cell_x * CELL_SIDE) + CELL_SIDE/2;
      y = (cell_y * CELL_SIDE) + CELL_SIDE/2;
      z = PARTICLE_RADIUS;
      nx = 0;
      ny = 0;
      nz = 1;

      for (int boundary_slot = MAX_NUMBER_OF_PARTICLES_PER_CELL - 1; boundary_slot >= 0; boundary_slot--)
        if (calGet3Di(ca, Q->ID[boundary_slot], cell_x, cell_y, cell_z) != BORDER_ID )
          {
            calInit3Dr(ca, Q->Fx[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fy[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fz[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->px[boundary_slot],cell_x,cell_y,cell_z,x);
            calInit3Dr(ca, Q->py[boundary_slot],cell_x,cell_y,cell_z,y);
            calInit3Dr(ca, Q->pz[boundary_slot],cell_x,cell_y,cell_z,z);
            calInit3Dr(ca, Q->vx[boundary_slot],cell_x,cell_y,cell_z,nx);
            calInit3Dr(ca, Q->vy[boundary_slot],cell_x,cell_y,cell_z,ny);
            calInit3Dr(ca, Q->vz[boundary_slot],cell_x,cell_y,cell_z,nz);
            calInit3Di(ca, Q->ID[boundary_slot],cell_x,cell_y,cell_z,BORDER_ID);
            break;
          }
    }

  if (cell_z == Z_CELLS-1)
    {
      x = (cell_x * CELL_SIDE) + CELL_SIDE/2;
      y = (cell_y * CELL_SIDE) + CELL_SIDE/2;
      z = Z_CELLS * CELL_SIDE - PARTICLE_RADIUS;
      nx = 0;
      ny = 0;
      nz = -1;

      for (int boundary_slot = MAX_NUMBER_OF_PARTICLES_PER_CELL - 1; boundary_slot >= 0; boundary_slot--)
        if (calGet3Di(ca, Q->ID[boundary_slot], cell_x, cell_y, cell_z) != BORDER_ID )
          {
            calInit3Dr(ca, Q->Fx[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fy[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->Fz[boundary_slot],cell_x,cell_y,cell_z,0.0);
            calInit3Dr(ca, Q->px[boundary_slot],cell_x,cell_y,cell_z,x);
            calInit3Dr(ca, Q->py[boundary_slot],cell_x,cell_y,cell_z,y);
            calInit3Dr(ca, Q->pz[boundary_slot],cell_x,cell_y,cell_z,z);
            calInit3Dr(ca, Q->vx[boundary_slot],cell_x,cell_y,cell_z,nx);
            calInit3Dr(ca, Q->vy[boundary_slot],cell_x,cell_y,cell_z,ny);
            calInit3Dr(ca, Q->vz[boundary_slot],cell_x,cell_y,cell_z,nz);
            calInit3Di(ca, Q->ID[boundary_slot],cell_x,cell_y,cell_z,BORDER_ID);
            break;
          }
    }
}


void boundaryCellsSerial(struct CALModel3D* ca)
{
  for (int x_cell = 0; x_cell < ca->rows; x_cell++)
    for (int y_cell = 0; y_cell < ca->columns; y_cell++)
      for (int z_cell = 0; z_cell < ca->slices; z_cell++)
        {
//#define Z_HOLE
#ifdef Z_HOLE
          if (z_cell==0 && !(x_cell < 5 || x_cell > 14) && !(y_cell < 5 || y_cell > 14))
            continue;
#endif
          mb(ca,&Q,x_cell,y_cell,z_cell);
        }
}
