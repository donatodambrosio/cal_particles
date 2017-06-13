#include <ep_boundary.h>

void mb(struct CALModel3D* ca, struct Substates *Q, int x_cell, int y_cell, int z_cell)
{
  CALreal x, y, z, nx, ny, nz;

  x  = PARTICLE_NODATA;
  y  = PARTICLE_NODATA;
  z  = PARTICLE_NODATA;
  nx = PARTICLE_NODATA;
  ny = PARTICLE_NODATA;
  nz = PARTICLE_NODATA;

  // x_cell = 0 (i.e. x = CELL_SIDE)
  if (x_cell == 0)
    {
      x = CELL_SIDE;
      y = (y_cell * CELL_SIDE) + CELL_SIDE/2;
      z = (z_cell * CELL_SIDE) + CELL_SIDE/2;
      nx = 1;
      ny = 0;
      nz = 0;
    }

  // x_cell = X_CELLS-1 (i.e. x = (X_CELLS-1)*CELL_SIDE)
  if (x_cell == X_CELLS-1)
    {
      x = (X_CELLS-1)*CELL_SIDE;
      y = (y_cell * CELL_SIDE) + CELL_SIDE/2;
      z = (z_cell * CELL_SIDE) + CELL_SIDE/2;
      nx = -1;
      ny = 0;
      nz = 0;
    }

  // y_cell = 0 (i.e. y = CELL_SIDE)
  if (y_cell == 0)
    {
      x = (x_cell * CELL_SIDE) + CELL_SIDE/2;
      y = CELL_SIDE;
      z = (z_cell * CELL_SIDE) + CELL_SIDE/2;
      nx = 0;
      ny = 1;
      nz = 0;
    }

  // y_cell = Y_CELLS-1 (i.e. x = (Y_CELLS-1)*CELL_SIDE)
  if (y_cell == Y_CELLS-1)
    {
      x = (x_cell * CELL_SIDE) + CELL_SIDE/2;
      y = (Y_CELLS-1)*CELL_SIDE;
      z = (z_cell * CELL_SIDE) + CELL_SIDE/2;
      nx = 0;
      ny = -1;
      nz = 0;
    }

  // z_cell = 0 (i.e. z = CELL_SIDE)
  if (z_cell == 0)
    {
      x = (x_cell * CELL_SIDE) + CELL_SIDE/2;
      y = (y_cell * CELL_SIDE) + CELL_SIDE/2;
      z = CELL_SIDE;
      nx = 0;
      ny = 0;
      nz = 1;
    }

  // z_cell = Z_CELLS-1 (i.e. z = (Z_CELLS-1)*CELL_SIDE)
  if (z_cell == Z_CELLS-1)
    {
      x = (x_cell * CELL_SIDE) + CELL_SIDE/2;
      y = (y_cell * CELL_SIDE) + CELL_SIDE/2;
      z = (Z_CELLS-1)*CELL_SIDE;
      nx = 0;
      ny = 0;
      nz = -1;
    }


  for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      calInit3Dr(ca, Q->px[slot],x_cell,y_cell,z_cell,x);
      calInit3Dr(ca, Q->py[slot],x_cell,y_cell,z_cell,y);
      calInit3Dr(ca, Q->pz[slot],x_cell,y_cell,z_cell,z);
      calInit3Dr(ca, Q->vx[slot],x_cell,y_cell,z_cell,nx);
      calInit3Dr(ca, Q->vy[slot],x_cell,y_cell,z_cell,ny);
      calInit3Dr(ca, Q->vz[slot],x_cell,y_cell,z_cell,nz);
      calInit3Di(ca, Q->imove[slot],x_cell,y_cell,z_cell,PARTICLE_BORDER);
    }
}

void boundary_cells(struct CALModel3D* ca, int x_cell, int y_cell, int z_cell)
{
  if (x_cell==0 || x_cell==X_CELLS-1 || y_cell==0 || y_cell==Y_CELLS-1 || z_cell==0 || z_cell==Z_CELLS-1)
    mb(ca,&Q,x_cell,y_cell,z_cell);
}
