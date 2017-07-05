#include <ep_boundary.h>

void mb(struct CALModel3D* ca, struct Substates *Q, int cell_x, int cell_y, int cell_z)
{
  CALreal x, y, z, nx, ny, nz;

  x  = 0.0;
  y  = 0.0;
  z  = 0.0;
  nx = 0.0;
  ny = 0.0;
  nz = 0.0;

  // x_cell = 0 (i.e. x = CELL_SIDE)
  if (cell_x == 0)
    {
      x = CELL_SIDE;
      y = (cell_y * CELL_SIDE) + CELL_SIDE/2;
      z = (cell_z * CELL_SIDE) + CELL_SIDE/2;
      nx = 1;
      ny = 0;
      nz = 0;
    }

  // x_cell = X_CELLS-1 (i.e. x = (X_CELLS-1)*CELL_SIDE)
  if (cell_x == X_CELLS-1)
    {
      x = (X_CELLS-1)*CELL_SIDE;
      y = (cell_y * CELL_SIDE) + CELL_SIDE/2;
      z = (cell_z * CELL_SIDE) + CELL_SIDE/2;
      nx = -1;
      ny = 0;
      nz = 0;
    }

  // y_cell = 0 (i.e. y = CELL_SIDE)
  if (cell_y == 0)
    {
      x = (cell_x * CELL_SIDE) + CELL_SIDE/2;
      y = CELL_SIDE;
      z = (cell_z * CELL_SIDE) + CELL_SIDE/2;
      nx = 0;
      ny = 1;
      nz = 0;
    }

  // y_cell = Y_CELLS-1 (i.e. x = (Y_CELLS-1)*CELL_SIDE)
  if (cell_y == Y_CELLS-1)
    {
      x = (cell_x * CELL_SIDE) + CELL_SIDE/2;
      y = (Y_CELLS-1)*CELL_SIDE;
      z = (cell_z * CELL_SIDE) + CELL_SIDE/2;
      nx = 0;
      ny = -1;
      nz = 0;
    }

  // z_cell = 0 (i.e. z = CELL_SIDE)
  if (cell_z == 0)
    {
      x = (cell_x * CELL_SIDE) + CELL_SIDE/2;
      y = (cell_y * CELL_SIDE) + CELL_SIDE/2;
      z = CELL_SIDE;
      nx = 0;
      ny = 0;
      nz = 1;
    }

  // z_cell = Z_CELLS-1 (i.e. z = (Z_CELLS-1)*CELL_SIDE)
  if (cell_z == Z_CELLS-1)
    {
      x = (cell_x * CELL_SIDE) + CELL_SIDE/2;
      y = (cell_y * CELL_SIDE) + CELL_SIDE/2;
      z = (Z_CELLS-1)*CELL_SIDE;
      nx = 0;
      ny = 0;
      nz = -1;
    }

  for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      calInit3Dr(ca, Q->Fx[slot],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q->Fy[slot],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q->Fz[slot],cell_x,cell_y,cell_z,0.0);
      calInit3Dr(ca, Q->px[slot],cell_x,cell_y,cell_z,x);
      calInit3Dr(ca, Q->py[slot],cell_x,cell_y,cell_z,y);
      calInit3Dr(ca, Q->pz[slot],cell_x,cell_y,cell_z,z);
      calInit3Dr(ca, Q->vx[slot],cell_x,cell_y,cell_z,nx);
      calInit3Dr(ca, Q->vy[slot],cell_x,cell_y,cell_z,ny);
      calInit3Dr(ca, Q->vz[slot],cell_x,cell_y,cell_z,nz);
      calInit3Di(ca, Q->ID[slot],cell_x,cell_y,cell_z,BORDER_ID);
    }
}

void boundary_cells(struct CALModel3D* ca, int x_cell, int y_cell, int z_cell)
{
  if (x_cell==0 || x_cell==X_CELLS-1 || y_cell==0 || y_cell==Y_CELLS-1 || (z_cell==0 && (x_cell < 5 || x_cell > 14)) || z_cell==Z_CELLS-1)
    mb(ca,&Q,x_cell,y_cell,z_cell);
}
