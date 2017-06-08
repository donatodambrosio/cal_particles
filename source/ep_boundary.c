#include <ep_boundary.h>

void mb(struct CALModel3D* ca, struct Substates *Q, int x_cell, int y_cell, int z_cell)
{
  for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      calSet3Dr(ca, Q->px[slot],x_cell,y_cell,z_cell,PARTICLE_NODATA);
      calSet3Dr(ca, Q->py[slot],x_cell,y_cell,z_cell,PARTICLE_NODATA);
      calSet3Dr(ca, Q->pz[slot],x_cell,y_cell,z_cell,PARTICLE_NODATA);
      calSet3Dr(ca, Q->vx[slot],x_cell,y_cell,z_cell,PARTICLE_NODATA);
      calSet3Dr(ca, Q->vy[slot],x_cell,y_cell,z_cell,PARTICLE_NODATA);
      calSet3Dr(ca, Q->vz[slot],x_cell,y_cell,z_cell,PARTICLE_NODATA);
      calSet3Di(ca, Q->imove[slot],x_cell,y_cell,z_cell,PARTICLE_BORDER);
    }
}

void boundary_cells(struct CALModel3D* ca, int x_cell, int y_cell, int z_cell)
{
  if (x_cell==0 || x_cell==X_CELLS-1 || y_cell==0 || y_cell==Y_CELLS-1 || z_cell==0 || z_cell==Z_CELLS-1)
    mb(ca,&Q,x_cell,y_cell,z_cell);
}
