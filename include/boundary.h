#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <model.h>

void boundary_cells(struct CALModel3D* ca, int i, int j, int k)
{
    if (i==0 || i==ROWS-1 || j==0 || j==COLS-1 || k==0 || k==SLICES-1)
        for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
        {
            calSet3Dr(ca, Q.px[slot],i,j,k,NODATA);
            calSet3Dr(ca, Q.py[slot],i,j,k,NODATA);
            calSet3Dr(ca, Q.pz[slot],i,j,k,NODATA);
            calSet3Dr(ca, Q.vx[slot],i,j,k,NODATA);
            calSet3Dr(ca, Q.vy[slot],i,j,k,NODATA);
            calSet3Dr(ca, Q.vz[slot],i,j,k,NODATA);
            calSet3Di(ca, Q.imove[slot],i,j,k,PARTICLE_EDGE);
        }
}



#endif
