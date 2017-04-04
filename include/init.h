#ifndef INIT_H
#define INIT_H

#include <model.h>

void mmiscali_nta_cella(struct CALModel3D* ca, int i, int j, int k)
{
    if ((SLICES-1-k) > TOP_LAYERS || calGet3Di(ca, Q.imove[0],i,j,k) == PARTICLE_BORDER)
        return;

    CALreal c;
    for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
        c = rand()/RAND_MAX; // 0 <= c <= 1
        if (c < CELL_FILL_RATE)
        {

            CALreal cx = (CALreal)rand()/(CALreal)(RAND_MAX); // 0 <= c <= 1
            CALreal cy = (CALreal)rand()/(CALreal)(RAND_MAX); // 0 <= c <= 1
            CALreal cz = (CALreal)rand()/(CALreal)(RAND_MAX); // 0 <= c <= 1

            CALreal x = CELL_SIDE * (j + cx);
            CALreal y = CELL_SIDE * (ROWS-1-i + cy);
            CALreal z = CELL_SIDE * (SLICES-1-k + cz);

            calSet3Dr(ca, Q.px[slot],i,j,k,x);
            calSet3Dr(ca, Q.py[slot],i,j,k,y);
            calSet3Dr(ca, Q.pz[slot],i,j,k,z);
            calSet3Dr(ca, Q.vx[slot],i,j,k,0.0);
            calSet3Dr(ca, Q.vy[slot],i,j,k,0.0);
            calSet3Dr(ca, Q.vz[slot],i,j,k,0.0);
            calSet3Di(ca, Q.imove[slot],i,j,k,PARTICLE_PRESENT);
        }
    }
}


#endif
