#ifndef INIT_H
#define INIT_H

#include <model.h>

void mmiscali_nta_cella(struct CALModel3D* ca, int x_cell, int y_cell, int z_cell)
{
    if (z_cell < TOP_LAYERS || calGet3Di(ca, Q.imove[0],x_cell,y_cell,z_cell) == PARTICLE_BORDER)
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

            CALreal px = CELL_SIDE * (x_cell + cx);
            CALreal py = CELL_SIDE * (y_cell + cy);
            CALreal pz = CELL_SIDE * (z_cell + cz);

            calSet3Dr(ca, Q.px[slot],x_cell,y_cell,z_cell,px);
            calSet3Dr(ca, Q.py[slot],x_cell,y_cell,z_cell,py);
            calSet3Dr(ca, Q.pz[slot],x_cell,y_cell,z_cell,pz);
            calSet3Dr(ca, Q.vx[slot],x_cell,y_cell,z_cell,0.0);
            calSet3Dr(ca, Q.vy[slot],x_cell,y_cell,z_cell,0.0);
            calSet3Dr(ca, Q.vz[slot],x_cell,y_cell,z_cell,0.0);
            calSet3Di(ca, Q.imove[slot],x_cell,y_cell,z_cell,PARTICLE_PRESENT);
        }
    }
}


#endif
