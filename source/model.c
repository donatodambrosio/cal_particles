#include <boundary.h>
#include <init.h>
#include <model.h>

struct CALModel3D* u_modellu;
struct CALRun3D* a_simulazioni;
struct Substates Q;


CALbyte ncestiArmenuNaParticella(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z, int n)
{
    for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
        if (calGetX3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
            return CAL_TRUE;

    return CAL_FALSE;
}

void movili(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
    if (!ncestiArmenuNaParticella(ca, cell_x, cell_y, cell_z, 0))
        return;

    CALreal delta_z = -CELL_SIDE/3;
    CALreal z;
    CALreal z_new;

    for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
        if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
        {
            z = calGet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z);
            z_new = z + delta_z;
            calSet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z,z_new);
        }
}

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

void sucala(int slot, int e, struct CALModel3D* ca, int cell_x, int cell_y, int cell_z, int n)
{
    calSet3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.px[e],cell_x,cell_y,cell_z,n));
    calSet3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.py[e],cell_x,cell_y,cell_z,n));
    calSet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.pz[e],cell_x,cell_y,cell_z,n));
    calSet3Dr(ca, Q.vx[slot],cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.vx[e],cell_x,cell_y,cell_z,n));
    calSet3Dr(ca, Q.vy[slot],cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.vy[e],cell_x,cell_y,cell_z,n));
    calSet3Dr(ca, Q.vz[slot],cell_x,cell_y,cell_z,calGetX3Dr(ca,Q.vz[e],cell_x,cell_y,cell_z,n));
    calSet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,calGetX3Di(ca,Q.imove[e],cell_x,cell_y,cell_z,n));
}

void moviliCazzu(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
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
        for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
            if (calGetX3Di(ca,Q.imove[slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
            {
                px = calGetX3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z,n);
                py = calGetX3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z,n);
                pz = calGetX3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z,n);

                new_cell_x = px/CELL_SIDE;
                new_cell_y = py/CELL_SIDE;
                new_cell_z = pz/CELL_SIDE;

                if ((cell_x == new_cell_x) && (cell_y == new_cell_y) && (cell_z == new_cell_z))
                {
                    int c;
                    for (c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
                        if (calGetNext3Di(ca,Q.imove[c],cell_x,cell_y,cell_z) == PARTICLE_ABSENT)
                        {
                            slot = c;
                            break;
                        }

                    if (c < MAX_NUMBER_OF_PARTICLES_PER_CELL)
                        sucala(c,slot,ca,cell_x,cell_y,cell_z,n);
                }
            }
}

CALbyte caminalu(struct CALModel3D* modello)
{
    CALint step = a_simulazioni->step;

    if (step <= STEPS)
        return CAL_FALSE;

    return CAL_TRUE;
}

void transizioniGlobali(struct CALModel3D* modello)
{
    calApplyElementaryProcess3D(modello,movili);
    calUpdate3D(modello);

    calApplyElementaryProcess3D(modello,moviliCazzu);
    calUpdate3D(modello);
}

void partilu()
{
    CALint NX = X_CELLS, NY=Y_CELLS, NZ=Z_CELLS;
    u_modellu = calCADef3D(X_CELLS,Y_CELLS,Z_CELLS,CAL_MOORE_NEIGHBORHOOD_3D,CAL_SPACE_TOROIDAL,CAL_NO_OPT);

    for(int slot=0;slot<MAX_NUMBER_OF_PARTICLES_PER_CELL;slot++)
    {
        Q.px[slot]    = calAddSubstate3Dr(u_modellu);
        Q.py[slot]    = calAddSubstate3Dr(u_modellu);
        Q.pz[slot]    = calAddSubstate3Dr(u_modellu);
        Q.vx[slot]    = calAddSubstate3Dr(u_modellu);
        Q.vy[slot]    = calAddSubstate3Dr(u_modellu);
        Q.vz[slot]    = calAddSubstate3Dr(u_modellu);
        Q.imove[slot] = calAddSubstate3Di(u_modellu);

        calInitSubstate3Dr(u_modellu,Q.px[slot],   PARTICLE_NODATA);
        calInitSubstate3Dr(u_modellu,Q.py[slot],   PARTICLE_NODATA);
        calInitSubstate3Dr(u_modellu,Q.pz[slot],   PARTICLE_NODATA);
        calInitSubstate3Dr(u_modellu,Q.vx[slot],   0);
        calInitSubstate3Dr(u_modellu,Q.vy[slot],   0);
        calInitSubstate3Dr(u_modellu,Q.vz[slot],   0);
        calInitSubstate3Di(u_modellu,Q.imove[slot],PARTICLE_ABSENT);
    }

    // Boundary
    calApplyElementaryProcess3D(u_modellu, boundary_cells);
    calUpdate3D(u_modellu);

    // Initial conditions
    calApplyElementaryProcess3D(u_modellu, mmiscali_nta_cella);
    calUpdate3D(u_modellu);

    // Simulation
    a_simulazioni = calRunDef3D(u_modellu,0,CAL_RUN_LOOP,CAL_UPDATE_IMPLICIT);
    calRunAddGlobalTransitionFunc3D(a_simulazioni, transizioniGlobali);
    calRunAddStopConditionFunc3D(a_simulazioni, caminalu);
}
