#include <boundary.h>
#include <init.h>
#include <math.h>
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

void applyForce(CALreal* p0, CALreal* v0, CALreal m, CALreal t, CALreal* pf, CALreal* vf)
{
    CALreal F[3], a[3];

    // F[0] =  0;
    // F[1] =  0;
    // F[2] = -m*G;

    F[0] =  m*0 - 6*M_PI*AIR_VISCOSITY*PARTICLE_RADIUS*v0[0];
    F[1] =  m*0 - 6*M_PI*AIR_VISCOSITY*PARTICLE_RADIUS*v0[1];
    F[2] = -m*G - 6*M_PI*AIR_VISCOSITY*PARTICLE_RADIUS*v0[2];


    a[0] = F[0]/m;
    a[1] = F[1]/m;
    a[2] = F[2]/m;

    for (int i=0; i<3; i++)
    {
        vf[i] = v0[i]+a[i]*t;
        pf[i] = p0[i] + v0[i]*t + 0.5*a[i]*t*t;

        if (fabs(pf[i]-p0[i]) >= CELL_SIDE)
        {
            pf[i] = p0[i];
            vf[i] = v0[i];
        }
    }
}

void movili(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
    if (!ncestiArmenuNaParticella(ca, cell_x, cell_y, cell_z, 0))
        return;

    CALreal p0[3];
    CALreal v0[3];
    CALreal pf[3];
    CALreal vf[3];

    for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
        if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
        {
            p0[0] = pf[0] = calGet3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z);
            p0[1] = pf[1] = calGet3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z);
            p0[2] = pf[2] = calGet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z);

            v0[0] = vf[0] = calGet3Dr(ca, Q.vx[slot],cell_x,cell_y,cell_z);
            v0[1] = vf[1] = calGet3Dr(ca, Q.vy[slot],cell_x,cell_y,cell_z);
            v0[2] = vf[2] = calGet3Dr(ca, Q.vz[slot],cell_x,cell_y,cell_z);

            applyForce(p0, v0, PARTICLE_MASS, DELTA_T, pf, vf);

            calSet3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z,pf[0]);
            calSet3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z,pf[1]);
            calSet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z,pf[2]);

            calSet3Dr(ca, Q.vx[slot],cell_x,cell_y,cell_z,vf[0]);
            calSet3Dr(ca, Q.vy[slot],cell_x,cell_y,cell_z,vf[1]);
            calSet3Dr(ca, Q.vz[slot],cell_x,cell_y,cell_z,vf[2]);
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

    Q.px = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
    Q.py = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
    Q.pz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
    Q.vx = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
    Q.vy = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
    Q.vz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
    Q.imove = (struct CALSubstate3Di**)malloc(sizeof(struct CALSubstate3Di*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);

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
