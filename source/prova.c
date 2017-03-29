#include "prova.h"
#include "parser.h"

struct CALModel3D* u_modellu;
struct CALRun3D* a_simulazioni;
struct Substates Q;

CALbyte ncestiArmenuNaParticella(struct CALModel3D* ca, int i, int j, int k, int n)
{
    for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
        if (calGetX3Di(ca, Q.imove[slot],i,j,k,n) == PARTICLE_PRESENT)
            return CAL_TRUE;

    return CAL_FALSE;
}

void movili(struct CALModel3D* ca, int i, int j, int k)
{
    if (!ncestiArmenuNaParticella(ca, i, j, k, 0))
        return;

    CALreal delta_z = -CL/3;
    CALreal z;
    CALreal z_new;

    for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
        if (calGet3Di(ca, Q.imove[c],i,j,k) == PARTICLE_PRESENT)
        {
            z = calGet3Dr(ca, Q.pz[c],i,j,k);
            z_new = z + delta_z;
            calSet3Dr(ca, Q.pz[c],i,j,k,z_new);

            //            CALreal z_new_test = calGetNext3Dr(ca,Q.pz[c],i,j,k);
            //            int a = 0;
        }
}

void pezziala(int slot, struct CALModel3D* ca, int i, int j, int k)
{
    calSet3Dr(ca, Q.px[slot],   i,j,k,NODATA);
    calSet3Dr(ca, Q.py[slot],   i,j,k,NODATA);
    calSet3Dr(ca, Q.pz[slot],   i,j,k,NODATA);
    calSet3Dr(ca, Q.vx[slot],   i,j,k,0);
    calSet3Dr(ca, Q.vy[slot],   i,j,k,0);
    calSet3Dr(ca, Q.vz[slot],   i,j,k,0);
    calSet3Di(ca, Q.imove[slot],i,j,k,PARTICLE_ABSENT);
}

void sucala(int slot, int e, struct CALModel3D* ca, int i, int j, int k, int n)
{
    calSet3Dr(ca, Q.px[slot],i,j,k,calGetX3Dr(ca,Q.px[e],i,j,k,n));
    calSet3Dr(ca, Q.py[slot],i,j,k,calGetX3Dr(ca,Q.py[e],i,j,k,n));
    calSet3Dr(ca, Q.pz[slot],i,j,k,calGetX3Dr(ca,Q.pz[e],i,j,k,n));
    calSet3Dr(ca, Q.vx[slot],i,j,k,calGetX3Dr(ca,Q.vx[e],i,j,k,n));
    calSet3Dr(ca, Q.vy[slot],i,j,k,calGetX3Dr(ca,Q.vy[e],i,j,k,n));
    calSet3Dr(ca, Q.vz[slot],i,j,k,calGetX3Dr(ca,Q.vz[e],i,j,k,n));
    calSet3Di(ca, Q.imove[slot],i,j,k,calGetX3Di(ca,Q.imove[e],i,j,k,n));
}

void moviliCazzu(struct CALModel3D* ca, int i, int j, int k)
{
    //    CALreal z_min = k*CL;
    //    CALreal z_max = (k+1)*CL;
    //    CALreal z_min = (SLICES-1 - (k+1))*CL;
    //    CALreal z_max = (SLICES-1 - k)*CL;

    CALreal x,  y,  z;
    CALint _i, _j, _k;

    if (ncestiArmenuNaParticella(ca, i, j, k, 0))
    {
        //pezziali
        for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
            if (calGet3Di(ca, Q.imove[slot],i,j,k) == PARTICLE_PRESENT)
            {
                x = calGet3Dr(ca, Q.px[slot],i,j,k);
                y = calGet3Dr(ca, Q.py[slot],i,j,k);
                z = calGet3Dr(ca, Q.pz[slot],i,j,k);

                _i = x/CL;
                _j = y/CL;
                _k = z/CL;

                if ((i != _i) || (j != _j) || (k != _k))
                    pezziala(slot, ca,i,j,k);
            }
    }

    //sucali
    for (int n=1; n<ca->sizeof_X; n++)
        for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
            if (calGetX3Di(ca,Q.imove[slot],i,j,k,n) == PARTICLE_PRESENT)
            {
                x = calGetX3Dr(ca, Q.px[slot],i,j,k,n);
                y = calGetX3Dr(ca, Q.py[slot],i,j,k,n);
                z = calGetX3Dr(ca, Q.pz[slot],i,j,k,n);

                _i = x/CL;
                _j = y/CL;
                _k = z/CL;

                int a;
                if (i==25 && j==25 && k==24 && n==18 && slot==0)
                    a = 0;

                if ((i == _i) && (j == _j) && (k == _k))
                {
                    int c;
                    for (c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
                        if (calGetNext3Di(ca,Q.imove[c],i,j,k) == PARTICLE_ABSENT)
                        {
                            slot = c;
                            break;
                        }

                    if (c < MAX_NUMBER_OF_PARTICLES_PER_CELL)
                        sucala(c,slot,ca,i,j,k,n);
                }
            }
}

void danciNaPosizioni(struct CALModel3D* ca, const CALreal x, const CALreal y, const CALreal z,const CALint imove)
{
    int i = x/CL;
    int j = y/CL;
    int k = z/CL;
    //int k = (SLICES -1) -  z/CL;

    CALint slot = -1;
    for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
        if (calGet3Di(ca,Q.imove[c],i,j,k) == PARTICLE_ABSENT)
        {
            slot = c;
            break;
        }

    if(slot != -1)
    {
        calInit3Dr(ca,Q.px[slot],   i,j,k,x    );
        calInit3Dr(ca,Q.py[slot],   i,j,k,y    );
        calInit3Dr(ca,Q.pz[slot],   i,j,k,z    );
        calInit3Di(ca,Q.imove[slot],i,j,k,imove);
    }
}

CALbyte caminalu(struct CALModel3D* modello)
{
    CALint step = a_simulazioni->step;

    if (step <= STEPS)
        return CAL_FALSE;

    return CAL_TRUE;
}

void stampaAPosizioni(struct CALModel3D* ca, int i, int j, int k)
{
    CALreal z = 0;
    if (k == 10)
        for (int c=0; c<MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
        {
            z = calGet3Dr(ca,Q.pz[0],i,j,k);
            if (z != NODATA)
                printf ("[(i,j,k),[c],z]=[(%d,%d,%d),[%d],%f]\n", i,j,k,c,z);
        }
}

void sbacanta(struct CALModel3D* ca, int i, int j, int k)
{
    if (ncestiArmenuNaParticella(ca, i, j, k, 0))
        for (int slot = 0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
            pezziala(slot, ca, i, j, k);
}

void transizioniGlobali(struct CALModel3D* modello)
{
    calApplyElementaryProcess3D(modello,movili);
    calUpdate3D(modello);

    //    printf("z(25,25,25) = %f\n", calGet3Dr(modello,Q.pz[0],25,25,25));
    //    for (int n=0; n<modello->sizeof_X; n++)
    //        printf("z(25,25,24; %d) = %f\n", n, calGetX3Dr(modello,Q.pz[0],25,25,24,n));

    calApplyElementaryProcess3D(modello,moviliCazzu);
    calUpdate3D(modello);
}


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

void mmiscali_nta_cella(struct CALModel3D* ca, int i, int j, int k)
{
    if (k > TOP_LAYERS && calGet3Di(ca, Q.imove[0],i,j,k) != PARTICLE_EDGE)
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

            CALreal x = CL * (j + cx);
            CALreal y = CL * (ROWS-1-i + cy);
            CALreal z = CL * (SLICES-1-k + cz);

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

void partilu()
{
    CALint NX = COLS, NY=ROWS, NZ=SLICES;
    u_modellu = calCADef3D(ROWS,COLS,SLICES,CAL_MOORE_NEIGHBORHOOD_3D,CAL_SPACE_TOROIDAL,CAL_NO_OPT);

    for(int i=0;i<MAX_NUMBER_OF_PARTICLES_PER_CELL;i++)
    {
        Q.px[i]    = calAddSubstate3Dr(u_modellu);
        Q.py[i]    = calAddSubstate3Dr(u_modellu);
        Q.pz[i]    = calAddSubstate3Dr(u_modellu);
        Q.vx[i]    = calAddSubstate3Dr(u_modellu);
        Q.vy[i]    = calAddSubstate3Dr(u_modellu);
        Q.vz[i]    = calAddSubstate3Dr(u_modellu);
        Q.imove[i] = calAddSubstate3Di(u_modellu);

        calInitSubstate3Dr(u_modellu,Q.px[i],   NODATA);
        calInitSubstate3Dr(u_modellu,Q.py[i],   NODATA);
        calInitSubstate3Dr(u_modellu,Q.pz[i],   NODATA);
        calInitSubstate3Dr(u_modellu,Q.vx[i],   0);
        calInitSubstate3Dr(u_modellu,Q.vy[i],   0);
        calInitSubstate3Dr(u_modellu,Q.vz[i],   0);
        calInitSubstate3Di(u_modellu,Q.imove[i],PARTICLE_ABSENT);
    }

    // Boundary
    calApplyElementaryProcess3D(u_modellu, boundary_cells);
    calUpdate3D(u_modellu);

    // Initial conditions
    calApplyElementaryProcess3D(u_modellu, mmiscali_nta_cella);
    calUpdate3D(u_modellu);

    //    particellaT particelle[NUMBER_OF_PARTICLES];
    //    parse(particelle);
    //    for(int i=0;i<NUMBER_OF_PARTICLES;i++)
    //        danciNaPosizioni(u_modellu, particelle[i].posizione[0],particelle[i].posizione[1],particelle[i].posizione[2],particelle[i].imove);

    // calApplyElementaryProcess3D(modello, printPos);

    // calApplyElementaryProcess3D(modello, sbacanta);
    // calUpdate3D(modello);

    // setPosition(modello, 0.000, 0.000, 0.000, PARTICLE_PRESENT);
    // setPosition(modello, 0.000, 0.000, 0.050, PARTICLE_PRESENT);
    // setPosition(modello, 0.025, 0.025, 0.025, PARTICLE_PRESENT);
    // setPosition(modello, 0.025, 0.025, 0.026, PARTICLE_PRESENT);
    calUpdate3D(u_modellu);

    a_simulazioni = calRunDef3D(u_modellu,1,CAL_RUN_LOOP,CAL_UPDATE_IMPLICIT);
    calRunAddGlobalTransitionFunc3D(a_simulazioni, transizioniGlobali);
    calRunAddStopConditionFunc3D(a_simulazioni, caminalu);
}
