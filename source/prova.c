#include "prova.h"
#include "parser.h"

struct CALModel3D* modello;
struct CALRun3D* simulazione;
struct Substates Q;

CALbyte estiDintra(struct CALModel3D* ca, int i, int j, int k, int n)
{
    for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
        if (calGetX3Di(ca, Q.imove[slot],i,j,k,n) == PARTICLE_PRESENT)
            return CAL_TRUE;

    return CAL_FALSE;
}

void movili(struct CALModel3D* ca, int i, int j, int k)
{
//    if (!estiDintra(ca, i, j, k, 0))
//        return;

    CALreal delta_z = CL;
    CALreal z;
    CALreal z_new;

    for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
        if (calGet3Di(ca, Q.imove[c],i,j,k) == PARTICLE_PRESENT)
        {
            z = calGet3Dr(ca, Q.pz[c],i,j,k);
            z_new = z + delta_z;
            calSet3Dr(ca, Q.pz[c],i,j,k,z_new);

            CALreal z_new_test = calGetNext3Dr(ca,Q.pz[c],i,j,k);
            int a = 0;
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
//    if (!estiDintra(ca, i, j, k, 0))
//        return;

    CALreal z_min = k*CL;
    CALreal z_max = (k+1)*CL;
    CALreal z;

    //pezziali
    for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
        if (calGet3Di(ca, Q.imove[c],i,j,k) == PARTICLE_PRESENT)
        {
            z = calGet3Dr(ca, Q.pz[c],i,j,k);
            if (!(z_min <= z && z < z_max))
                pezziala(c, ca,i,j,k);
        }

    //sucali
    for (int n=1; n<ca->sizeof_X; n++)
        //if (estiDintra(ca,i,j,k,n))
           for (int e = 0; e < MAX_NUMBER_OF_PARTICLES_PER_CELL; e++)
               if (calGetX3Dr(ca,Q.px[e],i,j,k,n) == NODATA)
                {
                    z = calGetX3Dr(ca, Q.pz[e],i,j,k,n);

                    if (z_min <= z && z < z_max)
                    {
                        CALint slot = -1;
                        for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
                            if (calGetNext3Di(ca,Q.imove[c],i,j,k) == PARTICLE_ABSENT)
                            {
                                slot = c;
                                break;
                            }

                        if (slot != -1)
                            sucala(slot,e,ca,i,j,k,n);
                    }
                }
}

void setPosition(struct CALModel3D* ca, const CALreal x, const CALreal y, const CALreal z,const CALint imove)
{
	int i = x/CL;
	int j = y/CL;
	int k = z/CL;

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
    CALint step = simulazione->step;

    if (step <= step)
        return CAL_TRUE;

    return CAL_FALSE;
}

void printPos(struct CALModel3D* ca, int i, int j, int k)
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
    if (estiDintra(ca, i, j, k, 0))
        for (int slot = 0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
            pezziala(slot, ca, i, j, k);
}

void partilu()
{
    modello = calCADef3D(ROWS,COLS,LAYERS,CAL_MOORE_NEIGHBORHOOD_3D,CAL_SPACE_TOROIDAL,CAL_NO_OPT);

    for(int i=0;i<MAX_NUMBER_OF_PARTICLES_PER_CELL;i++)
    {
        Q.px[i]    = calAddSubstate3Dr(modello);
        Q.py[i]    = calAddSubstate3Dr(modello);
        Q.pz[i]    = calAddSubstate3Dr(modello);
        Q.vx[i]    = calAddSubstate3Dr(modello);
        Q.vy[i]    = calAddSubstate3Dr(modello);
        Q.vz[i]    = calAddSubstate3Dr(modello);
        Q.imove[i] = calAddSubstate3Di(modello);

        calInitSubstate3Dr(modello,Q.px[i],   NODATA);
        calInitSubstate3Dr(modello,Q.py[i],   NODATA);
        calInitSubstate3Dr(modello,Q.pz[i],   NODATA);
        calInitSubstate3Dr(modello,Q.vx[i],   0);
        calInitSubstate3Dr(modello,Q.vy[i],   0);
        calInitSubstate3Dr(modello,Q.vz[i],   0);
        calInitSubstate3Di(modello,Q.imove[i],PARTICLE_ABSENT);
	}

	particellaT particelle[NUMBER_OF_PARTICLES];
	parse(particelle);
    for(int i=0;i<NUMBER_OF_PARTICLES;i++)
        setPosition(modello, particelle[i].posizione[0],particelle[i].posizione[1],particelle[i].posizione[2],particelle[i].imove);


    // calApplyElementaryProcess3D(modello, printPos);

    calApplyElementaryProcess3D(modello, sbacanta);
    calUpdate3D(modello);

    setPosition(modello, 0.025, 0.025, 0.025, PARTICLE_PRESENT);
    setPosition(modello, 0.025, 0.025, 0.026, PARTICLE_PRESENT);
    calUpdate3D(modello);

//    calAddElementaryProcess3D(modello,movili);
//    calAddElementaryProcess3D(modello,moviliCazzu);

    simulazione = calRunDef3D(modello,1,CAL_RUN_LOOP,CAL_UPDATE_IMPLICIT);
    calRunAddStopConditionFunc3D(simulazione, caminalu);
}
