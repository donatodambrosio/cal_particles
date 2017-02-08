/*
 * prova.c
 *
 *  Created on: 25 gen 2017
 *      Author: giuseppe
 */

#include "prova.h"
#include "parser.h"


struct CALModel3D* modello;
struct CALRun3D* simulazione;
struct Substates Q;

CALbyte estiDintra(struct CALModel3D* ca, int i, int j, int k, int n)
{
    for (int e = 0; e < MAX_NUMBER_OF_PARTICLES_PER_CELL; e++)
        if (calGet3Di(ca, Q.imove[e],i,j,k) != -3)
            return CAL_FALSE;

    return CAL_TRUE;
}

void movili(struct CALModel3D* ca, int i, int j, int k)
{
//    if (!estiDintra(ca, i, j, k, 0))
//        return;

    CALreal delta_z = CL;
    CALreal z;
    CALreal z_new;

    for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
        if (calGet3Di(ca, Q.imove[c],i,j,k) == 1)
        {
            z = calGet3Dr(ca, Q.pz[c],i,j,k);
            z_new = z + delta_z;
            calSet3Dr(ca, Q.pz[c],i,j,k,z_new);

            CALreal z_new_test = calGetNext3Dr(ca,Q.pz[c],i,j,k);
            int a = 0;
        }
}

void pezziala(int e, struct CALModel3D* ca, int i, int j, int k)
{
    calSet3Dr(ca, Q.px[e],i,j,k,NODATA);
    calSet3Dr(ca, Q.py[e],i,j,k,NODATA);
    calSet3Dr(ca, Q.pz[e],i,j,k,NODATA);
    calSet3Dr(ca, Q.vx[e],i,j,k,0);
    calSet3Dr(ca, Q.vy[e],i,j,k,0);
    calSet3Dr(ca, Q.vz[e],i,j,k,0);
    calSet3Di(ca, Q.imove[e],i,j,k,0);
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
        if (calGet3Di(ca, Q.imove[c],i,j,k) == 1)
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
                            if (calGetNext3Di(ca,Q.imove[c],i,j,k) == 0)
                            {
                                slot = c;
                                break;
                            }

                        if (slot != -1)
                            sucala(slot,e,ca,i,j,k,n);
                    }
                }
}



//void sucali(struct CALModel3D* ca, int i, int j, int k)
//{
//    CALreal z_min = k*CL;
//    CALreal z_max = (k+1)*CL;
//    CALreal z;

//    for (int n=0; n<ca->sizeof_X; n++)
//        if (estiDintra(ca,i,j,k,n))
//           for (int e = 0; e < MAX_NUMBER_OF_PARTICLES_PER_CELL; e++)
//               if (calGetX3Di(ca,Q.imove[e],i,j,k,n) == 1)
//                {
//                    z = calGetX3Dr(ca, Q.pz[e],i,j,k,n);
//                    CALint slot;
//                    if (z_min <= z && z < z_max)
//                    {
//                        for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
//                            if (calGetNext3Di(ca,Q.imove[c],i,j,k) == 0)
//                            {
//                                slot = c;
//                                break;
//                            }
//                        sucala(slot,e,ca,i,j,k,n);
//                    }
//                }
//}

void transition1(struct CALModel3D* ca ,int i,int j,int k)
{
    //printf("step = %d\n", simulazione->step);
	for(int w=0;w<MAX_NUMBER_OF_PARTICLES_PER_CELL;w++)
	{
		if((calGet3Dr(modello,Q.px[w],i,j,k)!=NODATA) && (calGet3Di(modello,Q.imove[w],i,j,k)==1)){
			CALreal posizione = calGet3Dr(modello,Q.pz[w],i,j,k);
			CALreal velocita = calGet3Dr(modello,Q.vz[w],i,j,k);
			CALreal deltaT = (simulazione->step - (simulazione->step-1))*0.005;
			CALreal nuovaPosizione =((0.5*-1*(deltaT*deltaT))+(velocita*deltaT)+posizione);
			CALreal nuovaVelocita = velocita+(-1*(deltaT));
			CALint nuovaK = nuovaPosizione/CL;

            nuovaK = k-1;

//            printf("%i,%i\n",k,nuovaK);

			if(nuovaK!=k){
				for(int q=0;q<MAX_NUMBER_OF_PARTICLES_PER_CELL;q++){
					if(calGet3Dr(modello,Q.px[q],i,j,nuovaK)==NODATA){
					calSet3Dr(modello,Q.pz[q],i,j,nuovaK,nuovaPosizione);
					calSet3Dr(modello,Q.px[q],i,j,nuovaK,calGet3Dr(modello,Q.px[w],i,j,k));
					calSet3Dr(modello,Q.py[q],i,j,nuovaK,calGet3Dr(modello,Q.py[w],i,j,k));
					calSet3Dr(modello,Q.vz[q],i,j,nuovaK,nuovaVelocita);
					calSet3Dr(modello,Q.vx[q],i,j,nuovaK,calGet3Dr(modello,Q.vx[w],i,j,k));
					calSet3Dr(modello,Q.vy[q],i,j,nuovaK,calGet3Dr(modello,Q.vy[w],i,j,k));
					calSet3Di(modello,Q.imove[q],i,j,nuovaK,1);
                    calSet3Dr(modello,Q.px[w],i,j,k,NODATA);
                    calSet3Dr(modello,Q.py[w],i,j,k,NODATA);
                    calSet3Dr(modello,Q.pz[w],i,j,k,NODATA);
                    calSet3Dr(modello,Q.vx[w],i,j,k,0);
                    calSet3Dr(modello,Q.vy[w],i,j,k,0);
                    calSet3Dr(modello,Q.vz[w],i,j,k,0);
                    calSet3Di(modello,Q.imove[w],i,j,k,0);
					break;
					}
				}
			}
		}
	}
}

void transition2(struct CALModel3D* ca ,int i,int j,int k)
{
	//	for(int p=1;p<27;p++)
	//	{
	//		for(int w=0;w<MAX_NUMBER_OF_PARTICLES_PER_CELL;w++)
	//		{
	//			if(k<=49){
	//				CALint checkX = calGetX3Dr(modello,Q.px[w],i,j,k,p)/CL;
	//				CALint checkY = calGetX3Dr(modello,Q.py[w],i,j,k,p)/CL;
	//				CALint checkZ = calGetX3Dr(modello,Q.pz[w],i,j,k,p)/CL;
	//				if(checkX == i && checkY == j && checkZ == k)
	//				{
	//					for(int z=0;z<MAX_NUMBER_OF_PARTICLES_PER_CELL;z++)
	//					{
	//						if(calGet3Dr(modello,Q.px[z],i,j,k)==NODATA)
	//						{
	//							calSet3Dr(modello,Q.px[z],i,j,k,calGetX3Dr(modello,Q.px[w],i,j,k,p));
	//							calSet3Dr(modello,Q.py[z],i,j,k,calGetX3Dr(modello,Q.py[w],i,j,k,p));
	//							calSet3Dr(modello,Q.pz[z],i,j,k,calGetX3Dr(modello,Q.pz[w],i,j,k,p));
	//							calSet3Di(modello,Q.imove[z],i,j,k,calGetX3Di(modello,Q.imove[w],i,j,k,p));
	//							break;
	//						}
	//					}
	//					calSetX3Dr(modello,Q.px[w],i,j,k,p,NODATA);
	//					calSetX3Dr(modello,Q.py[w],i,j,k,p,NODATA);
	//					calSetX3Dr(modello,Q.pz[w],i,j,k,p,NODATA);
	//					calSetX3Di(modello,Q.imove[w],i,j,k,p,0);
	//				}
	//			}
	//		}
	//	}

}
void setPosition(const CALreal x, const CALreal y, const CALreal z,const CALint imove, const int p)
{
	int i = x/CL;
	int j = y/CL;
	int k = z/CL;

	//	if(p<=14000){
	//		//printf("%f,%f,%f,%i,%i,%i\n",x,y,z,i,j,k);
	//	}

	int w=0;

	for (w=0; w<MAX_NUMBER_OF_PARTICLES_PER_CELL; w++)
		if( calGet3Dr(modello,Q.px[w],i,j,k) == NODATA ){
			calInit3Dr(modello,Q.px[w],i,j,k,x);
			calInit3Dr(modello,Q.py[w],i,j,k,y);
			calInit3Dr(modello,Q.pz[w],i,j,k,z);

			calInit3Di(modello,Q.imove[w],i,j,k,imove);
			// CALint numP = calGet3Di(modello,cella.numPresenti,i,j,k);
			// calInit3Di(modello,cella.numPresenti,i,j,k,numP+1);
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


void partilu()
{
    modello = calCADef3D(ROWS,COLS,LAYERS,CAL_MOORE_NEIGHBORHOOD_3D,CAL_SPACE_TOROIDAL,CAL_NO_OPT);

	for(int i=0;i<MAX_NUMBER_OF_PARTICLES_PER_CELL;i++){
		Q.px[i] = calAddSubstate3Dr(modello);
		Q.py[i] = calAddSubstate3Dr(modello);
		Q.pz[i] = calAddSubstate3Dr(modello);

		calInitSubstate3Dr(modello,Q.px[i],NODATA);
		calInitSubstate3Dr(modello,Q.py[i],NODATA);
		calInitSubstate3Dr(modello,Q.pz[i],NODATA);

		Q.vx[i] = calAddSubstate3Dr(modello);
		Q.vy[i] = calAddSubstate3Dr(modello);
		Q.vz[i] = calAddSubstate3Dr(modello);

		calInitSubstate3Dr(modello,Q.vx[i],0);
		calInitSubstate3Dr(modello,Q.vy[i],0);
		calInitSubstate3Dr(modello,Q.vz[i],0);

		Q.imove[i]=calAddSubstate3Di(modello);
        calInitSubstate3Di(modello,Q.imove[i],0);
	}
	//    cella.numPresenti = calAddSubstate3Di(modello);
	//    calInitSubstate3Di(modello,cella.numPresenti,0);


	particellaT particelle[NUMBER_OF_PARTICLES];
	parse(particelle);
	for(int i=0;i<NUMBER_OF_PARTICLES;i++){
		setPosition(particelle[i].posizione[0],particelle[i].posizione[1],particelle[i].posizione[2],particelle[i].imove,i);
	}
	//
	//	setPosition(0.005,0.005,0.005,1,0);
	//	setPosition(0.003,0.003,0.003,1,1);

//	calAddElementaryProcess3D(modello,transition1);
//	calAddElementaryProcess3D(modello,transition2);


    calApplyElementaryProcess3D(modello, printPos);


    calAddElementaryProcess3D(modello,movili);
    calAddElementaryProcess3D(modello,moviliCazzu);
//    calAddElementaryProcess3D(modello,sucali);




    simulazione = calRunDef3D(modello,1,CAL_RUN_LOOP,CAL_UPDATE_IMPLICIT);
    calRunAddStopConditionFunc3D(simulazione, caminalu);

}

