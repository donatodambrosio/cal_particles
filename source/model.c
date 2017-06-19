#include <ep_boundary.h>
#include <ep_collision.h>
#include <ep_init.h>
#include <ep_movili.h>
#include <ep_movili_cazzu.h>
#include <ep_physics.h>
#include <ep_utils.h>
#include <sim_stop.h>
#include <model.h>
#include <stdlib.h>

struct CALModel3D* u_modellu;
struct CALRun3D* a_simulazioni;
struct Substates Q;
CALreal elapsed_time;

void transizioniGlobali(struct CALModel3D* modello)
{
  calApplyElementaryProcess3D(modello, resetF);
  calUpdate3D(modello);

  calApplyElementaryProcess3D(modello,collision);
  calUpdate3D(modello);

  calApplyElementaryProcess3D(modello,movili);
  calUpdate3D(modello);

  calApplyElementaryProcess3D(modello,moviliCazzu);
  calUpdate3D(modello);

  elapsed_time += DELTA_T;

#ifdef VERBOSE
  printSummary(modello);
#endif
}

void partilu()
{
  u_modellu = calCADef3D(X_CELLS,Y_CELLS,Z_CELLS,CAL_MOORE_NEIGHBORHOOD_3D,CAL_SPACE_TOROIDAL,CAL_NO_OPT);

  Q.Fx = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.Fy = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.Fz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.rx = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.ry = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.rz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.vx = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.vy = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.vz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.imove = (struct CALSubstate3Di**)malloc(sizeof(struct CALSubstate3Di*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);

  for(int slot=0;slot<MAX_NUMBER_OF_PARTICLES_PER_CELL;slot++)
    {
      Q.Fx[slot]    = calAddSubstate3Dr(u_modellu);
      Q.Fy[slot]    = calAddSubstate3Dr(u_modellu);
      Q.Fz[slot]    = calAddSubstate3Dr(u_modellu);
      Q.rx[slot]    = calAddSubstate3Dr(u_modellu);
      Q.ry[slot]    = calAddSubstate3Dr(u_modellu);
      Q.rz[slot]    = calAddSubstate3Dr(u_modellu);
      Q.vx[slot]    = calAddSubstate3Dr(u_modellu);
      Q.vy[slot]    = calAddSubstate3Dr(u_modellu);
      Q.vz[slot]    = calAddSubstate3Dr(u_modellu);
      Q.imove[slot] = calAddSubstate3Di(u_modellu);

      calInitSubstate3Dr(u_modellu,Q.Fx[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.Fy[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.Fz[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.rx[slot],PARTICLE_NODATA);
      calInitSubstate3Dr(u_modellu,Q.ry[slot],PARTICLE_NODATA);
      calInitSubstate3Dr(u_modellu,Q.rz[slot],PARTICLE_NODATA);
      calInitSubstate3Dr(u_modellu,Q.vx[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.vy[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.vz[slot],0.0);
      calInitSubstate3Di(u_modellu,Q.imove[slot],PARTICLE_ABSENT);
    }

  // Boundary
  calApplyElementaryProcess3D(u_modellu, boundary_cells);

  // Initial conditions
  elapsed_time = 0.0;
  calApplyElementaryProcess3D(u_modellu, mmiscali_nta_cella);
  // ATTENZIONE, QUESTO PROCESSO ELEMENTARE DEVE AVVENIRE IN SEQUENZIALE
  calApplyElementaryProcess3D(u_modellu, cancella_particelle_in_urto);

  // Simulation
  a_simulazioni = calRunDef3D(u_modellu,0,CAL_RUN_LOOP,CAL_UPDATE_IMPLICIT);
  calRunAddGlobalTransitionFunc3D(a_simulazioni, transizioniGlobali);
  calRunAddStopConditionFunc3D(a_simulazioni, caminalu);
}
