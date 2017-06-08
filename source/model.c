#include <ep_boundary.h>
#include <ep_init.h>
#include <ep_movili.h>
#include <ep_movili_cazzu.h>
#include <ep_utils.h>
#include <sim_stop.h>
#include <model.h>

struct CALModel3D* u_modellu;
struct CALRun3D* a_simulazioni;
struct Substates Q;


void transizioniGlobali(struct CALModel3D* modello)
{
  calApplyElementaryProcess3D(modello,movili);
  calUpdate3D(modello);

  calApplyElementaryProcess3D(modello,moviliCazzu);
  calUpdate3D(modello);
}

void partilu()
{
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
