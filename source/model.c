#include <boundary.h>
#include <ep_collision.h>
#include <ep_movili.h>
#include <ep_movili_cazzu.h>
#include <ep_physics.h>
#include <ep_utils.h>
#include <init.h>
#include <model.h>
#include <utils_io.h>
#include <sim_stop.h>
#include <stdlib.h>

struct CALModel3D* u_modellu;
struct CALRun3D* a_simulazioni;
struct Substates Q;
CALint initial_nummber_of_particles;
CALreal elapsed_time;

void updateF(struct CALModel3D* ca)
{
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      calUpdateSubstate3Dr(ca,Q.Fx[slot]);
      calUpdateSubstate3Dr(ca,Q.Fy[slot]);
      calUpdateSubstate3Dr(ca,Q.Fz[slot]);
    }
}

void updateP(struct CALModel3D* ca)
{
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      calUpdateSubstate3Dr(ca,Q.px[slot]);
      calUpdateSubstate3Dr(ca,Q.py[slot]);
      calUpdateSubstate3Dr(ca,Q.pz[slot]);
    }
}

void updateV(struct CALModel3D* ca)
{
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      calUpdateSubstate3Dr(ca,Q.vx[slot]);
      calUpdateSubstate3Dr(ca,Q.vy[slot]);
      calUpdateSubstate3Dr(ca,Q.vz[slot]);
    }
}

void transizioniGlobali(struct CALModel3D* modello)
{
  calApplyElementaryProcess3D(modello, resetF);
  updateF(modello);

  calApplyElementaryProcess3D(modello,inner_collision);
  calApplyElementaryProcess3D(modello,outer_collision);
  updateF(modello);

  calApplyElementaryProcess3D(modello,movili);
  updateP(modello);
  updateV(modello);


  calApplyElementaryProcess3D(modello,moviliCazzu);
  calUpdate3D(modello);

  elapsed_time += DELTA_T;

#ifdef VERBOSE
  printSummary(modello);
#endif

  CALint S = INTEGRITY_CHECK_STEPS;
  if (a_simulazioni->step % S == 0)
    {
      CALint missing_particle = findMissingParticle(modello);
      if (missing_particle)
        {
          printf("ERROR: missing particle with ID %d\n", missing_particle);
          exit(EXIT_FAILURE);
        }
    }
}

void partilu()
{
  u_modellu = calCADef3D(X_CELLS,Y_CELLS,Z_CELLS,CAL_MOORE_NEIGHBORHOOD_3D,CAL_SPACE_TOROIDAL,CAL_NO_OPT);

  Q.Fx = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.Fy = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.Fz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.px = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.py = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.pz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.vx = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.vy = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.vz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.ID = (struct CALSubstate3Di**)malloc(sizeof(struct CALSubstate3Di*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);

  for(int slot=0;slot<MAX_NUMBER_OF_PARTICLES_PER_CELL;slot++)
    {
      Q.Fx[slot] = calAddSubstate3Dr(u_modellu);
      Q.Fy[slot] = calAddSubstate3Dr(u_modellu);
      Q.Fz[slot] = calAddSubstate3Dr(u_modellu);
      Q.px[slot] = calAddSubstate3Dr(u_modellu);
      Q.py[slot] = calAddSubstate3Dr(u_modellu);
      Q.pz[slot] = calAddSubstate3Dr(u_modellu);
      Q.vx[slot] = calAddSubstate3Dr(u_modellu);
      Q.vy[slot] = calAddSubstate3Dr(u_modellu);
      Q.vz[slot] = calAddSubstate3Dr(u_modellu);
      Q.ID[slot] = calAddSubstate3Di(u_modellu);

      calInitSubstate3Dr(u_modellu,Q.Fx[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.Fy[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.Fz[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.px[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.py[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.pz[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.vx[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.vy[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.vz[slot],0.0);
      calInitSubstate3Di(u_modellu,Q.ID[slot],NULL_ID);
    }

  // Boundary
  boundaryCellsSerial(u_modellu);

  // Initial conditions
  initial_nummber_of_particles = 0;
  elapsed_time = 0.0;
  mmiscali_nta_cella_seriale(u_modellu);
  cancella_particelle_in_urto(u_modellu);

  // Simulation
  a_simulazioni = calRunDef3D(u_modellu,0,CAL_RUN_LOOP,CAL_UPDATE_IMPLICIT);
  calRunAddGlobalTransitionFunc3D(a_simulazioni, transizioniGlobali);
  calRunAddStopConditionFunc3D(a_simulazioni, caminalu);

#ifdef VERBOSE
  printf("The 3D particles computational model\n");
#ifdef OMP
  printf("OpenMP parallel execution enabled!\n");
#endif
  printf("A system of %d particles will be simulated for %f s, subdivided in %d steps, each one corresponding to %f s\n", initial_nummber_of_particles, TOTAL_SIMULATION_TIME, STEPS, DELTA_T);
#endif
}
