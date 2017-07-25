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
struct Substates Q;
CALreal** Q_Fx_current;
CALreal** Q_Fx_next;
CALreal** Q_Fy_current;
CALreal** Q_Fy_next;
CALreal** Q_Fz_current;
CALreal** Q_Fz_next;
CALreal** Q_px_current;
CALreal** Q_px_next;
CALreal** Q_py_current;
CALreal** Q_py_next;
CALreal** Q_pz_current;
CALreal** Q_pz_next;
CALreal** Q_vx_current;
CALreal** Q_vx_next;
CALreal** Q_vy_current;
CALreal** Q_vy_next;
CALreal** Q_vz_current;
CALreal** Q_vz_next;
CALint** Q_ID_current;
CALint** Q_ID_next;
CALint* Xi = NULL;
CALint* Xj = NULL;
CALint* Xk = NULL;
CALint step;
CALint initial_nummber_of_particles;
CALreal elapsed_time;

void updateF()
{
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      calCopyBuffer3Dr(Q_Fx_next[slot],Q_Fx_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      calCopyBuffer3Dr(Q_Fy_next[slot],Q_Fy_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      calCopyBuffer3Dr(Q_Fz_next[slot],Q_Fz_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
    }
}

void updateP()
{
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      calCopyBuffer3Dr(Q_px_next[slot],Q_px_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      calCopyBuffer3Dr(Q_py_next[slot],Q_py_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      calCopyBuffer3Dr(Q_pz_next[slot],Q_pz_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
    }
}

void updateV()
{
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      calCopyBuffer3Dr(Q_vx_next[slot],Q_vx_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      calCopyBuffer3Dr(Q_vy_next[slot],Q_vy_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      calCopyBuffer3Dr(Q_vz_next[slot],Q_vz_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
    }
}

void updateID()
{
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
      calCopyBuffer3Di(Q_ID_next[slot],Q_ID_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
}

void transizioniGlobali(struct CALModel3D* modello)
{
  for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        resetF(modello, i, j, k);
  updateF();

  for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        inner_collision(modello, i, j, k);

  for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        outer_collision(modello, i, j, k);
  updateF();

  for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        movili(modello, i, j, k);
  updateP();
  updateV();

  for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        moviliCazzu(modello, i, j, k);
  updateF();
  updateP();
  updateV();
  updateID();

  elapsed_time += DELTA_T;

#ifdef VERBOSE
  printSummary(modello);
#endif

  CALint S = INTEGRITY_CHECK_STEPS;
  if (step % S == 0)
    {
      CALint missing_particle = findMissingParticle(modello);
      if (missing_particle)
        {
          printf("ERROR: missing particle with ID %d\n", missing_particle);
          exit(EXIT_FAILURE);
        }
    }
}

CALbyte runCAStep3D(struct CALModel3D* modello)
{
  transizioniGlobali(modello);

  if (caminalu(modello))
    return CAL_FALSE;
  return CAL_TRUE;
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

  Q_Fx_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fx_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fy_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fy_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fz_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fz_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_px_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_px_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_py_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_py_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_pz_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_pz_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vx_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vx_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vy_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vy_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vz_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vz_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_ID_current = (CALint**)malloc(sizeof(CALint*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_ID_next    = (CALint**)malloc(sizeof(CALint*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);

  for(int slot=0;slot<MAX_NUMBER_OF_PARTICLES_PER_CELL;slot++)
    {
      Q_Fx_current[slot] = Q.Fx[slot]->current;
      Q_Fx_next[slot]    = Q.Fx[slot]->next;
      Q_Fy_current[slot] = Q.Fy[slot]->current;
      Q_Fy_next[slot]    = Q.Fy[slot]->next;
      Q_Fz_current[slot] = Q.Fz[slot]->current;
      Q_Fz_next[slot]    = Q.Fz[slot]->next;
      Q_px_current[slot] = Q.px[slot]->current;
      Q_px_next[slot]    = Q.px[slot]->next;
      Q_py_current[slot] = Q.py[slot]->current;
      Q_py_next[slot]    = Q.py[slot]->next;
      Q_pz_current[slot] = Q.pz[slot]->current;
      Q_pz_next[slot]    = Q.pz[slot]->next;
      Q_vx_current[slot] = Q.vx[slot]->current;
      Q_vx_next[slot]    = Q.vx[slot]->next;
      Q_vy_current[slot] = Q.vy[slot]->current;
      Q_vy_next[slot]    = Q.vy[slot]->next;
      Q_vz_current[slot] = Q.vz[slot]->current;
      Q_vz_next[slot]    = Q.vz[slot]->next;
      Q_ID_current[slot] = Q.ID[slot]->current;
      Q_ID_next[slot]    = Q.ID[slot]->next;
    }

  // Boundary
  boundaryCellsSerial(u_modellu);

  // Initial conditions
  initial_nummber_of_particles = 0;
  elapsed_time = 0.0;
  mmiscali_nta_cella_seriale(u_modellu);
  cancella_particelle_in_urto(u_modellu);

  // Simulation step
  step = 1;

  // allocate and build neighborhood arrays
  Xi = (CALint*)malloc(sizeof(CALint)*u_modellu->sizeof_X);
  Xj = (CALint*)malloc(sizeof(CALint)*u_modellu->sizeof_X);
  Xk = (CALint*)malloc(sizeof(CALint)*u_modellu->sizeof_X);
  for (int n=0; n<u_modellu->sizeof_X; n++)
    {
      Xi[n] = u_modellu->X[n].i;
      Xj[n] = u_modellu->X[n].j;
      Xk[n] = u_modellu->X[n].k;
    }


#ifdef VERBOSE
  printf("The 3D particles computational model\n");
#ifdef OMP
  printf("OpenMP parallel execution enabled!\n");
#endif
  printf("A system of %d particles will be simulated for %f s, subdivided in %d steps, each one corresponding to %f s\n", initial_nummber_of_particles, TOTAL_SIMULATION_TIME, STEPS, DELTA_T);
#endif
}
