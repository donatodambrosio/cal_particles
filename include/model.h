#ifndef MODEL_H
#define MODEL_H

#ifdef OMP
#include <OpenCAL-OMP/cal3D.h>
#include <OpenCAL-OMP/cal3DRun.h>
#include <OpenCAL-OMP/cal3DIO.h>
#include <OpenCAL-OMP/cal3DUnsafe.h>
#include <OpenCAL-OMP/cal3DBuffer.h>
#else
#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DRun.h>
#include <OpenCAL/cal3DIO.h>
#include <OpenCAL/cal3DUnsafe.h>
#include <OpenCAL/cal3DBuffer.h>
#endif
#include <math.h>

// PHYSICAL CONSTANTS AND FLAGS
#define KN 100//50000
#define GRAVITY
#define G 9.81
#define VISCOELASTIC
#define ETHA 0.0015

// Particle mass, radius and volume
#define PI 3.14159265358979
#define PARTICLE_MASS 0.0001 //1.0471975511966E-06 // -> density = 2000
#define PARTICLE_RADIUS 0.0005
#define PARTICLE_VOLUME ((4.0/3.0)*PI*PARTICLE_RADIUS*PARTICLE_RADIUS*PARTICLE_RADIUS)

// Cell side [m], volume [m^3] and max occupancy volume [m^3] according to Kepler's conjecture
#define CELL_SIDE (0.002)
#define CELL_VOLUME (CELL_SIDE*CELL_SIDE*CELL_SIDE)
#define KEPLER_OCCUPANCY_FACTOR (0.74)
#define MAX_OCCUPANCY_VOLUME ((KEPLER_OCCUPANCY_FACTOR)*(CELL_VOLUME))

// Max number of particles per cell (slots) according to Kepler's conjecture
#define MAX_NUMBER_OF_PARTICLES_PER_CELL  (int)(((MAX_OCCUPANCY_VOLUME)/(PARTICLE_VOLUME))+1)

// SLOT (or PARTICLE) IDs
#define BORDER_ID -1            // BORDER ID
#define NULL_ID 0               // NO PARTICLE IN SLOT
#define DEFAULT_PARTICLE_ID 1   // DEFAULT PARTICLE ID

// Domain dimensions in m and in cells along x, y and z directions
#define X_SIZE 0.04
#define Y_SIZE 0.04
#define Z_SIZE 0.04
#define X_CELLS (int)((X_SIZE)/(CELL_SIDE))
#define Y_CELLS (int)((Y_SIZE)/(CELL_SIDE))
#define Z_CELLS (int)((Z_SIZE)/(CELL_SIDE))
#define SIZE X_CELLS*Y_CELLS*Z_CELLS

// Particles are randomly distributed on the CELL_FILL_RATE*MAX_NUMBER_OF_PARTICLES_PER_CELL top layers
#define TOP_LAYERS      (Z_CELLS) - 0.8 * (Z_CELLS)
#define CELL_FILL_RATE  0.75 // 0.59 // 1.0/(MAX_NUMBER_OF_PARTICLES_PER_CELL)

// PHYSICAL TIME AND COMPUTATIONAL STEPS
#define TOTAL_SIMULATION_TIME 1.0 //0.1 //[s]
#define DELTA_T (0.1 * sqrt(PARTICLE_MASS/KN)) //[s]
#define STEPS (int)((double)(TOTAL_SIMULATION_TIME)/(double)(DELTA_T))
#define INTEGRITY_CHECK_STEPS STEPS

//SUBSTATES
struct Substates
{
  struct CALSubstate3Dr **Fx;
  struct CALSubstate3Dr **Fy;
  struct CALSubstate3Dr **Fz;
  struct CALSubstate3Dr **px;
  struct CALSubstate3Dr **py;
  struct CALSubstate3Dr **pz;
  struct CALSubstate3Dr **vx;
  struct CALSubstate3Dr **vy;
  struct CALSubstate3Dr **vz;
  struct CALSubstate3Di **ID;
};
/*extern CALreal** Q_Fx_current;
extern CALreal** Q_Fx_next;
extern CALreal** Q_Fy_current;
extern CALreal** Q_Fy_next;
extern CALreal** Q_Fz_current;
extern CALreal** Q_Fz_next;
extern CALreal** Q_px_current;
extern CALreal** Q_px_next;
extern CALreal** Q_py_current;
extern CALreal** Q_py_next;
extern CALreal** Q_pz_current;
extern CALreal** Q_pz_next;
extern CALreal** Q_vx_current;
extern CALreal** Q_vx_next;
extern CALreal** Q_vy_current;
extern CALreal** Q_vy_next;
extern CALreal** Q_vz_current;
extern CALreal** Q_vz_next;
extern CALint** Q_ID_current;
extern CALint** Q_ID_next;
*/
#define REAL_SUBSTATES_NUMBER 9
#define FX 0
#define FY 1
#define FZ 2
#define PX 3
#define PY 4
#define PZ 5
#define VX 6
#define VY 7
#define VZ 8
#define INT_SUBSTATES_NUMBER 1
#define PID 0
extern CALreal* Q_current;
extern CALreal* Q_next;
extern CALint* ID_current;
extern CALint* ID_next;
#define REAL_SUBSTATE(S, slot, substate) ((S)+((slot)*(REAL_SUBSTATES_NUMBER*SIZE)+((substate)*(SIZE))))
#define INT_SUBSTATE(S, slot, substate) ((S)+((slot)*(INT_SUBSTATES_NUMBER*SIZE)+((substate)*(SIZE))))

// Main objcts
extern struct CALModel3D* u_modellu;
extern struct Substates Q;
extern CALint* Xi;
extern CALint* Xj;
extern CALint* Xk;
extern CALint step;
extern CALint initial_nummber_of_particles;
extern CALreal elapsed_time;

// Functions
void partilu();
CALbyte runCAStep3D(struct CALModel3D* modello);
void mapperToSubstates3D(struct CALModel3D *model, CALreal * realSubstate_current_OUT, CALint* intSubstate_current_OUT);
void defaultMapperToSubstates3D(struct CALModel3D *model, CALreal * realSubstate_current_OUT, CALint* intSubstate_current_OUT);

#endif /* MODEL_H */
