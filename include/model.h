#ifndef MODEL_H
#define MODEL_H

#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DRun.h>
#include <OpenCAL/cal3DIO.h>
#include <OpenCAL/cal3DUnsafe.h>
#include <math.h>

// PHYSICAL FLAGS AND CONSTANTS
// #define AIR_VISCOSITY 1.81e-5
#define GRAVITY
#define G 9.81
#define VISCOELASTIC
#define KN 1000
#define ETHA 0.01

// Particle mass, radius and volume
#define PI 3.14159265358979
#define PARTICLE_MASS 0.001
#define PARTICLE_RADIUS 0.00025
#define PARTICLE_VOLUME ((4.0/3.0)*PI*PARTICLE_RADIUS*PARTICLE_RADIUS*PARTICLE_RADIUS)

// PHYSICAL TIME
//#define DELTA_T 0.001 //[s]
#define DELTA_T (0.05 * sqrt(PARTICLE_MASS/KN)) //[s]

// Cell side [m], volume [m^3] and max occupancy volume [m^3] according to Kepler's conjecture
#define CELL_SIDE (0.002)
#define CELL_VOLUME (CELL_SIDE*CELL_SIDE*CELL_SIDE)
#define KEPLER_OCCUPANCY_FACTOR (0.74)
#define MAX_OCCUPANCY_VOLUME ((KEPLER_OCCUPANCY_FACTOR)*(CELL_VOLUME))

// Max number of particles per cell according to Kepler's conjecture
#define MAX_NUMBER_OF_PARTICLES_PER_CELL  (int)(((MAX_OCCUPANCY_VOLUME)/(PARTICLE_VOLUME))+1)

// Domain dimensions in m and in cells along x, y and z directions
#define X 0.02
#define Y 0.02
#define Z 0.02
#define X_CELLS (int)((X)/(CELL_SIDE))
#define Y_CELLS (int)((Y)/(CELL_SIDE))
#define Z_CELLS (int)((Z)/(CELL_SIDE))

// SLOT (or PARTICLE) IDs
#define BORDER_ID -1            // BORDER ID
#define NULL_ID 0               // NO PARTICLE IN SLOT
#define DEFAULT_PARTICLE_ID 1   // DEFAULT PARTICLE ID

// Particles are randomly distributed on the CELL_FILL_RATE*MAX_NUMBER_OF_PARTICLES_PER_CELL top layers
#define TOP_LAYERS      (Z_CELLS) - 0.4 * (Z_CELLS)
#define CELL_FILL_RATE  0.1 // 0.59 // 1.0/(MAX_NUMBER_OF_PARTICLES_PER_CELL)

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

// Main objcts
extern struct CALModel3D* u_modellu;
extern struct Substates Q;
extern struct CALRun3D* a_simulazioni;
extern CALint initial_nummber_of_particles;
extern CALreal elapsed_time;

// Computational steps
#define STEPS 10

// Verbose mode
#define VERBOSE

// Functions
void partilu();

#endif /* MODEL_H */
