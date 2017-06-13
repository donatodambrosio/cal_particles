#ifndef MODEL_H
#define MODEL_H

#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DRun.h>
#include <OpenCAL/cal3DIO.h>
#include <OpenCAL/cal3DUnsafe.h>
#include <GL/glut.h>
#include <stdlib.h>

// ACCELERATION OF GRAVITY
#define G 9.81
#define AIR_VISCOSITY 1.81e-5

// PHYSICAL TIME
//#define DELTA_T 0.001 //[s]
#define DELTA_T 0.0005 //[s]


// Particle mass, radius and volume
#define PI 3.14159265358979
#define PARTICLE_MASS 0.001
#define PARTICLE_RADIUS (0.0005)
#define PARTICLE_VOLUME ((4.0/3.0)*PI*PARTICLE_RADIUS*PARTICLE_RADIUS*PARTICLE_RADIUS)

// Cell side [m], volume [m^3] and max occupancy volume [m^3] according to Kepler's conjecture
#define CELL_SIDE (0.002)
#define CELL_VOLUME (CELL_SIDE*CELL_SIDE*CELL_SIDE)
#define KEPLER_OCCUPANCY_FACTOR (0.74)
#define MAX_OCCUPANCY_VOLUME ((KEPLER_OCCUPANCY_FACTOR)*(CELL_VOLUME))

// max allowed velocity
// #define V_MAX 0.9*CELL_SIDE/DELTA_T

// Max number of particles per cell according to Kepler's conjecture
#define MAX_NUMBER_OF_PARTICLES_PER_CELL  (int)(((MAX_OCCUPANCY_VOLUME)/(PARTICLE_VOLUME))+1)

// Domain dimensions in m
#define X 0.02
#define Y 0.02
#define Z 0.02

// Domain dimensions in cells along x, y and z directions
#define Y_CELLS (int)((Y)/(CELL_SIDE))
#define X_CELLS (int)((X)/(CELL_SIDE))
#define Z_CELLS (int)((Z)/(CELL_SIDE))

//#define MAX_NUMBER_OF_PARTICLES_PER_CELL 10
#define PARTICLE_NODATA -9999    // No particle condition (used in px, py and pz)
#define PARTICLE_BORDER -1
#define PARTICLE_ABSENT  0
#define PARTICLE_PRESENT 1


// Particles are randomly distributed on the 20% top layers
#define TOP_LAYERS      (Z_CELLS) - 0.4 * (Z_CELLS)
#define CELL_FILL_RATE  1.0 // 0.59 // 1.0/(MAX_NUMBER_OF_PARTICLES_PER_CELL)

//Sottostati
struct Substates
{
  //  struct CALSubstate3Dr *px[MAX_NUMBER_OF_PARTICLES_PER_CELL];
  //	struct CALSubstate3Dr *py[MAX_NUMBER_OF_PARTICLES_PER_CELL];
  //	struct CALSubstate3Dr *pz[MAX_NUMBER_OF_PARTICLES_PER_CELL];
  //	struct CALSubstate3Dr *vx[MAX_NUMBER_OF_PARTICLES_PER_CELL];
  //	struct CALSubstate3Dr *vy[MAX_NUMBER_OF_PARTICLES_PER_CELL];
  //	struct CALSubstate3Dr *vz[MAX_NUMBER_OF_PARTICLES_PER_CELL];
  //	struct CALSubstate3Di *imove[MAX_NUMBER_OF_PARTICLES_PER_CELL];

  struct CALSubstate3Dr **Fx;
  struct CALSubstate3Dr **Fy;
  struct CALSubstate3Dr **Fz;
  struct CALSubstate3Dr **px;
  struct CALSubstate3Dr **py;
  struct CALSubstate3Dr **pz;
  struct CALSubstate3Dr **vx;
  struct CALSubstate3Dr **vy;
  struct CALSubstate3Dr **vz;
  struct CALSubstate3Di **imove;
};

// Main objcts
extern struct CALModel3D* u_modellu;
extern struct Substates Q;
extern struct CALRun3D* a_simulazioni;

// Computational steps
#define STEPS 10000

// Verbose mode
#define VERBOSE

// Functions
void partilu();

#endif /* MODEL_H */
