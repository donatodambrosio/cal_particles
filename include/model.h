#ifndef MODEL_H
#define MODEL_H

#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DRun.h>
#include <OpenCAL/cal3DIO.h>
#include <OpenCAL/cal3DUnsafe.h>
#include <GL/glut.h>
#include <stdlib.h>


// Domain dimensions in m
#define X 0.02
#define Y 0.02
#define Z 0.04

// Cell side  in m
#define CELL_SIDE 0.002

// Domain dimensions in cells along x, y and z directions
#define Y_CELLS (int)((Y)/(CELL_SIDE))
#define X_CELLS (int)((X)/(CELL_SIDE))
#define Z_CELLS (int)((Z)/(CELL_SIDE))


#define MAX_NUMBER_OF_PARTICLES_PER_CELL 10
#define PARTICLE_NODATA -9999    // No particle condition (used in px, py and pz)
#define PARTICLE_BORDER -1
#define PARTICLE_ABSENT  0
#define PARTICLE_PRESENT 1


// Particles are randomly distributed on the 10% top layers
#define TOP_LAYERS      (Z_CELLS) - 0.1 * (Z_CELLS)
#define CELL_FILL_RATE  0.5 * (MAX_NUMBER_OF_PARTICLES_PER_CELL)

//Sottostati
struct Substates
{
	struct CALSubstate3Dr *px[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *py[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *pz[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *vx[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *vy[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *vz[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Di *imove[MAX_NUMBER_OF_PARTICLES_PER_CELL];
};

// Main objcts
extern struct CALModel3D* u_modellu;
extern struct Substates Q;
extern struct CALRun3D* a_simulazioni;

// Computational steps
#define STEPS 10

// Verbose mode
#define VERBOSE

// Functions
void partilu();

#endif /* MODEL_H */
