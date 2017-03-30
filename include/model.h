#ifndef MODEL_H
#define MODEL_H

#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DRun.h>
#include <OpenCAL/cal3DIO.h>
#include <OpenCAL/cal3DUnsafe.h>
#include <GL/glut.h>
#include <stdlib.h>


#define NUMBER_OF_PARTICLES 14988

// Domain dimensions in m
#define X 0.02
#define Y 0.02
#define Z 0.02

// Cell side  in m
#define CL 0.002

//Cell side for divisions
#define CLD 2000

#define MAX_NUMBER_OF_PARTICLES_PER_CELL 10
#define NODATA -9999    // No particle condition (used in px, py and pz)
#define PARTICLE_EDGE   -1
#define PARTICLE_ABSENT  0
#define PARTICLE_PRESENT 1

// Domain dimensions in rows, columns and layers
#define ROWS    (int)((X)/(CL))
#define COLS    (int)((Y)/(CL))
#define SLICES  (int)((Z)/(CL))

// Particles are randomly distributed on the 20% top layers
#define TOP_LAYERS      0.2 * (SLICES)
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
#define STEPS 30
#define VERBOSE

// Functions
void partilu();

#endif /* PROVA_H_ */
