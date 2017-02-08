#ifndef PROVA_H_
#define PROVA_H_

#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DRun.h>
#include <OpenCAL/cal3DIO.h>
#include <OpenCAL/cal3DUnsafe.h>
#include <GL/glut.h>
#include <stdlib.h>


#define NUMBER_OF_PARTICLES 14988

// Domain dimensions in m
#define X 0.10
#define Y 0.10
#define Z 0.05

// Cell side  in m
#define CL 0.001

//Cell side for divisions
#define CLD 1000

// Domain dimensions in rows, columns and layers
#define ROWS    (int)((X)*(CLD))+1
#define COLS    (int)((Y)*(CLD))+1
#define LAYERS  (int)((Z)*(CLD))+1

//Sottostati
#define MAX_NUMBER_OF_PARTICLES_PER_CELL 3
#define NODATA -10 // No particle condition (used in px, py and pz)

struct Substates
{
	struct CALSubstate3Dr *px[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *py[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *pz[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *vx[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *vy[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *vz[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Di *imove[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	//struct CALSubstate3Di *numPresenti;
};

// Main objcts
extern struct CALModel3D* modello;
extern struct Substates Q;
extern struct CALRun3D* simulazione;

// Computational steps
#define STEPS 1
#define VERBOSE

// Functions
void transition(struct CALModel3D*,int,int,int);
void startModello();
void initFunction();
void finalizeModel();
void setPosition(const double, const double, const double,const CALint,const int);
void partilu();

#endif /* PROVA_H_ */
