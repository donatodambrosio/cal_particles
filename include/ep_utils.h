#ifndef EP_UTILS_H
#define EP_UTILS_H

#include <model.h>

void printSummary(struct CALModel3D* ca);

CALbyte ncestiArmenuNaParticella(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z, int n);

CALreal distance (CALreal* p0, CALreal* p1);

#endif
