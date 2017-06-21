#ifndef UTILS_IO_H
#define UTILS_IO_H

#include <model.h>

void saveParticles(struct CALModel3D* ca, CALint step, CALreal elapsed_time, double CPU_time, char* path);

#endif
