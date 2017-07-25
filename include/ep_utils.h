#ifndef EP_UTILS_H
#define EP_UTILS_H

#include <model.h>

CALreal distance (CALreal* p0, CALreal* p1);

CALreal pointPlaneDistance(CALreal* P0, CALreal* Pp, CALreal* n);
//void reflect(CALreal* v, CALreal* n);
void orthogonalProjectedPointToPlane(CALreal* Pi, CALreal* Pp, CALreal* n, CALreal* ri);

CALint findMissingParticle(struct CALModel3D* ca);

#endif
