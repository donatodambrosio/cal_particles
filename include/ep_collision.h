#ifndef EP_COLLISION_H
#define EP_COLLISION_H

#include <model.h>

void inner_collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z);
void outer_collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z);

#endif
