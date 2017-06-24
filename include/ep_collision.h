#ifndef EP_COLLISION_H
#define EP_COLLISION_H

#include <model.h>

void inner_collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z);
void outer_collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z);
void boundary_collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z);
//void collision(struct CALModel3D* ca, int x_cell, int y_cell, int z_cell);

#endif
