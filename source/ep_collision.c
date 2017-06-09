#include <ep_collision.h>
#include <ep_utils.h>
#include <math.h>

#include <model.h>

CALreal distance (CALreal* p0, CALreal* p1)
{
  return sqrt((p0[0]-p1[0])*(p0[0]-p1[0]) +
              (p0[1]-p1[1])*(p0[1]-p1[1]) +
              (p0[2]-p1[2])*(p0[2]-p1[2]));
}

CALreal distancePointToPlane (CALreal* p0, CALreal* p1)
{
  // z = CELL_SIDE plane
  if (p1[2] == CELL_SIDE)
    return sqrt((p0[2]-p1[2])*(p0[2]-p1[2]));

  return CELL_SIDE;
}

void collision(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  if (!ncestiArmenuNaParticella(ca, cell_x, cell_y, cell_z, 0))
    return;

  CALreal d = CELL_SIDE;
  CALreal p0[3], p1[3];

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      p0[0] = calGet3Dr(ca, Q.px[slot],cell_x,cell_y,cell_z);
      p0[1] = calGet3Dr(ca, Q.py[slot],cell_x,cell_y,cell_z);
      p0[2] = calGet3Dr(ca, Q.pz[slot],cell_x,cell_y,cell_z);

      for (int n=1; n<ca->sizeof_X; n++)
        {
          // particle-plane collision

          if (calGetX3Di(ca, Q.imove[0], cell_x,cell_y, cell_z, n) == PARTICLE_BORDER)
            {
              if (calGetX3Dr(ca, Q.pz[0],cell_x,cell_y,cell_z,n) == CELL_SIDE)  // z = CELL_SIDE plane
                {
                  p1[0] = calGetX3Dr(ca, Q.px[0],cell_x,cell_y,cell_z,n);
                  p1[1] = calGetX3Dr(ca, Q.py[0],cell_x,cell_y,cell_z,n);
                  p1[2] = calGetX3Dr(ca, Q.pz[0],cell_x,cell_y,cell_z,n);

                  d = distancePointToPlane (p0, p1);
                  //d = fabs(p0[2]-p1[2]);

                  CALreal particle_radious = PARTICLE_RADIUS;

                  if (d < particle_radious)
                    {
                      CALreal vz = calGet3Dr(ca, Q.vz[slot],cell_x,cell_y,cell_z);
                      calSet3Dr(ca, Q.vz[slot],cell_x,cell_y,cell_z, -vz);
                      //break;
                    }
                }

            }

          // particle-particle collision

        }
    }

}

