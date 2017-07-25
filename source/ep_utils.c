#include <ep_utils.h>
#include <math.h>

CALreal distance (CALreal* p0, CALreal* p1)
{
  return sqrt((p0[0]-p1[0])*(p0[0]-p1[0]) +
              (p0[1]-p1[1])*(p0[1]-p1[1]) +
              (p0[2]-p1[2])*(p0[2]-p1[2]));
}

CALreal pointPlaneDistance(CALreal* P0, CALreal* Pp, CALreal* n)
{
  return fabs(n[0]*P0[0] + n[1]*P0[1] + n[2]*P0[2] - n[0]*Pp[0] - n[1]*Pp[1] - n[2]*Pp[2]);
}

CALreal scalar(CALreal *v, CALreal *n)
{
  return v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
}

//void reflect(CALreal *v, CALreal *n)
//{
//  CALreal s = scalar(v, n);
//  for (int i=0; i<3; i++)
//    v[i] = v[i] - 2*s*n[i];
//}

void orthogonalProjectedPointToPlane(CALreal* Pi, CALreal* Pp, CALreal* n, CALreal* Pj)
{
  CALreal d = - scalar(n, Pp);
  CALreal t = - d - scalar(n, Pi);
  for (int i=0; i<3; i++)
    Pj[i] = t * n[i] + Pi[i];
}


CALint findMissingParticle(struct CALModel3D* ca)
{
  CALbyte particle_found = CAL_FALSE;

  for (int id = 0; id < initial_nummber_of_particles; id ++)
    {
      particle_found = CAL_FALSE;
      // find the particle with id = id
      for (int cell_x=0; cell_x<ca->rows; cell_x++)
        for (int cell_y=0; cell_y<ca->columns; cell_y++)
          for (int cell_z = 0; cell_z<ca->slices; cell_z++)
            for (int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
              if (calGet3Di(ca,Q.ID[slot],cell_x,cell_y,cell_z) == id)
                  particle_found = CAL_TRUE;

      if (!particle_found)
        return id;
    }

  return 0; // no missing particles
}
