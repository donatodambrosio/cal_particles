#include <ep_utils.h>
#include <math.h>

CALint number_of_particles = 0;
CALreal total_energy = 0.0;
CALreal max_velocity = 0.0;
CALreal max_displacement = 0.0;

void summary(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  CALreal velocity[3];
  CALreal v = 0.0;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
      {
        number_of_particles++;
        velocity[0] = calGet3Dr(ca, Q.vx[slot], cell_x, cell_y, cell_z);
        velocity[1] = calGet3Dr(ca, Q.vy[slot], cell_x, cell_y, cell_z);
        velocity[2] = calGet3Dr(ca, Q.vz[slot], cell_x, cell_y, cell_z);

        total_energy += PARTICLE_MASS*G*calGet3Dr(ca, Q.rz[slot], cell_x, cell_y, cell_z);
        total_energy += 0.5*PARTICLE_MASS*calGet3Dr(ca, Q.vx[slot], cell_x, cell_y, cell_z)*calGet3Dr(ca, Q.vx[slot], cell_x, cell_y, cell_z);
        total_energy += 0.5*PARTICLE_MASS*calGet3Dr(ca, Q.vy[slot], cell_x, cell_y, cell_z)*calGet3Dr(ca, Q.vy[slot], cell_x, cell_y, cell_z);
        total_energy += 0.5*PARTICLE_MASS*calGet3Dr(ca, Q.vz[slot], cell_x, cell_y, cell_z)*calGet3Dr(ca, Q.vz[slot], cell_x, cell_y, cell_z);

        v = sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);

        if (max_velocity < v)
          {
            max_velocity = v;
            max_displacement = v*DELTA_T;
          }
      }
}

void printSummary(struct CALModel3D* ca)
{
  number_of_particles = 0;
  total_energy = 0.0;
  max_velocity = 0.0;
  max_displacement = 0.0;

  calApplyElementaryProcess3D(ca,summary);
  printf("Step %d, DELTA_T: %.6f, elapsed_time: %.3f s, number_of_particles: %d, totoal_energy: %.9f, max_v: %.6f, max_displacement: %.6f\n", a_simulazioni->step, DELTA_T, elapsed_time, number_of_particles, total_energy, max_velocity, max_displacement);
}

CALbyte ncestiArmenuNaParticella(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z, int n)
{
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGetX3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z,n) == PARTICLE_PRESENT)
      return CAL_TRUE;

  return CAL_FALSE;
}

CALreal distance (CALreal* p0, CALreal* p1)
{
  return sqrt((p0[0]-p1[0])*(p0[0]-p1[0]) +
              (p0[1]-p1[1])*(p0[1]-p1[1]) +
              (p0[2]-p1[2])*(p0[2]-p1[2]));
}

CALreal pointPlaneDistance(CALreal* P0, CALreal* n, CALreal* Pp)
{
  return fabs(n[0]*P0[0] + n[1]*P0[1] + n[2]*P0[2] - n[0]*Pp[0] - n[1]*Pp[1] - n[2]*Pp[2]);
}

CALreal scalar(CALreal *v, CALreal *n)
{
  return v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
}

void reflect(CALreal *v, CALreal *n)
{
  CALreal s = scalar(v, n);
  for (int i=0; i<3; i++)
    v[i] = v[i] - 2*s*n[i];
}

void orthogonalProjectedPointToPlane(CALreal* Pi, CALreal* Pp, CALreal* n, CALreal* Pj)
{
  CALreal d = - scalar(n, Pp);
  CALreal t = - d - scalar(n, Pi);
  for (int i=0; i<3; i++)
    Pj[i] = t * n[i] + Pi[i];
}
