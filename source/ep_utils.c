#include <ep_utils.h>
#include <math.h>

CALint number_of_particles = 0;
CALreal total_energy = 0.0;

void summary(struct CALModel3D* ca, int cell_x, int cell_y, int cell_z)
{
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    if (calGet3Di(ca, Q.imove[slot],cell_x,cell_y,cell_z) == PARTICLE_PRESENT)
      {
        number_of_particles++;
        total_energy += PARTICLE_MASS*G*calGet3Dr(ca, Q.rz[slot], cell_x, cell_y, cell_z);
        total_energy += 0.5*PARTICLE_MASS*calGet3Dr(ca, Q.vx[slot], cell_x, cell_y, cell_z)*calGet3Dr(ca, Q.vx[slot], cell_x, cell_y, cell_z);
        total_energy += 0.5*PARTICLE_MASS*calGet3Dr(ca, Q.vy[slot], cell_x, cell_y, cell_z)*calGet3Dr(ca, Q.vy[slot], cell_x, cell_y, cell_z);
        total_energy += 0.5*PARTICLE_MASS*calGet3Dr(ca, Q.vz[slot], cell_x, cell_y, cell_z)*calGet3Dr(ca, Q.vz[slot], cell_x, cell_y, cell_z);
      }
}

void printSummary(struct CALModel3D* ca)
{
  number_of_particles = 0;
  total_energy = 0.0;

  calApplyElementaryProcess3D(ca,summary);
  printf("Step %d, elapsed_time %.3f s, number of particles: %d, totoal energy: %.3f\n", a_simulazioni->step, elapsed_time, number_of_particles, total_energy);
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

