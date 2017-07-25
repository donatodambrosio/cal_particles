#include <model.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <utils_io.h>
#include <omp.h>

CALbyte simulationStep()
{
  CALbyte again;

  //exectutes the global transition function, the steering function and check for the stop condition.
  //again = calRunCAStep3D(a_simulazioni);
  again = runCAStep3D(u_modellu);

  //simulation main loop
  //a_simulazioni->step++;
  step++;

  return again;
}

int main(int argc, char** argv)
{
  double time_spent = 0.0;

  partilu();

  char t0_path[2048], tf_path[2048];
  strcpy(t0_path, argv[0]);
  t0_path[strlen(t0_path)-7] = '\0';
  strcpy(tf_path, t0_path);
  strcat(t0_path, "data/particles_t0.txt");
  strcat(tf_path, "data/particles_tf.txt");
#ifdef VERBOSE
  printf("argv[0] = %s; t0_path = %s\n", argv[0], t0_path);
#endif

  saveParticles(u_modellu, step, elapsed_time, time_spent, t0_path);

#ifdef OMP
  double begin, end;
  begin = omp_get_wtime();
#else
  clock_t begin = clock();
  clock_t end = begin;
#endif

  //calRun3D(a_simulazioni);
  CALbyte again;
  do
    again = simulationStep();
  while (again);

#ifdef OMP
  end = omp_get_wtime();
  time_spent = end - begin;
#else
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
#endif

  mapperToSubstates3D(u_modellu, Q_current, ID_current);
  saveParticles(u_modellu, step, elapsed_time, time_spent, tf_path);

  return 0;
}
