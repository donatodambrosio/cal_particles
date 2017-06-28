#include <model.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <utils_io.h>

CALbyte simulationStep()
{
  CALbyte again;

  //exectutes the global transition function, the steering function and check for the stop condition.
  again = calRunCAStep3D(a_simulazioni);

  //simulation main loop
  a_simulazioni->step++;

  return again;
}

int main(int argc, char** argv)
{
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
  clock_t begin = clock();
  clock_t end = begin;
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  saveParticles(u_modellu, a_simulazioni->step, elapsed_time, time_spent, t0_path);

  //calRun3D(a_simulazioni);
  CALbyte again; CALreal x = 0;
  do
    again = simulationStep();
  while (again);

  end =  clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  saveParticles(u_modellu, a_simulazioni->step, elapsed_time, time_spent, tf_path);

  return 0;
}
