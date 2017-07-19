#include <h5part_io.h>
#include <stdlib.h>

void test_h5part()
{
  int timeStep = 0;
  long long nparticles = 2;
  double *x = (double*)malloc(nparticles*sizeof(double));
  x[0] = 0.0;
  x[1] = 1.0;

  for (int i=0; i<2; i++)
    printf("x[%d]=%f\n", i, x[i]);

  MPI_Comm comm = 0;
  h5_file_t *f = H5OpenFile ("provola.h5part", H5_O_WRONLY, comm);

  H5SetStep(f,timeStep);
  H5PartSetNumParticles(f, nparticles);
  H5PartWriteDataFloat64(f,"x",x);

  H5CloseFile(f);
}


void read_test_h5part()
{
  int timeStep = 0;
  long long nparticles;
  double *x = (double*)malloc(2*sizeof(double));

  MPI_Comm comm = 0;
  h5_file_t *f = H5OpenFile ("./provola.h5part", H5_O_RDONLY, comm);

//  char name[64];
//  int index, nds;
//  nds=H5PartGetNumDatasets(f);
//    for(index=0;index<nds;index++){
//      H5PartGetDatasetName(f,index,name,64);
//      printf("\tDataset[%u] name=[%s]\n", index,name);
//    }
//    return;

  H5SetStep(f,timeStep);
  nparticles = H5PartGetNumParticles(f);
  H5PartReadDataFloat64(f,"x",x);

  H5CloseFile(f);

  for (int i=0; i<2; i++)
    printf("x[%d]=%f\n", i, x[i]);
}
