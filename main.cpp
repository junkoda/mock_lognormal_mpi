#include <cstdio>
#include <mpi.h>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(myrank == 0) {
    printf("Hello World!\n");
  }
  
  MPI_Finalize();
  
  return 0;
}
