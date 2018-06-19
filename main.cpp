#include <cstdio>
#include "comm.h"

int main(int argc, char* argv[])
{
  comm_init(&argc, &argv);

  //int myrank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  //if(myrank == 0) {
  //  printf("Hello World!\n");
  //}

  comm_finalise();
  
  return 0;
}
