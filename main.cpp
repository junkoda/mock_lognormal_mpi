#include <cstdio>
#include "comm.h"
#include "input_power.h"

int main(int argc, char* argv[])
{
  comm_init(&argc, &argv);

  InputPower ps= InputPower("planck_matterpower.dat");

  //int myrank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  //if(myrank == 0) {
  //  printf("Hello World!\n");
  //}

  comm_finalise();
  
  return 0;
}
