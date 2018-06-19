//
// Communication
//
#include <mpi.h>
#include "comm.h"

static int this_node, n_node;

void comm_init(int* p_argc, char*** p_argv)
{
  MPI_Init(p_argc, p_argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &this_node);
  MPI_Comm_size(MPI_COMM_WORLD, &n_node);
}

void comm_finalise()
{
  MPI_Finalize();
}

int comm_this_node()
{
  return this_node;
}

int comm_nnode()
{
  return n_node;
}

/*
int mpi_share_int(int x, MPI_Op op)
{
  int x_global;
  MPI_Reduce(&x, &x_global, 1, MPI_INT, op, 0, MPI_COMM_WORLD);

  MPI_Bcast(&x_global, 1, MPI_INT, 0, MPI_COMM_WORLD); 

  return x_global;
}

int mpi_share_int(int x, MPI_Op op)
{
  int x_global;
  MPI_Reduce(&x, &x_global, 1, MPI_INT, op, 0, MPI_COMM_WORLD);

  MPI_Bcast(&x_global, 1, MPI_INT, 0, MPI_COMM_WORLD); 

  return x_global;
}
*/


