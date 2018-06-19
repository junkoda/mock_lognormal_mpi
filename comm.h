#ifndef COMM_H
#define COMM_H 1

#include <typeinfo> 
#include <mpi.h>

#include "error.h"

void comm_init(int* p_argc, char*** p_argv);
void comm_finalise();

int comm_this_node();
int comm_n_nodes();

template <class T> MPI_Datatype comm_data_type(T x)
{
  if(typeid(x) == typeid(int))
    return MPI_INT;
  else if(typeid(x) == typeid(float))
    return MPI_FLOAT;
  else if(typeid(x) == typeid(double))
    return MPI_DOUBLE;
  else
    throw TypeError();

  return MPI_BYTE;
}

template <class T> T comm_bcast(T x)
{
  MPI_Datatype data_type= comm_data_type(x);
  
  MPI_Bcast(&x, 1, data_type, 0, MPI_COMM_WORLD);

  return x;
}

#endif
