#ifndef COMM_H
#define COMM_H 1

#include <typeinfo> 
#include <mpi.h>

#include "msg.h"
#include "error.h"

void comm_init(int* p_argc, char*** p_argv);
void comm_finalise();

int comm_this_node();
int comm_n_nodes();

static inline MPI_Datatype mpi_datatype(const std::type_info& type_id)
{
  if(type_id == typeid(int))
    return MPI_INT;
  else if(type_id == typeid(long))
    return MPI_LONG;
  else if(type_id == typeid(long long))
    return MPI_LONG_LONG;
  else if(type_id == typeid(float))
    return MPI_FLOAT;
  else if(type_id == typeid(double))
    return MPI_DOUBLE;

  msg_printf(msg_fatal, "Error: unknown data type\n");
  throw RuntimeError();

  return 0;
}

/*
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
*/

template <class T> T comm_bcast(T& x)
{
  MPI_Bcast(&x, 1, mpi_datatype(typeid(T)), 0, MPI_COMM_WORLD);

  return x;
}

template<class T> T comm_partial_sum(T x)
{
  T x_reduced;
  MPI_Scan(&x, &x_reduced, 1, mpi_datatype(typeid(T)),
	   MPI_SUM, MPI_COMM_WORLD);

  return x_reduced;
}

template<class T> T comm_sum(T x)
{
  T x_reduced;
  MPI_Reduce(&x, &x_reduced, 1, mpi_datatype(typeid(T)),
	     MPI_SUM, 0, MPI_COMM_WORLD);

  return x_reduced;
}


#endif
