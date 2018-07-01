#include <cassert>
#include <hdf5.h>

#include "comm.h"
#include "msg.h"
#include "error.h"
#include "grid.h"
#include "config.h"
#include "hdf5_write.h"

namespace {
void write_data_table(hid_t loc, const char name[], 
		      const hsize_t nrow, const hsize_t ncol,
		      const hsize_t stride,
		      const hid_t mem_type, const hid_t save_type,
		      void const * const data);
}
    
void hdf5_write_grid_real(const char filename[],
			     Grid const * const grid)
{
  assert(grid->mode == grid_mode_x);
  const int nc= grid->nc;
  const size_t ncz= 2*(nc/2 + 1);
  double const * const fx= grid->fx;
  
  //H5Eset_auto2(H5E_DEFAULT, NULL, 0);

  // Parallel file access
  hid_t plist= H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);

  hid_t file= H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  if(file < 0) {
    msg_printf(msg_error, "Error: unable to create HDF5 file: %s\n",
	       filename);
    throw IOError();
  }    
  msg_printf(msg_debug, "Created a new HDF5 file: %s\n", filename);

  /*
  hid_t file= H5Fopen(filename, H5F_ACC_RDWR, plist);
  if(file < 0) {
  }
  else {
    msg_printf(msg_debug, "Opened HDF5 file, %s\n", filename);
  }
  */


  // Data structure in memory
  const hsize_t data_size_mem[]= {grid->local_nx, nc, ncz};
  hid_t memspace= H5Screate_simple(3, data_size_mem, 0);

  const hsize_t offset_mem[]= {0, 0, 0};
  const hsize_t count_mem[]= {grid->local_nx, nc, nc};

  H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
  offset_mem, NULL, count_mem, NULL);

  // Data structure in file
  const hsize_t dim= 3;
  const hsize_t data_size_file[]= {nc, nc, nc};
  hid_t filespace= H5Screate_simple(dim, data_size_file, NULL);

  // local subset of data for this node
  const hsize_t offset_file[]= {grid->local_x0, 0, 0, 0};
  const hsize_t count_file[]= {grid->local_nx, nc, nc};


  H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
		      offset_file, NULL, count_file, NULL);


  hid_t dataset= H5Dcreate(file, "fx", H5T_IEEE_F64LE, filespace,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(dataset < 0) {
    throw IOError();
  }


  hid_t plist_data= H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_data, H5FD_MPIO_COLLECTIVE);
    
  const herr_t status= H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
				plist_data, fx);



  H5Pclose(plist_data);
  
  H5Dclose(dataset);

  H5Sclose(filespace);
  H5Sclose(memspace);
  
  assert(status >= 0);

  H5Pclose(plist);
  H5Fclose(file);

  msg_printf(msg_info, "%s written.\n", filename);
}



void hdf5_write_grid_complex(const char filename[],
			     Grid const * const grid)
{
  assert(grid->mode == grid_mode_k);
  const int nc= grid->nc;
  const size_t nckz= nc/2 + 1;
  double const * const fx= grid->fx;
  
  //H5Eset_auto2(H5E_DEFAULT, NULL, 0);

  // Parallel file access
  hid_t plist= H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);

  hid_t file= H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  if(file < 0) {
    msg_printf(msg_error, "Error: unable to create HDF5 file: %s\n",
	       filename);
    throw IOError();
  }    
  msg_printf(msg_debug, "Created a new HDF5 file: %s\n", filename);

  // Data structure in memory
  const hsize_t data_size_mem[]= {grid->local_nx, nc, nckz, 2};
  hid_t memspace= H5Screate_simple(4, data_size_mem, 0);

  // Data structure in file
  const hsize_t dim= 4;
  const hsize_t data_size_file[]= {nc, nc, nckz, 2};
  hid_t filespace= H5Screate_simple(dim, data_size_file, NULL);

  // local subset of data for this node
  const hsize_t offset_file[]= {grid->local_x0, 0, 0, 0};
  const hsize_t count_file[]= {grid->local_nx, nc, nckz, 2};


  H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
		      offset_file, NULL, count_file, NULL);


  hid_t dataset= H5Dcreate(file, "fk", H5T_IEEE_F64LE, filespace,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(dataset < 0) {
    throw IOError();
  }


  hid_t plist_data= H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_data, H5FD_MPIO_COLLECTIVE);
    
  const herr_t status= H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
				plist_data, fx);
  assert(status >= 0);
  
  H5Pclose(plist_data);  
  H5Dclose(dataset);
  H5Sclose(filespace);
  H5Sclose(memspace);

  H5Pclose(plist);
  H5Fclose(file);

  msg_printf(msg_info, "%s written.\n", filename);
}

void hdf5_write_grid_particles(const char filename[],
			       const Particles& particles)
{
  //H5Eset_auto2(H5E_DEFAULT, NULL, 0);

  // Parallel file access
  hid_t plist= H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);

  hid_t file= H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  if(file < 0) {
    msg_printf(msg_error, "Error: unable to create HDF5 file: %s\n",
	       filename);
    throw IOError();
  }    
  msg_printf(msg_debug, "Created a new HDF5 file: %s\n", filename);

  // collect particle information across nodes
  long long np_local= particles.size();
  assert(np_local > 0);

  // Data structure in memory
  assert(sizeof(Particle) % sizeof(double) == 0);
  const int stride= sizeof(Particle)/sizeof(double);

  write_data_table(file, "x", np_local, 3, stride,
		   H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE,
		   particles.front().x);

  write_data_table(file, "v", np_local, 3, stride,
		   H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE,
		   particles.front().v);

  H5Pclose(plist);
  H5Fclose(file);

  msg_printf(msg_info, "particles written to %s.\n", filename);
}

namespace {
void write_data_table(hid_t loc, const char name[], 
		      const hsize_t nrow, const hsize_t ncol,
		      const hsize_t stride,
		      const hid_t mem_type, const hid_t save_type,
		      void const * const data)
{
  // Write an array of a struct to a 2-dimensional table
  
  // Float    FLOAT_MEM_TYPE    FLOAT_SAVE_TYPE (config.h)
  // uint64_t H5T_NATIVE_UINT64 H5T_STD_U64LE 
  
  // Gather inter-node information
  long long offset_ll= comm_partial_sum<long long>(nrow);
  offset_ll -= nrow;

  long long nrow_total= comm_sum<long long>(nrow);
  if(nrow_total == 0) {
    msg_printf(msg_warn, "Warning: zero data given to write_data_table\n"); 
    return;
  }


  // Data structure in memory
  const hsize_t data_size_mem= nrow*stride;
  hid_t memspace= H5Screate_simple(1, &data_size_mem, 0);

  const hsize_t offset_mem= 0;
  const hsize_t size_mem= ncol;
  const hsize_t count_mem= nrow;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
		      &offset_mem, &stride, &count_mem, &size_mem);

  // Data structure in file
  const hsize_t dim= ncol == 1 ? 1 : 2;
  const hsize_t data_size_file[]= {hsize_t(nrow_total), ncol};
  hid_t filespace= H5Screate_simple(dim, data_size_file, NULL);

  // local subset of data for this node
  const hsize_t offset_file[]= {hsize_t(offset_ll), 0};
  const hsize_t count_file[]= {nrow, ncol};

  H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
		      offset_file, NULL, count_file, NULL);

  hid_t dataset= H5Dcreate(loc, name, save_type, filespace,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(dataset < 0)
    return;

  hid_t plist= H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
    
  const herr_t status = H5Dwrite(dataset, mem_type, memspace, filespace,
				 plist, data);
  assert(status >= 0);
  
  H5Pclose(plist);
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dataset);
  

}
} // unnamed namespace
