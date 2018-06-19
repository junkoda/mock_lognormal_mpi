#ifndef HDF5_WRITE
#define HDF5_WRITE 1

void hdf5_write_grid_real(const char filename[],
			  double const * const fx,
			  const int nc);

void hdf5_write_grid_complex(const char filename[],
			     Grid const * const grid);

#endif
