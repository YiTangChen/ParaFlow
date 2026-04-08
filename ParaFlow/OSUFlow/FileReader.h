#ifdef _MPI
#include <mpi.h>
#endif

#ifndef _FILEREADER_H
#define _FILEREADER_H

typedef enum {
  RAW = 0,
  RAW_HEADER,
  NETCDF,
  HDF_double,
  HDF_DOUBLE,
} DataMode;

void DatasetFiles(char **dataset_files, int num_dataset_files);

double* ReadStaticDataRaw(char* fname, int* dimension); 

double* ReadStaticDataRaw(char *fname, int* dimension, 
			 double* sMin, double* sMax); 

double* ReadStaticDataRaw(char *fname, double *sMin, double* sMax, int* dim);

double** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
			       int *dimension); 

double** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
			       int *dimension, 
			       double *minB, double *maxB, 
			       int min_t, int max_t); 

double** ReadTimeVaryingDataRaw(char *fname, double* sMin, double* sMax, 
			       int* dim, int bt_max, int t_min, int t_max);

double** ReadData(char *fname, double *dim, double *minB, 
		 double *maxB, int min_t, int max_t,
		 DataMode data_mode);

#ifdef _MPI
void Mpi_ioReadDataRaw(MPI_File f, double *dim, double *minB, 
		       double *maxB, double *p, DataMode dm);
#endif

void PosixReadDataRaw(FILE *f, double *dim, double *minB, 
		      double *maxB, double *p, DataMode dm);
void swap4(char *n);

#endif
