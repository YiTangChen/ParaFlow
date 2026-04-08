#include "CurvilinearGrid.h" //added by lijie


#include "OSUFlow.h"
#include "Plot3DReader.h" //added by lijie
#include "MPASOReader.h"


#pragma warning(disable : 4251 4100 4244 4101)

OSUFlow::OSUFlow()
{
  bUseRandomSeeds = false;
  flowName = NULL;
  meshName = NULL;
  dataset_files = NULL;
  num_dataset_files = 0;
  flowField = NULL;
  bStaticFlow = true;
  seedPtr = NULL; 
  seedTimeArray = NULL; 
  nSeeds = 0; 
  pStreakLine = NULL; 
  has_data = false; 
  deferred_load_case = -1; 
  order = FOURTH;
}

OSUFlow::~OSUFlow()
{
  delete[] flowName;
  if (has_data)
    delete flowField;
//   if (seedPtr!=NULL) delete[] seedPtr; 
  if (seedTimeArray!=NULL) delete [] seedTimeArray; 
  flowName = NULL;
  flowField = NULL;
}

void OSUFlow::DeleteData(void)
{
	delete flowField; 
	has_data = false; 
}

////////////////////////////////////////////////////////////////////////
// 
// Load the whole static or time-varying data set
//
void OSUFlow::LoadData(const char* fname, bool bStatic, bool deferred)
{
	flowName = new char[255];
	strcpy(flowName, fname);
	bStaticFlow = bStatic;
	
	has_data = false; 

	if(bStaticFlow) {
		numTimesteps = 1; 
		MinT = MaxT = 0; 
		if (deferred == true) {
		  deferred_load_case = 0; 
		  return; 
		}
		InitStaticFlowField();
	}
	else {
		if (deferred == true) {
		  deferred_load_case = 0; 
		  return; 
		}
		InitTimeVaryingFlowField();
	}
	has_data = true; 
}

/////////////////////////////////////////////////////////////
//
//   Read a partial data set 
//   sMin/sMax are local min and max range of the data that are held 
//
void OSUFlow::LoadData(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax, bool deferred)  
{
	flowName = new char[255];
	strcpy(flowName, fname);
	bStaticFlow = bStatic;
	lMin = sMin; lMax = sMax; 
	has_data = false; 

	if(bStaticFlow) {
	  numTimesteps = 1; 
	  MinT = MaxT = 0; 
	  if (deferred == true) {
	    deferred_load_case = 1; 
	    return; 
	  }
	  InitStaticFlowField(sMin, sMax);
	}
	else {
	  if (deferred == true) {
	    deferred_load_case = 1; 
	    return; 
	  }
	  InitTimeVaryingFlowField(sMin, sMax); 
	}
	has_data = true; 
}

/////////////////////////////////////////////////////////////
//
//  Load a partial time-varying data set 
//  sMin/sMax are local min and max range of the data that are held
//  t_min/t_max are the time range (for time-varying field) 
//
void OSUFlow::LoadData(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax, int min_t, int max_t, 
		       bool deferred)  
{
	flowName = new char[255];
	strcpy(flowName, fname);

	bStaticFlow = bStatic;
	lMin = sMin; lMax = sMax; 
	has_data = false; 

	if (max_t >= min_t) {
	  numTimesteps = max_t-min_t+1; 
	  MinT = min_t; MaxT = max_t; 
	}
	else {   //exception. goes back to default 
	  numTimesteps = 1; 
	  MinT = MaxT = min_t; 
	}
	  
	if(bStaticFlow) {  // ignore the time range 
	  numTimesteps = 1; 
	  MinT = MaxT = 0; 
	  if (deferred == true) {
	    deferred_load_case = 2; 
	    return; 
	  }
	  InitStaticFlowField(sMin, sMax);
	}
	else {
	  if (deferred == true) {
	    deferred_load_case = 2; 
	    return; 
	  }
	  InitTimeVaryingFlowField(sMin, sMax, min_t, max_t); 
	}
	has_data = true; 
}

//--------------------------------------------------------------------------
//
// Load a static or time-varying data set 
//
// fname: dataset file name
// sMin, sMax: corners of subdomain (inclusive, node-centered)
// ie, the range [0, 10] contains 10 voxels
// dim: size of total domain (leave uninitialized of mode == 1)
// t_min, t_max: temporal extents (inclusive, node-centered)
// ie, first timestep has extent [0, 1]
// t_min + 1 = tmax for static, tmin + 1 < tmax for time-varying
// mode: 0 = raw with no header data
//       1 = raw with dimensions at start
//       2 = netCDF
//       3 = HDF5
// data: data that has already been read in (defaults to NULL)
// currently only modes 0, 1, and 2 are implemented
//
void OSUFlow::LoadData(const char* fname, double *sMin, double *sMax, 
		       double *dim, int min_t, int max_t, DataMode mode,
		       double **data) {
  
  // init
  flowName = new char[255];
  strcpy(flowName, fname);
  bStaticFlow = false;
  lMin.Set(sMin[0], sMin[1], sMin[2]);
  lMax.Set(sMax[0], sMax[1], sMax[2]);
  has_data = false; 

  if (max_t >= min_t) {
    numTimesteps = max_t - min_t + 1; 
    MinT = min_t; MaxT = max_t; 
  }
  else { // defaults to 1 time step
    numTimesteps = 1; 
    MinT = MaxT = min_t; 
  }
	  
  switch (mode) {

  case 0:
  case 1:
  case 2:
    InitFlowField(sMin, sMax, dim, min_t, max_t, mode, data); 
    break;
  default:
    fprintf(stderr, "Error: LoadData() currently does not support modes other than 0, 1, and 2\n");
    break;

  }

  has_data = true; 

}
//--------------------------------------------------------------------------
//
// Load a static or time-varying data set 
//
// dataset_files: dataset timestep files
// num_dataset_files: number of timestep files
// sMin, sMax: corners of subdomain (inclusive, node-centered)
// ie, the range [0, 10] contains 10 voxels
// dim: size of total domain (leave uninitialized of mode == 1)
// t_min, t_max: temporal extents (inclusive, node-centered)
// ie, first timestep has extent [0, 1]
// t_min + 1 = tmax for static, tmin + 1 < tmax for time-varying
// mode: 0 = raw with no header data
//       1 = raw with dimensions at start
//       2 = netCDF
//       3 = HDF5
// data: data that has already been read in (defaults to NULL)
// currently only modes 0, 1, and 2 are implemented
//
void OSUFlow::LoadData(char **dataset_files, int num_dataset_files,
		       double *sMin, double *sMax, double *dim, int min_t, 
		       int max_t, DataMode mode, double **data) {
  
  // init
  this->dataset_files = dataset_files;
  this->num_dataset_files = num_dataset_files;
  bStaticFlow = false;
  lMin.Set(sMin[0], sMin[1], sMin[2]);
  lMax.Set(sMax[0], sMax[1], sMax[2]);
  has_data = false; 

  if (max_t >= min_t) {
    numTimesteps = max_t - min_t + 1; 
    MinT = min_t; MaxT = max_t; 
  }
  else { // defaults to 1 time step
    numTimesteps = 1; 
    MinT = MaxT = min_t; 
  }
	  
  switch (mode) {

  case 0:
  case 1:
  case 2:
    InitFlowField(sMin, sMax, dim, min_t, max_t, mode, data); 
    break;
  default:
    fprintf(stderr, "Error: LoadData() currently does not support modes other than 0, 1, and 2\n");
    break;

  }

  has_data = true; 

}
//--------------------------------------------------------------------------

bool OSUFlow::DeferredLoadData() 
{
  if (deferred_load_case == -1) return(false); 

  switch(deferred_load_case) {
  case 0: 
    if(bStaticFlow) {
      InitStaticFlowField();
    }
    else
      InitTimeVaryingFlowField();
    has_data = true; 
    break; 
  case 1: 
    if(bStaticFlow) {
      InitStaticFlowField(lMin, lMax);
    }
    else
      InitTimeVaryingFlowField(lMin, lMax); 
    has_data = true; 
    break; 
  case 2: 
    if(bStaticFlow) {  // ignore the time range 
      InitStaticFlowField(lMin, lMax);
    }
    else
      InitTimeVaryingFlowField(lMin, lMax, MinT, MaxT); 
    has_data = true; 
    break; 
  }
  ScaleField(10.0); 
  return(true); 
}

/////////////////////////////////////////////////////////////////
//
//Read the whole datafile and create a vectorfield object based on that 
//
void OSUFlow::InitStaticFlowField(void)
{
  int dimension[3]; 
  double* pData = NULL;

  printf("Initialize Static Flow Field\n");
	
  pData = ReadStaticDataRaw(flowName, dimension); 
	
  // update the domain bounds 
  gMin.Set(0.0, 0.0, 0.0);
  gMax.Set((double)(dimension[0]-1), (double)(dimension[1]-1), (double)(dimension[2]-1));
  
  lMin.Set(0.0, 0.0, 0.0);
  lMax.Set((double)(dimension[0]-1), (double)(dimension[1]-1), (double)(dimension[2]-1));

  double minB[3], maxB[3]; 
	
  minB[0] = minB[1] = minB[2] = 0; 
  maxB[0] = dimension[0]-1; 
  maxB[1] = dimension[1]-1; 
  maxB[2] = dimension[2]-1; 

  flowField = CreateStaticFlowField(pData, dimension[0], dimension[1], 
				    dimension[2],  minB, maxB); 
  
  if(!flowField->IsNormalized())
    flowField->NormalizeField(true);
}

/////////////////////////////////////////////////////////////////
//
// Read only a subdomain and create a vectorfield object based on that 
//
void OSUFlow::InitStaticFlowField(VECTOR3 sMin, VECTOR3 sMax)
{
	int dimension[3], totalNum;
	double* pData = NULL;
	int lxdim, lydim, lzdim; 
	double minB[3], maxB[3]; 

	minB[0] = sMin[0]; minB[1] = sMin[1]; minB[2] = sMin[2];
	maxB[0] = sMax[0]; maxB[1] = sMax[1]; maxB[2] = sMax[2];

	pData = ReadStaticDataRaw(flowName, dimension, minB, maxB); 

	gMin.Set(0.0,0.0,0.0); 
	gMax.Set((double)(dimension[0]-1), (double)(dimension[1]-1), (double)(dimension[2]-1));
	lMin = sMin; lMax = sMax; //local data min/max range

	dimension[0] = (int)(maxB[0]-minB[0]+1); 
	dimension[1] = (int)(maxB[1]-minB[1]+1); 
	dimension[2] = (int)(maxB[2]-minB[2]+1); 
	
	flowField = CreateStaticFlowField(pData, dimension[0], dimension[1], dimension[2], 
			      minB,maxB); 

	// 2023-05-07 Add by Yi-Tang Chen
	if(!flowField->IsNormalized())
    	flowField->NormalizeField(true);
}


///////////////////////////////////////////////////////////////////////
//
// read the whole time sequence and create a time-varying vector field 
void OSUFlow:: InitTimeVaryingFlowField(void)
{
	int n_timesteps;
	int dimension[3]; 
	double** ppData = NULL;
	double minB[3], maxB[3]; 

	ppData = ReadTimeVaryingDataRaw(flowName, n_timesteps, dimension); 

	// update the global bounds of the field. Assuming the same for all time steps 
	gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dimensions 
	gMax.Set((double)(dimension[0]-1), (double)(dimension[1]-1), (double)(dimension[2]-1));
        lMin = gMin; lMax = gMax; 

	minB[0]=minB[1]=minB[2] = 0; 

	maxB[0]= dimension[0]-1; 
	maxB[1]= dimension[1]-1; 
	maxB[2]= dimension[2]-1; 
       
	flowField = CreateTimeVaryingFlowField(ppData, dimension[0], dimension[1], 
					       dimension[2], minB, maxB, 
					       0, n_timesteps-1);  
}

/////////////////////////////////////////////////////////////////
//
// read the whole time sequence within the subdomain 
//
void OSUFlow:: InitTimeVaryingFlowField(VECTOR3 minB, VECTOR3 maxB)
{
  	int n_timesteps;
	int dimension[3]; 
	double** ppData = NULL;
	double min_B[3], max_B[3]; 

	min_B[0]= minB[0]; 
	min_B[1]= minB[1]; 
	min_B[2]= minB[2]; 

	max_B[0]= maxB[0]; 
	max_B[1]= maxB[1]; 
	max_B[2]= maxB[2]; 

	ppData = ReadTimeVaryingDataRaw(flowName, n_timesteps, dimension, 
				     min_B, max_B, -1, -1); 

	// update the global bounds of the field. Assuming the same for all time steps 
	gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dimensions 
	gMax.Set((double)(dimension[0]-1), (double)(dimension[1]-1), (double)(dimension[2]-1));
        lMin = minB; lMax = maxB; 

	dimension[0] = (int)(maxB[0]-minB[0]+1); 
	dimension[1] = (int)(maxB[1]-minB[1]+1); 
	dimension[2] = (int)(maxB[2]-minB[2]+1); 
	
	flowField = CreateTimeVaryingFlowField(ppData, dimension[0], dimension[1], 
					       dimension[2], min_B, max_B, 
					       0, n_timesteps-1);  
}


void OSUFlow:: InitTimeVaryingFlowField(VECTOR3 minB, VECTOR3 maxB, int min_t, int max_t)
{
  	int n_timesteps;
	int dimension[3]; 
	double** ppData = NULL;
	double min_B[3], max_B[3]; 

	min_B[0]= minB[0]; 
	min_B[1]= minB[1]; 
	min_B[2]= minB[2]; 

	max_B[0]= maxB[0]; 
	max_B[1]= maxB[1]; 
	max_B[2]= maxB[2]; 

	ppData = ReadTimeVaryingDataRaw(flowName, n_timesteps, dimension, 
				     min_B, max_B, min_t, max_t); 

	// update the global bounds of the field. Assuming the same for all time steps 
	gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dimensions 
	gMax.Set((double)(dimension[0]-1), (double)(dimension[1]-1), (double)(dimension[2]-1));
        lMin = minB; lMax = maxB; 

	dimension[0] = (int)(maxB[0]-minB[0]+1); 
	dimension[1] = (int)(maxB[1]-minB[1]+1); 
	dimension[2] = (int)(maxB[2]-minB[2]+1); 
	
	flowField = CreateTimeVaryingFlowField(ppData, dimension[0], dimension[1], 
					       dimension[2], min_B, max_B, 
					       min_t, max_t);  

}

//--------------------------------------------------------------------------
//
// Initialize a static or time-varing flow field object
//
// sMin, sMax: corners of subdomain (inclusive, node-centered)
// ie, the range [0, 10] contains 10 voxels
// dim: size of total domain (leave uninitialized if read_dims == true)
// t_min, t_max: temporal extents (inclusive, node-centered)
// ie, first timestep has extent [0, 1]
// t_min + 1 = tmax for static, tmin + 1 < tmax for time-varying
// dm: data mode
// data: data that has already been loaded (defaults to NULL)
//
void OSUFlow::InitFlowField(double *sMin, double *sMax, 
				       double *dim, int t_min, int t_max,
				       DataMode dm, double **data) {

  double** ppData = NULL;

  if (dataset_files != NULL)
    DatasetFiles(dataset_files, num_dataset_files);
  if (data == NULL)
    ppData = ReadData(flowName, dim, sMin, sMax, t_min, t_max, dm); 
  else
    ppData = (double **)data;

  // update the global bounds of the field (same for all time steps)
  gMin.Set(0.0,0.0,0.0);
  gMax.Set(dim[0], dim[1], dim[2]);
  lMin.Set(sMin[0], sMin[1], sMin[2]);
  lMax.Set(sMax[0], sMax[1], sMax[2]);
  MinT = t_min;
  MaxT = t_max; 

  // create the flow field	
  if (t_max == t_min)
    flowField = CreateStaticFlowField(ppData[0], 
				      (int)(sMax[0] - sMin[0] + 1),
				      (int)(sMax[1] - sMin[1] + 1),
				      (int)(sMax[2] - sMin[2] + 1),
				      sMin, sMax);  
  else
    flowField = CreateTimeVaryingFlowField(ppData, 
					   (int)(sMax[0] - sMin[0] + 1),
					   (int)(sMax[1] - sMin[1] + 1),
					   (int)(sMax[2] - sMin[2] + 1),
					   sMin, sMax, t_min, t_max);  

}
//--------------------------------------------------------------------------

///////////////////////////////////////////////
//
//   Create a static flow field 
//   Input: (double) vector data, minB, maxB
//   note this flow field can be a subfield so 
//   minB[0/1/2] may not be zero, and maxB[0/1/2]
//   may not be the max of the entire field 
//   
//
void OSUFlow::CreateStaticFlowField(double* pData, VECTOR3 minB, 
				  VECTOR3 maxB)
{
	int dimension[3], totalNum;

	dimension[0] = (int)(maxB[0]-minB[0]+1); 
	dimension[1] = (int)(maxB[1]-minB[1]+1); 
	dimension[2] = (int)(maxB[2]-minB[2]+1); 

	totalNum = dimension[0] * dimension[1] * dimension[2];

	// create field object
	Solution* pSolution;
	RegularCartesianGrid* pRegularCGrid;
	VECTOR3* pVector;
	VECTOR3** ppVector;
	pVector = new VECTOR3[totalNum];

	for(int iFor = 0; iFor < totalNum; iFor++)
		pVector[iFor].Set(pData[iFor*3], pData[iFor*3+1], pData[iFor*3+2]);
	delete [] pData; 
	ppVector = new VECTOR3*[1];
	ppVector[0] = pVector;
	pSolution = new Solution(ppVector, totalNum, 1);
	pRegularCGrid = new RegularCartesianGrid(dimension[0], dimension[1], dimension[2]);
	lMin = minB; lMax = maxB; //local data min/max range
	pRegularCGrid->SetBoundary(lMin, lMax);
	assert(pSolution != NULL && pRegularCGrid != NULL);
	flowField = new CVectorField(pRegularCGrid, pSolution, 1);
	if(!flowField->IsNormalized())
		flowField->NormalizeField(true);
	delete []pVector; 
	delete[] ppVector; 
}

///////////////////////////////////////////////////////////////////
//
//   pData is assumed to contain the vector data
//
//   xdim, ydim, zdim are the resolutions of the data block
//   minB and maxB specify the actual physical bounds of the block 
//
CVectorField* OSUFlow::CreateStaticFlowField(double *pData, 
				    int xdim, int ydim, int zdim, 
				    double* minB, double* maxB) 
{
  CVectorField* field; 
  Solution* pSolution;
  RegularCartesianGrid* pRegularCGrid;
  VECTOR3* pVector;
  VECTOR3** ppVector;
  VECTOR3 min_b, max_b; 

  int totalNum = xdim*ydim*zdim; 
  pVector = new VECTOR3[totalNum]; 
  ppVector = new VECTOR3*[1]; 

  for(int i=0; i<totalNum; i++) {
    pVector[i].Set(pData[i*3], pData[i*3+1], pData[i*3+2]); 
  }
  delete [] pData; 
  ppVector[0] = pVector; 

  pSolution = new Solution(ppVector, totalNum, 1);
  pRegularCGrid = new RegularCartesianGrid(xdim, ydim, zdim); 

  min_b[0] = minB[0]; min_b[1] = minB[1]; min_b[2] = minB[2]; 
  max_b[0] = maxB[0]; max_b[1] = maxB[1]; max_b[2] = maxB[2]; 

  pRegularCGrid->SetBoundary(min_b, max_b);

  assert(pSolution != NULL && pRegularCGrid != NULL);
  
  field = new CVectorField(pRegularCGrid, pSolution, 1);

  delete [] ppVector; 

  flowField = field; 
  return(field); 
}

//////////////////////////////////////////////////////////////////
//
//   ppData are assumed to contain the vector data of 
//   max_t-min_t+1 time steps 
//
//   xdim, ydim, zdim are the resolutions of the data block
//   minB and maxB specify the actual physical bounds of the block 
//   min_t and max_t are the time interval that this block represents 
//
CVectorField* OSUFlow::CreateTimeVaryingFlowField(double** ppData, 
						  int xdim, int ydim, int zdim, 
						  double* minB, double* maxB, 
						  int min_t, int max_t) 
{
  CVectorField* field; 
  Solution* pSolution;
  RegularCartesianGrid* pRegularCGrid;
  VECTOR3* pVector;
  VECTOR3** ppVector;
  VECTOR3 min_b, max_b; 
  
  int totalNum = xdim*ydim*zdim; 

  numTimesteps = max_t-min_t+1;  
  ppVector = new VECTOR3 *[numTimesteps]; 
  for (int i=0; i<numTimesteps; i++) {
    pVector = new VECTOR3[totalNum]; 
    if (ppData[i] == NULL) 
      printf("panic: null data pointer in CreateTimeVaryingFlowField\n"); 
    for (int j = 0; j<totalNum; j++) {
      pVector[j].Set(ppData[i][j*3], ppData[i][j*3+1], ppData[i][j*3+2]); 
    }
#ifdef USE_BIL
	// BIL is using malloc to allocate the blocks, so use free here.
    free(ppData[i]); 
#else
    delete[] ppData[i];
#endif
    ppVector[i] = pVector; 
  }
#ifdef USE_BIL
  // When BIL is being used, the blocks being passed to OSUFlow are being
  // allocated with malloc, so use free instead of delete.
  free(ppData);
#else
  delete[] ppData;
#endif
  min_b[0] = minB[0]; min_b[1] = minB[1]; min_b[2] = minB[2]; 
  max_b[0] = maxB[0]; max_b[1] = maxB[1]; max_b[2] = maxB[2]; 
  // create the flow field now 
  pSolution = new Solution(ppVector, totalNum, numTimesteps, min_t, max_t);
  pRegularCGrid = new RegularCartesianGrid(xdim, ydim, zdim); 			

  pRegularCGrid->SetBoundary(min_b, max_b);

  assert(pSolution != NULL && pRegularCGrid != NULL);
  field = new CVectorField(pRegularCGrid, pSolution, numTimesteps, min_t);

  for (int i=0; i<numTimesteps; i++) {
    delete [] ppVector[i]; 
  }
  delete[] ppVector; 

  flowField = field; 
  return(field); 
}


//////////////////////////////////////////////////////////////////////////
// specify a set of seed points randomly generated over the specified
// spatial interval. Points can be in axis aligned dimension 0, 1, 2, 3
//////////////////////////////////////////////////////////////////////////
void OSUFlow::SetRandomSeedPoints(const double min[3], 
				  const double max[3], 
				  int num)
{
	for(int iFor = 0; iFor < 3; iFor++)
	{
		minRakeExt[iFor] = min[iFor];
		maxRakeExt[iFor] = max[iFor];
	}
	this->numSeeds[0] = num;
	this->numSeeds[1] = 1;
	this->numSeeds[2] = 1;

	bUseRandomSeeds = true;

	// generate seeds
	nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];

	if (seedPtr!=NULL) delete[]seedPtr; 
	seedPtr = new VECTOR3[nSeeds];

	std::size_t SeedSize[3]; 
	SeedSize[0] = numSeeds[0]; 
	SeedSize[1] = numSeeds[1]; 
	SeedSize[2] = numSeeds[2]; 

	SeedGenerator* pSeedGenerator = new SeedGenerator((const double*)minRakeExt, 
							  (const double*)maxRakeExt, 
							  (const std::size_t*)SeedSize); 
	pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
	delete pSeedGenerator;

}

//////////////////////////////////////////////////////////////////////////
// specify a set of seed points regularly generated over the specified
// spatial interval. Points can be in axis aligned dimension 0, 1, 2, 3
//////////////////////////////////////////////////////////////////////////
void OSUFlow::SetRegularSeedPoints(const double min[3], 
								   const double max[3], 
				   const std::size_t numSeeds[3])
{
	for(int iFor = 0; iFor < 3; iFor++)
	{
		minRakeExt[iFor] = min[iFor];
		maxRakeExt[iFor] = max[iFor];
	}
	this->numSeeds[0] = (unsigned int)numSeeds[0];
	this->numSeeds[1] = (unsigned int)numSeeds[1];
	this->numSeeds[2] = (unsigned int)numSeeds[2];

	bUseRandomSeeds = false;

	// generate seeds
	nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
	if (seedPtr!=NULL) delete[]seedPtr; 
	seedPtr = new VECTOR3[nSeeds];

	SeedGenerator* pSeedGenerator = new SeedGenerator((const double*)minRakeExt, 
							  (const double*)maxRakeExt, 
							  (const std::size_t*)numSeeds);
	pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
	delete pSeedGenerator;
}

//////////////////////////////////////////////////////////////////////////
// explicitly set the values of seeds
//////////////////////////////////////////////////////////////////////////
void OSUFlow::SetSeedPoints(VECTOR3* seeds, int num_seeds)
{
  this->numSeeds[0] = num_seeds;
  this->numSeeds[1] = 1;
  this->numSeeds[2] = 1;
  nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];

  if (seedPtr != NULL)
  {
    delete [] seedPtr; 
  }
  seedPtr = new VECTOR3[nSeeds];

  for(int i=0; i<num_seeds; i++)
  {
    seedPtr[i].Set(seeds[i][0], seeds[i][1], seeds[i][2]);
  }
}

void OSUFlow::SetIntegrationParams(double initStepSize, double maxStepSize)
{
	initialStepSize = initStepSize;
	this->maxStepSize = maxStepSize;
}

void OSUFlow::SetIntegrationOrder(INTEG_ORD inOrder)
{
	this->order = inOrder;
}

//////////////////////////////////////////////////////////////////////////
// generate streamlines
// input
// listSeedTraces: STL list to keep the advection result
// maxPoints: how many advection points each streamline
// randomSeed: seed for random number generator
//////////////////////////////////////////////////////////////////////////
bool OSUFlow::GenStreamLines(list<vtListSeedTrace*>& listSeedTraces, 
			     TRACE_DIR traceDir,
			     int maxPoints,
			     unsigned int randomSeed)
{
	if (has_data == false) DeferredLoadData(); 

	// first generate seeds
    if (seedPtr==NULL)  {
		nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
		seedPtr = new VECTOR3[nSeeds];
		SeedGenerator* pSeedGenerator = new SeedGenerator((const double*)minRakeExt, 
							    (const double*)maxRakeExt, 
							    (const std::size_t*)numSeeds);
		pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
		delete pSeedGenerator;
	}
	// or, seeds have already been generated previously, do nothing. 
	listSeedTraces.clear();

	// execute streamline
	vtCStreamLine* pStreamLine;
	double currentT = 0.0;
	pStreamLine = new vtCStreamLine(flowField);
	switch(traceDir)
	{
	case BACKWARD_DIR:
		pStreamLine->setForwardTracing(false);
		break;
	case FORWARD_DIR:
		pStreamLine->setBackwardTracing(false);
		break;
	case BACKWARD_AND_FORWARD:
		break;
	default:
		break;
	}
	pStreamLine->SetLowerUpperAngle(3.0, 15.0);
	pStreamLine->setMaxPoints(maxPoints);
	pStreamLine->setSeedPoints(seedPtr, nSeeds, currentT);
	pStreamLine->SetInitStepSize(initialStepSize);
	pStreamLine->SetMaxStepSize(maxStepSize);
	pStreamLine->setIntegrationOrder(this->order);
	pStreamLine->execute((void *)&currentT, listSeedTraces);
	// release resource
	delete pStreamLine;
	return true;
}




/////////////////////////////////////////////////////////
// Generate streamlines from the given seeds
bool OSUFlow::GenStreamLines(VECTOR3* seeds, 
			     TRACE_DIR traceDir,
			     const int seedNum,
			     const int maxPoints, 
			     list<vtListSeedTrace*>& listSeedTraces,
           int64_t *seedIds,
           list<int64_t> *listSeedIds)
{
  if (has_data == false) DeferredLoadData(); 

//   if (seedPtr!=NULL) delete [] seedPtr; 
        nSeeds = seedNum; 
	seedPtr = seeds; 

	listSeedTraces.clear();
	if (listSeedIds != NULL)
		(*listSeedIds).clear();

	// execute streamline
	vtCStreamLine* pStreamLine;
	double currentT = 0.0;  // always starts from the default time 0 
	pStreamLine = new vtCStreamLine(flowField);
	switch(traceDir)
	{
	case BACKWARD_DIR:
		pStreamLine->setForwardTracing(false);
		break;
	case FORWARD_DIR:
		pStreamLine->setBackwardTracing(false);
		break;
	case BACKWARD_AND_FORWARD:
		break;
	default:
		break;
	}
	pStreamLine->SetLowerUpperAngle(3.0, 15.0);
	pStreamLine->setMaxPoints(maxPoints);
	pStreamLine->setSeedPoints(seedPtr, nSeeds, currentT, seedIds);
	pStreamLine->SetInitStepSize(initialStepSize);
	pStreamLine->SetMaxStepSize(maxStepSize);
	pStreamLine->setIntegrationOrder(this->order);
	pStreamLine->execute((void *)&currentT, listSeedTraces, listSeedIds);
	// release resource
	delete pStreamLine;
	return true;
}

//////////////////////////////////////////////////////////////////////////////
//   all particles start at the same time: currentT 
bool OSUFlow::GenPathLines(list<vtListTimeSeedTrace*>& listSeedTraces, 
			   TIME_DIR  dir, 
			   int maxPoints,
			   double currentT)
{

  if (has_data == false) DeferredLoadData(); 

	// first generate seeds if not exist before 
        if (seedPtr==NULL)  {
	  nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
	  seedPtr = new VECTOR3[nSeeds];
	  SeedGenerator* pSeedGenerator = new SeedGenerator((const double*)minRakeExt, 
							    (const double*)maxRakeExt, 
							    (const std::size_t*)numSeeds);
	  pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
	  delete pSeedGenerator;
	}

	listSeedTraces.clear();

	// execute streamline
	vtCPathLine* pPathLine;

	pPathLine = new vtCPathLine(flowField);

	pPathLine->SetTimeDir(dir); 

	pPathLine->SetLowerUpperAngle(3.0, 15.0);
	pPathLine->setMaxPoints(maxPoints);
	pPathLine->setSeedPoints(seedPtr, nSeeds, currentT);
	pPathLine->SetInitStepSize(initialStepSize);
	pPathLine->SetMaxStepSize(maxStepSize);
	pPathLine->setIntegrationOrder(FOURTH);
	pPathLine->execute(listSeedTraces);
	// release resource
	delete pPathLine;
	return true;
}

////////////////////////////////////////////////////////////////////////////
//
// Take an input list of seeds, which all start at the same time currentT 
bool OSUFlow::GenPathLines(VECTOR3* seeds, list<vtListTimeSeedTrace*>& listSeedTraces, 
			   TIME_DIR  dir, 
			   int num_seeds, 
			   int maxPoints,
			   double currentT)
{

  if (has_data == false) DeferredLoadData(); 
  if (seedPtr!=NULL) delete [] seedPtr; 

        seedPtr = seeds; 
	nSeeds = num_seeds; 

	listSeedTraces.clear();

	// execute streamline
	vtCPathLine* pPathLine;

	pPathLine = new vtCPathLine(flowField);

	pPathLine->SetTimeDir(dir); 

	pPathLine->SetLowerUpperAngle(3.0, 15.0);
	pPathLine->setMaxPoints(maxPoints);
	pPathLine->setSeedPoints(seedPtr, nSeeds, currentT);
	pPathLine->SetInitStepSize(initialStepSize);
	pPathLine->SetMaxStepSize(maxStepSize);
	pPathLine->setIntegrationOrder(FOURTH);
	pPathLine->execute(listSeedTraces);
	// release resource
	delete pPathLine;
	return true;
}

////////////////////////////////////////////////////////////////////////////
//
//    Take an input list of seeds, which 
//    can start from different times (tarray) 
//
bool OSUFlow::GenPathLines(VECTOR3* seeds, list<vtListTimeSeedTrace*>& listSeedTraces, 
			   TIME_DIR  dir, 
			   int num_seeds, 
			   int maxPoints,
			   double* tarray)
{

  if (has_data == false) DeferredLoadData(); 

	// MOD-BY-LEETEN 02/07/2011-FROM:
		// if (seedPtr!=NULL) delete[] seedPtr; 
	// TO:
	if ( seedPtr != seeds )  
		delete[] seedPtr; 
	// MOD-BY-LEETEN 02/07/2011-END

        nSeeds = num_seeds; 
	seedPtr = seeds; 
	seedTimeArray = tarray; 

	listSeedTraces.clear();

	// execute streamline
	vtCPathLine* pPathLine;

	pPathLine = new vtCPathLine(flowField);

	pPathLine->SetTimeDir(dir); 
	pPathLine->SetLowerUpperAngle(3.0, 15.0);
	pPathLine->setMaxPoints(maxPoints);
	pPathLine->setSeedPoints(seedPtr, nSeeds, tarray);
	pPathLine->SetInitStepSize(initialStepSize);
	pPathLine->SetMaxStepSize(maxStepSize);
	pPathLine->setIntegrationOrder(FOURTH);
	pPathLine->execute(listSeedTraces);
	// release resource
	delete pPathLine;
	return true;
}

////////////////////////////////////////////////////////////////////////////
//
//    Take an input list of seeds (VECTOR4) consisting of 
//    positions and times 
//
bool OSUFlow::GenPathLines(VECTOR4* seeds, 
			   list<vtListTimeSeedTrace*>& listSeedTraces, 
			   TIME_DIR  dir, 
			   int num_seeds, 
			   int maxPoints,
				 int64_t *seedIds,
				 list<int64_t> *listSeedIds)
{

	assert((seedIds == NULL && listSeedIds == NULL) ||
			(seedIds != NULL && listSeedIds != NULL));

  if (has_data == false)
    DeferredLoadData(); 

  nSeeds = num_seeds; 
  if (seedPtr!=NULL) delete [] seedPtr; 
  seedPtr = new VECTOR3[nSeeds]; 
  if (seedTimeArray!=NULL) delete [] seedTimeArray; 
  seedTimeArray = new double[nSeeds]; 
  for (int i=0; i<nSeeds; i++) {
    seedPtr[i][0] = seeds[i][0]; 
    seedPtr[i][1] = seeds[i][1]; 
    seedPtr[i][2] = seeds[i][2]; 
    seedTimeArray[i] = seeds[i][3]; 
  }

  listSeedTraces.clear();

  // execute streamline
  vtCPathLine* pPathLine;

  pPathLine = new vtCPathLine(flowField);
  pPathLine->SetTimeDir(dir); 
  pPathLine->SetLowerUpperAngle(3.0, 15.0);
  pPathLine->setMaxPoints(maxPoints);
  pPathLine->setSeedPoints(seeds, nSeeds, seedIds);
  pPathLine->SetInitStepSize(initialStepSize);
  pPathLine->SetMaxStepSize(maxStepSize);
  pPathLine->setIntegrationOrder(FOURTH);
  pPathLine->execute(listSeedTraces, listSeedIds);
  // release resource
  delete pPathLine;
  return true;
}

///////////////////////////////////////////////////////////////
//
// Use preset streakline seeds, all starting from current_time 
//
bool OSUFlow::GenStreakLines(vtStreakTraces& streakTraces, TIME_DIR dir, 
			     double current_time)
{

  if (has_data == false) DeferredLoadData(); 

  if (seedPtr==NULL)  {
    nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
    seedPtr = new VECTOR3[nSeeds];
    SeedGenerator* pSeedGenerator = new SeedGenerator((const double*)minRakeExt, 
						      (const double*)maxRakeExt, 
						      (const std::size_t*)numSeeds);
    pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
    delete pSeedGenerator;
  }

  pStreakLine = new vtCStreakLine(flowField); 

  double currentT = current_time; 
  
  pStreakLine->SetTimeDir(dir); 
  pStreakLine->SetLowerUpperAngle(3.0, 15.0);
  pStreakLine->setSeedPoints(seedPtr, nSeeds, currentT);
  pStreakLine->SetInitStepSize(initialStepSize);
  pStreakLine->SetMaxStepSize(maxStepSize);
  pStreakLine->setIntegrationOrder(FOURTH);
  pStreakLine->execute((void*) &currentT, streakTraces); 
  delete pStreakLine; 
  return true; 
} 


bool OSUFlow::GenStreakLines(VECTOR3* seeds, vtStreakTraces& streakTraces, TIME_DIR dir,
			     int num_seeds, double current_time)
{

  if (has_data == false) DeferredLoadData(); 

  nSeeds = num_seeds; 
  seedPtr = seeds; 

  pStreakLine = new vtCStreakLine(flowField); 

  //otherwise one sterakline has already been created before 
  double currentT = current_time; 
  
  pStreakLine->SetTimeDir(dir); 
  pStreakLine->SetLowerUpperAngle(3.0, 15.0);
  pStreakLine->setSeedPoints(seedPtr, nSeeds, currentT);
  pStreakLine->SetInitStepSize(initialStepSize);
  pStreakLine->SetMaxStepSize(maxStepSize);
  pStreakLine->setIntegrationOrder(FOURTH);

  pStreakLine->execute((void*) &currentT, streakTraces); 
  delete pStreakLine; 
  return true; 
} 


//------------------------------------------------------------------------
//
// Error()
// mpi error handler
//
void Error(const char *fmt, ...){

  va_list argp;
  vfprintf(stderr, fmt, argp);
  sleep(5);
  exit(0);

}

//code added for curvilinear grid
//added by lijie for curvilinear grid
void OSUFlow::LoadDataCurvilinear(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax)  
{
	flowName = new char[255];
	strcpy(flowName, fname);

	bStaticFlow = bStatic;

	if(bStaticFlow)
	  InitStaticCurvilinearFlowField(sMin, sMax);
	else
	  InitTimeVaryingCurvilinearFlowField(); // to be implemented 
}

//added by lijie for irregular grid
void OSUFlow::LoadDataIrregular(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax)  
{
	flowName = new char[255];
	strcpy(flowName, fname);

	bStaticFlow = bStatic;

	if(bStaticFlow)
	  InitStaticIrregularFlowField(sMin, sMax);
//	else
//	  InitTimeVaryingCurvilinearFlowField(); // to be implemented 
}


//added by lijie to read static curvilinear grid
void OSUFlow::InitStaticCurvilinearFlowField(VECTOR3 sMin, VECTOR3 sMax)
{

	int dimension[3], totalNum;
	int lxdim, lydim, lzdim; 

	char fsz[255],gsz[255];
	sprintf(fsz,"%s.soln",flowName);
	sprintf(gsz,"%s.grid",flowName);
	CPlot3DReader* reader=new CPlot3DReader(fsz,gsz,"");

	Solution* pSol=reader->CreateSolution();
	CurvilinearGrid* pGrid=reader->CreateCurvilinearGrid();

        InitStaticFlowField(pSol, pGrid, sMin, sMax); 
	delete reader;
}
//added by lijie to read static irregular grid
void OSUFlow::InitStaticIrregularFlowField(VECTOR3 sMin, VECTOR3 sMax)
{
        int dimension[3], totalNum;
	int lxdim, lydim, lzdim; 

	char fsz[255],gsz[255],hsz[255];
	sprintf(fsz,"%s.soln",flowName);
	sprintf(gsz,"%s.grid",flowName);
	sprintf(hsz,"%s.tetra",flowName);
	CPlot3DReader* reader=new CPlot3DReader(fsz,gsz,hsz);

	Solution* pSol=reader->CreateSolution();
	IrregularGrid* pGrid=reader->CreateIrregularGrid(false,false);

       InitStaticFlowField(pSol, pGrid, sMin, sMax); 
	delete reader;
}
//added by lijie
void OSUFlow::InitStaticFlowField(Solution* pSolution, Grid* pGrid, VECTOR3 minB, 
				  VECTOR3 maxB)
{
	
	int dimension[3], totalNum;

	dimension[0] = maxB[0]-minB[0]+1; 
	dimension[1] = maxB[1]-minB[1]+1; 
	dimension[2] = maxB[2]-minB[2]+1; 

	totalNum = dimension[0] * dimension[1] * dimension[2];

	// create field object

	lMin = minB; lMax = maxB; //local data min/max range
	//pGrid->SetBoundary(lMin, lMax);
	assert(pSolution != NULL && pGrid != NULL);

	flowField = new CVectorField(pGrid, pSolution, 1);
	if(!flowField->IsNormalized())
		flowField->NormalizeField(true);
}
//added by lijie 
void OSUFlow:: InitTimeVaryingCurvilinearFlowField(void)
{
	FILE *fIn;
	FILE *fVecIn;
	int timesteps;
	char* filename = new char[100];
	int dimension[3], totalNum, tStep;
	double** ppData = NULL;
	double* pData = NULL;
	VECTOR3** ppVector;
	Solution* pSolution;
	CurvilinearGrid* pCurvilinearGrid;
	
	char fsz[255],gsz[255];
	sprintf(fsz,"%s",flowName);//solution info file
	sprintf(gsz,"%s.grid",flowName);
	CPlot3DReader* reader=new CPlot3DReader(fsz,gsz,"");
	reader->m_nTimevarying=1;
	pSolution = reader->CreateSolution();
	pCurvilinearGrid = reader->CreateCurvilinearGrid();
	// set the boundary of physical grid
	VECTOR3 minB, maxB;
	reader->GetDims(dimension);
	minB.Set(0.0, 0.0, 0.0);
	maxB.Set((double)dimension[0], (double)dimension[1], (double)dimension[2]);
	//pRegularCGrid->SetBoundary(minB, maxB);
	assert(pSolution != NULL && pCurvilinearGrid != NULL);
	timesteps=reader->m_nTimeSteps;
	flowField = new CVectorField(pCurvilinearGrid, pSolution, timesteps);
	delete reader;
}
//added by lijie
void OSUFlow::UserPickedSeedPoints(VECTOR3* ptr, int num)
{
	if (seedPtr!=NULL)  {
			delete [] seedPtr;
		}
	seedPtr=new VECTOR3 [num];
	for(int i=0; i<num; i++)
	{
		seedPtr[i].Set(ptr[i].x(),ptr[i].y(),ptr[i].z());
	}

	nSeeds=num;
	this->numSeeds[0] = num;
	this->numSeeds[1] = 1;
	this->numSeeds[2] = 1;
}


void OSUFlow::LoadMPASOData(const char* fname)
{
	flowName = new char[255];
	strcpy(flowName, fname);

	bStaticFlow = false;
	InitMPASOFlowField();
}

void OSUFlow::InitMPASOFlowField(void)
{
	mpaso_reader = new MPASOReader(this->flowName);
	MPASOGrid* grid = mpaso_reader->CreateMPASOGrid();
	Solution *pSolution = nullptr;
	Solution *vSolution = nullptr;
	mpaso_reader->InitSolutions(grid, pSolution, vSolution);

	int timesteps = grid->getNTimestepsLoaded();
	assert(grid != NULL);
	assert(pSolution != NULL && vSolution != NULL);
	flowField = new CVectorField(grid, pSolution, vSolution, timesteps);
}

void OSUFlow::LoadMPASOData(const char* fname, int areaId, int nTotalCell, int* areaIndices, 
							std::vector<int> &neighborInd, int* &LC2GC, int &nLC, int &nVLvl, int* &GC2LC, int &nGC, TimeVaryingDataConfig &cfg)
{
	flowName = new char[255];
	strcpy(flowName, fname);

	bStaticFlow = false;
	InitMPASOFlowField(areaId, nTotalCell, areaIndices, neighborInd, LC2GC, nLC, nVLvl, GC2LC, nGC, cfg);
}

void OSUFlow::LoadMPASOData(const char* meshfilename, const char* datafilename, int areaId, int nTotalCell, int* areaIndices, 
							std::vector<int> &neighborInd, int* &LC2GC, int &nLC, int &nVLvl, int* &GC2LC, int &nGC, TimeVaryingDataConfig &cfg)
{
	this->flowName = new char[255];
	strcpy(this->flowName, datafilename);

	this->meshName = new char[255];
	strcpy(this->meshName, meshfilename);

	bStaticFlow = false;
	InitMPASOFlowField(areaId, nTotalCell, areaIndices, neighborInd, LC2GC, nLC, nVLvl, GC2LC, nGC, cfg);
}

void OSUFlow::InitMPASOFlowField(int areaId, int nTotalCell, int* areaIndices, std::vector<int> &neighborInd, 
								int* &LC2GC, int &nLC, int &nVLvl, int* &GC2LC, int &nGC, TimeVaryingDataConfig &cfg)
{
	if (this->meshName == NULL || *this->meshName == '\0')
		mpaso_reader = new MPASOReader(this->flowName, areaId, nTotalCell, areaIndices, cfg);
	else
		mpaso_reader = new MPASOReader(this->meshName, this->flowName, areaId, nTotalCell, areaIndices, cfg);

	MPASOGrid* grid;
	grid = mpaso_reader->CreateMPASOGrid();
	std::cout << "[OSUFlow::InitMPASOFlowField]: Created MPASOGrid. for areaId = " << areaId << std::endl;
	mpaso_reader->GetNeighborIndices(neighborInd);
	mpaso_reader->GetLocalCell2GlobalCell(LC2GC, nLC);
	mpaso_reader->GetGlobalCell2LocalCell(GC2LC, nGC);
	std::cout << "[OSUFlow::InitMPASOFlowField]: Global to Local Cell mapping retrieved." << std::endl;
	
	Solution *pSolution = nullptr;
	Solution *vSolution = nullptr;
	mpaso_reader->InitSolutions(grid, pSolution, vSolution);
	int timesteps = grid->getNTimestepsLoaded();
	std::cout << "[OSUFlow::InitMPASOFlowField]: Solutions initialized." << std::endl;
	mpaso_reader->GetnVertLevels(nVLvl);
	assert(grid != NULL);
	assert(pSolution != NULL && vSolution != NULL);
	flowField = new CVectorField(grid, pSolution, vSolution, timesteps);
	flowField->SetTimeRange(pSolution->GetRealMinTime(), pSolution->GetRealMaxTime());
	std::cout << "[OSUFlow::InitMPASOFlowField]: CVectorField created." << std::endl;
}

bool OSUFlow::HasMoreTimesteps() const
{
    if (!mpaso_reader) return false;
    return mpaso_reader->getTimestepOffset() < mpaso_reader->getTotalTimesteps();
}

int OSUFlow::GetTotalTimesteps() const
{
    if (!mpaso_reader) return 0;
    return mpaso_reader->getTotalTimesteps();
}

int OSUFlow::GetTimestepOffset() const
{
    if (!mpaso_reader) return 0;
    return mpaso_reader->getTimestepOffset();
}

void OSUFlow::AdvanceTimestep()
{
    assert(mpaso_reader && flowField);
    MPASOGrid* grid = static_cast<MPASOGrid*>(flowField->GetGrid());
    Solution*  pSol = flowField->GetHorizontalSolution();
    Solution*  vSol = flowField->GetVerticalSolution();
    mpaso_reader->UpdateSolutions(grid, pSol, vSol);
    flowField->setHorizontalSolution(pSol);
    flowField->setVerticalSolution(vSol);
    flowField->SetTimeRange(pSol->GetRealMinTime(), pSol->GetRealMaxTime());
}

bool OSUFlow::inBlock(VECTOR3& pos)
{
	// 1. Check if flowField is non-nullptr
    assert(flowField != nullptr && "flowField is non-nullptr");

    // 2. Check if a point is in the block
    VECTOR3 vec;
    int res = flowField->at_phys_truelocal(pos, 0, vec);
    return (res == 1);
}

bool OSUFlow::inBlock_WithGhost(VECTOR3& pos)
{
	// 1. Check if flowField is non-nullptr
    assert(flowField != nullptr && "flowField is non-nullptr");

    // 2. Check if a point is in the block
    VECTOR3 vec;
    int res = flowField->at_phys(pos, 0, vec);
    return (res == 1);
}

//////////////////////////////////////////////////////////////////////////
// generate streamlines
// input
// listSeedTraces: STL list to keep the advection result
// maxPoints: how many advection points each streamline
// randomSeed: seed for random number generator
//////////////////////////////////////////////////////////////////////////
bool OSUFlow::GenStreamLines(list<vtListSeedTrace*>& listSeedTraces, 
			     vector<int>& fromCells,
				 vector<int>& toCells,
				 TRACE_DIR traceDir,
			     int maxPoints,
			     unsigned int randomSeed)
{
	if (has_data == false) DeferredLoadData(); 

	// first generate seeds
    if (seedPtr==NULL)  {
		nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
		seedPtr = new VECTOR3[nSeeds];
		SeedGenerator* pSeedGenerator = new SeedGenerator((const double*)minRakeExt, 
							    (const double*)maxRakeExt, 
							    (const std::size_t*)numSeeds);
		pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
		delete pSeedGenerator;
	}
	// or, seeds have already been generated previously, do nothing. 
	listSeedTraces.clear();

	// execute streamline
	vtCStreamLine* pStreamLine;
	double currentT = 0.0;
	pStreamLine = new vtCStreamLine(flowField);
	switch(traceDir)
	{
	case BACKWARD_DIR:
		pStreamLine->setForwardTracing(false);
		break;
	case FORWARD_DIR:
		pStreamLine->setBackwardTracing(false);
		break;
	case BACKWARD_AND_FORWARD:
		break;
	default:
		break;
	}
	pStreamLine->SetLowerUpperAngle(3.0, 15.0);
	pStreamLine->setMaxPoints(maxPoints);
	pStreamLine->setSeedPoints(seedPtr, nSeeds, currentT);
	pStreamLine->SetInitStepSize(initialStepSize);
	pStreamLine->SetMaxStepSize(maxStepSize);
	pStreamLine->setIntegrationOrder(this->order);
	pStreamLine->setSaveInterval(m_nSaveInterval);
	pStreamLine->execute((void *)&currentT, listSeedTraces, fromCells, toCells, mpaso_reader);
	m_lastTraceStepCounts = pStreamLine->getTraceActualStepCounts();
	// release resource
	delete pStreamLine;
	return true;
}

bool OSUFlow::GenPathLines(VECTOR3* seeds, 
							list<vtListTimeSeedTrace*>& listSeedTraces, 
							TIME_DIR  dir,
							vector<int>& fromCells,
							vector<int>& toCells, 
							int num_seeds, 
							int maxPoints,
							double* tarray)
{

  if (has_data == false) DeferredLoadData(); 

	// MOD-BY-LEETEN 02/07/2011-FROM:
		// if (seedPtr!=NULL) delete[] seedPtr; 
	// TO:
	if ( seedPtr != seeds )  
		delete[] seedPtr; 
	// MOD-BY-LEETEN 02/07/2011-END

        nSeeds = num_seeds; 
	seedPtr = seeds; 
	seedTimeArray = tarray; 

	listSeedTraces.clear();

	// execute streamline
	vtCPathLine* pPathLine;

	pPathLine = new vtCPathLine(flowField);

	pPathLine->SetTimeDir(dir);
	pPathLine->SetLowerUpperAngle(3.0, 15.0);
	pPathLine->setMaxPoints(maxPoints);
	pPathLine->setSeedPoints(seedPtr, nSeeds, tarray);
	pPathLine->SetInitStepSize(initialStepSize);
	pPathLine->SetMaxStepSize(maxStepSize);
	pPathLine->setIntegrationOrder(this->order);
	pPathLine->setSaveInterval(m_nSaveInterval);
	pPathLine->execute(listSeedTraces, fromCells, toCells);
	m_lastTraceStepCounts = pPathLine->getTraceActualStepCounts();
	// We do not own seedPtr or seedTimeArray (both point into the caller's vectors).
	// Null them out so SetSeedPoints / the destructor do not attempt to free them.
	seedPtr = nullptr;
	seedTimeArray = nullptr;
	// release resource
	delete pPathLine;
	return true;
}

void OSUFlow::dumpData(string &prefix)
{
	if (flowField == nullptr) {
		std::cerr << "Flow field is not initialized." << std::endl;
		return;
	}

	// Dump the flow field data to a file or console
	flowField->dumpData(prefix);
}

#ifdef _MPI
//-----------------------------------------------------------------------

#endif
