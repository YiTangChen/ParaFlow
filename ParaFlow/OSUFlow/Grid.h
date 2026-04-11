/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Streamlines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef _GRID_H_
#define _GRID_H_

#include "header.h"
#include "Element.h"
#include "Interpolator.h"

//Silence an annoying and unnecessary compiler warning
//#pragma warning(disable : 4251 4100 4244)

enum CellType
{
	TRIANGLE,
	CUBE,
	POLYGONE,
	TETRAHEDRON,
	VORONOI
};

// define the cell type
enum CellTopoType
{
	T0_CELL,					// vertex
	T1_CELL,					// edge
	T2_CELL,					// triangle, quarilateral
	T3_CELL,					// tetrahedra, cube
	T4_CELL,					// hetrahedra, added by lijie
	T5_CELL						// MPAS Ocean mesh
};

enum SliceType
{
	X_ALIGNED,
	Y_ALIGNED,
	Z_ALIGNED
};

//////////////////////////////////////////////////////////////////////////
//
// base class for grid
//
//////////////////////////////////////////////////////////////////////////
class Grid
{
public:
        Grid() {}; 
	virtual ~Grid(){}; 
	// get the dimension
	virtual void GetDimension(int& xdim, int& ydim, int& zdim) = 0;
	// physical coordinate of vertex verIdx
	virtual bool at_vertex(int verIdx, VECTOR3& pos) = 0;
	// whether the physical point is in the boundary
	virtual bool at_phys(VECTOR3& pos) = 0;			
	// get vertex list of a cell
	virtual int getCellVertices(int cellId, CellTopoType cellType, vector<int>& vVertices) = 0;
	// get the cell id and also interpolating coefficients for the given physical position
	virtual int phys_to_cell(PointInfo& pInfo, double t, int* cachedLowT) = 0;
	// interpolation
	virtual void interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, double* coeff) = 0;
	// the volume of cell
	virtual double cellVolume(int cellId) = 0;
	// type of cell
	virtual CellType GetCellType(void) = 0;
	// get min and maximal boundary
	virtual void Boundary(VECTOR3& minB, VECTOR3& maxB) = 0;
	// set bounding box
	virtual void SetBoundary(VECTOR3& minB, VECTOR3& maxB) = 0;
	// get grid spacing in x,y,z dimensions
	virtual void GetGridSpacing(int cellId, double& xspace, double& yspace, double& zspace) = 0;
	// boundary intersection
	virtual void BoundaryIntersection(VECTOR3& intersectP, VECTOR3& startP, 
					  VECTOR3& endP,double* stepSize, double oldStepSize) = 0;
	virtual void dumpData(string &prefix) {};

protected:
	// reset parameters
	virtual void Reset(void) = 0;
	// compute bounding box
	virtual void ComputeBBox(void) = 0;
	// whether the point is in the bounding box
	virtual bool isInBBox(VECTOR3& pos) = 0;
	// whether in a cell
	virtual bool isInCell(PointInfo& pInfo, const int cellId) = 0;
};

//////////////////////////////////////////////////////////////////////////
//
// Cartesian Grid (Regular and Irregular)
//
//////////////////////////////////////////////////////////////////////////
class CartesianGrid : public Grid
{
public:
	// constructor and destructor
       CartesianGrid(int xdim, int ydim, int zdim);
	CartesianGrid();
	~CartesianGrid();
	inline void GetDimension(int& xdim, int& ydim, int& zdim)
	       {xdim = m_nDimension[0]; ydim = m_nDimension[1]; zdim = m_nDimension[2];}

	// physical coordinate of vertex verIdx
	virtual bool at_vertex(int verIdx, VECTOR3& pos) =0; 
	// whether the physical point is in the boundary
	virtual bool at_phys(VECTOR3& pos) =0; 
	// get vertex list of a cell
	virtual int getCellVertices(int cellId, CellTopoType cellType, vector<int>& vVertices) =0; 
	// get the cell id and also interpolating coefficients for the given physical position
	virtual int phys_to_cell(PointInfo& pInfo, double t, int* cachedLowT) =0;
	
	// interpolation
	virtual void interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, double* coeff) =0; 
	// the volume of cell
	virtual double cellVolume(int cellId) = 0; 
	// type of cell
	virtual CellType GetCellType(void) = 0; 
	// get min and maximal boundary
	virtual void Boundary(VECTOR3& minB, VECTOR3& maxB) = 0; 
	// set bounding box
	virtual void SetBoundary(VECTOR3& minB, VECTOR3& maxB) = 0; 
	// get grid spacing in x,y,z dimensions
	virtual void GetGridSpacing(int cellId, double& xspace, double& yspace, double& zspace) = 0; 
	// boundary intersection
	virtual void BoundaryIntersection(VECTOR3& intersectP, VECTOR3& startP, 
					  VECTOR3& endP,double* stepSize, double oldStepSize) = 0; 


protected:
	// reset parameters
	void Reset(void);
	// dimension related
	inline int xdim(void) { return m_nDimension[0];}
	inline int ydim(void) { return m_nDimension[1];}
	inline int zdim(void) { return m_nDimension[2];}
	inline int xcelldim(void) {return (m_nDimension[0] - 1);}
	inline int ycelldim(void) {return (m_nDimension[1] - 1);}
	inline int zcelldim(void) {return (m_nDimension[2] - 1);}

	int m_nDimension[3];				// dimension
	VECTOR3 m_vMinBound, m_vMaxBound;	// min and maximal boundary
};

//////////////////////////////////////////////////////////////////////////
//
// regular cartesian grid
//
//////////////////////////////////////////////////////////////////////////
// map coordinates in computational space to physical space
#define UCGridPhy2Comp(x, y, f) (((x) - (y))*(f))

class RegularCartesianGrid : public CartesianGrid
{
private:
	double mappingFactorX;				// mapping from physical space to computational space
	double mappingFactorY;
	double mappingFactorZ;
	double oneOvermappingFactorX;
	double oneOvermappingFactorY;
	double oneOvermappingFactorZ;
	double gridSpacing;			        // the minimal grid spacing of all dimensions

public:
	RegularCartesianGrid(int xdim, int ydim, int zdim);
	RegularCartesianGrid();
	~RegularCartesianGrid();
	// physical coordinate of vertex verIdx
	bool at_vertex(int verIdx, VECTOR3& pos);
	// whether the physical point is in the boundary
	bool at_phys(VECTOR3& pos);			
	// get vertex list of a cell
	int getCellVertices(int cellId, CellTopoType cellType, vector<int>& vVertices);
	// get the cell id and also interpolating coefficients for the given physical position
	int phys_to_cell(PointInfo& pInfo, double t, int* cachedLowT);
	// interpolation
	  void interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, double* coeff); 
	// the volume of cell
	double cellVolume(int cellId);
	// cell type
	CellType GetCellType(void) {return CUBE;}
	// set bounding box
	void SetBoundary(VECTOR3& minB, VECTOR3& maxB);
	// get min and maximal boundary
	void Boundary(VECTOR3& minB, VECTOR3& maxB);
	// get grid spacing in x,y,z dimensions
	void GetGridSpacing(int cellId, double& xspace, double& yspace, double& zspace) 
	{ xspace = oneOvermappingFactorX; yspace = oneOvermappingFactorY; zspace = oneOvermappingFactorZ; }
	void BoundaryIntersection(VECTOR3&, VECTOR3&, VECTOR3&, double*, double);

protected:
	void Reset(void);
	// compute bounding box
	void ComputeBBox(void);
	// whether the point is in the bounding box
	bool isInBBox(VECTOR3& pos);
	// whether in a cell
	bool isInCell(PointInfo& pInfo, const int cellId);
};

/* 
//Comment the following code out since they have not been implemented
// 
//////////////////////////////////////////////////////////////////////////
//
//	irregular cartesian grid
//
//////////////////////////////////////////////////////////////////////////
class IrregularCartesianGrid : public CartesianGrid
{
private:
	float* m_pXSpacing;			// space array for x, y, z dimension
	float* m_pYSpacing;
	float* m_pZSpacing;

public:
	IrregularCartesianGrid(int xdim, int ydim, int zdim);
	IrregularCartesianGrid();
	~IrregularCartesianGrid();

protected:
	void Reset(void);
};

//////////////////////////////////////////////////////////////////////////
//
// curvilinear grid
//
//////////////////////////////////////////////////////////////////////////
class CurvilinearGrid : public Grid
{
private:
	int m_nDimension[3];				// dimension

public:
	// constructor and deconstructor
	CurvilinearGrid(int xdim, int ydim, int zdim);
	CurvilinearGrid();
	~CurvilinearGrid();
};

*/ 

//////////////////////////////////////////////////////////////////////////
//
// irregular grid
//
//////////////////////////////////////////////////////////////////////////
class IrregularGrid : public Grid
{
private:
	int m_nNodeNum;						// number of nodes
	int m_nTetraNum;					// number of tetras
	CVertex* m_pVertexGeom;				// geometry of all vertices
	CTetra* m_pTetra;					// tetra
	TetraInfo* m_pTetraInfo;			// pre-computed tetra information
	TVertex* m_pVertexTopo;				// vertex topology
	VECTOR3 m_vMinBound, m_vMaxBound;	// min and maximal boundary
	bool m_bTetraInfoInit;				// whether the tetra information is pre-computed

public:
	// constructor and deconstructor
	IrregularGrid();
	IrregularGrid(int nodeNum, int tetraNum, CVertex* pVertexGeom, CTetra* pTetra, TVertex* pVertexTopo);
	~IrregularGrid();

	// from virtual functions
	void Reset(void);
	void GetDimension(int& xdim, int& ydim, int& zdim) {xdim = m_nNodeNum; ydim = m_nTetraNum; zdim = 0;}
	bool at_vertex(int verIdx, VECTOR3& pos);
	bool at_phys(VECTOR3& pos);
	int getCellVertices(int cellId, CellTopoType cellType, vector<int>& vVertices);
	int phys_to_cell(PointInfo& pInfo, double t, int* cachedLowT);
	void interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, double* coeff);
	double cellVolume(int cellId);
	bool isInCell(PointInfo& pInfo, const int cellId);
	CellType GetCellType(void) {return TETRAHEDRON;}

	void ComputeBBox(void);
	void SetBoundary(VECTOR3& minB, VECTOR3& maxB);
	void Boundary(VECTOR3& minB, VECTOR3& maxB);
	bool isInBBox(VECTOR3& pos);

	// irregular specific functions
	void SetTetraInfoInit(bool bInit);
	bool GetTetraInfoInit(void);
	int nextTetra(PointInfo& pInfo, int tetraId);
	void PreGetP2NMatrix(MATRIX3& m, int cellId);
	bool Physical2NaturalCoord(VECTOR3& nCoord, VECTOR3& pCoord, int cellId);

	void GetGridSpacing(int cellId, double& xspace, double& yspace, double& zspace) {}; 
	virtual void BoundaryIntersection(VECTOR3& intersectP, VECTOR3& startP, 
					  VECTOR3& endP,double* stepSize, double oldStepSize){}; 

};

double getStepSize(VECTOR3& p, VECTOR3& p1, VECTOR3& p2, double oldStepSize);


//////////////////////////////////////////////////////////////////////////
//
// MPAS Ocean grid
//
//////////////////////////////////////////////////////////////////////////
class MPASOGrid : public Grid
{
private:
	// Time independent
	// Store all of the cell and vertex coordinates
    VECTOR3* cellCoord;
    VECTOR3* vertexCoord;
    int* verticesOnCell;
    int* cellsOnVertex;
    int* cellsOnCell;
    int* numVerticesOnCell;
    int* maxLevelCell;
	// Time-varying
	double* zTop;
	std::vector<double> zTopTimestamps;

	// Parameters
	int nCells;
	int nTrueLocalCells;  // local cells that are not ghost cells
	int nLocalVertices;
    int nMaxEdges;
	int nTimestepsLoaded;  // timesteps currently loaded in memory
	int nVertLevels;
	VECTOR3 m_MinBound, m_MaxBound;	// min and maximal boundary
	double lat_min, lat_max;
    double lon_center;      // in [0, 2pi)
    double lon_half_width;  // in [0, pi]
	double earth_radius;

public:
	MPASOGrid();
	~MPASOGrid();

	// from virtual functions
	
	// Might be useless in MPAS-Ocean data
	void GetDimension(int& xdim, int& ydim, int& zdim) {xdim = 0; ydim = 0; zdim = 0;}
	// get the physical coordinate of the vertex
	bool at_vertex(int verIdx, VECTOR3& pos);
	// difficult to determine pos is in the grid or not
	bool at_phys(VECTOR3& pos);
	// What is cellType?
	// Need to know how to map cellId to cell index in distributed version
	int getCellVertices(int nodeId, CellTopoType cellType, vector<int>& vVertices);
	// Stack-array overload: avoids heap allocation on hot path. Returns nVert or -1.
	int getCellVertices(int nodeId, CellTopoType cellType, int* vVertices);
	// get cellId for pInfo.phyCoord
	int phys_to_cell(PointInfo& pInfo, double t, int* cachedLowT);
	int phys_to_truelocalcell(PointInfo& pInfo, double t, int* cachedLowT);
	void interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, double* coeff);
	// Stack-array overload: avoids heap allocation on hot path.
	void interpolate(VECTOR3& nodeData, VECTOR3* vData, int nV, double* coeff);
	// Read the data from MPAS Ocean data may help us compute this
	double cellVolume(int cellId);
	bool isInCell(PointInfo& pInfo, const int cellId);
	CellType GetCellType(void) {return VORONOI;}

	// get min and maximal boundary
	void Boundary(VECTOR3& minB, VECTOR3& maxB);
	// set bounding box
	void SetBoundary(VECTOR3& minB, VECTOR3& maxB);
	// get grid spacing in x,y,z dimensions
	void GetGridSpacing(int cellId, double& xspace, double& yspace, double& zspace);
	// boundary intersection
	void BoundaryIntersection(VECTOR3& intersectP, VECTOR3& startP, VECTOR3& endP,double* stepSize, double oldStepSize);

	// Set MPASO grid data
	void setCellCoord(VECTOR3* inCellCoord);
	void setVertexCoord(VECTOR3* inVertexCoord);
	void setVertexIndexOnCell(int* inVerticesOnCell);
	void setCellIndexOnVertex(int* inCellsOnVertex);
	void setCellIndexOnCell(int* inCellsOnCell);
	void setNumVertexOnCell(int* inNumVerticesOnCell);
	void setMaxLevelCell(int* inMaxLevelCell);

	void setZTop(double* inZTop, int nTimesteps = 1);
	void setZTopTimestamps(const std::vector<double>& timestamps);
	const double* getZTopTimestep(int timestep) const;

	// set and get parameters
	void setNumCell(int numOfCell) { this->nCells = numOfCell; }
	void setNumTrueLocalCell(int numOfTrueLocalCells) { this->nTrueLocalCells = numOfTrueLocalCells; }
	void setNumLocalVertex(int numOfVertices) { this->nLocalVertices = numOfVertices; }
	void setNMaxEdges(int numOfMaxEdges) { this->nMaxEdges = numOfMaxEdges; }
	void setNTimestepsLoaded(int n) { this->nTimestepsLoaded = n; }
	void setNVertLevels(int numOfVertLevels) { this->nVertLevels = numOfVertLevels; }
	void setLatRange(double lat_min, double lat_max) { this->lat_min = lat_min; this->lat_max = lat_max; }
	void setLonRange(double lon_center, double lon_half_width) { this->lon_center = lon_center; this->lon_half_width = lon_half_width; }
	void setEarthRadius(double radius) { this->earth_radius = radius; }

	int getNumCell() { return this->nCells; }
	int getNumTrueLocalCell() { return this->nTrueLocalCells; }
	int getNumLocalVert() { return this->nLocalVertices; }
	int getNMaxEdges() { return this->nMaxEdges; }
	int getNTimestepsLoaded() { return this->nTimestepsLoaded; }
	int getNVertLevels() { return this->nVertLevels; }
	double getLatMin() { return this->lat_min; }
	double getLatMax() { return this->lat_max; }
	double getLonCenter() { return this->lon_center; }
	double getLonHalfWidth() { return this->lon_half_width; }
	double getEarthRadius() { return this->earth_radius; }

	// Memory accounting (analytical formula from dimension variables)
	// Returns bytes occupied by all static topology + coordinate arrays in MPASOGrid.
	size_t getGridMemBytes() const {
		size_t b = 0;
		b += (size_t)nCells         * sizeof(VECTOR3);          // cellCoord
		b += (size_t)nLocalVertices * sizeof(VECTOR3);          // vertexCoord
		b += (size_t)nCells         * nMaxEdges * sizeof(int);  // verticesOnCell
		b += (size_t)nLocalVertices * 3         * sizeof(int);  // cellsOnVertex
		b += (size_t)nCells         * nMaxEdges * sizeof(int);  // cellsOnCell
		b += (size_t)nCells         * sizeof(int);              // numVerticesOnCell
		b += (size_t)nCells         * sizeof(int);              // maxLevelCell
		b += (size_t)nTimestepsLoaded * nLocalVertices * nVertLevels * sizeof(double);  // zTop
		return b;
	}

	// Returns bytes for both horizontal + vertical Solution arrays.
	// nTimestepsLoaded=1 for streamline, =loadNTimeSteps for pathline.
	size_t getSolutionMemBytes() const {
		// pSolution + vSolution, each [nTimestepsLoaded][nCells * nVertLevels] VECTOR3
		return (size_t)2 * nTimestepsLoaded * nCells * nVertLevels * sizeof(VECTOR3);
	}

	// dump grid data in the block
	void dumpData(string &prefix);

protected:
	void Reset();
	int getZTopTimeWindow(double t, int cachedLowT, int& lowT, int& highT, double& ratio) const;
	// compute bounding box
	void ComputeBBox(void);
	// whether the point is in the bounding box
	bool isInBBox(VECTOR3& pos);
};

double triangle_area_by_length(double a, double b, double c);
double triangle_area(VECTOR3 &v1, VECTOR3 &v2, VECTOR3 &v3);
// double Wachspress_coeff(VECTOR3* vertices, VECTOR3 &p, int nVertices, int idx);

static void xyz2latlon(double& lat, double& lon, VECTOR3& xyz)
{
    lon = atan2(xyz[1],xyz[0]);
    lat = atan2(xyz[2], hypot(xyz[0], xyz[1]));
    lon = std::fmod(lon, 2*M_PI);
    if (lon < 0) lon += 2*M_PI;
}

static inline double wrapTo2Pi(double x)
{
    x = std::fmod(x, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x;
}

static inline double wrapToPi(double x)
{
    x = std::fmod(x + M_PI, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x - M_PI;
}

#endif
