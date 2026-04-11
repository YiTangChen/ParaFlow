/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Vector Field: 3D Static or Time-Varying
//
///////////////////////////////////////////////////////////////////////////////


#ifndef _VECTOR_FIELD_H_
#define _VECTOR_FIELD_H_

#include "header.h"
#include "VectorMatrix.h"
#include "Grid.h"
#include "Solution.h"
#include "MPASOReader.h"
#include <typeinfo>

//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251 4100 4244)

//////////////////////////////////////////////////////////////////////////
// vector field class
//////////////////////////////////////////////////////////////////////////

class CVectorField
{
private:
	Grid* m_pGrid;						// grid
	Solution* m_pSolution;				// vector data
	Solution* m_vSolution;				// vertical vector data: only MPAS-Ocean
	int m_nTimeSteps;
	bool m_bIsNormalized;				// whether the solution is normalized or not
	double m_MinT, m_MaxT; // the min and max time step of the data field

public:
	// constructor and destructor
	CVectorField();
	CVectorField(Grid* pGrid, Solution* pSolution, int timesteps, double min_t=0.0);
	CVectorField(Grid* pGrid, Solution* pSolution, Solution* vSolution, int timesteps, double min_t=0.0);
	CVectorField(Grid* pGrid, int timesteps, double min_t=0.0);
	~CVectorField();

	int lerp_phys_coord(int cellId, CellTopoType eCellTopoType, double* coeff, VECTOR3& pos);
	int at_cell(int cellId, CellTopoType eCellTopoType, const double t, vector<VECTOR3>& vNodeData, int* cachedLowT = NULL);
	int at_vertcell(int cellId, CellTopoType eCellTopoType, const double t, vector<VECTOR3>& vNodeData, int* cachedLowT = NULL);
	int at_slice(int slice, SliceType eSliceType, const double t, vector<VECTOR3>&vSliceData);
	int at_vert(const int i, const int j, const int k, const double t, VECTOR3& dataValue);
	int at_phys(VECTOR3 pos, double t, VECTOR3& vecData, int* cachedLowT = NULL);
	int at_phys_truelocal(VECTOR3 pos, double t, VECTOR3& vecData);
	int at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const double t, VECTOR3& nodeData, int* cachedLowT = NULL);
	int at_phys(VECTOR3 pos, double t, VECTOR4& vecData, int* cachedLowT = NULL);
	int at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const double t, VECTOR4& nodeData, int* cachedLowT = NULL);
	int at_comp(const int i, const int j, const int k, const double t, VECTOR3& dataValue);
	double volume_of_cell(int cellId);
	void NormalizeField(bool bLocal);
	void ScaleField(double scale);
	void setHorizontalSolution(Solution* pSolution) { this->m_pSolution = pSolution; }
	void setVerticalSolution(Solution* vSolution) { this->m_vSolution = vSolution; }
	Grid* GetGrid() { return m_pGrid; }
	Solution* GetHorizontalSolution() { return m_pSolution; }
	Solution* GetVerticalSolution()   { return m_vSolution; }
	void SetTimeRange(double minT, double maxT) { m_MinT = minT; m_MaxT = maxT; }
	void dumpData(string &prefix);


	bool IsNormalized(void);
	void getDimension(int& xdim, int& ydim, int& zdim);
	CellType GetCellType(void) { return m_pGrid->GetCellType(); }
	int GetTimeSteps(void) {return m_nTimeSteps;}
	double GetMinTimeStep(void) {return m_MinT;}
	double GetMaxTimeStep(void) {return m_MaxT;}
	void GetInflowRegion(vector<VECTOR3>& inflowVerts, const double t);
	void GetOutflowRegion(vector<VECTOR3>& outflowVerts, const double t);
	void GetTangentialflowRegion(vector<VECTOR3>& tanflowVerts, const double t);
	void GetInflowSlice(vector<VECTOR3>& inflowVerts, const double t, const int slice, const SliceType eSliceType);
	void GetOutflowSlice(vector<VECTOR3>& outflowVerts, const double t, const int slice, const SliceType eSliceType);
	void GetTangentialflowSlice(vector<VECTOR3>& tanflowVerts, const double t, const int slice, const SliceType eSliceType);
	void Boundary(VECTOR3& minB, VECTOR3& maxB) { m_pGrid->Boundary(minB, maxB); };
	void SetBoundary(VECTOR3 minB, VECTOR3 maxB) {m_pGrid->SetBoundary(minB, maxB); 
	}
	void at_curl(int, VECTOR3&, VECTOR3&);
	void BoundaryIntersection(VECTOR3& intersectP,VECTOR3& startP,VECTOR3& endP, double* stepSize,double oldStepSize) 
	{ m_pGrid->BoundaryIntersection(intersectP, startP, endP, stepSize, oldStepSize); }
	void GenerateVortField(int t, bool bToNormalize, VECTOR3* pVort);
	void GenerateLapField(int t, bool bToNormalize, VECTOR3* pLap);
	// Resolve pInfo.inCell for the given position without computing velocity.
	// Use instead of at_phys when only cell mapping is needed.
	int resolve_cell(int fromCell, VECTOR3& pos, PointInfo& pInfo);
	

protected:
	// reset
	void Reset(void);
	// field functions
	bool isTimeVarying(void);
	// curl
	void curl(double, double, double, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&);
};

#endif
