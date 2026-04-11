
/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Grid: Irregular, Curvilinear, Cartesian grid
//
///////////////////////////////////////////////////////////////////////////////


#include "Grid.h"
#include <algorithm>
#include <cstring>

#pragma warning(disable : 4251 4100 4244 4101)

// Set to 1 to enable phys_to_cell debug output, 0 to disable
#define PHYS_TO_CELL_DEBUG 0

#if PHYS_TO_CELL_DEBUG
#define DBG(fmt, ...) do { fprintf(stderr, fmt, ##__VA_ARGS__); fflush(stderr); } while(0)
#else
#define DBG(fmt, ...) do {} while(0)
#endif


//////////////////////////////////////////////////////////////////////////
//
//	definition of Cartesian Grid Class
//
//////////////////////////////////////////////////////////////////////////

CartesianGrid::CartesianGrid(int xdim, int ydim, int zdim)
{
  m_nDimension[0] = xdim;  //the grid dimensions in C space 
  m_nDimension[1] = ydim;
  m_nDimension[2] = zdim;
}

CartesianGrid::CartesianGrid()
{
  Reset();
}

CartesianGrid::~CartesianGrid()
{
}

void CartesianGrid::Reset()
{
  m_nDimension[0] = m_nDimension[1] = m_nDimension[2] = 0;
  m_vMinBound.Zero();
  m_vMaxBound.Zero();
}


//////////////////////////////////////////////////////////////////////////
// set bounding box
//////////////////////////////////////////////////////////////////////////
void RegularCartesianGrid::SetBoundary(VECTOR3& minB, VECTOR3& maxB)
{
  m_vMinBound = minB;
  m_vMaxBound = maxB;
  mappingFactorX = (double)(xdim()-1)/(m_vMaxBound[0] - m_vMinBound[0]);
  mappingFactorY = (double)(ydim()-1)/(m_vMaxBound[1] - m_vMinBound[1]);
  mappingFactorZ = (double)(zdim()-1)/(m_vMaxBound[2] - m_vMinBound[2]);
  oneOvermappingFactorX = (m_vMaxBound[0] - m_vMinBound[0])/(double)(xdim()-1);
  oneOvermappingFactorY = (m_vMaxBound[1] - m_vMinBound[1])/(double)(ydim()-1);
  oneOvermappingFactorZ = (m_vMaxBound[2] - m_vMinBound[2])/(double)(zdim()-1);

  /*
  printf(" 1/p ****** %f %f %f *****\n", 
	 oneOvermappingFactorX,	 oneOvermappingFactorY, 
	 oneOvermappingFactorZ); 
  */
 
  // grid spacing
  gridSpacing = min(min(oneOvermappingFactorX, oneOvermappingFactorY), 
		    oneOvermappingFactorZ);
}


//////////////////////////////////////////////////////////////////////////
//
//
//	definition of Regular Cartesian Grid Class
//
//
//////////////////////////////////////////////////////////////////////////
// constructor and deconstructor
RegularCartesianGrid::RegularCartesianGrid():CartesianGrid()
{
  Reset();
}

RegularCartesianGrid::RegularCartesianGrid(int xdim, int ydim, int zdim):CartesianGrid(xdim, ydim, zdim)
{
  Reset();
  VECTOR3 a = VECTOR3(0,0,0);   // the default is from 0 to xdim-1, etc. 
  VECTOR3 b = VECTOR3(xdim-1, ydim-1, zdim-1); 
  SetBoundary(a, b); 
}

RegularCartesianGrid::~RegularCartesianGrid()
{
}

void RegularCartesianGrid::Reset(void)
{
  mappingFactorX = mappingFactorY = mappingFactorZ = 0.0;
  oneOvermappingFactorX = oneOvermappingFactorY = oneOvermappingFactorZ = 0.0;
  gridSpacing = 1.0;
}

void RegularCartesianGrid::Boundary(VECTOR3& minB, VECTOR3& maxB)
{
  minB = m_vMinBound;
  maxB = m_vMaxBound;
}

//////////////////////////////////////////////////////////////////////////
// whether the physical point pos is in bounding box
//////////////////////////////////////////////////////////////////////////
bool RegularCartesianGrid::isInBBox(VECTOR3& pos)
{



  if( (pos[0] >= m_vMinBound[0]) && (pos[0] <= m_vMaxBound[0]) &&
      (pos[1] >= m_vMinBound[1]) && (pos[1] <= m_vMaxBound[1]) &&
      (pos[2] >= m_vMinBound[2]) && (pos[2] <= m_vMaxBound[2]))
    return true;
  else {
    return false;
  }
}

// compute a default boundary 
void RegularCartesianGrid::ComputeBBox(void)
{
  VECTOR3 minB, maxB;
  
  minB.Set(0, 0, 0);  // default is from zero to xdim-1, etc. 
  maxB.Set((double)(xdim()-1), (double)(ydim()-1), (double)(zdim()-1));

  SetBoundary(minB, maxB);
}

//////////////////////////////////////////////////////////////////////////
// for Cartesian grid, this funcion means whether the physical point is 
// in the boundary
//////////////////////////////////////////////////////////////////////////
bool RegularCartesianGrid::at_phys(VECTOR3& pos)
{
	// whether in the bounding box
  if(!isInBBox(pos))
    return false;
  
  return true;
}

//////////////////////////////////////////////////////////////////////////
// get vertex list of a cell
// input
//		cellId:		cell Id
//		cellType:	cell type
// output
//		vVertices: the vertex lis of the cell
//////////////////////////////////////////////////////////////////////////
int RegularCartesianGrid::getCellVertices(int cellId, 
					  CellTopoType cellType, 
					  vector<int>& vVertices)
{
  int totalCell = xcelldim() * ycelldim() * zcelldim();
  int xidx, yidx, zidx, index;

  if((cellId < 0) || (cellId >= totalCell))
    //    return 0;
    return -1;

  vVertices.clear();
  zidx = cellId / (xcelldim() * ycelldim());
  yidx = cellId % (xcelldim() * ycelldim());
  yidx = yidx / xcelldim();
  xidx = cellId - zidx * xcelldim() * ycelldim() - yidx * xcelldim();

  for(int kFor = 0; kFor < 2; kFor++)
    for(int jFor = 0; jFor < 2; jFor++)
      for(int iFor = 0; iFor < 2; iFor++)
	{
	  index = (zidx+kFor) * ydim() * xdim() + (yidx + jFor) * xdim() + (xidx + iFor);
	  vVertices.push_back(index);
	}
  return 1;
}

//////////////////////////////////////////////////////////////////////////
// get the physical coordinate of the vertex
//
// input:
// verIdx: index of vertex
// output:
// pos: physical coordinate of vertex
//////////////////////////////////////////////////////////////////////////
bool RegularCartesianGrid::at_vertex(int verIdx, VECTOR3& pos)
{
  int xidx, yidx, zidx;
  int totalVer = xdim() * ydim() * zdim();
  if((verIdx < 0) || (verIdx >= totalVer))
    return false;

  zidx = verIdx / (xdim() * ydim());
  yidx = verIdx % (xdim() * ydim());
  yidx = verIdx / xdim();
  xidx = verIdx - zidx * xdim() * ydim() - yidx * xdim();

  double xpos = m_vMinBound[0] + xidx*oneOvermappingFactorX; 
  double ypos = m_vMinBound[1] + yidx*oneOvermappingFactorY; 
  double zpos = m_vMinBound[2] + zidx*oneOvermappingFactorZ; 

  // pos.Set((float)xidx, (float)yidx, (float)zidx);
  pos.Set((double)xpos, (double)ypos, (double)zpos);
  return true;
}

//////////////////////////////////////////////////////////////////////////
// whether the point in the physical position is in the cell
//
// input:
// phyCoord:	physical position
// interpolant:	interpolation coefficients
// output:
// pInfo.interpolant: interpolation coefficient
// return:		returns 1 if in cell
//////////////////////////////////////////////////////////////////////////
bool RegularCartesianGrid::isInCell(PointInfo& pInfo, const int cellId)
{
  if(!isInBBox(pInfo.phyCoord))
    return false;

  double cx, cy, cz; // computatnoial space x, y, and z 
  cx = (pInfo.phyCoord[0]-m_vMinBound[0])/oneOvermappingFactorX; 
  cy = (pInfo.phyCoord[1]-m_vMinBound[1])/oneOvermappingFactorY; 
  cz = (pInfo.phyCoord[2]-m_vMinBound[2])/oneOvermappingFactorZ; 

  int xidx, yidx, zidx;
  xidx = (int)floor(cx); 
  yidx = (int)floor(cy); 
  zidx = (int)floor(cz); 

  int inCell = zidx * ycelldim() * xcelldim() + yidx * xcelldim() + xidx;
  if(cellId == inCell)
    {
      pInfo.interpolant.Set(cx - (double)xidx, cy - (double)yidx, cz - (double)zidx);
      return true;
    }
  else
    return true;
}

//////////////////////////////////////////////////////////////////////////
// get the cell id and also interpolating coefficients for the given
// physical position
//
// input:
// phyCoord:	physical position
// interpolant:	interpolation coefficients
// fromCell:	if -1, initial point; otherwise this point is advected from others
// output:
// pInfo.inCell: in which cell
// pInfo.interpolant: interpolation coefficients
// return:	returns 1 if successful; otherwise returns -1
//////////////////////////////////////////////////////////////////////////
int RegularCartesianGrid::phys_to_cell(PointInfo& pInfo, double, int*)
{
  if(!isInBBox(pInfo.phyCoord))
    return -1;

  double cx, cy, cz; // computatnoial space x, y, and z 
  cx = (pInfo.phyCoord[0]-m_vMinBound[0])/oneOvermappingFactorX; 
  cy = (pInfo.phyCoord[1]-m_vMinBound[1])/oneOvermappingFactorY; 
  cz = (pInfo.phyCoord[2]-m_vMinBound[2])/oneOvermappingFactorZ; 

  int xidx, yidx, zidx;
  xidx = (int)floor(cx); 
  yidx = (int)floor(cy); 
  zidx = (int)floor(cz); 

  int inCell = zidx * ycelldim() * xcelldim() + yidx * xcelldim() + xidx;

  pInfo.inCell = inCell;
  pInfo.interpolant.Set(cx - (double)xidx, cy - (double)yidx, cz - (double)zidx);
  return 1;
}

//////////////////////////////////////////////////////////////////////////
// trilinear interpolation
// input
// nodeData:	8 corners of cube cell
// coeff:	bilinear interpolation coefficents
// output
// vData:	output
//////////////////////////////////////////////////////////////////////////
void RegularCartesianGrid::interpolate(VECTOR3& nodeData, 
				       vector<VECTOR3>& vData,
				       double* coeff)
{
  double fCoeff[3];
  fCoeff[0] = coeff[0];
  fCoeff[1] = coeff[1];
  fCoeff[2] = coeff[2];

  for(int iFor = 0; iFor < 3; iFor++)
    {
      if (vData.size() == 0) printf(" vData panic.\n"); 
	  if (vData.size() == 4) {
		  double bCoeff[2];
		  bCoeff[0] = coeff[0];
		  bCoeff[1] = coeff[1];
		  nodeData[iFor] = BiLerp(vData[0][iFor], vData[1][iFor], vData[2][iFor], vData[3][iFor],
			       bCoeff);
	  }
	  else {
		  nodeData[iFor] = TriLerp(vData[0][iFor], vData[1][iFor], vData[2][iFor], vData[3][iFor],
			       vData[4][iFor], vData[5][iFor], vData[6][iFor], vData[7][iFor],
			       fCoeff);
	  }
    }
	// if(vData.size() > 4) {
	// 	for(int iFor = 0; iFor < vData.size(); iFor++) {
	// 		std::cout << vData[iFor][0] << "," << vData[iFor][1] << "," << vData[iFor][2] << std::endl;
	// 	}
	// }
	// std::cout << "-------------------------------" << std::endl;
}

//////////////////////////////////////////////////////////////////////////
// get tetra volume
// input
// cellId:	which cell
// return the volume of this cell
//////////////////////////////////////////////////////////////////////////
double RegularCartesianGrid::cellVolume(int cellId)
{
  double volume = (double)oneOvermappingFactorX*oneOvermappingFactorY*oneOvermappingFactorZ;
  return volume;
}

void RegularCartesianGrid::BoundaryIntersection(VECTOR3& intersectP,
						VECTOR3& startP,
						VECTOR3& endP,
						double* stepSize, 
						double oldStepSize)
{
	VECTOR3 hitPoint;
	double thit;
	double xMin, xMax, yMin, yMax, zMin, zMax;

	intersectP.Set(-1, -1, -1);

	xMax = (double)(m_nDimension[0] - 1);
	yMax = (double)(m_nDimension[1] - 1);
	zMax = (double)(m_nDimension[2] - 1);
	xMin = yMin = zMin = 0;

	VECTOR3 planeN[6];
	planeN[0].Set(1, 0, 0);			//right
	planeN[1].Set(-1, 0, 0);		//left
	planeN[2].Set(0, 1, 0);			//top
	planeN[3].Set(0,-1, 0);			//bottom	
	planeN[4].Set(0, 0, 1);			//front
	planeN[5].Set(0, 0,-1);			//back

	VECTOR3 rayNormal;
	rayNormal.minus(endP, startP);
	rayNormal.Normalize();

	double vertex[8][3] = {	{xMin, yMin, zMin},		// rear, left, bottom
				{xMin, yMax, zMin},		// rear, left, top
				{xMax, yMin, zMin},		// rear, right, bottom
				{xMax, yMax, zMin},		// rear, right, top
				{xMin, yMin, zMax},		// front, left, bottom
				{xMin, yMax, zMax},		// front, left, top
				{xMax, yMin, zMax},		// front, right, bottom
				{xMax, yMax, zMax}};	// front, right, toop
	
	VECTOR3 A;
	A = startP;

	// a point on the plane
	VECTOR3 B[6];
	B[0].Set(xMax, (yMin + yMax)/2, (zMin + zMax)/2);			// right
	B[1].Set(xMin, (yMin + yMax)/2, (zMin + zMax)/2);			// left
	B[2].Set((xMax + xMin)/2, yMax, (zMin + zMax)/2);			// top
	B[3].Set((xMax + xMin)/2, yMin, (zMin + zMax)/2);			// bottom
	B[4].Set((xMax + xMin)/2, (yMin + yMax)/2, zMax);			// front
	B[5].Set((xMax + xMin)/2, (yMin + yMax)/2, zMin);			// rear

	VECTOR3 bMinusA[6];
	bMinusA[0].minus(B[0], A);
	bMinusA[1].minus(B[1], A);
	bMinusA[2].minus(B[2], A);
	bMinusA[3].minus(B[3], A);
	bMinusA[4].minus(B[4], A);
	bMinusA[5].minus(B[5], A);

	double nDotBMinusA[6];
	nDotBMinusA[0] = dot(planeN[0], bMinusA[0]);
	nDotBMinusA[1] = dot(planeN[1], bMinusA[1]);
	nDotBMinusA[2] = dot(planeN[2], bMinusA[2]);
	nDotBMinusA[3] = dot(planeN[3], bMinusA[3]);
	nDotBMinusA[4] = dot(planeN[4], bMinusA[4]);
	nDotBMinusA[5] = dot(planeN[5], bMinusA[5]);

	double nDotc[6];
	nDotc[0] = dot(rayNormal, planeN[0]);
	nDotc[1] = dot(rayNormal, planeN[1]);
	nDotc[2] = dot(rayNormal, planeN[2]);
	nDotc[3] = dot(rayNormal, planeN[3]);
	nDotc[4] = dot(rayNormal, planeN[4]);
	nDotc[5] = dot(rayNormal, planeN[5]);

	// intersect with the right plan?
	if(nDotc[0] != 0)							// ray does not parallel with the right plane
	{
		// get hit point
		thit = nDotBMinusA[0]/nDotc[0];
		if(thit >= 0)
		{
			hitPoint.Set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if((hitPoint[0] >= (xMax-EPS)) && (hitPoint[0] <= (xMax+EPS)) &&
				(hitPoint[1] >= (yMin-EPS))&& (hitPoint[1] <= (yMax+EPS))		&&
				(hitPoint[2] >= (zMin-EPS))&& (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[0] = xMax;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the left plan?
	if(nDotc[1] != 0)							// ray does not parallel with the left plane
	{
		// get hit point
		thit = nDotBMinusA[1]/nDotc[1];
		if(thit >= 0)
		{
			hitPoint.Set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if((hitPoint[0] >= (xMin-EPS)) && (hitPoint[0] <= (xMin+EPS)) &&
				(hitPoint[1] >= (yMin-EPS))&& (hitPoint[1] <= (yMax+EPS))		&&
				(hitPoint[2] >= (zMin-EPS))&& (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[0] = xMin;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the top plan?
	if(nDotc[2] != 0)							// ray does not parallel with the top plane
	{
		// get hit point
		thit = nDotBMinusA[2]/nDotc[2];
		if(thit >= 0)
		{
			hitPoint.Set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if((hitPoint[0] >= (xMin-EPS))	&& (hitPoint[0] <= (xMax+EPS))		&&
				(hitPoint[1] >= (yMax-EPS))	&& (hitPoint[1] <= (yMax+EPS)) &&
				(hitPoint[2] >= (zMin-EPS))	&& (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[1] = yMax;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the bottom plan?
	if(nDotc[3] != 0)							// ray does not parallel with the bottom plane
	{
		// get hit point
		thit = nDotBMinusA[3]/nDotc[3];
		if(thit >= 0)
		{
			hitPoint.Set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if( (hitPoint[0] >= (xMin-EPS))&& (hitPoint[0] <= (xMax+EPS)) &&
				(hitPoint[1] >= (yMin-EPS))&& (hitPoint[1] <= (yMin+EPS)) &&
				(hitPoint[2] >= (zMin-EPS))&& (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[1] = yMin;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the front plan?
	if(nDotc[4] != 0)							// ray does not parallel with the front plane
	{
		// get hit point
		thit = nDotBMinusA[4]/nDotc[4];
		if(thit >= 0)
		{
			hitPoint.Set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if( (hitPoint[0] >= (xMin-EPS)) && (hitPoint[0] <= (xMax+EPS)) &&
				(hitPoint[1] >= (yMin-EPS)) && (hitPoint[1] <= (yMax+EPS)) &&
				(hitPoint[2] >= (zMax-EPS)) && (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[2] = zMax;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the rear plan?
	if(nDotc[5] != 0)							// ray does not parallel with the rear plane
	{
		// get hit point
		thit = nDotBMinusA[5]/nDotc[5];
		if(thit >= 0)
		{
			hitPoint.Set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if( (hitPoint[0] >= (xMin-EPS)) && (hitPoint[0] <= (xMax+EPS)) &&
				(hitPoint[1] >= (yMin-EPS)) && (hitPoint[1] <= (yMax+EPS)) &&
				(hitPoint[2] >= (zMin-EPS)) && (hitPoint[2] <= (zMin+EPS)))
			{
				// in plane
				hitPoint[2] = zMin;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////
//																		//
//	definition of Irregular CartesianGrid Class							//
//																		//
//////////////////////////////////////////////////////////////////////////

/*
// constructor and deconstructor
IrregularCartesianGrid::IrregularCartesianGrid():CartesianGrid()
{
	Reset();
}

IrregularCartesianGrid::IrregularCartesianGrid(int xdim, int ydim, int zdim):CartesianGrid(xdim, ydim, zdim)
{
	Reset();
}

IrregularCartesianGrid::~IrregularCartesianGrid()
{
	delete[] m_pXSpacing;
	delete[] m_pYSpacing;
	delete[] m_pZSpacing;
	Reset();
}

void IrregularCartesianGrid::Reset(void)
{
	m_pXSpacing = m_pYSpacing = m_pZSpacing = NULL;
}
*/

//////////////////////////////////////////////////////////////////////////
//																		//
//	definition of Irregular Grid Class									//
//																		//
//////////////////////////////////////////////////////////////////////////
// constructor and deconstructor
IrregularGrid::IrregularGrid()
{
	Reset();
}

IrregularGrid::IrregularGrid(int nodeNum, int tetraNum, 
							 CVertex* pVertexGeom, 
							 CTetra* pTetra, 
							 TVertex* pVertexTopo)
{
	Reset();
	m_nNodeNum = nodeNum;
	m_nTetraNum = tetraNum;
	m_pVertexGeom = pVertexGeom;
	m_pTetra = pTetra;
	m_pVertexTopo = pVertexTopo;
}

IrregularGrid::~IrregularGrid()
{
	if(m_pVertexGeom != NULL)
		delete[] m_pVertexGeom;

	if(m_pTetra != NULL)
		delete[] m_pTetra;

	if(m_pTetraInfo != NULL)
		delete[] m_pTetraInfo;

	if(m_pVertexTopo != NULL)
		delete[] m_pVertexTopo;
}

void IrregularGrid::Reset(void)
{
	m_nNodeNum = m_nTetraNum = 0;
	m_bTetraInfoInit = false;
	m_pVertexGeom = NULL;
	m_pTetra = NULL;
	m_pTetraInfo = NULL;
	m_pVertexTopo = NULL;
}

void IrregularGrid::SetTetraInfoInit(bool bInit)
{
	m_bTetraInfoInit = bInit;
}

bool IrregularGrid::GetTetraInfoInit(void)
{
	return m_bTetraInfoInit;
}

//////////////////////////////////////////////////////////////////////////
// set bounding box
//////////////////////////////////////////////////////////////////////////
void IrregularGrid::SetBoundary(VECTOR3& minB, VECTOR3& maxB)
{
	m_vMinBound = minB;
	m_vMaxBound = maxB;
}

void IrregularGrid::Boundary(VECTOR3& minB, VECTOR3& maxB)
{
	minB = m_vMinBound;
	maxB = m_vMaxBound;
}

//////////////////////////////////////////////////////////////////////////
// whether the physical point pos is in bounding box
//////////////////////////////////////////////////////////////////////////
bool IrregularGrid::isInBBox(VECTOR3& pos)
{



	if( (pos[0] >= m_vMinBound[0]) && (pos[0] <= m_vMaxBound[0]) &&
		(pos[1] >= m_vMinBound[1]) && (pos[1] <= m_vMaxBound[1]) &&
		(pos[2] >= m_vMinBound[2]) && (pos[2] <= m_vMaxBound[2]))
		return true;
	else{
		return false;
	}
}

void IrregularGrid::ComputeBBox(void)
{
	VECTOR3 minB, maxB;

	minB.Set(FLT_MAX, FLT_MAX, FLT_MAX);
	maxB.Set(-FLT_MAX, -FLT_MAX, -FLT_MAX);

	for(int iFor = 0; iFor < m_nNodeNum; iFor++)
	{
		// x
		if(m_pVertexGeom[iFor].position[0] < minB[0])
			minB[0] = m_pVertexGeom[iFor].position[0];
		if(m_pVertexGeom[iFor].position[0] > maxB[0])
			maxB[0] = m_pVertexGeom[iFor].position[0];

		// y
		if(m_pVertexGeom[iFor].position[1] < minB[1])
			minB[1] = m_pVertexGeom[iFor].position[1];
		if(m_pVertexGeom[iFor].position[1] > maxB[1])
			maxB[1] = m_pVertexGeom[iFor].position[1];

		// z
		if(m_pVertexGeom[iFor].position[2] < minB[2])
			minB[2] = m_pVertexGeom[iFor].position[2];
		if(m_pVertexGeom[iFor].position[2] > maxB[2])
			maxB[2] = m_pVertexGeom[iFor].position[2];
	}

	SetBoundary(minB, maxB);
}

//////////////////////////////////////////////////////////////////////////
// whether the physical point is in the boundary
//////////////////////////////////////////////////////////////////////////
bool IrregularGrid::at_phys(VECTOR3& pos)
{
	PointInfo pInfo;

	// whether in the bounding box
	if(!isInBBox(pos))
		return false;

	pInfo.Set(pos, pInfo.interpolant, -1, -1);
	if(phys_to_cell(pInfo, 0.0, NULL) == 1)
		return true;

	return false;
}

//////////////////////////////////////////////////////////////////////////
// get vertex list of a cell
//////////////////////////////////////////////////////////////////////////
int IrregularGrid::getCellVertices(int cellId, 
								   CellTopoType cellType,
								   vector<int>& vVertices)
{
	if((cellId < 0) || (cellId >= m_nTetraNum))
		return 0;

	vVertices.clear();
	for(int iFor = 0; iFor < 4; iFor++)
	{
		vVertices.push_back(m_pTetra[cellId].ver[iFor]);
	}

	return 1;
}

//////////////////////////////////////////////////////////////////////////
// from physical coordinate to natural coordinate
//
// input:
// pCoord:	physical position
// output:
// nCoord:	interpolation coefficients
// return:	return false if determinant is 0
//////////////////////////////////////////////////////////////////////////
bool IrregularGrid::Physical2NaturalCoord(VECTOR3& nCoord, VECTOR3& pCoord, int cellId)
{
	VECTOR3 B, v, AmulB;
	int idx, i;
	double determinant, oneOverDeter;

	try
	{
		i = 0;
		idx = m_pTetra[cellId].ver[i];
		v.Set(m_pVertexGeom[idx].position[0], m_pVertexGeom[idx].position[1], m_pVertexGeom[idx].position[2]);	

		B.Set(pCoord[0] - v[0], pCoord[1] - v[1], pCoord[2] - v[2]);

		if(GetTetraInfoInit() == false)
		{
			SetTetraInfoInit(true);
			if(m_pTetraInfo == NULL)
				m_pTetraInfo = new TetraInfo[m_nTetraNum];

			for(i = 0; i < m_nTetraNum; i++)
			{
				PreGetP2NMatrix(m_pTetraInfo[i].p2nMatrix, i);
				m_pTetraInfo[i].volume = cellVolume(i);
			}
		}

		AmulB = m_pTetraInfo[cellId].p2nMatrix * B;

		determinant = m_pTetraInfo[cellId].volume * (double)6.0;
		if(determinant == 0) return false;

		oneOverDeter = (double)1.0 / determinant;
		nCoord[0] = AmulB[0] * oneOverDeter;
		nCoord[1] = AmulB[1] * oneOverDeter;
		nCoord[2] = AmulB[2] * oneOverDeter;
		return true;
	}
	catch(...)
	{
		printf("Exceptions happen in function IrregularGrid::Physical2NaturalCoord()!\n");
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////
// get the physical coordinate of the vertex
//
// input:
// verIdx: index of vertex
// output:
// pos: physical coordinate of vertex
//////////////////////////////////////////////////////////////////////////
bool IrregularGrid::at_vertex(int verIdx, VECTOR3& pos)
{
	if((verIdx < 0) || (verIdx >= m_nNodeNum))
		return false;

	pos.Set(m_pVertexGeom[verIdx].position[0],
			m_pVertexGeom[verIdx].position[1],
			m_pVertexGeom[verIdx].position[2]);
	return true;
}

//////////////////////////////////////////////////////////////////////////
// whether the point in the physical position is in the cell
//
// input:
// phyCoord:	physical position
// interpolant:	interpolation coefficients
// output:
// pInfo.interpolant: natural coordinate
// return:		returns 1 if in cell
//////////////////////////////////////////////////////////////////////////
bool IrregularGrid::isInCell(PointInfo& pInfo, const int cellId)
{
	VECTOR3 nCoord;

	if(!Physical2NaturalCoord(nCoord, pInfo.phyCoord, cellId))			// determinant is 0
	{
		return false;
	}

	pInfo.interpolant.Set(nCoord[0], nCoord[1], nCoord[2]);
	// four conditions to be in a tetra
	// all three natural coordinates >= 0, and the sum of them <= 1
	if( (nCoord[0] >= 0.0) && (nCoord[1] >= 0.0) && (nCoord[2] >= 0.0) && ((1.0 - nCoord[0] - nCoord[1] - nCoord[2]) >= 0))	
		return true;
	else															// not in tetra
		return false;
}

//////////////////////////////////////////////////////////////////////////
//	find the next tetra, from which the point goes out
//
//	input: 
//	tetraId:	the tetra in check
//	pInfo:		the physical and natural coordinates of point
//
//	output:		the next tetra to search
//////////////////////////////////////////////////////////////////////////
//changed by lijie, find the tetra with the highest absolute value
int IrregularGrid::nextTetra(PointInfo& pInfo, int tetraId)
{
	
	int nextT;
	VECTOR3 nCoord;
/*
	nCoord.Set(pInfo.interpolant[0], pInfo.interpolant[1], pInfo.interpolant[2]);

	if(nCoord[0] < 0.0)
		nextT = m_pTetra[tetraId].tetra[1];

	if(nCoord[1] < 0.0)
		nextT = m_pTetra[tetraId].tetra[2];

	if(nCoord[2] < 0.0)
		nextT = m_pTetra[tetraId].tetra[3];

	if((1.0 - nCoord[0] - nCoord[1] - nCoord[2]) < 0.0)
		nextT = m_pTetra[tetraId].tetra[0];

	//printf("nextT=%d\n",tetraId);
	return nextT;
	*/
	
	 nextT=-1;
//	VECTOR3 nCoord;

	nCoord.Set(pInfo.interpolant[0], pInfo.interpolant[1], pInfo.interpolant[2]);

	double tmp_ncord[4];
	tmp_ncord[0]=1.0 - nCoord[0] - nCoord[1] - nCoord[2];
	tmp_ncord[1]=nCoord[0];
	tmp_ncord[2]=nCoord[1];
	tmp_ncord[3]=nCoord[2];

	bool found=false;
	int next_id;
	double maxi=0;
	for(int i=0; i<4; i++)
	{
		if(tmp_ncord[i]<0)
		{
			if(fabs(tmp_ncord[i])>=maxi)
			{
				maxi=fabs(tmp_ncord[i]);
				next_id=i;
				found=true;
				break;
			}
		}
	}
	
	if(found==true)
	nextT = m_pTetra[tetraId].tetra[next_id];
	
	//printf("ncord=%f %f %f %f\n",tmp_ncord[0],tmp_ncord[1],tmp_ncord[2],tmp_ncord[3]);
	//printf("nextT=%d\n",next_id);
	return nextT;

}

//////////////////////////////////////////////////////////////////////////
// get the cell id and also interpolating coefficients for the given
// physical position
//
// input:
// phyCoord:	physical position
// interpolant:	interpolation coefficients
// fromCell:	if -1, initial point; otherwise this point is advected from others
// output:
// pInfo.inCell: in which cell
// pInfo.interpolant: natural coordinate
// return:	returns 1 if successful; otherwise returns -1
//////////////////////////////////////////////////////////////////////////
int IrregularGrid::phys_to_cell(PointInfo& pInfo, double, int*)
{
	int iFor, tetraInCheck, nextT;
	bool bFind;
	bool* pNonCheckedCells;

	if(pInfo.fromCell == -1)		// no previous information
	{
		bFind = false;
		for(iFor = 0; iFor < m_nTetraNum; iFor++)
		{
			if(isInCell(pInfo, iFor))			// in this cell
			{
				pInfo.inCell = iFor;
				bFind = true;
				return 1;
			}
			if((bFind == false) && (iFor == m_nTetraNum))
				return -1;
		}
	}
	else							// generated from previous point
	{
		bFind = false;
		tetraInCheck = pInfo.fromCell;

		// in order to prevent circular check
		pNonCheckedCells = new bool[m_nTetraNum];
		for(iFor = 0; iFor < m_nTetraNum; iFor++)
			pNonCheckedCells[iFor] = false;

                while(!bFind && tetraInCheck>=0)
		{
			if(isInCell(pInfo, tetraInCheck))		// whether in this tetra
			{
				//pInfo.inCell = iFor;//wrong
				pInfo.inCell = 	tetraInCheck;		//changed by lijie
				delete pNonCheckedCells;
				return 1;
			}
			else									// not in tetraInCheck
			{
				// find which face to go out
				pNonCheckedCells[tetraInCheck] = true;
				nextT = nextTetra(pInfo, tetraInCheck);
				tetraInCheck = nextT;
                                if(tetraInCheck>=0 && pNonCheckedCells[tetraInCheck] == true)	// circular check
				{
					delete pNonCheckedCells;
					return -1;
				}
			}
		}
	}

	return -1;
}

//////////////////////////////////////////////////////////////////////////
// barycentric interpolation
// input
// nodeData:	four corners of tetra
// coeff:		barycentric coordinate
// output
// vData:		output
//////////////////////////////////////////////////////////////////////////
void IrregularGrid::interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, double* coeff)
{
	double fTemp[4];
	double fCoeff[3];

	fCoeff[0] = coeff[0];
	fCoeff[1] = coeff[1];
	fCoeff[2] = coeff[2];

	for(int iFor = 0; iFor < 3; iFor++)
	{
		fTemp[0] = vData[0][iFor];
		fTemp[1] = vData[1][iFor];
		fTemp[2] = vData[2][iFor];
		fTemp[3] = vData[3][iFor];

		nodeData[iFor] = BaryInterp(fTemp, fCoeff);
	}
}

//////////////////////////////////////////////////////////////////////////
// get tetra volume
// input
// cellId:	which cell
// return the volume of this cell
//////////////////////////////////////////////////////////////////////////
double IrregularGrid::cellVolume(int cellId)
{
	double volume;
    VECTOR3 v[4];
	int idx, i;
	
	for(i = 0; i < 4; i++)
	{
		idx = m_pTetra[cellId].ver[i];
		v[i].Set(m_pVertexGeom[idx].position[0], m_pVertexGeom[idx].position[1], m_pVertexGeom[idx].position[2]);	
	}

	volume = (v[1][0] - v[0][0])*((v[2][1] - v[0][1])*(v[3][2] - v[0][2]) - (v[2][2] - v[0][2])*(v[3][1] - v[0][1])) + 
			 (v[2][0] - v[0][0])*((v[0][1] - v[1][1])*(v[3][2] - v[0][2]) - (v[0][2] - v[1][2])*(v[3][1] - v[0][1])) +
			 (v[3][0] - v[0][0])*((v[1][1] - v[0][1])*(v[2][2] - v[0][2]) - (v[1][2] - v[0][2])*(v[2][1] - v[0][1]));

	volume /= 6.0;
	return volume;
}

//////////////////////////////////////////////////////////////////////////
// used to compute the natural coordinate of points
//////////////////////////////////////////////////////////////////////////


//new version, modified by lijie
void IrregularGrid::PreGetP2NMatrix(MATRIX3& m, int cellId)
{
	VECTOR3 v[5];
	int idx, i;
	
	//printf("cell %d vertices=\n",cellId);

	for(i = 1; i < 5; i++)
	{
		idx = m_pTetra[cellId].ver[i-1];
		v[i].Set(m_pVertexGeom[idx].position[0], m_pVertexGeom[idx].position[1], m_pVertexGeom[idx].position[2]);	
	//	printf("%f %f %f \n",v[i].x(),v[i].y(),v[i].z());
	}

	/*
	m[0][0] = (v[4].z()-v[1].z()) * (v[3].y() - v[4].y()) - (v[3].z()-v[4].z()) * (v[4].y()-v[1].y());
	m[0][1] = (v[4].z()-v[1].z()) * (v[1].y() - v[2].y()) - (v[1].z()-v[2].z()) * (v[4].y()-v[1].y());
	m[0][2] = (v[2].z()-v[3].z()) * (v[1].y() - v[2].y()) - (v[1].z()-v[2].z()) * (v[2].y()-v[3].y());

	m[1][0] = (v[4].x()-v[1].x()) * (v[3].z() - v[4].z()) - (v[3].x()-v[4].x()) * (v[4].z()-v[1].z());
	m[1][1] = (v[4].x()-v[1].x()) * (v[1].z() - v[2].z()) - (v[1].x()-v[2].x()) * (v[4].z()-v[1].z());
	m[1][2] = (v[2].x()-v[3].x()) * (v[1].z() - v[2].z()) - (v[1].x()-v[2].x()) * (v[2].z()-v[3].z());

	m[2][0] = (v[4].y()-v[1].y()) * (v[3].x() - v[4].x()) - (v[3].y()-v[4].y()) * (v[4].x()-v[1].x());
	m[2][1] = (v[4].y()-v[1].y()) * (v[1].x() - v[2].x()) - (v[1].y()-v[2].y()) * (v[4].x()-v[1].x());
	m[2][2] = (v[2].y()-v[3].y()) * (v[1].x() - v[2].x()) - (v[1].y()-v[2].y()) * (v[2].x()-v[3].x());
	*/

	m[0][0] = (v[4].z()-v[1].z()) * (v[3].y() - v[4].y()) - (v[3].z()-v[4].z()) * (v[4].y()-v[1].y());
	m[1][0] = (v[4].z()-v[1].z()) * (v[1].y() - v[2].y()) - (v[1].z()-v[2].z()) * (v[4].y()-v[1].y());
	m[2][0] = (v[2].z()-v[3].z()) * (v[1].y() - v[2].y()) - (v[1].z()-v[2].z()) * (v[2].y()-v[3].y());

	m[0][1] = (v[4].x()-v[1].x()) * (v[3].z() - v[4].z()) - (v[3].x()-v[4].x()) * (v[4].z()-v[1].z());
	m[1][1] = (v[4].x()-v[1].x()) * (v[1].z() - v[2].z()) - (v[1].x()-v[2].x()) * (v[4].z()-v[1].z());
	m[2][1] = (v[2].x()-v[3].x()) * (v[1].z() - v[2].z()) - (v[1].x()-v[2].x()) * (v[2].z()-v[3].z());

	m[0][2] = (v[4].y()-v[1].y()) * (v[3].x() - v[4].x()) - (v[3].y()-v[4].y()) * (v[4].x()-v[1].x());
	m[1][2] = (v[4].y()-v[1].y()) * (v[1].x() - v[2].x()) - (v[1].y()-v[2].y()) * (v[4].x()-v[1].x());
	m[2][2] = (v[2].y()-v[3].y()) * (v[1].x() - v[2].x()) - (v[1].y()-v[2].y()) * (v[2].x()-v[3].x());

}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// with the newly intersection point, get the interpolated step size for this new point
///////////////////////////////////////////////////////////////////////////////////////////////////////
double getStepSize(VECTOR3& p, VECTOR3& p1, VECTOR3& p2, double oldStepSize)
{
	double newStepSize;
	double p1p, p1p2;

	p1p = (double)sqrt((p[0] - p1[0]) * (p[0] - p1[0])+
					  (p[1] - p1[1]) * (p[1] - p1[1])+
					  (p[2] - p1[2]) * (p[2] - p1[2]));

	p1p2 = (double)sqrt((p2[0] - p1[0]) * (p2[0] - p1[0])+
					   (p2[1] - p1[1]) * (p2[1] - p1[1])+
					   (p2[2] - p1[2]) * (p2[2] - p1[2]));

	assert(p1p2 != 0);

	newStepSize = oldStepSize * (p1p / p1p2);

	return newStepSize;
}


//////////////////////////////////////////////////////////////
//
//    MPAS Ocean Grid
//    By Yi-Tang Chen
//
//////////////////////////////////////////////////////////////

double triangle_area_by_length(double a, double b, double c)
{
    double s = (a+b+c) * 0.5;
    return sqrt(s*(s-a)*(s-b)*(s-c));
}

double triangle_area(VECTOR3 &v1, VECTOR3 &v2, VECTOR3 &v3)
{
    VECTOR3 u = v2 - v1;
    VECTOR3 w = v3 - v1;
    return cross(u, w).GetMag() * 0.5;
}

MPASOGrid::MPASOGrid()
{
	this->Reset();
}

MPASOGrid::~MPASOGrid()
{
	delete[] this->cellCoord;
	delete[] this->vertexCoord;
	delete[] this->verticesOnCell;
	delete[] this->cellsOnVertex;
	delete[] this->cellsOnCell;
	delete[] this->numVerticesOnCell;
	delete[] this->maxLevelCell;
	delete[] this->zTop;
}

bool MPASOGrid::at_vertex(int verIdx, VECTOR3& pos)
{
	if(verIdx < 0 || verIdx >= this->nLocalVertices)
		return false;
	pos = this->vertexCoord[verIdx];
	return true;
}

bool MPASOGrid::at_phys(VECTOR3& pos)
{
	assert(this->cellCoord != nullptr);
	if(pos.GetMag() > this->earth_radius)
		return false;
	if(!isInBBox(pos))
		return false;
	return true;
}

int MPASOGrid::getCellVertices(int nodeId, CellTopoType cellType, vector<int>& vVertices)
{
	// nodeId is a local flat index (localCellIdx * nVertLevels + vLevel)
	int cellidx = nodeId / this->nVertLevels;
	if (cellidx < 0 || cellidx >= this->nCells) return -1;
	int vLevel = nodeId % this->nVertLevels;
	int nVert = this->numVerticesOnCell[cellidx];
	int offset = cellidx * this->nMaxEdges;

	for(int i = 0; i < nVert; i++) {
		int lVertIdx = this->verticesOnCell[offset + i]; // 0-based local
		if (lVertIdx < 0 || lVertIdx >= this->nLocalVertices) return -1;
		vVertices.push_back(lVertIdx * this->nVertLevels + vLevel);
	}
	for(int i = 0; i < nVert; i++) {
		int lVertIdx = this->verticesOnCell[offset + i]; // 0-based local
		if (lVertIdx < 0 || lVertIdx >= this->nLocalVertices) return -1;
		vVertices.push_back(lVertIdx * this->nVertLevels + vLevel + 1);
	}
	return 1;
}

// Stack-array overload: writes upper-layer indices to vVertices[0..nVert-1]
// and lower-layer indices to vVertices[nVert..2*nVert-1]. Returns nVert or -1.
int MPASOGrid::getCellVertices(int nodeId, CellTopoType cellType, int* vVertices)
{
	if (cellType != T5_CELL) return -1;
	int cellidx = nodeId / this->nVertLevels;
	if (cellidx < 0 || cellidx >= this->nCells) return -1;
	int vLevel = nodeId % this->nVertLevels;
	int nVert = this->numVerticesOnCell[cellidx];
	int offset = cellidx * this->nMaxEdges;
	for(int i = 0; i < nVert; i++) {
		int lVertIdx = this->verticesOnCell[offset + i]; // 0-based local
		if (lVertIdx < 0 || lVertIdx >= this->nLocalVertices) return -1;
		vVertices[i]         = lVertIdx * this->nVertLevels + vLevel;
		vVertices[i + nVert] = lVertIdx * this->nVertLevels + vLevel + 1;
	}
	return nVert;
}

// Stack-array overload: vData[0..nV-1] upper layer, vData[nV..2nV-1] lower layer.
void MPASOGrid::interpolate(VECTOR3& nodeData, VECTOR3* vData, int nV, double* coeff)
{
	nodeData.Zero();
	VECTOR3 upperData, lowerData;
	for(int vertIdx = 0; vertIdx < nV; vertIdx++) {
		for(int i = 0; i < 3; i++) {
			upperData[i] += coeff[vertIdx] * vData[vertIdx][i];
			lowerData[i] += coeff[vertIdx] * vData[vertIdx + nV][i];
		}
	}
	for(int i = 0; i < 3; i++)
		nodeData[i] = coeff[nV] * upperData[i] + (1.0 - coeff[nV]) * lowerData[i];
}

// double Wachspress_coeff(VECTOR3* vertices, VECTOR3 &p, int nVertices, int idx)
// {
// 	double B = triangle_area(vertices[idx],
// 							vertices[(idx + 1) % nVertices],
// 							vertices[(idx + 2) % nVertices]);
// 	for(int i = 0; i < nVertices-2; i++) {
// 		double A = triangle_area(vertices[(idx+2+i) % nVertices],
// 								vertices[(idx+3+i) % nVertices],
// 								p);
// 		B = B*A;
// 	}
// 	return B;
// }

int MPASOGrid::phys_to_cell(PointInfo& pInfo, double t, int* cachedLowT)
{
	VECTOR3 phys = pInfo.phyCoord;
	VECTOR3 phys_ontop = phys;
	phys_ontop.Normalize();
	phys_ontop.scale(this->earth_radius);

	if(!isInBBox(phys))
		return -1;

	if(phys.GetMag() > this->earth_radius) {
		return -1;
	}

	// pInfo.fromCell is a local flat index (localCellIdx * nVertLevels + vLevel)
	int nearest_cellIdx = -1;
	int vLevel = -1;
	int fromCellIdx = (pInfo.fromCell >= 0) ? (pInfo.fromCell / this->nVertLevels) : -1;
	if(fromCellIdx >= 0 && fromCellIdx < this->nCells) {
		int curr_nVOC = this->numVerticesOnCell[fromCellIdx];
		int voc_offset = fromCellIdx * this->nMaxEdges;

		VECTOR3 d0 = phys_ontop - this->cellCoord[fromCellIdx];
		double min_dist2 = dot(d0, d0);
		nearest_cellIdx = fromCellIdx;
		DBG("[DBG-A] fast-path: fromCellIdx=%d nCells=%d curr_nVOC=%d\n", fromCellIdx, this->nCells, curr_nVOC);
		for(int i = 0; i < curr_nVOC; i++) {
			int neighbor_idx = this->cellsOnCell[voc_offset+i]; // 0-based local, -1 for none
			if(neighbor_idx >= 0) {
				VECTOR3 d = phys_ontop - this->cellCoord[neighbor_idx];
				double curr_dist2 = dot(d, d);
				if(curr_dist2 < min_dist2) {
					min_dist2 = curr_dist2;
					nearest_cellIdx = neighbor_idx;
				}
			}
		}
	}
	else {
		double min_dist2 = this->earth_radius * this->earth_radius;
		for(int lCellIdx = 0; lCellIdx < this->nCells; lCellIdx++) {
			VECTOR3 d = phys_ontop - this->cellCoord[lCellIdx];
			double curr_dist2 = dot(d, d);
			if(curr_dist2 < min_dist2) {
				min_dist2 = curr_dist2;
				nearest_cellIdx = lCellIdx;
			}
		}
	}

	DBG("[DBG-C] nearest_cellIdx=%d nCells=%d\n", nearest_cellIdx, this->nCells);
	if (nearest_cellIdx < 0 || nearest_cellIdx >= this->nCells)
		return -1;

	// check if the phys is in the cell
	int curr_nVertices = numVerticesOnCell[nearest_cellIdx];
	int nvoc_offset = nearest_cellIdx * nMaxEdges;
	bool dot_sign = true;
	constexpr int MAX_EDGES = MPASO_MAX_EDGES;
	VECTOR3 vertices[MAX_EDGES];
	double  omegas[MAX_EDGES];
	int     vertexidx[MAX_EDGES];
	DBG("[DBG-D] curr_nVertices=%d nearest_cellIdx=%d nvoc_offset=%d\n", curr_nVertices, nearest_cellIdx, nvoc_offset);
	for(int i = 0; i < curr_nVertices; i++) {
		vertexidx[i] = this->verticesOnCell[nvoc_offset+i]; // 0-based local, -1 for none
		if (vertexidx[i] < 0 || vertexidx[i] >= this->nLocalVertices) {
			DBG("[DBG-D] OOB! i=%d vertexidx=%d nLocalVertices=%d\n", i, vertexidx[i], this->nLocalVertices);
			return -1;
		}
		vertices[i] = this->vertexCoord[vertexidx[i]];
	}
	for(int i = 0; i < curr_nVertices; i++) {
		VECTOR3 normal = cross(vertices[i], vertices[(i+1)%curr_nVertices]);
		double dot_res = dot(normal, phys);
		if(i == 0)
			dot_sign = (dot_res > 0.0);
		else if(dot_sign != (dot_res > 0.0))
			return -1;
	}

	// check which level
	double depth = phys.GetMag() - this->cellCoord[nearest_cellIdx].GetMag();
	int cellOffset = nearest_cellIdx * this->nVertLevels;

	constexpr int MAX_VERT_LEVELS = 120;
	double currztop[MAX_VERT_LEVELS];
	double Aj[MAX_EDGES], Bk[MAX_EDGES];
	for(int j = 0; j < curr_nVertices; j++)
		Aj[j] = triangle_area(vertices[j], vertices[(j+1)%curr_nVertices], phys_ontop);
	for(int k = 0; k < curr_nVertices; k++)
		Bk[k] = triangle_area(vertices[k], vertices[(k+1)%curr_nVertices], vertices[(k+2)%curr_nVertices]);

	double omega_sum = 0.0;
	for(int k = 0; k < curr_nVertices; k++) {
		double w = Bk[k];
		for(int i = 2; i < curr_nVertices; i++)
			w *= Aj[(k+i) % curr_nVertices];
		omegas[(k+1)%curr_nVertices] = w;
		omega_sum += w;
	}
	if (omega_sum == 0.0)
		omega_sum = 1e-6;
	for(int i = 0; i < curr_nVertices; i++)
		omegas[i] = omegas[i] / omega_sum;

	DBG("[DBG-E] checking vertex locality and precomputing zTop offsets\n");
	// vertexidx[i] is already a local vertex index; use directly for zTop
	int lvoffsets[MAX_EDGES];
	for(int i = 0; i < curr_nVertices; i++) {
		if (vertexidx[i] < 0 || vertexidx[i] >= this->nLocalVertices) {
			DBG("[DBG-E] OOB! i=%d lv=%d nLocalVertices=%d\n", i, vertexidx[i], this->nLocalVertices);
			return -1;
		}
		lvoffsets[i] = vertexidx[i] * this->nVertLevels;
	}

	DBG("[DBG-F] computing zTop (nVertLevels=%d nLocalVertices=%d)\n", this->nVertLevels, this->nLocalVertices);
	int lowT, highT;
	double timeRatio;
	if(this->getZTopTimeWindow(t, cachedLowT ? *cachedLowT : -1, lowT, highT, timeRatio) == -1)
		return -1;
	if(cachedLowT)
		*cachedLowT = lowT;
	const double* zTopLow = this->getZTopTimestep(lowT);
	const double* zTopHigh = this->getZTopTimestep(highT);
	for(int vl = 0; vl < this->nVertLevels; vl++)
		currztop[vl] = 0.0;
	for(int i = 0; i < curr_nVertices; i++) {
		for(int vl = 0; vl < this->nVertLevels; vl++) {
			double z = zTopLow[lvoffsets[i] + vl];
			if(highT != lowT)
				z += (zTopHigh[lvoffsets[i] + vl] - z) * timeRatio;
			currztop[vl] += z * omegas[i];
		}
	}

	auto it = std::upper_bound(currztop, currztop + this->nVertLevels, depth, std::greater<double>());
	if (it != currztop + this->nVertLevels) {
		int vl = (int)(it - currztop);
		vLevel = (vl > 0) ? vl - 1 : 0;
	}

	if (vLevel == -1) {
		return -1;
	}

	// pInfo.inCell is a local flat index (localCellIdx * nVertLevels + vLevel)
	pInfo.inCell = cellOffset + vLevel;
	pInfo.nVertices = curr_nVertices;
	for(int i = 0; i < curr_nVertices; i++)
		pInfo.MPASO_interpolant[i] = omegas[i];
	double dz = currztop[vLevel+1] - currztop[vLevel];
	pInfo.MPASO_interpolant[curr_nVertices] = (dz != 0.0) ? (currztop[vLevel+1] - depth) / dz : 0.5;

	return 1;
}

int MPASOGrid::phys_to_truelocalcell(PointInfo& pInfo, double t, int* cachedLowT)
{
	VECTOR3 phys = pInfo.phyCoord;
	VECTOR3 phys_ontop = phys;
	phys_ontop.Normalize();
	phys_ontop.scale(this->earth_radius);

	if(!isInBBox(phys))
		return -1;

	if(phys.GetMag() > this->earth_radius) {
		return -1;
	}

	// pInfo.fromCell is a local flat index (localCellIdx * nVertLevels + vLevel)
	int nearest_cellIdx = -1;
	int vLevel = -1;
	int fromCellIdx = (pInfo.fromCell >= 0) ? (pInfo.fromCell / this->nVertLevels) : -1;
	if(fromCellIdx >= 0 && fromCellIdx < this->nCells) {
		int curr_nVOC = this->numVerticesOnCell[fromCellIdx];
		int voc_offset = fromCellIdx * this->nMaxEdges;

		VECTOR3 d0 = phys_ontop - this->cellCoord[fromCellIdx];
		double min_dist2 = dot(d0, d0);
		nearest_cellIdx = fromCellIdx;
		DBG("[DBG-A] fast-path: fromCellIdx=%d nCells=%d curr_nVOC=%d\n", fromCellIdx, this->nCells, curr_nVOC);
		for(int i = 0; i < curr_nVOC; i++) {
			int neighbor_idx = this->cellsOnCell[voc_offset+i]; // 0-based local, -1 for none
			if(neighbor_idx >= 0) {
				VECTOR3 d = phys_ontop - this->cellCoord[neighbor_idx];
				double curr_dist2 = dot(d, d);
				if(curr_dist2 < min_dist2) {
					min_dist2 = curr_dist2;
					nearest_cellIdx = neighbor_idx;
				}
			}
		}
	}
	else {
		// slow-path: search only true local cells (not ghost cells)
		double min_dist2 = this->earth_radius * this->earth_radius;
		for(int lCellIdx = 0; lCellIdx < this->nTrueLocalCells; lCellIdx++) {
			VECTOR3 d = phys_ontop - this->cellCoord[lCellIdx];
			double curr_dist2 = dot(d, d);
			if(curr_dist2 < min_dist2) {
				min_dist2 = curr_dist2;
				nearest_cellIdx = lCellIdx;
			}
		}
	}

	DBG("[DBG-C] nearest_cellIdx=%d nCells=%d\n", nearest_cellIdx, this->nCells);
	if (nearest_cellIdx < 0 || nearest_cellIdx >= this->nCells)
		return -1;

	// check if the phys is in the cell
	int curr_nVertices = numVerticesOnCell[nearest_cellIdx];
	int nvoc_offset = nearest_cellIdx * nMaxEdges;
	bool dot_sign = true;
	constexpr int MAX_EDGES = MPASO_MAX_EDGES;
	VECTOR3 vertices[MAX_EDGES];
	double  omegas[MAX_EDGES];
	int     vertexidx[MAX_EDGES];
	DBG("[DBG-D] curr_nVertices=%d nearest_cellIdx=%d nvoc_offset=%d\n", curr_nVertices, nearest_cellIdx, nvoc_offset);
	for(int i = 0; i < curr_nVertices; i++) {
		vertexidx[i] = this->verticesOnCell[nvoc_offset+i]; // 0-based local, -1 for none
		if (vertexidx[i] < 0 || vertexidx[i] >= this->nLocalVertices) {
			DBG("[DBG-D] OOB! i=%d vertexidx=%d nLocalVertices=%d\n", i, vertexidx[i], this->nLocalVertices);
			return -1;
		}
		vertices[i] = this->vertexCoord[vertexidx[i]];
	}
	for(int i = 0; i < curr_nVertices; i++) {
		VECTOR3 normal = cross(vertices[i], vertices[(i+1)%curr_nVertices]);
		double dot_res = dot(normal, phys);
		if(i == 0)
			dot_sign = (dot_res > 0.0);
		else if(dot_sign != (dot_res > 0.0))
			return -1;
	}

	double depth = phys.GetMag() - this->cellCoord[nearest_cellIdx].GetMag();
	int cellOffset = nearest_cellIdx * this->nVertLevels;

	constexpr int MAX_VERT_LEVELS = 120;
	double currztop[MAX_VERT_LEVELS];
	double Aj[MAX_EDGES], Bk[MAX_EDGES];
	for(int j = 0; j < curr_nVertices; j++)
		Aj[j] = triangle_area(vertices[j], vertices[(j+1)%curr_nVertices], phys_ontop);
	for(int k = 0; k < curr_nVertices; k++)
		Bk[k] = triangle_area(vertices[k], vertices[(k+1)%curr_nVertices], vertices[(k+2)%curr_nVertices]);

	double omega_sum = 0.0;
	for(int k = 0; k < curr_nVertices; k++) {
		double w = Bk[k];
		for(int i = 2; i < curr_nVertices; i++)
			w *= Aj[(k+i) % curr_nVertices];
		omegas[(k+1)%curr_nVertices] = w;
		omega_sum += w;
	}
	if (omega_sum == 0.0)
		omega_sum = 1e-6;
	for(int i = 0; i < curr_nVertices; i++)
		omegas[i] = omegas[i] / omega_sum;

	DBG("[DBG-E] checking vertex locality and precomputing zTop offsets\n");
	// vertexidx[i] is already a local vertex index; use directly for zTop
	int lvoffsets[MAX_EDGES];
	for(int i = 0; i < curr_nVertices; i++) {
		if (vertexidx[i] < 0 || vertexidx[i] >= this->nLocalVertices) {
			DBG("[DBG-E] OOB! i=%d lv=%d nLocalVertices=%d\n", i, vertexidx[i], this->nLocalVertices);
			return -1;
		}
		lvoffsets[i] = vertexidx[i] * this->nVertLevels;
	}

	DBG("[DBG-F] computing zTop (nVertLevels=%d nLocalVertices=%d)\n", this->nVertLevels, this->nLocalVertices);
	int lowT, highT;
	double timeRatio;
	if(this->getZTopTimeWindow(t, cachedLowT ? *cachedLowT : -1, lowT, highT, timeRatio) == -1)
		return -1;
	if(cachedLowT)
		*cachedLowT = lowT;
	const double* zTopLow = this->getZTopTimestep(lowT);
	const double* zTopHigh = this->getZTopTimestep(highT);
	for(int vl = 0; vl < this->nVertLevels; vl++)
		currztop[vl] = 0.0;
	for(int i = 0; i < curr_nVertices; i++) {
		for(int vl = 0; vl < this->nVertLevels; vl++) {
			double z = zTopLow[lvoffsets[i] + vl];
			if(highT != lowT)
				z += (zTopHigh[lvoffsets[i] + vl] - z) * timeRatio;
			currztop[vl] += z * omegas[i];
		}
	}

	auto it = std::upper_bound(currztop, currztop + this->nVertLevels, depth, std::greater<double>());
	if (it != currztop + this->nVertLevels) {
		int vl = (int)(it - currztop);
		vLevel = (vl > 0) ? vl - 1 : 0;
	}

	if (vLevel == -1) {
		phys.Normalize();
		phys.scale(this->earth_radius + currztop[this->nVertLevels-1] + 0.1);
		vLevel = this->nVertLevels-2;
		pInfo.phyCoord = phys;
	}

	// pInfo.inCell is a local flat index (localCellIdx * nVertLevels + vLevel)
	pInfo.inCell = cellOffset + vLevel;
	pInfo.nVertices = curr_nVertices;
	for(int i = 0; i < curr_nVertices; i++)
		pInfo.MPASO_interpolant[i] = omegas[i];
	double dz = currztop[vLevel+1] - currztop[vLevel];
	pInfo.MPASO_interpolant[curr_nVertices] = (dz != 0.0) ? (currztop[vLevel+1] - depth) / dz : 0.5;

	return 1;
}

///////////////////////////////////////////////////////////////////////////////////
//
//    vData contains 2 layer's mesh
//    coeff last element is the vertical coefficient
//    other coeff elements are for planar Wachspress interpolation
//
///////////////////////////////////////////////////////////////////////////////////
void MPASOGrid::interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, double* coeff)
{
	nodeData.Zero();
	VECTOR3 upperData, lowerData;
	int nV = vData.size()/2;
	for(int i = 0; i < 3; i++) {
		for(int vertIdx = 0; vertIdx < nV; vertIdx++) {
			upperData[i] += coeff[vertIdx]*vData[vertIdx][i];
			lowerData[i] += coeff[vertIdx]*vData[vertIdx+nV][i];
		}
	}
	for(int i = 0; i < 3; i++) {
		nodeData[i] = coeff[nV]*upperData[i] + (1.0-coeff[nV])*lowerData[i];
	}
}

double MPASOGrid::cellVolume(int cellId)
{
	// cellId is a local flat index: localCellIdx * nVertLevels + vLevel.
	int localCellIdx = cellId / this->nVertLevels;
	int localVLevel = cellId % this->nVertLevels;
	if (localCellIdx < 0 || localCellIdx >= this->nCells)
		return 2.7e13;  // fallback: (30,000 m)^3
	int nVert = this->numVerticesOnCell[localCellIdx];
	if (nVert <= 0) return 2.7e13;
	int offset = localCellIdx * this->nMaxEdges;
	VECTOR3 center = this->cellCoord[localCellIdx];
	double height = 1.0;
	if (this->zTop && localVLevel < this->nVertLevels-1) {
		height = this->zTop[localCellIdx * this->nVertLevels + localVLevel] - this->zTop[localCellIdx * this->nVertLevels + localVLevel + 1];
	}
	double area = 0.0;
	for (int i = 0; i < nVert; i++) {
		int vidx = this->verticesOnCell[offset + i]; // 0-based local
		int vidx_next = this->verticesOnCell[offset + (i+1)%nVert];
		area += triangle_area(center, this->vertexCoord[vidx], this->vertexCoord[vidx_next]);
	}
	return area * height;
}

bool MPASOGrid::isInCell(PointInfo& pInfo, const int cellId)
{
	int res = this->phys_to_cell(pInfo, 0.0, NULL);
	if(res == -1) return false;
	return pInfo.inCell == cellId;
}

void MPASOGrid::Boundary(VECTOR3& minB, VECTOR3& maxB)
{
	minB = m_MinBound;
	maxB = m_MaxBound;
}

void MPASOGrid::SetBoundary(VECTOR3& minB, VECTOR3& maxB)
{
	m_MinBound = minB;
	m_MaxBound = maxB;
}

void MPASOGrid::GetGridSpacing(int cellId, double& xspace, double& yspace, double& zspace)
{
	// MPAS-Ocean is irregular mesh. This function is useless
}

void MPASOGrid::BoundaryIntersection(VECTOR3& intersectP, VECTOR3& startP, VECTOR3& endP,double* stepSize, double oldStepSize)
{

}


void MPASOGrid::setCellCoord(VECTOR3* inCellCoord)
{
	assert(this->nCells != 0);
	this->cellCoord = new VECTOR3[this->nCells];
	for(int lc = 0; lc < this->nCells; lc++) {
		this->cellCoord[lc] = inCellCoord[lc];
		this->earth_radius = max(this->earth_radius, this->cellCoord[lc].GetMag());
	}
}

void MPASOGrid::setVertexCoord(VECTOR3* inVertexCoord)
{
	assert(this->nLocalVertices != 0);
	this->vertexCoord = new VECTOR3[this->nLocalVertices];
	memcpy(this->vertexCoord, inVertexCoord,
	       this->nLocalVertices * sizeof(VECTOR3));
}

void MPASOGrid::setVertexIndexOnCell(int* inVerticesOnCell)
{
	assert(this->nCells != 0);
	assert(this->nMaxEdges != 0);
	this->verticesOnCell = new int[this->nCells * this->nMaxEdges];
	memcpy(this->verticesOnCell, inVerticesOnCell,
	       this->nCells * this->nMaxEdges * sizeof(int));
}

void MPASOGrid::setCellIndexOnVertex(int* inCellsOnVertex)
{
	assert(this->nLocalVertices != 0);
	this->cellsOnVertex = new int[this->nLocalVertices * 3];
	memcpy(this->cellsOnVertex, inCellsOnVertex,
	       this->nLocalVertices * 3 * sizeof(int));
}

void MPASOGrid::setCellIndexOnCell(int* inCellsOnCell)
{
	assert(this->nCells != 0);
	assert(this->nMaxEdges != 0);
	this->cellsOnCell = new int[this->nCells * this->nMaxEdges];
	memcpy(this->cellsOnCell, inCellsOnCell,
	       this->nCells * this->nMaxEdges * sizeof(int));
}

void MPASOGrid::setNumVertexOnCell(int* inNumVerticesOnCell)
{
	assert(this->nCells != 0);
	this->numVerticesOnCell = new int[this->nCells];
	memcpy(this->numVerticesOnCell, inNumVerticesOnCell,
	       this->nCells * sizeof(int));
}

void MPASOGrid::setMaxLevelCell(int* inMaxLevelCell)
{
	assert(this->nCells != 0);
	this->maxLevelCell = new int[this->nCells];
	memcpy(this->maxLevelCell, inMaxLevelCell,
	       this->nCells * sizeof(int));
}

void MPASOGrid::setZTop(double *inZTop, int nTimesteps)
{
	assert(this->nLocalVertices != 0);
	assert(this->nVertLevels != 0);
	assert(nTimesteps > 0);
	size_t nValues = static_cast<size_t>(nTimesteps) * this->nLocalVertices * this->nVertLevels;
	delete[] this->zTop;
	this->zTop = new double[nValues];
	for(size_t i = 0; i < nValues; i++)
		this->zTop[i] = inZTop[i];
}

void MPASOGrid::setZTopTimestamps(const std::vector<double>& timestamps)
{
	assert(timestamps.empty() || static_cast<int>(timestamps.size()) == this->nTimestepsLoaded);
	this->zTopTimestamps = timestamps;
}

const double* MPASOGrid::getZTopTimestep(int timestep) const
{
	assert(this->zTop != nullptr);
	assert(timestep >= 0 && timestep < this->nTimestepsLoaded);
	size_t timestepSize = static_cast<size_t>(this->nLocalVertices) * this->nVertLevels;
	return this->zTop + timestep * timestepSize;
}

int MPASOGrid::getZTopTimeWindow(double t, int cachedLowT,
								 int& lowT, int& highT, double& ratio) const
{
	lowT = 0;
	highT = 0;
	ratio = 0.0;
	if(this->nTimestepsLoaded <= 0)
		return -1;

	if(static_cast<int>(this->zTopTimestamps.size()) == this->nTimestepsLoaded) {
		if(t < this->zTopTimestamps.front() || t > this->zTopTimestamps.back())
			return -1;
		if(this->nTimestepsLoaded == 1)
			return 1;
		if(t == this->zTopTimestamps.back()) {
			lowT = this->nTimestepsLoaded - 2;
			highT = this->nTimestepsLoaded - 1;
			ratio = 1.0;
			return 1;
		}

		if(cachedLowT >= 0 && cachedLowT < this->nTimestepsLoaded - 1 &&
		   t >= this->zTopTimestamps[cachedLowT]) {
			lowT = cachedLowT;
			while(lowT + 1 < this->nTimestepsLoaded - 1 &&
				  t >= this->zTopTimestamps[lowT + 1])
				lowT++;
		}
		else {
			lowT = static_cast<int>(std::upper_bound(this->zTopTimestamps.begin(),
													 this->zTopTimestamps.end(), t)
									- this->zTopTimestamps.begin()) - 1;
			if(lowT < 0)
				lowT = 0;
		}
		highT = lowT + 1;
		double span = this->zTopTimestamps[highT] - this->zTopTimestamps[lowT];
		ratio = (span != 0.0) ? (t - this->zTopTimestamps[lowT]) / span : 0.0;
		return 1;
	}

	double adjusted_t = t;
	if(adjusted_t < 0.0 || adjusted_t > static_cast<double>(this->nTimestepsLoaded - 1))
		return -1;
	if(this->nTimestepsLoaded == 1)
		return 1;
	if(adjusted_t == static_cast<double>(this->nTimestepsLoaded - 1)) {
		lowT = this->nTimestepsLoaded - 2;
		highT = this->nTimestepsLoaded - 1;
		ratio = 1.0;
		return 1;
	}
	lowT = static_cast<int>(floor(adjusted_t));
	highT = lowT + 1;
	ratio = adjusted_t - static_cast<double>(lowT);
	return 1;
}

void MPASOGrid::dumpData(string &prefix){
	// cell coordinates
	string filename = prefix + "_CellCoord.bin";
	ofstream cellCoordOutFile(filename, ios::binary);
	if (!cellCoordOutFile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	// All arrays are now local-indexed; write directly.
	cellCoordOutFile.write(reinterpret_cast<const char*>(&this->nCells), sizeof(int));
	for(int i = 0; i < this->nCells; i++) {
		cellCoordOutFile.write(reinterpret_cast<const char*>(&this->cellCoord[i][0]), sizeof(double));
		cellCoordOutFile.write(reinterpret_cast<const char*>(&this->cellCoord[i][1]), sizeof(double));
		cellCoordOutFile.write(reinterpret_cast<const char*>(&this->cellCoord[i][2]), sizeof(double));
	}
	cellCoordOutFile.close();

	// vertex coordinates
	filename = prefix + "_VertexCoord.bin";
	ofstream vertexCoordOutFile(filename, ios::binary);
	if (!vertexCoordOutFile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	vertexCoordOutFile.write(reinterpret_cast<const char*>(&this->nLocalVertices), sizeof(int));
	for(int i = 0; i < this->nLocalVertices; i++) {
		vertexCoordOutFile.write(reinterpret_cast<const char*>(&this->vertexCoord[i][0]), sizeof(double));
		vertexCoordOutFile.write(reinterpret_cast<const char*>(&this->vertexCoord[i][1]), sizeof(double));
		vertexCoordOutFile.write(reinterpret_cast<const char*>(&this->vertexCoord[i][2]), sizeof(double));
	}
	vertexCoordOutFile.close();

	// vertex indices on cell (0-based local)
	filename = prefix + "_VerticesOnCell.bin";
	ofstream verticesOnCellOutFile(filename, ios::binary);
	if (!verticesOnCellOutFile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	verticesOnCellOutFile.write(reinterpret_cast<const char*>(&this->nCells), sizeof(int));
	verticesOnCellOutFile.write(reinterpret_cast<const char*>(&this->nMaxEdges), sizeof(int));
	for(int i = 0; i < this->nCells * this->nMaxEdges; i++)
		verticesOnCellOutFile.write(reinterpret_cast<const char*>(&this->verticesOnCell[i]), sizeof(int));
	verticesOnCellOutFile.close();

	// cell indices on vertex (0-based local)
	filename = prefix + "_CellsOnVertex.bin";
	ofstream cellsOnVertexOutFile(filename, ios::binary);
	if (!cellsOnVertexOutFile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	cellsOnVertexOutFile.write(reinterpret_cast<const char*>(&this->nLocalVertices), sizeof(int));
	for(int i = 0; i < this->nLocalVertices * 3; i++)
		cellsOnVertexOutFile.write(reinterpret_cast<const char*>(&this->cellsOnVertex[i]), sizeof(int));
	cellsOnVertexOutFile.close();

	// cell indices on cell (0-based local)
	filename = prefix + "_CellsOnCell.bin";
	ofstream cellsOnCellOutFile(filename, ios::binary);
	if (!cellsOnCellOutFile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	cellsOnCellOutFile.write(reinterpret_cast<const char*>(&this->nCells), sizeof(int));
	cellsOnCellOutFile.write(reinterpret_cast<const char*>(&this->nMaxEdges), sizeof(int));
	for(int i = 0; i < this->nCells * this->nMaxEdges; i++)
		cellsOnCellOutFile.write(reinterpret_cast<const char*>(&this->cellsOnCell[i]), sizeof(int));
	cellsOnCellOutFile.close();

	// num vertices on cell
	filename = prefix + "_NumVerticesOnCell.bin";
	ofstream numVerticesOnCellOutFile(filename, ios::binary);
	if (!numVerticesOnCellOutFile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	numVerticesOnCellOutFile.write(reinterpret_cast<const char*>(&this->nCells), sizeof(int));
	for(int i = 0; i < this->nCells; i++)
		numVerticesOnCellOutFile.write(reinterpret_cast<const char*>(&this->numVerticesOnCell[i]), sizeof(int));
	numVerticesOnCellOutFile.close();

	// max level on cell
	filename = prefix + "_MaxLevelCell.bin";
	ofstream maxLevelCellOutFile(filename, ios::binary);
	if (!maxLevelCellOutFile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	maxLevelCellOutFile.write(reinterpret_cast<const char*>(&this->nCells), sizeof(int));
	for(int i = 0; i < this->nCells; i++)
		maxLevelCellOutFile.write(reinterpret_cast<const char*>(&this->maxLevelCell[i]), sizeof(int));
	maxLevelCellOutFile.close();

	// zTop
	filename = prefix + "_ZTop.bin";
	ofstream zTopOutFile(filename, ios::binary);
	if (!zTopOutFile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	zTopOutFile.write(reinterpret_cast<const char*>(&this->nLocalVertices), sizeof(int));
	zTopOutFile.write(reinterpret_cast<const char*>(&this->nTimestepsLoaded), sizeof(int));
	zTopOutFile.write(reinterpret_cast<const char*>(&this->nVertLevels), sizeof(int));
	size_t nTimeVaryingValues = static_cast<size_t>(this->nTimestepsLoaded) * this->nLocalVertices * this->nVertLevels;
	for(size_t i = 0; i < nTimeVaryingValues; i++) {
		zTopOutFile.write(reinterpret_cast<const char*>(&this->zTop[i]), sizeof(double));
	}
	zTopOutFile.close();

}

void MPASOGrid::Reset()
{
	// Time independent
	this->cellCoord = NULL;
    this->vertexCoord = NULL;
    this->verticesOnCell = NULL;
	this->cellsOnVertex = NULL;
	this->cellsOnCell = NULL;
    this->numVerticesOnCell = NULL;
	this->maxLevelCell = NULL;

	// Time-varying
	this->zTop = NULL;
	this->zTopTimestamps.clear();

	// Parameters
    this->nCells = 0;
	this->nTrueLocalCells = 0;
    this->nLocalVertices = 0;
    this->nMaxEdges = 0;
	this->nTimestepsLoaded = 0;
	this->nVertLevels = 0;
	this->m_MinBound.Zero();
	this->m_MaxBound.Zero();
	this->earth_radius = 0.0;
	this->lat_min = 0.0;
	this->lat_max = 0.0;
    this->lon_center = 0.0;      // in [0, 2pi)
    this->lon_half_width = 0.0;  // in [0, pi]

}

void MPASOGrid::ComputeBBox(void)
{
	return;
}

bool MPASOGrid::isInBBox(VECTOR3& pos)
{
	double lat, lon;
    xyz2latlon(lat, lon, pos);
	if (lat < this->lat_min || lat > this->lat_max)
        return false;

    double dlon = wrapToPi(lon - this->lon_center);
    if (std::abs(dlon) > this->lon_half_width)
        return false;

	return true;
}
