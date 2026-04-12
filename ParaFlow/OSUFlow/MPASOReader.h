//////////////////////////////////////////////////////////////////////////
// File Reader for MPAS Ocean Data
//
// Created:	Yi-Tang Chen
// Date:	07/22/2023
//////////////////////////////////////////////////////////////////////////
#ifndef _MPASOREADER_H
#define _MPASOREADER_H

#include "MPASOheader.h"
#include "Field.h"
#include <unordered_set>
#include <algorithm>
#include <fstream>

static void convertENUVelocityToXYZ(const VECTOR3& xyzPoint, const double& Uzon, const double& Umer, const double& Uup, VECTOR3& xyzVel) 
{
    double Rxy, Rxyz, slon, clon, slat, clat;

    // Test for singularities at the poles
    if (xyzPoint[0] == 0.0 && xyzPoint[1] == 0.0) 
    {
        xyzVel[0] = 0.0;
        xyzVel[1] = 0.0;
        xyzVel[2] = Uup;  // Only the vertical component remains
        return;
    }

    // Compute geometric coordinate transform coefficients
    Rxy = std::sqrt(xyzPoint[0] * xyzPoint[0] + xyzPoint[1] * xyzPoint[1]);
    Rxyz = std::sqrt(xyzPoint[0] * xyzPoint[0] + xyzPoint[1] * xyzPoint[1] + xyzPoint[2] * xyzPoint[2]);
    slon = xyzPoint[1] / Rxy;
    clon = xyzPoint[0] / Rxy;
    slat = xyzPoint[2] / Rxyz;
    clat = Rxy / Rxyz;

    // Convert ENU to XYZ
    xyzVel[0] = -slon * Uzon - slat * clon * Umer + clon * clat * Uup;
    xyzVel[1] = clon * Uzon - slat * slon * Umer + slon * clat * Uup;
    xyzVel[2] = clat * Umer + slat * Uup;
}

class MPASOReader
{
private:
    // file info
    int ncDataid;
    int ncMeshid;
    int retval;
    int loadMesh;
    int areaId;
    // load grid info
    int nTotalCells;
    int nTotalVertices;
    int nTotalEdges;
    int verticesDegree;
    int nVertLevels;
    int nVertLevelsP1;
    int nMaxEdges;
    // compute grid info
    int nLocalCells;
    int nTrueLocalCells;
    int nLocalVertices;
    int nTimestepsInFile;  // total timesteps available in the NetCDF file

    int* areaIndices;
    std::vector<int> neighborAreaIndices;

    // Mapping
    int* globalCell2LocalCell;
    int* localCell2GlobalCell;
    int* globalVert2LocalVert;
    int* localVert2GlobalVert;
    
    // Time independent
    VECTOR3* cellCoord;
    VECTOR3* vertexCoord;
    int* verticesOnCell;
    int* cellsOnVertex;
    int* cellsOnCell;
    int* numVerticesOnCell;
    int* maxLevelCell;
    double* bottomDepth;
    int* verticesOnEdge;
    int* edgesOnCell;

    double lat_min, lat_max;
    double lon_center;      // in [0, 2pi)
    double lon_half_width;  // in [0, pi]

    // Vertex interpolation weights (cell-on-vertex barycentric weights)
    VECTOR3* cov_weights = nullptr;
    std::vector<int> requiredCellIndices;

    // Edge normal directions (computed once from mesh topology)
    VECTOR3* edgeNormal = nullptr;
    bool edgeNormalComputed = false;

    // Time varying
    VECTOR3** cellVelocity;    // normal velocity, [nTimestepsInMem][nCells * nVertLevels]
    double** cellVertVelocity; // vertical velocity, [nTimestepsInMem][nCells * nVertLevels]
    int nTimestepsInMem = 0;   // timesteps currently loaded for all time-varying arrays
    int timestepOffset  = 0;   // offset within the current file
    double* cellZTop;
    TimeVaryingDataConfig cfg;
    std::vector<double> m_fileTimestamps;  // real timestamps (seconds) across all files
    std::vector<int> m_timestepFileIdx;     // global timestep -> dataFiles index
    std::vector<int> m_timestepLocalIdx;    // global timestep -> timestep within that file
    std::vector<int> m_fileTimestepBase;    // dataFiles index -> first global timestep
    std::vector<int> m_fileTimestepCount;   // dataFiles index -> number of timesteps

    // Multi-file support
    std::vector<std::string> dataFiles;       // ordered list of all data files
    int currentFileIdx         = 0;           // index into dataFiles currently open
    int nTimestepsTotal        = 0;           // total timesteps across all files
    int currentFileTimestepBase = 0;          // global index of the first ts in current file

public:
    MPASOReader(const char* infile);
    MPASOReader(const char* infile, int inAreaId, int inNTotalCells, int* inAreaIndices, TimeVaryingDataConfig &incfg);
    MPASOReader(const char* meshfile, const char* datafile, int inAreaId, int inNTotalCells, int* inAreaIndices, TimeVaryingDataConfig &incfg);
    ~MPASOReader();

    void Reset();

    MPASOGrid* CreateMPASOGrid();
    // Time independent
    void readDimensions();
    void readCellCoord();
    void readVertexCoord();
    void readVertexOnCell();
    void readCellOnVertex();
    void readCellOnCell();
    void readNumVertexOnCell();
    void readMaxLevelCell();
    void readBottomDepth();
    void readVerticesOnEdge();
    void readEdgesOnCell();

    void readCellVelocity(int startTimestep, int N);
    void readZonMeridVelocity(int startTimestep, int N);
    void readNormalVelocity(int startTimestep, int N);
    void readVertVelocityTop(int startTimestep, int N);
    void readZTop(int timestep);
    void readLayerThickness(int timestep);
    void readAllTimestamps();  // parse xtime strings into m_fileTimestamps (seconds since t=0)

    // CreateSolution should be called after create MPASOGrid
    void InitSolutions(MPASOGrid* grid, Solution* &pSolution, Solution* &vSolution);
    void UpdateSolutions(MPASOGrid* grid, Solution* &pSolution, Solution* &vSolution);
    void GetNeighborIndices(std::vector<int> &neighborInd);
    void GetLocalCell2GlobalCell(int* &LC2GC, int &nLC);
    void GetGlobalCell2LocalCell(int* &GC2LC, int &nGC);
    void GetnVertLevels(int &nVLvl);
    int getTotalTimesteps() const { return this->nTimestepsTotal > 0 ? this->nTimestepsTotal : this->nTimestepsInFile; }
    int getTimestepOffset() const { return this->currentFileTimestepBase + this->timestepOffset; }

private:
    static int queryTimeDimFromNcid(int ncid);   // read "Time" dimension from an open NetCDF id
    void buildTimeIndex();                       // build global timestep -> (file, local timestep)
    void ensureDataFileOpen(int fileIdx);        // switch ncDataid to dataFiles[fileIdx] if needed
    void setGlobalTimestepOffset(int globalOffset);
    void loadGlobalTimestepIntoSlot(int globalTimestep, int outputSlot,
                                    int nVertNodes,
                                    VECTOR3** vertexVelocity,
                                    VECTOR3** vertexVertVelocity,
                                    double* vertexZTop);
    void appendTimestampsFromNcid(int ncid, double& t0, bool isFirst);  // parse xtime into m_fileTimestamps
    void buildRequiredCellSet();
    void computeCOVweights();
    void computeEdgeNormalDirection();
    void computeCellVelocity(double* normalVelocity, VECTOR3*& cellVelocity);
    void computeTimeVaryingVar(int localVertIdx,
                               VECTOR3* cellVelocity, double* cellVertVelocity,
                               double* cellZTop,
                               VECTOR3* pVelocityOut, VECTOR3* vVelocityOut,
                               double* ztopOut);
};

#endif 
