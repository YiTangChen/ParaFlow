#include "MPASOReader.h"
#include <cmath>
#include <cstring>

MPASOReader::MPASOReader(const char* infile)
{
    if ((retval = nc_open(infile, NC_NOWRITE, &ncDataid))) {
        fprintf(stderr, "Error: Cannot open %s (%s)\n", infile, nc_strerror(retval));
        throw std::runtime_error("MPASOReader: failed to open NetCDF file");
    }
    this->Reset();
}

MPASOReader::MPASOReader(const char* infile, int inAreaId, int inNTotalCells, int* inAreaIndices, TimeVaryingDataConfig &incfg)
{
    if ((retval = nc_open(infile, NC_NOWRITE, &ncDataid))) {
        fprintf(stderr, "Error: Cannot open %s (%s)\n", infile, nc_strerror(retval));
        throw std::runtime_error("MPASOReader: failed to open NetCDF file");
    }
    this->Reset();
    this->areaId = inAreaId;
    this->nTotalCells = inNTotalCells;
    this->areaIndices = inAreaIndices;
    this->cfg = incfg;
    this->dataFiles = incfg.dataFiles;
    this->buildTimeIndex();
}

MPASOReader::MPASOReader(const char* meshfile, const char* datafile, int inAreaId, int inNTotalCells, int* inAreaIndices, TimeVaryingDataConfig &incfg)
{
    if ((retval = nc_open(meshfile, NC_NOWRITE, &ncMeshid))) {
        fprintf(stderr, "Error: Cannot open %s (%s)\n", meshfile, nc_strerror(retval));
        throw std::runtime_error("MPASOReader: failed to open NetCDF file");
    }
    if ((retval = nc_open(datafile, NC_NOWRITE, &ncDataid))) {
        fprintf(stderr, "Error: Cannot open %s (%s)\n", datafile, nc_strerror(retval));
        throw std::runtime_error("MPASOReader: failed to open NetCDF file");
    }

    this->Reset();
    this->loadMesh = 1; // indicate that mesh is loaded
    this->areaId = inAreaId;
    this->nTotalCells = inNTotalCells;
    this->areaIndices = inAreaIndices;
    this->cfg = incfg;
    this->dataFiles = incfg.dataFiles;
    this->buildTimeIndex();
}

MPASOReader::~MPASOReader()
{
    if(this->loadMesh) {
        if ((retval = nc_close(ncMeshid)))
            fprintf(stderr, "Error: Cannot close NetCDF file (ncMeshid=%d): %s\n", ncMeshid, nc_strerror(retval));
    }
    if ((retval = nc_close(ncDataid)))
        fprintf(stderr, "Error: Cannot close NetCDF file (ncDataid=%d): %s\n", ncDataid, nc_strerror(retval));

    delete[] this->cellCoord;
    delete[] this->vertexCoord;
    delete[] this->globalCell2LocalCell;
    delete[] this->localCell2GlobalCell;
    delete[] this->globalVert2LocalVert;
    delete[] this->localVert2GlobalVert;
    delete[] this->verticesOnCell;
    delete[] this->cellsOnVertex;
    delete[] this->cellsOnCell;
    delete[] this->numVerticesOnCell;
    delete[] this->maxLevelCell;
    delete[] this->bottomDepth;
    delete[] this->edgesOnCell;
    delete[] this->verticesOnEdge;
    delete[] this->edgeNormal;
    delete[] this->cov_weights;
    if (this->cellVelocity != nullptr) {
        for (int t = 0; t < this->nTimestepsInMem; ++t)
            delete[] this->cellVelocity[t];
        delete[] this->cellVelocity;
    }
    if (this->cellVertVelocity != nullptr) {
        for (int t = 0; t < this->nTimestepsInMem; ++t)
            delete[] this->cellVertVelocity[t];
        delete[] this->cellVertVelocity;
    }
    delete[] this->cellZTop;
}

void MPASOReader::Reset()
{
    this->areaId = -1;
    this->nTotalCells = 0;
    this->nTotalVertices = 0;
    this->verticesDegree = 0;
    this->loadMesh = 0;
    this->areaIndices = nullptr;

    // Mapping
    this->globalCell2LocalCell = nullptr;
    this->localCell2GlobalCell = nullptr;
    this->globalVert2LocalVert = nullptr;
    this->localVert2GlobalVert = nullptr;
    
    // Time independent
    this->cellCoord = nullptr;
    this->vertexCoord = nullptr;
    this->verticesOnCell = nullptr;
    this->cellsOnVertex = nullptr;
    this->cellsOnVertex = nullptr;
    this->cellsOnCell = nullptr;
    this->numVerticesOnCell = nullptr;
    this->maxLevelCell = nullptr;
    this->bottomDepth = nullptr;
    this->verticesOnEdge = nullptr;
    this->edgesOnCell = nullptr;
    this->lat_min = 100.0;
    this->lat_max = -100.0;
    this->lon_center = 0.0;
    this->lon_half_width = 0.0;

    // Time varying
    this->cellVelocity    = nullptr;
    this->cellVertVelocity = nullptr;
    this->nTimestepsInMem  = 0;
    this->timestepOffset   = 0;
    this->cellZTop = nullptr;
    this->requiredCellIndices.clear();
    this->m_fileTimestamps.clear();
    this->m_timestepFileIdx.clear();
    this->m_timestepLocalIdx.clear();
    this->m_fileTimestepBase.clear();
    this->m_fileTimestepCount.clear();
    this->currentFileIdx          = 0;
    this->nTimestepsTotal         = 0;
    this->currentFileTimestepBase = 0;

    this->nLocalCells = 0;
    this->nTrueLocalCells = 0;
    this->nLocalVertices = 0;
    this->nTotalEdges = 0;
    this->nTimestepsInFile = 0;
    this->nVertLevels = 0;
    this->nMaxEdges = 0;
    this->nVertLevelsP1 = 0;
}

// ---------------------------------------------------------------------------
// Multi-file helpers
// ---------------------------------------------------------------------------

int MPASOReader::queryTimeDimFromNcid(int ncid)
{
    int dimid;
    size_t len = 0;
    if (nc_inq_dimid(ncid, "Time", &dimid) == NC_NOERR)
        nc_inq_dimlen(ncid, dimid, &len);
    return static_cast<int>(len);
}

void MPASOReader::buildTimeIndex()
{
    this->m_timestepFileIdx.clear();
    this->m_timestepLocalIdx.clear();
    this->m_fileTimestepBase.clear();
    this->m_fileTimestepCount.clear();

    int nFiles = this->dataFiles.empty() ? 1 : static_cast<int>(this->dataFiles.size());
    this->m_fileTimestepBase.assign(nFiles, -1);
    this->m_fileTimestepCount.assign(nFiles, 0);

    int total = 0;
    for (int fileIdx = 0; fileIdx < nFiles; ++fileIdx) {
        int count = 0;
        if (fileIdx == this->currentFileIdx) {
            count = queryTimeDimFromNcid(this->ncDataid);
        } else if (!this->dataFiles.empty()) {
            int tmpid;
            if (nc_open(this->dataFiles[fileIdx].c_str(), NC_NOWRITE, &tmpid) == NC_NOERR) {
                count = queryTimeDimFromNcid(tmpid);
                nc_close(tmpid);
            } else {
                std::cerr << "[MPASOReader::buildTimeIndex] Skipping unreadable file "
                          << fileIdx << " (" << this->dataFiles[fileIdx] << ").\n";
            }
        }

        this->m_fileTimestepBase[fileIdx] = total;
        this->m_fileTimestepCount[fileIdx] = count;
        for (int localTimestep = 0; localTimestep < count; ++localTimestep) {
            this->m_timestepFileIdx.push_back(fileIdx);
            this->m_timestepLocalIdx.push_back(localTimestep);
        }
        total += count;
    }

    this->nTimestepsTotal = total;
    if (this->currentFileIdx >= 0 && this->currentFileIdx < nFiles) {
        this->currentFileTimestepBase = this->m_fileTimestepBase[this->currentFileIdx];
        this->nTimestepsInFile = this->m_fileTimestepCount[this->currentFileIdx];
    }
}

void MPASOReader::ensureDataFileOpen(int fileIdx)
{
    if (fileIdx == this->currentFileIdx)
        return;

    if (this->dataFiles.empty()) {
        if (fileIdx == 0)
            return;
        throw std::runtime_error("MPASOReader: requested a non-existent data file");
    }
    if (fileIdx < 0 || fileIdx >= static_cast<int>(this->dataFiles.size()))
        throw std::runtime_error("MPASOReader: data file index out of range");

    nc_close(this->ncDataid);
    const char* nextFile = this->dataFiles[fileIdx].c_str();
    if ((this->retval = nc_open(nextFile, NC_NOWRITE, &this->ncDataid))) {
        fprintf(stderr, "Error: Cannot open %s (%s)\n", nextFile, nc_strerror(this->retval));
        throw std::runtime_error("MPASOReader: failed to open NetCDF file");
    }

    this->currentFileIdx = fileIdx;
    if (fileIdx >= static_cast<int>(this->m_fileTimestepCount.size()))
        this->buildTimeIndex();
    this->nTimestepsInFile = this->m_fileTimestepCount[fileIdx];
    this->currentFileTimestepBase = this->m_fileTimestepBase[fileIdx];
    this->timestepOffset = 0;
    this->nTimestepsInMem = 0;
}

void MPASOReader::setGlobalTimestepOffset(int globalOffset)
{
    if (this->m_timestepFileIdx.empty())
        this->buildTimeIndex();

    if (globalOffset <= 0) {
        int fileIdx = this->m_timestepFileIdx.empty() ? 0 : this->m_timestepFileIdx.front();
        this->ensureDataFileOpen(fileIdx);
        this->timestepOffset = 0;
        return;
    }

    if (globalOffset >= this->nTimestepsTotal) {
        if (!this->m_timestepFileIdx.empty()) {
            int lastGlobal = this->nTimestepsTotal - 1;
            int fileIdx = this->m_timestepFileIdx[lastGlobal];
            this->ensureDataFileOpen(fileIdx);
            this->timestepOffset = this->m_timestepLocalIdx[lastGlobal] + 1;
        }
        return;
    }

    int fileIdx = this->m_timestepFileIdx[globalOffset];
    this->ensureDataFileOpen(fileIdx);
    this->timestepOffset = this->m_timestepLocalIdx[globalOffset];
}

void MPASOReader::loadGlobalTimestepIntoSlot(int globalTimestep, int outputSlot,
                                             int nVertNodes,
                                             VECTOR3** vertexVelocity,
                                             VECTOR3** vertexVertVelocity,
                                             double* vertexZTop)
{
    if (this->m_timestepFileIdx.empty())
        this->buildTimeIndex();
    assert(globalTimestep >= 0 && globalTimestep < this->nTimestepsTotal);

    int fileIdx = this->m_timestepFileIdx[globalTimestep];
    int localTimestep = this->m_timestepLocalIdx[globalTimestep];
    this->ensureDataFileOpen(fileIdx);

    if (cfg.zTop)
        this->readZTop(localTimestep);
    else
        this->readLayerThickness(localTimestep);

    if (cfg.velocity.cartesian.X)
        this->readCellVelocity(localTimestep, 1);
    else if (cfg.velocity.spherical.zonal)
        this->readZonMeridVelocity(localTimestep, 1);
    else if (cfg.velocity.normal)
        this->readNormalVelocity(localTimestep, 1);
    else
        std::cerr << "[MPASOReader] Error: cfg.velocity is not properly set." << std::endl;
    this->readVertVelocityTop(localTimestep, 1);

    size_t gridBase = static_cast<size_t>(outputSlot) * nVertNodes;
    for (int localVertIdx = 0; localVertIdx < this->nLocalVertices; localVertIdx++) {
        int vertexOffset = localVertIdx * this->nVertLevels;
        this->computeTimeVaryingVar(localVertIdx,
                    this->cellVelocity[0], this->cellVertVelocity[0], this->cellZTop,
                    vertexVelocity[outputSlot] + vertexOffset,
                    vertexVertVelocity[outputSlot] + vertexOffset,
                    vertexZTop + gridBase + vertexOffset);
    }
}

void MPASOReader::appendTimestampsFromNcid(int ncid, double& t0, bool isFirst)
{
    if (!cfg.xtime) return;
    int varid;
    if (nc_inq_varid(ncid, cfg.xtime.value().c_str(), &varid) != NC_NOERR) return;

    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    if (nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dimids, nullptr) != NC_NOERR) return;
    if (ndims != 2) return;

    size_t nTime, strLen;
    nc_inq_dimlen(ncid, dimids[0], &nTime);
    nc_inq_dimlen(ncid, dimids[1], &strLen);

    std::vector<char> buf(nTime * strLen, '\0');
    size_t start[2] = {0, 0};
    size_t count[2] = {nTime, strLen};
    if (nc_get_vara_text(ncid, varid, start, count, buf.data()) != NC_NOERR) return;

    // Convert calendar date to Julian Day Number (proleptic Gregorian calendar).
    // Using integer arithmetic so there is no floating-point approximation error
    // from variable month/year lengths (avoids the ~30.4375 days/month hack).
    auto toJulianDay = [](int Y, int Mo, int D) -> long long {
        long long a = (14 - Mo) / 12;
        long long y = Y + 4800 - a;
        long long m = Mo + 12 * a - 3;
        return D + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045LL;
    };

    auto parseSeconds = [&](int idx) -> double {
        const char* s = buf.data() + idx * strLen;
        int Y=0, Mo=0, D=0, H=0, Mi=0, S=0;
        sscanf(s, "%d-%d-%d_%d:%d:%d", &Y, &Mo, &D, &H, &Mi, &S);
        long long jd = toJulianDay(Y, Mo, D);
        return (double)jd * 86400.0 + H * 3600.0 + Mi * 60.0 + (double)S;
    };

    if (isFirst)
        t0 = parseSeconds(0);

    for (size_t i = 0; i < nTime; i++)
        m_fileTimestamps.push_back(parseSeconds(static_cast<int>(i)) - t0);
}

// ---------------------------------------------------------------------------

MPASOGrid* MPASOReader::CreateMPASOGrid()
{
    MPASOGrid* grid;
    grid = new MPASOGrid();
    this->readDimensions();
    this->readCellCoord();
    this->readVertexCoord();
    this->readVertexOnCell();
    this->readCellOnVertex();
    this->readCellOnCell();
    this->readNumVertexOnCell();
    this->readMaxLevelCell();
    this->readBottomDepth();
    this->readVerticesOnEdge();
    this->readEdgesOnCell();

    if(this->areaId >= 0) {
        // create global cell index to local cell index mapping
        this->nLocalCells = 0;
        this->globalCell2LocalCell = new int[this->nTotalCells];
        for(int i = 0; i < this->nTotalCells; i++) {
            this->globalCell2LocalCell[i] = 0;
            if(this->areaIndices[i] == this->areaId) {
                this->globalCell2LocalCell[i] = this->nLocalCells + 1;
                this->nLocalCells++;
            }
        }

        this->nTrueLocalCells = this->nLocalCells;

        // add ghost cells to mapping
        std::unordered_set<int> neighbor_idx_set;
        for(int gCellIdx = 0; gCellIdx < this->nTotalCells; gCellIdx++) {
            if(this->areaIndices[gCellIdx] == this->areaId) {
                int nNeighbors = this->numVerticesOnCell[gCellIdx];
                int global_offset = gCellIdx * this->nMaxEdges;
                for(int i = 0; i < nNeighbors; i++) {
                    int neighGCellIdx = this->cellsOnCell[global_offset + i]-1;
                    if(neighGCellIdx >= 0 && this->globalCell2LocalCell[neighGCellIdx] == 0) {
                        this->globalCell2LocalCell[neighGCellIdx] = this->nLocalCells + 1;
                        this->nLocalCells++;
                        int neighbor_areaidx = this->areaIndices[neighGCellIdx];
                        neighbor_idx_set.insert(neighbor_areaidx);
                    }
                }
            }
        }

        // store neighbor block indices
        for (auto iter = neighbor_idx_set.begin(); iter != neighbor_idx_set.end(); ++iter) {
            this->neighborAreaIndices.push_back((*iter));
        }
        std::sort(this->neighborAreaIndices.begin(), this->neighborAreaIndices.end());

        // create local cell index mapping
        this->localCell2GlobalCell = new int[this->nLocalCells];
        for(int i = 0; i < this->nTotalCells; i++) {
            if(this->globalCell2LocalCell[i] > 0)
                this->localCell2GlobalCell[this->globalCell2LocalCell[i]-1] = i+1;
        }

        // create global vertex index to local vertex index mapping
        this->nLocalVertices = 0;
        this->globalVert2LocalVert = new int[this->nTotalVertices];
        for(int gCellIdx = 0; gCellIdx < this->nTotalVertices; gCellIdx++)
            this->globalVert2LocalVert[gCellIdx] = 0;
        for(int localCellIdx = 0; localCellIdx < this->nLocalCells; localCellIdx++) {
            int gCellIdx = this->localCell2GlobalCell[localCellIdx]-1;
            int nNeighbors = this->numVerticesOnCell[gCellIdx];
            int voc_offset = gCellIdx * this->nMaxEdges;
            for(int j = 0; j < nNeighbors; j++) {
                int neighVertIdx = this->verticesOnCell[voc_offset + j]-1;
                if(this->globalVert2LocalVert[neighVertIdx] == 0) {
                    this->globalVert2LocalVert[neighVertIdx] = this->nLocalVertices + 1;
                    this->nLocalVertices++;
                }
            }
        }

        // create local vertex index to global vertex index mapping
        this->localVert2GlobalVert = new int[this->nLocalVertices];
        for(int gVertIdx = 0; gVertIdx < this->nTotalVertices; gVertIdx++) {
            if(this->globalVert2LocalVert[gVertIdx] > 0)
                this->localVert2GlobalVert[this->globalVert2LocalVert[gVertIdx]-1] = gVertIdx+1;
        }

        // create bounding box for the area
        // 1. Transform vertex coordinates to lat/lon and find lat range and lon range
        double lat, lon;
        double pole_eps = 1e-2;
        std::vector<double> lons;
        lons.reserve(this->nLocalVertices);
        
        for (int i = 0; i < this->nLocalVertices; i++) {
            xyz2latlon(lat, lon, vertexCoord[this->localVert2GlobalVert[i] - 1]);
            lons.push_back(wrapTo2Pi(lon));
            if (lat < this->lat_min) this->lat_min = lat;
            if (lat > this->lat_max) this->lat_max = lat;
        }
        // 2. Clamp latitude to poles if very close
        if (this->lat_max > M_PI * 0.5 - pole_eps) this->lat_max = M_PI * 0.5;
        if (this->lat_min < -M_PI * 0.5 + pole_eps) this->lat_min = -M_PI * 0.5;

        std::sort(lons.begin(), lons.end());

        // 3. Find largest cyclic gap in longitude to determine the minimal covering arc
        double max_gap = -1.0;
        double curr_gap;
        int max_gap_idx = -1;

        for (int i = 0; i < this->nLocalVertices - 1; ++i) {
            curr_gap = lons[i + 1] - lons[i];
            if (curr_gap > max_gap) {
                max_gap = curr_gap;
                max_gap_idx = i;
            }
        }

        curr_gap = (lons[0] + 2.0 * M_PI) - lons[this->nLocalVertices - 1];
        if (curr_gap > max_gap) {
            max_gap = curr_gap;
            max_gap_idx = this->nLocalVertices - 1;
        }

        // 4. Minimal covering arc is everything except the largest gap
        double start = lons[(max_gap_idx + 1) % this->nLocalVertices];
        double end   = lons[max_gap_idx];

        if (end < start) end += 2.0 * M_PI;

        double width = end - start;
        double center = 0.5 * (start + end);

        this->lon_center = wrapTo2Pi(center);
        this->lon_half_width = 0.5 * width;
    }
    else {
        this->nLocalCells = this->nTotalCells;
        this->nTrueLocalCells = this->nTotalCells;
        this->nLocalVertices = this->nTotalVertices;
        this->lat_min = -M_PI_2;
        this->lat_max =  M_PI_2;
        this->lon_center     = M_PI;
        this->lon_half_width = M_PI;

        this->globalCell2LocalCell = new int[this->nTotalCells];
        this->localCell2GlobalCell = new int[this->nTotalCells];
        for(int i = 0; i < this->nLocalCells; i++) {
            this->globalCell2LocalCell[i] = i+1;
            this->localCell2GlobalCell[i] = i+1;
        }
        this->globalVert2LocalVert = new int[this->nTotalVertices];
        this->localVert2GlobalVert = new int[this->nTotalVertices];
        for(int i = 0; i < this->nTotalVertices; i++) {
            this->globalVert2LocalVert[i] = i+1;
            this->localVert2GlobalVert[i] = i+1;
        }
    }

    // Set grid parameters
    grid->setNumCell(this->nLocalCells);
    grid->setNumTrueLocalCell(this->nTrueLocalCells);
    grid->setNumLocalVertex(this->nLocalVertices);
    grid->setNMaxEdges(this->nMaxEdges);
    grid->setNVertLevels(this->nVertLevels);
    grid->setLatRange(this->lat_min, this->lat_max);
    grid->setLonRange(this->lon_center, this->lon_half_width);

    // Set grid data
// Build local cellCoord (size nLocalCells)
    VECTOR3* localCellCoord = new VECTOR3[this->nLocalCells];
    for(int lc = 0; lc < this->nLocalCells; lc++) {
        int gc = this->localCell2GlobalCell[lc] - 1;
        localCellCoord[lc] = this->cellCoord[gc];
    }
    grid->setCellCoord(localCellCoord);
    delete[] localCellCoord;
    // cellCoord (global) kept alive: computeCOVweights, computeTimeVaryingVar, computeEdgeNormalDirection need it

    // Build local vertexCoord (size nLocalVertices)
    VECTOR3* localVertexCoord = new VECTOR3[this->nLocalVertices];
    for(int lv = 0; lv < this->nLocalVertices; lv++) {
        int gv = this->localVert2GlobalVert[lv] - 1;
        localVertexCoord[lv] = this->vertexCoord[gv];
    }
    grid->setVertexCoord(localVertexCoord);
    delete[] localVertexCoord;
    // vertexCoord (global) kept alive: computeCOVweights and computeEdgeNormalDirection need it
    // Build local verticesOnCell (0-based local vertex indices, size nLocalCells × nMaxEdges)
    int* localVerticesOnCell = new int[this->nLocalCells * this->nMaxEdges];
    for(int lc = 0; lc < this->nLocalCells; lc++) {
        int gc = this->localCell2GlobalCell[lc] - 1;
        for(int j = 0; j < this->nMaxEdges; j++) {
            int gv = this->verticesOnCell[gc * this->nMaxEdges + j] - 1;
            int lv = (gv >= 0) ? this->globalVert2LocalVert[gv] - 1 : -1;
            localVerticesOnCell[lc * this->nMaxEdges + j] = lv;
        }
    }
    delete[] this->verticesOnCell;
    this->verticesOnCell = nullptr;
    grid->setVertexIndexOnCell(localVerticesOnCell);
    delete[] localVerticesOnCell;
    // Build local cellsOnVertex (0-based local cell indices, size nLocalVertices × verticesDegree)
    int* localCellsOnVertex = new int[this->nLocalVertices * this->verticesDegree];
    for(int lv = 0; lv < this->nLocalVertices; lv++) {
        int gv = this->localVert2GlobalVert[lv] - 1;
        for(int j = 0; j < this->verticesDegree; j++) {
            int gc = this->cellsOnVertex[gv * this->verticesDegree + j] - 1;
            int lc = (gc >= 0) ? this->globalCell2LocalCell[gc] - 1 : -1;
            localCellsOnVertex[lv * this->verticesDegree + j] = lc;
        }
    }
    grid->setCellIndexOnVertex(localCellsOnVertex);
    delete[] localCellsOnVertex;
    // cellsOnVertex (global) kept alive: computeCOVweights and computeTimeVaryingVar need it
    // Build local cellsOnCell (0-based local neighbor cell indices, size nLocalCells × nMaxEdges)
    int* localCellsOnCell = new int[this->nLocalCells * this->nMaxEdges];
    for(int lc = 0; lc < this->nLocalCells; lc++) {
        int gc = this->localCell2GlobalCell[lc] - 1;
        for(int j = 0; j < this->nMaxEdges; j++) {
            int gnbr = this->cellsOnCell[gc * this->nMaxEdges + j] - 1;
            int lnbr = (gnbr >= 0) ? this->globalCell2LocalCell[gnbr] - 1 : -1;
            localCellsOnCell[lc * this->nMaxEdges + j] = lnbr;
        }
    }
    delete[] this->cellsOnCell;
    this->cellsOnCell = nullptr;
    grid->setCellIndexOnCell(localCellsOnCell);
    delete[] localCellsOnCell;
    // Build local numVerticesOnCell (size nLocalCells)
    int* localNumVerticesOnCell = new int[this->nLocalCells];
    for(int lc = 0; lc < this->nLocalCells; lc++) {
        int gc = this->localCell2GlobalCell[lc] - 1;
        localNumVerticesOnCell[lc] = this->numVerticesOnCell[gc];
    }
    grid->setNumVertexOnCell(localNumVerticesOnCell);
    delete[] localNumVerticesOnCell;
    // numVerticesOnCell (global) kept alive: computeEdgeNormalDirection needs it
    // Build local maxLevelCell (size nLocalCells)
    int* localMaxLevelCell = new int[this->nLocalCells];
    for(int lc = 0; lc < this->nLocalCells; lc++) {
        int gc = this->localCell2GlobalCell[lc] - 1;
        localMaxLevelCell[lc] = this->maxLevelCell[gc];
    }
    grid->setMaxLevelCell(localMaxLevelCell);
    delete[] localMaxLevelCell;
    return grid;
}

void MPASOReader::InitSolutions(MPASOGrid* grid, Solution* &pSolution, Solution* &vSolution)
{
    std::cout << "[MPASOReader::InitSolutions] Initializing solutions..." << std::endl;
    if (this->m_timestepFileIdx.empty())
        this->buildTimeIndex();
    this->readAllTimestamps();
    assert(this->nVertLevels > 0);
    assert(this->nTimestepsTotal > 0);
    grid->setNVertLevels(this->nVertLevels);

    int N = (cfg.isTimeVarying && cfg.isTimeVarying.value() && cfg.loadNTimeSteps)
                ? cfg.loadNTimeSteps.value() : 1;
    N = std::min(N, this->nTimestepsTotal);
    grid->setNTimestepsLoaded(N);
    std::cout << "[MPASOReader::InitSolutions] Loading global timestep 0 - "
              << N - 1 << " of " << this->nTimestepsTotal
              << " timestep(s)." << std::endl;

    this->buildRequiredCellSet();
    this->computeCOVweights();
    std::cout << "[MPASOReader::InitSolutions] Finished computing COV weights." << std::endl;

    int nVertNodes = this->nLocalVertices * this->nVertLevels;
    VECTOR3** vertexVelocity     = new VECTOR3*[N];
    VECTOR3** vertexVertVelocity = new VECTOR3*[N];
    for (int ts = 0; ts < N; ts++) {
        vertexVelocity[ts]     = new VECTOR3[nVertNodes];
        vertexVertVelocity[ts] = new VECTOR3[nVertNodes];
    }

    double* vertexZTop           = new double[static_cast<size_t>(N) * nVertNodes];

    // Stream: read one timestep at a time, convert immediately, then free cell-space data
    for (int ts = 0; ts < N; ts++) {
        this->loadGlobalTimestepIntoSlot(ts, ts, nVertNodes,
                                         vertexVelocity,
                                         vertexVertVelocity,
                                         vertexZTop);
        std::cout << "[MPASOReader::InitSolutions] Finished vertex conversion for global timestep "
                  << "(0, " << ts << ", " << N - 1 << ")" << std::endl;
    }

    // Clean up cell-space data (last timestep's allocation from the streaming loop)
    for (int t = 0; t < this->nTimestepsInMem; ++t) delete[] this->cellVelocity[t];
    delete[] this->cellVelocity;
    this->cellVelocity = nullptr;
    for (int t = 0; t < this->nTimestepsInMem; ++t) delete[] this->cellVertVelocity[t];
    delete[] this->cellVertVelocity;
    this->cellVertVelocity = nullptr;
    this->nTimestepsInMem = 0;
    delete[] this->cellZTop;
    this->cellZTop = nullptr;

    // Slice real timestamps for [0, N-1].
    std::vector<double> windowTs;
    if (!m_fileTimestamps.empty() && N <= (int)m_fileTimestamps.size()) {
        windowTs.assign(m_fileTimestamps.begin(), m_fileTimestamps.begin() + N);
    }
    std::vector<double> gridTs = windowTs;
    if (gridTs.empty()) {
        gridTs.reserve(N);
        for (int i = 0; i < N; i++)
            gridTs.push_back(static_cast<double>(i));
    }

    pSolution = new Solution(vertexVelocity, nVertNodes, N);
    if (!windowTs.empty()) {
        pSolution->setTimestamps(windowTs);
        pSolution->SetMinMaxTime(windowTs.front(), windowTs.back());
    } else {
        pSolution->SetMinMaxTime(0.0, (double)(N - 1));
    }
    vSolution = new Solution(vertexVertVelocity, nVertNodes, N);
    if (!windowTs.empty()) {
        vSolution->setTimestamps(windowTs);
        vSolution->SetMinMaxTime(windowTs.front(), windowTs.back());
    } else {
        vSolution->SetMinMaxTime(0.0, (double)(N - 1));
    }
    grid->setZTop(vertexZTop, N);
    grid->setZTopTimestamps(gridTs);

    for (int ts = 0; ts < N; ts++) {
        delete[] vertexVelocity[ts];
        delete[] vertexVertVelocity[ts];
    }
    delete[] vertexVelocity;
    delete[] vertexVertVelocity;
    delete[] vertexZTop;
    this->setGlobalTimestepOffset(N);
}

void MPASOReader::UpdateSolutions(MPASOGrid* grid, Solution* &pSolution, Solution* &vSolution)
{
    if (this->m_timestepFileIdx.empty())
        this->buildTimeIndex();

    int N = (cfg.isTimeVarying && cfg.isTimeVarying.value() && cfg.loadNTimeSteps)
                ? cfg.loadNTimeSteps.value() : 1;
    int nextGlobal = this->getTimestepOffset();
    if (nextGlobal >= this->nTimestepsTotal) {
        std::cerr << "[MPASOReader::UpdateSolutions] No more timesteps to load.\n";
        return;
    }

    int globalStart = (nextGlobal > 0 && N > 1) ? nextGlobal - 1 : nextGlobal;
    N = std::min(N, this->nTimestepsTotal - globalStart);
    if (N <= 0) {
        std::cerr << "[MPASOReader::UpdateSolutions] No more timesteps to load.\n";
        return;
    }

    int nVertNodes = this->nLocalVertices * this->nVertLevels;
    bool hasBoundary = (globalStart < nextGlobal && pSolution && vSolution);

    // Copy boundary timestep (start == old window's last ts) from old solutions to avoid
    // recomputing the vertex-space conversion for data we already have.
    VECTOR3* boundary_p = nullptr;
    VECTOR3* boundary_v = nullptr;
    int old_last = -1;
    if (hasBoundary) {
        old_last = pSolution->GetNTimeSteps() - 1;
        boundary_p = new VECTOR3[nVertNodes];
        boundary_v = new VECTOR3[nVertNodes];
        std::memcpy(boundary_p, pSolution->GetTimestepData(old_last), nVertNodes * sizeof(VECTOR3));
        std::memcpy(boundary_v, vSolution->GetTimestepData(old_last), nVertNodes * sizeof(VECTOR3));
    }

    VECTOR3** vertexVelocity     = new VECTOR3*[N];
    VECTOR3** vertexVertVelocity = new VECTOR3*[N];
    for (int ts = 0; ts < N; ts++) {
        vertexVelocity[ts]     = new VECTOR3[nVertNodes];
        vertexVertVelocity[ts] = new VECTOR3[nVertNodes];
    }
    double* vertexZTop           = new double[static_cast<size_t>(N) * nVertNodes];

    // Slot 0: copy boundary from old solution (no recomputation)
    if (hasBoundary && boundary_p) {
        std::memcpy(vertexVelocity[0],     boundary_p, nVertNodes * sizeof(VECTOR3));
        std::memcpy(vertexVertVelocity[0], boundary_v, nVertNodes * sizeof(VECTOR3));
        std::memcpy(vertexZTop, grid->getZTopTimestep(old_last), nVertNodes * sizeof(double));
        delete[] boundary_p; boundary_p = nullptr;
        delete[] boundary_v; boundary_v = nullptr;
        std::cerr << "[MPASOReader::UpdateSolutions] Copied boundary timestep from previous window for global timestep "
                  << globalStart << "." << std::endl;
    }

    // Fill the remaining slots by global timestep; this can span any number of files.
    int ts_out_start = hasBoundary ? 1 : 0;
    for (int ts_out = ts_out_start; ts_out < N; ++ts_out) {
        int globalTimestep = globalStart + ts_out;
        this->loadGlobalTimestepIntoSlot(globalTimestep, ts_out, nVertNodes,
                                         vertexVelocity,
                                         vertexVertVelocity,
                                         vertexZTop);
        std::cout << "[MPASOReader::UpdateSolutions] Finished vertex conversion for global timestep "
                  << "(" << globalStart << ", " << globalTimestep << ", "
                  << globalStart + N - 1 << ")" << std::endl;
    }

    // Slice timestamps using global index.
    std::vector<double> windowTs;
    if (!m_fileTimestamps.empty() && globalStart >= 0 && globalStart + N <= (int)m_fileTimestamps.size()) {
        windowTs.assign(m_fileTimestamps.begin() + globalStart,
                        m_fileTimestamps.begin() + globalStart + N);
    }
    std::vector<double> gridTs = windowTs;
    if (gridTs.empty()) {
        gridTs.reserve(N);
        for (int i = 0; i < N; i++)
            gridTs.push_back(static_cast<double>(globalStart + i));
    }

    if (pSolution) { delete pSolution; pSolution = nullptr; }
    if (vSolution) { delete vSolution; vSolution = nullptr; }
    pSolution = new Solution(vertexVelocity, nVertNodes, N);
    if (!windowTs.empty()) {
        pSolution->setTimestamps(windowTs);
        pSolution->SetMinMaxTime(windowTs.front(), windowTs.back());
    } else {
        pSolution->SetMinMaxTime((double)globalStart, (double)(globalStart + N - 1));
    }
    vSolution = new Solution(vertexVertVelocity, nVertNodes, N);
    if (!windowTs.empty()) {
        vSolution->setTimestamps(windowTs);
        vSolution->SetMinMaxTime(windowTs.front(), windowTs.back());
    } else {
        vSolution->SetMinMaxTime((double)globalStart, (double)(globalStart + N - 1));
    }
    grid->setNTimestepsLoaded(N);
    grid->setZTop(vertexZTop, N);
    grid->setZTopTimestamps(gridTs);

    // Clean up cell-space data
    for (int t = 0; t < this->nTimestepsInMem; ++t) delete[] this->cellVelocity[t];
    delete[] this->cellVelocity;
    this->cellVelocity = nullptr;
    for (int t = 0; t < this->nTimestepsInMem; ++t) delete[] this->cellVertVelocity[t];
    delete[] this->cellVertVelocity;
    this->cellVertVelocity = nullptr;
    this->nTimestepsInMem = 0;
    delete[] this->cellZTop;
    this->cellZTop = nullptr;
    for (int ts = 0; ts < N; ts++) {
        delete[] vertexVelocity[ts];
        delete[] vertexVertVelocity[ts];
    }
    delete[] vertexVelocity;
    delete[] vertexVertVelocity;
    delete[] vertexZTop;

    this->setGlobalTimestepOffset(globalStart + N);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////
///////  Read time independent variables
///////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void MPASOReader::readDimensions()
{
    int ncid = ncDataid;
    if (this->loadMesh)
        ncid = ncMeshid;

    int dimid;
    size_t len;

    nc_inq_dimid(ncid, "nCells", &dimid);
    nc_inq_dimlen(ncid, dimid, &len);
    this->nTotalCells = static_cast<int>(len);

    nc_inq_dimid(ncid, "nVertices", &dimid);
    nc_inq_dimlen(ncid, dimid, &len);
    this->nTotalVertices = static_cast<int>(len);

    nc_inq_dimid(ncid, "nEdges", &dimid);
    nc_inq_dimlen(ncid, dimid, &len);
    this->nTotalEdges = static_cast<int>(len);

    nc_inq_dimid(ncid, "vertexDegree", &dimid);
    nc_inq_dimlen(ncid, dimid, &len);
    this->verticesDegree = static_cast<int>(len);

    nc_inq_dimid(ncid, "nVertLevels", &dimid);
    nc_inq_dimlen(ncid, dimid, &len);
    this->nVertLevels = static_cast<int>(len);

    nc_inq_dimid(ncid, "nVertLevelsP1", &dimid);
    nc_inq_dimlen(ncid, dimid, &len);
    this->nVertLevelsP1 = static_cast<int>(len);

    nc_inq_dimid(ncid, "maxEdges", &dimid);
    nc_inq_dimlen(ncid, dimid, &len);
    this->nMaxEdges = static_cast<int>(len);
}

void MPASOReader::readCellCoord()
{
    std::string xcell_str = "xCell";
    std::string ycell_str = "yCell";
    std::string zcell_str = "zCell";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;

    // 1. acquire varid for each variable
    int xcell_varid, ycell_varid, zcell_varid;
    nc_inq_varid(ncid, xcell_str.c_str(), &xcell_varid);
    nc_inq_varid(ncid, ycell_str.c_str(), &ycell_varid);
    nc_inq_varid(ncid, zcell_str.c_str(), &zcell_varid);

    // 2. Assert number of total cells matches dimension read in readDimensions()
    int ndims, dimids[NC_MAX_DIMS];
    nc_inq_var(ncid, xcell_varid, nullptr, nullptr, &ndims, dimids, nullptr);
    size_t len;
    nc_inq_dimlen(ncid, dimids[0], &len);
    assert(this->nTotalCells == static_cast<int>(len));

    // 3. Allocate memory for cell coordinates
    this->cellCoord = new VECTOR3[this->nTotalCells];
    double* cellCoord_tmp = new double[this->nTotalCells];

    nc_get_var_double(ncid, xcell_varid, cellCoord_tmp);
    for (int i = 0; i < this->nTotalCells; i++) {
        this->cellCoord[i][0] = cellCoord_tmp[i];
    }

    nc_get_var_double(ncid, ycell_varid, cellCoord_tmp);
    for (int i = 0; i < this->nTotalCells; i++) {
        this->cellCoord[i][1] = cellCoord_tmp[i];
    }

    nc_get_var_double(ncid, zcell_varid, cellCoord_tmp);
    for (int i = 0; i < this->nTotalCells; i++) {
        this->cellCoord[i][2] = cellCoord_tmp[i];
    }

    delete[] cellCoord_tmp;
}

void MPASOReader::readVertexCoord()
{
    std::string xvert_str = "xVertex";
    std::string yvert_str = "yVertex";
    std::string zvert_str = "zVertex";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;

    // 1. Acquire varid for each vertex coordinate variable
    int xvert_varid, yvert_varid, zvert_varid;
    nc_inq_varid(ncid, xvert_str.c_str(), &xvert_varid);
    nc_inq_varid(ncid, yvert_str.c_str(), &yvert_varid);
    nc_inq_varid(ncid, zvert_str.c_str(), &zvert_varid);

    // 2. Assert number of total vertices matches dimension read in readDimensions()
    int ndims, dimids[NC_MAX_DIMS];
    nc_inq_var(ncid, xvert_varid, nullptr, nullptr, &ndims, dimids, nullptr);
    size_t vert_len;
    nc_inq_dimlen(ncid, dimids[0], &vert_len);
    assert(this->nTotalVertices == static_cast<int>(vert_len));

    // 3. Allocate memory for vertex coordinates and a temporary buffer
    this->vertexCoord = new VECTOR3[this->nTotalVertices];
    double* vertexCoord_tmp = new double[this->nTotalVertices];

    // 4. Read xVertex data and store into vertexCoord[][0]
    nc_get_var_double(ncid, xvert_varid, vertexCoord_tmp);
    for (int i = 0; i < this->nTotalVertices; ++i) {
        this->vertexCoord[i][0] = vertexCoord_tmp[i];
    }

    // 5. Read yVertex data and store into vertexCoord[][1]
    nc_get_var_double(ncid, yvert_varid, vertexCoord_tmp);
    for (int i = 0; i < this->nTotalVertices; ++i) {
        this->vertexCoord[i][1] = vertexCoord_tmp[i];
    }

    // 6. Read zVertex data and store into vertexCoord[][2]
    nc_get_var_double(ncid, zvert_varid, vertexCoord_tmp);
    for (int i = 0; i < this->nTotalVertices; ++i) {
        this->vertexCoord[i][2] = vertexCoord_tmp[i];
    }

    delete[] vertexCoord_tmp;
}

void MPASOReader::readVertexOnCell()
{
    std::string voc_str = "verticesOnCell";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;
    
    // 1. Acquire varid for verticesOnCell
    int voc_varid;
    nc_inq_varid(ncid, voc_str.c_str(), &voc_varid);

    // 2. Assert dimensions of verticesOnCell match values read in readDimensions()
    int ndims, dimids[NC_MAX_DIMS];
    nc_inq_var(ncid, voc_varid, nullptr, nullptr, &ndims, dimids, nullptr);
    size_t len;
    nc_inq_dimlen(ncid, dimids[0], &len);
    assert(this->nTotalCells == static_cast<int>(len));

    size_t dim_len_max_edges;
    nc_inq_dimlen(ncid, dimids[1], &dim_len_max_edges);
    assert(this->nMaxEdges == static_cast<int>(dim_len_max_edges));

    // 4. Allocate array to hold vertex indices for all cells
    int total_size = this->nTotalCells * this->nMaxEdges;
    this->verticesOnCell = new int[total_size];

    // 5. Read data into the allocated array
    nc_get_var_int(ncid, voc_varid, this->verticesOnCell);
}

void MPASOReader::readCellOnVertex()
{
    std::string cov_str = "cellsOnVertex";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;

    // 1. Acquire varid for cellsOnVertex
    int cov_varid;
    nc_inq_varid(ncid, cov_str.c_str(), &cov_varid);

    // 2. Determine dimensions of cellsOnVertex
    int ndims, dimids[NC_MAX_DIMS];
    nc_inq_var(ncid, cov_varid, nullptr, nullptr, &ndims, dimids, nullptr);

    // 3. Assert dimensions of cellsOnVertex match values read in readDimensions()
    size_t nVert, vertDeg;
    nc_inq_dimlen(ncid, dimids[0], &nVert);
    nc_inq_dimlen(ncid, dimids[1], &vertDeg);
    assert(this->nTotalVertices == static_cast<int>(nVert));
    assert(this->verticesDegree == static_cast<int>(vertDeg));

    // 4. Allocate array to hold cell indices for all vertices
    int total_size = this->nTotalVertices*this->verticesDegree;
    this->cellsOnVertex = new int[total_size];

    // 5. Read data into the allocated array
    nc_get_var_int(ncid, cov_varid, this->cellsOnVertex);
}

void MPASOReader::readCellOnCell()
{
    std::string coc_str = "cellsOnCell";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;

    // 1. Acquire varid for cellsOnCell
    int coc_varid;
    nc_inq_varid(ncid, coc_str.c_str(), &coc_varid);

    // 2. Determine the dimensions of cellsOnCell
    int ndims, dimids[NC_MAX_DIMS];
    nc_inq_var(ncid, coc_varid, nullptr, nullptr, &ndims, dimids, nullptr);

    // 3. Query the lengths of the first two dimensions
    size_t len_cells, len_edges;
    nc_inq_dimlen(ncid, dimids[0], &len_cells);
    nc_inq_dimlen(ncid, dimids[1], &len_edges);

    // 4. Assert dimensions of cellsOnCell match values read in readDimensions()
    assert(this->nTotalCells == static_cast<int>(len_cells));
    assert(this->nMaxEdges == static_cast<int>(len_edges));

    // 5. Allocate array to hold neighboring cell indices for each cell
    int total_size = this->nTotalCells * this->nMaxEdges;
    this->cellsOnCell = new int[total_size];

    // 6. Read the data into the allocated array
    nc_get_var_int(ncid, coc_varid, this->cellsOnCell);
}

void MPASOReader::readNumVertexOnCell()
{
    std::string nvoc_str = "nEdgesOnCell";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;

    // 1. Acquire varid for nEdgesOnCell
    int nvoc_varid;
    nc_inq_varid(ncid, nvoc_str.c_str(), &nvoc_varid);
    
    int ndims, dimids[NC_MAX_DIMS];
    nc_inq_var(ncid, nvoc_varid, nullptr, nullptr, &ndims, dimids, nullptr);
    size_t len_cells;
    nc_inq_dimlen(ncid, dimids[0], &len_cells);

    assert(this->nTotalCells == static_cast<int>(len_cells));

    // 2. Allocate array to hold the number of vertices per cell
    this->numVerticesOnCell = new int[this->nTotalCells];

    // 3. Read data into the allocated array
    nc_get_var_int(ncid, nvoc_varid, this->numVerticesOnCell);
}

void MPASOReader::readMaxLevelCell()
{
    std::string mlc_str = "maxLevelCell";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;

    // 1. Acquire varid for maxLevelCell
    int mlc_varid;
    if ((retval = nc_inq_varid(ncid, mlc_str.c_str(), &mlc_varid))) {
        std::cerr << "[MPASOReader] nc_inq_varid maxLevelCell failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    int ndims, dimids[NC_MAX_DIMS];
    if ((retval = nc_inq_var(ncid, mlc_varid, nullptr, nullptr,
                             &ndims, dimids, nullptr))) {
        std::cerr << "[MPASOReader] nc_inq_var maxLevelCell failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    size_t len_cells;
    if ((retval = nc_inq_dimlen(ncid, dimids[0], &len_cells))) {
        std::cerr << "[MPASOReader] nc_inq_dimlen maxLevelCell failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    assert(this->nTotalCells == static_cast<int>(len_cells));

    // 2. Allocate array to hold max level for each cell
    delete[] this->maxLevelCell;
    this->maxLevelCell = new int[this->nTotalCells];

    // 3. Read data into the allocated array
    if ((retval = nc_get_var_int(ncid, mlc_varid, this->maxLevelCell))) {
        std::cerr << "[MPASOReader] nc_get_var_int maxLevelCell failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
}

void MPASOReader::readBottomDepth()
{
    std::string bd_str = "bottomDepth";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;

    // 1. Acquire varid for bottomDepth
    int bd_varid;
    if ((retval = nc_inq_varid(ncid, bd_str.c_str(), &bd_varid))) {
        std::cerr << "[MPASOReader] nc_inq_varid bottomDepth failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    int ndims, dimids[NC_MAX_DIMS];
    if ((retval = nc_inq_var(ncid, bd_varid, nullptr, nullptr,
                             &ndims, dimids, nullptr))) {
        std::cerr << "[MPASOReader] nc_inq_var bottomDepth failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    size_t len_cells;
    if ((retval = nc_inq_dimlen(ncid, dimids[0], &len_cells))) {
        std::cerr << "[MPASOReader] nc_inq_dimlen bottomDepth failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    assert(this->nTotalCells == static_cast<int>(len_cells));

    delete[] this->bottomDepth;
    this->bottomDepth = new double[this->nTotalCells];

    if ((retval = nc_get_var_double(ncid, bd_varid, this->bottomDepth))) {
        std::cerr << "[MPASOReader] nc_get_var_double bottomDepth failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
}

void MPASOReader::readVerticesOnEdge()
{
    std::string voe_str = "verticesOnEdge";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;

    // 1. Acquire varid for verticesOnEdge
    int voe_varid;
    nc_inq_varid(ncid, voe_str.c_str(), &voe_varid);

    // 2. Determine dimensions of verticesOnEdge
    int ndims, dimids[NC_MAX_DIMS];
    nc_inq_var(ncid, voe_varid, nullptr, nullptr, &ndims, dimids, nullptr);
    size_t len_edge;
    nc_inq_dimlen(ncid, dimids[0], &len_edge);

    assert(this->nTotalEdges == static_cast<int>(len_edge));

    size_t dim_len_2;
    nc_inq_dimlen(ncid, dimids[1], &dim_len_2);
    assert(dim_len_2 == 2);

    // 4. Allocate array to hold vertex indices for all edges
    int total_size = this->nTotalEdges * 2;
    this->verticesOnEdge = new int[total_size];

    // 5. Read data into the allocated array
    nc_get_var_int(ncid, voe_varid, this->verticesOnEdge);
}

void MPASOReader::readEdgesOnCell()
{
    std::string eoc_str = "edgesOnCell";
    int ncid = ncDataid; // default to data file id

    if(this->loadMesh)
        ncid = ncMeshid;

    // 1. Acquire varid for edgesOnCell
    int eoc_varid;
    nc_inq_varid(ncid, eoc_str.c_str(), &eoc_varid);

    // 2. Determine dimensions of edgesOnCell
    int ndims, dimids[NC_MAX_DIMS];
    nc_inq_var(ncid, eoc_varid, nullptr, nullptr, &ndims, dimids, nullptr);
    size_t len_cells;
    nc_inq_dimlen(ncid, dimids[0], &len_cells);

    assert(this->nTotalCells == static_cast<int>(len_cells));

    size_t dim_len_max_edges;
    nc_inq_dimlen(ncid, dimids[1], &dim_len_max_edges);
    assert(dim_len_max_edges == static_cast<size_t>(this->nMaxEdges));

    // 4. Allocate array to hold edge indices for all cells
    int total_size = this->nTotalCells * this->nMaxEdges;
    this->edgesOnCell = new int[total_size];

    // 5. Read data into the allocated array
    nc_get_var_int(ncid, eoc_varid, this->edgesOnCell);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////
///////  Read time varying variables
///////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void MPASOReader::readCellVelocity(int startTimestep, int N)
{
    int velX_id, velY_id, velZ_id;
    int dimids[NC_MAX_VAR_DIMS], ndims;
    size_t dimlen0, dimlen1, dimlen2;
    size_t start[3], count[3];
    int cellNodeNum;
    double* cellVel_tmp;

    assert(this->cfg.velocity.cartesian.X);
    assert(this->cfg.velocity.cartesian.Y);
    assert(this->cfg.velocity.cartesian.Z);
    assert(N > 0);

    // 1. Get variable IDs for velocity components
    if ((retval = nc_inq_varid(ncDataid, this->cfg.velocity.cartesian.X->c_str(), &velX_id))) {
        std::cerr << "nc_inq_varid velocityX failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_varid(ncDataid, this->cfg.velocity.cartesian.Y->c_str(), &velY_id))) {
        std::cerr << "nc_inq_varid velocityY failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_varid(ncDataid, this->cfg.velocity.cartesian.Z->c_str(), &velZ_id))) {
        std::cerr << "nc_inq_varid velocityZ failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 2. query the dimensions of velocityX variable
    if ((retval = nc_inq_var(ncDataid, velX_id,
                             /*name=*/ nullptr,
                             /*xtype=*/ nullptr,
                             &ndims,
                             dimids,
                             /*natts=*/ nullptr))) {
        std::cerr << "nc_inq_var failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if (ndims != 3) {
        std::cerr << "Unexpected number of dims for velocityX: " << ndims << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 3. read the dimension lengths
    if ((retval = nc_inq_dimlen(ncDataid, dimids[0], &dimlen0))) {
        std::cerr << "nc_inq_dimlen dim0 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[1], &dimlen1))) {
        std::cerr << "nc_inq_dimlen dim1 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[2], &dimlen2))) {
        std::cerr << "nc_inq_dimlen dim2 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 4. initialize nTimestepsInFile; assert cells and levels match readDimensions()
    if (this->nTimestepsInFile <= 0) this->nTimestepsInFile = static_cast<int>(dimlen0);
    assert(this->nTimestepsInFile  == static_cast<int>(dimlen0));
    assert(this->nTotalCells == static_cast<int>(dimlen1));
    assert(this->nVertLevels == static_cast<int>(dimlen2));
    assert(startTimestep + N <= this->nTimestepsInFile);

    cellNodeNum = this->nTotalCells * this->nVertLevels;

    if (this->requiredCellIndices.empty())
        this->buildRequiredCellSet();

    // 5. free any previously loaded data and allocate N slots
    if (this->cellVelocity != nullptr) {
        for (int t = 0; t < this->nTimestepsInMem; ++t)
            delete[] this->cellVelocity[t];
        delete[] this->cellVelocity;
    }
    this->cellVelocity  = new VECTOR3*[N];
    this->nTimestepsInMem = N;
    for (int t = 0; t < N; ++t)
        this->cellVelocity[t] = new VECTOR3[cellNodeNum];

    // 6. read N timesteps
    cellVel_tmp = new double[cellNodeNum];
    count[0] = 1;
    count[1] = this->nTotalCells;
    count[2] = this->nVertLevels;
    start[1]  = 0;
    start[2]  = 0;

    for (int t = 0; t < N; ++t) {
        start[0] = startTimestep + t;

        auto copyCartesianComponent = [&](int dim) {
            if (!this->requiredCellIndices.empty()) {
                for (int gCellIdx : this->requiredCellIndices) {
                    int offset = gCellIdx * this->nVertLevels;
                    for (int j = 0; j < this->nVertLevels; j++) {
                        double val = cellVel_tmp[offset+j];
                        this->cellVelocity[t][offset+j][dim] =
                            (val > 1.0e30 || val < -1.0e30) ? 0.0 : val;
                    }
                }
            } else {
                for (int i = 0; i < cellNodeNum; ++i) {
                    double val = cellVel_tmp[i];
                    this->cellVelocity[t][i][dim] =
                        (val > 1.0e30 || val < -1.0e30) ? 0.0 : val;
                }
            }
        };

        if ((retval = nc_get_vara_double(ncDataid, velX_id, start, count, cellVel_tmp))) {
            std::cerr << "nc_get_vara_double velocityX t=" << start[0] << " failed: " << nc_strerror(retval) << "\n";
            std::exit(EXIT_FAILURE);
        }
        copyCartesianComponent(0);

        if ((retval = nc_get_vara_double(ncDataid, velY_id, start, count, cellVel_tmp))) {
            std::cerr << "nc_get_vara_double velocityY t=" << start[0] << " failed: " << nc_strerror(retval) << "\n";
            std::exit(EXIT_FAILURE);
        }
        copyCartesianComponent(1);

        if ((retval = nc_get_vara_double(ncDataid, velZ_id, start, count, cellVel_tmp))) {
            std::cerr << "nc_get_vara_double velocityZ t=" << start[0] << " failed: " << nc_strerror(retval) << "\n";
            std::exit(EXIT_FAILURE);
        }
        copyCartesianComponent(2);
    }

    delete[] cellVel_tmp;
}

void MPASOReader::readZonMeridVelocity(int startTimestep, int N)
{
    // Variable names in C-string form
    const char* velZon_cstr = this->cfg.velocity.spherical.zonal->c_str();
    const char* velMer_cstr = this->cfg.velocity.spherical.meridional->c_str();
    int velZon_varid, velMer_varid;
    int dimids[NC_MAX_VAR_DIMS], ndims;
    size_t dimlen0, dimlen1, dimlen2;
    size_t start[3], count[3];
    int cellNodeNum;
    double* cellZonVel_tmp;
    double* cellMerVel_tmp;

    assert(this->cfg.velocity.spherical.zonal);
    assert(this->cfg.velocity.spherical.meridional);
    assert(N > 0);

    // 1. Acquire varid for each velocity component
    retval = nc_inq_varid(this->ncDataid, velZon_cstr, &velZon_varid);
    if (retval != NC_NOERR) {
        std::cerr << "[MPASOReader] Error acquiring varid for " << velZon_cstr << ": " << nc_strerror(retval) << std::endl;
        exit(1);
    }
    retval = nc_inq_varid(this->ncDataid, velMer_cstr, &velMer_varid);
    if (retval != NC_NOERR) {
        std::cerr << "[MPASOReader] Error acquiring varid for " << velMer_cstr << ": " << nc_strerror(retval) << std::endl;
        exit(1);
    }

    // 2. query the dimensions of velocityX variable
    if ((retval = nc_inq_var(ncDataid, velZon_varid,
                             /*name=*/ nullptr,
                             /*xtype=*/ nullptr,
                             &ndims,
                             dimids,
                             /*natts=*/ nullptr))) {
        std::cerr << "nc_inq_var failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if (ndims != 3) {
        std::cerr << "Unexpected number of dims for velocityX: " << ndims << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 3. read the dimension lengths
    if ((retval = nc_inq_dimlen(ncDataid, dimids[0], &dimlen0))) {
        std::cerr << "nc_inq_dimlen dim0 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[1], &dimlen1))) {
        std::cerr << "nc_inq_dimlen dim1 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[2], &dimlen2))) {
        std::cerr << "nc_inq_dimlen dim2 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 4. initialize nTimestepsInFile; assert cells and levels match readDimensions()
    if (this->nTimestepsInFile <= 0) this->nTimestepsInFile = static_cast<int>(dimlen0);
    assert(this->nTimestepsInFile  == static_cast<int>(dimlen0));
    assert(this->nTotalCells == static_cast<int>(dimlen1));
    assert(this->nVertLevels == static_cast<int>(dimlen2));
    assert(startTimestep + N <= this->nTimestepsInFile);

    cellNodeNum = this->nTotalCells * this->nVertLevels;

    if (this->requiredCellIndices.empty())
        this->buildRequiredCellSet();

    // 5. free any previously loaded data and allocate N slots
    if (this->cellVelocity != nullptr) {
        for (int t = 0; t < this->nTimestepsInMem; ++t)
            delete[] this->cellVelocity[t];
        delete[] this->cellVelocity;
    }
    this->cellVelocity  = new VECTOR3*[N];
    this->nTimestepsInMem = N;
    for (int t = 0; t < N; ++t) {
        this->cellVelocity[t] = new VECTOR3[cellNodeNum];
        for (int i = 0; i < cellNodeNum; ++i)
            this->cellVelocity[t][i].Zero();
    }

    // 6. read N timesteps
    cellZonVel_tmp = new double[cellNodeNum];
    cellMerVel_tmp = new double[cellNodeNum];
    count[0] = 1;
    count[1] = this->nTotalCells;
    count[2] = this->nVertLevels;
    start[1]  = 0;
    start[2]  = 0;

    for (int t = 0; t < N; ++t) {
        start[0] = startTimestep + t;

        if ((retval = nc_get_vara_double(ncDataid, velZon_varid, start, count, cellZonVel_tmp))) {
            std::cerr << "nc_get_vara_double zonal velocity t=" << start[0] << " failed: " << nc_strerror(retval) << "\n";
            std::exit(EXIT_FAILURE);
        }
        if ((retval = nc_get_vara_double(ncDataid, velMer_varid, start, count, cellMerVel_tmp))) {
            std::cerr << "nc_get_vara_double meridional velocity t=" << start[0] << " failed: " << nc_strerror(retval) << "\n";
            std::exit(EXIT_FAILURE);
        }
        for (int gCellIdx : this->requiredCellIndices) {
            int offset = gCellIdx * this->nVertLevels;
            for (int j = 0; j < this->nVertLevels; j++) {
                double uzon = cellZonVel_tmp[offset+j];
                double umer = cellMerVel_tmp[offset+j];
                // Mask out fill values (9.96921e+36) — land/dry cells have no valid velocity
                if (uzon > 1.0e30 || uzon < -1.0e30) uzon = 0.0;
                if (umer > 1.0e30 || umer < -1.0e30) umer = 0.0;
                convertENUVelocityToXYZ(this->cellCoord[gCellIdx], uzon, umer, 0.0, this->cellVelocity[t][offset+j]);
            }
        }
    }

    delete[] cellZonVel_tmp;
    delete[] cellMerVel_tmp;
}

void MPASOReader::readNormalVelocity(int startTimestep, int N) {
    // normal velocity: timestep, edge, level
    std::cout << "[MPASOReader] Reading normal velocity..." << std::endl;
    const char* normalVel_cstr = this->cfg.velocity.normal->c_str();
    int normalVel_varid;
    int dimids[NC_MAX_VAR_DIMS], ndims;
    size_t dimlen0, dimlen1, dimlen2;
    size_t start[3], count[3];
    int edgeNum;
    double* edgeNormalVel_tmp;

    assert(this->cfg.velocity.normal);
    assert(N > 0);

    // 1. Acquire varid for normal velocity
    retval = nc_inq_varid(this->ncDataid, normalVel_cstr, &normalVel_varid);
    if (retval != NC_NOERR) {
        std::cerr << "[MPASOReader] Error acquiring varid for " << normalVel_cstr << ": " << nc_strerror(retval) << std::endl;
        exit(1);
    }

    // 2. query the dimensions of normal velocity variable
    if ((retval = nc_inq_var(ncDataid, normalVel_varid,
                             /*name=*/ nullptr,
                             /*xtype=*/ nullptr,
                             &ndims,
                             dimids,
                             /*natts=*/ nullptr))) {
        std::cerr << "nc_inq_var failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if (ndims != 3) {
        std::cerr << "Unexpected number of dims for normal velocity: " << ndims << "\n";
        std::exit(EXIT_FAILURE);
    }
    // 3. read the dimension lengths
    if ((retval = nc_inq_dimlen(ncDataid, dimids[0], &dimlen0))) {
        std::cerr << "nc_inq_dimlen dim0 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[1], &dimlen1))) {
        std::cerr << "nc_inq_dimlen dim1 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[2], &dimlen2))) {
        std::cerr << "nc_inq_dimlen dim2 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    // 4. initialize nTimestepsInFile; assert edges and levels match readDimensions()
    if (this->nTimestepsInFile <= 0) this->nTimestepsInFile = static_cast<int>(dimlen0);
    assert(this->nTimestepsInFile  == static_cast<int>(dimlen0));
    assert(this->nTotalEdges == static_cast<int>(dimlen1));
    assert(this->nVertLevels == static_cast<int>(dimlen2));
    assert(startTimestep + N <= this->nTimestepsInFile);

    edgeNum = this->nTotalEdges * this->nVertLevels;

    if (this->requiredCellIndices.empty())
        this->buildRequiredCellSet();

    // 5. free any previously loaded data and allocate N slots
    if (this->cellVelocity != nullptr) {
        for (int t = 0; t < this->nTimestepsInMem; ++t)
            delete[] this->cellVelocity[t];
        delete[] this->cellVelocity;
    }
    this->cellVelocity  = new VECTOR3*[N];
    this->nTimestepsInMem = N;
    for (int t = 0; t < N; ++t)
        this->cellVelocity[t] = nullptr;   // computeCellVelocity asserts nullptr and allocates

    // 6. read N timesteps and convert each to cell-centered velocity
    edgeNormalVel_tmp = new double[edgeNum];
    count[0] = 1;
    count[1] = this->nTotalEdges;
    count[2] = this->nVertLevels;
    start[1]  = 0;
    start[2]  = 0;

    for (int t = 0; t < N; ++t) {
        start[0] = startTimestep + t;
        if ((retval = nc_get_vara_double(ncDataid, normalVel_varid, start, count, edgeNormalVel_tmp))) {
            std::cerr << "nc_get_vara_double normal velocity t=" << start[0] << " failed: " << nc_strerror(retval) << "\n";
            std::exit(EXIT_FAILURE);
        }
        std::cout << "[MPASOReader] Converting normal velocity t=" << start[0] << " to cell-centered velocity..." << std::endl;
        this->computeCellVelocity(edgeNormalVel_tmp, this->cellVelocity[t]);
    }

    delete[] edgeNormalVel_tmp;
}

void MPASOReader::readVertVelocityTop(int startTimestep, int N)
{
    int varid;
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    size_t dimlen0, dimlen1, dimlen2;
    size_t start[3], count[3];
    double* vertVelocity_tmp;

    assert(this->cfg.vertVelocityTop);
    assert(N > 0);

    // 1. Get the variable ID for vertVelocityTop
    if ((retval = nc_inq_varid(ncDataid, this->cfg.vertVelocityTop->c_str(), &varid))) {
        std::cerr << "nc_inq_varid vertVelocityTop failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 2. Query the dimensions of vertVelocityTop
    if ((retval = nc_inq_var(ncDataid, varid,
                             /*name*/    nullptr,
                             /*xtype*/   nullptr,
                             &ndims,
                             dimids,
                             /*natts*/   nullptr))) {
        std::cerr << "nc_inq_var failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if (ndims != 3) {
        std::cerr << "Unexpected number of dims for vertVelocityTop: "
                  << ndims << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 3. Query (time, cell, level+1)
    if ((retval = nc_inq_dimlen(ncDataid, dimids[0], &dimlen0))) {
        std::cerr << "nc_inq_dimlen dim0 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[1], &dimlen1))) {
        std::cerr << "nc_inq_dimlen dim1 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[2], &dimlen2))) {
        std::cerr << "nc_inq_dimlen dim2 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 4. initialize nTimestepsInFile; assert cells and levelsP1 match readDimensions()
    if (this->nTimestepsInFile <= 0) this->nTimestepsInFile = static_cast<int>(dimlen0);
    assert(this->nTimestepsInFile  == static_cast<int>(dimlen0));
    assert(this->nTotalCells       == static_cast<int>(dimlen1));
    assert(this->nVertLevelsP1     == static_cast<int>(dimlen2));
    assert(startTimestep + N <= this->nTimestepsInFile);

    // 5. free any previously loaded data and allocate N slots
    if (this->cellVertVelocity != nullptr) {
        for (int t = 0; t < this->nTimestepsInMem; ++t)
            delete[] this->cellVertVelocity[t];
        delete[] this->cellVertVelocity;
    }
    int totalBottomLevels = this->nVertLevelsP1 - 1;
    assert(totalBottomLevels == this->nVertLevels);
    if (this->requiredCellIndices.empty())
        this->buildRequiredCellSet();

    this->cellVertVelocity  = new double*[N];
    this->nTimestepsInMem = N;
    for (int t = 0; t < N; ++t)
        this->cellVertVelocity[t] = new double[this->nTotalCells * totalBottomLevels];

    // 6. read N timesteps
    vertVelocity_tmp = new double[this->nTotalCells * this->nVertLevelsP1];
    count[0] = 1;
    count[1] = this->nTotalCells;
    count[2] = this->nVertLevelsP1;
    start[1]  = 0;
    start[2]  = 0;

    for (int t = 0; t < N; ++t) {
        start[0] = startTimestep + t;
        if ((retval = nc_get_vara_double(ncDataid, varid, start, count, vertVelocity_tmp))) {
            std::cerr << "nc_get_vara_double vertVelocityTop t=" << start[0] << " failed: "
                      << nc_strerror(retval) << "\n";
            std::exit(EXIT_FAILURE);
        }
        auto copyVertVelocity = [&](int gCellIdx) {
            size_t offset   = static_cast<size_t>(gCellIdx) * totalBottomLevels;
            size_t offsetp1 = static_cast<size_t>(gCellIdx) * this->nVertLevelsP1;
            for (int v = 0; v < totalBottomLevels; ++v) {
                double val = vertVelocity_tmp[offsetp1 + v];
                this->cellVertVelocity[t][offset + v] = (val > 1.0e30 || val < -1.0e30) ? 0.0 : val;
            }
        };

        if (!this->requiredCellIndices.empty()) {
            for (int gCellIdx : this->requiredCellIndices)
                copyVertVelocity(gCellIdx);
        } else {
            for (int gCellIdx = 0; gCellIdx < this->nTotalCells; ++gCellIdx)
                copyVertVelocity(gCellIdx);
        }
    }

    delete[] vertVelocity_tmp;
}


void MPASOReader::readZTop(int timestep)
{
    int varid;
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    size_t dimlen0, dimlen1, dimlen2;
    size_t start[3], count[3];

    assert(this->cfg.zTop);

    // 1. Get the variable ID for zTop
    if ((retval = nc_inq_varid(ncDataid, this->cfg.zTop->c_str(), &varid))) {
        std::cerr << "nc_inq_varid zTop failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 2. Query the dimensions of zTop
    if ((retval = nc_inq_var(ncDataid, varid,
                             /*name*/    nullptr,
                             /*xtype*/   nullptr,
                             &ndims,
                             dimids,
                             /*natts*/   nullptr))) {
        std::cerr << "nc_inq_var failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if (ndims != 3) {
        std::cerr << "Unexpected number of dims for zTop: "
                  << ndims << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 3. Query (time, cell, level)
    if ((retval = nc_inq_dimlen(ncDataid, dimids[0], &dimlen0))) {
        std::cerr << "nc_inq_dimlen dim0 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[1], &dimlen1))) {
        std::cerr << "nc_inq_dimlen dim1 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[2], &dimlen2))) {
        std::cerr << "nc_inq_dimlen dim2 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 4. initialize nTimesteps; assert cells and levels match readDimensions()
    if (this->nTimestepsInFile <= 0) this->nTimestepsInFile = static_cast<int>(dimlen0);
    assert(this->nTimestepsInFile  == static_cast<int>(dimlen0));
    assert(this->nTotalCells == static_cast<int>(dimlen1));
    assert(this->nVertLevels == static_cast<int>(dimlen2));

    // 5. Set up start and count vectors for reading
    start[0] = timestep;               count[0] = 1;
    start[1] = 0;                      count[1] = this->nTotalCells;
    start[2] = 0;                      count[2] = this->nVertLevels;

    // 6. Read zTop directly into cellZTop
    size_t total = static_cast<size_t>(this->nTotalCells) * this->nVertLevels;
    delete[] this->cellZTop;
    this->cellZTop = new double[total];
    if ((retval = nc_get_vara_double(ncDataid, varid,
                                     start, count,
                                     this->cellZTop))) {
        std::cerr << "nc_get_vara_double zTop failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
}

void MPASOReader::readLayerThickness(int timestep)
{
    int varid;
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    size_t dimlen0, dimlen1, dimlen2;
    size_t start[3], count[3];
    double* layerThickness_tmp;

    assert(this->cfg.layerThickness);

    // 1. Get the variable ID for layerThickness
    if ((retval = nc_inq_varid(ncDataid, this->cfg.layerThickness->c_str(), &varid))) {
        std::cerr << "[MPASOReader] nc_inq_varid layerThickness failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 2. Query the dimensions of layerThickness
    if ((retval = nc_inq_var(ncDataid, varid,
                             /*name*/    nullptr,
                             /*xtype*/   nullptr,
                             &ndims,
                             dimids,
                             /*natts*/   nullptr))) {
        std::cerr << "nc_inq_var failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if (ndims != 3) {
        std::cerr << "Unexpected number of dims for layerThickness: "
                  << ndims << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 3. Query (time, cell, level)
    if ((retval = nc_inq_dimlen(ncDataid, dimids[0], &dimlen0))) {
        std::cerr << "nc_inq_dimlen dim0 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[1], &dimlen1))) {
        std::cerr << "nc_inq_dimlen dim1 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_dimlen(ncDataid, dimids[2], &dimlen2))) {
        std::cerr << "nc_inq_dimlen dim2 failed: " << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 4. Initialize the number of timesteps, total cells, and vertical levels
    if (this->nTimestepsInFile  <= 0) this->nTimestepsInFile  = static_cast<int>(dimlen0);
    if (this->nTotalCells <= 0) this->nTotalCells = static_cast<int>(dimlen1);
    if (this->nVertLevels <= 0) this->nVertLevels = static_cast<int>(dimlen2);

    assert(this->nTimestepsInFile  == static_cast<int>(dimlen0));
    assert(this->nTotalCells == static_cast<int>(dimlen1));
    assert(this->nVertLevels == static_cast<int>(dimlen2));

    // 5. Set up start and count vectors for reading
    start[0] = timestep;                   count[0] = 1;
    start[1] = 0;                          count[1] = this->nTotalCells;
    start[2] = 0;                          count[2] = this->nVertLevels;

    // 6. Allocate memory for temporary layer thickness storage
    size_t total = static_cast<size_t>(this->nTotalCells) * this->nVertLevels;
    layerThickness_tmp = new double[total];
    if ((retval = nc_get_vara_double(ncDataid, varid,
                                     start, count,
                                     layerThickness_tmp))) {
        std::cerr << "nc_get_vara_double layerThickness failed: "
                  << nc_strerror(retval) << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 7. Compute zTop from layerThickness when zTop is unavailable
    delete[] this->cellZTop;
    if(this->maxLevelCell == nullptr)
        this->readMaxLevelCell();
    if(this->bottomDepth == nullptr)
        this->readBottomDepth();
    this->cellZTop = new double[total];

    if (this->requiredCellIndices.empty() &&
        this->nLocalVertices > 0 &&
        this->localVert2GlobalVert != nullptr &&
        this->cellsOnVertex != nullptr) {
        this->buildRequiredCellSet();
    }

    auto computeCellZTop = [&](int gCellIdx) {
        int offset = gCellIdx * this->nVertLevels;
        int activeLevels = this->maxLevelCell[gCellIdx];
        if(activeLevels < 0 || activeLevels > this->nVertLevels) {
            std::cerr << "[MPASOReader] Invalid maxLevelCell[" << gCellIdx << "] = "
                      << activeLevels << "\n";
            std::exit(EXIT_FAILURE);
        }
        double surface = -this->bottomDepth[gCellIdx];
        for(int j = 0; j < activeLevels; j++)
            surface += layerThickness_tmp[offset+j];
        for(int j = 0; j < this->nVertLevels; j++) {
            if(j == 0)
                this->cellZTop[offset+j] = surface;
            else
                this->cellZTop[offset+j] = this->cellZTop[offset+j-1] - layerThickness_tmp[offset+j-1];
        }
    };

    if (!this->requiredCellIndices.empty()) {
        for (int gCellIdx : this->requiredCellIndices)
            computeCellZTop(gCellIdx);
    } else {
        for (int gCellIdx = 0; gCellIdx < this->nTotalCells; gCellIdx++)
            computeCellZTop(gCellIdx);
    }

    delete[] layerThickness_tmp;
}

void MPASOReader::readAllTimestamps()
{
    if (!this->cfg.xtime || this->cfg.xtime->empty()) {
        std::cerr << "[MPASOReader::readAllTimestamps] xtime variable not configured, skipping.\n";
        return;
    }
    if (!m_fileTimestamps.empty())
        return;  // already read
    if (this->m_timestepFileIdx.empty())
        this->buildTimeIndex();

    double t0 = 0.0;
    m_fileTimestamps.clear();
    m_fileTimestamps.reserve(this->nTimestepsTotal);

    int nFiles = this->dataFiles.empty() ? 1 : static_cast<int>(this->dataFiles.size());
    for (int fileIdx = 0; fileIdx < nFiles; fileIdx++) {
        if (fileIdx < static_cast<int>(this->m_fileTimestepCount.size()) &&
            this->m_fileTimestepCount[fileIdx] <= 0)
            continue;

        if (fileIdx == this->currentFileIdx) {
            this->appendTimestampsFromNcid(this->ncDataid, t0, m_fileTimestamps.empty());
        } else if (!this->dataFiles.empty()) {
            int tmpid;
            if (nc_open(this->dataFiles[fileIdx].c_str(), NC_NOWRITE, &tmpid) == NC_NOERR) {
                this->appendTimestampsFromNcid(tmpid, t0, m_fileTimestamps.empty());
                nc_close(tmpid);
            }
        }
    }

    if (!m_fileTimestamps.empty() && (int)m_fileTimestamps.size() != this->nTimestepsTotal) {
        std::cerr << "[MPASOReader::readAllTimestamps] Warning: read "
                  << m_fileTimestamps.size() << " timestamp(s), expected "
                  << this->nTimestepsTotal << ". Falling back to index time for incomplete windows.\n";
        m_fileTimestamps.clear();
    }

    std::cout << "[MPASOReader::readAllTimestamps] Read " << m_fileTimestamps.size()
              << " timestamps across " << nFiles << " file(s)."
              << " t[0]=0.0, t[last]=" << (m_fileTimestamps.empty() ? 0.0 : m_fileTimestamps.back()) << " s\n";
}

void MPASOReader::GetNeighborIndices(std::vector<int> &neighborInd)
{
    if (this->neighborAreaIndices.empty()) {
        std::cerr << "[Warning] neighborAreaIndices is empty! areaId = " << this->areaId << std::endl;
        return;
    }
    for(int neighboridx: this->neighborAreaIndices)
        neighborInd.push_back(neighboridx);
}

void MPASOReader::GetLocalCell2GlobalCell(int* &LC2GC, int &nLC)
{
    assert(LC2GC == nullptr);
    LC2GC = new int[this->nLocalCells];
    nLC = this->nLocalCells;
    for(int lcidx = 0; lcidx < nLC; lcidx++)
        LC2GC[lcidx] = this->localCell2GlobalCell[lcidx];
}

void MPASOReader::GetGlobalCell2LocalCell(int* &GC2LC, int &nGC)
{
    assert(GC2LC == nullptr);
    GC2LC = new int[this->nTotalCells];
    nGC = this->nTotalCells;
    for(int gcidx = 0; gcidx < nGC; gcidx++) {
        GC2LC[gcidx] = this->globalCell2LocalCell[gcidx];
    }
}

void MPASOReader::GetnVertLevels(int &nVLvl)
{
    nVLvl = this->nVertLevels;
}

void MPASOReader::buildRequiredCellSet()
{
    if (!this->requiredCellIndices.empty())
        return;

    assert(this->nLocalVertices > 0);
    assert(this->nTotalCells > 0);
    assert(this->verticesDegree > 0);
    assert(this->localVert2GlobalVert != nullptr);
    assert(this->cellsOnVertex != nullptr);

    std::vector<unsigned char> requiredCellMask(this->nTotalCells, 0);
    this->requiredCellIndices.reserve(this->nLocalVertices * this->verticesDegree);

    for (int localVertIdx = 0; localVertIdx < this->nLocalVertices; ++localVertIdx) {
        int gVertIdx = this->localVert2GlobalVert[localVertIdx] - 1;
        if (gVertIdx < 0 || gVertIdx >= this->nTotalVertices)
            throw std::runtime_error("localVert2GlobalVert produced out-of-range vertex index");

        int covOffset = gVertIdx * this->verticesDegree;
        for (int i = 0; i < this->verticesDegree; ++i) {
            int gCell1 = this->cellsOnVertex[covOffset + i];
            if (gCell1 <= 0)
                continue;

            int gCellIdx = gCell1 - 1;
            if (gCellIdx < 0 || gCellIdx >= this->nTotalCells)
                throw std::runtime_error("cellsOnVertex produced out-of-range cell index");

            if (!requiredCellMask[gCellIdx]) {
                requiredCellMask[gCellIdx] = 1;
                this->requiredCellIndices.push_back(gCellIdx);
            }
        }
    }

    std::sort(this->requiredCellIndices.begin(), this->requiredCellIndices.end());
    std::cout << "[MPASOReader] Built required cell set: "
              << this->requiredCellIndices.size() << " / " << this->nTotalCells
              << " global cells are needed by local vertices." << std::endl;
}

void MPASOReader::computeCOVweights()
{
    assert(this->cov_weights == nullptr);
    assert(this->nLocalVertices > 0);
    this->cov_weights = new VECTOR3[this->nLocalVertices];
    for (int localVertIdx = 0; localVertIdx < this->nLocalVertices; localVertIdx++) {
        int gVertIdx = this->localVert2GlobalVert[localVertIdx] - 1;
        int cellIndices[3];
        int nNeighbers = 0;
        for (int i = 0; i < 3; i++) {
            if (this->cellsOnVertex[gVertIdx*3 + i] > 0) {
                cellIndices[i] = this->cellsOnVertex[gVertIdx*3 + i] - 1;
                nNeighbers++;
            }
        }

        if (nNeighbers <= 2) {
            this->cov_weights[localVertIdx].Zero();
        } else if (nNeighbers == 3) {
            VECTOR3 diff1 = this->cellCoord[cellIndices[0]] - this->vertexCoord[gVertIdx];
            VECTOR3 diff2 = this->cellCoord[cellIndices[1]] - this->vertexCoord[gVertIdx];
            VECTOR3 diff3 = this->cellCoord[cellIndices[2]] - this->vertexCoord[gVertIdx];

            VECTOR3 diff12 = this->cellCoord[cellIndices[1]] - this->cellCoord[cellIndices[0]];
            VECTOR3 diff23 = this->cellCoord[cellIndices[2]] - this->cellCoord[cellIndices[1]];
            VECTOR3 diff31 = this->cellCoord[cellIndices[0]] - this->cellCoord[cellIndices[2]];

            double dist1 = diff1.GetMag();
            double dist2 = diff2.GetMag();
            double dist3 = diff3.GetMag();

            double dist12 = diff12.GetMag();
            double dist23 = diff23.GetMag();
            double dist31 = diff31.GetMag();

            double w1 = triangle_area_by_length(dist2, dist3, dist23);
            double w2 = triangle_area_by_length(dist1, dist3, dist31);
            double w3 = triangle_area_by_length(dist1, dist2, dist12);
            double wsum = w1 + w2 + w3;
            if (std::abs(wsum) <= 1.0e-30)
                this->cov_weights[localVertIdx].Zero();
            else
                this->cov_weights[localVertIdx].Set(w1/wsum, w2/wsum, w3/wsum);
        }
    }
}

void MPASOReader::computeEdgeNormalDirection()
{
    if (this->edgeNormalComputed) return;
    std::cout << "Compute edge normal direction\n";

    assert(this->numVerticesOnCell    != nullptr);
    assert(this->edgesOnCell          != nullptr);
    assert(this->verticesOnEdge       != nullptr);
    assert(this->vertexCoord          != nullptr);
    assert(this->cellCoord            != nullptr);
    assert(this->nTotalEdges > 0);
    assert(this->nTotalCells > 0);
    assert(this->nMaxEdges > 0);

    if (this->requiredCellIndices.empty())
        this->buildRequiredCellSet();

    // 1) Initialize
    this->edgeNormal = new VECTOR3[this->nTotalEdges];
    for (int e = 0; e < this->nTotalEdges; ++e) this->edgeNormal[e].Zero();

    // 2) Compute normals for the edges touched by required cells
    for (int gCellidx : this->requiredCellIndices) {
        if (gCellidx < 0 || gCellidx >= this->nTotalCells)
            throw std::runtime_error("requiredCellIndices contains out-of-range cell index: "
                                     + std::to_string(gCellidx));

        int nVOC = this->numVerticesOnCell[gCellidx];
        if (nVOC <= 0 || nVOC > this->nMaxEdges)
            throw std::runtime_error("numVerticesOnCell invalid at cell "
                                     + std::to_string(gCellidx) + ": " + std::to_string(nVOC));

        for (int j = 0; j < nVOC; ++j) {
            int edge1 = this->edgesOnCell[gCellidx * this->nMaxEdges + j];
            if (edge1 <= 0)
                throw std::runtime_error("edgesOnCell has non-positive id at cell="
                                         + std::to_string(gCellidx) + " j=" + std::to_string(j));
            int gEdgeIdx = edge1 - 1;
            if (gEdgeIdx < 0 || gEdgeIdx >= this->nTotalEdges)
                throw std::runtime_error("gEdgeIdx out of range: " + std::to_string(gEdgeIdx));

            // Only compute if not yet set
            if (this->edgeNormal[gEdgeIdx].GetMag() == 0.0) {
                int v1 = this->verticesOnEdge[gEdgeIdx * 2];
                int v2 = this->verticesOnEdge[gEdgeIdx * 2 + 1];
                if (v1 <= 0 || v2 <= 0)
                    throw std::runtime_error("verticesOnEdge has non-positive vertex id at edge "
                                             + std::to_string(gEdgeIdx));
                int gV1Idx = v1 - 1;
                int gV2Idx = v2 - 1;
                if (gV1Idx < 0 || gV1Idx >= this->nTotalVertices ||
                    gV2Idx < 0 || gV2Idx >= this->nTotalVertices)
                    throw std::runtime_error("vertex index out of range for edge " + std::to_string(gEdgeIdx));

                VECTOR3 edgeVec = this->vertexCoord[gV2Idx] - this->vertexCoord[gV1Idx];
                VECTOR3 tangent = this->cellCoord[gCellidx];
                VECTOR3 n = cross(edgeVec, tangent);
                double mag = n.GetMag();
                if (mag == 0.0)
                    continue;
                n.Normalize();
                this->edgeNormal[gEdgeIdx] = n;
            }
        }
    }

    this->edgeNormalComputed = true;
}

void MPASOReader::computeCellVelocity(double* normalVelocity, VECTOR3*& cellVelocity)
{
    std::cout << "Compute cell velocity from edge normal velocity\n";

    this->computeEdgeNormalDirection();

    assert(cellVelocity == nullptr);
    assert(this->nLocalCells > 0 && this->nVertLevels > 0);
    assert(this->nMaxEdges > 0);
    assert(this->localCell2GlobalCell != nullptr);
    assert(this->numVerticesOnCell    != nullptr);
    assert(this->edgesOnCell          != nullptr);
    assert(this->edgeNormal           != nullptr);
    assert(normalVelocity             != nullptr);
    assert(this->nTotalEdges > 0);

    cellVelocity = new VECTOR3[this->nTotalCells * this->nVertLevels];
    for (int i = 0; i < this->nTotalCells * this->nVertLevels; ++i)
        cellVelocity[i].Zero();

    constexpr double magic_number = 2.0; // from MPAS-Ocean doc

    auto accumulateCellVelocity = [&](int globalCellIdx) {
        int nVOC = this->numVerticesOnCell[globalCellIdx];
        if (nVOC <= 0 || nVOC > this->nMaxEdges)
            throw std::runtime_error("numVerticesOnCell invalid: nVOC=" + std::to_string(nVOC) +
                                     " for globalCellIdx=" + std::to_string(globalCellIdx));

        size_t cellOffset = size_t(globalCellIdx) * this->nVertLevels;

        for (int vl = 0; vl < this->nVertLevels; ++vl) {
            VECTOR3 acc{0,0,0};

            for (int j = 0; j < nVOC; ++j) {
                int edge1 = this->edgesOnCell[size_t(globalCellIdx) * this->nMaxEdges + j];
                if (edge1 <= 0)
                    throw std::runtime_error("edgesOnCell has non-positive edge id at cell=" +
                                             std::to_string(globalCellIdx) + " j=" + std::to_string(j));
                int gEdgeIdx = edge1 - 1;
                if (gEdgeIdx < 0 || gEdgeIdx >= this->nTotalEdges)
                    throw std::runtime_error("gEdgeIdx out of range: " + std::to_string(gEdgeIdx));

                size_t edgeOffset = size_t(gEdgeIdx) * this->nVertLevels;
                double vel = normalVelocity[edgeOffset + vl];
                if (vel > 1.0e30 || vel < -1.0e30) vel = 0.0; // mask fill values

                acc[0] += vel * this->edgeNormal[gEdgeIdx][0];
                acc[1] += vel * this->edgeNormal[gEdgeIdx][1];
                acc[2] += vel * this->edgeNormal[gEdgeIdx][2];
            }

            acc.scale(magic_number / double(nVOC));
            cellVelocity[cellOffset + vl] = acc;
        }
    };

    if (!this->requiredCellIndices.empty()) {
        for (int globalCellIdx : this->requiredCellIndices)
            accumulateCellVelocity(globalCellIdx);
    } else {
        for (int globalCellIdx = 0; globalCellIdx < this->nTotalCells; ++globalCellIdx)
            accumulateCellVelocity(globalCellIdx);
    }
}

void MPASOReader::computeTimeVaryingVar(int localVertIdx,
                                        VECTOR3* cellVelocity, double* cellVertVelocity,
                                        double* cellZTop,
                                        VECTOR3* pVelocityOut, VECTOR3* vVelocityOut,
                                        double* ztopOut)
{
    assert(this->cov_weights != nullptr);
    assert(pVelocityOut != nullptr);
    assert(vVelocityOut != nullptr);
    assert(ztopOut != nullptr);

    if (this->cov_weights[localVertIdx].GetMax() == 0.0) {
        for (int vl = 0; vl < this->nVertLevels; vl++) {
            pVelocityOut[vl].Zero();
            vVelocityOut[vl].Zero();
            ztopOut[vl] = 0.0;
        }
    } else {
        int gVertIdx = this->localVert2GlobalVert[localVertIdx] - 1;
        int cellIndices[3];
        for (int i = 0; i < 3; i++)
            cellIndices[i] = this->cellsOnVertex[gVertIdx*3 + i] - 1;

        int c1offset = cellIndices[0] * this->nVertLevels;
        int c2offset = cellIndices[1] * this->nVertLevels;
        int c3offset = cellIndices[2] * this->nVertLevels;
        VECTOR3 w = this->cov_weights[localVertIdx];
        for (int vl = 0; vl < this->nVertLevels; vl++) {
            vVelocityOut[vl].Zero();
            for (int dim = 0; dim < 3; dim++) {
                pVelocityOut[vl][dim] = w[0]*cellVelocity[c1offset+vl][dim]
                                       + w[1]*cellVelocity[c2offset+vl][dim]
                                       + w[2]*cellVelocity[c3offset+vl][dim];
            }
            vVelocityOut[vl][0]   = w[0]*cellVertVelocity[c1offset+vl]
                                   + w[1]*cellVertVelocity[c2offset+vl]
                                   + w[2]*cellVertVelocity[c3offset+vl];
            ztopOut[vl]           = w[0]*cellZTop[c1offset+vl]
                                   + w[1]*cellZTop[c2offset+vl]
                                   + w[2]*cellZTop[c3offset+vl];
        }
    }
}
