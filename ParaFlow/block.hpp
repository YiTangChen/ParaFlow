#ifndef _BLOCK_HPP
#define _BLOCK_HPP

#include <iostream>
#include <fstream>
#include <OSUFlow.h>
#include <LoadTiming.h>
#include <diy/master.hpp>
#include "utils.hpp"
#include "timing.hpp"
#ifdef OSUFLOW_ENABLE_CUDA
#include "GPU/CUDA/MPASOGPUTracer.h"
#endif

using namespace std;

typedef     diy::RegularGridLink        RGLink;

// the block structure
struct Block
{
    int rank;
    OSUFlow *osuflow;
    std::vector<VECTOR3> currentSeeds;
    std::vector<PtInfo> startPts;
    std::list<vtListSeedTrace *> sl_list;
    std::list<vtListTimeSeedTrace *> pl_list;
    std::vector<Segment> segs;
    std::vector<int> sl_stepcounts;  // actual integration steps per sl_list entry
    std::vector<int> pl_stepcounts;  // actual integration steps per pl_list entry
    std::vector<int> neighborIndices;
    std::vector<int> areaID2neighborID;
    std::vector<int> fromCells;
    std::vector<int> toCells;
    std::vector<PtInfo> temporalBoundaryPts;
    int* areaIndicesArr;
    int* LocalCell2GlobalCell;
    int* GlobalCell2LocalCell;
    int nLocalCells;
    int nGlobalCells;
    int nVertLevels;
    int areaId;
    double integrationDt;
    BlockTiming timing;
#ifdef OSUFLOW_ENABLE_CUDA
    mpaso_gpu_host::GPUBlockContext* gpuContext;
#endif

    // following is mandatory
    Block() : osuflow(nullptr), areaIndicesArr(nullptr), LocalCell2GlobalCell(nullptr), GlobalCell2LocalCell(nullptr)
#ifdef OSUFLOW_ENABLE_CUDA
            , gpuContext(nullptr)
#endif
    {}
    ~Block() {
#ifdef OSUFLOW_ENABLE_CUDA
        mpaso_gpu_host::DestroyGPUBlockContext(gpuContext);
        gpuContext = nullptr;
#endif
        delete[] areaIndicesArr;
        delete[] LocalCell2GlobalCell;
        delete[] GlobalCell2LocalCell;
        delete osuflow;
    }
    void            set_rank(int r)     { rank = r; }
    void            set_areaId(int g)      { this->areaId = g; }
    int             get_rank()          { return rank; }
    static void*    create()            { return new Block; }
    static void     destroy(void* b)    { delete static_cast<Block*>(b); }

    void set_data(int gid, const char* meshfile, const char* datafile, std::vector<VECTOR3>& allseeds, std::vector<int>& areaIndices, TimeVaryingDataConfig &cfg) {
        this->integrationDt = cfg.dt.value_or(1.0);
        int* areaIdxPointer = &(areaIndices[0]);
        areaIndicesArr = new int[areaIndices.size()];
        for(int i = 0; i < areaIndices.size(); i++)
            areaIndicesArr[i] = areaIndices[i];
        this->LocalCell2GlobalCell = nullptr;
        this->GlobalCell2LocalCell = nullptr;

        osuflow = new OSUFlow();
        if (meshfile == NULL || *meshfile == '\0') {
            std::cout << "Loading data without mesh file for gid: " << gid << std::endl;
            osuflow->LoadMPASOData(datafile, gid, areaIndices.size(), areaIdxPointer, this->neighborIndices, 
                                this->LocalCell2GlobalCell, this->nLocalCells, this->nVertLevels, this->GlobalCell2LocalCell, this->nGlobalCells, cfg);
        } else {
            std::cout << "Loading data with mesh file for gid: " << gid << std::endl;
            osuflow->LoadMPASOData(meshfile, datafile, gid, areaIndices.size(), areaIdxPointer, this->neighborIndices, 
                                this->LocalCell2GlobalCell, this->nLocalCells, this->nVertLevels, this->GlobalCell2LocalCell, this->nGlobalCells, cfg);
        }

        if (!this->neighborIndices.empty()) {
            int maxneighboridx = *std::max_element(this->neighborIndices.begin(), this->neighborIndices.end());
            this->areaID2neighborID.assign(maxneighboridx + 1, -1);
            for(int i = 0 ; i < (int)this->neighborIndices.size(); i++)
                this->areaID2neighborID[this->neighborIndices[i]] = i;
        }

        std::cerr << "[Block::set_data]: areaID2neighborID size = " << areaID2neighborID.size() << std::endl;
        assert(osuflow != nullptr && "osuflow should not be null!");

        if (allseeds.empty()) {
            std::cerr << "[Block::set_data] warning: allseeds is empty, skip\n";
            return;
        }


        currentSeeds.reserve(allseeds.size());
        startPts.reserve(allseeds.size());
        fromCells.reserve(allseeds.size());

        load_timing::clock::time_point _t_sf;
        if (load_timing::g.enabled) _t_sf = load_timing::tic();
        int seedCnt = 0;
        int checkpt = (allseeds.size() >= 10) ? (int)(allseeds.size() / 10) : 0;
        for (auto& seed : allseeds) {
            assert(seed.Dimension() == 3 && "VECTOR3 must have 3 dimensions");

            if (osuflow->inBlock(seed)) {
                VECTOR3 tmp;
                for (int i = 0; i < seed.Dimension(); ++i) {
                    tmp[i] = seed[i];
                }
                currentSeeds.push_back(tmp);

                PtInfo pt;
                for (int i = 0; i < seed.Dimension(); ++i) {
                    pt.coord[i] = seed[i];
                }
                pt.pid      = seedCnt;
                pt.gid      = gid;
                pt.sid      = 0;
                pt.nsteps   = 0;
                pt.fromCell = -1;
                pt.ts = 0.0;
                startPts.push_back(pt);

                fromCells.push_back(-1);
            }
            ++seedCnt;
            if(checkpt > 0 && seedCnt % checkpt == 0) {
                std::cerr << "[Block::set_data] gid: " << gid << ", processed " << seedCnt << "/" << allseeds.size() << " seeds" << std::endl;
            }
        }
        if (load_timing::g.enabled) load_timing::g.seed_filter_s += load_timing::toc(_t_sf);

        this->setStartPtsDone();
    }

    void set_data_with_ghost_cells(int gid, const char* meshfile, const char* datafile, std::vector<VECTOR3>& allseeds, std::vector<int>& areaIndices, TimeVaryingDataConfig &cfg) {
        this->integrationDt = cfg.dt.value_or(1.0);
        int* areaIdxPointer = &(areaIndices[0]);
        areaIndicesArr = new int[areaIndices.size()];
        for(int i = 0; i < areaIndices.size(); i++)
            areaIndicesArr[i] = areaIndices[i];
        this->LocalCell2GlobalCell = nullptr;
        this->GlobalCell2LocalCell = nullptr;

        osuflow = new OSUFlow();
        if (meshfile == NULL || *meshfile == '\0') {
            std::cout << "Loading data without mesh file for gid: " << gid << std::endl;
            osuflow->LoadMPASOData(datafile, gid, areaIndices.size(), areaIdxPointer, this->neighborIndices, 
                                this->LocalCell2GlobalCell, this->nLocalCells, this->nVertLevels, this->GlobalCell2LocalCell, this->nGlobalCells, cfg);
        } else {
            std::cout << "Loading data with mesh file for gid: " << gid << std::endl;
            osuflow->LoadMPASOData(meshfile, datafile, gid, areaIndices.size(), areaIdxPointer, this->neighborIndices, 
                                this->LocalCell2GlobalCell, this->nLocalCells, this->nVertLevels, this->GlobalCell2LocalCell, this->nGlobalCells, cfg);
        }

        if (!this->neighborIndices.empty()) {
            int maxneighboridx = *std::max_element(this->neighborIndices.begin(), this->neighborIndices.end());
            this->areaID2neighborID.assign(maxneighboridx + 1, -1);
            for(int i = 0 ; i < (int)this->neighborIndices.size(); i++)
                this->areaID2neighborID[this->neighborIndices[i]] = i;
        }

        std::cerr << "[Block::set_data]: areaID2neighborID size = " << areaID2neighborID.size() << std::endl;
        assert(osuflow != nullptr && "osuflow should not be null!");

        if (allseeds.empty()) {
            std::cerr << "[Block::set_data] warning: allseeds is empty, skip\n";
            return;
        }


        currentSeeds.reserve(allseeds.size());
        startPts.reserve(allseeds.size());
        fromCells.reserve(allseeds.size());

        int seedCnt = 0;
        int checkpt = (allseeds.size() >= 10) ? (int)(allseeds.size() / 10) : 0;
        for (auto& seed : allseeds) {
            assert(seed.Dimension() == 3 && "VECTOR3 must have 3 dimensions");

            if (osuflow->inBlock_WithGhost(seed)) {
                VECTOR3 tmp;
                for (int i = 0; i < seed.Dimension(); ++i) {
                    tmp[i] = seed[i];
                }
                currentSeeds.push_back(tmp);

                PtInfo pt;
                for (int i = 0; i < seed.Dimension(); ++i) {
                    pt.coord[i] = seed[i];
                }
                pt.pid      = seedCnt;
                pt.gid      = gid;
                pt.sid      = 0;
                pt.nsteps   = 0;
                pt.fromCell = -1;
                pt.ts = 0.0;
                startPts.push_back(pt);

                fromCells.push_back(-1);
            }
            ++seedCnt;
            if(checkpt > 0 && seedCnt % checkpt == 0) {
                std::cerr << "[Block::set_data] gid: " << gid << ", processed " << seedCnt << "/" << allseeds.size() << " seeds" << std::endl;
            }
        }

        this->setStartPtsDone();
    }

    void dumpData(string &prefix) {
        osuflow->dumpData(prefix);
    }

    void endTracing() {
        this->fromCells.clear();
        this->toCells.clear();
        this->startPts.clear();
        this->currentSeeds.clear();
        // temporalBoundaryPts is intentionally NOT cleared here;
        // GenPathLines re-injects them after advanceTimestep()
    }

    void addStartPts(PtInfo& pt) {
        this->startPts.push_back(pt);
        VECTOR3 tmp;
        for(int i = 0; i < 3; i++)
            tmp[i] = pt.coord[i];
        this->currentSeeds.push_back(tmp);
        // pt.fromCell is a global flat index; convert to local for phys_to_cell hint
        this->fromCells.push_back(GetLocalCellidx(pt.fromCell));
    }

    void setStartPtsDone() {
        assert(osuflow != nullptr);
        if (this->currentSeeds.empty()) return;  // no seeds in this block, skip
        VECTOR3* seedsPointer = this->currentSeeds.data();
        // VECTOR3* seedsPointer = &(this->currentSeeds[0]);
        osuflow->SetSeedPoints(seedsPointer, this->currentSeeds.size());
    }

    void GenStreamLineByOSUFlow(int maxSteps, int saveInterval = 1,
                                const std::vector<int>& perSeedMaxPoints = std::vector<int>()) {
        // Free previous trace objects to prevent memory leak:
        // each vtListSeedTrace* and its VECTOR3* elements are heap-allocated by OSUFlow
        // and sl_list does not own them; we must delete them before clearing.
        for (auto* trace : this->sl_list) {
            for (auto* pt : *trace)
                delete pt;
            delete trace;
        }
        this->sl_list.clear();
        osuflow->SetIntegrationParams(this->integrationDt, this->integrationDt);
        osuflow->SetIntegrationOrder(MPASO_FOURTH);
        osuflow->SetSaveInterval(saveInterval);
        osuflow->SetPerSeedMaxPoints(perSeedMaxPoints);
        osuflow->GenStreamLines(this->sl_list, this->fromCells, this->toCells, FORWARD_DIR, maxSteps, 0);
        this->sl_stepcounts = osuflow->GetLastTraceStepCounts();
        // NOTE: segs are built by trace_block in ParaFlow.cpp, not here
    }

    void GenPathLineByOSUFlow(int maxSteps, double* tarray, int saveInterval = 1) {
        // Free previous pathline trace objects (each VECTOR4* is heap-allocated by OSUFlow)
        for (auto* trace : this->pl_list) {
            for (auto* pt : *trace)
                delete pt;
            delete trace;
        }
        this->pl_list.clear();
        this->toCells.clear();
        osuflow->SetIntegrationParams(this->integrationDt, this->integrationDt);
        osuflow->SetIntegrationOrder(MPASO_FOURTH);
        osuflow->SetSaveInterval(saveInterval);
        osuflow->GenPathLines(this->currentSeeds.data(), this->pl_list, FORWARD,
                              this->fromCells, this->toCells,
                              (int)this->currentSeeds.size(), maxSteps, tarray);
        this->pl_stepcounts = osuflow->GetLastTraceStepCounts();
        // NOTE: segs are built by trace_block_pathline in ParaFlow.cpp, not here
    }

    int getTotalTimesteps() {
        assert(osuflow);
        return osuflow->GetTotalTimesteps();
    }

    int getGlobalTimestepOffset() {
        assert(osuflow);
        return osuflow->GetTimestepOffset();
    }

    double getWindowMaxTime() {
        assert(osuflow && osuflow->GetFlowField());
        return (double)osuflow->GetFlowField()->GetMaxTimeStep();
    }

    void advanceTimestep() {
        assert(osuflow);
        osuflow->AdvanceTimestep();
    }

    // cellidx is a local flat index (localCellIdx * nVertLevels + vLevel)
    // as stored in toCells by phys_to_cell
    int GetNeighborId(int cellidx) {
        if (cellidx < 0) return -1;  // dead/invalid cell
        int localCellId = cellidx / this->nVertLevels;
        if (localCellId < 0 || localCellId >= this->nLocalCells) return -1;
        int globalCellId = this->LocalCell2GlobalCell[localCellId] - 1;
        if (globalCellId < 0 || globalCellId >= this->nGlobalCells) return -1;
        int partitionId  = this->areaIndicesArr[globalCellId];
        if (partitionId < 0 || partitionId >= (int)this->areaID2neighborID.size())
            return -1;
        return this->areaID2neighborID[partitionId];
    }

    // Convert local flat index → global flat index (for cross-block communication)
    int GetGlobalCellidx(int cellidx) {
        if (cellidx < 0) return -1;
        int vLevel = cellidx % this->nVertLevels;
        int localCellId = cellidx / this->nVertLevels;
        if (localCellId < 0 || localCellId >= this->nLocalCells) return -1;
        int globalCellId = this->LocalCell2GlobalCell[localCellId] - 1;
        return globalCellId * this->nVertLevels + vLevel;
    }

    // Convert global flat index → local flat index (for at_phys hint)
    int GetLocalCellidx(int cellidx) {
        if (cellidx < 0) return -1;
        int vLevel = cellidx % this->nVertLevels;
        int globalCellId = cellidx / this->nVertLevels;
        int localCellId = this->GlobalCell2LocalCell[globalCellId] - 1;
        if (localCellId < 0) return -1;
        return localCellId * this->nVertLevels + vLevel;
    }

    // writeToDisk=false skips the file entirely -- output still gets every segment,
    // just with no bytes ever touching storage (used by streamline_storage: memory).
    void write_trajectory(const std::string& filename, std::list<std::vector<VECTOR3>>& output, int interval = 1, bool writeToDisk = true) {
        std::ofstream outfile;
        if (writeToDisk) outfile.open(filename, ios::out | ios::binary);
        for(int segidx = 0; segidx < (int)this->segs.size(); segidx++) {
            int pid  = segs[segidx].pid;
            int gid  = segs[segidx].gid;
            int sid  = segs[segidx].sid;
            int nPts = (int)segs[segidx].coords.size();
            // coords are already downsampled by interval during tracing
            if (writeToDisk) {
                outfile.write((char *) &pid,  sizeof(int));
                outfile.write((char *) &gid,  sizeof(int));
                outfile.write((char *) &sid,  sizeof(int));
                outfile.write((char *) &nPts, sizeof(int));
            }
            std::vector<VECTOR3> seg2vec;
            seg2vec.reserve(nPts);
            for(int ptidx = 0; ptidx < nPts; ptidx++) {
                seg2vec.push_back(segs[segidx].coords[ptidx]);
                if (writeToDisk) {
                    outfile.write((char *) &(segs[segidx].coords[ptidx][0]), sizeof(double));
                    outfile.write((char *) &(segs[segidx].coords[ptidx][1]), sizeof(double));
                    outfile.write((char *) &(segs[segidx].coords[ptidx][2]), sizeof(double));
                }
            }
            output.push_back(std::move(seg2vec));
        }
        if (writeToDisk) outfile.close();
    }

    // Serialize this block's segments into the SAME wire format write_trajectory
    // writes to a file (repeated [pid][gid][sid][nPts][xyz]*nPts), but appended
    // to an in-memory buffer instead. Used by LIC combine mode
    // (rendering/LIC/lic_render.hpp lic_gather_segments) to MPI-gather segments
    // straight to rank 0 without any rank ever touching disk.
    void append_trajectory_bytes(std::vector<char>& buf) const {
        auto put = [&buf](const void* p, size_t n) {
            const char* c = (const char*)p;
            buf.insert(buf.end(), c, c + n);
        };
        for (size_t segidx = 0; segidx < segs.size(); ++segidx) {
            const int pid  = segs[segidx].pid;
            const int gid  = segs[segidx].gid;
            const int sid  = segs[segidx].sid;
            const int nPts = (int)segs[segidx].coords.size();
            put(&pid, sizeof(int));
            put(&gid, sizeof(int));
            put(&sid, sizeof(int));
            put(&nPts, sizeof(int));
            for (int ptidx = 0; ptidx < nPts; ++ptidx) {
                put(&(segs[segidx].coords[ptidx][0]), sizeof(double));
                put(&(segs[segidx].coords[ptidx][1]), sizeof(double));
                put(&(segs[segidx].coords[ptidx][2]), sizeof(double));
            }
        }
    }
};

#endif
