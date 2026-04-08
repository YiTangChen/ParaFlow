#include "ParaFlow.hpp"

void ParaFlow::ReadRegularGridXYZ(const char* regularField, int &domain_x, int &domain_y, int &domain_z) 
{
    std::ifstream myFile;
    myFile.open(regularField, std::ifstream::binary);
    if(myFile) {
        myFile.read(reinterpret_cast<char *>(&domain_x), sizeof(int));
        myFile.read(reinterpret_cast<char *>(&domain_y), sizeof(int));
        myFile.read(reinterpret_cast<char *>(&domain_z), sizeof(int));
        std::cout << "x = "<< domain_x << ", y = " << domain_y << ", z = " << domain_z << std::endl;
        myFile.close();
    }
    else {
        throw "File doesn't exist!";
    }
}

void ParaFlow::ReadSeedFile(const char* seedFile, std::vector<VECTOR3>& seeds) 
{
    std::ifstream myFile;
    std::string line;
    std::string delimiter = ",";
    int maxstep;
    myFile.open(seedFile);
    while(getline(myFile, line)) {
        VECTOR3 tmp;
        int pos = line.find(delimiter);
        tmp[0] = stof(line.substr(0, pos));
        line.erase(0, pos + delimiter.length());
        pos = line.find(delimiter);
        tmp[1] = stof(line.substr(0, pos));
        line.erase(0, pos + delimiter.length());
        tmp[2] = stof(line.substr(0));
        seeds.push_back(tmp);
        // for only load 1 seed
        // break;
    }
    myFile.close();
}

void ParaFlow::ReadSeedBinFile(const char* seedFile, std::vector<VECTOR3>& seeds) 
{
    double x,y,z;
    std::ifstream inFile(seedFile, std::ios::in | std::ios::binary);
    while(1) {
        if(inFile.read(reinterpret_cast<char*>(&x), sizeof(x))) {
            inFile.read(reinterpret_cast<char*>(&y), sizeof(y));
            inFile.read(reinterpret_cast<char*>(&z), sizeof(z));
            VECTOR3 tmp(x, y, z);
            seeds.push_back(tmp);
        }
        else
            break;
    }
    inFile.close();
}

void ParaFlow::ReadAreaIndices(const char* graph_partition_file)
{
    ifstream f(graph_partition_file, ios::out);
    int areaIdx;
    while(f >> areaIdx) {
        this->areaIndices.push_back(areaIdx);
    }
    f.close();
}


ParaFlow::ParaFlow(int argc, char* argv[], const char* configFile)
    : env(argc, argv), world()
{
    YAML::Node config = YAML::LoadFile(configFile);

    // Read configs from YAML
    // Read timing flag first so it applies to seed I/O below
    this->enable_timing = config["enable_timing"] ? config["enable_timing"].as<bool>() : false;

    try {
        if(config["nproc"])
            this->size = config["nproc"].as<int>();
        else
            throw "The number of processes isn't specified correctly.\n";
        if(config["nblocks"])
            this->nblocks = config["nblocks"].as<int>();
        else
            throw "The number of nblocks isn't specified correctly.\n";
        if(config["maxsteps"])
            this->max_steps = config["maxsteps"].as<int>();
        else
            throw "The number of steps isn't specified correctly.\n";
        this->interval = config["record_interval"] ? config["record_interval"].as<int>() : 1;
        if(this->interval < 1) this->interval = 1;
        if(!config["datafile_src"])
            throw "Error: no vector file(s).\n";
    }
    catch(const char* msg) {
        std::cout << msg;
    }

    this->outputDir = (config["outputDir"] && !config["outputDir"].IsNull())
                      ? config["outputDir"].as<std::string>() : "output";
    if (world.rank() == 0)
        std::filesystem::create_directories(this->outputDir);
    world.barrier();  // all ranks wait until directory exists before any writes

    if(config["graph_partition_indices"]) {
        std::string gp_src = config["graph_partition_indices"].as<string>();
        this->ReadAreaIndices(gp_src.c_str());
    }

    if(config["meshfile_src"] && !config["meshfile_src"].IsNull())
        this->meshFile = config["meshfile_src"].as<string>();
    else
        this->meshFile = "";
    if (config["datafile_src"].IsSequence()) {
        for (const auto& node : config["datafile_src"])
            this->vectorFiles.push_back(node.as<string>());
    } else {
        this->vectorFiles.push_back(config["datafile_src"].as<string>());
    }

    if(config["seedFile"]) {
        this->seedFile = config["seedFile"].as<string>();
        if(config["isBinarySeed"].as<bool>()) 
            seedReadMode = SeedReadMode::Binary;
        else 
            seedReadMode = SeedReadMode::Text;
    }

    // Time-varying variables
    cfg.velocity.cartesian.X  = read_opt_str(config["velocity"]["cartesian"], "X");
    cfg.velocity.cartesian.Y  = read_opt_str(config["velocity"]["cartesian"], "Y");
    cfg.velocity.cartesian.Z  = read_opt_str(config["velocity"]["cartesian"], "Z");
    cfg.velocity.spherical.zonal      = read_opt_str(config["velocity"]["spherical"], "zonal");
    cfg.velocity.spherical.meridional = read_opt_str(config["velocity"]["spherical"], "meridional");
    cfg.velocity.normal = read_opt_str(config["velocity"], "normal");
    cfg.zTop            = read_opt_str(config, "zTop");
    cfg.layerThickness  = read_opt_str(config, "layerThickness");
    cfg.vertVelocityTop = read_opt_str(config, "vertVelocityTop");
    cfg.xtime           = read_opt_str(config, "xtime");
    cfg.isTimeVarying   = config["isTimeVarying"] ? config["isTimeVarying"].as<bool>() : false;
    if (cfg.isTimeVarying) {
        cfg.loadNTimeSteps = config["loadNTimeSteps"] ? config["loadNTimeSteps"].as<int>() : 1;
    }
    cfg.dt = config["dt"] ? std::optional<double>(config["dt"].as<double>()) : std::nullopt;
    cfg.dataFiles = this->vectorFiles;  // pass full file list so MPASOReader can switch files

    // Read seeds on rank 0, broadcast to all — done once here for all Gen* calls
    double t_seed_read = 0.0;
    if (world.rank() == 0) {
        double _t = pf_now(enable_timing);
        if (seedReadMode == SeedReadMode::Binary)
            ReadSeedBinFile(seedFile.c_str(), seeds);
        else if (seedReadMode == SeedReadMode::Text)
            ReadSeedFile(seedFile.c_str(), seeds);
        pf_accum(t_seed_read, _t, enable_timing);
    }

    double t_seed_bcast = 0.0;
    double _t_bcast = pf_now(enable_timing);
    diy::mpi::broadcast(world, seeds, 0);
    pf_accum(t_seed_bcast, _t_bcast, enable_timing);

    if (enable_timing) {
        if (world.rank() == 0)
            fprintf(stderr, "TIMING phase=seed_read  rank=0 n_seeds=%zu t=%.6f\n",
                    seeds.size(), t_seed_read);
        fprintf(stderr, "TIMING phase=seed_bcast rank=%d t=%.6f\n",
                world.rank(), t_seed_bcast);
    }
}

ParaFlow::~ParaFlow() {

}

void ParaFlow::deq_incoming_iexchange(Block* b, const diy::Master::ProxyWithLink& cp)
{
    diy::Link *l = cp.link();
    int received = 0;
    for (size_t i = 0; i < l->size(); ++i)
    {
        int nbr_gid = l->target(i).gid;
        while (cp.incoming(nbr_gid))
        {
            PtInfo incoming_pt;
            cp.dequeue(nbr_gid, incoming_pt);
            b->addStartPts(incoming_pt);
            received++;
        }
        b->setStartPtsDone();
    }
    if (enable_timing) b->timing.n_particles_received += received;
}

void ParaFlow::trace_block(Block*                               b,
                            const diy::Master::ProxyWithLink&   cp)
{
    diy::Link* l = cp.link();

    do
    {
        double _t_comm = pf_now(enable_timing);
        this->deq_incoming_iexchange(b, cp);
        pf_accum(b->timing.t_trace_comm, _t_comm, enable_timing);

        if(b->currentSeeds.size() == 0)
            continue;
        // std::cerr << "[rank " << b->rank << ", gid " << cp.gid() << "] Tracing " << b->currentSeeds.size() << " seeds\n";
        int maxRemaining = 0;
        for (auto& pt : b->startPts)
            maxRemaining = std::max(maxRemaining, this->max_steps - pt.nsteps);

        double _t_comp = pf_now(enable_timing);
        b->GenStreamLineByOSUFlow(maxRemaining, this->interval);
        pf_accum(b->timing.t_trace_compute, _t_comp, enable_timing);
        if (enable_timing)
            for (auto s : b->sl_stepcounts) b->timing.n_steps_total += s;
        // deal with current seeds
        int sent = 0;
        int seedCnt = 0;
        for(auto sl = b->sl_list.begin(); sl != b->sl_list.end(); sl++) {
            vtListSeedTrace * trace = *sl;
            bool finished = false;

            Segment currseg;
            currseg.gid = b->startPts[seedCnt].gid;
            currseg.pid = b->startPts[seedCnt].pid;
            currseg.sid = b->startPts[seedCnt].sid;
            currseg.nsteps = b->startPts[seedCnt].nsteps;
            // Collect already-downsampled coords from OSUFlow
            for(auto pIt = trace->begin(); pIt != trace->end(); pIt++) {
                VECTOR3 tmp;
                for (size_t i = 0; i < 3; i++)
                    tmp[i] = (**pIt)[i];
                currseg.coords.push_back(tmp);
            }
            // Use actual step count reported by OSUFlow (interval-aware)
            if (seedCnt < (int)b->sl_stepcounts.size())
                currseg.nsteps += b->sl_stepcounts[seedCnt];

            if(currseg.nsteps >= this->max_steps)
                finished = true;

            // skip seeds that were out of boundary from the start
            if(currseg.coords.empty()) {
                seedCnt++;
                continue;
            }

            int nextNeighborId = b->GetNeighborId(b->toCells[seedCnt]);

            b->segs.push_back(currseg);

            if(!finished) {
                PtInfo endPt;
                for(int i = 0; i < 3; i++)
                    endPt.coord[i] = currseg.coords.back()[i];
                endPt.gid = currseg.gid;
                endPt.pid = currseg.pid;
                endPt.sid = currseg.sid + 1;
                endPt.nsteps = currseg.nsteps;
                endPt.fromCell = b->GetGlobalCellidx(b->toCells[seedCnt]);

                if(nextNeighborId >= 0) {
                    diy::BlockID bid = l->target(nextNeighborId); // in case of multiple dests, send to first dest only
                    cp.enqueue(bid, endPt);
                    sent++;
                }
            }
            seedCnt++;
        }
        b->endTracing();
    } while (cp.fill_incoming());
}

bool ParaFlow::trace_block_iexchange(Block*                              b,
                                    const diy::Master::ProxyWithLink&    cp)
{
    this->trace_block(b, cp);
    return true;
}

void ParaFlow::GenStreamLines(std::list<std::vector<VECTOR3>>& streamlines)
{
    diy::ContiguousAssigner   assigner(this->size, this->nblocks);
    int rank = world.rank();
    diy::Master master(world,
                       1,                              // one thread
                       -1,                             // all blocks in memory
                       &Block::create,
                       &Block::destroy);

    std::vector<int> gids;                     // global ids of local blocks
    assigner.local_gids(world.rank(), gids);   // get the gids of local blocks
    for (size_t i = 0; i < gids.size(); ++i)   // for the local blocks in this processor
    {
        int gid = gids[i];

        Block *b = new Block;
        diy::Link*   link = new diy::Link;   // link is this block's neighborhood
        b->set_rank(rank);
        if (enable_timing) b->timing.mem_vmrss_before_kb = pf_read_vmrss_kb();
        double _t_load = pf_now(enable_timing);
        b->set_data(gid, this->meshFile.c_str(), this->vectorFiles[0].c_str(), this->seeds, this->areaIndices, this->cfg);
        pf_accum(b->timing.t_block_load, _t_load, enable_timing);
        if (enable_timing) {
            b->timing.n_seeds_initial      = (int)b->currentSeeds.size();
            b->timing.mem_vmrss_after_kb   = pf_read_vmrss_kb();
            MPASOGrid* grid = dynamic_cast<MPASOGrid*>(b->osuflow->GetFlowField()->GetGrid());
            if (grid) {
                b->timing.mem_grid_bytes     = grid->getGridMemBytes();
                b->timing.mem_solution_bytes = grid->getSolutionMemBytes();
                // Add Block-owned index arrays (not visible to MPASOGrid)
                // b->timing.mem_grid_bytes    += (size_t)b->nLocalCells  * sizeof(int); // LocalCell2GlobalCell
                // b->timing.mem_grid_bytes    += (size_t)b->nGlobalCells * sizeof(int); // GlobalCell2LocalCell
                // b->timing.mem_grid_bytes    += (size_t)b->nGlobalCells * sizeof(int); // areaIndicesArr
            }
        }

        diy::BlockID neighbor;               // one neighbor in the neighborhood
        for(int neighborIdx: b->neighborIndices) {
            neighbor.gid  = neighborIdx;                    // gid of the neighbor block
            neighbor.proc = assigner.rank(neighbor.gid);    // process of the neighbor block
            link->add_neighbor(neighbor);
        }
        master.add(gid, b, link); // add block to the master (mandatory)
    }

    double t_iexchange = 0.0;
    double _t_iex = pf_now(enable_timing);
    int ncalls = 0;
    master.iexchange([&](Block* b, const diy::Master::ProxyWithLink& icp) -> bool
    {
        ncalls++;
        bool val = this->trace_block_iexchange(b, icp);
        return val;
    });
    world.barrier();
    pf_accum(t_iexchange, _t_iex, enable_timing);

    // Write output — always runs; timing is gated inside
    for (size_t i = 0; i < gids.size(); ++i)   // for the local blocks in this processor
    {
        int gid = gids[i];
        std::string filename = this->outputDir + "/" + std::to_string(gid) + ".bin";
        double _t_write = pf_now(enable_timing);
        Block* blk = (Block*)master.block(master.lid(gid));
        blk->write_trajectory(filename, streamlines, this->interval);
        pf_accum(blk->timing.t_output_write, _t_write, enable_timing);
    }

    // All measurement output — fully gated: zero overhead when enable_timing=false
    if (enable_timing) {
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp) {
            b->timing.mem_peak_vmhwm_kb = pf_read_vmhwm_kb();
            fprintf(stderr,
                "TIMING phase=block_load    rank=%d gid=%d t=%.6f nseeds=%d\n",
                b->rank, cp.gid(), b->timing.t_block_load, b->timing.n_seeds_initial);
            fprintf(stderr,
                "TIMING phase=trace_compute rank=%d gid=%d t=%.6f nsteps=%ld nrecv=%d\n",
                b->rank, cp.gid(), b->timing.t_trace_compute,
                b->timing.n_steps_total, b->timing.n_particles_received);
            fprintf(stderr,
                "TIMING phase=trace_comm    rank=%d gid=%d t=%.6f\n",
                b->rank, cp.gid(), b->timing.t_trace_comm);
            fprintf(stderr,
                "TIMING phase=output_write  rank=%d gid=%d t=%.6f\n",
                b->rank, cp.gid(), b->timing.t_output_write);
            fprintf(stderr,
                "MEM_ANALYTICAL rank=%d gid=%d grid_bytes=%zu solution_bytes=%zu total_bytes=%zu\n",
                b->rank, cp.gid(),
                b->timing.mem_grid_bytes, b->timing.mem_solution_bytes,
                b->timing.mem_grid_bytes + b->timing.mem_solution_bytes);
            fprintf(stderr,
                "MEM_DELTA      rank=%d gid=%d before_kb=%ld after_kb=%ld delta_kb=%ld\n",
                b->rank, cp.gid(),
                b->timing.mem_vmrss_before_kb, b->timing.mem_vmrss_after_kb,
                b->timing.mem_vmrss_after_kb - b->timing.mem_vmrss_before_kb);
            fprintf(stderr,
                "MEM_PEAK       rank=%d gid=%d vmhwm_kb=%ld\n",
                b->rank, cp.gid(), b->timing.mem_peak_vmhwm_kb);
        });
        fprintf(stderr, "TIMING phase=iexchange_total rank=%d t=%.6f\n", rank, t_iexchange);
        double iex_min, iex_max, iex_sum;
        MPI_Reduce(&t_iexchange, &iex_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_iexchange, &iex_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_iexchange, &iex_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
            fprintf(stderr,
                "TIMING phase=iexchange_global iex_min=%.6f iex_max=%.6f iex_avg=%.6f iex_imbalance=%.3f\n",
                iex_min, iex_max, iex_sum / this->size,
                (iex_min > 0.0) ? iex_max / iex_min : 0.0);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PathLine support
// ─────────────────────────────────────────────────────────────────────────────

void ParaFlow::trace_block_pathline(Block*                              b,
                                    const diy::Master::ProxyWithLink&   cp)
{
    diy::Link* l = cp.link();

    do
    {
        double _t_comm = pf_now(enable_timing);
        this->deq_incoming_iexchange(b, cp);
        pf_accum(b->timing.t_trace_comm, _t_comm, enable_timing);

        if (b->currentSeeds.size() == 0)
            continue;

        // Build per-seed start times from PtInfo.ts (0.0 for initial seeds,
        // non-zero for particles that were forwarded from a neighbour block).
        std::vector<double> tarray(b->startPts.size());
        for (size_t i = 0; i < b->startPts.size(); i++)
            tarray[i] = b->startPts[i].ts;

        // ── Advect particles through the time-varying field ──────────────────
        int maxRemaining = 0;
        for (auto& pt : b->startPts)
            maxRemaining = std::max(maxRemaining, this->max_steps - pt.nsteps);

        double _t_comp = pf_now(enable_timing);
        b->GenPathLineByOSUFlow(maxRemaining, tarray.data(), this->interval);
        pf_accum(b->timing.t_trace_compute, _t_comp, enable_timing);
        if (enable_timing)
            for (auto s : b->pl_stepcounts) b->timing.n_steps_total += s;

        int sent = 0;
        int seedCnt = 0;
        double maxT = b->getWindowMaxTime();
        for (auto pl = b->pl_list.begin(); pl != b->pl_list.end(); pl++, seedCnt++)
        {
            vtListTimeSeedTrace* trace = *pl;

            Segment currseg;
            currseg.gid    = b->startPts[seedCnt].gid;
            currseg.pid    = b->startPts[seedCnt].pid;
            currseg.sid    = b->startPts[seedCnt].sid;
            currseg.nsteps = b->startPts[seedCnt].nsteps;

            // Each VECTOR4 entry: [x, y, z, t]; collect already-downsampled coords
            double exitTime = tarray[seedCnt];
            for (auto pIt = trace->begin(); pIt != trace->end(); pIt++)
            {
                VECTOR3 tmp;
                for (int i = 0; i < 3; i++)
                    tmp[i] = (**pIt)[i];
                currseg.coords.push_back(tmp);
                exitTime = (**pIt)[3];  // time of last stored point
            }
            // Use actual step count reported by OSUFlow (interval-aware)
            if (seedCnt < (int)b->pl_stepcounts.size())
                currseg.nsteps += b->pl_stepcounts[seedCnt];

            if (currseg.coords.empty())
                continue;

            bool finished = (currseg.nsteps >= this->max_steps);

            int toCell         = (seedCnt < (int)b->toCells.size()) ? b->toCells[seedCnt] : -1;
            int nextNeighborId = b->GetNeighborId(toCell);
            bool isTemporal    = (!finished && nextNeighborId < 0 && exitTime >= maxT - b->integrationDt);

            b->segs.push_back(currseg);

            if (!finished)
            {
                PtInfo endPt;
                for (int i = 0; i < 3; i++)
                    endPt.coord[i] = currseg.coords.back()[i];
                endPt.gid      = currseg.gid;
                endPt.pid      = currseg.pid;
                endPt.sid      = currseg.sid + 1;
                endPt.nsteps   = currseg.nsteps;
                endPt.ts = exitTime;

                if (isTemporal)
                {
                    endPt.fromCell = (toCell >= 0) ? b->GetGlobalCellidx(toCell) : -1;
                    b->temporalBoundaryPts.push_back(endPt);
                }
                else if (nextNeighborId >= 0)
                {
                    endPt.fromCell  = b->GetGlobalCellidx(toCell);
                    diy::BlockID bid = l->target(nextNeighborId);
                    cp.enqueue(bid, endPt);
                    sent++;
                }
            }
        }

        int temporal = (int)b->temporalBoundaryPts.size();
        b->endTracing();
    } while (cp.fill_incoming());
}

bool ParaFlow::trace_block_pathline_iexchange(Block*                              b,
                                               const diy::Master::ProxyWithLink&   cp)
{
    this->trace_block_pathline(b, cp);
    return true;
}

void ParaFlow::GenPathLines(std::list<std::vector<VECTOR3>>& pathlines)
{
    diy::ContiguousAssigner   assigner(this->size, this->nblocks);
    int rank = world.rank();
    diy::Master master(world,
                       1,                      // one thread
                       -1,                     // all blocks in memory
                       &Block::create,
                       &Block::destroy);

    std::vector<int> gids;
    assigner.local_gids(world.rank(), gids);
    for (size_t i = 0; i < gids.size(); ++i)
    {
        int gid = gids[i];

        Block*     b    = new Block;
        diy::Link* link = new diy::Link;
        b->set_rank(rank);
        if (enable_timing) b->timing.mem_vmrss_before_kb = pf_read_vmrss_kb();
        double _t_load = pf_now(enable_timing);
        b->set_data(gid, this->meshFile.c_str(), this->vectorFiles[0].c_str(),
                    this->seeds, this->areaIndices, this->cfg);
        pf_accum(b->timing.t_block_load, _t_load, enable_timing);
        if (enable_timing) {
            b->timing.n_seeds_initial      = (int)b->currentSeeds.size();
            b->timing.mem_vmrss_after_kb   = pf_read_vmrss_kb();
            MPASOGrid* grid = dynamic_cast<MPASOGrid*>(b->osuflow->GetFlowField()->GetGrid());
            if (grid) {
                b->timing.mem_grid_bytes     = grid->getGridMemBytes();
                b->timing.mem_solution_bytes = grid->getSolutionMemBytes();
            }
        }

        diy::BlockID neighbor;
        for (int neighborIdx : b->neighborIndices) {
            neighbor.gid  = neighborIdx;
            neighbor.proc = assigner.rank(neighbor.gid);
            link->add_neighbor(neighbor);
        }
        master.add(gid, b, link);
    }

    // Get total timesteps in file from first local block (all blocks share the same file)
    int totalTimestepsInFile = 0;
    if (!gids.empty())
        totalTimestepsInFile = ((Block*)master.block(master.lid(gids[0])))->getTotalTimesteps();
    diy::mpi::broadcast(world, totalTimestepsInFile, 0);

    double t_iexchange = 0.0;   // accumulated across all temporal windows
    int windowIdx = 0;
    while (true) {
        double _t_iex = pf_now(enable_timing);
        master.iexchange([&](Block* b, const diy::Master::ProxyWithLink& icp) -> bool {
            return this->trace_block_pathline_iexchange(b, icp);
        });
        world.barrier();
        double t_window = 0.0;
        pf_accum(t_window, _t_iex, enable_timing);
        t_iexchange += t_window;
        if (enable_timing)
            fprintf(stderr, "TIMING phase=window_iexchange rank=%d window=%d t=%.6f\n",
                    rank, windowIdx, t_window);
        windowIdx++;

        // Use the reader's actual global offset to determine if more data exists.
        // This correctly handles file-boundary clamping that would cause windowStart to drift.
        int globalOffset = 0;
        if (!gids.empty())
            globalOffset = ((Block*)master.block(master.lid(gids[0])))->getGlobalTimestepOffset();
        diy::mpi::broadcast(world, globalOffset, 0);
        if (globalOffset >= totalTimestepsInFile)
            break;

        // Load next window and re-inject particles that hit the temporal boundary
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp) {
            b->advanceTimestep();
            for (auto& pt : b->temporalBoundaryPts)
                b->addStartPts(pt);
            if (!b->temporalBoundaryPts.empty()) {
                b->setStartPtsDone();
                b->temporalBoundaryPts.clear();
            }
        });
        world.barrier();
    }

    // Write output — always runs; timing is gated inside
    for (size_t i = 0; i < gids.size(); ++i)
    {
        int gid = gids[i];
        std::string filename = this->outputDir + "/" + std::to_string(gid) + ".bin";
        double _t_write = pf_now(enable_timing);
        Block* blk = (Block*)master.block(master.lid(gid));
        blk->write_trajectory(filename, pathlines, this->interval);
        pf_accum(blk->timing.t_output_write, _t_write, enable_timing);
    }

    // All measurement output — fully gated: zero overhead when enable_timing=false
    if (enable_timing) {
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp) {
            b->timing.mem_peak_vmhwm_kb = pf_read_vmhwm_kb();
            fprintf(stderr,
                "TIMING phase=block_load    rank=%d gid=%d t=%.6f nseeds=%d\n",
                b->rank, cp.gid(), b->timing.t_block_load, b->timing.n_seeds_initial);
            fprintf(stderr,
                "TIMING phase=trace_compute rank=%d gid=%d t=%.6f nsteps=%ld nrecv=%d\n",
                b->rank, cp.gid(), b->timing.t_trace_compute,
                b->timing.n_steps_total, b->timing.n_particles_received);
            fprintf(stderr,
                "TIMING phase=trace_comm    rank=%d gid=%d t=%.6f\n",
                b->rank, cp.gid(), b->timing.t_trace_comm);
            fprintf(stderr,
                "TIMING phase=output_write  rank=%d gid=%d t=%.6f\n",
                b->rank, cp.gid(), b->timing.t_output_write);
            fprintf(stderr,
                "MEM_ANALYTICAL rank=%d gid=%d grid_bytes=%zu solution_bytes=%zu total_bytes=%zu\n",
                b->rank, cp.gid(),
                b->timing.mem_grid_bytes, b->timing.mem_solution_bytes,
                b->timing.mem_grid_bytes + b->timing.mem_solution_bytes);
            fprintf(stderr,
                "MEM_DELTA      rank=%d gid=%d before_kb=%ld after_kb=%ld delta_kb=%ld\n",
                b->rank, cp.gid(),
                b->timing.mem_vmrss_before_kb, b->timing.mem_vmrss_after_kb,
                b->timing.mem_vmrss_after_kb - b->timing.mem_vmrss_before_kb);
            fprintf(stderr,
                "MEM_PEAK       rank=%d gid=%d vmhwm_kb=%ld\n",
                b->rank, cp.gid(), b->timing.mem_peak_vmhwm_kb);
        });
        fprintf(stderr, "TIMING phase=iexchange_total rank=%d t=%.6f\n", rank, t_iexchange);
        double iex_min, iex_max, iex_sum;
        MPI_Reduce(&t_iexchange, &iex_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_iexchange, &iex_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_iexchange, &iex_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
            fprintf(stderr,
                "TIMING phase=iexchange_global iex_min=%.6f iex_max=%.6f iex_avg=%.6f iex_imbalance=%.3f\n",
                iex_min, iex_max, iex_sum / this->size,
                (iex_min > 0.0) ? iex_max / iex_min : 0.0);
    }
}

void merge_traces(void* b_, const diy::ReduceProxy& rp, const diy::RegularMergePartners&)
{
    Block* b = static_cast<Block*>(b_);

    // dequeue and merge
    for (unsigned i = 0; i < rp.in_link().size(); ++i)
    {
        int nbr_gid = rp.in_link().target(i).gid;
        if (nbr_gid == rp.gid())                    // skip self
            continue;

        std::vector<Segment> in_traces;
        rp.dequeue(nbr_gid, in_traces);

        // append in_traces to segments, leaving trajectories segmented and disorganized
        // eventually could sort into continuous long trajectories, but not necessary at this time
        b->segs.insert(b->segs.end(), in_traces.begin(), in_traces.end());
    }

    // enqueue
    if (rp.out_link().size())
    {
        int nbr_gid = rp.out_link().target(0).gid;  // for a merge, the out_link size is 1; ie, there is only one target
        if (nbr_gid != rp.gid())                    // skip self
            rp.enqueue(rp.out_link().target(0), b->segs);
    }
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///// For Draw Subdomain (more like a test function for now) //////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

void ParaFlow::CheckPointsInBlock(std::function<void(int gid, const std::vector<VECTOR3>&)> onBlock)
{
    diy::ContiguousAssigner   assigner(this->size, this->nblocks);
    int rank = world.rank();
    diy::Master master(world,
                       1,                              // one thread
                       -1,                             // all blocks in memory
                       &Block::create,
                       &Block::destroy);

    std::vector<int> gids;                     // global ids of local blocks
    assigner.local_gids(world.rank(), gids);   // get the gids of local blocks
    for (size_t i = 0; i < gids.size(); ++i)   // for the local blocks in this processor
    {
        int gid = gids[i];
        std::cout << "gid: " << gid << ", seeds size: " << this->seeds.size() << std::endl;

        Block *b = new Block;
        diy::Link*   link = new diy::Link;   // link is this block's neighborhood
        b->set_rank(rank);
        b->set_data_with_ghost_cells(gid, this->meshFile.c_str(), this->vectorFiles[0].c_str(), this->seeds, this->areaIndices, this->cfg);
        std::cout << "set data done for gid: " << gid << std::endl;

        diy::BlockID neighbor;               // one neighbor in the neighborhood
        for(int neighborIdx: b->neighborIndices) {
            neighbor.gid  = neighborIdx;                    // gid of the neighbor block
            neighbor.proc = assigner.rank(neighbor.gid);    // process of the neighbor block
            link->add_neighbor(neighbor);
        }
        master.add(gid, b, link); // add block to the master (mandatory)

        if (onBlock) {
            onBlock(gid, b->currentSeeds);
        } else {
            // fallback: write seeds to file
            std::string filename = "drawSubdomain/seeds_" + std::to_string(gid) + ".bin";
            std::ofstream outFile(filename, std::ios::out | std::ios::binary);
            for (size_t j = 0; j < b->currentSeeds.size(); j++) {
                VECTOR3 tmp = b->currentSeeds[j];
                outFile.write(reinterpret_cast<const char*>(&tmp[0]), sizeof(double));
                outFile.write(reinterpret_cast<const char*>(&tmp[1]), sizeof(double));
                outFile.write(reinterpret_cast<const char*>(&tmp[2]), sizeof(double));
            }
        }
    }

    master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp) {
        std::ifstream f("/proc/self/status");
        std::string line;
        while (std::getline(f, line)) {
            if (line.rfind("VmHWM:", 0) == 0) {
                fprintf(stderr, "MEM_PEAK rank=%d gid=%d %s\n",
                        b->rank, cp.gid(), line.c_str());
                break;
            }
        }
    });
}
