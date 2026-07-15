#include "ParaFlow.hpp"

#include <limits>

#ifdef OSUFLOW_ENABLE_CUDA
#include "GPU/CUDA/MPASOGPUTracer.h"
#include <cuda_runtime.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>

// CPU/GPU parity test for the MPAS-O RK4 tracer. Enabled at runtime with:
//   OSUFLOW_GPU_SMOKE=1 <paraflow-command>
//
// Reference trajectory is built by walking CPU CVectorField::at_phys() + a
// host copy of vtCFieldLine::geodesic_step (FieldLine.C:173-195), mirroring
// vtCFieldLine::MPASO_rk4 step-for-step. The GPU tracer is then launched with
// the same seeds and compared position-by-position. Target (option b,
// "acceptably close"): trajectory rel err < 1e-5 after 10 RK4 steps.

static int gpu_smoke_env_int(const char* name, int fallback)
{
    const char* value = std::getenv(name);
    if (value == nullptr) return fallback;
    int parsed = std::atoi(value);
    return parsed > 0 ? parsed : fallback;
}

static double gpu_smoke_env_double(const char* name, double fallback)
{
    const char* value = std::getenv(name);
    if (value == nullptr) return fallback;
    double parsed = std::atof(value);
    return parsed > 0.0 ? parsed : fallback;
}

static int gpu_streamline_chunk_steps()
{
    return gpu_smoke_env_int("OSUFLOW_GPU_STREAMLINE_CHUNK_STEPS", 43200);
}

static void print_gpu_rank_mapping(int rank, int gid)
{
    const char* visible = std::getenv("CUDA_VISIBLE_DEVICES");
    int device_count = 0;
    cudaError_t count_err = cudaGetDeviceCount(&device_count);
    if (count_err != cudaSuccess || device_count <= 0) {
        std::fprintf(stderr,
            "GPU_RANK_MAP rank=%d gid=%d CUDA_VISIBLE_DEVICES=%s selected_device=-1 device_count=0 gpu_name=\"unavailable\" cuda_error=\"%s\"\n",
            rank, gid, visible ? visible : "(unset)",
            cudaGetErrorString(count_err));
        (void)cudaGetLastError();
        return;
    }

    int selected = -1;
    cudaError_t dev_err = cudaGetDevice(&selected);
    if (dev_err != cudaSuccess || selected < 0 || selected >= device_count) {
        std::fprintf(stderr,
            "GPU_RANK_MAP rank=%d gid=%d CUDA_VISIBLE_DEVICES=%s selected_device=-1 device_count=%d gpu_name=\"unavailable\" cuda_error=\"%s\"\n",
            rank, gid, visible ? visible : "(unset)", device_count,
            cudaGetErrorString(dev_err));
        (void)cudaGetLastError();
        return;
    }

    cudaDeviceProp prop;
    cudaError_t prop_err = cudaGetDeviceProperties(&prop, selected);
    std::fprintf(stderr,
        "GPU_RANK_MAP rank=%d gid=%d CUDA_VISIBLE_DEVICES=%s selected_device=%d device_count=%d gpu_name=\"%s\"\n",
        rank, gid, visible ? visible : "(unset)", selected, device_count,
        (prop_err == cudaSuccess) ? prop.name : "unavailable");
    if (prop_err != cudaSuccess) (void)cudaGetLastError();
}

static void accumulate_gpu_timing(BlockTiming& timing,
                                  const mpaso_gpu_host::GPUTimingBreakdown& gpu)
{
    constexpr double ms_to_s = 0.001;
    timing.t_trace_integrate_gpu   += gpu.kernel_ms * ms_to_s;

    // Transfer drill-down: 3 meaningful stages kept + everything sub-millisecond
    // (per-launch alloc / particle H2D / result D2H / free / release) folded into misc.
    const double flatten  = gpu.host_prepare_ms;
    const double topo     = gpu.upload_topology_ms;
    const double velocity = gpu.upload_velocity_ms;
    const double misc     = gpu.alloc_ms + gpu.upload_particles_ms
                          + gpu.download_results_ms + gpu.free_ms + gpu.field_release_ms;
    timing.t_transfer_flatten         += flatten  * ms_to_s;
    timing.t_transfer_upload_topology += topo     * ms_to_s;
    timing.t_transfer_upload_velocity += velocity * ms_to_s;
    timing.t_transfer_misc            += misc     * ms_to_s;

    // trace_transfer is exactly the sum of its children (single axis, drill-downable),
    // and mirrors the CPU path where this bucket is 0 — so both devices decompose
    // trace_local_wall into the same {prepare, transfer, integrate, postprocess, enqueue}.
    timing.t_trace_transfer           += (flatten + topo + velocity + misc) * ms_to_s;
}

#include "gpu_parity_tests.hpp"  // CPU-vs-GPU parity/smoke tests (dormant unless OSUFLOW_GPU_SMOKE set)
#endif // OSUFLOW_ENABLE_CUDA

namespace {
constexpr int PF_UNLIMITED_STEPS = -1;
constexpr int PF_STEP_LIMIT_SENTINEL = std::numeric_limits<int>::max() / 4;

inline bool pf_unlimited_steps(int max_steps)
{
    return max_steps == PF_UNLIMITED_STEPS;
}

inline int pf_remaining_steps(int max_steps, int nsteps_done)
{
    return pf_unlimited_steps(max_steps)
        ? PF_STEP_LIMIT_SENTINEL
        : std::max(0, max_steps - nsteps_done);
}

inline bool pf_finished_by_steps(int max_steps, int nsteps_done)
{
    return !pf_unlimited_steps(max_steps) && nsteps_done >= max_steps;
}
} // namespace


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

    // GPU execution flag. Resolve against runtime availability below once we
    // know whether the binary was built with CUDA and a device is present.
    bool requested_gpu = config["useGPU"] ? config["useGPU"].as<bool>() : false;
#ifdef OSUFLOW_ENABLE_CUDA
    if (requested_gpu) {
        if (mpaso_gpu_host::isAvailable()) {
            this->useGPU = true;
            if (world.rank() == 0)
                std::fprintf(stderr, "[ParaFlow] useGPU=true — CUDA device available, using GPU tracer\n");
        } else {
            this->useGPU = false;
            if (world.rank() == 0)
                std::fprintf(stderr, "[ParaFlow] useGPU=true but no CUDA device visible — falling back to CPU tracer\n");
        }
    } else {
        this->useGPU = false;
    }
#else
    if (requested_gpu && world.rank() == 0)
        std::fprintf(stderr, "[ParaFlow] useGPU=true but binary built without CUDA — falling back to CPU tracer\n");
    this->useGPU = false;
#endif

    try {
        if(config["nproc"])
            this->size = config["nproc"].as<int>();
        else
            throw "The number of processes isn't specified correctly.\n";
        if(config["nblocks"])
            this->nblocks = config["nblocks"].as<int>();
        else
            throw "The number of nblocks isn't specified correctly.\n";
        bool requested_time_varying =
            config["isTimeVarying"] ? config["isTimeVarying"].as<bool>() : false;
        if(config["maxsteps"] && !config["maxsteps"].IsNull())
            this->max_steps = config["maxsteps"].as<int>();
        else if (requested_time_varying)
            this->max_steps = PF_UNLIMITED_STEPS;
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
            fprintf(stderr, "TIMING phase=run_seed_read rank=0 n_seeds=%zu t=%.6f\n",
                    seeds.size(), t_seed_read);
        fprintf(stderr, "TIMING phase=run_seed_bcast rank=%d t=%.6f\n",
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
        pf_accum(b->timing.t_trace_dequeue, _t_comm, enable_timing);

        if(b->currentSeeds.size() == 0)
            continue;
        double _t_local = pf_now(enable_timing);
        double _t_prepare = pf_now(enable_timing);
        // std::cerr << "[rank " << b->rank << ", gid " << cp.gid() << "] Tracing " << b->currentSeeds.size() << " seeds\n";
        int maxRemaining = 0;
        std::vector<int> perSeedMax;
        perSeedMax.reserve(b->startPts.size());
        for (auto& pt : b->startPts) {
            int rem = pf_remaining_steps(this->max_steps, pt.nsteps);
            maxRemaining = std::max(maxRemaining, rem);
            // Per-seed budget: each streamline gets its OWN remaining budget so a
            // re-injected seed is not re-granted a fresh (batch-max) budget. This
            // matches the GPU's per-seed cap (trace_block_gpu) and makes the total
            // step count partition-invariant.
            perSeedMax.push_back(rem);
        }
        pf_accum(b->timing.t_trace_prepare, _t_prepare, enable_timing);

        double _t_comp = pf_now(enable_timing);
        b->GenStreamLineByOSUFlow(maxRemaining, this->interval, perSeedMax);
        pf_accum(b->timing.t_trace_integrate_cpu, _t_comp, enable_timing);
        if (enable_timing)
            for (auto s : b->sl_stepcounts) b->timing.n_steps_total += s;
        // deal with current seeds
        int sent = 0;
        int seedCnt = 0;
        double _t_post = pf_now(enable_timing);
        double enqueue_local = 0.0;
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

            // A seed that exhausts its per-seed step budget is finished even though
            // currseg.nsteps tops out at max_steps-1 (the seed counts as one stored
            // point, so the integrator takes at most maxPoints-1 RK steps). Mirrors
            // the GPU's budget-exhausted check in trace_block_gpu.
            if(pf_finished_by_steps(this->max_steps, currseg.nsteps) ||
               (!pf_unlimited_steps(this->max_steps) &&
                currseg.nsteps >= this->max_steps - 1))
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
                    double _t_enqueue = pf_now(enable_timing);
                    cp.enqueue(bid, endPt);
                    if (enable_timing) {
                        double elapsed = MPI_Wtime() - _t_enqueue;
                        enqueue_local += elapsed;
                        b->timing.t_trace_enqueue += elapsed;
                    }
                    sent++;
                }
            }
            seedCnt++;
        }
        if (enable_timing) {
            double post_elapsed = MPI_Wtime() - _t_post;
            b->timing.t_trace_postprocess += std::max(0.0, post_elapsed - enqueue_local);
        }
        b->endTracing();
        pf_accum(b->timing.t_trace_local_wall, _t_local, enable_timing);
    } while (cp.fill_incoming());
}

#ifdef OSUFLOW_ENABLE_CUDA
void ParaFlow::trace_block_gpu(Block*                             b,
                               const diy::Master::ProxyWithLink&  cp)
{
    diy::Link* l = cp.link();
    const int chunkSteps = gpu_streamline_chunk_steps();

    do
    {
        double _t_comm = pf_now(enable_timing);
        this->deq_incoming_iexchange(b, cp);
        pf_accum(b->timing.t_trace_dequeue, _t_comm, enable_timing);

        if (b->currentSeeds.empty())
            continue;
        double _t_local = pf_now(enable_timing);
        double _t_prepare = pf_now(enable_timing);

        CVectorField* field = b->osuflow->GetFlowField();
        MPASOGrid* grid = dynamic_cast<MPASOGrid*>(field->GetGrid());
        Solution*  pSol = field->GetHorizontalSolution();
        Solution*  vSol = field->GetVerticalSolution();
        if (!grid || !pSol || !vSol) {
            std::fprintf(stderr,
                "[ParaFlow] GPU streamline skipped for gid=%d: grid/solution unavailable, using CPU\n",
                cp.gid());
            this->trace_block(b, cp);
            return;
        }

        const int n_particles = (int)b->currentSeeds.size();
        std::vector<int> remaining_steps(n_particles, 0);
        for (int i = 0; i < n_particles && i < (int)b->startPts.size(); i++) {
            int remaining = pf_remaining_steps(this->max_steps, b->startPts[i].nsteps);
            // CPU streamline path (OSUFlow::computeFieldLine with m_nMaxsize) advances
            // at most (maxsteps - 1) RK steps because the seed counts as one stored point.
            // Keep GPU launch budget aligned with that legacy CPU semantics.
            if (!pf_unlimited_steps(this->max_steps))
                remaining = std::max(0, remaining - 1);
            remaining_steps[i] = remaining;
        }

        std::vector<int> current_cells(n_particles, -1);
        for (int i = 0; i < n_particles; i++) {
            int fromCell = (i < (int)b->fromCells.size()) ? b->fromCells[i] : -1;
            VECTOR4 v4;
            PointInfo ci;
            ci.phyCoord = b->currentSeeds[i];
            ci.fromCell = fromCell;
            ci.inCell = fromCell;
            int ok = field->at_phys(fromCell, b->currentSeeds[i], ci, 0.0, v4, nullptr);
            if (ok == 1 && ci.inCell >= 0)
                current_cells[i] = ci.inCell / b->nVertLevels;
        }

        std::vector<VECTOR3> current_pos = b->currentSeeds;
        std::vector<Segment> segments(n_particles);
        std::vector<int> total_steps_taken(n_particles, 0);
        std::vector<unsigned char> active(n_particles, 0);
        b->toCells.assign(n_particles, -1);

        int active_count = 0;
        for (int i = 0; i < n_particles && i < (int)b->startPts.size(); i++) {
            Segment& seg = segments[i];
            seg.gid    = b->startPts[i].gid;
            seg.pid    = b->startPts[i].pid;
            seg.sid    = b->startPts[i].sid;
            seg.nsteps = b->startPts[i].nsteps;
            if (remaining_steps[i] > 0) {
                active[i] = 1;
                active_count++;
            }
        }
        pf_accum(b->timing.t_trace_prepare, _t_prepare, enable_timing);

        mpaso_gpu_host::GPUTimingBreakdown gpu_timing;
        if (b->gpuContext == nullptr)
            b->gpuContext = mpaso_gpu_host::CreateGPUBlockContext(
                grid, enable_timing ? &gpu_timing : nullptr);
        if (b->gpuContext == nullptr) {
            std::fprintf(stderr,
                "[ParaFlow] GPU streamline skipped for gid=%d: could not create GPU context, using CPU\n",
                cp.gid());
            this->trace_block(b, cp);
            return;
        }
        mpaso_gpu_host::UploadGPUVelocityWindow(
            b->gpuContext, grid, pSol, vSol,
            /*use_real_timestamps=*/false,
            enable_timing ? &gpu_timing : nullptr);
        while (active_count > 0)
        {
            double _t_batch_prepare = pf_now(enable_timing);
            std::vector<int> active_indices;
            active_indices.reserve(active_count);
            int launchSteps = 0;
            for (int i = 0; i < n_particles; i++) {
                if (!active[i]) continue;
                int perSeedSteps = std::min(remaining_steps[i], chunkSteps);
                if (perSeedSteps <= 0) continue;
                active_indices.push_back(i);
                launchSteps = std::max(launchSteps, perSeedSteps);
            }
            if (active_indices.empty() || launchSteps <= 0)
                break;

            const int n_launch = (int)active_indices.size();
            std::vector<VECTOR3> launch_seeds(n_launch);
            std::vector<int> launch_cells(n_launch, -1);
            std::vector<int> launch_max_steps(n_launch, 0);
            std::vector<int> launch_step_offsets(n_launch, 0);
            for (int j = 0; j < n_launch; j++) {
                int seedCnt = active_indices[j];
                launch_seeds[j] = current_pos[seedCnt];
                launch_cells[j] = current_cells[seedCnt];
                launch_max_steps[j] = std::min(remaining_steps[seedCnt], launchSteps);
                launch_step_offsets[j] = segments[seedCnt].nsteps;
            }

            const int saveInterval = std::max(1, this->interval);
            const int maxSaved = (launchSteps + saveInterval - 1) / saveInterval + 2;
            std::vector<VECTOR3> gpu_trace((size_t)n_launch * maxSaved);
            std::vector<int> gpu_steps_taken(n_launch, 0);
            std::vector<int> gpu_final_cell(n_launch, -1);
            std::vector<int> gpu_saved_counts(n_launch, 0);
            pf_accum(b->timing.t_trace_prepare, _t_batch_prepare, enable_timing);

            mpaso_gpu_host::TraceParticlesOnGPUContext(
                b->gpuContext,
                launch_seeds.data(),
                launch_cells.data(),
                n_launch,
                launchSteps,
                b->integrationDt,
                /*t_start=*/0.0,
                /*use_euler=*/false,
                gpu_trace.data(),
                launch_max_steps.data(),
                launch_step_offsets.data(),
                gpu_steps_taken.data(),
                gpu_final_cell.data(),
                saveInterval,
                maxSaved,
                gpu_saved_counts.data(),
                nullptr,
                enable_timing ? &gpu_timing : nullptr);

            double _t_post = pf_now(enable_timing);
            double enqueue_local = 0.0;
            for (int j = 0; j < n_launch; j++)
            {
                int seedCnt = active_indices[j];
                Segment& seg = segments[seedCnt];
                int steps = gpu_steps_taken[j];
                int savedCount = gpu_saved_counts[j];
                int row = j * maxSaved;

                VECTOR3 lastPos = (savedCount > 0)
                    ? gpu_trace[row + savedCount - 1]
                    : current_pos[seedCnt];
                int oldStep = seg.nsteps;
                int toCell = (gpu_final_cell[j] >= 0)
                    ? gpu_final_cell[j] * b->nVertLevels
                    : -1;
                b->toCells[seedCnt] = toCell;
                current_pos[seedCnt] = lastPos;
                current_cells[seedCnt] = gpu_final_cell[j];
                seg.nsteps += steps;
                total_steps_taken[seedCnt] += steps;
                remaining_steps[seedCnt] = std::max(0, remaining_steps[seedCnt] - steps);

                // A seed that exhausts its step budget is finished even though
                // seg.nsteps tops out at maxsteps-1 (the seed counts as one point,
                // see the remaining-1 above) and so never reaches maxsteps.
                bool finished = pf_finished_by_steps(this->max_steps, seg.nsteps) ||
                                (!pf_unlimited_steps(this->max_steps) &&
                                 remaining_steps[seedCnt] == 0);
                bool stoppedEarly = steps < launch_max_steps[j] || toCell < 0;
                int nextNeighborId = b->GetNeighborId(toCell);
                bool handoff = !finished && nextNeighborId >= 0;
                bool dead = !finished && stoppedEarly && nextNeighborId < 0;
                bool terminal = finished || handoff || dead;

                int appendCount = savedCount;
                bool forcedChunkFinal = steps > 0 &&
                    ((oldStep + steps) % saveInterval) != 0;
                if (!terminal && forcedChunkFinal && appendCount > 1)
                    appendCount--;

                if (appendCount > 0) {
                    if (seg.coords.empty())
                        seg.coords.push_back(gpu_trace[row]);
                    for (int k = 1; k < appendCount; k++)
                        seg.coords.push_back(gpu_trace[row + k]);
                }

                if (terminal)
                {
                    bool firstStepFailed = (total_steps_taken[seedCnt] == 0 && toCell < 0);
                    if (!firstStepFailed) {
                        bool duplicateLast =
                            !seg.coords.empty() &&
                            seg.coords.back()[0] == lastPos[0] &&
                            seg.coords.back()[1] == lastPos[1] &&
                            seg.coords.back()[2] == lastPos[2];
                        if (!duplicateLast)
                            seg.coords.push_back(lastPos);
                        b->segs.push_back(seg);
                    }

                    if (handoff)
                    {
                        PtInfo endPt;
                        for (int i = 0; i < 3; i++)
                            endPt.coord[i] = lastPos[i];
                        endPt.gid      = seg.gid;
                        endPt.pid      = seg.pid;
                        endPt.sid      = seg.sid + 1;
                        endPt.nsteps   = seg.nsteps;
                        endPt.fromCell = b->GetGlobalCellidx(toCell);

                        diy::BlockID bid = l->target(nextNeighborId);
                        double _t_enqueue = pf_now(enable_timing);
                        cp.enqueue(bid, endPt);
                        if (enable_timing) {
                            double elapsed = MPI_Wtime() - _t_enqueue;
                            enqueue_local += elapsed;
                            b->timing.t_trace_enqueue += elapsed;
                        }
                    }

                    active[seedCnt] = 0;
                    active_count--;
                }
            }
            if (enable_timing) {
                double post_elapsed = MPI_Wtime() - _t_post;
                b->timing.t_trace_postprocess += std::max(0.0, post_elapsed - enqueue_local);
            }
        }
        if (enable_timing) {
            accumulate_gpu_timing(b->timing, gpu_timing);
        }

        b->sl_stepcounts = total_steps_taken;
        if (enable_timing)
            for (auto s : total_steps_taken) b->timing.n_steps_total += s;
        b->endTracing();
        pf_accum(b->timing.t_trace_local_wall, _t_local, enable_timing);
    } while (cp.fill_incoming());
}
#endif

bool ParaFlow::trace_block_iexchange(Block*                              b,
	                                    const diy::Master::ProxyWithLink&    cp)
{
#ifdef OSUFLOW_ENABLE_CUDA
    if (this->useGPU)
        this->trace_block_gpu(b, cp);
    else
        this->trace_block(b, cp);
#else
    this->trace_block(b, cp);
#endif
    return true;
}

void ParaFlow::GenStreamLines(std::list<std::vector<VECTOR3>>& streamlines, double dt_sign, bool writeToDisk,
                              std::vector<char>* segBytesOut)
{
    diy::ContiguousAssigner   assigner(this->size, this->nblocks);
    int rank = world.rank();
    double _t_run = pf_now(enable_timing);   // end-to-end driver wall: block_load + exchange + write
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
        if (enable_timing) { load_timing::g.enabled = true; load_timing::g.reset(); }
        double _t_load = pf_now(enable_timing);
        b->set_data(gid, this->meshFile.c_str(), this->vectorFiles[0].c_str(), this->seeds, this->areaIndices, this->cfg);
        // Direction control (same method, opposite integration sign). set_data set
        // integrationDt = cfg.dt (> 0); dt_sign = -1 integrates the identical RK4 +
        // cross-block handoff backward from each seed. dt_sign = +1 leaves it as-is,
        // so existing callers (ParaFlow_streamline) are unaffected.
        if (dt_sign < 0.0) b->integrationDt = -b->integrationDt;
        pf_accum(b->timing.t_block_load, _t_load, enable_timing);
        if (enable_timing) {
            b->timing.n_seeds_initial      = (int)b->currentSeeds.size();
            b->timing.t_load_grid_build    = load_timing::g.grid_build_s;
            b->timing.t_load_read          = load_timing::g.read_s;
            b->timing.t_load_convert       = load_timing::g.convert_s;
            b->timing.t_load_seed_filter   = load_timing::g.seed_filter_s;
            load_timing::g.enabled         = false;
            b->timing.mem_vmrss_after_kb   = pf_read_vmrss_kb();
            b->timing.n_local_cells        = b->nLocalCells;
            b->timing.n_global_cells       = b->nGlobalCells;
            MPASOGrid* grid = dynamic_cast<MPASOGrid*>(b->osuflow->GetFlowField()->GetGrid());
            if (grid) {
                b->timing.mem_grid_bytes     = grid->getGridMemBytes();
                b->timing.mem_solution_bytes = grid->getSolutionMemBytes();
                // Block-owned index arrays (not visible to MPASOGrid)
                b->timing.mem_grid_bytes    += (size_t)b->nLocalCells  * sizeof(int); // LocalCell2GlobalCell
                b->timing.mem_grid_bytes    += (size_t)b->nGlobalCells * sizeof(int); // GlobalCell2LocalCell
                b->timing.mem_grid_bytes    += (size_t)b->nGlobalCells * sizeof(int); // areaIndicesArr
                b->timing.mem_grid_bytes    += b->areaID2neighborID.size() * sizeof(int); // areaID2neighborID
            }
        }

        diy::BlockID neighbor;               // one neighbor in the neighborhood
        for(int neighborIdx: b->neighborIndices) {
            neighbor.gid  = neighborIdx;                    // gid of the neighbor block
            neighbor.proc = assigner.rank(neighbor.gid);    // process of the neighbor block
            link->add_neighbor(neighbor);
        }
        master.add(gid, b, link); // add block to the master (mandatory)

#ifdef OSUFLOW_ENABLE_CUDA
        if (this->useGPU)
            print_gpu_rank_mapping(rank, gid);
        run_gpu_smoke_test(b, rank);
#endif
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
        // dt_sign < 0 (backward pass) gets its own file so a forward+backward
        // pair with writeToDisk=true doesn't have the second write clobber the
        // first. Forward keeps the plain name so single-direction callers
        // (ParaFlow_streamline) are unaffected.
        std::string filename = this->outputDir + "/" + std::to_string(gid)
                              + (dt_sign < 0.0 ? "_bwd" : "") + ".bin";
        double _t_write = pf_now(enable_timing);
        Block* blk = (Block*)master.block(master.lid(gid));
        blk->write_trajectory(filename, streamlines, this->interval, writeToDisk);
        if (segBytesOut) blk->append_trajectory_bytes(*segBytesOut);
        pf_accum(blk->timing.t_output_write, _t_write, enable_timing);
    }

    // End-to-end driver wall clock (block_load + exchange + write), reduced across ranks.
    if (enable_timing) {
        double t_run = MPI_Wtime() - _t_run;
        fprintf(stderr, "TIMING phase=run_streamline_rank rank=%d t=%.6f\n", rank, t_run);
        double run_min, run_max, run_sum;
        MPI_Reduce(&t_run, &run_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_run, &run_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_run, &run_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
            fprintf(stderr,
                "TIMING phase=run_streamline_summary min=%.6f max=%.6f avg=%.6f imbalance=%.3f\n",
                run_min, run_max, run_sum / this->size,
                (run_min > 0.0) ? run_max / run_min : 0.0);
    }

    // All measurement output — fully gated: zero overhead when enable_timing=false
    if (enable_timing) {
        double rank_block_load = 0.0;   // per-rank sum of block_load, for the reduce below
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp) {
            b->timing.mem_peak_vmhwm_kb = pf_read_vmhwm_kb();
            fprintf(stderr,
                "TIMING phase=block_load    rank=%d gid=%d t=%.6f nseeds=%d\n",
                b->rank, cp.gid(), b->timing.t_block_load, b->timing.n_seeds_initial);
            if (b->timing.t_load_grid_build > 0.0)
                fprintf(stderr, "TIMING phase=load_grid_build rank=%d gid=%d t=%.6f\n",
                        b->rank, cp.gid(), b->timing.t_load_grid_build);
            if (b->timing.t_load_read > 0.0)
                fprintf(stderr, "TIMING phase=load_read rank=%d gid=%d t=%.6f\n",
                        b->rank, cp.gid(), b->timing.t_load_read);
            if (b->timing.t_load_convert > 0.0)
                fprintf(stderr, "TIMING phase=load_convert rank=%d gid=%d t=%.6f\n",
                        b->rank, cp.gid(), b->timing.t_load_convert);
            if (b->timing.t_load_seed_filter > 0.0)
                fprintf(stderr, "TIMING phase=load_seed_filter rank=%d gid=%d t=%.6f\n",
                        b->rank, cp.gid(), b->timing.t_load_seed_filter);
            fprintf(stderr,
                "TIMING phase=trace_local_wall rank=%d gid=%d t=%.6f nsteps=%ld nrecv=%d\n",
                b->rank, cp.gid(), b->timing.t_trace_local_wall,
                b->timing.n_steps_total, b->timing.n_particles_received);
            if (b->timing.t_trace_prepare > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_prepare rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_prepare);
            if (b->timing.t_trace_integrate_cpu > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_integrate_cpu rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_integrate_cpu);
            if (b->timing.t_trace_integrate_gpu > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_integrate_gpu rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_integrate_gpu);
            if (b->timing.t_trace_transfer > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_transfer rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_transfer);
            if (b->timing.t_trace_postprocess > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_postprocess rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_postprocess);
            if (b->timing.t_trace_enqueue > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_enqueue rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_enqueue);
            // Drill-down of trace_transfer (children sum to it; GPU only).
            if (b->timing.t_trace_transfer > 0.0) {
                fprintf(stderr,
                    "TIMING phase=transfer_flatten rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_transfer_flatten);
                fprintf(stderr,
                    "TIMING phase=transfer_upload_topology rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_transfer_upload_topology);
                fprintf(stderr,
                    "TIMING phase=transfer_upload_velocity rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_transfer_upload_velocity);
                fprintf(stderr,
                    "TIMING phase=transfer_misc rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_transfer_misc);
            }
            fprintf(stderr,
                "TIMING phase=trace_dequeue rank=%d gid=%d t=%.6f\n",
                b->rank, cp.gid(), b->timing.t_trace_dequeue);
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
            fprintf(stderr,
                "MEM_CELLCOUNT  rank=%d gid=%d n_local_cells=%d n_global_cells=%d\n",
                b->rank, cp.gid(), b->timing.n_local_cells, b->timing.n_global_cells);
            rank_block_load += b->timing.t_block_load;
        });
        double bl_min, bl_max, bl_sum;
        MPI_Reduce(&rank_block_load, &bl_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&rank_block_load, &bl_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&rank_block_load, &bl_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
            fprintf(stderr,
                "TIMING phase=block_load_summary min=%.6f max=%.6f avg=%.6f imbalance=%.3f\n",
                bl_min, bl_max, bl_sum / this->size, (bl_min > 0.0) ? bl_max / bl_min : 0.0);
        fprintf(stderr, "TIMING phase=trace_exchange_rank rank=%d t=%.6f\n", rank, t_iexchange);
        double iex_min, iex_max, iex_sum;
        MPI_Reduce(&t_iexchange, &iex_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_iexchange, &iex_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_iexchange, &iex_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
            fprintf(stderr,
                "TIMING phase=trace_exchange_summary min=%.6f max=%.6f avg=%.6f imbalance=%.3f\n",
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
        pf_accum(b->timing.t_trace_dequeue, _t_comm, enable_timing);

        if (b->currentSeeds.size() == 0)
            continue;
        double _t_local = pf_now(enable_timing);
        double _t_prepare = pf_now(enable_timing);

        // Build per-seed start times from PtInfo.ts (0.0 for initial seeds,
        // non-zero for particles that were forwarded from a neighbour block).
        std::vector<double> tarray(b->startPts.size());
        for (size_t i = 0; i < b->startPts.size(); i++)
            tarray[i] = b->startPts[i].ts;

        // ── Advect particles through the time-varying field ──────────────────
        int maxRemaining = 0;
        for (auto& pt : b->startPts)
            maxRemaining = std::max(maxRemaining, pf_remaining_steps(this->max_steps, pt.nsteps));
        pf_accum(b->timing.t_trace_prepare, _t_prepare, enable_timing);

        double _t_comp = pf_now(enable_timing);
        b->GenPathLineByOSUFlow(maxRemaining, tarray.data(), this->interval);
        pf_accum(b->timing.t_trace_integrate_cpu, _t_comp, enable_timing);
        if (enable_timing)
            for (auto s : b->pl_stepcounts) b->timing.n_steps_total += s;

        int sent = 0;
        int seedCnt = 0;
        double maxT = b->getWindowMaxTime();
        double _t_post = pf_now(enable_timing);
        double enqueue_local = 0.0;
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

            bool finished = pf_finished_by_steps(this->max_steps, currseg.nsteps);

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
                    double _t_enqueue = pf_now(enable_timing);
                    cp.enqueue(bid, endPt);
                    if (enable_timing) {
                        double elapsed = MPI_Wtime() - _t_enqueue;
                        enqueue_local += elapsed;
                        b->timing.t_trace_enqueue += elapsed;
                    }
                    sent++;
                }
            }
        }
        if (enable_timing) {
            double post_elapsed = MPI_Wtime() - _t_post;
            b->timing.t_trace_postprocess += std::max(0.0, post_elapsed - enqueue_local);
        }

        int temporal = (int)b->temporalBoundaryPts.size();
        b->endTracing();
        pf_accum(b->timing.t_trace_local_wall, _t_local, enable_timing);
    } while (cp.fill_incoming());
}

#ifdef OSUFLOW_ENABLE_CUDA
void ParaFlow::trace_block_pathline_gpu(Block*                              b,
                                        const diy::Master::ProxyWithLink&   cp)
{
    diy::Link* l = cp.link();

    do
    {
        double _t_comm = pf_now(enable_timing);
        this->deq_incoming_iexchange(b, cp);
        pf_accum(b->timing.t_trace_dequeue, _t_comm, enable_timing);

        if (b->currentSeeds.empty())
            continue;
        double _t_local = pf_now(enable_timing);
        double _t_prepare = pf_now(enable_timing);

        CVectorField* field = b->osuflow->GetFlowField();
        MPASOGrid* grid = dynamic_cast<MPASOGrid*>(field->GetGrid());
        Solution*  pSol = field->GetHorizontalSolution();
        Solution*  vSol = field->GetVerticalSolution();
        if (!grid || !pSol || !vSol) {
            std::fprintf(stderr,
                "[ParaFlow] GPU pathline skipped for gid=%d: grid/solution unavailable, using CPU\n",
                cp.gid());
            this->trace_block_pathline(b, cp);
            return;
        }

        const int n_particles = (int)b->currentSeeds.size();
        std::vector<double> tarray(n_particles, 0.0);
        std::vector<int> seed_max_steps(n_particles, 0);
        for (int i = 0; i < n_particles && i < (int)b->startPts.size(); i++)
            tarray[i] = b->startPts[i].ts;

        int maxRemaining = 0;
        int maxWindowSteps = 0;
        const double dt = b->integrationDt;
        const double maxT = b->getWindowMaxTime();
        for (int i = 0; i < n_particles && i < (int)b->startPts.size(); i++) {
            int remaining = pf_remaining_steps(this->max_steps, b->startPts[i].nsteps);
            // Match CPU pathline budget: vtCPathLine counts the seed as one stored
            // point (count starts at 1), so it advances at most (maxsteps - 1) RK
            // steps. Without this the GPU takes one extra step at termination,
            // shifting each pathline's final point. Mirrors the streamline GPU path.
            if (!pf_unlimited_steps(this->max_steps))
                remaining = std::max(0, remaining - 1);
            seed_max_steps[i] = remaining;
            maxRemaining = std::max(maxRemaining, remaining);
            double windowRemaining = maxT - tarray[i];
            int windowSteps = (windowRemaining > 0.0 && dt > 0.0)
                ? (int)std::floor((windowRemaining + 1.0e-9) / dt)
                : 0;
            maxWindowSteps = std::max(maxWindowSteps, std::min(remaining, windowSteps));
        }

        int n_steps = std::min(maxRemaining, maxWindowSteps);
        if (n_steps <= 0) {
            pf_accum(b->timing.t_trace_prepare, _t_prepare, enable_timing);
            double _t_post = pf_now(enable_timing);
            for (int i = 0; i < n_particles && i < (int)b->startPts.size(); i++) {
                if (pf_finished_by_steps(this->max_steps, b->startPts[i].nsteps))
                    continue;
                b->temporalBoundaryPts.push_back(b->startPts[i]);
            }
            pf_accum(b->timing.t_trace_postprocess, _t_post, enable_timing);
            b->endTracing();
            pf_accum(b->timing.t_trace_local_wall, _t_local, enable_timing);
            continue;
        }

        std::vector<int> seed_cells(n_particles, -1);
        for (int i = 0; i < n_particles; i++) {
            int fromCell = (i < (int)b->fromCells.size()) ? b->fromCells[i] : -1;
            VECTOR4 v4;
            PointInfo ci;
            ci.phyCoord = b->currentSeeds[i];
            ci.fromCell = fromCell;
            ci.inCell = fromCell;
            int ok = field->at_phys(fromCell, b->currentSeeds[i], ci, tarray[i], v4, nullptr);
            if (ok == 1 && ci.inCell >= 0)
                seed_cells[i] = ci.inCell / b->nVertLevels;
        }

        const int max_saved_points = n_steps / this->interval + 2;
        std::vector<VECTOR3> gpu_trace((size_t)n_particles * max_saved_points);
        std::vector<double>  gpu_final_time(n_particles, 0.0);
        std::vector<int>     gpu_steps_taken(n_particles, 0);
        std::vector<int>     gpu_final_cell(n_particles, -1);
        std::vector<int>     gpu_saved_counts(n_particles, 0);
        mpaso_gpu_host::GPUTimingBreakdown gpu_timing;
        if (b->gpuContext == nullptr)
            b->gpuContext = mpaso_gpu_host::CreateGPUBlockContext(
                grid, enable_timing ? &gpu_timing : nullptr);
        if (b->gpuContext == nullptr) {
            std::fprintf(stderr,
                "[ParaFlow] GPU pathline skipped for gid=%d: could not create GPU context, using CPU\n",
                cp.gid());
            this->trace_block_pathline(b, cp);
            return;
        }
        mpaso_gpu_host::UploadGPUVelocityWindow(
            b->gpuContext, grid, pSol, vSol,
            /*use_real_timestamps=*/true,
            enable_timing ? &gpu_timing : nullptr);
        pf_accum(b->timing.t_trace_prepare, _t_prepare, enable_timing);

        mpaso_gpu_host::TracePathlineBatchOnGPUContext(
            b->gpuContext,
            b->currentSeeds.data(),
            seed_cells.data(),
            tarray.data(),
            seed_max_steps.data(),
            n_particles,
            n_steps,
            dt,
            gpu_trace.data(),
            gpu_final_time.data(),
            gpu_steps_taken.data(),
            gpu_final_cell.data(),
            this->interval,
            max_saved_points,
            gpu_saved_counts.data(),
            nullptr,
            enable_timing ? &gpu_timing : nullptr);
        if (enable_timing) {
            accumulate_gpu_timing(b->timing, gpu_timing);
        }

        b->pl_stepcounts = gpu_steps_taken;
        b->toCells.assign(n_particles, -1);
        if (enable_timing)
            for (auto s : gpu_steps_taken) b->timing.n_steps_total += s;

        double _t_post = pf_now(enable_timing);
        double enqueue_local = 0.0;
        for (int seedCnt = 0; seedCnt < n_particles; seedCnt++)
        {
            Segment currseg;
            currseg.gid    = b->startPts[seedCnt].gid;
            currseg.pid    = b->startPts[seedCnt].pid;
            currseg.sid    = b->startPts[seedCnt].sid;
            currseg.nsteps = b->startPts[seedCnt].nsteps + gpu_steps_taken[seedCnt];

            int saved = gpu_saved_counts[seedCnt];
            int row = seedCnt * max_saved_points;
            for (int s = 0; s < saved; s++)
                currseg.coords.push_back(gpu_trace[row + s]);

            if (currseg.coords.empty())
                continue;

            int toCell = (gpu_final_cell[seedCnt] >= 0)
                ? gpu_final_cell[seedCnt] * b->nVertLevels
                : -1;
            b->toCells[seedCnt] = toCell;

            bool finished = pf_finished_by_steps(this->max_steps, currseg.nsteps);
            int nextNeighborId = b->GetNeighborId(toCell);
            double exitTime = gpu_final_time[seedCnt];
            bool isTemporal = (!finished && nextNeighborId < 0 && exitTime >= maxT - b->integrationDt);

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
                endPt.ts       = exitTime;

                if (isTemporal)
                {
                    endPt.fromCell = (toCell >= 0) ? b->GetGlobalCellidx(toCell) : -1;
                    b->temporalBoundaryPts.push_back(endPt);
                }
                else if (nextNeighborId >= 0)
                {
                    endPt.fromCell = b->GetGlobalCellidx(toCell);
                    diy::BlockID bid = l->target(nextNeighborId);
                    double _t_enqueue = pf_now(enable_timing);
                    cp.enqueue(bid, endPt);
                    if (enable_timing) {
                        double elapsed = MPI_Wtime() - _t_enqueue;
                        enqueue_local += elapsed;
                        b->timing.t_trace_enqueue += elapsed;
                    }
                }
            }
        }
        if (enable_timing) {
            double post_elapsed = MPI_Wtime() - _t_post;
            b->timing.t_trace_postprocess += std::max(0.0, post_elapsed - enqueue_local);
        }

        b->endTracing();
        pf_accum(b->timing.t_trace_local_wall, _t_local, enable_timing);
    } while (cp.fill_incoming());
}
#endif

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
    double _t_run = pf_now(enable_timing);   // end-to-end driver wall: block_load + windows + write
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
        if (enable_timing) { load_timing::g.enabled = true; load_timing::g.reset(); }
        double _t_load = pf_now(enable_timing);
        b->set_data(gid, this->meshFile.c_str(), this->vectorFiles[0].c_str(),
                    this->seeds, this->areaIndices, this->cfg);
        pf_accum(b->timing.t_block_load, _t_load, enable_timing);
        if (enable_timing) {
            b->timing.n_seeds_initial      = (int)b->currentSeeds.size();
            b->timing.t_load_grid_build    = load_timing::g.grid_build_s;
            b->timing.t_load_read          = load_timing::g.read_s;
            b->timing.t_load_convert       = load_timing::g.convert_s;
            b->timing.t_load_seed_filter   = load_timing::g.seed_filter_s;
            load_timing::g.enabled         = false;
            b->timing.mem_vmrss_after_kb   = pf_read_vmrss_kb();
            b->timing.n_local_cells        = b->nLocalCells;
            b->timing.n_global_cells       = b->nGlobalCells;
            MPASOGrid* grid = dynamic_cast<MPASOGrid*>(b->osuflow->GetFlowField()->GetGrid());
            if (grid) {
                b->timing.mem_grid_bytes     = grid->getGridMemBytes();
                b->timing.mem_solution_bytes = grid->getSolutionMemBytes();
                // Block-owned index arrays (not visible to MPASOGrid)
                b->timing.mem_grid_bytes    += (size_t)b->nLocalCells  * sizeof(int); // LocalCell2GlobalCell
                b->timing.mem_grid_bytes    += (size_t)b->nGlobalCells * sizeof(int); // GlobalCell2LocalCell
                b->timing.mem_grid_bytes    += (size_t)b->nGlobalCells * sizeof(int); // areaIndicesArr
                b->timing.mem_grid_bytes    += b->areaID2neighborID.size() * sizeof(int); // areaID2neighborID
            }
        }

        diy::BlockID neighbor;
        for (int neighborIdx : b->neighborIndices) {
            neighbor.gid  = neighborIdx;
            neighbor.proc = assigner.rank(neighbor.gid);
            link->add_neighbor(neighbor);
        }
        master.add(gid, b, link);

#ifdef OSUFLOW_ENABLE_CUDA
        if (this->useGPU)
            print_gpu_rank_mapping(rank, gid);
        run_gpu_pathline_parity_test(b, rank);
#endif
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
#ifdef OSUFLOW_ENABLE_CUDA
            if (this->useGPU) {
                this->trace_block_pathline_gpu(b, icp);
                return true;
            }
#endif
            return this->trace_block_pathline_iexchange(b, icp);
        });
        world.barrier();
        double t_window = 0.0;
        pf_accum(t_window, _t_iex, enable_timing);
        t_iexchange += t_window;
        if (enable_timing)
            fprintf(stderr, "TIMING phase=dist_trace_window rank=%d window=%d t=%.6f\n",
                    rank, windowIdx, t_window);
        windowIdx++;

        int local_active = 0;
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink&) {
            if (!b->temporalBoundaryPts.empty())
                local_active = 1;
        });
        int active_sum = 0;
        diy::mpi::all_reduce(world, local_active, active_sum, std::plus<int>());
        if (active_sum == 0)
            break;

        // Use the reader's actual global offset to determine if more data exists.
        // This correctly handles file-boundary clamping that would cause windowStart to drift.
        int globalOffset = 0;
        if (!gids.empty())
            globalOffset = ((Block*)master.block(master.lid(gids[0])))->getGlobalTimestepOffset();
        diy::mpi::broadcast(world, globalOffset, 0);
        if (globalOffset >= totalTimestepsInFile)
            break;

        // Load next window and re-inject particles that hit the temporal boundary
        double t_window_load = 0.0;
        double t_window_reinject = 0.0;
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp) {
            double _t_load = pf_now(enable_timing);
            b->advanceTimestep();
            pf_accum(t_window_load, _t_load, enable_timing);

            double _t_reinject = pf_now(enable_timing);
            for (auto& pt : b->temporalBoundaryPts)
                b->addStartPts(pt);
            if (!b->temporalBoundaryPts.empty()) {
                b->setStartPtsDone();
                b->temporalBoundaryPts.clear();
            }
            pf_accum(t_window_reinject, _t_reinject, enable_timing);
        });
        world.barrier();
        if (enable_timing) {
            fprintf(stderr, "TIMING phase=window_load rank=%d window=%d t=%.6f\n",
                    rank, windowIdx - 1, t_window_load);
            fprintf(stderr, "TIMING phase=window_reinject rank=%d window=%d t=%.6f\n",
                    rank, windowIdx - 1, t_window_reinject);
        }
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

    // End-to-end driver wall clock (block_load + windows + write), reduced across ranks.
    if (enable_timing) {
        double t_run = MPI_Wtime() - _t_run;
        fprintf(stderr, "TIMING phase=run_pathline_rank rank=%d t=%.6f\n", rank, t_run);
        double run_min, run_max, run_sum;
        MPI_Reduce(&t_run, &run_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_run, &run_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_run, &run_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
            fprintf(stderr,
                "TIMING phase=run_pathline_summary min=%.6f max=%.6f avg=%.6f imbalance=%.3f\n",
                run_min, run_max, run_sum / this->size,
                (run_min > 0.0) ? run_max / run_min : 0.0);
    }

    // All measurement output — fully gated: zero overhead when enable_timing=false
    if (enable_timing) {
        double rank_block_load = 0.0;   // per-rank sum of block_load, for the reduce below
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp) {
            b->timing.mem_peak_vmhwm_kb = pf_read_vmhwm_kb();
            fprintf(stderr,
                "TIMING phase=block_load    rank=%d gid=%d t=%.6f nseeds=%d\n",
                b->rank, cp.gid(), b->timing.t_block_load, b->timing.n_seeds_initial);
            if (b->timing.t_load_grid_build > 0.0)
                fprintf(stderr, "TIMING phase=load_grid_build rank=%d gid=%d t=%.6f\n",
                        b->rank, cp.gid(), b->timing.t_load_grid_build);
            if (b->timing.t_load_read > 0.0)
                fprintf(stderr, "TIMING phase=load_read rank=%d gid=%d t=%.6f\n",
                        b->rank, cp.gid(), b->timing.t_load_read);
            if (b->timing.t_load_convert > 0.0)
                fprintf(stderr, "TIMING phase=load_convert rank=%d gid=%d t=%.6f\n",
                        b->rank, cp.gid(), b->timing.t_load_convert);
            if (b->timing.t_load_seed_filter > 0.0)
                fprintf(stderr, "TIMING phase=load_seed_filter rank=%d gid=%d t=%.6f\n",
                        b->rank, cp.gid(), b->timing.t_load_seed_filter);
            fprintf(stderr,
                "TIMING phase=trace_local_wall rank=%d gid=%d t=%.6f nsteps=%ld nrecv=%d\n",
                b->rank, cp.gid(), b->timing.t_trace_local_wall,
                b->timing.n_steps_total, b->timing.n_particles_received);
            if (b->timing.t_trace_prepare > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_prepare rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_prepare);
            if (b->timing.t_trace_integrate_cpu > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_integrate_cpu rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_integrate_cpu);
            if (b->timing.t_trace_integrate_gpu > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_integrate_gpu rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_integrate_gpu);
            if (b->timing.t_trace_transfer > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_transfer rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_transfer);
            if (b->timing.t_trace_postprocess > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_postprocess rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_postprocess);
            if (b->timing.t_trace_enqueue > 0.0)
                fprintf(stderr,
                    "TIMING phase=trace_enqueue rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_trace_enqueue);
            // Drill-down of trace_transfer (children sum to it; GPU only).
            if (b->timing.t_trace_transfer > 0.0) {
                fprintf(stderr,
                    "TIMING phase=transfer_flatten rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_transfer_flatten);
                fprintf(stderr,
                    "TIMING phase=transfer_upload_topology rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_transfer_upload_topology);
                fprintf(stderr,
                    "TIMING phase=transfer_upload_velocity rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_transfer_upload_velocity);
                fprintf(stderr,
                    "TIMING phase=transfer_misc rank=%d gid=%d t=%.6f\n",
                    b->rank, cp.gid(), b->timing.t_transfer_misc);
            }
            fprintf(stderr,
                "TIMING phase=trace_dequeue rank=%d gid=%d t=%.6f\n",
                b->rank, cp.gid(), b->timing.t_trace_dequeue);
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
            fprintf(stderr,
                "MEM_CELLCOUNT  rank=%d gid=%d n_local_cells=%d n_global_cells=%d\n",
                b->rank, cp.gid(), b->timing.n_local_cells, b->timing.n_global_cells);
            rank_block_load += b->timing.t_block_load;
        });
        double bl_min, bl_max, bl_sum;
        MPI_Reduce(&rank_block_load, &bl_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&rank_block_load, &bl_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&rank_block_load, &bl_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
            fprintf(stderr,
                "TIMING phase=block_load_summary min=%.6f max=%.6f avg=%.6f imbalance=%.3f\n",
                bl_min, bl_max, bl_sum / this->size, (bl_min > 0.0) ? bl_max / bl_min : 0.0);
        fprintf(stderr, "TIMING phase=trace_exchange_rank rank=%d t=%.6f\n", rank, t_iexchange);
        double iex_min, iex_max, iex_sum;
        MPI_Reduce(&t_iexchange, &iex_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_iexchange, &iex_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_iexchange, &iex_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
            fprintf(stderr,
                "TIMING phase=trace_exchange_summary min=%.6f max=%.6f avg=%.6f imbalance=%.3f\n",
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

    // Peak-RSS probe: gated by enable_timing like all other profiling output, and
    // emitted in the current MEM_PEAK format (vmhwm_kb=) rather than the legacy string.
    if (enable_timing) {
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp) {
            fprintf(stderr, "MEM_PEAK       rank=%d gid=%d vmhwm_kb=%ld\n",
                    b->rank, cp.gid(), pf_read_vmhwm_kb());
        });
    }
}
