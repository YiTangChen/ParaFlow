#include "MPASOGPUTracer.h"

#include "Grid.h"
#include "Solution.h"
#include "MPASODeviceField.h"
#include "Kernel/MPASOTracerKernels.h"

#include <cuda_runtime.h>
#include <chrono>
#include <cstring>
#include <vector>

// VECTOR3 stores a single `double vec[3]` and no virtuals, so it is
// standard-layout with size == sizeof(mpaso_vec3) == 24 bytes. We exploit that
// to reinterpret_cast arrays of VECTOR3 to arrays of mpaso_vec3 without any
// per-element copying. Any future change to VECTOR3's layout would break this
// assumption — the static_assert below catches it at compile time.
static_assert(sizeof(VECTOR3) == sizeof(mpaso_vec3),
              "VECTOR3 and mpaso_vec3 must be bit-compatible for reinterpret_cast");

namespace {

using Clock = std::chrono::steady_clock;

inline Clock::time_point gpu_timer_now()
{
    return Clock::now();
}

inline float gpu_timer_ms_since(Clock::time_point start)
{
    return std::chrono::duration<float, std::milli>(Clock::now() - start).count();
}

// Pack Solution's jagged [ts][node] layout into a single flat [ts * node]
// buffer of mpaso_vec3 that uploadVelocityWindow expects.
std::vector<mpaso_vec3> flattenSolution(Solution* sol, int n_vert_nodes, int n_timesteps)
{
    std::vector<mpaso_vec3> out(static_cast<size_t>(n_timesteps) * n_vert_nodes);
    for (int ts = 0; ts < n_timesteps; ++ts) {
        const VECTOR3* src = sol->GetTimestepData(ts);
        std::memcpy(&out[static_cast<size_t>(ts) * n_vert_nodes],
                    src,
                    static_cast<size_t>(n_vert_nodes) * sizeof(mpaso_vec3));
    }
    return out;
}

} // namespace

namespace mpaso_gpu_host {

struct GPUBlockContext {
    MPASODeviceField field;
    bool topology_uploaded = false;
};

bool isAvailable()
{
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err != cudaSuccess) {
        // Clear the sticky error so later cudaGetLastError calls aren't confused.
        (void)cudaGetLastError();
        return false;
    }
    return count > 0;
}

GPUBlockContext* CreateGPUBlockContext(MPASOGrid* grid,
                                       GPUTimingBreakdown* timing_out)
{
    if (grid == nullptr) return nullptr;

    const bool collect_timing = (timing_out != nullptr);
    auto t_pipeline = collect_timing ? gpu_timer_now() : Clock::time_point{};
    auto t_host_prepare = collect_timing ? gpu_timer_now() : Clock::time_point{};

    GPUBlockContext* context = new GPUBlockContext;
    MPASODeviceField& field = context->field;
    field.n_cells            = grid->getNumCell();
    field.n_local_vertices   = grid->getNumLocalVert();
    field.n_max_edges        = grid->getNMaxEdges();
    field.n_vert_levels      = grid->getNVertLevels();
    field.n_vert_levels_p1   = field.n_vert_levels + 1;
    field.n_timesteps_loaded = grid->getNTimestepsLoaded();
    field.earth_radius       = grid->getEarthRadius();

    if (timing_out) timing_out->host_prepare_ms += gpu_timer_ms_since(t_host_prepare);

    auto t_upload_topology = collect_timing ? gpu_timer_now() : Clock::time_point{};
    field.uploadTopology(
        reinterpret_cast<const mpaso_vec3*>(grid->getCellCoord()),
        reinterpret_cast<const mpaso_vec3*>(grid->getVertexCoord()),
        grid->getVerticesOnCell(),
        grid->getCellsOnCell(),
        grid->getNumVerticesOnCell(),
        grid->getMaxLevelCell());
    context->topology_uploaded = true;

    if (timing_out) {
        timing_out->upload_topology_ms += gpu_timer_ms_since(t_upload_topology);
        timing_out->pipeline_wall_ms += gpu_timer_ms_since(t_pipeline);
    }
    return context;
}

void UploadGPUVelocityWindow(GPUBlockContext* context,
                             MPASOGrid*       grid,
                             Solution*        pSolution,
                             Solution*        vSolution,
                             bool             use_real_timestamps,
                             GPUTimingBreakdown* timing_out)
{
    if (context == nullptr || grid == nullptr || pSolution == nullptr || vSolution == nullptr) return;
    if (!context->topology_uploaded) return;

    const bool collect_timing = (timing_out != nullptr);
    auto t_pipeline = collect_timing ? gpu_timer_now() : Clock::time_point{};
    auto t_host_prepare = collect_timing ? gpu_timer_now() : Clock::time_point{};

    MPASODeviceField& field = context->field;
    field.n_timesteps_loaded = grid->getNTimestepsLoaded();

    const int n_vert_nodes = field.n_local_vertices * field.n_vert_levels;
    const int n_ts         = field.n_timesteps_loaded;

    std::vector<mpaso_vec3> h_vel  = flattenSolution(pSolution, n_vert_nodes, n_ts);
    std::vector<mpaso_vec3> h_vvel = flattenSolution(vSolution, n_vert_nodes, n_ts);
    const double* h_ztop = grid->getZTopAll();

    std::vector<double> h_timestamps(n_ts, 0.0);
    if (use_real_timestamps) {
        const std::vector<double>& ts = grid->getZTopTimestamps();
        for (int i = 0; i < n_ts && i < (int)ts.size(); ++i) h_timestamps[i] = ts[i];
    }

    if (timing_out) timing_out->host_prepare_ms += gpu_timer_ms_since(t_host_prepare);

    auto t_upload_velocity = collect_timing ? gpu_timer_now() : Clock::time_point{};
    field.uploadVelocityWindow(h_vel.data(), h_vvel.data(), h_ztop, h_timestamps.data());

    if (timing_out) {
        timing_out->upload_velocity_ms += gpu_timer_ms_since(t_upload_velocity);
        timing_out->pipeline_wall_ms += gpu_timer_ms_since(t_pipeline);
    }
}

void TraceParticlesOnGPUContext(GPUBlockContext* context,
                                const VECTOR3*  seeds,
                                const int*      seed_cell_id,
                                int             n_particles,
                                int             n_steps,
                                double          dt,
                                double          t_start,
                                bool            use_euler,
                                VECTOR3*        traces_out,
                                const int*      seed_max_steps,
                                const int*      seed_step_offset,
                                int*            steps_taken_out,
                                int*            final_cell_out,
                                int             save_interval,
                                int             max_saved_points,
                                int*            saved_counts_out,
                                float*          kernel_ms_out,
                                GPUTimingBreakdown* timing_out)
{
    if (context == nullptr || n_particles <= 0 || n_steps <= 0) return;

    const bool collect_timing = (timing_out != nullptr) || (kernel_ms_out != nullptr);
    auto t_pipeline = collect_timing ? gpu_timer_now() : Clock::time_point{};
    float kernel_ms = 0.0f;
    mpaso_gpu::LaunchTracer(
        context->field,
        reinterpret_cast<const mpaso_vec3*>(seeds),
        seed_cell_id,
        n_particles,
        n_steps,
        dt,
        t_start,
        use_euler,
        reinterpret_cast<mpaso_vec3*>(traces_out),
        seed_max_steps,
        seed_step_offset,
        steps_taken_out,
        final_cell_out,
        save_interval,
        max_saved_points,
        saved_counts_out,
        collect_timing ? &kernel_ms : nullptr,
        timing_out ? &timing_out->alloc_ms : nullptr,
        timing_out ? &timing_out->upload_particles_ms : nullptr,
        timing_out ? &timing_out->download_results_ms : nullptr,
        timing_out ? &timing_out->free_ms : nullptr);
    if (timing_out) {
        timing_out->kernel_ms += kernel_ms;
        timing_out->pipeline_wall_ms += gpu_timer_ms_since(t_pipeline);
    }
    if (kernel_ms_out) *kernel_ms_out += kernel_ms;
}

void TracePathlineBatchOnGPUContext(GPUBlockContext* context,
                                    const VECTOR3*  seeds,
                                    const int*      seed_cell_id,
                                    const double*   seed_t_start,
                                    const int*      seed_max_steps,
                                    int             n_particles,
                                    int             n_steps,
                                    double          dt,
                                    VECTOR3*        traces_out,
                                    double*         final_time_out,
                                    int*            steps_taken_out,
                                    int*            final_cell_out,
                                    int             save_interval,
                                    int             max_saved_points,
                                    int*            saved_counts_out,
                                    float*          kernel_ms_out,
                                    GPUTimingBreakdown* timing_out)
{
    if (context == nullptr || n_particles <= 0 || n_steps <= 0) return;

    const bool collect_timing = (timing_out != nullptr) || (kernel_ms_out != nullptr);
    auto t_pipeline = collect_timing ? gpu_timer_now() : Clock::time_point{};
    float kernel_ms = 0.0f;
    mpaso_gpu::LaunchPathlineTracer(
        context->field,
        reinterpret_cast<const mpaso_vec3*>(seeds),
        seed_cell_id,
        seed_t_start,
        seed_max_steps,
        n_particles,
        n_steps,
        dt,
        reinterpret_cast<mpaso_vec3*>(traces_out),
        final_time_out,
        steps_taken_out,
        final_cell_out,
        save_interval,
        max_saved_points,
        saved_counts_out,
        collect_timing ? &kernel_ms : nullptr,
        timing_out ? &timing_out->alloc_ms : nullptr,
        timing_out ? &timing_out->upload_particles_ms : nullptr,
        timing_out ? &timing_out->download_results_ms : nullptr,
        timing_out ? &timing_out->free_ms : nullptr);
    if (timing_out) {
        timing_out->kernel_ms += kernel_ms;
        timing_out->pipeline_wall_ms += gpu_timer_ms_since(t_pipeline);
    }
    if (kernel_ms_out) *kernel_ms_out += kernel_ms;
}

void DestroyGPUBlockContext(GPUBlockContext* context,
                            GPUTimingBreakdown* timing_out)
{
    if (context == nullptr) return;
    const bool collect_timing = (timing_out != nullptr);
    auto t_pipeline = collect_timing ? gpu_timer_now() : Clock::time_point{};
    auto t_release = collect_timing ? gpu_timer_now() : Clock::time_point{};
    context->field.release();
    if (timing_out) {
        timing_out->field_release_ms += gpu_timer_ms_since(t_release);
        timing_out->pipeline_wall_ms += gpu_timer_ms_since(t_pipeline);
    }
    delete context;
}

} // namespace mpaso_gpu_host
