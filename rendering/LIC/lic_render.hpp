// lic_render.hpp
// -----------------------------------------------------------------------------
// The LIC method, and nothing else. Given streamlines that were ALREADY generated
// OUTSIDE this folder -- by ParaFlow::GenStreamLines, called from the repo-root
// driver (ParaFlow_lic.cpp / ParaFlow_lic_gpu.cpp) exactly the way
// ParaFlow_streamline.cpp calls it -- this folds them into a Line-Integral-
// Convolution image and writes ONE grayscale PNG on rank 0.
//
// This header does NO streamline generation and does NOT init MPI. It is purely
// the LIC rendering step, so everything under rendering/LIC/ stays "LIC only".
//
// Include AFTER <ParaFlow.hpp> (it needs VECTOR3) and after MPI is initialized
// (the ParaFlow ctor does that in the drivers).
// -----------------------------------------------------------------------------
#ifndef LIC_RENDER_HPP
#define LIC_RENDER_HPP

#include <list>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <random>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <mpi.h>
#include <yaml-cpp/yaml.h>

#ifndef LIC_WITH_MPI
#define LIC_WITH_MPI
#endif
#include "lic_core.hpp"     // LicImage, lic_add_streamline, lic_reduce, PNG (pure engine)

// Which integration directions to run, read from config key `lic_direction`.
// Since a scientific streamline is traced one way from its seed but a LIC line
// should run THROUGH the seed:
//     forward  -> {+1}   backward -> {-1}   both (default) -> {+1, -1}
// The caller runs pf.GenStreamLines(sls, sign) for each returned sign -- that call
// (the streamline "run") lives in the driver, OUTSIDE this folder.
inline std::vector<double> lic_directions(const std::string& configPath, int rank = 0)
{
    YAML::Node config = YAML::LoadFile(configPath);
    std::string d = config["lic_direction"] ? config["lic_direction"].as<std::string>() : "both";
    std::transform(d.begin(), d.end(), d.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    if (d == "forward")  { if (rank == 0) std::cerr << "[lic] direction = forward [+dt]\n";  return {+1.0}; }
    if (d == "backward") { if (rank == 0) std::cerr << "[lic] direction = backward [-dt]\n"; return {-1.0}; }
    if (d != "both" && rank == 0)
        std::cerr << "[lic] unknown lic_direction='" << d << "', using 'both'\n";
    if (rank == 0) std::cerr << "[lic] direction = both [+dt][-dt]\n";
    return {+1.0, -1.0};
}

// Number of seeds, read straight from the seed file named in the config -- so the
// LIC image auto-sizing needs NO access to ParaFlow's private seed list. Binary
// seed files are 3 float64 per seed (24 bytes); text files are one seed per line.
// This is the raw grid count (e.g. 1024*1024), matching what ParaFlow reads in.
inline size_t lic_count_seeds(const std::string& configPath)
{
    YAML::Node config = YAML::LoadFile(configPath);
    if (!config["seedFile"] || config["seedFile"].IsNull()) return 0;
    const std::string seedFile = config["seedFile"].as<std::string>();
    const bool isBin = config["isBinarySeed"] && config["isBinarySeed"].as<bool>();

    if (isBin) {
        std::error_code ec;
        const auto sz = std::filesystem::file_size(seedFile, ec);
        if (ec) return 0;
        return (size_t)sz / (3 * sizeof(double));      // 3 doubles (x,y,z) per seed
    }
    std::ifstream in(seedFile);
    size_t n = 0;
    std::string line;
    while (std::getline(in, line))
        if (line.find(',') != std::string::npos) ++n;  // "x,y,z" lines, like the seed parser
    return n;
}

// Name for this run's output folder: "<mode>_<direction>_<datafile>", e.g.
// "gpu_both_...mpaso.hist.am.timeSeriesStatsMonthly.0015-01-01". `mode` is "gpu"
// or "cpu" (from the driver); the direction and datafile name come from the config
// (datafile = basename of datafile_src without directory or trailing .nc). Two runs
// with the same mode/direction/data reuse the same folder.
inline std::string lic_run_label(const std::string& mode, const std::string& configPath)
{
    YAML::Node config = YAML::LoadFile(configPath);
    std::string d = config["lic_direction"] ? config["lic_direction"].as<std::string>() : "both";
    std::transform(d.begin(), d.end(), d.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    if (d != "forward" && d != "backward") d = "both";

    // datafile_src may be a scalar or a sequence -- use the first entry.
    std::string dataPath;
    if (config["datafile_src"]) {
        if (config["datafile_src"].IsSequence() && config["datafile_src"].size() > 0)
            dataPath = config["datafile_src"][0].as<std::string>();
        else if (config["datafile_src"].IsScalar())
            dataPath = config["datafile_src"].as<std::string>();
    }
    std::string dataName = std::filesystem::path(dataPath).filename().string();
    if (dataName.size() > 3 && dataName.substr(dataName.size() - 3) == ".nc")
        dataName = dataName.substr(0, dataName.size() - 3);        // strip trailing .nc
    if (dataName.empty()) dataName = "data";

    return mode + "_" + d + "_" + dataName;
}

// Which of the two folding strategies the driver should run (config key
// `lic_combine`), so both can be A/B compared with the same binary:
//   false = METHOD 2 (no combine): each rank folds its OWN per-block streamline
//           SEGMENTS straight into a local image (lic_add_streamline), then the
//           per-rank IMAGES are summed (lic_reduce). Distributed; each fragment's
//           footprint is averaged on its own.
//   true  = METHOD 1 (combine): rank 0 reads the per-block <gid>.bin files and
//           REASSEMBLES the segments (by pid/sid) into WHOLE streamlines -- the
//           same thing scripts/plot/read_traces.py does -- then folds those whole
//           lines. Each complete streamline's footprint is averaged as one.
// The two give DIFFERENT images (fragment averages vs whole-line averages); that
// is the point of the comparison. NOTE: <gid>.bin is truncated per direction, so a
// combine run reflects the LAST traced direction only -- use lic_direction:
// forward or backward for a clean combine comparison.
inline bool lic_combine(const std::string& configPath)
{
    YAML::Node config = YAML::LoadFile(configPath);
    return config["lic_combine"] && config["lic_combine"].as<bool>();
}

// Where the traced streamlines end up, from config key `streamline_storage`:
//   memory = fold the streamlines straight out of RAM and DELETE the per-block
//            <gid>.bin files as soon as the run is done. The only artifact left is
//            the PNG. Each direction's traces are also freed right after folding
//            (see the drivers), so peak memory is one direction, not both.
//   disk   = KEEP the <gid>.bin files (moved into the run folder) so the raw
//            trajectories can be reused later (read_traces.py, plotting, ...).
//
// CORE LIMITATION (why `memory` cannot skip the write): ParaFlow::GenStreamLines
// always calls block.hpp's write_trajectory, which fills the in-memory list and
// writes <gid>.bin in the SAME pass -- there is no fill-only path. So `memory`
// guarantees the files are gone afterwards, but the transient write still happens.
// To keep even that write off real disk, point outputDir at a tmpfs (RAM) path:
//     outputDir: /dev/shm/lic_run        # tmpfs = RAM; copy the PNG out afterwards
//
// Back-compat: if `streamline_storage` is absent, fall back to the older boolean
// `keep_streamline_bin` (true -> disk, false -> memory).
inline bool lic_keep_bin(const std::string& configPath, int rank = 0)
{
    YAML::Node config = YAML::LoadFile(configPath);
    if (config["streamline_storage"] && !config["streamline_storage"].IsNull()) {
        std::string s = config["streamline_storage"].as<std::string>();
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c){ return (char)std::tolower(c); });
        if (s == "disk" || s == "storage") {
            if (rank == 0) std::cerr << "[lic] storage = disk (<gid>.bin kept in the run folder)\n";
            return true;
        }
        if (s != "memory" && rank == 0)
            std::cerr << "[lic] unknown streamline_storage='" << s << "', using 'memory'\n";
        if (rank == 0)
            std::cerr << "[lic] storage = memory (folded from RAM, <gid>.bin deleted, traces freed)\n";
        return false;
    }
    return config["keep_streamline_bin"] && config["keep_streamline_bin"].as<bool>();
}

// Whether GenStreamLines should write <gid>.bin to disk in the first place. Only
// meaningful for METHOD 2 (no combine) -- METHOD 1 (combine) always needs the
// on-disk file to reassemble whole streamlines, so callers must ignore this in
// combine mode. Absent key -> true (old behavior: always write, unaffected).
inline bool lic_write_bin_to_disk(const std::string& configPath)
{
    YAML::Node config = YAML::LoadFile(configPath);
    if (!config["streamline_storage"] || config["streamline_storage"].IsNull())
        return true;
    std::string s = config["streamline_storage"].as<std::string>();
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    return (s == "disk" || s == "storage");   // memory (or unrecognized) -> skip write
}

// A LIC run in progress: the accumulating image plus what the output step needs.
struct LicRun {
    LicImage    img;
    uint64_t    seed = 0;
    std::string mode = "cpu";
    long        folded = 0;        // traces folded so far (for logging)
    LicRun(int w, int h, uint64_t s, std::string m)
        : img(w, h, s), seed(s), mode(std::move(m)) {}
};

// Start a run: build the noise image. Reads lic_width/lic_height/lic_seed; auto-
// sizes from nSeeds when width/height are 0. The noise map is IDENTICAL on every
// rank (config seed, or a fresh seed drawn by rank 0 and broadcast). `mode`
// ("gpu"/"cpu") is remembered only to name the output folder.
inline LicRun lic_begin(const std::string& configPath, size_t nSeeds,
                        const std::string& mode = "cpu", MPI_Comm comm = MPI_COMM_WORLD)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    YAML::Node config = YAML::LoadFile(configPath);

    int licW = config["lic_width"]  ? config["lic_width"].as<int>()  : 0;
    int licH = config["lic_height"] ? config["lic_height"].as<int>() : 0;

    uint64_t licSeed = 0;
    if (config["lic_seed"] && !config["lic_seed"].IsNull()) {
        licSeed = config["lic_seed"].as<uint64_t>();
    } else {
        if (rank == 0) {
            std::random_device rd;
            licSeed = ((uint64_t)rd() << 32) ^ (uint64_t)rd()
                    ^ (uint64_t)std::chrono::high_resolution_clock::now()
                          .time_since_epoch().count();
        }
        MPI_Bcast(&licSeed, 1, MPI_UINT64_T, 0, comm);
    }
    if (rank == 0)
        std::cerr << "[lic] noise seed = " << licSeed
                  << (config["lic_seed"] ? " (from config)" : " (random per run)") << "\n";

    // Resolution: one pixel per seed unless the config forced a size. Seed grids
    // are N*N (square) or (2N)*N (2:1).
    if (licW <= 0 || licH <= 0) {
        const size_t n = nSeeds;
        const size_t s = size_t(std::llround(std::sqrt(double(n))));
        if (s > 0 && s * s == n) { licW = int(s); licH = int(s); }                    // N x N
        else {
            const size_t s2 = size_t(std::llround(std::sqrt(double(n) / 2.0)));
            if (s2 > 0 && 2 * s2 * s2 == n) { licW = int(2 * s2); licH = int(s2); }    // 2N x N
            else { licW = 2048; licH = 1024; }                                        // fallback
        }
        if (rank == 0)
            std::cerr << "[lic] auto resolution " << licW << "x" << licH
                      << " from " << n << " seeds\n";
    }
    return LicRun(licW, licH, licSeed, mode);
}

// Fold one batch of streamlines into the run image (used by METHOD 2, once per
// direction). A VECTOR3 is three contiguous doubles, so a vector<VECTOR3>
// reinterpret-casts straight to the flat [x,y,z,...] array lic_add_streamline wants.
inline void lic_fold(LicRun& run, const std::list<std::vector<VECTOR3>>& streamlines)
{
    for (const auto& tr : streamlines) {
        if (tr.empty()) continue;
        lic_add_streamline(run.img, reinterpret_cast<const double*>(tr.data()), tr.size());
    }
    run.folded += (long)streamlines.size();
}

// Rank-0 output step, shared by both methods: create the per-run folder
//     <outputDir>/<mode>_<direction>_<datafile>/
// write lic.png + run_info.txt (a copy of the config = the run's features) there,
// then move-or-delete the per-block <gid>.bin files GenStreamLines left in the base
// outputDir (config key keep_streamline_bin -- applied here so the core
// GenStreamLines is untouched). Assumes run.img already holds the FINAL image on
// rank 0 (after lic_reduce for METHOD 2, or after the gather-fold for METHOD 1).
inline void lic_write_outputs(LicRun& run, const std::string& configPath,
                              MPI_Comm comm = MPI_COMM_WORLD)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    YAML::Node config = YAML::LoadFile(configPath);

    const std::string outputDir = (config["outputDir"] && !config["outputDir"].IsNull())
                                  ? config["outputDir"].as<std::string>() : "lic_output";
    const bool keepBin = lic_keep_bin(configPath, rank);   // streamline_storage: memory|disk
    const int  nblocks = config["nblocks"] ? config["nblocks"].as<int>() : 0;

    if (rank == 0) {
        namespace fs = std::filesystem;
        const int licW = run.img.w, licH = run.img.h;

        const std::string runDir = outputDir + "/" + lic_run_label(run.mode, configPath);
        fs::create_directories(runDir);

        // 1) the LIC image
        const std::string pngPath = runDir + "/lic.png";
        std::vector<uint8_t> gray = lic_to_gray(run.img);
        if (lic_write_png_gray(pngPath, licW, licH, gray.data()))
            std::cerr << "[lic] wrote " << pngPath << " (" << licW << "x" << licH << ")\n";
        else
            std::cerr << "[lic] ERROR writing " << pngPath << "\n";

        // 2) run_info.txt: a short header + a verbatim copy of the config (features)
        {
            std::ofstream info(runDir + "/run_info.txt");
            info << "# LIC run\n";
            info << "mode:        " << run.mode << "\n";
            info << "resolution:  " << licW << "x" << licH << "\n";
            info << "noise_seed:  " << run.seed << "\n";
            info << "keep_bin:    " << (keepBin ? "true" : "false") << "\n";
            info << "config_file: " << configPath << "\n";
            info << "\n# ---- config (" << configPath << ") ----\n";
            std::ifstream in(configPath);
            if (in) info << in.rdbuf();
        }
        std::cerr << "[lic] wrote " << runDir << "/run_info.txt\n";

        // 3) streamline .bin: move into the run folder (keep) or delete. Assumes a
        // shared filesystem -- the .bin workflow (all ranks writing into one
        // outputDir) already requires that.
        int handled = 0;
        for (int gid = 0; gid < nblocks; ++gid) {
            const std::string bin = outputDir + "/" + std::to_string(gid) + ".bin";
            std::error_code ec;
            if (!fs::exists(bin, ec)) continue;
            if (keepBin) fs::rename(bin, runDir + "/" + std::to_string(gid) + ".bin", ec);
            else         fs::remove(bin, ec);
            if (!ec) ++handled;
        }
        std::cerr << "[lic] streamline .bin " << (keepBin ? "moved into run folder" : "deleted")
                  << " (" << handled << "/" << nblocks << " blocks)\n";
    }

    MPI_Barrier(comm);   // no rank reaches MPI_Finalize until rank 0 finished writing
}

// METHOD 2 (no combine): reduce the per-rank images onto rank 0, then write.
inline void lic_finish(LicRun& run, const std::string& configPath,
                       MPI_Comm comm = MPI_COMM_WORLD)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) std::cerr << "[lic] folded " << run.folded << " local traces\n";
    lic_reduce(run.img, comm, 0);          // combine per-rank IMAGES
    lic_write_outputs(run, configPath, comm);
}

// One streamline segment read back from a <gid>.bin file (coords flattened xyz).
struct LicSeg { int pid = 0, sid = 0; std::vector<double> xyz; };

// Read one <gid>.bin (format written by block.hpp write_trajectory):
//   repeated: [int32 pid][int32 gid][int32 sid][int32 nPts][float64 xyz]*nPts
// Appends each segment to `out`. (Little-endian host, matching the writer.)
inline void lic_read_bin_segments(const std::string& path, std::vector<LicSeg>& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) return;
    while (true) {
        int hdr[4];
        f.read(reinterpret_cast<char*>(hdr), sizeof(hdr));
        if (f.gcount() < (std::streamsize)sizeof(hdr)) break;      // clean EOF
        const int pid = hdr[0], sid = hdr[2], nPts = hdr[3];
        if (nPts < 0) break;
        LicSeg s; s.pid = pid; s.sid = sid; s.xyz.resize((size_t)nPts * 3);
        const std::streamsize want = (std::streamsize)nPts * 3 * (std::streamsize)sizeof(double);
        f.read(reinterpret_cast<char*>(s.xyz.data()), want);
        if (f.gcount() < want) break;                              // truncated tail
        out.push_back(std::move(s));
    }
}

// METHOD 1 (combine): rank 0 reads every <gid>.bin from the base outputDir,
// reassembles the segments into whole streamlines (group by pid, order by sid,
// drop the duplicated boundary point -- same as scripts/plot/read_traces.py), folds
// those whole lines into the image, then writes (no image reduce -- rank 0 has the
// whole picture). Rank 0 holds all the streamline data at once.
inline void lic_combine_finish(LicRun& run, const std::string& configPath,
                               MPI_Comm comm = MPI_COMM_WORLD)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    YAML::Node config = YAML::LoadFile(configPath);
    const std::string outputDir = (config["outputDir"] && !config["outputDir"].IsNull())
                                  ? config["outputDir"].as<std::string>() : "lic_output";
    const int nblocks = config["nblocks"] ? config["nblocks"].as<int>() : 0;

    if (rank == 0) {
        // 1) read all per-block segment files
        std::vector<LicSeg> segs;
        for (int gid = 0; gid < nblocks; ++gid)
            lic_read_bin_segments(outputDir + "/" + std::to_string(gid) + ".bin", segs);
        std::cerr << "[lic] read " << segs.size() << " segments from " << nblocks << " block files\n";

        // 2) group segment indices by pid
        std::unordered_map<int, std::vector<size_t>> byPid;
        for (size_t i = 0; i < segs.size(); ++i)
            byPid[segs[i].pid].push_back(i);

        // 3) per pid: order by sid, concatenate (drop each later segment's first
        //    point -- it repeats the previous segment's last), fold the whole line
        std::vector<double> line;
        long nLines = 0;
        for (auto& kv : byPid) {
            auto& idxs = kv.second;
            std::sort(idxs.begin(), idxs.end(),
                      [&](size_t a, size_t b){ return segs[a].sid < segs[b].sid; });
            line.clear();
            for (size_t k = 0; k < idxs.size(); ++k) {
                const std::vector<double>& xyz = segs[idxs[k]].xyz;
                const size_t startDbl = (k == 0) ? 0 : 3;   // skip 1 duplicated point
                if (xyz.size() > startDbl)
                    line.insert(line.end(), xyz.begin() + startDbl, xyz.end());
            }
            if (!line.empty()) {
                lic_add_streamline(run.img, line.data(), line.size() / 3);
                ++nLines;
            }
        }
        run.folded = nLines;
        std::cerr << "[lic] reassembled + folded " << nLines << " whole streamlines on rank 0\n";
    }

    lic_write_outputs(run, configPath, comm);   // no reduce: rank 0 has the full image
}

#endif // LIC_RENDER_HPP
