// lic_render.hpp
// -----------------------------------------------------------------------------
// The LIC method, and nothing else. Folds streamlines already generated OUTSIDE
// this folder (ParaFlow::GenStreamLines, called from ParaFlow_lic.cpp /
// ParaFlow_lic_gpu.cpp) into a Line-Integral-Convolution image and writes ONE
// grayscale PNG on rank 0. No streamline generation, no MPI init here.
//
// Include AFTER <ParaFlow.hpp> (needs VECTOR3) and after MPI is initialized
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
#include <cstring>
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

// Integration directions to run (config key `lic_direction`): a LIC line must
// run THROUGH its seed, unlike a one-way scientific streamline.
//     forward -> {+1}   backward -> {-1}   both (default) -> {+1, -1}
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

// Seed count read from the config's seed file, so image auto-sizing needs no
// access to ParaFlow's private seed list. Binary: 3 float64/seed; text: 1/line.
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

// Output folder name: "<mode>_<direction>_<datafile>" (datafile = basename of
// datafile_src, .nc stripped). Same mode/direction/data reuses the same folder.
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

// Where traces end up (config `streamline_storage`): rank 0 reassembles whole
// streamlines from <gid>.bin (GenStreamLines always writes it -- no fill-only
// path), so the file exists regardless; this only controls what happens to it
// afterward.
//   memory = DELETE <gid>.bin once folded -- only the PNG remains.
//   disk   = KEEP <gid>.bin in the run folder for later reuse (read_traces.py, ...).
// Point outputDir at tmpfs (e.g. /dev/shm/lic_run) to keep the transient write
// off real disk too.
// Back-compat: falls back to bool `keep_streamline_bin` if the key is absent.
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

// A run in progress: the accumulating image plus output-step state. img folds at
// the NATIVE (seed-matched) resolution to avoid black stripes; stretchW/H scale
// the finished image at write time (lic_write_outputs) to the output aspect ratio.
struct LicRun {
    LicImage    img;
    uint64_t    seed = 0;
    std::string mode = "cpu";
    long        folded = 0;        // traces folded so far (for logging)
    int         stretchW = 1;
    int         stretchH = 1;
    LicRun(int w, int h, uint64_t s, std::string m)
        : img(w, h, s), seed(s), mode(std::move(m)) {}
};

// Start a run: build the noise image, picking the NATIVE fold resolution:
//   1) seed_width/height (recommended) -> fold matches the seed grid's density
//      (every pixel reachable); lic_ratio_w/h then just STRETCHES the output at
//      write time (lic_write_outputs) -- a plain resize, so no black stripes.
//      e.g. seed_width/height 256 + ratio 2:1 folds 256x256, outputs 512x256.
//   2) lic_width/height (legacy) -> used directly, no stretch. Must stay <= the
//      seed grid's resolution per dimension or stripes reappear.
//   3) neither -> auto-size from nSeeds (assumes an N*N or 2N*N seed grid).
// Noise map is IDENTICAL on every rank (config seed, or rank 0 draws + broadcasts).
// `mode` ("gpu"/"cpu") is only remembered to name the output folder.
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

    const bool haveSeedRes = config["seed_width"]  && !config["seed_width"].IsNull()
                           && config["seed_height"] && !config["seed_height"].IsNull();

    int stretchW = config["lic_ratio_w"] && !config["lic_ratio_w"].IsNull()
                 ? config["lic_ratio_w"].as<int>() : 1;
    int stretchH = config["lic_ratio_h"] && !config["lic_ratio_h"].IsNull()
                 ? config["lic_ratio_h"].as<int>() : 1;
    if (stretchW < 1) stretchW = 1;
    if (stretchH < 1) stretchH = 1;

    if (haveSeedRes) {
        licW = config["seed_width"].as<int>();
        licH = config["seed_height"].as<int>();
        if (rank == 0)
            std::cerr << "[lic] fold resolution " << licW << "x" << licH
                      << " (seed grid), output stretched to "
                      << licW * stretchW << "x" << licH * stretchH << "\n";
    } else {
        stretchW = stretchH = 1;    // no seed grid given -> no stretch, legacy behavior
        if (licW <= 0 || licH <= 0) {
            // Resolution: one pixel per seed unless the config forced a size. Seed
            // grids are N*N (square) or (2N)*N (2:1).
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
    }
    LicRun run(licW, licH, licSeed, mode);
    run.stretchW = stretchW;
    run.stretchH = stretchH;
    return run;
}

// One streamline segment (coords flattened xyz), gathered directly from a
// rank's in-memory buffer (see lic_gather_segments) -- never read from a file.
struct LicSeg { int pid = 0, gid = 0, sid = 0; std::vector<double> xyz; };

// Parse the wire format Block::append_trajectory_bytes / write_trajectory both
// use: repeated [int32 pid][int32 gid][int32 sid][int32 nPts][float64 xyz]*nPts.
// Appends each segment to `out`. (Little-endian host, matching the writer.)
inline void lic_parse_segments(const char* data, size_t len, std::vector<LicSeg>& out)
{
    size_t off = 0;
    while (off + 4 * sizeof(int) <= len) {
        int hdr[4];
        std::memcpy(hdr, data + off, sizeof(hdr));
        off += sizeof(hdr);
        const int pid = hdr[0], gid = hdr[1], sid = hdr[2], nPts = hdr[3];
        if (nPts < 0) break;
        const size_t want = (size_t)nPts * 3 * sizeof(double);
        if (off + want > len) break;                              // truncated tail
        LicSeg s; s.pid = pid; s.gid = gid; s.sid = sid;
        s.xyz.resize((size_t)nPts * 3);
        std::memcpy(s.xyz.data(), data + off, want);
        off += want;
        out.push_back(std::move(s));
    }
}

// MPI-gather one direction's segments (already serialized locally by
// GenStreamLines's segBytesOut) straight to rank 0 and parse them there -- no
// rank ever writes a file. `out` is only populated on rank 0.
inline void lic_gather_segments(const std::vector<char>& localBytes,
                                std::vector<LicSeg>& out,
                                MPI_Comm comm = MPI_COMM_WORLD)
{
    int rank = 0, nranks = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    const int localSize = (int)localBytes.size();
    std::vector<int> sizes(rank == 0 ? nranks : 0);
    MPI_Gather(&localSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, comm);

    std::vector<int> displs;
    std::vector<char> allBytes;
    if (rank == 0) {
        displs.resize(nranks);
        int total = 0;
        for (int r = 0; r < nranks; ++r) { displs[r] = total; total += sizes[r]; }
        allBytes.resize(total);
    }
    MPI_Gatherv(localBytes.data(), localSize, MPI_CHAR,
                rank == 0 ? allBytes.data() : nullptr,
                rank == 0 ? sizes.data()    : nullptr,
                rank == 0 ? displs.data()   : nullptr,
                MPI_CHAR, 0, comm);

    if (rank == 0) lic_parse_segments(allBytes.data(), allBytes.size(), out);
}

// Rank-0 output step: create <outputDir>/<mode>_<direction>_<datafile>/, write
// lic.png + run_info.txt (a config copy), then -- only if streamline_storage:
// disk asked for it -- write the gathered segments out as <gid>.bin (both
// directions merged per gid, same wire format read_traces.py already expects,
// so it needs no changes). With streamline_storage: memory (default) the
// segments only ever existed in RAM (moved rank-to-rank via MPI, never a
// file), so there's nothing to clean up -- they're just dropped here.
// Assumes run.img already holds the FINAL image on rank 0.
inline void lic_write_outputs(LicRun& run, const std::string& configPath,
                              const std::vector<LicSeg>& fwdSegs,
                              const std::vector<LicSeg>& bwdSegs,
                              MPI_Comm comm = MPI_COMM_WORLD)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    YAML::Node config = YAML::LoadFile(configPath);

    const std::string outputDir = (config["outputDir"] && !config["outputDir"].IsNull())
                                  ? config["outputDir"].as<std::string>() : "lic_output";
    const bool keepBin = lic_keep_bin(configPath, rank);   // streamline_storage: memory|disk

    if (rank == 0) {
        namespace fs = std::filesystem;
        const int foldW = run.img.w, foldH = run.img.h;
        const int licW = foldW * run.stretchW, licH = foldH * run.stretchH;

        const std::string runDir = outputDir + "/" + lic_run_label(run.mode, configPath);
        fs::create_directories(runDir);

        // 1) the LIC image: fold at native resolution, then stretch to the
        // requested output size (lic_ratio_w/h) -- a plain resize, no stripes.
        const std::string pngPath = runDir + "/lic.png";
        std::vector<uint8_t> gray = lic_to_gray(run.img);
        if (run.stretchW != 1 || run.stretchH != 1)
            gray = lic_stretch_gray(gray, foldW, foldH, licW, licH);
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

        // 3) streamline .bin: only written if streamline_storage: disk asked for
        // it. Segments never touched a filesystem before this point (they moved
        // rank-to-rank over MPI -- see lic_gather_segments), so `memory` mode
        // means literally zero streamline bytes ever hit storage.
        if (keepBin) {
            std::unordered_map<int, std::vector<const LicSeg*>> byGid;
            for (const auto& s : fwdSegs) byGid[s.gid].push_back(&s);
            for (const auto& s : bwdSegs) byGid[s.gid].push_back(&s);
            for (auto& kv : byGid) {
                std::ofstream out(runDir + "/" + std::to_string(kv.first) + ".bin", std::ios::binary);
                for (const LicSeg* s : kv.second) {
                    const int nPts = (int)(s->xyz.size() / 3);
                    out.write((const char*)&s->pid, sizeof(int));
                    out.write((const char*)&s->gid, sizeof(int));
                    out.write((const char*)&s->sid, sizeof(int));
                    out.write((const char*)&nPts,   sizeof(int));
                    out.write((const char*)s->xyz.data(), (std::streamsize)(s->xyz.size() * sizeof(double)));
                }
            }
            std::cerr << "[lic] streamline .bin written to run folder (" << byGid.size() << " blocks)\n";
        } else {
            std::cerr << "[lic] streamline_storage: memory -- .bin never touched disk\n";
        }
    }

    MPI_Barrier(comm);   // no rank reaches MPI_Finalize until rank 0 finished writing
}

// Concatenate one direction's pid-ordered segments (sorted by sid ascending)
// into a single polyline: segment 0 starts at the seed, each later segment
// drops its first point (it repeats the previous segment's last -- the
// block-boundary crossing point). Same rule read_traces.py uses.
inline std::vector<double> lic_build_line(const std::vector<LicSeg>& segs,
                                          const std::vector<size_t>& idxs)
{
    std::vector<double> line;
    for (size_t k = 0; k < idxs.size(); ++k) {
        const std::vector<double>& xyz = segs[idxs[k]].xyz;
        const size_t startDbl = (k == 0) ? 0 : 3;   // skip 1 duplicated point
        if (xyz.size() > startDbl)
            line.insert(line.end(), xyz.begin() + startDbl, xyz.end());
    }
    return line;
}

// Reassemble whole streamlines from already-gathered segments (see
// lic_gather_segments -- rank 0 only, everyone else's vectors are empty),
// fold them, then write outputs (no reduce -- rank 0 already has it all).
inline void lic_finish(LicRun& run, const std::string& configPath,
                       const std::vector<LicSeg>& fwdSegs,
                       const std::vector<LicSeg>& bwdSegs,
                       MPI_Comm comm = MPI_COMM_WORLD)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        std::cerr << "[lic] gathered " << fwdSegs.size() << " forward + " << bwdSegs.size()
                  << " backward segments (in RAM, no disk)\n";

        // Group each direction's segment indices by pid, ordered by sid. Each
        // direction's sid restarts at 0 from the shared seed, so a fwd segment
        // and a bwd segment can share the same sid -- fwdSegs/bwdSegs must stay
        // separate groupings, never merged into one sorted list, or reassembly
        // would interleave two unrelated segments.
        auto groupByPid = [](const std::vector<LicSeg>& segs) {
            std::unordered_map<int, std::vector<size_t>> byPid;
            for (size_t i = 0; i < segs.size(); ++i) byPid[segs[i].pid].push_back(i);
            for (auto& kv : byPid)
                std::sort(kv.second.begin(), kv.second.end(),
                          [&](size_t a, size_t b){ return segs[a].sid < segs[b].sid; });
            return byPid;
        };
        std::unordered_map<int, std::vector<size_t>> fwdByPid = groupByPid(fwdSegs);
        std::unordered_map<int, std::vector<size_t>> bwdByPid = groupByPid(bwdSegs);

        std::vector<int> pids;
        for (auto& kv : fwdByPid) pids.push_back(kv.first);
        for (auto& kv : bwdByPid)
            if (!fwdByPid.count(kv.first)) pids.push_back(kv.first);

        // Per pid: build each direction's line separately (both start at the
        // shared seed point), then stitch reversed-backward (far -> seed) +
        // forward (seed -> far) so the LIC line runs through the seed, dropping
        // the duplicated seed point at the join. Single-direction pids (only
        // fwd or only bwd present) fold that one line as-is.
        long nLines = 0;
        for (int pid : pids) {
            std::vector<double> fwdLine, bwdLine;
            auto fit = fwdByPid.find(pid);
            if (fit != fwdByPid.end()) fwdLine = lic_build_line(fwdSegs, fit->second);
            auto bit = bwdByPid.find(pid);
            if (bit != bwdByPid.end()) bwdLine = lic_build_line(bwdSegs, bit->second);

            std::vector<double> line;
            if (!fwdLine.empty() && !bwdLine.empty()) {
                const size_t n = bwdLine.size() / 3;
                for (size_t i = 0; i < n; ++i)
                    line.insert(line.end(), bwdLine.end() - 3 * (i + 1), bwdLine.end() - 3 * i);
                line.insert(line.end(), fwdLine.begin() + 3, fwdLine.end());   // skip dup seed point
            } else if (!fwdLine.empty()) {
                line = std::move(fwdLine);
            } else {
                line = std::move(bwdLine);
            }

            if (!line.empty()) {
                lic_add_streamline(run.img, line.data(), line.size() / 3);
                ++nLines;
            }
        }
        run.folded = nLines;
        std::cerr << "[lic] reassembled + folded " << nLines << " whole streamlines on rank 0\n";
    }

    lic_write_outputs(run, configPath, fwdSegs, bwdSegs, comm);   // no reduce: rank 0 has it all
}

#endif // LIC_RENDER_HPP
