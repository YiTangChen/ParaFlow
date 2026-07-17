// lic_core.hpp
// -----------------------------------------------------------------------------
// Shared streamline-LIC engine: equirectangular projection, white-noise map,
// per-streamline accumulation, an 8-bit grayscale PNG writer, regular sphere
// seeding, and (optionally) an MPI image reduce.
//
// Design goals
//   * Pure: depends only on the C++ standard library. Streamline points are
//     passed as a FLAT [x0,y0,z0, x1,y1,z1, ...] double array, so this header
//     needs no VECTOR3 and is usable by both the standalone tool and the
//     ParaFlow driver. (A std::vector<VECTOR3> is layout-compatible — just
//     reinterpret_cast its .data() to const double*.)
//   * MPI is opt-in: #define LIC_WITH_MPI before including to get lic_reduce().
//
// Pipeline (parallel, no trajectories stored):
//   LicImage img(w, h, seed);                       // 1. noise map (same seed on every rank)
//   for each finished trace: lic_add_streamline(img, xyz, n);  // 2. fold in, then free trace
//   lic_reduce(img, comm, 0);                       // 3. combine per-rank images
//   if (rank==0) lic_write_png_gray(path, w, h, lic_to_gray(img).data());  // 4. only output
// -----------------------------------------------------------------------------
#ifndef LIC_CORE_HPP
#define LIC_CORE_HPP

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>

#ifndef LIC_PI
#define LIC_PI 3.14159265358979323846
#endif

// Equirectangular projection: a 3-D point on the Earth sphere -> pixel index.
// lon = atan2(y,x), lat = asin(z/r) using each point's OWN radius r so the
// latitude is depth-independent. Returns row*w + col.
inline int64_t lic_to_pixel_index(double x, double y, double z, int w, int h)
{
    double r   = std::sqrt(x * x + y * y + z * z);
    double lon = std::atan2(y, x);
    double zr  = (r != 0.0) ? z / r : 0.0;
    zr = zr > 1.0 ? 1.0 : (zr < -1.0 ? -1.0 : zr);
    double lat = std::asin(zr);
    double u = (lon + LIC_PI) / (2.0 * LIC_PI);
    double v = (LIC_PI / 2.0 - lat) / LIC_PI;
    int col = int(u * w); col = col < 0 ? 0 : (col >= w ? w - 1 : col);
    int row = int(v * h); row = row < 0 ? 0 : (row >= h ? h - 1 : row);
    return int64_t(row) * w + col;
}

// Minimal grayscale PNG writer: signature + IHDR + IDAT(zlib "stored") + IEND.
// Needs no external libraries.
namespace lic_detail {

inline uint32_t crc32_buf(const uint8_t* buf, size_t len)
{
    static uint32_t table[256];
    static bool init = false;
    if (!init) {
        for (uint32_t n = 0; n < 256; ++n) {
            uint32_t c = n;
            for (int k = 0; k < 8; ++k)
                c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            table[n] = c;
        }
        init = true;
    }
    uint32_t crc = 0xFFFFFFFFu;
    for (size_t i = 0; i < len; ++i)
        crc = table[(crc ^ buf[i]) & 0xFF] ^ (crc >> 8);
    return crc ^ 0xFFFFFFFFu;
}

inline uint32_t adler32(const uint8_t* data, size_t len)
{
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < len; ++i) { a = (a + data[i]) % 65521; b = (b + a) % 65521; }
    return (b << 16) | a;
}

inline void put_be32(std::vector<uint8_t>& out, uint32_t v)
{
    out.push_back(uint8_t(v >> 24)); out.push_back(uint8_t(v >> 16));
    out.push_back(uint8_t(v >> 8));  out.push_back(uint8_t(v));
}

// Wrap raw bytes in a zlib stream using uncompressed "stored" deflate blocks.
inline std::vector<uint8_t> zlib_store(const std::vector<uint8_t>& raw)
{
    std::vector<uint8_t> out;
    out.push_back(0x78); out.push_back(0x01);          // zlib header
    size_t pos = 0, n = raw.size();
    do {
        size_t block = std::min<size_t>(65535, n - pos);
        bool last = (pos + block >= n);
        out.push_back(last ? 1 : 0);                   // BFINAL, BTYPE=00 (stored)
        out.push_back(uint8_t(block));      out.push_back(uint8_t(block >> 8));     // LEN
        uint16_t nlen = uint16_t(~block);
        out.push_back(uint8_t(nlen));       out.push_back(uint8_t(nlen >> 8));      // NLEN
        out.insert(out.end(), raw.begin() + pos, raw.begin() + pos + block);
        pos += block;
    } while (pos < n);
    put_be32(out, adler32(raw.data(), raw.size()));
    return out;
}

inline void write_chunk(std::ofstream& f, const char type[4], const std::vector<uint8_t>& data)
{
    std::vector<uint8_t> len; put_be32(len, uint32_t(data.size()));
    f.write(reinterpret_cast<char*>(len.data()), 4);
    std::vector<uint8_t> crcbuf(type, type + 4);
    crcbuf.insert(crcbuf.end(), data.begin(), data.end());
    uint32_t crc = crc32_buf(crcbuf.data(), crcbuf.size());
    f.write(type, 4);
    if (!data.empty()) f.write(reinterpret_cast<const char*>(data.data()), data.size());
    std::vector<uint8_t> c; put_be32(c, crc);
    f.write(reinterpret_cast<char*>(c.data()), 4);
}

} // namespace lic_detail

inline bool lic_write_png_gray(const std::string& path, int w, int h, const uint8_t* pix)
{
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    const uint8_t sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    f.write(reinterpret_cast<const char*>(sig), 8);

    std::vector<uint8_t> ihdr;
    lic_detail::put_be32(ihdr, uint32_t(w)); lic_detail::put_be32(ihdr, uint32_t(h));
    ihdr.push_back(8);   // bit depth
    ihdr.push_back(0);   // color type 0 = grayscale
    ihdr.push_back(0); ihdr.push_back(0); ihdr.push_back(0); // compression, filter, interlace
    lic_detail::write_chunk(f, "IHDR", ihdr);

    std::vector<uint8_t> raw;
    raw.reserve(size_t(h) * (size_t(w) + 1));
    for (int y = 0; y < h; ++y) {                      // each scanline: filter byte 0 + pixels
        raw.push_back(0);
        raw.insert(raw.end(), pix + size_t(y) * w, pix + size_t(y) * w + w);
    }
    lic_detail::write_chunk(f, "IDAT", lic_detail::zlib_store(raw));
    lic_detail::write_chunk(f, "IEND", {});
    return bool(f);
}

// LIC image state: a white-noise map plus the two accumulators. The noise map
// MUST be identical on every rank, so it is seeded deterministically — pass the
// SAME (broadcast) seed on all ranks.
struct LicImage {
    int w, h;
    std::vector<double> noise;      // w*h, identical on every rank
    std::vector<double> noise_sum;  // w*h, accumulator (partial sums until reduce)
    std::vector<double> count;      // w*h, accumulator

    LicImage(int w_, int h_, uint64_t seed)
        : w(w_), h(h_),
          noise(size_t(w_) * size_t(h_)),
          noise_sum(size_t(w_) * size_t(h_), 0.0),
          count(size_t(w_) * size_t(h_), 0.0)
    {
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (auto& n : noise) n = dist(rng);
    }
};

// Fold ONE streamline/segment into the image, then the caller frees the trace.
// pts3 is a flat [x,y,z, x,y,z, ...] array of length npts*3.
// Stores PARTIAL sums; final pixels are computed in lic_to_gray AFTER reduce.
inline void lic_add_streamline(LicImage& img, const double* pts3, size_t npts)
{
    if (npts == 0) return;

    std::vector<int64_t> idx;                   // pixels this line covers
    idx.reserve(npts);
    for (size_t i = 0; i < npts; ++i)
        idx.push_back(lic_to_pixel_index(pts3[i*3], pts3[i*3+1], pts3[i*3+2], img.w, img.h));
    std::sort(idx.begin(), idx.end());
    idx.erase(std::unique(idx.begin(), idx.end()), idx.end());   // distinct pixels

    double sum = 0.0;
    for (int64_t q : idx) sum += img.noise[q];
    double mean = sum / double(idx.size());                      // this line's gray value
    for (int64_t q : idx) { img.noise_sum[q] += mean; img.count[q] += 1.0; }
}

// Finalize: produce an 8-bit grayscale image (call on the root rank AFTER the
// reduce). Pixels untouched by any streamline stay black (0).
inline std::vector<uint8_t> lic_to_gray(const LicImage& img)
{
    const size_t N = size_t(img.w) * size_t(img.h);
    std::vector<uint8_t> g(N, 0);
    for (size_t i = 0; i < N; ++i) {
        if (img.count[i] > 0.0) {
            double v = img.noise_sum[i] / img.count[i];
            v = v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
            g[i] = uint8_t(std::lround(v * 255.0));
        }
    }
    return g;
}

// Nearest-neighbor stretch of an 8-bit grayscale image from srcW x srcH to
// dstW x dstH. Used to change the OUTPUT picture's aspect ratio (e.g. widen a
// square LIC into a 2:1 world map) WITHOUT touching the resolution streamlines
// were folded at -- that fold resolution must stay matched to the seed grid's
// own density, or pixels between seeds go untouched and render black (see
// lic_add_streamline / lic_to_gray above). Stretching the finished image is a
// plain resize, so it never introduces that artifact.
inline std::vector<uint8_t> lic_stretch_gray(const std::vector<uint8_t>& src,
                                             int srcW, int srcH, int dstW, int dstH)
{
    std::vector<uint8_t> out(size_t(dstW) * size_t(dstH));
    for (int y = 0; y < dstH; ++y) {
        int sy = int(size_t(y) * srcH / dstH); sy = sy >= srcH ? srcH - 1 : sy;
        const uint8_t* srow = src.data() + size_t(sy) * srcW;
        uint8_t* drow = out.data() + size_t(y) * dstW;
        for (int x = 0; x < dstW; ++x) {
            int sx = int(size_t(x) * srcW / dstW); sx = sx >= srcW ? srcW - 1 : sx;
            drow[x] = srow[sx];
        }
    }
    return out;
}

// Regular lat/lon seeding on a sphere of the given radius (for "seeding: regular").
// Returns a flat [x,y,z, ...] array of lon_count*lat_count seed points. Kept here
// (VECTOR3-free) so the caller wraps these into its own seed type.
inline std::vector<double> lic_regular_sphere_seeds(
    int lon_count, int lat_count, double radius,
    double lon_min = -LIC_PI, double lon_max = LIC_PI,
    double lat_min = -LIC_PI / 2.0, double lat_max = LIC_PI / 2.0)
{
    std::vector<double> seeds;
    if (lon_count < 1 || lat_count < 1) return seeds;
    seeds.reserve(size_t(lon_count) * size_t(lat_count) * 3);
    const double dlon = (lon_count > 1) ? (lon_max - lon_min) / (lon_count - 1) : 0.0;
    const double dlat = (lat_count > 1) ? (lat_max - lat_min) / (lat_count - 1) : 0.0;
    for (int j = 0; j < lat_count; ++j) {
        double lat = lat_min + j * dlat;
        double cl = std::cos(lat), sl = std::sin(lat);
        for (int i = 0; i < lon_count; ++i) {
            double lon = lon_min + i * dlon;
            seeds.push_back(radius * cl * std::cos(lon));   // x
            seeds.push_back(radius * cl * std::sin(lon));   // y
            seeds.push_back(radius * sl);                   // z
        }
    }
    return seeds;
}

// Optional MPI reduce of the accumulators onto `root`. Opt-in to keep this
// header MPI-free for the standalone tool: #define LIC_WITH_MPI before include.
#ifdef LIC_WITH_MPI
#include <mpi.h>
inline void lic_reduce(LicImage& img, MPI_Comm comm, int root)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    const int N = int(size_t(img.w) * size_t(img.h));   // assumes w*h < 2^31
    if (rank == root) {
        MPI_Reduce(MPI_IN_PLACE, img.noise_sum.data(), N, MPI_DOUBLE, MPI_SUM, root, comm);
        MPI_Reduce(MPI_IN_PLACE, img.count.data(),     N, MPI_DOUBLE, MPI_SUM, root, comm);
    } else {
        MPI_Reduce(img.noise_sum.data(), nullptr, N, MPI_DOUBLE, MPI_SUM, root, comm);
        MPI_Reduce(img.count.data(),     nullptr, N, MPI_DOUBLE, MPI_SUM, root, comm);
    }
}
#endif // LIC_WITH_MPI

#endif // LIC_CORE_HPP
