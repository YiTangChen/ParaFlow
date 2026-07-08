#ifndef _PARAFLOW_HPP
#define _PARAFLOW_HPP

#include <iostream>
#include <fstream>
#include <functional>
#include <filesystem>

#include <diy/serialization.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>

#include <OSUFlow.h>
#include <yaml-cpp/yaml.h>
#include "block.hpp"
// #include <api/MOPS.h>

typedef     diy::RegularGridLink        RGLink;
typedef     diy::DiscreteBounds         Bounds;
typedef     diy::RegularDecomposer<Bounds> Decomposer;
enum class SeedReadMode {Binary, Text, NONDEF};

static std::optional<std::string> read_opt_str(const YAML::Node& node, const char* key) {
    if (!node[key] || node[key].IsNull())
        return std::nullopt;
    return node[key].as<std::string>();
}

class ParaFlow
{
public:
    ParaFlow(int argc, char* argv[], const char* configFile);
    ~ParaFlow();

    void GenStreamLines(std::list<std::vector<VECTOR3>>& sl_list);
    void GenPathLines(std::list<std::vector<VECTOR3>>& pl_list);
    void CheckPointsInBlock(std::function<void(int gid, const std::vector<VECTOR3>&)> onBlock = nullptr);
    void setSeeds(std::vector<VECTOR3> s) { seeds = std::move(s); }
    void allreduceSum(std::vector<double>& data)
    {
        std::vector<double> result(data.size());
        diy::mpi::all_reduce(world, data, result, std::plus<double>());
        data = std::move(result);
    }
private:
    void ReadRegularGridXYZ(const char* regularField, int &domain_x, int &domain_y, int &domain_z);
    void ReadSeedFile(const char* seedFile, std::vector<VECTOR3>& seeds);
    void ReadSeedBinFile(const char* seedFile, std::vector<VECTOR3>& seeds);
    void deq_incoming_iexchange(Block* b, const diy::Master::ProxyWithLink& cp);
    void trace_block(Block*                             b,
                    const diy::Master::ProxyWithLink&   cp);
    void trace_block_gpu(Block*                         b,
                    const diy::Master::ProxyWithLink&   cp);
    bool trace_block_iexchange(Block*                                b,
                                const diy::Master::ProxyWithLink&    cp);
    void trace_block_pathline(Block*                            b,
                              const diy::Master::ProxyWithLink& cp);
    void trace_block_pathline_gpu(Block*                        b,
                                  const diy::Master::ProxyWithLink& cp);
    bool trace_block_pathline_iexchange(Block*                               b,
                                        const diy::Master::ProxyWithLink&    cp);

    void ReadAreaIndices(const char* graph_partition_file);

    // DIY / MPI — env must be declared first so it is destroyed last
    diy::mpi::environment     env;
    diy::mpi::communicator    world;

    int                       size;                 // total number of MPI processes
    int                       nblocks;              // total number of blocks in global domain
    int                       max_steps;
    int                       interval;             // output every N-th point (1 = all points)
    int                       domain_x, domain_y, domain_z;
    bool                      isTimeVarying;
    bool                      enable_timing = false; // controlled by YAML enable_timing key
    bool                      useGPU = false;        // YAML useGPU, resolved against runtime availability
    SeedReadMode              seedReadMode = SeedReadMode::NONDEF;
    std::vector<std::string>  vectorFiles;
    std::string               meshFile;
    std::string               seedFile;
    std::string               outputDir;
    std::vector<int>          areaIndices;
    std::vector<VECTOR3>      seeds;
    TimeVaryingDataConfig     cfg;
};

#endif
