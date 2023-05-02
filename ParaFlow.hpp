#ifndef _PARAFLOW_HPP
#define _PARAFLOW_HPP

#include <iostream>
#include <fstream>

#include <diy/serialization.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>

#include <OSUFlow.h>
#include <yaml-cpp/yaml.h>
#include "block.hpp"

typedef     diy::RegularGridLink        RGLink;
typedef     diy::DiscreteBounds         Bounds;
typedef     diy::RegularDecomposer<Bounds> Decomposer;

class ParaFlow
{
public:
    ParaFlow(const char* configFile);
    ~ParaFlow();

    void GenStreamLines(int argc, char* argv[], std::list<std::list<VECTOR3>>& sl_list);
private:
    void ReadRegularGridXYZ(const char* regularField, int &domain_x, int &domain_y, int &domain_z);
    void ReadSeedFile(const char* seedFile, std::list<VECTOR3>& seeds);
    void deq_incoming_iexchange(Block* b, const diy::Master::ProxyWithLink& cp);
    void trace_block(Block*                             b,
                    const diy::Master::ProxyWithLink&   cp,
                    const Decomposer&                   decomposer);
    bool trace_block_iexchange(Block*                                b,
                                const diy::Master::ProxyWithLink&    cp,
                                const Decomposer&                    decomposer);
    void write_traces(diy::Master&                  master,
                    const Decomposer&               decomposer,
                    const diy::Assigner&            assigner,
                    std::list<std::list<VECTOR3>>&  sl_list);

    int                       size;                 // total number of MPI processes
    int                       nblocks;              // total number of blocks in global domain
    int                       max_steps;
    int                       domain_x, domain_y, domain_z;
    std::string               vectorFile;
    std::list<VECTOR3>        seeds;
};

#endif