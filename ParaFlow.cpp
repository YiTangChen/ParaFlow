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

void ParaFlow::ReadSeedFile(const char* seedFile, std::list<VECTOR3>& seeds) 
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
    }
    myFile.close();
}

ParaFlow::ParaFlow(const char* configFile) 
{
    YAML::Node config = YAML::LoadFile(configFile);

    // Read configs from YAML
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
        if(!config["vecFile"])
            throw "Error: no vector file(s).\n";
    }
    catch(const char* msg) {
        std::cout << msg;
    }

    if(config["vecFile"].IsScalar()) {
        this->vectorFile = config["vecFile"].as<string>();
        this->ReadRegularGridXYZ(this->vectorFile.c_str(), this->domain_x, this->domain_y, this->domain_z);
    }
    else
        throw "We don't have unstable vector field now.\n";
    
    if(config["seedFile"])
        this->ReadSeedFile(config["seedFile"].as<string>().c_str(), this->seeds);
}

ParaFlow::~ParaFlow() {

}

void ParaFlow::deq_incoming_iexchange(Block* b, const diy::Master::ProxyWithLink& cp)
{
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds>*>(cp.link());
    for (size_t i = 0; i < l->size(); ++i)
    {
        int nbr_gid = l->target(i).gid;
        while (cp.incoming(nbr_gid))
        {
            PtInfo incoming_pt;
            cp.dequeue(nbr_gid, incoming_pt);
            b->addStartPts(incoming_pt);
        }
    }
}

void ParaFlow::trace_block(Block*                               b,
                            const diy::Master::ProxyWithLink&   cp,
                            const Decomposer&                   decomposer)
{
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds>*>(cp.link());

    do
    {
        string info = "trace block " + to_string(cp.gid()) + "\n";
        b->printInfo(info);
        this->deq_incoming_iexchange(b, cp);
        b->GenStreamLineByOSUFlow(this->max_steps);
        info = "There are " + to_string(b->sl_list.size()) + " particles in this block\n";
        b->printInfo(info);
        // deal with current seeds
        int seedCnt = 0;
        for(auto sl = b->sl_list.begin(); sl != b->sl_list.end(); sl++) {
            info = "current seed cnt = " + to_string(seedCnt) + "\n";
            b->printInfo(info);
            vtListSeedTrace * trace = *sl;
            bool finished = false;
            Segment currseg;
            currseg.gid = b->startPts[seedCnt].gid;
            currseg.pid = b->startPts[seedCnt].pid;
            currseg.sid = b->startPts[seedCnt].sid;
            currseg.nsteps = b->startPts[seedCnt].nsteps;
            info = "before go through all positions in trace\n";
            b->printInfo(info);
            for(auto pIt = trace->begin(); pIt != trace->end() && currseg.nsteps < this->max_steps; pIt++) {
                VECTOR3 tmp;
                for (size_t i = 0; i < 3; i++)
                    tmp[i] = (**pIt)[i];
                info = "(" + to_string(tmp[0]) + ", " + to_string(tmp[1]) + ", " + to_string(tmp[2]) + ")\n";
                b->printInfo(info);
                currseg.coords.push_back(tmp);
                if (pIt != trace->begin())
                    currseg.nsteps++;
            }
            info = "after go through all positions in trace\n";
            b->printInfo(info);
            b->segs.push_back(currseg);
            if(currseg.nsteps >= this->max_steps) {
                info = "current trace steps > max steps\n";
                b->printInfo(info);
                finished = true;
            }
                
            if(!inDomainVec(currseg.coords.back(), decomposer.domain)) {
                info = "the particle goes outside the domain\n";
                b->printInfo(info);
                finished = true;
            }

            if(!finished) {
                info = "The particle will go to other block\n";
                b->printInfo(info);
                Bounds neigh_bounds {0}; // neighbor block bounds
                PtInfo endPt;
                int next_block = -1;
                for(int i = 0; i < 3; i++)
                    endPt.coord[i] = currseg.coords.back()[i];
                endPt.gid = currseg.gid;
                endPt.pid = currseg.pid;
                endPt.sid = currseg.sid;
                endPt.nsteps = currseg.nsteps;

                info = "Find the block...\n";
                b->printInfo(info);
                for(int n = 0; n < l->size(); n++) {
                    neigh_bounds = l->core(n);
                    bool inBlock = inDomainPt(endPt, neigh_bounds);
                    if(inBlock) {
                        next_block = n; // directly use n for target?
                        break;
                    }
                }
                info = "The block index is " + to_string(next_block) + "\n";
                b->printInfo(info);

                if(next_block >= 0) {
                    diy::BlockID bid = l->target(next_block); // in case of multiple dests, send to first dest only
                    info = "The block is " + to_string(bid.gid) + "\n";
                    b->printInfo(info);
                    cp.enqueue(bid, endPt);
                    info = "Enqueue!\n";
                    b->printInfo(info);
                }
                else {
                    info = "Cannot find the block.\n";
                    b->printInfo(info);
                }
            }
            seedCnt++;
        }
        b->currentSeeds.clear();
        b->startPts.clear();
        b->sl_list.clear();
    } while (cp.fill_incoming());
}

bool ParaFlow::trace_block_iexchange(Block*                              b,
                                    const diy::Master::ProxyWithLink&    cp,
                                    const Decomposer&                    decomposer)
{
    this->trace_block(b, cp, decomposer);
    return true;
}

void ParaFlow::GenStreamLines(int argc, char* argv[], std::list<std::list<VECTOR3>>& sl_list)
{
    // Init DIY
    diy::mpi::environment     env(argc, argv);         // diy equivalent of MPI_Init
    diy::mpi::communicator    world;                   // diy equivalent of MPI communicator
    diy::ContiguousAssigner   assigner(this->size, this->nblocks);

    Bounds domain(3);                                   // global data size
    domain.min[0] = domain.min[1] = domain.min[2] = 0;
    domain.max[0] = domain_x-1;
    domain.max[1] = domain_y-1;
    domain.max[2] = domain_z-1;
    
    int rank = world.rank();                           // MPI rank of this process
    diy::Master master(world,
                       1,                              // one thread
                       -1,                             // all blocks in memory
                       &Block::create,
                       &Block::destroy);

    // share_face is an n-dim (size 3 in this example) vector of bools
    // indicating whether faces are shared in each dimension
    // uninitialized values default to false
    diy::RegularDecomposer<Bounds>::BoolVector          share_face;
    share_face.push_back(true); share_face.push_back(true); share_face.push_back(true);

    // wrap is an n-dim (size 3 in this example) vector of bools
    // indicating whether boundary conditions are periodic in each dimension
    // uninitialized values default to false
    diy::RegularDecomposer<Bounds>::BoolVector          wrap;

    // ghosts is an n-dim (size 3 in this example) vector of ints
    // indicating number of ghost cells per side in each dimension
    // uninitialized values default to 0
    diy::RegularDecomposer<Bounds>::CoordinateVector    ghosts;
    ghosts.push_back(1); ghosts.push_back(1); ghosts.push_back(1);

    // either create the regular decomposer and call its decompose function
    // (having the decomposer available is useful for its other member functions
    diy::RegularDecomposer<Bounds> decomposer(3,
                                              domain,
                                              nblocks,
                                              share_face,
                                              wrap,
                                              ghosts);

    //////////////// Record file ///////////////////
    std::ofstream outfile;
    outfile.open(std::to_string(rank)+".txt");
    outfile << "";
    outfile.close();
    ////////////////////////////////////////////////

    decomposer.decompose(rank, assigner,
                         [&](int gid,                          // block global id
                             const Bounds& core,               // block bounds without any ghost added
                             const Bounds& bounds,             // block bounds including ghost region
                             const Bounds& domain,             // global data bounds
                             const RGLink& link)               // neighborhood
                         {
                             Block *b = new Block;
                             RGLink *l = new RGLink(link);
                             // Init block
                             b->set_rank(rank);
                             b->set_data(gid, this->vectorFile.c_str(), this->seeds, core, bounds);
                             master.add(gid, b, l); // add block to the master (mandatory)
                         });
    int ncalls = 0;
    master.iexchange([&](Block* b, const diy::Master::ProxyWithLink& icp) -> bool
    {
        ncalls++;
        bool val = this->trace_block_iexchange(b, icp, decomposer);
        return val;
    });
    world.barrier();
    this->write_traces(master, decomposer, assigner, sl_list);
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

void ParaFlow::write_traces(diy::Master&                    master,
                            const Decomposer&               decomposer,
                            const diy::Assigner&            assigner,
                            std::list<std::list<VECTOR3>>&  sl_list) {
    // merge-reduce traces to one block
    int k = 2;                               // the radix of the k-ary reduction tree
    diy::RegularMergePartners  partners(decomposer, k);
    diy::reduce(master, assigner, partners, &merge_traces);

    if (master.communicator().rank() == 0)
    {
        std::string filename = "traces.txt";
        ((Block*)master.block(0))->write_trajectory(filename, sl_list);
    }
}