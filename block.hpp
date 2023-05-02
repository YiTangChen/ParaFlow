#ifndef _BLOCK_HPP
#define _BLOCK_HPP

#include <iostream>
#include <fstream>
#include <OSUFlow.h>
#include <diy/master.hpp>
#include "utils.hpp"

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
    std::vector<Segment> segs;
    
    // following is mandatory
    Block()                             {}
    void            set_rank(int r)     { rank = r; }
    int             get_rank()          { return rank; }
    static void*    create()            { return new Block; }
    static void     destroy(void* b)    { delete static_cast<Block*>(b); }

    void set_data(int gid, const char* fname, std::list<VECTOR3>& allseeds, const Bounds& core, const Bounds &bounds) {
        osuflow = new OSUFlow();
        
        VECTOR3 boundaryMin, boundaryMax;
        for(int i = 0; i < 3; i++) {
            boundaryMin[i] = float(bounds.min[i]);
            boundaryMax[i] = float(bounds.max[i]);
        }
        osuflow->LoadData(fname, true, boundaryMin, boundaryMax, false);
        
        int seedCnt = 0;
        for (std::list<VECTOR3>::iterator it = allseeds.begin(); it != allseeds.end(); ++it) {
            if(inDomainVec((*it), core)) {
                VECTOR3 tmp;
                for(int i = 0; i < 3; i++)
                    tmp[i] = (*it)[i];
                currentSeeds.push_back(tmp);
                PtInfo ptInDomain;
                for(int i = 0; i < 3; i++)
                    ptInDomain.coord[i] = (*it)[i];
                ptInDomain.pid = seedCnt;
                ptInDomain.gid = gid;
                ptInDomain.sid = 0;
                ptInDomain.nsteps = 0;
                startPts.push_back(ptInDomain);
            }
            seedCnt++;
        }
        this->setStartPtsDone();
    }

    void addStartPts(PtInfo& pt) {
        this->startPts.push_back(pt);
        VECTOR3 tmp;
        for(int i = 0; i < 3; i++)
            tmp[i] = pt.coord[i];
        this->currentSeeds.push_back(tmp);
    }

    void setStartPtsDone() {
        VECTOR3* seedsPointer = &(this->currentSeeds[0]);
        osuflow->SetSeedPoints(seedsPointer, this->currentSeeds.size());
    }

    void GenStreamLineByOSUFlow(int maxSteps) {
        this->sl_list.clear();
        osuflow->SetIntegrationParams(1, 5);
	    osuflow->GenStreamLines(this->sl_list, FORWARD_DIR, maxSteps, 0);
    }

    void printInfo(string &str) {
        std::ofstream outfile;
        outfile.open(std::to_string(this->get_rank())+".txt", std::ios_base::app);
        outfile << str;
        outfile.close();
    }

    void write_trajectory(std::string filename, std::list<std::list<VECTOR3>>& sl_list) {
        // collect segments by pid
        std::map<int, std::vector<std::pair<int, int>>> pididx;
        for(int i = 0; i < this->segs.size(); i++) {
            int pid = this->segs[i].pid;
            if(pididx.find(pid) == pididx.end()) {
                pididx[pid] = vector<std::pair<int, int>>();
            }
            pididx[pid].push_back(std::pair<int ,int>(this->segs[i].sid, i));
            // std::cout <<  "find trajectory pid = " << to_string(pid) << '\n';
        }
        // sort segments by sid
        for(auto it = pididx.begin(); it != pididx.end(); it++) {
            int pid = it->first;
            std::sort(pididx[pid].begin(), pididx[pid].end(), [](std::pair<int, int> &lhs, std::pair<int, int> &rhs)
                      { return lhs.first < rhs.first; });
        }

        std::ofstream outfile;
        outfile.open(filename);
        for(auto it = pididx.begin(); it != pididx.end(); it++) {
            std::list<VECTOR3> streamline;
            int pid = it->first;
            outfile << "Trajectory\n";
            for(int i = 0; i < pididx[pid].size(); i++) {
                int segidx = pididx[pid][i].second;
                if(i == 0) {
                    outfile << this->segs[segidx].coords[0][0] << "," << this->segs[segidx].coords[0][1] << "," << this->segs[segidx].coords[0][2] << "\n";
                    streamline.push_back(VECTOR3(this->segs[segidx].coords[0]));
                }
                for(int j = 1; j < this->segs[segidx].coords.size(); j++) {
                    outfile << this->segs[segidx].coords[j][0] << "," << this->segs[segidx].coords[j][1] << "," << this->segs[segidx].coords[j][2] << "\n";
                    streamline.push_back(VECTOR3(this->segs[segidx].coords[j]));
                }
            }
        }
        outfile.close();
    }
};

#endif